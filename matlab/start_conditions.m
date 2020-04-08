%% simulation parameters (cf. table in assignement)
clear all; close all; 
T_cel = 25;
nu_u = 20.8 / 100;
nu_v = 0.04 / 100;


% model parameters

% Radial and axial diffusivity of oxygen in pear tissue:
sigma_ur = 2.8e-10;
sigma_uz = 1.1e-9;
% Radial and axial diffusivity of carbon dioxide in pear tissue:
sigma_vr = 2.32e-9;
sigma_vz = 6.97e-9;

% universal gas constant
global R_g
R_g = 8.314;
% Maximum oxygen consumption rate:
global T
T = T_cel + 273.15;           % actual temperature in Kelvin
T_ref = 293.15;               % reference temperature
V_mu = 2.39e-4 * exp( 80200/R_g * (1/T_ref - 1/T) );
% Maximum fermentative carbon dioxide production rate:
V_mfv = 1.61e-4 * exp( 56700/R_g * (1/T_ref - 1/T) );
% Michaelis-Menten constants and respiration quotient:
K_mu = 0.4103;          % oxygen consumption
K_mv = 27.2438;         % non-competitive carbon dioxide inhibition
K_mfu = 0.1149;         % oxygen inhibition on fermentative carbon dioxide production
r_q = 0.97;             % respiration quotient

% Convective mass transfer coefficients:
rho_u = 7.0e-7;
rho_v = 7.5e-7;

% Ambient conditions:
global p_atm
p_atm = 101300;
C_uamb = p_atm * nu_u / (R_g * T);
C_vamb = p_atm * nu_v / (R_g * T);

save('var.mat');
vars = struct('V_mu',num2cell(V_mu),'K_mu',num2cell(K_mu),'K_mv',num2cell(K_mv),...
    'r_q',num2cell(r_q),'V_mfv',num2cell(V_mfv),'K_mfu',num2cell(K_mfu));

% generate mesh
global model
model = createpde;

load('pear_data.mat');% obtained via getPearShape.m
x = x/7; x = x - x(1);
y = y/7; y(end) = y(end-1); y = y-y(end);
pgon = polyshape(x,y);      % create polygon from (x,y) points
tr = triangulation(pgon);

tnodes = [x; y];
telements = tr.ConnectivityList';   


geometryFromMesh(model,tnodes,telements);   % create 
clear tnodes telements tr pgon;
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',0.005,'Hmin',0.002, 'Hgrad', 1.5);

global nodes
nodes = mesh.Nodes;

figure(1); clf;
subplot(121); pdegplot(model,'EdgeLabels','off'); ylim([0 0.2]); axis off;
subplot(122); pdemesh(model,'NodeLabels','on'); ylim([0 0.2]); axis off;



%% export mesh to text file

[p,e,t] = meshToPet(mesh);

Points = p;
edgeLabels = e(1,:);
triangleLabels = zeros(1, 3*length(t(1,:)));
for i=1:length(t(1,:)) 
    for j=1:3
        triangleLabels(3*(i-1)+j) = t(j,i);
    end
end 
edge1Labels = [1 3:33 2];		% vertical edge, top to bottom
edge2Labels = [1 80:-1:34 2];	% round edge, top to bottom
sizes = [length(Points(1,:)), length(triangleLabels), length(edge1Labels), length(edge2Labels) ];

writematrix(Points(1,:),'../c++/mesh1/Xpoints.txt','Delimiter',' ')  ;


writematrix(Points(2,:),'../c++/mesh1/Ypoints.txt','Delimiter',' ')  ;


writematrix(edge1Labels-1,'../c++/mesh1/edge1Labels.txt','Delimiter',' ');  


writematrix(edge2Labels-1,'../c++/mesh1/edge2Labels.txt','Delimiter',' ');  


writematrix(triangleLabels-1,'../c++/mesh1/triangleLabels.txt','Delimiter',' ');  
  

writematrix(sizes,'../c++/mesh1/sizes.txt','Delimiter',' ');  
  
        
%% Get initial solution (linearisation) & stiffness matrix


% Get stiffness matrix K and constant term f:
[K, K_lin, f, f_lin] = create_stiffness(mesh);

% First solution: 
C_0 = (K) \ -(f+f_lin);
%C_0 = K\-(f+f_lin_gross);

%C_0(C_0<0) = 0;

% Shift O2 upwards: 
%maximum = max(abs(C_0));
%C_0(1:length(nodes)) = maximum + C_0(1:length(nodes));



% Plot initial solution found by linearisation: 
C_0_plot = C_0 * R_g * T / p_atm;
figure(1); clf;
subplot(221); hold on;
pdeplot(model,'XYData',C_0_plot(1:length(nodes)),'Contour','on','ColorMap','jet');
title('O_2 concentration');
subplot(223);
pdeplot(model,'XYData',C_0_plot(1:length(nodes)),'ZData',C_0_plot(1:length(nodes)),'ColorMap','jet');
% scatter3(nodes(1,:), nodes(2,:), C_0(1:length(nodes)));
% trisurf(triangleLabels', nodes(1,:), nodes(2,:), C_0(1:length(nodes)));
% shading interp
% colorbar();



subplot(222);
hold on;
pdeplot(model,'XYData',C_0_plot(length(nodes)+1:end),'Contour','on','ColorMap','jet');
title('CO_2 concentration');
subplot(224);
pdeplot(model,'XYData',C_0_plot(length(nodes)+1:end),'ZData',C_0_plot(length(nodes)+1:end),'ColorMap','jet');
% scatter3(nodes(1,:), nodes(2,:), C_0(length(nodes)+1:end))
% trisurf(triangleLabels', nodes(1,:), nodes(2,:), C_0(length(nodes)+1:end))
% shading interp
% colorbar();

%% Test eval_nonlinear: 

%H = eval_nonlinear(mesh, [ones(length(nodes),1)*C_uamb; ones(length(nodes),1)*C_vamb], vars);
H = eval_nonlinear(mesh, C_0, vars);
figure(1); clf;
subplot(221); hold on;
pdeplot(model,'XYData',H(1:length(nodes)),'Contour','on','ColorMap','jet');
title('O_2 concentration');
subplot(223);
pdeplot(model,'XYData',H(1:length(nodes)),'ZData',H(1:length(nodes)),'ColorMap','jet');
% scatter3(nodes(1,:), nodes(2,:), C_0(1:length(nodes)));
% trisurf(triangleLabels', nodes(1,:), nodes(2,:), C_0(1:length(nodes)));
% shading interp
% colorbar();



subplot(222);
hold on;
pdeplot(model,'XYData',H(length(nodes)+1:end),'Contour','on','ColorMap','jet');
title('CO_2 concentration');
subplot(224);
pdeplot(model,'XYData',H(length(nodes)+1:end),'ZData',H(length(nodes)+1:end),'ColorMap','jet');

%% Solve nonlinear system (with intermediate plots) 
options = optimoptions('fsolve',...
    'Display','iter','FunctionTolerance',1e-20, 'OptimalityTolerance',1e-10,'UseParallel', true, 'OutputFcn', @outfun);

[C,fval,exitflag,output, J] = fsolve(fun,C_0,options);
% J is jacobian at solution of solver...

%% Plot result: 
Cplot = C * R_g * T / p_atm;
figure(1); clf;
subplot(221); hold on;
pdeplot(model,'XYData',Cplot(1:length(nodes)),'Contour','on','ColorMap','jet');
title('O_2 concentration: end result');
subplot(223);
pdeplot(model,'XYData',Cplot(1:length(nodes)),'ZData',Cplot(1:length(nodes)),'ColorMap','jet');
% scatter3(nodes(1,:), nodes(2,:), C(1:length(nodes)));
% trisurf(triangleLabels', nodes(1,:), nodes(2,:), C(1:length(nodes)));
% shading interp
% colorbar();



subplot(222);
hold on;
pdeplot(model,'XYData',Cplot(length(nodes)+1:end),'Contour','on','ColorMap','jet');
title('CO_2 concentration: end result');
subplot(224);
pdeplot(model,'XYData',Cplot(length(nodes)+1:end),'ZData',Cplot(length(nodes)+1:end),'ColorMap','jet');
% scatter3(nodes(1,:), nodes(2,:), C(length(nodes)+1:end))
% trisurf(triangleLabels', nodes(1,:), nodes(2,:), C(length(nodes)+1:end))
% shading interp
% colorbar();

%% Functions (load this before the rest...) 

fun = @(C) 10e4*( K*C + f + eval_nonlinear(mesh, C, vars));

function stop = outfun(C_, optimValues, stats)
    
    global model nodes
    figure(2); clf;
        global R_g T p_atm
        Cplot = C_ * (R_g * T / p_atm);
        
    
        subplot(121); hold on;
        pdeplot(model,'XYData',Cplot(1:length(nodes)));
        title('O_2 concentration: nonlinear');
        %scatter3(nodes(1,:), nodes(2,:), C_0(1:length(nodes)));
        %trisurf(triangleLabels', nodes(1,:), nodes(2,:), C_(1:length(nodes)));
        colorbar();


        subplot(122);
        hold on;
        pdeplot(model,'XYData',Cplot(length(nodes)+1:end));
        title('CO_2 concentration: nonlinear');
        %scatter3(nodes(1,:), nodes(2,:), C_0(length(nodes)+1:end))
        %trisurf(triangleLabels', nodes(1,:), nodes(2,:), C_(length(nodes)+1:end))
        colorbar();

    
    stop = false; 
     if (optimValues.iteration == 50)
         stop = true;
     end

end

