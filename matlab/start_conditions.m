%% simulation parameters (cf. table in assignement)
clear all;  close all; 
T_cel = 25;
nu_u =  20.8/100;
nu_v = 0.04/100;


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
T_ref = 293.15;  % reference temperature
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
    'r_q',num2cell(r_q),'V_mfv',num2cell(V_mfv),'K_mfu',num2cell(K_mfu), 'R_g', num2cell(R_g), 'T', num2cell(T), 'p_atm', num2cell(p_atm));

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
mesh = generateMesh(model,'GeometricOrder','linear','Hmax',0.004,'Hmin',0.001, 'Hgrad', 1.5);

global nodes
nodes = mesh.Nodes;

figure(1); clf;
subplot(121); pdegplot(model,'EdgeLabels','off'); ylim([-0.02 0.12]); axis off;
subplot(122); pdemesh(model,'NodeLabels','on'); ylim([-0.02 0.12]); axis off;

save('mesh.mat', 'mesh')

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
edge1Labels = [1 3:31 2];		% vertical edge, top to bottom
edge2Labels = [1 75:-1:32 2];	% round edge, top to bottom
sizes = [length(Points(1,:)), length(triangleLabels), length(edge1Labels), length(edge2Labels) ];

writematrix(Points(1,:),'../c++/mesh1/Xpoints.txt','Delimiter',' ')  ;
writematrix(Points(2,:),'../c++/mesh1/Ypoints.txt','Delimiter',' ')  ;
writematrix(edge1Labels-1,'../c++/mesh1/edge1Labels.txt','Delimiter',' ');  
writematrix(edge2Labels-1,'../c++/mesh1/edge2Labels.txt','Delimiter',' ');  
writematrix(triangleLabels-1,'../c++/mesh1/triangleLabels.txt','Delimiter',' ');  
writematrix(sizes,'../c++/mesh1/sizes.txt','Delimiter',' ');  
        
%% Get initial solution (linearisation) & stiffness matrix

a = 0.117644615194776;
b = 61.175199901283264;

% Get stiffness matrix K and constant term f:
[K, K_lin, f, f_lin] = create_stiffness(mesh, vars);
%[K, K_lin, f, f_lin] = create_stiffness_2(mesh);
% f = [a*f(1:length(nodes)); b*f(length(nodes)+1:end)];


% Third integral: 
%[fu,Ku,fv,Kv] = third_integral(edge2Labels, mesh.Nodes(1,:), mesh.Nodes(2,:), T, nu_u, nu_v);
% [fu;fv]-f is zero --> check
% norm(K(1:length(fu), 1:length(fu))-Ku) is zero --> check
% norm(K(length(fv)+1:end, length(fv)+1:end)-Kv) is zero --> check

%C_0_u =( K(1:length(fu), 1:length(fu))+Ku) \ (fu - f_lin(1:length(fu)));
%C_0_v = (K(length(fv)+1:end, length(fv)+1:end)+Kv )\ (fv - f_lin(length(fu)+1:end));
%C_0 = [C_0_u; C_0_v]; 

% First solution: 
C_0 = (K+K_lin) \ -(f+f_lin);
%C_0 = K\-(f+f_lin_gross);

CAMB = [C_uamb*ones(length(f)/2, 1); C_vamb*ones(length(f)/2, 1)];

% err = K * [(1/a)*ones(length(nodes),1); (1/b)*ones(length(nodes),1)];
% test = err - f;
% test1 = norm(err)
% test2 = norm(err - f)

% Shift O2 upwards: 
%  maximum = max(abs(C_0));
%  C_0 = maximum*ones(length(C_0),1)+C_0;
% for elem = 1:length(C_0)/2
%     C_0(elem) = 15*C_0(elem);
% end 
% for elem = length(C_0)/2+1:length(C_0)
%     C_0(elem) = 3*C_0(elem);
% end 



% Plot initial solution found by linearisation: 
C_0_plot = 100*C_0 * R_g * T / p_atm;
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



%% Solve nonlinear system (with intermediate plots) 
options = optimoptions('fsolve',...
    'Display','iter','FunctionTolerance',1e-20, 'OptimalityTolerance',1e-10,'UseParallel', true, 'OutputFcn', @outfun);

[C,fval,exitflag,output, J] = fsolve(fun,C_0,options);
% J is jacobian at solution of solver...

%% Plot result: 
Cplot = 100 * C * R_g * T / p_atm;
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


%% Functions (load this before the rest...) 

fun = @(C) 10e4*( K*C + f + eval_nonlinear(mesh, C, vars));

function stop = outfun(C_, optimValues, stats)
    
    global model nodes
    figure(2); clf;
        global R_g T p_atm
        Cplot = 100* C_ * (R_g * T / p_atm);
        
    
        subplot(121); hold on;
        ylim([0 0.12]); axis equal;
        pdeplot(model,'XYData',Cplot(1:length(nodes)),'Contour','on','ColorMap','jet');
        title('O_2 concentration: nonlinear');
        %scatter3(nodes(1,:), nodes(2,:), C_0(1:length(nodes)));
        %trisurf(triangleLabels', nodes(1,:), nodes(2,:), C_(1:length(nodes)));
        colorbar();


        subplot(122);
        hold on;
        ylim([0 0.12]); axis equal;
        pdeplot(model,'XYData',Cplot(length(nodes)+1:end),'Contour','on','ColorMap','jet');
        title('CO_2 concentration: nonlinear');
        %scatter3(nodes(1,:), nodes(2,:), C_0(length(nodes)+1:end))
        %trisurf(triangleLabels', nodes(1,:), nodes(2,:), C_(length(nodes)+1:end))
        colorbar();

    
    stop = false; 
     if (optimValues.iteration == 50)
         stop = true;
     end

end


