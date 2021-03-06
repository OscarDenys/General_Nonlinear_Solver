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
C_0 = (K+K_lin) \ -(f+f_lin/100);
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

%% Plot linear output C0 C++: 
C_0 = [ -16.1892 -21.8251 -52.8162 -86.5367 -117.618 -146.388 -173.213 -198.254 -221.578 -243.184 -263.03 -281.035 -297.09 -311.066 -322.82 -332.197 -339.037 -343.174 -344.442 -342.676 -337.711 -329.383 -317.528 -301.979 -282.567 -259.116 -231.44 -199.341 -162.594 -120.953 -74.1103 -21.5299 -21.0062 -20.6272 -19.8718 -18.4652 -16.5529 -15.119 -13.9063 -12.4264 -10.747 -10.115 -9.08997 -7.40388 -6.72598 -6.11963 -5.41298 -4.86355 -4.66446 -4.7035 -5.03581 -5.38915 -5.79111 -6.1938 -7.54231 -9.25837 -10.1995 -10.9545 -11.8217 -12.9759 -14.5892 -12.3331 -10.6574 -9.30859 -8.82238 -9.47902 -9.79719 -10.178 -10.0032 -9.33969 -8.99845 -9.50953 -12.5244 -12.9757 -14.824 -227.587 -232.64 -230.491 -223.237 -52.976 -57.2799 -60.5803 -46.6449 -61.8098 -37.2424 -64.8321 -67.4529 -44.2907 -71.3329 -45.1786 -73.8311 -73.9374 -69.0917 -64.3771 -101.184 -61.5677 -128.847 -58.2872 -156.277 -55.6648 -182.094 -53.4012 -206.239 -51.3329 -228.739 -50.3568 -249.533 -50.5103 -268.572 -49.3682 -285.741 -47.3529 -300.914 -47.7243 -313.951 -48.8398 -324.699 -48.3441 -332.995 -48.7048 -338.674 -48.6388 -341.567 -47.2435 -341.509 -48.1994 -338.332 -51.3164 -331.874 -51.4842 -321.971 -45.586 -308.461 -55.7303 -291.178 -56.616 -269.954 -50.458 -244.611 -60.3719 -214.957 -63.2496 -180.089 -63.0122 -140.886 -66.4734 -98.73 -65.5975 -58.2027 -75.8228 -82.2574 -91.7242 -100.483 -104.265 -68.301 -109.094 -114.575 -120.094 -125.286 -130.053 -127.859 -118.705 -111.206 -103.756 -105.818 -135.281 -100.99 -160.179 -95.2607 -185.127 -92.0527 -207.453 -91.1489 -230.145 -88.3637 -250.303 -87.3893 -268.712 -87.094 -285.194 -79.2041 -299.658 -87.1272 -311.937 -86.4624 -321.859 -84.7732 -329.254 -84.6587 -333.948 -84.6202 -335.769 -78.0299 -334.552 -86.3437 -330.132 -87.8067 -322.347 -84.292 -311.038 -296.048 -88.8111 -277.216 -91.5973 -254.377 -226.832 -99.3683 -195.833 -102.567 -157.972 -98.9847 -114.881 -73.7068 -110.859 -131.244 -137.318 -144.839 -151.501 -161.495 -166.444 -171.89 -183.734 -169.912 -160.009 -151.91 -145.085 -138.954 -159.899 -126.819 -176.55 -131.05 -203.988 -127.685 -223.6 -118.24 -243.254 -123.659 -263.101 -236.336 -279.177 -119.166 -293.19 -122.07 -304.933 -120.661 -314.252 -114.288 -320.964 -120.584 -324.868 -116.517 -325.81 -323.624 -122.039 -318.146 -122.581 -309.217 -116.781 -296.688 -280.413 -126.44 -260.244 -129.363 -235.028 -205.474 -136.148 -175.566 -139.287 -175.463 -191.029 -206.218 -200.427 -214.597 -214.301 -203.29 -196.195 -187.27 -180.757 -169.084 -166.478 -157.503 -233.069 -156.317 -250.926 -153.23 -266.985 -281.755 -293.24 -152.308 -301.866 -308.442 -151.973 -311.736 -151.63 -311.999 -309.039 -154.438 -302.664 -153.182 -292.753 -279.174 -261.799 -159.518 -239.193 -159.066 -211.81 -175.038 -155.471 -265.596 -239.884 -221.109 -220.882 -213.214 -205.005 -196.299 -191.371 -185.871 -252.043 -275.21 -183.576 -286.18 -293.059 -175.895 -294.968 -180.995 -294.426 -291.111 -184.581 -283.977 -187.567 -273.201 -258.705 -239.511 -187.581 -215.625 -180.602 -190.052 -254.208 -247.376 -243.905 -228.006 -224.895 -216.845 -203.077 -272.157 -273.919 -273.698 -270.048 -214.07 -263.097 -251.474 -235.506 -207.052 -212.826 -271.938 -209.326 -252.899 -248.098 -248.54 -239.562 10.8046 16.4294 15.0937 19.1373 22.9852 26.637 30.1186 33.4384 36.5932 39.573 42.3631 44.9457 47.2993 49.4008 51.2257 52.7487 53.9445 54.7879 55.2545 55.3206 54.9637 54.1614 52.8921 51.1342 48.8655 46.0628 42.7015 38.7547 34.1907 28.9734 23.0607 16.2611 15.9506 15.5734 14.9654 14.0456 12.9525 12.0395 11.2241 10.335 9.46631 8.90064 8.20964 7.38112 6.87833 6.46026 6.05735 5.75676 5.63185 5.66224 5.84074 6.07772 6.36783 6.74904 7.4914 8.41919 9.1382 9.72742 10.307 10.9371 11.3988 10.659 9.76435 8.99397 8.58667 8.60199 8.64278 8.66944 8.5077 8.18667 8.0244 8.22832 9.11971 9.62645 10.3094 40.064 38.9787 39.8663 38.9627 13.97 14.6521 15.2047 13.1279 15.4568 12.0692 15.9446 16.4373 13.391 17.2324 13.9581 17.9674 18.4614 18.1987 17.5414 20.8386 16.8988 24.3227 16.0996 27.8615 15.3563 31.2622 14.605 34.5055 13.8519 37.5854 13.3373 40.4844 13.0843 43.1895 12.7193 45.6782 12.2921 47.9276 12.3043 49.9138 12.5141 51.6112 12.5594 52.9939 12.8188 54.0363 13.0725 54.7133 13.1927 55.0009 13.7558 54.8761 14.6506 54.3165 15.1205 53.3006 14.7568 51.8069 16.7302 49.8142 17.4549 47.3003 17.1762 44.242 19.1261 40.6133 20.2125 36.2945 20.7562 31.3929 21.6121 26.0971 21.7662 20.9815 17.7206 17.9115 19.1318 20.3489 20.9569 16.2389 21.7031 22.557 23.4613 24.4005 25.3667 25.4637 24.5655 23.6185 20.9433 22.7538 25.0222 21.8651 28.3048 20.7896 31.645 20.0217 34.6895 19.6064 37.8309 18.9429 40.6723 18.5955 43.3133 18.4274 45.7251 17.1624 47.8905 18.3784 49.7813 18.3483 51.3702 18.229 52.6308 18.4182 53.5369 18.681 54.0638 18.0009 54.1881 19.5668 53.8874 20.2228 53.1402 20.2792 51.9257 50.2235 21.6092 48.013 22.7275 45.2728 41.9109 24.5366 38.0908 25.5787 33.3501 25.548 27.9103 17.2466 21.7215 24.4035 25.2626 26.3519 27.3591 28.8289 29.6574 30.6363 32.5152 31.056 29.9796 28.9699 27.9853 27.0036 28.2479 25.1232 30.5244 25.5434 34.2894 24.8733 37.0462 23.3201 39.8388 24.0134 42.6927 39.8634 45.0608 23.2619 47.1702 23.7364 48.9889 23.6317 50.492 22.8163 51.6516 23.979 52.4394 23.7269 52.8335 52.8115 24.9854 52.3513 25.526 51.4325 25.1944 50.036 48.1432 27.0876 45.7351 28.1166 42.6605 39.0087 29.647 35.3145 30.6209 30.4776 32.6653 34.8532 34.2741 36.3571 36.7756 35.5793 34.7981 33.6831 32.7867 31.0736 30.5575 29.1315 38.6245 28.8774 41.2234 28.4199 43.5908 45.7995 47.5741 28.5513 48.9675 50.0905 28.8037 50.775 29.1201 51.0537 50.9015 29.9492 50.2946 30.3217 49.2177 47.6542 45.5881 31.8204 42.8289 32.2899 39.4349 34.8353 28.8478 43.8116 40.5414 38.1809 38.2353 37.2256 36.0672 34.7931 34.0585 33.2686 41.7759 45.3158 33.3666 46.9871 48.1264 32.5091 48.5855 33.5776 48.715 48.4871 34.4838 47.7439 35.3942 46.5183 44.8001 42.4595 35.9387 39.4776 32.5514 34.0915 42.6597 41.8823 41.5041 39.315 38.9094 37.8018 36.5002 45.3634 45.7465 45.8911 45.5788 38.8026 44.8415 43.4709 41.5132 36.456 37.645 45.1882 37.7332 42.8636 42.3165 42.5748 41.5139 ];
C_0 = 100 * C_0 * R_g * T / p_atm;
figure(1); clf;
subplot(221); hold on;
pdeplot(model,'XYData',C_0(1:length(nodes)),'Contour','on','ColorMap','jet');
title('O_2 concentration: end result');
subplot(223);
pdeplot(model,'XYData',C_0(1:length(nodes)),'ZData',C_0(1:length(nodes)),'ColorMap','jet');
% scatter3(nodes(1,:), nodes(2,:), C(1:length(nodes)));
% trisurf(triangleLabels', nodes(1,:), nodes(2,:), C(1:length(nodes)));
% shading interp
% colorbar();



subplot(222);
hold on;
pdeplot(model,'XYData',C_0(length(nodes)+1:end),'Contour','on','ColorMap','jet');
title('CO_2 concentration: end result');
subplot(224);
pdeplot(model,'XYData',C_0(length(nodes)+1:end),'ZData',C_0(length(nodes)+1:end),'ColorMap','jet');
% scatter3(nodes(1,:), nodes(2,:), C(length(nodes)+1:end))
% trisurf(triangleLabels', nodes(1,:), nodes(2,:), C(length(nodes)+1:end))
% shading interp
% colorbar();

%% Plot nonlinear output C++: 
C_nonlin = [ 3.08719 -21.8251 -52.8162 -86.5367 -117.618 -146.388 -173.213 -198.254 -221.578 -243.184 -263.03 -281.035 -297.09 -311.066 -322.82 -332.197 -339.037 -343.174 -344.442 -342.676 -337.711 -329.383 -317.528 -301.979 -282.567 -259.116 -231.44 -199.341 -162.594 -120.953 -74.1103 -21.5299 -21.0062 -20.6272 -19.8718 -18.4652 -16.5529 -15.119 -13.9063 -12.4264 -10.747 -10.115 -9.08997 -7.40388 -6.72598 -6.11963 -5.41298 -4.86355 -4.66446 -4.7035 -5.03581 -5.38915 -5.79111 -6.1938 -7.54231 -9.25837 -10.1995 -10.9545 -11.8217 -12.9759 -14.5892 -12.3331 -10.6574 -9.30859 -8.82238 -9.47902 -9.79719 -10.178 -10.0032 -9.33969 -8.99845 -9.50953 -12.5244 -12.9757 -14.824 -227.587 -232.64 -230.491 -223.237 -52.976 -57.2799 -60.5803 -46.6449 -61.8098 -37.2424 -64.8321 -67.4529 -44.2907 -71.3329 -45.1786 -73.8311 -73.9374 -69.0917 -64.3771 -101.184 -61.5677 -128.847 -58.2872 -156.277 -55.6648 -182.094 -53.4012 -206.239 -51.3329 -228.739 -50.3568 -249.533 -50.5103 -268.572 -49.3682 -285.741 -47.3529 -300.914 -47.7243 -313.951 -48.8398 -324.699 -48.3441 -332.995 -48.7048 -338.674 -48.6388 -341.567 -47.2435 -341.509 -48.1994 -338.332 -51.3164 -331.874 -51.4842 -321.971 -45.586 -308.461 -55.7303 -291.178 -56.616 -269.954 -50.458 -244.611 -60.3719 -214.957 -63.2496 -180.089 -63.0122 -140.886 -66.4734 -98.73 -65.5975 -58.2027 -75.8228 -82.2574 -91.7242 -100.483 -104.265 -68.301 -109.094 -114.575 -120.094 -125.286 -130.053 -127.859 -118.705 -111.206 -103.756 -105.818 -135.281 -100.99 -160.179 -95.2607 -185.127 -92.0527 -207.453 -91.1489 -230.145 -88.3637 -250.303 -87.3893 -268.712 -87.094 -285.194 -79.2041 -299.658 -87.1272 -311.937 -86.4624 -321.859 -84.7732 -329.254 -84.6587 -333.948 -84.6202 -335.769 -78.0299 -334.552 -86.3437 -330.132 -87.8067 -322.347 -84.292 -311.038 -296.048 -88.8111 -277.216 -91.5973 -254.377 -226.832 -99.3683 -195.833 -102.567 -157.972 -98.9847 -114.881 -73.7068 -110.859 -131.244 -137.318 -144.839 -151.501 -161.495 -166.444 -171.89 -183.734 -169.912 -160.009 -151.91 -145.085 -138.954 -159.899 -126.819 -176.55 -131.05 -203.988 -127.685 -223.6 -118.24 -243.254 -123.659 -263.101 -236.336 -279.177 -119.166 -293.19 -122.07 -304.933 -120.661 -314.252 -114.288 -320.964 -120.584 -324.868 -116.517 -325.81 -323.624 -122.039 -318.146 -122.581 -309.217 -116.781 -296.688 -280.413 -126.44 -260.244 -129.363 -235.028 -205.474 -136.148 -175.566 -139.287 -175.463 -191.029 -206.218 -200.427 -214.597 -214.301 -203.29 -196.195 -187.27 -180.757 -169.084 -166.478 -157.503 -233.069 -156.317 -250.926 -153.23 -266.985 -281.755 -293.24 -152.308 -301.866 -308.442 -151.973 -311.736 -151.63 -311.999 -309.039 -154.438 -302.664 -153.182 -292.753 -279.174 -261.799 -159.518 -239.193 -159.066 -211.81 -175.038 -155.471 -265.596 -239.884 -221.109 -220.882 -213.214 -205.005 -196.299 -191.371 -185.871 -252.043 -275.21 -183.576 -286.18 -293.059 -175.895 -294.968 -180.995 -294.426 -291.111 -184.581 -283.977 -187.567 -273.201 -258.705 -239.511 -187.581 -215.625 -180.602 -190.052 -254.208 -247.376 -243.905 -228.006 -224.895 -216.845 -203.077 -272.157 -273.919 -273.698 -270.048 -214.07 -263.097 -251.474 -235.506 -207.052 -212.826 -271.938 -209.326 -252.899 -248.098 -248.54 -239.562 10.8046 16.4294 15.0937 19.1373 22.9852 26.637 30.1186 33.4384 36.5932 39.573 42.3631 44.9457 47.2993 49.4008 51.2257 52.7487 53.9445 54.7879 55.2545 55.3206 54.9637 54.1614 52.8921 51.1342 48.8655 46.0628 42.7015 38.7547 34.1907 28.9734 23.0607 16.2611 15.9506 15.5734 14.9654 14.0456 12.9525 12.0395 11.2241 10.335 9.46631 8.90064 8.20964 7.38112 6.87833 6.46026 6.05735 5.75676 5.63185 5.66224 5.84074 6.07772 6.36783 6.74904 7.4914 8.41919 9.1382 9.72742 10.307 10.9371 11.3988 10.659 9.76435 8.99397 8.58667 8.60199 8.64278 8.66944 8.5077 8.18667 8.0244 8.22832 9.11971 9.62645 10.3094 40.064 38.9787 39.8663 38.9627 13.97 14.6521 15.2047 13.1279 15.4568 12.0692 15.9446 16.4373 13.391 17.2324 13.9581 17.9674 18.4614 18.1987 17.5414 20.8386 16.8988 24.3227 16.0996 27.8615 15.3563 31.2622 14.605 34.5055 13.8519 37.5854 13.3373 40.4844 13.0843 43.1895 12.7193 45.6782 12.2921 47.9276 12.3043 49.9138 12.5141 51.6112 12.5594 52.9939 12.8188 54.0363 13.0725 54.7133 13.1927 55.0009 13.7558 54.8761 14.6506 54.3165 15.1205 53.3006 14.7568 51.8069 16.7302 49.8142 17.4549 47.3003 17.1762 44.242 19.1261 40.6133 20.2125 36.2945 20.7562 31.3929 21.6121 26.0971 21.7662 20.9815 17.7206 17.9115 19.1318 20.3489 20.9569 16.2389 21.7031 22.557 23.4613 24.4005 25.3667 25.4637 24.5655 23.6185 20.9433 22.7538 25.0222 21.8651 28.3048 20.7896 31.645 20.0217 34.6895 19.6064 37.8309 18.9429 40.6723 18.5955 43.3133 18.4274 45.7251 17.1624 47.8905 18.3784 49.7813 18.3483 51.3702 18.229 52.6308 18.4182 53.5369 18.681 54.0638 18.0009 54.1881 19.5668 53.8874 20.2228 53.1402 20.2792 51.9257 50.2235 21.6092 48.013 22.7275 45.2728 41.9109 24.5366 38.0908 25.5787 33.3501 25.548 27.9103 17.2466 21.7215 24.4035 25.2626 26.3519 27.3591 28.8289 29.6574 30.6363 32.5152 31.056 29.9796 28.9699 27.9853 27.0036 28.2479 25.1232 30.5244 25.5434 34.2894 24.8733 37.0462 23.3201 39.8388 24.0134 42.6927 39.8634 45.0608 23.2619 47.1702 23.7364 48.9889 23.6317 50.492 22.8163 51.6516 23.979 52.4394 23.7269 52.8335 52.8115 24.9854 52.3513 25.526 51.4325 25.1944 50.036 48.1432 27.0876 45.7351 28.1166 42.6605 39.0087 29.647 35.3145 30.6209 30.4776 32.6653 34.8532 34.2741 36.3571 36.7756 35.5793 34.7981 33.6831 32.7867 31.0736 30.5575 29.1315 38.6245 28.8774 41.2234 28.4199 43.5908 45.7995 47.5741 28.5513 48.9675 50.0905 28.8037 50.775 29.1201 51.0537 50.9015 29.9492 50.2946 30.3217 49.2177 47.6542 45.5881 31.8204 42.8289 32.2899 39.4349 34.8353 28.8478 43.8116 40.5414 38.1809 38.2353 37.2256 36.0672 34.7931 34.0585 33.2686 41.7759 45.3158 33.3666 46.9871 48.1264 32.5091 48.5855 33.5776 48.715 48.4871 34.4838 47.7439 35.3942 46.5183 44.8001 42.4595 35.9387 39.4776 32.5514 34.0915 42.6597 41.8823 41.5041 39.315 38.9094 37.8018 36.5002 45.3634 45.7465 45.8911 45.5788 38.8026 44.8415 43.4709 41.5132 36.456 37.645 45.1882 37.7332 42.8636 42.3165 42.5748 41.5139 ];
C_nonlin = 100 * C_nonlin * R_g * T / p_atm;
figure(1); clf;
subplot(221); hold on;
pdeplot(model,'XYData',C_nonlin(1:length(nodes)),'Contour','on','ColorMap','jet');
title('O_2 concentration: end result');
subplot(223);
pdeplot(model,'XYData',C_nonlin(1:length(nodes)),'ZData',C_nonlin(1:length(nodes)),'ColorMap','jet');
% scatter3(nodes(1,:), nodes(2,:), C(1:length(nodes)));
% trisurf(triangleLabels', nodes(1,:), nodes(2,:), C(1:length(nodes)));
% shading interp
% colorbar();



subplot(222);
hold on;
pdeplot(model,'XYData',C_nonlin(length(nodes)+1:end),'Contour','on','ColorMap','jet');
title('CO_2 concentration: end result');
subplot(224);
pdeplot(model,'XYData',C_nonlin(length(nodes)+1:end),'ZData',C_nonlin(length(nodes)+1:end),'ColorMap','jet');
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


