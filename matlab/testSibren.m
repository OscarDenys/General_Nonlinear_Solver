clear all; close all; 
load('Sibren');

pdemesh(mesh,'NodeLabels','on')

load('var.mat');
vars = struct('V_mu',num2cell(V_mu),'K_mu',num2cell(K_mu),'K_mv',num2cell(K_mv),...
    'r_q',num2cell(r_q),'V_mfv',num2cell(V_mfv),'K_mfu',num2cell(K_mfu));

C_0 = [ones(length(Nodes),1)*C_uamb; ones(length(Nodes),1)*C_vamb];
nodes = mesh.Nodes;
%% Test eval_nonlinear: 

%H = eval_nonlinear(mesh, [ones(length(nodes),1)*C_uamb; ones(length(nodes),1)*C_vamb], vars);
H = eval_nonlinear(mesh, C_0, vars);
figure(1); clf;
subplot(221); hold on;
pdeplot(mesh,'XYData',H(1:length(nodes)),'Contour','on','ColorMap','jet');
title('O_2 concentration');
subplot(223);
pdeplot(mesh,'XYData',H(1:length(nodes)),'ZData',H(1:length(nodes)),'ColorMap','jet');
% scatter3(nodes(1,:), nodes(2,:), C_0(1:length(nodes)));
% trisurf(triangleLabels', nodes(1,:), nodes(2,:), C_0(1:length(nodes)));
% shading interp
% colorbar();



subplot(222);
hold on;
pdeplot(mesh,'XYData',H(length(nodes)+1:end),'Contour','on','ColorMap','jet');
title('CO_2 concentration');
subplot(224);
pdeplot(mesh,'XYData',H(length(nodes)+1:end),'ZData',H(length(nodes)+1:end),'ColorMap','jet');

%% Agathe
figure(2); clf;
subplot(221); hold on;
pdeplot(mesh,'XYData',Sibren(1:length(nodes)),'Contour','on','ColorMap','jet');
title('O_2 concentration');
subplot(223);
pdeplot(mesh,'XYData',Sibren(1:length(nodes)),'ZData',Sibren(1:length(nodes)),'ColorMap','jet');
% scatter3(nodes(1,:), nodes(2,:), C_0(1:length(nodes)));
% trisurf(triangleLabels', nodes(1,:), nodes(2,:), C_0(1:length(nodes)));
% shading interp
% colorbar();



subplot(222);
hold on;
pdeplot(mesh,'XYData',Sibren(length(nodes)+1:end),'Contour','on','ColorMap','jet');
title('CO_2 concentration');
subplot(224);
pdeplot(mesh,'XYData',Sibren(length(nodes)+1:end),'ZData',Sibren(length(nodes)+1:end),'ColorMap','jet');


