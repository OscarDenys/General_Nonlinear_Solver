%% simulation parameters (cf. table in assignement)

T_cel = 25;
nu_u = 20.8 / 100;
nu_v = 0.04 / 100;


%% model parameters

% Radial and axial diffusivity of oxygen in pear tissue:
sigma_ur = 2.8e-10;
sigma_uz = 1.1e-9;
% Radial and axial diffusivity of carbon dioxide in pear tissue:
sigma_vr = 2.32e-9;
sigma_vz = 6.97e-9;

% universal gas constant
R_g = 8.314;
% Maximum oxygen consumption rate:
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
p_atm = 101300;
C_uamb = p_atm * nu_u / (R_g * T);
C_vamb = p_atm * nu_v / (R_g * T);


%% generate mesh

model = createpde;

load('pear_data.mat');      % obtained via getPearShape.m
pgon = polyshape(x,y);      % create polygon from (x,y) points
tr = triangulation(pgon);

tnodes = [x; y];
telements = tr.ConnectivityList';      

geometryFromMesh(model,tnodes,telements);   % create geometry
clear tnodes telements tr pgon;
mesh = generateMesh(model, 'GeometricOrder', 'linear');

figure(1);
subplot(121); pdegplot(model,'EdgeLabels','on'); ylim([0 1]); axis off;
subplot(122); pdemesh(model); ylim([0 1]); axis off;

%% finite elements
% TODO: change mesh to piecewise linear <-> quadratic

nodes = mesh.Nodes;
nb_nodes = length(nodes);
nb_elements_total = length(mesh.Elements(1,:));

%% Iterate over all test functions
stiffness_matrix = zeros(2*nb_nodes, 2*nb_nodes);
b = zeros(2*nb_nodes, 1);
for elem_index = 1:nb_elements_total  
    
    % For each element
    
    elem_stiff = zeros(3); 
    element = mesh.Elements(:,elem_index);     % node indexes of this element
    n1 = element(1);           % node index
    n2 = element(2);
    n3 = element(3);
    P1 = nodes(:,n1);          % node coordinates
    P2 = nodes(:,n2);     
    P3 = nodes(:,n3);

    dr_dy = P3(1)-P1(1);   
    dr_dksi = P2(1)-P1(1);
    dz_dy = P3(2) - P1(2);
    dz_dksi = P2(2) - P1(2);
    Jac = [[dr_dy, dr_dksi];[dz_dy, dz_dksi]];
    det_jac = det(Jac);
    
    % =====================   integraal 1
    temp = (r(0,0,P1,P2,P3) + r(0,1,P1,P2,P3)) / 4 / det_jac;     

    % node 1 -> test phi_1
    elem_stiff(1,1) = temp * (sigma_ur*power(P3(1)-P2(1),2) + sigma_uz*power(P3(2)-P2(2),2) ) ;  
    elem_stiff(1,2) = temp * (-sigma_ur*(P3(1)-P2(1))*(P3(1)-P1(1)) + sigma_uz*(P2(2)-P3(2))*(P3(2)-P1(2)) ); 
    elem_stiff(1,3) = temp * (sigma_ur*(P3(1)-P2(1))*(P2(1)-P1(1)) - sigma_uz*(P2(2)-P3(2))*(P2(2)-P1(2)) ); 

    % node 2 -> test phi_2
    elem_stiff(2,1) = temp * (-sigma_ur*(P3(1)-P2(1))*(P3(1)-P1(1)) + sigma_uz*(P2(2)-P3(2))*(P3(2)-P1(2)) );
    elem_stiff(2,2) = temp * (sigma_ur*power(P3(1)-P1(1),2) + sigma_uz*power(P3(2)-P1(2),2) ) ; 
    elem_stiff(2,3) = temp * (-sigma_ur*(P3(1)-P1(1))*(P2(1)-P1(1)) - sigma_uz*(P3(2)-P1(2))*(P2(2)-P1(2)) );

    % node 3 -> test phi_3
    elem_stiff(3,1) = temp * (sigma_ur*(P3(1)-P2(1))*(P2(1)-P1(1)) - sigma_uz*(P2(2)-P3(2))*(P2(2)-P1(2)) );
    elem_stiff(3,2) = temp * (-sigma_ur*(P3(1)-P1(1))*(P2(1)-P1(1)) - sigma_uz*(P3(2)-P1(2))*(P2(2)-P1(2)) );
    elem_stiff(3,3) = temp * (sigma_ur*power(P2(1)-P1(1),2) + sigma_uz*power(P2(2)-P1(2),2) ) ; 
     
    stiffness_matrix(n1,n1) = stiffness_matrix(n1,n1) + elem_stiff(1,1);
    stiffness_matrix(n1,n2) = stiffness_matrix(n1,n2) + elem_stiff(1,2);
    stiffness_matrix(n1,n3) = stiffness_matrix(n1,n3) + elem_stiff(1,3);
    stiffness_matrix(n2,n1) = stiffness_matrix(n2,n1) + elem_stiff(2,1);
    stiffness_matrix(n2,n2) = stiffness_matrix(n2,n2) + elem_stiff(2,2);
    stiffness_matrix(n2,n3) = stiffness_matrix(n2,n3) + elem_stiff(2,3);
    stiffness_matrix(n3,n1) = stiffness_matrix(n3,n1) + elem_stiff(3,1);
    stiffness_matrix(n3,n2) = stiffness_matrix(n3,n2) + elem_stiff(3,2);
    stiffness_matrix(n3,n3) = stiffness_matrix(n3,n3) + elem_stiff(3,3);
     
    
    % =====================   integraal 2 - lineair
    temp = det_jac * V_mu / 4;
    F1 = temp * ( r(0,0,P1,P2,P3));
    F2 = temp * ( r(1,0,P1,P2,P3));
%     F3 = 0;
    
    b(n1) = b(n1) + F1;
    b(n2) = b(n2) + F2;
    
end
   




function [r_] = r(ksi, y, P1, P2, P3)
    r_ = P1(1) * (1-ksi - y) + P2(1) * (ksi) + P3(1) * (y);
end 




