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
mesh = generateMesh(model);

figure(1);
subplot(121); pdegplot(model,'EdgeLabels','on'); ylim([0 1]); axis off;
subplot(122); pdemesh(model); ylim([0 1]); axis off;

%% finite elements
% TODO: change mesh to piecewise linear <-> quadratic

nodes = mesh.Nodes;
nb_nodes = length(nodes);
nb_elements_total = length(mesh.Elements(1,:));
%% Iterate over all test functions
A = zeros(2*nb_nodes, 2*nb_nodes);
b = zeros(2*nb_nodes, 1);
for testfunction_i = 1:nb_nodes
     % Get the node where the testfunction takes value 1 (zero in all other
     % nodes).
     centralpoint = nodes(testfunction_i);
     
     % Find elements that have centralpoint as one of its nodes.
     nonzero_elements = [];
     for element = 1:nb_elements_total
         if any(mesh.Elements(:,element) == testfunction_i)
            nonzero_elements = [nonzero_elements, mesh.Elements(element)];
         end 
     end
     
     % For each element in nonzero_elements: 
     % 1) Introduce local coordinate system
     %      Make sure that X1 is the node associated with the testfunction
     %      (testfunction_i). 
     % 2) Apply quadrature formula (Trapezoidal rule) and solve to obtain
     % equation in ci --> see mathematical derivation. 
     % 
     for j = 1:length(nonzero_elements(1,:))
         element = nonzero_elements(:,j);
         index_test_function = find(element == testfunction_i);
         
         node1 = element(index_test_function); % Should equal testfunction_i
         node2 = element(mod(index_test_function+1,3));
         node3 = element(mod(index_test_function+2,3));
         X1 = nodes(node1); % (x1,y1) = (r1,z1)
         X2 = nodes(node2);
         X3 = nodes(node3);
        
         % Formula for area of triangle given three vertices:
         area = abs((X1(1)*X2(2) + X2(1)*X3(2) + X3(1)*X1(2) - X1(1)*X3(2) - X2(1)*X1(2) - X3(1)*X2(2))/2);
         
         % Calculations of Jacobian:
         A = (X3(1) - X1(1))*(X2(2) - X1(2))/((X2(1)-X1(1))*(X3(2)-X1(2)));
         dr_dy = 2*area*(X3(1)-X1(1))/((1-A)*(X2(1)-X3(1))*(X3(2)-X1(2)));
         dr_dksi = 2*area/((1-A)*(X3(2)-X1(2)));
         dz_dy = 2*area/((1-A)*(X3(2)-X1(2)));
         dz_dksi = 2*area*(X2(2)-X1(2))/((1-A)*(X3(2)-X1(2))*(X2(1)-X1(1)));
         Jac = [[dr_dy, dr_dksi];[dz_dy, dz_dksi]];
         det_jac = det(Jac);
         
         % FIRST INTEGRAL of (5) --------------------------------------------
         
         % Calculations of D (see mathematical derivation):
         D = (r(0,0,X1,X2,X3)+r(0,1,X1,X2,X3))/4;
         
         % Coefficients of Ci:
         c1 = D*det_jac*(sigma_ur+sigma_uz);
         c2 = -D*det_jac*sigma_ur;
         c3 = -D*det_jac*sigma_uz;
         
         % SECOND INTEGRAL of (5) ----------------------------------------
         E = r(0,0,X1,X2,X3)*det_jac*V_mu;
         
         % NONLINEAR:
         %nonlinear_term = E*c(node1)/(K_mu+K_mu/K_mv*c(node1+nb_nodes)+c(node1)+c(node1)*c(node1+nb_nodes)/K_mv);
         
         % LINEAR:
         int_5_2_lin = -V_mu*det_jac*r(0,0,X1,X2,X3)/4;
         
         % FIRST INTEGRAL of (6) ----------------------------------------
         c1_hat = D*det_jac*(sigma_vr+sigma_vz);
         c2_hat = -D*det_jac*sigma_vr;
         c3_hat = -D*det_jac*sigma_vz;
         
         % SECOND INTEGRAL of (6) ----------------------------------------
         % NONLINEAR:
         %nonlinear_term = ...
         
         % LINEAR:
         int_6_2_lin = -r_q*V_mu*det_jac*r(0,0,X1,X2,X3)/4;
         
         % Alle termen toevoegen aan sparse system of equations
         % TODO: termen van lijnintegraal Sem nog toevoegen...
         % Ik denk dat het += is omdat we integraal van 1 testfunctie
         % hebben opgesplitst in meerdere bijdrages over verschillende
         % elementen. 
         A(testfunction_i,node1) = A(testfunction_i,node1) + c1;
         A(testfunction_i,node2) = A(testfunction_i,node2) + c2;
         A(testfunction_i,node3) = A(testfunction_i,node3) + c3;
         b(testfunction_i) = b(testfunction_i) + int_5_2_lin;
         A(nb_nodes+testfunction_i,node1) = A(nb_nodes+testfunction_i,node1) + c1_hat;
         A(nb_nodes+testfunction_i,node2) = A(nb_nodes+testfunction_i,node2) + c2_hat;
         A(nb_nodes+testfunction_i,node3) = A(nb_nodes+testfunction_i,node3) + c3_hat;
         b(nb_nodes+testfunction_i) = b(nd_nodes+testfunction_i) + int_6_2_lin;
     end
    
end

function [r] = r(ksi, y, X1, X2, X3)
    area = abs((X1(1)*X2(2) + X2(1)*X3(2) + X3(1)*X1(2) - X1(1)*X3(2) ...
                    - X2(1)*X1(2) - X3(1)*X2(2))/2);
    A = (X3(1) - X1(1))*(X2(2) - X1(2))/((X2(1)-X1(1))*(X3(2)-X1(2)));
    r = 1/((1-A)*(X3(2)-X1(2)))*(ksi*2*area+(X3(1)-X1(1))*...
        (y*2*area-(X2(2)-X1(2))*X1(1))/(X2(1)-X3(1)))+X1(1)/(1-A);
end 

function [z] = z(ksi, y, X1, X2, X3)
    area = abs((X1(1)*X2(2) + X2(1)*X3(2) + X3(1)*X1(2) - X1(1)*X3(2) ...
                    - X2(1)*X1(2) - X3(1)*X2(2))/2);
    A = (X3(1) - X1(1))*(X2(2) - X1(2))/((X2(1)-X1(1))*(X3(2)-X1(2)));
    z = 1/((1-A)*(X2(1)-X1(1)))*(y*2*area+(X2(2)-X1(2))*...
        (ksi*2*area+(X3(1)-X1(1))*X1(2))/(X3(2)-X1(2)))+X1(2)/(1-A);
end 


