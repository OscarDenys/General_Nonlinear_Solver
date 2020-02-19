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
mesh = generateMesh(model,'GeometricOrder','linear','Hmin',0.25);

figure(1);
subplot(121); pdegplot(model,'EdgeLabels','on'); ylim([0 1]); axis off;
subplot(122); pdemesh(model,'NodeLabels','on'); ylim([0 1]); axis off;

%% export mesh to text file

[p,e,t] = meshToPet(mesh);

Points = p;
edgeLabels = e(1,:);
triangleLabels = t(1:3,:);

writematrix(Points,'points.txt','Delimiter',',')  
type points.txt

writematrix(edgeLabels,'edgeLabels.txt','Delimiter',',')  
type edgeLabels.txt

writematrix(triangleLabels,'triangleLabels.txt','Delimiter',',')  
type triangleLabels.txt


%% finite elements


% opstellen van lijnintegraal

nodes = mesh.Nodes;
nb_nodes = length(nodes);       % aantal knopen = aantal basisfuncties

A = zeros(2*nb_nodes);
gamma = zeros(2*nb_nodes,1);

%boundary_nodes = [1 58:-1:25 2];    % hard coded, andere methode nodig in C++

boundary_nodes = [1 8:-1:5 2];

boundary_length = 0;
for i = 1:length(boundary_nodes)-1
    node = boundary_nodes(i);
    next_node = boundary_nodes(i+1);
    boundary_length = boundary_length + sqrt( (nodes(1,node)-nodes(1,next_node))^2 + (nodes(2,node)-nodes(2,next_node))^2);
end


for a = 0:1
    C_amb = C_uamb;
    if (a==1)
        C_amb = C_vamb;
    end
    Dt = 1/length(boundary_nodes);
%     k = boundary_length;
    for i = 1:length(boundary_nodes)
        node = boundary_nodes(i);
        
        if (i~=1)
            prev_node = boundary_nodes(i-1);

            Dr = (nodes(1,node)-nodes(1,prev_node));
            Dz = (nodes(2,node)-nodes(2,prev_node));
            k = sqrt(Dr^2 + Dz^2) / Dt^2;
    %         Dt = sqrt(Dr^2 + Dz^2) / boundary_length;

            A(a*nb_nodes+node,a*nb_nodes+prev_node) = A(a*nb_nodes+node,a*nb_nodes+prev_node) + ...
                       ( -Dt^4/4*Dr + ...
                          Dt^3/3*(nodes(1,node)-2*nodes(1,prev_node)) + ...
                          Dt^2/2*nodes(1,prev_node) ) * k;
            A(a*nb_nodes+node,a*nb_nodes+node) = A(a*nb_nodes+node,a*nb_nodes+node) + ...
                       ( Dt^4/4*Dr + Dt^3/3*nodes(1,prev_node) ) * k;
            gamma(a*nb_nodes+node) = - Dt^2*C_amb*k * (Dt/3*Dr - nodes(1,prev_node)/2);
        end
        
        if (i~=length(boundary_nodes))
            next_node = boundary_nodes(i+1);
        
            Dr = (nodes(1,next_node)-nodes(1,node));
            Dz = (nodes(2,next_node)-nodes(2,node));
    %         Dt = sqrt(Dr^2 + Dz^2) / boundary_length;        
            k = sqrt(Dr^2 + Dz^2) / Dt^2;
            
            A(a*nb_nodes+node,a*nb_nodes+node) = A(a*nb_nodes+node,a*nb_nodes+node) + ...
                           k * (-Dt^4*Dr/12 + ...
                                Dt^3*( nodes(1,next_node)/2 ...
                                       -2*nodes(1,node)/3 ...
                                       -1/3 ) + ...
                                Dt^2*nodes(1,node)/2 );
            A(a*nb_nodes+node,a*nb_nodes+next_node) = A(a*nb_nodes+node,a*nb_nodes+next_node) + ...
                                k * ( Dt^4*Dr/12 + Dt^3*nodes(1,node)/6 );
            gamma(a*nb_nodes+node) = gamma(a*nb_nodes+node) + ...
                           k * C_amb * (Dt^3*(1/3 - Dr/2) - Dt^2*nodes(1,node)/2);
        end
    end
end






