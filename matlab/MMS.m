%% simulation parameters (cf. table in assignement)
clear all; % close all; 
T_cel = 25;
nu_u = 20.8 /100;
nu_v = 0.04 /100;


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

% a = C_uamb / (2*sigma_ur/rho_u + 1);
% b = C_vamb / (2*sigma_vr/rho_v + 1);

save('var.mat');
vars = struct('V_mu',num2cell(V_mu),'K_mu',num2cell(K_mu),'K_mv',num2cell(K_mv),...
    'r_q',num2cell(r_q),'V_mfv',num2cell(V_mfv),'K_mfu',num2cell(K_mfu),...
    'C_uamb',num2cell(C_uamb),'C_vamb',num2cell(C_vamb));

%% generate circular mesh
semi_circ;


%%

M = length(p);
mesh = struct('Nodes', p, 'Elements', t(1:3,:));

global nodes
nodes = p;

[K_1, K_3, f, ~] = create_stiffness(mesh, vars);
K = K_1+K_3;
C_amb_vect = [C_uamb*ones(length(nodes),1);C_vamb*ones(length(nodes),1)];
f_lin = eval_nonlinear(mesh, C_amb_vect, vars); 

% First solution: 
C_0 = (K) \ -(f+0.85*f_lin);

%% gen ManSol

C_man = zeros(2*M, 1);
for i = 1:M
    R_sq = p(1,i)^2 + p(2,i)^2;
    C_man(i) = C_uamb * (-R_sq + 2*sqrt(R_sq));
    C_man(i+M) = C_vamb * (-R_sq + 2*sqrt(R_sq));
end

Cplot = C_0 - C_man;

norm(Cplot(1:length(p)))
norm(Cplot(length(p)+1:end))
figure(10);

subplot(321); scatter3(p(1,:), p(2,:), C_0(1:length(p))); title('C_0');
subplot(322); scatter3(p(1,:), p(2,:), C_0(length(p)+1:end));
subplot(323); scatter3(p(1,:), p(2,:), C_man(1:length(p))); title('manufactured');
subplot(324); scatter3(p(1,:), p(2,:), C_man(length(p)+1:end));
subplot(325); scatter3(p(1,:), p(2,:), Cplot(1:length(p))); title('C_0 - C_{man}');
subplot(326); scatter3(p(1,:), p(2,:), Cplot(length(p)+1:end));

%%
norm(K*C_man + f + eval_nonlinear(mesh, C_man, vars))

%% Solve nonlinear system (with intermediate plots) 
options = optimoptions('fsolve',...
    'Display','iter','FunctionTolerance',1e-20, 'OptimalityTolerance',1e-10,'UseParallel', true, 'OutputFcn', @outfun);

[C,fval,exitflag,output, J] = fsolve(fun,C_man,options);

%% Functions (load this before the rest...) 

fun = @(C) 10e4 * ( K*C + f + 0.01*eval_nonlinear(mesh, C, vars));
% fun = @(C) 10e4*( K*C + f + second_integral(nodes(1,:),nodes(2,:),tri,length(nodes),C,vars));

function stop = outfun(C_, optimValues, stats)
    
    global model nodes
    figure(2); clf;
    global R_g T p_atm 
    Cplot = 100* C_ * (R_g * T / p_atm);

    subplot(121); scatter3(nodes(1,:), nodes(2,:), C_(1:length(nodes)));
    subplot(122); scatter3(nodes(1,:), nodes(2,:), C_(length(nodes)+1:end));

    stop = false; 
    if (optimValues.iteration == 50)
        stop = true;
    end

end