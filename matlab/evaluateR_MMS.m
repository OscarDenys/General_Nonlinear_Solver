function [R_u,R_v] = evaluateR_MMS(r,z, vars)
    
    R = sqrt(r^2+z^2);
    % Convective mass transfer coefficients:
rho_u = 7.0e-7;
rho_v = 7.5e-7;
    % Radial and axial diffusivity of oxygen in pear tissue:
    sigma_ur = 2.8e-10;
%     sigma_uz = 1.1e-9;
    % Radial and axial diffusivity of carbon dioxide in pear tissue:
    sigma_vr = 2.32e-9;
%     sigma_vz = 6.97e-9;

    R_u = vars.C_uamb * 4*sigma_ur* rho_u / (2*sigma_ur+rho_u);
    R_v = - vars.C_vamb * 4*sigma_vr* rho_v / (2*sigma_vr+rho_v);
    
end