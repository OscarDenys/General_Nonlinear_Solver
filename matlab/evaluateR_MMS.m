function [R_u,R_v] = evaluateR_MMS(r,z, vars)
    
    R = sqrt(r^2+z^2);
    
    % Radial and axial diffusivity of oxygen in pear tissue:
    sigma_ur = 2.8e-10;
    sigma_uz = 1.1e-9;
    % Radial and axial diffusivity of carbon dioxide in pear tissue:
    sigma_vr = 2.32e-9;
    sigma_vz = 6.97e-9;

    R_u = vars.C_uamb * (sigma_ur*(-4+4/R-2*r^2/R^3) + sigma_uz*(-4+4/R-2*z^2/R^3))/2;
    R_v = - vars.C_vamb * (sigma_vr*(-4+4/R-2*r^2/R^3) + sigma_vz*(-4+4/R-2*z^2/R^3))/2;
    
end