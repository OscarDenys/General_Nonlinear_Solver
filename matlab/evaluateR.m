function [R_u,R_v] = evaluateR(C_u,C_v, vars)
    
    R_u = vars.V_mu * C_u / ((vars.K_mu+C_u)*(1+ C_v/vars.K_mv));
    R_v = vars.r_q * R_u + vars.V_mfv / (1+ C_u/vars.K_mfu);
    
end