function [R_u,R_v] = respiration(index1,index2,prev_sol, vars)
% Function retuning formulae (3) on the point halfway between the nodes
% with index1 and index2. prev_sol is  the solution of last iteration.
    M = length(prev_sol)/2;
     assert(M>=index1);
     assert(M>=index2);
   
    C_u = (prev_sol(index1) + prev_sol(index2)) / 2;
    C_v = (prev_sol(M+index1) + prev_sol(M+index2)) / 2;
 
    
    
    R_u = vars.V_mu * C_u / ((vars.K_mu+C_u)*(1+ C_v/vars.K_mv));
    R_v = vars.r_q * R_u + vars.V_mfv / (1+ C_u/vars.K_mfu);
 
end

