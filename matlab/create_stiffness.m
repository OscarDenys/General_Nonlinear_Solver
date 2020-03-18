function [K, K_lin, f, f_lin] = create_stiffness(mesh)
    load('var.mat');
    % finite elements
    nodes = mesh.Nodes;
    nb_nodes = length(nodes);
    nb_elements_total = length(mesh.Elements(1,:));

    %boundary_nodes = [1 58:-1:25 2]; % Hardcode Hmin = 0.001;
    boundary_nodes = [1 80:-1:34 2]; %'Hmax',0.005,'Hmin',0.0005
    %boundary_nodes = [1 8:-1:5 2]; % Hardcode Hmin = 0.25
    
    % all constants with prefix lin_ are used for the linearised second
    % integral
    lin_k = V_mu/(K_mu+C_uamb)*(1+C_vamb/K_mv) /120;
%     lin_alpha = 5*C_uamb*(K_mv*C_uamb + ...
%                           2*C_uamb*C_vamb + ...
%                           K_mu*C_vamb);
%     lin_beta = K_mu*(K_mv+C_vamb);
%     lin_gamma = - C_uamb*(K_mu+C_uamb);
%     lin_kappa = K_mfu*V_mfv/(120*(K_mfu+C_uamb)^2);
%     lin_delta = 5*K_mfu + 10*C_uamb;
    

    % Iterate over all elements
    stiffness_matrix = zeros(2*nb_nodes, 2*nb_nodes);
    stiffness_matrix_lin = zeros(2*nb_nodes, 2*nb_nodes);
    b = zeros(2*nb_nodes, 1);
    f_lin = zeros(2*nb_nodes, 1);
    for elem_index = 1:nb_elements_total  

        % For each element
        elem_stiff = zeros(3); 
        element = mesh.Elements(:,elem_index);     % node indexes of this element
        n1 = element(1);           % node indices
        n2 = element(2);
        n3 = element(3);
        n1_ = element(1) + nb_nodes;  % phi indices of CO_2 concentrations
        n2_ = element(2) + nb_nodes;
        n3_ = element(3) + nb_nodes;
        P1 = nodes(:,n1);          % node coordinates
        P2 = nodes(:,n2);     
        P3 = nodes(:,n3);

        dr_dy = P3(1)-P1(1);   
        dr_dksi = P2(1)-P1(1);
        dz_dy = P3(2) - P1(2);
        dz_dksi = P2(2) - P1(2);
        Jac = [[dr_dy, dr_dksi];[dz_dy, dz_dksi]];
        det_jac = abs(det(Jac));

        % =====================   integraal 1 - (5)
        temp = (P1(1) + P2(1) + P3(1)) / 6 / det_jac;     

        % node 1 -> test phi_1
        elem_stiff(1,1) = temp * (sigma_uz*power(P3(1)-P2(1),2) + sigma_ur*power(P3(2)-P2(2),2) ) ;  
        elem_stiff(1,2) = temp * (-sigma_uz*(P3(1)-P2(1))*(P3(1)-P1(1)) + sigma_ur*(P2(2)-P3(2))*(P3(2)-P1(2)) ); 
        elem_stiff(1,3) = temp * (sigma_uz*(P3(1)-P2(1))*(P2(1)-P1(1)) - sigma_ur*(P2(2)-P3(2))*(P2(2)-P1(2)) ); 

        % node 2 -> test phi_2
        elem_stiff(2,1) = temp * (-sigma_uz*(P3(1)-P2(1))*(P3(1)-P1(1)) + sigma_ur*(P2(2)-P3(2))*(P3(2)-P1(2)) );
        elem_stiff(2,2) = temp * (sigma_uz*power(P3(1)-P1(1),2) + sigma_ur*power(P3(2)-P1(2),2) ) ; 
        elem_stiff(2,3) = temp * (-sigma_uz*(P3(1)-P1(1))*(P2(1)-P1(1)) - sigma_ur*(P3(2)-P1(2))*(P2(2)-P1(2)) );

        % node 3 -> test phi_3
        elem_stiff(3,1) = temp * (sigma_uz*(P3(1)-P2(1))*(P2(1)-P1(1)) - sigma_ur*(P2(2)-P3(2))*(P2(2)-P1(2)) );
        elem_stiff(3,2) = temp * (-sigma_uz*(P3(1)-P1(1))*(P2(1)-P1(1)) - sigma_ur*(P3(2)-P1(2))*(P2(2)-P1(2)) );
        elem_stiff(3,3) = temp * (sigma_uz*power(P2(1)-P1(1),2) + sigma_ur*power(P2(2)-P1(2),2) ) ; 

        stiffness_matrix(n1,n1) = stiffness_matrix(n1,n1) + elem_stiff(1,1);
        stiffness_matrix(n1,n2) = stiffness_matrix(n1,n2) + elem_stiff(1,2);
        stiffness_matrix(n1,n3) = stiffness_matrix(n1,n3) + elem_stiff(1,3);
        stiffness_matrix(n2,n1) = stiffness_matrix(n2,n1) + elem_stiff(2,1);
        stiffness_matrix(n2,n2) = stiffness_matrix(n2,n2) + elem_stiff(2,2);
        stiffness_matrix(n2,n3) = stiffness_matrix(n2,n3) + elem_stiff(2,3);
        stiffness_matrix(n3,n1) = stiffness_matrix(n3,n1) + elem_stiff(3,1);
        stiffness_matrix(n3,n2) = stiffness_matrix(n3,n2) + elem_stiff(3,2);
        stiffness_matrix(n3,n3) = stiffness_matrix(n3,n3) + elem_stiff(3,3);


   
        % =====================   integraal 2 - lineair (5)
        % constant part
        A1 = det_jac/24* 120*lin_k*C_uamb *(2*P1(1) + P2(1) + P3(1));
        A2 = det_jac/24* 120*lin_k*C_uamb *(P1(1) + 2*P2(1) + P3(1));
        A3 = det_jac/24* 120*lin_k*C_uamb *(P1(1) + P2(1) + 2*P3(1));
        
        f_lin(n1) = f_lin(n1) + A1;
        f_lin(n2) = f_lin(n2) + A2;
        f_lin(n3) = f_lin(n3) + A3;
        
        % part linear in c
%         elem_stiff(1,1) = (6*P1(1) + 2*P2(1) + 2*P3(1)) *lin_k*det_jac;
%         elem_stiff(1,2) = (2*P1(1) + 2*P2(1) +   P3(1)) *lin_k*det_jac;
%         elem_stiff(1,3) = (2*P1(1) +   P2(1) + 2*P3(1)) *lin_k*det_jac;
%         elem_stiff(2,1) = (2*P1(1) + 2*P2(1) +   P3(1)) *lin_k*det_jac;
%         elem_stiff(2,2) = (2*P1(1) + 6*P2(1) + 2*P3(1)) *lin_k*det_jac;
%         elem_stiff(2,3) = (  P1(1) + 2*P2(1) + 2*P3(1)) *lin_k*det_jac;
%         elem_stiff(3,1) = (2*P1(1) +   P2(1) + 2*P3(1)) *lin_k*det_jac;
%         elem_stiff(3,2) = (  P1(1) + 2*P2(1) + 2*P3(1)) *lin_k*det_jac;
%         elem_stiff(3,3) = (2*P1(1) + 2*P2(1) + 6*P3(1)) *lin_k*det_jac;
%         
%         stiffness_matrix_lin(n1,n1) = stiffness_matrix_lin(n1,n1) + elem_stiff(1,1);
%         stiffness_matrix_lin(n1,n2) = stiffness_matrix_lin(n1,n2) + elem_stiff(1,2);
%         stiffness_matrix_lin(n1,n3) = stiffness_matrix_lin(n1,n3) + elem_stiff(1,3);
%         stiffness_matrix_lin(n2,n1) = stiffness_matrix_lin(n2,n1) + elem_stiff(2,1);
%         stiffness_matrix_lin(n2,n2) = stiffness_matrix_lin(n2,n2) + elem_stiff(2,2);
%         stiffness_matrix_lin(n2,n3) = stiffness_matrix_lin(n2,n3) + elem_stiff(2,3);
%         stiffness_matrix_lin(n3,n1) = stiffness_matrix_lin(n3,n1) + elem_stiff(3,1);
%         stiffness_matrix_lin(n3,n2) = stiffness_matrix_lin(n3,n2) + elem_stiff(3,2);
%         stiffness_matrix_lin(n3,n3) = stiffness_matrix_lin(n3,n3) + elem_stiff(3,3);
       
        
        % =====================   integraal 2 - lineair (6)
        % constant part
        F1 = det_jac/24* (V_mfv/(1+C_uamb/K_mfu))* (2*P1(1) + P2(1) + P3(1));
        F2 = det_jac/24* (V_mfv/(1+C_uamb/K_mfu))* (P1(1) + 2*P2(1) + P3(1));
        F3 = det_jac/24* (V_mfv/(1+C_uamb/K_mfu))* (P1(1) + P2(1) + 2*P3(1));

        A1 = det_jac/24* (120*lin_k*C_uamb)* r_q *(2*P1(1) + P2(1) + P3(1));
        A2 = det_jac/24* (120*lin_k*C_uamb)* r_q *(P1(1) + 2*P2(1) + P3(1));
        A3 = det_jac/24* (120*lin_k*C_uamb)* r_q *(P1(1) + P2(1) + 2*P3(1));
        
        f_lin(n1_) = f_lin(n1_) - F1 - A1; % negatief vanwege vgl (6)
        f_lin(n2_) = f_lin(n2_) - F2 - A2;
        f_lin(n3_) = f_lin(n3_) - F3 - A3;
        
        % part linear in c
%         elem_stiff(1,1) = (6*P1(1) + 2*P2(1) + 2*P3(1)) *lin_k*det_jac*r_q;
%         elem_stiff(1,2) = (2*P1(1) + 2*P2(1) +   P3(1)) *lin_k*det_jac*r_q;
%         elem_stiff(1,3) = (2*P1(1) +   P2(1) + 2*P3(1)) *lin_k*det_jac*r_q;
%         elem_stiff(2,1) = (2*P1(1) + 2*P2(1) +   P3(1)) *lin_k*det_jac*r_q;
%         elem_stiff(2,2) = (2*P1(1) + 6*P2(1) + 2*P3(1)) *lin_k*det_jac*r_q;
%         elem_stiff(2,3) = (  P1(1) + 2*P2(1) + 2*P3(1)) *lin_k*det_jac*r_q;
%         elem_stiff(3,1) = (2*P1(1) +   P2(1) + 2*P3(1)) *lin_k*det_jac*r_q;
%         elem_stiff(3,2) = (  P1(1) + 2*P2(1) + 2*P3(1)) *lin_k*det_jac*r_q;
%         elem_stiff(3,3) = (2*P1(1) + 2*P2(1) + 6*P3(1)) *lin_k*det_jac*r_q;
%         
%         stiffness_matrix_lin(n1_,n1) = stiffness_matrix_lin(n1_,n1) + elem_stiff(1,1);
%         stiffness_matrix_lin(n1_,n2) = stiffness_matrix_lin(n1_,n2) + elem_stiff(1,2);
%         stiffness_matrix_lin(n1_,n3) = stiffness_matrix_lin(n1_,n3) + elem_stiff(1,3);
%         stiffness_matrix_lin(n2_,n1) = stiffness_matrix_lin(n2_,n1) + elem_stiff(2,1);
%         stiffness_matrix_lin(n2_,n2) = stiffness_matrix_lin(n2_,n2) + elem_stiff(2,2);
%         stiffness_matrix_lin(n2_,n3) = stiffness_matrix_lin(n2_,n3) + elem_stiff(2,3);
%         stiffness_matrix_lin(n3_,n1) = stiffness_matrix_lin(n3_,n1) + elem_stiff(3,1);
%         stiffness_matrix_lin(n3_,n2) = stiffness_matrix_lin(n3_,n2) + elem_stiff(3,2);
%         stiffness_matrix_lin(n3_,n3) = stiffness_matrix_lin(n3_,n3) + elem_stiff(3,3);        
        

        % =====================   integraal 1 - (6)
        temp = (P1(1) + P2(1) + P3(1)) / 6 / det_jac;     

        % node 1 -> test phi_1
        elem_stiff(1,1) = temp * (sigma_vz*power(P3(1)-P2(1),2) + sigma_vr*power(P3(2)-P2(2),2) ) ;  
        elem_stiff(1,2) = temp * (-sigma_vz*(P3(1)-P2(1))*(P3(1)-P1(1)) + sigma_vr*(P2(2)-P3(2))*(P3(2)-P1(2)) ); 
        elem_stiff(1,3) = temp * (sigma_vz*(P3(1)-P2(1))*(P2(1)-P1(1)) - sigma_vr*(P2(2)-P3(2))*(P2(2)-P1(2)) ); 

        % node 2 -> test phi_2
        elem_stiff(2,1) = temp * (-sigma_vz*(P3(1)-P2(1))*(P3(1)-P1(1)) + sigma_vr*(P2(2)-P3(2))*(P3(2)-P1(2)) );
        elem_stiff(2,2) = temp * (sigma_vz*power(P3(1)-P1(1),2) + sigma_vr*power(P3(2)-P1(2),2) ) ; 
        elem_stiff(2,3) = temp * (-sigma_vz*(P3(1)-P1(1))*(P2(1)-P1(1)) - sigma_vr*(P3(2)-P1(2))*(P2(2)-P1(2)) );

        % node 3 -> test phi_3
        elem_stiff(3,1) = temp * (sigma_vz*(P3(1)-P2(1))*(P2(1)-P1(1)) - sigma_vr*(P2(2)-P3(2))*(P2(2)-P1(2)) );
        elem_stiff(3,2) = temp * (-sigma_vz*(P3(1)-P1(1))*(P2(1)-P1(1)) - sigma_vr*(P3(2)-P1(2))*(P2(2)-P1(2)) );
        elem_stiff(3,3) = temp * (sigma_vz*power(P2(1)-P1(1),2) + sigma_vr*power(P2(2)-P1(2),2) ) ; 

        stiffness_matrix(n1_,n1_) = stiffness_matrix(n1_,n1_) + elem_stiff(1,1);
        stiffness_matrix(n1_,n2_) = stiffness_matrix(n1_,n2_) + elem_stiff(1,2);
        stiffness_matrix(n1_,n3_) = stiffness_matrix(n1_,n3_) + elem_stiff(1,3);
        stiffness_matrix(n2_,n1_) = stiffness_matrix(n2_,n1_) + elem_stiff(2,1);
        stiffness_matrix(n2_,n2_) = stiffness_matrix(n2_,n2_) + elem_stiff(2,2);
        stiffness_matrix(n2_,n3_) = stiffness_matrix(n2_,n3_) + elem_stiff(2,3);
        stiffness_matrix(n3_,n1_) = stiffness_matrix(n3_,n1_) + elem_stiff(3,1);
        stiffness_matrix(n3_,n2_) = stiffness_matrix(n3_,n2_) + elem_stiff(3,2);
        stiffness_matrix(n3_,n3_) = stiffness_matrix(n3_,n3_) + elem_stiff(3,3);        
    end


    % Integral 3 of (5) and (6) - boundary conditions: 
    for a = 0:1
        C_amb = C_uamb;
        rho = rho_u;
        if (a==1)
            C_amb = C_vamb;
            rho = rho_v;
        end
        Dt = 1/length(boundary_nodes);
        for i = 3:length(boundary_nodes)-2
            node = boundary_nodes(i);
%             if (i ~= 2)
                prev_node = boundary_nodes(i-1);

                Dr = (nodes(1,node)-nodes(1,prev_node));
                Dz = (nodes(2,node)-nodes(2,prev_node));
                k = sqrt(Dr^2 + Dz^2);

                stiffness_matrix(a*nb_nodes+node,a*nb_nodes+prev_node) = stiffness_matrix(a*nb_nodes+node,a*nb_nodes+prev_node) + ...
                      rho* ( nodes(1,node) + nodes(1,prev_node) ) * k / 12;
                stiffness_matrix(a*nb_nodes+node,a*nb_nodes+node) = stiffness_matrix(a*nb_nodes+node,a*nb_nodes+node) + ...
                      rho* ( nodes(1,node)/4 + nodes(1,prev_node)/12 ) * k;
                b(a*nb_nodes+node) = - rho* k * C_amb * (nodes(1,node)/3 + nodes(1,prev_node)/6);
%             end
%             if (i ~= length(boundary_nodes)-1)
                next_node = boundary_nodes(i+1);

                Dr = (nodes(1,next_node)-nodes(1,node));
                Dz = (nodes(2,next_node)-nodes(2,node));
                k = sqrt(Dr^2 + Dz^2);

                stiffness_matrix(a*nb_nodes+node,a*nb_nodes+node) = stiffness_matrix(a*nb_nodes+node,a*nb_nodes+node) + ...
                         rho*  k * ( nodes(1,node)/4 + nodes(1,next_node)/12);
                stiffness_matrix(a*nb_nodes+node,a*nb_nodes+next_node) = stiffness_matrix(a*nb_nodes+node,a*nb_nodes+next_node) + ...
                          rho*      k * ( nodes(1,node)/12 + nodes(1,next_node)/12 );
                b(a*nb_nodes+node) = b(a*nb_nodes+node) + ...
                          rho* k * C_amb * (- nodes(1,node)/3 - nodes(1,next_node)/6 );
%             end
        end
    end
    K = stiffness_matrix;
    K_lin = stiffness_matrix_lin;
    f = b;
    f_lin = f_lin;
end

