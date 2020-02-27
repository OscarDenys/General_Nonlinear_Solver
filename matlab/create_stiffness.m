function [K, f, f_lin] = create_stiffness(mesh)
    load('var.mat');
    % finite elements
    nodes = mesh.Nodes;
    nb_nodes = length(nodes);
    nb_elements_total = length(mesh.Elements(1,:));

    boundary_nodes = [1 58:-1:25 2]; % Hardcode Hmin = 0.05

    boundary_length = 0;
    for i = 1:length(boundary_nodes)-1
        node = boundary_nodes(i);
        next_node = boundary_nodes(i+1);
        boundary_length = boundary_length + sqrt( (nodes(1,node)-nodes(1,next_node))^2 + (nodes(2,node)-nodes(2,next_node))^2);
    end

    % Iterate over all elements
    stiffness_matrix = zeros(2*nb_nodes, 2*nb_nodes);
    b = zeros(2*nb_nodes, 1);
    f_lin = zeros(2*nb_nodes, 1);
    for elem_index = 1:nb_elements_total  

        % For each element
        elem_stiff = zeros(3); 
        element = mesh.Elements(:,elem_index);     % node indexes of this element
        n1 = element(1);           % node indices
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

        % =====================   integraal 1 - (5)
        temp = (P1(1) + P2(1) + P3(1)) / 6 / det_jac;     

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


   
        % =====================   integraal 2 - lineair (5-)
        F1 = V_mu * det_jac * (2*P1(1) + P2(1) + P3(1)) / 24;
        F2 = V_mu * det_jac * (2*P2(1) + P3(1) + P1(1)) / 24;
        F3 = V_mu * det_jac * (2*P3(1) + P1(1) + P2(1)) / 24;

        f_lin(n1) = f_lin(n1) + F1;
        f_lin(n2) = f_lin(n2) + F2;
        f_lin(n3) = f_lin(n3) + F3;
        
        
        % =====================   integraal 2 - lineair (6)
        n1 = element(1) + nb_nodes;           % node index
        n2 = element(2) + nb_nodes;
        n3 = element(3) + nb_nodes;
        
        
        F1 = - (r_q * V_mu + V_mfv) * det_jac * (2*P1(1) + P2(1) + P3(1)) / 24;
        F2 = - (r_q * V_mu + V_mfv) * det_jac * (2*P2(1) + P3(1) + P1(1)) / 24;
        F3 = - (r_q * V_mu + V_mfv) * det_jac * (2*P3(1) + P1(1) + P2(1)) / 24;

        f_lin(n1) = f_lin(n1) + F1;
        f_lin(n2) = f_lin(n2) + F2;
        f_lin(n3) = f_lin(n3) + F3;
        
        

        % =====================   integraal 1 - (6)
        temp = (P1(1) + P2(1) + P3(1)) / 6 / det_jac;     

        % node 1 -> test phi_1
        elem_stiff(1,1) = temp * (sigma_vr*power(P3(1)-P2(1),2) + sigma_vz*power(P3(2)-P2(2),2) ) ;  
        elem_stiff(1,2) = temp * (-sigma_vr*(P3(1)-P2(1))*(P3(1)-P1(1)) + sigma_vz*(P2(2)-P3(2))*(P3(2)-P1(2)) ); 
        elem_stiff(1,3) = temp * (sigma_vr*(P3(1)-P2(1))*(P2(1)-P1(1)) - sigma_vz*(P2(2)-P3(2))*(P2(2)-P1(2)) ); 

        % node 2 -> test phi_2
        elem_stiff(2,1) = temp * (-sigma_vr*(P3(1)-P2(1))*(P3(1)-P1(1)) + sigma_vz*(P2(2)-P3(2))*(P3(2)-P1(2)) );
        elem_stiff(2,2) = temp * (sigma_vr*power(P3(1)-P1(1),2) + sigma_vz*power(P3(2)-P1(2),2) ) ; 
        elem_stiff(2,3) = temp * (-sigma_vr*(P3(1)-P1(1))*(P2(1)-P1(1)) - sigma_vz*(P3(2)-P1(2))*(P2(2)-P1(2)) );

        % node 3 -> test phi_3
        elem_stiff(3,1) = temp * (sigma_vr*(P3(1)-P2(1))*(P2(1)-P1(1)) - sigma_vz*(P2(2)-P3(2))*(P2(2)-P1(2)) );
        elem_stiff(3,2) = temp * (-sigma_vr*(P3(1)-P1(1))*(P2(1)-P1(1)) - sigma_vz*(P3(2)-P1(2))*(P2(2)-P1(2)) );
        elem_stiff(3,3) = temp * (sigma_vr*power(P2(1)-P1(1),2) + sigma_vz*power(P2(2)-P1(2),2) ) ; 

        stiffness_matrix(n1,n1) = stiffness_matrix(n1,n1) + elem_stiff(1,1);
        stiffness_matrix(n1,n2) = stiffness_matrix(n1,n2) + elem_stiff(1,2);
        stiffness_matrix(n1,n3) = stiffness_matrix(n1,n3) + elem_stiff(1,3);
        stiffness_matrix(n2,n1) = stiffness_matrix(n2,n1) + elem_stiff(2,1);
        stiffness_matrix(n2,n2) = stiffness_matrix(n2,n2) + elem_stiff(2,2);
        stiffness_matrix(n2,n3) = stiffness_matrix(n2,n3) + elem_stiff(2,3);
        stiffness_matrix(n3,n1) = stiffness_matrix(n3,n1) + elem_stiff(3,1);
        stiffness_matrix(n3,n2) = stiffness_matrix(n3,n2) + elem_stiff(3,2);
        stiffness_matrix(n3,n3) = stiffness_matrix(n3,n3) + elem_stiff(3,3);        
    end


    % Integral 3 of (5) and (6) - boundary conditions: 
    for a = 0:1
        C_amb = C_uamb;
        rho = rho_u;
        if (a==1)
            C_amb = C_vamb;
            rho = rho_v;
        end
    %     Dt = 1/length(boundary_nodes);
        k = boundary_length;
        for i = 1:length(boundary_nodes)
            node = boundary_nodes(i);

            if (i~=1)
                prev_node = boundary_nodes(i-1);

                Dr = (nodes(1,node)-nodes(1,prev_node));
                Dz = (nodes(2,node)-nodes(2,prev_node));
    %             k = sqrt(Dr^2 + Dz^2) / Dt^2;
                Dt = sqrt(Dr^2 + Dz^2) / boundary_length;

                stiffness_matrix(a*nb_nodes+node,a*nb_nodes+prev_node) = stiffness_matrix(a*nb_nodes+node,a*nb_nodes+prev_node) + ...
                      rho* ( -Dt^4/4*Dr + ...
                              Dt^3/3*(nodes(1,node)-2*nodes(1,prev_node)) + ...
                              Dt^2/2*nodes(1,prev_node) ) * k;
                stiffness_matrix(a*nb_nodes+node,a*nb_nodes+node) = stiffness_matrix(a*nb_nodes+node,a*nb_nodes+node) + ...
                      rho* ( Dt^4/4*Dr + Dt^3/3*nodes(1,prev_node) ) * k;
                b(a*nb_nodes+node) = - rho* Dt^2*C_amb*k * (Dt/3*Dr - nodes(1,prev_node)/2);
            end

            if (i~=length(boundary_nodes))
                next_node = boundary_nodes(i+1);

                Dr = (nodes(1,next_node)-nodes(1,node));
                Dz = (nodes(2,next_node)-nodes(2,node));
                Dt = sqrt(Dr^2 + Dz^2) / boundary_length;        
    %             k = sqrt(Dr^2 + Dz^2) / Dt^2;

                stiffness_matrix(a*nb_nodes+node,a*nb_nodes+node) = stiffness_matrix(a*nb_nodes+node,a*nb_nodes+node) + ...
                         rho*  k * (-Dt^4*Dr/12 + ...
                                    Dt^3*( nodes(1,next_node)/2 ...
                                           -2*nodes(1,node)/3 ...
                                           -1/3 ) + ...
                                    Dt^2*nodes(1,node)/2 );
                stiffness_matrix(a*nb_nodes+node,a*nb_nodes+next_node) = stiffness_matrix(a*nb_nodes+node,a*nb_nodes+next_node) + ...
                          rho*      k * ( Dt^4*Dr/12 + Dt^3*nodes(1,node)/6 );
                b(a*nb_nodes+node) = b(a*nb_nodes+node) + ...
                          rho* k * C_amb * (Dt^3*(1/3 - Dr/2) - Dt^2*nodes(1,node)/2);
            end
        end
    end
    K = stiffness_matrix;
    f = b;
    f_lin = f + f_lin;
end

