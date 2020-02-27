function [f_lin] = lin_start_conditions(mesh)
% finite elements
    nodes = mesh.Nodes;
    nb_nodes = length(nodes);
    nb_elements_total = length(mesh.Elements(1,:));

    boundary_nodes = [1 8:-1:5 2]; % Hardcode 

    boundary_length = 0;
    for i = 1:length(boundary_nodes)-1
        node = boundary_nodes(i);
        next_node = boundary_nodes(i+1);
        boundary_length = boundary_length + sqrt( (nodes(1,node)-nodes(1,next_node))^2 + (nodes(2,node)-nodes(2,next_node))^2);
    end

    % Iterate over all elements
    b = zeros(2*nb_nodes, 1);
    for elem_index = 1:nb_elements_total  

        % For each element
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

        % =====================   integraal 2 - lineair (5-)
        F1 = V_mu * det_jac * (2*P1(1) + P2(1) + P3(1)) / 24;
        F2 = V_mu * det_jac * (2*P2(1) + P3(1) + P1(1)) / 24;
        F3 = V_mu * det_jac * (2*P3(1) + P1(1) + P2(1)) / 24;

        b(n1) = b(n1) + F1;
        b(n2) = b(n2) + F2;
        b(n3) = b(n3) + F3;
        
        
        % =====================   integraal 2 - lineair (6)
        n1 = element(1) + nb_nodes;           % node index
        n2 = element(2) + nb_nodes;
        n3 = element(3) + nb_nodes;
        
        F1 = - (r_q * V_mu + V_mfv) * det_jac * (2*P1(1) + P2(1) + P3(1)) / 24;
        F2 = - (r_q * V_mu + V_mfv) * det_jac * (2*P2(1) + P3(1) + P1(1)) / 24;
        F3 = - (r_q * V_mu + V_mfv) * det_jac * (2*P3(1) + P1(1) + P2(1)) / 24;

        b(n1) = b(n1) + F1;
        b(n2) = b(n2) + F2;
        b(n3) = b(n3) + F3;


    end
    
    f_lin = b;

end

