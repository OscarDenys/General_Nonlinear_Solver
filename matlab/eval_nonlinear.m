function [b] = eval_nonlinear(mesh, C, vars)
    
    % finite elements:
    M = length(mesh.Nodes(1,:));
    nb_elements_total = length(mesh.Elements(1,:));

    % Iterate over all elements
    b = zeros(2*M, 1);

    for elem_index = 1:nb_elements_total  

        % For each element
        element = mesh.Elements(:,elem_index);     % node indexes of this element
        n1 = element(1);           % node indices
        n2 = element(2);
        n3 = element(3);
        P1 = mesh.Nodes(:,n1);          % node coordinates
        P2 = mesh.Nodes(:,n2);     
        P3 = mesh.Nodes(:,n3);
       

        dr_dy = P3(1)-P1(1);   
        dr_dksi = P2(1)-P1(1);
        dz_dy = P3(2) - P1(2);
        dz_dksi = P2(2) - P1(2);
        Jac = [[dr_dy, dr_dksi];[dz_dy, dz_dksi]];
        det_jac = abs(det(Jac)) ;
        
        [resp12u, resp12v] = evaluateR((C(n1)+C(n2))/2,(C(n1+M)+C(n2+M))/2,vars);
        [resp13u, resp13v] = evaluateR((C(n1)+C(n3))/2,(C(n1+M)+C(n3+M))/2,vars);
        [resp23u, resp23v] = evaluateR((C(n2)+C(n3))/2,(C(n2+M)+C(n3+M))/2,vars);

%         [resp12u, resp12v] = evaluateR_MMS((P1(1)+P2(1))/2,(P1(2)+P2(2))/2,vars);
%         [resp13u, resp13v] = evaluateR_MMS((P1(1)+P3(1))/2,(P1(2)+P3(2))/2,vars);
%         [resp23u, resp23v] = evaluateR_MMS((P2(1)+P3(1))/2,(P2(2)+P3(2))/2,vars);

        % =====================   integraal 2 - lineair (5)
        factor = 24;
        b(n1) = b(n1) + det_jac * ((P1(1) + P2(1))*resp12u + (P1(1)+P3(1))*resp13u)/factor;
        b(n2) = b(n2) + det_jac * ((P1(1) + P2(1))*resp12u + (P2(1)+P3(1))*resp23u)/factor;
        b(n3) = b(n3) + det_jac * ((P1(1) + P3(1))*resp13u + (P2(1)+P3(1))*resp23u)/factor;
        
        b(n1) = b(n1) + det_jac * ((P1(1) + P2(1))*resp12u + (P1(1)+P3(1))*resp13u) /24;
        b(n2) = b(n2) + det_jac * ((P1(1) + P2(1))*resp12u + (P2(1)+P3(1))*resp23u) /24;
        b(n3) = b(n3) + det_jac * ((P1(1) + P3(1))*resp13u + (P2(1)+P3(1))*resp23u) /24;
      
        % =====================   integraal 2 - lineair (6)
        n1 = element(1) + M;           % node index
        n2 = element(2) + M;
        n3 = element(3) + M;

        b(n1) = b(n1) - det_jac * ((P1(1) + P2(1))*resp12v + (P1(1)+P3(1))*resp13v)/factor;
        b(n2) = b(n2) - det_jac * ((P1(1) + P2(1))*resp12v + (P2(1)+P3(1))*resp23v)/factor;
        b(n3) = b(n3) - det_jac * ((P1(1) + P3(1))*resp13v + (P2(1)+P3(1))*resp23v)/factor;


    end
   

end

