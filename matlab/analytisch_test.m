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
