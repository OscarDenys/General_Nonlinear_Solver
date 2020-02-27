function [r_] = r(ksi, y, P1, P2, P3)
    r_ = P1(1) * (1-ksi - y) + P2(1) * (ksi) + P3(1) * (y);
end 
