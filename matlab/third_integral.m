function [fu,Ku,fv,Kv] = third_integral(edges, x, y, T, eta_u, eta_v)
%edges: M by 2 matrix of index pairs for the points representing the edge
%x,y: coordinates of the mesh points

rho_u = 7 * 10^-7;
rho_v = 7.5 * 10^-7;
Cu_amb = ambientOx(T, eta_u);
Cv_amb = ambientCarb(T, eta_v);
M = length(x);
Ku = zeros(M);
fu = zeros(M,1);
Kv = zeros(M);
fv = zeros(M,1);

% NOTE: t1 = 1-t2
%gausskwadratuur punten voor 2 puntsformule (ook functiewaarden in deze punten voor de
%rechten y=x en y = 1-x)(dus de phi's)
t1= 1/2 - 1/(2*sqrt(3));
t2= 1/2 + 1/(2*sqrt(3));
fprod = t1*t2; %phi1*phi2 symmetrisch tov t= 0.5, in beide punten dus dezelfde waarde

for i = 1:length(edges)-1 
    %coordinatentransformatie van parametrisatie recht lijnstuk
    g_deriv_norm = sqrt( ( x(edges(i+1))-x(edges(i)))^2 + (y(edges(i+1))-y(edges(i)) )^2); 
    %x van stijgende rechte (0-->1) in de gausskwadratuurpunten
    %de andere basisfunctie heeft dezelfde waarden, maar omgekeerd, want
    %dalend van 1 naar 0
    x1 = x(edges(i)) * t2 + x(edges(i+1))*t1;
    x2 = x(edges(i)) * t1 + x(edges(i+1))*t2;
    
    %bijdragen aan K matrix: op rand 2 basisfuncties niet 0 --> 4 bijdragen
    %diagonaalbijdragen i = j (basisfctie indices)
    %1e basisfctie dalend (1->0), simpele gausskwadratuur voor interval [0-1]: (f(t1) +
    %f(t2))/2, met t1 en t2 de kwadratuurpunten
    Ku(edges(i),edges(i)) = Ku(edges(i),edges(i)) + rho_u *g_deriv_norm/2*(t2^2 *x1+ t1^2*x2);
    %2e basisfctie stijgend (0->1)
    Ku(edges(i+1),edges(i+1)) = Ku(edges(i+1),edges(i+1)) + rho_u * g_deriv_norm/2*(t1^2 *x1+ t2^2*x2);
    %bijdragen kruistermen. Functiewaarden (fprod) van phi_i*phi_j is hier gelijk
    %in gausspunten wegens symmetrie rond 0.5 
    Ku(edges(i),edges(i+1)) = Ku(edges(i),edges(i+1)) + rho_u * g_deriv_norm/2*fprod*(x1+x2);
    Ku(edges(i+1),edges(i)) = Ku(edges(i+1),edges(i)) + rho_u * g_deriv_norm/2*fprod*(x1+x2);
    
    %analoog voor 2de deel vgl
    Kv(edges(i),edges(i)) = Kv(edges(i),edges(i)) + rho_v * g_deriv_norm/2*(t2^2 *x1+ t1^2*x2);
    Kv(edges(i+1),edges(i+1)) = Kv(edges(i+1),edges(i+1)) + rho_v * g_deriv_norm/2*(t1^2 *x1+ t2^2*x2);
    Kv(edges(i),edges(i+1)) = Kv(edges(i),edges(i+1)) + rho_v * g_deriv_norm/2*fprod*(x1+x2);
    Kv(edges(i+1),edges(i)) = Kv(edges(i+1),edges(i)) + rho_v * g_deriv_norm/2*fprod*(x1+x2);
    
    %bijdrage van -Camb deel van derde integraal, ook gausskwadratuur
    fu(edges(i),1) = fu(edges(i),1) - rho_u * g_deriv_norm/2*Cu_amb*(t2*x1+t1*x2);
    fu(edges(i+1),1) = fu(edges(i+1),1) - rho_u * g_deriv_norm/2*Cu_amb*(t1*x1+t2*x2);
    
    fv(edges(i),1) = fv(edges(i),1) - rho_v * g_deriv_norm/2*Cv_amb*(t2*x1+t1*x2);
    fv(edges(i+1),1) = fv(edges(i+1),1) - rho_v * g_deriv_norm/2*Cv_amb*(t1*x1+t2*x2);

end
end

% Returns ambient oxygen concentration given the ambient temperature 
% and the fraction of oxygen in the atmosphere.
function Cua = ambientOx(T,etau)
    patm = 101300;
    Rg = 8.314;
    Cua = patm*etau/(Rg*T);
end

% Returns ambient carbon dioxide concentration given the ambient temperature 
% and the fraction of carbon dioxide in the atmosphere.
function Cva = ambientCarb(T,etav)
    patm = 101300;
    Rg = 8.314;
    Cva = patm*etav/(Rg*T);
end

