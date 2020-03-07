#ifndef constants_hpp
#define constants_hpp
#include <cmath>
#include <vector>

namespace std {

    // Setup --------------------------------------------------------------------
    float T_cel = 25.0;
    float nu_u = 20.8 / 100;
    float nu_v = 0.04 / 100;
    std::vector<float> setup = {T_cel, nu_u, nu_v};


    // --- MODEL PARAMETERS:
    // Radial and axial diffusivity of oxygen in pear tissue:
    const float sigma_ur = 2.8e-10;
    const float sigma_uz = 1.1e-9;
    // Radial and axial diffusivity of carbon dioxide in pear tissue:
    const float sigma_vr = 2.32e-9;
    const float sigma_vz = 6.97e-9;

    // Universal gas constant
    const float R_g = 8.314;
    // Maximum oxygen consumption rate:
    const float T = T_cel + 273.15;           // actual temperature in Kelvin
    const float T_ref = 293.15;               // reference temperature
    const float V_mu = 2.39e-4 * exp( 80200/R_g * (1/T_ref - 1/T) );
    // Maximum fermentative carbon dioxide production rate:
    const float V_mfv = 1.61e-4 * exp( 56700/R_g * (1/T_ref - 1/T) );
    // Michaelis-Menten constants and respiration quotient:
    const float K_mu = 0.4103;          // oxygen consumption
    const float K_mv = 27.2438;         // non-competitive carbon dioxide inhibition
    const float K_mfu = 0.1149;         // oxygen inhibition on fermentative carbon dioxide production
    const float r_q = 0.97;             // respiration quotient

    // Convective mass transfer coefficients:
    const float rho_u = 7.0e-7;
    const float rho_v = 7.5e-7;

    // Ambient conditions:
    const float p_atm = 101300;
    const float C_uamb = p_atm * nu_u / (R_g * T);
    const float C_vamb = p_atm * nu_v / (R_g * T);
    // -----------------------------------------------------------------------------

    // constants for evaluation the linearised second integral (cf. math derivation sec. 1.3-1.4)
    // all these constants use the prefix lin_
    // double lin_k = K_mv*V_mu / (120* pow(K_mu+C_uamb,2)* pow(K_mv+C_vamb,2));
    // double lin_alpha = 5* C_uamb* (K_mv*C_uamb + 2*C_uamb*C_vamb + K_mu*C_vamb );
    // double lin_beta = K_mu* (K_mv + C_vamb);
    // double lin_gamma = -C_uamb* (K_mu + C_uamb);
    // double lin_kappa = K_mfu*V_mfv/ (120* pow(K_mfu+C_uamb,2));
    // double lin_delta = 5* (K_mfu + 2*C_uamb);

} // namespace std

#endif