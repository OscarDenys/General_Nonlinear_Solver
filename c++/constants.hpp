#ifndef constants
#define constants

namespace std {

    // Setup --------------------------------------------------------------------
    float T_cel = 25;
    float nu_u = 20.8 / 100;
    float nu_v = 0.04 / 100;
    std::vector<float> setup = {T_cel, nu_u, nu_v};

    // TODO: nadenken waar we deze constanten best definieeren...
    // Model parameters:


    // Radial and axial diffusivity of oxygen in pear tissue:
    float sigma_ur = 2.8e-10;
    float sigma_uz = 1.1e-9;
    // Radial and axial diffusivity of carbon dioxide in pear tissue:
    float sigma_vr = 2.32e-9;
    float sigma_vz = 6.97e-9;

    // Universal gas constant
    float R_g = 8.314;
    // Maximum oxygen consumption rate:
    float T = T_cel + 273.15;           // actual temperature in Kelvin
    float T_ref = 293.15;               // reference temperature
    float V_mu = 2.39e-4 * exp( 80200/R_g * (1/T_ref - 1/T) );
    // Maximum fermentative carbon dioxide production rate:
    float V_mfv = 1.61e-4 * exp( 56700/R_g * (1/T_ref - 1/T) );
    // Michaelis-Menten constants and respiration quotient:
    float K_mu = 0.4103;          // oxygen consumption
    float K_mv = 27.2438;         // non-competitive carbon dioxide inhibition
    float K_mfu = 0.1149;         // oxygen inhibition on fermentative carbon dioxide production
    float r_q = 0.97;             // respiration quotient

    // Convective mass transfer coefficients:
    float rho_u = 7.0e-7;
    float rho_v = 7.5e-7;

    // Ambient conditions:
    float p_atm = 101300;
    float C_uamb = p_atm * nu_u / (R_g * T);
    float C_vamb = p_atm * nu_v / (R_g * T);
    // -----------------------------------------------------------------------------

    // constants for evaluation the linearised second integral (cf. math derivation sec. 1.3-1.4)
    // all these constants use the prefix lin_
    double lin_k = K_mv*V_mu / (120* pow(K_mu+C_uamb,2)* pow(K_mv+C_vamb,2));
    double lin_alpha = 5* C_uamb* (K_mv*C_uamb + 2*C_uamb*C_vamb + K_mu*C_vamb );
    double lin_beta = K_mu* (K_mv + C_vamb);
    double lin_gamma = -C_uamb* (K_mu + C_uamb);
    double lin_kappa = K_mfu*V_mfv/ (120* pow(K_mfu+C_uamb,2));
    double lin_delta = 5* (K_mfu + 2*C_uamb);

} // namespace std

#endif