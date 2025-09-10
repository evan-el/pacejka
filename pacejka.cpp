/** @file pacejka.cpp
 * @brief Library for computing force from the Pacejka model given a tir file.
 * 
 * Reference: Tyre and Vehicle Dynamics 2nd. ed., Hans B. Pacejka, p. 172-191. 
 * Initializes with a full path to a tir file, then calculates forces.
 * 
 * 
*/

#include "pacejka.hpp"

// Public

Pacejka::Pacejka(std::string tir_file) {
    readTIR(tir_file);
}

Pacejka::~Pacejka() {

}

double Pacejka::calcLongForceCombined(double fz, double kappa, double gamma, double alpha, double Vcx, double Vcy, double Vsx, double Vsy) {
    double gamma_star = calcGammaStar(gamma);
    double alpha_star = calcAlphaStar(alpha, Vcx);

    double FX_o = calcLongForcePure(fz, kappa, gamma, Vcx, Vcy, Vsx, Vsy);

    double dfz = calcDfz(fz);

    double SHX_alpha = tir_params["RHX1"];
    double EX_alpha = tir_params["REX1"] + tir_params["REX2"]*dfz;
    double CX_alpha = tir_params["RCX1"];
    double BX_alpha = ((tir_params["RBX1"] + tir_params["RBX3"]*pow(gamma_star,2)) \
        *cos(atan(tir_params["RBX2"]*kappa)))*tir_params["LXAL"];
    double alpha_s = alpha_star + SHX_alpha;
    double GX_alpha_o = cos(CX_alpha*atan(BX_alpha*SHX_alpha - \
        EX_alpha*(BX_alpha*SHX_alpha - atan(BX_alpha*SHX_alpha))));
    double GX_alpha = cos(CX_alpha*(atan(BX_alpha*alpha_s - EX_alpha*(BX_alpha*alpha_s - \
        atan(BX_alpha*alpha_s)))))/GX_alpha_o;
    double Fx = GX_alpha*FX_o;
    return Fx;
}

double Pacejka::calcLongForcePure(double fz, double kappa, double gamma, double Vcx, double Vcy, double Vsx, double Vsy) {
    double dfz = calcDfz(fz);

    double SHX = (tir_params["PHX1"] + tir_params["PHX2"]*dfz)*tir_params["LHX"];
    double kappa_x = kappa + tir_params["SHX"];
    double Cx = tir_params["PCX1"]*tir_params["LCX"];

    double Vs = sqrt(pow(Vsx, 2) + pow(Vsy, 2));

    double LMUX_star = calcLMUStar(tir_params["LMUX"], Vs, tir_params["LONGVL"]);

    double mu_x = tir_params["PDX1"] + tir_params["PDX2"]*LMUX_star;
    double Dx = mu_x*fz*zeta_1;
    double Ex = (tir_params["PEX1"] + tir_params["PEX2"]*dfz + tir_params["PEX3"]*pow(dfz, 2))* \
        (1.0 - tir_params["PEX4"]*sgn(kappa_x))*tir_params["LEX"];
    
    double KXK = fz*(tir_params["PKX1"] + tir_params["PKX2"]*dfz)*exp(tir_params["PKX3"]*dfz)*tir_params["LKX"];
    double Bx = KXK/(Cx*Dx + epsilon_Vx);
    double SVX = fz*(tir_params["PVX1"] + tir_params["PVX2"]*dfz)*(std::abs(Vcx)/(epsilon_Vx + std::abs(Vcx)))\
        *tir_params["LVX"]*tir_params["LMUX"]*zeta_1;
    double Fx_o = Dx*sin(Cx*atan(Bx*kappa_x - Ex*(Bx*kappa_x - atan(Bx*kappa_x)))) + SVX;
    return Fx_o;
}

double Pacejka::calcLatForceCombined(double fz, double kappa, double gamma, double alpha, double Vcx, double Vcy, double Vsx, double Vsy) {
    double dfz = calcDfz(fz);
    double Fyo = calcLatForcePure(fz, alpha, kappa, gamma, Vcx, Vcy, Vsx, Vsy);
    double Cy_kappa = tir_params["RCY1"];
    double Ey_kappa = tir_params["REY1"] + tir_params["REY2"]*dfz;
    double Vs = sqrt(pow(Vsx, 2) + pow(Vsy, 2));
    double Vc = sqrt(pow(Vcx, 2) + pow(Vcy, 2));
    double LMUY_star = calcLMUStar(tir_params["LMUY"], Vs, Vc);
    double gamma_star = calcGammaStar(gamma);
    double alpha_star = calcAlphaStar(alpha, Vcx);
    double MUy = ((tir_params["PDY1"] + tir_params["PDY2"]*dfz)/(1.0 + tir_params["PDY3"]*pow(gamma_star, 2)))*LMUY_star;
    double DVY_kappa = MUy*fz*(tir_params["RVY1"] + tir_params["RVY2"]*dfz + tir_params["RVY3"]*gamma_star)*\
        cos(atan(tir_params["RVY4"]*alpha_star))*zeta_2;
    double By_kappa = (tir_params["RBY1"] + tir_params["RBY4"]*pow(gamma, 2))*cos(atan(tir_params["RBY2"]*\
        (alpha_star - tir_params["RBY3"])))*tir_params["LYK"];
    double SVY_kappa = DVY_kappa*sin(tir_params["RVY5"]*atan(tir_params["RVY6"]*kappa))*tir_params["LVYKA"];
    double SHY_kappa = tir_params["RHY1"] + tir_params["RHY2"]*dfz;
    double kappa_s = kappa + SHY_kappa;
    double Gy_kappa_o = cos(Cy_kappa*atan(By_kappa*SHY_kappa - Ey_kappa*(By_kappa*SHY_kappa - atan(By_kappa*SHY_kappa))));
    double Gy_kappa = cos(Cy_kappa*atan(By_kappa*kappa_s - Ey_kappa*(By_kappa*kappa_s - atan(By_kappa*kappa_s))))/Gy_kappa_o;
    double Fy = Gy_kappa*Fyo + SVY_kappa;
    return Fy;
}

double Pacejka::calcLatForcePure(double fz, double kappa, double gamma, double alpha, double Vcx, double Vcy, double Vsx, double Vsy) {
    double dfz = calcDfz(fz);
    double Vs = sqrt(pow(Vsx, 2) + pow(Vsy, 2));
    double Vc = sqrt(pow(Vcx, 2) + pow(Vcy, 2));
    double gamma_star = calcGammaStar(gamma);
    double alpha_star = calcAlphaStar(alpha, Vcx);
    double Fzo_prime = calcFZOPrime(tir_params["FNOMIN"]);
    double LMUY_star = calcLMUStar(tir_params["LMUY"], Vs, Vc);
    double K_y_gamma_o = fz*(tir_params["PKY6"] + tir_params["PKY7"]*dfz)*tir_params["LKYC"];
    double SVY_gamma = fz*(tir_params["PVY3"] + tir_params["PVY4"]*dfz)*gamma_star*tir_params["LKYC"]*tir_params["LMUY_star"]*zeta_2;
    double KY_alpha = tir_params["PKY1"]*Fzo_prime*sin(tir_params["PKY4"]*atan(fz/((tir_params["PKY2"] \
        + tir_params["PKY5"]*pow(gamma_star, 2))*Fzo_prime)))/(1.0 + tir_params["PKY3"]*pow(gamma_star, 2))*zeta_3*tir_params["LYKA"];
    double SHY = (tir_params["PHY1"] + tir_params["PHY2"]*dfz)*tir_params["LHY"] + \
        (K_y_gamma_o*gamma_star - SVY_gamma)*zeta_0/(KY_alpha + epsilon_K) + zeta_4 - 1;
    double alpha_y = alpha_star + SHY;
    double SVY = fz*(tir_params["PVY1"] + tir_params["PVY2"]*dfz)*tir_params["LVY"]*zeta_2 + SVY_gamma;
    double Cy = tir_params["PCY1"]*tir_params["LCY"];
    double MUy = ((tir_params["PDY1"] + tir_params["PDY2"]*dfz)/(1.0 + tir_params["PDY3"]*pow(gamma_star, 2)))*LMUY_star;
    double Dy = MUy*fz*zeta_2;
    double Ey = (tir_params["PEY1"] + tir_params["PEY2"]*dfz)*(1.0 + tir_params["PEY5"]*pow(gamma_star, 2) - \
        (tir_params["PEY3"] + tir_params["PEY4"]*gamma_star)*sgn(alpha_y))*tir_params["LEY"];
    double By = KY_alpha/(Cy*Dy + epsilon_y);
    double Fy_o = Dy*sin(Cy*atan(By*alpha_y - Ey*(By*alpha_y - atan(By*alpha_y)))) + SVY;
    return Fy_o;
}

// Private

void Pacejka::readTIR(std::string tir_file) {
    std::cout << tir_file << std::endl;
    std::ifstream file(tir_file);
    std::stringstream buffer;
    buffer << file.rdbuf();
    file.close();

    // parse tir
    std::string line;
    std::vector<double> values;
    std::vector<std::string> param_names;
    int line_num = 1;
    while (getline(buffer, line, '\n')) {

        if (isalpha(line[0]) && (line_num<321)) {
            std::stringstream line_stream(line);
            std::string value;
            int element_index = 0;
            bool got_param_value = false;
            while (getline(line_stream, value, ' ')) {
                if (element_index==0) {
                    param_names.push_back(value);
                }
                if (!got_param_value && (element_index>1) && ((value[0]!=' ') && \
                        (value.length()!=0) && (value[0]!='='))) {
                            got_param_value = true;
                            if (value[0]!='\'') {
                                tir_params[param_names.back()] = stod(value);
                            }
                }
                element_index++;
            }
        }
        line_num++;
    }

    std::cout << "RHX1: " << tir_params["RHX1"] << std::endl;

}

double Pacejka::calcGammaStar(double gamma) {
    double gamma_star = sin(gamma);
    return gamma_star;
}

double Pacejka::calcAlphaStar(double alpha, double Vcx) {
    double alpha_star = tan(alpha)*sgn(Vcx);
    return alpha_star;
}

double Pacejka::calcDfz(double fz) {
    double dfz = (fz - tir_params["FNOMIN"])/tir_params["FNOMIN"];
    return dfz;
}

double Pacejka::calcLMUStar(double LMU, double Vs, double Vo) {
    double LMU_star = LMU/(1.0 + LMUV*(Vs/Vo));
    return LMU_star;
}

double Pacejka::calcFZOPrime(double Fzo) {
    double Fzo_prime = tir_params["LFZO"]*Fzo;
    return Fzo_prime;
}

double Pacejka::sgn(double num) {
    if (num < 0.0) {
        return -1.0;
    }
    else {
        return 1.0;
    }
}
