/** @file pacejka.hpp
 * @brief Library for computing force from the Pacejka model given a tir file.
 * 
 * Reference: Tyre and Vehicle Dynamics 2nd. ed., Hans B. Pacejka, p. 172-191. 
 * Initializes with a full path to a tir file, then calculates forces.
 * 
 * 
*/

#ifndef PACEJKA_HPP
#define PACEJKA_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>

class Pacejka {
    public:
        /** @brief Constructor for Pacejka class.
         * 
         * @param tir_file Full path to tir file.
        */
        Pacejka(std::string tir_file);
        ~Pacejka();
        
        /** @brief Calculates longitudinal force at the tire under combined (long. and lat.) slip conditions.
         * 
         * @param fz Vertical force on tire (Newtons).
         * @param kappa Longitudinal slip ratio.
         * @param gamma Camber (rad).
         * @param alpha Lateral slip angle (rad).
         * @param Vcx Velocity x of contact center (Wheel velocity).
         * @param Vcy Velocity y of contact center (Wheel velocity).
         * @param Vsx Slip velocity x.
         * @param Vsy Slip velocity y.
         * 
         * @return Fx, longitudinal force (N).
        */
        double calcLongForceCombined(double fz, double kappa, double gamma, double alpha, double Vcx, double Vcy, double Vsx, double Vsy);
        
        /** @brief Calculates longitudinal force at the tire under pure longitudinal slip conditions.
         * 
         * @param fz Vertical force on tire (Newtons).
         * @param kappa Longitudinal slip ratio.
         * @param gamma Camber (rad).
         * @param alpha Lateral slip angle (rad).
         * @param Vcx Velocity x of contact center (Wheel velocity).
         * @param Vcy Velocity y of contact center (Wheel velocity).
         * @param Vsx Slip velocity x.
         * @param Vsy Slip velocity y.
         * 
         * @return Fx_o, Longitudinal force under pure long. slip (N).
        */
        double calcLongForcePure(double fz, double kappa, double gamma, double Vcx, double Vcy, double Vsx, double Vsy);
        
        /** @brief Calculates lateral force at the tire under combined (long. and lat.) slip conditions.
         * 
         * @param fz Vertical force on tire (Newtons).
         * @param kappa Longitudinal slip ratio.
         * @param gamma Camber (rad).
         * @param alpha Lateral slip angle (rad).
         * @param Vcx Velocity x of contact center (Wheel velocity).
         * @param Vcy Velocity y of contact center (Wheel velocity).
         * @param Vsx Slip velocity x.
         * @param Vsy Slip velocity y.
         * 
         * @return Fy, lateral force (N).
        */
        double calcLatForceCombined(double fz, double kappa, double gamma, double alpha, double Vcx, double Vcy, double Vsx, double Vsy);
        
        /** @brief Calculates lateral force at the tire under pure lateral slip conditions.
         * 
         * @param fz Vertical force on tire (Newtons).
         * @param kappa Longitudinal slip ratio.
         * @param gamma Camber (rad).
         * @param alpha Lateral slip angle (rad).
         * @param Vcx Velocity x of contact center (Wheel velocity).
         * @param Vcy Velocity y of contact center (Wheel velocity).
         * @param Vsx Slip velocity x.
         * @param Vsy Slip velocity y.
         * 
         * @return FY_o, Lateral force under pure lat. slip (N).
        */
        double calcLatForcePure(double fz, double kappa, double gamma, double alpha, double Vcx, double Vcy, double Vsx, double Vsy);

    private:

        /** @brief Reads tir file and saves numeric values to a map.
         * 
         * @param tir_file Full path to tir file.
         * @return void
        */
        void readTIR(std::string tir_file);

        /** @brief Calculates gamma_star for use in Pacejka equations.
        *   @param gamma Camber angle (rad).
        *   @return gamma_star
        */
        double calcGammaStar(double gamma);

        /**
         * @brief Calculates alpha_star for use in Pacejka equations.
         * @param alpha Lateral slip angle (rad) 
         * @param Vcx Longitudinal velocity of contact patch (wheel speed) (m/s).
        */
        double calcAlphaStar(double alpha, double Vcx);

        /** @brief Calculates dfz (deviation from nominal fz).
         *  @param fz Vertical force (N).
         *  @return dfz
        */
        double calcDfz(double fz);

        /** @brief Calculates lambda mu star.
         *  @param LMU lambda mu, scaling factor from tir file.
         *  @param Vs Slip velocity total.
         *  @param Vo Nominal long. velocity from tir file.
         * 
         *  @return LMU_star
         * 
        */
        double calcLMUStar(double LMU, double Vs, double Vo);

        /** @brief Calculates modified FZO.
         *  @param Fzo Nominal vertical force from tir file
         *  @return Fzo_prime
        */
        double calcFZOPrime(double Fzo);
        
        /** @brief Returns the sign of the input.
         *  @param num Input number.
         *  @return Sign of num as double.
        */
        double sgn(double num);

        const double g = 9.81;

        std::map<std::string, double> tir_params;
        std::map<std::string, double> scaling_factors;

        // scaling factors and epsilons to prevent division by 0
        double epsilon_Vx = 0.01;
        double epsilon_K = 0.01;
        double epsilon_y = 0.01;
        double zeta_0 = 1.0;
        double zeta_1 = 1.0;
        double zeta_2 = 1.0;
        double zeta_3 = 1.0;
        double zeta_4 = 1.0;

        // can make this value larger than 0 to reflect a wet surface
        double LMUV = 0.0;
};

#endif