#include <cmath>
#include "matplot/matplot.h"

#include "pacejka.hpp"

using namespace matplot;

int main(int argc, char* argv[]) {
    Pacejka tire("example.tir");
    double fz = (1200*9.8/4.0);
    double gamma = -0.01; // camber angle
    double Vcx = 20.0; // long vel of contact patch (long. wheel speed)
    std::vector<double> tan_alpha;
    std::vector<double> alpha;
    std::vector<double> fy;
    std::vector<double> fx;
    std::vector<double> kappa;
    double Vsy, Vcy, Vsx;
    for (int i = 0; i < 200; i++) {
        kappa.push_back(-1.0 + 0.01*i);
        tan_alpha.push_back(-1.0 + 0.01*i);
        alpha.push_back(atan(tan_alpha[i]));
        Vsy = -Vcx*tan_alpha[i];
        Vcy = Vsy;
        Vsx = -Vcx*kappa[i];
        fx.push_back(tire.calcLongForceCombined(fz, kappa[i], gamma, alpha[i], Vcx, Vcy, Vsx, Vsy));
        fy.push_back(tire.calcLatForceCombined(fz, alpha[i], kappa[i], gamma, Vcx, Vcy, Vsx, Vsy));
        std::cout << fy[i] << std::endl;
        alpha[i] = alpha[i]*180.0/M_PI;
    }

    auto fig1 = figure();
    fig1->width(1000);
    fig1->height(600);
    plot(alpha, fy)->line_width(3.0);
    title("Tire Fy vs Lateral Slip Angle");
    xlabel("[Rad]");
    ylabel("[N]");
    show();

    auto fig2 = figure();
    fig2->width(1000);
    fig2->height(600);
    plot(kappa, fx)->line_width(3.0);
    title("Tire Fx vs Longitudinal Slip Ratio");
    ylabel("[N]");
    show();
    return 0;
}