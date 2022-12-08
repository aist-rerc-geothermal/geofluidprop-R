
#include "IAPWS-Viscosity-85.h"

#include <math.h>

static const double ViH1[4] = {
    1.0, 0.978197, 0.579829, -0.202354
};

static const double ViH2[6][7] = {
    { 0.5132047,  0.2151778, -0.2818107,  0.1778064, -0.04176610,  0.0,         0.0},
    { 0.3205656,  0.7317883, -1.070786,   0.4605040,  0.0,        -0.01578386,  0.0},
    { 0.0,        1.241044,  -1.263184,   0.2340379,  0.0,         0.0,         0.0},
    { 0.0,        1.476783,   0.0,       -0.4924179,  0.1600435,   0.0,        -0.003629481},
    {-0.7782567,  0.0,        0.0,        0.0,        0.0,         0.0,         0.0},
    { 0.1885447,  0.0,        0.0,        0.0,        0.0,         0.0,         0.0}
};

double iapws85_viscosity_rhoT(double rho, double tmpK)
{
    double tempAst = 647.226, densAst = 317.763, viscAst = 5.5071E-5;

    double tauR = tempAst / tmpK;
    double tauR1 = tauR - 1.0;
    double rx = rho / densAst;
    double rx1 = rx - 1.0;

    double vis00 = ViH1[0] + ViH1[1] * tauR + ViH1[2] * tauR * tauR + ViH1[3] * tauR * tauR * tauR;
    double vis0 = sqrt(tmpK / tempAst) / vis00;

    double xx = ViH2[0][0];

    double p2 = 1.0;
    for (int i = 1; i <= 6; i++)
    {
        p2 *= rx1;
        xx += ViH2[0][i] * p2;
    }

    double p1 = 1.0;
    for (int j = 1; j <= 5; j++)
    {
        p1 *= tauR1;
        xx += ViH2[j][0] * p1;

        p2 = 1.0;
        for (int i = 1; i <= 6; i++)
        {
            p2 *= rx1;
            xx += ViH2[j][i] * p1 * p2;
        }
    }

    return viscAst * vis0 * exp(rx * xx);
}
