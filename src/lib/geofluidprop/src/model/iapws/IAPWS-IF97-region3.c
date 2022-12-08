// Region 3
//  623.15 <= T <= T(p)
//  p(T) <= p <= 100MPa

#include "IAPWS-IF97.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

#include "IAPWS-IF97-const.h"

//double R3I1[41], R3J1[41], R3An[41];

static const double R3I1[41] = {
    0,
    0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,6,6,6,7,8,9,9,10,10,11
};

static const double R3J1[41] = {
    0,
    0,0,1,2,7,10,12,23,2,6,15,17,0,2,6,7,22,26,0,2,4,16,26,0,2,4,26,1,3,26,0,2,26,2,26,2,26,0,1,26
};

//double R3J1[41], R3An[41];
static const double R3An[41] = {
    0,
    0.10658070028513E+01,-0.15732845290239E+02,0.20944396974307E+02,-0.76867707878716E+01,
    0.26185947787954E+01,-0.28080781148620E+01,0.12053369696517E+01,-0.84566812812502E-02,
    -0.12654315477714E+01,-0.11524407806681E+01,0.88521043984318E+00,-0.64207765181607E+00,
    0.38493460186671E+00,-0.85214708824206E+00,0.48972281541877E+01,-0.30502617256965E+01,
    0.39420536879154E-01,0.12558408424308E+00,-0.27999329698710E+00,0.13899799569460E+01,
    -0.20189915023570E+01,-0.82147637173963E-02,-0.47596035734923E+00,0.43984074473500E-01,
    -0.44476435428739E+00,0.90572070719733E+00,0.70522450087967E+00,0.10770512626332E+00,
    -0.32913623258954E+00,-0.50871062041158E+00,-0.22175400873096E-01,0.94260751665092E-01,
    0.16436278447961E+00,-0.13503372241348E-01,-0.14834345352472E-01,0.57922953628084E-03,
    0.32308904703711E-02,0.80964802996215E-04,-0.16557679795037E-03,-0.44923899061815E-04
};

void if97_init_region3()
{
    #if 0
    R3I1[1] = 0; R3J1[1] = 0; R3An[1] = 0.10658070028513E+01;
    R3I1[2] = 0; R3J1[2] = 0; R3An[2] = -0.15732845290239E+02;
    R3I1[3] = 0; R3J1[3] = 1; R3An[3] = 0.20944396974307E+02;
    R3I1[4] = 0; R3J1[4] = 2; R3An[4] = -0.76867707878716E+01;
    R3I1[5] = 0; R3J1[5] = 7; R3An[5] = 0.26185947787954E+01;
    R3I1[6] = 0; R3J1[6] = 10; R3An[6] = -0.28080781148620E+01;
    R3I1[7] = 0; R3J1[7] = 12; R3An[7] = 0.12053369696517E+01;
    R3I1[8] = 0; R3J1[8] = 23; R3An[8] = -0.84566812812502E-02;
    R3I1[9] = 1; R3J1[9] = 2; R3An[9] = -0.12654315477714E+01;
    R3I1[10] = 1; R3J1[10] = 6; R3An[10] = -0.11524407806681E+01;
    R3I1[11] = 1; R3J1[11] = 15; R3An[11] = 0.88521043984318E+00;
    R3I1[12] = 1; R3J1[12] = 17; R3An[12] = -0.64207765181607E+00;
    R3I1[13] = 2; R3J1[13] = 0; R3An[13] = 0.38493460186671E+00;
    R3I1[14] = 2; R3J1[14] = 2; R3An[14] = -0.85214708824206E+00;
    R3I1[15] = 2; R3J1[15] = 6; R3An[15] = 0.48972281541877E+01;
    R3I1[16] = 2; R3J1[16] = 7; R3An[16] = -0.30502617256965E+01;
    R3I1[17] = 2; R3J1[17] = 22; R3An[17] = 0.39420536879154E-01;
    R3I1[18] = 2; R3J1[18] = 26; R3An[18] = 0.12558408424308E+00;
    R3I1[19] = 3; R3J1[19] = 0; R3An[19] = -0.27999329698710E+00;
    R3I1[20] = 3; R3J1[20] = 2; R3An[20] = 0.13899799569460E+01;
    R3I1[21] = 3; R3J1[21] = 4; R3An[21] = -0.20189915023570E+01;
    R3I1[22] = 3; R3J1[22] = 16; R3An[22] = -0.82147637173963E-02;
    R3I1[23] = 3; R3J1[23] = 26; R3An[23] = -0.47596035734923E+00;
    R3I1[24] = 4; R3J1[24] = 0; R3An[24] = 0.43984074473500E-01;
    R3I1[25] = 4; R3J1[25] = 2; R3An[25] = -0.44476435428739E+00;
    R3I1[26] = 4; R3J1[26] = 4; R3An[26] = 0.90572070719733E+00;
    R3I1[27] = 4; R3J1[27] = 26; R3An[27] = 0.70522450087967E+00;
    R3I1[28] = 5; R3J1[28] = 1; R3An[28] = 0.10770512626332E+00;
    R3I1[29] = 5; R3J1[29] = 3; R3An[29] = -0.32913623258954E+00;
    R3I1[30] = 5; R3J1[30] = 26; R3An[30] = -0.50871062041158E+00;
    R3I1[31] = 6; R3J1[31] = 0; R3An[31] = -0.22175400873096E-01;
    R3I1[32] = 6; R3J1[32] = 2; R3An[32] = 0.94260751665092E-01;
    R3I1[33] = 6; R3J1[33] = 26; R3An[33] = 0.16436278447961E+00;
    R3I1[34] = 7; R3J1[34] = 2; R3An[34] = -0.13503372241348E-01;
    R3I1[35] = 8; R3J1[35] = 26; R3An[35] = -0.14834345352472E-01;
    R3I1[36] = 9; R3J1[36] = 2; R3An[36] = 0.57922953628084E-03;
    R3I1[37] = 9; R3J1[37] = 26; R3An[37] = 0.32308904703711E-02;
    R3I1[38] = 10; R3J1[38] = 0; R3An[38] = 0.80964802996215E-04;
    R3I1[39] = 10; R3J1[39] = 1; R3An[39] = -0.16557679795037E-03;
    R3I1[40] = 11; R3J1[40] = 26; R3An[40] = -0.44923899061815E-04;
    #endif
}

double if97_region3_delta(double dens, double temK)
{
    return dens / Rhocrt;
}

double if97_region3_tau(double dens, double temK)
{
    return TKcrt / temK;
}

#define USE_OPTIMIZED
#ifdef USE_OPTIMIZED
double if97_region3_phi_tau(double delta, double tau)
{
    double dlt00 = 1.0;
    double dlt01 = delta;
    double dlt02 = dlt01 * dlt01;
    double dlt03 = dlt01 * dlt02;
    double dlt04 = dlt01 * dlt03;
    double dlt05 = dlt01 * dlt04;
    double dlt06 = dlt01 * dlt05;
    double dlt07 = dlt01 * dlt06;
    double dlt08 = dlt01 * dlt07;
    double dlt09 = dlt01 * dlt08;
    double dlt10 = dlt01 * dlt09;
    double dlt11 = dlt01 * dlt10;

    double tau00 = 1.0;
    double tau01 = tau;
    double tau02 = tau01 * tau01;
    double tau03 = tau01 * tau02;
    double tau04 = tau01 * tau03;
    double tau05 = tau01 * tau04;
    double tau06 = tau01 * tau05;
    double tau07 = tau01 * tau06;
    double tau09 = tau02 * tau07;
    double tau10 = tau01 * tau09;
    double tau11 = tau01 * tau10;
    double tau12 = tau01 * tau11;
    double tau14 = tau02 * tau12;
    double tau15 = tau01 * tau14;
    double tau16 = tau01 * tau15;
    double tau17 = tau01 * tau16;
    double tau21 = tau04 * tau17;
    double tau22 = tau01 * tau21;
    double tau23 = tau01 * tau22;
    double tau25 = tau02 * tau23;
    double tau26 = tau01 * tau25;
    double taum1 = 1.0 / tau01;

    double ft = R3An[3]  * dlt00 * R3J1[3]  * tau00 + R3An[4]  * dlt00 * R3J1[4]  * tau01
        + R3An[5]  * dlt00 * R3J1[5]  * tau06 + R3An[6]  * dlt00 * R3J1[6]  * tau09
        + R3An[7]  * dlt00 * R3J1[7]  * tau11 + R3An[8]  * dlt00 * R3J1[8]  * tau22
        + R3An[9]  * dlt01 * R3J1[9]  * tau01 + R3An[10] * dlt01 * R3J1[10] * tau05
        + R3An[11] * dlt01 * R3J1[11] * tau14 + R3An[12] * dlt01 * R3J1[12] * tau16
        + R3An[13] * dlt02 * R3J1[13] * taum1 + R3An[14] * dlt02 * R3J1[14] * tau01
        + R3An[15] * dlt02 * R3J1[15] * tau05 + R3An[16] * dlt02 * R3J1[16] * tau06
        + R3An[17] * dlt02 * R3J1[17] * tau21 + R3An[18] * dlt02 * R3J1[18] * tau25
        + R3An[19] * dlt03 * R3J1[19] * taum1 + R3An[20] * dlt03 * R3J1[20] * tau01
        + R3An[21] * dlt03 * R3J1[21] * tau03 + R3An[22] * dlt03 * R3J1[22] * tau15
        + R3An[23] * dlt03 * R3J1[23] * tau25 + R3An[24] * dlt04 * R3J1[24] * taum1
        + R3An[25] * dlt04 * R3J1[25] * tau01 + R3An[26] * dlt04 * R3J1[26] * tau03
        + R3An[27] * dlt04 * R3J1[27] * tau25 + R3An[28] * dlt05 * R3J1[28] * tau00
        + R3An[29] * dlt05 * R3J1[29] * tau02 + R3An[30] * dlt05 * R3J1[30] * tau25
        + R3An[31] * dlt06 * R3J1[31] * taum1 + R3An[32] * dlt06 * R3J1[32] * tau01
        + R3An[33] * dlt06 * R3J1[33] * tau25 + R3An[34] * dlt07 * R3J1[34] * tau01
        + R3An[35] * dlt08 * R3J1[35] * tau25 + R3An[36] * dlt09 * R3J1[36] * tau01
        + R3An[37] * dlt09 * R3J1[37] * tau25 + R3An[38] * dlt10 * R3J1[38] * taum1
        + R3An[39] * dlt10 * R3J1[39] * tau00 + R3An[40] * dlt11 * R3J1[40] * tau25;

    return ft;
}

double if97_region3_phi_delta(double delta, double tau)
{
    double dlt00 = 1.0;
    double dlt01 = delta;
    double dlt02 = dlt01 * dlt01;
    double dlt03 = dlt01 * dlt02;
    double dlt04 = dlt01 * dlt03;
    double dlt05 = dlt01 * dlt04;
    double dlt06 = dlt01 * dlt05;
    double dlt07 = dlt01 * dlt06;
    double dlt08 = dlt01 * dlt07;
    double dlt09 = dlt01 * dlt08;
    double dlt10 = dlt01 * dlt09;
    double dlt11 = dlt01 * dlt10;

    double tau00 = 1.0;
    double tau01 = tau;
    double tau02 = tau01 * tau01;
    double tau03 = tau01 * tau02;
    double tau04 = tau01 * tau03;
    double tau05 = tau01 * tau04;
    double tau06 = tau01 * tau05;
    double tau07 = tau01 * tau06;
    double tau09 = tau02 * tau07;
    double tau10 = tau01 * tau09;
    double tau11 = tau01 * tau10;
    double tau12 = tau01 * tau11;
    double tau14 = tau02 * tau12;
    double tau15 = tau01 * tau14;
    double tau16 = tau01 * tau15;
    double tau17 = tau01 * tau16;
    double tau21 = tau04 * tau17;
    double tau22 = tau01 * tau21;
    double tau23 = tau01 * tau22;
    double tau25 = tau02 * tau23;
    double tau26 = tau01 * tau25;
    double taum1 = 1.0 / tau01;

    double fd = R3An[1] / dlt01;
    fd += R3An[9]  * R3I1[9]  * dlt00 * tau02 + R3An[10] * R3I1[10] * dlt00 * tau06
        + R3An[11] * R3I1[11] * dlt00 * tau15 + R3An[12] * R3I1[12] * dlt00 * tau17
        + R3An[13] * R3I1[13] * dlt01 * tau00 + R3An[14] * R3I1[14] * dlt01 * tau02
        + R3An[15] * R3I1[15] * dlt01 * tau06 + R3An[16] * R3I1[16] * dlt01 * tau07
        + R3An[17] * R3I1[17] * dlt01 * tau22 + R3An[18] * R3I1[18] * dlt01 * tau26
        + R3An[19] * R3I1[19] * dlt02 * tau00 + R3An[20] * R3I1[20] * dlt02 * tau02
        + R3An[21] * R3I1[21] * dlt02 * tau04 + R3An[22] * R3I1[22] * dlt02 * tau16
        + R3An[23] * R3I1[23] * dlt02 * tau26 + R3An[24] * R3I1[24] * dlt03 * tau00
        + R3An[25] * R3I1[25] * dlt03 * tau02 + R3An[26] * R3I1[26] * dlt03 * tau04
        + R3An[27] * R3I1[27] * dlt03 * tau26 + R3An[28] * R3I1[28] * dlt04 * tau01
        + R3An[29] * R3I1[29] * dlt04 * tau03 + R3An[30] * R3I1[30] * dlt04 * tau26
        + R3An[31] * R3I1[31] * dlt05 * tau00 + R3An[32] * R3I1[32] * dlt05 * tau02
        + R3An[33] * R3I1[33] * dlt05 * tau26 + R3An[34] * R3I1[34] * dlt06 * tau02
        + R3An[35] * R3I1[35] * dlt07 * tau26 + R3An[36] * R3I1[36] * dlt08 * tau02
        + R3An[37] * R3I1[37] * dlt08 * tau26 + R3An[38] * R3I1[38] * dlt09 * tau00
        + R3An[39] * R3I1[39] * dlt09 * tau01 + R3An[40] * R3I1[40] * dlt10 * tau26;

    return fd;
}

#else

double if97_region3_phi_tau(double delta, double tau)
{
    // sum_{i=2}^{40} n*delta^I*J*tau^(J-1)
    double val = 0.0;
    for (int i=2; i<=40; i++)
        val += R3An[i]*pow(delta, R3I1[i])*R3J1[i]*pow(tau, R3J1[i]-1.);
    return val;
}

double if97_region3_phi_delta(double delta, double tau)
{
    // n1/delta + sum_{i=2}^{40} n*I*delta^(I-1)*tau^J
    double val = R3An[1]/delta;
    for (int i=2; i<=40; i++)
        val += R3An[i]*R3I1[i]*pow(delta, R3I1[i]-1.)*pow(tau, R3J1[i]);
    return val;
}

#endif

double if97_region3_phi_tautau(double delta, double tau)
{
    // sum_{i=2}^{40} n*delta^I*J*(J-1)*tau^(J-2)
    double val = 0.0;
    for (int i=2; i<=40; i++)
    {
        val += R3An[i]*pow(delta, R3I1[i])*R3J1[i]*(R3J1[i]-1.)*pow(tau, R3J1[i]-2.);
    }
    return val;
}

double if97_region3_phi_deltadelta(double delta, double tau)
{
    // -n1/delta^2 + sum_{i=2}^{40} n*I*(I-1)*delta^(I-2)*tau^J
    double val = - R3An[1]/(delta*delta);
    for (int i=2; i<=40; i++)
    {
        val += R3An[i]*R3I1[i]*(R3I1[i]-1.)*pow(delta, R3I1[i]-2.)*pow(tau, R3J1[i]);
    }
    return val;
}

double if97_region3_phi_deltatau(double delta, double tau)
{
    // sum_{i=2}^{40} n*I*delta^(I-1)*J*tau^(J-1)
    double val = 0.0;
    for (int i=2; i<=40; i++)
    {
        val += R3An[i]*R3I1[i]*pow(delta, R3I1[i]-1.)*R3J1[i]*pow(tau, R3J1[i]-1.);
    }
    return val;
}


double if97_h_rhoT_region3(double dens, double temK)
{
    double delta = if97_region3_delta(dens, temK);
    double tau = if97_region3_tau(dens, temK);
    double phit = if97_region3_phi_tau(delta, tau);
    double phid = if97_region3_phi_delta(delta, tau);
    return (tau * phit + delta * phid) * ConstR * temK;
}

double if97_cp_rhoT_region3(double dens, double temK)
{
    // double cp = -tau^2*phi_tt+(delta*phid-deulta*tau*phi_dt)^2/(2*delta*phid+delta^2*phi_dd)
    double delta = if97_region3_delta(dens, temK);
    double tau = if97_region3_tau(dens, temK);
    double phid = if97_region3_phi_delta(delta, tau);
    double phitt = if97_region3_phi_tautau(delta, tau);
    double phidd = if97_region3_phi_deltadelta(delta, tau);
    double phidt = if97_region3_phi_deltatau(delta, tau);
    double cp_R = -tau*tau*phitt + pow(delta*phid - delta*tau*phidt, 2.)/(2.*delta*phid + delta*delta*phidd);
    return cp_R * ConstR;
}

double if97_cv_rhoT_region3(double dens, double temK)
{
    double delta = if97_region3_delta(dens, temK);
    double tau = if97_region3_tau(dens, temK);
    double phitt = if97_region3_phi_tautau(delta, tau);
    double cv = - tau*tau*phitt*ConstR;
    return cv;
}

double if97_p_vT_region3(double volm, double temK)
{
    double dens = 1./volm;
    double delta = if97_region3_delta(dens, temK);
    double tau = if97_region3_tau(dens, temK);
    double phid = if97_region3_phi_delta(delta, tau);
    double p = delta * phid * dens * ConstR * temK;
    return p * 1.e-3;
}


double if97_dpdv_vT_region3(double volm, double temK)
{
    double dltx01 = 1.0 / volm / Rhocrt;
    double dltx02 = dltx01 * dltx01;
    double dltx03 = dltx01 * dltx02;
    double dltx04 = dltx01 * dltx03;
    double dltx05 = dltx01 * dltx04;
    double dltx06 = dltx01 * dltx05;
    double dltx07 = dltx01 * dltx06;
    double dltx08 = dltx01 * dltx07;
    double dltx09 = dltx01 * dltx08;
    double dltx10 = dltx01 * dltx09;
    double dltx11 = dltx01 * dltx10;

    double taux00 = 1.0;
    double taux01 = TKcrt / temK;
    double taux02 = taux01 * taux01;
    double taux03 = taux01 * taux02;
    double taux04 = taux01 * taux03;
    double taux06 = taux02 * taux04;
    double taux07 = taux01 * taux06;
    double taux15 = taux01 * taux07 * taux07;
    double taux16 = taux01 * taux15;
    double taux17 = taux01 * taux16;
    double taux22 = taux06 * taux16;
    double taux26 = taux04 * taux22;

    double dPdV = R3An[1]
        + R3An[9]  * R3I1[9]  * dltx01 * taux02
        + R3An[10] * R3I1[10] * dltx01 * taux06
        + R3An[11] * R3I1[11] * dltx01 * taux15
        + R3An[12] * R3I1[12] * dltx01 * taux17
        + R3An[13] * R3I1[13] * dltx02 * taux00
        + R3An[14] * R3I1[14] * dltx02 * taux02
        + R3An[15] * R3I1[15] * dltx02 * taux06
        + R3An[16] * R3I1[16] * dltx02 * taux07
        + R3An[17] * R3I1[17] * dltx02 * taux22
        + R3An[18] * R3I1[18] * dltx02 * taux26
        + R3An[19] * R3I1[19] * dltx03 * taux00
        + R3An[20] * R3I1[20] * dltx03 * taux02
        + R3An[21] * R3I1[21] * dltx03 * taux04
        + R3An[22] * R3I1[22] * dltx03 * taux16
        + R3An[23] * R3I1[23] * dltx03 * taux26
        + R3An[24] * R3I1[24] * dltx04 * taux00
        + R3An[25] * R3I1[25] * dltx04 * taux02
        + R3An[26] * R3I1[26] * dltx04 * taux04
        + R3An[27] * R3I1[27] * dltx04 * taux26
        + R3An[28] * R3I1[28] * dltx05 * taux01
        + R3An[29] * R3I1[29] * dltx05 * taux03
        + R3An[30] * R3I1[30] * dltx05 * taux26
        + R3An[31] * R3I1[31] * dltx06 * taux00
        + R3An[32] * R3I1[32] * dltx06 * taux02
        + R3An[33] * R3I1[33] * dltx06 * taux26
        + R3An[34] * R3I1[34] * dltx07 * taux02
        + R3An[35] * R3I1[35] * dltx08 * taux26
        + R3An[36] * R3I1[36] * dltx09 * taux02
        + R3An[37] * R3I1[37] * dltx09 * taux26
        + R3An[38] * R3I1[38] * dltx10 * taux00
        + R3An[39] * R3I1[39] * dltx10 * taux01
        + R3An[40] * R3I1[40] * dltx11 * taux26;

    dPdV *= -temK * 1.0E-3 * ConstR / volm / volm;

    return dPdV;
}


double if97_spin_vL_T_region3(double temK)
{
    double vMin = 1.3E-3, vCrt = 1.0 / Rhocrt;
    double v1 = vCrt;
    double v2 = vMin;

    double vm = 0.5 * (v1 + v2);
    for (int n = 1; n <= 30; n++)
    {
        double dPdV = if97_dpdv_vT_region3(vm, temK);

        if (dPdV >= 0.0)
            v1 = vm;
        else
            v2 = vm;
        vm = 0.5 * (v1 + v2);
        if (fabs(v1-v2)<vm*1e-8)
            return vm;
    }
    assert(false);
    return -1;
}


double if97_spin_vG_T_region3(double temK)
{
    double vCrt = 1.0 / Rhocrt, vMax = 8.9E-3;
    double v1 = vMax;
    double v2 = vCrt;

    double vm = 0.5 * (v1 + v2);
    for (int n = 1; n <= 30; n++)
    {
        double dPdV = if97_dpdv_vT_region3(vm, temK);

        if (dPdV <= 0.0)
            v1 = vm;
        else
            v2 = vm;
        vm = 0.5 * (v1 + v2);
        if (fabs(v1-v2)<vm*1e-8)
            return vm;
    }

    assert(false);
    return -1;
}


double if97_rho_pT_region3(double pres, double temK)
{
    const double vMin = 1.3E-3, vMax = 8.9E-3;

    double v1 = vMax;
    double v2 = vMin;
    if (temK < TKcrt)
    {
        const double psat = if97_sat_p_T(temK);
        if (pres <= psat) /**********  VPT_32(P_in, T_in, V_out)  **********/
        {
            v1 = vMax;
            v2 = if97_spin_vG_T_region3(temK);
        }
        else /**********  VPT_31(P_in, T_in, V_out)  **********/
        {
            v1 = if97_spin_vL_T_region3(temK);
            v2 = vMin;
        }
    }

#if 1
    for (int n = 0; n < 5; n++)
    {
        double vm = 0.5 * (v1 + v2);
        double pre0 = if97_p_vT_region3(vm, temK);
        if (pre0 <= pres)
            v1 = vm;
        else
            v2 = vm;
    }

    double vm = 0.5 * (v1 + v2);
    const double dv = vm * 1e-3;
    for (int i=0; i<100; i++)
    {
        double p1 = if97_p_vT_region3(vm, temK);
        double r = p1 - pres;
        if (fabs(r) < pres*1e-8)
            return 1/vm;
        double p2 = if97_p_vT_region3(vm + dv, temK);
        double dr_dv = (p2-p1)/dv;
        vm -= r/dr_dv;
    }
    assert(false);
    return -1;
#else
    for (int n = 0; n < 40; n++)
    {
        double vm = 0.5 * (v1 + v2);
        double pre0 = if97_p_vT_region3(vm, temK);
        if (pre0 <= pres)
            v1 = vm;
        else
            v2 = vm;
    }
    double volume = 0.5 * (v1 + v2);
    double rho = 1./volume;
    return rho;
#endif
}


double if97_sat_vL_T_region3(double temK)
{
    double vMin = 1.3E-3;
    double pSat = if97_sat_p_T(temK);
    double volm = if97_spin_vL_T_region3(temK);

    double v1 = volm;
    double v2 = vMin;

    double vm = 0.5 * (v1 + v2);
    for (int n = 1; n <= 40; n++)
    {
        double pres = if97_p_vT_region3(vm, temK);

        if (pres <= pSat)
            v1 = vm;
        else
            v2 = vm;
        vm = 0.5 * (v1 + v2);
        if (fabs(v1-v2)<vm*1e-8)
            return vm;
    }

    assert(false);
    return -1;
}

double if97_sat_vG_T_region3(double temK)
{
    double vMax = 8.9E-3;
    double pSat = if97_sat_p_T(temK);
    double volm = if97_spin_vG_T_region3(temK);

    double v1 = vMax;
    double v2 = volm;

    double vm = 0.5 * (v1 + v2);
    for (int n = 1; n <= 40; n++)
    {
        double pres = if97_p_vT_region3(vm, temK);

        if (pres <= pSat)
            v1 = vm;
        else
            v2 = vm;
        vm = 0.5 * (v1 + v2);
        if (fabs(v1-v2)<vm*1e-8)
            return vm;
    }

    assert(false);
    return -1;
}



double if97_T_ph_region30(double pres, double enth)
{
#if 1
    double t1 = if97_boundary23_T_p(pres);
    double t2 = TKref;

    for (int i = 1; i <= 5; i++)
    {
        double tm = 0.5 * (t1 + t2);
        double dd = if97_rho_pT_region3(pres, tm);
        double hh = if97_h_rhoT_region3(dd, tm);

        if (hh >= enth)
            t1 = tm;
        else
            t2 = tm;
    }

    double tm = 0.5 * (t1 + t2);
    double dt = tm * 1e-3;
    for (int i=0; i<100; i++)
    {
        double h1 = if97_h_rhoT_region3(if97_rho_pT_region3(pres, tm), tm);
        double res = h1 - enth;
        if (fabs(res) < enth*1e-8)
            return tm;
        double h2 = if97_h_rhoT_region3(if97_rho_pT_region3(pres, tm+dt), tm+dt);
        double drdt = (h2-h1)/dt;
        tm -= res/drdt;
    }
    assert(false);
    return -1;

#else
    double t1 = if97_boundary23_T_p(pres);
    double t2 = TKref;

    for (int i = 1; i <= 40; i++)
    {
        double tm = 0.5 * (t1 + t2);
        double dd = if97_rho_pT_region3(pres, tm);
        double hh = if97_h_rhoT_region3(dd, tm);

        if (hh >= enth)
            t1 = tm;
        else
            t2 = tm;
    }

    return 0.5 * (t1 + t2);
#endif
}


double if97_T_ph_region31(double pres, double enth)
{
    double t1 = if97_sat_T_p(pres);
    double t2 = TKref;

    double tm = 0.5 * (t1 + t2);
    for (int i = 1; i <= 40; i++)
    {
        double dd = if97_rho_pT_region3(pres, tm);
        double hh = if97_h_rhoT_region3(dd, tm);

        if (hh >= enth)
            t1 = tm;
        else
            t2 = tm;
        tm = 0.5 * (t1 + t2);
        if (fabs(t1-t2)<tm*1e-8)
            return tm;
    }

    assert(false);
    return -1;
}

double if97_T_ph_region32(double pres, double enth)
{
    double t1 = if97_boundary23_T_p(pres);
    double t2 = if97_sat_T_p(pres);

    double tm = 0.5 * (t1 + t2);
    for (int i = 1; i <= 40; i++)
    {
        double dd = if97_rho_pT_region3(pres, tm);
        double hh = if97_h_rhoT_region3(dd, tm);

        if (hh >= enth)
            t1 = tm;
        else
            t2 = tm;
        tm = 0.5 * (t1 + t2);
        if (fabs(t1-t2)<tm*1e-8)
            return tm;
    }

    assert(false);
    return -1;
}


double if97_h_pT_region3(double pres, double temK)
{
    double dd = if97_rho_pT_region3(pres, temK);
    return if97_h_rhoT_region3(dd, temK);
}

