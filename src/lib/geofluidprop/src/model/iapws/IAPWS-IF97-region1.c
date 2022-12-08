
#include <math.h>
#include <stdio.h>

#include "IAPWS-IF97-const.h"

double R1I1[35], R1J1[35], R1An[35];
double R1I3[21], R1J3[21], R1Cn[21];


void if97_init_region1()
{
    R1I1[1] = 0; R1J1[1] = -2; R1An[1] = 0.14632971213167E+00;
    R1I1[2] = 0; R1J1[2] = -1; R1An[2] = -0.84548187169114E+00;
    R1I1[3] = 0; R1J1[3] = 0; R1An[3] = -0.37563603672040E+01;
    R1I1[4] = 0; R1J1[4] = 1; R1An[4] = 0.33855169168385E+01;
    R1I1[5] = 0; R1J1[5] = 2; R1An[5] = -0.95791963387872E+00;
    R1I1[6] = 0; R1J1[6] = 3; R1An[6] = 0.15772038513228E+00;
    R1I1[7] = 0; R1J1[7] = 4; R1An[7] = -0.16616417199501E-01;
    R1I1[8] = 0; R1J1[8] = 5; R1An[8] = 0.81214629983568E-03;
    R1I1[9] = 1; R1J1[9] = -9; R1An[9] = 0.28319080123804E-03;
    R1I1[10] = 1; R1J1[10] = -7; R1An[10] = -0.60706301565874E-03;
    R1I1[11] = 1; R1J1[11] = -1; R1An[11] = -0.18990068218419E-01;
    R1I1[12] = 1; R1J1[12] = 0; R1An[12] = -0.32529748770505E-01;
    R1I1[13] = 1; R1J1[13] = 1; R1An[13] = -0.21841717175414E-01;
    R1I1[14] = 1; R1J1[14] = 3; R1An[14] = -0.52838357969930E-04;
    R1I1[15] = 2; R1J1[15] = -3; R1An[15] = -0.47184321073267E-03;
    R1I1[16] = 2; R1J1[16] = 0; R1An[16] = -0.30001780793026E-03;
    R1I1[17] = 2; R1J1[17] = 1; R1An[17] = 0.47661393906987E-04;
    R1I1[18] = 2; R1J1[18] = 3; R1An[18] = -0.44141845330846E-05;
    R1I1[19] = 2; R1J1[19] = 17; R1An[19] = -0.72694996297594E-15;
    R1I1[20] = 3; R1J1[20] = -4; R1An[20] = -0.31679644845054E-04;
    R1I1[21] = 3; R1J1[21] = 0; R1An[21] = -0.28270797985312E-05;
    R1I1[22] = 3; R1J1[22] = 6; R1An[22] = -0.85205128120103E-09;
    R1I1[23] = 4; R1J1[23] = -5; R1An[23] = -0.22425281908000E-05;
    R1I1[24] = 4; R1J1[24] = -2; R1An[24] = -0.65171222895601E-06;
    R1I1[25] = 4; R1J1[25] = 10; R1An[25] = -0.14341729937924E-12;
    R1I1[26] = 5; R1J1[26] = -8; R1An[26] = -0.40516996860117E-06;
    R1I1[27] = 8; R1J1[27] = -11; R1An[27] = -0.12734301741641E-08;
    R1I1[28] = 8; R1J1[28] = -6; R1An[28] = -0.17424871230634E-09;
    R1I1[29] = 21; R1J1[29] = -29; R1An[29] = -0.68762131295531E-18;
    R1I1[30] = 23; R1J1[30] = -31; R1An[30] = 0.14478307828521E-19;
    R1I1[31] = 29; R1J1[31] = -38; R1An[31] = 0.26335781662795E-22;
    R1I1[32] = 30; R1J1[32] = -39; R1An[32] = -0.11947622640071E-22;
    R1I1[33] = 31; R1J1[33] = -40; R1An[33] = 0.18228094581404E-23;
    R1I1[34] = 32; R1J1[34] = -41; R1An[34] = -0.93537087292458E-25;

    R1I3[1] = 0; R1J3[1] = 0; R1Cn[1] = -0.23872489924521E+03;
    R1I3[2] = 0; R1J3[2] = 1; R1Cn[2] = 0.40421188637945E+03;
    R1I3[3] = 0; R1J3[3] = 2; R1Cn[3] = 0.11349746881718E+03;
    R1I3[4] = 0; R1J3[4] = 6; R1Cn[4] = -0.58457616048039E+01;
    R1I3[5] = 0; R1J3[5] = 22; R1Cn[5] = -0.15285482413140E-03;
    R1I3[6] = 0; R1J3[6] = 32; R1Cn[6] = -0.10866707695377E-05;
    R1I3[7] = 1; R1J3[7] = 0; R1Cn[7] = -0.13391744872602E+02;
    R1I3[8] = 1; R1J3[8] = 1; R1Cn[8] = 0.43211039183559E+02;
    R1I3[9] = 1; R1J3[9] = 2; R1Cn[9] = -0.54010067170506E+02;
    R1I3[10] = 1; R1J3[10] = 3; R1Cn[10] = 0.30535892203916E+02;
    R1I3[11] = 1; R1J3[11] = 4; R1Cn[11] = -0.65964749423638E+01;
    R1I3[12] = 1; R1J3[12] = 10; R1Cn[12] = 0.93965400878363E-02;
    R1I3[13] = 1; R1J3[13] = 32; R1Cn[13] = 0.11573647505340E-06;
    R1I3[14] = 2; R1J3[14] = 10; R1Cn[14] = -0.25858641282073E-04;
    R1I3[15] = 2; R1J3[15] = 32; R1Cn[15] = -0.40644363084799E-08;
    R1I3[16] = 3; R1J3[16] = 10; R1Cn[16] = 0.66456186191635E-07;
    R1I3[17] = 3; R1J3[17] = 32; R1Cn[17] = 0.80670734103027E-10;
    R1I3[18] = 4; R1J3[18] = 32; R1Cn[18] = -0.93477771213947E-12;
    R1I3[19] = 5; R1J3[19] = 32; R1Cn[19] = 0.58265442020601E-14;
    R1I3[20] = 6; R1J3[20] = 32; R1Cn[20] = -0.15020185953503E-16;
}


double if97_region1_pi(double pres)
{
    return pres / 16.53;
}

double if97_region1_tau(double temK)
{
    return 1386.0 / temK;
}

double if97_region1_gamma_tau(double pi, double tau)
{
    double pn000 = 1.0;
    double pn001 = 7.1 - pi;
    double pn002 = pn001 * pn001;
    double pn003 = pn001 * pn002;
    double pn004 = pn001 * pn003;
    double pn005 = pn001 * pn004;
    double pn008 = pn003 * pn005;
    double pn021 = pn005 * pn008 * pn008;
    double pn023 = pn002 * pn021;
    double pn029 = pn001 * pn005 * pn023;
    double pn030 = pn001 * pn029;
    double pn031 = pn001 * pn030;
    double pn032 = pn001 * pn031;

    double tn000 = 1.0;
    double tn001 = tau - 1.222;
    double tn002 = tn001 * tn001;
    double tn003 = tn001 * tn002;
    double tn004 = tn001 * tn003;
    double tn005 = tn001 * tn004;
    double tn009 = tn004 * tn005;
    double tn016 = tn002 * tn005 * tn009;
    double tnm01 = 1.0 / tn001;
    double tnm02 = tnm01 * tnm01;
    double tnm03 = tnm01 * tnm02;
    double tnm04 = tnm01 * tnm03;
    double tnm05 = tnm01 * tnm04;
    double tnm06 = tnm01 * tnm05;
    double tnm07 = tnm01 * tnm06;
    double tnm08 = tnm01 * tnm07;
    double tnm09 = tnm01 * tnm08;
    double tnm10 = tnm01 * tnm09;
    double tnm12 = tnm02 * tnm10;
    double tnm30 = tnm08 * tnm10 * tnm12;
    double tnm32 = tnm02 * tnm30;
    double tnm39 = tnm07 * tnm32;
    double tnm40 = tnm01 * tnm39;
    double tnm41 = tnm01 * tnm40;
    double tnm42 = tnm01 * tnm41;

    double gt  = R1An[1]  * pn000 * R1J1[1]  * tnm03 + R1An[2]  * pn000 * R1J1[2]  * tnm02
        + R1An[3]  * pn000 * R1J1[3]  * tnm01 + R1An[4]  * pn000 * R1J1[4]  * tn000
        + R1An[5]  * pn000 * R1J1[5]  * tn001 + R1An[6]  * pn000 * R1J1[6]  * tn002
        + R1An[7]  * pn000 * R1J1[7]  * tn003 + R1An[8]  * pn000 * R1J1[8]  * tn004
        + R1An[9]  * pn001 * R1J1[9]  * tnm10 + R1An[10] * pn001 * R1J1[10] * tnm08
        + R1An[11] * pn001 * R1J1[11] * tnm02 + R1An[12] * pn001 * R1J1[12] * tnm01
        + R1An[13] * pn001 * R1J1[13] * tn000 + R1An[14] * pn001 * R1J1[14] * tn002
        + R1An[15] * pn002 * R1J1[15] * tnm04 + R1An[16] * pn002 * R1J1[16] * tnm01
        + R1An[17] * pn002 * R1J1[17] * tn000 + R1An[18] * pn002 * R1J1[18] * tn002
        + R1An[19] * pn002 * R1J1[19] * tn016 + R1An[20] * pn003 * R1J1[20] * tnm05
        + R1An[21] * pn003 * R1J1[21] * tnm01 + R1An[22] * pn003 * R1J1[22] * tn005
        + R1An[23] * pn004 * R1J1[23] * tnm06 + R1An[24] * pn004 * R1J1[24] * tnm03
        + R1An[25] * pn004 * R1J1[25] * tn009 + R1An[26] * pn005 * R1J1[26] * tnm09
        + R1An[27] * pn008 * R1J1[27] * tnm12 + R1An[28] * pn008 * R1J1[28] * tnm07
        + R1An[29] * pn021 * R1J1[29] * tnm30 + R1An[30] * pn023 * R1J1[30] * tnm32
        + R1An[31] * pn029 * R1J1[31] * tnm39 + R1An[32] * pn030 * R1J1[32] * tnm40
        + R1An[33] * pn031 * R1J1[33] * tnm41 + R1An[34] * pn032 * R1J1[34] * tnm42;

    return gt;
}

double if97_region1_gamma_pi(double pi, double tau)
{
    double pn00 = 1.0;
    double pn01 = 7.1 - pi;
    double pn02 = pn01 * pn01;
    double pn03 = pn01 * pn02;
    double pn04 = pn01 * pn03;
    double pn07 = pn03 * pn04;
    double pn20 = pn02 * pn04 * pn07 * pn07;
    double pn22 = pn02 * pn20;
    double pn28 = pn02 * pn04 * pn22;
    double pn29 = pn01 * pn28;
    double pn30 = pn01 * pn29;
    double pn31 = pn01 * pn30;

    double tn000 = 1.0;
    double tn001 = tau - 1.222;
    double tn003 = tn001 * tn001 * tn001;
    double tn006 = tn003 * tn003;
    double tn010 = tn001 * tn003 * tn006;
    double tn017 = tn001 * tn006 * tn010;
    double tnm01 = 1.0 / tn001;
    double tnm02 = tnm01 * tnm01;
    double tnm03 = tnm01 * tnm02;
    double tnm04 = tnm01 * tnm03;
    double tnm05 = tnm01 * tnm04;
    double tnm06 = tnm01 * tnm05;
    double tnm07 = tnm01 * tnm06;
    double tnm08 = tnm01 * tnm07;
    double tnm09 = tnm01 * tnm08;
    double tnm11 = tnm02 * tnm09;
    double tnm29 = tnm07 * tnm11 * tnm11;
    double tnm31 = tnm02 * tnm29;
    double tnm38 = tnm07 * tnm31;
    double tnm39 = tnm01 * tnm38;
    double tnm40 = tnm01 * tnm39;
    double tnm41 = tnm01 * tnm40;

    double gp = -R1An[9]  * R1I1[9]  * pn00 * tnm09 - R1An[10] * R1I1[10] * pn00 * tnm07
        - R1An[11] * R1I1[11] * pn00 * tnm01 - R1An[12] * R1I1[12] * pn00 * tn000
        - R1An[13] * R1I1[13] * pn00 * tn001 - R1An[14] * R1I1[14] * pn00 * tn003
        - R1An[15] * R1I1[15] * pn01 * tnm03 - R1An[16] * R1I1[16] * pn01 * tn000
        - R1An[17] * R1I1[17] * pn01 * tn001 - R1An[18] * R1I1[18] * pn01 * tn003
        - R1An[19] * R1I1[19] * pn01 * tn017 - R1An[20] * R1I1[20] * pn02 * tnm04
        - R1An[21] * R1I1[21] * pn02 * tn000 - R1An[22] * R1I1[22] * pn02 * tn006
        - R1An[23] * R1I1[23] * pn03 * tnm05 - R1An[24] * R1I1[24] * pn03 * tnm02
        - R1An[25] * R1I1[25] * pn03 * tn010 - R1An[26] * R1I1[26] * pn04 * tnm08
        - R1An[27] * R1I1[27] * pn07 * tnm11 - R1An[28] * R1I1[28] * pn07 * tnm06
        - R1An[29] * R1I1[29] * pn20 * tnm29 - R1An[30] * R1I1[30] * pn22 * tnm31
        - R1An[31] * R1I1[31] * pn28 * tnm38 - R1An[32] * R1I1[32] * pn29 * tnm39
        - R1An[33] * R1I1[33] * pn30 * tnm40 - R1An[34] * R1I1[34] * pn31 * tnm41;
    return gp;
}

double if97_region1_gamma_tautau(double pi, double tau)
{
    double pn00 = 1.0;
    double pn01 = 7.1 - pi;
    double pn02 = pn01 * pn01;
    double pn03 = pn01 * pn02;
    double pn04 = pn01 * pn03;
    double pn05 = pn01 * pn04;
    double pn08 = pn03 * pn05;
    double pn21 = pn05 * pn08 * pn08;
    double pn23 = pn02 * pn21;
    double pn29 = pn01 * pn05 * pn23;
    double pn30 = pn01 * pn29;
    double pn31 = pn01 * pn30;
    double pn32 = pn01 * pn31;

    double tn000 = 1.0;
    double tn001 = tau - 1.222;
    double tn002 = tn001 * tn001;
    double tn003 = tn001 * tn002;
    double tn004 = tn001 * tn003;
    double tn008 = tn004 * tn004;
    double tn015 = tn003 * tn004 * tn008;
    double tnm01 = 1.0 / tn001;
    double tnm02 = tnm01 * tnm01;
    double tnm03 = tnm01 * tnm02;
    double tnm04 = tnm01 * tnm03;
    double tnm05 = tnm01 * tnm04;
    double tnm06 = tnm01 * tnm05;
    double tnm07 = tnm01 * tnm06;
    double tnm08 = tnm01 * tnm07;
    double tnm09 = tnm01 * tnm08;
    double tnm10 = tnm01 * tnm09;
    double tnm11 = tnm01 * tnm10;
    double tnm13 = tnm02 * tnm11;
    double tnm31 = tnm05 * tnm13 * tnm13;
    double tnm33 = tnm02 * tnm31;
    double tnm40 = tnm07 * tnm33;
    double tnm41 = tnm01 * tnm40;
    double tnm42 = tnm01 * tnm41;
    double tnm43 = tnm01 * tnm42;

    double gtt = R1An[1]  * pn00 * R1J1[1]  * (R1J1[1]  - 1.0) * tnm04
        + R1An[2]  * pn00 * R1J1[2]  * (R1J1[2]  - 1.0) * tnm03
        + R1An[3]  * pn00 * R1J1[3]  * (R1J1[3]  - 1.0) * tnm02
        + R1An[4]  * pn00 * R1J1[4]  * (R1J1[4]  - 1.0) * tnm01
        + R1An[5]  * pn00 * R1J1[5]  * (R1J1[5]  - 1.0) * tn000
        + R1An[6]  * pn00 * R1J1[6]  * (R1J1[6]  - 1.0) * tn001
        + R1An[7]  * pn00 * R1J1[7]  * (R1J1[7]  - 1.0) * tn002
        + R1An[8]  * pn00 * R1J1[8]  * (R1J1[8]  - 1.0) * tn003
        + R1An[9]  * pn01 * R1J1[9]  * (R1J1[9]  - 1.0) * tnm11
        + R1An[10] * pn01 * R1J1[10] * (R1J1[10] - 1.0) * tnm09
        + R1An[11] * pn01 * R1J1[11] * (R1J1[11] - 1.0) * tnm03
        + R1An[12] * pn01 * R1J1[12] * (R1J1[12] - 1.0) * tnm02
        + R1An[13] * pn01 * R1J1[13] * (R1J1[13] - 1.0) * tnm01
        + R1An[14] * pn01 * R1J1[14] * (R1J1[14] - 1.0) * tn001
        + R1An[15] * pn02 * R1J1[15] * (R1J1[15] - 1.0) * tnm05
        + R1An[16] * pn02 * R1J1[16] * (R1J1[16] - 1.0) * tnm02
        + R1An[17] * pn02 * R1J1[17] * (R1J1[17] - 1.0) * tnm01
        + R1An[18] * pn02 * R1J1[18] * (R1J1[18] - 1.0) * tn001
        + R1An[19] * pn02 * R1J1[19] * (R1J1[19] - 1.0) * tn015
        + R1An[20] * pn03 * R1J1[20] * (R1J1[20] - 1.0) * tnm06
        + R1An[21] * pn03 * R1J1[21] * (R1J1[21] - 1.0) * tnm02
        + R1An[22] * pn03 * R1J1[22] * (R1J1[22] - 1.0) * tn004
        + R1An[23] * pn04 * R1J1[23] * (R1J1[23] - 1.0) * tnm07
        + R1An[24] * pn04 * R1J1[24] * (R1J1[24] - 1.0) * tnm04
        + R1An[25] * pn04 * R1J1[25] * (R1J1[25] - 1.0) * tn008
        + R1An[26] * pn05 * R1J1[26] * (R1J1[26] - 1.0) * tnm10
        + R1An[27] * pn08 * R1J1[27] * (R1J1[27] - 1.0) * tnm13
        + R1An[28] * pn08 * R1J1[28] * (R1J1[28] - 1.0) * tnm08
        + R1An[29] * pn21 * R1J1[29] * (R1J1[29] - 1.0) * tnm31
        + R1An[30] * pn23 * R1J1[30] * (R1J1[30] - 1.0) * tnm33
        + R1An[31] * pn29 * R1J1[31] * (R1J1[31] - 1.0) * tnm40
        + R1An[32] * pn30 * R1J1[32] * (R1J1[32] - 1.0) * tnm41
        + R1An[33] * pn31 * R1J1[33] * (R1J1[33] - 1.0) * tnm42
        + R1An[34] * pn32 * R1J1[34] * (R1J1[34] - 1.0) * tnm43;
    return gtt;
}

double if97_h_pT_region1(double pres, double temK)
{
    double pi = if97_region1_pi(pres);
    double tau = if97_region1_tau(temK);
    double gt  = if97_region1_gamma_tau(pi, tau);
    return tau * gt * ConstR * temK;
}


double if97_rho_pT_region1(double pres, double temK)
{
    double pi = if97_region1_pi(pres);
    double tau = if97_region1_tau(temK);
    double gp  = if97_region1_gamma_pi(pi, tau);
    return 1.0E+3 * pres / (pi * gp * ConstR * temK);
}


double if97_cp_pT_region1(double pres, double temK)
{
    double pi = if97_region1_pi(pres);
    double tau = if97_region1_tau(temK);
    double gtt = if97_region1_gamma_tautau(pi, tau);
    return - tau * tau * gtt * ConstR;
}

double if97_T_ph_region1(double pres, double enth)
{
    double pn0 = 1.0;
    double pn1 = pres;
    double pn2 = pn1 * pn1;
    double pn3 = pn1 * pn2;
    double pn4 = pn1 * pn3;
    double pn5 = pn1 * pn4;
    double pn6 = pn1 * pn5;

    double eta00 = 1.0;
    double eta01 = enth / 2500.0 + 1.0;
    double eta02 = eta01 * eta01;
    double eta03 = eta01 * eta02;
    double eta04 = eta01 * eta03;
    double eta06 = eta02 * eta04;
    double eta10 = eta04 * eta06;
    double eta22 = eta02 * eta10 * eta10;
    double eta32 = eta10 * eta22;

    double tt = R1Cn[1]   * pn0 * eta00
        + R1Cn[2]  * pn0 * eta01
        + R1Cn[3]  * pn0 * eta02
        + R1Cn[4]  * pn0 * eta06
        + R1Cn[5]  * pn0 * eta22
        + R1Cn[6]  * pn0 * eta32
        + R1Cn[7]  * pn1 * eta00
        + R1Cn[8]  * pn1 * eta01
        + R1Cn[9]  * pn1 * eta02
        + R1Cn[10] * pn1 * eta03
        + R1Cn[11] * pn1 * eta04
        + R1Cn[12] * pn1 * eta10
        + R1Cn[13] * pn1 * eta32
        + R1Cn[14] * pn2 * eta10
        + R1Cn[15] * pn2 * eta32
        + R1Cn[16] * pn3 * eta10
        + R1Cn[17] * pn3 * eta32
        + R1Cn[18] * pn4 * eta32
        + R1Cn[19] * pn5 * eta32
        + R1Cn[20] * pn6 * eta32;

    return tt;
}

// JSME Adjustment
double if97_T_ph_region1_jsme(double pres, double enth)
{
    const double eps = 1.0E-6;
    int converged = 0;
    double tt = if97_T_ph_region1(pres, enth);
    for (int i = 1; i <= 10; i++)
    {
        double hh = if97_h_pT_region1(pres, tt);
        double dlt = enth - hh;
        if (fabs(dlt) <= eps)
        {
            converged = 1;
            break;
        }

        double cp = if97_cp_pT_region1(pres, tt);
        tt += dlt / cp;
    }

    if (!converged)
    {
        printf("%s not converged! pres = %g, h = %g\n", __FUNCTION__, pres, enth);
        return -1;
    }

    return tt;
}
