
#include <math.h>
#include <stdio.h>

#include "IAPWS-IF97-const.h"


double R2J1[10], R2An[10];
double R2I2[44], R2J2[44], R2Bn[44];
double R2I3[35], R2J3[35], R2A3[35];
double R2I4[39], R2J4[39], R2A4[39];
double R2I5[24], R2J5[24], R2A5[24];


void if97_init_region2()
{
    R2J1[1] = 0; R2An[1] = -0.96927686500217E+01;
    R2J1[2] = 1; R2An[2] = 0.10086655968018E+02;
    R2J1[3] = -5; R2An[3] = -0.56087911283020E-02;
    R2J1[4] = -4; R2An[4] = 0.71452738081455E-01;
    R2J1[5] = -3; R2An[5] = -0.40710498223928E+00;
    R2J1[6] = -2; R2An[6] = 0.14240819171444E+01;
    R2J1[7] = -1; R2An[7] = -0.43839511319450E+01;
    R2J1[8] = 2; R2An[8] = -0.28408632460772E+00;
    R2J1[9] = 3; R2An[9] = 0.21268463753307E-01;

    R2I2[1] = 1; R2J2[1] = 0; R2Bn[1] = -0.17731742473213E-02;
    R2I2[2] = 1; R2J2[2] = 1; R2Bn[2] = -0.17834862292358E-01;
    R2I2[3] = 1; R2J2[3] = 2; R2Bn[3] = -0.45996013696365E-01;
    R2I2[4] = 1; R2J2[4] = 3; R2Bn[4] = -0.57581259083432E-01;
    R2I2[5] = 1; R2J2[5] = 6; R2Bn[5] = -0.50325278727930E-01;
    R2I2[6] = 2; R2J2[6] = 1; R2Bn[6] = -0.33032641670203E-04;
    R2I2[7] = 2; R2J2[7] = 2; R2Bn[7] = -0.18948987516315E-03;
    R2I2[8] = 2; R2J2[8] = 4; R2Bn[8] = -0.39392777243355E-02;
    R2I2[9] = 2; R2J2[9] = 7; R2Bn[9] = -0.43797295650573E-01;
    R2I2[10] = 2; R2J2[10] = 36; R2Bn[10] = -0.26674547914087E-04;
    R2I2[11] = 3; R2J2[11] = 0; R2Bn[11] = 0.20481737692309E-07;
    R2I2[12] = 3; R2J2[12] = 1; R2Bn[12] = 0.43870667284435E-06;
    R2I2[13] = 3; R2J2[13] = 3; R2Bn[13] = -0.32277677238570E-04;
    R2I2[14] = 3; R2J2[14] = 6; R2Bn[14] = -0.15033924542148E-02;
    R2I2[15] = 3; R2J2[15] = 35; R2Bn[15] = -0.40668253562649E-01;
    R2I2[16] = 4; R2J2[16] = 1; R2Bn[16] = -0.78847309559367E-09;
    R2I2[17] = 4; R2J2[17] = 2; R2Bn[17] = 0.12790717852285E-07;
    R2I2[18] = 4; R2J2[18] = 3; R2Bn[18] = 0.48225372718507E-06;
    R2I2[19] = 5; R2J2[19] = 7; R2Bn[19] = 0.22922076337661E-05;
    R2I2[20] = 6; R2J2[20] = 3; R2Bn[20] = -0.16714766451061E-10;
    R2I2[21] = 6; R2J2[21] = 16; R2Bn[21] = -0.21171472321355E-02;
    R2I2[22] = 6; R2J2[22] = 35; R2Bn[22] = -0.23895741934104E+02;
    R2I2[23] = 7; R2J2[23] = 0; R2Bn[23] = -0.59059564324270E-17;
    R2I2[24] = 7; R2J2[24] = 11; R2Bn[24] = -0.12621808899101E-05;
    R2I2[25] = 7; R2J2[25] = 25; R2Bn[25] = -0.38946842435739E-01;
    R2I2[26] = 8; R2J2[26] = 8; R2Bn[26] = 0.11256211360459E-10;
    R2I2[27] = 8; R2J2[27] = 36; R2Bn[27] = -0.82311340897998E+01;
    R2I2[28] = 9; R2J2[28] = 13; R2Bn[28] = 0.19809712802088E-07;
    R2I2[29] = 10; R2J2[29] = 4; R2Bn[29] = 0.10406965210174E-18;
    R2I2[30] = 10; R2J2[30] = 10; R2Bn[30] = -0.10234747095929E-12;
    R2I2[31] = 10; R2J2[31] = 14; R2Bn[31] = -0.10018179379511E-08;
    R2I2[32] = 16; R2J2[32] = 29; R2Bn[32] = -0.80882908646985E-10;
    R2I2[33] = 16; R2J2[33] = 50; R2Bn[33] = 0.10693031879409E+00;
    R2I2[34] = 18; R2J2[34] = 57; R2Bn[34] = -0.33662250574171E+00;
    R2I2[35] = 20; R2J2[35] = 20; R2Bn[35] = 0.89185845355421E-24;
    R2I2[36] = 20; R2J2[36] = 35; R2Bn[36] = 0.30629316876232E-12;
    R2I2[37] = 20; R2J2[37] = 48; R2Bn[37] = -0.42002467698208E-05;
    R2I2[38] = 21; R2J2[38] = 21; R2Bn[38] = -0.59056029685639E-25;
    R2I2[39] = 22; R2J2[39] = 53; R2Bn[39] = 0.37826947613457E-05;
    R2I2[40] = 23; R2J2[40] = 39; R2Bn[40] = -0.12768608934681E-14;
    R2I2[41] = 24; R2J2[41] = 26; R2Bn[41] = 0.73087610595061E-28;
    R2I2[42] = 24; R2J2[42] = 40; R2Bn[42] = 0.55414715350778E-16;
    R2I2[43] = 24; R2J2[43] = 58; R2Bn[43] = -0.94369707241210E-06;

    R2I3[1] = 0; R2J3[1] = 0; R2A3[1] = 0.10898952318288E+04;
    R2I3[2] = 0; R2J3[2] = 1; R2A3[2] = 0.84951654495535E+03;
    R2I3[3] = 0; R2J3[3] = 2; R2A3[3] = -0.10781748091826E+03;
    R2I3[4] = 0; R2J3[4] = 3; R2A3[4] = 0.33153654801263E+02;
    R2I3[5] = 0; R2J3[5] = 7; R2A3[5] = -0.74232016790248E+01;
    R2I3[6] = 0; R2J3[6] = 20; R2A3[6] = 0.11765048724356E+02;
    R2I3[7] = 1; R2J3[7] = 0; R2A3[7] = 0.18445749355790E+01;
    R2I3[8] = 1; R2J3[8] = 1; R2A3[8] = -0.41792700549624E+01;
    R2I3[9] = 1; R2J3[9] = 2; R2A3[9] = 0.62478196935812E+01;
    R2I3[10] = 1; R2J3[10] = 3; R2A3[10] = -0.17344563108114E+02;
    R2I3[11] = 1; R2J3[11] = 7; R2A3[11] = -0.20058176862096E+03;
    R2I3[12] = 1; R2J3[12] = 9; R2A3[12] = 0.27196065473796E+03;
    R2I3[13] = 1; R2J3[13] = 11; R2A3[13] = -0.45511318285818E+03;
    R2I3[14] = 1; R2J3[14] = 18; R2A3[14] = 0.30919688604755E+04;
    R2I3[15] = 1; R2J3[15] = 44; R2A3[15] = 0.25226640357872E+06;
    R2I3[16] = 2; R2J3[16] = 0; R2A3[16] = -0.61707422868339E-02;
    R2I3[17] = 2; R2J3[17] = 2; R2A3[17] = -0.31078046629583E+00;
    R2I3[18] = 2; R2J3[18] = 7; R2A3[18] = 0.11670873077107E+02;
    R2I3[19] = 2; R2J3[19] = 36; R2A3[19] = 0.12812798404046E+09;
    R2I3[20] = 2; R2J3[20] = 38; R2A3[20] = -0.98554909623276E+09;
    R2I3[21] = 2; R2J3[21] = 40; R2A3[21] = 0.28224546973002E+10;
    R2I3[22] = 2; R2J3[22] = 42; R2A3[22] = -0.35948971410703E+10;
    R2I3[23] = 2; R2J3[23] = 44; R2A3[23] = 0.17227349913197E+10;
    R2I3[24] = 3; R2J3[24] = 24; R2A3[24] = -0.13551334240775E+05;
    R2I3[25] = 3; R2J3[25] = 44; R2A3[25] = 0.12848734664650E+08;
    R2I3[26] = 4; R2J3[26] = 12; R2A3[26] = 0.13865724283226E+01;
    R2I3[27] = 4; R2J3[27] = 32; R2A3[27] = 0.23598832556514E+06;
    R2I3[28] = 4; R2J3[28] = 44; R2A3[28] = -0.13105236545054E+08;
    R2I3[29] = 5; R2J3[29] = 32; R2A3[29] = 0.73999835474766E+04;
    R2I3[30] = 5; R2J3[30] = 36; R2A3[30] = -0.55196697030060E+06;
    R2I3[31] = 5; R2J3[31] = 42; R2A3[31] = 0.37154085996233E+07;
    R2I3[32] = 6; R2J3[32] = 34; R2A3[32] = 0.19127729239660E+05;
    R2I3[33] = 6; R2J3[33] = 44; R2A3[33] = -0.41535164835634E+06;
    R2I3[34] = 7; R2J3[34] = 28; R2A3[34] = -0.62459855192507E+02;

    R2I4[1] = 0; R2J4[1] = 0; R2A4[1] = 0.14895041079516E+04;
    R2I4[2] = 0; R2J4[2] = 1; R2A4[2] = 0.74307798314034E+03;
    R2I4[3] = 0; R2J4[3] = 2; R2A4[3] = -0.97708318797837E+02;
    R2I4[4] = 0; R2J4[4] = 12; R2A4[4] = 0.24742464705674E+01;
    R2I4[5] = 0; R2J4[5] = 18; R2A4[5] = -0.63281320016026E+00;
    R2I4[6] = 0; R2J4[6] = 24; R2A4[6] = 0.11385952129658E+01;
    R2I4[7] = 0; R2J4[7] = 28; R2A4[7] = -0.47811863648625E+00;
    R2I4[8] = 0; R2J4[8] = 40; R2A4[8] = 0.85208123431544E-02;
    R2I4[9] = 1; R2J4[9] = 0; R2A4[9] = 0.93747147377932E+00;
    R2I4[10] = 1; R2J4[10] = 2; R2A4[10] = 0.33593118604916E+01;
    R2I4[11] = 1; R2J4[11] = 6; R2A4[11] = 0.33809355601454E+01;
    R2I4[12] = 1; R2J4[12] = 12; R2A4[12] = 0.16844539671904E+00;
    R2I4[13] = 1; R2J4[13] = 18; R2A4[13] = 0.73875745236695E+00;
    R2I4[14] = 1; R2J4[14] = 24; R2A4[14] = -0.47128737436186E+00;
    R2I4[15] = 1; R2J4[15] = 28; R2A4[15] = 0.15020273139707E+00;
    R2I4[16] = 1; R2J4[16] = 40; R2A4[16] = -0.21764114219750E-02;
    R2I4[17] = 2; R2J4[17] = 2; R2A4[17] = -0.21810755324761E-01;
    R2I4[18] = 2; R2J4[18] = 8; R2A4[18] = -0.10829784403677E+00;
    R2I4[19] = 2; R2J4[19] = 18; R2A4[19] = -0.46333324635812E-01;
    R2I4[20] = 2; R2J4[20] = 40; R2A4[20] = 0.71280351959551E-04;
    R2I4[21] = 3; R2J4[21] = 1; R2A4[21] = 0.11032831789999E-03;
    R2I4[22] = 3; R2J4[22] = 2; R2A4[22] = 0.18955248387902E-03;
    R2I4[23] = 3; R2J4[23] = 12; R2A4[23] = 0.30891541160537E-02;
    R2I4[24] = 3; R2J4[24] = 24; R2A4[24] = 0.13555504554949E-02;
    R2I4[25] = 4; R2J4[25] = 2; R2A4[25] = 0.28640237477456E-06;
    R2I4[26] = 4; R2J4[26] = 12; R2A4[26] = -0.10779857357512E-04;
    R2I4[27] = 4; R2J4[27] = 18; R2A4[27] = -0.76462712454814E-04;
    R2I4[28] = 4; R2J4[28] = 24; R2A4[28] = 0.14052392818316E-04;
    R2I4[29] = 4; R2J4[29] = 28; R2A4[29] = -0.31083814331434E-04;
    R2I4[30] = 4; R2J4[30] = 40; R2A4[30] = -0.10302738212103E-05;
    R2I4[31] = 5; R2J4[31] = 18; R2A4[31] = 0.28217281635040E-06;
    R2I4[32] = 5; R2J4[32] = 24; R2A4[32] = 0.12704902271945E-05;
    R2I4[33] = 5; R2J4[33] = 40; R2A4[33] = 0.73803353468292E-07;
    R2I4[34] = 6; R2J4[34] = 28; R2A4[34] = -0.11030139238909E-07;
    R2I4[35] = 7; R2J4[35] = 2; R2A4[35] = -0.81456365207833E-13;
    R2I4[36] = 7; R2J4[36] = 28; R2A4[36] = -0.25180545682962E-10;
    R2I4[37] = 9; R2J4[37] = 1; R2A4[37] = -0.17565233969407E-17;
    R2I4[38] = 9; R2J4[38] = 40; R2A4[38] = 0.86934156344163E-14;

    R2I5[1] = -7; R2J5[1] = 0; R2A5[1] = -0.32368398555242E+13;
    R2I5[2] = -7; R2J5[2] = 4; R2A5[2] = 0.73263350902181E+13;
    R2I5[3] = -6; R2J5[3] = 0; R2A5[3] = 0.35825089945447E+12;
    R2I5[4] = -6; R2J5[4] = 2; R2A5[4] = -0.58340131851590E+12;
    R2I5[5] = -5; R2J5[5] = 0; R2A5[5] = -0.10783068217470E+11;
    R2I5[6] = -5; R2J5[6] = 2; R2A5[6] = 0.20825544563171E+11;
    R2I5[7] = -2; R2J5[7] = 0; R2A5[7] = 0.61074783564516E+06;
    R2I5[8] = -2; R2J5[8] = 1; R2A5[8] = 0.85977722535580E+06;
    R2I5[9] = -1; R2J5[9] = 0; R2A5[9] = -0.25745723604170E+05;
    R2I5[10] = -1; R2J5[10] = 2; R2A5[10] = 0.31081088422714E+05;
    R2I5[11] = 0; R2J5[11] = 0; R2A5[11] = 0.12082315865936E+04;
    R2I5[12] = 0; R2J5[12] = 1; R2A5[12] = 0.48219755109255E+03;
    R2I5[13] = 1; R2J5[13] = 4; R2A5[13] = 0.37966001272486E+01;
    R2I5[14] = 1; R2J5[14] = 8; R2A5[14] = -0.10842984880077E+02;
    R2I5[15] = 2; R2J5[15] = 4; R2A5[15] = -0.45364172676660E-01;
    R2I5[16] = 6; R2J5[16] = 0; R2A5[16] = 0.14559115658698E-12;
    R2I5[17] = 6; R2J5[17] = 1; R2A5[17] = 0.11261597407230E-11;
    R2I5[18] = 6; R2J5[18] = 4; R2A5[18] = -0.17804982240686E-10;
    R2I5[19] = 6; R2J5[19] = 10; R2A5[19] = 0.12324579690832E-06;
    R2I5[20] = 6; R2J5[20] = 12; R2A5[20] = -0.11606921130984E-05;
    R2I5[21] = 6; R2J5[21] = 16; R2A5[21] = 0.27846367088554E-04;
    R2I5[22] = 6; R2J5[22] = 20; R2A5[22] = -0.59270038474176E-03;
    R2I5[23] = 6; R2J5[23] = 22; R2A5[23] = 0.12918582991878E-02;
}

double if97_region2_pi(double p)
{
    return p/1.0;
}

double if97_region2_tau(double T)
{
    return 540.0/T;
}

double if97_region2_gamma_ideal_tau(double pi, double tau)
{
    double tau00 = 1.0;
    double tau01 = tau;
    double tau02 = tau01 * tau01;
    double taum1 = 1.0 / tau01;
    double taum2 = taum1 * taum1;
    double taum3 = taum1 * taum2;
    double taum4 = taum1 * taum3;
    double taum5 = taum1 * taum4;
    double taum6 = taum1 * taum5;

    double gt = R2An[1] * R2J1[1] * taum1 + R2An[2] * R2J1[2] * tau00 + R2An[3] * R2J1[3] * taum6
        + R2An[4] * R2J1[4] * taum5 + R2An[5] * R2J1[5] * taum4 + R2An[6] * R2J1[6] * taum3
        + R2An[7] * R2J1[7] * taum2 + R2An[8] * R2J1[8] * tau01 + R2An[9] * R2J1[9] * tau02;

    return gt;
}

double if97_region2_gamma_ideal_pi(double pi, double tau)
{
    return 1.0 / pi;
}

double if97_region2_gamma_ideal_tautau(double pi, double tau)
{
    double tau00 = 1.0;
    double tau01 = tau;
    double taum1 = 1.0 / tau01;
    double taum2 = taum1 * taum1;
    double taum3 = taum1 * taum2;
    double taum4 = taum1 * taum3;
    double taum5 = taum1 * taum4;
    double taum6 = taum1 * taum5;
    double taum7 = taum1 * taum6;
    double gtt = R2An[1] * R2J1[1] * (R2J1[1] - 1.0) * taum2
         + R2An[2] * R2J1[2] * (R2J1[2] - 1.0) * taum1
         + R2An[3] * R2J1[3] * (R2J1[3] - 1.0) * taum7
         + R2An[4] * R2J1[4] * (R2J1[4] - 1.0) * taum6
         + R2An[5] * R2J1[5] * (R2J1[5] - 1.0) * taum5
         + R2An[6] * R2J1[6] * (R2J1[6] - 1.0) * taum4
         + R2An[7] * R2J1[7] * (R2J1[7] - 1.0) * taum3
         + R2An[8] * R2J1[8] * (R2J1[8] - 1.0) * tau00
         + R2An[9] * R2J1[9] * (R2J1[9] - 1.0) * tau01;
    return gtt;
}


double if97_region2_gamma_real_tau(double pi, double tau)
{
    double pai01 = pi;
    double pai02 = pai01 * pai01;
    double pai03 = pai01 * pai02;
    double pai04 = pai01 * pai03;
    double pai05 = pai01 * pai04;
    double pai06 = pai01 * pai05;
    double pai07 = pai01 * pai06;
    double pai08 = pai01 * pai07;
    double pai09 = pai01 * pai08;
    double pai10 = pai01 * pai09;
    double pai16 = pai06 * pai10;
    double pai18 = pai02 * pai16;
    double pai20 = pai02 * pai18;
    double pai21 = pai01 * pai20;
    double pai22 = pai01 * pai21;
    double pai23 = pai01 * pai22;
    double pai24 = pai01 * pai23;

    double tn00 = 1.0;
    double tn01 = tau - 0.5;
    double tn02 = tn01 * tn01;
    double tn03 = tn01 * tn02;
    double tn05 = tn02 * tn03;
    double tn06 = tn01 * tn05;
    double tn07 = tn01 * tn06;
    double tn09 = tn02 * tn07;
    double tn10 = tn01 * tn09;
    double tn12 = tn02 * tn10;
    double tn13 = tn01 * tn12;
    double tn15 = tn02 * tn13;
    double tn19 = tn06 * tn13;
    double tn20 = tn01 * tn19;
    double tn24 = tn05 * tn19;
    double tn25 = tn01 * tn24;
    double tn28 = tn03 * tn25;
    double tn34 = tn06 * tn28;
    double tn35 = tn01 * tn34;
    double tn38 = tn03 * tn35;
    double tn39 = tn01 * tn38;
    double tn47 = tn09 * tn38;
    double tn49 = tn02 * tn47;
    double tn52 = tn03 * tn49;
    double tn56 = tn07 * tn49;
    double tn57 = tn01 * tn56;
    double tnm1 = 1.0 / tn01;

    double gt = R2Bn[2]  * pai01 * R2J2[2]  * tn00 + R2Bn[3]  * pai01 * R2J2[3]  * tn01
        + R2Bn[4]  * pai01 * R2J2[4]  * tn02 + R2Bn[5]  * pai01 * R2J2[5]  * tn05
        + R2Bn[6]  * pai02 * R2J2[6]  * tn00 + R2Bn[7]  * pai02 * R2J2[7]  * tn01
        + R2Bn[8]  * pai02 * R2J2[8]  * tn03 + R2Bn[9]  * pai02 * R2J2[9]  * tn06
        + R2Bn[10] * pai02 * R2J2[10] * tn35 + R2Bn[11] * pai03 * R2J2[11] * tnm1
        + R2Bn[12] * pai03 * R2J2[12] * tn00 + R2Bn[13] * pai03 * R2J2[13] * tn02
        + R2Bn[14] * pai03 * R2J2[14] * tn05 + R2Bn[15] * pai03 * R2J2[15] * tn34
        + R2Bn[16] * pai04 * R2J2[16] * tn00 + R2Bn[17] * pai04 * R2J2[17] * tn01
        + R2Bn[18] * pai04 * R2J2[18] * tn02 + R2Bn[19] * pai05 * R2J2[19] * tn06
        + R2Bn[20] * pai06 * R2J2[20] * tn02 + R2Bn[21] * pai06 * R2J2[21] * tn15
        + R2Bn[22] * pai06 * R2J2[22] * tn34 + R2Bn[23] * pai07 * R2J2[23] * tnm1
        + R2Bn[24] * pai07 * R2J2[24] * tn10 + R2Bn[25] * pai07 * R2J2[25] * tn24
        + R2Bn[26] * pai08 * R2J2[26] * tn07 + R2Bn[27] * pai08 * R2J2[27] * tn35
        + R2Bn[28] * pai09 * R2J2[28] * tn12 + R2Bn[29] * pai10 * R2J2[29] * tn03
        + R2Bn[30] * pai10 * R2J2[30] * tn09 + R2Bn[31] * pai10 * R2J2[31] * tn13
        + R2Bn[32] * pai16 * R2J2[32] * tn28 + R2Bn[33] * pai16 * R2J2[33] * tn49
        + R2Bn[34] * pai18 * R2J2[34] * tn56 + R2Bn[35] * pai20 * R2J2[35] * tn19
        + R2Bn[36] * pai20 * R2J2[36] * tn34 + R2Bn[37] * pai20 * R2J2[37] * tn47
        + R2Bn[38] * pai21 * R2J2[38] * tn20 + R2Bn[39] * pai22 * R2J2[39] * tn52
        + R2Bn[40] * pai23 * R2J2[40] * tn38 + R2Bn[41] * pai24 * R2J2[41] * tn25
        + R2Bn[42] * pai24 * R2J2[42] * tn39 + R2Bn[43] * pai24 * R2J2[43] * tn57;
    return gt;
}

double if97_region2_gamma_real_pi(double pi, double tau)
{
    double pai00 = 1.0;
    double pai01 = pi;
    double pai02 = pai01 * pai01;
    double pai03 = pai01 * pai02;
    double pai04 = pai01 * pai03;
    double pai05 = pai01 * pai04;
    double pai06 = pai01 * pai05;
    double pai07 = pai01 * pai06;
    double pai08 = pai01 * pai07;
    double pai09 = pai01 * pai08;
    double pai15 = pai06 * pai09;
    double pai17 = pai02 * pai15;
    double pai19 = pai02 * pai17;
    double pai20 = pai01 * pai19;
    double pai21 = pai01 * pai20;
    double pai22 = pai01 * pai21;
    double pai23 = pai01 * pai22;

    double tn00 = 1.0;
    double tn01 = tau - 0.5;
    double tn02 = tn01 * tn01;
    double tn03 = tn01 * tn02;
    double tn04 = tn01 * tn03;
    double tn06 = tn02 * tn04;
    double tn07 = tn01 * tn06;
    double tn08 = tn01 * tn07;
    double tn10 = tn02 * tn08;
    double tn11 = tn01 * tn10;
    double tn13 = tn02 * tn11;
    double tn14 = tn01 * tn13;
    double tn16 = tn02 * tn14;
    double tn20 = tn04 * tn16;
    double tn21 = tn01 * tn20;
    double tn25 = tn04 * tn21;
    double tn26 = tn01 * tn25;
    double tn29 = tn03 * tn26;
    double tn35 = tn06 * tn29;
    double tn36 = tn01 * tn35;
    double tn39 = tn03 * tn36;
    double tn40 = tn01 * tn39;
    double tn48 = tn08 * tn40;
    double tn50 = tn02 * tn48;
    double tn53 = tn03 * tn50;
    double tn57 = tn04 * tn53;
    double tn58 = tn01 * tn57;

    double gp =
        R2Bn[1]  * R2I2[1]  * pai00 * tn00 + R2Bn[2]  * R2I2[2]  * pai00 * tn01
        + R2Bn[3]  * R2I2[3]  * pai00 * tn02 + R2Bn[4]  * R2I2[4]  * pai00 * tn03
        + R2Bn[5]  * R2I2[5]  * pai00 * tn06 + R2Bn[6]  * R2I2[6]  * pai01 * tn01
        + R2Bn[7]  * R2I2[7]  * pai01 * tn02 + R2Bn[8]  * R2I2[8]  * pai01 * tn04
        + R2Bn[9]  * R2I2[9]  * pai01 * tn07 + R2Bn[10] * R2I2[10] * pai01 * tn36
        + R2Bn[11] * R2I2[11] * pai02 * tn00 + R2Bn[12] * R2I2[12] * pai02 * tn01
        + R2Bn[13] * R2I2[13] * pai02 * tn03 + R2Bn[14] * R2I2[14] * pai02 * tn06
        + R2Bn[15] * R2I2[15] * pai02 * tn35 + R2Bn[16] * R2I2[16] * pai03 * tn01
        + R2Bn[17] * R2I2[17] * pai03 * tn02 + R2Bn[18] * R2I2[18] * pai03 * tn03
        + R2Bn[19] * R2I2[19] * pai04 * tn07 + R2Bn[20] * R2I2[20] * pai05 * tn03
        + R2Bn[21] * R2I2[21] * pai05 * tn16 + R2Bn[22] * R2I2[22] * pai05 * tn35
        + R2Bn[23] * R2I2[23] * pai06 * tn00 + R2Bn[24] * R2I2[24] * pai06 * tn11
        + R2Bn[25] * R2I2[25] * pai06 * tn25 + R2Bn[26] * R2I2[26] * pai07 * tn08
        + R2Bn[27] * R2I2[27] * pai07 * tn36 + R2Bn[28] * R2I2[28] * pai08 * tn13
        + R2Bn[29] * R2I2[29] * pai09 * tn04 + R2Bn[30] * R2I2[30] * pai09 * tn10
        + R2Bn[31] * R2I2[31] * pai09 * tn14 + R2Bn[32] * R2I2[32] * pai15 * tn29
        + R2Bn[33] * R2I2[33] * pai15 * tn50 + R2Bn[34] * R2I2[34] * pai17 * tn57
        + R2Bn[35] * R2I2[35] * pai19 * tn20 + R2Bn[36] * R2I2[36] * pai19 * tn35
        + R2Bn[37] * R2I2[37] * pai19 * tn48 + R2Bn[38] * R2I2[38] * pai20 * tn21
        + R2Bn[39] * R2I2[39] * pai21 * tn53 + R2Bn[40] * R2I2[40] * pai22 * tn39
        + R2Bn[41] * R2I2[41] * pai23 * tn26 + R2Bn[42] * R2I2[42] * pai23 * tn40
        + R2Bn[43] * R2I2[43] * pai23 * tn58;
    return gp;
}

double if97_region2_gamma_real_tautau(double pi, double tau)
{
    double pai01 = pi;
    double pai02 = pai01 * pai01;
    double pai03 = pai01 * pai02;
    double pai04 = pai01 * pai03;
    double pai05 = pai01 * pai04;
    double pai06 = pai01 * pai05;
    double pai07 = pai01 * pai06;
    double pai08 = pai01 * pai07;
    double pai09 = pai01 * pai08;
    double pai10 = pai01 * pai09;
    double pai16 = pai06 * pai10;
    double pai18 = pai02 * pai16;
    double pai20 = pai02 * pai18;
    double pai21 = pai01 * pai20;
    double pai22 = pai01 * pai21;
    double pai23 = pai01 * pai22;
    double pai24 = pai01 * pai23;

    double tn00 = 1.0;
    double tn01 = tau - 0.5;
    double tn02 = tn01 * tn01;
    double tn04 = tn02 * tn02;
    double tn05 = tn01 * tn04;
    double tn06 = tn01 * tn05;
    double tn08 = tn02 * tn06;
    double tn09 = tn01 * tn08;
    double tn11 = tn02 * tn09;
    double tn12 = tn01 * tn11;
    double tn14 = tn02 * tn12;
    double tn18 = tn04 * tn14;
    double tn19 = tn01 * tn18;
    double tn23 = tn04 * tn19;
    double tn24 = tn01 * tn23;
    double tn27 = tn04 * tn23;
    double tn33 = tn06 * tn27;
    double tn34 = tn01 * tn33;
    double tn37 = tn04 * tn33;
    double tn38 = tn01 * tn37;
    double tn46 = tn08 * tn38;
    double tn48 = tn02 * tn46;
    double tn51 = tn05 * tn46;
    double tn55 = tn04 * tn51;
    double tn56 = tn01 * tn55;
    double tnm1 = 1.0 / tn01;
    double tnm2 = tnm1 * tnm1;

    double gtt = R2Bn[3]  * pai01 * R2J2[3]  * (R2J2[3]  - 1.0) * tn00
         + R2Bn[4]  * pai01 * R2J2[4]  * (R2J2[4]  - 1.0) * tn01
         + R2Bn[5]  * pai01 * R2J2[5]  * (R2J2[5]  - 1.0) * tn04
         + R2Bn[6]  * pai02 * R2J2[6]  * (R2J2[6]  - 1.0) * tnm1
         + R2Bn[7]  * pai02 * R2J2[7]  * (R2J2[7]  - 1.0) * tn00
         + R2Bn[8]  * pai02 * R2J2[8]  * (R2J2[8]  - 1.0) * tn02
         + R2Bn[9]  * pai02 * R2J2[9]  * (R2J2[9]  - 1.0) * tn05
         + R2Bn[10] * pai02 * R2J2[10] * (R2J2[10] - 1.0) * tn34
         + R2Bn[11] * pai03 * R2J2[11] * (R2J2[11] - 1.0) * tnm2
         + R2Bn[12] * pai03 * R2J2[12] * (R2J2[12] - 1.0) * tnm1
         + R2Bn[13] * pai03 * R2J2[13] * (R2J2[13] - 1.0) * tn01
         + R2Bn[14] * pai03 * R2J2[14] * (R2J2[14] - 1.0) * tn04
         + R2Bn[15] * pai03 * R2J2[15] * (R2J2[15] - 1.0) * tn33
         + R2Bn[16] * pai04 * R2J2[16] * (R2J2[16] - 1.0) * tnm1
         + R2Bn[17] * pai04 * R2J2[17] * (R2J2[17] - 1.0) * tn00
         + R2Bn[18] * pai04 * R2J2[18] * (R2J2[18] - 1.0) * tn01
         + R2Bn[19] * pai05 * R2J2[19] * (R2J2[19] - 1.0) * tn05
         + R2Bn[20] * pai06 * R2J2[20] * (R2J2[20] - 1.0) * tn01
         + R2Bn[21] * pai06 * R2J2[21] * (R2J2[21] - 1.0) * tn14
         + R2Bn[22] * pai06 * R2J2[22] * (R2J2[22] - 1.0) * tn33
         + R2Bn[23] * pai07 * R2J2[23] * (R2J2[23] - 1.0) * tnm2
         + R2Bn[24] * pai07 * R2J2[24] * (R2J2[24] - 1.0) * tn09
         + R2Bn[25] * pai07 * R2J2[25] * (R2J2[25] - 1.0) * tn23
         + R2Bn[26] * pai08 * R2J2[26] * (R2J2[26] - 1.0) * tn06
         + R2Bn[27] * pai08 * R2J2[27] * (R2J2[27] - 1.0) * tn34
         + R2Bn[28] * pai09 * R2J2[28] * (R2J2[28] - 1.0) * tn11
         + R2Bn[29] * pai10 * R2J2[29] * (R2J2[29] - 1.0) * tn02
         + R2Bn[30] * pai10 * R2J2[30] * (R2J2[30] - 1.0) * tn08
         + R2Bn[31] * pai10 * R2J2[31] * (R2J2[31] - 1.0) * tn12
         + R2Bn[32] * pai16 * R2J2[32] * (R2J2[32] - 1.0) * tn27
         + R2Bn[33] * pai16 * R2J2[33] * (R2J2[33] - 1.0) * tn48
         + R2Bn[34] * pai18 * R2J2[34] * (R2J2[34] - 1.0) * tn55
         + R2Bn[35] * pai20 * R2J2[35] * (R2J2[35] - 1.0) * tn18
         + R2Bn[36] * pai20 * R2J2[36] * (R2J2[36] - 1.0) * tn33
         + R2Bn[37] * pai20 * R2J2[37] * (R2J2[37] - 1.0) * tn46
         + R2Bn[38] * pai21 * R2J2[38] * (R2J2[38] - 1.0) * tn19
         + R2Bn[39] * pai22 * R2J2[39] * (R2J2[39] - 1.0) * tn51
         + R2Bn[40] * pai23 * R2J2[40] * (R2J2[40] - 1.0) * tn37
         + R2Bn[41] * pai24 * R2J2[41] * (R2J2[41] - 1.0) * tn24
         + R2Bn[42] * pai24 * R2J2[42] * (R2J2[42] - 1.0) * tn38
         + R2Bn[43] * pai24 * R2J2[43] * (R2J2[43] - 1.0) * tn56;
    return gtt;
}

double if97_h_pT_region2(double pres, double temK)
{
    double pi = if97_region2_pi(pres);
    double tau = if97_region2_tau(temK);
    double gt_ideal = if97_region2_gamma_ideal_tau(pi, tau);
    double gt_real = if97_region2_gamma_real_tau(pi, tau);
    return tau * (gt_ideal + gt_real) * ConstR * temK;
}


double if97_rho_pT_region2(double pres, double temK)
{
    double pi = if97_region2_pi(pres);
    double tau = if97_region2_tau(temK);
    double gp_ideal = if97_region2_gamma_ideal_pi(pi, tau);
    double gp_real = if97_region2_gamma_real_pi(pi, tau);
    return 1.0E+3 * pres / (pi * (gp_ideal + gp_real) * ConstR * temK);
}


double if97_cp_pT_region2(double pres, double temK)
{
    double pi = if97_region2_pi(pres);
    double tau = if97_region2_tau(temK);
    double gtt_i = if97_region2_gamma_ideal_tautau(pi, tau);
    double gtt_r = if97_region2_gamma_real_tautau(pi, tau);
    return - tau * tau * (gtt_i + gtt_r) * ConstR;
}


double if97_T_ph_region2(double pres, double enth)
{
    double pnA0 = 1.0;
    double pnA1 = pres;
    double pnA2 = pnA1 * pnA1;
    double pnA3 = pnA1 * pnA2;
    double pnA4 = pnA1 * pnA3;
    double pnA5 = pnA1 * pnA4;
    double pnA6 = pnA1 * pnA5;
    double pnA7 = pnA1 * pnA6;

    double etaA00 = 1.0;
    double etaA01 = 0.0005 * enth - 2.1;
    double etaA02 = etaA01 * etaA01;
    double etaA03 = etaA01 * etaA02;
    double etaA04 = etaA01 * etaA03;
    double etaA07 = etaA03 * etaA04;
    double etaA09 = etaA02 * etaA07;
    double etaA11 = etaA02 * etaA09;
    double etaA12 = etaA01 * etaA11;
    double etaA18 = etaA07 * etaA11;
    double etaA20 = etaA02 * etaA18;
    double etaA24 = etaA04 * etaA20;
    double etaA28 = etaA04 * etaA24;
    double etaA32 = etaA04 * etaA28;
    double etaA34 = etaA02 * etaA32;
    double etaA36 = etaA02 * etaA34;
    double etaA38 = etaA02 * etaA36;
    double etaA40 = etaA02 * etaA38;
    double etaA42 = etaA02 * etaA40;
    double etaA44 = etaA02 * etaA42;

    double pnB0 = 1.0;
    double pnB1 = pres - 2.0;
    double pnB2 = pnB1 * pnB1;
    double pnB3 = pnB1 * pnB2;
    double pnB4 = pnB1 * pnB3;
    double pnB5 = pnB1 * pnB4;
    double pnB6 = pnB1 * pnB5;
    double pnB7 = pnB1 * pnB6;
    double pnB9 = pnB2 * pnB7;

    double etaB00 = 1.0;
    double etaB01 = 0.0005 * enth - 2.6;
    double etaB02 = etaB01 * etaB01;
    double etaB04 = etaB02 * etaB02;
    double etaB06 = etaB02 * etaB04;
    double etaB08 = etaB02 * etaB06;
    double etaB12 = etaB04 * etaB08;
    double etaB18 = etaB06 * etaB12;
    double etaB24 = etaB06 * etaB18;
    double etaB28 = etaB04 * etaB24;
    double etaB40 = etaB12 * etaB28;

    double pnC00 = 1.0;
    double pnC01 = pres + 25.0;
    double pnC02 = pnC01 * pnC01;
    double pnC06 = pnC02 * pnC02 * pnC02;
    double pnCm1 = 1.0 / pnC01;
    double pnCm2 = pnCm1 * pnCm1;
    double pnCm5 = pnCm1 * pnCm2 * pnCm2;
    double pnCm6 = pnCm1 * pnCm5;
    double pnCm7 = pnCm1 * pnCm6;

    double etaC00 = 1.0;
    double etaC01 = 0.0005 * enth - 1.8;
    double etaC02 = etaC01 * etaC01;
    double etaC04 = etaC02 * etaC02;
    double etaC08 = etaC04 * etaC04;
    double etaC10 = etaC02 * etaC08;
    double etaC12 = etaC02 * etaC10;
    double etaC16 = etaC04 * etaC12;
    double etaC20 = etaC04 * etaC16;
    double etaC22 = etaC02 * etaC20;

    double aa = 0.12809002730136E-3;
    double bb = 0.67955786399241E+0;
    double cc = 0.90584278514723E+3;

    double p585 = (aa * enth - bb) * enth + cc;

    double tt = 0;
    if (pres <= 4.0) //********** Tph2A **********
    {
        tt = R2A3[1]   * pnA0 * etaA00
            + R2A3[2]  * pnA0 * etaA01
            + R2A3[3]  * pnA0 * etaA02
            + R2A3[4]  * pnA0 * etaA03
            + R2A3[5]  * pnA0 * etaA07
            + R2A3[6]  * pnA0 * etaA20
            + R2A3[7]  * pnA1 * etaA00
            + R2A3[8]  * pnA1 * etaA01
            + R2A3[9]  * pnA1 * etaA02
            + R2A3[10] * pnA1 * etaA03
            + R2A3[11] * pnA1 * etaA07
            + R2A3[12] * pnA1 * etaA09
            + R2A3[13] * pnA1 * etaA11
            + R2A3[14] * pnA1 * etaA18
            + R2A3[15] * pnA1 * etaA44
            + R2A3[16] * pnA2 * etaA00
            + R2A3[17] * pnA2 * etaA02
            + R2A3[18] * pnA2 * etaA07
            + R2A3[19] * pnA2 * etaA36
            + R2A3[20] * pnA2 * etaA38
            + R2A3[21] * pnA2 * etaA40
            + R2A3[22] * pnA2 * etaA42
            + R2A3[23] * pnA2 * etaA44
            + R2A3[24] * pnA3 * etaA24
            + R2A3[25] * pnA3 * etaA44
            + R2A3[26] * pnA4 * etaA12
            + R2A3[27] * pnA4 * etaA32
            + R2A3[28] * pnA4 * etaA44
            + R2A3[29] * pnA5 * etaA32
            + R2A3[30] * pnA5 * etaA36
            + R2A3[31] * pnA5 * etaA42
            + R2A3[32] * pnA6 * etaA34
            + R2A3[33] * pnA6 * etaA44
            + R2A3[34] * pnA7 * etaA28;
    }
    else if (pres <= p585) //********** Tph2B **********
    {
        tt = R2A4[1]   * pnB0 * etaB00
            + R2A4[2]  * pnB0 * etaB01
            + R2A4[3]  * pnB0 * etaB02
            + R2A4[4]  * pnB0 * etaB12
            + R2A4[5]  * pnB0 * etaB18
            + R2A4[6]  * pnB0 * etaB24
            + R2A4[7]  * pnB0 * etaB28
            + R2A4[8]  * pnB0 * etaB40
            + R2A4[9]  * pnB1 * etaB00
            + R2A4[10] * pnB1 * etaB02
            + R2A4[11] * pnB1 * etaB06
            + R2A4[12] * pnB1 * etaB12
            + R2A4[13] * pnB1 * etaB18
            + R2A4[14] * pnB1 * etaB24
            + R2A4[15] * pnB1 * etaB28
            + R2A4[16] * pnB1 * etaB40
            + R2A4[17] * pnB2 * etaB02
            + R2A4[18] * pnB2 * etaB08
            + R2A4[19] * pnB2 * etaB18
            + R2A4[20] * pnB2 * etaB40
            + R2A4[21] * pnB3 * etaB01
            + R2A4[22] * pnB3 * etaB02
            + R2A4[23] * pnB3 * etaB12
            + R2A4[24] * pnB3 * etaB24
            + R2A4[25] * pnB4 * etaB02
            + R2A4[26] * pnB4 * etaB12
            + R2A4[27] * pnB4 * etaB18
            + R2A4[28] * pnB4 * etaB24
            + R2A4[29] * pnB4 * etaB28
            + R2A4[30] * pnB4 * etaB40
            + R2A4[31] * pnB5 * etaB18
            + R2A4[32] * pnB5 * etaB24
            + R2A4[33] * pnB5 * etaB40
            + R2A4[34] * pnB6 * etaB28
            + R2A4[35] * pnB7 * etaB02
            + R2A4[36] * pnB7 * etaB28
            + R2A4[37] * pnB9 * etaB01
            + R2A4[38] * pnB9 * etaB40;
    }
    else //********** Tph2C **********
    {
        tt = R2A5[1]   * pnCm7 * etaC00
            + R2A5[2]  * pnCm7 * etaC04
            + R2A5[3]  * pnCm6 * etaC00
            + R2A5[4]  * pnCm6 * etaC02
            + R2A5[5]  * pnCm5 * etaC00
            + R2A5[6]  * pnCm5 * etaC02
            + R2A5[7]  * pnCm2 * etaC00
            + R2A5[8]  * pnCm2 * etaC01
            + R2A5[9]  * pnCm1 * etaC00
            + R2A5[10] * pnCm1 * etaC02
            + R2A5[11] * pnC00 * etaC00
            + R2A5[12] * pnC00 * etaC01
            + R2A5[13] * pnC01 * etaC04
            + R2A5[14] * pnC01 * etaC08
            + R2A5[15] * pnC02 * etaC04
            + R2A5[16] * pnC06 * etaC00
            + R2A5[17] * pnC06 * etaC01
            + R2A5[18] * pnC06 * etaC04
            + R2A5[19] * pnC06 * etaC10
            + R2A5[20] * pnC06 * etaC12
            + R2A5[21] * pnC06 * etaC16
            + R2A5[22] * pnC06 * etaC20
            + R2A5[23] * pnC06 * etaC22;
    }

    return tt;
}


double if97_T_ph_region2_jsme(double pres, double enth)
{
    const double eps = 1.0E-6;
    int converged = 0;
    double tt = if97_T_ph_region2(pres, enth);
    for (int i = 1; i <= 10; i++)
    {
        double hh = if97_h_pT_region2(pres, tt);
        double dlt = enth - hh;

        if (fabs(dlt) <= eps)
        {
            converged = 1;
            break;
        }

        double cp = if97_cp_pT_region2(pres, tt);
        tt += dlt / cp;
    }

    if (!converged)
    {
        printf("%s not converged! pres = %g, h=%g\n", __FUNCTION__, pres, enth);
        return -1;
    }

    return tt;
}

