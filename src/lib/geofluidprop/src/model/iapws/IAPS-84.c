
#include "IAPS-84.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

#define GAS_CONST 0.461522
#define TEMP_ZERO 647.073
#define UREF -4328.455039
#define SREF 7.6180802

#define TP_TEMP 273.16     /***** Triple point temperature [k] *****/
#define TP_PRES 0.00061173 /***** Triple point pressure [MPa]  *****/

#define CP_TEMP 647.126    /***** Critical point temperature [k] *****/
#define CP_PRES 22.055     /***** Critical point pressure [MPa]  *****/

#define MAX_TEMP 2500.0
#define MAX_PRES 3000.0
// #define MAX_TEMP 1273.15
// #define MAX_PRES 1500.0

#define DCNV_GCM3_TO_KGM3 1000.0 /***** Convert density unit from g/cm^3 to kg/m^3 *****/
#define DCNV_KGM3_TO_GCM3 0.001  /***** Convert density unit from kg/m^3 to g/cm^3 *****/

#define JCNV_KJKG_TO_JKG 1000.0 /***** Convert specific energy unit from kJ/kg to j/kg *****/
#define JCNV_JKG_TO_KJKG 0.001  /***** Convert specific energy unit from J/kg to kJ/kg *****/


void calc_bb(double temp, double* b1, double* b2)
{
	static const double BP[10] = {0.7478629, -0.3540782, 0.0, 0.0, 0.007159876, 0.0, -0.003528426, 0.0, 0.0, 0.0};
	static const double BQ[10] = {1.1278334, 0.0, -0.5944001, -5.010996, 0.0, 0.63684256, 0.0, 0.0, 0.0, 0.0};


	double tzr = TEMP_ZERO / temp;
	double tr1 = 1.0 / temp;
	double tr2 = 1.0 / (temp * temp);

	double ri = 0.0;

	double vv[10] = { 1.0 };


	for (int i = 1; i < 10; i++)
	{
		vv[i] = vv[i - 1] * tzr;
	}

	b1[0] = BP[0] + BP[1] * log(1.0 / vv[1]);
	b2[0] = BQ[0];

	b1[1] = BP[1] * vv[1] / TEMP_ZERO;
	b2[1] = 0.0;

	b1[2] = 0.0;
	b2[2] = 0.0;

	for (int i = 2; i < 10; i++)
	{
		ri = (double)i - 1.0;

		b1[0] += BP[i] * vv[i - 1];
		b2[0] += BQ[i] * vv[i - 1];

		b1[1] -= ri * BP[i] * vv[i - 1] * tr1;
		b2[1] -= ri * BQ[i] * vv[i - 1] * tr1;

		b1[2] += BP[i] * ri * ri * vv[i - 1] * tr2;
		b2[2] += BQ[i] * ri * ri * vv[i - 1] * tr2;
	}

	b1[2] -= b1[1] * tr1;
	b2[2] -= b2[1] * tr1;
}

void calc_base(double dens, double temp, double* b1, double* b2, double* basef)
{
	static const double GCPZ = GAS_CONST / 0.101325;  //***** GasConst / P0 *****

	static const double G1 = 11.0, G2 = 44.333333333333, GF = 3.5, G4 = 28.16666667, G5 = 15.166666667;


	double b20b10g = b2[0] / b1[0] - GF;
	double b11b10t = b1[1] / b1[0] * temp;

	double xx = 0.0, yy = 0.0, zz = 0.0, dz = 0.0, z0 = 0.0, d0 = 0.0;


	yy = 0.25 * b1[0] * dens;

	xx = 1.0 - yy;

	z0 = (1.0 + G1 * yy + G2 * yy * yy) / (xx * xx * xx);
	zz = z0 + 4.0 * yy * b20b10g;

	d0 = (G1 + 2.0 * G2 * yy) / (xx * xx * xx) + 3.0 * (1.0 + G1 * yy + G2 * yy * yy) / (xx * xx * xx * xx);
	dz = d0 + 4.0 * b20b10g;

	basef[0] = -log(xx) - (G2 - 1.0) / xx + G4 / (xx * xx) + 4.0 * yy * b20b10g + G5 + log(dens * temp * GCPZ);

	basef[1] = basef[0] + zz;

	basef[3] = -b11b10t * (zz - 1.0 - dens * b2[0]) - dens * temp * b2[1];

	basef[4] = basef[3] + zz;

	basef[2] = basef[3] - basef[0];

	basef[5] = 2.0 * basef[3] + (z0 - 1.0) * (b11b10t * b11b10t - temp * temp * b1[2] / b1[0])
			 - dens * temp * temp * (b2[2] - GF * b1[2]) - b11b10t * b11b10t * yy * d0;

	basef[6] = zz / temp + dens * (0.25 * dz * b1[1] + b2[1] - b2[0] / b1[0] * b1[1]);

	basef[7] = yy;
	basef[8] = zz;
	basef[9] = dz;
}

void calc_residual(double dens, double temp, double* resif)
{
	static const double GG[40] = { -5.3062968529023e+2,  2.2744901424408e+3,  7.8779333020687e+2, -6.9830527374994e+1,
									1.7863832875422e+4, -3.9514731563338e+4,  3.3803884280753e+4, -1.3855050202703e+4,
								   -2.5637436613260e+5,  4.8212575981415e+5, -3.4183016969660e+5,  1.2223156417448e+5,
									1.1797433655832e+6, -2.1734810110373e+6,  1.0829952168620e+6, -2.5441998064049e+5,
								   -3.1377774947767e+6,  5.2911910757704e+6, -1.3802577177877e+6, -2.5109914369001e+5,
									4.6561826115608e+6, -7.2752773275387e+6,  4.1774246148294e+5,  1.4016358244614e+6,
								   -3.1555231392127e+6,  4.7929666384584e+6,  4.0912664781209e+5, -1.3626369388386e+6,
									6.9625220862664e+5, -1.0834900096447e+6, -2.2722827401688e+5,  3.8365486000660e+5,
									6.8833257944332e+3,  2.1757245522644e+4, -2.6627944829770e+3, -7.0730418082074e+4,
								   -0.225, -1.68, 0.055, -93.0 };

	static const int II[40] = { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4,
								5, 5, 5, 5, 6, 6, 6, 6, 8, 8, 8, 8, 2, 2, 0, 4, 2, 2, 2, 4 };

	static const int JJ[40] = { 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7,
								2, 3, 5, 7, 2, 3, 5, 7, 2, 3, 5, 7, 1, 4, 4, 4, 0, 2, 0, 0 };

	static const double ATZ[4] = { 640.0, 640.0, 641.6, 270.0 };

	static const double ADZ[4] = { 0.319, 0.319, 0.319, 1.55 };

	static const double AAT[4] = { 20000.0, 20000.0, 40000.0, 25.0 };

	static const double AAD[4] = { 34.0, 40.0, 30.0, 1050.0 };


	int iii = 0, jjj = 0;

	double rt = 0.0, dadt = 0.0, ee = 0.0, q10 = 0.0, q20 = 0.0, vv = 0.0, qp = 0.0, ri = 0.0, rj = 0.0;

	double dfdt = 0.0, d2f = 0.0, dpt = 0.0;

	double ddz = 0.0, del = 0.0, dd = 0.0, ex1 = 0.0, dex = 0.0, att = 0.0, tx = 0.0, tau = 0.0, ex2 = 0.0, tex = 0.0;

	double qm = 0.0, fct = 0.0, q5t = 0.0, q2a = 0.0;

	double qr[11] = { 0.0 }, qt[10] = { 0.0 }, qzr[9] = { 0.0 }, qrt[9] = { 0.0 };


	rt = GAS_CONST * temp;

	ee = exp(-dens);

	q10 = dens * dens * ee;
	q20 = 1.0 - ee;

	qr[1] = q10;

	vv = TEMP_ZERO / temp;

	qt[0] = temp / TEMP_ZERO;

	for (int i = 1; i < 10; i++)
	{
		qr[i + 1] = qr[i] * q20;

		qt[i] = qt[i - 1] * vv;
	}

	for (int i = 0; i < 9; i++)	resif[i] = 0.0;

	for (int i = 0; i < 36; i++)
	{
		iii = II[i] + 1;
		jjj = JJ[i];

		ri = (double)iii;
		rj = (double)jjj;

		qp = GG[i] * qr[iii] * qt[jjj];

		resif[7] += qp;

		resif[8] += (2.0 / dens - (1.0 - ee * (ri - 1.0) / q20)) * qp;

		resif[0] += GG[i] * qr[iii + 1] * qt[jjj] / (q10 * ri * rt);

		dfdt = pow(q20, ri) * (1.0 - rj) * qt[jjj + 1] / (TEMP_ZERO * ri);

		d2f = rj * dfdt;

		dpt = dfdt * q10 * ri / q20;

		dadt += GG[i] * dfdt;

		resif[6] += GG[i] * dpt;

		resif[5] += GG[i] * d2f / GAS_CONST;
	}

	qp = 0.0;

	for (int j = 36; j < 40; j++)
	{
		iii = II[j];
		jjj = JJ[j];

		ri = (double)iii;
		rj = (double)jjj;

		ddz = ADZ[j - 36];

		del = dens / ddz - 1.0;

		if (fabs(del) < 1.0e-10) del = 1.0e-10;

		dd = del * del;

		ex1 = -AAD[j - 36] * pow(del, ri);

		dex = exp(ex1) * pow(del, rj);

		att = AAT[j - 36];

		tx = ATZ[j - 36];

		tau = temp / tx - 1.0;

		ex2 = -att * tau * tau;

		tex = exp(ex2);

		q10 = dex * tex;

		qm = rj / del - ri * AAD[j - 36] * pow(del, ri - 1.0);

		fct = qm * dens * dens * q10 / ddz;

		q5t = fct * (2.0 / dens + qm / ddz) - (dens / ddz) * (dens / ddz) * q10
			* (rj / (del * del) + ri * (ri - 1.0) * AAD[j - 36] * pow(del, ri - 2.0));

		resif[8] += q5t * GG[j];

		qp += GG[j] * fct;

		dadt -= 2.0 * GG[j] * att * tau * q10 / tx;

		resif[6] -= 2.0 * GG[j] * att * tau * fct / tx;

		q2a += temp * GG[j] * (4.0 * att * ex2 + 2.0 * att) * q10 / (tx * tx);

		resif[0] += q10 * GG[j] / rt;
	}

	resif[2] = -dadt / GAS_CONST;

	resif[3] = resif[0] + resif[2];

	resif[5] += q2a / GAS_CONST;

	resif[7] += qp;
}

void calc_ideal(double temp, double* igasf)
{
	static const double CC[18] = { 1.97302710180e+01,  2.09662681977e+01, -4.83429455355e-01,  6.05743189245e+00,
								   2.25602388500e+01, -9.87532442000e+00, -4.31355385130e+00,  4.58155781000e-01,
								  -4.77549018830e-02,  4.12384606330e-03, -2.79290528520e-04,  1.44816952610e-05,
								  -5.64736587480e-07,  1.62004460000e-08, -3.30382279600e-10,  4.51916067368e-12,
								  -3.70734122708e-14,  1.37546068238e-16 };


	double tt = 0.0, tl = 0.0, ri = 0.0, rw = 0.0;


	tt = 0.01 * temp;
	tl = log(tt);

	igasf[1] = -(CC[0] / tt + CC[1]) * tl;

	igasf[4] = CC[1] + CC[0] * (1.0 - tl) / tt;

	igasf[6] = CC[1] - CC[0] / tt;

	for (int i = 2; i < 18; i++)
	{
		ri = (double)i;
		rw = pow(tt, ri - 5.0);

		igasf[1] -= CC[i] * rw;

		igasf[4] += CC[i] * (ri - 5.0) * rw;

		igasf[6] += CC[i] * (ri - 5.0) * (ri - 4.0) * rw;
	}

	igasf[0] = igasf[1] - 1.0;

	igasf[3] = igasf[4] - 1.0;

	igasf[5] = igasf[6] - 1.0;

	igasf[2] = igasf[3] - igasf[0];
}

void dfind1(double pres, double temp, double dess, double* dens, double* dpdd)
{
	double dd = 0.0, rt = 0.0, pp = 0.0, dpdx = 0.0, dp = 0.0, xx = 0.0;

	double  b1[3] = { 0.0 }, b2[3] = { 0.0 }, basef[10] = { 0.0 }, resif[9] = { 0.0 };


	dd = dess;

	if (dd <= 0.0) dd = 1.0e-8;

	if (dd > 1.9) dd = 1.9;

	rt = GAS_CONST * temp;

	calc_bb(temp, b1, b2);

	for (int i = 1; i <= 31; i++)
	{
		if (dd <= 0.0) dd = 1.0e-8;

		if (dd > 1.9) dd = 1.9;

		calc_base(dd, temp, b1, b2, basef);

		calc_residual(dd, temp, resif);

		pp = rt * dd * basef[8] + resif[7];

		*dpdd = rt * (basef[8] + basef[7] * basef[9]) + resif[8];

		if (*dpdd <= 0.0)
		{
			if (dess >= 0.2967)
			{
				dd *= 1.02;
			}
			else
			{
				dd *= 0.98;
			}

			if (i <= 10) continue;
		}

		dpdx = *dpdd * 1.1;

		if (dpdx < 0.1) dpdx = 0.1;

		dp = fabs(1.0 - pp / pres);

		if (dp < 1.0e-10) break;

		if (dess > 0.3 && dp < 1.0e-7) break;

		if (dess > 0.7 && dp < 1.0e-6) break;

		xx = (pres - pp) / dpdx;

		if (fabs(xx) > 0.1) xx *= 0.1 / fabs(xx);

		dd += xx;

		if (dd <= 0.0) dd = 1.0e-8;
	}

	*dens = dd;
}

void dfind2(double pres, double temp, double dess, double* b1, double* b2, double* basef, double* resif, double* dens, double* dpdd)
{
	double dd = 0.0, rt = 0.0, pp = 0.0, dpdx = 0.0, dp = 0.0, xx = 0.0;


	dd = dess;

	if (dd <= 0.0) dd = 1.0e-8;

	if (dd > 1.9) dd = 1.9;

	rt = GAS_CONST * temp;

	for (int i = 1; i <= 31; i++)
	{
		if (dd <= 0.0) dd = 1.0e-8;

		if (dd > 1.9) dd = 1.9;

		calc_base(dd, temp, b1, b2, basef);

		calc_residual(dd, temp, resif);

		pp = rt * dd * basef[8] + resif[7];

		*dpdd = rt * (basef[8] + basef[7] * basef[9]) + resif[8];

		if (*dpdd <= 0.0)
		{
			if (dess >= 0.2967)
			{
				dd *= 1.02;
			}
			else
			{
				dd *= 0.98;
			}

			if (i <= 10) continue;
		}

		dpdx = *dpdd * 1.1;

		if (dpdx < 0.1) dpdx = 0.1;

		dp = fabs(1.0 - pp / pres);

		if (dp < 1.0e-10) break;

		if (dess > 0.3 && dp < 1.0e-7) break;

		if (dess > 0.7 && dp < 1.0e-6) break;

		xx = (pres - pp) / dpdx;

		if (fabs(xx) > 0.1) xx *= 0.1 / fabs(xx);

		dd += xx;

		if (dd <= 0.0) dd = 1.0e-8;
	}

	*dens = dd;
}

static double psat(double temp)
{
	static const double AA[8] = { -7.8889166e+0,  2.5514255e+0, -6.7161690e+0,  3.3239495e+1,
								  -1.0538479e+2,  1.7435319e+2, -1.4839348e+2,  4.8631602e+1 };


	double pl = 0.0, ps = 0.0, vv = 0.0, ww = 0.0, bb = 0.0, zz = 0.0, qq = 0.0;


	if (temp <= 314.0)
	{
		pl = 6.3573118 - 8858.843 / temp + 607.56335 * pow(temp, -0.6);

		ps = 0.1 * exp(pl);

		return ps;
	}

	vv = temp / 647.25;

	ww = fabs(1.0 - vv);

	for (int i = 0; i < 8; i++)
	{
		zz = ((double)i + 2.0) / 2.0;

		bb += AA[i] * pow(ww, zz);
	}

	qq = bb / vv;

	ps = 22.093 * exp(qq);


	return ps;
}

double tdpsdt(double temp)
{
	static const double AA[8] = { -7.8889166e+0,  2.5514255e+0, -6.7161690e+0,  3.3239495e+1,
								  -1.0538479e+2,  1.7435319e+2, -1.4839348e+2,  4.8631602e+1 };


	double tdpsdt = 0.0, vv = 0.0, ww = 0.0, bb = 0.0, cc = 0.0, z1 = 0.0, z2 = 0.0, yy = 0.0, qq = 0.0;


	vv = temp / 647.25;

	ww = 1.0 - vv;

	for (int i = 0; i < 8; i++)
	{
		z1 = (double)i + 1.0;

		z2 = (z1 + 1.0) / 2.0;

		yy = AA[i] * pow(ww, z2);

		bb += yy;

		cc += yy / ww * (0.5 - 0.5 * z1 - 1.0 / vv);
	}

	qq = bb / vv;

	tdpsdt = 22.093 * exp(qq) * cc;


	return tdpsdt;
}

static double tsat(double pres)
{
	double pp = 0.0, tg = 0.0, pl = 0.0, dp = 0.0;


	if (pres > CP_PRES) return 0.0;

	pl = 2.302585 + log(pres);

	tg = 372.83 + pl * (27.7589 + pl * (2.3819 + pl * (0.24834 + pl * 0.0193855)));

	for (int i = 0; i < 8; i++)
	{
		if (tg < 273.15) tg = 273.15;
		if (tg > 647.0)  tg = 647.0;

		pp = psat(tg);

		dp = tdpsdt(tg);

		if (fabs(1.0 - pp / pres) < 1.0e-5) break;

		tg *= 1.0 + (pres - pp) / dp;
	}


	return tg;
}

static void corr(double temp, double* pres, double* dl, double* dv, double* delg)
{
	double rt = 0.0, dliq = 0.0, dvap = 0.0, dpd = 0.0, gl = 0.0, gv = 0.0, tau = 0.0;

	double b1[3] = { 0.0 }, b2[3] = { 0.0 }, basef[10] = { 0.0 }, resif[9] = { 0.0 }, igasf[7] = { 0.0 }, fcts[11] = { 0.0 };


	rt = GAS_CONST * temp;

	if (temp <= 646.3)
	{
		calc_bb(temp, b1, b2);

		calc_ideal(temp, igasf);

		if(*dl > 0.0)
		{
			dliq = *dl;
		}
		else
		{
			dliq = 1.11 - 0.0004 * temp;
		}

		dfind2(*pres, temp, dliq, b1, b2, basef, resif, dl, &dpd);

		gl = basef[0] + resif[0] + igasf[0] - UREF / temp + SREF + basef[8] + resif[7] / (rt * (*dl));

		if(*dv > 0.0)
		{
			dvap = *dv;
		}
		else
		{
			dvap = *pres / rt;
		}

		dfind2(*pres, temp, dvap, b1, b2, basef, resif, dv, &dpd);

		if (*dv < 5.0e-7) *dv = 5.0e-7;

		gv = basef[0] + resif[0] + igasf[0] - UREF / temp + SREF + basef[8] + resif[7] / (rt * (*dv));

		*delg = gl - gv;


		return;
	}


	*pres = 0.0;

	if (temp > 647.126) return;

	*delg = 0.0;

	calc_bb(temp, b1, b2);

	tau = 0.657128 * pow((1.0 - temp / 647.126), 0.325);

	*dl = 0.32189 + tau;
	*dv = 0.32189 - tau;

	calc_base(*dv, temp, b1, b2, basef);

	calc_residual(*dv, temp, resif);

	*pres = rt * (*dv) * basef[8] + resif[7];
}

static void pcorr(double temp, double* pres, double* dl, double* dv)
{
	double delg = 0.0, dp = 0.0;


	*pres = psat(temp);

	for (;;)
	{
		corr(temp, pres, dl, dv, &delg);

		dp = delg * GAS_CONST * temp / (1.0 / *dv - 1.0 / *dl);

		*pres += dp;

		if (fabs(delg) < 1.0e-4) return;
	}
}

static void tcorr(double* pres, double* temp, double* dl, double* dv)
{
	double delg = 0.0, dp = 0.0;


	*temp = tsat(*pres);

	if (*temp <= 0.0) return;

	for (;;)
	{
		corr(*temp, pres, dl, dv, &delg);

		dp = delg * GAS_CONST * (*temp) / (1.0 / *dv - 1.0 / *dl);

		*temp *= 1.0 - dp / tdpsdt(*temp);

		if (fabs(delg) < 1.0e-4) return;
	}
}


static void baseParam_1(double dens, double temp, double* basef, double* resif)
{
	double b1[3] = { 0.0 }, b2[3] = { 0.0 };

	calc_bb(temp, b1, b2);

	calc_base(dens, temp, b1, b2, basef);

	calc_residual(dens, temp, resif);
}

static void baseParam_2(double dens, double temp, double* basef, double* resif, double* igasf)
{
	double b1[3] = { 0.0 }, b2[3] = { 0.0 };

	calc_bb(temp, b1, b2);

	calc_base(dens, temp, b1, b2, basef);

	calc_residual(dens, temp, resif);

	calc_ideal(temp, igasf);
}


double calc_A_for_DT(double dens, double temp)
{
	/*****  Calculate Helmholtz Function  ****/

	double basef[10] = { 0.0 }, resif[9] = { 0.0 }, igasf[7] = { 0.0 };

	baseParam_2(dens, temp, basef, resif, igasf);

	return GAS_CONST * temp * (basef[0] + resif[0] + igasf[0] - UREF / temp + SREF);
}

double calc_G_for_DT(double dens, double temp)
{
	/*****  Calculate Gibbs Function  ****/

	double basef[10] = { 0.0 }, resif[9] = { 0.0 }, igasf[7] = { 0.0 };

	baseParam_2(dens, temp, basef, resif, igasf);

	return GAS_CONST * temp * (basef[0] + resif[0] + igasf[0] - UREF / temp + SREF + basef[8]) + resif[7] / dens;
}

double calc_U_for_DT(double dens, double temp)
{
	double basef[10] = { 0.0 }, resif[9] = { 0.0 }, igasf[7] = { 0.0 };

	baseParam_2(dens, temp, basef, resif, igasf);

	return GAS_CONST * temp * (basef[3] + resif[3] + igasf[3] - UREF / temp);
}

double calc_H_for_DT(double dens, double temp)
{
	double basef[10] = { 0.0 }, resif[9] = { 0.0 }, igasf[7] = { 0.0 };

	baseParam_2(dens, temp, basef, resif, igasf);

	return GAS_CONST * temp * (basef[3] + resif[3] + igasf[3] - UREF / temp + basef[8]) + resif[7] / dens;
}

double calc_S_for_DT(double dens, double temp)
{
	double basef[10] = { 0.0 }, resif[9] = { 0.0 }, igasf[7] = { 0.0 };

	baseParam_2(dens, temp, basef, resif, igasf);

	return GAS_CONST * (basef[2] + resif[2] + igasf[2] - SREF);
}

double calc_P_for_DT(double dens, double temp)
{
	double basef[10] = { 0.0 }, resif[9] = { 0.0 };

	baseParam_1(dens, temp, basef, resif);

	return GAS_CONST * temp * dens * basef[8] + resif[7];
}

double calc_dPdD_for_DT(double dens, double temp)
{
	double basef[10] = { 0.0 }, resif[9] = { 0.0 };

	baseParam_1(dens, temp, basef, resif);

	return GAS_CONST * temp * (basef[8] + basef[7] * basef[9]) + resif[8];
}

double calc_dPdT_for_DT(double dens, double temp)
{
	double basef[10] = { 0.0 }, resif[9] = { 0.0 };

	baseParam_1(dens, temp, basef, resif);

	return GAS_CONST * temp * dens * basef[6] + resif[6];
}

double calc_Cp_for_DT(double dens, double temp)
{
	double basef[10] = { 0.0 }, resif[9] = { 0.0 }, igasf[7] = { 0.0 };

	double dpdd = 0.0, dpdt = 0.0;

	baseParam_2(dens, temp, basef, resif, igasf);

	dpdd = GAS_CONST * temp * (basef[8] + basef[7] * basef[9]) + resif[8];

	dpdt = GAS_CONST * temp * dens * basef[6] + resif[6];

	return GAS_CONST * (basef[5] + resif[5] + igasf[5]) + temp * dpdt * dpdt / (dens * dens * dpdd);
}

double calc_Cv_for_DT(double dens, double temp)
{
	double basef[10] = { 0.0 }, resif[9] = { 0.0 }, igasf[7] = { 0.0 };

	baseParam_2(dens, temp, basef, resif, igasf);

	return GAS_CONST * (basef[5] + resif[5] + igasf[5]);
}


double calc_D_for_PT(double pres, double temp)
{
	double dens = 0.0, dpdd = 0.0;

	double dgss = pres / (0.4 * temp);

	double psat = 2.0e+4, dll = 0.0, dvv = 0.0;


	if (temp < TEMP_ZERO) pcorr(temp, &psat, &dll, &dvv);

	if (pres > psat) dgss = dll;

	dfind1(pres, temp, dgss, &dens, &dpdd);

	return dens;
}

#define RTOL 1e-8

double calc_T_for_DP(double dens, double pres)
{
	// ad-hoc sampling for initial guess
	double T0 = 0;
	{
		double temp1 = 0.01 + 273.15, temp2 = 999.99 + 273.15;
		double tempM = 0.5 * (temp1 + temp2);
		for (int n = 0; n < 5; n++)
		{
			double presM = calc_P_for_DT(dens, tempM);
			double error = presM - pres;
			if (error < 0.0)
				temp1 = tempM;
			else
				temp2 = tempM;
			tempM = 0.5 * (temp1 + temp2);
		}
		T0 = tempM;
	}

	double temp = T0;
	const double dT = 100*1e-3;
	for (int n=0; n<100; n++)
	{
		double p1 = calc_P_for_DT(dens, temp);
		double res = p1 - pres;
		// printf("%d: T=%g, r=%g\n", n, temp, res);
		if (fabs(res)<RTOL*pres)
			return temp;
		double p2 = calc_P_for_DT(dens, temp + dT);
		double drdT = (p2-p1)/dT;
		temp += -res / drdT;
		temp = temp >= 0 ? temp : iaps84_Tt(); //TODO
	}
	assert(false);
	return -1;
}

double calc_D_for_PH(double pres, double enth)
{
	if (pres < iaps84_pc())
	{
		// check two phase
		double hl = iaps84_sat_hl_p(pres);
		double hv = iaps84_sat_hv_p(pres);
		if (hl < enth && enth < hv)
		{
			double wv = (enth - hl) / (hv - hl);
			double rhol = iaps84_sat_rhol_p(pres);
			double rhov = iaps84_sat_rhov_p(pres);
			double v = 1./rhol + wv * (1./rhov - 1./rhol);
			return 1./v;
		}
	}

	// get initial guess
	double densM = 0;
	{
		double temp1 = 0.01 + 273.15, temp2 = 999.99 + 273.15, tempM = 0.0;
		double dens1 = 0.0, dens2 = 0.0, enthM = 0.0, error = 0.0;

		dens1 = calc_D_for_PT(pres, temp1);
		dens2 = calc_D_for_PT(pres, temp2);
		densM = 0.5 * (dens1 + dens2);

		for (int n = 0; n < 30; n++)
		{
			tempM = calc_T_for_DP(densM, pres);
			enthM = calc_H_for_DT(densM, tempM);
			error = enthM - enth;
			if (error < 0.0)
				dens1 = densM;
			else
				dens2 = densM;

			densM = 0.5 * (dens1 + dens2);
		}
	}

	double dens = densM;
	const double dRho = dens*1e-3;
	for (int n=0; n<100; n++)
	{
		double h1 = calc_H_for_DT(dens, calc_T_for_DP(dens, pres));
		double res = h1 - enth;
		if (fabs(res)<RTOL*enth)
			return dens;
		double h2 = calc_H_for_DT(dens + dRho, calc_T_for_DP(dens + dRho, pres));
		double drdrho = (h2-h1)/dRho;
		dens += -res / drdrho;
	}
	assert(false);
	return -1;
}

double calc_T_for_PH(double pres, double enth)
{
	if (pres < iaps84_pc())
	{
		// check two phase
		double hl = iaps84_sat_hl_p(pres);
		double hv = iaps84_sat_hv_p(pres);
		if (hl < enth && enth < hv)
			return iaps84_sat_T_p(pres);
	}

	double temp = 273.15 + 400; //TODO
	const double dT = 100*1e-3;
	for (int n=0; n<100; n++)
	{
		double h1 = calc_H_for_DT(calc_D_for_PT(pres, temp), temp);
		double res = h1 - enth;
		if (fabs(res)<RTOL*enth)
			return temp;
		double h2 = calc_H_for_DT(calc_D_for_PT(pres, temp + dT), temp + dT);
		double drdT = (h2-h1)/dT;
		temp += -res / drdT;
	}
	assert(false);
	return -1;
}


double iaps84_Tt()
{
	/*** Return the temperature[K] at triple point ***/

	return TP_TEMP;
}

double iaps84_pt()
{
	/*** Return the pressure[MPa] at triple point ***/

	return TP_PRES;
}

double iaps84_rholt()
{
	/*** Return the density[kg/m^3] of liquid phase at triple point ***/

	return 999.78;
}

double iaps84_rhovt()
{
	/*** Return the density[kg/m^3] of vapor phase at triple point ***/

	return 0.004855;
}

double iaps84_hlt()
{
	/*** Return the specific enthalpy[J/kg] of liquid phase at triple point ***/

	return 0.0;
}

double iaps84_hvt()
{
	/*** Return the specific enthalpy[J/kg] of vapor phase at triple point ***/

	return 2500500.0;
}


double iaps84_Tc()
{
	/*** Return the temperature[K] at critical point ***/

	return CP_TEMP;
}

double iaps84_pc()
{
	/*** Return the pressure[MPa] at critical point ***/

	return CP_PRES;
}

double iaps84_rhoc()
{
	/*** Return the density[kg/m^3] at critical point ***/

	return 322.0;
}

double iaps84_hc()
{
	/*** Return the specific enthalpy[J/kg] at critical point ***/

	return 2086000.0;
}


double iaps84_f_rhoT(double rhoS, double temp)
{
	/*** Return the Helmholtz function for given density[kg/m^3] and temperature[K] ***/

	if (TP_TEMP <= temp && temp <= MAX_TEMP)
	{
		return calc_A_for_DT(DCNV_KGM3_TO_GCM3 * rhoS, temp);
	}
	else
	{
		assert(false);
		return -1;
	}
}

double iaps84_g_rhoT(double rhoS, double temp)
{
	/*** Return the Gibbs function for given density[kg/m^3] and temperature[K] ***/

	if (TP_TEMP <= temp && temp <= MAX_TEMP)
	{
		return calc_G_for_DT(DCNV_KGM3_TO_GCM3 * rhoS, temp);
	}
	else
	{
		assert(false);
		return -1;
	}
}

double iaps84_u_rhoT(double rhoS, double temp)
{
	/*** Return the specific internal energy[J/kg] for given density[kg/m^3] and temperature[K] ***/

	if (TP_TEMP <= temp && temp <= MAX_TEMP)
	{
		return JCNV_KJKG_TO_JKG * calc_U_for_DT(DCNV_KGM3_TO_GCM3 * rhoS, temp);
	}
	else
	{
		assert(false);
		return -1;
	}
}

double iaps84_h_rhoT(double rhoS, double temp)
{
	/*** Return the specific enthalpy[J/kg] for given density[kg/m^3] and temperature[K] ***/

	if (TP_TEMP <= temp && temp <= MAX_TEMP)
	{
		return JCNV_KJKG_TO_JKG * calc_H_for_DT(DCNV_KGM3_TO_GCM3 * rhoS, temp);
	}
	else
	{
		assert(false);
		return -1;
	}
}

double iaps84_s_rhoT(double rhoS, double temp)
{
	/*** Return the specific entropy[J/(kg*K)] for given density[kg/m^3] and temperature[K] ***/

	if (TP_TEMP <= temp && temp <= MAX_TEMP)
	{
		return JCNV_KJKG_TO_JKG * calc_S_for_DT(DCNV_KGM3_TO_GCM3 * rhoS, temp);
	}
	else
	{
		assert(false);
		return -1;
	}
}

double iaps84_p_rhoT(double rhoS, double temp)
{
	/*** Return the pressure[MPa] for given density[kg/m^3] and temperature[K] ***/

	if (TP_TEMP <= temp && temp <= MAX_TEMP)
	{
		return calc_P_for_DT(DCNV_KGM3_TO_GCM3 * rhoS, temp);
	}
	else
	{
		assert(false);
		return -1;
	}
}

double iaps84_cp_rhoT(double rhoS, double temp)
{
	/*** Return the specific heat capacity[J/(kg*K)] at constant pressure for given density[kg/m^3] and temperature[K] ***/

	if (TP_TEMP <= temp && temp <= MAX_TEMP)
	{
		return JCNV_KJKG_TO_JKG * calc_Cp_for_DT(DCNV_KGM3_TO_GCM3 * rhoS, temp);
	}
	else
	{
		assert(false);
		return -1;
	}
}

double iaps84_cv_rhoT(double rhoS, double temp)
{
	/*** Return the specific heat capacity[J/(kg*K)] at constant volume for given density[kg/m^3] and temperature[K] ***/

	if (TP_TEMP <= temp && temp <= MAX_TEMP)
	{
		return JCNV_KJKG_TO_JKG * calc_Cv_for_DT(DCNV_KGM3_TO_GCM3 * rhoS, temp);
	}
	else
	{
		assert(false);
		return -1;
	}
}


double iaps84_rho_pT(double pres, double temp)
{
	/*** Return the density[kg/m^3] for given pressure[MPa] and temperature[K] ***/

	if (TP_PRES <= pres && pres <= MAX_PRES && TP_TEMP <= temp && temp <= MAX_TEMP)
	{
		return DCNV_GCM3_TO_KGM3 * calc_D_for_PT(pres, temp);
	}
	else
	{
		assert(false);
		return -1;
	}
}

double iaps84_rho_ph(double pres, double hJkg)
{
	/*** Return the density[kg/m^3] for given pressure[MPa] and specific enthalpy[J/kg] ***/

	if (TP_PRES <= pres && pres <= MAX_PRES)
	{
		return DCNV_GCM3_TO_KGM3 * calc_D_for_PH(pres, JCNV_JKG_TO_KJKG * hJkg);
	}
	else
	{
		assert(false);
		return -1;
	}
}

double iaps84_drho_dp_ph(double pres, double hJkg)
{
	/*** Return the drho/dp[(kg/m^3)/MPa] for given pressure[MPa] and specific enthalpy[J/kg] ***/

	static const double DELTA_PRES = 1.0e-7;

	if (TP_PRES <= pres && pres <= MAX_PRES)
	{
		double dens1 = DCNV_GCM3_TO_KGM3 * calc_D_for_PH(pres + DELTA_PRES, JCNV_JKG_TO_KJKG * hJkg);

		double dens2 = DCNV_GCM3_TO_KGM3 * calc_D_for_PH(pres, JCNV_JKG_TO_KJKG * hJkg);

		return (dens1 - dens2) / DELTA_PRES;
	}
	else
	{
		assert(false);
		return -1;
	}
}

double iaps84_drho_dh_ph(double pres, double hJkg)
{
	/*** Return the drho/dh[(kg/m^3)/(J/kg)] for given pressure[MPa] and specific enthalpy[J/kg] ***/

	static const double DELTA_ENTH = 1.0e-2;

	if (TP_PRES <= pres && pres <= MAX_PRES)
	{
		double dens1 = DCNV_GCM3_TO_KGM3 * calc_D_for_PH(pres, JCNV_JKG_TO_KJKG * hJkg + DELTA_ENTH);

		double dens2 = DCNV_GCM3_TO_KGM3 * calc_D_for_PH(pres, JCNV_JKG_TO_KJKG * hJkg);

		return (dens1 - dens1) / DELTA_ENTH;
	}
	else
	{
		assert(false);
		return -1;
	}
}


double iaps84_T_rhop(double rhoS, double pres)
{
	/*** Return the temperature[K] for given density[kg/m^3] and pressure[MPa] ***/

	if (TP_PRES <= pres && pres <= MAX_PRES)
	{
		return calc_T_for_DP(DCNV_KGM3_TO_GCM3 * rhoS, pres);
	}
	else
	{
		assert(false);
		return -1;
	}
}

double iaps84_T_ph(double pres, double hJkg)
{
	/*** Return the temperature[K] for given pressure[MPa] and specific enthalpy[J/kg] ***/

	if (TP_PRES <= pres && pres <= MAX_PRES)
	{
		return calc_T_for_PH(pres, JCNV_JKG_TO_KJKG * hJkg);
	}
	else
	{
		assert(false);
		return -1;
	}
}

double iaps84_dT_dp_ph(double pres, double hJkg)
{
	/*** Return the dT/dp[K/MPa] for given pressure[MPa] and specific enthalpy[J/kg] ***/

	static const double DELTA_PRES = 1.0e-7;

	if (TP_PRES <= pres && pres <= MAX_PRES)
	{
		double temp1 = calc_T_for_PH(pres + DELTA_PRES, JCNV_JKG_TO_KJKG * hJkg);

		double temp2 = calc_T_for_PH(pres, JCNV_JKG_TO_KJKG * hJkg);

		return (temp1 - temp2) / DELTA_PRES;
	}
	else
	{
		assert(false);
		return -1;
	}
}

double iaps84_dT_dh_ph(double pres, double hJkg)
{
	/*** Return the dT/dh[K/(J/kg)] for given pressure[MPa] and specific enthalpy[J/kg] ***/

	static const double DELTA_ENTH = 1.0e-2;

	if (TP_PRES <= pres && pres <= MAX_PRES)
	{
		double temp1 = calc_T_for_PH(pres, JCNV_JKG_TO_KJKG * hJkg + DELTA_ENTH);

		double temp2 = calc_T_for_PH(pres, JCNV_JKG_TO_KJKG * hJkg);

		return (temp1 - temp1) / DELTA_ENTH;
	}
	else
	{
		assert(false);
		return -1;
	}
}


double iaps84_viscosity_rhoT(double rhoS, double temp)
{
	/*** Return the viscosity[Pa*s] for given density[kg/m^3] and temperature[K] ***/
	/*** Density unit is kg/m^3 in this function ***/

	static const double REF_TEMP = 647.27, REF_DENS = 317.763;

	static const double AA[4] = { 0.0181583, 0.0177624, 0.0105287, -0.0036744 };

	static const double BB[6][5] = { {  0.501938,  0.235622,  -0.274637,  0.145831, -0.0270448 },
									 {  0.162888,  0.789393,  -0.743539,  0.263129, -0.0253093 },
									 { -0.130356,  0.673665,  -0.959456,  0.347247, -0.0267758 },
									 {  0.907919,  1.207552,  -0.687343,  0.213486, -0.0822904 },
									 { -0.551119,  0.0670665, -0.497089,  0.100754,  0.0602253 },
									 {  0.146543, -0.0843370,  0.195286, -0.032932, -0.0202595 } };

	if (TP_TEMP <= temp && temp <= MAX_TEMP)
	{
		double nTemp = temp / REF_TEMP;
		double rTemp = REF_TEMP / temp;
		double nDens = rhoS / REF_DENS;

		double sum = 0.0, mul1 = 0.0, mul2 = 0.0, v0 = 0.0, v1 = 0.0;


		sum = 0.0;
		for (int k = 0; k < 4; k++)
		{
			mul1 = 1.0;
			for (int n = 0; n < k; n++) mul1 *= rTemp;

			sum += AA[k] * mul1;
		}

		v0 = 1.0e-6 * sqrt(nTemp) / sum;


		sum = 0.0;
		for (int i = 0; i < 6; i++)
		{
			mul1 = 1.0;
			for (int n = 0; n < i; n++) mul1 *= rTemp - 1.0;

			for (int j = 0; j < 5; j++)
			{
				mul2 = 1.0;
				for (int n = 0; n < j; n++) mul2 *= nDens - 1.0;

				sum += BB[i][j] * mul1 * mul2;
			}
		}

		v1 = exp(nDens * sum);


		return v0 * v1;
	}
	else
	{
		assert(false);
		return -1;
	}
}

double iaps84_thermal_conductivity_rhoT(double rhoS, double temp)
{
	/*** Return the thermal conductivity[W/(K*m)] for given density[kg/m^3] and temperature[K] ***/
	/*** Density unit is kg/m^3 in this function ***/

	static const double REF_TEMP = 647.27, REF_DENS = 317.763, REF_PRES = 22.115;

	static const double PARAM_C = 3.7711e-8, OMEGA = 0.4678, PARAM_A = 18.66;

	static const double AA[4] = { 2.02223, 14.11166, 5.25597, -2.01870 };

	static const double BB[5][6] = { {  1.3293046, -0.40452437, 0.2440949,  0.018660751, -0.12961068,  0.044809953 },
									 {  1.7018363, -2.2156845,  1.6511057, -0.76736002,   0.37283344, -0.11203160  },
									 {  5.2246158, -10.124111,  4.9874687, -0.27297694,  -0.43083393,  0.13333849  },
									 {  8.7127675, -9.5000611,  4.3786606, -0.91783782,   0.0,         0.0         },
									 { -1.8525999,  0.9340469,  0.0,        0.0,          0.0,         0.0         } };


	if (TP_TEMP <= temp && temp <= MAX_TEMP)
	{
		double nTemp = temp / REF_TEMP;
		double rTemp = REF_TEMP / temp;
		double nDens = rhoS / REF_DENS;
		double rDens = REF_DENS / rhoS;

		double dpdt = 0.0, dddp = 0.0, visc = 0.0;

		double chi = 0.0, tcD = 0.0, tc0 = 0.0, tc1 = 0.0;

		double sum = 0.0, mul1 = 0.0, mul2 = 0.0, ww = 0.0, wt = 0.0, wd = 0.0;


		dpdt = (REF_TEMP / REF_PRES) * calc_dPdT_for_DT(DCNV_KGM3_TO_GCM3 * rhoS, temp);

		dddp = (REF_PRES / REF_DENS) / (calc_dPdD_for_DT(DCNV_KGM3_TO_GCM3 * rhoS, temp) / DCNV_GCM3_TO_KGM3);


		visc = iaps84_viscosity_rhoT(rhoS, temp);


		chi = nDens * dddp;

		ww = nTemp * rDens * dpdt;
		wt = rTemp - 1.0;
		wd = nDens - 1.0;

		tcD = PARAM_C / visc * ww * ww * pow(chi, OMEGA) * sqrt(nDens) * exp(-PARAM_A * wt * wt - wd * wd * wd * wd);


		sum = 0.0;

		for (int k = 0; k < 4; k++)
		{
			mul1 = 1.0;
			for (int n = 0; n < k; n++) mul1 *= rTemp;

			sum += AA[k] * mul1;
		}

		tc0 = sqrt(nTemp) / sum;


		sum = 0.0;

		for (int i = 0; i < 5; i++)
		{
			mul1 = 1.0;
			for (int n = 0; n < i; n++) mul1 *= rTemp - 1.0;

			for (int j = 0; j < 6; j++)
			{
				mul2 = 1.0;
				for (int n = 0; n < j; n++) mul2 *= nDens - 1.0;

				sum += BB[i][j] * mul1 * mul2;
			}
		}

		tc1 = exp(nDens * sum);


		return tc0 * tc1 + tcD;
	}
	else
	{
		assert(false);
		return -1;
	}
}


double iaps84_wv_ph(double pres, double hJkg)
{
	/*** Return the vapor quality[mass fraction] for given pressure[MPa] and specific enthalpy[J/kg] ***/

	if (TP_PRES <= pres && pres <= MAX_PRES)
	{
		if (pres >= iaps84_pc())
			return 0.0;

		double temp = 0.0, dl = 0.0, dv = 0.0, hl = 0.0, hv = 0.0;

		tcorr(&pres, &temp, &dl, &dv);

		hl = calc_H_for_DT(dl, temp);
		hv = calc_H_for_DT(dv, temp);

		return (JCNV_JKG_TO_KJKG * hJkg - hl) / (hv - hl);
	}
	else
	{
		assert(false);
		return -1;
	}
}

double iaps84_satv_ph(double pres, double hJkg)
{
	/*** Return the vapor saturation[volume fraction] for given pressure[MPa] and specific enthalpy[J/kg] ***/

	if (TP_PRES <= pres && pres <= MAX_PRES)
	{
		double temp = 0.0, dl = 0.0, dv = 0.0, hl = 0.0, hv = 0.0, wv = 0.0, vl = 0.0, vv = 0.0;

		tcorr(&pres, &temp, &dl, &dv);

		hl = calc_H_for_DT(dl, temp);
		hv = calc_H_for_DT(dv, temp);

		wv = (JCNV_JKG_TO_KJKG * hJkg - hl) / (hv - hl);

		vl = (1.0 - wv) / dl;
		vv = wv / dv;

		return vv / (vl + vv);
	}
	else
	{
		assert(false);
		return -1;
	}
}


double iaps84_sat_p_T(double temp)
{
	/*** Return the saturation pressure[MPa] for a given temperature[K] ***/

	if (TP_TEMP <= temp && temp < CP_TEMP)
	{
		double pres = 0.0, dl = 0.0, dv = 0.0;

		pcorr(temp, &pres, &dl, &dv);

		return pres;
	}
	else
	{
		assert(false);
		return -1.0;
	}
}

double iaps84_sat_T_p(double pres)
{
	/*** Return the saturation temperature[k] for a given pressure[MPa] ***/

	if (TP_PRES <= pres && pres < CP_PRES)
	{
		double temp = 0.0, dl = 0.0, dv = 0.0;

		tcorr(&pres, &temp, &dl, &dv);

		return temp;
	}
	else
	{
		assert(false);
		return -1.0;
	}
}

double iaps84_sat_rhol_T(double temp)
{
	/*** Return the liquid phase density[kg/m^3] for a given temperature[K] ***/

	if (TP_TEMP <= temp && temp < CP_TEMP)
	{
		double pres = 0.0, dl = 0.0, dv = 0.0;

		pcorr(temp, &pres, &dl, &dv);

		return DCNV_GCM3_TO_KGM3 * dl;
	}
	else
	{
		assert(false);
		return -1.0;
	}
}

double iaps84_sat_rhol_p(double pres)
{
	/*** Return the liquid phase density[kg/m^3] for a given pressure[MPa] ***/

	if (TP_PRES <= pres && pres < CP_PRES)
	{
		double temp = 0.0, dl = 0.0, dv = 0.0;

		tcorr(&pres, &temp, &dl, &dv);

		return DCNV_GCM3_TO_KGM3 * dl;
	}
	else
	{
		assert(false);
		return -1.0;
	}
}

double iaps84_sat_rhov_T(double temp)
{
	/*** Return the vapor phase density[kg/m^3] for a given temperature[K] ***/

	if (TP_TEMP <= temp && temp < CP_TEMP)
	{
		double pres = 0.0, dl = 0.0, dv = 0.0;

		pcorr(temp, &pres, &dl, &dv);

		return DCNV_GCM3_TO_KGM3 * dv;
	}
	else
	{
		assert(false);
		return -1.0;
	}
}

double iaps84_sat_rhov_p(double pres)
{
	/*** Return the vapor phase density[kg/m^3] for a given pressure[MPa] ***/

	if (TP_PRES <= pres && pres < CP_PRES)
	{
		double temp = 0.0, dl = 0.0, dv = 0.0;

		tcorr(&pres, &temp, &dl, &dv);

		return DCNV_GCM3_TO_KGM3 * dv;
	}
	else
	{
		assert(false);
		return -1.0;
	}
}

double iaps84_sat_hl_T(double temp)
{
	/*** Return the specific enthalpy[J/kg] of liquid phase for a given temperature[K] ***/

	if (TP_TEMP <= temp && temp < CP_TEMP)
	{
		double pres = 0.0, dl = 0.0, dv = 0.0;

		pcorr(temp, &pres, &dl, &dv);

		return JCNV_KJKG_TO_JKG * calc_H_for_DT(dl, temp);
	}
	else
	{
		assert(false);
		return -1.0;
	}
}

double iaps84_sat_hl_p(double pres)
{
	/*** Return the specific enthalpy[J/kg] of liquid phase for a given pressure[MPa] ***/

	if (TP_PRES <= pres && pres < CP_PRES)
	{
		double temp = 0.0, dl = 0.0, dv = 0.0;

		tcorr(&pres, &temp, &dl, &dv);

		return JCNV_KJKG_TO_JKG * calc_H_for_DT(dl, temp);
	}
	else
	{
		assert(false);
		return -1.0;
	}
}

double iaps84_sat_hv_T(double temp)
{
	/*** Return the specific enthalpy[J/kg] of vapor phase for a given temperature[K] ***/

	if (TP_TEMP <= temp && temp < CP_TEMP)
	{
		double pres = 0.0, dl = 0.0, dv = 0.0;

		pcorr(temp, &pres, &dl, &dv);

		return JCNV_KJKG_TO_JKG * calc_H_for_DT(dv, temp);
	}
	else
	{
		assert(false);
		return -1.0;
	}
}

double iaps84_sat_hv_p(double pres)
{
	/*** Return the specific enthalpy[J/kg] of vapor phase for a given pressure[MPa] ***/

	if (TP_PRES <= pres && pres < CP_PRES)
	{
		double temp = 0.0, dl = 0.0, dv = 0.0;

		tcorr(&pres, &temp, &dl, &dv);

		return JCNV_KJKG_TO_JKG * calc_H_for_DT(dv, temp);
	}
	else
	{
		assert(false);
		return -1.0;
	}
}

int iaps84_sat_prholv_T(double temp, double* pres, double* rhoL, double* rhoV)
{
	/*** Calculate pressure[MPa] and densities[kg/m^3] at saturation for a given temperature[K] ***/

	if (TP_TEMP <= temp && temp < CP_TEMP)
	{
		double dl = 0.0, dv = 0.0;

		pcorr(temp, pres, &dl, &dv);

		*rhoL = DCNV_GCM3_TO_KGM3 * dl;
		*rhoV = DCNV_GCM3_TO_KGM3 * dv;

		return 1;
	}
	else
	{
		assert(false);
		return 0;
	}
}

int iaps84_sat_Trholv_p(double pres, double* temp, double* rhoL, double* rhoV)
{
	/*** Calculate temperature[K] and densities[kg/m^3] at saturation for a given pressure[MPa] ***/

	if (TP_PRES <= pres && pres < CP_PRES)
	{
		double pp = pres, dl = 0.0, dv = 0.0;

		tcorr(&pp, temp, &dl, &dv);

		*rhoL = DCNV_GCM3_TO_KGM3 * dl;
		*rhoV = DCNV_GCM3_TO_KGM3 * dv;

		return 1;
	}
	else
	{
		assert(false);
		return 0;
	}
}
