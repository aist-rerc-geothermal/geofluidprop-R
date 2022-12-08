
#include "IAPWS-Viscosity-08.h"

#include <math.h>

#include "util/utility.h"

#include "IAPWS-95.h"

// reference constants
#define Tr   647.096
#define pr   22.064e6
#define rhor 322.0
#define mur  1.e-6 // [Pa s]

// static const double Tr   = 647.096;
// static const double pr   = 22.064e6;
// static const double rhor = 322.0;
// static const double mur  = 1.e-6; // [Pa s]
// // // //
#define H vis2008_H
// H_ij
static const double vis2008_H[42] = {
     5.20094e-1,  2.22531e-1, -2.81378e-1,  1.61913e-1, -3.25372e-2,           0,           0,
     8.50895e-2,  9.99115e-1, -9.06851e-1,  2.57399e-1,           0,           0,           0,
       -1.08374,     1.88797, -7.72479e-1,           0,           0,           0,           0,
    -2.89555e-1,     1.26613, -4.89837e-1,           0,  6.98452e-2,           0, -4.35673e-3,
              0,           0, -2.57040e-1,           0,           0,  8.72102e-3,           0,
              0,  1.20573e-1,           0,           0,           0,           0, -5.93264e-4
};

#define xmu 0.068
#define qc 1./1.9
#define qd 1./1.1
#define v 0.630
#define gamma_ 1.239
#define xi0 0.13
#define Gamma0 0.06
#define TbarR 1.5
// static const double xmu = 0.068;
// static const double qc = 1./1.9;
// static const double qd = 1./1.1;
// static const double v = 0.630;
// static const double gamma_ = 1.239;
// static const double xi0 = 0.13;
// static const double Gamma0 = 0.06;
// static const double TbarR = 1.5;
// // // // // // // //

// viscosity in the dilute gas limit
double iapws08_vis_mu0bar(double Tbar)
{
    double a =   1.67752    / pow(Tbar,0)
               + 2.20462    / pow(Tbar,1)
               + 0.6366564  / pow(Tbar,2)
               - 0.241605   / pow(Tbar,3);
    return 100.*sqrt(Tbar)/a;
}

// contribution due to finite density
double iapws08_vis_mu1bar(double Tbar, double rhobar)
{
    double a = 0;
    for (int i=0; i<6; i++)
    {
        double b = 0;
        for (int j=0; j<7; j++)
            b += H[i*7+j] * pow(rhobar-1, j);
        a += pow(1/Tbar-1,i)*b;
    }
    a *= rhobar;
    return exp(a);
}

static double L(double xi, double w)
{
    if (qc*xi>1)
        return log((1+w)/(1-w));
    else
        return 2*atan(fabs(w));
}

static double sigma(double Tbar, double rhobar)
{
    // (drhobar/dpbar)_T = pr/rhor*(dp/drho)^-1
    return pr/rhor/iapws95_dp_drho_rhoT(rhobar*rhor, Tbar*Tr);
}

//critical enhancement
double iapws08_vis_mu2bar(double Tbar, double rhobar)
{
    double dchi = rhobar*(sigma(Tbar, rhobar) - sigma(TbarR,rhobar)*TbarR/Tbar);
    dchi = max(dchi, 0);
    double xi = xi0*pow(dchi/Gamma0, v/gamma_);

    double Y = 0;
    if (0<=xi && xi<=0.3817016416) {
        Y = 1./5*qc*xi*pow(qd*xi,5)*(1-qc*xi+pow(qc*xi,2)-765./504*pow(qd*xi,2));
    } else  if (xi>0.3817016416) {
        double psid = acos(pow(1+pow(qd*xi,2),-0.5));
        double qcxi2 = pow(qc*xi,2);
        double w = pow(fabs((qc*xi-1)/(qc*xi+1)),0.5)*tan(psid*0.5);
        Y = 1./12*sin(3*psid) -1/(4*qc*xi)*sin(2*psid)
            + 1./qcxi2*(1-5./4*qcxi2)*sin(psid)
            - 1./pow(qc*xi,3)*(
                (1-1.5*qcxi2)*psid - pow(fabs(qcxi2-1),1.5)*L(xi,w)
            );
    } else {
        //error
        //printf("***ERROR: this line should be called in mu2bar()\n");
        return -1;
    }


    return exp(xmu*Y);
}

double iapws08_viscosity_rhoT(double rho, double T)
{
    double Tbar = T/Tr;
    double rhobar = rho/rhor;
    return mur*iapws08_vis_mu0bar(Tbar)*iapws08_vis_mu1bar(Tbar,rhobar)*iapws08_vis_mu2bar(Tbar,rhobar);
}

double iapws08_viscosity_rhoT_simplified(double rho, double T)
{
    double Tbar = T/Tr;
    double rhobar = rho/rhor;
    if (645.91 < T && T < 650.77 && 245.8 < rho && rho < 405.3)
        return mur*iapws08_vis_mu0bar(Tbar)*iapws08_vis_mu1bar(Tbar,rhobar)*iapws08_vis_mu2bar(Tbar,rhobar);

    return mur*iapws08_vis_mu0bar(Tbar)*iapws08_vis_mu1bar(Tbar,rhobar);
}
