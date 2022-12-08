
#include <stdlib.h>
#include "utest.h"

#include "util/gauss3.h"

UTEST(gauss3, double)
{
    const double c = 1.2345678901234567890;
    double A[9];
    srand(0);
    for (int i=0; i<9; i++)
        A[i] = (double)rand()/(double)RAND_MAX;
    double x[3];
    for (int i=0; i<3; i++)
        x[i] = c * pow(0.1, i);
    double b[3];
    for (int i=0; i<3; i++)
        b[i] = A[i*3+0]*x[0] + A[i*3+1]*x[1] + A[i*3+2]*x[2];
    gauss3(A, b, x);
    EXPECT_NEAR(c,      x[0]);
    EXPECT_NEAR(c*0.1,  x[1]);
    EXPECT_NEAR(c*0.01, x[2]);
    for (int i=0; i<3; i++) {
        double e = c * pow(0.1, i);
        // printf("expected: %.20e\n", e);
        // printf("actual  : %.20e\n", x[i]);
        printf("error   : %.20e\n", fabs(x[i]-e));
    }
}

#ifdef USE_LONGDOUBLE
UTEST(gauss3, long)
{
    const double c = 1.2345678901234567890;
    long double A[9];
    srand(0);
    for (int i=0; i<9; i++)
        A[i] = (double)rand()/(double)RAND_MAX;
    long double x[3];
    for (int i=0; i<3; i++)
        x[i] = c * pow(0.1, i);
    long double b[3];
    for (int i=0; i<3; i++)
        b[i] = A[i*3+0]*x[0] + A[i*3+1]*x[1] + A[i*3+2]*x[2];
    gauss3_long(A, b, x);
    EXPECT_NEAR(c,      x[0]);
    EXPECT_NEAR(c*0.1,  x[1]);
    EXPECT_NEAR(c*0.01, x[2]);
    for (int i=0; i<3; i++) {
        double e = c * pow(0.1, i);
        // printf("expected: %.20e\n", e);
        // printf("actual  : %.20e\n", x[i]);
        printf("error   : %.20e\n", fabs(x[i]-e));
    }
}
#endif

UTEST(cramer3, double)
{
    const double c = 1.2345678901234567890;
    double A[9];
    srand(0);
    for (int i=0; i<9; i++)
        A[i] = (double)rand()/(double)RAND_MAX;
    double x[3];
    for (int i=0; i<3; i++)
        x[i] = c * pow(0.1, i);
    double b[3];
    for (int i=0; i<3; i++)
        b[i] = A[i*3+0]*x[0] + A[i*3+1]*x[1] + A[i*3+2]*x[2];
    cramer3(A, b, x);
    EXPECT_NEAR(c,      x[0]);
    EXPECT_NEAR(c*0.1,  x[1]);
    EXPECT_NEAR(c*0.01, x[2]);
    for (int i=0; i<3; i++) {
        double e = c * pow(0.1, i);
        // printf("expected: %.20e\n", e);
        // printf("actual  : %.20e\n", x[i]);
        printf("error   : %.20e\n", fabs(x[i]-e));
    }
}

#ifdef USE_LONGDOUBLE
UTEST(cramer3, long)
{
    const double c = 1.2345678901234567890;
    long double A[9];
    srand(0);
    for (int i=0; i<9; i++)
        A[i] = (double)rand()/(double)RAND_MAX;
    long double x[3];
    for (int i=0; i<3; i++)
        x[i] = c * pow(0.1, i);
    long double b[3];
    for (int i=0; i<3; i++)
        b[i] = A[i*3+0]*x[0] + A[i*3+1]*x[1] + A[i*3+2]*x[2];
    cramer3_long(A, b, x);
    EXPECT_NEAR(c,      x[0]);
    EXPECT_NEAR(c*0.1,  x[1]);
    EXPECT_NEAR(c*0.01, x[2]);
    for (int i=0; i<3; i++) {
        double e = c * pow(0.1, i);
        // printf("expected: %.20e\n", e);
        // printf("actual  : %.20e\n", x[i]);
        printf("error   : %.20e\n", fabs(x[i]-e));
    }
}
#endif
