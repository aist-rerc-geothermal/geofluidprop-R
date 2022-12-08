
#include <stdio.h>
#include <math.h>

#include "model/driesner07/Driesner2007_H2ONaCl.h"
#include "model/klyukinetal17/KlyukinEtAl2017.h"

#include "util/utility.h"
#include "utest.h"

double to_mole_fraction(double mass_frac);
double molality_to_mole_fraction(double molality);
double to_mass_frac(double mole_frac);

UTEST(klyukinetal17, test1)
{
    EXPECT_NEAR0(1137.56861587956, klyukinetal2017_viscosity(999.102002811697, 15+273.15, 0)*1e6, 1e-2);
    EXPECT_NEAR0(141.088472897387, klyukinetal2017_viscosity(878.877573449578, 200+273.15, 1.*1e-2)*1e6, 3e-2);
    EXPECT_NEAR0(73.9318857035372, klyukinetal2017_viscosity(621.255502889116, 350+273.15, 1.*1e-2)*1e6, 4e-2);
    EXPECT_NEAR0(93.9076585322679, klyukinetal2017_viscosity(751.534089680878, 350+273.15, 10.*1e-2)*1e6, 3e-2);
}
