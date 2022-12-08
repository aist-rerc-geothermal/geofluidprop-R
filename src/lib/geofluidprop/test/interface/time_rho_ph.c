
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef _MSC_VER
#define POSIX
#endif

#ifdef POSIX
#include <sys/time.h>
#endif


#include "eos.h"

#include "util/utility.h"

int main(int argc, char **argv)
{
    int print_only_time = 0;
    if (argc > 2) {
        sscanf(argv[1], "%d", &print_only_time);
    }
    if (!print_only_time) {
        printf("-------------------------------------------------\n");
        printf("Measure run time\n");
        printf("-------------------------------------------------\n");
    }
    int eos_type = EOS_TYPE_WATER_IF97;
    if (argc > 1) {
        sscanf(argv[1], "%d", &eos_type);
    }
    char eosname[30];
    eos_get_eos_type_name(eos_type, eosname);
    if (!print_only_time)
    printf("-> EOS type = %s\n", eosname);
    EOS* eos = eos_create(eos_type);
    if (!eos)
        return 0;

    EOS_ARGS args;
    args.p = 0.992418352e-1*1e6;
    args.h = 112.652982e3;

#ifdef POSIX
    struct timeval s, e;
    gettimeofday(&s, NULL);
#endif

//    printf("rho=%g\n", eos_rho_ph(eos, &args));
    const int n = 10000;
    for (int i=0; i<n; i++)
//        eos_T_ph(eos, &args);
        eos_rho_ph(eos, &args);


#ifdef POSIX
    gettimeofday(&e, NULL);
    if (!print_only_time) printf("Elapsed time:");
    printf("%lf\n", (e.tv_sec - s.tv_sec) + (e.tv_usec - s.tv_usec)*1.0E-6);
#endif

    // finalize
    eos_free(eos);

    return 0;
}
