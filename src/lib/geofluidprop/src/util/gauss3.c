
#include "gauss3.h"

#define REAL_AS_DOUBLE
#include "gauss3-template.h"
#undef REAL_AS_DOUBLE

#define REAL_AS_LONG_DOUBLE
#include "gauss3-template.h"
#undef REAL_AS_LONG_DOUBLE

#ifdef USE_QUAD
#define REAL_AS_QUAD
#include "gauss3-template.h"
#undef REAL_AS_QUAD
#endif
