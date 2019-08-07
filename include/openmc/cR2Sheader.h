#ifndef OPENMC_CR2SHEADER_H
#define OPENMC_CR2SHEADER_H

//**** External libraries ***************************************************//


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <sys/types.h>
#include <ctype.h>
#include <dirent.h>
#include <sys/timeb.h>
#include <signal.h>
#include <unistd.h>

#ifndef NO_GFX_MODE
#include <gd.h>
#endif

namespace openmc{

#ifdef DEBUG

void CheckPointer(char *, char *, long, long);
#else
#define CheckPointer(a, b, c, d)
#endif

//**** Function prototypes **************************************************//

void MCR2SSrc(long, double *, double *, double *, double *,  double *, double *,
 double *, double *, double *, long);
void ReadCDGS();


} // namespace openmc

#endif // OPENMC_CR2SHEADER_H