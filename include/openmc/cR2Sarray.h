#ifndef OPENMC_CR2SARRAY_H
#define OPENMC_CR2SARRAY_H

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

namespace openmc {
// Maximum array sizes //
#define MAX_STR                  256

///*** Global arrays and variables ******************************************//
//                                                                           //
//                                                                           //
// WDB    Writable main data block that contains nuclear data, geometry,     //
//        pointters and calculation parameters. Must be OpenMP protected     //
//        when written during threaded routines.                             //
//                                                                           //
// RDB    Pointer to WDB array, but type-cast to constant double. This is to //
//        induce a compiler warning when trying to write in the main data    //
//        at the wrong place.                                                //
//                                                                           //
// PRIVA  Data block for storing OpenMP private values. Divided into         //
//        segments and can be accessed by special routines without atomic or //
//        critical pragmas.                                                  //
//                                                                           //
// BUF    Buffer for storing cycle/batch wise data for statistics. Divided   //
//        into segments or accessed with atomic pragmas by special routines. //
//                                                                           //
// RES1   First results array, used for storing statistics. Not accessed by  //
//        OpenMP threads.                                                    //
//                                                                           //
// RES2   Second results array, containing large tables of integral data for //
//        which the statistics is not needed. Divided into segments or       //
//        accessed with atomic pragmas by special routines.                  //
//                                                                           //
// RES3   Same as RES2, but used for storing material-wise data for burnup   //
//        calculation when domain decomposition is in use. Always shared     //
//        between OpenMP threads.                                            //
//                                                                           //
// ASCII  Data array for storing text strings.                               //
//                                                                           //
// ACE    Data array for storing cross sections and nuclide data before      //
//        processing. Freed before transport cycle.                          //
//                                                                           //
//***************************************************************************//

// Arrays //

double *ACE;
double *WDB;
double *PRIVA;
double *BUF;
double *RES1;
double *RES2;
double *RES3;

const double *RDB;

char *ASCII;

// This size must be larger than cache line width //

#define RNG_SZ 100

unsigned long *SEED;
unsigned long *SEED0;
} //namespace openmc
#endif // OPENMC_CR2SARRAY_H