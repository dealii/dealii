// $Id$

#ifndef __base_types_h
#define __base_types_h

#ifndef _stdlib_h
#include <stdlib.h>
#endif

/// Maximum Number of boundary conditions (divided by 8 )

#define BOUNDARY_CONDITION_BYTES 1

/// Maximum number of multi purpose flags per cell ( / 8 )

#define CELL_FLAG_BYTES 1

/// Maximum number of multi purpose flags per vertex ( / 8 )

#define VERTEX_FLAG_BYTES 1

#define UNSIGNED unsigned long

typedef double COORD;
#define DOUBLECOORD

typedef int SHORT;

#define ABS(x) ( ((x)>=0) ? (x) : -(x) )
#define MAX(a,b) ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b) ( ((a)<(b)) ? (a) : (b) )
#define SIGN(a) ( ((a)<0) ? (-1) : 1 )
#ifndef NULL
#define NULL 0
#endif

#endif
