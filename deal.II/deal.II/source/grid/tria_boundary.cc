/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */


#include <grid/tria_boundary.h>


// explicit instantiations
template class Boundary<deal_II_dimension>;
template class StraightBoundary<deal_II_dimension>;
template class HyperBallBoundary<deal_II_dimension>;
