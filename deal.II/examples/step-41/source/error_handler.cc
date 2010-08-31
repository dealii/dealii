#include "../include/error_handler.templates.h"
#include <lac/block_vector.h>
#include <lac/petsc_vector.h>
#include <lac/vector.h>

template class ErrorHandler<deal_II_dimension, Vector<double> >;
template class ErrorHandler<deal_II_dimension, BlockVector<double> >;
#ifdef CONTRIB_USE_PETSC
template class ErrorHandler<deal_II_dimension, PETScWrappers::Vector >;
#endif
//template class ErrorHandler<deal_II_dimension, PETScWrappers::BlockVector >;
