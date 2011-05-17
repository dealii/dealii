#include "../include/error_handler.templates.h"
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/vector.h>

template class ErrorHandler<deal_II_dimension, Vector<double> >;
template class ErrorHandler<deal_II_dimension, BlockVector<double> >;
#ifdef CONTRIB_USE_PETSC
template class ErrorHandler<deal_II_dimension, PETScWrappers::Vector >;
#endif
//template class ErrorHandler<deal_II_dimension, PETScWrappers::BlockVector >;
