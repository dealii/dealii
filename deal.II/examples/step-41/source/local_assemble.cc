#include "../include/local_assemble.templates.h"
#include <dofs/dof_handler.h>
#include <multigrid/mg_dof_handler.h>

template class LocalAssemble<deal_II_dimension, DoFHandler<deal_II_dimension> >;
template class LocalAssemble<deal_II_dimension, MGDoFHandler<deal_II_dimension> >;
