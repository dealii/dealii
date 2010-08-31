#include "../include/local_assemble_base.templates.h"

template class LocalAssembleBase<deal_II_dimension, DoFHandler<deal_II_dimension> >;
template class LocalAssembleBase<deal_II_dimension, MGDoFHandler<deal_II_dimension> >;

