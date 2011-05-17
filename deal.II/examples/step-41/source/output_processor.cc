#include "../include/output_processor.templates.h"
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/vector.h>

template <>
void OutputProcessor<deal_II_dimension,Vector<double> >::dump_vector (const Vector<double> &rhs, 
								    const std::string &filename )
{
    deallog.push("Dump");
    deallog << "Writing: " << filename << std::endl;
    std::ofstream out ( (filename).c_str());
    rhs.block_write(out);
    deallog.pop();
}

template class FilteredDataOut<deal_II_dimension, DoFHandler<deal_II_dimension> >;

template class OutputProcessor<deal_II_dimension, Vector<double> >;
template class OutputProcessor<deal_II_dimension, BlockVector<double> >;

#ifdef DEAL_II_USE_PETSC
template class OutputProcessor<deal_II_dimension, PETScWrappers::Vector >;
#endif

// template class OutputProcessor<deal_II_dimension, PETScWrappers::BlockVector >;
