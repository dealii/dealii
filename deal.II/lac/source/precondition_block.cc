/*----------------------------   precondition_block.cc     ---------------------------*/
/*      $Id$                 */
/*                Ralf Hartmann, University of Heidelberg                            */



#include <lac/precondition_block.templates.h>




// explicit instantiations for "float" PreconditionBlock
template class PreconditionBlock<float, float>;

// the instantiation for class PreconditionBlock<float, double> is skipped
// because it does not make sense to have inverse block matrices with
// higher precision than the matrix itself


// explicit instantiations for "double" PreconditionBlock
template class PreconditionBlock<double, float>;

template class PreconditionBlock<double, double>;


/*--------------------- PreconditionBlockGaussSeidel -----------------------*/


// explicit instantiations for "float" PreconditionBlock
template class PreconditionBlockSOR<float, float>;

template void PreconditionBlockSOR<float, float>::operator() (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<float, float>::operator() (
  Vector<double> &, const Vector<double> &) const;


// the instantiation for class PreconditionBlockSOR<float, double> is skipped
// because it does not make sense to have inverse block matrices with
// higher precision than the matrix itself



// explicit instantiations for "double" PreconditionBlockSOR
template class PreconditionBlockSOR<double, float>;


template void PreconditionBlockSOR<double, float>::operator() (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<double, float>::operator() (
  Vector<double> &, const Vector<double> &) const;

template class PreconditionBlockSOR<double, double>;

template void PreconditionBlockSOR<double, double>::operator() (
  Vector<float> &, const Vector<float> &) const;
template void PreconditionBlockSOR<double, double>::operator() (
  Vector<double> &, const Vector<double> &) const;






/*----------------------------   precondition_block.cc     ---------------------------*/
