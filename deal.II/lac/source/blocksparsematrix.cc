/*----------------------------   blocksparsematrix.cc     ---------------------------*/
/*      $Id$                 */
/*                Ralf Hartmann, University of Heidelberg                            */
/*----------------------------   blocksparsematrix.cc     ---------------------------*/



#include <lac/blocksparsematrix.h>
#include <lac/blocksparsematrix.templates.h>
//#include <lac/vector.h>


template <typename number, typename blocknumber>
bool datatype_compatible ();


// explicit instantiations for "float" BlockSparseMatrix
template class BlockSparseMatrix<float, float>;

template void BlockSparseMatrix<float, float>::precondition_BlockSOR (
  Vector<float> &, const Vector<float> &, float) const;
template void BlockSparseMatrix<float, float>::precondition_BlockSOR (
  Vector<double> &, const Vector<double> &, float) const;

// the instantiation for class BlockSparseMatrix<float, double> is skipped
// because it does not make sence to have inverse block matrices with
// higher precision than the matrix itself

// explicit instantiations for "double" BlockSparseMatrix
template class BlockSparseMatrix<double, float>;

template void BlockSparseMatrix<double, float>::precondition_BlockSOR (
  Vector<float> &, const Vector<float> &, double) const;
template void BlockSparseMatrix<double, float>::precondition_BlockSOR (
  Vector<double> &, const Vector<double> &, double) const;

template class BlockSparseMatrix<double, double>;

template void BlockSparseMatrix<double, double>::precondition_BlockSOR (
  Vector<float> &, const Vector<float> &, double) const;
template void BlockSparseMatrix<double, double>::precondition_BlockSOR (
  Vector<double> &, const Vector<double> &, double) const;




/*----------------------------   blocksparsematrix.cc     ---------------------------*/
