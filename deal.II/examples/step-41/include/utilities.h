#ifndef UTILITIES_HLT
#define UTILITIES_HLT

#include <base/utilities.h>
#include <grid/grid_tools.h>
#include <lac/vector_memory.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/block_vector.h>
#include <lac/block_sparse_matrix.h>

using namespace dealii;

template <class Matrix>
class InverseMatrix : public Subscriptor
{
 public:
  InverseMatrix (const Matrix &m);
  
  void vmult (Vector<double>       &dst,
	      const Vector<double> &src) const;
  
  void Tvmult (Vector<double>       &dst,
	       const Vector<double> &src) const;
 
 private:
  const SmartPointer<const Matrix> matrix;
 
  mutable GrowingVectorMemory<> vector_memory;    
};
 
  
class BBt : public Subscriptor
{
 public:
  void reinit (const BlockSparseMatrix<double> &A);
  
  
  void vmult (Vector<double>       &dst,
              const Vector<double> &src) const;

  void Tvmult (Vector<double>       &dst,
              const Vector<double> &src) const;

 private:
  SmartPointer<const BlockSparseMatrix<double> > system_matrix;

  mutable BlockVector<double> u, v;

  unsigned int dim;
};

template <typename TYPE>
void smart_delete (SmartPointer<TYPE> &sp) {
  if(sp) {
    TYPE * p = sp;
    sp = 0;
    delete p;
  }
}

#endif
