// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/block_matrix_array.h>
#include <deal.II/multigrid/mg_block_smoother.h>
#include <fstream>
#include <iomanip>
#include <iomanip>
#include <algorithm>

using namespace std;

template <typename number>
class ScalingMatrix : public Subscriptor
{
public:
  /**
   * Constructor setting the
   * scaling factor. Default is
   * constructing the identity
   * matrix.
   */
  ScalingMatrix(number scaling_factor = 1.);
  /**
  * Apply preconditioner.
  */
  template<class VECTOR>
  void vmult (VECTOR &, const VECTOR &) const;

  /**
   * Apply transpose
   * preconditioner. Since this is
   * the identity, this function is
   * the same as
   * vmult().
   */
  template<class VECTOR>
  void Tvmult (VECTOR &, const VECTOR &) const;
  /**
   * Apply preconditioner, adding to the previous value.
   */
  template<class VECTOR>
  void vmult_add (VECTOR &, const VECTOR &) const;

  /**
   * Apply transpose
   * preconditioner, adding. Since this is
   * the identity, this function is
   * the same as
   * vmult_add().
   */
  template<class VECTOR>
  void Tvmult_add (VECTOR &, const VECTOR &) const;

private:
  number factor;
};


//----------------------------------------------------------------------//

template<typename number>
ScalingMatrix<number>::ScalingMatrix(number factor)
  :
  factor(factor)
{}


template<typename number>
template<class VECTOR>
inline void
ScalingMatrix<number>::vmult (VECTOR &dst, const VECTOR &src) const
{
  dst.equ(factor, src);
}

template<typename number>
template<class VECTOR>
inline void
ScalingMatrix<number>::Tvmult (VECTOR &dst, const VECTOR &src) const
{
  dst.equ(factor, src);
}

template<typename number>
template<class VECTOR>
inline void
ScalingMatrix<number>::vmult_add (VECTOR &dst, const VECTOR &src) const
{
  dst.add(factor, src);
}



template<typename number>
template<class VECTOR>
inline void
ScalingMatrix<number>::Tvmult_add (VECTOR &dst, const VECTOR &src) const
{
  dst.add(factor, src);
}

//----------------------------------------------------------------------//

template<class MATRIX, class RELAX>
void check_smoother(const MGLevelObject<MATRIX> &m,
                    const MGLevelObject<RELAX> &r)
{
  GrowingVectorMemory<BlockVector<double> > mem;
  MGSmootherBlock<MATRIX, RELAX, double> smoother(mem);

  smoother.initialize(m, r);

  for (unsigned int l=m.min_level(); l<= m.max_level(); ++l)
    {
      deallog << "Level " << l << std::endl;

      BlockVector<double> &u = *mem.alloc();
      BlockVector<double> &f = *mem.alloc();
      u.reinit(m[l].n_block_rows(), 3);
      f.reinit(u);
      for (unsigned int b=0; b<f.n_blocks(); ++b)
        for (unsigned int i=0; i<f.block(b).size(); ++i)
          f.block(b)(i) = (b+1)*(i+l);

      deallog << "First step" << std::endl;
      smoother.set_steps(1);
      smoother.smooth(l, u, f);

      for (unsigned int b=0; b<u.n_blocks(); ++b)
        {
          for (unsigned int i=0; i<u.block(b).size(); ++i)
            deallog << '\t' << (int) (u.block(b)(i)+.5);
          deallog << std::endl;
        }

      deallog << "Second step" << std::endl;
      smoother.smooth(l, u, f);

      for (unsigned int b=0; b<u.n_blocks(); ++b)
        {
          for (unsigned int i=0; i<u.block(b).size(); ++i)
            deallog << '\t' << (int) (u.block(b)(i)+.5);
          deallog << std::endl;
        }

      deallog << "Two steps" << std::endl;
      u = 0.;
      smoother.set_steps(2);
      smoother.smooth(l, u, f);

      for (unsigned int b=0; b<u.n_blocks(); ++b)
        {
          for (unsigned int i=0; i<u.block(b).size(); ++i)
            deallog << '\t' << (int) (u.block(b)(i)+.5);
          deallog << std::endl;
        }

      mem.free(&u);
      mem.free(&f);
    }
}

void check()
{
  ScalingMatrix<double> s1(-1.);
  ScalingMatrix<double> s2(2.);
  ScalingMatrix<double> s8(8.);

  GrowingVectorMemory<Vector<double> > mem;
  MGLevelObject<BlockMatrixArray<double> > A (2,4);
  MGLevelObject<BlockTrianglePrecondition<double> > P(2,4);

  for (unsigned int l=A.min_level(); l<= A.max_level(); ++l)
    {
      A[l].initialize(3, 3);
      P[l].reinit(3);
      for (unsigned int b=0; b<A[l].n_block_rows(); ++b)
        {
          P[l].enter(s2, b, b, A[l].n_block_rows()-b);
          A[l].enter(s8, b, b, 1);
          for (unsigned int b2=0; b2<A[l].n_block_rows(); ++b2)
            A[l].enter(s1,b,b2,1.);
        }
    }

  check_smoother(A, P);
}


int main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check();
}

