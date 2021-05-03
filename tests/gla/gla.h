// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/generic_linear_algebra.h>

class LA_PETSc
{
public:
  class MPI
  {
  public:
    using Vector            = LinearAlgebraPETSc::MPI::Vector;
    using SparseMatrix      = LinearAlgebraPETSc::MPI::SparseMatrix;
    using BlockVector       = LinearAlgebraPETSc::MPI::BlockVector;
    using BlockSparseMatrix = LinearAlgebraPETSc::MPI::BlockSparseMatrix;
  };
};

class LA_Trilinos
{
public:
  class MPI
  {
  public:
    using Vector            = LinearAlgebraTrilinos::MPI::Vector;
    using SparseMatrix      = LinearAlgebraTrilinos::MPI::SparseMatrix;
    using BlockVector       = LinearAlgebraTrilinos::MPI::BlockVector;
    using BlockSparseMatrix = LinearAlgebraTrilinos::MPI::BlockSparseMatrix;
  };
};

class LA_Dummy
{
public:
  class MPI
  {
  public:
    class Vector
    {
    public:
      Vector()
      {}

      Vector(const IndexSet local, const MPI_Comm &comm)
      {}

      Vector(const IndexSet &local, const IndexSet &ghost, const MPI_Comm &comm)
      {}

      void
      reinit(const IndexSet local, const MPI_Comm &comm)
      {}

      void
      reinit(const IndexSet local, const IndexSet &ghost, const MPI_Comm &comm)
      {}

      void
      compress(VectorOperation::values op)
      {}

      bool
      all_zero()
      {
        return false;
      }

      bool
      has_ghost_elements()
      {
        return false;
      }

      unsigned int
      size()
      {
        return 0;
      }


      const Vector &
      operator=(const double number)
      {
        return *this;
      }


      const Vector &
      operator*=(const double factor)
      {
        return *this;
      }

      double &
      operator()(unsigned int)
      {
        static double d;
        return d;
      }

      const double &
      operator()(unsigned int) const
      {
        static double d;
        return d;
      }
    };

    class SparseMatrix
    {
    public:
      template <typename SP>
      SparseMatrix(const IndexSet &local,
                   const IndexSet &,
                   SP &            sp,
                   const MPI_Comm &comm = MPI_COMM_WORLD)
      {}

      void
      set(unsigned int, unsigned int, double);

      const double &
      operator()(unsigned int, unsigned int) const
      {
        static double d;
        return d;
      }
    };
  };
};
