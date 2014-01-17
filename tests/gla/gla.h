
#include <deal.II/lac/generic_linear_algebra.h>

class LA_PETSc
{
public:
  class MPI
  {
  public:
    typedef LinearAlgebraPETSc::MPI::Vector Vector;
    typedef LinearAlgebraPETSc::MPI::SparseMatrix SparseMatrix;
    typedef LinearAlgebraPETSc::MPI::BlockVector BlockVector;
    typedef LinearAlgebraPETSc::MPI::BlockSparseMatrix BlockSparseMatrix;
  };
};

class LA_Trilinos
{
public:
  class MPI
  {
  public:
    typedef LinearAlgebraTrilinos::MPI::Vector Vector;
    typedef LinearAlgebraTrilinos::MPI::SparseMatrix SparseMatrix;
    typedef LinearAlgebraTrilinos::MPI::BlockVector BlockVector;
    typedef LinearAlgebraTrilinos::MPI::BlockSparseMatrix BlockSparseMatrix;
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

      void reinit(const IndexSet local, const MPI_Comm &comm)
      {}

      void reinit(const IndexSet local, const IndexSet &ghost, const MPI_Comm &comm)
      {}

      void compress(VectorOperation::values op)
      {}

	bool all_zero()
	  {
	    return false;
	  }	
	
      bool has_ghost_elements()
      {
        return false;
      }

      unsigned int size()
      {
        return 0;
      }


      const Vector &operator=(const double number)
      {
        return *this;
      }


      const Vector &operator*=(const double factor)
      {
        return *this;
      }

      double &operator()(unsigned int)
      {
        static double d;
        return d;
      }

      const double &operator()(unsigned int) const
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
                   SP &sp,
                   const MPI_Comm &comm=MPI_COMM_WORLD)
      {}

      void set(unsigned int, unsigned int, double);

      const double &operator()(unsigned int, unsigned int) const
      {
        static double d;
        return d;
      }

    };


  };
};
