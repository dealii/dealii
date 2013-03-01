
#include <deal.II/lac/abstract_linear_algebra.h>

class LA_PETSc
{
  public:
    class MPI
    {
      public:
	typedef LinearAlgebraPETSc::MPI::Vector Vector;
    };
};

class LA_Trilinos
{
  public:
    class MPI
    {
      public:
	typedef LinearAlgebraPETSc::MPI::Vector Vector;
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
	    
	    Vector(const IndexSet local, const MPI_Comm &comm=MPI_COMM_WORLD)
	      {}
	    
	    Vector(const IndexSet &local, const IndexSet &ghost, const MPI_Comm &comm=MPI_COMM_WORLD) 
	      {}
	    
	    Vector(const MPI_Comm &comm, const IndexSet local) 
	      {}
	    
	    Vector(const MPI_Comm &commm, const IndexSet &local, const IndexSet &ghost)
	      {}
	    
	    void compress(VectorOperation::values op)
	      {}
	    
	    bool has_ghost_elements()
	      {
		return false;}

	    unsigned int size()
	      {return 0;}
	    
	    
	    const Vector & operator=(const double number)
	      {
		return *this;
	      }
	    

	      const Vector & operator*=(const double factor)
	      {
		return *this;
	      }
	    
	    double & operator()(unsigned int)
	      {
		static double d;
		return d;
	      }
	    
	    const double & operator()(unsigned int) const
	      {
		static double d;
		return d;
	      }
	    
	    
	};
	
	
    };
};
