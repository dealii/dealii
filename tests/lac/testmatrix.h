// $Id$

#include <lac/forward-declarations.h>
#include <base/exceptions.h>
#include <lac/mgbase.h>

/**
 * Finite difference matrix on uniform grid.
 * Generator for simple 5-point discretization of Laplace problem.
 */

class FDMatrix
{
  public:
				     /**
				      * Constructor specifying grid resolution.
				      */
    FDMatrix(unsigned int nx, unsigned int ny);
    
				     /**
				      * Generate the matrix structure.
				      */
    void build_structure(SparseMatrixStruct& structure) const;
    
				     /**
				      * Fill the matrix with values.
				      */
    template <typename number>
    void laplacian(SparseMatrix<number>&) const;

    template <typename number>
    void gnuplot_print(ostream&, const Vector<number>&) const;
    
  private:
				     /**
				      * Number of gridpoints in x-direction.
				      */
    unsigned int nx;

				     /**
				      * Number of gridpoints in y-direction.
				      */
    unsigned int ny;
};

/**
 * Grid transfer for finite differences on uniform grid.
 */

class FDMGTransfer
  :
  public MGTransferBase
{
  public:
				     /**
				      * Constructor. Prepares grid
				      * transfer matrices for #nlevels#
				      * levels on an #nx# times #ny#
				      * grid.
				      */
    FDMGTransfer(unsigned int nx, unsigned int ny,
		 unsigned int nlevels);

				     /**
				      * Implementation of abstract
				      * function in #MGTranferBase#.
				      */
    virtual void prolongate (const unsigned int   to_level,
			     Vector<double>       &dst,
			     const Vector<double> &src) const;

				     /**
				      * Implementation of abstract
				      * function in #MGTranferBase#.
				      */
    virtual void restrict (const unsigned int   from_level,
			   Vector<double>       &dst,
			   const Vector<double> &src) const;

				     /**
				      * Exception.
				      */
    DeclException2(ExcDivide, unsigned int, unsigned int,
		   << "Cannot divide " << arg1 << " by " << arg2);
    
  private:
				     /**
				      * Prolongation matrix structures.
				      */
    vector<SparseMatrixStruct > structures;
				     /**
				      * Prolongation matrices.
				      */
    vector<SparseMatrix<double> > matrices;

				     /**
				      * Matrix generator.
				      * The arguments #nx# and #ny#
				      * are the numbers on the COARSE level.
				      * Builds a transfer matrix from
				      * fine to coarse (#vmult#).
				      */
    void build_matrix(unsigned int nx, unsigned int ny,
		      SparseMatrixStruct&, SparseMatrix<double>&);
};
