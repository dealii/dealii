// $Id$

#include <lac/forward-declarations.h>

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
