// testmatrix.h,v 1.15 2002/10/03 22:41:22 wolf Exp

#include "../tests.h"
#include <base/exceptions.h>
#include <iostream>

template <typename number> class Vector;
class SparsityPattern;

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
    void five_point_structure(SparsityPattern& structure) const;
    
				     /**
				      * Generate the matrix structure.
				      */
    void nine_point_structure(SparsityPattern& structure) const;
    
				     /**
				      * Fill the matrix with values.
				      */
    template <typename MATRIX>
    void five_point(MATRIX&, bool nonsymmetric = false) const;

				     /**
				      * Fill the matrix with values.
				      */
    template <typename MATRIX>
    void nine_point(MATRIX&, bool nonsymmetric = false) const;

				     /**
				      * Fill the matrix with values.
				      */
    template <typename MATRIX>
    void upwind(MATRIX&, bool back = false) const;

    template <typename number>
    void gnuplot_print(std::ostream&, const Vector<number>&) const;
    
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
