
#include "../tests.h"
#include <deal.II/base/exceptions.h>
#include <iostream>
#include <deal.II/lac/petsc_matrix_free.h>
#include <deal.II/lac/vector.h>

// A variant of the tests/lac/testmatrix.h file that uses the PETSc
// matrix free interface


/**
 * Finite difference PETSc matrix-free object on uniform grid.
 * Generator for simple 5-point discretization of Laplace problem.
 */

class PetscFDMatrix : public dealii::PETScWrappers::MatrixFree
{
public:
  /**
   * Constructor specifying grid resolution.
   */
  PetscFDMatrix(unsigned int size, unsigned int dim);

  /**
   * Matrix-vector multiplication:
   * let <i>dst = M*src</i> with
   * <i>M</i> being this matrix.
   *
   * Source and destination must
   * not be the same vector.
   */
  void vmult (dealii::PETScWrappers::VectorBase       &dst,
              const dealii::PETScWrappers::VectorBase &src) const;

  /**
   * Matrix-vector multiplication: let
   * <i>dst = M<sup>T</sup>*src</i> with
   * <i>M</i> being this matrix. This
   * function does the same as @p vmult()
   * but takes the transposed matrix.
   *
   * Source and destination must
   * not be the same vector.
   */
  void Tvmult (dealii::PETScWrappers::VectorBase       &dst,
               const dealii::PETScWrappers::VectorBase &src) const;

  /**
   * Adding Matrix-vector
   * multiplication. Add
   * <i>M*src</i> on <i>dst</i>
   * with <i>M</i> being this
   * matrix.
   *
   * Source and destination must
   * not be the same vector.
   */
  void vmult_add (dealii::PETScWrappers::VectorBase       &dst,
                  const dealii::PETScWrappers::VectorBase &src) const;

  /**
   * Adding Matrix-vector
   * multiplication. Add
   * <i>M<sup>T</sup>*src</i> to
   * <i>dst</i> with <i>M</i> being
   * this matrix. This function
   * does the same as @p vmult_add()
   * but takes the transposed
   * matrix.
   *
   * Source and destination must
   * not be the same vector.
   */
  void Tvmult_add (dealii::PETScWrappers::VectorBase       &dst,
                   const dealii::PETScWrappers::VectorBase &src) const;

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


// --------------- inline and template functions -----------------

inline
PetscFDMatrix::PetscFDMatrix(unsigned int size, unsigned int dim)
  : PETScWrappers::MatrixFree (dim, dim, dim, dim),
    nx (size), ny (size)
{}



inline
void
PetscFDMatrix::vmult_add (dealii::PETScWrappers::VectorBase       &dst,
                          const dealii::PETScWrappers::VectorBase &src) const
{
  for (unsigned int i=0; i<=ny-2; i++)
    {
      for (unsigned int j=0; j<=nx-2; j++)
        {
          // Number of the row to be entered
          unsigned int row = j+(nx-1)*i;

          dst(row) += 4. * src(row);              // A.set(row, row, 4.);

          if (j>0)
            {
              dst(row-1) += -1. * src(row);       // A.set(row-1, row, -1.);
              dst(row) += -1. * src(row-1);       // A.set(row, row-1, -1.);
            }
          if (i>0)
            {
              dst(row-(nx-1)) += -1. * src(row);  // A.set(row-(nx-1), row, -1.);
              dst(row) += -1. * src(row-(nx-1));  // A.set(row, row-(nx-1), -1.);
            }
        }
    }
}



inline
void
PetscFDMatrix::vmult (dealii::PETScWrappers::VectorBase       &dst,
                      const dealii::PETScWrappers::VectorBase &src) const
{
  dst = 0;
  vmult_add (dst, src);
}



inline
void
PetscFDMatrix::Tvmult (dealii::PETScWrappers::VectorBase       &dst,
                       const dealii::PETScWrappers::VectorBase &src) const
{
  dst = 0;
  vmult_add (dst, src);
}



inline
void
PetscFDMatrix::Tvmult_add (dealii::PETScWrappers::VectorBase       &dst,
                           const dealii::PETScWrappers::VectorBase &src) const
{
  vmult_add (dst, src);
}


