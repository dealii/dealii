
#include "../tests.h"
#include <deal.II/base/exceptions.h>
#include <iostream>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/vector.h>


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
  template <typename SP>
  void five_point_structure(SP &structure) const;

  /**
   * Generate the matrix structure.
   */
  template <typename SP>
  void nine_point_structure(SP &structure) const;

  /**
   * Fill the matrix with values.
   */
  template <typename MATRIX>
  void five_point(MATRIX &, bool nonsymmetric = false) const;

  /**
   * Fill the matrix with values.
   */
  template <typename MATRIX>
  void nine_point(MATRIX &, bool nonsymmetric = false) const;

  /**
   * Fill the matrix with values.
   */
  template <typename MATRIX>
  void upwind(MATRIX &, bool back = false) const;

  template <typename number>
  void gnuplot_print(std::ostream &, const Vector<number> &) const;

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
FDMatrix::FDMatrix(unsigned int nx, unsigned int ny)
  :
  nx(nx), ny(ny)
{}



template <typename SP>
inline
void
FDMatrix::five_point_structure(SP &structure) const
{
  for (unsigned int i=0; i<=ny-2; i++)
    {
      for (unsigned int j=0; j<=nx-2; j++)
        {
          // Number of the row to be entered
          unsigned int row = j+(nx-1)*i;
          structure.add(row, row);
          if (j>0)
            {
              structure.add(row-1, row);
              structure.add(row, row-1);
            }
          if (j<nx-2)
            {
              structure.add(row+1, row);
              structure.add(row, row+1);
            }
          if (i>0)
            {
              structure.add(row-(nx-1), row);
              structure.add(row, row-(nx-1));
            }
          if (i<ny-2)
            {
              structure.add(row+(nx-1), row);
              structure.add(row, row+(nx-1));
            }
        }
    }
}



template <typename SP>
inline
void
FDMatrix::nine_point_structure(SP &structure) const
{
  for (unsigned int i=0; i<=ny-2; i++)
    {
      for (unsigned int j=0; j<=nx-2; j++)
        {
          // Number of the row to be entered
          unsigned int row = j+(nx-1)*i;
          structure.add(row, row);
          if (j>0)
            {
              structure.add(row-1, row);
              structure.add(row, row-1);
              if (i>0)
                {
                  structure.add(row-1, row-(nx-1));
                  structure.add(row-(nx-1), row-1);
                }
              if (i<ny-2)
                {
                  structure.add(row-1, row+(nx-1));
                  structure.add(row+(nx-1), row-1);
                }
            }
          if (j<nx-2)
            {
              structure.add(row+1, row);
              structure.add(row, row+1);
              if (i>0)
                {
                  structure.add(row+1, row-(nx-1));
                  structure.add(row-(nx-1), row+1);
                }
              if (i<ny-2)
                {
                  structure.add(row+1, row+(nx-1));
                  structure.add(row+(nx-1), row+1);
                }
            }
          if (i>0)
            {
              structure.add(row-(nx-1), row);
              structure.add(row, row-(nx-1));
            }
          if (i<ny-2)
            {
              structure.add(row+(nx-1), row);
              structure.add(row, row+(nx-1));
            }
        }
    }
}



template<typename MATRIX>
void
FDMatrix::nine_point(MATRIX &A, bool) const
{
  for (unsigned int i=0; i<=ny-2; i++)
    {
      for (unsigned int j=0; j<=nx-2; j++)
        {
          // Number of the row to be entered
          unsigned int row = j+(nx-1)*i;

          A.set(row, row, 20.);
          if (j>0)
            {
              A.set(row-1, row, -4.);
              A.set(row, row-1, -4.);
              if (i>0)
                {
                  A.set(row-1, row-(nx-1), -1.);
                  A.set(row-(nx-1), row-1, -1.);
                }
              if (i<ny-2)
                {
                  A.set(row-1, row+(nx-1), -1.);
                  A.set(row+(nx-1), row-1, -1.);
                }
            }
          if (j<nx-2)
            {
              A.set(row+1, row, -4.);
              A.set(row, row+1, -4.);
              if (i>0)
                {
                  A.set(row+1, row-(nx-1), -1.);
                  A.set(row-(nx-1), row+1, -1.);
                }
              if (i<ny-2)
                {
                  A.set(row+1, row+(nx-1), -1.);
                  A.set(row+(nx-1), row+1, -1.);
                }
            }
          if (i>0)
            {
              A.set(row-(nx-1), row, -4.);
              A.set(row, row-(nx-1), -4.);
            }
          if (i<ny-2)
            {
              A.set(row+(nx-1), row, -4.);
              A.set(row, row+(nx-1), -4.);
            }
        }
    }
}

template<typename MATRIX>
inline
void
FDMatrix::five_point(MATRIX &A, bool nonsymmetric) const
{
  for (unsigned int i=0; i<=ny-2; i++)
    {
      for (unsigned int j=0; j<=nx-2; j++)
        {
          // Number of the row to be entered
          unsigned int row = j+(nx-1)*i;
          if (nonsymmetric)
            A.set(row, row, 5.);
          else
            A.set(row, row, 4.);
          if (j>0)
            {
              if (nonsymmetric)
                A.set(row-1, row, -2.);
              else
                A.set(row-1, row, -1.);
              A.set(row, row-1, -1.);
            }
          if (j<nx-2)
            {
              A.set(row+1, row, -1.);
              A.set(row, row+1, -1.);
            }
          if (i>0)
            {
              A.set(row-(nx-1), row, -1.);
              A.set(row, row-(nx-1), -1.);
            }
          if (i<ny-2)
            {
              A.set(row+(nx-1), row, -1.);
              A.set(row, row+(nx-1), -1.);
            }
        }
    }
}



template<typename MATRIX>
inline
void
FDMatrix::upwind(MATRIX &A, bool back) const
{
  for (unsigned int i=0; i<=ny-2; i++)
    {
      for (unsigned int j=0; j<=nx-2; j++)
        {
          // Number of the row to be entered
          unsigned int row = j+(nx-1)*i;
          A.set(row, row, 3.);

          if (j>0 && !back)
            A.set(row, row-1, -1.);

          if (j<nx-2 && back)
            A.set(row, row+1, -1.);
        }
    }
}

template<typename number>
inline
void
FDMatrix::gnuplot_print(std::ostream &s, const Vector<number> &V) const
{
  for (unsigned int i=0; i<=ny-2; i++)
    {
      for (unsigned int j=0; j<=nx-2; j++)
        {
          // Number of the row to be entered
          unsigned int row = j+(nx-1)*i;
          s << (j+1)/(float)nx << '\t' << (i+1)/(float)ny << '\t' << V(row) << std::endl;
        }
      s << std::endl;
    }
  s << std::endl;
}
