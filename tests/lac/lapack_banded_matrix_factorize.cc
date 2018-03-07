#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_banded_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <numeric>

#include "../tests.h"

// Verify that LAPACKBandedMatrix::compute_lu_factorization and friends are
// implemented correctly by comparing their results to that of FullMatrix and
// SparseMatrix.

int
main()
{
  initlog();

  // test 1: badly conditioned tridiagonal matrix. The backward error should
  // be much smaller than the forward error.
  {
    const std::size_t          n_rows = 100;
    LAPACKBandedMatrix<double> banded(n_rows, n_rows, 1, 1);
    deallog << "banded shape: " << banded.m() << ", " << banded.n() << std::endl
            << "banded subdiagonals: " << banded.n_stored_subdiagonals()
            << std::endl
            << "banded superdiagonals: " << banded.n_stored_superdiagonals()
            << std::endl;

    // This is a well-known badly conditioned tridiagonal matrix. Citation:
    // F.W. Dorr, An example of ill-conditioning in the numerical solution of
    // singular perturbation problems, Math. Comp., 25 (1971), pp. 271-283.
    const double epsilon = 0.01;
    const double h       = 1.0 / n_rows;

    // while we are here, test the set function as well
    for (std::size_t row_n = 0; row_n < n_rows / 2; ++row_n)
      {
        const double position    = double(row_n + 1) / double(n_rows);
        banded(row_n, row_n + 1) = -epsilon / (h * h) - (0.5 - position) / h;

        if (row_n == 0)
          {
            banded.set(
              row_n, row_n, -banded(row_n, row_n + 1) + epsilon / (h * h));
          }
        else
          {
            banded.set(row_n, row_n - 1, -epsilon / (h * h));
            banded.set(row_n,
                       row_n,
                       -(banded(row_n, row_n - 1) + banded(row_n, row_n + 1)));
          }
      }
    for (std::size_t row_n = n_rows / 2; row_n < n_rows; ++row_n)
      {
        const double position    = double(row_n + 1) / double(n_rows);
        banded(row_n, row_n - 1) = -epsilon / (h * h) + (0.5 - position) / h;

        if (row_n == n_rows - 1)
          {
            banded(row_n, row_n) =
              -banded(row_n, row_n - 1) + epsilon / (h * h);
          }
        else
          {
            banded(row_n, row_n + 1) = -epsilon / (h * h);
            banded(row_n, row_n) =
              -(banded(row_n, row_n - 1) + banded(row_n, row_n + 1));
          }
      }
    const LAPACKBandedMatrix<double> copy(banded);
    banded.compute_lu_factorization();
    Vector<double> rhs(n_rows);
    std::fill(rhs.begin(), rhs.end(), 1.0);
    Vector<double> solution = rhs;

    const std::pair<double, double> errors =
      banded.solve(solution, false, true);

    Vector<double> other_rhs(n_rows);
    copy.vmult(other_rhs, solution);
    other_rhs -= rhs;
    deallog << "difference less than 1e-11: " << (other_rhs.l2_norm() < 1e-14)
            << std::endl;
    deallog << "backward error less than 1e-14: " << (errors.first < 1e-14)
            << std::endl;
    deallog << "forward error less than 1e-11: " << (errors.second < 1e-11)
            << std::endl;
  }

  // test 2: something a bit bigger: check with a 2D Laplacian
  {
    Triangulation<2> triangulation;
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(2);
    FE_Q<2>       fe(1);
    DoFHandler<2> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe);

    SparsityPattern sparsity_pattern(dof_handler.n_dofs(),
                                     dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
    sparsity_pattern.compress();
    SparseMatrix<double> mass_matrix(sparsity_pattern);

    MatrixTools::create_mass_matrix(
      dof_handler, QGauss<2>(fe.degree + 2), mass_matrix);
    Vector<double> rhs(dof_handler.n_dofs());
    VectorTools::create_right_hand_side(
      dof_handler, QGauss<2>(fe.degree + 2), ConstantFunction<2>(42.0), rhs);

    LAPACKBandedMatrix<double> banded(mass_matrix);
    // TODO 'const LAPACKBandedMatrix<double> copy = banded;' doesn't work
    const LAPACKBandedMatrix<double> copy(banded);
    deallog << "banded shape: " << banded.m() << ", " << banded.n() << std::endl
            << "banded subdiagonals: " << banded.n_stored_subdiagonals()
            << std::endl
            << "banded superdiagonals: " << banded.n_stored_superdiagonals()
            << std::endl;

    banded.compute_lu_factorization();
    Vector<double> solution = rhs;
    banded.solve(solution, false, true);
    Vector<double> other_rhs(banded.m());
    copy.vmult(other_rhs, solution);
    other_rhs -= rhs;
    other_rhs.print(deallog.get_file_stream(), 2, false);
    deallog << "difference less than 1e-14: " << (other_rhs.l2_norm() < 1.0e-14)
            << std::endl;
  }
}
