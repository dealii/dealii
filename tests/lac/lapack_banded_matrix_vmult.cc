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

#include <numeric>

#include "../tests.h"

// Verify that LAPACKBandedMatrix::vmult and friends are implemented correctly
// by comparing their results to that of FullMatrix and SparseMatrix.

int
main()
{
  initlog();

  // test 1: tridiagonal matrix copied from a FullMatrix
  {
    const unsigned int n_rows = 8;
    FullMatrix<double> full_matrix(n_rows, n_rows);
    for (unsigned int row_n = 0; row_n < n_rows; ++row_n)
      {
        full_matrix(row_n, row_n) = 2.0 + row_n / 2;
        if (0 < row_n)
          full_matrix(row_n, row_n - 1) = -1.0 - row_n;
        if (row_n < n_rows - 1)
          full_matrix(row_n, row_n + 1) = -1.0 - row_n / 2;
      }

    LAPACKBandedMatrix<double> banded(full_matrix);
    deallog << "banded shape: " << banded.m() << ", " << banded.n() << std::endl
            << "banded subdiagonals: " << banded.n_stored_subdiagonals()
            << std::endl
            << "banded superdiagonals: " << banded.n_stored_superdiagonals()
            << std::endl;

    {
      Vector<double> in(banded.m());
      Vector<double> out(banded.m());
      std::iota(in.begin(), in.end(), 2.0);
      deallog << "vmult results should be equal:" << std::endl;
      banded.vmult(out, in);
      out.print(deallog);
      full_matrix.vmult(out, in);
      out.print(deallog);
    }

    {
      Vector<double> in1(banded.m());
      Vector<double> in2(banded.m());
      Vector<double> out(banded.m());
      std::iota(in1.begin(), in1.end(), 2.0);
      std::iota(in2.begin(), in2.end(), -1.0);

      deallog << "vmult(x, y, true) results should be equal:" << std::endl;
      banded.vmult(out, in1, false);
      banded.vmult(out, in2, true);
      out.print(deallog);

      full_matrix.vmult(out, in1, false);
      full_matrix.vmult(out, in2, true);
      out.print(deallog);
    }

    {
      Vector<double> in1(banded.m());
      Vector<double> in2(banded.m());
      Vector<double> out(banded.m());
      std::iota(in1.begin(), in1.end(), 2.0);
      std::iota(in2.begin(), in2.end(), -1.0);

      deallog << "vmult_add results should be equal:" << std::endl;
      banded.vmult(out, in1);
      banded.vmult_add(out, in2);
      out.print(deallog);

      full_matrix.vmult(out, in1);
      full_matrix.vmult_add(out, in2);
      out.print(deallog);
    }

    {
      Vector<double> in(banded.m());
      Vector<double> out(banded.m());
      std::iota(in.begin(), in.end(), 2.0);
      deallog << "Tvmult results should be equal:" << std::endl;
      banded.Tvmult(out, in);
      out.print(deallog);
      full_matrix.Tvmult(out, in);
      out.print(deallog);
    }

    {
      Vector<double> in1(banded.m());
      Vector<double> in2(banded.m());
      Vector<double> out(banded.m());
      std::iota(in1.begin(), in1.end(), 2.0);
      std::iota(in2.begin(), in2.end(), -1.0);

      deallog << "Tvmult(x, y, true) results should be equal:" << std::endl;
      banded.Tvmult(out, in1, false);
      banded.Tvmult(out, in2, true);
      out.print(deallog);

      full_matrix.Tvmult(out, in1, false);
      full_matrix.Tvmult(out, in2, true);
      out.print(deallog);
    }

    {
      Vector<double> in1(banded.m());
      Vector<double> in2(banded.m());
      Vector<double> out(banded.m());
      std::iota(in1.begin(), in1.end(), 2.0);
      std::iota(in2.begin(), in2.end(), -1.0);

      deallog << "Tvmult_add results should be equal:" << std::endl;
      banded.Tvmult(out, in1);
      banded.Tvmult_add(out, in2);
      out.print(deallog);

      full_matrix.Tvmult(out, in1);
      full_matrix.Tvmult_add(out, in2);
      out.print(deallog);
    }
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
    SparseMatrix<double> laplace_matrix(sparsity_pattern);

    MatrixTools::create_laplace_matrix(
      dof_handler, QGauss<2>(fe.degree + 2), laplace_matrix);
    LAPACKBandedMatrix<double> banded(laplace_matrix);
    deallog << "banded shape: " << banded.m() << ", " << banded.n() << std::endl
            << "banded subdiagonals: " << banded.n_stored_subdiagonals()
            << std::endl
            << "banded superdiagonals: " << banded.n_stored_superdiagonals()
            << std::endl;

    Vector<double> in(banded.m());
    std::iota(in.begin(), in.end(), -double(in.size() / 2));
    Vector<double> out(banded.m());
    Vector<double> out2(banded.m());
    banded.vmult(out, in);
    out.print(deallog.get_file_stream(), 0, false);
    laplace_matrix.vmult(out2, in);
    out.print(deallog.get_file_stream(), 0, false);
    out -= out2;
    AssertThrow(out.l2_norm() < 1e-14, ExcInternalError());
  }
}
