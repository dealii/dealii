
#include "../tests.h"
#include <fstream>
#include <iomanip>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_fire.h>

using namespace dealii;


// Test to verify correctness of SolverFIRE::solve()
// The objective function is f(x,y) = x^2 + y^2.


using vector_t = typename dealii::Vector<double>;


double compute (vector_t &G, const vector_t &X)
{
  AssertThrow (X.size() == 2 && G.size() == 2,
               ExcInternalError());

  G(0) = 2*X(0);
  G(1) = 2*X(1);

  return X.norm_sqr();
}



void check_value (const double x,
                  const double y,
                  const double tol)
{
  vector_t X;

  X.reinit(2, true);

  // Use this to initialize DiagonalMatrix
  X = 1.;

  // Create inverse diagonal matrix.
  DiagonalMatrix<vector_t> inv_mass;
  inv_mass.reinit(X);

  // Set initial iterate.
  X(0) = x;
  X(1) = y;

  auto additional_data =
    SolverFIRE<vector_t>::AdditionalData(0.1, 1., 1);

  SolverControl solver_control (200, tol);

  SolverFIRE<vector_t> fire (solver_control, additional_data);

  fire.solve(compute, X, inv_mass);

  deallog << "FIRE::Solution vector: ";

  X.print(deallog);
}

int main ()
{
  std::ofstream logfile("output");
//  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  check_value(  10,  -2, 1e-15);
  check_value(-0.1, 0.1, 1e-15);
  check_value( 9.1,-6.1, 1e-15);

}
