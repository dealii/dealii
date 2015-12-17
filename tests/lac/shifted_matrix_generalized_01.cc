#include <deal.II/base/logstream.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/shifted_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>

int main()
{
  using namespace dealii;

  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(4);
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  FullMatrix<double> m1(2, 2);
  FullMatrix<double> m2(2, 2);

  m1(0, 0) = 2.0;
  m1(0, 1) = 1.0;
  m1(1, 0) = 1.0;
  m1(1, 1) = 2.0;

  m2(0, 0) = 1.0;
  m2(0, 1) = 1.0;
  m2(1, 0) = 1.0;
  m2(1, 1) = 1.0;

  ShiftedMatrixGeneralized<FullMatrix<double>, FullMatrix<double>, Vector<double> >
  shifted_matrix(m1, m2, 2.0);

  Vector<double> src(2);
  src[0] = 3.0;
  src[1] = 7.0;
  Vector<double> dst(2);
  shifted_matrix.vmult(dst, src);

  deallog << dst[0] << ", ";
  deallog << dst[1] << std::endl;
}
