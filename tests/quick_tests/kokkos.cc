
#include <deal.II/base/config.h>

#include <Kokkos_Core.hpp>

struct FillFunctor
{
  KOKKOS_FUNCTION void
  operator()(const int i) const
  {
    m_values(i) = i;
  }
  Kokkos::View<int *> m_values;
};

struct SumFunctor
{
  KOKKOS_FUNCTION void
  operator()(const int i, int &sum) const
  {
    sum += m_values(i);
  }
  Kokkos::View<int *> m_values;
};

int
main(int argc, char *argv[])
{
  const long n = 100;
  int        sum;
  Kokkos::initialize(argc, argv);
  {
    Kokkos::View<int *> values("values", n);
    Kokkos::parallel_for(n, FillFunctor{values});
    Kokkos::parallel_reduce(n, SumFunctor{values}, sum);
  }
  Kokkos::finalize();
  return (sum != (n * (n - 1)) / 2);
}
