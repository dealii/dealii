#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/vector_memory.templates.h>

#include <vector>

#include "../tests.h"

// Check that we can create VectorMemory objects with vector-like classes in
// the standard library.

template <typename VectorType>
void
test_std_vector_pointer()
{
  GrowingVectorMemory<VectorType> mem;

  std::vector<typename VectorMemory<VectorType>::Pointer> va;
  va.push_back(typename VectorMemory<VectorType>::Pointer(mem));
  va.emplace_back(mem);
}

int
main()
{
  initlog();
  PrimitiveVectorMemory<std::vector<double>> primitive_memory;
  test_std_vector_pointer<std::vector<double>>();

  deallog << "OK" << std::endl;
}
