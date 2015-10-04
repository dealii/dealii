// Test the sizes of fundamental types, their pointers, and vectors comprised
// of either.
#include <deal.II/base/config.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/memory_consumption.h>

#include <fstream>
#include <iostream>
#include <vector>

#define CONCAT(a, b) a##b
#define MAKE_LABEL(a, b) CONCAT(a, b)
#define EQUALITY_TEST(a) EqualityWithSizeofTest<a> MAKE_LABEL(b_, __LINE__);   \
  EqualityWithSizeofTest<a *> MAKE_LABEL(c_, __LINE__);                        \
  EqualityWithSizeofTest<a **> MAKE_LABEL(d_, __LINE__);

using namespace dealii;

template<typename T>
struct EqualityWithSizeofTest
{
  EqualityWithSizeofTest()
  {
    T t = 0;
    deallog << (sizeof(t) == dealii::MemoryConsumption::memory_consumption(t))
            << std::endl;

    std::vector<T> vector_test(42);
    deallog << ((sizeof(std::vector<T>) + vector_test.capacity()*sizeof(T))
                == dealii::MemoryConsumption::memory_consumption(vector_test))
            << std::endl;
  }
};

int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  // do not test value type (just pointers) with void
  EqualityWithSizeofTest<void *> a;
  EqualityWithSizeofTest<void **> b;

  // std::vector<bool> is overloaded to work as a bit field, so the vector
  // version should be unequal
  EQUALITY_TEST(bool)
  EQUALITY_TEST(signed char)
  EQUALITY_TEST(unsigned char)
  // char* is overloaded to measure the length of a C string, so the pointer
  // version should be unequal
  EQUALITY_TEST(char)
  EQUALITY_TEST(wchar_t)

  EQUALITY_TEST(short)
  EQUALITY_TEST(short int)
  EQUALITY_TEST(signed short)
  EQUALITY_TEST(signed short int)
  EQUALITY_TEST(unsigned short)
  EQUALITY_TEST(unsigned short int)
  EQUALITY_TEST(int)
  EQUALITY_TEST(signed)
  EQUALITY_TEST(signed int)
  EQUALITY_TEST(unsigned)
  EQUALITY_TEST(unsigned int)
  EQUALITY_TEST(long)
  EQUALITY_TEST(long int)
  EQUALITY_TEST(signed long)
  EQUALITY_TEST(signed long int)
  EQUALITY_TEST(unsigned long int)

  EQUALITY_TEST(float)
  EQUALITY_TEST(double)
  EQUALITY_TEST(long double)
}
