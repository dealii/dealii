// tests column padding for different sizes

#include <RFAStDFT/block_csr_matrix.h>

#include <fstream>
#include <iostream>

using namespace dealii;
using namespace RealFAStDFT;


template <typename NumberType = double>
void test(const std::vector<unsigned int> sizes)
{
  const NumberType dummy = NumberType(0.);
  const unsigned int CL = 64;
  std::cout << "Testing  " << dealii::Utilities::type_to_string(dummy)
            << std::endl
            << "  size:  " << sizeof(NumberType) << std::endl
            << "  SIMD:  " << VectorizedArray<NumberType>::n_array_elements
            << std::endl
            << "  CL:    " << CL << std::endl
            << "  el/CL: " << CL / sizeof(NumberType) << std::endl
            << std::endl;

  for (const auto s : sizes)
    {
      const auto upper = RealFAStDFT::internal::padded_size<NumberType>(s);
      Assert(RealFAStDFT::internal::ceil_divisible_by(
               s, VectorizedArray<NumberType>::n_array_elements) <= upper,
             ExcInternalError());
      Assert(RealFAStDFT::internal::ceil_divisible_by(
               s, CL / sizeof(NumberType)) == upper,
             ExcInternalError());
      std::cout << s << " -> " << upper << std::endl;
    }
}

int main()
{

  test<double>({3,4,7,8,12,15,16,24,32});
  test<float>({3,4,7,8,12,15,16,24,32});

  return 0;
}
