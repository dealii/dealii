// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


// Test that all of the low-level math operations function as expected
// and that their derivatives can be computed using the low-level
// SymEngineWrapper operations.

#include <deal.II/base/logstream.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/differentiation/sd.h>

#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iomanip>

#include "../tests.h"

using namespace dealii;
namespace SD = Differentiation::SD;

template <typename NumberType>
struct FunctionAdd
{
  static std::string
  name()
  {
    return "Add";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return (s1 + s2) * (s1 + s2);
  };
};

template <typename NumberType>
struct FunctionSub
{
  static std::string
  name()
  {
    return "Subtract";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return (s1 - s2) * (s1 - s2);
  };
};

template <typename NumberType>
struct FunctionMul
{
  static std::string
  name()
  {
    return "Multiply";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return (s1 * s2) * (s1 * s2);
  };
};

template <typename NumberType>
struct FunctionDiv
{
  static std::string
  name()
  {
    return "Divide";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return (s1 / s2) * (s1 / s2);
  };
};

template <typename NumberType>
struct FunctionPow
{
  static std::string
  name()
  {
    return "Power";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return NumberType(2.0) * std::pow(s1, s2);
  };
};

template <typename NumberType>
struct FunctionSqrt
{
  static std::string
  name()
  {
    return "Sqrt";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return NumberType(2.0) * std::sqrt(s1 * s2);
  };
};

template <typename NumberType>
struct FunctionExp
{
  static std::string
  name()
  {
    return "Exponential";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return s1 * s1 * std::exp(s2);
  };
};

template <typename NumberType>
struct FunctionLog
{
  static std::string
  name()
  {
    return "Logarithm";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return s1 * s1 * std::log(s2);
  };
};

template <typename NumberType>
struct FunctionLogBase
{
  static std::string
  name()
  {
    return "Logarithm with Base";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return std::log(s2) / std::log(s1);
  };
};

template <>
struct FunctionLogBase<SD::Expression>
{
  static std::string
  name()
  {
    return "Logarithm with Base";
  }

  static SD::Expression
  psi(const SD::Expression &s1, const SD::Expression &s2)
  {
    return SD::log(s2, s1);
  };
};

template <typename NumberType>
struct FunctionLog10
{
  static std::string
  name()
  {
    return "Logarithm base 10";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return s1 * s1 * std::log10(s2);
  };
};

template <typename NumberType>
struct FunctionSin
{
  static std::string
  name()
  {
    return "Sin";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return std::sin(s2 / s1);
  };
};

template <typename NumberType>
struct FunctionCos
{
  static std::string
  name()
  {
    return "Cos";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return std::cos(s2 / s1);
  };
};

template <typename NumberType>
struct FunctionTan
{
  static std::string
  name()
  {
    return "Tan";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return std::tan(s2 / s1);
  };
};

template <typename NumberType>
struct FunctionASin
{
  static std::string
  name()
  {
    return "ASin";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return std::asin(NumberType(1.0) / (s1 * s2));
  };
};

template <typename NumberType>
struct FunctionACos
{
  static std::string
  name()
  {
    return "ACos";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return std::acos(NumberType(1.0) / (s1 * s2));
  };
};

template <typename NumberType>
struct FunctionATan
{
  static std::string
  name()
  {
    return "ATan";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return std::atan(NumberType(1.0) / (s1 * s2));
  };
};

template <typename NumberType>
struct FunctionATan2
{
  static std::string
  name()
  {
    return "ATan2";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return NumberType(3.0) * std::atan2(s2, s1);
  };
};

// TODO: cosec; sec; cot; acosec; asec; acot

template <typename NumberType>
struct FunctionSinh
{
  static std::string
  name()
  {
    return "Sinh";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return std::sinh(s2 / s1);
  };
};

template <typename NumberType>
struct FunctionCosh
{
  static std::string
  name()
  {
    return "Cosh";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return std::cosh(s2 / s1);
  };
};

template <typename NumberType>
struct FunctionTanh
{
  static std::string
  name()
  {
    return "Tanh";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return std::tanh(s2 / s1);
  };
};

template <typename NumberType>
struct FunctionErf
{
  static std::string
  name()
  {
    return "Erf";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return s1 * s1 * std::erf(s2);
  };
};

template <typename NumberType>
struct FunctionErfc
{
  static std::string
  name()
  {
    return "Erfc";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return s1 * s1 * std::erfc(s2);
  };
};

template <typename NumberType>
struct FunctionASinh
{
  static std::string
  name()
  {
    return "ASinh";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return std::asinh(NumberType(1.0) / (s1 * s2));
  };
};

template <typename NumberType>
struct FunctionACosh
{
  static std::string
  name()
  {
    return "ACosh";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return std::acosh(s1 * s2);
  };
};

template <typename NumberType>
struct FunctionATanh
{
  static std::string
  name()
  {
    return "ATanh";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return std::atanh(NumberType(1.0) / (s1 * s2));
  };
};

template <typename NumberType>
struct FunctionAbs
{
  static std::string
  name()
  {
    return "Abs";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return std::abs(s1 * s2) * std::abs(-1.0 * s1 / s2);
  };
};

template <typename NumberType>
struct FunctionMax
{
  static std::string
  name()
  {
    return "Max";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return 2.5 * std::max(s1, s2) * std::max(s1, s2);
  };
};

template <typename NumberType>
struct FunctionMin
{
  static std::string
  name()
  {
    return "Min";
  }

  static NumberType
  psi(const NumberType &s1, const NumberType &s2)
  {
    return 1.5 * std::min(s1, s2) * std::min(s1, s2);
  };
};

template <typename number_t, template <typename> class FunctionStruct>
void
test_functions()
{
  deallog.push(FunctionStruct<number_t>::name());

  // Evaluation point
  const double       s1                      = 3.0;
  const double       s2                      = 2.0;
  const unsigned int n_independent_variables = 2;

  // Pure SD differentiation
  double             psi_sd_sd = 0.0;
  Vector<double>     Dpsi_sd_sd(n_independent_variables);
  FullMatrix<double> D2psi_sd_sd(n_independent_variables,
                                 n_independent_variables);
  {
    typedef SD::Expression SDNumberType;

    // Function and its derivatives
    typedef FunctionStruct<SDNumberType> func_sd;

    deallog.push("Symbolic computation: Function initialisation");
    const SDNumberType symb_s0("s1");
    const SDNumberType symb_s1("s2");
    const SDNumberType symb_psi = func_sd::psi(symb_s0, symb_s1);
    deallog << "symb_s0: " << symb_s0 << std::endl;
    deallog << "symb_s1: " << symb_s1 << std::endl;
    deallog << "symb_psi: " << symb_psi << std::endl;
    deallog.pop();

    deallog.push("Symbolic computation: Differentiation");
    const SDNumberType symb_dpsi_ds0 = symb_psi.differentiate(symb_s0);
    const SDNumberType symb_dpsi_ds1 = symb_psi.differentiate(symb_s1);
    const SDNumberType symb_d2psi_ds0_ds0 =
      symb_dpsi_ds0.differentiate(symb_s0);
    const SDNumberType symb_d2psi_ds1_ds0 =
      symb_dpsi_ds0.differentiate(symb_s1);
    const SDNumberType symb_d2psi_ds0_ds1 =
      symb_dpsi_ds1.differentiate(symb_s0);
    const SDNumberType symb_d2psi_ds1_ds1 =
      symb_dpsi_ds1.differentiate(symb_s1);
    deallog << "symb_dpsi_ds0: " << symb_dpsi_ds0 << std::endl;
    deallog << "symb_dpsi_ds1: " << symb_dpsi_ds1 << std::endl;
    deallog << "symb_d2psi_ds0_ds0: " << symb_d2psi_ds0_ds0 << std::endl;
    deallog << "symb_d2psi_ds1_ds0: " << symb_d2psi_ds1_ds0 << std::endl;
    deallog << "symb_d2psi_ds0_ds1: " << symb_d2psi_ds0_ds1 << std::endl;
    deallog << "symb_d2psi_ds1_ds1: " << symb_d2psi_ds1_ds1 << std::endl;
    deallog.pop();

    deallog.push("Symbolic computation: Substitution");
    SD::types::substitution_map sub_vals;
    SD::add_to_substitution_map(sub_vals, symb_s0, s1);
    SD::add_to_substitution_map(sub_vals, symb_s1, s2);

    const SDNumberType subs_psi      = symb_psi.substitute(sub_vals);
    const SDNumberType subs_dpsi_ds0 = symb_dpsi_ds0.substitute(sub_vals);
    const SDNumberType subs_dpsi_ds1 = symb_dpsi_ds1.substitute(sub_vals);
    const SDNumberType subs_d2psi_ds0_ds0 =
      symb_d2psi_ds0_ds0.substitute(sub_vals);
    const SDNumberType subs_d2psi_ds1_ds0 =
      symb_d2psi_ds1_ds0.substitute(sub_vals);
    const SDNumberType subs_d2psi_ds0_ds1 =
      symb_d2psi_ds0_ds1.substitute(sub_vals);
    const SDNumberType subs_d2psi_ds1_ds1 =
      symb_d2psi_ds1_ds1.substitute(sub_vals);
    deallog << "subs_psi: " << subs_psi << std::endl;
    deallog << "subs_dpsi_ds0: " << subs_dpsi_ds0 << std::endl;
    deallog << "subs_dpsi_ds1: " << subs_dpsi_ds1 << std::endl;
    deallog << "subs_d2psi_ds0_ds0: " << subs_d2psi_ds0_ds0 << std::endl;
    deallog << "subs_d2psi_ds1_ds0: " << subs_d2psi_ds1_ds0 << std::endl;
    deallog << "subs_d2psi_ds0_ds1: " << subs_d2psi_ds0_ds1 << std::endl;
    deallog << "subs_d2psi_ds1_ds1: " << subs_d2psi_ds1_ds1 << std::endl;

    psi_sd_sd         = static_cast<double>(subs_psi);
    Dpsi_sd_sd[0]     = static_cast<double>(subs_dpsi_ds0);
    Dpsi_sd_sd[1]     = static_cast<double>(subs_dpsi_ds1);
    D2psi_sd_sd[0][0] = static_cast<double>(subs_d2psi_ds0_ds0);
    D2psi_sd_sd[0][1] = static_cast<double>(subs_d2psi_ds1_ds0);
    D2psi_sd_sd[1][0] = static_cast<double>(subs_d2psi_ds0_ds1);
    D2psi_sd_sd[1][1] = static_cast<double>(subs_d2psi_ds1_ds1);

    deallog << "psi_sd_sd: " << psi_sd_sd << std::endl;
    deallog << "Dpsi_sd_sd[0]: " << Dpsi_sd_sd[0] << std::endl;
    deallog << "Dpsi_sd_sd[1]: " << Dpsi_sd_sd[1] << std::endl;
    deallog << "D2psi_sd_sd[0][0]: " << D2psi_sd_sd[0][0] << std::endl;
    deallog << "D2psi_sd_sd[0][1]: " << D2psi_sd_sd[0][1] << std::endl;
    deallog << "D2psi_sd_sd[1][0]: " << D2psi_sd_sd[1][0] << std::endl;
    deallog << "D2psi_sd_sd[1][1]: " << D2psi_sd_sd[1][1] << std::endl;
    deallog.pop();
  }

  deallog.pop();
}

template <typename number_t>
void
test_all_functions()
{
  test_functions<number_t, FunctionAdd>();
  test_functions<number_t, FunctionSub>();
  test_functions<number_t, FunctionMul>();
  test_functions<number_t, FunctionDiv>();
  test_functions<number_t, FunctionPow>();
  test_functions<number_t, FunctionSqrt>();
  test_functions<number_t, FunctionExp>();
  test_functions<number_t, FunctionLog>();
  test_functions<number_t, FunctionLogBase>();
  test_functions<number_t, FunctionLog10>();
  test_functions<number_t, FunctionSin>();
  test_functions<number_t, FunctionCos>();
  test_functions<number_t, FunctionTan>();
  test_functions<number_t, FunctionASin>();
  test_functions<number_t, FunctionACos>();
  test_functions<number_t, FunctionATan>();
  test_functions<number_t, FunctionATan2>();
  test_functions<number_t, FunctionSinh>();
  test_functions<number_t, FunctionCosh>();
  test_functions<number_t, FunctionTanh>();
  //  test_functions<number_t,FunctionAbs>(); // TODO: Error message: "Cannot
  //  perform substitution on derivative"
  //  test_functions<number_t,FunctionMax>(); // TODO: Error message: "Cannot
  //  perform substitution on derivative"
  //  test_functions<number_t,FunctionMin>(); // TODO: Error message: "Cannot
  //  perform substitution on derivative"

  test_functions<number_t, FunctionASinh>();
  test_functions<number_t, FunctionACosh>();
  test_functions<number_t, FunctionATanh>();
  test_functions<number_t, FunctionErf>();
  test_functions<number_t, FunctionErfc>();
}

int
main()
{
  initlog();

  deallog.push("Double");
  {
    test_all_functions<double>();
    deallog << "OK" << std::endl;
  }
  deallog.pop();

  deallog << "OK" << std::endl;
}
