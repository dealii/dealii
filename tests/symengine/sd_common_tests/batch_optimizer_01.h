// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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


// Header file:
// Test that all of the low-level math operations function as expected and
// that their values and derivatives can be computed using the optimizer.

#include <deal.II/differentiation/ad.h>
#include <deal.II/differentiation/sd.h>

#include <iostream>

using namespace dealii;
namespace AD = Differentiation::AD;
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

template <typename number_t,
          enum SD::OptimizerType     opt_method,
          enum SD::OptimizationFlags opt_flags,
          template <typename> class FunctionStruct>
void
test_functions()
{
  deallog.push(FunctionStruct<number_t>::name());

  // Evaluation point
  const number_t     s1                      = 3.0;
  const number_t     s2                      = 2.0;
  const unsigned int n_independent_variables = 2;

  // Pure SD differentiation
  number_t             psi_sd_sd = 0.0;
  Vector<number_t>     Dpsi_sd_sd(n_independent_variables);
  FullMatrix<number_t> D2psi_sd_sd(n_independent_variables,
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

    deallog.push("Symbolic computation: Optimization");
    SD::BatchOptimizer<number_t> optimizer;
    optimizer.set_optimization_method(opt_method, opt_flags);
    SD::types::substitution_map sub_vals_optim;
    SD::add_to_symbol_map(sub_vals_optim, symb_s0, symb_s1);
    optimizer.register_symbols(sub_vals_optim); // Independent symbols
    optimizer.register_functions(
      symb_psi,
      symb_dpsi_ds0,
      symb_dpsi_ds1,
      symb_d2psi_ds0_ds0,
      symb_d2psi_ds1_ds0,
      symb_d2psi_ds0_ds1,
      symb_d2psi_ds1_ds1); // Dependent symbolic functions
    std::cout << FunctionStruct<number_t>::name() << ": About to optimise... "
              << std::flush;
    optimizer.optimize();
    std::cout << "done." << std::endl;
    deallog.pop();

    deallog.push("Symbolic computation: Substitution");
    SD::types::substitution_map sub_vals;
    SD::add_to_substitution_map(sub_vals, symb_s0, s1);
    SD::add_to_substitution_map(sub_vals, symb_s1, s2);

    // Perform substitution of symbols
    std::cout << FunctionStruct<number_t>::name() << ": About to substitute... "
              << std::flush;
    optimizer.substitute(sub_vals);
    std::cout << "done." << std::endl;
    // optimizer.print(std::cout);

    // Evaluate
    deallog << "subs_psi: " << optimizer.evaluate(symb_psi) << std::endl;
    deallog << "subs_dpsi_ds0: " << optimizer.evaluate(symb_dpsi_ds0)
            << std::endl;
    deallog << "subs_dpsi_ds1: " << optimizer.evaluate(symb_dpsi_ds1)
            << std::endl;
    deallog << "subs_d2psi_ds0_ds0: " << optimizer.evaluate(symb_d2psi_ds0_ds0)
            << std::endl;
    deallog << "subs_d2psi_ds1_ds0: " << optimizer.evaluate(symb_d2psi_ds1_ds0)
            << std::endl;
    deallog << "subs_d2psi_ds0_ds1: " << optimizer.evaluate(symb_d2psi_ds0_ds1)
            << std::endl;
    deallog << "subs_d2psi_ds1_ds1: " << optimizer.evaluate(symb_d2psi_ds1_ds1)
            << std::endl;
    deallog.pop();
  }

  deallog.pop();
}

template <typename number_t,
          enum SD::OptimizerType     opt_method,
          enum SD::OptimizationFlags opt_flags>
void
test_all_functions()
{
  const bool skip_non_implemented_functions =
    (numbers::NumberTraits<number_t>::is_complex &&
     (opt_method == SD::OptimizerType::dictionary ||
      opt_method == SD::OptimizerType::lambda));

  std::cout << "Optimisation method: " << opt_method << "\n"
            << "Optimisation flags: " << opt_flags << "\n"
            << "Skip non-implemented functions: "
            << skip_non_implemented_functions << std::endl;

  test_functions<number_t, opt_method, opt_flags, FunctionAdd>();
  test_functions<number_t, opt_method, opt_flags, FunctionSub>();
  test_functions<number_t, opt_method, opt_flags, FunctionMul>();
  test_functions<number_t, opt_method, opt_flags, FunctionDiv>();
  test_functions<number_t, opt_method, opt_flags, FunctionPow>();
  test_functions<number_t, opt_method, opt_flags, FunctionSqrt>();
  test_functions<number_t, opt_method, opt_flags, FunctionExp>();
  test_functions<number_t, opt_method, opt_flags, FunctionLog>();
  test_functions<number_t, opt_method, opt_flags, FunctionLogBase>();
  test_functions<number_t, opt_method, opt_flags, FunctionLog10>();
  test_functions<number_t, opt_method, opt_flags, FunctionSin>();
  test_functions<number_t, opt_method, opt_flags, FunctionCos>();
  test_functions<number_t, opt_method, opt_flags, FunctionTan>();
  test_functions<number_t, opt_method, opt_flags, FunctionASin>();
  test_functions<number_t, opt_method, opt_flags, FunctionACos>();
  test_functions<number_t, opt_method, opt_flags, FunctionATan>();
  if (!skip_non_implemented_functions)
    test_functions<number_t, opt_method, opt_flags, FunctionATan2>();
  test_functions<number_t, opt_method, opt_flags, FunctionSinh>();
  test_functions<number_t, opt_method, opt_flags, FunctionCosh>();
  test_functions<number_t, opt_method, opt_flags, FunctionTanh>();
  //  test_functions<number_t,opt_method,opt_flags,FunctionAbs>(); // Cannot
  //  perform substitution on derivative
  //  test_functions<number_t,opt_method,opt_flags,FunctionMax>(); // Cannot
  //  perform substitution on derivative
  //  test_functions<number_t,opt_method,opt_flags,FunctionMin>(); // Cannot
  //  perform substitution on derivative

  test_functions<number_t, opt_method, opt_flags, FunctionASinh>();
  test_functions<number_t, opt_method, opt_flags, FunctionACosh>();
  test_functions<number_t, opt_method, opt_flags, FunctionATanh>();

  if (!skip_non_implemented_functions)
    test_functions<number_t, opt_method, opt_flags, FunctionErf>();
  if (!skip_non_implemented_functions)
    test_functions<number_t, opt_method, opt_flags, FunctionErfc>();

  deallog << "OK" << std::endl;
}

template <enum SD::OptimizerType     opt_method,
          enum SD::OptimizationFlags opt_flags>
void
run_tests(const int n_runs = 1)
{
  deallog.push("Float");
  test_all_functions<float, opt_method, opt_flags>();
  deallog.pop();

  deallog.push("Double");
  test_all_functions<double, opt_method, opt_flags>();
  deallog.pop();

  // The LLVM optimizer does not currently support complex numbers.
  if (opt_method != SD::OptimizerType::llvm)
    {
      deallog.push("Complex float");
      test_all_functions<std::complex<float>, opt_method, opt_flags>();
      deallog.pop();

      deallog.push("Complex double");
      test_all_functions<std::complex<double>, opt_method, opt_flags>();
      deallog.pop();
    }
}
