// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2019 by the deal.II authors
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
// Evaluation of a coupled system (tensor + vector + scalar components)
// using a helper class.
// This test is based off of helper_scalar_coupled_3_components_01.h, and checks
// that everything still works in tapeless only mode (i.e. when the
// start_recording_operations and stop_recording_operations calls are
// removed).

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <iostream>

#include "../tests.h"

namespace AD = dealii::Differentiation::AD;

// Function and its derivatives
template <int dim, typename NumberType>
struct FunctionsTestTensorVectorScalarCoupled
{
  static NumberType
  det_t(const Tensor<2, dim, NumberType> &t)
  {
    return determinant(t);
  }

  static Tensor<2, dim, NumberType>
  ddet_t_dt(const Tensor<2, dim, NumberType> &t)
  {
    return det_t(t) * transpose(invert(t));
  }

  static Tensor<4, dim, NumberType>
  d2det_t_dt_dt(const Tensor<2, dim, NumberType> &t)
  {
    const Tensor<2, dim, NumberType> t_inv = invert(t);
    Tensor<4, dim, NumberType>       dt_inv_trans_dt;
    // https://en.wikiversity.org/wiki/Introduction_to_Elasticity/Tensors#Derivative_of_the_determinant_of_a_tensor
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = 0; l < dim; ++l)
            dt_inv_trans_dt[i][j][k][l] = -t_inv[l][i] * t_inv[j][k];

    return det_t(t) * outer_product(transpose(t_inv), transpose(t_inv)) +
           det_t(t) * dt_inv_trans_dt;
  }

  static NumberType
  v_squ(const Tensor<1, dim, NumberType> &v)
  {
    return v * v;
  }

  static Tensor<1, dim, NumberType>
  dv_squ_dv(const Tensor<1, dim, NumberType> &v)
  {
    return 2.0 * v;
  }

  static Tensor<2, dim, NumberType>
  d2v_squ_dv_dv(const Tensor<1, dim, NumberType> &v)
  {
    static const Tensor<2, dim, NumberType> I(
      unit_symmetric_tensor<dim, NumberType>());
    return 2.0 * I;
  }

  // --------

  static const double sf;

  static NumberType
  psi(const Tensor<2, dim, NumberType> &t,
      const Tensor<1, dim, NumberType> &v,
      const NumberType &                s)
  {
    return std::pow(det_t(t), 2) * std::pow(v_squ(v), 3) * std::pow(s, sf);
  };

  static Tensor<2, dim, NumberType>
  dpsi_dt(const Tensor<2, dim, NumberType> &t,
          const Tensor<1, dim, NumberType> &v,
          const NumberType &                s)
  {
    return 2.0 * std::pow(det_t(t), 1) * ddet_t_dt(t) * std::pow(v_squ(v), 3) *
           std::pow(s, sf);
  };

  static Tensor<1, dim, NumberType>
  dpsi_dv(const Tensor<2, dim, NumberType> &t,
          const Tensor<1, dim, NumberType> &v,
          const NumberType &                s)
  {
    return std::pow(det_t(t), 2) * 3.0 * std::pow(v_squ(v), 2) * dv_squ_dv(v) *
           std::pow(s, sf);
  };

  static NumberType
  dpsi_ds(const Tensor<2, dim, NumberType> &t,
          const Tensor<1, dim, NumberType> &v,
          const NumberType &                s)
  {
    return std::pow(det_t(t), 2) * std::pow(v_squ(v), 3) * sf *
           std::pow(s, sf - 1.0);
  };

  static Tensor<4, dim, NumberType>
  d2psi_dt_dt(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType &                s)
  {
    return 2.0 * std::pow(v_squ(v), 3) *
           (pow(det_t(t), 0) * outer_product(ddet_t_dt(t), ddet_t_dt(t)) +
            std::pow(det_t(t), 1) * d2det_t_dt_dt(t)) *
           std::pow(s, sf);
  };

  static Tensor<3, dim, NumberType>
  d2psi_dv_dt(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType &                s)
  {
    return 2.0 * std::pow(det_t(t), 1) * 3.0 * std::pow(v_squ(v), 2) *
           outer_product(ddet_t_dt(t), dv_squ_dv(v)) * std::pow(s, sf);
  };

  static Tensor<2, dim, NumberType>
  d2psi_ds_dt(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType &                s)
  {
    return 2.0 * std::pow(det_t(t), 1) * ddet_t_dt(t) * std::pow(v_squ(v), 3) *
           sf * std::pow(s, sf - 1.0);
  };

  static Tensor<3, dim, NumberType>
  d2psi_dt_dv(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType &                s)
  {
    return 2.0 * std::pow(det_t(t), 1) * 3.0 * std::pow(v_squ(v), 2) *
           outer_product(dv_squ_dv(v), ddet_t_dt(t)) * std::pow(s, sf);
  };

  static Tensor<2, dim, NumberType>
  d2psi_dv_dv(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType &                s)
  {
    return std::pow(det_t(t), 2) * 3.0 *
           (2.0 * std::pow(v_squ(v), 1) *
              outer_product(dv_squ_dv(v), dv_squ_dv(v)) +
            std::pow(v_squ(v), 2) * d2v_squ_dv_dv(v)) *
           std::pow(s, sf);
  };

  static Tensor<1, dim, NumberType>
  d2psi_ds_dv(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType &                s)
  {
    return std::pow(det_t(t), 2) * 3.0 * std::pow(v_squ(v), 2) * dv_squ_dv(v) *
           sf * std::pow(s, sf - 1.0);
  };

  static Tensor<2, dim, NumberType>
  d2psi_dt_ds(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType &                s)
  {
    return 2.0 * std::pow(det_t(t), 1) * ddet_t_dt(t) * std::pow(v_squ(v), 3) *
           sf * std::pow(s, sf - 1.0);
  };

  static Tensor<1, dim, NumberType>
  d2psi_dv_ds(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType &                s)
  {
    return std::pow(det_t(t), 2) * 3.0 * std::pow(v_squ(v), 2) * dv_squ_dv(v) *
           sf * std::pow(s, sf - 1.0);
  };

  static NumberType
  d2psi_ds_ds(const Tensor<2, dim, NumberType> &t,
              const Tensor<1, dim, NumberType> &v,
              const NumberType &                s)
  {
    return std::pow(det_t(t), 2) * std::pow(v_squ(v), 3) * sf * (sf - 1.0) *
           std::pow(s, sf - 2.0);
  };
};

template <int dim, typename NumberType>
const double FunctionsTestTensorVectorScalarCoupled<dim, NumberType>::sf = 2.2;

template <int dim, typename number_t, enum AD::NumberTypes ad_type_code>
void
test_tensor_vector_scalar_coupled()
{
  typedef AD::ScalarFunction<dim, ad_type_code, number_t> ADHelper;
  typedef typename ADHelper::ad_type                      ADNumberType;
  typedef typename ADHelper::scalar_type                  ScalarNumberType;

  const unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  ConditionalOStream pcout(deallog.get_console(), this_mpi_process == 0);

  pcout << "*** Test variables: Tensor + Vector + Scalar (coupled), "
        << "dim = " << Utilities::to_string(dim) << ", "
        << "Type code: " << static_cast<int>(ad_type_code) << std::endl;

  // Values computed from the AD energy function
  ScalarNumberType             psi;
  Vector<ScalarNumberType>     Dpsi;
  FullMatrix<ScalarNumberType> D2psi;

  // Function and its derivatives
  typedef FunctionsTestTensorVectorScalarCoupled<dim, ADNumberType> func_ad;

  // Setup the variable components and choose a value at which to
  // evaluate the tape
  const FEValuesExtractors::Tensor<2> t_dof(0);
  const FEValuesExtractors::Vector    v_dof(
    Tensor<2, dim>::n_independent_components);
  const FEValuesExtractors::Scalar s_dof(
    Tensor<2, dim>::n_independent_components +
    Tensor<1, dim>::n_independent_components);
  const unsigned int n_AD_components =
    Tensor<2, dim>::n_independent_components +
    Tensor<1, dim>::n_independent_components + 1;
  ADHelper ad_helper(n_AD_components);
  ad_helper.set_tape_buffer_sizes(); // Increase the buffer size from the
                                     // default values

  Tensor<2, dim, ScalarNumberType> t =
    unit_symmetric_tensor<dim, ScalarNumberType>();
  for (unsigned int i = 0; i < t.n_independent_components; ++i)
    t[t.unrolled_to_component_indices(i)] += 0.11 * (i + 0.125);
  Tensor<1, dim, ScalarNumberType> v;
  ScalarNumberType                 s = 0.57;
  for (unsigned int i = 0; i < dim; ++i)
    v[i] = 0.275 * (1.0 + i);

  // Perform the differentiation operations
  {
    ad_helper.register_independent_variable(t, t_dof);
    ad_helper.register_independent_variable(v, v_dof);
    ad_helper.register_independent_variable(s, s_dof);

    const Tensor<2, dim, ADNumberType> t_ad =
      ad_helper.get_sensitive_variables(t_dof);
    const Tensor<1, dim, ADNumberType> v_ad =
      ad_helper.get_sensitive_variables(v_dof);
    ADNumberType s_ad = ad_helper.get_sensitive_variables(s_dof);

    const ADNumberType psi(func_ad::psi(t_ad, v_ad, s_ad));

    ad_helper.register_dependent_variable(psi);

    pcout << "Recorded data..." << std::endl;
    pcout << "independent variable values: " << std::flush;
    if (this_mpi_process == 0)
      ad_helper.print_values(pcout.get_stream());
    pcout << "t_ad: " << t_ad << std::endl;
    pcout << "v_ad: " << v_ad << std::endl;
    pcout << "s_ad: " << s_ad << std::endl;
    pcout << "psi: " << psi << std::endl;
    pcout << std::endl;
  }

  pcout << "independent variable values: " << std::flush;
  if (this_mpi_process == 0)
    ad_helper.print_values(pcout.get_stream());

  // Compute the function value, gradient and hessian for the new evaluation
  // point
  psi = ad_helper.compute_value();
  ad_helper.compute_gradient(Dpsi);
  if (AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2)
    {
      ad_helper.compute_hessian(D2psi);
    }

  // Output the full stored function, gradient vector and hessian matrix
  pcout << "psi: " << psi << std::endl;
  pcout << "Dpsi: \n";
  if (this_mpi_process == 0)
    Dpsi.print(pcout.get_stream());
  if (AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2)
    {
      pcout << "D2psi: \n";
      if (this_mpi_process == 0)
        D2psi.print_formatted(pcout.get_stream(), 3, true, 0, "0.0");
    }

  // Extract components of the solution
  const Tensor<2, dim, ScalarNumberType> dpsi_dt =
    ad_helper.extract_gradient_component(Dpsi, t_dof);
  const Tensor<1, dim, ScalarNumberType> dpsi_dv =
    ad_helper.extract_gradient_component(Dpsi, v_dof);
  const Tensor<0, dim, ScalarNumberType> dpsi_ds =
    ad_helper.extract_gradient_component(Dpsi, s_dof);
  pcout << "extracted Dpsi (t): " << dpsi_dt << "\n"
        << "extracted Dpsi (v): " << dpsi_dv << "\n"
        << "extracted Dpsi (s): " << dpsi_ds << "\n";

  // Verify the result
  typedef FunctionsTestTensorVectorScalarCoupled<dim, ScalarNumberType> func;
  static const ScalarNumberType                                         tol =
    1e5 * std::numeric_limits<ScalarNumberType>::epsilon();
  Assert(std::abs(psi - func::psi(t, v, s)) < tol,
         ExcMessage("No match for function value."));
  Assert(std::abs((dpsi_dt - func::dpsi_dt(t, v, s)).norm()) < tol,
         ExcMessage("No match for first derivative."));
  Assert(std::abs((dpsi_dv - func::dpsi_dv(t, v, s)).norm()) < tol,
         ExcMessage("No match for first derivative."));
  Assert(std::abs(dpsi_ds - func::dpsi_ds(t, v, s)) < tol,
         ExcMessage("No match for first derivative."));
  if (AD::ADNumberTraits<ADNumberType>::n_supported_derivative_levels >= 2)
    {
      const Tensor<4, dim, ScalarNumberType> d2psi_dt_dt =
        ad_helper.extract_hessian_component(D2psi, t_dof, t_dof);
      const Tensor<3, dim, ScalarNumberType> d2psi_dv_dt =
        ad_helper.extract_hessian_component(D2psi, t_dof, v_dof);
      const Tensor<2, dim, ScalarNumberType> d2psi_ds_dt =
        ad_helper.extract_hessian_component(D2psi, t_dof, s_dof);
      const Tensor<3, dim, ScalarNumberType> d2psi_dt_dv =
        ad_helper.extract_hessian_component(D2psi, v_dof, t_dof);
      const Tensor<2, dim, ScalarNumberType> d2psi_dv_dv =
        ad_helper.extract_hessian_component(D2psi, v_dof, v_dof);
      const Tensor<1, dim, ScalarNumberType> d2psi_ds_dv =
        ad_helper.extract_hessian_component(D2psi, v_dof, s_dof);
      const Tensor<2, dim, ScalarNumberType> d2psi_dt_ds =
        ad_helper.extract_hessian_component(D2psi, s_dof, t_dof);
      const Tensor<1, dim, ScalarNumberType> d2psi_dv_ds =
        ad_helper.extract_hessian_component(D2psi, s_dof, v_dof);
      const Tensor<0, dim, ScalarNumberType> d2psi_ds_ds =
        ad_helper.extract_hessian_component(D2psi, s_dof, s_dof);
      pcout << "extracted Dpsi (t): " << dpsi_dt << "\n"
            << "extracted Dpsi (v): " << dpsi_dv << "\n"
            << "extracted Dpsi (s): " << dpsi_ds << "\n"
            << "extracted D2psi (t,t): " << d2psi_dt_dt << "\n"
            << "extracted D2psi (t,v): " << d2psi_dv_dt << "\n"
            << "extracted D2psi (t,s): " << d2psi_ds_dt << "\n"
            << "extracted D2psi (v,t): " << d2psi_dt_dv << "\n"
            << "extracted D2psi (v,v): " << d2psi_dv_dv << "\n"
            << "extracted D2psi (v,s): " << d2psi_ds_dv << "\n"
            << "extracted D2psi (s,t): " << d2psi_dt_ds << "\n"
            << "extracted D2psi (s,v): " << d2psi_dv_ds << "\n"
            << "extracted D2psi (s,s): " << d2psi_ds_ds << "\n"
            << std::endl;
      Assert(std::abs((d2psi_dt_dt - func::d2psi_dt_dt(t, v, s)).norm()) < tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dv_dt - func::d2psi_dv_dt(t, v, s)).norm()) < tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_ds_dt - func::d2psi_ds_dt(t, v, s)).norm()) < tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dt_dv - func::d2psi_dt_dv(t, v, s)).norm()) < tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dv_dv - func::d2psi_dv_dv(t, v, s)).norm()) < tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_ds_dv - func::d2psi_ds_dv(t, v, s)).norm()) < tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dt_ds - func::d2psi_dt_ds(t, v, s)).norm()) < tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dv_ds - func::d2psi_dv_ds(t, v, s)).norm()) < tol,
             ExcMessage("No match for second derivative."));
      Assert(std::abs(d2psi_ds_ds - func::d2psi_ds_ds(t, v, s)) < tol,
             ExcMessage("No match for second derivative."));
    }
}
