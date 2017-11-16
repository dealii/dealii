// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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


// Evaluation of a coupled system (tensor + vector + scalar components)
// using a helper class

#include "../tests.h"
#include <deal.II/differentiation/ad.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iomanip>

using namespace dealii;
namespace AD = dealii::Differentiation::AD;

// Function and its derivatives
template<int dim, typename NumberType>
struct FunctionsTestTensorVectorScalarCoupled
{
  static NumberType
  det_t(const Tensor<2,dim,NumberType> &t)
  {
    return determinant(t);
  }

  static Tensor<2,dim,NumberType>
  ddet_t_dt(const Tensor<2,dim,NumberType> &t)
  {
    return det_t(t)*transpose(invert(t));
  }

  static Tensor<4,dim,NumberType>
  d2det_t_dt_dt(const Tensor<2,dim,NumberType> &t)
  {
    const Tensor<2,dim,NumberType> t_inv = invert(t);
    Tensor<4,dim,NumberType> dt_inv_trans_dt;
    // https://en.wikiversity.org/wiki/Introduction_to_Elasticity/Tensors#Derivative_of_the_determinant_of_a_tensor
    for (unsigned int i=0; i<dim; ++i)
      for (unsigned int j=0; j<dim; ++j)
        for (unsigned int k=0; k<dim; ++k)
          for (unsigned int l=0; l<dim; ++l)
            dt_inv_trans_dt[i][j][k][l] = -t_inv[l][i]*t_inv[j][k];

    return det_t(t)*outer_product(transpose(t_inv),transpose(t_inv))
           + det_t(t)*dt_inv_trans_dt;
  }

  static NumberType
  v_squ(const Tensor<1,dim,NumberType> &v)
  {
    return v*v;
  }

  static Tensor<1,dim,NumberType>
  dv_squ_dv(const Tensor<1,dim,NumberType> &v)
  {
    return 2.0*v;
  }

  static Tensor<2,dim,NumberType>
  d2v_squ_dv_dv(const Tensor<1,dim,NumberType> &v)
  {
    static const Tensor<2,dim,NumberType> I (unit_symmetric_tensor<dim>());
    return 2.0*I;
  }

  // --------

  static const double sf;

  static NumberType
  psi (const Tensor<2,dim,NumberType> &t, const Tensor<1,dim,NumberType> &v, const NumberType &s)
  {
    return pow(det_t(t),2)*pow(v_squ(v),3)*pow(s,sf);
  };

  static Tensor<2,dim,NumberType>
  dpsi_dt (const Tensor<2,dim,NumberType> &t, const Tensor<1,dim,NumberType> &v, const NumberType &s)
  {
    return 2.0*pow(det_t(t),1)*ddet_t_dt(t)*pow(v_squ(v),3)*pow(s,sf);
  };

  static Tensor<1,dim,NumberType>
  dpsi_dv (const Tensor<2,dim,NumberType> &t, const Tensor<1,dim,NumberType> &v, const NumberType &s)
  {
    return pow(det_t(t),2)*3.0*pow(v_squ(v),2)*dv_squ_dv(v)*pow(s,sf);
  };

  static NumberType
  dpsi_ds (const Tensor<2,dim,NumberType> &t, const Tensor<1,dim,NumberType> &v, const NumberType &s)
  {
    return pow(det_t(t),2)*pow(v_squ(v),3)*sf*pow(s,sf-1.0);
  };

  static Tensor<4,dim,NumberType>
  d2psi_dt_dt (const Tensor<2,dim,NumberType> &t, const Tensor<1,dim,NumberType> &v, const NumberType &s)
  {
    return 2.0*pow(v_squ(v),3)*
           (pow(det_t(t),0)*outer_product(ddet_t_dt(t),ddet_t_dt(t))
            + pow(det_t(t),1)*d2det_t_dt_dt(t))*pow(s,sf);
  };

  static Tensor<3,dim,NumberType>
  d2psi_dv_dt (const Tensor<2,dim,NumberType> &t, const Tensor<1,dim,NumberType> &v, const NumberType &s)
  {
    return 2.0*pow(det_t(t),1)*3.0*pow(v_squ(v),2)*outer_product(ddet_t_dt(t),dv_squ_dv(v))*pow(s,sf);
  };

  static Tensor<2,dim,NumberType>
  d2psi_ds_dt (const Tensor<2,dim,NumberType> &t, const Tensor<1,dim,NumberType> &v, const NumberType &s)
  {
    return 2.0*pow(det_t(t),1)*ddet_t_dt(t)*pow(v_squ(v),3)*sf*pow(s,sf-1.0);
  };

  static Tensor<3,dim,NumberType>
  d2psi_dt_dv (const Tensor<2,dim,NumberType> &t, const Tensor<1,dim,NumberType> &v, const NumberType &s)
  {
    return 2.0*pow(det_t(t),1)*3.0*pow(v_squ(v),2)*outer_product(dv_squ_dv(v),ddet_t_dt(t))*pow(s,sf);
  };

  static Tensor<2,dim,NumberType>
  d2psi_dv_dv (const Tensor<2,dim,NumberType> &t, const Tensor<1,dim,NumberType> &v, const NumberType &s)
  {
    return pow(det_t(t),2)*3.0*
           (2.0*pow(v_squ(v),1)*outer_product(dv_squ_dv(v),dv_squ_dv(v))
            + pow(v_squ(v),2)*d2v_squ_dv_dv(v))*pow(s,sf);
  };

  static Tensor<1,dim,NumberType>
  d2psi_ds_dv (const Tensor<2,dim,NumberType> &t, const Tensor<1,dim,NumberType> &v, const NumberType &s)
  {
    return pow(det_t(t),2)*3.0*pow(v_squ(v),2)*dv_squ_dv(v)*sf*pow(s,sf-1.0);
  };

  static Tensor<2,dim,NumberType>
  d2psi_dt_ds (const Tensor<2,dim,NumberType> &t, const Tensor<1,dim,NumberType> &v, const NumberType &s)
  {
    return 2.0*pow(det_t(t),1)*ddet_t_dt(t)*pow(v_squ(v),3)*sf*pow(s,sf-1.0);
  };

  static Tensor<1,dim,NumberType>
  d2psi_dv_ds (const Tensor<2,dim,NumberType> &t, const Tensor<1,dim,NumberType> &v, const NumberType &s)
  {
    return pow(det_t(t),2)*3.0*pow(v_squ(v),2)*dv_squ_dv(v)*sf*pow(s,sf-1.0);
  };

  static NumberType
  d2psi_ds_ds (const Tensor<2,dim,NumberType> &t, const Tensor<1,dim,NumberType> &v, const NumberType &s)
  {
    return pow(det_t(t),2)*pow(v_squ(v),3)*sf*(sf-1.0)*pow(s,sf-2.0);
  };
};

template<int dim, typename NumberType>
const double
FunctionsTestTensorVectorScalarCoupled<dim,NumberType>::sf = 2.2;

template<int dim, typename number_t, enum AD::NumberTypes ad_type_code>
void test_tensor_vector_scalar_coupled ()
{
  typedef AD::ADHelperScalarFunction<dim,ad_type_code,number_t> ADHelper;
  typedef typename ADHelper::ad_type ADNumberType;

  const unsigned int this_mpi_process = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  ConditionalOStream pcout (deallog.get_console(), this_mpi_process==0);

  pcout
      << "*** Test variables: Tensor + Vector + Scalar (coupled), "
      << (AD::ADNumberTraits<ADNumberType>::is_taped == true ? "Taped" : "Tapeless")
      << std::endl;

  // Values computed from the AD energy function
  double psi;
  Vector<double> Dpsi;
  FullMatrix<double> D2psi;

  // Function and its derivatives
  typedef FunctionsTestTensorVectorScalarCoupled<dim,ADNumberType> func_ad;

  // Setup the variable components and choose a value at which to
  // evaluate the tape
  const FEValuesExtractors::Tensor<2> t_dof (0);
  const FEValuesExtractors::Vector    v_dof (Tensor<2,dim>::n_independent_components);
  const FEValuesExtractors::Scalar    s_dof (Tensor<2,dim>::n_independent_components
                                             + Tensor<1,dim>::n_independent_components);
  const unsigned int n_AD_components = Tensor<2,dim>::n_independent_components
                                       + Tensor<1,dim>::n_independent_components
                                       + 1;
  ADHelper ad_helper (n_AD_components);
  Tensor<2,dim> t = unit_symmetric_tensor<dim>();
  for (unsigned int i=0; i<t.n_independent_components; ++i)
    t[t.unrolled_to_component_indices(i)] += 0.11*(i+0.125);
  Tensor<1,dim> v;
  double s = 0.57;
  for (unsigned int i=0; i<dim; ++i)
    v[i] = 0.275*(1.0+i);

  const int tape_no = 1;
  const bool is_recording = ad_helper.enable_record_sensitivities(tape_no /*material_id*/,
                            true /*overwrite_tape*/,
                            true /*keep*/);
  if (is_recording == true)
    {
      ad_helper.register_independent_variable(t, t_dof);
      ad_helper.register_independent_variable(v, v_dof);
      ad_helper.register_independent_variable(s, s_dof);

      const Tensor<2,dim,ADNumberType> t_ad = ad_helper.get_sensitive_variables(t_dof);
      const Tensor<1,dim,ADNumberType> v_ad = ad_helper.get_sensitive_variables(v_dof);
      ADNumberType s_ad = ad_helper.get_sensitive_variables(s_dof);

      const ADNumberType psi (func_ad::psi(t_ad,v_ad,s_ad));

      ad_helper.register_dependent_variable(psi);
      ad_helper.disable_record_sensitivities(false /*write_tapes_to_file*/);

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
  else
    {
      Assert(is_recording==true, ExcInternalError());
    }

  // Do some work :-)
  // Set a new evaluation point
  if (AD::ADNumberTraits<ADNumberType>::is_taped == true)
    {
      pcout << "Using tape with different values for independent variables..." << std::endl;
      ad_helper.activate_tape(tape_no);
      t *= 0.9;
      v *= 0.63;
      s *= 1.21;
      ad_helper.set_independent_variable(t, t_dof);
      ad_helper.set_independent_variable(v, v_dof);
      ad_helper.set_independent_variable(s, s_dof);
    }

  pcout << "independent variable values: " << std::flush;
  if (this_mpi_process == 0)
    ad_helper.print_values(pcout.get_stream());

  // Compute the function value, gradient and hessian for the new evaluation point
  psi = ad_helper.compute_value();
  Dpsi = ad_helper.compute_gradient();
  if (AD::ADNumberTraits<ADNumberType>::is_taped == true)
    {
      D2psi = ad_helper.compute_hessian();
    }

  // Output the full stored function, gradient vector and hessian matrix
  pcout << "psi: " << psi << std::endl;
  pcout << "Dpsi: \n";
  if (this_mpi_process == 0)
    Dpsi.print(pcout.get_stream());
  if (AD::ADNumberTraits<ADNumberType>::is_taped == true)
    {
      pcout << "D2psi: \n";
      if (this_mpi_process == 0)
        D2psi.print_formatted(pcout.get_stream(),3,true,0,"0.0");
    }

  // Extract components of the solution
  const Tensor<2,dim,double> dpsi_dt = ad_helper.extract_gradient_component(Dpsi,t_dof);
  const Tensor<1,dim,double> dpsi_dv = ad_helper.extract_gradient_component(Dpsi,v_dof);
  const Tensor<0,dim,double> dpsi_ds = ad_helper.extract_gradient_component(Dpsi,s_dof);
  pcout
      << "extracted Dpsi (t): " << dpsi_dt << "\n"
      << "extracted Dpsi (v): " << dpsi_dv << "\n"
      << "extracted Dpsi (s): " << dpsi_ds << "\n";

  // Verify the result
  typedef FunctionsTestTensorVectorScalarCoupled<dim,double> func;
  static const double tol = 1e-12;
  Assert(std::abs(psi - func::psi(t,v,s)) < tol, ExcMessage("No match for function value."));
  Assert(std::abs((dpsi_dt - func::dpsi_dt(t,v,s)).norm()) < tol, ExcMessage("No match for first derivative."));
  Assert(std::abs((dpsi_dv - func::dpsi_dv(t,v,s)).norm()) < tol, ExcMessage("No match for first derivative."));
  Assert(std::abs( dpsi_ds - func::dpsi_ds(t,v,s))         < tol, ExcMessage("No match for first derivative."));
  if (AD::ADNumberTraits<ADNumberType>::is_taped == true)
    {
      const Tensor<4,dim,double> d2psi_dt_dt = ad_helper.extract_hessian_component(D2psi,t_dof,t_dof);
      const Tensor<3,dim,double> d2psi_dv_dt = ad_helper.extract_hessian_component(D2psi,t_dof,v_dof);
      const Tensor<2,dim,double> d2psi_ds_dt = ad_helper.extract_hessian_component(D2psi,t_dof,s_dof);
      const Tensor<3,dim,double> d2psi_dt_dv = ad_helper.extract_hessian_component(D2psi,v_dof,t_dof);
      const Tensor<2,dim,double> d2psi_dv_dv = ad_helper.extract_hessian_component(D2psi,v_dof,v_dof);
      const Tensor<1,dim,double> d2psi_ds_dv = ad_helper.extract_hessian_component(D2psi,v_dof,s_dof);
      const Tensor<2,dim,double> d2psi_dt_ds = ad_helper.extract_hessian_component(D2psi,s_dof,t_dof);
      const Tensor<1,dim,double> d2psi_dv_ds = ad_helper.extract_hessian_component(D2psi,s_dof,v_dof);
      const Tensor<0,dim,double> d2psi_ds_ds = ad_helper.extract_hessian_component(D2psi,s_dof,s_dof);
      pcout
          << "extracted Dpsi (t): " << dpsi_dt << "\n"
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
      Assert(std::abs((d2psi_dt_dt - func::d2psi_dt_dt(t,v,s)).norm()) < tol, ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dv_dt - func::d2psi_dv_dt(t,v,s)).norm()) < tol, ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_ds_dt - func::d2psi_ds_dt(t,v,s)).norm()) < tol, ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dt_dv - func::d2psi_dt_dv(t,v,s)).norm()) < tol, ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dv_dv - func::d2psi_dv_dv(t,v,s)).norm()) < tol, ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_ds_dv - func::d2psi_ds_dv(t,v,s)).norm()) < tol, ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dt_ds - func::d2psi_dt_ds(t,v,s)).norm()) < tol, ExcMessage("No match for second derivative."));
      Assert(std::abs((d2psi_dv_ds - func::d2psi_dv_ds(t,v,s)).norm()) < tol, ExcMessage("No match for second derivative."));
      Assert(std::abs( d2psi_ds_ds - func::d2psi_ds_ds(t,v,s))         < tol, ExcMessage("No match for second derivative."));
    }
}

int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

  const unsigned int this_mpi_process = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  if (this_mpi_process == 0)
    {
      initlog();

      // --- Taped ---
      test_tensor_vector_scalar_coupled<2,double,AD::NumberTypes::adolc_taped>();
      test_tensor_vector_scalar_coupled<3,double,AD::NumberTypes::adolc_taped>();
      deallog << "Taped OK" << std::endl;

      // --- Tapeless ---
      test_tensor_vector_scalar_coupled<2,double,AD::NumberTypes::adolc_tapeless>();
      test_tensor_vector_scalar_coupled<3,double,AD::NumberTypes::adolc_tapeless>();
      deallog << "Tapeless OK" << std::endl;

      deallog << "OK" << std::endl;
    }
  else
    {
      // --- Taped ---
      test_tensor_vector_scalar_coupled<2,double,AD::NumberTypes::adolc_taped>();
      test_tensor_vector_scalar_coupled<3,double,AD::NumberTypes::adolc_taped>();

      // --- Tapeless ---
      test_tensor_vector_scalar_coupled<2,double,AD::NumberTypes::adolc_tapeless>();
      test_tensor_vector_scalar_coupled<3,double,AD::NumberTypes::adolc_tapeless>();
    }
}
