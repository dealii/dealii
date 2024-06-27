// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// This is a modified version of step-44, which tests the implementation of
// QP-level symbolic-differentiation.
// The differentiation is performed inline with the update step.
// Compared to the run-time for step-44-sd-quadrature_level-[02,03].cc this
// serves only to gauge how quick/slow symbolic tensor differentiation is.

#include "../tests.h"

#include "sd_common_tests/step-44-sd-quadrature_level.h"

namespace Step44
{
  namespace SD = Differentiation::SD;

  template <int dim>
  class PointHistory
  {
  public:
    PointHistory()
      : F_inv(StandardTensors<dim>::I)
      , tau(SymmetricTensor<2, dim>())
      , d2Psi_vol_dJ2(0.0)
      , dPsi_vol_dJ(0.0)
      , Jc(SymmetricTensor<4, dim>())
    {}
    virtual ~PointHistory()
    {}
    void
    setup_lqp(const Parameters::AllParameters &parameters)
    {
      material.reset(
        new Material_Compressible_Neo_Hook_Three_Field<dim>(parameters.mu,
                                                            parameters.nu));
      update_values(Tensor<2, dim>(), 0.0, 1.0);
    }
    void
    update_values(const Tensor<2, dim> &Grad_u_n,
                  const double          p_tilde,
                  const double          J_tilde)
    {
      const Tensor<2, dim> F =
        (Tensor<2, dim>(StandardTensors<dim>::I) + Grad_u_n);
      material->update_material_data(F, p_tilde, J_tilde);
      F_inv = invert(F);

      namespace SD       = Differentiation::SD;
      using SDNumberType = SD::Expression;
      const bool debug   = false;

      // NOTE: We could actually perform all of the symbolic computations ONCE
      // up front and then reuse them as necessary. It is out of sheer laziness
      // that I'm not doing this here. Doing this might lead to a good
      // performance gain since we don't have to do the symbolic tensor
      // differentiation many times, but rather only once. The alternate
      // approach is adopted in step-44-sd_02.cc

      // Define some symbols
      const SymmetricTensor<2, dim, SDNumberType> C_SD(
        SD::make_symmetric_tensor_of_symbols<2, dim>("C"));
      const SDNumberType p_tilde_SD("p");
      const SDNumberType J_tilde_SD(
        SD::make_symbol("J")); // Check an unnecessary complication

      // Build the symbol substitution map
      // NOTE: I would typically leave this to the end (after symbol
      // calculations have been performed) but since it has to be used twice for
      // this problem so I do it up front
      const SymmetricTensor<2, dim> C = symmetrize(transpose(F) * F);
      SD::types::substitution_map   sub_vals_unresolved;
      SD::add_to_substitution_map(sub_vals_unresolved,
                                  SD::make_substitution_map(C_SD, C));
      SD::add_to_substitution_map(sub_vals_unresolved,
                                  SD::make_substitution_map(p_tilde_SD,
                                                            p_tilde));
      SD::add_to_substitution_map(sub_vals_unresolved,
                                  SD::make_substitution_map(J_tilde_SD,
                                                            J_tilde));
      // NOTE: The recursive substitution is not really required in this case,
      // but good to use in practise in case a more complex energy function is
      // employed later
      const SD::types::substitution_map sub_vals =
        SD::resolve_explicit_dependencies(sub_vals_unresolved);

      if (debug)
        {
          SD::Utilities::print_substitution_map(std::cout, sub_vals);
        }

      // Step 1: Update stress and material tangent
      {
        // 1a: Perform some symbolic calculations
        // Compute additional kinematic quantities
        // const SDNumberType det_F_SD = sqrt(determinant(C_SD));
        // const SymmetricTensor<2,dim,SDNumberType> C_bar_SD (pow(det_F_SD,
        // -2.0/dim) * C_SD);
        const SDNumberType det_C_SD = determinant(C_SD);
        const SymmetricTensor<2, dim, SDNumberType> C_bar_SD(
          pow(det_C_SD, -1.0 / dim) * C_SD);
        // Compute energy
        const SDNumberType det_F_SD     = sqrt(det_C_SD);
        const SDNumberType symbolic_psi = material->get_Psi_iso(C_bar_SD) +
                                          p_tilde_SD * (det_F_SD - J_tilde_SD);
        // Compute partial derivatives of energy function
        const SymmetricTensor<2, dim, SDNumberType> symbolic_S =
          2.0 * SD::differentiate(symbolic_psi, C_SD); // S = 2*dpsi_dC
        const SymmetricTensor<4, dim, SDNumberType> symbolic_H =
          2.0 *
          SD::differentiate(symbolic_S, C_SD); // H = 2*dS_dC = 4*d2psi_dC_dC

        if (debug)
          {
            std::cout << "symbolic_S: " << symbolic_S << std::endl;
            std::cout << "symbolic_H: " << symbolic_H << std::endl;
          }

        // 1b: Perform substitution of symbols
        const SymmetricTensor<2, dim, SDNumberType> substitution_S =
          SD::substitute(symbolic_S, sub_vals);
        const SymmetricTensor<4, dim, SDNumberType> substitution_H =
          SD::substitute(symbolic_H, sub_vals);

        if (debug)
          {
            std::cout << "substitution_S: " << substitution_S << std::endl;
            std::cout << "substitution_H: " << substitution_H << std::endl;
          }

        // Compute real-valued deal.II tensors
        const SymmetricTensor<2, dim> S =
          static_cast<SymmetricTensor<2, dim>>(substitution_S);
        const SymmetricTensor<4, dim> H =
          static_cast<SymmetricTensor<4, dim>>(substitution_H);

        // Push-forward
        tau = Physics::Transformations::Contravariant::push_forward(S, F);
        Jc  = Physics::Transformations::Contravariant::push_forward(H, F);
      }

      // Step 2: Update volumetric penalty terms
      {
        // 1a: Perform some symbolic calculations
        // Compute energy
        const SDNumberType symbolic_psi_vol = material->get_Psi_vol(J_tilde_SD);
        // Compute partial derivatives of energy function
        const SDNumberType symbolic_dPsi_vol_dJ =
          SD::differentiate(symbolic_psi_vol, J_tilde_SD);
        const SDNumberType symbolic_d2Psi_vol_dJ2 =
          SD::differentiate(symbolic_dPsi_vol_dJ, J_tilde_SD);

        if (debug)
          {
            std::cout << "symbolic_dPsi_vol_dJ: " << symbolic_dPsi_vol_dJ
                      << std::endl;
            std::cout << "symbolic_d2Psi_vol_dJ2: " << symbolic_d2Psi_vol_dJ2
                      << std::endl;
          }

        // 1b: Perform substitution of symbols
        const SDNumberType substitution_dPsi_vol_dJ =
          SD::substitute(symbolic_dPsi_vol_dJ, sub_vals);
        const SDNumberType substitution_d2Psi_vol_dJ2 =
          SD::substitute(symbolic_d2Psi_vol_dJ2, sub_vals);

        if (debug)
          {
            std::cout << "substitution_dPsi_vol_dJ: "
                      << substitution_dPsi_vol_dJ << std::endl;
            std::cout << "substitution_d2Psi_vol_dJ2: "
                      << substitution_d2Psi_vol_dJ2 << std::endl;
          }

        // Compute real values
        dPsi_vol_dJ   = static_cast<double>(substitution_dPsi_vol_dJ);
        d2Psi_vol_dJ2 = static_cast<double>(substitution_d2Psi_vol_dJ2);
      }

      if (debug)
        {
          std::cout << "dPsi_vol_dJ:                   " << dPsi_vol_dJ
                    << std::endl;
          std::cout << "material->get_dPsi_vol_dJ():   "
                    << material->get_dPsi_vol_dJ() << std::endl;
          std::cout << "d2Psi_vol_dJ2:                 " << d2Psi_vol_dJ2
                    << std::endl;
          std::cout << "material->get_d2Psi_vol_dJ2(): "
                    << material->get_d2Psi_vol_dJ2() << std::endl;

          std::cout << "tau:                 " << tau << std::endl;
          std::cout << "material->get_tau(): " << material->get_tau()
                    << std::endl;
          std::cout << "Jc:                  " << Jc << std::endl;
          std::cout << "material->get_Jc():  " << material->get_Jc()
                    << std::endl;
        }

      {
        static const double tol =
          1e-6; // Minor computation error due to order of operations

        // Zero strain --> zero stress
        if (std::abs(determinant(F) - 1.0) > 1e-9)
          {
            const SymmetricTensor<2, dim> tau_ref = material->get_tau();
            const SymmetricTensor<4, dim> Jc_ref  = material->get_Jc();

            Assert((tau - tau_ref).norm() / tau_ref.norm() < tol,
                   ExcMessage("SD computed stress is incorrect."));
            Assert((Jc - Jc_ref).norm() / Jc_ref.norm() < tol,
                   ExcMessage("SD computed tangent is incorrect."));
            Assert(std::abs((dPsi_vol_dJ - material->get_dPsi_vol_dJ()) /
                            material->get_dPsi_vol_dJ()) < tol,
                   ExcMessage("SD computed dPsi_vol_dJ is incorrect."));
            Assert(std::abs((d2Psi_vol_dJ2 - material->get_d2Psi_vol_dJ2()) /
                            material->get_d2Psi_vol_dJ2()) < tol,
                   ExcMessage("SD computed d2Psi_vol_dJ2 is incorrect."));
          }
        else
          {
            Assert((tau - material->get_tau()).norm() < tol,
                   ExcMessage("SD computed stress is incorrect."));
            Assert((Jc - material->get_Jc()).norm() < tol,
                   ExcMessage("SD computed tangent is incorrect."));
            Assert(std::abs(dPsi_vol_dJ - material->get_dPsi_vol_dJ()) < tol,
                   ExcMessage("SD computed dPsi_vol_dJ is incorrect."));
            Assert(std::abs(d2Psi_vol_dJ2 - material->get_d2Psi_vol_dJ2()) <
                     tol,
                   ExcMessage("SD computed d2Psi_vol_dJ2 is incorrect."));
          }
      }
    }
    double
    get_J_tilde() const
    {
      return material->get_J_tilde();
    }
    double
    get_det_F() const
    {
      return material->get_det_F();
    }
    const Tensor<2, dim> &
    get_F_inv() const
    {
      return F_inv;
    }
    double
    get_p_tilde() const
    {
      return material->get_p_tilde();
    }
    const SymmetricTensor<2, dim> &
    get_tau() const
    {
      return tau;
    }
    double
    get_dPsi_vol_dJ() const
    {
      return dPsi_vol_dJ;
    }
    double
    get_d2Psi_vol_dJ2() const
    {
      return d2Psi_vol_dJ2;
    }
    const SymmetricTensor<4, dim> &
    get_Jc() const
    {
      return Jc;
    }

  private:
    std::shared_ptr<Material_Compressible_Neo_Hook_Three_Field<dim>> material;
    Tensor<2, dim>                                                   F_inv;
    SymmetricTensor<2, dim>                                          tau;
    double                  d2Psi_vol_dJ2;
    double                  dPsi_vol_dJ;
    SymmetricTensor<4, dim> Jc;
  };
} // namespace Step44
int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  using namespace Step44;
  try
    {
      const unsigned int dim = 3;
      Solid<dim>         solid(SOURCE_DIR "/prm/parameters-step-44.prm");
      solid.run();
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
