// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/vector.h>

#include <deal.II/sundials/ida.h>

#include "../tests.h"

int
main()
{
  initlog();

  using VectorType = Vector<double>;
  SUNDIALS::IDA<VectorType>::AdditionalData data;
  data.initial_time      = 0.0;
  data.final_time        = 0.2;
  data.initial_step_size = 0.01;
  data.minimum_step_size = 0.02;
  data.output_period     = 0.1;
  data.ic_type = SUNDIALS::IDA<VectorType>::AdditionalData::use_y_diff;
  data.maximum_number_of_steps           = 100;
  data.maximum_step_size                 = 0.05;
  data.nonlinear_convergence_coefficient = 0.33;
  data.maximum_error_test_failures       = 7;

  ParameterHandler prm;
  data.add_parameters(prm);

  VectorType y(1), yp(1), constraints(1);
  y[0]           = 1.0;
  yp[0]          = 1.0;
  constraints[0] = 0.0;
  SUNDIALS::IDA<VectorType> ida(data);
  Assert(ida.get_ida_memory() == nullptr, ExcInternalError());
#if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
  Assert(ida.get_sun_context() != nullptr, ExcInternalError());
#endif

  ida.reinit_vector        = [](VectorType &v) { v.reinit(1); };
  unsigned int setup_count = 0;
  ida.residual             = [&](const double,
                     const VectorType &,
                     const VectorType &ydot,
                     VectorType       &r) {
    Assert(setup_count > 0, ExcInternalError());
    r = ydot;
  };
  double jacobian_scale = 1.0;
  ida.setup_jacobian    = [&](const double,
                           const VectorType &,
                           const VectorType &,
                           const double alpha) { jacobian_scale = alpha; };
  ida.solve_with_jacobian =
    [&](const VectorType &rhs, VectorType &dst, const double) {
      dst = rhs;
      dst /= jacobian_scale;
    };
  ida.get_constraint_vector = [&]() -> VectorType & { return constraints; };

  ida.custom_setup = [&](void *memory) {
    ++setup_count;
    void *user_data = nullptr;
    Assert(IDAGetUserData(memory, &user_data) == IDA_SUCCESS,
           ExcInternalError());
    Assert(user_data == &ida, ExcInternalError());
    Assert(memory == ida.get_ida_memory(), ExcInternalError());
  };

  ida.output_step = [](const double,
                       const VectorType &,
                       const VectorType &,
                       const unsigned int) {};
  ida.solve_dae(y, yp);
  Assert(setup_count == 1, ExcInternalError());
  Assert(ida.get_ida_memory() != nullptr, ExcInternalError());
  long int steps = -1;
  Assert(IDAGetNumSteps(ida.get_ida_memory(), &steps) == IDA_SUCCESS,
         ExcInternalError());
  Assert(steps > 0, ExcInternalError());
#if DEAL_II_SUNDIALS_VERSION_GTE(6, 2, 0)
  SUNDIALS::realtype last_step = 0.0;
  Assert(IDAGetLastStep(ida.get_ida_memory(), &last_step) == IDA_SUCCESS,
         ExcInternalError());
  Assert(last_step >= data.minimum_step_size, ExcInternalError());
#endif

  ida.reset(0.1, 0.01, y, yp);
  Assert(setup_count == 2, ExcInternalError());
  Assert(ida.get_ida_memory() != nullptr, ExcInternalError());
#if DEAL_II_SUNDIALS_VERSION_GTE(6, 0, 0)
  Assert(ida.get_sun_context() != nullptr, ExcInternalError());
#endif

  auto constrained_data       = data;
  constrained_data.ic_type    = SUNDIALS::IDA<VectorType>::AdditionalData::none;
  constrained_data.final_time = 0.8;
  VectorType y_constrained(1), yp_constrained(1), constraint_vector(1);
  y_constrained[0]     = 1.0;
  yp_constrained[0]    = -1.0;
  constraint_vector[0] = 2.0;
  SUNDIALS::IDA<VectorType> constrained_ida(constrained_data);
  double                    constrained_jacobian_scale = 1.0;
  constrained_ida.reinit_vector = [](VectorType &v) { v.reinit(1); };
  constrained_ida.residual      = [](const double,
                                const VectorType &,
                                const VectorType &ydot,
                                VectorType       &r) {
    r = ydot;
    r[0] += 1.0;
  };
  constrained_ida.setup_jacobian = [&](const double,
                                       const VectorType &,
                                       const VectorType &,
                                       const double alpha) {
    constrained_jacobian_scale = alpha;
  };
  constrained_ida.solve_with_jacobian =
    [&](const VectorType &rhs, VectorType &dst, const double) {
      dst = rhs;
      dst /= constrained_jacobian_scale;
    };
  constrained_ida.get_constraint_vector = [&]() -> VectorType & {
    return constraint_vector;
  };
  constrained_ida.output_step = [](const double,
                                   const VectorType &solution,
                                   const VectorType &,
                                   const unsigned int) {
    Assert(solution[0] > 0.0, ExcInternalError());
  };
  constrained_ida.solve_dae(y_constrained, yp_constrained);
  deallog << "IDA optional inputs and expert setup OK" << std::endl;
}
