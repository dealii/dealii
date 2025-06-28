// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_theta_timestepping_templates_h
#define dealii_theta_timestepping_templates_h


#include <deal.II/base/config.h>

#include <deal.II/algorithms/theta_timestepping.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/signaling_nan.h>

#include <deal.II/lac/vector_memory.h>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  template <typename VectorType>
  ThetaTimestepping<VectorType>::ThetaTimestepping(OperatorBase &e,
                                                   OperatorBase &i)
    : vtheta(0.5)
    , adaptive(false)
    , op_explicit(&e)
    , op_implicit(&i)
  {
    d_explicit.step = numbers::signaling_nan<double>();
    d_explicit.time = numbers::signaling_nan<double>();

    d_implicit.step = numbers::signaling_nan<double>();
    d_implicit.time = numbers::signaling_nan<double>();
  }


  template <typename VectorType>
  void
  ThetaTimestepping<VectorType>::notify(const Event &e)
  {
    op_explicit->notify(e);
    op_implicit->notify(e);
  }

  template <typename VectorType>
  void
  ThetaTimestepping<VectorType>::declare_parameters(ParameterHandler &param)
  {
    param.enter_subsection("ThetaTimestepping");
    TimestepControl::declare_parameters(param);
    param.declare_entry("Theta", ".5", Patterns::Double(0., 1.));
    param.declare_entry("Adaptive", "false", Patterns::Bool());
    param.leave_subsection();
  }

  template <typename VectorType>
  void
  ThetaTimestepping<VectorType>::parse_parameters(ParameterHandler &param)
  {
    param.enter_subsection("ThetaTimestepping");
    control.parse_parameters(param);
    vtheta   = param.get_double("Theta");
    adaptive = param.get_bool("Adaptive");
    param.leave_subsection();
  }


  template <typename VectorType>
  void
  ThetaTimestepping<VectorType>::operator()(AnyData &out, const AnyData &in)
  {
    Assert(!adaptive, ExcNotImplemented());

    LogStream::Prefix prefix("Theta");

    VectorType                     &solution = *out.entry<VectorType *>(0);
    GrowingVectorMemory<VectorType> mem;
    typename VectorMemory<VectorType>::Pointer aux(mem);
    aux->reinit(solution);

    control.restart();

    d_explicit.time = control.now();

    // The data used to compute the
    // vector associated with the old
    // timestep
    AnyData src1;
    src1.add<const VectorType *>(&solution, "Previous iterate");
    src1.add<const double *>(&d_explicit.time, "Time");
    src1.add<const double *>(&d_explicit.step, "Timestep");
    src1.add<const double *>(&vtheta, "Theta");
    src1.merge(in);

    AnyData src2;

    AnyData out1;
    out1.add<VectorType *>(aux.get(), "Solution");
    // The data provided to the inner solver
    src2.add<const VectorType *>(aux.get(), "Previous time");
    src2.add<const VectorType *>(&solution, "Previous iterate");
    src2.add<const double *>(&d_implicit.time, "Time");
    src2.add<const double *>(&d_implicit.step, "Timestep");
    src2.add<const double *>(&vtheta, "Theta");
    src2.merge(in);

    if (output != nullptr)
      (*output) << 0U << out;

    // When using nvcc 12.6 with C++20, the compiler has trouble determining the
    // type of d_explicit.time. We need to use a static_cast to help it.
    for (unsigned int count = 1;
         static_cast<double>(d_explicit.time) < control.final();
         ++count)
      {
        const bool step_change = control.advance();
        d_implicit.time        = control.now();
        d_explicit.step        = (1. - vtheta) * control.step();
        d_implicit.step        = vtheta * control.step();
        deallog << "Time step:" << d_implicit.time << std::endl;

        op_explicit->notify(Events::new_time);
        op_implicit->notify(Events::new_time);
        if (step_change)
          {
            op_explicit->notify(Events::new_timestep_size);
            op_implicit->notify(Events::new_timestep_size);
          }

        // Compute
        // (I + (1-theta)dt A) u
        (*op_explicit)(out1, src1);
        (*op_implicit)(out, src2);

        if (output != nullptr && control.print())
          (*output) << count << out;

        d_explicit.time = control.now();
      }
  }
} // namespace Algorithms

DEAL_II_NAMESPACE_CLOSE

#endif
