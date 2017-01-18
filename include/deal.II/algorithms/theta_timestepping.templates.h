// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2016 by the deal.II authors
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

#ifndef dealii__theta_timestepping_templates_h
#define dealii__theta_timestepping_templates_h


#include <deal.II/algorithms/theta_timestepping.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/vector_memory.h>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  template <typename VectorType>
  ThetaTimestepping<VectorType>::ThetaTimestepping (OperatorBase &e,
                                                    OperatorBase &i)
    :
    vtheta(0.5),
    adaptive(false),
    op_explicit(&e),
    op_implicit(&i)
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
    TimestepControl::declare_parameters (param);
    param.declare_entry("Theta", ".5", Patterns::Double());
    param.declare_entry("Adaptive", "false", Patterns::Bool());
    param.leave_subsection();
  }

  template <typename VectorType>
  void
  ThetaTimestepping<VectorType>::parse_parameters (ParameterHandler &param)
  {
    param.enter_subsection("ThetaTimestepping");
    control.parse_parameters (param);
    vtheta = param.get_double("Theta");
    adaptive = param.get_bool("Adaptive");
    param.leave_subsection ();
  }


  template <typename VectorType>
  void
  ThetaTimestepping<VectorType>::operator() (AnyData &out, const AnyData &in)
  {
    Assert(!adaptive, ExcNotImplemented());

    deallog.push ("Theta");

    VectorType &solution = *out.entry<VectorType *>(0);
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
    out1.add<VectorType *>(aux, "Solution");
    // The data provided to the inner solver
    src2.add<const VectorType *>(aux, "Previous time");
    src2.add<const VectorType *>(&solution, "Previous iterate");
    src2.add<const double *>(&d_implicit.time, "Time");
    src2.add<const double *>(&d_implicit.step, "Timestep");
    src2.add<const double *>(&vtheta, "Theta");
    src2.merge(in);

    if (output != 0)
      (*output) << 0U << out;

    for (unsigned int count = 1; d_explicit.time < control.final(); ++count)
      {
        const bool step_change = control.advance();
        d_implicit.time = control.now();
        d_explicit.step = (1.-vtheta)*control.step();
        d_implicit.step = vtheta*control.step();
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

        if (output != 0 && control.print())
          (*output) << count << out;

        d_explicit.time = control.now();
      }
    deallog.pop();
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
