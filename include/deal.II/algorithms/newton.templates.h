// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_newton_templates_h
#define dealii_newton_templates_h


#include <deal.II/base/config.h>

#include <deal.II/algorithms/newton.h>

#include <deal.II/base/logstream.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/vector_memory.h>

#include <iomanip>


DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  template <typename VectorType>
  Newton<VectorType>::Newton(OperatorBase &residual,
                             OperatorBase &inverse_derivative)
    : residual(&residual)
    , inverse_derivative(&inverse_derivative)
    , assemble_now(false)
    , n_stepsize_iterations(21)
    , assemble_threshold(0.)
    , debug_vectors(false)
    , debug(0)
  {}


  template <typename VectorType>
  void
  Newton<VectorType>::declare_parameters(ParameterHandler &param)
  {
    param.enter_subsection("Newton");
    ReductionControl::declare_parameters(param);
    param.declare_entry("Assemble threshold", "0.", Patterns::Double(0.));
    param.declare_entry("Stepsize iterations", "21", Patterns::Integer(0));
    param.declare_entry("Debug level", "0", Patterns::Integer(0));
    param.declare_entry("Debug vectors", "false", Patterns::Bool());
    param.leave_subsection();
  }

  template <typename VectorType>
  void
  Newton<VectorType>::parse_parameters(ParameterHandler &param)
  {
    param.enter_subsection("Newton");
    control.parse_parameters(param);
    assemble_threshold    = param.get_double("Assemble threshold");
    n_stepsize_iterations = param.get_integer("Stepsize iterations");
    debug                 = param.get_integer("Debug level");
    debug_vectors         = param.get_bool("Debug vectors");
    param.leave_subsection();
  }

  template <typename VectorType>
  void
  Newton<VectorType>::initialize(OutputOperator<VectorType> &output)
  {
    data_out = &output;
  }

  template <typename VectorType>
  void
  Newton<VectorType>::notify(const Event &e)
  {
    residual->notify(e);
    inverse_derivative->notify(e);
  }


  template <typename VectorType>
  double
  Newton<VectorType>::threshold(const double thr)
  {
    const double t     = assemble_threshold;
    assemble_threshold = thr;
    return t;
  }


  template <typename VectorType>
  void
  Newton<VectorType>::operator()(AnyData &out, const AnyData &in)
  {
    Assert(out.size() == 1, ExcNotImplemented());
    LogStream::Prefix prefix("Newton");

    VectorType &u = *out.entry<VectorType *>(0);

    if (debug > 2)
      deallog << "u: " << u.l2_norm() << std::endl;

    GrowingVectorMemory<VectorType>            mem;
    typename VectorMemory<VectorType>::Pointer Du(mem);
    typename VectorMemory<VectorType>::Pointer res(mem);

    Du->reinit(u);
    res->reinit(u);
    AnyData src1;
    AnyData src2;
    src1.add<const VectorType *>(&u, "Newton iterate");
    src1.merge(in);
    src2.add<const VectorType *>(res.get(), "Newton residual");
    src2.merge(src1);
    AnyData out1;
    out1.add<VectorType *>(res.get(), "Residual");
    AnyData out2;
    out2.add<VectorType *>(Du.get(), "Update");

    unsigned int step = 0;
    // fill res with (f(u), v)
    (*residual)(out1, src1);
    double resnorm      = res->l2_norm();
    double old_residual = 0.;

    if (debug_vectors)
      {
        AnyData     tmp;
        VectorType *p = &u;
        tmp.add<const VectorType *>(p, "solution");
        p = Du.get();
        tmp.add<const VectorType *>(p, "update");
        p = res.get();
        tmp.add<const VectorType *>(p, "residual");
        *data_out << step;
        *data_out << tmp;
      }

    while (control.check(step++, resnorm) == SolverControl::iterate)
      {
        // assemble (Df(u), v)
        if ((step > 1) && (resnorm / old_residual >= assemble_threshold))
          inverse_derivative->notify(Events::bad_derivative);

        Du->reinit(u);
        try
          {
            (*inverse_derivative)(out2, src2);
          }
        catch (const SolverControl::NoConvergence &e)
          {
            deallog << "Inner iteration failed after " << e.last_step
                    << " steps with residual " << e.last_residual << std::endl;
          }

        if (debug_vectors)
          {
            AnyData     tmp;
            VectorType *p = &u;
            tmp.add<const VectorType *>(p, "solution");
            p = Du.get();
            tmp.add<const VectorType *>(p, "update");
            p = res.get();
            tmp.add<const VectorType *>(p, "residual");
            *data_out << step;
            *data_out << tmp;
          }

        u.add(-1., *Du);
        old_residual = resnorm;
        (*residual)(out1, src1);
        resnorm = res->l2_norm();

        // Step size control
        unsigned int step_size = 0;
        while (resnorm >= old_residual)
          {
            ++step_size;
            if (step_size > n_stepsize_iterations)
              {
                deallog << "No smaller stepsize allowed!";
                break;
              }
            if (control.log_history())
              deallog << "Trying step size: 1/" << (1 << step_size)
                      << " since residual was " << resnorm << std::endl;
            u.add(1. / (1 << step_size), *Du);
            (*residual)(out1, src1);
            resnorm = res->l2_norm();
          }
      }

    // in case of failure: throw exception
    if (control.last_check() != SolverControl::success)
      AssertThrow(false,
                  SolverControl::NoConvergence(control.last_step(),
                                               control.last_value()));
    // otherwise exit as normal
  }
} // namespace Algorithms


DEAL_II_NAMESPACE_CLOSE

#endif
