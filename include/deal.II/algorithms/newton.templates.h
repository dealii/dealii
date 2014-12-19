// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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


#include <deal.II/algorithms/newton.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector_memory.h>

#include <iomanip>


DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
  template <class VECTOR>
  Newton<VECTOR>::Newton(Operator<VECTOR> &residual, Operator<VECTOR> &inverse_derivative)
    :
    residual(&residual), inverse_derivative(&inverse_derivative),
    assemble_now(false),
    n_stepsize_iterations(21),
    assemble_threshold(0.),
    debug_vectors(false),
    debug(0)
  {}


  template <class VECTOR>
  void
  Newton<VECTOR>::declare_parameters(ParameterHandler &param)
  {
    param.enter_subsection("Newton");
    ReductionControl::declare_parameters (param);
    param.declare_entry("Assemble threshold", "0.", Patterns::Double());
    param.declare_entry("Stepsize iterations", "21", Patterns::Integer());
    param.declare_entry("Debug level", "0", Patterns::Integer());
    param.declare_entry("Debug vectors", "false", Patterns::Bool());
    param.leave_subsection();
  }

  template <class VECTOR>
  void
  Newton<VECTOR>::parse_parameters (ParameterHandler &param)
  {
    param.enter_subsection("Newton");
    control.parse_parameters (param);
    assemble_threshold = param.get_double("Assemble threshold");
    n_stepsize_iterations = param.get_integer("Stepsize iterations");
    debug_vectors = param.get_bool("Debug vectors");
    param.leave_subsection ();
  }

  template <class VECTOR>
  void
  Newton<VECTOR>::initialize (ParameterHandler &param)
  {
    parse_parameters(param);
  }

  template <class VECTOR>
  void
  Newton<VECTOR>::initialize (OutputOperator<VECTOR> &output)
  {
    data_out = &output;
  }

  template <class VECTOR>
  void
  Newton<VECTOR>::notify(const Event &e)
  {
    residual->notify(e);
    inverse_derivative->notify(e);
  }


  template <class VECTOR>
  double
  Newton<VECTOR>::threshold(const double thr)
  {
    const double t = assemble_threshold;
    assemble_threshold = thr;
    return t;
  }


  template <class VECTOR>
  void
  Newton<VECTOR>::operator() (NamedData<VECTOR *> &out, const NamedData<VECTOR *> &in)
  {
    Operator<VECTOR>::operator() (out, in);
  }

  template <class VECTOR>
  void
  Newton<VECTOR>::operator() (AnyData &out, const AnyData &in)
  {
    Assert (out.size() == 1, ExcNotImplemented());
    deallog.push ("Newton");

    VECTOR &u = *out.entry<VECTOR *>(0);

    if (debug>2)
      deallog << "u: " << u.l2_norm() << std::endl;

    GrowingVectorMemory<VECTOR> mem;
    typename VectorMemory<VECTOR>::Pointer Du(mem);
    typename VectorMemory<VECTOR>::Pointer res(mem);

    res->reinit(u);
    AnyData src1;
    AnyData src2;
    src1.add<const VECTOR *>(&u, "Newton iterate");
    src1.merge(in);
    src2.add<const VECTOR *>(res, "Newton residual");
    src2.merge(src1);
    AnyData out1;
    out1.add<VECTOR *>(res, "Residual");
    AnyData out2;
    out2.add<VECTOR *>(Du, "Update");

    unsigned int step = 0;
    // fill res with (f(u), v)
    (*residual)(out1, src1);
    double resnorm = res->l2_norm();
    double old_residual = resnorm / assemble_threshold + 1;

    if (debug_vectors)
      {
        NamedData<VECTOR *> out;
        VECTOR *p = &u;
        out.add(p, "solution");
        p = Du;
        out.add(p, "update");
        p = res;
        out.add(p, "residual");
        *data_out << step;
        *data_out << out;
      }

    while (control.check(step++, resnorm) == SolverControl::iterate)
      {
        // assemble (Df(u), v)
        if (resnorm/old_residual >= assemble_threshold)
          inverse_derivative->notify (Events::bad_derivative);

        Du->reinit(u);
        try
          {
            (*inverse_derivative)(out2, src2);
          }
        catch (SolverControl::NoConvergence &e)
          {
            deallog << "Inner iteration failed after "
                    << e.last_step << " steps with residual "
                    << e.last_residual << std::endl;
          }

        if (debug_vectors)
          {
            NamedData<VECTOR *> out;
            VECTOR *p = &u;
            out.add(p, "solution");
            p = Du;
            out.add(p, "update");
            p = res;
            out.add(p, "residual");
            *data_out << step;
            *data_out << out;
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
              deallog << "Trying step size: 1/" << (1<<step_size)
                      << " since residual was " << resnorm << std::endl;
            u.add(1./(1<<step_size), *Du);
            (*residual)(out1, src1);
            resnorm = res->l2_norm();
          }
      }
    deallog.pop();

    // in case of failure: throw exception
    if (control.last_check() != SolverControl::success)
      AssertThrow(false, SolverControl::NoConvergence (control.last_step(),
                                                       control.last_value()));
    // otherwise exit as normal
  }
}


DEAL_II_NAMESPACE_CLOSE
