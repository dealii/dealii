// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2000 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_dof_print_solver_step_h
#define dealii_dof_print_solver_step_h

#include <deal.II/base/config.h>

#include <deal.II/base/logstream.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iomanip>
#include <sstream>

DEAL_II_NAMESPACE_OPEN


/**
 * Print intermediate solutions in solvers.  This is derived from a solver
 * class provided as template argument.  It implements the @p print_vector
 * function of the solver using a DoFHandler. This way, the intermediate
 * vectors can be viewed as finite element functions. This class might be used
 * first to understand how solvers work (for example to visualize the
 * smoothing properties of various solvers, e.g. in a multigrid context), and
 * second to investigate why and how a solver fails to solve certain classes
 * of problems.
 *
 * Objects of this class are provided with a solver class through a template
 * argument, and with a file name (as a string), with which a new file is
 * constructed in each iteration (named <tt>basename.[step].[suffix]</tt>) and
 * into which the solution is written as a finite element field using the
 * DataOut class. Please note that this class may produce enormous amounts of
 * data!
 *
 * @ingroup output
 */
template <int dim, typename SolverType, typename VectorType = Vector<double>>
class DoFPrintSolverStep : public SolverType
{
public:
  /**
   * Constructor.  First, we take the arguments needed for the solver. @p
   * data_out is the object doing the output as a finite element function.
   *
   * One output file with the name <tt>basename.[step].[suffix]</tt> will be
   * produced for each iteration step.
   */
  DoFPrintSolverStep(SolverControl            &control,
                     VectorMemory<VectorType> &mem,
                     DataOut<dim>             &data_out,
                     const std::string        &basename);

  /**
   * Call-back function for the iterative method.
   */
  virtual void
  print_vectors(const unsigned int step,
                const VectorType  &x,
                const VectorType  &r,
                const VectorType  &d) const;

private:
  /**
   * Output object.
   */
  DataOut<dim> &out;

  /**
   * Base of filenames.
   */
  const std::string basename;
};


/* ----------------------- template functions --------------- */

template <int dim, typename SolverType, typename VectorType>
DoFPrintSolverStep<dim, SolverType, VectorType>::DoFPrintSolverStep(
  SolverControl            &control,
  VectorMemory<VectorType> &mem,
  DataOut<dim>             &data_out,
  const std::string        &basename)
  : SolverType(control, mem)
  , out(data_out)
  , basename(basename)
{}


template <int dim, typename SolverType, typename VectorType>
void
DoFPrintSolverStep<dim, SolverType, VectorType>::print_vectors(
  const unsigned int step,
  const VectorType  &x,
  const VectorType  &r,
  const VectorType  &d) const
{
  out.clear_data_vectors();
  out.add_data_vector(x, "solution");
  out.add_data_vector(r, "residual");
  out.add_data_vector(d, "update");

  std::ostringstream filename;
  filename << basename << std::setw(3) << std::setfill('0') << step
           << out.default_suffix();

  const std::string fname = filename.str();

  deallog << "Writing file:" << fname << std::endl;

  out.build_patches();
  std::ofstream of(fname);
  out.write(of);
}

DEAL_II_NAMESPACE_CLOSE

#endif
