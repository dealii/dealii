//----------------------------  solver_cg.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  solver_cg.h  ---------------------------
#ifndef __deal2__dof_print_solver_step_h
#define __deal2__dof_print_solver_step_h

#include <base/config.h>
#include <base/logstream.h>
#include <numerics/data_out.h>

#include <strstream>
#include <iomanip>
#include <fstream>



/**
 * Print intermediate solutions in solvers.  This is derived from a
 * solver class provided as template argument.  It implements the
 * @p{print_vector} function of the solver using a
 * @ref{DoFHandler}. This way, the intermediate vectors can be viewed
 * as finite element functions. This class might be used first to
 * understand how solvers work (for example to visualize the smoothing
 * properties of various solvers, e.g. in a multigrid context), and
 * second to investigate why and how a solver fails to solve certain
 * classes of problems.
 *
 * Objects of this class are provided with a solver class through a
 * template argument, and with a file name (as a string), with which a
 * new file is constructed in each iteration (named
 * @p{basename.[step].[suffix]}) and into which the solution is
 * written as a finite element field using the @ref{DataOut} class.
 * Please note that this class may produce enormous amounts of data!
 *
 * @author Guido Kanschat, 2000
 */
template<int dim, class SOLVER, class VECTOR = Vector<double> >
class DoFPrintSolverStep : public SOLVER
{
  public:
				     /**
				      * Constructor.  First, we take
				      * the arguments needed for the
				      * solver. @p{data_out} is the
				      * object doing the output as a
				      * finite element function.
				      *
				      * One output file with the name
				      * @p{basename.[step].[suffix]}
				      * will be produced for each
				      * iteration step.
				      */
    DoFPrintSolverStep (SolverControl& control,
			VectorMemory<VECTOR>& mem,
			DataOut<dim>& data_out,
			const std::string& basename);

				     /**
				      * Call-back function for the
				      * iterative method.
				      */
    virtual void print_vectors (const unsigned int step,
				const VECTOR& x,
				const VECTOR& r,
				const VECTOR& d) const;
  private:
				     /**
				      * Output object.
				      */
    DataOut<dim>& out;

				     /**
				      * Base of filenames.
				      */
    const std::string basename;
};


/* ----------------------- template functions --------------- */

template<int dim, class SOLVER, class VECTOR>
DoFPrintSolverStep<dim, SOLVER, VECTOR>::DoFPrintSolverStep (SolverControl& control,
							     VectorMemory<VECTOR>& mem,
							     DataOut<dim>& data_out,
							     const std::string& basename)
		: SOLVER (control, mem),
		  out (data_out),
		  basename (basename)
{}


template<int dim, class SOLVER, class VECTOR>
void
DoFPrintSolverStep<dim, SOLVER, VECTOR>::print_vectors (const unsigned int step,
							const VECTOR& x,
							const VECTOR& r,
							const VECTOR& d) const
{
  out.clear_data_vectors();
  out.add_data_vector(x, "solution");
  out.add_data_vector(r, "residual");
  out.add_data_vector(d, "update");

  ostrstream filename;
  filename << basename
	   << setw(3) << setfill('0') << step
	   << out.default_suffix() << ends;

  const std::string fname = filename.str();

  deallog << "Writing file:" << fname << std::endl;

  out.build_patches();
  std::ofstream of (fname.c_str());
  out.write (of);
}

#endif
