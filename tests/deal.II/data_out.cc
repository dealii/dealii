//----------------------------  data_out.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  data_out.cc  ---------------------------


/* Author: Wolfgang Bangerth, University of Heidelberg, 1999 */
/* adapted from step-4. */


#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/continuous.h>
#include <fe/mapping_q1.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>

#include <numerics/data_out.h>
#include <numerics/data_out_rotation.h>
#include <numerics/data_out_faces.h>
#include <fstream>

#include <base/logstream.h>


ofstream logfile("data_out.output");


template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem ();
    void run ();
    
  private:
    void make_grid_and_dofs ();
    void solve ();
    void output_results () const;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
};


template <int dim>
class RightHandSide : public Function<dim> 
{
  public:
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};



template <int dim>
class BoundaryValues : public Function<dim> 
{
  public:
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};




template <int dim>
double RightHandSide<dim>::value (const Point<dim> &p,
				  const unsigned int) const 
{
  double return_value = 0;
  for (unsigned int i=0; i<dim; ++i)
    return_value += 4*pow(p(i), 4);

  return return_value;
};


template <int dim>
double BoundaryValues<dim>::value (const Point<dim> &p,
				   const unsigned int) const 
{
  return p.square();
};




template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
  fe (1), dof_handler (triangulation)
{};



template <int dim>
void LaplaceProblem<dim>::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, 0, 1);
  triangulation.refine_global (dim==2 ? 2 : 1);
  for (unsigned int i=0; i<2; ++i)
    {
      triangulation.begin_active()->set_refine_flag ();
      triangulation.execute_coarsening_and_refinement ();
    };
  
  
  logfile << "   Number of active cells: "
	  << triangulation.n_active_cells()
	  << endl
	  << "   Total number of cells: "
	  << triangulation.n_cells()
	  << endl;

  dof_handler.distribute_dofs (fe);

  logfile << "   Number of degrees of freedom: "
	  << dof_handler.n_dofs()
	  << endl;

  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
};




template <int dim>
void LaplaceProblem<dim>::solve () 
{
				   // dummy solve. just insert some
				   // arbitrary values
  for (unsigned int i=0; i<solution.size(); ++i)
    solution(i) = i;
};



template <>
void LaplaceProblem<2>::output_results () const
{
  const unsigned int dim = 2;

				   // test regular output in 2d
  if (true)
    {
      DataOut<dim> data_out;
      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (solution, "solution");
      data_out.build_patches ();
      data_out.write_gmv (logfile);
      data_out.write_gnuplot (logfile);
      data_out.write_ucd (logfile);
      data_out.write_povray (logfile);
      data_out.write_eps (logfile);
    };

				   // test DataOutRotation in 2d
  if (true)
    {
      DataOutRotation<dim> data_out;
      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (solution, "solution");
      data_out.build_patches (10);
      data_out.write_gmv (logfile);
      data_out.write_gnuplot (logfile);
      data_out.write_ucd (logfile);
    };
};



template <>
void LaplaceProblem<3>::output_results () const
{
  const unsigned int dim = 3;

				   // test regular output in 3d
  if (true)
    {
      DataOut<dim> data_out;
      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (solution, "solution");
      data_out.build_patches ();
      data_out.write_gmv (logfile);
      data_out.write_gnuplot (logfile);
      data_out.write_ucd (logfile);
    };

				   // test DataOutFaces in 3d. note:
				   // not all output formats support
				   // this
  if (true)
    {
      DataOutFaces<dim> data_out;
      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (solution, "solution");
      data_out.build_patches (10);
      data_out.write_gmv (logfile);
      data_out.write_gnuplot (logfile);
      data_out.write_ucd (logfile);
    };
};



template <int dim>
void LaplaceProblem<dim>::run () 
{
  make_grid_and_dofs();
  solve ();
  output_results ();
};

    

int main () 
{
  deallog.depth_console (0);

  LaplaceProblem<2> laplace_problem_2d;
  laplace_problem_2d.run ();

  LaplaceProblem<3> laplace_problem_3d;
  laplace_problem_3d.run ();
  
  return 0;
};
