/* $Id$ */
/* Authors: Andrea Bonito, Sebastian Pauletti. */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2010, 2011 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

				 // @sect3{Include files}

				 // If you've read through step-4 and step-7,
				 // you will recognize that we have used all
				 // of the following include files there
				 // already. Consequently, we will not explain
				 // their meaning here again.
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_generator.h>
#include <grid/grid_tools.h>
#include <lac/full_matrix.h>
#include <lac/vector.h>
#include <lac/solver_control.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/sparse_matrix.h>
#include <lac/compressed_sparsity_pattern.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <numerics/data_out.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>

#include <fstream>
#include <iostream>

using namespace dealii;

				 // @sect3{The <code>LaplaceBeltramiProblem</code> class template}

				 // This class is almost exactly similar to
				 // the <code>LaplaceProblem</code> class in
				 // step-4.

				 // The essential differences are these:
				 //
				 // - The template parameter now denotes the
				 //   dimensionality of the embedding space,
				 //   which is no longer the same as the
				 //   dimensionality of the domain and the
				 //   triangulation on which we compute. We
				 //   indicate this by calling the parameter
				 //   @p spacedim , and introducing a constant
				 //   @p dim equal to the dimensionality of
				 //   the domain -- here equal to
				 //   <code>spacedim-1</code>.
				 // - All member variables that have geometric
				 //   aspects now need to know about both
				 //   their own dimensionality as well as that
				 //   of the embedding space. Consequently, we
				 //   need to specify both of their template
				 //   parameters one for the dimension of the
				 //   mesh @p dim, and the other for the
				 //   dimension of the embedding space,
				 //   @p spacedim. This is exactly what we
				 //   did in step-34, take a look there for
				 //   a deeper explanation.

				 // - We need an object that describes which
				 //   kind of mapping to use from the
				 //   reference cell to the cells that the
				 //   triangulation is composed of. The
				 //   classes derived from the Mapping base
				 //   class do exactly this. Throughout most
				 //   of deal.II, if you don't do anything at
				 //   all, the library assumes that you want
				 //   an object of kind MappingQ1 that uses a
				 //   (bi-, tri-)linear mapping. In many
				 //   cases, this is quite sufficient, which
				 //   is why the use of these objects is
				 //   mostly optional: for example, if you
				 //   have a polygonal two-dimensional domain
				 //   in two-dimensional space, a bilinear
				 //   mapping of the reference cell to the
				 //   cells of the triangulation yields an
				 //   exact representation of the domain. If
				 //   you have a curved domain, one may want
				 //   to use a higher order mapping for those
				 //   cells that lie at the boundary of the
				 //   domain -- this is what we did in
				 //   step-11, for example. However, here we
				 //   have a curved domain, not just a curved
				 //   boundary, and while we can approximate
				 //   it with bilinearly mapped cells, it is
				 //   really only prodent to use a higher
				 //   order mapping for all
				 //   cells. Consequently, this class has a
				 //   member variable of type MappingQ; we
				 //   will choose the polynomial degree of the
				 //   mapping equal to the polynomial degree
				 //   of the finite element used in the
				 //   computations, though this
				 //   iso-parametricity is not necessary.
template <int spacedim>
class LaplaceBeltramiProblem 
{
  public:
    LaplaceBeltramiProblem (const unsigned degree = 2);
    void run ();
  
  private:
    static const unsigned int dim = spacedim-1;

    void make_grid_and_dofs ();
    void assemble_system ();
    void solve ();
    void output_results () const;
    void compute_error () const;
    
  
    Triangulation<dim,spacedim>   triangulation;
    FE_Q<dim,spacedim>            fe;
    DoFHandler<dim,spacedim>      dof_handler;
    MappingQ<dim, spacedim>       mapping;

    SparsityPattern               sparsity_pattern;
    SparseMatrix<double>          system_matrix;
  
    Vector<double>                solution;
    Vector<double>                system_rhs;
};


				 // @sect3{Equation data}

                                 // Next, let us define the classes that
                                 // describe the exact solution and the right
                                 // hand sides of the problem. This is in
                                 // analogy to step-4 and step-7 where we also
                                 // defined such objects. Given the discussion
                                 // in the introduction, the actual formulas
                                 // should be self-explanatory. A point of
                                 // interest may be how we define the value
                                 // and gradient functions for the 2d and 3d
                                 // cases separately, using explicit
                                 // specializations of the general
                                 // template. An alternative to doing it this
                                 // way might have been to define the general
                                 // template and have a <code>switch</code>
                                 // statement (or a sequence of
                                 // <code>if</code>s) for each possible value
                                 // of the spatial dimension.
template <int dim>
class Solution  : public Function<dim>
{
  public:
    Solution () : Function<dim>() {}
  
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
  
    virtual Tensor<1,dim> gradient (const Point<dim>   &p,
				    const unsigned int  component = 0) const;

};


template <>
double
Solution<2>::value (const Point<2> &p,
		    const unsigned int) const 
{
  return ( -2. * p(0) * p(1) );
}


template <>
Tensor<1,2>
Solution<2>::gradient (const Point<2>   &p,
		       const unsigned int) const
{
  Tensor<1,2> return_value;
  return_value[0] = -2. * p(1) * (1 - 2. * p(0) * p(0));
  return_value[1] = -2. * p(0) * (1 - 2. * p(1) * p(1));

  return return_value;
}


template <>
double
Solution<3>::value (const Point<3> &p,
		    const unsigned int) const 
{
  return (std::sin(numbers::PI * p(0)) *
	  std::cos(numbers::PI * p(1))*exp(p(2)));
}


template <>
Tensor<1,3>
Solution<3>::gradient (const Point<3>   &p,
		       const unsigned int) const
{
  using numbers::PI;

  Tensor<1,3> return_value;

  return_value[0] = PI *cos(PI * p(0))*cos(PI * p(1))*exp(p(2));
  return_value[1] = -PI *sin(PI * p(0))*sin(PI * p(1))*exp(p(2));
  return_value[2] = sin(PI * p(0))*cos(PI * p(1))*exp(p(2));
  
  return return_value;
}



template <int dim>
class RightHandSide : public Function<dim>
{
  public:
    RightHandSide () : Function<dim>() {}
  
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};

template <>
double
RightHandSide<2>::value (const Point<2> &p,
			 const unsigned int comp) const 
{
  return ( -8. * p(0) * p(1) ); 
}


template <>
double
RightHandSide<3>::value (const Point<3> &p,
			 const unsigned int comp) const 
{
  using numbers::PI;
  
  Tensor<2,3> hessian;

  hessian[0][0] = -PI*PI*sin(PI*p(0))*cos(PI*p(1))*exp(p(2));
  hessian[1][1] = -PI*PI*sin(PI*p(0))*cos(PI*p(1))*exp(p(2));
  hessian[2][2] = sin(PI*p(0))*cos(PI*p(1))*exp(p(2));

  hessian[0][1] = -PI*PI*cos(PI*p(0))*sin(PI*p(1))*exp(p(2));
  hessian[1][0] = -PI*PI*cos(PI*p(0))*sin(PI*p(1))*exp(p(2));

  hessian[0][2] = PI*cos(PI*p(0))*cos(PI*p(1))*exp(p(2));
  hessian[2][0] = PI*cos(PI*p(0))*cos(PI*p(1))*exp(p(2));

  hessian[1][2] = -PI*sin(PI*p(0))*sin(PI*p(1))*exp(p(2));
  hessian[2][1] = -PI*sin(PI*p(0))*sin(PI*p(1))*exp(p(2));

  Tensor<1,3> gradient;
  gradient[0] = PI * cos(PI*p(0))*cos(PI*p(1))*exp(p(2));
  gradient[1] = - PI * sin(PI*p(0))*sin(PI*p(1))*exp(p(2));
  gradient[2] = sin(PI*p(0))*cos(PI*p(1))*exp(p(2));

  Point<3> normal = p;
  normal /= p.norm();
   
  return (- trace(hessian)
	  - (2-3-1) * (gradient * normal)
	  + (hessian * normal) * normal);
}


                                 // @sect3{Implementation of the <code>LaplaceBeltramiProblem</code> class}

template <int spacedim>
LaplaceBeltramiProblem<spacedim>::
LaplaceBeltramiProblem (const unsigned degree)
		:
		fe (degree),
		dof_handler(triangulation),
		mapping(degree)
{}



template <int spacedim>
void LaplaceBeltramiProblem<spacedim>::make_grid_and_dofs ()
{
  Triangulation<spacedim> volume_mesh;
  GridGenerator::half_hyper_ball(volume_mesh);
  
  static HyperBallBoundary<dim,spacedim> surface_description;
  triangulation.set_boundary (0, surface_description);
  
  std::set<unsigned char> boundary_ids;
  boundary_ids.insert(0);
  
  GridTools::extract_boundary_mesh (volume_mesh, triangulation,
				    boundary_ids);
  triangulation.refine_global(4);

  std::cout << "Surface mesh has " << triangulation.n_active_cells()
	    << " cells."
	    << std::endl;

  dof_handler.distribute_dofs (fe);

  std::cout << "Surface mesh has " << dof_handler.n_dofs()
	    << " degrees of freedom."
	    << std::endl;
  
  CompressedSparsityPattern csp (dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, csp);
  sparsity_pattern.copy_from (csp);

  system_matrix.reinit (sparsity_pattern);
  
  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}


template <int spacedim>
void LaplaceBeltramiProblem<spacedim>::assemble_system () 
{
  system_matrix = 0;
  system_rhs = 0;
  
  QGauss<dim>  quadrature_formula(2);

  FEValues<dim,spacedim> fe_values (mapping, fe, quadrature_formula, 
				    update_values              |
				    update_gradients           |
				    update_quadrature_points   |
				    update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector< double > rhs_values(n_q_points);
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  const RightHandSide<spacedim> rhs;
  
  for (typename DoFHandler<dim,spacedim>::active_cell_iterator
	 cell = dof_handler.begin_active(),
	 endc = dof_handler.end(); cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit (cell);

      rhs.value_list (fe_values.get_quadrature_points(), rhs_values); 

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j) 
	  for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	    cell_matrix(i,j) 
	      += fe_values.shape_grad(i,q_point) *
	      fe_values.shape_grad(j,q_point) *
	      fe_values.JxW(q_point);
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	  cell_rhs(i) += fe_values.shape_value(i,q_point) *
			 rhs_values[q_point]*
			 fe_values.JxW(q_point);

      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    system_matrix.add (local_dof_indices[i],
			       local_dof_indices[j],
			       cell_matrix(i,j));
	
	  system_rhs(local_dof_indices[i]) += cell_rhs(i);
	}
    }


  std::map<unsigned int,double> boundary_values; 
  VectorTools::interpolate_boundary_values (mapping,
					    dof_handler,
					    0,
					    Solution<spacedim>(),
					    boundary_values);
  
  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      solution,
				      system_rhs,false);
}


template <int spacedim>
void LaplaceBeltramiProblem<spacedim>::solve () 
{
  SolverControl solver_control (solution.size(), 1e-7);
  SolverCG<>    cg (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);
}



template <int spacedim>
void LaplaceBeltramiProblem<spacedim>::output_results () const
{
  DataOut<dim,DoFHandler<dim,spacedim> > data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution",
			    DataOut<dim,DoFHandler<dim,spacedim> >::type_dof_data);
  data_out.build_patches (mapping,
			  mapping.get_degree());

  std::string filename ("solution-");
  filename += ('0'+spacedim);filename += "d.vtk"; 
  std::ofstream output (filename.c_str());
  data_out.write_vtk (output);
}



template <int spacedim>
void LaplaceBeltramiProblem<spacedim>::compute_error () const
{  
  Vector<float> difference_per_cell (triangulation.n_active_cells());
  VectorTools::integrate_difference (mapping, dof_handler, solution,
				     Solution<spacedim>(),
				     difference_per_cell,
				     QGauss<dim>(2*fe.degree+1),
				     VectorTools::H1_norm);  
  
  std::cout << "H1 error = "
	    << difference_per_cell.l2_norm()
	    << std::endl;
}



template <int spacedim>
void LaplaceBeltramiProblem<spacedim>::run () 
{
  make_grid_and_dofs();
  assemble_system ();
  solve ();
  output_results ();
  compute_error ();
}


                                 // @sect3{The main() function}

				 // The remainder of the program is taken up
				 // by the <code>main()</code> function. It
				 // follows exactly the general layout first
				 // introduced in step-6 and used in all
				 // following tutorial programs:
int main ()
{
  try
    {
      deallog.depth_console (0);

      LaplaceBeltramiProblem<3> laplace_beltrami;
      laplace_beltrami.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
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
