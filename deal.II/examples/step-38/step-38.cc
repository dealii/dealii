/* $Id$ */
/* Authors: Andrea Bonito, Sebastian Pauletti. */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2010 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


#include <base/quadrature_lib.h>
#include <base/function.h>

#include <lac/full_matrix.h>
#include <lac/vector.h>
#include <lac/solver_control.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/sparse_matrix.h>
#include <lac/compressed_sparsity_pattern.h>

#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_generator.h>
#include <grid/grid_tools.h>

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
RightHandSide<3>::value (const Point<3> &p,
			 const unsigned int comp) const 
{
  using numbers::PI;

				   // LB: u = Delta u - nu D2 u nu - (Grad u nu ) div (nu)

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

  double curvature;
  Point<3> normal;
  double dLength;  

  curvature = 3-1; // dim-1
  dLength = sqrt(p(0)*p(0)+p(1)*p(1)+p(2)*p(2));
  
  normal[0] = p(0)/dLength;
  normal[1] = p(1)/dLength;
  normal[2] = p(2)/dLength;
  
  return (-trace(hessian) +
	  (hessian * normal) * normal +
	  (gradient * normal) * curvature);
}






template <int dim>
class LaplaceBeltrami 
{
  public:
    LaplaceBeltrami ();
				     // arguments are:
				     // triangulation
				     // right-hand side
				     // fe_degree for solution
				     // fe_degree for mapping
				     // exact solution is known
  
    ~LaplaceBeltrami ();
    void run ();
  
  private:

    void make_mesh ();
    void setup_system ();
    void assemble_system ();
    void solve ();
    void output_results () const;
    void compute_error () const;
  
  
  
    Triangulation<dim-1,dim>   triangulation;
    FE_Q<dim-1,dim>            fe;
    DoFHandler<dim-1,dim>      dof_handler;
    MappingQ<dim-1, dim>       mapping;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
  
    Vector<double>       solution;
    Vector<double>       system_rhs;
};


template <int dim>
LaplaceBeltrami<dim>::LaplaceBeltrami ()
		:
		fe (2),
		dof_handler(triangulation),
		mapping(2)
{}


template <int dim>
LaplaceBeltrami<dim>::~LaplaceBeltrami () 
{
  dof_handler.clear ();
}


template <int dim>
void LaplaceBeltrami<dim>::make_mesh ()
{
  HyperBallBoundary<dim> boundary_description;
  Triangulation<dim> volume_mesh;
  GridGenerator::half_hyper_ball(volume_mesh);
  
  volume_mesh.set_boundary (1, boundary_description);
  volume_mesh.set_boundary (0, boundary_description);
  volume_mesh.refine_global (1);
  
  static HyperBallBoundary<dim-1,dim> surface_description;
  triangulation.set_boundary (1, surface_description);
  triangulation.set_boundary (0, surface_description);
  
  std::set<unsigned char> boundary_ids;
  boundary_ids.insert(0);
  
  GridTools::extract_boundary_mesh (volume_mesh, triangulation,
				    boundary_ids);
  triangulation.refine_global(3);

  std::cout << "Surface mesh has " << triangulation.n_active_cells()
	    << " cells."
	    << std::endl;
}


template <int dim>
void LaplaceBeltrami<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  std::cout << "Surface mesh has " << dof_handler.n_dofs()
	    << " degrees of freedom."
	    << std::endl;
  
  CompressedSparsityPattern csp (dof_handler.n_dofs(),
                                 dof_handler.n_dofs());

  DoFTools::make_sparsity_pattern (dof_handler, csp);

  sparsity_pattern.copy_from (csp);

  system_matrix.reinit (sparsity_pattern);
  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}

template <int dim>
void LaplaceBeltrami<dim>::assemble_system () 
{  

  system_matrix = 0;
  system_rhs = 0;
  
  QGauss<dim-1>  quadrature_formula(2);

  FEValues<dim-1,dim> fe_values (mapping,fe, quadrature_formula, 
				 update_values   | update_cell_normal_vectors |
				 update_gradients |
				 update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector< double > rhs_values(n_q_points);
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  const RightHandSide<dim> rhs;
  
  typename DoFHandler<dim-1,dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit (cell);

      rhs.value_list (fe_values.get_quadrature_points(), rhs_values); 

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{	  
	  for (unsigned int j=0; j<dofs_per_cell; ++j) 
	    {	      
	      for (unsigned int q_point=0; q_point<n_q_points;
		   ++q_point)
		{
		  cell_matrix(i,j) 
		    += fe_values.shape_grad(i,q_point) *
		    fe_values.shape_grad(j,q_point) *
		    fe_values.JxW(q_point);
		}
	    }
	}

      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  
	  for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	    cell_rhs(i) += fe_values.shape_value(i,q_point) *
			   rhs_values[q_point]*
			   fe_values.JxW(q_point);
	}

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


  std::map<unsigned int,double> bdy_values; 
  VectorTools::interpolate_boundary_values (mapping,dof_handler,
					    0,
					    Solution<dim>(),
					    bdy_values
  );

  MatrixTools::apply_boundary_values (bdy_values,
				      system_matrix,
				      solution,
				      system_rhs,false);
}


template <int dim>
void LaplaceBeltrami<dim>::solve () 
{
  SolverControl           solver_control (solution.size(), 1e-7);
  SolverCG<>              cg (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);
}



template <int dim>
void LaplaceBeltrami<dim>::output_results () const
{
  
  std::string filename = "solution.vtk";

  std::ofstream output (filename.c_str());

  DataOut<dim-1,DoFHandler<dim-1,dim> > data_out;
  data_out.attach_dof_handler (dof_handler);

 

  data_out.add_data_vector (solution, "solution",DataOut_DoFData<DoFHandler<dim-1,dim>,dim-1,dim>::type_dof_data);
  data_out.build_patches (mapping,
			  mapping.get_degree());
  data_out.write_vtk (output);
}



template <int dim>
void LaplaceBeltrami<dim>::compute_error () const
{  
  Vector<float> difference_per_cell (triangulation.n_active_cells());
  VectorTools::integrate_difference (mapping, dof_handler, solution,
				     Solution<dim>(),
				     difference_per_cell,
				     QGauss<dim-1>(2*fe.degree+1),
				     VectorTools::H1_norm);
  
  
  std::cout << "H1 error = "
	    << difference_per_cell.l2_norm()
	    << std::endl;
}



template <int dim>
void LaplaceBeltrami<dim>::run () 
{
  make_mesh ();
  setup_system ();
  assemble_system ();
  solve ();
  output_results ();
  compute_error ();
}



int main ( int argc, char **argv )
{
  try
    {
      deallog.depth_console (0);

      LaplaceBeltrami<3> laplace_beltrami_2d;
      laplace_beltrami_2d.run();
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
