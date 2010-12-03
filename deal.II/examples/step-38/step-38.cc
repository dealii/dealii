#include <base/polynomial_space.h>
#include <base/parsed_function.h>
#include <base/smartpointer.h>
#include <base/convergence_table.h>
#include <base/quadrature_lib.h>
#include <base/quadrature_selector.h>
#include <base/utilities.h>

#include <lac/full_matrix.h>
#include <lac/vector.h>
#include <lac/solver_control.h>
#include <lac/solver_gmres.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/constraint_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/compressed_sparsity_pattern.h>

#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_out.h>
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

#include <sys/types.h>
#include <sys/stat.h>

#define deal_II_dimension 3

using std::cout;
using std::endl;
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
 
template <int dim>
double Solution<dim>::value (const Point<dim> &p,
			  const unsigned int) const 
{
  return sin(numbers::PI * p(0))*cos(numbers::PI * p(1))*exp(p(2)); 
}

template <int dim>
Tensor<1,dim> Solution<dim>::gradient (const Point<dim>   &p,
				       const unsigned int) const
{
  double dPi = numbers::PI;

  Tensor<1,dim> return_value;

  return_value[0] = dPi *cos(dPi * p(0))*cos(dPi * p(1))*exp(p(2));
  return_value[1] = -dPi *sin(dPi * p(0))*sin(dPi * p(1))*exp(p(2));
  return_value[2] = sin(dPi * p(0))*cos(dPi * p(1))*exp(p(2));
  
  // tangential gradient: nabla u - (nabla u nu)nu
  Point<dim> normal;
  double dLength;  

  dLength = sqrt(p(0)*p(0)+p(1)*p(1)+p(2)*p(2));
  
  normal[0] = p(0)/dLength;
  normal[1] = p(1)/dLength;
  normal[2] = p(2)/dLength;
  
  return return_value - (return_value*normal)*normal;
}

template <int dim>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide () : Function<dim>() {}
  
  virtual double value (const Point<dim>   &p,
			const unsigned int  component = 0) const;
};

template <int dim>
double RightHandSide<dim>::value (const Point<dim> &p,
				  const unsigned int comp) const 
{
  Assert(dim == 3,  ExcInternalError());
  
  double dPi = numbers::PI;

  // LB: u = Delta u - nu D2 u nu - (Grad u nu ) div (nu)

  Tensor<2,dim> hessian;

  hessian[0][0] = -dPi*dPi*sin(dPi*p(0))*cos(dPi*p(1))*exp(p(2));
  hessian[1][1] = -dPi*dPi*sin(dPi*p(0))*cos(dPi*p(1))*exp(p(2));
  hessian[2][2] = sin(dPi*p(0))*cos(dPi*p(1))*exp(p(2));

  hessian[0][1] = -dPi*dPi*cos(dPi*p(0))*sin(dPi*p(1))*exp(p(2));
  hessian[1][0] = -dPi*dPi*cos(dPi*p(0))*sin(dPi*p(1))*exp(p(2));

  hessian[0][2] = dPi*cos(dPi*p(0))*cos(dPi*p(1))*exp(p(2));
  hessian[2][0] = dPi*cos(dPi*p(0))*cos(dPi*p(1))*exp(p(2));

  hessian[1][2] = -dPi*sin(dPi*p(0))*sin(dPi*p(1))*exp(p(2));
  hessian[2][1] = -dPi*sin(dPi*p(0))*sin(dPi*p(1))*exp(p(2));

  Tensor<1,dim> gradient;
  gradient[0] = dPi * cos(dPi*p(0))*cos(dPi*p(1))*exp(p(2));
  gradient[1] = - dPi * sin(dPi*p(0))*sin(dPi*p(1))*exp(p(2));
  gradient[2] = sin(dPi*p(0))*cos(dPi*p(1))*exp(p(2));

  double curvature;
  Point<dim> normal;
  double dLength;  

  curvature = dim-1;
  dLength = sqrt(p(0)*p(0)+p(1)*p(1)+p(2)*p(2));
  
  normal[0] = p(0)/dLength;
  normal[1] = p(1)/dLength;
  normal[2] = p(2)/dLength;
  
  return -trace(hessian) + (hessian * normal) * normal + (gradient * normal)*curvature;
}






template <int dim>
class LaplaceBeltrami 
{
  public:
  LaplaceBeltrami (Triangulation<dim-1,dim> *tria, Function<dim> &func_data, 
		   unsigned int fe_degree = 1, unsigned int mapping_degree = 1,
		   Function<dim> *pExact = 0);
  // arguments are:
  // triangulation
  // right-hand side
  // fe_degree for solution
  // fe_degree for mapping
  // exact solution is known
  
  ~LaplaceBeltrami ();
  void run ();
  double compute_error (VectorTools::NormType norm_type) const;
  
  private:
  
  void setup_system ();
  void assemble_system ();
  void solve ();
  void output_results () const;
  
  
  
  Triangulation<dim-1,dim>   *pTria;
  FE_Q<dim-1,dim>            fe;
  DoFHandler<dim-1,dim>      dh;
  MappingQ<dim-1, dim>       mapping;

  ConstraintMatrix     matrix_constraints;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;
  
  Vector<double>       solution;
  Vector<double>       system_rhs;

  Function<dim> &rhs_func; // function data

  Function<dim> *pExact; // exact solution if provided

};


template <int dim>
LaplaceBeltrami<dim>::LaplaceBeltrami (Triangulation<dim-1,dim> *tria,
				       Function<dim> &func_data,
				       unsigned int fe_degree,
				       unsigned int mapping_degree,
				       Function<dim> *pExact
				       )
		:
  fe (fe_degree),
  dh(*tria),
  mapping(mapping_degree),
  rhs_func(func_data)
{
  pTria = tria;
  this->pExact = pExact;

}

template <int dim>
LaplaceBeltrami<dim>::~LaplaceBeltrami () 
{
  dh.clear ();
}

template <int dim>
void LaplaceBeltrami<dim>::setup_system ()
{
  dh.distribute_dofs (fe);

  matrix_constraints.clear ();
  matrix_constraints.close ();

  CompressedSparsityPattern csp (dh.n_dofs(),
                                 dh.n_dofs());

  DoFTools::make_sparsity_pattern (dh, csp);
  matrix_constraints.condense (csp);

  sparsity_pattern.copy_from (csp);

  system_matrix.reinit (sparsity_pattern);
  solution.reinit (dh.n_dofs());
  system_rhs.reinit (dh.n_dofs());
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

  typename DoFHandler<dim-1,dim>::active_cell_iterator
    cell = dh.begin_active(),
    endc = dh.end();

  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit (cell);

      rhs_func.value_list (fe_values.get_quadrature_points(), rhs_values); 

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
   VectorTools::interpolate_boundary_values (mapping,dh,
					     0,
					     *pExact,
					     bdy_values
					     );

   MatrixTools::apply_boundary_values (bdy_values,
				       system_matrix,
				       solution,
				       system_rhs,false);

  // condense matrices
  matrix_constraints.condense (system_matrix);
  matrix_constraints.condense (system_rhs);
}


template <int dim>
void LaplaceBeltrami<dim>::solve () 
{
  SolverControl           solver_control (1000, 1e-7);
  SolverCG<>              cg (solver_control);

  cg.solve (system_matrix, solution, system_rhs,
	    PreconditionIdentity());


  matrix_constraints.distribute (solution);

}



template <int dim>
void LaplaceBeltrami<dim>::output_results () const
{
  
  std::string filename = "solution.vtk";

  std::ofstream output (filename.c_str());

  DataOut<dim-1,DoFHandler<dim-1,dim> > data_out;
  data_out.attach_dof_handler (dh);

 

  data_out.add_data_vector (solution, "solution",DataOut_DoFData<DoFHandler<dim-1,dim>,dim-1,dim>::type_dof_data);
  data_out.build_patches (mapping,
			  mapping.get_degree());
  data_out.write_vtk (output);
}



template <int dim>
void LaplaceBeltrami<dim>::run () 
{

  setup_system();
  assemble_system();
  solve();
  output_results();

}

//################################################################################//

template <int dim>
double LaplaceBeltrami<dim>::compute_error (VectorTools::NormType norm_type) const
{  
  Assert(pExact != 0, ExcInternalError());
  
  Vector<float> difference_per_cell (pTria->n_active_cells());
  VectorTools::integrate_difference (mapping, dh, solution,
				     *pExact, difference_per_cell,
				     QGauss<(dim-1)>(2*fe.degree+1),
				     norm_type);

  return difference_per_cell.l2_norm();
}




int main ( int argc, char **argv )
{
  std::cout<<std::endl;
  std::cout<<"================================"<<std::endl;
  std::cout<<"          LB PROBLEM "<<std::endl;
  std::cout<<"================================"<<std::endl;
  std::cout<<std::endl;
  
  
  Triangulation<deal_II_dimension-1,deal_II_dimension> tria;
  
  // create a mesh consisting on the boundary of the half sphere... thx Seba
  std::map< Triangulation<deal_II_dimension-1,deal_II_dimension>::cell_iterator,
    Triangulation<deal_II_dimension,deal_II_dimension>::face_iterator> surface_to_volume_mapping;

  HyperBallBoundary<deal_II_dimension> boundary_description;
  Triangulation<deal_II_dimension> volume_mesh;
  GridGenerator::half_hyper_ball(volume_mesh);
  
  volume_mesh.set_boundary (1, boundary_description);
  volume_mesh.set_boundary (0, boundary_description);
  volume_mesh.refine_global (1);
  
  static HyperBallBoundary<deal_II_dimension-1,deal_II_dimension> surface_description;
  tria.set_boundary (1, surface_description);
  tria.set_boundary (0, surface_description);
  
  std::set<unsigned char> boundary_ids;
  boundary_ids.insert(0);
  
  GridTools::extract_boundary_mesh (volume_mesh, tria,
				    surface_to_volume_mapping,
				    boundary_ids);
 

  tria.refine_global(4);
      

  RightHandSide<deal_II_dimension> rhs;
  Solution<deal_II_dimension> exact;
  
  LaplaceBeltrami<deal_II_dimension> laplace_beltrami_2d(&tria,rhs,2,2,&exact);
  
  laplace_beltrami_2d.run();
  
  std::cout<<laplace_beltrami_2d.compute_error(VectorTools::H1_norm)<<std::endl;

  tria.set_boundary(0);
  return 0;
}











