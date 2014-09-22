// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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


// test by Denis Davydov, see bug report 76: can't create hanging node
// constraints for a particular combination of finite elements

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>

#include <deal.II/base/quadrature_lib.h>   //needed for assembling the matrix using quadrature on each cell
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_direct.h>     //direct solver
#include <deal.II/lac/sparse_matrix.h>     //definition of the sparse matrix
#include <deal.II/lac/compressed_sparsity_pattern.h> //for the intermediate sparsity pattern structure
#include <deal.II/lac/solver_cg.h>         //CG solver
#include <deal.II/lac/precondition.h>      //and a preconditioner
#include <deal.II/lac/constraint_matrix.h> //conform hanging nodes DoF to certain constrains to make the solution continuous
#include <deal.II/lac/iterative_inverse.h> //need for Schur matrix
#include <deal.II/grid/tria.h>             //triangulation class
#include <deal.II/grid/grid_generator.h>   //standard functions to generate grid
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>    //these two to loop over cells/faces
#include <deal.II/grid/tria_iterator.h>    // ^
#include <deal.II/grid/tria_boundary_lib.h>//to use not straightline boundaries
#include <deal.II/dofs/dof_handler.h>      //associate DoF to cells/vertices/lines
#include <deal.II/dofs/dof_accessor.h>     //provides information about the degrees of freedom local to a cell
#include <deal.II/dofs/dof_tools.h>        //needed for the creation of sparsity patterns of sparse matrices
#include <deal.II/dofs/dof_renumbering.h>  //use to renumber DoF to have a better sparsity pattern
#include <deal.II/fe/fe_values.h>          //needed for assembling the matrix using quadrature on each cell
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/fe_field_function.h> //for evaluating solution at a given point
#include <deal.II/fe/fe_system.h>          /*for vector-valued FE; vector-valude FE will be composed from the regular Q1 included below*/
#include <deal.II/fe/fe_q.h>               /*Q1 FE: 1 DoF in each vertex, none on sides...)*/
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/q_collection.h>
#include <fstream>
#include <iostream>
#include <cmath>                           //sqrt and fabs functions
#include <memory> //smart pointers
#include <algorithm> //min,max


template <int dim>
class ElasticProblem
{
  public:
  ElasticProblem();
  ~ElasticProblem ();
  void run();
  private:
  enum {
	u_block = 0,
	lambda_block = 1
  };
    
  void make_grid();
  void setup_system ();

  const unsigned int n_blocks;//total number of blocks
  const unsigned int n_components;// total number of components
  const unsigned int first_u_comp;//where displacements starts
  const unsigned int first_lamda_comp;//where lagrage multipliers start
  const unsigned int degree;//of shape functions

  std::vector<types::global_dof_index> dofs_per_block;

  Triangulation<dim>   triangulation; //a triangulation object of the "dim"-dimensional domain;
  hp::DoFHandler<dim>      dof_handler;   //is assotiated with triangulation;

  FESystem<dim>         elasticity_fe;
  FESystem<dim>         elasticity_w_lagrange_fe;
  hp::FECollection<dim>        fe; //used to stack several other elements together to form one vector-valued finite element.

  ConstraintMatrix     hanging_node_constraints; //object to hold hanging node constraints after refinement
  ConstraintMatrix     constraints;

  BlockSparsityPattern      sparsity_pattern; //store sparsity pattern

  SymmetricTensor<4,dim> stress_strain_tensor; // elastic constitutive equations

  hp::QCollection<dim>  quadrature_formula;//a Gauss quadrature formula to be used for integral evaluation

  const FEValuesExtractors::Vector u_fe;
  const FEValuesExtractors::Vector lamda_fe;

  unsigned int id_of_lagrange_mult;
  double beta1;

  SparsityPattern      sparsity_pattern_nb;//needed when switching between block system and non-block.
};

///
template <int dim>
SymmetricTensor<4,dim>
get_stress_strain_tensor (const double lambda, const double mu)
{
  SymmetricTensor<4,dim> tmp;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	for (unsigned int l=0; l<dim; ++l)
	  tmp[i][j][k][l] = (
	    ((i==k) && (j==l) ? mu : 0.0) +
	    ((i==l) && (j==k) ? mu : 0.0) +
	    ((i==j) && (k==l) ? lambda : 0.0)
	  );
  return tmp;
}

///
/*vector-valued BC function */
template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  BoundaryValues () : Function<dim>(dim) {}

  virtual void vector_value (const Point<dim> &p,
			     Vector<double>   &values) const;
  virtual double value (const Point<dim> &p,
			unsigned int compoment) const;
};

template <int dim>
inline
void BoundaryValues<dim>::vector_value (const Point<dim> &p,
					Vector<double>   &values) const
{
  values(0) = -0.001;
}

template <int dim>
inline
double BoundaryValues<dim>::value (const Point<dim> &p,
				   unsigned int compoment) const
{
  return -0.001;
}

///
/* vector valued constrained displacements*/
/*vector-valued RHS function */
template <int dim>
class ConstrainValues :  public Function<dim>
{
public:
  ConstrainValues (); //need to pass down to the base class of how many components the function consists; default - 1
  virtual void vector_value (const Point<dim> &p,
			     Vector<double>   &values) const;//returns calculated values in the second argument;
  virtual void vector_value_list (const std::vector<Point<dim> > &points,
				  std::vector<Vector<double> >   &value_list) const;//values at several points at once
  //prevent from calling virtual function "vector_value" to frequently
};

template <int dim>
ConstrainValues<dim>::ConstrainValues ()
:
Function<dim> (dim) /*pass down to the main class the number of components of the function*/
{}

template <int dim>
inline
void ConstrainValues<dim>::vector_value (const Point<dim> &p,
					 Vector<double>   &values) const
{
  Assert (values.size() == dim,
	  ExcDimensionMismatch (values.size(), dim));//check is the size of "values" is correct
  //Assert (dim >= 2, ExcNotImplemented());//not implemented for 1d
  values[0] = 0.0;
}

template <int dim>
void ConstrainValues<dim>::vector_value_list (const std::vector<Point<dim> > &points,
					      std::vector<Vector<double> >   &value_list) const
{
  Assert (value_list.size() == points.size(),
	  ExcDimensionMismatch (value_list.size(), points.size()));//check if input-output is consistent

  const unsigned int n_points = points.size();//number of points at which the function is to be evaluated

  for (unsigned int p=0; p<n_points; ++p)
    ConstrainValues<dim>::vector_value (points[p],
                                        value_list[p]);
}
  

/* constructor*/
template <int dim>
ElasticProblem<dim>::ElasticProblem ()
:
n_blocks(2),
n_components(dim*n_blocks),
first_u_comp(0),
first_lamda_comp(dim),
degree(1),
dofs_per_block(n_blocks),
dof_handler (triangulation),  /*assotiate dof_handler to the triangulation */
elasticity_fe (FE_Q<dim>(degree), dim, //use dim FE_Q of a given degree to represent displacements
	       FE_Nothing<dim>(), dim  // zero extension of lagrange multipliers elsewhere in the domain
),
elasticity_w_lagrange_fe(
  FE_Q<dim>(degree), dim,
  FE_Q<dim>(degree), dim  // same for lagrange multipliers
),
u_fe(first_u_comp),
lamda_fe(first_lamda_comp),
id_of_lagrange_mult(1), 
beta1(1.0)            
{
  fe.push_back (elasticity_fe);//FE index 0
  fe.push_back (elasticity_w_lagrange_fe);//FE index 1

  quadrature_formula.push_back(QGauss<dim>(degree+2));
  quadrature_formula.push_back(QGauss<dim>(degree+2));

  stress_strain_tensor
    = get_stress_strain_tensor<dim> (/*lambda = */ 9.695e10,
				     /*mu     = */ 7.617e10);
}

template <int dim>
ElasticProblem<dim>::~ElasticProblem ()
{
  dof_handler.clear (); 
}


template <int dim>
void ElasticProblem<dim>::make_grid ()
{

  Triangulation<dim>   triangulationL;
  Triangulation<dim>   triangulationR;
  GridGenerator::hyper_cube (triangulationL, -1,0); //create a square [-1,1]^d domain
  GridGenerator::hyper_cube (triangulationR, -1,0); //create a square [-1,1]^d domain
  Point<dim> shift_vector;
  shift_vector[0] = 1.0;
  GridTools::shift(shift_vector,triangulationR);

  const unsigned int n_faces_per_cell = GeometryInfo< dim >::faces_per_cell;
    
  triangulationL.begin_active()->set_material_id(0);
  triangulationR.begin_active()->set_material_id(id_of_lagrange_mult);

  for (unsigned int i=0;i<n_faces_per_cell;i++) {
    triangulationL.begin_active()->face(i)->set_boundary_indicator(i);
    triangulationR.begin_active()->face(i)->set_boundary_indicator(n_faces_per_cell+i);
  }
    
  GridGenerator::merge_triangulations (triangulationL, triangulationR, triangulation);
 
  //triangulation.refine_global (2); //refine twice globally before solving
}


template <int dim>
void ElasticProblem<dim>::setup_system ()
{
  std::vector<unsigned int> block_component(n_components, u_block);// init to represent u everywhere
  for (unsigned int i = 0; i < dim; i++)
    block_component[i+dim] = lambda_block;

  //(1) set active FE indices based in material id...
  typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
						     endc = dof_handler.end();
  unsigned int n_lagrange_cells = 0;
  unsigned int n_elasticity_cells = 0;
  for (; cell!=endc; ++cell) //loop over all cells
    {
      if (cell->material_id() == id_of_lagrange_mult ) {
	cell->set_active_fe_index(1); // for lagrangian cells..
	n_lagrange_cells++;
      } else {
	cell->set_active_fe_index(0);
	n_elasticity_cells++;
      }
    }
  deallog<<" number of cells (L/E): "<< n_lagrange_cells << "; "<<n_elasticity_cells <<std::endl;
  Assert (n_lagrange_cells > 0, ExcInternalError()); //there should be at least 1 cell! Otherwise DoFHanlder crashes with 0 dofs for block 2!
  //
  //(2) distribute DoFs
  dof_handler.distribute_dofs (fe); /*enumerate DoF based on fe object (which knows about shape-functions used, their order and dimension*/
  DoFRenumbering::Cuthill_McKee (dof_handler);//do renumberring; must be done right after distributing DoF !!!
  DoFRenumbering::component_wise(dof_handler,block_component);
  DoFTools::count_dofs_per_block(dof_handler,dofs_per_block,block_component);

  deallog<<"dofs per block:  U="<<dofs_per_block[u_block]<<" L="<<dofs_per_block[lambda_block]<<std::endl;

  //related to hanging nodes and refinement:
  //(3) create hanging nodes constrains (also add normal diricle BC here?
  hanging_node_constraints.clear ();//clear all content that may left from previous calculations
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);//fill the object representing constraints
  hanging_node_constraints.close ();//close the object - sort and rearrange constraints for effective calculations

  //(4) create block sparsity pattern (define which elements in sparse matrix are non-zero; prescribe coupling between blocks)
  // following step-22 use of simple compressed block sparsity pattern for efficiency
  BlockCompressedSparsityPattern compressed_sparsity_pattern(n_blocks,
							     n_blocks);

  compressed_sparsity_pattern.block(u_block     ,u_block     ).reinit(dofs_per_block[u_block]     ,dofs_per_block[u_block]);
  compressed_sparsity_pattern.block(u_block     ,lambda_block).reinit(dofs_per_block[u_block]     ,dofs_per_block[lambda_block]);
  compressed_sparsity_pattern.block(lambda_block,u_block     ).reinit(dofs_per_block[lambda_block],dofs_per_block[u_block]);
  compressed_sparsity_pattern.block(lambda_block,lambda_block).reinit(dofs_per_block[lambda_block],dofs_per_block[lambda_block]);

  compressed_sparsity_pattern.collect_sizes();

  Table<2, DoFTools::Coupling> coupling(n_components, n_components);
  for (unsigned int ii = 0; ii< n_components;++ii)
    for (unsigned int jj=0; jj< n_components; ++jj)
      {
	if ( (block_component[ii]==lambda_block) &&
	     (block_component[jj]==lambda_block)
	)
	  coupling[ii][jj] = DoFTools::none;//diagonal = 0
	else
	  coupling[ii][jj] = DoFTools::always;//full coupling (u,u), (u,lambda)
      }

  hanging_node_constraints.condense (compressed_sparsity_pattern);

  DoFTools::make_sparsity_pattern (dof_handler,
				   coupling,
				   compressed_sparsity_pattern,
				   hanging_node_constraints,
				   false);

  compressed_sparsity_pattern.print (deallog.get_file_stream());
  hanging_node_constraints.print (deallog.get_file_stream());
} 
  
  
  
template <int dim>
void ElasticProblem<dim>::run ()
{
  make_grid ();
  setup_system ();
}  
  

int main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  ElasticProblem <2> elastic_problem;
  elastic_problem.run();

  deallog << "OK" << std::endl;
}



