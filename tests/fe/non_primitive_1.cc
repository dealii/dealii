//----------------------------  non_primitive_1.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  non_primitive_1.cc  ---------------------------

// assemble a matrix for the stokes equation in different ways for
// primitive elements. the name of the program is only to indicate
// that we use the way to assemble matrices necessary for
// non-primitive elements, and compare the resulting matrix with what
// we get if we use the usual way which we can take if the FE is
// primitive. then use a third way where we use some optimizations,
// and compare again whether the matrices are the same


#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <lac/sparse_matrix.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>

#include <fstream>


				 // create the matrix in the plain old
				 // way, which is valid if the finite
				 // element shape functions are
				 // primitive
template <int dim>
void
create_stokes_matrix_1 (const DoFHandler<dim> &dof_handler,
			SparseMatrix<double>  &A)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  FullMatrix<double> local_matrix (dofs_per_cell, dofs_per_cell);

  QGauss<dim> quadrature (3);
  const unsigned int n_q_points = quadrature.n_quadrature_points;
  
  FEValues<dim> fe_values (fe, quadrature,
			   update_values | update_gradients |
			   update_JxW_values);

  const double nu = 3.14159265358e-2;
  
  for (; cell!=endc; ++cell)
    {
      local_matrix.clear ();
      fe_values.reinit (cell);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  for (unsigned int q=0; q<n_q_points; ++q)
	    {
	      const unsigned int
		comp_i = fe.system_to_component_index(i).first,
		comp_j = fe.system_to_component_index(j).first;

					       // velocity-velocity coupling?
	      if ((comp_i<dim) && (comp_j<dim))
		if (comp_i == comp_j)
		  local_matrix(i,j) += (nu *
				       (fe_values.shape_grad(i,q) * 
					fe_values.shape_grad(j,q)  ) *
				       fe_values.JxW(q)); 

					       // velocity-pressure coupling
	      if ((comp_i<dim) && (comp_j==dim))
		local_matrix(i,j) += (-fe_values.shape_grad(i,q)[comp_i] * 
				     fe_values.shape_value(j,q)        *
				     fe_values.JxW(q)); 

					       // pressure-velocity coupling
	      if ((comp_i==dim) && (comp_j<dim))
		local_matrix(i,j) += (fe_values.shape_value(i,q) * 
				     fe_values.shape_grad(j,q)[comp_j] *
				     fe_values.JxW(q)); 
	    };
      
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  A.add (local_dof_indices[i],
		 local_dof_indices[j],
		 local_matrix(i,j));
    };
};



				 // create the matrix in the simple
				 // way that is necessary when you
				 // want to use non-primitive shape
				 // functions
template <int dim>
void
create_stokes_matrix_2 (const DoFHandler<dim> &dof_handler,
			SparseMatrix<double>  &A)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  FullMatrix<double> local_matrix (dofs_per_cell, dofs_per_cell);

  QGauss<dim> quadrature (3);
  const unsigned int n_q_points = quadrature.n_quadrature_points;
  
  FEValues<dim> fe_values (fe, quadrature,
			   update_values | update_gradients |
			   update_JxW_values);

  const double nu = 3.14159265358e-2;
  
  for (; cell!=endc; ++cell)
    {
      local_matrix.clear ();
      fe_values.reinit (cell);
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	for (unsigned int comp_i=0; comp_i<fe.n_components(); ++comp_i)
	  for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
	    for (unsigned int comp_j=0; comp_j<fe.n_components(); ++comp_j)
	      for (unsigned int q=0; q<n_q_points; ++q) 
		{
						   // velocity-velocity coupling?
		  if ((comp_i<dim) && (comp_j<dim))
		    if (comp_i == comp_j)
		      local_matrix(i,j)
			+= (nu *
			    (fe_values.shape_grad_component(i,q,comp_i) * 
			     fe_values.shape_grad_component(j,q,comp_j)  ) *
			    fe_values.JxW(q)); 

						   // velocity-pressure coupling
		  if ((comp_i<dim) && (comp_j==dim))
		    local_matrix(i,j)
		      += (-fe_values.shape_grad_component(i,q,comp_i)[comp_i] * 
			  fe_values.shape_value_component(j,q,comp_j)        *
			  fe_values.JxW(q)); 

						   // pressure-velocity coupling
		  if ((comp_i==dim) && (comp_j<dim))
		    local_matrix(i,j)
		      += (fe_values.shape_value_component(i,q,comp_i) * 
			  fe_values.shape_grad_component(j,q,comp_j)[comp_j] *
			  fe_values.JxW(q)); 
		};

      
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  A.add (local_dof_indices[i],
		 local_dof_indices[j],
		 local_matrix(i,j));
    };
};



				 // create the matrix in a way that is
				 // necessary when you want to use
				 // non-primitive shape
				 // functions. compared to the second
				 // possibility used above, use some
				 // optimizations
template <int dim>
void
create_stokes_matrix_3 (const DoFHandler<dim> &dof_handler,
			SparseMatrix<double>  &A)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  FullMatrix<double> local_matrix (dofs_per_cell, dofs_per_cell);

  QGauss<dim> quadrature (3);
  const unsigned int n_q_points = quadrature.n_quadrature_points;
  
  FEValues<dim> fe_values (fe, quadrature,
			   update_values | update_gradients |
			   update_JxW_values);

  const double nu = 3.14159265358e-2;
  
  for (; cell!=endc; ++cell)
    {
      local_matrix.clear ();
      fe_values.reinit (cell);
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	for (unsigned int comp_i=0; comp_i<fe.n_components(); ++comp_i)
	  if (fe.get_nonzero_components(i)[comp_i] == true)
	    for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
	      for (unsigned int comp_j=0; comp_j<fe.n_components(); ++comp_j)
		if (fe.get_nonzero_components(j)[comp_j] == true)		
		  for (unsigned int q=0; q<n_q_points; ++q) 
		    {
						       // velocity-velocity coupling?
		      if ((comp_i<dim) && (comp_j<dim))
			if (comp_i == comp_j)
			  local_matrix(i,j)
			    += (nu *
				(fe_values.shape_grad_component(i,q,comp_i) * 
				 fe_values.shape_grad_component(j,q,comp_j)  ) *
				fe_values.JxW(q)); 
		      
						       // velocity-pressure coupling
		      if ((comp_i<dim) && (comp_j==dim))
			local_matrix(i,j)
			  += (-fe_values.shape_grad_component(i,q,comp_i)[comp_i] * 
			      fe_values.shape_value_component(j,q,comp_j)        *
			      fe_values.JxW(q)); 
		      
						       // pressure-velocity coupling
		      if ((comp_i==dim) && (comp_j<dim))
			local_matrix(i,j)
			  += (fe_values.shape_value_component(i,q,comp_i) * 
			      fe_values.shape_grad_component(j,q,comp_j)[comp_j] *
			      fe_values.JxW(q)); 
		    };
      
      
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  A.add (local_dof_indices[i],
		 local_dof_indices[j],
		 local_matrix(i,j));
    };
};




template <int dim>
void
test ()
{
  Triangulation<dim> triangulation;
  if (dim != 2)
    GridGenerator::hyper_cube (triangulation);
  else
    GridGenerator::hyper_cube_slit (triangulation);    
  triangulation.refine_global (1);
  triangulation.begin_active()->set_refine_flag();
  triangulation.execute_coarsening_and_refinement();
  triangulation.last_active()->set_refine_flag ();
  triangulation.execute_coarsening_and_refinement();

  if (dim < 3)
    triangulation.refine_global (1);

  deallog << "dim=" << dim
	  << ", n_cells=" << triangulation.n_active_cells()
	  << std::endl;
  
  FESystem<dim> fe (FE_Q<dim>(2), dim,
		    FE_Q<dim>(1), 1);
  DoFHandler<dim> dof_handler (triangulation);
  dof_handler.distribute_dofs (fe);
  
  SparsityPattern sparsity(dof_handler.n_dofs(),
			   dof_handler.n_dofs());
  
				   // have a mask that indicates that
				   // for the stokes equation the
				   // pressure does not couple to
				   // itself
  std::vector<std::vector<bool> > mask (dim+1, std::vector<bool> (dim+1, true));
  mask[dim][dim] = false;
  
  DoFTools::make_sparsity_pattern (dof_handler, mask, sparsity);
  sparsity.compress ();
  
  SparseMatrix<double> A1 (sparsity);
  SparseMatrix<double> A2 (sparsity);
  SparseMatrix<double> A3 (sparsity);

  create_stokes_matrix_1 (dof_handler, A1);
  create_stokes_matrix_2 (dof_handler, A2);
  create_stokes_matrix_3 (dof_handler, A3);

				   // write out the contents of the
				   // matrix and compare for equality
				   // with the other matrices. to
				   // reduce the amount of data
				   // written out a little bit, only
				   // write every so-many-th element
  for (unsigned int i=0; i<A1.n_nonzero_elements(); ++i)
    {
      if (i % (dim*dim*dim) == 0)
	deallog << i << ' ' << A1.global_entry(i) << std::endl;
      Assert (A1.global_entry(i) == A2.global_entry(i),
	      ExcInternalError());
      Assert (A1.global_entry(i) == A3.global_entry(i),
	      ExcInternalError());
    };
};



int main ()
{
  std::ofstream logfile("non_primitive_1.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<1> ();
  test<2> ();
  test<3> ();
};



	 
