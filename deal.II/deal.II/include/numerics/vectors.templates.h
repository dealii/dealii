//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name:  $
//
//    Copyright (C) 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef _deal2__vectors_templates_h
#define _deal2__vectors_templates_h

#include <base/function.h>
#include <base/quadrature.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/petsc_vector.h>
#include <lac/petsc_block_vector.h>
#include <lac/trilinos_vector.h>
#include <lac/trilinos_block_vector.h>
#include <lac/sparse_matrix.h>
#include <lac/compressed_simple_sparsity_pattern.h>
#include <lac/precondition.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/filtered_matrix.h>
#include <lac/constraint_matrix.h>
#include <grid/tria_iterator.h>
#include <grid/grid_tools.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe.h>
#include <fe/fe_tools.h>
#include <hp/fe_values.h>
#include <fe/mapping_q1.h>
#include <hp/mapping_collection.h>
#include <hp/q_collection.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>

#include <numeric>
#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>

DEAL_II_NAMESPACE_OPEN

template <class VECTOR, class DH>
void VectorTools::interpolate (const Mapping<DH::dimension,DH::space_dimension>    &mapping,
			       const DH              &dof,
			       const Function<DH::space_dimension>   &function,
			       VECTOR                &vec)
{
  const unsigned int dim=DH::dimension;

  Assert (dof.get_fe().n_components() == function.n_components,
	  ExcDimensionMismatch(dof.get_fe().n_components(),
			       function.n_components));
  
  const hp::FECollection<DH::dimension,DH::space_dimension> fe (dof.get_fe());
  const unsigned int          n_components = fe.n_components();
  const bool                  fe_is_system = (n_components != 1);
  
  typename DH::active_cell_iterator cell = dof.begin_active(),
				    endc = dof.end();

				   // For FESystems many of the
				   // unit_support_points will appear
				   // multiply, as a point may be
				   // unit_support_point for several of the
				   // components of the system.  The following
				   // is rather complicated, but at least
				   // attempts to avoid evaluating the
				   // vectorfunction multiple times at the
				   // same point on a cell.
				   //
				   // note that we have to set up all of the
				   // following arrays for each of the
				   // elements in the FECollection (which
				   // means only once if this is for a regular
				   // DoFHandler)
  std::vector<std::vector<Point<dim> > > unit_support_points (fe.size());
  for (unsigned int fe_index=0; fe_index<fe.size(); ++fe_index)
    {
      unit_support_points[fe_index] = fe[fe_index].get_unit_support_points();
      Assert (unit_support_points[fe_index].size() != 0,
	      ExcNonInterpolatingFE());
    }
  

				   // Find the support points 
				   // on a cell that
				   // are multiply mentioned in 
				   // unit_support_points.
				   // Mark the first representative
				   // of each multiply mentioned
				   // support point by appending its
				   // dof index to dofs_of_rep_points.
				   // Each multiple point gets to know
				   // the dof index of its representative
				   // point by the dof_to_rep_dof_table.

				   // the following vector collects all dofs i,
				   // 0<=i<fe.dofs_per_cell, for that
				   // unit_support_points[i] 
				   // is a representative one. i.e.
				   // the following vector collects all rep dofs.
				   // the position of a rep dof within this vector
				   // is called rep index.
  std::vector<std::vector<unsigned int> > dofs_of_rep_points(fe.size());
				   // the following table converts a dof i
				   // to the rep index.
  std::vector<std::vector<unsigned int> > dof_to_rep_index_table(fe.size());

  std::vector<unsigned int> n_rep_points (fe.size(), 0);
  
  for (unsigned int fe_index=0; fe_index<fe.size(); ++fe_index)
    {
      for (unsigned int i=0; i<fe[fe_index].dofs_per_cell; ++i)
	{
	  bool representative=true;
					   // the following loop is looped
					   // the other way round to get
					   // the minimal effort of
					   // O(fe.dofs_per_cell) for multiple
					   // support points that are placed
					   // one after the other.
	  for (unsigned int j=dofs_of_rep_points[fe_index].size(); j>0; --j)
	    if (unit_support_points[fe_index][i] 
		== unit_support_points[fe_index][dofs_of_rep_points[fe_index][j-1]])
	      {
		dof_to_rep_index_table[fe_index].push_back(j-1);
		representative=false;
		break;
	      }
      
	  if (representative)
	    {
					       // rep_index=dofs_of_rep_points.size()
	      dof_to_rep_index_table[fe_index].push_back(dofs_of_rep_points[fe_index].size());
					       // dofs_of_rep_points[rep_index]=i
	      dofs_of_rep_points[fe_index].push_back(i);
	      ++n_rep_points[fe_index];
	    }
	}

      Assert(dofs_of_rep_points[fe_index].size()==n_rep_points[fe_index],
	     ExcInternalError());
      Assert(dof_to_rep_index_table[fe_index].size()==fe[fe_index].dofs_per_cell,
	     ExcInternalError());
    }
  
  const unsigned int max_rep_points = *std::max_element (n_rep_points.begin(),
							 n_rep_points.end());
  std::vector<unsigned int> dofs_on_cell (fe.max_dofs_per_cell());
  std::vector<Point<DH::space_dimension> >  rep_points (max_rep_points);
  
				   // get space for the values of the
				   // function at the rep support points.
				   //
				   // have two versions, one for system fe
				   // and one for scalar ones, to take the
				   // more efficient one respectively
  std::vector<std::vector<double> >         function_values_scalar(fe.size());
  std::vector<std::vector<Vector<double> > > function_values_system(fe.size());

				   // Make a quadrature rule from support points
				   // to feed it into FEValues
  hp::QCollection<dim> support_quadrature;
  for (unsigned int fe_index=0; fe_index<fe.size(); ++fe_index)
    support_quadrature.push_back (Quadrature<dim>(unit_support_points[fe_index]));

				   // Transformed support points are computed by
				   // FEValues
  hp::MappingCollection<dim,DH::space_dimension> mapping_collection (mapping);
  
  hp::FEValues<dim, DH::space_dimension> fe_values (mapping_collection,
			       fe, support_quadrature, update_quadrature_points);
  
  for (; cell!=endc; ++cell)
    {
      const unsigned int fe_index = cell->active_fe_index();
      
				       // for each cell:
				       // get location of finite element
				       // support_points
      fe_values.reinit(cell);
      const std::vector<Point<DH::space_dimension> >& support_points =
	fe_values.get_present_fe_values().get_quadrature_points();
      
				       // pick out the representative
				       // support points
      rep_points.resize (dofs_of_rep_points[fe_index].size());
      for (unsigned int j=0; j<dofs_of_rep_points[fe_index].size(); ++j)
	rep_points[j] = support_points[dofs_of_rep_points[fe_index][j]];

				       // get indices of the dofs on this cell
      dofs_on_cell.resize (fe[fe_index].dofs_per_cell);      
      cell->get_dof_indices (dofs_on_cell);


      if (fe_is_system)
	{
					   // get function values at
					   // these points. Here: get
					   // all components
	  function_values_system[fe_index]
	    .resize (n_rep_points[fe_index],
		     Vector<double> (fe[fe_index].n_components()));
	  function.vector_value_list (rep_points,
				      function_values_system[fe_index]);
					   // distribute the function
					   // values to the global
					   // vector
	  for (unsigned int i=0; i<fe[fe_index].dofs_per_cell; ++i)
	    {
	      const unsigned int component
		= fe[fe_index].system_to_component_index(i).first;
	      const unsigned int rep_dof=dof_to_rep_index_table[fe_index][i];
	      vec(dofs_on_cell[i])
		= function_values_system[fe_index][rep_dof](component);
	    }
	}
      else
	{
					   // get first component only,
					   // which is the only component
					   // in the function anyway
	  function_values_scalar[fe_index].resize (n_rep_points[fe_index]);
	  function.value_list (rep_points,
			       function_values_scalar[fe_index],
			       0);
					   // distribute the function
					   // values to the global
					   // vector
	  for (unsigned int i=0; i<fe[fe_index].dofs_per_cell; ++i)
	    vec(dofs_on_cell[i]) 
	      = function_values_scalar[fe_index][dof_to_rep_index_table[fe_index][i]];
	}
    }
}


template <class VECTOR, class DH>
void VectorTools::interpolate (const DH              &dof,
			       const Function<DH::space_dimension>   &function,
			       VECTOR                &vec)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  interpolate(StaticMappingQ1<DH::dimension, DH::space_dimension>::mapping, 
	      dof, function, vec);
}




template <int dim, class InVector, class OutVector, int spacedim>
void
VectorTools::interpolate (const DoFHandler<dim,spacedim>           &dof_1,
			  const DoFHandler<dim,spacedim>           &dof_2,
			  const FullMatrix<double>        &transfer,
			  const InVector                  &data_1,
			  OutVector                       &data_2)
{
  Vector<double> cell_data_1(dof_1.get_fe().dofs_per_cell);
  Vector<double> cell_data_2(dof_2.get_fe().dofs_per_cell);

  std::vector<short unsigned int> touch_count (dof_2.n_dofs(), 0);
  std::vector<unsigned int>       local_dof_indices (dof_2.get_fe().dofs_per_cell);
  
  typename DoFHandler<dim,spacedim>::active_cell_iterator h = dof_1.begin_active();
  typename DoFHandler<dim,spacedim>::active_cell_iterator l = dof_2.begin_active();
  const typename DoFHandler<dim,spacedim>::cell_iterator endh = dof_1.end();
  
  for(; h != endh; ++h, ++l)
  {
    h->get_dof_values(data_1, cell_data_1);
    transfer.vmult(cell_data_2, cell_data_1);

    l->get_dof_indices (local_dof_indices);
  
				   // distribute cell vector
    for (unsigned int j=0; j<dof_2.get_fe().dofs_per_cell; ++j) 
      {
	data_2(local_dof_indices[j]) += cell_data_2(j);

					 // count, how often we have
					 // added to this dof
	Assert (touch_count[local_dof_indices[j]] < 255,
		ExcInternalError());	
	++touch_count[local_dof_indices[j]];
      };
  };

				   // compute the mean value of the
				   // sum which we have placed in each
				   // entry of the output vector
  for (unsigned int i=0; i<dof_2.n_dofs(); ++i)
    {
      Assert (touch_count[i] != 0,
	      ExcInternalError());
      
      data_2(i) /= touch_count[i];
    };
}


namespace internal
{
  namespace VectorTools
  {
#if deal_II_dimension == 1

    template <int dim>
    void
    interpolate_zero_boundary_values (const dealii::DoFHandler<dim>         &dof_handler,
                                      std::map<unsigned int,double> &boundary_values)
    {
                                       // we only need to find the
                                       // left-most and right-most
                                       // vertex and query its vertex
                                       // dof indices. that's easy :-)
      for (unsigned direction=0; direction<2; ++direction)
        {
          typename dealii::DoFHandler<dim>::cell_iterator
              cell = dof_handler.begin(0);
          while (!cell->at_boundary(direction))
            cell = cell->neighbor(direction);

          for (unsigned int i=0; i<dof_handler.get_fe().dofs_per_vertex; ++i)
            boundary_values[cell->vertex_dof_index (direction, i)] = 0.;
        }
    }


//codimension 1
    void
      interpolate_zero_boundary_values (const dealii::DoFHandler<1,2>         &dof_handler,
                                      std::map<unsigned int,double> &boundary_values)
    {
                                       // we only need to find the
                                       // left-most and right-most
                                       // vertex and query its vertex
                                       // dof indices. that's easy :-)
      for (unsigned direction=0; direction<2; ++direction)
        {
          dealii::DoFHandler<1,2>::cell_iterator
              cell = dof_handler.begin(0);
          while (!cell->at_boundary(direction))
            cell = cell->neighbor(direction);

          for (unsigned int i=0; i<dof_handler.get_fe().dofs_per_vertex; ++i)
            boundary_values[cell->vertex_dof_index (direction, i)] = 0.;
        }
    }


#else
    
    template <int dim, int spacedim>
    void
    interpolate_zero_boundary_values (const dealii::DoFHandler<dim,spacedim>       &dof_handler,
                                      std::map<unsigned int,double> &boundary_values)
    {
      const FiniteElement<dim,spacedim> &fe = dof_handler.get_fe();

				       // loop over all boundary faces
				       // to get all dof indices of
				       // dofs on the boundary. note
				       // that in 3d there are cases
				       // where a face is not at the
				       // boundary, yet one of its
				       // lines is, and we should
				       // consider the degrees of
				       // freedom on it as boundary
				       // nodes. likewise, in 2d and
				       // 3d there are cases where a
				       // cell is only at the boundary
				       // by one vertex. nevertheless,
				       // since we do not support
				       // boundaries with dimension
				       // less or equal to dim-2, each
				       // such boundary dof is also
				       // found from some other face
				       // that is actually wholly on
				       // the boundary, not only by
				       // one line or one vertex
      typename dealii::DoFHandler<dim,spacedim>::active_face_iterator
        face = dof_handler.begin_active_face(),
        endf = dof_handler.end_face();
      std::vector<unsigned int> face_dof_indices (fe.dofs_per_face);
      for (; face!=endf; ++face)
	if (face->at_boundary())
	  {
	    face->get_dof_indices (face_dof_indices);
	    for (unsigned int i=0; i<fe.dofs_per_face; ++i)
					       // enter zero boundary values
					       // for all boundary nodes
					       //
					       // we need not care about
					       // vector valued elements here,
					       // since we set all components
	      boundary_values[face_dof_indices[i]] = 0.;
	  }
    }
    
#endif
  }
}



template <int dim, class VECTOR, int spacedim>
void VectorTools::project (const Mapping<dim, spacedim>       &mapping,
			   const DoFHandler<dim,spacedim>    &dof,
			   const ConstraintMatrix   &constraints,
			   const Quadrature<dim>    &quadrature,
			   const Function<spacedim>      &function,
			   VECTOR                   &vec_result,
			   const bool                enforce_zero_boundary,
			   const Quadrature<dim-1>  &q_boundary,
			   const bool                project_to_boundary_first)
{
  Assert (dof.get_fe().n_components() == function.n_components,
	  ExcDimensionMismatch(dof.get_fe().n_components(),
			       function.n_components));

  Assert (vec_result.size() == dof.n_dofs(),
          ExcDimensionMismatch (vec_result.size(), dof.n_dofs()));
  
				   // make up boundary values
  std::map<unsigned int,double> boundary_values;
  
  if (enforce_zero_boundary == true) 
				     // no need to project boundary
				     // values, but enforce
				     // homogeneous boundary values
				     // anyway
    internal::VectorTools::
      interpolate_zero_boundary_values (dof, boundary_values);
  
  else
				     // no homogeneous boundary values
    if (project_to_boundary_first == true)
				       // boundary projection required
      {
					 // set up a list of boundary
					 // functions for the
					 // different boundary
					 // parts. We want the
					 // function to hold on
					 // all parts of the boundary
	typename FunctionMap<spacedim>::type boundary_functions;
	for (unsigned char c=0; c<255; ++c)
	  boundary_functions[c] = &function;
	project_boundary_values (dof, boundary_functions, q_boundary,
				 boundary_values);
      }

				   // set up mass matrix and right hand side
  Vector<double> vec (dof.n_dofs());
  SparsityPattern sparsity;

				   // use csp to consume less memory and to
				   // still be fast
  {
    CompressedSimpleSparsityPattern csp (dof.n_dofs(), dof.n_dofs());
    DoFTools::make_sparsity_pattern (dof, csp, constraints);

    sparsity.copy_from (csp);
  }

  SparseMatrix<double> mass_matrix (sparsity);
  Vector<double> tmp (mass_matrix.n());

				   // create mass matrix and rhs at once,
				   // which is faster.
  MatrixCreator::create_mass_matrix (mapping, dof, quadrature, mass_matrix,
				     function, tmp);

  constraints.condense (mass_matrix);
  constraints.condense (tmp);
  if (boundary_values.size() != 0)
    MatrixTools::apply_boundary_values (boundary_values,
					mass_matrix, vec, tmp,
					true);
				   // Allow for a maximum of 5*n
				   // steps to reduce the residual by
				   // 10^-12. n steps may not be
				   // sufficient, since roundoff
				   // errors may accumulate for badly
				   // conditioned matrices
  ReductionControl        control(5*tmp.size(), 0., 1e-12, false, false);
  GrowingVectorMemory<> memory;
  SolverCG<>              cg(control,memory);

  PreconditionSSOR<> prec;
  prec.initialize(mass_matrix, 1.2);
				   // solve
  cg.solve (mass_matrix, vec, tmp, prec);
  
				   // distribute solution
  constraints.distribute (vec);

                                   // copy vec into vec_result. we
                                   // can't use ve_result itself
                                   // above, since it may be of
                                   // another type than Vector<double>
                                   // and that wouldn't necessarily go
                                   // together with the matrix and
                                   // other functions
  for (unsigned int i=0; i<vec.size(); ++i)
    vec_result(i) = vec(i);
}


template <int dim, class VECTOR, int spacedim>
void VectorTools::project (const DoFHandler<dim,spacedim>    &dof,
			   const ConstraintMatrix   &constraints,
			   const Quadrature<dim>    &quadrature,
			   const Function<spacedim>      &function,
			   VECTOR                   &vec,
			   const bool                enforce_zero_boundary,
			   const Quadrature<dim-1>  &q_boundary,
			   const bool                project_to_boundary_first)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  project(StaticMappingQ1<dim,spacedim>::mapping, dof, constraints, quadrature, function, vec,
	  enforce_zero_boundary, q_boundary, project_to_boundary_first);
}




template <int dim, int spacedim>
void VectorTools::create_right_hand_side (const Mapping<dim, spacedim>    &mapping,
					  const DoFHandler<dim,spacedim> &dof_handler,
					  const Quadrature<dim> &quadrature,
					  const Function<spacedim>   &rhs_function,
					  Vector<double>        &rhs_vector)
{
  const FiniteElement<dim,spacedim> &fe  = dof_handler.get_fe();
  Assert (fe.n_components() == rhs_function.n_components,
	  ExcDimensionMismatch(fe.n_components(), rhs_function.n_components));
  Assert (rhs_vector.size() == dof_handler.n_dofs(),
	  ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
  rhs_vector = 0;
  
  UpdateFlags update_flags = UpdateFlags(update_values   |
					 update_quadrature_points |
					 update_JxW_values);
  FEValues<dim,spacedim> fe_values (mapping, fe, quadrature, update_flags);

  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points,
		     n_components  = fe.n_components();
  
  std::vector<unsigned int> dofs (dofs_per_cell);
  Vector<double> cell_vector (dofs_per_cell);

  typename DoFHandler<dim,spacedim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  if (n_components==1)
    {
      std::vector<double> rhs_values(n_q_points);
      
      for (; cell!=endc; ++cell) 
	{
	  fe_values.reinit(cell);
	  
	  const std::vector<double> &weights   = fe_values.get_JxW_values ();
	  rhs_function.value_list (fe_values.get_quadrature_points(),
				   rhs_values);
	  
	  cell_vector = 0;
	  for (unsigned int point=0; point<n_q_points; ++point)
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      cell_vector(i) += rhs_values[point] *
				fe_values.shape_value(i,point) *
				weights[point];
	
	  cell->get_dof_indices (dofs);
	  
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    rhs_vector(dofs[i]) += cell_vector(i);
	}
      
    }
  else
    {
      std::vector<Vector<double> > rhs_values(n_q_points,
					      Vector<double>(n_components));
      
      for (; cell!=endc; ++cell) 
	{
	  fe_values.reinit(cell);
	  
	  const std::vector<double> &weights   = fe_values.get_JxW_values ();
	  rhs_function.vector_value_list (fe_values.get_quadrature_points(),
					  rhs_values);
	      
	  cell_vector = 0;
					   // Use the faster code if the
					   // FiniteElement is primitive
	  if (fe.is_primitive ())
	    {
	      for (unsigned int point=0; point<n_q_points; ++point)
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  {
		    const unsigned int component
		      = fe.system_to_component_index(i).first;
		    
		    cell_vector(i) += rhs_values[point](component) *
		                      fe_values.shape_value(i,point) *
		                      weights[point];
		  }
	    }
	  else
	    {
					       // Otherwise do it the way
					       // proposed for vector valued
					       // elements
	      for (unsigned int point=0; point<n_q_points; ++point)
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
		    if (fe.get_nonzero_components(i)[comp_i])
		      {
			cell_vector(i) += rhs_values[point](comp_i) *
			                  fe_values.shape_value_component(i,point,comp_i) *
			                  weights[point];
		      }
	    }
      
	  cell->get_dof_indices (dofs);
	  
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    rhs_vector(dofs[i]) += cell_vector(i);
	}
    }
}



template <int dim, int spacedim>
void VectorTools::create_right_hand_side (const DoFHandler<dim,spacedim>    &dof_handler,
					  const Quadrature<dim>    &quadrature,
					  const Function<spacedim>      &rhs_function,
					  Vector<double>           &rhs_vector)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_right_hand_side(StaticMappingQ1<dim,spacedim>::mapping, dof_handler, quadrature,
			 rhs_function, rhs_vector);
}




template <int dim, int spacedim>
void VectorTools::create_right_hand_side (const hp::MappingCollection<dim,spacedim>    &mapping,
					  const hp::DoFHandler<dim,spacedim> &dof_handler,
					  const hp::QCollection<dim> &quadrature,
					  const Function<spacedim>   &rhs_function,
					  Vector<double>        &rhs_vector)
{
  const hp::FECollection<dim,spacedim> &fe  = dof_handler.get_fe();
  Assert (fe.n_components() == rhs_function.n_components,
	  ExcDimensionMismatch(fe.n_components(), rhs_function.n_components));
  Assert (rhs_vector.size() == dof_handler.n_dofs(),
	  ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
  rhs_vector = 0;
  
  UpdateFlags update_flags = UpdateFlags(update_values   |
					 update_quadrature_points |
					 update_JxW_values);
  hp::FEValues<dim,spacedim> x_fe_values (mapping, fe, quadrature, update_flags);

  const unsigned int n_components  = fe.n_components();
  
  std::vector<unsigned int> dofs (fe.max_dofs_per_cell());
  Vector<double> cell_vector (fe.max_dofs_per_cell());

  typename hp::DoFHandler<dim,spacedim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  if (n_components==1)
    {
      std::vector<double> rhs_values;
      
      for (; cell!=endc; ++cell) 
	{
	  x_fe_values.reinit(cell);

	  const FEValues<dim,spacedim> &fe_values = x_fe_values.get_present_fe_values();
	  
	  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
			     n_q_points    = fe_values.n_quadrature_points;
	  rhs_values.resize (n_q_points);
	  dofs.resize (dofs_per_cell);
	  cell_vector.reinit (dofs_per_cell);
	  
	  const std::vector<double> &weights   = fe_values.get_JxW_values ();
	  rhs_function.value_list (fe_values.get_quadrature_points(),
				   rhs_values);
	  
	  cell_vector = 0;
	  for (unsigned int point=0; point<n_q_points; ++point)
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      cell_vector(i) += rhs_values[point] *
				fe_values.shape_value(i,point) *
				weights[point];
	
	  cell->get_dof_indices (dofs);
	  
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    rhs_vector(dofs[i]) += cell_vector(i);
	}
      
    }
  else
    {
      std::vector<Vector<double> > rhs_values;
      
      for (; cell!=endc; ++cell)
	{
	  x_fe_values.reinit(cell);

	  const FEValues<dim,spacedim> &fe_values = x_fe_values.get_present_fe_values();
	  
	  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
			     n_q_points    = fe_values.n_quadrature_points;
	  rhs_values.resize (n_q_points,
			     Vector<double>(n_components));
	  dofs.resize (dofs_per_cell);
	  cell_vector.reinit (dofs_per_cell);
	      
	  const std::vector<double> &weights   = fe_values.get_JxW_values ();
	  rhs_function.vector_value_list (fe_values.get_quadrature_points(),
					  rhs_values);
	      
	  cell_vector = 0;
	  
					   // Use the faster code if the
					   // FiniteElement is primitive
	  if (cell->get_fe().is_primitive ())
	    {
	      for (unsigned int point=0; point<n_q_points; ++point)
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  {
		    const unsigned int component
		      = cell->get_fe().system_to_component_index(i).first;
		    
		    cell_vector(i) += rhs_values[point](component) *
				      fe_values.shape_value(i,point) *
				      weights[point];
		  }
	    }
	  else
	    {
					       // Otherwise do it the way proposed
					       // for vector valued elements
	      for (unsigned int point=0; point<n_q_points; ++point)
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
		    if (cell->get_fe().get_nonzero_components(i)[comp_i])
		      {
			cell_vector(i) += rhs_values[point](comp_i) *
					  fe_values.shape_value_component(i,point,comp_i) *
					  weights[point];
		      }
	    }
	      
	  cell->get_dof_indices (dofs);
	  
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    rhs_vector(dofs[i]) += cell_vector(i);
	}
    }
}



template <int dim, int spacedim>
void VectorTools::create_right_hand_side (const hp::DoFHandler<dim,spacedim>    &dof_handler,
					  const hp::QCollection<dim>    &quadrature,
					  const Function<spacedim>      &rhs_function,
					  Vector<double>           &rhs_vector)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_right_hand_side(hp::StaticMappingQ1<dim,spacedim>::mapping_collection, 
			 dof_handler, quadrature,
			 rhs_function, rhs_vector);
}




template <int dim, int spacedim>
void VectorTools::create_point_source_vector (const Mapping<dim, spacedim>       &mapping,
                                              const DoFHandler<dim,spacedim>    &dof_handler,
                                              const Point<spacedim>         &p,
                                              Vector<double>           &rhs_vector)
{
   Assert (rhs_vector.size() == dof_handler.n_dofs(),
           ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
   Assert (dof_handler.get_fe().n_components() == 1,
	   ExcMessage ("This function only works for scalar finite elements"));
   
   rhs_vector = 0;

   std::pair<typename DoFHandler<dim,spacedim>::active_cell_iterator, Point<spacedim> >
      cell_point =
      GridTools::find_active_cell_around_point (mapping, dof_handler, p);

   Quadrature<dim> q(GeometryInfo<dim>::project_to_unit_cell(cell_point.second));

   FEValues<dim,spacedim> fe_values(mapping, dof_handler.get_fe(),
			   q, UpdateFlags(update_values));
   fe_values.reinit(cell_point.first);

   const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   cell_point.first->get_dof_indices (local_dof_indices);

   for(unsigned int i=0; i<dofs_per_cell; i++)
      rhs_vector(local_dof_indices[i]) =  fe_values.shape_value(i,0);
}



template <int dim, int spacedim>
void VectorTools::create_point_source_vector (const DoFHandler<dim,spacedim>    &dof_handler,
                                              const Point<spacedim>         &p,
                                              Vector<double>           &rhs_vector)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_point_source_vector(StaticMappingQ1<dim,spacedim>::mapping, dof_handler,
                             p, rhs_vector);
}


template <int dim, int spacedim>
void VectorTools::create_point_source_vector (const hp::MappingCollection<dim,spacedim>       &mapping,
                                              const hp::DoFHandler<dim,spacedim>    &dof_handler,
                                              const Point<spacedim>         &p,
                                              Vector<double>           &rhs_vector)
{
   Assert (rhs_vector.size() == dof_handler.n_dofs(),
           ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
   Assert (dof_handler.get_fe().n_components() == 1,
	   ExcMessage ("This function only works for scalar finite elements"));
   
   rhs_vector = 0;

   std::pair<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator, Point<spacedim> >
      cell_point =
      GridTools::find_active_cell_around_point (mapping, dof_handler, p);

   Quadrature<dim> q(GeometryInfo<dim>::project_to_unit_cell(cell_point.second));

   FEValues<dim> fe_values(mapping[cell_point.first->active_fe_index()],
			   cell_point.first->get_fe(), q, UpdateFlags(update_values));
   fe_values.reinit(cell_point.first);

   const unsigned int dofs_per_cell = cell_point.first->get_fe().dofs_per_cell;

   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   cell_point.first->get_dof_indices (local_dof_indices);

   for(unsigned int i=0; i<dofs_per_cell; i++)
      rhs_vector(local_dof_indices[i]) =  fe_values.shape_value(i,0);
}



template <int dim, int spacedim>
void VectorTools::create_point_source_vector (const hp::DoFHandler<dim,spacedim>    &dof_handler,
                                              const Point<spacedim>         &p,
                                              Vector<double>           &rhs_vector)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_point_source_vector(hp::StaticMappingQ1<dim>::mapping_collection,
			     dof_handler,
                             p, rhs_vector);
}


#if deal_II_dimension != 1

template <int dim, int spacedim>
void
VectorTools::create_boundary_right_hand_side (const Mapping<dim, spacedim>      &mapping,
					      const DoFHandler<dim,spacedim>   &dof_handler,
					      const Quadrature<dim-1> &quadrature,
					      const Function<spacedim>     &rhs_function,
					      Vector<double>          &rhs_vector,
					      const std::set<unsigned char> &boundary_indicators)
{
  const FiniteElement<dim> &fe  = dof_handler.get_fe();
  Assert (fe.n_components() == rhs_function.n_components,
	  ExcDimensionMismatch(fe.n_components(), rhs_function.n_components));
  Assert (rhs_vector.size() == dof_handler.n_dofs(),
	  ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
  
  rhs_vector = 0;
  
  UpdateFlags update_flags = UpdateFlags(update_values   |
					 update_quadrature_points |
					 update_JxW_values);
  FEFaceValues<dim> fe_values (mapping, fe, quadrature, update_flags);

  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points,
		     n_components  = fe.n_components();
  
  std::vector<unsigned int> dofs (dofs_per_cell);
  Vector<double> cell_vector (dofs_per_cell);

  typename DoFHandler<dim,spacedim>::active_cell_iterator cell = dof_handler.begin_active(),
						 endc = dof_handler.end();

  if (n_components==1)
    {
      std::vector<double> rhs_values(n_q_points);
      
      for (; cell!=endc; ++cell)
	for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	  if (cell->face(face)->at_boundary () &&
	      (boundary_indicators.find (cell->face(face)->boundary_indicator())
	       !=
	       boundary_indicators.end()))
	    {
	      fe_values.reinit(cell, face);
	  
	      const std::vector<double> &weights   = fe_values.get_JxW_values ();
	      rhs_function.value_list (fe_values.get_quadrature_points(), rhs_values);
	      
	      cell_vector = 0;
	      for (unsigned int point=0; point<n_q_points; ++point)
		for (unsigned int i=0; i<dofs_per_cell; ++i) 
		  cell_vector(i) += rhs_values[point] *
				    fe_values.shape_value(i,point) *
				    weights[point];
	
	      cell->get_dof_indices (dofs);
	  
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		rhs_vector(dofs[i]) += cell_vector(i);
	    }
    }
  else
    {
      std::vector<Vector<double> > rhs_values(n_q_points, Vector<double>(n_components));
      
      for (; cell!=endc; ++cell) 
	for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	  if (cell->face(face)->at_boundary () &&
	      (boundary_indicators.find (cell->face(face)->boundary_indicator())
	       !=
	       boundary_indicators.end()))
	    {
	      fe_values.reinit(cell, face);
	      
	      const std::vector<double> &weights   = fe_values.get_JxW_values ();
	      rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);
	      
	      cell_vector = 0;
	      
					       // Use the faster code if the
					       // FiniteElement is primitive
	      if (fe.is_primitive ())
		{		  
		  for (unsigned int point=0; point<n_q_points; ++point)
		    for (unsigned int i=0; i<dofs_per_cell; ++i)
		      {
			const unsigned int component
			  = fe.system_to_component_index(i).first;
			
			cell_vector(i) += rhs_values[point](component) *
				          fe_values.shape_value(i,point) *
				          weights[point];
		      }
		}
	      else
		{
						   // And the full featured
						   // code, if vector valued
						   // FEs are used
		  for (unsigned int point=0; point<n_q_points; ++point)
		    for (unsigned int i=0; i<dofs_per_cell; ++i)
		      for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
			if (fe.get_nonzero_components(i)[comp_i])
			  {
			    cell_vector(i)
			      += rhs_values[point](comp_i) *
			      fe_values.shape_value_component(i,point,comp_i) *
			      weights[point];
			  }
		}
		  
	      cell->get_dof_indices (dofs);
	      
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		rhs_vector(dofs[i]) += cell_vector(i);
	    }
    }
}

#else

// Implementation for 1D
template <int dim, int spacedim>
void
VectorTools::create_boundary_right_hand_side (const Mapping<dim, spacedim>    &,
					      const DoFHandler<dim,spacedim> &,
					      const Quadrature<dim-1> &,
					      const Function<spacedim>   &,
					      Vector<double>      &,
					      const std::set<unsigned char> &)
{
  Assert (false, ExcImpossibleInDim(dim));
}

#endif

template <int dim, int spacedim>
void
VectorTools::create_boundary_right_hand_side (const DoFHandler<dim,spacedim>   &dof_handler,
					      const Quadrature<dim-1> &quadrature,
					      const Function<spacedim>     &rhs_function,
					      Vector<double>          &rhs_vector,
					      const std::set<unsigned char> &boundary_indicators)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));

  create_boundary_right_hand_side(StaticMappingQ1<dim>::mapping, dof_handler,
				  quadrature,
				  rhs_function, rhs_vector,
				  boundary_indicators);
}



#if deal_II_dimension != 1

template <int dim, int spacedim>
void
VectorTools::create_boundary_right_hand_side (const hp::MappingCollection<dim,spacedim>      &mapping,
					      const hp::DoFHandler<dim,spacedim>   &dof_handler,
					      const hp::QCollection<dim-1> &quadrature,
					      const Function<spacedim>     &rhs_function,
					      Vector<double>          &rhs_vector,
					      const std::set<unsigned char> &boundary_indicators)
{
  const hp::FECollection<dim> &fe  = dof_handler.get_fe();
  Assert (fe.n_components() == rhs_function.n_components,
	  ExcDimensionMismatch(fe.n_components(), rhs_function.n_components));
  Assert (rhs_vector.size() == dof_handler.n_dofs(),
	  ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
  
  rhs_vector = 0;
  
  UpdateFlags update_flags = UpdateFlags(update_values   |
					 update_quadrature_points |
					 update_JxW_values);
  hp::FEFaceValues<dim> x_fe_values (mapping, fe, quadrature, update_flags);

  const unsigned int n_components  = fe.n_components();
  
  std::vector<unsigned int> dofs (fe.max_dofs_per_cell());
  Vector<double> cell_vector (fe.max_dofs_per_cell());

  typename hp::DoFHandler<dim,spacedim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  if (n_components==1)
    {
      std::vector<double> rhs_values;
      
      for (; cell!=endc; ++cell)
	for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	  if (cell->face(face)->at_boundary () &&
	      (boundary_indicators.find (cell->face(face)->boundary_indicator())
	       !=
	       boundary_indicators.end()))
	    {
	      x_fe_values.reinit(cell, face);

	      const FEFaceValues<dim> &fe_values = x_fe_values.get_present_fe_values();

	      const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
				 n_q_points    = fe_values.n_quadrature_points;
	      rhs_values.resize (n_q_points);
	      
	      const std::vector<double> &weights   = fe_values.get_JxW_values ();
	      rhs_function.value_list (fe_values.get_quadrature_points(), rhs_values);
	      
	      cell_vector = 0;
	      for (unsigned int point=0; point<n_q_points; ++point)
		for (unsigned int i=0; i<dofs_per_cell; ++i) 
		  cell_vector(i) += rhs_values[point] *
				    fe_values.shape_value(i,point) *
				    weights[point];
	
	      cell->get_dof_indices (dofs);
	  
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		rhs_vector(dofs[i]) += cell_vector(i);
	    }
    }
  else
    {
      std::vector<Vector<double> > rhs_values;
      
      for (; cell!=endc; ++cell) 
	for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	  if (cell->face(face)->at_boundary () &&
	      (boundary_indicators.find (cell->face(face)->boundary_indicator())
	       !=
	       boundary_indicators.end()))
	    {
	      x_fe_values.reinit(cell, face);

	      const FEFaceValues<dim> &fe_values = x_fe_values.get_present_fe_values();

	      const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
				 n_q_points    = fe_values.n_quadrature_points;
	      rhs_values.resize (n_q_points, Vector<double>(n_components));
	      
	      const std::vector<double> &weights   = fe_values.get_JxW_values ();
	      rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);
	      
	      cell_vector = 0;
	      
					       // Use the faster code if the
					       // FiniteElement is primitive
	      if (cell->get_fe().is_primitive ())
		{		  
		  for (unsigned int point=0; point<n_q_points; ++point)
		    for (unsigned int i=0; i<dofs_per_cell; ++i)
		      {
			const unsigned int component
			  = cell->get_fe().system_to_component_index(i).first;
			
			cell_vector(i) += rhs_values[point](component) *
				          fe_values.shape_value(i,point) *
				          weights[point];
		      }
		}
	      else
		{
						   // And the full featured
						   // code, if vector valued
						   // FEs are used
		  for (unsigned int point=0; point<n_q_points; ++point)
		    for (unsigned int i=0; i<dofs_per_cell; ++i)
		      for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
			if (cell->get_fe().get_nonzero_components(i)[comp_i])
			  {
			    cell_vector(i)
			      += rhs_values[point](comp_i) *
			      fe_values.shape_value_component(i,point,comp_i) *
			      weights[point];
			  }
		}
		  
	      cell->get_dof_indices (dofs);
	      
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		rhs_vector(dofs[i]) += cell_vector(i);
	    }
    }
}

#else

// Implementation for 1D
template <int dim, int spacedim>
void
VectorTools::create_boundary_right_hand_side (const hp::MappingCollection<dim,spacedim>    &,
					      const hp::DoFHandler<dim,spacedim> &,
					      const hp::QCollection<dim-1> &,
					      const Function<spacedim>   &,
					      Vector<double>      &,
					      const std::set<unsigned char> &)
{
  Assert (false, ExcImpossibleInDim(1));
}

#endif

template <int dim, int spacedim>
void
VectorTools::create_boundary_right_hand_side (const hp::DoFHandler<dim,spacedim>   &dof_handler,
					      const hp::QCollection<dim-1> &quadrature,
					      const Function<spacedim>     &rhs_function,
					      Vector<double>          &rhs_vector,
					      const std::set<unsigned char> &boundary_indicators)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_boundary_right_hand_side(hp::StaticMappingQ1<dim>::mapping_collection,
				  dof_handler, quadrature,
				  rhs_function, rhs_vector,
				  boundary_indicators);
}



// ----------- interpolate_boundary_values for std::map --------------------

#if deal_II_dimension == 1

//TODO[?] Actually the Mapping object should be a MappingCollection object for the
// hp::DoFHandler.

//template <int dim, template <int, int> class DH, int spacedim>

template <class DH>
void
VectorTools::interpolate_boundary_values (const Mapping<DH::dimension, DH::space_dimension> &,
					  const DH              &dof,
					  const unsigned char    boundary_component,
					  const Function<DH::space_dimension> &boundary_function,
					  std::map<unsigned int,double> &boundary_values,
					  const std::vector<bool>       &component_mask_)
{
  const unsigned int dim=DH::dimension;
  const unsigned int spacedim=DH::space_dimension;

  Assert (boundary_component != 255,
	  ExcInvalidBoundaryIndicator());
  Assert ((component_mask_.size() == 0) ||
	  (component_mask_.size() == dof.get_fe().n_components()),
	  ExcMessage ("The number of components in the mask has to be either "
		      "zero or equal to the number of components in the finite "
		      "element."));
    
				   // check whether boundary values at
				   // the left or right boundary of
				   // the line are
				   // requested. direction denotes
				   // the neighboring direction in
				   // which we seek the boundary,
				   // i.e. 0 is left boundary and 1 is
				   // right.
  const unsigned int direction = boundary_component;
  Assert (direction < 2, ExcInvalidBoundaryIndicator());
  
				   // first find the outermost active
				   // cell by first traversing the coarse
				   // grid to its end and then going
				   // to the children
  typename DH::cell_iterator outermost_cell = dof.begin(0);
  while (outermost_cell->neighbor(direction).state() == IteratorState::valid)
    outermost_cell = outermost_cell->neighbor(direction);
  
  while (outermost_cell->has_children())
    outermost_cell = outermost_cell->child(direction);

                                   // get the FE corresponding to this
                                   // cell
  const FiniteElement<dim,spacedim> &fe = outermost_cell->get_fe();
  Assert (fe.n_components() == boundary_function.n_components,
	  ExcDimensionMismatch(fe.n_components(), boundary_function.n_components));

				   // set the component mask to either
				   // the original value or a vector
				   // of trues
  const std::vector<bool> component_mask ((component_mask_.size() == 0) ?
					  std::vector<bool> (fe.n_components(), true) :
					  component_mask_);
  Assert (std::count(component_mask.begin(), component_mask.end(), true) > 0,
	  ExcNoComponentSelected());

				   // now set the value of the
				   // outermost degree of
				   // freedom. setting also
				   // creates the entry in the map
				   // if it did not exist
				   // beforehand
				   //
				   // save some time by requesting
				   // values only once for each point,
				   // irrespective of the number of
				   // components of the function
  Vector<double> function_values (fe.n_components());
  if (fe.n_components() == 1)
    function_values(0)
      = boundary_function.value (outermost_cell->vertex(direction));
  else
    boundary_function.vector_value (outermost_cell->vertex(direction),
				    function_values);
  
  for (unsigned int i=0; i<fe.dofs_per_vertex; ++i)
    if (component_mask[fe.face_system_to_component_index(i).first])
      boundary_values[outermost_cell->vertex_dof_index(direction,i)]
	= function_values(fe.face_system_to_component_index(i).first);
}


//TODO[?] Actually the Mapping object should be a MappingCollection object for the
// hp::DoFHandler.

// Implementation for 1D
template <class DH>
void
VectorTools::interpolate_boundary_values (const Mapping<DH::dimension, DH::space_dimension>              &mapping,
					  const DH         &dof,
					  const typename FunctionMap<DH::space_dimension>::type    &function_map,
					  std::map<unsigned int,double> &boundary_values,
					  const std::vector<bool>       &component_mask)
{
  for (typename FunctionMap<DH::space_dimension>::type::const_iterator i=function_map.begin();
       i!=function_map.end(); ++i)
    interpolate_boundary_values (mapping, dof, i->first, *i->second,
				 boundary_values, component_mask);
}


//TODO[?] Actually the Mapping object should be a MappingCollection object for the
// hp::DoFHandler.


#else


template <class DH>
void
VectorTools::
interpolate_boundary_values (const Mapping<DH::dimension, DH::space_dimension>            &mapping,
                             const DH                 &dof,
                             const typename FunctionMap<DH::space_dimension>::type &function_map,
                             std::map<unsigned int,double> &boundary_values,
                             const std::vector<bool>       &component_mask_)
{
  const unsigned int dim=DH::dimension;

  Assert ((component_mask_.size() == 0) ||
	  (component_mask_.size() == dof.get_fe().n_components()),
	  ExcMessage ("The number of components in the mask has to be either "
		      "zero or equal to the number of components in the finite "
		      "element."));


				   // if for whatever reason we were
				   // passed an empty map, return
				   // immediately
  if (function_map.size() == 0)
    return;
  
  Assert (function_map.find(255) == function_map.end(),
	  ExcInvalidBoundaryIndicator());

  const unsigned int        n_components = DoFTools::n_components(dof);
  const bool                fe_is_system = (n_components != 1);

  for (typename FunctionMap<DH::space_dimension>::type::const_iterator i=function_map.begin();
       i!=function_map.end(); ++i)
    Assert (n_components == i->second->n_components,
	    ExcDimensionMismatch(n_components, i->second->n_components));

				   // set the component mask to either
				   // the original value or a vector
				   // of trues
  const std::vector<bool> component_mask ((component_mask_.size() == 0) ?
					  std::vector<bool> (n_components, true) :
					  component_mask_);
  Assert (std::count(component_mask.begin(), component_mask.end(), true) > 0,
	  ExcNoComponentSelected());

				   // field to store the indices
  std::vector<unsigned int> face_dofs;
  face_dofs.reserve (DoFTools::max_dofs_per_face(dof));

  std::vector<Point<DH::space_dimension> >  dof_locations;
  dof_locations.reserve (DoFTools::max_dofs_per_face(dof));

				   // array to store the values of
				   // the boundary function at the
				   // boundary points. have to arrays
				   // for scalar and vector functions
				   // to use the more efficient one
				   // respectively
  std::vector<double>          dof_values_scalar;
  std::vector<Vector<double> > dof_values_system;
  dof_values_scalar.reserve (DoFTools::max_dofs_per_face (dof));
  dof_values_system.reserve (DoFTools::max_dofs_per_face (dof));

				   // before we start with the loop
				   // over all cells create an
				   // hp::FEValues object that holds
				   // the interpolation points of all
				   // finite elements that may ever be
				   // in use
  hp::FECollection<dim> finite_elements (dof.get_fe());
  hp::QCollection<dim-1>  q_collection;
  for (unsigned int f=0; f<finite_elements.size(); ++f)
    {
      const FiniteElement<dim> &fe = finite_elements[f];
      
				       // generate a quadrature rule
				       // on the face from the unit
				       // support points. this will be
				       // used to obtain the
				       // quadrature points on the
				       // real cell's face
				       //
				       // to do this, we check whether
				       // the FE has support points on
				       // the face at all:
      if (fe.has_face_support_points())
	q_collection.push_back (Quadrature<dim-1>(fe.get_unit_face_support_points()));
      else
	{
					   // if not, then we should
					   // try a more clever
					   // way. the idea is that a
					   // finite element may not
					   // offer support points for
					   // all its shape functions,
					   // but maybe only some. if
					   // it offers support points
					   // for the components we
					   // are interested in in
					   // this function, then
					   // that's fine. if not, the
					   // function we call in the
					   // finite element will
					   // raise an exception. the
					   // support points for the
					   // other shape functions
					   // are left uninitialized
					   // (well, initialized by
					   // the default
					   // constructor), since we
					   // don't need them anyway.
					   //
					   // As a detour, we must
					   // make sure we only query
					   // face_system_to_component_index
					   // if the index corresponds
					   // to a primitive shape
					   // function. since we know
					   // that all the components
					   // we are interested in are
					   // primitive (by the above
					   // check), we can safely
					   // put such a check in
					   // front
	  std::vector<Point<dim-1> > unit_support_points (fe.dofs_per_face);

	  for (unsigned int i=0; i<fe.dofs_per_face; ++i)
	    if (fe.is_primitive (fe.face_to_equivalent_cell_index(i)))
	      if (component_mask[fe.face_system_to_component_index(i).first]
		  == true)
		unit_support_points[i] = fe.unit_face_support_point(i);
	
	  q_collection.push_back (Quadrature<dim-1>(unit_support_points));
        }    
    }
				   // now that we have a q_collection
				   // object with all the right
				   // quadrature points, create an
				   // hp::FEFaceValues object that we
				   // can use to evaluate the boundary
				   // values at
  hp::MappingCollection<dim> mapping_collection (mapping);
  hp::FEFaceValues<dim> x_fe_values (mapping_collection, finite_elements, q_collection,
				     update_quadrature_points);
  
  typename DH::active_cell_iterator cell = dof.begin_active(),
				    endc = dof.end();
  for (; cell!=endc; ++cell)
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell;
	 ++face_no)
      {
        const FiniteElement<dim,DH::space_dimension> &fe = cell->get_fe();

					 // we can presently deal only with
					 // primitive elements for boundary
					 // values. this does not preclude
					 // us using non-primitive elements
					 // in components that we aren't
					 // interested in, however. make
					 // sure that all shape functions
					 // that are non-zero for the
					 // components we are interested in,
					 // are in fact primitive
	for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
	  {
	    const std::vector<bool> &nonzero_component_array
	      = cell->get_fe().get_nonzero_components (i);
	    for (unsigned int c=0; c<n_components; ++c)
	      if ((nonzero_component_array[c] == true)
		  &&
		  (component_mask[c] == true))
		Assert (cell->get_fe().is_primitive (i),
			ExcMessage ("This function can only deal with requested boundary "
				    "values that correspond to primitive (scalar) base "
				    "elements"));
	  }
	
	typename DH::face_iterator face = cell->face(face_no);
	const unsigned char boundary_component = face->boundary_indicator();
	if (function_map.find(boundary_component) != function_map.end()) 
	  {
					     // face is of the right component
	    x_fe_values.reinit(cell, face_no);
	    const FEFaceValues<dim> &fe_values = x_fe_values.get_present_fe_values();

					     // get indices, physical location and
					     // boundary values of dofs on this
					     // face
            face_dofs.resize (fe.dofs_per_face);
	    face->get_dof_indices (face_dofs, cell->active_fe_index());
	    const std::vector<Point<DH::space_dimension> > &dof_locations
              = fe_values.get_quadrature_points ();
	    
	    if (fe_is_system)
	      {
                                                 // resize
                                                 // array. avoid
                                                 // construction of a
                                                 // memory allocating
                                                 // temporary if
                                                 // possible
                if (dof_values_system.size() < fe.dofs_per_face)
                  dof_values_system.resize (fe.dofs_per_face,
                                            Vector<double>(fe.n_components()));
                else
                  dof_values_system.resize (fe.dofs_per_face);
                
		function_map.find(boundary_component)->second
                  ->vector_value_list (dof_locations, dof_values_system);
		
						 // enter those dofs
						 // into the list that
						 // match the
						 // component
						 // signature. avoid
						 // the usual
						 // complication that
						 // we can't just use
						 // *_system_to_component_index
						 // for non-primitive
						 // FEs
		for (unsigned int i=0; i<face_dofs.size(); ++i)
                  {
                    unsigned int component;
                    if (fe.is_primitive())
                      component = fe.face_system_to_component_index(i).first;
                    else
                      {
                                                         // non-primitive
                                                         // case. make
                                                         // sure that
                                                         // this
                                                         // particular
                                                         // shape
                                                         // function
                                                         // _is_
                                                         // primitive,
                                                         // and get at
                                                         // it's
                                                         // component. use
                                                         // usual
                                                         // trick to
                                                         // transfer
                                                         // face dof
                                                         // index to
                                                         // cell dof
                                                         // index
                        const unsigned int cell_i
                          = (dim == 1 ?
                             i
                             :
                             (dim == 2 ?
                              (i<2*fe.dofs_per_vertex ? i : i+2*fe.dofs_per_vertex)
                              :
                              (dim == 3 ?
                               (i<4*fe.dofs_per_vertex ?
                                i
                                :
                                (i<4*fe.dofs_per_vertex+4*fe.dofs_per_line ?
                                 i+4*fe.dofs_per_vertex
                                 :
                                 i+4*fe.dofs_per_vertex+8*fe.dofs_per_line))
                               :
                               numbers::invalid_unsigned_int)));
                        Assert (cell_i < fe.dofs_per_cell, ExcInternalError());

                                                         // make sure
                                                         // that if
                                                         // this is
                                                         // not a
                                                         // primitive
                                                         // shape function,
                                                         // then all
                                                         // the
                                                         // corresponding
                                                         // components
                                                         // in the
                                                         // mask are
                                                         // not set
                        if (!fe.is_primitive(cell_i))
                          for (unsigned int c=0; c<n_components; ++c)
                            if (fe.get_nonzero_components(cell_i)[c])
                              Assert (component_mask[c] == false,
                                      FETools::ExcFENotPrimitive());

                                                         // let's pick
                                                         // the first
                                                         // of
                                                         // possibly
                                                         // more than
                                                         // one
                                                         // non-zero
                                                         // components. if
                                                         // shape
                                                         // function
                                                         // is
                                                         // non-primitive,
                                                         // then we
                                                         // will
                                                         // ignore the
                                                         // result in
                                                         // the
                                                         // following
                                                         // anyway,
                                                         // otherwise
                                                         // there's
                                                         // only one
                                                         // non-zero
                                                         // component
                                                         // which we
                                                         // will use
                        component = (std::find (fe.get_nonzero_components(cell_i).begin(),
                                                fe.get_nonzero_components(cell_i).end(),
                                                true)
                                     -
                                     fe.get_nonzero_components(cell_i).begin());
                      }
                    
                    if (component_mask[component] == true)
                      boundary_values[face_dofs[i]] = dof_values_system[i](component);
                  } 
	      }
	    else
					       // fe has only one component,
					       // so save some computations
	      {
						 // get only the one component that
						 // this function has
                dof_values_scalar.resize (fe.dofs_per_face);
		function_map.find(boundary_component)->second
                  ->value_list (dof_locations, dof_values_scalar, 0);
		
						 // enter into list
		
		for (unsigned int i=0; i<face_dofs.size(); ++i)
		  boundary_values[face_dofs[i]] = dof_values_scalar[i];
	      }
	  }
      }
}



template <class DH>
void
VectorTools::interpolate_boundary_values (const Mapping<DH::dimension, DH::space_dimension>            &mapping,
					  const DH                 &dof,
					  const unsigned char            boundary_component,
					  const Function<DH::space_dimension>           &boundary_function,
					  std::map<unsigned int,double> &boundary_values,
					  const std::vector<bool>       &component_mask)
{
  typename FunctionMap<DH::space_dimension>::type function_map;
  function_map[boundary_component] = &boundary_function;
  interpolate_boundary_values (mapping, dof, function_map, boundary_values,
			       component_mask);
}

#endif
  

template <class DH>
void
VectorTools::interpolate_boundary_values (const DH                 &dof,
					  const unsigned char            boundary_component,
					  const Function<DH::space_dimension>           &boundary_function,
					  std::map<unsigned int,double> &boundary_values,
					  const std::vector<bool>       &component_mask)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  interpolate_boundary_values(StaticMappingQ1<DH::dimension,DH::space_dimension>::mapping, 
			      dof, boundary_component,
			      boundary_function, boundary_values, component_mask);
}



template <class DH>
void
VectorTools::interpolate_boundary_values (const DH                 &dof,
					  const typename FunctionMap<DH::space_dimension>::type &function_map,
					  std::map<unsigned int,double> &boundary_values,
					  const std::vector<bool>       &component_mask)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  interpolate_boundary_values(StaticMappingQ1<DH::dimension,DH::space_dimension>::mapping, 
			      dof, function_map,
			      boundary_values, component_mask);
}




// ----------- interpolate_boundary_values for ConstraintMatrix --------------

// TODO (M.K.): There is a lot of duplicated code with the above
// function. We should unify all these interpolate_boundary_values functions
// in one way or the other.

#if deal_II_dimension == 1

//TODO[?] Actually the Mapping object should be a MappingCollection object for the
// hp::DoFHandler.

//template <int dim, template <int, int> class DH, int spacedim>

template <class DH>
void
VectorTools::interpolate_boundary_values (const Mapping<DH::dimension, DH::space_dimension> &,
					  const DH                &dof,
					  const unsigned char      boundary_component,
					  const Function<DH::space_dimension> &boundary_function,
					  ConstraintMatrix        &constraints,
					  const std::vector<bool> &component_mask_)
{
  const unsigned int dim=DH::dimension;
  const unsigned int spacedim=DH::space_dimension;

  Assert (boundary_component != 255,
	  ExcInvalidBoundaryIndicator());
  Assert ((component_mask_.size() == 0) ||
	  (component_mask_.size() == dof.get_fe().n_components()),
	  ExcMessage ("The number of components in the mask has to be either "
		      "zero or equal to the number of components in the finite "
		      "element."));
    
				   // check whether boundary values at
				   // the left or right boundary of
				   // the line are
				   // requested. direction denotes
				   // the neighboring direction in
				   // which we seek the boundary,
				   // i.e. 0 is left boundary and 1 is
				   // right.
  const unsigned int direction = boundary_component;
  Assert (direction < 2, ExcInvalidBoundaryIndicator());
  
				   // first find the outermost active
				   // cell by first traversing the coarse
				   // grid to its end and then going
				   // to the children
  typename DH::cell_iterator outermost_cell = dof.begin(0);
  while (outermost_cell->neighbor(direction).state() == IteratorState::valid)
    outermost_cell = outermost_cell->neighbor(direction);
  
  while (outermost_cell->has_children())
    outermost_cell = outermost_cell->child(direction);

                                   // get the FE corresponding to this
                                   // cell
  const FiniteElement<dim,spacedim> &fe = outermost_cell->get_fe();
  Assert (fe.n_components() == boundary_function.n_components,
	  ExcDimensionMismatch(fe.n_components(), boundary_function.n_components));

				   // set the component mask to either
				   // the original value or a vector
				   // of trues
  const std::vector<bool> component_mask ((component_mask_.size() == 0) ?
					  std::vector<bool> (fe.n_components(), true) :
					  component_mask_);
  Assert (std::count(component_mask.begin(), component_mask.end(), true) > 0,
	  ExcNoComponentSelected());

				   // now set the value of the
				   // outermost degree of
				   // freedom. setting also
				   // creates the entry in the map
				   // if it did not exist
				   // beforehand
				   //
				   // save some time by requesting
				   // values only once for each point,
				   // irrespective of the number of
				   // components of the function
  Vector<double> function_values (fe.n_components());
  if (fe.n_components() == 1)
    function_values(0)
      = boundary_function.value (outermost_cell->vertex(direction));
  else
    boundary_function.vector_value (outermost_cell->vertex(direction),
				    function_values);
  
  for (unsigned int i=0; i<fe.dofs_per_vertex; ++i)
    if (component_mask[fe.face_system_to_component_index(i).first])
      {
				   // TODO: should we clear the other
				   // entries in the line here?
	const unsigned int row = outermost_cell->vertex_dof_index(direction,i);
	constraints.add_line (row);
	constraints.set_inhomogeneity (row, 
		    function_values(fe.face_system_to_component_index(i).first));
      }
}


//TODO[?] Actually the Mapping object should be a MappingCollection object for the
// hp::DoFHandler.

// Implementation for 1D
template <class DH>
void
VectorTools::interpolate_boundary_values 
  (const Mapping<DH::dimension, DH::space_dimension>     &mapping,
   const DH                                              &dof,
   const typename FunctionMap<DH::space_dimension>::type &function_map,
   ConstraintMatrix                                      &constraints,
   const std::vector<bool>                               &component_mask)
{
  for (typename FunctionMap<DH::space_dimension>::type::const_iterator i=function_map.begin();
       i!=function_map.end(); ++i)
    interpolate_boundary_values (mapping, dof, i->first, *i->second,
				 constraints, component_mask);
}


//TODO[?] Actually the Mapping object should be a MappingCollection object for the
// hp::DoFHandler.


#else


template <class DH>
void
VectorTools::interpolate_boundary_values 
 (const Mapping<DH::dimension, DH::space_dimension>     &mapping,
  const DH                                              &dof,
  const typename FunctionMap<DH::space_dimension>::type &function_map,
  ConstraintMatrix                                      &constraints,
  const std::vector<bool>                               &component_mask_)
{
  const unsigned int dim=DH::dimension;

  Assert ((component_mask_.size() == 0) ||
	  (component_mask_.size() == dof.get_fe().n_components()),
	  ExcMessage ("The number of components in the mask has to be either "
		      "zero or equal to the number of components in the finite "
		      "element."));


				   // if for whatever reason we were
				   // passed an empty map, return
				   // immediately
  if (function_map.size() == 0)
    return;
  
  Assert (function_map.find(255) == function_map.end(),
	  ExcInvalidBoundaryIndicator());

  const unsigned int        n_components = DoFTools::n_components(dof);
  const bool                fe_is_system = (n_components != 1);

  for (typename FunctionMap<DH::space_dimension>::type::const_iterator i=function_map.begin();
       i!=function_map.end(); ++i)
    Assert (n_components == i->second->n_components,
	    ExcDimensionMismatch(n_components, i->second->n_components));

				   // set the component mask to either
				   // the original value or a vector
				   // of trues
  const std::vector<bool> component_mask ((component_mask_.size() == 0) ?
					  std::vector<bool> (n_components, true) :
					  component_mask_);
  Assert (std::count(component_mask.begin(), component_mask.end(), true) > 0,
	  ExcNoComponentSelected());

				   // field to store the indices
  std::vector<unsigned int> face_dofs;
  face_dofs.reserve (DoFTools::max_dofs_per_face(dof));

  std::vector<Point<DH::space_dimension> >  dof_locations;
  dof_locations.reserve (DoFTools::max_dofs_per_face(dof));

				   // array to store the values of the
				   // boundary function at the boundary
				   // points. have to arrays for scalar and
				   // vector functions to use the more
				   // efficient one respectively
  std::vector<double>          dof_values_scalar;
  std::vector<Vector<double> > dof_values_system;
  dof_values_scalar.reserve (DoFTools::max_dofs_per_face (dof));
  dof_values_system.reserve (DoFTools::max_dofs_per_face (dof));

				   // before we start with the loop over all
				   // cells create an hp::FEValues object
				   // that holds the interpolation points of
				   // all finite elements that may ever be
				   // in use
  hp::FECollection<dim> finite_elements (dof.get_fe());
  hp::QCollection<dim-1>  q_collection;
  for (unsigned int f=0; f<finite_elements.size(); ++f)
    {
      const FiniteElement<dim> &fe = finite_elements[f];
      
				       // generate a quadrature rule on the
				       // face from the unit support
				       // points. this will be used to
				       // obtain the quadrature points on
				       // the real cell's face
				       //
				       // to do this, we check whether
				       // the FE has support points on
				       // the face at all:
      if (fe.has_face_support_points())
	q_collection.push_back (Quadrature<dim-1>(fe.get_unit_face_support_points()));
      else
	{
					   // if not, then we should try a
					   // more clever way. the idea is
					   // that a finite element may not
					   // offer support points for all
					   // its shape functions, but maybe
					   // only some. if it offers
					   // support points for the
					   // components we are interested
					   // in in this function, then
					   // that's fine. if not, the
					   // function we call in the finite
					   // element will raise an
					   // exception. the support points
					   // for the other shape functions
					   // are left uninitialized (well,
					   // initialized by the default
					   // constructor), since we don't
					   // need them anyway.
					   //
					   // As a detour, we must make sure
					   // we only query
					   // face_system_to_component_index
					   // if the index corresponds to a
					   // primitive shape
					   // function. since we know that
					   // all the components we are
					   // interested in are primitive
					   // (by the above check), we can
					   // safely put such a check in
					   // front
	  std::vector<Point<dim-1> > unit_support_points (fe.dofs_per_face);

	  for (unsigned int i=0; i<fe.dofs_per_face; ++i)
	    if (fe.is_primitive (fe.face_to_equivalent_cell_index(i)))
	      if (component_mask[fe.face_system_to_component_index(i).first]
		  == true)
		unit_support_points[i] = fe.unit_face_support_point(i);
	
	  q_collection.push_back (Quadrature<dim-1>(unit_support_points));
        }    
    }
				   // now that we have a q_collection object
				   // with all the right quadrature points,
				   // create an hp::FEFaceValues object that
				   // we can use to evaluate the boundary
				   // values at
  hp::MappingCollection<dim> mapping_collection (mapping);
  hp::FEFaceValues<dim> x_fe_values (mapping_collection, finite_elements, q_collection,
				     update_quadrature_points);
  
  typename DH::active_cell_iterator cell = dof.begin_active(),
				    endc = dof.end();
  for (; cell!=endc; ++cell)
    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell;
	 ++face_no)
      {
        const FiniteElement<dim,DH::space_dimension> &fe = cell->get_fe();

					 // we can presently deal only with
					 // primitive elements for boundary
					 // values. this does not preclude
					 // us using non-primitive elements
					 // in components that we aren't
					 // interested in, however. make
					 // sure that all shape functions
					 // that are non-zero for the
					 // components we are interested in,
					 // are in fact primitive
	for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
	  {
	    const std::vector<bool> &nonzero_component_array
	      = cell->get_fe().get_nonzero_components (i);
	    for (unsigned int c=0; c<n_components; ++c)
	      if ((nonzero_component_array[c] == true)
		  &&
		  (component_mask[c] == true))
		Assert (cell->get_fe().is_primitive (i),
			ExcMessage ("This function can only deal with requested boundary "
				    "values that correspond to primitive (scalar) base "
				    "elements"));
	  }
	
	typename DH::face_iterator face = cell->face(face_no);
	const unsigned char boundary_component = face->boundary_indicator();
	if (function_map.find(boundary_component) != function_map.end()) 
	  {
					     // face is of the right
					     // component
	    x_fe_values.reinit(cell, face_no);
	    const FEFaceValues<dim> &fe_values = x_fe_values.get_present_fe_values();

					     // get indices, physical
					     // location and boundary values
					     // of dofs on this face
            face_dofs.resize (fe.dofs_per_face);
	    face->get_dof_indices (face_dofs, cell->active_fe_index());
	    const std::vector<Point<DH::space_dimension> > &dof_locations
              = fe_values.get_quadrature_points ();
	    
	    if (fe_is_system)
	      {
                                                 // resize array. avoid
                                                 // construction of a memory
                                                 // allocating temporary if
                                                 // possible
                if (dof_values_system.size() < fe.dofs_per_face)
                  dof_values_system.resize (fe.dofs_per_face,
                                            Vector<double>(fe.n_components()));
                else
                  dof_values_system.resize (fe.dofs_per_face);
                
		function_map.find(boundary_component)->second
                  ->vector_value_list (dof_locations, dof_values_system);
		
						 // enter those dofs into
						 // the list that match the
						 // component
						 // signature. avoid the
						 // usual complication that
						 // we can't just use
						 // *_system_to_component_index
						 // for non-primitive FEs
		for (unsigned int i=0; i<face_dofs.size(); ++i)
                  {
                    unsigned int component;
                    if (fe.is_primitive())
                      component = fe.face_system_to_component_index(i).first;
                    else
                      {
                                                         // non-primitive
                                                         // case. make sure
                                                         // that this
                                                         // particular shape
                                                         // function _is_
                                                         // primitive, and
                                                         // get at it's
                                                         // component. use
                                                         // usual trick to
                                                         // transfer face
                                                         // dof index to
                                                         // cell dof index
                        const unsigned int cell_i
                          = (dim == 1 ?
                             i
                             :
                             (dim == 2 ?
                              (i<2*fe.dofs_per_vertex ? i : i+2*fe.dofs_per_vertex)
                              :
                              (dim == 3 ?
                               (i<4*fe.dofs_per_vertex ?
                                i
                                :
                                (i<4*fe.dofs_per_vertex+4*fe.dofs_per_line ?
                                 i+4*fe.dofs_per_vertex
                                 :
                                 i+4*fe.dofs_per_vertex+8*fe.dofs_per_line))
                               :
                               numbers::invalid_unsigned_int)));
                        Assert (cell_i < fe.dofs_per_cell, ExcInternalError());

                                                         // make sure that
                                                         // if this is not a
                                                         // primitive shape
                                                         // function, then
                                                         // all the
                                                         // corresponding
                                                         // components in
                                                         // the mask are not
                                                         // set
                        if (!fe.is_primitive(cell_i))
                          for (unsigned int c=0; c<n_components; ++c)
                            if (fe.get_nonzero_components(cell_i)[c])
                              Assert (component_mask[c] == false,
                                      FETools::ExcFENotPrimitive());

                                                         // let's pick the
                                                         // first of
                                                         // possibly more
                                                         // than one
                                                         // non-zero
                                                         // components. if
                                                         // shape function
                                                         // is
                                                         // non-primitive,
                                                         // then we will
                                                         // ignore the
                                                         // result in the
                                                         // following
                                                         // anyway,
                                                         // otherwise
                                                         // there's only one
                                                         // non-zero
                                                         // component which
                                                         // we will use
                        component = (std::find (fe.get_nonzero_components(cell_i).begin(),
                                                fe.get_nonzero_components(cell_i).end(),
                                                true)
                                     -
                                     fe.get_nonzero_components(cell_i).begin());
                      }
                    
                    if (component_mask[component] == true)
		      {
				   // TODO: check whether we should clear
				   // the current line first...
			constraints.add_line (face_dofs[i]);
			constraints.set_inhomogeneity (face_dofs[i],
						       dof_values_system[i](component));
		      }
                  } 
	      }
	    else
					       // fe has only one component,
					       // so save some computations
	      {
						 // get only the one
						 // component that this
						 // function has
                dof_values_scalar.resize (fe.dofs_per_face);
		function_map.find(boundary_component)->second
                  ->value_list (dof_locations, dof_values_scalar, 0);
		
						 // enter into list
		
		for (unsigned int i=0; i<face_dofs.size(); ++i)
				   // TODO: check whether we should clear
				   // the current line first...
		  {
		    constraints.add_line (face_dofs[i]);
		    constraints.set_inhomogeneity (face_dofs[i],
						   dof_values_scalar[i]);
		  }
	      }
	  }
      }
}



template <class DH>
void
VectorTools::interpolate_boundary_values 
  (const Mapping<DH::dimension, DH::space_dimension> &mapping,
   const DH                                          &dof,
   const unsigned char                                boundary_component,
   const Function<DH::space_dimension>               &boundary_function,
   ConstraintMatrix                                  &constraints,
   const std::vector<bool>                           &component_mask)
{
  typename FunctionMap<DH::space_dimension>::type function_map;
  function_map[boundary_component] = &boundary_function;
  interpolate_boundary_values (mapping, dof, function_map, constraints,
			       component_mask);
}

#endif
  

template <class DH>
void
VectorTools::interpolate_boundary_values 
  (const DH                            &dof,
   const unsigned char                  boundary_component,
   const Function<DH::space_dimension> &boundary_function,
   ConstraintMatrix                    &constraints,
   const std::vector<bool>             &component_mask)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  interpolate_boundary_values(StaticMappingQ1<DH::dimension,DH::space_dimension>::mapping, 
			      dof, boundary_component,
			      boundary_function, constraints, component_mask);
}



template <class DH>
void
VectorTools::interpolate_boundary_values 
  (const DH                                              &dof,
   const typename FunctionMap<DH::space_dimension>::type &function_map,
   ConstraintMatrix                                      &constraints,
   const std::vector<bool>                               &component_mask)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  interpolate_boundary_values(StaticMappingQ1<DH::dimension,DH::space_dimension>::mapping, 
			      dof, function_map,
			      constraints, component_mask);
}




// -------- implementation for project_boundary_values with std::map --------

#if deal_II_dimension == 1

// Implementation for 1D
template <int dim, int spacedim>
void
VectorTools::project_boundary_values (const Mapping<dim, spacedim>       &mapping,
				      const DoFHandler<dim,spacedim>    &dof,
				      const typename FunctionMap<spacedim>::type &boundary_functions,
				      const Quadrature<dim-1>  &,
				      std::map<unsigned int,double> &boundary_values,
				      std::vector<unsigned int> component_mapping)
{
  Assert (component_mapping.size() == 0, ExcNotImplemented());
				   // projection in 1d is equivalent
				   // to interpolation
  interpolate_boundary_values (mapping, dof, boundary_functions,
			       boundary_values, std::vector<bool>());
}

#else


template <int dim, int spacedim>
void
VectorTools::project_boundary_values (const Mapping<dim, spacedim>       &mapping,
				      const DoFHandler<dim,spacedim>    &dof,
				      const typename FunctionMap<spacedim>::type &boundary_functions,
				      const Quadrature<dim-1>  &q,
				      std::map<unsigned int,double> &boundary_values,
				      std::vector<unsigned int> component_mapping)
{
//TODO:[?] In VectorTools::project_boundary_values, no condensation of sparsity
//    structures, matrices and right hand sides or distribution of
//    solution vectors is performed. This is ok for dim<3 because then
//    there are no constrained nodes on the boundary, but is not
//    acceptable for higher dimensions. Fix this.

  if (component_mapping.size() == 0)
    {
      AssertDimension (dof.get_fe().n_components(), boundary_functions.begin()->second->n_components);
				       // I still do not see why i
				       // should create another copy
				       // here
      component_mapping.resize(dof.get_fe().n_components());
      for (unsigned int i= 0 ;i < component_mapping.size() ; ++i)
	component_mapping[i] = i;
    }
  else
    AssertDimension (dof.get_fe().n_components(), component_mapping.size());
  
  std::vector<unsigned int> dof_to_boundary_mapping;
  std::set<unsigned char> selected_boundary_components;
  for (typename FunctionMap<spacedim>::type::const_iterator i=boundary_functions.begin();
       i!=boundary_functions.end(); ++i)
    selected_boundary_components.insert (i->first);
  
  DoFTools::map_dof_to_boundary_indices (dof, selected_boundary_components,
					 dof_to_boundary_mapping);
  
				   // Done if no degrees of freedom on
				   // the boundary
  if (dof.n_boundary_dofs (boundary_functions) == 0)
    return;
				   // set up sparsity structure
  SparsityPattern sparsity(dof.n_boundary_dofs (boundary_functions),
			   dof.max_couplings_between_boundary_dofs());
  DoFTools::make_boundary_sparsity_pattern (dof,
					    boundary_functions,
					    dof_to_boundary_mapping,
					    sparsity);

				   // note: for three or more dimensions, there
				   // may be constrained nodes on the boundary
				   // in this case the boundary mass matrix has
				   // to be condensed and the solution is to
				   // be distributed afterwards, which is not
				   // yet implemented. The reason for this is
				   // that we cannot simply use the condense
				   // family of functions, since the matrices
				   // and vectors do not use the global
				   // numbering but rather the boundary
				   // numbering, i.e. the condense function
				   // needs to use another indirection. There
				   // should be not many technical problems,
				   // but it needs to be implemented
  if (dim>=3)
    {
#ifdef DEBUG
// Assert that there are no hanging nodes at the boundary      
      int level = -1;
      for (typename DoFHandler<dim,spacedim>::active_cell_iterator cell = dof.begin_active();
	   cell != dof.end(); ++cell)
	for (unsigned int f=0;f<GeometryInfo<dim>::faces_per_cell;++f)
	  {
	    if (cell->at_boundary(f))
	      {
		if (level == -1)
		  level = cell->level();
		else
		  {		          
		    Assert (level == cell->level(), ExcNotImplemented());
		  }
	      }
	  }
#endif
    }
  sparsity.compress();
  

				   // make mass matrix and right hand side
  SparseMatrix<double> mass_matrix(sparsity);
  Vector<double>       rhs(sparsity.n_rows());


  MatrixCreator::create_boundary_mass_matrix (mapping, dof, q, 
					      mass_matrix, boundary_functions,
					      rhs, dof_to_boundary_mapping, (const Function<spacedim>*) 0,
					      component_mapping);

				   // For certain weird elements,
				   // there might be degrees of
				   // freedom on the boundary, but
				   // their shape functions do not
				   // have support there. Let's
				   // eliminate them here.

//TODO: Maybe we should figure out if the element really needs this
  
  FilteredMatrix<Vector<double> > filtered_mass_matrix(mass_matrix, true);
  FilteredMatrix<Vector<double> > filtered_precondition;
  std::vector<bool> excluded_dofs(mass_matrix.m(), false);

  double max_element = 0.;
  for (unsigned int i=0;i<mass_matrix.m();++i)
    if (mass_matrix.diag_element(i) > max_element)
      max_element = mass_matrix.diag_element(i);
  
  for (unsigned int i=0;i<mass_matrix.m();++i)
    if (mass_matrix.diag_element(i) < 1.e-8 * max_element)
      {
	filtered_mass_matrix.add_constraint(i, 0.);
	filtered_precondition.add_constraint(i, 0.);
	mass_matrix.diag_element(i) = 1.;
	excluded_dofs[i] = true;
      }
  
  Vector<double> boundary_projection (rhs.size());

				   // Allow for a maximum of 5*n
				   // steps to reduce the residual by
				   // 10^-12. n steps may not be
				   // sufficient, since roundoff
				   // errors may accumulate for badly
				   // conditioned matrices
  ReductionControl        control(5*rhs.size(), 0., 1.e-12, false, false);
  GrowingVectorMemory<> memory;
  SolverCG<>              cg(control,memory);

  PreconditionSSOR<> prec;
  prec.initialize(mass_matrix, 1.2);
  filtered_precondition.initialize(prec, true);
				   // solve
  cg.solve (filtered_mass_matrix, boundary_projection, rhs, filtered_precondition);
  filtered_precondition.apply_constraints(boundary_projection, true);
  filtered_precondition.clear();  
				   // fill in boundary values
  for (unsigned int i=0; i<dof_to_boundary_mapping.size(); ++i)
    if (dof_to_boundary_mapping[i] != DoFHandler<dim,spacedim>::invalid_dof_index
    && ! excluded_dofs[dof_to_boundary_mapping[i]])
				       // this dof is on one of the
				       // interesting boundary parts
				       //
				       // remember: i is the global dof
				       // number, dof_to_boundary_mapping[i]
				       // is the number on the boundary and
				       // thus in the solution vector
      boundary_values[i] = boundary_projection(dof_to_boundary_mapping[i]);
}

#endif

template <int dim, int spacedim>
void
VectorTools::project_boundary_values (const DoFHandler<dim,spacedim>    &dof,
				      const typename FunctionMap<spacedim>::type &boundary_functions,
				      const Quadrature<dim-1>  &q,
				      std::map<unsigned int,double> &boundary_values,
				      std::vector<unsigned int> component_mapping)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  project_boundary_values(StaticMappingQ1<dim,spacedim>::mapping, dof, boundary_functions, q,
			  boundary_values, component_mapping);
}



// ----- implementation for project_boundary_values with ConstraintMatrix -----

#if deal_II_dimension == 1

// Implementation for 1D
template <int dim, int spacedim>
void
VectorTools::project_boundary_values (const Mapping<dim, spacedim>       &mapping,
				      const DoFHandler<dim,spacedim>    &dof,
				      const typename FunctionMap<spacedim>::type &boundary_functions,
				      const Quadrature<dim-1>  &,
				      ConstraintMatrix &constraints,
				      std::vector<unsigned int> component_mapping)
{
  Assert (component_mapping.size() == 0, ExcNotImplemented());
				   // projection in 1d is equivalent
				   // to interpolation
  interpolate_boundary_values (mapping, dof, boundary_functions,
			       constraints, std::vector<bool>());
}

#else


template <int dim, int spacedim>
void
VectorTools::project_boundary_values (const Mapping<dim, spacedim>       &mapping,
				      const DoFHandler<dim,spacedim>    &dof,
				      const typename FunctionMap<spacedim>::type &boundary_functions,
				      const Quadrature<dim-1>  &q,
				      ConstraintMatrix &constraints,
				      std::vector<unsigned int> component_mapping)
{
//TODO:[?] In VectorTools::project_boundary_values, no condensation of sparsity
//    structures, matrices and right hand sides or distribution of
//    solution vectors is performed. This is ok for dim<3 because then
//    there are no constrained nodes on the boundary, but is not
//    acceptable for higher dimensions. Fix this.

  if (component_mapping.size() == 0)
    {
      AssertDimension (dof.get_fe().n_components(), boundary_functions.begin()->second->n_components);
				       // I still do not see why i
				       // should create another copy
				       // here
      component_mapping.resize(dof.get_fe().n_components());
      for (unsigned int i= 0 ;i < component_mapping.size() ; ++i)
	component_mapping[i] = i;
    }
  else
    AssertDimension (dof.get_fe().n_components(), component_mapping.size());
  
  std::vector<unsigned int> dof_to_boundary_mapping;
  std::set<unsigned char> selected_boundary_components;
  for (typename FunctionMap<spacedim>::type::const_iterator i=boundary_functions.begin();
       i!=boundary_functions.end(); ++i)
    selected_boundary_components.insert (i->first);
  
  DoFTools::map_dof_to_boundary_indices (dof, selected_boundary_components,
					 dof_to_boundary_mapping);
  
				   // Done if no degrees of freedom on
				   // the boundary
  if (dof.n_boundary_dofs (boundary_functions) == 0)
    return;
				   // set up sparsity structure
  SparsityPattern sparsity(dof.n_boundary_dofs (boundary_functions),
			   dof.max_couplings_between_boundary_dofs());
  DoFTools::make_boundary_sparsity_pattern (dof,
					    boundary_functions,
					    dof_to_boundary_mapping,
					    sparsity);

				   // note: for three or more dimensions, there
				   // may be constrained nodes on the boundary
				   // in this case the boundary mass matrix has
				   // to be condensed and the solution is to
				   // be distributed afterwards, which is not
				   // yet implemented. The reason for this is
				   // that we cannot simply use the condense
				   // family of functions, since the matrices
				   // and vectors do not use the global
				   // numbering but rather the boundary
				   // numbering, i.e. the condense function
				   // needs to use another indirection. There
				   // should be not many technical problems,
				   // but it needs to be implemented
  if (dim>=3)
    {
#ifdef DEBUG
// Assert that there are no hanging nodes at the boundary      
      int level = -1;
      for (typename DoFHandler<dim,spacedim>::active_cell_iterator cell = dof.begin_active();
	   cell != dof.end(); ++cell)
	for (unsigned int f=0;f<GeometryInfo<dim>::faces_per_cell;++f)
	  {
	    if (cell->at_boundary(f))
	      {
		if (level == -1)
		  level = cell->level();
		else
		  {		          
		    Assert (level == cell->level(), ExcNotImplemented());
		  }
	      }
	  }
#endif
    }
  sparsity.compress();
  

				   // make mass matrix and right hand side
  SparseMatrix<double> mass_matrix(sparsity);
  Vector<double>       rhs(sparsity.n_rows());


  MatrixCreator::create_boundary_mass_matrix (mapping, dof, q, 
					      mass_matrix, boundary_functions,
					      rhs, dof_to_boundary_mapping, (const Function<spacedim>*) 0,
					      component_mapping);

				   // For certain weird elements,
				   // there might be degrees of
				   // freedom on the boundary, but
				   // their shape functions do not
				   // have support there. Let's
				   // eliminate them here.

//TODO: Maybe we should figure out if the element really needs this
  
  FilteredMatrix<Vector<double> > filtered_mass_matrix(mass_matrix, true);
  FilteredMatrix<Vector<double> > filtered_precondition;
  std::vector<bool> excluded_dofs(mass_matrix.m(), false);

  double max_element = 0.;
  for (unsigned int i=0;i<mass_matrix.m();++i)
    if (mass_matrix.diag_element(i) > max_element)
      max_element = mass_matrix.diag_element(i);
  
  for (unsigned int i=0;i<mass_matrix.m();++i)
    if (mass_matrix.diag_element(i) < 1.e-8 * max_element)
      {
	filtered_mass_matrix.add_constraint(i, 0.);
	filtered_precondition.add_constraint(i, 0.);
	mass_matrix.diag_element(i) = 1.;
	excluded_dofs[i] = true;
      }
  
  Vector<double> boundary_projection (rhs.size());

				   // Allow for a maximum of 5*n
				   // steps to reduce the residual by
				   // 10^-12. n steps may not be
				   // sufficient, since roundoff
				   // errors may accumulate for badly
				   // conditioned matrices
  ReductionControl        control(5*rhs.size(), 0., 1.e-12, false, false);
  GrowingVectorMemory<> memory;
  SolverCG<>              cg(control,memory);

  PreconditionSSOR<> prec;
  prec.initialize(mass_matrix, 1.2);
  filtered_precondition.initialize(prec, true);
				   // solve
  cg.solve (filtered_mass_matrix, boundary_projection, rhs, filtered_precondition);
  filtered_precondition.apply_constraints(boundary_projection, true);
  filtered_precondition.clear();  
				   // fill in boundary values
  for (unsigned int i=0; i<dof_to_boundary_mapping.size(); ++i)
    if (dof_to_boundary_mapping[i] != DoFHandler<dim,spacedim>::invalid_dof_index
    && ! excluded_dofs[dof_to_boundary_mapping[i]])
				       // this dof is on one of the
				       // interesting boundary parts
				       //
				       // remember: i is the global dof
				       // number, dof_to_boundary_mapping[i]
				       // is the number on the boundary and
				       // thus in the solution vector
      {
				   // TODO: check whether we should clear
				   // the entries in this line first.
	constraints.add_line (i);
	constraints.set_inhomogeneity (i, boundary_projection(dof_to_boundary_mapping[i]));
      }
}

#endif

template <int dim, int spacedim>
void
VectorTools::project_boundary_values (const DoFHandler<dim,spacedim>    &dof,
				      const typename FunctionMap<spacedim>::type &boundary_functions,
				      const Quadrature<dim-1>  &q,
				      ConstraintMatrix &constraints,
				      std::vector<unsigned int> component_mapping)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  project_boundary_values(StaticMappingQ1<dim,spacedim>::mapping, dof, boundary_functions, q,
			  constraints, component_mapping);
}




namespace internal
{
  namespace VectorTools
  {
				     /**
				      * A structure that stores the dim DoF
				      * indices that correspond to a
				      * vector-valued quantity at a single
				      * support point.
				      */
    template <int dim>
    struct VectorDoFTuple
    {
	unsigned int dof_indices[dim];

	VectorDoFTuple ()
	  {
	    for (unsigned int i=0; i<dim; ++i)
	      dof_indices[i] = numbers::invalid_unsigned_int;
	  }
	
	
	bool operator < (const VectorDoFTuple<dim> &other) const
	  {
	    for (unsigned int i=0; i<dim; ++i)
	      if (dof_indices[i] < other.dof_indices[i])
		return true;
	      else
		if (dof_indices[i] > other.dof_indices[i])
		  return false;
	    return false;
	  }

	bool operator == (const VectorDoFTuple<dim> &other) const
	  {
	    for (unsigned int i=0; i<dim; ++i)
	      if (dof_indices[i] != other.dof_indices[i])
		return false;

	    return true;
	  }

	bool operator != (const VectorDoFTuple<dim> &other) const
	  {
	    return ! (*this == other);
	  }
    };



				     /**
				      * Add the constraint
				      * $\vec n \cdot \vec u = 0$
				      * to the list of constraints.
				      *
				      * Here, $\vec u$ is represented
				      * by the set of given DoF
				      * indices, and $\vec n$ by the
				      * vector specified as the second
				      * argument.
				      */
    template <int dim>
    void
    add_constraint (const VectorDoFTuple<dim> &dof_indices,
		    const Tensor<1,dim>       &constraining_vector,
		    ConstraintMatrix          &constraints)
    {	    
				       // choose the DoF that has the
				       // largest component in the
				       // constraining_vector as the
				       // one to be constrained as
				       // this makes the process
				       // stable in cases where the
				       // constraining_vector has the
				       // form n=(1,0) or n=(0,1)
				       //
				       // we get constraints of the form
				       //   x0 = a_1*x1 + a2*x2 + ...
				       // if one of the weights is
				       // essentially zero then skip
				       // this part. the ConstraintMatrix
				       // can also deal with cases like
				       //   x0 = 0
				       // if necessary
      switch (dim)
	{
	  case 2:
	  {
	    if (std::fabs(constraining_vector[0]) > std::fabs(constraining_vector[1]))
	      {
		constraints.add_line (dof_indices.dof_indices[0]);

		if (std::fabs (constraining_vector[1]/constraining_vector[0])
		    > std::numeric_limits<double>::epsilon())
		  constraints.add_entry (dof_indices.dof_indices[0],
					 dof_indices.dof_indices[1],
					 -constraining_vector[1]/constraining_vector[0]);
	      }
	    else
	      {
		constraints.add_line (dof_indices.dof_indices[1]);

		if (std::fabs (constraining_vector[0]/constraining_vector[1])
		    > std::numeric_limits<double>::epsilon())
		  constraints.add_entry (dof_indices.dof_indices[1],
					 dof_indices.dof_indices[0],
					 -constraining_vector[0]/constraining_vector[1]);
	      }
	    break;
	  }

	  case 3:
	  {
	    if ((std::fabs(constraining_vector[0]) >= std::fabs(constraining_vector[1]))
		&&
		(std::fabs(constraining_vector[0]) >= std::fabs(constraining_vector[2])))
	      {
		constraints.add_line (dof_indices.dof_indices[0]);
		
		if (std::fabs (constraining_vector[1]/constraining_vector[0])
		    > std::numeric_limits<double>::epsilon())
		  constraints.add_entry (dof_indices.dof_indices[0],
					 dof_indices.dof_indices[1],
					 -constraining_vector[1]/constraining_vector[0]);

		if (std::fabs (constraining_vector[2]/constraining_vector[0])
		    > std::numeric_limits<double>::epsilon())
		  constraints.add_entry (dof_indices.dof_indices[0],
					 dof_indices.dof_indices[2],
					 -constraining_vector[2]/constraining_vector[0]);
	      }
	    else
	      if ((std::fabs(constraining_vector[1]) >= std::fabs(constraining_vector[0]))
		  &&
		  (std::fabs(constraining_vector[1]) >= std::fabs(constraining_vector[2])))
		{
		  constraints.add_line (dof_indices.dof_indices[1]);

		  if (std::fabs (constraining_vector[0]/constraining_vector[1])
		      > std::numeric_limits<double>::epsilon())
		    constraints.add_entry (dof_indices.dof_indices[1],
					   dof_indices.dof_indices[0],
					   -constraining_vector[0]/constraining_vector[1]);

		  if (std::fabs (constraining_vector[2]/constraining_vector[1])
		      > std::numeric_limits<double>::epsilon())
		    constraints.add_entry (dof_indices.dof_indices[1],
					   dof_indices.dof_indices[2],
					   -constraining_vector[2]/constraining_vector[1]);
		}
	      else
		{
		  constraints.add_line (dof_indices.dof_indices[2]);

		  if (std::fabs (constraining_vector[0]/constraining_vector[2])
		      > std::numeric_limits<double>::epsilon())
		    constraints.add_entry (dof_indices.dof_indices[2],
					   dof_indices.dof_indices[0],
					   -constraining_vector[0]/constraining_vector[2]);

		  if (std::fabs (constraining_vector[1]/constraining_vector[2])
		      > std::numeric_limits<double>::epsilon())
		    constraints.add_entry (dof_indices.dof_indices[2],
					   dof_indices.dof_indices[1],
					   -constraining_vector[1]/constraining_vector[2]);
		}
	    
	    break;
	  }

	  default:
		Assert (false, ExcNotImplemented());
	}
      
    }



				     /**
				      * Given a vector, compute a set
				      * of dim-1 vectors that are
				      * orthogonal to the first one
				      * and mutually orthonormal as
				      * well.
				      */
    template <int dim>
    void
    compute_orthonormal_vectors (const Tensor<1,dim> &vector,
				 Tensor<1,dim> (&orthonormals)[dim-1])
    {	    
      switch (dim)
	{
	  case 3:
	  {
					     // to do this in 3d, take
					     // one vector that is
					     // guaranteed to be not
					     // aligned with the
					     // average tangent and
					     // form the cross
					     // product. this yields
					     // one vector that is
					     // certainly
					     // perpendicular to the
					     // tangent; then take the
					     // cross product between
					     // this vector and the
					     // tangent and get one
					     // vector that is
					     // perpendicular to both

					     // construct a
					     // temporary vector
					     // by swapping the
					     // larger two
					     // components and
					     // flipping one
					     // sign; this can
					     // not be collinear
					     // with the average
					     // tangent
	    Tensor<1,dim> tmp = vector;
	    if ((std::fabs(tmp[0]) > std::fabs(tmp[1]))
		&&
		(std::fabs(tmp[0]) > std::fabs(tmp[2])))
	      {
						 // entry zero
						 // is the
						 // largest
		if ((std::fabs(tmp[1]) > std::fabs(tmp[2])))
		  std::swap (tmp[0], tmp[1]);
		else
		  std::swap (tmp[0], tmp[2]);

		tmp[0] *= -1;
	      }
	    else if ((std::fabs(tmp[1]) > std::fabs(tmp[0]))
		     &&
		     (std::fabs(tmp[1]) > std::fabs(tmp[2])))
	      {
						 // entry one
						 // is the
						 // largest
		if ((std::fabs(tmp[0]) > std::fabs(tmp[2])))
		  std::swap (tmp[1], tmp[0]);
		else
		  std::swap (tmp[1], tmp[2]);

		tmp[1] *= -1;
	      }
	    else
	      {
						 // entry two
						 // is the
						 // largest
		if ((std::fabs(tmp[0]) > std::fabs(tmp[1])))
		  std::swap (tmp[2], tmp[0]);
		else
		  std::swap (tmp[2], tmp[1]);

		tmp[2] *= -1;
	      }

	    Assert (std::fabs(vector * tmp) < 1e-12,
		    ExcInternalError());

					     // now compute the
					     // two normals
	    cross_product (orthonormals[0], vector, tmp);
	    cross_product (orthonormals[1], vector, orthonormals[0]);

	    break;
	  }

	  default:
		Assert (false, ExcNotImplemented());
	}      
    }
    
  }
}



template <int dim, template <int, int> class DH, int spacedim>
void
VectorTools::compute_no_normal_flux_constraints (const DH<dim,spacedim>         &dof_handler,
						 const unsigned int     first_vector_component,
						 const std::set<unsigned char> &boundary_ids,
						 ConstraintMatrix      &constraints,
						 const Mapping<dim, spacedim>    &mapping)
{
  Assert (dim > 1,
	  ExcMessage ("This function is not useful in 1d because it amounts "
		      "to imposing Dirichlet values on the vector-valued "
		      "quantity."));
	  
  const FiniteElement<dim> &fe = dof_handler.get_fe();

  std::vector<unsigned int> face_dofs (fe.dofs_per_face);
  std::vector<Point<spacedim> >  dof_locations  (fe.dofs_per_face);

				   // have a map that stores normal vectors
				   // for each vector-dof tuple we want to
				   // constrain. since we can get at the same
				   // vector dof tuple more than once (for
				   // example if it is located at a vertex
				   // that we visit from all adjacent cells),
				   // we will want to average later on the
				   // normal vectors computed on different
				   // cells as described in the documentation
				   // of this function. however, we can only
				   // average if the contributions came from
				   // different cells, whereas we want to
				   // constrain twice or more in case the
				   // contributions came from different faces
				   // of the same cell. consequently, we also
				   // have to store which cell a normal vector
				   // was computed on
  typedef 
    std::multimap<internal::VectorTools::VectorDoFTuple<dim>,
    std::pair<Tensor<1,dim>, typename DH<dim,spacedim>::active_cell_iterator> >
    DoFToNormalsMap;

  DoFToNormalsMap dof_to_normals_map;

				   // now loop over all cells and all faces
  typename DH<dim,spacedim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    for (unsigned int face_no=0; face_no < GeometryInfo<dim>::faces_per_cell;
	 ++face_no)
      if (boundary_ids.find(cell->face(face_no)->boundary_indicator())
	  != boundary_ids.end())
	{
	  typename DH<dim,spacedim>::face_iterator face = cell->face(face_no);

					   // get the indices of the
					   // dofs on this cell...
	  face->get_dof_indices (face_dofs, cell->active_fe_index());

					   // ...and the normal
					   // vectors at the locations
					   // where they are defined:
	  const std::vector<Point<dim-1> > &
	    unit_support_points = fe.get_unit_face_support_points();  
	  Quadrature<dim-1> aux_quad (unit_support_points);
	  FEFaceValues<dim> fe_values (mapping, fe, aux_quad,
				       update_normal_vectors);
	  fe_values.reinit(cell, face_no);

					   // then identify which of
					   // them correspond to the
					   // selected set of vector
					   // components
	  for (unsigned int i=0; i<face_dofs.size(); ++i)
	    if (fe.face_system_to_component_index(i).first ==
		first_vector_component)
	      {
						 // find corresponding other
						 // components of vector
		internal::VectorTools::VectorDoFTuple<dim> vector_dofs;
		vector_dofs.dof_indices[0] = face_dofs[i];
		
		for (unsigned int k=0; k<fe.dofs_per_face; ++k)
		  if ((k != i)
		      &&
		      (unit_support_points[k] == unit_support_points[i])
		      &&
		      (fe.face_system_to_component_index(k).first >=
		      first_vector_component)
		      &&
		      (fe.face_system_to_component_index(k).first <
		       first_vector_component + dim))
		    vector_dofs.dof_indices[fe.face_system_to_component_index(k).first -
					    first_vector_component]
		      = face_dofs[k];
	  
		for (unsigned int d=0; d<dim; ++d)
		  Assert (vector_dofs.dof_indices[d] < dof_handler.n_dofs(),
			  ExcInternalError());

						 // and enter the
						 // (dofs,(normal_vector,cell))
						 // entry into the map
		dof_to_normals_map
		  .insert (std::make_pair (vector_dofs,
					   std::make_pair (fe_values.normal_vector(i),
							   cell)));
	      }
	}

				   // Now do something with the
				   // collected information. To this
				   // end, loop through all sets of
				   // pairs (dofs,normal_vector) and
				   // identify which entries belong to
				   // the same set of dofs and then do
				   // as described in the
				   // documentation, i.e. either
				   // average the normal vector or
				   // don't for this particular set of
				   // dofs
  typename DoFToNormalsMap::const_iterator
    p = dof_to_normals_map.begin();
  
  while (p != dof_to_normals_map.end())
    {
				       // first find the range of entries in
				       // the multimap that corresponds to the
				       // same vector-dof tuple. as usual, we
				       // define the range half-open. the
				       // first entry of course is 'p'
      typename DoFToNormalsMap::const_iterator same_dof_range[2]
	= { p };
      for (++p; p != dof_to_normals_map.end(); ++p)
	if (p->first != same_dof_range[0]->first)
	  {
	    same_dof_range[1] = p;
	    break;
	  }
      if (p == dof_to_normals_map.end())
	same_dof_range[1] = dof_to_normals_map.end();

				       // now compute the reverse mapping: for
				       // each of the cells that contributed
				       // to the current set of vector dofs,
				       // add up the normal vectors. the
				       // values of the map are pairs of
				       // normal vectors and number of cells
				       // that have contributed
      typedef
	std::map
	<typename DH<dim,spacedim>::active_cell_iterator,
	std::pair<Tensor<1,dim>, unsigned int> >
	CellToNormalsMap;

      CellToNormalsMap cell_to_normals_map;
      for (typename DoFToNormalsMap::const_iterator
	     q = same_dof_range[0];
	   q != same_dof_range[1]; ++q)
	if (cell_to_normals_map.find (q->second.second)
	    == cell_to_normals_map.end())
	  cell_to_normals_map[q->second.second]
	    = std::make_pair (q->second.first, 1U);
	else
	  {
	    const Tensor<1,dim> old_normal
	      = cell_to_normals_map[q->second.second].first;
	    const unsigned int old_count
	      = cell_to_normals_map[q->second.second].second;

	    Assert (old_count > 0, ExcInternalError());

					     // in the same entry,
					     // store again the now
					     // averaged normal vector
					     // and the new count
	    cell_to_normals_map[q->second.second]
	      = std::make_pair ((old_normal * old_count + q->second.first) / (old_count + 1),
				old_count + 1);
	  }

      Assert (cell_to_normals_map.size() >= 1, ExcInternalError());

				       // count the maximum number of
				       // contributions from each cell
      unsigned int max_n_contributions_per_cell = 1;
      for (typename CellToNormalsMap::const_iterator
	     x = cell_to_normals_map.begin();
	   x != cell_to_normals_map.end(); ++x)
	max_n_contributions_per_cell
	  = std::max (max_n_contributions_per_cell,
		      x->second.second);

				       // verify that each cell can have only
				       // contributed at most dim times, since
				       // that is the maximum number of faces
				       // that come together at a single place
      Assert (max_n_contributions_per_cell <= dim, ExcInternalError());

      switch (max_n_contributions_per_cell)
	{
					   // first deal with the case that a
					   // number of cells all have
					   // registered that they have a
					   // normal vector defined at the
					   // location of a given vector dof,
					   // and that each of them have
					   // encountered this vector dof
					   // exactly once while looping over
					   // all their faces. as stated in
					   // the documentation, this is the
					   // case where we want to simply
					   // average over all normal vectors
	  case 1:
	  {
	    
					     // compute the average
					     // normal vector from all
					     // the ones that have the
					     // same set of dofs. we
					     // could add them up and
					     // divide them by the
					     // number of additions,
					     // or simply normalize
					     // them right away since
					     // we want them to have
					     // unit length anyway
	    Tensor<1,dim> normal;
	    for (typename CellToNormalsMap::const_iterator
		   x = cell_to_normals_map.begin();
		 x != cell_to_normals_map.end(); ++x)
	      normal += x->second.first;
	    normal /= normal.norm();	    
	    
					     // then construct constraints
					     // from this:
	    const internal::VectorTools::VectorDoFTuple<dim> &
	      dof_indices = same_dof_range[0]->first;	    
	    internal::VectorTools::add_constraint (dof_indices, normal,
						   constraints);

	    break;
	  }


					   // this is the slightly
					   // more complicated case
					   // that a single cell has
					   // contributed with exactly
					   // DIM normal vectors to
					   // the same set of vector
					   // dofs. this is what
					   // happens in a corner in
					   // 2d and 3d (but not on an
					   // edge in 3d, where we
					   // have only 2, i.e. <DIM,
					   // contributions. Here we
					   // do not want to average
					   // the normal
					   // vectors. Since we have
					   // DIM contributions, let's
					   // assume (and verify) that
					   // they are in fact all
					   // linearly independent; in
					   // that case, all vector
					   // components are
					   // constrained and we need
					   // to set them to zero
	  case dim:
	  {
					     // assert that indeed
					     // only a single cell has
					     // contributed
	    Assert (cell_to_normals_map.size() == 1,
		    ExcInternalError());

					     // check linear
					     // independence by
					     // computing the
					     // determinant of the
					     // matrix created from
					     // all the normal
					     // vectors. if they are
					     // linearly independent,
					     // then the determinant
					     // is nonzero. if they
					     // are orthogonal, then
					     // the matrix is in fact
					     // equal to 1 (since they
					     // are all unit vectors);
					     // make sure the
					     // determinant is larger
					     // than 1e-3 to avoid
					     // cases where cells are
					     // degenerate
	    {
	      Tensor<2,dim> t;

	      typename DoFToNormalsMap::const_iterator x = same_dof_range[0];
	      for (unsigned int i=0; i<dim; ++i, ++x)
		for (unsigned int j=0; j<dim; ++j)
		  t[i][j] = x->second.first[j];
	      
	      Assert (std::fabs(determinant (t)) > 1e-3,
		      ExcMessage("Found a set of normal vectors that are nearly collinear."));
	    }
	    
					     // so all components of
					     // this vector dof are
					     // constrained. enter
					     // this into the
					     // constraint matrix
	    for (unsigned int i=0; i<dim; ++i)
	      {
		constraints.add_line (same_dof_range[0]->first.dof_indices[i]);
						 // no add_entries here
	      }

	    break;
	  }
	  

					   // this is the case of an
					   // edge contribution in 3d,
					   // i.e. the vector is
					   // constrained in two
					   // directions but not the
					   // third.
	  default:
	  {
	    Assert (dim >= 3, ExcNotImplemented());
	    Assert (max_n_contributions_per_cell == 2, ExcInternalError());

					     // as described in the
					     // documentation, let us
					     // first collect what
					     // each of the cells
					     // contributed at the
					     // current point. we use
					     // a std::list instead of
					     // a std::set (which
					     // would be more natural)
					     // because std::set
					     // requires that the
					     // stored elements are
					     // comparable with
					     // operator<
	    typedef
	      std::map<typename DH<dim,spacedim>::active_cell_iterator, std::list<Tensor<1,dim> > >
	      CellContributions;
	    CellContributions cell_contributions;

	    for (typename DoFToNormalsMap::const_iterator
		   q = same_dof_range[0];
		 q != same_dof_range[1]; ++q)
	      cell_contributions[q->second.second].push_back (q->second.first);
	    Assert (cell_contributions.size() >= 1, ExcInternalError());
	    
					     // now for each cell that
					     // has contributed
					     // determine the number
					     // of normal vectors it
					     // has contributed. we
					     // currently only
					     // implement if this is
					     // dim-1 for all cells
					     // (if a single cell has
					     // contributed dim, or if
					     // all adjacent cells
					     // have contributed 1
					     // normal vector, this is
					     // already handled above)
					     //
					     // for each contributing
					     // cell compute the
					     // tangential vector that
					     // remains unconstrained
	    std::list<Tensor<1,dim> > tangential_vectors;
	    for (typename CellContributions::const_iterator
		   contribution = cell_contributions.begin();
		 contribution != cell_contributions.end();
		 ++contribution)
	      {
		Assert (contribution->second.size() == dim-1, ExcNotImplemented());

		Tensor<1,dim> normals[dim-1];
		{
		  unsigned int index=0;
		  for (typename std::list<Tensor<1,dim> >::const_iterator
			 t = contribution->second.begin();
		       t != contribution->second.end();
		       ++t, ++index)
		    normals[index] = *t;
		  Assert (index == dim-1, ExcInternalError());
		}

						 // calculate the
						 // tangent as the
						 // outer product of
						 // the normal vectors
		Tensor<1,dim> tangent;
		switch (dim)
		  {
		    case 3:
							   // take
							   // cross
							   // product
							   // between
							   // normals[0]
							   // and
							   // normals[1]. write
							   // it in
							   // the
							   // current
							   // form to
							   // make
							   // sure
							   // that
							   // compilers
							   // don't
							   // warn
							   // about
							   // out-of-bounds
							   // accesses
							   // -- the
							   // warnings
							   // are
							   // bogus
							   // since we
							   // get here
							   // only for
							   // dim==3,
							   // but at
							   // least
							   // one
							   // isn't
							   // quite
							   // smart
							   // enough
							   // to
							   // notice
							   // this and
							   // warns
							   // when
							   // compiling
							   // the
							   // function
							   // in 2d
			  cross_product (tangent, normals[0], normals[dim-2]);
			  break;
		    default:
			  Assert (false, ExcNotImplemented());
		  }

		Assert (std::fabs (tangent.norm()-1) < 1e-12,
			ExcInternalError());
		
		tangential_vectors.push_back (tangent);
	      }

					     // go through the list of
					     // tangents and make sure
					     // that they all roughly
					     // point in the same
					     // direction as the first
					     // one (i.e. have an
					     // angle less than 90
					     // degrees); if they
					     // don't then flip their
					     // sign
	    {
	      const Tensor<1,dim> first_tangent = tangential_vectors.front();
	      typename std::list<Tensor<1,dim> >::iterator
		t = tangential_vectors.begin();
	      ++t;
	      for (; t != tangential_vectors.end(); ++t)
		if (*t * first_tangent < 0)
		  *t *= -1;
	    }
	    
					     // now compute the
					     // average tangent and
					     // normalize it
	    Tensor<1,dim> average_tangent;
	    for (typename std::list<Tensor<1,dim> >::const_iterator
		   t = tangential_vectors.begin();
		 t != tangential_vectors.end();
		 ++t)
	      average_tangent += *t;
	    average_tangent /= average_tangent.norm();

					     // from the tangent
					     // vector we now need to
					     // again reconstruct dim-1
					     // normal directions in
					     // which the vector field
					     // is to be constrained
	    Tensor<1,dim> constraining_normals[dim-1];
	    internal::VectorTools::
	      compute_orthonormal_vectors<dim> (average_tangent,
						constraining_normals);

					     // now all that is left
					     // is that we add the
					     // constraints for these
					     // dim-1 vectors
	    const internal::VectorTools::VectorDoFTuple<dim> &
		dof_indices = same_dof_range[0]->first;
	    for (unsigned int c=0; c<dim-1; ++c)
	      internal::VectorTools::add_constraint (dof_indices,
						     constraining_normals[c],
						     constraints);
	  }
	}
    }
}



namespace internal
{
  namespace VectorTools
  {
    template <int dim, class InVector, class OutVector, class DH, int spacedim>
    static
    void
    do_integrate_difference (const dealii::hp::MappingCollection<dim,spacedim>    &mapping,
                             const DH              &dof,
                             const InVector        &fe_function,
                             const Function<spacedim>   &exact_solution,
                             OutVector             &difference,
                             const dealii::hp::QCollection<dim> &q,
                             const dealii::VectorTools::NormType &norm,
                             const Function<spacedim>   *weight,
                             const double           exponent_1)
    {
                                       // we mark the "exponent" parameter
                                       // to this function "const" since
                                       // it is strictly incoming, but we
                                       // need to set it to something
                                       // different later on, if
                                       // necessary, so have a read-write
                                       // version of it:
      double exponent = exponent_1;
  
      const unsigned int        n_components = dof.get_fe().n_components();
      const bool                fe_is_system = (n_components != 1);

      if (weight!=0)
        {
          Assert ((weight->n_components==1) || (weight->n_components==n_components),
                  ExcDimensionMismatch(weight->n_components, n_components));
        }

      difference.reinit (dof.get_tria().n_active_cells());
  
      switch (norm)
        {
          case dealii::VectorTools::L2_norm:
          case dealii::VectorTools::H1_seminorm:
          case dealii::VectorTools::H1_norm:
                exponent = 2.;
                break;
          case dealii::VectorTools::L1_norm:
                exponent = 1.;
                break;
          default:
                break;
        }
  
      UpdateFlags update_flags = UpdateFlags (update_quadrature_points  |
                                              update_JxW_values);
      switch (norm)
        {
          case dealii::VectorTools::H1_seminorm:
          case dealii::VectorTools::W1p_seminorm:
          case dealii::VectorTools::W1infty_seminorm:
                update_flags |= UpdateFlags (update_gradients);
                break;
          case dealii::VectorTools::H1_norm:
          case dealii::VectorTools::W1p_norm:
          case dealii::VectorTools::W1infty_norm:
                update_flags |= UpdateFlags (update_gradients);
                                                 // no break!
          default:
                update_flags |= UpdateFlags (update_values);
                break;
        }  

      dealii::hp::FECollection<dim,spacedim> fe_collection (dof.get_fe());
      dealii::hp::FEValues<dim,spacedim> x_fe_values(mapping, fe_collection, q, update_flags);

      const unsigned int max_n_q_points = q.max_n_quadrature_points ();
      
      std::vector< dealii::Vector<double> >
        function_values (max_n_q_points, dealii::Vector<double>(n_components));
      std::vector<std::vector<Tensor<1,spacedim> > >
        function_grads (max_n_q_points, std::vector<Tensor<1,spacedim> >(n_components));

      std::vector<double>
        weight_values (max_n_q_points);
      std::vector<dealii::Vector<double> >
        weight_vectors (max_n_q_points, dealii::Vector<double>(n_components));

      std::vector<dealii::Vector<double> >
        psi_values (max_n_q_points, dealii::Vector<double>(n_components));
      std::vector<std::vector<Tensor<1,spacedim> > >
        psi_grads (max_n_q_points, std::vector<Tensor<1,spacedim> >(n_components));
      std::vector<double>
        psi_scalar (max_n_q_points);

                                       // tmp vector when we use the
                                       // Function<spacedim> functions for
                                       // scalar functions
      std::vector<double>         tmp_values (max_n_q_points);
      std::vector<Tensor<1,spacedim> > tmp_gradients (max_n_q_points);
  
                                       // loop over all cells
      typename DH::active_cell_iterator cell = dof.begin_active(),
                                        endc = dof.end();
      for (unsigned int index=0; cell != endc; ++cell, ++index)
        {
          double diff=0;
                                           // initialize for this cell
          x_fe_values.reinit (cell);

          const dealii::FEValues<dim, spacedim> &fe_values  = x_fe_values.get_present_fe_values ();
          const unsigned int   n_q_points = fe_values.n_quadrature_points;
          
                                           // resize all out scratch
                                           // arrays to the number of
                                           // quadrature points we use
                                           // for the present cell
          function_values.resize (n_q_points,
                                  dealii::Vector<double>(n_components));
          function_grads.resize (n_q_points,
                                 std::vector<Tensor<1,spacedim> >(n_components));

          weight_values.resize (n_q_points);
          weight_vectors.resize (n_q_points,
                                 dealii::Vector<double>(n_components));

          psi_values.resize (n_q_points,
                             dealii::Vector<double>(n_components));
          psi_grads.resize (n_q_points,
                            std::vector<Tensor<1,spacedim> >(n_components));
          psi_scalar.resize (n_q_points);

          tmp_values.resize (n_q_points);
          tmp_gradients.resize (n_q_points);

          if (weight!=0)
            {
              if (weight->n_components>1)
                weight->vector_value_list (fe_values.get_quadrature_points(),
                                           weight_vectors);
              else
                {
                  weight->value_list (fe_values.get_quadrature_points(),
                                      weight_values);
                  for (unsigned int k=0;k<n_q_points;++k)
                    weight_vectors[k] = weight_values[k];
                }
            }
          else
            {
              for (unsigned int k=0;k<n_q_points;++k)
                weight_vectors[k] = 1.;
            }
      
      
          if (update_flags & update_values)
            {
                                               // first compute the exact solution
                                               // (vectors) at the quadrature points
                                               // try to do this as efficient as
                                               // possible by avoiding a second
                                               // virtual function call in case
                                               // the function really has only
                                               // one component
              if (fe_is_system)
                exact_solution.vector_value_list (fe_values.get_quadrature_points(),
                                                  psi_values);
              else
                {
                  exact_solution.value_list (fe_values.get_quadrature_points(),
                                             tmp_values);
                  for (unsigned int i=0; i<n_q_points; ++i)
                    psi_values[i](0) = tmp_values[i];
                }
	  
                                               // then subtract finite element
                                               // fe_function
              fe_values.get_function_values (fe_function, function_values);
              for (unsigned int q=0; q<n_q_points; ++q)
                psi_values[q] -= function_values[q];
            }

                                           // Do the same for gradients, if required
          if (update_flags & update_gradients)
            {
                                               // try to be a little clever
                                               // to avoid recursive virtual
                                               // function calls when calling
                                               // gradient_list for functions
                                               // that are really scalar
                                               // functions
              if (fe_is_system)
                exact_solution.vector_gradient_list (fe_values.get_quadrature_points(),
                                                     psi_grads);
              else
                {
                  exact_solution.gradient_list (fe_values.get_quadrature_points(),
                                                tmp_gradients);
                  for (unsigned int i=0; i<n_q_points; ++i)
                    psi_grads[i][0] = tmp_gradients[i];
                }      
	  
                                               // then subtract finite element
                                               // function_grads
              fe_values.get_function_grads (fe_function, function_grads);
              for (unsigned int k=0; k<n_components; ++k)
                for (unsigned int q=0; q<n_q_points; ++q)
                  psi_grads[q][k] -= function_grads[q][k];
            }
      
          switch (norm)
            {
              case dealii::VectorTools::mean:
                    std::fill_n (psi_scalar.begin(), n_q_points, 0.0);
                                                     // Compute values in
                                                     // quadrature points
                    for (unsigned int k=0; k<n_components; ++k)
                      for (unsigned int q=0; q<n_q_points; ++q)
                        psi_scalar[q] += psi_values[q](k)
                                         * weight_vectors[q](k);

                                                     // Integrate
                    diff = std::inner_product (psi_scalar.begin(), psi_scalar.end(),
                                               fe_values.get_JxW_values().begin(),
                                               0.0);
                    break;
              case dealii::VectorTools::Lp_norm:
              case dealii::VectorTools::L1_norm:
              case dealii::VectorTools::W1p_norm:
                    std::fill_n (psi_scalar.begin(), n_q_points, 0.0);
                                                     // Compute values in
                                                     // quadrature points
                    for (unsigned int k=0; k<n_components; ++k)
                      for (unsigned int q=0; q<n_q_points; ++q)
                        psi_scalar[q] += std::pow(psi_values[q](k)*psi_values[q](k),
                                                  exponent/2.)
                                         * weight_vectors[q](k);
	    
                                                     // Integrate
                    diff = std::inner_product (psi_scalar.begin(), psi_scalar.end(),
                                               fe_values.get_JxW_values().begin(),
                                               0.0);
                                                     // Compute the root only,
                                                     // if no derivative
                                                     // values are added later
                    if (!(update_flags & update_gradients))
                      diff = std::pow(diff, 1./exponent);
                    break;
              case dealii::VectorTools::L2_norm:
              case dealii::VectorTools::H1_norm:
                    std::fill_n (psi_scalar.begin(), n_q_points, 0.0);
                                                     // Compute values in
                                                     // quadrature points
                    for (unsigned int k=0; k<n_components; ++k)
                      for (unsigned int q=0; q<n_q_points; ++q)
                        psi_scalar[q] += psi_values[q](k)*psi_values[q](k)
                                         * weight_vectors[q](k);

                                                     // Integrate
                    diff = std::inner_product (psi_scalar.begin(), psi_scalar.end(),
                                               fe_values.get_JxW_values().begin(),
                                               0.0);
                                                     // Compute the root only,
                                                     // if no derivative
                                                     // values are added later
                    if (norm == dealii::VectorTools::L2_norm)
                      diff=std::sqrt(diff);
                    break;
              case dealii::VectorTools::Linfty_norm:
              case dealii::VectorTools::W1infty_norm:
                    std::fill_n (psi_scalar.begin(), n_q_points, 0.0);
                    for (unsigned int k=0; k<n_components; ++k)
                      for (unsigned int q=0; q<n_q_points; ++q)
                        {
                          double newval = std::fabs(psi_values[q](k))
                                          * weight_vectors[q](k);
                          if (psi_scalar[q]<newval)
                            psi_scalar[q] = newval;
                        }
                                                     // Maximum on one cell
                    diff = *std::max_element (psi_scalar.begin(), psi_scalar.end());
                    break;
              case dealii::VectorTools::H1_seminorm:
              case dealii::VectorTools::W1p_seminorm:
              case dealii::VectorTools::W1infty_seminorm:
                    break;
              default:
                    Assert (false, ExcNotImplemented());
                    break;
            }

          switch (norm)
            {
              case dealii::VectorTools::W1p_seminorm:
              case dealii::VectorTools::W1p_norm:
                    std::fill_n (psi_scalar.begin(), n_q_points, 0.0);
                    for (unsigned int k=0; k<n_components; ++k)
                      for (unsigned int q=0; q<n_q_points; ++q)
                        psi_scalar[q] += std::pow(psi_grads[q][k] * psi_grads[q][k],
                                                  exponent/2.)
                                         * weight_vectors[q](k);
	    
                    diff += std::inner_product (psi_scalar.begin(), psi_scalar.end(),
                                                fe_values.get_JxW_values().begin(),
                                                0.0);
                    diff = std::pow(diff, 1./exponent);
                    break;
              case dealii::VectorTools::H1_seminorm:
              case dealii::VectorTools::H1_norm:
                                                     // take square of integrand
                    std::fill_n (psi_scalar.begin(), n_q_points, 0.0);
                    for (unsigned int k=0; k<n_components; ++k)
                      for (unsigned int q=0; q<n_q_points; ++q)
                        psi_scalar[q] += (psi_grads[q][k] * psi_grads[q][k])
                                         * weight_vectors[q](k);

                                                     // add seminorm to L_2 norm or
                                                     // to zero
                    diff += std::inner_product (psi_scalar.begin(), psi_scalar.end(),
                                                fe_values.get_JxW_values().begin(),
                                                0.0);
                    diff = std::sqrt(diff);
                    break;
              case dealii::VectorTools::W1infty_seminorm:
              case dealii::VectorTools::W1infty_norm:
                    Assert(false, ExcNotImplemented());
                    std::fill_n (psi_scalar.begin(), n_q_points, 0.0);
                    for (unsigned int k=0; k<n_components; ++k)
                      for (unsigned int q=0; q<n_q_points; ++q)
                        {
                          double t = 0.;
                          for (unsigned int d=0;d<dim;++d)
                            t = std::max(t,std::fabs(psi_grads[q][k][d])
                                         * weight_vectors[q](k));
		  
                          psi_scalar[q] = std::max(psi_scalar[q],t);
                        }

                    for (unsigned int i=0;i<psi_scalar.size();++i)
                      diff = std::max (diff, psi_scalar[i]);
                    break;
              default:
                    break;
            }
                                           // append result of this cell
                                           // to the end of the vector
          Assert (numbers::is_finite(diff), ExcInternalError());
          difference(index) = diff;
        }
    }

  } //namespace VectorTools
} // namespace internal




template <int dim, class InVector, class OutVector, int spacedim>
void
VectorTools::integrate_difference (const Mapping<dim, spacedim>    &mapping,
				   const DoFHandler<dim,spacedim> &dof,
				   const InVector        &fe_function,
				   const Function<spacedim>   &exact_solution,
				   OutVector             &difference,
				   const Quadrature<dim> &q,
				   const NormType        &norm,
				   const Function<spacedim>   *weight,
				   const double           exponent)
{
  internal::VectorTools
    ::do_integrate_difference (hp::MappingCollection<dim,spacedim>(mapping),
                               dof, fe_function, exact_solution,
                               difference, hp::QCollection<dim>(q),
                               norm, weight, exponent);
}


template <int dim, class InVector, class OutVector, int spacedim>
void
VectorTools::integrate_difference (const DoFHandler<dim,spacedim>    &dof,
				   const InVector           &fe_function,
				   const Function<spacedim>      &exact_solution,
				   OutVector                &difference,
				   const Quadrature<dim>    &q,
				   const NormType           &norm,
				   const Function<spacedim>      *weight,
				   const double              exponent)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  internal::VectorTools
      ::do_integrate_difference(hp::StaticMappingQ1<dim,spacedim>::mapping_collection,
                              dof, fe_function, exact_solution,
                              difference, hp::QCollection<dim>(q),
                              norm, weight, exponent);
}



template <int dim, class InVector, class OutVector, int spacedim>
void
VectorTools::integrate_difference (const dealii::hp::MappingCollection<dim,spacedim>    &mapping,
				   const dealii::hp::DoFHandler<dim,spacedim> &dof,
				   const InVector        &fe_function,
				   const Function<spacedim>   &exact_solution,
				   OutVector             &difference,
				   const dealii::hp::QCollection<dim> &q,
				   const NormType        &norm,
				   const Function<spacedim>   *weight,
				   const double           exponent)
{
  internal::VectorTools
    ::do_integrate_difference (hp::MappingCollection<dim,spacedim>(mapping),
                               dof, fe_function, exact_solution,
                               difference, q,
                               norm, weight, exponent);
}


template <int dim, class InVector, class OutVector, int spacedim>
void
VectorTools::integrate_difference (const dealii::hp::DoFHandler<dim,spacedim>    &dof,
				   const InVector           &fe_function,
				   const Function<spacedim>      &exact_solution,
				   OutVector                &difference,
				   const dealii::hp::QCollection<dim>    &q,
				   const NormType           &norm,
				   const Function<spacedim>      *weight,
				   const double              exponent)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  internal::VectorTools
    ::do_integrate_difference(hp::StaticMappingQ1<dim>::mapping_collection,
                              dof, fe_function, exact_solution,
                              difference, q,
                              norm, weight, exponent);
}



template <int dim, class InVector, int spacedim>
void
VectorTools::point_difference (const DoFHandler<dim,spacedim> &dof,
			       const InVector        &fe_function,
			       const Function<spacedim>   &exact_function,
			       Vector<double>        &difference,
			       const Point<spacedim>      &point)
{
   point_difference(StaticMappingQ1<dim>::mapping,
                    dof,
                    fe_function,
                    exact_function,
                    difference,
                    point);
}


template <int dim, class InVector, int spacedim>
void
VectorTools::point_difference (const Mapping<dim, spacedim>    &mapping,
                               const DoFHandler<dim,spacedim> &dof,
			       const InVector        &fe_function,
			       const Function<spacedim>   &exact_function,
			       Vector<double>        &difference,
			       const Point<spacedim>      &point)
{
  const FiniteElement<dim>& fe = dof.get_fe();

  Assert(difference.size() == fe.n_components(),
	 ExcDimensionMismatch(difference.size(), fe.n_components()));

                                   // first find the cell in which this point
                                   // is, initialize a quadrature rule with
                                   // it, and then a FEValues object
  const std::pair<typename DoFHandler<dim,spacedim>::active_cell_iterator, Point<spacedim> >
    cell_point = GridTools::find_active_cell_around_point (mapping, dof, point);

  Assert(GeometryInfo<dim>::distance_to_unit_cell(cell_point.second) < 1e-10,
         ExcInternalError());
  
  const Quadrature<dim>
    quadrature (GeometryInfo<dim>::project_to_unit_cell(cell_point.second));
  FEValues<dim> fe_values(mapping, fe, quadrature, update_values);
  fe_values.reinit(cell_point.first);

                                   // then use this to get at the values of
                                   // the given fe_function at this point
  std::vector<Vector<double> > u_value(1, Vector<double> (fe.n_components()));
  fe_values.get_function_values(fe_function, u_value);

  if (fe.n_components() == 1)
    difference(0) = exact_function.value(point);
  else
    exact_function.vector_value(point, difference);
    
  for (unsigned int i=0; i<difference.size(); ++i)
    difference(i) -= u_value[0](i);
}


template <int dim, class InVector, int spacedim>
void
VectorTools::point_value (const DoFHandler<dim,spacedim> &dof,
			  const InVector        &fe_function,
			  const Point<spacedim>      &point,
			  Vector<double>        &value)
{

  point_value (StaticMappingQ1<dim,spacedim>::mapping,
               dof,
               fe_function,
               point,
               value);
}



template <int dim, class InVector, int spacedim>
double
VectorTools::point_value (const DoFHandler<dim,spacedim> &dof,
			  const InVector        &fe_function,
			  const Point<spacedim>      &point)
{
  return point_value (StaticMappingQ1<dim,spacedim>::mapping,
                      dof,
                      fe_function,
                      point);
}

template <int dim, class InVector, int spacedim>
void
VectorTools::point_value (const Mapping<dim, spacedim>    &mapping,
                          const DoFHandler<dim,spacedim> &dof,
			  const InVector        &fe_function,
			  const Point<spacedim>      &point,
			  Vector<double>        &value)
{
  const FiniteElement<dim>& fe = dof.get_fe();

  Assert(value.size() == fe.n_components(),
	 ExcDimensionMismatch(value.size(), fe.n_components()));

                                   // first find the cell in which this point
                                   // is, initialize a quadrature rule with
                                   // it, and then a FEValues object
  const std::pair<typename DoFHandler<dim,spacedim>::active_cell_iterator, Point<spacedim> >
    cell_point
    = GridTools::find_active_cell_around_point (mapping, dof, point);

  Assert(GeometryInfo<dim>::distance_to_unit_cell(cell_point.second) < 1e-10,
         ExcInternalError());

  const Quadrature<dim>
    quadrature (GeometryInfo<dim>::project_to_unit_cell(cell_point.second));

  FEValues<dim> fe_values(mapping, fe, quadrature, update_values);
  fe_values.reinit(cell_point.first);

                                   // then use this to get at the values of
                                   // the given fe_function at this point
  std::vector<Vector<double> > u_value(1, Vector<double> (fe.n_components()));
  fe_values.get_function_values(fe_function, u_value);

  value = u_value[0];
}



template <int dim, class InVector, int spacedim>
double
VectorTools::point_value (const Mapping<dim, spacedim>    &mapping,
                          const DoFHandler<dim,spacedim> &dof,
			  const InVector        &fe_function,
			  const Point<spacedim>      &point)
{
  const FiniteElement<dim>& fe = dof.get_fe();

  Assert(fe.n_components() == 1,
	 ExcMessage ("Finite element is not scalar as is necessary for this function"));

                                   // first find the cell in which this point
                                   // is, initialize a quadrature rule with
                                   // it, and then a FEValues object
  const std::pair<typename DoFHandler<dim,spacedim>::active_cell_iterator, Point<spacedim> >
    cell_point = GridTools::find_active_cell_around_point (mapping, dof, point);

  Assert(GeometryInfo<dim>::distance_to_unit_cell(cell_point.second) < 1e-10,
         ExcInternalError());
  
  const Quadrature<dim>
    quadrature (GeometryInfo<dim>::project_to_unit_cell(cell_point.second));
  FEValues<dim> fe_values(mapping, fe, quadrature, update_values);
  fe_values.reinit(cell_point.first);

                                   // then use this to get at the values of
                                   // the given fe_function at this point
  std::vector<double> u_value(1);
  fe_values.get_function_values(fe_function, u_value);

  return u_value[0];
}


template <int dim, class InVector, int spacedim>
double
VectorTools::compute_mean_value (const Mapping<dim, spacedim>    &mapping,
				 const DoFHandler<dim,spacedim> &dof,
				 const Quadrature<dim> &quadrature,
				 const InVector        &v,
				 const unsigned int     component)
{
  Assert (v.size() == dof.n_dofs(),
	  ExcDimensionMismatch (v.size(), dof.n_dofs()));
  Assert (component < dof.get_fe().n_components(),
	  ExcIndexRange(component, 0, dof.get_fe().n_components()));
  
  FEValues<dim> fe(mapping, dof.get_fe(), quadrature,
		   UpdateFlags(update_JxW_values
			       | update_values));

  typename DoFHandler<dim,spacedim>::active_cell_iterator c;
  std::vector<Vector<double> > values(quadrature.size(),
				      Vector<double> (dof.get_fe().n_components()));
  
  double mean = 0.;
  double area = 0.;
				   // Compute mean value
  for (c = dof.begin_active(); c != dof.end(); ++c)
    {
      fe.reinit (c);
      fe.get_function_values(v, values);
      for (unsigned int k=0; k< quadrature.size(); ++k)
	{
	  mean += fe.JxW(k) * values[k](component);
	  area += fe.JxW(k);
	};
    };
  
  return (mean/area);
}


template <int dim, class InVector, int spacedim>
double
VectorTools::compute_mean_value (const DoFHandler<dim,spacedim> &dof,
				 const Quadrature<dim> &quadrature,
				 const InVector        &v,
				 const unsigned int     component)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  return compute_mean_value(StaticMappingQ1<dim>::mapping, dof, quadrature, v, component);
}

DEAL_II_NAMESPACE_CLOSE

#endif
