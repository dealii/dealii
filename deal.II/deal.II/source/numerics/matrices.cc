//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/function.h>
#include <base/quadrature.h>
#include <base/thread_management.h>
#include <base/work_stream.h>
#include <base/multithread_info.h>
#include <grid/tria_iterator.h>
#include <grid/geometry_info.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <grid/tria_iterator.h>
#include <base/geometry_info.h>
#include <base/quadrature.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <hp/fe_values.h>
#include <hp/mapping_collection.h>
#include <numerics/matrices.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <fe/mapping_q1.h>

#include <algorithm>
#include <set>
#include <cmath>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MatrixCreator
  {
    namespace AssemblerData
    {
      template <int dim,
		int spacedim>
      struct Scratch
      {
	  Scratch (const FiniteElement<dim,spacedim> &fe,
		   const UpdateFlags         update_flags,
		   const Function<spacedim>      *coefficient,
		   const Function<spacedim>      *rhs_function,
		   const Quadrature<dim>    &quadrature,
		   const Mapping<dim,spacedim>       &mapping)
			  :
			  fe_collection (fe),
			  quadrature_collection (quadrature),
			  mapping_collection (mapping),
			  x_fe_values (mapping_collection,
				       fe_collection,
				       quadrature_collection,
				       update_flags),
			  coefficient_values(quadrature_collection.max_n_quadrature_points()),
			  coefficient_vector_values (quadrature_collection.max_n_quadrature_points(),
						     dealii::Vector<double> (fe_collection.n_components())),
			  rhs_values(quadrature_collection.max_n_quadrature_points()),
			  rhs_vector_values (quadrature_collection.max_n_quadrature_points(),
					     dealii::Vector<double> (fe_collection.n_components())),
			  coefficient (coefficient),
			  rhs_function (rhs_function),
			  update_flags (update_flags)
	    {}

	  Scratch (const ::dealii::hp::FECollection<dim,spacedim> &fe,
		   const UpdateFlags         update_flags,
		   const Function<spacedim>      *coefficient,
		   const Function<spacedim>      *rhs_function,
		   const ::dealii::hp::QCollection<dim>    &quadrature,
		   const ::dealii::hp::MappingCollection<dim,spacedim>       &mapping)
			  :
			  fe_collection (fe),
			  quadrature_collection (quadrature),
			  mapping_collection (mapping),
			  x_fe_values (mapping_collection,
				       fe_collection,
				       quadrature_collection,
				       update_flags),
			  coefficient_values(quadrature_collection.max_n_quadrature_points()),
			  coefficient_vector_values (quadrature_collection.max_n_quadrature_points(),
						     dealii::Vector<double> (fe_collection.n_components())),
			  rhs_values(quadrature_collection.max_n_quadrature_points()),
			  rhs_vector_values(quadrature_collection.max_n_quadrature_points(),
					    dealii::Vector<double> (fe_collection.n_components())),
			  coefficient (coefficient),
			  rhs_function (rhs_function),
			  update_flags (update_flags)
	    {}
	
	  Scratch (const Scratch &data)
			  :
			  fe_collection (data.fe_collection),
			  quadrature_collection (data.quadrature_collection),
			  mapping_collection (data.mapping_collection),
			  x_fe_values (mapping_collection,
				       fe_collection,
				       quadrature_collection,
				       data.update_flags),
			  coefficient_values (data.coefficient_values),
			  coefficient_vector_values (data.coefficient_vector_values),
			  rhs_values (data.rhs_values),
			  rhs_vector_values (data.rhs_vector_values),
			  coefficient (data.coefficient),
			  rhs_function (data.rhs_function),
			  update_flags (data.update_flags)
	    {}
	
	  const ::dealii::hp::FECollection<dim,spacedim> fe_collection;
	  const ::dealii::hp::QCollection<dim> quadrature_collection;
	  const ::dealii::hp::MappingCollection<dim,spacedim> mapping_collection;
	
	  ::dealii::hp::FEValues<dim,spacedim> x_fe_values;
	
	  std::vector<double>                  coefficient_values;
	  std::vector<dealii::Vector<double> > coefficient_vector_values;
	  std::vector<double>                  rhs_values;
	  std::vector<dealii::Vector<double> > rhs_vector_values;

	  const Function<spacedim>   *coefficient;
	  const Function<spacedim>   *rhs_function;

	  const UpdateFlags update_flags;
      };


      struct CopyData
      {
	  std::vector<unsigned int> dof_indices;
	  FullMatrix<double>        cell_matrix;
	  dealii::Vector<double>    cell_rhs;
      };
      

      
    }
    

    template <int dim,
	      int spacedim,
	      typename CellIterator>
    void mass_assembler (const CellIterator &cell,
			 internal::MatrixCreator::AssemblerData::Scratch<dim,spacedim> &data,
			 internal::MatrixCreator::AssemblerData::CopyData              &copy_data)
    {
      data.x_fe_values.reinit (cell);
      const FEValues<dim,spacedim> &fe_values = data.x_fe_values.get_present_fe_values ();
      
      const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
			 n_q_points    = fe_values.n_quadrature_points;
      const FiniteElement<dim,spacedim>    &fe  = fe_values.get_fe();
      const unsigned int n_components  = fe.n_components();

      Assert(data.rhs_function == 0 ||
	     data.rhs_function->n_components==1 ||
	     data.rhs_function->n_components==n_components,
	     ::dealii::MatrixCreator::ExcComponentMismatch());
      Assert(data.coefficient == 0 ||
	     data.coefficient->n_components==1 ||
	     data.coefficient->n_components==n_components,
	     ::dealii::MatrixCreator::ExcComponentMismatch());

      copy_data.cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
      copy_data.cell_matrix = 0;

      copy_data.cell_rhs.reinit (dofs_per_cell);
      copy_data.cell_rhs = 0;
      
      copy_data.dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (copy_data.dof_indices);

      if (data.rhs_function != 0)
	{
	  if (data.rhs_function->n_components==1)
	    {
	      data.rhs_values.resize (n_q_points);
	      data.rhs_function->value_list (fe_values.get_quadrature_points(),
					     data.rhs_values);
	    }
	  else
	    {
	      data.rhs_vector_values.resize (n_q_points, 
					     dealii::Vector<double>(n_components));
	      data.rhs_function->vector_value_list (fe_values.get_quadrature_points(),
						    data.rhs_vector_values);
	    }
	}

      if (data.coefficient != 0)
	{
	  if (data.coefficient->n_components==1)
	    {
	      data.coefficient_values.resize (n_q_points);
 	      data.coefficient->value_list (fe_values.get_quadrature_points(),
					    data.coefficient_values);
	    }
	  else
	    {
	      data.coefficient_vector_values.resize (n_q_points, 
						     dealii::Vector<double>(n_components));
	      data.coefficient->vector_value_list (fe_values.get_quadrature_points(),
						   data.coefficient_vector_values);
	    }
	}


      if (data.coefficient != 0)
	{
	  if (data.coefficient->n_components==1)
	    {
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
		  const unsigned int component_i =
		    fe.system_to_component_index(i).first;
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    if ((n_components==1) ||
			(fe.system_to_component_index(j).first ==
			 component_i))
		      for (unsigned int point=0; point<n_q_points; ++point)
			copy_data.cell_matrix(i,j)
			  += (fe_values.shape_value(i,point) * 
			      fe_values.shape_value(j,point) * 
			      fe_values.JxW(point) *
			      data.coefficient_values[point]);

		  if (data.rhs_function != 0)
		    {
		      if (data.rhs_function->n_components==1)
			for (unsigned int point=0; point<n_q_points; ++point)
			  copy_data.cell_rhs(i) += fe_values.shape_value(i, point) *
			    data.rhs_values[point] * fe_values.JxW(point);
		      else
			for (unsigned int point=0; point<n_q_points; ++point)
			  copy_data.cell_rhs(i) += fe_values.shape_value(i, point) *
			    data.rhs_vector_values[point](component_i) * 
			    fe_values.JxW(point);
		    }
		}
	    }
	  else
	    {
	      if (fe.is_primitive ())
		{
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    {
		      const unsigned int component_i =
			fe.system_to_component_index(i).first;
		      for (unsigned int j=0; j<dofs_per_cell; ++j)
			if ((n_components==1) ||
			    (fe.system_to_component_index(j).first ==
			     component_i))
			  for (unsigned int point=0; point<n_q_points; ++point)
			    copy_data.cell_matrix(i,j) +=
			      (fe_values.shape_value(i,point) * 
			       fe_values.shape_value(j,point) * 
			       fe_values.JxW(point) *
			       data.coefficient_vector_values[point](component_i));

		      if (data.rhs_function != 0)
			{
			  if (data.rhs_function->n_components==1)
			    for (unsigned int point=0; point<n_q_points; ++point)
			      copy_data.cell_rhs(i) += fe_values.shape_value(i, point) *
				data.rhs_values[point] * fe_values.JxW(point);
			  else
			    for (unsigned int point=0; point<n_q_points; ++point)
			      copy_data.cell_rhs(i) += fe_values.shape_value(i, point) *
				data.rhs_vector_values[point](component_i) * 
				fe_values.JxW(point);
			}
		    }
		}
	      else
						 // non-primitive vector-valued FE
		{
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
		      if (fe.get_nonzero_components(i)[comp_i])
			{
			  for (unsigned int j=0; j<dofs_per_cell; ++j)
			    for (unsigned int comp_j = 0; comp_j < n_components; ++comp_j)
			      if (fe.get_nonzero_components(j)[comp_j])
				if (comp_i == comp_j)
				  for (unsigned int point=0; point<n_q_points; ++point)
				    copy_data.cell_matrix(i,j) += 
				      (fe_values.shape_value_component(i,point,comp_i) * 
				       fe_values.shape_value_component(j,point,comp_j) * 
				       fe_values.JxW(point) *
				       data.coefficient_vector_values[point](comp_i));
			      
			  if (data.rhs_function != 0)
			    {
			      if (data.rhs_function->n_components==1)
				for (unsigned int point=0; point<n_q_points; ++point)
				  copy_data.cell_rhs(i) += 
				    fe_values.shape_value_component(i,point,comp_i) *
				    data.rhs_values[point] * fe_values.JxW(point);
			      else
				for (unsigned int point=0; point<n_q_points; ++point)
				  copy_data.cell_rhs(i) += 
				    fe_values.shape_value_component(i,point,comp_i) *
				    data.rhs_vector_values[point](comp_i) * 
				    fe_values.JxW(point);
			    }
			}
		}
	    }
	}
      else
					 // no coefficient
	{
	  if (fe.is_primitive ())
	    {
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
		  const unsigned int component_i =
		    fe.system_to_component_index(i).first;
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    if ((n_components==1) ||
			(fe.system_to_component_index(j).first ==
			 component_i))
		      for (unsigned int point=0; point<n_q_points; ++point)
			copy_data.cell_matrix(i,j) +=
			  (fe_values.shape_value(i,point) * 
			   fe_values.shape_value(j,point) * 
			   fe_values.JxW(point));

		  if (data.rhs_function != 0)
		    {
		      if (data.rhs_function->n_components==1)
			for (unsigned int point=0; point<n_q_points; ++point)
			  copy_data.cell_rhs(i) += fe_values.shape_value(i, point) *
			    data.rhs_values[point] * fe_values.JxW(point);
		      else
			for (unsigned int point=0; point<n_q_points; ++point)
			  copy_data.cell_rhs(i) += fe_values.shape_value(i, point) *
			    data.rhs_vector_values[point](component_i) * 
			    fe_values.JxW(point);
		    }
		}
	    }
	  else
					     // non-primitive FE, no coefficient
	    {
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
		  if (fe.get_nonzero_components(i)[comp_i])
		    {
		      for (unsigned int j=0; j<dofs_per_cell; ++j)
			for (unsigned int comp_j = 0; comp_j < n_components; ++comp_j)
			  if (fe.get_nonzero_components(j)[comp_j])
			    if (comp_i == comp_j)
			      for (unsigned int point=0; point<n_q_points; ++point)
				copy_data.cell_matrix(i,j) += 
				  (fe_values.shape_value_component(i,point,comp_i) * 
				   fe_values.shape_value_component(j,point,comp_j) * 
				   fe_values.JxW(point));

		      if (data.rhs_function != 0)
			{
			  if (data.rhs_function->n_components==1)
			    for (unsigned int point=0; point<n_q_points; ++point)
			      copy_data.cell_rhs(i) += 
				fe_values.shape_value_component(i,point,comp_i) *
				data.rhs_values[point] * fe_values.JxW(point);
			  else
			    for (unsigned int point=0; point<n_q_points; ++point)
			      copy_data.cell_rhs(i) += 
				fe_values.shape_value_component(i,point,comp_i) *
				data.rhs_vector_values[point](comp_i) * 
				fe_values.JxW(point);
			}
		    }
	    }
	}
    }



    template <int dim,
	      int spacedim,
	      typename CellIterator>
    void laplace_assembler (const CellIterator &cell,
			    internal::MatrixCreator::AssemblerData::Scratch<dim,spacedim> &data,
			    internal::MatrixCreator::AssemblerData::CopyData              &copy_data)
    {
      data.x_fe_values.reinit (cell);
      const FEValues<dim,spacedim> &fe_values = data.x_fe_values.get_present_fe_values ();      

      const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
			 n_q_points    = fe_values.n_quadrature_points;
      const FiniteElement<dim,spacedim>    &fe  = fe_values.get_fe();
      const unsigned int n_components  = fe.n_components();

      Assert(data.rhs_function == 0 ||
	     data.rhs_function->n_components==1 ||
	     data.rhs_function->n_components==n_components,
	     ::dealii::MatrixCreator::ExcComponentMismatch());
      Assert(data.coefficient == 0 ||
	     data.coefficient->n_components==1 ||
	     data.coefficient->n_components==n_components,
	     ::dealii::MatrixCreator::ExcComponentMismatch());

      copy_data.cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
      copy_data.cell_matrix = 0;

      copy_data.cell_rhs.reinit (dofs_per_cell);
      copy_data.cell_rhs = 0;
      
      copy_data.dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (copy_data.dof_indices);


      if (data.rhs_function != 0)
	{
	  if (data.rhs_function->n_components==1)
	    {
	      data.rhs_values.resize (n_q_points);
	      data.rhs_function->value_list (fe_values.get_quadrature_points(),
					     data.rhs_values);
	    }
	  else
	    {
	      data.rhs_vector_values.resize (n_q_points, 
					     dealii::Vector<double>(n_components));
	      data.rhs_function->vector_value_list (fe_values.get_quadrature_points(),
						    data.rhs_vector_values);
	    }
	}

      if (data.coefficient != 0)
	{
	  if (data.coefficient->n_components==1)
	    {
	      data.coefficient_values.resize (n_q_points);
 	      data.coefficient->value_list (fe_values.get_quadrature_points(),
					    data.coefficient_values);
	    }
	  else
	    {
	      data.coefficient_vector_values.resize (n_q_points, 
						     dealii::Vector<double>(n_components));
	      data.coefficient->vector_value_list (fe_values.get_quadrature_points(),
						   data.coefficient_vector_values);
	    }
	}

      
      if (data.coefficient != 0)
	{
	  if (data.coefficient->n_components==1)
	    {
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
		  const unsigned int component_i =
		    fe.system_to_component_index(i).first;
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    if ((n_components==1) ||
			(fe.system_to_component_index(j).first ==
			 component_i))
		      for (unsigned int point=0; point<n_q_points; ++point)
			copy_data.cell_matrix(i,j)
			  += (fe_values.shape_grad(i,point) * 
			      fe_values.shape_grad(j,point) * 
			      fe_values.JxW(point) *
			      data.coefficient_values[point]);

		  if (data.rhs_function != 0)
		    {
		      if (data.rhs_function->n_components==1)
			for (unsigned int point=0; point<n_q_points; ++point)
			  copy_data.cell_rhs(i) += fe_values.shape_value(i, point) *
			    data.rhs_values[point] * fe_values.JxW(point);
		      else
			for (unsigned int point=0; point<n_q_points; ++point)
			  copy_data.cell_rhs(i) += fe_values.shape_value(i, point) *
			    data.rhs_vector_values[point](component_i) * 
			    fe_values.JxW(point);
		    }
		}
	    }
	  else
	    {
	      if (fe.is_primitive ())
		{
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    {
		      const unsigned int component_i =
			fe.system_to_component_index(i).first;
		      for (unsigned int j=0; j<dofs_per_cell; ++j)
			if ((n_components==1) ||
			    (fe.system_to_component_index(j).first ==
			     component_i))
			  for (unsigned int point=0; point<n_q_points; ++point)
			    copy_data.cell_matrix(i,j) +=
			      (fe_values.shape_grad(i,point) * 
			       fe_values.shape_grad(j,point) * 
			       fe_values.JxW(point) *
			       data.coefficient_vector_values[point](component_i));

		      if (data.rhs_function != 0)
			{
			  if (data.rhs_function->n_components==1)
			    for (unsigned int point=0; point<n_q_points; ++point)
			      copy_data.cell_rhs(i) += fe_values.shape_value(i, point) *
				data.rhs_values[point] * fe_values.JxW(point);
			  else
			    for (unsigned int point=0; point<n_q_points; ++point)
			      copy_data.cell_rhs(i) += fe_values.shape_value(i, point) *
				data.rhs_vector_values[point](component_i) * 
				fe_values.JxW(point);
			}
		    }
		}
	      else
						 // non-primitive vector-valued FE
		{
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
		      if (fe.get_nonzero_components(i)[comp_i])
			{
			  for (unsigned int j=0; j<dofs_per_cell; ++j)
			    for (unsigned int comp_j = 0; comp_j < n_components; ++comp_j)
			      if (fe.get_nonzero_components(j)[comp_j])
				if (comp_i == comp_j)
				  for (unsigned int point=0; point<n_q_points; ++point)
				    copy_data.cell_matrix(i,j) += 
				      (fe_values.shape_grad_component(i,point,comp_i) * 
				       fe_values.shape_grad_component(j,point,comp_j) * 
				       fe_values.JxW(point) *
				       data.coefficient_vector_values[point](comp_i));
			      
			  if (data.rhs_function != 0)
			    {
			      if (data.rhs_function->n_components==1)
				for (unsigned int point=0; point<n_q_points; ++point)
				  copy_data.cell_rhs(i) += 
				    fe_values.shape_value_component(i,point,comp_i) *
				    data.rhs_values[point] * fe_values.JxW(point);
			      else
				for (unsigned int point=0; point<n_q_points; ++point)
				  copy_data.cell_rhs(i) += 
				    fe_values.shape_value_component(i,point,comp_i) *
				    data.rhs_vector_values[point](comp_i) * 
				    fe_values.JxW(point);
			    }
			}
		}
	    }
	}
      else
					 // no coefficient
	{
	  if (fe.is_primitive ())
	    {
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
		  const unsigned int component_i =
		    fe.system_to_component_index(i).first;
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    if ((n_components==1) ||
			(fe.system_to_component_index(j).first ==
			 component_i))
		      for (unsigned int point=0; point<n_q_points; ++point)
			copy_data.cell_matrix(i,j) +=
			  (fe_values.shape_grad(i,point) * 
			   fe_values.shape_grad(j,point) * 
			   fe_values.JxW(point));

		  if (data.rhs_function != 0)
		    {
		      if (data.rhs_function->n_components==1)
			for (unsigned int point=0; point<n_q_points; ++point)
			  copy_data.cell_rhs(i) += fe_values.shape_value(i, point) *
			    data.rhs_values[point] * fe_values.JxW(point);
		      else
			for (unsigned int point=0; point<n_q_points; ++point)
			  copy_data.cell_rhs(i) += fe_values.shape_value(i, point) *
			    data.rhs_vector_values[point](component_i) * 
			    fe_values.JxW(point);
		    }
		}
	    }
	  else
					     // non-primitive FE, no coefficient
	    {
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
		  if (fe.get_nonzero_components(i)[comp_i])
		    {
		      for (unsigned int j=0; j<dofs_per_cell; ++j)
			for (unsigned int comp_j = 0; comp_j < n_components; ++comp_j)
			  if (fe.get_nonzero_components(j)[comp_j])
			    if (comp_i == comp_j)
			      for (unsigned int point=0; point<n_q_points; ++point)
				copy_data.cell_matrix(i,j) += 
				  (fe_values.shape_grad_component(i,point,comp_i) * 
				   fe_values.shape_grad_component(j,point,comp_j) * 
				   fe_values.JxW(point));

		      if (data.rhs_function != 0)
			{
			  if (data.rhs_function->n_components==1)
			    for (unsigned int point=0; point<n_q_points; ++point)
			      copy_data.cell_rhs(i) += 
				fe_values.shape_value_component(i,point,comp_i) *
				data.rhs_values[point] * fe_values.JxW(point);
			  else
			    for (unsigned int point=0; point<n_q_points; ++point)
			      copy_data.cell_rhs(i) += 
				fe_values.shape_value_component(i,point,comp_i) *
				data.rhs_vector_values[point](comp_i) * 
				fe_values.JxW(point);
			}
		    }
	    }
	}
    }
    

    
    template <typename MatrixType,
	      typename VectorType>
    void copy_local_to_global (const AssemblerData::CopyData &data,
			       MatrixType *matrix,
			       VectorType *right_hand_side)
    {
      const unsigned int dofs_per_cell = data.dof_indices.size();

      Assert (data.cell_matrix.m() == dofs_per_cell,
	      ExcInternalError());
      Assert (data.cell_matrix.n() == dofs_per_cell,
	      ExcInternalError());
      Assert ((right_hand_side == 0)
	      ||
	      (data.cell_rhs.size() == dofs_per_cell),
	      ExcInternalError());

      matrix->add (data.dof_indices, data.cell_matrix);

      if (right_hand_side != 0)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  (*right_hand_side)(data.dof_indices[i]) += data.cell_rhs(i);
    }
  }
}



template <typename DH>
inline
MatrixCreator::IteratorRange<DH>::
IteratorRange (const active_cell_iterator &first,
	       const active_cell_iterator &second)
		:
		first (first),
		second (second)
{}



template <typename DH>
inline
MatrixCreator::IteratorRange<DH>::IteratorRange (const iterator_pair &ip)
		:
		first (ip.first),
		second (ip.second)
{}




template <int dim, typename number, int spacedim>
void MatrixCreator::create_mass_matrix (const Mapping<dim,spacedim>       &mapping,
					const DoFHandler<dim,spacedim>    &dof,
					const Quadrature<dim>    &q,
					SparseMatrix<number>     &matrix,
					const Function<spacedim> * const coefficient)
{
  Assert (matrix.m() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
  Assert (matrix.n() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

  internal::MatrixCreator::AssemblerData::Scratch<dim, spacedim>
    assembler_data (dof.get_fe(),
		    update_values | update_JxW_values |
		    (coefficient != 0 ? update_quadrature_points : UpdateFlags(0)),
		    coefficient, /*rhs_function=*/0,
		    q, mapping);
  
  internal::MatrixCreator::AssemblerData::CopyData copy_data;
  copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
			 assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());

  WorkStream::run (dof.begin_active(),
		   static_cast<typename DoFHandler<dim,spacedim>::active_cell_iterator>(dof.end()),
		   &internal::MatrixCreator::mass_assembler<dim, spacedim, typename DoFHandler<dim,spacedim>::active_cell_iterator>,
		   std_cxx1x::bind (&internal::MatrixCreator::
				    copy_local_to_global<SparseMatrix<number>, Vector<double> >,
				    _1, &matrix, (Vector<double>*)0),
		   assembler_data,
		   copy_data);
}



template <int dim, typename number, int spacedim>
void MatrixCreator::create_mass_matrix (const DoFHandler<dim,spacedim>    &dof,
					const Quadrature<dim>    &q,
					SparseMatrix<number>     &matrix,
					const Function<spacedim> * const coefficient)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_mass_matrix(StaticMappingQ1<dim,spacedim>::mapping, dof,
		     q, matrix, coefficient);
}



template <int dim, typename number, int spacedim>
void MatrixCreator::create_mass_matrix (const Mapping<dim,spacedim>       &mapping,
					const DoFHandler<dim,spacedim>    &dof,
					const Quadrature<dim>    &q,
					SparseMatrix<number>     &matrix,
					const Function<spacedim>      &rhs,
					Vector<double>           &rhs_vector,
					const Function<spacedim> * const coefficient)
{
  Assert (matrix.m() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
  Assert (matrix.n() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

  internal::MatrixCreator::AssemblerData::Scratch<dim, spacedim>
    assembler_data (dof.get_fe(),
		    update_values |
		    update_JxW_values | update_quadrature_points,
		    coefficient, &rhs,
		    q, mapping);
  internal::MatrixCreator::AssemblerData::CopyData copy_data;
  copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
			 assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());  

  WorkStream::run (dof.begin_active(),
		   static_cast<typename DoFHandler<dim,spacedim>::active_cell_iterator>(dof.end()),
		   &internal::MatrixCreator::mass_assembler<dim, spacedim, typename DoFHandler<dim,spacedim>::active_cell_iterator>,
		   std_cxx1x::bind(&internal::MatrixCreator::
				   copy_local_to_global<SparseMatrix<number>, Vector<double> >,
				   _1, &matrix, &rhs_vector),
		   assembler_data,
		   copy_data);
}



template <int dim, typename number, int spacedim>
void MatrixCreator::create_mass_matrix (const DoFHandler<dim,spacedim>    &dof,
					const Quadrature<dim>    &q,
					SparseMatrix<number>     &matrix,
					const Function<spacedim>      &rhs,
					Vector<double>           &rhs_vector,
					const Function<spacedim> * const coefficient)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_mass_matrix(StaticMappingQ1<dim,spacedim>::mapping,
		     dof, q, matrix, rhs, rhs_vector, coefficient);
}



template <int dim, typename number, int spacedim>
void MatrixCreator::create_mass_matrix (const hp::MappingCollection<dim,spacedim> &mapping,
					const hp::DoFHandler<dim,spacedim>    &dof,
					const hp::QCollection<dim>    &q,
					SparseMatrix<number>     &matrix,
					const Function<spacedim> * const coefficient)
{
  Assert (matrix.m() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
  Assert (matrix.n() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

  internal::MatrixCreator::AssemblerData::Scratch<dim, spacedim>
    assembler_data (dof.get_fe(),
		    update_values | update_JxW_values |
		    (coefficient != 0 ? update_quadrature_points : UpdateFlags(0)),
		    coefficient, /*rhs_function=*/0,
		    q, mapping);
  internal::MatrixCreator::AssemblerData::CopyData copy_data;
  copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
			 assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());
  
  WorkStream::run (dof.begin_active(),
		   static_cast<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator>(dof.end()),
		   &internal::MatrixCreator::mass_assembler<dim, spacedim, typename hp::DoFHandler<dim,spacedim>::active_cell_iterator>,
		   std_cxx1x::bind (&internal::MatrixCreator::
				    copy_local_to_global<SparseMatrix<number>, Vector<double> >,
				    _1, &matrix, (Vector<double>*)0),
		   assembler_data,
		   copy_data);
}



template <int dim, typename number, int spacedim>
void MatrixCreator::create_mass_matrix (const hp::DoFHandler<dim,spacedim>    &dof,
					const hp::QCollection<dim>    &q,
					SparseMatrix<number>     &matrix,
					const Function<spacedim> * const coefficient)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_mass_matrix(hp::StaticMappingQ1<dim,spacedim>::mapping_collection, dof, q, matrix, coefficient);
}



template <int dim, typename number, int spacedim>
void MatrixCreator::create_mass_matrix (const hp::MappingCollection<dim,spacedim>       &mapping,
					const hp::DoFHandler<dim,spacedim>    &dof,
					const hp::QCollection<dim>    &q,
					SparseMatrix<number>     &matrix,
					const Function<spacedim>      &rhs,
					Vector<double>           &rhs_vector,
					const Function<spacedim> * const coefficient)
{
  Assert (matrix.m() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
  Assert (matrix.n() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

  internal::MatrixCreator::AssemblerData::Scratch<dim, spacedim>
    assembler_data (dof.get_fe(),
		    update_values |
		    update_JxW_values | update_quadrature_points,
		    coefficient, &rhs,
		    q, mapping);
  internal::MatrixCreator::AssemblerData::CopyData copy_data;
  copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
			 assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());

  WorkStream::run (dof.begin_active(),
		   static_cast<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator>(dof.end()),
		   &internal::MatrixCreator::mass_assembler<dim, spacedim, typename hp::DoFHandler<dim,spacedim>::active_cell_iterator>,
		   std_cxx1x::bind (&internal::MatrixCreator::
				    copy_local_to_global<SparseMatrix<number>, Vector<double> >,
				    _1, &matrix, &rhs_vector),
		   assembler_data,
		   copy_data);
}



template <int dim, typename number, int spacedim>
void MatrixCreator::create_mass_matrix (const hp::DoFHandler<dim,spacedim>    &dof,
					const hp::QCollection<dim>    &q,
					SparseMatrix<number>     &matrix,
					const Function<spacedim>      &rhs,
					Vector<double>           &rhs_vector,
					const Function<spacedim> * const coefficient)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_mass_matrix(hp::StaticMappingQ1<dim,spacedim>::mapping_collection, dof, q,
		     matrix, rhs, rhs_vector, coefficient);
}




#if deal_II_dimension == 1

template <>
void MatrixCreator::create_boundary_mass_matrix (const Mapping<1,1>&,
						 const DoFHandler<1,1>&,
						 const Quadrature<0>&,
						 SparseMatrix<double>&,
						 const FunctionMap<1>::type&,
						 Vector<double>&,
						 std::vector<unsigned int>&,
						 const Function<1>* const,
						 std::vector<unsigned int>)
{
				   // what would that be in 1d? the
				   // identity matrix on the boundary
				   // dofs?
  Assert (false, ExcNotImplemented());
}

template <>
void MatrixCreator::create_boundary_mass_matrix (const Mapping<1,2>&,
						 const DoFHandler<1,2>&,
						 const Quadrature<0>&,
						 SparseMatrix<double>&,
						 const FunctionMap<2>::type&,
						 Vector<double>&,
						 std::vector<unsigned int>&,
						 const Function<2>* const,
						 std::vector<unsigned int>)
{
				   // what would that be in 1d? the
				   // identity matrix on the boundary
				   // dofs?
  Assert (false, ExcNotImplemented());
}


#else



template <int dim, int spacedim>
void
MatrixCreator::create_boundary_mass_matrix (const Mapping<dim, spacedim>  &mapping,
					    const DoFHandler<dim,spacedim> &dof,
					    const Quadrature<dim-1>  &q,
					    SparseMatrix<double>  &matrix,
					    const typename FunctionMap<spacedim>::type  &boundary_functions,
					    Vector<double>            &rhs_vector,
					    std::vector<unsigned int> &dof_to_boundary_mapping,
					    const Function<spacedim> * const coefficient,
					    std::vector<unsigned int> component_mapping)
{
  const FiniteElement<dim,spacedim> &fe = dof.get_fe();
  const unsigned int n_components  = fe.n_components();
  
  Assert (matrix.n() == dof.n_boundary_dofs(boundary_functions),
          ExcInternalError());
  Assert (matrix.n() == matrix.m(), ExcInternalError());
  Assert (matrix.n() == rhs_vector.size(), ExcInternalError());
  Assert (boundary_functions.size() != 0, ExcInternalError());
  Assert (dof_to_boundary_mapping.size() == dof.n_dofs(),
	  ExcInternalError());  
  Assert (coefficient ==0 ||
	  coefficient->n_components==1 ||
	  coefficient->n_components==n_components, ExcComponentMismatch());

  if (component_mapping.size() == 0)
    {
      AssertDimension (n_components, boundary_functions.begin()->second->n_components);
      for (unsigned int i=0;i<n_components;++i)
	component_mapping.push_back(i);
    }
  else
    AssertDimension (n_components, component_mapping.size());
  
  const unsigned int n_threads = multithread_info.n_default_threads;
  Threads::ThreadGroup<> threads;

				   // define starting and end point
				   // for each thread
  typedef typename DoFHandler<dim,spacedim>::active_cell_iterator active_cell_iterator;  
  std::vector<std::pair<active_cell_iterator,active_cell_iterator> > thread_ranges
    = Threads::split_range<active_cell_iterator> (dof.begin_active(),
						  dof.end(), n_threads);

				   // mutex to synchronise access to
				   // the matrix
  Threads::ThreadMutex mutex;

  typedef std_cxx1x::tuple<const Mapping<dim,spacedim>&,
    const DoFHandler<dim,spacedim>&,
    const Quadrature<dim-1>&> Commons;
  
				   // then assemble in parallel
  typedef void (*create_boundary_mass_matrix_1_t)
      (Commons,
       SparseMatrix<double>      &matrix,
       const typename FunctionMap<spacedim>::type &boundary_functions,
       Vector<double>            &rhs_vector,
       std::vector<unsigned int> &dof_to_boundary_mapping,
       const Function<spacedim> * const coefficient,
       const std::vector<unsigned int>& component_mapping,
       const IteratorRange<DoFHandler<dim,spacedim> >   range,
       Threads::ThreadMutex      &mutex);
  create_boundary_mass_matrix_1_t p = &MatrixCreator::template create_boundary_mass_matrix_1<dim,spacedim>;

//TODO: Use WorkStream here  
  for (unsigned int thread=0; thread<n_threads; ++thread)
    threads += Threads::new_thread (p,
				    Commons(mapping, dof, q), matrix,
				    boundary_functions, rhs_vector,
				    dof_to_boundary_mapping, coefficient,
				    component_mapping,
				    thread_ranges[thread], mutex);
  threads.join_all ();  
}


#if deal_II_dimension == 2

template <>
void
MatrixCreator::
create_boundary_mass_matrix_1<2,3> (std_cxx1x::tuple<const Mapping<2,3> &,
				    const DoFHandler<2,3> &,
			            const Quadrature<1> & > ,
			            SparseMatrix<double>  &,
			            const FunctionMap<3>::type &,
			            Vector<double> &,
			            std::vector<unsigned int> &,
			            const Function<3> * const ,
			            const std::vector<unsigned int> &,
				    const IteratorRange<DoFHandler<2,3> > ,
			            Threads::ThreadMutex &)
{
  Assert(false,ExcNotImplemented());
}

#endif


template <int dim, int spacedim>
void
MatrixCreator::
create_boundary_mass_matrix_1 (std_cxx1x::tuple<const Mapping<dim, spacedim> &,
			       const DoFHandler<dim,spacedim> &,
			       const Quadrature<dim-1> & >  commons,
			       SparseMatrix<double>      &matrix,
			       const typename FunctionMap<spacedim>::type &boundary_functions,
			       Vector<double>            &rhs_vector,
			       std::vector<unsigned int> &dof_to_boundary_mapping,
			       const Function<spacedim> * const coefficient,
			       const std::vector<unsigned int>& component_mapping,
			       const IteratorRange<DoFHandler<dim,spacedim> >   range,
			       Threads::ThreadMutex      &mutex)
{
				   // All assertions for this function
				   // are in the calling function
				   // before creating threads.
  const Mapping<dim,spacedim>& mapping = std_cxx1x::get<0>(commons);
  const DoFHandler<dim,spacedim>& dof = std_cxx1x::get<1>(commons);
  const Quadrature<dim-1>& q = std_cxx1x::get<2>(commons);
  
  const FiniteElement<dim,spacedim> &fe = dof.get_fe();
  const unsigned int n_components  = fe.n_components();
  const unsigned int n_function_components = boundary_functions.begin()->second->n_components;
  const bool         fe_is_system  = (n_components != 1);
  const bool         fe_is_primitive = fe.is_primitive();  
  
  const unsigned int dofs_per_cell = fe.dofs_per_cell,
		     dofs_per_face = fe.dofs_per_face;
  
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_vector(dofs_per_cell);


  UpdateFlags update_flags = UpdateFlags (update_values     |
					  update_JxW_values |
					  update_normal_vectors |
					  update_quadrature_points);
  FEFaceValues<dim,spacedim> fe_values (mapping, fe, q, update_flags);

				   // two variables for the coefficient,
				   // one for the two cases indicated in
				   // the name
  std::vector<double>          coefficient_values (fe_values.n_quadrature_points, 1.);
  std::vector<Vector<double> > coefficient_vector_values (fe_values.n_quadrature_points,
							  Vector<double>(n_components));
  const bool coefficient_is_vector = (coefficient != 0 && coefficient->n_components != 1)
				     ? true : false;
  
  std::vector<double>          rhs_values_scalar (fe_values.n_quadrature_points);
  std::vector<Vector<double> > rhs_values_system (fe_values.n_quadrature_points,
						  Vector<double>(n_function_components));

  std::vector<unsigned int> dofs (dofs_per_cell);
  std::vector<unsigned int> dofs_on_face_vector (dofs_per_face);
  
				   // for each dof on the cell, have a
				   // flag whether it is on the face
  std::vector<bool>         dof_is_on_face(dofs_per_cell);
  
  typename DoFHandler<dim,spacedim>::active_cell_iterator cell = range.first;
  for (; cell!=range.second; ++cell)
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
				       // check if this face is on that part of
				       // the boundary we are interested in
      if (boundary_functions.find(cell->face(face)->boundary_indicator()) !=
	  boundary_functions.end())
	{
	  cell_matrix = 0;
	  cell_vector = 0;
	  
	  fe_values.reinit (cell, face);
	  
	  if (fe_is_system)
					     // FE has several components
	    {
	      boundary_functions.find(cell->face(face)->boundary_indicator())
		->second->vector_value_list (fe_values.get_quadrature_points(),
					     rhs_values_system);

	      if (coefficient_is_vector)
						 // If coefficient is
						 // vector valued, fill
						 // all components
		coefficient->vector_value_list (fe_values.get_quadrature_points(),
						coefficient_vector_values);
	      else
		{
						   // If a scalar
						   // function is
						   // geiven, update
						   // the values, if
						   // not, use the
						   // default one set
						   // in the
						   // constructor above
		  if (coefficient != 0)
		    coefficient->value_list (fe_values.get_quadrature_points(),
					     coefficient_values);
						   // Copy scalar
						   // values into vector
		  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
		    coefficient_vector_values[point] = coefficient_values[point];
		}

					       // Special treatment
					       // for Hdiv and Hcurl
					       // elements, where only
					       // the normal or
					       // tangential component
					       // should be projected.
	      std::vector<std::vector<double> > normal_adjustment(fe_values.n_quadrature_points,
								  std::vector<double>(n_components, 1.));
	      
	      for (unsigned int comp = 0;comp<n_components;++comp)
		{
		  const FiniteElement<dim,spacedim>& base = fe.base_element(fe.component_to_base_index(comp).first);
		  const unsigned int bcomp = fe.component_to_base_index(comp).second;
		  
		  if (!base.conforms(FiniteElementData<dim>::H1) &&
		      base.conforms(FiniteElementData<dim>::Hdiv))
		    for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
		      normal_adjustment[point][comp] = fe_values.normal_vector(point)(bcomp)
						       * fe_values.normal_vector(point)(bcomp);
		}
	      
	      for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
		{
		  const double weight = fe_values.JxW(point);
		  for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i)
		    if (fe_is_primitive)
		      {
			for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
			  {
			    if (fe.system_to_component_index(j).first
				== fe.system_to_component_index(i).first)
			      {
				cell_matrix(i,j)
				  += weight
				  * fe_values.shape_value(j,point)
				  * fe_values.shape_value(i,point)
				  * coefficient_vector_values[point](fe.system_to_component_index(i).first);
			      }
			  }
			cell_vector(i) += fe_values.shape_value(i,point)
					  * rhs_values_system[point](component_mapping[fe.system_to_component_index(i).first])
					  * weight;
		      }
		    else
		      {
			for (unsigned int comp=0;comp<n_components;++comp)
			  {
			    for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
			      cell_matrix(i,j)
				+= fe_values.shape_value_component(j,point,comp)
				* fe_values.shape_value_component(i,point,comp)
				* normal_adjustment[point][comp]
				* weight * coefficient_vector_values[point](comp);
			    cell_vector(i) += fe_values.shape_value_component(i,point,comp) *
					      rhs_values_system[point](component_mapping[comp])
					      * normal_adjustment[point][comp]
					      * weight;
			  }
		      }
		}
	    }
	  else
					     // FE is a scalar one
	    {
	      boundary_functions.find(cell->face(face)->boundary_indicator())
		->second->value_list (fe_values.get_quadrature_points(), rhs_values_scalar);
	      
	      if (coefficient != 0)
		coefficient->value_list (fe_values.get_quadrature_points(),
					 coefficient_values);
	      for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
		{
		  const double weight = fe_values.JxW(point);
		  for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i) 
		    {
		      const double v = fe_values.shape_value(i,point);
		      for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
			{
			  const double u = fe_values.shape_value(j,point);
			  cell_matrix(i,j) += (u * v * weight * coefficient_values[point]);
			}
		      cell_vector(i) += v * rhs_values_scalar[point] *weight;
		    }
		}
	    }
					   // now transfer cell matrix and vector
					   // to the whole boundary matrix
					   //
					   // in the following: dof[i] holds the
					   // global index of the i-th degree of
					   // freedom on the present cell. If it
					   // is also a dof on the boundary, it
					   // must be a nonzero entry in the
					   // dof_to_boundary_mapping and then
					   // the boundary index of this dof is
					   // dof_to_boundary_mapping[dof[i]].
					   //
					   // if dof[i] is not on the boundary,
					   // it should be zero on the boundary
					   // therefore on all quadrature
					   // points and finally all of its
					   // entries in the cell matrix and
					   // vector should be zero. If not, we
					   // throw an error (note: because of
					   // the evaluation of the shape
					   // functions only up to machine
					   // precision, the term "must be zero"
					   // really should mean: "should be
					   // very small". since this is only an
					   // assertion and not part of the
					   // code, we may choose "very small"
					   // quite arbitrarily)
					   //
					   // the main problem here is that the
					   // matrix or vector entry should also
					   // be zero if the degree of freedom
					   // dof[i] is on the boundary, but not
					   // on the present face, i.e. on
					   // another face of the same cell also
					   // on the boundary. We can therefore
					   // not rely on the
					   // dof_to_boundary_mapping[dof[i]]
					   // being !=-1, we really have to
					   // determine whether dof[i] is a
					   // dof on the present face. We do so
					   // by getting the dofs on the
					   // face into @p{dofs_on_face_vector},
					   // a vector as always. Usually,
					   // searching in a vector is
					   // inefficient, so we copy the dofs
					   // into a set, which enables binary
					   // searches.
	  cell->get_dof_indices (dofs);
	  cell->face(face)->get_dof_indices (dofs_on_face_vector);
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    dof_is_on_face[i] = (std::find(dofs_on_face_vector.begin(),
					   dofs_on_face_vector.end(),
					   dofs[i])
				 !=
				 dofs_on_face_vector.end());
	  
                                           // lock the matrix
          Threads::ThreadMutex::ScopedLock lock (mutex);
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
	      if (dof_is_on_face[i] && dof_to_boundary_mapping[dofs[i]] != numbers::invalid_unsigned_int)
		{
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    if (dof_is_on_face[j] && dof_to_boundary_mapping[dofs[j]] != numbers::invalid_unsigned_int)
		      matrix.add(dof_to_boundary_mapping[dofs[i]],
				 dof_to_boundary_mapping[dofs[j]],
				 cell_matrix(i,j));
		  
		  rhs_vector(dof_to_boundary_mapping[dofs[i]]) += cell_vector(i);
		}
	    }
	}
}

#endif

template <int dim, int spacedim>
void MatrixCreator::create_boundary_mass_matrix (const DoFHandler<dim,spacedim>     &dof,
						 const Quadrature<dim-1>   &q,
						 SparseMatrix<double>      &matrix,
						 const typename FunctionMap<spacedim>::type &rhs,
						 Vector<double>            &rhs_vector,
						 std::vector<unsigned int> &dof_to_boundary_mapping,
						 const Function<spacedim> * const a,
						 std::vector<unsigned int> component_mapping)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_boundary_mass_matrix(StaticMappingQ1<dim,spacedim>::mapping, dof, q,
			      matrix,rhs, rhs_vector, dof_to_boundary_mapping, a, component_mapping);
}



#if deal_II_dimension == 1

template <int dim, int spacedim>
void
MatrixCreator::create_boundary_mass_matrix (const hp::MappingCollection<dim,spacedim>&,
					    const hp::DoFHandler<dim,spacedim>&,
					    const hp::QCollection<dim-1>&,
					    SparseMatrix<double>&,
					    const typename FunctionMap<spacedim>::type&,
					    Vector<double>&,
					    std::vector<unsigned int>&,
					    const Function<spacedim>* const,
					    std::vector<unsigned int>)
{
				   // what would that be in 1d? the
				   // identity matrix on the boundary
				   // dofs?
  Assert (false, ExcNotImplemented());
}


#else


template <int dim, int spacedim>
void
MatrixCreator::create_boundary_mass_matrix (const hp::MappingCollection<dim,spacedim>        &mapping,
					    const hp::DoFHandler<dim,spacedim>     &dof,
					    const hp::QCollection<dim-1>   &q,
					    SparseMatrix<double>      &matrix,
					    const typename FunctionMap<spacedim>::type         &boundary_functions,
					    Vector<double>            &rhs_vector,
					    std::vector<unsigned int> &dof_to_boundary_mapping,
					    const Function<spacedim> * const coefficient,
					    std::vector<unsigned int> component_mapping)
{
  const hp::FECollection<dim,spacedim> &fe_collection = dof.get_fe();
  const unsigned int n_components  = fe_collection.n_components();
  
  Assert (matrix.n() == dof.n_boundary_dofs(boundary_functions),
          ExcInternalError());
  Assert (matrix.n() == matrix.m(), ExcInternalError());
  Assert (matrix.n() == rhs_vector.size(), ExcInternalError());
  Assert (boundary_functions.size() != 0, ExcInternalError());
  Assert (dof_to_boundary_mapping.size() == dof.n_dofs(),
	  ExcInternalError());
  Assert (coefficient ==0 ||
	  coefficient->n_components==1 ||
	  coefficient->n_components==n_components, ExcComponentMismatch());

  if (component_mapping.size() == 0)
    {
      AssertDimension (n_components, boundary_functions.begin()->second->n_components);
      for (unsigned int i=0;i<n_components;++i)
	component_mapping.push_back(i);
    }
  else
    AssertDimension (n_components, component_mapping.size());

  const unsigned int n_threads = multithread_info.n_default_threads;
  Threads::ThreadGroup<> threads;

				   // define starting and end point
				   // for each thread
  typedef typename hp::DoFHandler<dim,spacedim>::active_cell_iterator active_cell_iterator;  
  std::vector<std::pair<active_cell_iterator,active_cell_iterator> > thread_ranges
    = Threads::split_range<active_cell_iterator> (dof.begin_active(),
						  dof.end(), n_threads);

  typedef std_cxx1x::tuple<const hp::MappingCollection<dim,spacedim>&,
    const hp::DoFHandler<dim,spacedim>&,
    const hp::QCollection<dim-1>&> Commons;
  
				   // mutex to synchronise access to
				   // the matrix
  Threads::ThreadMutex mutex;
  
				   // then assemble in parallel
  typedef void (*create_boundary_mass_matrix_1_t)
      (Commons,
       SparseMatrix<double>      &matrix,
       const typename FunctionMap<spacedim>::type &boundary_functions,
       Vector<double>            &rhs_vector,
       std::vector<unsigned int> &dof_to_boundary_mapping,
       const Function<spacedim> * const coefficient,
       const std::vector<unsigned int>& component_mapping,
       const IteratorRange<hp::DoFHandler<dim,spacedim> >   range,
       Threads::ThreadMutex      &mutex);
  create_boundary_mass_matrix_1_t p = &MatrixCreator::template create_boundary_mass_matrix_1<dim,spacedim>;

//TODO: Use WorkStream here  
  for (unsigned int thread=0; thread<n_threads; ++thread)
    threads += Threads::new_thread (p,
				    Commons(mapping, dof, q), matrix,
				    boundary_functions, rhs_vector,
				    dof_to_boundary_mapping, coefficient,
				    component_mapping,
				    thread_ranges[thread], mutex);
  threads.join_all ();  
}



template <int dim, int spacedim>
void
MatrixCreator::
create_boundary_mass_matrix_1 (std_cxx1x::tuple<const hp::MappingCollection<dim,spacedim> &,
			       const hp::DoFHandler<dim,spacedim> &,
			       const hp::QCollection<dim-1> &> commons,
			       SparseMatrix<double>      &matrix,
			       const typename FunctionMap<spacedim>::type &boundary_functions,
			       Vector<double>            &rhs_vector,
			       std::vector<unsigned int> &dof_to_boundary_mapping,
			       const Function<spacedim> * const coefficient,
			       const std::vector<unsigned int>& component_mapping,
			       const IteratorRange<hp::DoFHandler<dim,spacedim> >   range,
			       Threads::ThreadMutex      &mutex)
{
  const hp::MappingCollection<dim,spacedim>& mapping = std_cxx1x::get<0>(commons);
  const hp::DoFHandler<dim,spacedim>& dof = std_cxx1x::get<1>(commons);
  const hp::QCollection<dim-1>& q = std_cxx1x::get<2>(commons);
  const hp::FECollection<dim,spacedim> &fe_collection = dof.get_fe();
  const unsigned int n_components  = fe_collection.n_components();
  const unsigned int n_function_components = boundary_functions.begin()->second->n_components;
  const bool         fe_is_system  = (n_components != 1);
#ifdef DEBUG
  if (true)
    {
      unsigned int max_element = 0;
      for (std::vector<unsigned int>::const_iterator i=dof_to_boundary_mapping.begin();
	   i!=dof_to_boundary_mapping.end(); ++i)
	if ((*i != hp::DoFHandler<dim,spacedim>::invalid_dof_index) &&
	    (*i > max_element))
	  max_element = *i;
      Assert (max_element  == matrix.n()-1, ExcInternalError());
    };
#endif
  
  const unsigned int max_dofs_per_cell = fe_collection.max_dofs_per_cell(),
		     max_dofs_per_face = fe_collection.max_dofs_per_face();
  
  FullMatrix<double> cell_matrix(max_dofs_per_cell, max_dofs_per_cell);
  Vector<double>     cell_vector(max_dofs_per_cell);


  UpdateFlags update_flags = UpdateFlags (update_values     |
					  update_JxW_values |
					  update_quadrature_points);
  hp::FEFaceValues<dim> x_fe_values (mapping, fe_collection, q, update_flags);

				   // two variables for the coefficient,
				   // one for the two cases indicated in
				   // the name
  std::vector<double>          coefficient_values;
  std::vector<Vector<double> > coefficient_vector_values;

  std::vector<double>          rhs_values_scalar;
  std::vector<Vector<double> > rhs_values_system;

  std::vector<unsigned int> dofs (max_dofs_per_cell);
  std::vector<unsigned int> dofs_on_face_vector (max_dofs_per_face);

				   // for each dof on the cell, have a
				   // flag whether it is on the face
  std::vector<bool>         dof_is_on_face(max_dofs_per_cell);


  typename hp::DoFHandler<dim,spacedim>::active_cell_iterator cell = range.first;
  for (; cell!=range.second; ++cell)
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
				       // check if this face is on that part of
				       // the boundary we are interested in
      if (boundary_functions.find(cell->face(face)->boundary_indicator()) !=
	  boundary_functions.end())
	{
	  x_fe_values.reinit (cell, face);

	  const FEFaceValues<dim,spacedim> &fe_values = x_fe_values.get_present_fe_values ();

	  const FiniteElement<dim,spacedim> &fe = cell->get_fe();
	  const unsigned int dofs_per_cell = fe.dofs_per_cell;
	  const unsigned int dofs_per_face = fe.dofs_per_face;
	  
	  cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
	  cell_vector.reinit (dofs_per_cell);
	  cell_matrix = 0;
	  cell_vector = 0;	  
	  
	  if (fe_is_system)
					     // FE has several components
	    {
	      rhs_values_system.resize (fe_values.n_quadrature_points,
					Vector<double>(n_function_components));
	      boundary_functions.find(cell->face(face)->boundary_indicator())
		->second->vector_value_list (fe_values.get_quadrature_points(),
					     rhs_values_system);
	      
	      if (coefficient != 0)
		{
		  if (coefficient->n_components==1)
		    {
		      coefficient_values.resize (fe_values.n_quadrature_points);
		      coefficient->value_list (fe_values.get_quadrature_points(),
					       coefficient_values);
		      for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
			{
			  const double weight = fe_values.JxW(point);
			  for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i) 
			    {
			      const double v = fe_values.shape_value(i,point);
			      for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
				if (fe.system_to_component_index(i).first ==
				    fe.system_to_component_index(j).first)
				  {
				    const double u = fe_values.shape_value(j,point);
				    cell_matrix(i,j)
				      += (u * v * weight * coefficient_values[point]);
				  }
			      
			      cell_vector(i) += v *
						rhs_values_system[point](
						  component_mapping[fe.system_to_component_index(i).first]) * weight;
			    }
			}
		    }
		  else
		    {
		      coefficient_vector_values.resize (fe_values.n_quadrature_points,
							Vector<double>(n_components));
		      coefficient->vector_value_list (fe_values.get_quadrature_points(),
						      coefficient_vector_values);
		      for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
			{
			  const double weight = fe_values.JxW(point);
			  for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i) 
			    {
			      const double v = fe_values.shape_value(i,point);
			      const unsigned int component_i=
				fe.system_to_component_index(i).first;
			      for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
				if (fe.system_to_component_index(j).first ==
				    component_i)
				  {
				    const double u = fe_values.shape_value(j,point);
				    cell_matrix(i,j) +=
				      (u * v * weight * coefficient_vector_values[point](component_i));
				  }
			      cell_vector(i) += v * rhs_values_system[point](component_mapping[component_i]) * weight;
			    }
			}
		    }
		}
	      else	//      if (coefficient == 0)
		for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
		  {
		    const double weight = fe_values.JxW(point);
		    for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i) 
		      {
			const double v = fe_values.shape_value(i,point);
			for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
			  if (fe.system_to_component_index(i).first ==
			      fe.system_to_component_index(j).first)
			    {
			      const double u = fe_values.shape_value(j,point);
			      cell_matrix(i,j) += (u * v * weight);
			    }
			cell_vector(i) += v *
					  rhs_values_system[point](
					    fe.system_to_component_index(i).first) *
					  weight;
		      }
		  }
	    }
	  else
					     // FE is a scalar one
	    {
	      rhs_values_scalar.resize (fe_values.n_quadrature_points);
	      boundary_functions.find(cell->face(face)->boundary_indicator())
		->second->value_list (fe_values.get_quadrature_points(), rhs_values_scalar);
	      
	      if (coefficient != 0)
		{
		  coefficient_values.resize (fe_values.n_quadrature_points);
		  coefficient->value_list (fe_values.get_quadrature_points(),
					   coefficient_values);
		  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
		    {
		      const double weight = fe_values.JxW(point);
		      for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i) 
			{
			  const double v = fe_values.shape_value(i,point);
			  for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
			    {
			      const double u = fe_values.shape_value(j,point);
			      cell_matrix(i,j) += (u * v * weight * coefficient_values[point]);
			    }
			  cell_vector(i) += v * rhs_values_scalar[point] *weight;
			}
		    }
		}
	      else
		for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
		  {
		    const double weight = fe_values.JxW(point);
		    for (unsigned int i=0; i<fe_values.dofs_per_cell; ++i) 
		      {
			const double v = fe_values.shape_value(i,point);
			for (unsigned int j=0; j<fe_values.dofs_per_cell; ++j)
			  {
			    const double u = fe_values.shape_value(j,point);
			    cell_matrix(i,j) += (u * v * weight);
			  }
			cell_vector(i) += v * rhs_values_scalar[point] * weight;
		      }
		  }
	    }

					   // now transfer cell matrix and vector
					   // to the whole boundary matrix
					   //
					   // in the following: dof[i] holds the
					   // global index of the i-th degree of
					   // freedom on the present cell. If it
					   // is also a dof on the boundary, it
					   // must be a nonzero entry in the
					   // dof_to_boundary_mapping and then
					   // the boundary index of this dof is
					   // dof_to_boundary_mapping[dof[i]].
					   //
					   // if dof[i] is not on the boundary,
					   // it should be zero on the boundary
					   // therefore on all quadrature
					   // points and finally all of its
					   // entries in the cell matrix and
					   // vector should be zero. If not, we
					   // throw an error (note: because of
					   // the evaluation of the shape
					   // functions only up to machine
					   // precision, the term "must be zero"
					   // really should mean: "should be
					   // very small". since this is only an
					   // assertion and not part of the
					   // code, we may choose "very small"
					   // quite arbitrarily)
					   //
					   // the main problem here is that the
					   // matrix or vector entry should also
					   // be zero if the degree of freedom
					   // dof[i] is on the boundary, but not
					   // on the present face, i.e. on
					   // another face of the same cell also
					   // on the boundary. We can therefore
					   // not rely on the
					   // dof_to_boundary_mapping[dof[i]]
					   // being !=-1, we really have to
					   // determine whether dof[i] is a
					   // dof on the present face. We do so
					   // by getting the dofs on the
					   // face into @p{dofs_on_face_vector},
					   // a vector as always. Usually,
					   // searching in a vector is
					   // inefficient, so we copy the dofs
					   // into a set, which enables binary
					   // searches.
	  dofs.resize (dofs_per_cell);
	  dofs_on_face_vector.resize (dofs_per_face);
	  dof_is_on_face.resize (dofs_per_cell);

	  cell->get_dof_indices (dofs);
	  cell->face(face)->get_dof_indices (dofs_on_face_vector,
					     cell->active_fe_index());

					   // check for each of the
					   // dofs on this cell
					   // whether it is on the
					   // face
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    dof_is_on_face[i] = (std::find(dofs_on_face_vector.begin(),
					   dofs_on_face_vector.end(),
					   dofs[i])
				 !=
				 dofs_on_face_vector.end());
	  
					   // in debug mode: compute an element
					   // in the matrix which is
					   // guaranteed to belong to a boundary
					   // dof. We do this to check that the
					   // entries in the cell matrix are
					   // guaranteed to be zero if the
					   // respective dof is not on the
					   // boundary. Since because of
					   // round-off, the actual
					   // value of the matrix entry may be
					   // only close to zero, we assert that
					   // it is small relative to an element
					   // which is guaranteed to be nonzero.
					   // (absolute smallness does not
					   // suffice since the size of the
					   // domain scales in here)
					   //
					   // for this purpose we seek the
					   // diagonal of the matrix, where there
					   // must be an element belonging to
					   // the boundary. we take the maximum
					   // diagonal entry.
#ifdef DEBUG
	  double max_diag_entry = 0;
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    if (std::fabs(cell_matrix(i,i)) > max_diag_entry)
	      max_diag_entry = std::fabs(cell_matrix(i,i));
#endif  

                                           // lock the matrix
          Threads::ThreadMutex::ScopedLock lock (mutex);
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      if (dof_is_on_face[i] && dof_is_on_face[j])
		matrix.add(dof_to_boundary_mapping[dofs[i]],
			   dof_to_boundary_mapping[dofs[j]],
			   cell_matrix(i,j));
	      else
		{
						   // assume that all
						   // shape functions
						   // that are nonzero
						   // on the boundary
						   // are also listed
						   // in the
						   // @p{dof_to_boundary}
						   // mapping. if that
						   // is not the case,
						   // then the
						   // boundary mass
						   // matrix does not
						   // make that much
						   // sense anyway, as
						   // it only contains
						   // entries for
						   // parts of the
						   // functions living
						   // on the boundary
						   //
						   // these, we may
						   // compare here for
						   // relative
						   // smallness of all
						   // entries in the
						   // local matrix
						   // which are not
						   // taken over to
						   // the global one
		  Assert (std::fabs(cell_matrix(i,j)) <= 1e-10 * max_diag_entry,
			  ExcInternalError ());
		};
	  
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    if (dof_is_on_face[j])
	      rhs_vector(dof_to_boundary_mapping[dofs[j]]) += cell_vector(j);
	    else
	      {
						 // compare here for relative
						 // smallness
		Assert (std::fabs(cell_vector(j)) <= 1e-10 * max_diag_entry,
			ExcInternalError());
	      }
	}
}

#endif

template <int dim, int spacedim>
void MatrixCreator::create_boundary_mass_matrix (const hp::DoFHandler<dim,spacedim>     &dof,
						 const hp::QCollection<dim-1>   &q,
						 SparseMatrix<double>      &matrix,
						 const typename FunctionMap<spacedim>::type &rhs,
						 Vector<double>            &rhs_vector,
						 std::vector<unsigned int> &dof_to_boundary_mapping,
						 const Function<spacedim> * const a,
						 std::vector<unsigned int> component_mapping)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_boundary_mass_matrix(hp::StaticMappingQ1<dim,spacedim>::mapping_collection, dof, q,
			      matrix,rhs, rhs_vector, dof_to_boundary_mapping, a, component_mapping);
}



template <int dim, int spacedim>
void MatrixCreator::create_laplace_matrix (const Mapping<dim, spacedim>       &mapping,
					   const DoFHandler<dim,spacedim>    &dof,
					   const Quadrature<dim>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<spacedim> * const coefficient)
{
  Assert (matrix.m() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
  Assert (matrix.n() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

  internal::MatrixCreator::AssemblerData::Scratch<dim, spacedim>
    assembler_data (dof.get_fe(),
		    update_gradients  | update_JxW_values |
		    (coefficient != 0 ? update_quadrature_points : UpdateFlags(0)),
		    coefficient, /*rhs_function=*/0,
		    q, mapping);
  internal::MatrixCreator::AssemblerData::CopyData copy_data;
  copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
			 assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());  

  WorkStream::run (dof.begin_active(),
		   static_cast<typename DoFHandler<dim,spacedim>::active_cell_iterator>(dof.end()),
		   &internal::MatrixCreator::laplace_assembler<dim, spacedim, typename DoFHandler<dim,spacedim>::active_cell_iterator>,
		   std_cxx1x::bind (&internal::MatrixCreator::
				    copy_local_to_global<SparseMatrix<double>, Vector<double> >,
				    _1, &matrix, (Vector<double>*)0),
		   assembler_data,
		   copy_data);
}



template <int dim, int spacedim>
void MatrixCreator::create_laplace_matrix (const DoFHandler<dim,spacedim>    &dof,
					   const Quadrature<dim>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<spacedim> * const coefficient)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_laplace_matrix(StaticMappingQ1<dim,spacedim>::mapping, dof, q, matrix, coefficient);
}



template <int dim, int spacedim>
void MatrixCreator::create_laplace_matrix (const Mapping<dim, spacedim>       &mapping,
					   const DoFHandler<dim,spacedim>    &dof,
					   const Quadrature<dim>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<spacedim>      &rhs,
					   Vector<double>           &rhs_vector,
					   const Function<spacedim> * const coefficient)
{
  Assert (matrix.m() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
  Assert (matrix.n() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

  internal::MatrixCreator::AssemblerData::Scratch<dim, spacedim>
    assembler_data (dof.get_fe(),
		    update_gradients  | update_values |
		    update_JxW_values | update_quadrature_points,
		    coefficient, &rhs,
		    q, mapping);
  internal::MatrixCreator::AssemblerData::CopyData copy_data;
  copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
			 assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());
  WorkStream::run (dof.begin_active(),
		   static_cast<typename DoFHandler<dim,spacedim>::active_cell_iterator>(dof.end()),
		   &internal::MatrixCreator::laplace_assembler<dim, spacedim, typename DoFHandler<dim,spacedim>::active_cell_iterator>,
		   std_cxx1x::bind (&internal::MatrixCreator::
				    copy_local_to_global<SparseMatrix<double>, Vector<double> >,
				    _1, &matrix, &rhs_vector),
		   assembler_data,
		   copy_data);
}



template <int dim, int spacedim>
void MatrixCreator::create_laplace_matrix (const DoFHandler<dim,spacedim>    &dof,
					   const Quadrature<dim>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<spacedim>      &rhs,
					   Vector<double>           &rhs_vector,
					   const Function<spacedim> * const coefficient)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_laplace_matrix(StaticMappingQ1<dim,spacedim>::mapping, dof, q,
			matrix, rhs, rhs_vector, coefficient);
}



template <int dim, int spacedim>
void MatrixCreator::create_laplace_matrix (const hp::MappingCollection<dim,spacedim>       &mapping,
					   const hp::DoFHandler<dim,spacedim>    &dof,
					   const hp::QCollection<dim>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<spacedim> * const coefficient)
{
  Assert (matrix.m() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
  Assert (matrix.n() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

  internal::MatrixCreator::AssemblerData::Scratch<dim, spacedim>
    assembler_data (dof.get_fe(),
		    update_gradients  | update_JxW_values |
		    (coefficient != 0 ? update_quadrature_points : UpdateFlags(0)),
		    coefficient, /*rhs_function=*/0,
		    q, mapping);
  internal::MatrixCreator::AssemblerData::CopyData copy_data;
  copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
			 assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());  

  WorkStream::run (dof.begin_active(),
		   static_cast<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator>(dof.end()),
		   &internal::MatrixCreator::laplace_assembler<dim, spacedim, typename hp::DoFHandler<dim,spacedim>::active_cell_iterator>,
		   std_cxx1x::bind (&internal::MatrixCreator::
				    copy_local_to_global<SparseMatrix<double>, Vector<double> >,
				    _1, &matrix, (Vector<double>*)0),
		   assembler_data,
		   copy_data);
}



template <int dim, int spacedim>
void MatrixCreator::create_laplace_matrix (const hp::DoFHandler<dim,spacedim>    &dof,
					   const hp::QCollection<dim>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<spacedim> * const coefficient)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_laplace_matrix(hp::StaticMappingQ1<dim,spacedim>::mapping_collection, dof, q, matrix, coefficient);
}



template <int dim, int spacedim>
void MatrixCreator::create_laplace_matrix (const hp::MappingCollection<dim,spacedim>       &mapping,
					   const hp::DoFHandler<dim,spacedim>    &dof,
					   const hp::QCollection<dim>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<spacedim>      &rhs,
					   Vector<double>           &rhs_vector,
					   const Function<spacedim> * const coefficient)
{
  Assert (matrix.m() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
  Assert (matrix.n() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

  internal::MatrixCreator::AssemblerData::Scratch<dim, spacedim>
    assembler_data (dof.get_fe(),
		    update_gradients  | update_values |
		    update_JxW_values | update_quadrature_points,
		    coefficient, &rhs,
		    q, mapping);
  internal::MatrixCreator::AssemblerData::CopyData copy_data;
  copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
			 assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.cell_rhs.reinit (assembler_data.fe_collection.max_dofs_per_cell());
  copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());
  
  WorkStream::run (dof.begin_active(),
		   static_cast<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator>(dof.end()),
		   &internal::MatrixCreator::laplace_assembler<dim, spacedim, typename hp::DoFHandler<dim,spacedim>::active_cell_iterator>,
		   std_cxx1x::bind (&internal::MatrixCreator::
				    copy_local_to_global<SparseMatrix<double>, Vector<double> >,
				    _1, &matrix, &rhs_vector),
		   assembler_data,
		   copy_data);
}



template <int dim, int spacedim>
void MatrixCreator::create_laplace_matrix (const hp::DoFHandler<dim,spacedim>    &dof,
					   const hp::QCollection<dim>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<spacedim>      &rhs,
					   Vector<double>           &rhs_vector,
					   const Function<spacedim> * const coefficient)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_laplace_matrix(hp::StaticMappingQ1<dim,spacedim>::mapping_collection, dof, q,
			matrix, rhs, rhs_vector, coefficient);
}






// explicit instantiations

// non-hp version of create_mass_matrix
template
void MatrixCreator::create_mass_matrix<deal_II_dimension>
(const Mapping<deal_II_dimension>       &mapping,
 const DoFHandler<deal_II_dimension>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension>
(const DoFHandler<deal_II_dimension>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension>
(const Mapping<deal_II_dimension>       &mapping,
 const DoFHandler<deal_II_dimension>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension>      &rhs,
 Vector<double>           &rhs_vector,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension>
(const DoFHandler<deal_II_dimension>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension>      &rhs,
 Vector<double>           &rhs_vector,
 const Function<deal_II_dimension> * const coefficient);


template
void MatrixCreator::create_mass_matrix<deal_II_dimension>
(const Mapping<deal_II_dimension>       &mapping,
 const DoFHandler<deal_II_dimension>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<float>     &matrix,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension>
(const DoFHandler<deal_II_dimension>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<float>     &matrix,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension>
(const Mapping<deal_II_dimension>       &mapping,
 const DoFHandler<deal_II_dimension>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<float>     &matrix,
 const Function<deal_II_dimension>      &rhs,
 Vector<double>           &rhs_vector,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension>
(const DoFHandler<deal_II_dimension>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<float>     &matrix,
 const Function<deal_II_dimension>      &rhs,
 Vector<double>           &rhs_vector,
 const Function<deal_II_dimension> * const coefficient);

// hp versions of functions
template
void MatrixCreator::create_mass_matrix<deal_II_dimension>
(const hp::MappingCollection<deal_II_dimension>       &mapping,
 const hp::DoFHandler<deal_II_dimension>    &dof,
 const hp::QCollection<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension>
(const hp::DoFHandler<deal_II_dimension>    &dof,
 const hp::QCollection<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension>
(const hp::MappingCollection<deal_II_dimension>       &mapping,
 const hp::DoFHandler<deal_II_dimension>    &dof,
 const hp::QCollection<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension>      &rhs,
 Vector<double>           &rhs_vector,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension>
(const hp::DoFHandler<deal_II_dimension>    &dof,
 const hp::QCollection<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension>      &rhs,
 Vector<double>           &rhs_vector,
 const Function<deal_II_dimension> * const coefficient);


template
void MatrixCreator::create_mass_matrix<deal_II_dimension>
(const hp::MappingCollection<deal_II_dimension>       &mapping,
 const hp::DoFHandler<deal_II_dimension>    &dof,
 const hp::QCollection<deal_II_dimension>    &q,
 SparseMatrix<float>     &matrix,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension>
(const hp::DoFHandler<deal_II_dimension>    &dof,
 const hp::QCollection<deal_II_dimension>    &q,
 SparseMatrix<float>     &matrix,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension>
(const hp::MappingCollection<deal_II_dimension>       &mapping,
 const hp::DoFHandler<deal_II_dimension>    &dof,
 const hp::QCollection<deal_II_dimension>    &q,
 SparseMatrix<float>     &matrix,
 const Function<deal_II_dimension>      &rhs,
 Vector<double>           &rhs_vector,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension>
(const hp::DoFHandler<deal_II_dimension>    &dof,
 const hp::QCollection<deal_II_dimension>    &q,
 SparseMatrix<float>     &matrix,
 const Function<deal_II_dimension>      &rhs,
 Vector<double>           &rhs_vector,
 const Function<deal_II_dimension> * const coefficient);



template
void MatrixCreator::create_boundary_mass_matrix<deal_II_dimension>
(const DoFHandler<deal_II_dimension>     &dof,
 const Quadrature<deal_II_dimension-1>   &q,
 SparseMatrix<double>      &matrix,
 const FunctionMap<deal_II_dimension>::type &rhs,
 Vector<double>            &rhs_vector,
 std::vector<unsigned int> &dof_to_boundary_mapping,
 const Function<deal_II_dimension> * const a,
 std::vector<unsigned int>);

template
void MatrixCreator::create_boundary_mass_matrix<deal_II_dimension>
(const Mapping<deal_II_dimension> &,
 const DoFHandler<deal_II_dimension>     &dof,
 const Quadrature<deal_II_dimension-1>   &q,
 SparseMatrix<double>      &matrix,
 const FunctionMap<deal_II_dimension>::type &rhs,
 Vector<double>            &rhs_vector,
 std::vector<unsigned int> &dof_to_boundary_mapping,
 const Function<deal_II_dimension> * const a,
 std::vector<unsigned int>);

template
void
MatrixCreator::create_boundary_mass_matrix<deal_II_dimension>
(const hp::MappingCollection<deal_II_dimension>&,
 const hp::DoFHandler<deal_II_dimension>&,
 const hp::QCollection<deal_II_dimension-1>&,
 SparseMatrix<double>&,
 const FunctionMap<deal_II_dimension>::type&,
 Vector<double>&,
 std::vector<unsigned int>&,
 const Function<deal_II_dimension> * const,
 std::vector<unsigned int>);


template
void MatrixCreator::create_boundary_mass_matrix<deal_II_dimension>
(const hp::DoFHandler<deal_II_dimension>&,
 const hp::QCollection<deal_II_dimension-1>&,
 SparseMatrix<double>&,
 const FunctionMap<deal_II_dimension>::type&,
 Vector<double>&,
 std::vector<unsigned int>&,
 const Function<deal_II_dimension> * const,
 std::vector<unsigned int>);



// non-hp versions of create_laplace_matrix
template
void MatrixCreator::create_laplace_matrix<deal_II_dimension>
(const DoFHandler<deal_II_dimension>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_laplace_matrix<deal_II_dimension>
(const Mapping<deal_II_dimension>       &mapping,
 const DoFHandler<deal_II_dimension>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_laplace_matrix<deal_II_dimension>
(const Mapping<deal_II_dimension>       &mapping,
 const DoFHandler<deal_II_dimension>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension>      &rhs,
 Vector<double>           &rhs_vector,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_laplace_matrix<deal_II_dimension>
(const DoFHandler<deal_II_dimension>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension>      &rhs,
 Vector<double>           &rhs_vector,
 const Function<deal_II_dimension> * const coefficient);

// hp versions of create_laplace_matrix
template
void MatrixCreator::create_laplace_matrix<deal_II_dimension>
(const hp::DoFHandler<deal_II_dimension>    &dof,
 const hp::QCollection<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_laplace_matrix<deal_II_dimension>
(const hp::MappingCollection<deal_II_dimension>       &mapping,
 const hp::DoFHandler<deal_II_dimension>    &dof,
 const hp::QCollection<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_laplace_matrix<deal_II_dimension>
(const hp::MappingCollection<deal_II_dimension>       &mapping,
 const hp::DoFHandler<deal_II_dimension>    &dof,
 const hp::QCollection<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension>      &rhs,
 Vector<double>           &rhs_vector,
 const Function<deal_II_dimension> * const coefficient);
template
void MatrixCreator::create_laplace_matrix<deal_II_dimension>
(const hp::DoFHandler<deal_II_dimension>    &dof,
 const hp::QCollection<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension>      &rhs,
 Vector<double>           &rhs_vector,
 const Function<deal_II_dimension> * const coefficient);



#if deal_II_dimension != 3

// non-hp version of create_mass_matrix
template
void MatrixCreator::create_mass_matrix<deal_II_dimension,double,deal_II_dimension+1>
(const Mapping<deal_II_dimension,deal_II_dimension+1>       &mapping,
 const DoFHandler<deal_II_dimension,deal_II_dimension+1>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension+1> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension,double,deal_II_dimension+1>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension+1> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension,double,deal_II_dimension+1>
(const Mapping<deal_II_dimension,deal_II_dimension+1>       &mapping,
 const DoFHandler<deal_II_dimension,deal_II_dimension+1>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension+1>      &rhs,
 Vector<double>           &rhs_vector,
 const Function<deal_II_dimension+1> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension,double,deal_II_dimension+1>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<double>     &matrix,
 const Function<deal_II_dimension+1>      &rhs,
 Vector<double>           &rhs_vector,
 const Function<deal_II_dimension+1> * const coefficient);


template
void MatrixCreator::create_mass_matrix<deal_II_dimension,float,deal_II_dimension+1>
(const Mapping<deal_II_dimension,deal_II_dimension+1>       &mapping,
 const DoFHandler<deal_II_dimension,deal_II_dimension+1>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<float>     &matrix,
 const Function<deal_II_dimension+1> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension,float,deal_II_dimension+1>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<float>     &matrix,
 const Function<deal_II_dimension+1> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension,float,deal_II_dimension+1>
(const Mapping<deal_II_dimension,deal_II_dimension+1>       &mapping,
 const DoFHandler<deal_II_dimension,deal_II_dimension+1>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<float>     &matrix,
 const Function<deal_II_dimension+1>      &rhs,
 Vector<double>           &rhs_vector,
 const Function<deal_II_dimension+1> * const coefficient);
template
void MatrixCreator::create_mass_matrix<deal_II_dimension,float,deal_II_dimension+1>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>    &dof,
 const Quadrature<deal_II_dimension>    &q,
 SparseMatrix<float>     &matrix,
 const Function<deal_II_dimension+1>      &rhs,
 Vector<double>           &rhs_vector,
 const Function<deal_II_dimension+1> * const coefficient);

template
void MatrixCreator::create_boundary_mass_matrix<deal_II_dimension,deal_II_dimension+1>
(const DoFHandler<deal_II_dimension,deal_II_dimension+1>     &dof,
 const Quadrature<deal_II_dimension-1>   &q,
 SparseMatrix<double>      &matrix,
 const FunctionMap<deal_II_dimension+1>::type &rhs,
 Vector<double>            &rhs_vector,
 std::vector<unsigned int> &dof_to_boundary_mapping,
 const Function<deal_II_dimension+1> * const a,
 std::vector<unsigned int>);

template
void MatrixCreator::create_boundary_mass_matrix<deal_II_dimension,deal_II_dimension+1>
(const Mapping<deal_II_dimension,deal_II_dimension+1> &,
 const DoFHandler<deal_II_dimension,deal_II_dimension+1>     &dof,
 const Quadrature<deal_II_dimension-1>   &q,
 SparseMatrix<double>      &matrix,
 const FunctionMap<deal_II_dimension+1>::type &rhs,
 Vector<double>            &rhs_vector,
 std::vector<unsigned int> &dof_to_boundary_mapping,
 const Function<deal_II_dimension+1> * const a,
 std::vector<unsigned int>);

 

// #if deal_II_dimension != 1
// template
// void
// MatrixCreator::create_boundary_mass_matrix<deal_II_dimension,deal_II_dimension+1>
// (const Mapping<deal_II_dimension,deal_II_dimension+1>        &mapping,
//  const DoFHandler<deal_II_dimension,deal_II_dimension+1>     &dof,
//  const Quadrature<deal_II_dimension-1>   &q,
//  SparseMatrix<double>      &matrix,
//  const FunctionMap<deal_II_dimension+1>::type         &boundary_functions,
//  Vector<double>            &rhs_vector,
//  std::vector<unsigned int> &dof_to_boundary_mapping,
//  const Function<deal_II_dimension+1> * const a,
//  std::vector<unsigned int> &component_mapping);
// #endif

// template
// void MatrixCreator::create_boundary_mass_matrix<deal_II_dimension,deal_II_dimension+1>
// (const DoFHandler<deal_II_dimension,deal_II_dimension+1>     &dof,
//  const Quadrature<deal_II_dimension-1>   &q,
//  SparseMatrix<double>      &matrix,
//  const FunctionMap<deal_II_dimension+1>::type &rhs,
//  Vector<double>            &rhs_vector,
//  std::vector<unsigned int> &dof_to_boundary_mapping,
//  const Function<deal_II_dimension+1> * const a,
//  std::vector<unsigned int> &component_mapping);


// // non-hp version of create_mass_matrix
// template
// void MatrixCreator::create_mass_matrix
// (const Mapping<deal_II_dimension,deal_II_dimension+1>       &mapping,
//  const DoFHandler<deal_II_dimension,deal_II_dimension+1>    &dof,
//  const Quadrature<deal_II_dimension>    &q,
//  SparseMatrix<double>     &matrix,
//  const Function<deal_II_dimension+1> * const coefficient);
// template
// void MatrixCreator::create_mass_matrix
// (const DoFHandler<deal_II_dimension,deal_II_dimension+1>    &dof,
//  const Quadrature<deal_II_dimension>    &q,
//  SparseMatrix<double>     &matrix,
//  const Function<deal_II_dimension+1> * const coefficient);
// template
// void MatrixCreator::create_mass_matrix
// (const Mapping<deal_II_dimension,deal_II_dimension+1>       &mapping,
//  const DoFHandler<deal_II_dimension,deal_II_dimension+1>    &dof,
//  const Quadrature<deal_II_dimension>    &q,
//  SparseMatrix<double>     &matrix,
//  const Function<deal_II_dimension+1>      &rhs,
//  Vector<double>           &rhs_vector,
//  const Function<deal_II_dimension+1> * const coefficient);
// template
// void MatrixCreator::create_mass_matrix
// (const DoFHandler<deal_II_dimension,deal_II_dimension+1>    &dof,
//  const Quadrature<deal_II_dimension>    &q,
//  SparseMatrix<double>     &matrix,
//  const Function<deal_II_dimension+1>      &rhs,
//  Vector<double>           &rhs_vector,
//  const Function<deal_II_dimension+1> * const coefficient);


// template
// void MatrixCreator::create_mass_matrix
// (const Mapping<deal_II_dimension,deal_II_dimension+1>       &mapping,
//  const DoFHandler<deal_II_dimension,deal_II_dimension+1>    &dof,
//  const Quadrature<deal_II_dimension>    &q,
//  SparseMatrix<float>     &matrix,
//  const Function<deal_II_dimension+1> * const coefficient);
// template
// void MatrixCreator::create_mass_matrix
// (const DoFHandler<deal_II_dimension,deal_II_dimension+1>    &dof,
//  const Quadrature<deal_II_dimension>    &q,
//  SparseMatrix<float>     &matrix,
//  const Function<deal_II_dimension+1> * const coefficient);
// template
// void MatrixCreator::create_mass_matrix
// (const Mapping<deal_II_dimension,deal_II_dimension+1>       &mapping,
//  const DoFHandler<deal_II_dimension,deal_II_dimension+1>    &dof,
//  const Quadrature<deal_II_dimension>    &q,
//  SparseMatrix<float>     &matrix,
//  const Function<deal_II_dimension+1>      &rhs,
//  Vector<double>           &rhs_vector,
//  const Function<deal_II_dimension+1> * const coefficient);
// template
// void MatrixCreator::create_mass_matrix
// (const DoFHandler<deal_II_dimension,deal_II_dimension+1>    &dof,
//  const Quadrature<deal_II_dimension>    &q,
//  SparseMatrix<float>     &matrix,
//  const Function<deal_II_dimension+1>      &rhs,
//  Vector<double>           &rhs_vector,
//  const Function<deal_II_dimension+1> * const coefficient);


#endif

DEAL_II_NAMESPACE_CLOSE
