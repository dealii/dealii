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
void MatrixCreator::create_mass_matrix (const Mapping<dim, spacedim>       &mapping,
					const DoFHandler<dim,spacedim>    &dof,
					const Quadrature<dim>    &q,
					SparseMatrix<number>     &matrix,
					const Function<spacedim> * const coefficient)
{
  Assert (matrix.m() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
  Assert (matrix.n() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

  const unsigned int n_threads = multithread_info.n_default_threads;
  Threads::ThreadGroup<> threads;

				   // define starting and end point
				   // for each thread
  typedef typename DoFHandler<dim,spacedim>::active_cell_iterator active_cell_iterator;  
  const std::vector<std::pair<active_cell_iterator,active_cell_iterator> >
    thread_ranges = Threads::split_range<active_cell_iterator> (dof.begin_active(),
								 dof.end(), n_threads);

				   // mutex to synchronise access to
				   // the matrix
  Threads::ThreadMutex mutex;
  
				   // then assemble in parallel
  typedef void (*create_mass_matrix_1_t) (const Mapping<dim, spacedim>       &mapping,
					  const DoFHandler<dim,spacedim>    &dof,
					  const Quadrature<dim>    &q,
					  SparseMatrix<number>     &matrix,
					  const Function<spacedim> * const coefficient,
					  const IteratorRange<DoFHandler<dim,spacedim> >  range,
					  Threads::ThreadMutex     &mutex);
  create_mass_matrix_1_t p = &MatrixCreator::template create_mass_matrix_1<dim,number,spacedim>;
  for (unsigned int thread=0; thread<n_threads; ++thread)
    threads += Threads::spawn (p)(mapping, dof, q, matrix, coefficient,
                                  thread_ranges[thread], mutex);
  threads.join_all ();  
}



template <int dim, typename number, int spacedim>
void MatrixCreator::create_mass_matrix_1 (const Mapping<dim, spacedim>       &mapping,
					  const DoFHandler<dim,spacedim>    &dof,
					  const Quadrature<dim>    &q,
					  SparseMatrix<number>     &matrix,
					  const Function<spacedim> * const coefficient,
					  const IteratorRange<DoFHandler<dim,spacedim> >  range,
					  Threads::ThreadMutex     &mutex)
{
  UpdateFlags update_flags = UpdateFlags(update_values | update_JxW_values);
  if (coefficient != 0)
    update_flags = UpdateFlags (update_flags | update_quadrature_points);
  
  FEValues<dim,spacedim> fe_values (mapping, dof.get_fe(), q, update_flags);
    
  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points;
  const FiniteElement<dim,spacedim>    &fe  = fe_values.get_fe();
  const unsigned int n_components  = fe.n_components();

  Assert(coefficient == 0 ||
	 coefficient->n_components==1 ||
	 coefficient->n_components==n_components, ExcComponentMismatch());
  
  FullMatrix<number>  cell_matrix (dofs_per_cell, dofs_per_cell);
  std::vector<double> coefficient_values (n_q_points);
  std::vector<Vector<double> > coefficient_vector_values (n_q_points,
							  Vector<double> (n_components));
  
  std::vector<unsigned int> dof_indices (dofs_per_cell);
  
  typename DoFHandler<dim,spacedim>::active_cell_iterator cell = range.first;
  for (; cell!=range.second; ++cell)
    {
      fe_values.reinit (cell);
      
      cell_matrix = 0;
      cell->get_dof_indices (dof_indices);
      
      if (coefficient != 0)
	{
	  if (coefficient->n_components==1)
	    {
	      // Version for variable coefficient with 1 component
	      coefficient->value_list (fe_values.get_quadrature_points(),
				       coefficient_values);
	      for (unsigned int point=0; point<n_q_points; ++point)
		{
		  const double weight = fe_values.JxW(point);
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    {
		      const double v = fe_values.shape_value(i,point);
		      const unsigned int component_i = 
			fe.system_to_component_index(i).first;
		      for (unsigned int j=0; j<dofs_per_cell; ++j)
			if ((n_components==1) ||
			    (fe.system_to_component_index(j).first ==
			     component_i))
			  {
			    const double u = fe_values.shape_value(j,point);
			    cell_matrix(i,j) +=
			      (u * v * weight * coefficient_values[point]);
			  }
		    }
		}
	    }
	  else
	    {
	      coefficient->vector_value_list (fe_values.get_quadrature_points(),
					      coefficient_vector_values);
	      if (fe.is_primitive ())
		{
		  // Version for variable coefficient with multiple components
		  for (unsigned int point=0; point<n_q_points; ++point)
		    {
		      const double weight = fe_values.JxW(point);
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
			  const double v = fe_values.shape_value(i,point);
			  const unsigned int component_i=
			    fe.system_to_component_index(i).first;
			  for (unsigned int j=0; j<dofs_per_cell; ++j)
			    if ((n_components==1) ||
				(fe.system_to_component_index(j).first == 
				 component_i))
			      {
				const double u = fe_values.shape_value(j,point);
				cell_matrix(i,j) +=
				  (u * v * weight *
				   coefficient_vector_values[point](component_i));
			      }
			}
		    }
		}
	      else
		{
		  // Version for variable coefficient with multiple components and
		  // vector values FE.
		  for (unsigned int point=0; point<n_q_points; ++point)
		    {
		      const double weight = fe_values.JxW(point);
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
			  if (fe.get_nonzero_components(i)[comp_i])
			    {
			      const double v = fe_values.shape_value_component(i,point,comp_i);
			      for (unsigned int j=0; j<dofs_per_cell; ++j)
				for (unsigned int comp_j = 0; comp_j < n_components; ++comp_j)
				  if (fe.get_nonzero_components(j)[comp_j])
				    {
				      const double u = fe_values.shape_value_component(j,point,comp_j);
				      if (comp_i == comp_j)
					cell_matrix(i,j) +=
					  (u * v * weight *
					   coefficient_vector_values[point](comp_i));
				    }
			    }
		    }
		}
	    }
	}
      else
	{
	  if (fe.is_primitive ())
	    {
	      // Version for primitive FEs
	      for (unsigned int point=0; point<n_q_points; ++point)
		{
		  const double weight = fe_values.JxW(point);
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    {
		      const double v = fe_values.shape_value(i,point); 
		      for (unsigned int j=0; j<dofs_per_cell; ++j)
			if ((n_components==1) ||
			    (fe.system_to_component_index(i).first ==
			     fe.system_to_component_index(j).first))
			  {
			    const double u = fe_values.shape_value(j,point);
			    cell_matrix(i,j) += (u * v * weight);
			  }
		    }
		}
	    }
	  else
	    {
	      // Version for vector valued FEs
	      for (unsigned int point=0; point<n_q_points; ++point)
		{
		  const double weight = fe_values.JxW(point);
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
		      if (fe.get_nonzero_components(i)[comp_i])
			{
			  const double v = fe_values.shape_value_component(i,point,comp_i); 
			  for (unsigned int j=0; j<dofs_per_cell; ++j)
			    for (unsigned int comp_j = 0; comp_j < n_components; ++comp_j)
			      if (fe.get_nonzero_components(j)[comp_j])
				{
				  const double u = fe_values.shape_value_component(j,point,comp_j);
				  if (comp_i == comp_j)
				    cell_matrix(i,j) += (u * v * weight);
				}
			}
		}
	    }
	}
				       // transfer everything into the
				       // global object
      mutex.acquire ();
      matrix.add(dof_indices, cell_matrix);
      mutex.release ();
    };
}


template <int dim, typename number, int spacedim>
void MatrixCreator::create_mass_matrix (const DoFHandler<dim,spacedim>    &dof,
					const Quadrature<dim>    &q,
					SparseMatrix<number>     &matrix,
					const Function<spacedim> * const coefficient)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_mass_matrix(StaticMappingQ1<dim,spacedim>::mapping, dof, q, matrix, coefficient);
}


template <int dim, typename number, int spacedim>
void MatrixCreator::create_mass_matrix (const Mapping<dim, spacedim>       &mapping,
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
  
				   // then assemble in parallel
  typedef void (*create_mass_matrix_2_t) (const Mapping<dim, spacedim>       &mapping,
					  const DoFHandler<dim,spacedim>    &dof,
					  const Quadrature<dim>    &q,
					  SparseMatrix<number>     &matrix,
					  const Function<spacedim>      &rhs,
					  Vector<double>           &rhs_vector,
					  const Function<spacedim> * const coefficient,
					  const IteratorRange<DoFHandler<dim,spacedim> >  range,
					  Threads::ThreadMutex     &mutex);
  create_mass_matrix_2_t p = &MatrixCreator::template create_mass_matrix_2<dim,number,spacedim>;
  for (unsigned int thread=0; thread<n_threads; ++thread)
    threads += Threads::spawn (p)(mapping, dof, q, matrix, rhs,
                                  rhs_vector, coefficient,
                                  thread_ranges[thread], mutex);
  threads.join_all ();  
}



template <int dim, typename number, int spacedim>
void
MatrixCreator::create_mass_matrix_2 (const Mapping<dim, spacedim>       &mapping,
				     const DoFHandler<dim,spacedim>    &dof,
				     const Quadrature<dim>    &q,
				     SparseMatrix<number>     &matrix,
				     const Function<spacedim> &rhs,
				     Vector<double>           &rhs_vector,
				     const Function<spacedim> * const coefficient,
				     const IteratorRange<DoFHandler<dim,spacedim> >  range,
				     Threads::ThreadMutex     &mutex)
{
  UpdateFlags update_flags = UpdateFlags(update_values |
					 update_quadrature_points |
					 update_JxW_values);
  if (coefficient != 0)
    update_flags = UpdateFlags (update_flags | update_quadrature_points);

  FEValues<dim,spacedim> fe_values (mapping, dof.get_fe(), q, update_flags);
    
  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points;
  const FiniteElement<dim,spacedim>    &fe  = fe_values.get_fe();
  const unsigned int n_components  = fe.n_components();

  Assert(coefficient == 0 ||
	 coefficient->n_components==1 ||
	 coefficient->n_components==n_components, ExcComponentMismatch());
  Assert (rhs.n_components == 1 ||
	  rhs.n_components == n_components,ExcComponentMismatch());
  Assert (rhs_vector.size() == dof.n_dofs(),
	  ExcDimensionMismatch(rhs_vector.size(), dof.n_dofs()));
  
  FullMatrix<number>  cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>      local_rhs (dofs_per_cell);
  std::vector<double> coefficient_values (n_q_points);
  std::vector<Vector<double> > coefficient_vector_values (n_q_points,
						 Vector<double> (n_components));
  std::vector<double> rhs_values(n_q_points);
  std::vector<Vector<double> > rhs_vector_values(n_q_points,
						 Vector<double> (n_components));

  std::vector<unsigned int> dof_indices (dofs_per_cell);

  typename DoFHandler<dim,spacedim>::active_cell_iterator cell = range.first;

  for (; cell!=range.second; ++cell)
    {
      fe_values.reinit (cell);

      cell_matrix = 0;
      local_rhs = 0;
      cell->get_dof_indices (dof_indices);

				   // value_list for one component rhs,
				   // vector_value_list otherwise
      if (rhs.n_components==1)
	rhs.value_list (fe_values.get_quadrature_points(), rhs_values);
      else
	rhs.vector_value_list (fe_values.get_quadrature_points(),
			       rhs_vector_values);

				   // Case with coefficient
      if (coefficient != 0)
	{
	  if (coefficient->n_components == 1)
	    {
	      // Version for variable coefficient with 1 component
	      coefficient->value_list (fe_values.get_quadrature_points(),
				       coefficient_values);
	      for (unsigned int point=0; point<n_q_points; ++point)
		{
		  const double weight = fe_values.JxW(point);
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    {
		      const double v = fe_values.shape_value(i,point);
		      const unsigned int component_i =
			fe.system_to_component_index(i).first;
		      for (unsigned int j=0; j<dofs_per_cell; ++j)
			if ((n_components==1) ||
			    (fe.system_to_component_index(j).first ==
			     component_i))
			  {
			    const double u = fe_values.shape_value(j,point);
			    cell_matrix(i,j) +=
			      (u * v * weight * coefficient_values[point]);
			  }

		      if (rhs.n_components == 1)
			local_rhs(i) += v * rhs_values[point] * weight;
		      else
			local_rhs(i) += v * rhs_vector_values[point](component_i) 
			                  * weight;
		    }
		}
	    }
	  else
	    {
	      coefficient->vector_value_list (fe_values.get_quadrature_points(),
					      coefficient_vector_values);
	      if (fe.is_primitive ())
		{
		  // Version for variable coefficient with multiple components
		  for (unsigned int point=0; point<n_q_points; ++point)
		    {
		      const double weight = fe_values.JxW(point);
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			{
			  const double v = fe_values.shape_value(i,point);
			  const unsigned int component_i=
			    fe.system_to_component_index(i).first;
			  for (unsigned int j=0; j<dofs_per_cell; ++j)
			    if ((n_components==1) ||
				(fe.system_to_component_index(j).first == 
				 component_i))
			      {
				const double u = fe_values.shape_value(j,point);
				cell_matrix(i,j) +=
				  (u * v * weight *
				   coefficient_vector_values[point](component_i));
			      }

			  if (rhs.n_components == 1)
			    local_rhs(i) += v * rhs_values[point] * weight;
			  else
			    local_rhs(i) += v * rhs_vector_values[point](component_i) 
			                      * weight;
			}
		    }
		}
	      else
		{
		  // Version for variable coefficient with multiple components and
		  // vector valued FE.
		  for (unsigned int point=0; point<n_q_points; ++point)
		    {
		      const double weight = fe_values.JxW(point);
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
			  if (fe.get_nonzero_components(i)[comp_i])
			    {
			      const double v = fe_values.shape_value_component(i,point,comp_i);
			      for (unsigned int j=0; j<dofs_per_cell; ++j)
				for (unsigned int comp_j = 0; comp_j < n_components; ++comp_j)
				  if (fe.get_nonzero_components(j)[comp_j])
				    {
				      const double u = fe_values.shape_value_component(j,point,comp_j);
				      if (comp_i == comp_j)
					cell_matrix(i,j) +=
					  (u * v * weight *
					   coefficient_vector_values[point](comp_i));
				    }

			      if (rhs.n_components == 1)
				local_rhs(i) += v * rhs_values[point] * weight;
			      else
				local_rhs(i) += v * rhs_vector_values[point](comp_i) 
				                  * weight;
			    }
		    }
		}
	    }
	}
      else
	{
	  if (fe.is_primitive ())
	    {
	      // Version for primitive FEs
	      for (unsigned int point=0; point<n_q_points; ++point)
		{
		  const double weight = fe_values.JxW(point);
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    {
		      const double v = fe_values.shape_value(i,point);
		      const unsigned int component_i = 
			fe.system_to_component_index(i).first;
		      for (unsigned int j=0; j<dofs_per_cell; ++j)
			if ((n_components==1) ||
			    (fe.system_to_component_index(j).first ==
			     component_i))
			  {
			    const double u = fe_values.shape_value(j,point);
			    cell_matrix(i,j) += (u * v * weight);
			  }

		      if (rhs.n_components == 1)
			local_rhs(i) += v * rhs_values[point] * weight;
		      else
			local_rhs(i) += v * rhs_vector_values[point](component_i) 
			                  * weight;
		    }
		}
	    }
	  else
	    {
	      // Version for vector valued FEs
	      for (unsigned int point=0; point<n_q_points; ++point)
		{
		  const double weight = fe_values.JxW(point);
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
		      if (fe.get_nonzero_components(i)[comp_i])
			{
			  const double v = fe_values.shape_value_component(i,point,comp_i); 
			  for (unsigned int j=0; j<dofs_per_cell; ++j)
			    for (unsigned int comp_j = 0; comp_j < n_components; ++comp_j)
			      if (fe.get_nonzero_components(j)[comp_j])
				{
				  const double u = fe_values.shape_value_component(j,point,comp_j);
				  if (comp_i == comp_j)
				    cell_matrix(i,j) += (u * v * weight);
				}

			  if (rhs.n_components == 1)
			    local_rhs(i) += v * rhs_values[point] * weight;
			  else
			    local_rhs(i) += v * rhs_vector_values[point](comp_i) 
				              * weight;
			}
		}
	    }
	}

				       // transfer everything into the
				       // global object. lock the
				       // matrix meanwhile
      Threads::ThreadMutex::ScopedLock lock (mutex);
      matrix.add(dof_indices, cell_matrix);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	rhs_vector(dof_indices[i]) += local_rhs(i);
    }
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
void MatrixCreator::create_mass_matrix (const hp::MappingCollection<dim,spacedim>       &mapping,
					const hp::DoFHandler<dim,spacedim>    &dof,
					const hp::QCollection<dim>    &q,
					SparseMatrix<number>     &matrix,
					const Function<spacedim> * const coefficient)
{
  Assert (matrix.m() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
  Assert (matrix.n() == dof.n_dofs(),
	  ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

  const unsigned int n_threads = multithread_info.n_default_threads;
  Threads::ThreadGroup<> threads;

				   // define starting and end point
				   // for each thread
  typedef typename hp::DoFHandler<dim,spacedim>::active_cell_iterator active_cell_iterator;  
  std::vector<std::pair<active_cell_iterator,active_cell_iterator> > thread_ranges
    = Threads::split_range<active_cell_iterator> (dof.begin_active(),
						  dof.end(), n_threads);

				   // mutex to synchronise access to
				   // the matrix
  Threads::ThreadMutex mutex;
  
				   // then assemble in parallel
  typedef void (*create_mass_matrix_1_t) (const hp::MappingCollection<dim,spacedim>       &mapping,
					  const hp::DoFHandler<dim,spacedim>    &dof,
					  const hp::QCollection<dim>    &q,
					  SparseMatrix<number>     &matrix,
					  const Function<spacedim> * const coefficient,
					  const IteratorRange<hp::DoFHandler<dim,spacedim> >  range,
					  Threads::ThreadMutex     &mutex);
  create_mass_matrix_1_t p = &MatrixCreator::template create_mass_matrix_1<dim>;
  for (unsigned int thread=0; thread<n_threads; ++thread)
    threads += Threads::spawn (p)(mapping, dof, q, matrix, coefficient,
                                  thread_ranges[thread], mutex);
  threads.join_all ();  
}



template <int dim, typename number, int spacedim>
void
MatrixCreator::create_mass_matrix_1 (const hp::MappingCollection<dim,spacedim>       &mapping,
				     const hp::DoFHandler<dim,spacedim>    &dof,
				     const hp::QCollection<dim>    &q,
				     SparseMatrix<number>     &matrix,
				     const Function<spacedim> * const coefficient,
				     const IteratorRange<hp::DoFHandler<dim,spacedim> >  range,
				     Threads::ThreadMutex     &mutex)
{
  UpdateFlags update_flags = UpdateFlags(update_values |
					 update_JxW_values);
  if (coefficient != 0)
    update_flags = UpdateFlags (update_flags | update_quadrature_points);

  hp::FEValues<dim,spacedim> x_fe_values (mapping, dof.get_fe(), q, update_flags);
    
  const unsigned int n_components  = dof.get_fe().n_components();

  Assert(coefficient == 0 ||
	 coefficient->n_components==1 ||
	 coefficient->n_components==n_components, ExcComponentMismatch());

  FullMatrix<double>  cell_matrix;
  std::vector<double> coefficient_values;
  std::vector<Vector<double> > coefficient_vector_values;
  
  std::vector<unsigned int> dof_indices;
  
  typename hp::DoFHandler<dim,spacedim>::active_cell_iterator cell = range.first;
  for (; cell!=range.second; ++cell)
    {
      x_fe_values.reinit (cell);
      const FEValues<dim,spacedim> &fe_values = x_fe_values.get_present_fe_values ();

      const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
			 n_q_points    = fe_values.n_quadrature_points;
      const FiniteElement<dim,spacedim>    &fe  = fe_values.get_fe();

      cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
      coefficient_values.resize (n_q_points);
      coefficient_vector_values.resize (n_q_points,
					Vector<double> (n_components));
      dof_indices.resize (dofs_per_cell);

      
      cell_matrix = 0;
      cell->get_dof_indices (dof_indices);
      
      if (coefficient != 0)
	{
	  if (coefficient->n_components==1)
	    {
	      coefficient->value_list (fe_values.get_quadrature_points(),
				       coefficient_values);
	      for (unsigned int point=0; point<n_q_points; ++point)
		{
		  const double weight = fe_values.JxW(point);
		  for (unsigned int i=0; i<dofs_per_cell; ++i) 
		    {
		      for (unsigned int j=0; j<dofs_per_cell; ++j)
			{
			  if ((n_components==1) ||
			      (fe.system_to_component_index(i).first ==
			       fe.system_to_component_index(j).first))
			    cell_matrix(i,j) += (fe_values.shape_value(i,point) *
						 fe_values.shape_value(j,point) *
						 weight *
						 coefficient_values[point]);
			}
		    }
		}
	    }
	  else
	    {
	      coefficient->vector_value_list (fe_values.get_quadrature_points(),
					      coefficient_vector_values);
	      for (unsigned int point=0; point<n_q_points; ++point)
		{
		  const double weight = fe_values.JxW(point);
		  for (unsigned int i=0; i<dofs_per_cell; ++i) 
		    {
		      const unsigned int component_i=
			fe.system_to_component_index(i).first;		    
		      for (unsigned int j=0; j<dofs_per_cell; ++j)
			{
			  if ((n_components==1) ||
			      (fe.system_to_component_index(j).first == component_i))
			    cell_matrix(i,j) += (fe_values.shape_value(i,point) *
						 fe_values.shape_value(j,point) *
						 weight *
						 coefficient_vector_values[point](component_i));
			}
		    }
		}
	    }
	}
      else
	for (unsigned int point=0; point<n_q_points; ++point)
	  {
	    const double weight = fe_values.JxW(point);
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      {
		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  {
		    if ((n_components==1) ||
			(fe.system_to_component_index(i).first ==
			 fe.system_to_component_index(j).first))
		      cell_matrix(i,j) += (fe_values.shape_value(i,point) *
					   fe_values.shape_value(j,point) *
					   weight);
		  }
	      }
	  }

				       // transfer everything into the
				       // global object. lock the
				       // matrix meanwhile
      Threads::ThreadMutex::ScopedLock lock (mutex);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  if ((n_components==1) ||
	      (fe.system_to_component_index(i).first ==
	       fe.system_to_component_index(j).first))
	    matrix.add (dof_indices[i], dof_indices[j],
			cell_matrix(i,j));
    }
}



template <int dim, typename number, int spacedim>
void MatrixCreator::create_mass_matrix (const hp::DoFHandler<dim,spacedim>    &dof,
					const hp::QCollection<dim>    &q,
					SparseMatrix<number>     &matrix,
					const Function<spacedim> * const coefficient)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_mass_matrix(hp::StaticMappingQ1<dim>::mapping_collection, dof, q, matrix, coefficient);
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

  const unsigned int n_threads = multithread_info.n_default_threads;
  Threads::ThreadGroup<> threads;

				   // define starting and end point
				   // for each thread
  typedef typename hp::DoFHandler<dim,spacedim>::active_cell_iterator active_cell_iterator;  
  std::vector<std::pair<active_cell_iterator,active_cell_iterator> > thread_ranges
    = Threads::split_range<active_cell_iterator> (dof.begin_active(),
						  dof.end(), n_threads);

				   // mutex to synchronise access to
				   // the matrix
  Threads::ThreadMutex mutex;
  
				   // then assemble in parallel
  typedef void (*create_mass_matrix_2_t) (const hp::MappingCollection<dim,spacedim>       &mapping,
					  const hp::DoFHandler<dim,spacedim>    &dof,
					  const hp::QCollection<dim>    &q,
					  SparseMatrix<number>     &matrix,
					  const Function<spacedim>      &rhs,
					  Vector<double>           &rhs_vector,
					  const Function<spacedim> * const coefficient,
					  const IteratorRange<hp::DoFHandler<dim,spacedim> >  range,
					  Threads::ThreadMutex     &mutex);
  create_mass_matrix_2_t p = &MatrixCreator::template create_mass_matrix_2<dim>;
  for (unsigned int thread=0; thread<n_threads; ++thread)
    threads += Threads::spawn (p)(mapping, dof, q, matrix, rhs,
                                  rhs_vector, coefficient,
                                  thread_ranges[thread], mutex);
  threads.join_all ();  
}



template <int dim, typename number, int spacedim>
void
MatrixCreator::create_mass_matrix_2 (const hp::MappingCollection<dim,spacedim>       &mapping,
				     const hp::DoFHandler<dim,spacedim>    &dof,
				     const hp::QCollection<dim>    &q,
				     SparseMatrix<number>     &matrix,
				     const Function<spacedim>      &rhs,
				     Vector<double>           &rhs_vector,
				     const Function<spacedim> * const coefficient,
				     const IteratorRange<hp::DoFHandler<dim,spacedim> >  range,
				     Threads::ThreadMutex     &mutex)
{
  UpdateFlags update_flags = UpdateFlags(update_values    |
					 update_quadrature_points  |
					 update_JxW_values);
  if (coefficient != 0)
    update_flags = UpdateFlags (update_flags | update_quadrature_points);

  hp::FEValues<dim,spacedim> x_fe_values (mapping, dof.get_fe(), q, update_flags);
    
  const unsigned int n_components  = dof.get_fe().n_components();

  Assert(coefficient == 0 ||
	 coefficient->n_components==1 ||
	 coefficient->n_components==n_components, ExcComponentMismatch());

  FullMatrix<double>  cell_matrix;
  Vector<double>      local_rhs;
  std::vector<double> rhs_values;
  std::vector<double> coefficient_values;
  std::vector<Vector<double> > coefficient_vector_values;
  
  std::vector<unsigned int> dof_indices;
  
  typename hp::DoFHandler<dim,spacedim>::active_cell_iterator cell = range.first;
  for (; cell!=range.second; ++cell)
    {
      x_fe_values.reinit (cell);
      const FEValues<dim,spacedim> &fe_values = x_fe_values.get_present_fe_values ();

      const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
			 n_q_points    = fe_values.n_quadrature_points;
      const FiniteElement<dim,spacedim>    &fe  = fe_values.get_fe();

      cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
      local_rhs.reinit (dofs_per_cell);
      rhs_values.resize (fe_values.n_quadrature_points);
      coefficient_values.resize (n_q_points);
      coefficient_vector_values.resize (n_q_points,
					Vector<double> (n_components));
      dof_indices.resize (dofs_per_cell);

      
      cell_matrix = 0;
      local_rhs = 0;
      cell->get_dof_indices (dof_indices);
      
      rhs.value_list (fe_values.get_quadrature_points(), rhs_values);
      
      if (coefficient != 0)
	{
	  if (coefficient->n_components==1)
	    {
	      coefficient->value_list (fe_values.get_quadrature_points(),
				       coefficient_values);
	      for (unsigned int point=0; point<n_q_points; ++point)
		{
		  const double weight = fe_values.JxW(point);
		  for (unsigned int i=0; i<dofs_per_cell; ++i) 
		    {
		      for (unsigned int j=0; j<dofs_per_cell; ++j)
			{
			  if ((n_components==1) ||
			      (fe.system_to_component_index(i).first ==
			       fe.system_to_component_index(j).first))
			    cell_matrix(i,j) += (fe_values.shape_value(i,point) *
						 fe_values.shape_value(j,point) *
						 weight *
						 coefficient_values[point]);
			}
		      local_rhs(i) += fe_values.shape_value(i,point) *
				      rhs_values[point] * weight;
		    }
		}
	    }
	  else
	    {
	      coefficient->vector_value_list (fe_values.get_quadrature_points(),
					      coefficient_vector_values);
	      for (unsigned int point=0; point<n_q_points; ++point)
		{
		  const double weight = fe_values.JxW(point);
		  for (unsigned int i=0; i<dofs_per_cell; ++i) 
		    {
		      const unsigned int component_i=
			fe.system_to_component_index(i).first;		    
		      for (unsigned int j=0; j<dofs_per_cell; ++j)
			{
			  if ((n_components==1) ||
			      (fe.system_to_component_index(j).first == component_i))
			    cell_matrix(i,j) += (fe_values.shape_value(i,point) *
						 fe_values.shape_value(j,point) *
						 weight *
						 coefficient_vector_values[point](component_i));
			}			  
		      local_rhs(i) += fe_values.shape_value(i,point) *
				      rhs_values[point] * weight;
		    }
		}
	    }
	}
      else
	for (unsigned int point=0; point<n_q_points; ++point)
	  {
	    const double weight = fe_values.JxW(point);
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      {
		const double v = fe_values.shape_value(i,point);
		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  {
		    if ((n_components==1) ||
			(fe.system_to_component_index(i).first ==
			 fe.system_to_component_index(j).first))
		      cell_matrix(i,j) += (fe_values.shape_value(i,point) *
					   fe_values.shape_value(j,point) *
					   weight);
		  }
		local_rhs(i) += v * rhs_values[point] * weight;
	      }
	  }

				       // transfer everything into the
				       // global object. lock the
				       // matrix meanwhile
      Threads::ThreadMutex::ScopedLock lock (mutex);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  if ((n_components==1) ||
	      (fe.system_to_component_index(i).first ==
	       fe.system_to_component_index(j).first))
	    matrix.add (dof_indices[i], dof_indices[j],
			cell_matrix(i,j));
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	rhs_vector(dof_indices[i]) += local_rhs(i);
    }
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
  create_mass_matrix(hp::StaticMappingQ1<dim>::mapping_collection, dof, q,
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
  for (unsigned int thread=0; thread<n_threads; ++thread)
    threads += Threads::spawn (p)(Commons(mapping, dof, q), matrix,
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
  FEFaceValues<dim> fe_values (mapping, fe, q, update_flags);

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
  const hp::FECollection<dim> &fe_collection = dof.get_fe();
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

  typedef std_cxx1x::tuple<const hp::MappingCollection<dim>&,
    const hp::DoFHandler<dim>&,
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
  for (unsigned int thread=0; thread<n_threads; ++thread)
    threads += Threads::spawn (p)(Commons(mapping, dof, q), matrix,
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
  const hp::MappingCollection<dim>& mapping = std_cxx1x::get<0>(commons);
  const hp::DoFHandler<dim>& dof = std_cxx1x::get<1>(commons);
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

	  const FEFaceValues<dim> &fe_values = x_fe_values.get_present_fe_values ();

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
				{
				  const double u = fe_values.shape_value(j,point);
				  if (fe.system_to_component_index(i).first ==
				      fe.system_to_component_index(j).first)
				    {
				      cell_matrix(i,j)
					+= (u * v * weight * coefficient_values[point]);
				    }
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
				{
				  const double u = fe_values.shape_value(j,point);
				  if (fe.system_to_component_index(j).first ==
				      component_i)
				    {
				      cell_matrix(i,j) +=
					(u * v * weight * coefficient_vector_values[point](component_i));
				    }
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
			  {
			    const double u = fe_values.shape_value(j,point);
			    if (fe.system_to_component_index(i).first ==
				fe.system_to_component_index(j).first)
			      {
				cell_matrix(i,j) += (u * v * weight);
			      }
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
  
				   // then assemble in parallel
  typedef void (*create_laplace_matrix_1_t) (const Mapping<dim, spacedim>       &mapping,
					     const DoFHandler<dim,spacedim>    &dof,
					     const Quadrature<dim>    &q,
					     SparseMatrix<double>     &matrix,
					     const Function<spacedim> * const coefficient,
					     const IteratorRange<DoFHandler<dim,spacedim> >  range,
					     Threads::ThreadMutex     &mutex);
  create_laplace_matrix_1_t p = &MatrixCreator::template create_laplace_matrix_1<dim>;
  for (unsigned int thread=0; thread<n_threads; ++thread)
    threads += Threads::spawn (p)(mapping, dof, q, matrix, coefficient,
                                  thread_ranges[thread], mutex);
  threads.join_all ();  
}



template <int dim, int spacedim>
void MatrixCreator::create_laplace_matrix_1 (const Mapping<dim, spacedim>       &mapping,
					     const DoFHandler<dim,spacedim>    &dof,
					     const Quadrature<dim>    &q,
					     SparseMatrix<double>     &matrix,
					     const Function<spacedim> * const coefficient,
					     const IteratorRange<DoFHandler<dim,spacedim> >  range,
					     Threads::ThreadMutex     &mutex)
{
  UpdateFlags update_flags = UpdateFlags(update_JxW_values |
					 update_gradients);
  if (coefficient != 0)
    update_flags = UpdateFlags (update_flags | update_quadrature_points);

  FEValues<dim,spacedim> fe_values (mapping, dof.get_fe(), q, update_flags);
    
  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points;
  const FiniteElement<dim,spacedim>    &fe  = fe_values.get_fe();
  const unsigned int n_components  = fe.n_components();

  Assert(coefficient == 0 ||
	 coefficient->n_components==1 ||
	 coefficient->n_components==n_components, ExcComponentMismatch());

  FullMatrix<double>  cell_matrix (dofs_per_cell, dofs_per_cell);
  std::vector<double> coefficient_values (n_q_points);
  std::vector<Vector<double> > coefficient_vector_values (n_q_points,
							  Vector<double> (n_components));
  
  std::vector<unsigned int> dof_indices (dofs_per_cell);
  
  typename DoFHandler<dim,spacedim>::active_cell_iterator cell = range.first;
  for (; cell!=range.second; ++cell)
    {
      fe_values.reinit (cell);
      
      cell_matrix = 0;
      cell->get_dof_indices (dof_indices);
      
      if (coefficient != 0)
	{
	  if (coefficient->n_components==1)
	    {
	      coefficient->value_list (fe_values.get_quadrature_points(),
				       coefficient_values);
	      for (unsigned int point=0; point<n_q_points; ++point)
		{
		  const double weight = fe_values.JxW(point);
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    {
		      const Tensor<1,dim>& Dv = fe_values.shape_grad(i,point);
		      for (unsigned int j=0; j<dofs_per_cell; ++j)
			{
			  const Tensor<1,dim>& Du = fe_values.shape_grad(j,point);
			  if ((n_components==1) ||
			      (fe.system_to_component_index(i).first ==
			       fe.system_to_component_index(j).first))
			    cell_matrix(i,j) += (Du * Dv * weight *
						 coefficient_values[point]);
			}
		    }
		}
	    }
	  else
	    {
	      coefficient->vector_value_list (fe_values.get_quadrature_points(),
					      coefficient_vector_values);
	      for (unsigned int point=0; point<n_q_points; ++point)
		{
		  const double weight = fe_values.JxW(point);
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    {
		      const Tensor<1,dim>& Dv = fe_values.shape_grad(i,point);
		      const unsigned int component_i=
			fe.system_to_component_index(i).first;
		      for (unsigned int j=0; j<dofs_per_cell; ++j)
			{
			  const Tensor<1,dim>& Du = fe_values.shape_grad(j,point);
			  if ((n_components==1) ||
			      (fe.system_to_component_index(j).first == component_i))
			    cell_matrix(i,j) += (Du * Dv * weight *
						 coefficient_vector_values[point](component_i));
			  
			}
		    }
		}
	    }
	}
      else
	for (unsigned int point=0; point<n_q_points; ++point)
	  {
	    const double weight = fe_values.JxW(point);
	    for (unsigned int i=0; i<dofs_per_cell; ++i)
	      {
		const Tensor<1,dim>& Dv = fe_values.shape_grad(i,point);
		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  {
		    const Tensor<1,dim>& Du = fe_values.shape_grad(j,point);
		    if ((n_components==1) ||
			(fe.system_to_component_index(i).first ==
			 fe.system_to_component_index(j).first))
		      cell_matrix(i,j) += (Du * Dv * weight);
		  }
	      }
	  }
    
				       // transfer everything into the
				       // global object. lock the
				       // matrix meanwhile
      Threads::ThreadMutex::ScopedLock lock (mutex);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  if ((n_components==1) ||
	      (fe.system_to_component_index(i).first ==
	       fe.system_to_component_index(j).first))
	    matrix.add (dof_indices[i], dof_indices[j],
			cell_matrix(i,j));
    }
}



template <int dim, int spacedim>
void MatrixCreator::create_laplace_matrix (const DoFHandler<dim,spacedim>    &dof,
					   const Quadrature<dim>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<spacedim> * const coefficient)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_laplace_matrix(StaticMappingQ1<dim>::mapping, dof, q, matrix, coefficient);
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
  
				   // then assemble in parallel
  typedef void (*create_laplace_matrix_2_t) (const Mapping<dim, spacedim>       &mapping,
					     const DoFHandler<dim,spacedim>    &dof,
					     const Quadrature<dim>    &q,
					     SparseMatrix<double>     &matrix,
					     const Function<spacedim>      &rhs,
					     Vector<double>           &rhs_vector,
					     const Function<spacedim> * const coefficient,
					     const IteratorRange<DoFHandler<dim,spacedim> >  range,
					     Threads::ThreadMutex     &mutex);
  create_laplace_matrix_2_t p = &MatrixCreator::template create_laplace_matrix_2<dim>;
  for (unsigned int thread=0; thread<n_threads; ++thread)
    threads += Threads::spawn (p)(mapping, dof, q, matrix, rhs,
                                  rhs_vector, coefficient,
                                  thread_ranges[thread], mutex);
  threads.join_all ();  
}



template <int dim, int spacedim>
void
MatrixCreator::create_laplace_matrix_2 (const Mapping<dim, spacedim>       &mapping,
					const DoFHandler<dim,spacedim>    &dof,
					const Quadrature<dim>    &q,
					SparseMatrix<double>     &matrix,
					const Function<spacedim>      &rhs,
					Vector<double>           &rhs_vector,
					const Function<spacedim> * const coefficient,
					const IteratorRange<DoFHandler<dim,spacedim> >  range,
					Threads::ThreadMutex     &mutex)
{
  UpdateFlags update_flags = UpdateFlags(update_values    |
					 update_gradients |
					 update_quadrature_points  |
					 update_JxW_values);
  if (coefficient != 0)
    update_flags = UpdateFlags (update_flags | update_quadrature_points);

  FEValues<dim,spacedim> fe_values (mapping, dof.get_fe(), q, update_flags);
    
  const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
		     n_q_points    = fe_values.n_quadrature_points;
  const FiniteElement<dim,spacedim>    &fe  = fe_values.get_fe();
  const unsigned int n_components  = fe.n_components();

  Assert(coefficient == 0 ||
	 coefficient->n_components==1 ||
	 coefficient->n_components==n_components, ExcComponentMismatch());

  FullMatrix<double>  cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>      local_rhs (dofs_per_cell);
  std::vector<double> rhs_values (fe_values.n_quadrature_points);
  std::vector<double> coefficient_values (n_q_points);
  std::vector<Vector<double> > coefficient_vector_values (n_q_points,
							  Vector<double> (n_components));
  
  std::vector<unsigned int> dof_indices (dofs_per_cell);
  
  typename DoFHandler<dim,spacedim>::active_cell_iterator cell = range.first;
  for (; cell!=range.second; ++cell)
    {
      fe_values.reinit (cell);
      
      cell_matrix = 0;
      local_rhs = 0;
      cell->get_dof_indices (dof_indices);
      
      rhs.value_list (fe_values.get_quadrature_points(), rhs_values);
      
      if (coefficient != 0)
	{
	  if (coefficient->n_components==1)
	    {
	      coefficient->value_list (fe_values.get_quadrature_points(),
				       coefficient_values);
	      for (unsigned int point=0; point<n_q_points; ++point)
		  {
		    const double weight = fe_values.JxW(point);
		    for (unsigned int i=0; i<dofs_per_cell; ++i) 
		      {
			const double v = fe_values.shape_value(i,point);
			const Tensor<1,dim>& Dv = fe_values.shape_grad(i,point);
			for (unsigned int j=0; j<dofs_per_cell; ++j)
			  {
			    const Tensor<1,dim>& Du = fe_values.shape_grad(j,point);
			    if ((n_components==1) ||
				(fe.system_to_component_index(i).first ==
				 fe.system_to_component_index(j).first))
			      cell_matrix(i,j) += (Du * Dv * weight *
						   coefficient_values[point]);
			  }
			local_rhs(i) += v * rhs_values[point] * weight;  
		      }
		  }
	    }
	  else
	    {
	      coefficient->vector_value_list (fe_values.get_quadrature_points(),
					      coefficient_vector_values);
	      for (unsigned int point=0; point<n_q_points; ++point)
		  {
		    const double weight = fe_values.JxW(point);
		    for (unsigned int i=0; i<dofs_per_cell; ++i) 
		      {
			const double v = fe_values.shape_value(i,point);
			const Tensor<1,dim>& Dv = fe_values.shape_grad(i,point);
			const unsigned int component_i=
			  fe.system_to_component_index(i).first;		    
			for (unsigned int j=0; j<dofs_per_cell; ++j)
			  {
			    const Tensor<1,dim>& Du = fe_values.shape_grad(j,point);
			    if ((n_components==1) ||
				(fe.system_to_component_index(j).first == component_i))
			      cell_matrix(i,j) += (Du * Dv * weight *
						   coefficient_vector_values[point](component_i));
			  }
			local_rhs(i) += v * rhs_values[point] * weight;
		      }
		  }
	    }
	}
      else
	for (unsigned int point=0; point<n_q_points; ++point)
	  {
	    const double weight = fe_values.JxW(point);
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      {
		const double v = fe_values.shape_value(i,point);
		const Tensor<1,dim>& Dv = fe_values.shape_grad(i,point);
		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  {
		    const Tensor<1,dim>& Du = fe_values.shape_grad(j,point);
		    if ((n_components==1) ||
			(fe.system_to_component_index(i).first ==
			 fe.system_to_component_index(j).first))
		      cell_matrix(i,j) += (Du * Dv * weight);
		  }
		local_rhs(i) += v * rhs_values[point] * weight;
	      }
	  }

				       // transfer everything into the
				       // global object. lock the
				       // matrix meanwhile
      Threads::ThreadMutex::ScopedLock lock (mutex);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  if ((n_components==1) ||
	      (fe.system_to_component_index(i).first ==
	       fe.system_to_component_index(j).first))
	    matrix.add (dof_indices[i], dof_indices[j],
			cell_matrix(i,j));
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	rhs_vector(dof_indices[i]) += local_rhs(i);
    };
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
  create_laplace_matrix(StaticMappingQ1<dim>::mapping, dof, q,
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

  const unsigned int n_threads = multithread_info.n_default_threads;
  Threads::ThreadGroup<> threads;

				   // define starting and end point
				   // for each thread
  typedef typename hp::DoFHandler<dim,spacedim>::active_cell_iterator active_cell_iterator;  
  std::vector<std::pair<active_cell_iterator,active_cell_iterator> > thread_ranges
    = Threads::split_range<active_cell_iterator> (dof.begin_active(),
						  dof.end(), n_threads);

				   // mutex to synchronise access to
				   // the matrix
  Threads::ThreadMutex mutex;
  
				   // then assemble in parallel
  typedef void (*create_laplace_matrix_1_t) (const hp::MappingCollection<dim,spacedim>       &mapping,
					     const hp::DoFHandler<dim,spacedim>    &dof,
					     const hp::QCollection<dim>    &q,
					     SparseMatrix<double>     &matrix,
					     const Function<spacedim> * const coefficient,
					     const IteratorRange<hp::DoFHandler<dim,spacedim> >  range,
					     Threads::ThreadMutex     &mutex);
  create_laplace_matrix_1_t p = &MatrixCreator::template create_laplace_matrix_1<dim>;
  for (unsigned int thread=0; thread<n_threads; ++thread)
    threads += Threads::spawn (p)(mapping, dof, q, matrix, coefficient,
                                  thread_ranges[thread], mutex);
  threads.join_all ();  
}



template <int dim, int spacedim>
void
MatrixCreator::create_laplace_matrix_1 (const hp::MappingCollection<dim,spacedim>       &mapping,
					const hp::DoFHandler<dim,spacedim>    &dof,
					const hp::QCollection<dim>    &q,
					SparseMatrix<double>     &matrix,
					const Function<spacedim> * const coefficient,
					const IteratorRange<hp::DoFHandler<dim,spacedim> >  range,
					Threads::ThreadMutex     &mutex)
{
  UpdateFlags update_flags = UpdateFlags(update_gradients |
					 update_JxW_values);
  if (coefficient != 0)
    update_flags = UpdateFlags (update_flags | update_quadrature_points);

  hp::FEValues<dim,spacedim> x_fe_values (mapping, dof.get_fe(), q, update_flags);
    
  const unsigned int n_components  = dof.get_fe().n_components();

  Assert(coefficient == 0 ||
	 coefficient->n_components==1 ||
	 coefficient->n_components==n_components, ExcComponentMismatch());

  FullMatrix<double>  cell_matrix;
  std::vector<double> coefficient_values;
  std::vector<Vector<double> > coefficient_vector_values;
  
  std::vector<unsigned int> dof_indices;
  
  typename hp::DoFHandler<dim,spacedim>::active_cell_iterator cell = range.first;
  for (; cell!=range.second; ++cell)
    {
      x_fe_values.reinit (cell);
      const FEValues<dim,spacedim> &fe_values = x_fe_values.get_present_fe_values ();

      const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
			 n_q_points    = fe_values.n_quadrature_points;
      const FiniteElement<dim,spacedim>    &fe  = fe_values.get_fe();

      cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
      coefficient_values.resize (n_q_points);
      coefficient_vector_values.resize (n_q_points,
					Vector<double> (n_components));
      dof_indices.resize (dofs_per_cell);

      
      cell_matrix = 0;
      cell->get_dof_indices (dof_indices);
      
      if (coefficient != 0)
	{
	  if (coefficient->n_components==1)
	    {
	      coefficient->value_list (fe_values.get_quadrature_points(),
				       coefficient_values);
	      for (unsigned int point=0; point<n_q_points; ++point)
		  {
		    const double weight = fe_values.JxW(point);
		    for (unsigned int i=0; i<dofs_per_cell; ++i) 
		      {
			const Tensor<1,dim>& Dv = fe_values.shape_grad(i,point);
			for (unsigned int j=0; j<dofs_per_cell; ++j)
			  {
			    const Tensor<1,dim>& Du = fe_values.shape_grad(j,point);
			    if ((n_components==1) ||
				(fe.system_to_component_index(i).first ==
				 fe.system_to_component_index(j).first))
			      cell_matrix(i,j) += (Du * Dv * weight *
						   coefficient_values[point]);
			  }
		      }
		  }
	    }
	  else
	    {
	      coefficient->vector_value_list (fe_values.get_quadrature_points(),
					      coefficient_vector_values);
	      for (unsigned int point=0; point<n_q_points; ++point)
		  {
		    const double weight = fe_values.JxW(point);
		    for (unsigned int i=0; i<dofs_per_cell; ++i) 
		      {
			const Tensor<1,dim>& Dv = fe_values.shape_grad(i,point);
			const unsigned int component_i=
			  fe.system_to_component_index(i).first;		    
			for (unsigned int j=0; j<dofs_per_cell; ++j)
			  {
			    const Tensor<1,dim>& Du = fe_values.shape_grad(j,point);
			    if ((n_components==1) ||
				(fe.system_to_component_index(j).first == component_i))
			      cell_matrix(i,j) += (Du * Dv * weight *
						   coefficient_vector_values[point](component_i));
			  }
		      }
		  }
	    }
	}
      else
	for (unsigned int point=0; point<n_q_points; ++point)
	  {
	    const double weight = fe_values.JxW(point);
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      {
		const Tensor<1,dim>& Dv = fe_values.shape_grad(i,point);
		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  {
		    const Tensor<1,dim>& Du = fe_values.shape_grad(j,point);
		    if ((n_components==1) ||
			(fe.system_to_component_index(i).first ==
			 fe.system_to_component_index(j).first))
		      cell_matrix(i,j) += (Du * Dv * weight);
		  }
	      }
	  }

				       // transfer everything into the
				       // global object. lock the
				       // matrix meanwhile
      Threads::ThreadMutex::ScopedLock lock (mutex);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  if ((n_components==1) ||
	      (fe.system_to_component_index(i).first ==
	       fe.system_to_component_index(j).first))
	    matrix.add (dof_indices[i], dof_indices[j],
			cell_matrix(i,j));
    }
}



template <int dim, int spacedim>
void MatrixCreator::create_laplace_matrix (const hp::DoFHandler<dim,spacedim>    &dof,
					   const hp::QCollection<dim>    &q,
					   SparseMatrix<double>     &matrix,
					   const Function<spacedim> * const coefficient)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  create_laplace_matrix(hp::StaticMappingQ1<dim>::mapping_collection, dof, q, matrix, coefficient);
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

  const unsigned int n_threads = multithread_info.n_default_threads;
  Threads::ThreadGroup<> threads;

				   // define starting and end point
				   // for each thread
  typedef typename hp::DoFHandler<dim,spacedim>::active_cell_iterator active_cell_iterator;  
  std::vector<std::pair<active_cell_iterator,active_cell_iterator> > thread_ranges
    = Threads::split_range<active_cell_iterator> (dof.begin_active(),
						  dof.end(), n_threads);

				   // mutex to synchronise access to
				   // the matrix
  Threads::ThreadMutex mutex;
  
				   // then assemble in parallel
  typedef void (*create_laplace_matrix_2_t) (const hp::MappingCollection<dim,spacedim>       &mapping,
					     const hp::DoFHandler<dim,spacedim>    &dof,
					     const hp::QCollection<dim>    &q,
					     SparseMatrix<double>     &matrix,
					     const Function<spacedim>      &rhs,
					     Vector<double>           &rhs_vector,
					     const Function<spacedim> * const coefficient,
					     const IteratorRange<hp::DoFHandler<dim,spacedim> >  range,
					     Threads::ThreadMutex     &mutex);
  create_laplace_matrix_2_t p = &MatrixCreator::template create_laplace_matrix_2<dim>;
  for (unsigned int thread=0; thread<n_threads; ++thread)
    threads += Threads::spawn (p)(mapping, dof, q, matrix, rhs,
                                  rhs_vector, coefficient,
                                  thread_ranges[thread], mutex);
  threads.join_all ();  
}



template <int dim, int spacedim>
void
MatrixCreator::create_laplace_matrix_2 (const hp::MappingCollection<dim,spacedim>       &mapping,
					const hp::DoFHandler<dim,spacedim>    &dof,
					const hp::QCollection<dim>    &q,
					SparseMatrix<double>     &matrix,
					const Function<spacedim>      &rhs,
					Vector<double>           &rhs_vector,
					const Function<spacedim> * const coefficient,
					const IteratorRange<hp::DoFHandler<dim,spacedim> >  range,
					Threads::ThreadMutex     &mutex)
{
  UpdateFlags update_flags = UpdateFlags(update_values    |
					 update_gradients |
					 update_quadrature_points  |
					 update_JxW_values);
  if (coefficient != 0)
    update_flags = UpdateFlags (update_flags | update_quadrature_points);

  hp::FEValues<dim,spacedim> x_fe_values (mapping, dof.get_fe(), q, update_flags);
    
  const unsigned int n_components  = dof.get_fe().n_components();

  Assert(coefficient == 0 ||
	 coefficient->n_components==1 ||
	 coefficient->n_components==n_components, ExcComponentMismatch());

  FullMatrix<double>  cell_matrix;
  Vector<double>      local_rhs;
  std::vector<double> rhs_values;
  std::vector<double> coefficient_values;
  std::vector<Vector<double> > coefficient_vector_values;
  
  std::vector<unsigned int> dof_indices;
  
  typename hp::DoFHandler<dim,spacedim>::active_cell_iterator cell = range.first;
  for (; cell!=range.second; ++cell)
    {
      x_fe_values.reinit (cell);
      const FEValues<dim,spacedim> &fe_values = x_fe_values.get_present_fe_values ();

      const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
			 n_q_points    = fe_values.n_quadrature_points;
      const FiniteElement<dim,spacedim>    &fe  = fe_values.get_fe();

      cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
      local_rhs.reinit (dofs_per_cell);
      rhs_values.resize (fe_values.n_quadrature_points);
      coefficient_values.resize (n_q_points);
      coefficient_vector_values.resize (n_q_points,
					Vector<double> (n_components));
      dof_indices.resize (dofs_per_cell);

      
      cell_matrix = 0;
      local_rhs = 0;
      cell->get_dof_indices (dof_indices);
      
      rhs.value_list (fe_values.get_quadrature_points(), rhs_values);
      
      if (coefficient != 0)
	{
	  if (coefficient->n_components==1)
	    {
	      coefficient->value_list (fe_values.get_quadrature_points(),
				       coefficient_values);
	      for (unsigned int point=0; point<n_q_points; ++point)
		  {
		    const double weight = fe_values.JxW(point);
		    for (unsigned int i=0; i<dofs_per_cell; ++i) 
		      {
			const double v = fe_values.shape_value(i,point);
			const Tensor<1,dim>& Dv = fe_values.shape_grad(i,point);
			for (unsigned int j=0; j<dofs_per_cell; ++j)
			  {
			    const Tensor<1,dim>& Du = fe_values.shape_grad(j,point);
			    if ((n_components==1) ||
				(fe.system_to_component_index(i).first ==
				 fe.system_to_component_index(j).first))
			      cell_matrix(i,j) += (Du * Dv * weight *
						   coefficient_values[point]);
			  }
			local_rhs(i) += v * rhs_values[point] * weight;
		      }
		  }
	    }
	  else
	    {
	      coefficient->vector_value_list (fe_values.get_quadrature_points(),
					      coefficient_vector_values);
	      for (unsigned int point=0; point<n_q_points; ++point)
		  {
		    const double weight = fe_values.JxW(point);
		    for (unsigned int i=0; i<dofs_per_cell; ++i) 
		      {
			const double v = fe_values.shape_value(i,point);
			const Tensor<1,dim>& Dv = fe_values.shape_grad(i,point);
			const unsigned int component_i=
			  fe.system_to_component_index(i).first;		    
			for (unsigned int j=0; j<dofs_per_cell; ++j)
			  {
			    const Tensor<1,dim>& Du = fe_values.shape_grad(j,point);
			    if ((n_components==1) ||
				(fe.system_to_component_index(j).first == component_i))
			      cell_matrix(i,j) += (Du * Dv * weight *
						   coefficient_vector_values[point](component_i));
			  }
			local_rhs(i) += v * rhs_values[point] * weight;
		      }
		  }
	    }
	}
      else
	for (unsigned int point=0; point<n_q_points; ++point)
	  {
	    const double weight = fe_values.JxW(point);
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      {
		const double v = fe_values.shape_value(i,point);
		const Tensor<1,dim>& Dv = fe_values.shape_grad(i,point);
		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  {
		    const Tensor<1,dim>& Du = fe_values.shape_grad(j,point);
		    if ((n_components==1) ||
			(fe.system_to_component_index(i).first ==
			 fe.system_to_component_index(j).first))
		      cell_matrix(i,j) += (Du * Dv * weight);
		  }
		local_rhs(i) += v * rhs_values[point] * weight;
	      }
	  }

				       // transfer everything into the
				       // global object. lock the
				       // matrix meanwhile
      Threads::ThreadMutex::ScopedLock lock (mutex);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  if ((n_components==1) ||
	      (fe.system_to_component_index(i).first ==
	       fe.system_to_component_index(j).first))
	    matrix.add (dof_indices[i], dof_indices[j],
			cell_matrix(i,j));
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	rhs_vector(dof_indices[i]) += local_rhs(i);
    }
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
  create_laplace_matrix(hp::StaticMappingQ1<dim>::mapping_collection, dof, q,
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
