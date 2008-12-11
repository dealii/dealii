//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/memory_consumption.h>
#include <base/quadrature.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/petsc_vector.h>
#include <lac/petsc_block_vector.h>
#include <lac/trilinos_vector.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_boundary.h>
#include <dofs/dof_accessor.h>
#include <fe/mapping_q1.h>
#include <fe/fe_values.h>
#include <fe/fe.h>

#include <iomanip>

DEAL_II_NAMESPACE_OPEN


namespace
{
  template <int dim, int spacedim>
  inline
  std::vector<unsigned int>
  make_shape_function_to_row_table (const FiniteElement<dim,spacedim> &fe)
  {
    std::vector<unsigned int> shape_function_to_row_table (fe.dofs_per_cell);
    unsigned int row = 0;
    for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
      {
	shape_function_to_row_table[i] = row;
	row += fe.n_nonzero_components (i);
      }

    return shape_function_to_row_table;
  }
}



namespace FEValuesViews
{
  template <int dim, int spacedim>  
  Scalar<dim,spacedim>::Scalar (const FEValuesBase<dim,spacedim> &fe_values,
		       const unsigned int       component)
		  :
		  fe_values (fe_values),
		  component (component),
		  is_nonzero_shape_function_component (fe_values.fe->dofs_per_cell),
		  row_index (fe_values.fe->dofs_per_cell)
  {
    Assert (component < fe_values.fe->n_components(),
	    ExcIndexRange(component, 0, fe_values.fe->n_components()));

    const std::vector<unsigned int> shape_function_to_row_table
      = make_shape_function_to_row_table (*fe_values.fe);
    
    for (unsigned int i=0; i<fe_values.fe->dofs_per_cell; ++i)
      {
	const bool is_primitive = (fe_values.fe->is_primitive() ||
				   fe_values.fe->is_primitive(i));

	if (is_primitive == true)
	  is_nonzero_shape_function_component[i]
	    = (component ==
	       fe_values.fe->system_to_component_index(i).first);
	else
	  is_nonzero_shape_function_component[i]
	    = (fe_values.fe->get_nonzero_components(i)[component]
	       == true);

	if (is_nonzero_shape_function_component[i] == true)
	  {
	    if (is_primitive == true)
	      row_index[i] = shape_function_to_row_table[i];
	    else
	      row_index[i] = (shape_function_to_row_table[i]
			      +
			      std::count (fe_values.fe->get_nonzero_components(i).begin(),
					  fe_values.fe->get_nonzero_components(i).begin()+
					  component,
					  true));
	  }
	else
	  row_index[i] = numbers::invalid_unsigned_int;
      }
  }



  template <int dim, int spacedim>
  Scalar<dim,spacedim>::Scalar ()
		  :
		  fe_values (*static_cast<dealii::FEValuesBase<dim,spacedim>*>(0)),
		  component (numbers::invalid_unsigned_int)
  {}


  template <int dim, int spacedim>
  Scalar<dim,spacedim> &
  Scalar<dim,spacedim>::operator= (const Scalar<dim,spacedim> &)
  {
				     // we shouldn't be copying these objects
    Assert (false, ExcInternalError());
    return *this;
  }
  


  template <int dim, int spacedim>  
  Vector<dim,spacedim>::Vector (const FEValuesBase<dim,spacedim> &fe_values,
				const unsigned int       first_vector_component)
		  :
		  fe_values (fe_values),
		  first_vector_component (first_vector_component),
		  is_nonzero_shape_function_component (fe_values.fe->dofs_per_cell,
						       dim),
		  row_index (fe_values.fe->dofs_per_cell,
			     dim)
  {
    Assert (first_vector_component+dim-1 < fe_values.fe->n_components(),
	    ExcIndexRange(first_vector_component+dim-1, 0,
			  fe_values.fe->n_components()));

    const std::vector<unsigned int> shape_function_to_row_table
      = make_shape_function_to_row_table (*fe_values.fe);
    
    for (unsigned int d=0; d<dim; ++d)
      {
	const unsigned int component = first_vector_component + d;
	
	for (unsigned int i=0; i<fe_values.fe->dofs_per_cell; ++i)
	  {
	    const bool is_primitive = (fe_values.fe->is_primitive() ||
				       fe_values.fe->is_primitive(i));

	    if (is_primitive == true)
	      is_nonzero_shape_function_component[i][d]
		= (component ==
		   fe_values.fe->system_to_component_index(i).first);
	    else
	      is_nonzero_shape_function_component[i][d]
		= (fe_values.fe->get_nonzero_components(i)[component]
		   == true);

	    if (is_nonzero_shape_function_component[i][d] == true)
	      {
		if (is_primitive == true)
		  row_index[i][d] = shape_function_to_row_table[i];
		else
		  row_index[i][d] = (shape_function_to_row_table[i]
				     +
				     std::count (fe_values.fe->get_nonzero_components(i).begin(),
						 fe_values.fe->get_nonzero_components(i).begin()+
						 component,
						 true));
	      }
	    else
	      row_index[i][d] = numbers::invalid_unsigned_int;
	  }
      }
  }
  

  template <int dim, int spacedim>
  Vector<dim,spacedim>::Vector ()
		  :
		  fe_values (*static_cast<dealii::FEValuesBase<dim,spacedim>*>(0)),
		  first_vector_component (numbers::invalid_unsigned_int)
  {}



  template <int dim, int spacedim>
  Vector<dim,spacedim> &
  Vector<dim,spacedim>::operator= (const Vector<dim,spacedim> &)
  {
				     // we shouldn't be copying these objects
    Assert (false, ExcInternalError());
    return *this;
  }
}


namespace internal
{
  namespace FEValuesViews
  {
    template <int dim, int spacedim>
    Cache<dim,spacedim>::Cache (const FEValuesBase<dim,spacedim> &fe_values)
    {
      const FiniteElement<dim,spacedim> &fe = fe_values.get_fe();

				       // create the views objects. allocate a
				       // bunch of default-constructed ones
				       // then destroy them again and do
				       // in-place construction of those we
				       // actually want to use (copying stuff
				       // is wasteful and we can't do that
				       // anyway because the class has
				       // reference members)
      const unsigned int n_scalars = fe.n_components();
      scalars.resize (n_scalars);
      for (unsigned int component=0; component<n_scalars; ++component)
	{
	  scalars[component].
	    dealii::FEValuesViews::Scalar<dim,spacedim>::~Scalar ();
	  new (&scalars[component])
	    dealii::FEValuesViews::Scalar<dim,spacedim>(fe_values,
							component);
	}
      
      const unsigned int n_vectors = (fe.n_components() >= dim ?
				      fe.n_components()-dim+1 :
				      0);
      vectors.resize (n_vectors);
      for (unsigned int component=0; component<n_vectors; ++component)
	{
	  vectors[component].
	    dealii::FEValuesViews::Vector<dim,spacedim>::~Vector ();
	  new (&vectors[component])
	    dealii::FEValuesViews::Vector<dim,spacedim>(fe_values,
							component);
	}	
    }
  }
}


/* ---------------- FEValuesBase<dim,spacedim>::CellIteratorBase --------- */

template <int dim, int spacedim>
FEValuesBase<dim,spacedim>::CellIteratorBase::~CellIteratorBase ()
{}



/* ---------------- FEValuesBase<dim,spacedim>::CellIterator<CI> --------- */


template <int dim, int spacedim>
template <typename CI>
FEValuesBase<dim,spacedim>::CellIterator<CI>::CellIterator (const CI &cell)
                :
                cell(cell)
{}



template <int dim, int spacedim>
template <typename CI>
FEValuesBase<dim,spacedim>::CellIterator<CI>::
operator const typename Triangulation<dim,spacedim>::cell_iterator () const
{
  return cell;
}



template <int dim, int spacedim>
template <typename CI>
unsigned int
FEValuesBase<dim,spacedim>::CellIterator<CI>::n_dofs_for_dof_handler () const
{
  return cell->get_dof_handler().n_dofs();
}



template <int dim, int spacedim>
template <typename CI>
void
FEValuesBase<dim,spacedim>::CellIterator<CI>::
get_interpolated_dof_values (const Vector<double> &in,
                             Vector<double>       &out) const
{
  cell->get_interpolated_dof_values (in, out);
}


template <int dim, int spacedim>
template <typename CI>
void
FEValuesBase<dim,spacedim>::CellIterator<CI>::
get_interpolated_dof_values (const Vector<float> &in,
                             Vector<float>       &out) const
{
  cell->get_interpolated_dof_values (in, out);
}


template <int dim, int spacedim>
template <typename CI>
void
FEValuesBase<dim,spacedim>::CellIterator<CI>::
get_interpolated_dof_values (const Vector<long double> &in,
                             Vector<long double>       &out) const
{
  cell->get_interpolated_dof_values (in, out);
}


template <int dim, int spacedim>
template <typename CI>
void
FEValuesBase<dim,spacedim>::CellIterator<CI>::
get_interpolated_dof_values (const BlockVector<double> &in,
                             Vector<double>            &out) const
{
  cell->get_interpolated_dof_values (in, out);
}


template <int dim, int spacedim>
template <typename CI>
void
FEValuesBase<dim,spacedim>::CellIterator<CI>::
get_interpolated_dof_values (const BlockVector<float> &in,
                             Vector<float>            &out) const
{
  cell->get_interpolated_dof_values (in, out);
}


template <int dim, int spacedim>
template <typename CI>
void
FEValuesBase<dim,spacedim>::CellIterator<CI>::
get_interpolated_dof_values (const BlockVector<long double> &in,
                             Vector<long double>            &out) const
{
  cell->get_interpolated_dof_values (in, out);
}



#ifdef DEAL_II_USE_PETSC

template <int dim, int spacedim>
template <typename CI>
void
FEValuesBase<dim,spacedim>::CellIterator<CI>::
get_interpolated_dof_values (const PETScWrappers::Vector &in,
                             Vector<PetscScalar>         &out) const
{
  cell->get_interpolated_dof_values (in, out);
}



template <int dim, int spacedim>
template <typename CI>
void
FEValuesBase<dim,spacedim>::CellIterator<CI>::
get_interpolated_dof_values (const PETScWrappers::BlockVector &in,
                             Vector<PetscScalar>              &out) const
{
  cell->get_interpolated_dof_values (in, out);
}

#endif

#ifdef DEAL_II_USE_TRILINOS


template <int dim, int spacedim>
template <typename CI>
void
FEValuesBase<dim,spacedim>::CellIterator<CI>::
get_interpolated_dof_values (const TrilinosWrappers::Vector &in,
                             Vector<TrilinosScalar>         &out) const
{
  cell->get_interpolated_dof_values (in, out);
}



template <int dim, int spacedim>
template <typename CI>
void
FEValuesBase<dim,spacedim>::CellIterator<CI>::
get_interpolated_dof_values (const TrilinosWrappers::BlockVector &in,
                             Vector<TrilinosScalar>              &out) const
{
  cell->get_interpolated_dof_values (in, out);
}



template <int dim, int spacedim>
template <typename CI>
void
FEValuesBase<dim, spacedim>::CellIterator<CI>::
get_interpolated_dof_values (const TrilinosWrappers::MPI::Vector &in,
                             Vector<TrilinosScalar>              &out) const
{
  cell->get_interpolated_dof_values (in, out);
}



template <int dim, int spacedim>
template <typename CI>
void
FEValuesBase<dim, spacedim>::CellIterator<CI>::
get_interpolated_dof_values (const TrilinosWrappers::MPI::BlockVector &in,
                             Vector<TrilinosScalar>                   &out) const
{
  cell->get_interpolated_dof_values (in, out);
}

#endif


/* ---------------- FEValuesBase<dim,spacedim>::TriaCellIterator --------- */

template <int dim, int spacedim>
const char * const
FEValuesBase<dim,spacedim>::TriaCellIterator::message_string
= ("You have previously called the FEValues::reinit function with a\n"
   "cell iterator of type Triangulation<dim,spacedim>::cell_iterator. However,\n"
   "when you do this, you cannot call some functions in the FEValues\n"
   "class, such as the get_function_values/gradients/hessians\n"
   "functions. If you need these functions, then you need to call\n"
   "FEValues::reinit with an iterator type that allows to extract\n"
   "degrees of freedom, such as DoFHandler<dim,spacedim>::cell_iterator.");


template <int dim, int spacedim>
FEValuesBase<dim,spacedim>::TriaCellIterator::
TriaCellIterator (const typename Triangulation<dim,spacedim>::cell_iterator &cell)
                :
                cell(cell)
{}



template <int dim, int spacedim>
FEValuesBase<dim,spacedim>::TriaCellIterator::
operator const typename Triangulation<dim,spacedim>::cell_iterator () const
{
  return cell;
}



template <int dim, int spacedim>
unsigned int
FEValuesBase<dim,spacedim>::TriaCellIterator::n_dofs_for_dof_handler () const
{
  Assert (false, ExcMessage (message_string));
  return 0;
}


template <int dim, int spacedim>
void
FEValuesBase<dim,spacedim>::TriaCellIterator::
get_interpolated_dof_values (const Vector<double> &,
                             Vector<double>       &) const
{
  Assert (false, ExcMessage (message_string));
}


template <int dim, int spacedim>
void
FEValuesBase<dim,spacedim>::TriaCellIterator::
get_interpolated_dof_values (const Vector<float> &,
                             Vector<float>       &) const
{
  Assert (false, ExcMessage (message_string));
}


template <int dim, int spacedim>
void
FEValuesBase<dim,spacedim>::TriaCellIterator::
get_interpolated_dof_values (const Vector<long double> &,
                             Vector<long double>       &) const
{
  Assert (false, ExcMessage (message_string));
}


template <int dim, int spacedim>
void
FEValuesBase<dim,spacedim>::TriaCellIterator::
get_interpolated_dof_values (const BlockVector<double> &,
                             Vector<double>            &) const
{
  Assert (false, ExcMessage (message_string));
}



template <int dim, int spacedim>
void
FEValuesBase<dim,spacedim>::TriaCellIterator::
get_interpolated_dof_values (const BlockVector<float> &,
                             Vector<float>            &) const
{
  Assert (false, ExcMessage (message_string));
}


template <int dim, int spacedim>
void
FEValuesBase<dim,spacedim>::TriaCellIterator::
get_interpolated_dof_values (const BlockVector<long double> &,
                             Vector<long double>            &) const
{
  Assert (false, ExcMessage (message_string));
}



#ifdef DEAL_II_USE_PETSC

template <int dim, int spacedim>
void
FEValuesBase<dim,spacedim>::TriaCellIterator::
get_interpolated_dof_values (const PETScWrappers::Vector &,
                             Vector<PetscScalar>         &) const
{
  Assert (false, ExcMessage (message_string));
}



template <int dim, int spacedim>
void
FEValuesBase<dim,spacedim>::TriaCellIterator::
get_interpolated_dof_values (const PETScWrappers::BlockVector &,
                             Vector<PetscScalar>              &) const
{
  Assert (false, ExcMessage (message_string));
}

#endif

#ifdef DEAL_II_USE_TRILINOS

template <int dim, int spacedim>
void
FEValuesBase<dim,spacedim>::TriaCellIterator::
get_interpolated_dof_values (const TrilinosWrappers::Vector &,
                             Vector<TrilinosScalar>         &) const
{
  Assert (false, ExcMessage (message_string));
}



template <int dim, int spacedim>
void
FEValuesBase<dim,spacedim>::TriaCellIterator::
get_interpolated_dof_values (const TrilinosWrappers::BlockVector &,
                             Vector<TrilinosScalar>              &) const
{
  Assert (false, ExcMessage (message_string));
}



template <int dim, int spacedim>
void
FEValuesBase<dim, spacedim>::TriaCellIterator::
get_interpolated_dof_values (const TrilinosWrappers::MPI::Vector &,
                             Vector<TrilinosScalar>              &) const
{
  Assert (false, ExcMessage (message_string));
}



template <int dim, int spacedim>
void
FEValuesBase<dim,spacedim>::TriaCellIterator::
get_interpolated_dof_values (const TrilinosWrappers::MPI::BlockVector &,
                             Vector<TrilinosScalar>                   &) const
{
  Assert (false, ExcMessage (message_string));
}

#endif




/* --------------------- FEValuesData ----------------- */


template <int dim, int spacedim>
void
FEValuesData<dim,spacedim>::initialize (const unsigned int        n_quadrature_points,
			       const FiniteElement<dim,spacedim> &fe,
			       const UpdateFlags         flags)
{
  this->update_flags = flags;

				   // initialize the table mapping
				   // from shape function number to
				   // the rows in the tables denoting
				   // its first non-zero
				   // component
  this->shape_function_to_row_table
    = make_shape_function_to_row_table (fe);

				   // count the total number of non-zero
				   // components accumulated over all shape
				   // functions
  unsigned int n_nonzero_shape_components = 0;
  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    n_nonzero_shape_components += fe.n_nonzero_components (i);	 
  Assert (n_nonzero_shape_components >= fe.dofs_per_cell,
	  ExcInternalError());

				   // with the number of rows now
				   // known, initialize those fields
				   // that we will need to their
				   // correct size
  if (flags & update_values)
    this->shape_values.reinit(n_nonzero_shape_components,
			      n_quadrature_points);

  if (flags & update_gradients)
    this->shape_gradients.resize (n_nonzero_shape_components,
                                  std::vector<Tensor<1,spacedim> > (n_quadrature_points));

  if (flags & update_hessians)
    this->shape_hessians.resize (n_nonzero_shape_components,
                                        std::vector<Tensor<2,spacedim> > (n_quadrature_points));
  
  if (flags & update_quadrature_points)
    this->quadrature_points.resize(n_quadrature_points);

  if (flags & update_JxW_values)
    this->JxW_values.resize(n_quadrature_points);

  if (flags & update_jacobians)
    this->jacobians.resize(n_quadrature_points);

  if (flags & update_jacobian_grads)
    this->jacobian_grads.resize(n_quadrature_points);

  if (flags & update_inverse_jacobians)
    this->inverse_jacobians.resize(n_quadrature_points);

  if (flags & update_boundary_forms)
    this->boundary_forms.resize(n_quadrature_points);

  if (flags & update_normal_vectors)
    this->normal_vectors.resize(n_quadrature_points);

  if (flags & update_cell_JxW_values)
    this->cell_JxW_values.resize(n_quadrature_points);

  if (flags & update_cell_normal_vectors)
    this->cell_normal_vectors.resize(n_quadrature_points);
}


/*------------------------------- FEValuesBase ---------------------------*/


template <int dim, int spacedim>
FEValuesBase<dim,spacedim>::FEValuesBase (const unsigned int n_q_points,
				 const unsigned int dofs_per_cell,
				 const UpdateFlags flags,
				 const Mapping<dim,spacedim>       &mapping,
				 const FiniteElement<dim,spacedim> &fe)
		:
		n_quadrature_points (n_q_points),
		dofs_per_cell (dofs_per_cell),
		mapping(&mapping),
		fe(&fe),
		mapping_data(0),
		fe_data(0),
		fe_values_views_cache (*this)
{
  this->update_flags = flags;
}



template <int dim, int spacedim>
FEValuesBase<dim,spacedim>::~FEValuesBase ()
{
				   // delete those fields that were
				   // created by the mapping and
				   // finite element objects,
				   // respectively, but of which we
				   // have assumed ownership
  if (fe_data != 0)
    {
      typename Mapping<dim,spacedim>::InternalDataBase *tmp1=fe_data;
      fe_data=0;
      delete tmp1;
    }

  if (mapping_data != 0)
    {
      typename Mapping<dim,spacedim>::InternalDataBase *tmp1=mapping_data;
      mapping_data=0;
      delete tmp1;
    }
}



template <int dim, int spacedim>
template <class InputVector, typename number>
void FEValuesBase<dim,spacedim>::get_function_values (
  const InputVector            &fe_function,
  std::vector<number> &values) const
{
  Assert (this->update_flags & update_values, ExcAccessToUninitializedField());
  Assert (fe->n_components() == 1,
	  ExcDimensionMismatch(fe->n_components(), 1));
  Assert (values.size() == n_quadrature_points,
	  ExcDimensionMismatch(values.size(), n_quadrature_points));
  Assert (present_cell.get() != 0,
	  ExcMessage ("FEValues object is not reinit'ed to any cell"));
  Assert (fe_function.size() == present_cell->n_dofs_for_dof_handler(),
	  ExcDimensionMismatch(fe_function.size(),
			       present_cell->n_dofs_for_dof_handler()));

				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  present_cell->get_interpolated_dof_values(fe_function, dof_values);

				   // initialize with zero
  std::fill_n (values.begin(), n_quadrature_points, 0);

				   // add up contributions of trial
				   // functions. note that here we
				   // deal with scalar finite
				   // elements, so no need to check
				   // for non-primitivity of shape
				   // functions
  for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
    {
      const double value = dof_values(shape_func);
      const double *shape_value_ptr = &this->shape_values(shape_func, 0);
      for (unsigned int point=0; point<n_quadrature_points; ++point)
	values[point] += value * *shape_value_ptr++;
    }
}



template <int dim, int spacedim>
template <class InputVector, typename number>
void FEValuesBase<dim,spacedim>::get_function_values (
  const InputVector& fe_function,
  const VectorSlice<const std::vector<unsigned int> >& indices,
  std::vector<number> &values) const
{
  Assert (this->update_flags & update_values, ExcAccessToUninitializedField());
				   // This function fills a single
				   // component only
  Assert (fe->n_components() == 1,
	  ExcDimensionMismatch(fe->n_components(), 1));
				   // One index for each dof
  Assert (indices.size() == dofs_per_cell,
	  ExcDimensionMismatch(indices.size(), dofs_per_cell));
				   // This vector has one entry for
				   // each quadrature point
  Assert (values.size() == n_quadrature_points,
	  ExcDimensionMismatch(values.size(), n_quadrature_points));
  
				   // initialize with zero
  std::fill_n (values.begin(), n_quadrature_points, 0);
  
				   // add up contributions of trial
				   // functions. note that here we
				   // deal with scalar finite
				   // elements, so no need to check
				   // for non-primitivity of shape
				   // functions
  for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
    {
      const double value = fe_function(indices[shape_func]);
      const double *shape_value_ptr = &this->shape_values(shape_func, 0);
      for (unsigned int point=0; point<n_quadrature_points; ++point)
	values[point] += value * *shape_value_ptr++;
    }
}



template <int dim, int spacedim>
template <class InputVector, typename number>
void FEValuesBase<dim,spacedim>::get_function_values (
  const InputVector&            fe_function,
  std::vector<Vector<number> >& values) const
{
//TODO: Find out how to do this assertion.  
				   // This vector must correspond to a
				   // complete discretization
//  Assert (fe_function.size() == present_cell->get_dof_handler().n_dofs(),
//	  ExcDimensionMismatch(fe_function.size(),
//			     present_cell->get_dof_handler().n_dofs()));
				   // One entry per quadrature point
  Assert (present_cell.get() != 0,
	  ExcMessage ("FEValues object is not reinit'ed to any cell"));
  Assert (values.size() == n_quadrature_points,
	  ExcDimensionMismatch(values.size(), n_quadrature_points));

  const unsigned int n_components = fe->n_components();
				   // Assert that we can write all
				   // components into the result
				   // vectors
  for (unsigned i=0;i<values.size();++i)
    Assert (values[i].size() == n_components,
	    ExcDimensionMismatch(values[i].size(), n_components));
  
  Assert (this->update_flags & update_values, ExcAccessToUninitializedField());
  Assert (fe_function.size() == present_cell->n_dofs_for_dof_handler(),
	  ExcDimensionMismatch(fe_function.size(), present_cell->n_dofs_for_dof_handler()));
    
				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  present_cell->get_interpolated_dof_values(fe_function, dof_values);
  
				   // initialize with zero
  for (unsigned i=0;i<values.size();++i)
    std::fill_n (values[i].begin(), values[i].size(), 0);

				   // add up contributions of trial
				   // functions. now check whether the
				   // shape function is primitive or
				   // not. if it is, then set its only
				   // non-zero component, otherwise
				   // loop over components
  for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
    {
      const double value = dof_values(shape_func);

      if (fe->is_primitive(shape_func))
	{
	  const unsigned int 
	    row = fe->is_primitive() ? 
	          shape_func : this->shape_function_to_row_table[shape_func];

	  const double *shape_value_ptr = &this->shape_values(row, 0);
	  const unsigned int comp = fe->system_to_component_index(shape_func).first;
	  for (unsigned int point=0; point<n_quadrature_points; ++point)
	    values[point](comp) += value * *shape_value_ptr++;
	}
      else
	for (unsigned int c=0; c<n_components; ++c)
	  {
	    if (fe->get_nonzero_components(shape_func)[c] == false)
	      continue;

	    const unsigned int
	      row = (this->shape_function_to_row_table[shape_func]
		     +
		     std::count (fe->get_nonzero_components(shape_func).begin(),
				 fe->get_nonzero_components(shape_func).begin()+c,
				 true));

	    const double *shape_value_ptr = &this->shape_values(row, 0);

	    for (unsigned int point=0; point<n_quadrature_points; ++point)
	      values[point](c) += value * *shape_value_ptr++;
	  }
    }
}



template <int dim, int spacedim>
template <class InputVector, typename number>
void FEValuesBase<dim,spacedim>::get_function_values (
  const InputVector& fe_function,
  const VectorSlice<const std::vector<unsigned int> >& indices,
  std::vector<Vector<number> >& values) const
{
				   // One value per quadrature point
  Assert (n_quadrature_points == values.size(),
	  ExcDimensionMismatch(values.size(), n_quadrature_points));
  
  const unsigned int n_components = fe->n_components();
  
				   // Size of indices must be a
				   // multiple of dofs_per_cell such
				   // that an integer number of
				   // function values is generated in
				   // each point.
  Assert (indices.size() % dofs_per_cell == 0,
	  ExcNotMultiple(indices.size(), dofs_per_cell));

				   // The number of components of the
				   // result may be a multiple of the
				   // number of components of the
				   // finite element
  const unsigned int result_components = indices.size() * n_components / dofs_per_cell;
  
  for (unsigned i=0;i<values.size();++i)
    Assert (values[i].size() == result_components,
	    ExcDimensionMismatch(values[i].size(), result_components));

				   // If the result has more
				   // components than the finite
				   // element, we need this number for
				   // loops filling all components
  const unsigned int component_multiple = result_components / n_components;
  
  Assert (this->update_flags & update_values, ExcAccessToUninitializedField());
    
				   // initialize with zero
  for (unsigned i=0;i<values.size();++i)
    std::fill_n (values[i].begin(), values[i].size(), 0);

				   // add up contributions of trial
				   // functions. now check whether the
				   // shape function is primitive or
				   // not. if it is, then set its only
				   // non-zero component, otherwise
				   // loop over components
  for (unsigned int mc = 0; mc < component_multiple; ++mc)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      {
	const double value = fe_function(indices[shape_func+mc*dofs_per_cell]);

	if (fe->is_primitive(shape_func))
	  {
	    const unsigned int 
	      row = fe->is_primitive() ? 
	            shape_func : this->shape_function_to_row_table[shape_func];

	    const double *shape_value_ptr = &this->shape_values(row, 0);
	    const unsigned int comp = fe->system_to_component_index(shape_func).first 
	                               + mc * n_components;
	    for (unsigned int point=0; point<n_quadrature_points; ++point)
	      values[point](comp) += value * *shape_value_ptr++;
	  }
	else
	  for (unsigned int c=0; c<n_components; ++c)
	    {
	      if (fe->get_nonzero_components(shape_func)[c] == false)
		continue;

	      const unsigned int
		row = (this->shape_function_to_row_table[shape_func]
		       +
		       std::count (fe->get_nonzero_components(shape_func).begin(),
				   fe->get_nonzero_components(shape_func).begin()+c,
				   true));

	      const double *shape_value_ptr = &this->shape_values(row, 0);
	      const unsigned int comp = c + mc * n_components;
	    
	      for (unsigned int point=0; point<n_quadrature_points; ++point)
		values[point](comp) += value * *shape_value_ptr++;
	    }
      }
}



template <int dim, int spacedim>
template <class InputVector, typename number>
void FEValuesBase<dim,spacedim>::get_function_values (
  const InputVector& fe_function,
  const VectorSlice<const std::vector<unsigned int> >& indices,
  std::vector<std::vector<number> >& values,
  bool quadrature_points_fastest) const
{
  const unsigned int n_components = fe->n_components();
  
				   // Size of indices must be a
				   // multiple of dofs_per_cell such
				   // that an integer number of
				   // function values is generated in
				   // each point.
  Assert (indices.size() % dofs_per_cell == 0,
	  ExcNotMultiple(indices.size(), dofs_per_cell));

				   // The number of components of the
				   // result may be a multiple of the
				   // number of components of the
				   // finite element
  const unsigned int result_components = indices.size() * n_components / dofs_per_cell;

				   // Check if the value argument is
				   // initialized to the correct sizes
  if (quadrature_points_fastest)
    {
      Assert (values.size() == result_components,
	      ExcDimensionMismatch(values.size(), result_components));
      for (unsigned i=0;i<values.size();++i)
	Assert (values[i].size() == n_quadrature_points,
		ExcDimensionMismatch(values[i].size(), n_quadrature_points));
    }
  else
    {
      Assert(values.size() == n_quadrature_points,
	     ExcDimensionMismatch(values.size(), n_quadrature_points));
      for (unsigned i=0;i<values.size();++i)
	Assert (values[i].size() == result_components,
		ExcDimensionMismatch(values[i].size(), result_components));
    }
  
				   // If the result has more
				   // components than the finite
				   // element, we need this number for
				   // loops filling all components
  const unsigned int component_multiple = result_components / n_components;
  
  Assert (this->update_flags & update_values, ExcAccessToUninitializedField());
    
				   // initialize with zero
  for (unsigned i=0;i<values.size();++i)
    std::fill_n (values[i].begin(), values[i].size(), 0);

				   // add up contributions of trial
				   // functions. now check whether the
				   // shape function is primitive or
				   // not. if it is, then set its only
				   // non-zero component, otherwise
				   // loop over components
  for (unsigned int mc = 0; mc < component_multiple; ++mc)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      {
	const double value = fe_function(indices[shape_func+mc*dofs_per_cell]);

	if (fe->is_primitive(shape_func))
	  {
	    const unsigned int 
	      row = fe->is_primitive() ? 
	            shape_func : this->shape_function_to_row_table[shape_func];

	    const double *shape_value_ptr = &this->shape_values(row, 0);
	    const unsigned int comp = fe->system_to_component_index(shape_func).first 
	                               + mc * n_components;
	    if (quadrature_points_fastest)
	      for (unsigned int point=0; point<n_quadrature_points; ++point)
		values[comp][point] += value * *shape_value_ptr++;
	    else
	      for (unsigned int point=0; point<n_quadrature_points; ++point)
		values[point][comp] += value * *shape_value_ptr++;
	  }
	else
	  for (unsigned int c=0; c<n_components; ++c)
	    {
	      if (fe->get_nonzero_components(shape_func)[c] == false)
		continue;

	      const unsigned int
		row = (this->shape_function_to_row_table[shape_func]
		       +
		       std::count (fe->get_nonzero_components(shape_func).begin(),
				   fe->get_nonzero_components(shape_func).begin()+c,
				   true));

	      const double *shape_value_ptr = &this->shape_values(row, 0);
	      const unsigned int comp = c + mc * n_components;
	    
	      if (quadrature_points_fastest)
		for (unsigned int point=0; point<n_quadrature_points; ++point)
		  values[comp][point] += value * *shape_value_ptr++;
	      else
		for (unsigned int point=0; point<n_quadrature_points; ++point)
		  values[point][comp] += value * *shape_value_ptr++;
	    }
      }
}



template <int dim, int spacedim>
template <class InputVector>
void
FEValuesBase<dim,spacedim>::get_function_gradients (
  const InputVector           &fe_function,
  std::vector<Tensor<1,spacedim> > &gradients) const
{
  Assert (this->update_flags & update_gradients, ExcAccessToUninitializedField());

  Assert (fe->n_components() == 1,
	  ExcDimensionMismatch(fe->n_components(), 1));
  Assert (gradients.size() == n_quadrature_points,
	  ExcDimensionMismatch(gradients.size(), n_quadrature_points));
  Assert (present_cell.get() != 0,
	  ExcMessage ("FEValues object is not reinit'ed to any cell"));
  Assert (fe_function.size() == present_cell->n_dofs_for_dof_handler(),
	  ExcDimensionMismatch(fe_function.size(),
			       present_cell->n_dofs_for_dof_handler()));

				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  present_cell->get_interpolated_dof_values(fe_function, dof_values);

				   // initialize with zero
  std::fill_n (gradients.begin(), n_quadrature_points, Tensor<1,spacedim>());

				   // add up contributions of trial
				   // functions. note that here we
				   // deal with scalar finite
				   // elements, so no need to check
				   // for non-primitivity of shape
				   // functions
  for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
    {
      const double value = dof_values(shape_func);
      const Tensor<1,spacedim> *shape_gradient_ptr = &this->shape_gradients[shape_func][0];
      for (unsigned int point=0; point<n_quadrature_points; ++point)
	gradients[point] += value * *shape_gradient_ptr++;
    }
}



template <int dim, int spacedim>
template <class InputVector>
void FEValuesBase<dim,spacedim>::get_function_gradients (
  const InputVector& fe_function,
  const VectorSlice<const std::vector<unsigned int> >& indices,
  std::vector<Tensor<1,spacedim> > &gradients) const
{
  Assert (this->update_flags & update_gradients, ExcAccessToUninitializedField());
				   // This function fills a single
				   // component only
  Assert (fe->n_components() == 1,
	  ExcDimensionMismatch(fe->n_components(), 1));
				   // One index for each dof
  Assert (indices.size() == dofs_per_cell,
	  ExcDimensionMismatch(indices.size(), dofs_per_cell));
				   // This vector has one entry for
				   // each quadrature point
  Assert (gradients.size() == n_quadrature_points,
	  ExcDimensionMismatch(gradients.size(), n_quadrature_points));
  
				   // initialize with zero
  std::fill_n (gradients.begin(), n_quadrature_points, Tensor<1,spacedim>());
  
				   // add up contributions of trial
				   // functions. note that here we
				   // deal with scalar finite
				   // elements, so no need to check
				   // for non-primitivity of shape
				   // functions
  for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
    {
      const double value = fe_function(indices[shape_func]);
      const Tensor<1,spacedim> *shape_gradient_ptr = &this->shape_gradients[shape_func][0];
      for (unsigned int point=0; point<n_quadrature_points; ++point)
	gradients[point] += value * *shape_gradient_ptr++;
    }
}




template <int dim, int spacedim>
template <class InputVector>
void
FEValuesBase<dim,spacedim>::get_function_gradients (
  const InputVector                         &fe_function,
  std::vector<std::vector<Tensor<1,spacedim> > > &gradients) const
{
  Assert (gradients.size() == n_quadrature_points,
	  ExcDimensionMismatch(gradients.size(), n_quadrature_points));

  const unsigned int n_components = fe->n_components();
  for (unsigned i=0; i<gradients.size(); ++i)
    Assert (gradients[i].size() == n_components,
	    ExcDimensionMismatch(gradients[i].size(), n_components));

  Assert (this->update_flags & update_gradients, ExcAccessToUninitializedField());
  Assert (present_cell.get() != 0,
	  ExcMessage ("FEValues object is not reinit'ed to any cell"));
  Assert (fe_function.size() == present_cell->n_dofs_for_dof_handler(),
	  ExcDimensionMismatch(fe_function.size(),
			       present_cell->n_dofs_for_dof_handler()));

				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  present_cell->get_interpolated_dof_values(fe_function, dof_values);

				   // initialize with zero
  for (unsigned i=0;i<gradients.size();++i)
    std::fill_n (gradients[i].begin(), gradients[i].size(), Tensor<1,spacedim>());

				   // add up contributions of trial
				   // functions. now check whether the
				   // shape function is primitive or
				   // not. if it is, then set its only
				   // non-zero component, otherwise
				   // loop over components
  for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
    {
      const double value = dof_values(shape_func);

      if (fe->is_primitive(shape_func))
	{
	  const unsigned int 
	    row = fe->is_primitive() ? 
	          shape_func : this->shape_function_to_row_table[shape_func];
	  const Tensor<1,spacedim> *shape_gradient_ptr 
	    = &this->shape_gradients[row][0];
	  const unsigned int comp = fe->system_to_component_index(shape_func).first;
	  for (unsigned int point=0; point<n_quadrature_points; ++point)
	    gradients[point][comp] += value * *shape_gradient_ptr++;
	}
      else
	for (unsigned int c=0; c<n_components; ++c)
	  {
	    if (fe->get_nonzero_components(shape_func)[c] == false)
	      continue;

	    const unsigned int
	      row = (this->shape_function_to_row_table[shape_func]
		     +
		     std::count (fe->get_nonzero_components(shape_func).begin(),
				 fe->get_nonzero_components(shape_func).begin()+c,
				 true));

	    const Tensor<1,spacedim> *shape_gradient_ptr 
	      = &this->shape_gradients[row][0];

	    for (unsigned int point=0; point<n_quadrature_points; ++point)
	      gradients[point][c] += value * *shape_gradient_ptr++;
	  }
    }
}



template <int dim, int spacedim>
template <class InputVector>
void FEValuesBase<dim,spacedim>::get_function_gradients (
  const InputVector& fe_function,
  const VectorSlice<const std::vector<unsigned int> >& indices,
  std::vector<std::vector<Tensor<1,spacedim> > >& gradients,
  bool quadrature_points_fastest) const
{
  const unsigned int n_components = fe->n_components();
  
				   // Size of indices must be a
				   // multiple of dofs_per_cell such
				   // that an integer number of
				   // function values is generated in
				   // each point.
  Assert (indices.size() % dofs_per_cell == 0,
	  ExcNotMultiple(indices.size(), dofs_per_cell));

				   // The number of components of the
				   // result may be a multiple of the
				   // number of components of the
				   // finite element
  const unsigned int result_components = indices.size() * n_components / dofs_per_cell;
  
				   // Check if the value argument is
				   // initialized to the correct sizes
  if (quadrature_points_fastest)
    {
      Assert (gradients.size() == result_components,
	      ExcDimensionMismatch(gradients.size(), result_components));
      for (unsigned i=0;i<gradients.size();++i)
	Assert (gradients[i].size() == n_quadrature_points,
		ExcDimensionMismatch(gradients[i].size(), n_quadrature_points));
    }
  else
    {
      Assert(gradients.size() == n_quadrature_points,
	     ExcDimensionMismatch(gradients.size(), n_quadrature_points));
      for (unsigned i=0;i<gradients.size();++i)
	Assert (gradients[i].size() == result_components,
		ExcDimensionMismatch(gradients[i].size(), result_components));
    }

				   // If the result has more
				   // components than the finite
				   // element, we need this number for
				   // loops filling all components
  const unsigned int component_multiple = result_components / n_components;
  
  Assert (this->update_flags & update_values, ExcAccessToUninitializedField());
    
				   // initialize with zero
  for (unsigned i=0;i<gradients.size();++i)
    std::fill_n (gradients[i].begin(), gradients[i].size(), Tensor<1,spacedim>());

				   // add up contributions of trial
				   // functions. now check whether the
				   // shape function is primitive or
				   // not. if it is, then set its only
				   // non-zero component, otherwise
				   // loop over components
  for (unsigned int mc = 0; mc < component_multiple; ++mc)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      {
	const double value = fe_function(indices[shape_func+mc*dofs_per_cell]);

	if (fe->is_primitive(shape_func))
	  {
	    const unsigned int 
	      row = fe->is_primitive() ? 
	            shape_func : this->shape_function_to_row_table[shape_func];
	    const Tensor<1,spacedim> *shape_gradient_ptr 
	      = &this->shape_gradients[row][0];
	    const unsigned int comp = fe->system_to_component_index(shape_func).first
	                               + mc * n_components;

	  if (quadrature_points_fastest)
	    for (unsigned int point=0; point<n_quadrature_points; ++point)
	      gradients[comp][point] += value * *shape_gradient_ptr++;
	  else
	    for (unsigned int point=0; point<n_quadrature_points; ++point)
	      gradients[point][comp] += value * *shape_gradient_ptr++;
	  }
	else
	  for (unsigned int c=0; c<n_components; ++c)
	    {
	      if (fe->get_nonzero_components(shape_func)[c] == false)
		continue;

	      const unsigned int
		row = (this->shape_function_to_row_table[shape_func]
		       +
		       std::count (fe->get_nonzero_components(shape_func).begin(),
				   fe->get_nonzero_components(shape_func).begin()+c,
				   true));

	      const Tensor<1,spacedim> *shape_gradient_ptr 
		= &this->shape_gradients[row][0];
	      const unsigned int comp = c + mc * n_components;

	      if (quadrature_points_fastest)
		for (unsigned int point=0; point<n_quadrature_points; ++point)
		  gradients[comp][point] += value * *shape_gradient_ptr++;
	      else
		for (unsigned int point=0; point<n_quadrature_points; ++point)
		  gradients[point][comp] += value * *shape_gradient_ptr++;
	    }
      }
}



template <int dim, int spacedim>
template <class InputVector>
void
FEValuesBase<dim,spacedim>::
get_function_hessians (const InputVector           &fe_function,
		       std::vector<Tensor<2,spacedim> > &hessians) const
{
  Assert (fe->n_components() == 1,
	  ExcDimensionMismatch(fe->n_components(), 1));
  Assert (hessians.size() == n_quadrature_points,
	  ExcDimensionMismatch(hessians.size(), n_quadrature_points));
  Assert (this->update_flags & update_hessians, ExcAccessToUninitializedField());
  Assert (present_cell.get() != 0,
	  ExcMessage ("FEValues object is not reinit'ed to any cell"));
  Assert (fe_function.size() == present_cell->n_dofs_for_dof_handler(),
	  ExcDimensionMismatch(fe_function.size(),
			       present_cell->n_dofs_for_dof_handler()));

				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  present_cell->get_interpolated_dof_values(fe_function, dof_values);

				   // initialize with zero
  std::fill_n (hessians.begin(), n_quadrature_points, Tensor<2,spacedim>());

				   // add up contributions of trial
				   // functions. note that here we
				   // deal with scalar finite
				   // elements, so no need to check
				   // for non-primitivity of shape
				   // functions
  for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
    {
      const double value = dof_values(shape_func);
      const Tensor<2,spacedim> *shape_hessians_ptr = &this->shape_hessians[shape_func][0];
      for (unsigned int point=0; point<n_quadrature_points; ++point)
	hessians[point] += value * *shape_hessians_ptr++;
    }
}



template <int dim, int spacedim>
template <class InputVector>
void FEValuesBase<dim,spacedim>::get_function_hessians (
  const InputVector& fe_function,
  const VectorSlice<const std::vector<unsigned int> >& indices,
  std::vector<Tensor<2,spacedim> > &hessians) const
{
  Assert (this->update_flags & update_second_derivatives, ExcAccessToUninitializedField());
				   // This function fills a single
				   // component only
  Assert (fe->n_components() == 1,
	  ExcDimensionMismatch(fe->n_components(), 1));
				   // One index for each dof
  Assert (indices.size() == dofs_per_cell,
	  ExcDimensionMismatch(indices.size(), dofs_per_cell));
				   // This vector has one entry for
				   // each quadrature point
  Assert (hessians.size() == n_quadrature_points,
	  ExcDimensionMismatch(hessians.size(), n_quadrature_points));
  
				   // initialize with zero
  std::fill_n (hessians.begin(), n_quadrature_points, Tensor<2,spacedim>());
  
				   // add up contributions of trial
				   // functions. note that here we
				   // deal with scalar finite
				   // elements, so no need to check
				   // for non-primitivity of shape
				   // functions
  for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
    {
      const double value = fe_function(indices[shape_func]);
      const Tensor<2,spacedim> *shape_hessians_ptr = &this->shape_hessians[shape_func][0];
      for (unsigned int point=0; point<n_quadrature_points; ++point)
	hessians[point] += value * *shape_hessians_ptr++;
    }
}




template <int dim, int spacedim>
template <class InputVector>
void
FEValuesBase<dim,spacedim>::
get_function_hessians (const InputVector                         &fe_function,
		       std::vector<std::vector<Tensor<2,spacedim> > > &hessians,
		       bool quadrature_points_fastest) const
{
  Assert (n_quadrature_points == hessians.size(),
	  ExcDimensionMismatch(hessians.size(), n_quadrature_points));

  const unsigned int n_components = fe->n_components();
  for (unsigned i=0;i<hessians.size();++i)
    Assert (hessians[i].size() == n_components,
	    ExcDimensionMismatch(hessians[i].size(), n_components));

  Assert (this->update_flags & update_hessians, ExcAccessToUninitializedField());
  Assert (present_cell.get() != 0,
	  ExcMessage ("FEValues object is not reinit'ed to any cell"));
  Assert (fe_function.size() == present_cell->n_dofs_for_dof_handler(),
	  ExcDimensionMismatch(fe_function.size(),
			       present_cell->n_dofs_for_dof_handler()));

				   // get function values of dofs
				   // on this cell
  Vector<typename InputVector::value_type> dof_values (dofs_per_cell);
  present_cell->get_interpolated_dof_values(fe_function, dof_values);

				   // initialize with zero
  for (unsigned i=0;i<hessians.size();++i)
    std::fill_n (hessians[i].begin(), hessians[i].size(), Tensor<2,spacedim>());

				   // add up contributions of trial
				   // functions
  for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
    {
      const double value = dof_values(shape_func);

      if (fe->is_primitive(shape_func))
	{
	  const unsigned int 
	    row = fe->is_primitive() ? 
	          shape_func : this->shape_function_to_row_table[shape_func];

	  const Tensor<2,spacedim> *shape_hessian_ptr 
	    = &this->shape_hessians[row][0];
	  const unsigned int comp = fe->system_to_component_index(shape_func).first;

	  if (quadrature_points_fastest)
	    for (unsigned int point=0; point<n_quadrature_points; ++point)
	      hessians[comp][point] += value * *shape_hessian_ptr++;
	  else
	    for (unsigned int point=0; point<n_quadrature_points; ++point)
	      hessians[point][comp] += value * *shape_hessian_ptr++;
	}
      else
	for (unsigned int c=0; c<n_components; ++c)
	  {
	    if (fe->get_nonzero_components(shape_func)[c] == false)
	      continue;

	    const unsigned int
	      row = (this->shape_function_to_row_table[shape_func]
		     +
		     std::count (fe->get_nonzero_components(shape_func).begin(),
				 fe->get_nonzero_components(shape_func).begin()+c,
				 true));

	    const Tensor<2,spacedim> *shape_hessian_ptr 
	      = &this->shape_hessians[row][0];

	    if (quadrature_points_fastest)
	      for (unsigned int point=0; point<n_quadrature_points; ++point)
		hessians[c][point] += value * *shape_hessian_ptr++;
	    else
	      for (unsigned int point=0; point<n_quadrature_points; ++point)
		hessians[point][c] += value * *shape_hessian_ptr++;
	  }
    }
}



template <int dim, int spacedim>
template <class InputVector>
void FEValuesBase<dim, spacedim>::get_function_hessians (
  const InputVector& fe_function,
  const VectorSlice<const std::vector<unsigned int> >& indices,
  std::vector<std::vector<Tensor<2,spacedim> > >& hessians,
  bool quadrature_points_fastest) const
{
  Assert (this->update_flags & update_second_derivatives, ExcAccessToUninitializedField());
    
  const unsigned int n_components = fe->n_components();
  
				   // Size of indices must be a
				   // multiple of dofs_per_cell such
				   // that an integer number of
				   // function values is generated in
				   // each point.
  Assert (indices.size() % dofs_per_cell == 0,
	  ExcNotMultiple(indices.size(), dofs_per_cell));

				   // The number of components of the
				   // result may be a multiple of the
				   // number of components of the
				   // finite element
  const unsigned int result_components = indices.size() * n_components / dofs_per_cell;
  
				   // Check if the value argument is
				   // initialized to the correct sizes
  if (quadrature_points_fastest)
    {
      Assert (hessians.size() == result_components,
	      ExcDimensionMismatch(hessians.size(), result_components));
      for (unsigned i=0;i<hessians.size();++i)
	Assert (hessians[i].size() == n_quadrature_points,
		ExcDimensionMismatch(hessians[i].size(), n_quadrature_points));
    }
  else
    {
      Assert(hessians.size() == n_quadrature_points,
	     ExcDimensionMismatch(hessians.size(), n_quadrature_points));
      for (unsigned i=0;i<hessians.size();++i)
	Assert (hessians[i].size() == result_components,
		ExcDimensionMismatch(hessians[i].size(), result_components));
    }

				   // If the result has more
				   // components than the finite
				   // element, we need this number for
				   // loops filling all components
  const unsigned int component_multiple = result_components / n_components;
  
				   // initialize with zero
  for (unsigned i=0;i<hessians.size();++i)
    std::fill_n (hessians[i].begin(), hessians[i].size(), Tensor<2,spacedim>());

				   // add up contributions of trial
				   // functions. now check whether the
				   // shape function is primitive or
				   // not. if it is, then set its only
				   // non-zero component, otherwise
				   // loop over components
  for (unsigned int mc = 0; mc < component_multiple; ++mc)
    for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
      {
	const double value = fe_function(indices[shape_func+mc*dofs_per_cell]);

	if (fe->is_primitive(shape_func))
	  {
	    const unsigned int 
	      row = fe->is_primitive() ? 
	            shape_func : this->shape_function_to_row_table[shape_func];

	    const Tensor<2,spacedim> *shape_hessian_ptr 
	      = &this->shape_hessians[row][0];
	    const unsigned int comp = fe->system_to_component_index(shape_func).first
	                               + mc * n_components;

	  if (quadrature_points_fastest)
	    for (unsigned int point=0; point<n_quadrature_points; ++point)
	      hessians[comp][point] += value * *shape_hessian_ptr++;
	  else
	    for (unsigned int point=0; point<n_quadrature_points; ++point)
	      hessians[point][comp] += value * *shape_hessian_ptr++;
	  }
	else
	  for (unsigned int c=0; c<n_components; ++c)
	    {
	      if (fe->get_nonzero_components(shape_func)[c] == false)
		continue;

	      const unsigned int
		row = (this->shape_function_to_row_table[shape_func]
		       +
		       std::count (fe->get_nonzero_components(shape_func).begin(),
				   fe->get_nonzero_components(shape_func).begin()+c,
				   true));

	      const Tensor<2,spacedim> *shape_hessian_ptr 
		= &this->shape_hessians[row][0];
	      const unsigned int comp = c + mc * n_components;

	      if (quadrature_points_fastest)
		for (unsigned int point=0; point<n_quadrature_points; ++point)
		  hessians[comp][point] += value * *shape_hessian_ptr++;
	      else
		for (unsigned int point=0; point<n_quadrature_points; ++point)
		  hessians[point][comp] += value * *shape_hessian_ptr++;
	    }
      }
}



template <int dim, int spacedim>
const std::vector<Point<spacedim> > &
FEValuesBase<dim,spacedim>::get_cell_normal_vectors () const
{
  typedef FEValuesBase<dim,spacedim> FEVB;
  Assert (this->update_flags & update_cell_normal_vectors,
	  typename FEVB::ExcAccessToUninitializedField());
  return this->cell_normal_vectors;
}



template <int dim, int spacedim>
unsigned int
FEValuesBase<dim,spacedim>::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (this->shape_values) +
	  MemoryConsumption::memory_consumption (this->shape_gradients) +
	  MemoryConsumption::memory_consumption (this->shape_hessians) +
	  MemoryConsumption::memory_consumption (this->JxW_values) +
	  MemoryConsumption::memory_consumption (this->jacobians) +
	  MemoryConsumption::memory_consumption (this->jacobian_grads) +
	  MemoryConsumption::memory_consumption (this->inverse_jacobians) +
	  MemoryConsumption::memory_consumption (this->quadrature_points) +
	  MemoryConsumption::memory_consumption (this->normal_vectors) +
	  MemoryConsumption::memory_consumption (this->boundary_forms) +
	  MemoryConsumption::memory_consumption (this->cell_JxW_values) +
	  sizeof(this->update_flags) +
	  MemoryConsumption::memory_consumption (n_quadrature_points) +
	  MemoryConsumption::memory_consumption (dofs_per_cell) +
	  MemoryConsumption::memory_consumption (mapping) +
	  MemoryConsumption::memory_consumption (fe) +
	  MemoryConsumption::memory_consumption (mapping_data) +
	  MemoryConsumption::memory_consumption (*mapping_data) +
	  MemoryConsumption::memory_consumption (fe_data) +
	  MemoryConsumption::memory_consumption (*fe_data) +
	  MemoryConsumption::memory_consumption (this->shape_function_to_row_table));
}



template <int dim, int spacedim>
UpdateFlags
FEValuesBase<dim,spacedim>::compute_update_flags (const UpdateFlags update_flags) const
{
  
				   // first find out which objects
				   // need to be recomputed on each
				   // cell we visit. this we have to
				   // ask the finite element and mapping.
				   // elements are first since they
				   // might require update in mapping
  UpdateFlags flags = update_flags
		      | fe->update_once (update_flags)
		      | fe->update_each (update_flags);
  flags |= mapping->update_once (flags)
	   | mapping->update_each (flags);

  return flags;
}



/*------------------------------- FEValues -------------------------------*/


template <int dim, int spacedim>
const unsigned int FEValues<dim,spacedim>::dimension;

template <int dim, int spacedim>
const unsigned int FEValues<dim,spacedim>::integral_dimension;



template <int dim, int spacedim>
FEValues<dim,spacedim>::FEValues (const Mapping<dim,spacedim>       &mapping,
			 const FiniteElement<dim,spacedim> &fe,
			 const Quadrature<dim>    &q,
			 const UpdateFlags         update_flags)
		:
		FEValuesBase<dim,spacedim> (q.size(),
					    fe.dofs_per_cell,
					    update_default,
					    mapping,
					    fe),
		quadrature (q)
{
  initialize (update_flags);
}



template <int dim, int spacedim>
FEValues<dim,spacedim>::FEValues (const FiniteElement<dim,spacedim> &fe,
			 const Quadrature<dim>    &q,
			 const UpdateFlags         update_flags)
		:
		FEValuesBase<dim,spacedim> (q.size(),
				   fe.dofs_per_cell,
				   update_default,
					    StaticMappingQ1<dim,spacedim>::mapping,
				   fe),
  quadrature (q)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  initialize (update_flags);
}



template <int dim, int spacedim>
void
FEValues<dim,spacedim>::initialize (const UpdateFlags update_flags)
{
				   // you can't compute normal vectors
				   // on cells, only on faces
  typedef FEValuesBase<dim,spacedim> FEVB;
  Assert ((update_flags & update_normal_vectors) == false,
	  typename FEVB::ExcInvalidUpdateFlag());

				   // You can compute normal vectors
				   // to the cells only in the
				   // codimension one case.
  if(dim == spacedim)
    Assert ( (update_flags & update_cell_normal_vectors) == false,
	     typename FEVB::ExcInvalidUpdateFlag());


  const UpdateFlags flags = this->compute_update_flags (update_flags);
  
				   // then get objects into which the
				   // FE and the Mapping can store
				   // intermediate data used across
				   // calls to reinit
  this->mapping_data = this->mapping->get_data(flags, quadrature);
  this->fe_data      = this->fe->get_data(flags, *this->mapping, quadrature);

				   // set up objects within this class
  FEValuesData<dim,spacedim>::initialize (this->n_quadrature_points, *this->fe, flags);
}



template <int dim, int spacedim>
void
FEValues<dim,spacedim>::reinit (const typename DoFHandler<dim,spacedim>::cell_iterator &cell)
{
				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the DoFHandler used by
				   // this cell, are the same

  typedef FEValuesBase<dim,spacedim> FEVB;
  Assert (static_cast<const FiniteElementData<dim>&>(*this->fe) ==
	  static_cast<const FiniteElementData<dim>&>(cell->get_fe()),
	  typename FEVB::ExcFEDontMatch());

                                   // set new cell. auto_ptr will take
                                   // care that old object gets
                                   // destroyed and also that this
                                   // object gets destroyed in the
                                   // destruction of this class
  this->present_cell.reset 
    (new typename FEValuesBase<dim,spacedim>::template
     CellIterator<typename DoFHandler<dim,spacedim>::cell_iterator> (cell));

                                   // this was the part of the work
                                   // that is dependent on the actual
                                   // data type of the iterator. now
                                   // pass on to the function doing
                                   // the real work.
  do_reinit ();
}



template <int dim, int spacedim>
void
FEValues<dim,spacedim>::reinit (const typename hp::DoFHandler<dim,spacedim>::cell_iterator &cell)
{
				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the DoFHandler used by
				   // this cell, are the same
    typedef FEValuesBase<dim,spacedim> FEVB;
  Assert (static_cast<const FiniteElementData<dim>&>(*this->fe) ==
	  static_cast<const FiniteElementData<dim>&>(cell->get_fe()),
	  typename FEVB::ExcFEDontMatch());

                                   // set new cell. auto_ptr will take
                                   // care that old object gets
                                   // destroyed and also that this
                                   // object gets destroyed in the
                                   // destruction of this class
  this->present_cell.reset 
    (new typename FEValuesBase<dim,spacedim>::template
     CellIterator<typename hp::DoFHandler<dim,spacedim>::cell_iterator> (cell));

                                   // this was the part of the work
                                   // that is dependent on the actual
                                   // data type of the iterator. now
                                   // pass on to the function doing
                                   // the real work.
  do_reinit ();
}


#if deal_II_dimension == 1

template <>
void
FEValues<1,2>::reinit (const MGDoFHandler<1,2>::cell_iterator &)
{
  Assert(false, ExcNotImplemented());
}

#endif

#if deal_II_dimension == 2

template <>
void
FEValues<2,3>::reinit (const MGDoFHandler<2,3>::cell_iterator &)
{
  Assert(false, ExcNotImplemented());
}

#endif

template <int dim, int spacedim>
void
FEValues<dim,spacedim>::reinit (const typename MGDoFHandler<dim,spacedim>::cell_iterator &cell)
{
				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the DoFHandler used by
				   // this cell, are the same

//TODO: This was documented out ith the repository. Why?
  
//  Assert (static_cast<const FiniteElementData<dim>&>(*this->fe) ==
//	  static_cast<const FiniteElementData<dim>&>(cell->get_fe()),
//	  typename FEValuesBase<dim,spacedim>::ExcFEDontMatch());

                                   // set new cell. auto_ptr will take
                                   // care that old object gets
                                   // destroyed and also that this
                                   // object gets destroyed in the
                                   // destruction of this class
  this->present_cell.reset 
    (new typename FEValuesBase<dim,spacedim>::template
     CellIterator<typename MGDoFHandler<dim,spacedim>::cell_iterator> (cell));

                                   // this was the part of the work
                                   // that is dependent on the actual
                                   // data type of the iterator. now
                                   // pass on to the function doing
                                   // the real work.
  do_reinit ();
}



template <int dim, int spacedim>
void FEValues<dim,spacedim>::reinit (const typename Triangulation<dim,spacedim>::cell_iterator &cell)
{
                                   // no FE in this cell, so no check
                                   // necessary here

                                   // set new cell. auto_ptr will take
                                   // care that old object gets
                                   // destroyed and also that this
                                   // object gets destroyed in the
                                   // destruction of this class
  this->present_cell.reset 
    (new typename FEValuesBase<dim,spacedim>::TriaCellIterator (cell));
                                   // this was the part of the work
                                   // that is dependent on the actual
                                   // data type of the iterator. now
                                   // pass on to the function doing
                                   // the real work.
  do_reinit ();
}



template <int dim, int spacedim>
void FEValues<dim,spacedim>::do_reinit ()
{
  this->get_mapping().fill_fe_values(*this->present_cell,
				     quadrature,
				     *this->mapping_data,
				     this->quadrature_points,
				     this->JxW_values,
				     this->jacobians,
				     this->jacobian_grads,
				     this->inverse_jacobians,
				     this->cell_normal_vectors);
  
  this->get_fe().fill_fe_values(this->get_mapping(),
				*this->present_cell,
				quadrature,
				*this->mapping_data,
				*this->fe_data,
				*this);

  this->fe_data->clear_first_cell ();
  this->mapping_data->clear_first_cell ();
}



template <int dim, int spacedim>
unsigned int
FEValues<dim,spacedim>::memory_consumption () const
{
  return (FEValuesBase<dim,spacedim>::memory_consumption () +
	  MemoryConsumption::memory_consumption (quadrature));
}


/*------------------------------- FEFaceValuesBase --------------------------*/


template <int dim, int spacedim>
FEFaceValuesBase<dim,spacedim>::FEFaceValuesBase (const unsigned int n_q_points,
					 const unsigned int dofs_per_cell,
					 const UpdateFlags,
					 const Mapping<dim,spacedim> &mapping,      
					 const FiniteElement<dim,spacedim> &fe,
					 const Quadrature<dim-1>& quadrature)
		:
		FEValuesBase<dim,spacedim> (n_q_points,
				   dofs_per_cell,
				   update_default,
				   mapping,
				   fe),
                quadrature(quadrature)
{}



template <int dim, int spacedim>
const std::vector<Point<dim> > &
FEFaceValuesBase<dim,spacedim>::get_normal_vectors () const
{
  typedef FEValuesBase<dim,spacedim> FEVB;
  Assert (this->update_flags & update_normal_vectors,
	  typename FEVB::ExcAccessToUninitializedField());
  return this->normal_vectors;
}



template <int dim, int spacedim>
const std::vector<Tensor<1,spacedim> > &
FEFaceValuesBase<dim,spacedim>::get_boundary_forms () const
{
  typedef FEValuesBase<dim,spacedim> FEVB;
  Assert (this->update_flags & update_boundary_forms,
	  typename FEVB::ExcAccessToUninitializedField());
  return this->boundary_forms;
}



template <int dim, int spacedim>
unsigned int
FEFaceValuesBase<dim,spacedim>::memory_consumption () const
{
  return (FEValuesBase<dim,spacedim>::memory_consumption () +
	  MemoryConsumption::memory_consumption (quadrature));
}


/*------------------------------- FEFaceValues -------------------------------*/

template <int dim, int spacedim>
const unsigned int FEFaceValues<dim,spacedim>::dimension;

template <int dim, int spacedim>
const unsigned int FEFaceValues<dim,spacedim>::integral_dimension;


template <int dim, int spacedim>
FEFaceValues<dim,spacedim>::FEFaceValues (const Mapping<dim,spacedim>       &mapping,
				 const FiniteElement<dim,spacedim> &fe,
				 const Quadrature<dim-1>  &quadrature,
				 const UpdateFlags         update_flags)
		:
		FEFaceValuesBase<dim,spacedim> (quadrature.size(),
				       fe.dofs_per_cell,
				       update_flags,
				       mapping,
				       fe, quadrature)
{
  initialize (update_flags);
}



template <int dim, int spacedim>
FEFaceValues<dim,spacedim>::FEFaceValues (const FiniteElement<dim,spacedim> &fe,
				 const Quadrature<dim-1>  &quadrature,
				 const UpdateFlags         update_flags)
		:
		FEFaceValuesBase<dim,spacedim> (quadrature.size(),
				       fe.dofs_per_cell,
				       update_flags,
				       StaticMappingQ1<dim>::mapping,
				       fe, quadrature)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  initialize (update_flags);
}



template <int dim, int spacedim>
void
FEFaceValues<dim,spacedim>::initialize (const UpdateFlags update_flags)
{
  const UpdateFlags flags = this->compute_update_flags (update_flags);
  
				   // then get objects into which the
				   // FE and the Mapping can store
				   // intermediate data used across
				   // calls to reinit
  this->mapping_data = this->mapping->get_face_data(flags, this->quadrature);
  this->fe_data      = this->fe->get_face_data(flags, *this->mapping, this->quadrature);

				   // set up objects within this class
  FEValuesData<dim,spacedim>::initialize(this->n_quadrature_points, *this->fe, flags);
}



template <int dim, int spacedim>
void FEFaceValues<dim,spacedim>::reinit (const typename DoFHandler<dim,spacedim>::cell_iterator &cell,
				const unsigned int              face_no)
{
				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the DoFHandler used by
				   // this cell, are the same
//   Assert (static_cast<const FiniteElementData<dim>&>(*this->fe) ==
// 	  static_cast<const FiniteElementData<dim>&>(cell->get_dof_handler().get_fe()),
// 	  typename FEValuesBase<dim,spacedim>::ExcFEDontMatch());

  Assert (face_no < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_no, 0, GeometryInfo<dim>::faces_per_cell));
  
                                   // set new cell. auto_ptr will take
                                   // care that old object gets
                                   // destroyed and also that this
                                   // object gets destroyed in the
                                   // destruction of this class
  this->present_cell.reset 
    (new typename FEValuesBase<dim,spacedim>::template
     CellIterator<typename DoFHandler<dim,spacedim>::cell_iterator> (cell));
  
                                   // this was the part of the work
                                   // that is dependent on the actual
                                   // data type of the iterator. now
                                   // pass on to the function doing
                                   // the real work.
  do_reinit (face_no);
}



template <int dim, int spacedim>
void FEFaceValues<dim,spacedim>::reinit (const typename hp::DoFHandler<dim,spacedim>::cell_iterator &cell,
				const unsigned int              face_no)
{
				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the DoFHandler used by
				   // this cell, are the same
  typedef FEValuesBase<dim,spacedim> FEVB;
  Assert (static_cast<const FiniteElementData<dim>&>(*this->fe) ==
	  static_cast<const FiniteElementData<dim>&>(
	    cell->get_dof_handler().get_fe()[cell->active_fe_index ()]),
	  typename FEVB::ExcFEDontMatch());

  Assert (face_no < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_no, 0, GeometryInfo<dim>::faces_per_cell));
  
                                   // set new cell. auto_ptr will take
                                   // care that old object gets
                                   // destroyed and also that this
                                   // object gets destroyed in the
                                   // destruction of this class
  this->present_cell.reset 
    (new typename FEValuesBase<dim,spacedim>::template
     CellIterator<typename hp::DoFHandler<dim,spacedim>::cell_iterator> (cell));

                                   // this was the part of the work
                                   // that is dependent on the actual
                                   // data type of the iterator. now
                                   // pass on to the function doing
                                   // the real work.
  do_reinit (face_no);
}


#if deal_II_dimension == 1

template <>
void FEFaceValues<1,2>::reinit (const MGDoFHandler<1,2>::cell_iterator &,
				const unsigned int)
{
  Assert(false,ExcNotImplemented());
}

#endif

#if deal_II_dimension == 2

template <>
void FEFaceValues<2,3>::reinit (const MGDoFHandler<2,3>::cell_iterator &,
				const unsigned int)
{
  Assert(false,ExcNotImplemented());
}

#endif

template <int dim, int spacedim>
void FEFaceValues<dim,spacedim>::reinit (const typename MGDoFHandler<dim,spacedim>::cell_iterator &cell,
				const unsigned int              face_no)
{
				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the DoFHandler used by
				   // this cell, are the same
  typedef FEValuesBase<dim,spacedim> FEVB;
  Assert (static_cast<const FiniteElementData<dim>&>(*this->fe) ==
	  static_cast<const FiniteElementData<dim>&>(cell->get_dof_handler().get_fe()),
	  typename FEVB::ExcFEDontMatch());

  Assert (face_no < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_no, 0, GeometryInfo<dim>::faces_per_cell));
  
                                   // set new cell. auto_ptr will take
                                   // care that old object gets
                                   // destroyed and also that this
                                   // object gets destroyed in the
                                   // destruction of this class
  this->present_cell.reset 
    (new typename FEValuesBase<dim,spacedim>::template
     CellIterator<typename MGDoFHandler<dim,spacedim>::cell_iterator> (cell));

                                   // this was the part of the work
                                   // that is dependent on the actual
                                   // data type of the iterator. now
                                   // pass on to the function doing
                                   // the real work.
  do_reinit (face_no);
}


template <int dim, int spacedim>
void FEFaceValues<dim,spacedim>::reinit (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
				const unsigned int              face_no)
{
  Assert (face_no < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_no, 0, GeometryInfo<dim>::faces_per_cell));

                                   // set new cell. auto_ptr will take
                                   // care that old object gets
                                   // destroyed and also that this
                                   // object gets destroyed in the
                                   // destruction of this class
  this->present_cell.reset 
    (new typename FEValuesBase<dim,spacedim>::TriaCellIterator (cell));

                                   // this was the part of the work
                                   // that is dependent on the actual
                                   // data type of the iterator. now
                                   // pass on to the function doing
                                   // the real work.
  do_reinit (face_no);
}



template <int dim, int spacedim>
void FEFaceValues<dim,spacedim>::do_reinit (const unsigned int face_no)
{
				   // first of all, set the
				   // present_face_index (if
				   // available)
  const typename Triangulation<dim,spacedim>::cell_iterator cell=*this->present_cell;
  this->present_face_index=cell->face_index(face_no);
  
  this->get_mapping().fill_fe_face_values(*this->present_cell, face_no,
					  this->quadrature,
					  *this->mapping_data,
					  this->quadrature_points,
					  this->JxW_values,
					  this->boundary_forms,
					  this->normal_vectors,
					  this->cell_JxW_values);
  
  this->get_fe().fill_fe_face_values(this->get_mapping(),
				     *this->present_cell, face_no,
				     this->quadrature,
				     *this->mapping_data,
				     *this->fe_data,
				     *this);

  this->fe_data->clear_first_cell ();
  this->mapping_data->clear_first_cell ();
}


/*------------------------------- FESubFaceValues -------------------------------*/


template <int dim, int spacedim>
const unsigned int FESubfaceValues<dim,spacedim>::dimension;

template <int dim, int spacedim>
const unsigned int FESubfaceValues<dim,spacedim>::integral_dimension;



template <int dim, int spacedim>
FESubfaceValues<dim,spacedim>::FESubfaceValues (const Mapping<dim,spacedim>       &mapping,
				       const FiniteElement<dim,spacedim> &fe,
				       const Quadrature<dim-1>  &quadrature,
				       const UpdateFlags         update_flags)
		:
		FEFaceValuesBase<dim,spacedim> (quadrature.size(),
				       fe.dofs_per_cell,
				       update_flags,
				       mapping,
				       fe, quadrature)
{
  initialize (update_flags);
}



template <int dim, int spacedim>
FESubfaceValues<dim,spacedim>::FESubfaceValues (const FiniteElement<dim,spacedim> &fe,
				       const Quadrature<dim-1>  &quadrature,
				       const UpdateFlags         update_flags)
		:
		FEFaceValuesBase<dim,spacedim> (quadrature.size(),
				       fe.dofs_per_cell,
				       update_flags,
				       StaticMappingQ1<dim>::mapping,
				       fe, quadrature)
{
  Assert (DEAL_II_COMPAT_MAPPING, ExcCompatibility("mapping"));
  initialize (update_flags);
}



template <int dim, int spacedim>
void
FESubfaceValues<dim,spacedim>::initialize (const UpdateFlags update_flags)
{
  const UpdateFlags flags = this->compute_update_flags (update_flags);
  
				   // then get objects into which the
				   // FE and the Mapping can store
				   // intermediate data used across
				   // calls to reinit
  this->mapping_data = this->mapping->get_subface_data(flags, this->quadrature);
  this->fe_data      = this->fe->get_subface_data(flags,
						  *this->mapping,
						  this->quadrature);

				   // set up objects within this class
  FEValuesData<dim,spacedim>::initialize(this->n_quadrature_points, *this->fe, flags);
}



template <int dim, int spacedim>
void FESubfaceValues<dim,spacedim>::reinit (const typename DoFHandler<dim,spacedim>::cell_iterator &cell,
				   const unsigned int         face_no,
				   const unsigned int         subface_no)
{
				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the DoFHandler used by
				   // this cell, are the same
//   Assert (static_cast<const FiniteElementData<dim>&>(*this->fe) ==
// 	  static_cast<const FiniteElementData<dim>&>(cell->get_dof_handler().get_fe()),
// 	  typename FEValuesBase<dim,spacedim>::ExcFEDontMatch());
  Assert (face_no < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_no, 0, GeometryInfo<dim>::faces_per_cell));
				   // We would like to check for
				   // subface_no < cell->face(face_no)->n_children(),
				   // but unfortunately the current
				   // function is also called for
				   // faces without children (see
				   // tests/fe/mapping.cc). Therefore,
				   // we must use following workaround
				   // of two separate assertions
  Assert (cell->face(face_no)->has_children() ||
	  subface_no < GeometryInfo<dim>::max_children_per_face,
	  ExcIndexRange (subface_no, 0, GeometryInfo<dim>::max_children_per_face));
  Assert (!cell->face(face_no)->has_children() ||
	  subface_no < cell->face(face_no)->number_of_children(),
	  ExcIndexRange (subface_no, 0, cell->face(face_no)->number_of_children()));

  Assert (cell->has_children() == false,
          ExcMessage ("You can't use subface data for cells that are "
                      "already refined. Iterate over their children "
                      "instead in these cases."));

                                   // set new cell. auto_ptr will take
                                   // care that old object gets
                                   // destroyed and also that this
                                   // object gets destroyed in the
                                   // destruction of this class
  this->present_cell.reset 
    (new typename FEValuesBase<dim,spacedim>::template
     CellIterator<typename DoFHandler<dim,spacedim>::cell_iterator> (cell));

                                   // this was the part of the work
                                   // that is dependent on the actual
                                   // data type of the iterator. now
                                   // pass on to the function doing
                                   // the real work.
  do_reinit (face_no, subface_no);
}



template <int dim, int spacedim>
void FESubfaceValues<dim,spacedim>::reinit (const typename hp::DoFHandler<dim,spacedim>::cell_iterator &cell,
				   const unsigned int         face_no,
				   const unsigned int         subface_no)
{
				   // assert that the finite elements
				   // passed to the constructor and
				   // used by the hp::DoFHandler used by
				   // this cell, are the same
  typedef FEValuesBase<dim,spacedim> FEVB;
  Assert (static_cast<const FiniteElementData<dim>&>(*this->fe) ==
	  static_cast<const FiniteElementData<dim>&>(
	    cell->get_dof_handler().get_fe()[cell->active_fe_index ()]),
	  typename FEVB::ExcFEDontMatch());
  Assert (face_no < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_no, 0, GeometryInfo<dim>::faces_per_cell));
  Assert (subface_no < cell->face(face_no)->number_of_children(),
	  ExcIndexRange (subface_no, 0, cell->face(face_no)->number_of_children()));
  Assert (cell->has_children() == false,
          ExcMessage ("You can't use subface data for cells that are "
                      "already refined. Iterate over their children "
                      "instead in these cases."));

                                   // set new cell. auto_ptr will take
                                   // care that old object gets
                                   // destroyed and also that this
                                   // object gets destroyed in the
                                   // destruction of this class
  this->present_cell.reset
    (new typename FEValuesBase<dim,spacedim>::template
     CellIterator<typename hp::DoFHandler<dim,spacedim>::cell_iterator> (cell));

                                   // this was the part of the work
                                   // that is dependent on the actual
                                   // data type of the iterator. now
                                   // pass on to the function doing
                                   // the real work.
  do_reinit (face_no, subface_no);
}

  
template <int dim, int spacedim>
void FESubfaceValues<dim,spacedim>::reinit (const typename MGDoFHandler<dim,spacedim>::cell_iterator &cell,
				   const unsigned int         face_no,
				   const unsigned int         subface_no)
{
  typedef FEValuesBase<dim,spacedim> FEVB;
  Assert (static_cast<const FiniteElementData<dim>&>(*this->fe) ==
	   static_cast<const FiniteElementData<dim>&>(cell->get_dof_handler().get_fe()),
	   typename FEVB::ExcFEDontMatch());
  Assert (face_no < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_no, 0, GeometryInfo<dim>::faces_per_cell));
  Assert (subface_no < cell->face(face_no)->number_of_children(),
	  ExcIndexRange (subface_no, 0, cell->face(face_no)->number_of_children()));
  Assert (cell->has_children() == false,
          ExcMessage ("You can't use subface data for cells that are "
                      "already refined. Iterate over their children "
                      "instead in these cases."));

                                   // set new cell. auto_ptr will take
                                   // care that old object gets
                                   // destroyed and also that this
                                   // object gets destroyed in the
                                   // destruction of this class
  this->present_cell.reset 
    (new typename FEValuesBase<dim,spacedim>::template
     CellIterator<typename MGDoFHandler<dim,spacedim>::cell_iterator> (cell));

                                   // this was the part of the work
                                   // that is dependent on the actual
                                   // data type of the iterator. now
                                   // pass on to the function doing
                                   // the real work.
  do_reinit (face_no, subface_no);
}



template <int dim, int spacedim>
void FESubfaceValues<dim,spacedim>::reinit (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
				   const unsigned int         face_no,
				   const unsigned int         subface_no)
{
  Assert (face_no < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_no, 0, GeometryInfo<dim>::faces_per_cell));
  Assert (subface_no < cell->face(face_no)->n_children(),
	  ExcIndexRange (subface_no, 0, cell->face(face_no)->n_children()));
  
                                   // set new cell. auto_ptr will take
                                   // care that old object gets
                                   // destroyed and also that this
                                   // object gets destroyed in the
                                   // destruction of this class
  this->present_cell.reset 
    (new typename FEValuesBase<dim,spacedim>::TriaCellIterator (cell));

                                   // this was the part of the work
                                   // that is dependent on the actual
                                   // data type of the iterator. now
                                   // pass on to the function doing
                                   // the real work.
  do_reinit (face_no, subface_no);
}



template <int dim, int spacedim>
void FESubfaceValues<dim,spacedim>::do_reinit (const unsigned int face_no,
                                               const unsigned int subface_no)
{
				   // first of all, set the present_face_index
				   // (if available)
  const typename Triangulation<dim,spacedim>::cell_iterator cell=*this->present_cell;

  if (!cell->face(face_no)->has_children())
				     // no subfaces at all, so set
				     // present_face_index to this face rather
				     // than any subface
    this->present_face_index=cell->face_index(face_no);
  else
    if (dim!=3)
      this->present_face_index=cell->face(face_no)->child_index(subface_no);
    else
      {
					 // this is the same logic we use in
					 // cell->neighbor_child_on_subface(). See
					 // there for an explanation of the
					 // different cases
	unsigned int subface_index=numbers::invalid_unsigned_int;
	switch (cell->subface_case(face_no))
	  {
	    case internal::SubfaceCase<3>::case_x:
	    case internal::SubfaceCase<3>::case_y:
	    case internal::SubfaceCase<3>::case_xy:
		  subface_index=cell->face(face_no)->child_index(subface_no);
		  break;
	    case internal::SubfaceCase<3>::case_x1y2y:
	    case internal::SubfaceCase<3>::case_y1x2x:
		  subface_index=cell->face(face_no)->child(subface_no/2)->child_index(subface_no%2);
		  break;
	    case internal::SubfaceCase<3>::case_x1y:
	    case internal::SubfaceCase<3>::case_y1x:
		  switch (subface_no)
		    {
		      case 0:
		      case 1:
			    subface_index=cell->face(face_no)->child(0)->child_index(subface_no);
			    break;
		      case 2:
			    subface_index=cell->face(face_no)->child_index(1);
			    break;
		      default:
			    Assert(false, ExcInternalError());
		    }
		  break;
	    case internal::SubfaceCase<3>::case_x2y:
	    case internal::SubfaceCase<3>::case_y2x:
		  switch (subface_no)
		    {
		      case 0:
			    subface_index=cell->face(face_no)->child_index(0);
			    break;
		      case 1:
		      case 2:
			    subface_index=cell->face(face_no)->child(1)->child_index(subface_no-1);
			    break;
		      default:
			    Assert(false, ExcInternalError());
		    }
		  break;
	    default:
		  Assert(false, ExcInternalError());
		  break;
	  }
	Assert(subface_index!=numbers::invalid_unsigned_int,
	       ExcInternalError());
	this->present_face_index=subface_index;
      }
  
				   // now ask the mapping and the finite element
				   // to do the actual work
  this->get_mapping().fill_fe_subface_values(*this->present_cell,
                                             face_no, subface_no,
					     this->quadrature,
					     *this->mapping_data,
					     this->quadrature_points,
					     this->JxW_values,
					     this->boundary_forms,
					     this->normal_vectors,
					     this->cell_JxW_values);
  
  this->get_fe().fill_fe_subface_values(this->get_mapping(),
					*this->present_cell,
                                        face_no, subface_no,
					this->quadrature,
					*this->mapping_data,
					*this->fe_data,
					*this);

  this->fe_data->clear_first_cell ();
  this->mapping_data->clear_first_cell ();
}


/*------------------------------- Explicit Instantiations -------------*/

template class FEValuesData<deal_II_dimension>;
template class FEValuesBase<deal_II_dimension>;
template class FEValues<deal_II_dimension>;
template class FEValuesBase<deal_II_dimension>::
  CellIterator<DoFHandler<deal_II_dimension>::cell_iterator>;
template class FEValuesBase<deal_II_dimension>::
  CellIterator<MGDoFHandler<deal_II_dimension>::cell_iterator>;

#if deal_II_dimension >= 2
template class FEFaceValuesBase<deal_II_dimension>;
template class FEFaceValues<deal_II_dimension>;
template class FESubfaceValues<deal_II_dimension>;
#endif

#if deal_II_dimension != 3
template class FEValuesData<deal_II_dimension,deal_II_dimension+1>;
template class FEValuesBase<deal_II_dimension,deal_II_dimension+1>;
template class FEValues<deal_II_dimension,deal_II_dimension+1>;
template class FEValuesBase<deal_II_dimension,deal_II_dimension+1>:: 
  CellIterator<DoFHandler<deal_II_dimension,deal_II_dimension+1>::cell_iterator>; 
//template class FEValuesBase<deal_II_dimension,deal_II_dimension+1>::
//  CellIterator<MGDoFHandler<deal_II_dimension,deal_II_dimension+1>::cell_iterator>;

// #if deal_II_dimension == 2
// template class FEFaceValuesBase<deal_II_dimension,deal_II_dimension+1>;
// template class FEFaceValues<deal_II_dimension,deal_II_dimension+1>;
// template class FESubfaceValues<deal_II_dimension,deal_II_dimension+1>;
// #endif



#endif

//---------------------------------------------------------------------------
// Instantiations are in a different file using the macro IN for the vector type.
// This way, we avoid code reduplication

#include "fe_values.inst"

DEAL_II_NAMESPACE_CLOSE
