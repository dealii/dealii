//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <fe/fe_nothing.h>

DEAL_II_NAMESPACE_OPEN

namespace
{
  const char*
  zero_dof_message = "This element has no shape functions.";
}




template <int dim>
FE_Nothing<dim>::FE_Nothing ()
                :
                FiniteElement<dim>
		(FiniteElementData<dim>(std::vector<unsigned>(dim+1,0),
					1, 0,
					FiniteElementData<dim>::unknown),
		 std::vector<bool>(),
		 std::vector<std::vector<bool> >() )
{}


template <int dim>
FiniteElement<dim> *
FE_Nothing<dim>::clone() const
{
  return new FE_Nothing<dim>(*this);
}



template <int dim>
std::string
FE_Nothing<dim>::get_name () const
{
  std::ostringstream namebuf;
  namebuf << "FE_Nothing<" << dim << ">()";
  return namebuf.str();
}



template <int dim>
unsigned int
FE_Nothing<dim>::n_base_elements () const
{
  return 1;
}



template <int dim>
const FiniteElement<dim> &
FE_Nothing<dim>::base_element (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return *this;
}



template <int dim>
unsigned int
FE_Nothing<dim>::element_multiplicity (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return 1;
}



template <int dim>
UpdateFlags
FE_Nothing<dim>::update_once (const UpdateFlags /*flags*/) const
{
  return update_default;
}



template <int dim>
UpdateFlags
FE_Nothing<dim>::update_each (const UpdateFlags /*flags*/) const
{
  return update_default;
}



template <int dim>
double
FE_Nothing<dim>::shape_value (const unsigned int /*i*/,
			      const Point<dim> & /*p*/) const
{
  Assert(false,ExcMessage(zero_dof_message));
  return 0;
}



template <int dim>
typename Mapping<dim>::InternalDataBase *
FE_Nothing<dim>::get_data (const UpdateFlags  /*flags*/,
			   const Mapping<dim> & /*mapping*/,
			   const Quadrature<dim> & /*quadrature*/) const
{
  Assert(false,ExcMessage(zero_dof_message));
  return NULL;
}



template <int dim>
void
FE_Nothing<dim>::
fill_fe_values (const Mapping<dim>                      & /*mapping*/,
		const typename Triangulation<dim>::cell_iterator & /*cell*/,
		const Quadrature<dim> & /*quadrature*/,
		typename Mapping<dim>::InternalDataBase & /*mapping_data*/,
		typename Mapping<dim>::InternalDataBase & /*fedata*/,
		FEValuesData<dim,dim> & /*data*/,
		CellSimilarity::Similarity & /*cell_similarity*/) const
{
  Assert(false,ExcMessage(zero_dof_message));
}



template <int dim>
void
FE_Nothing<dim>::
fill_fe_face_values (const Mapping<dim> & /*mapping*/,
		     const typename Triangulation<dim>::cell_iterator & /*cell*/,
		     const unsigned int /*face*/,
		     const Quadrature<dim-1> & /*quadrature*/,
		     typename Mapping<dim>::InternalDataBase & /*mapping_data*/,
		     typename Mapping<dim>::InternalDataBase & /*fedata*/,
		     FEValuesData<dim,dim> & /*data*/) const
{
  Assert(false,ExcMessage(zero_dof_message));
}



template <int dim>
void
FE_Nothing<dim>::
fill_fe_subface_values (const Mapping<dim> & /*mapping*/,
			const typename Triangulation<dim>::cell_iterator & /*cell*/,
			const unsigned int /*face*/,
			const unsigned int /*subface*/,
			const Quadrature<dim-1> & /*quadrature*/,
			typename Mapping<dim>::InternalDataBase & /*mapping_data*/,
			typename Mapping<dim>::InternalDataBase & /*fedata*/,
			FEValuesData<dim,dim> & /*data*/) const
{
  Assert(false,ExcMessage(zero_dof_message));
}



template class FE_Nothing<deal_II_dimension>;

DEAL_II_NAMESPACE_CLOSE

