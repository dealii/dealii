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
#ifndef __deal2__fe_nothing_h
#define __deal2__fe_nothing_h

#include <base/config.h>
#include <fe/fe.h>

DEAL_II_NAMESPACE_OPEN


/*!@addtogroup fe */
/*@{*/

/**
 * Definition of a finite element with zero degrees of freedom.  This
 * class is useful (in the context of an hp method) to represent empty
 * cells in the triangulation on which no degrees of freedom should
 * be allocated.  Thus a triangulation may be divided into two regions:
 * an active region where normal elements are used, and an inactive
 * region where FE_Nothing elements are used.  The hp::DoFHandler will
 * therefore assign no degrees of freedom to the FE_Nothing cells, and
 * this subregion is therefore implicitly deleted from the computation.
 *
 * @author Joshua White
 */
template <int dim>
class FE_Nothing : public FiniteElement<dim>
{
  public:

    FE_Nothing ();

    virtual
    FiniteElement<dim> *
    clone() const;

    virtual
    std::string
    get_name() const;

    virtual
    unsigned int
    n_base_elements () const;

    virtual
    const FiniteElement<dim> &
    base_element (const unsigned int index) const;

    virtual
    unsigned int
    element_multiplicity (const unsigned int index) const;

    virtual
    UpdateFlags
    update_once (const UpdateFlags flags) const;

    virtual
    UpdateFlags
    update_each (const UpdateFlags flags) const;

    virtual
    double
    shape_value (const unsigned int i, const Point<dim> &p) const;

    virtual
    void
    fill_fe_values (const Mapping<dim> & mapping,
		    const typename Triangulation<dim>::cell_iterator & cell,
		    const Quadrature<dim> & quadrature,
		    typename Mapping<dim>::InternalDataBase & mapping_data,
		    typename Mapping<dim>::InternalDataBase & fedata,
		    FEValuesData<dim,dim> & data,
		    CellSimilarity::Similarity & cell_similarity) const;

    virtual
    void
    fill_fe_face_values (const Mapping<dim> & mapping,
			 const typename Triangulation<dim> :: cell_iterator & cell,
			 const unsigned int face,
			 const Quadrature<dim-1> & quadrature,
			 typename Mapping<dim> :: InternalDataBase & mapping_data,
			 typename Mapping<dim> :: InternalDataBase & fedata,
			 FEValuesData<dim,dim> & data) const;

    virtual
    void
    fill_fe_subface_values (const Mapping<dim> & mapping,
			    const typename Triangulation<dim>::cell_iterator & cell,
			    const unsigned int face,
			    const unsigned int subface,
			    const Quadrature<dim-1> & quadrature,
			    typename Mapping<dim>::InternalDataBase & mapping_data,
			    typename Mapping<dim>::InternalDataBase & fedata,
			    FEValuesData<dim,dim> & data) const;

    virtual
    typename Mapping<dim>::InternalDataBase *
    get_data (const UpdateFlags     update_flags,
	      const Mapping<dim>    & mapping,
	      const Quadrature<dim> & quadrature) const;
};


/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif

