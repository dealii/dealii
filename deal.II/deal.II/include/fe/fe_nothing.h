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
 * Note that some care must be taken that the resulting mesh topology
 * continues to make sense when FE_Nothing elements are introduced.
 * This is particularly true when dealing with hanging node constraints,
 * because the library makes some basic assumptions about the nature
 * of those constraints.  The following geometries are acceptable:
 *
 * +---------+----+----+
 * |         | 0  |    |
 * |    1    +----+----+
 * |         | 0  |    |
 * +---------+----+----+
 *
 * +---------+----+----+
 * |         | 1  |    |
 * |    0    +----+----+
 * |         | 1  |    |
 * +---------+----+----+
 *
 * Here, 0 denotes an FE_Nothing cell, and 1 denotes some other
 * element type.  The library has no difficulty computing the necessary
 * hanging node constraints in these cases (i.e. no constraint).  
 * However, the following geometry is NOT acceptable (at least 
 * in the current implementation):
 *
 * +---------+----+----+
 * |         | 0  |    |
 * |    1    +----+----+
 * |         | 1  |    |
 * +---------+----+----+
 * 
 * The distinction lies in the mixed nature of the child faces,
 * a case we have not implemented as of yet.
 *
 * @author Joshua White
 */
template <int dim>
class FE_Nothing : public FiniteElement<dim>
{
  public:

                                    /**
                                      * Constructor. Argument denotes the
                                      * number of components to give this
                                      * finite element (default = 1).  
                                      */
    FE_Nothing (unsigned int n_components = 1);
    
                                     /**
                                      * A sort of virtual copy
                                      * constructor. Some places in
                                      * the library, for example the
                                      * constructors of FESystem as
                                      * well as the hp::FECollection
                                      * class, need to make copied of
                                      * finite elements without
                                      * knowing their exact type. They
                                      * do so through this function.
                                      */
    virtual
    FiniteElement<dim> *
    clone() const;

                                     /**
                                      * Return a string that uniquely
                                      * identifies a finite
                                      * element. In this case it is
                                      * <code>FE_Nothing@<dim@></code>.
                                      */
    virtual
    std::string
    get_name() const;

                                     /**
                                      * Number of base elements in a
                                      * mixed discretization. In this case
                                      * we only have one.
                                      */
    virtual
    unsigned int
    n_base_elements () const;
    
                                     /**
                                      * Access to base element
                                      * objects. Since the element is
                                      * scalar, then
                                      * <code>base_element(0)</code> is
                                      * @p this.
                                      */
    virtual
    const FiniteElement<dim> &
    base_element (const unsigned int index) const;

                                     /**
                                      * This index denotes how often
                                      * the base element @p index is
                                      * used in a composed element. If
                                      * the element is scalar, then
                                      * the result is always equal to
                                      * one. See the documentation for
                                      * the n_base_elements()
                                      * function for more details.
                                      */
    virtual
    unsigned int
    element_multiplicity (const unsigned int index) const;

                                     /**
                                      * Determine the values a finite
                                      * element should compute on
                                      * initialization of data for
                                      * FEValues.
                                      *
                                      * Given a set of flags
                                      * indicating what quantities are
                                      * requested from a FEValues
                                      * object, update_once() and
                                      * update_each() compute which
                                      * values must really be
                                      * computed. Then, the
                                      * <tt>fill_*_values</tt> functions
                                      * are called with the result of
                                      * these.
                                      *
                                      * In this case, since the element 
                                      * has zero degrees of freedom and
                                      * no information can be computed on
                                      * it, this function simply returns
                                      * the default (empty) set of update
                                      * flags.
                                      */

    virtual
    UpdateFlags
    update_once (const UpdateFlags flags) const;

                                     /**
                                      * Complementary function for
                                      * update_once().
                                      *
                                      * While update_once() returns
                                      * the values to be computed on
                                      * the unit cell for yielding the
                                      * required data, this function
                                      * determines the values that
                                      * must be recomputed on each
                                      * cell.
                                      *
                                      * Refer to update_once() for
                                      * more details.
                                      */
    virtual
    UpdateFlags
    update_each (const UpdateFlags flags) const;

                                     /**
                                      * Return the value of the
                                      * @p ith shape function at the
                                      * point @p p. @p p is a point
                                      * on the reference element. Because the
                                      * current element has no degrees of freedom,
                                      * this function should obviously not be
                                      * called in practice.  All this function
                                      * really does, therefore, is trigger an
                                      * exception.
                                      */
    virtual
    double
    shape_value (const unsigned int i, const Point<dim> &p) const;

                                     /**
                                      * Fill the fields of
                                      * FEValues. This function
                                      * performs all the operations
                                      * needed to compute the data of an
                                      * FEValues object.
                                      *
                                      * In the current case, this function
                                      * returns no meaningful information, 
                                      * since the element has no degrees of
                                      * freedom.
                                      */
    virtual
    void
    fill_fe_values (const Mapping<dim> & mapping,
		    const typename Triangulation<dim>::cell_iterator & cell,
		    const Quadrature<dim> & quadrature,
		    typename Mapping<dim>::InternalDataBase & mapping_data,
		    typename Mapping<dim>::InternalDataBase & fedata,
		    FEValuesData<dim,dim> & data,
		    CellSimilarity::Similarity & cell_similarity) const;

                                     /**
                                      * Fill the fields of
                                      * FEFaceValues. This function
                                      * performs all the operations
                                      * needed to compute the data of an
                                      * FEFaceValues object.
                                      *
                                      * In the current case, this function
                                      * returns no meaningful information, 
                                      * since the element has no degrees of
                                      * freedom.
                                      */
    virtual
    void
    fill_fe_face_values (const Mapping<dim> & mapping,
			 const typename Triangulation<dim> :: cell_iterator & cell,
			 const unsigned int face,
			 const Quadrature<dim-1> & quadrature,
			 typename Mapping<dim> :: InternalDataBase & mapping_data,
			 typename Mapping<dim> :: InternalDataBase & fedata,
			 FEValuesData<dim,dim> & data) const;

                                     /**
                                      * Fill the fields of
                                      * FESubFaceValues. This function
                                      * performs all the operations
                                      * needed to compute the data of an
                                      * FESubFaceValues object.
                                      *
                                      * In the current case, this function
                                      * returns no meaningful information, 
                                      * since the element has no degrees of
                                      * freedom.
                                      */
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

                                     /**
                                      * Prepare internal data
                                      * structures and fill in values
                                      * independent of the
                                      * cell. Returns a pointer to an
                                      * object of which the caller of
                                      * this function then has to
                                      * assume ownership (which
                                      * includes destruction when it
                                      * is no more needed).
                                      *
                                      * In the current case, this function 
                                      * just returns a default pointer, since
                                      * no meaningful data exists for this 
                                      * element.
                                      */
    virtual
    typename Mapping<dim>::InternalDataBase *
    get_data (const UpdateFlags     update_flags,
	      const Mapping<dim>    & mapping,
	      const Quadrature<dim> & quadrature) const;
        
                                     /**
                                      * Return whether this element dominates
                                      * the one given as argument when they
                                      * meet at a common face,
                                      * whether it is the other way around,
                                      * whether neither dominates, or if
                                      * either could dominate.
                                      *
                                      * For a definition of domination, see
                                      * FiniteElementBase::Domination and in
                                      * particular the @ref hp_paper "hp paper".
                                      *
                                      * In the current case, this element
                                      * is always assumed to dominate, unless 
                                      * it is also of type FE_Nothing().  In
                                      * that situation, either element can
                                      * dominate.
                                      */
    virtual
    FiniteElementDomination::Domination
    compare_for_face_domination (const FiniteElement<dim> & fe_other) const;
    
    
    
    virtual
    std::vector<std::pair<unsigned int, unsigned int> >
    hp_vertex_dof_identities (const FiniteElement<dim> &fe_other) const;
    
    virtual
    std::vector<std::pair<unsigned int, unsigned int> >
    hp_line_dof_identities (const FiniteElement<dim> &fe_other) const;
    
    virtual
    std::vector<std::pair<unsigned int, unsigned int> >
    hp_quad_dof_identities (const FiniteElement<dim> &fe_other) const;
    
    virtual
    bool
    hp_constraints_are_implemented () const;

                                      /**
                                      * Return the matrix
                                      * interpolating from a face of
                                      * of one element to the face of
                                      * the neighboring element. 
                                      * The size of the matrix is
                                      * then <tt>source.#dofs_per_face</tt> times
                                      * <tt>this->#dofs_per_face</tt>.
                                      *
                                      * Since the current finite element has no
                                      * degrees of freedom, the interpolation
                                      * matrix is necessarily empty.
                                      */

    virtual
    void
    get_face_interpolation_matrix (const FiniteElement<dim> &source_fe,
                                   FullMatrix<double>       &interpolation_matrix) const;
                                   
                                   
                                     /**
                                      * Return the matrix
                                      * interpolating from a face of
                                      * of one element to the subface of
                                      * the neighboring element. 
                                      * The size of the matrix is
                                      * then <tt>source.#dofs_per_face</tt> times
                                      * <tt>this->#dofs_per_face</tt>.
                                      *
                                      * Since the current finite element has no
                                      * degrees of freedom, the interpolation
                                      * matrix is necessarily empty.
                                      */
                                      
    virtual
    void
    get_subface_interpolation_matrix (const FiniteElement<dim> & source_fe,
                                      const unsigned int index,
                                      FullMatrix<double>  &interpolation_matrix) const;

    
};


/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif

