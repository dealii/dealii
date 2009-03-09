//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__mapping_cartesian_h
#define __deal2__mapping_cartesian_h


#include <base/config.h>
#include <base/table.h>
#include <cmath>
#include <fe/mapping.h>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mapping */
/*@{*/

/**
 * Mapping of an axis-parallel cell.
 *
 * This class maps the unit cell to a grid cell with surfaces parallel
 * to the coordinate lines/planes. The mapping is therefore a scaling
 * along the coordinate directions. It is specifically developed for
 * cartesian meshes. Apply this mapping to a general mesh to get
 * strange results.
 * 
 * For more information about the <tt>spacedim</tt> template parameter
 * check the documentation of FiniteElement or the one of
 * Triangulation.
 *
 * @author Guido Kanschat, 2001; Ralf Hartmann, 2005
 */
template <int dim, int spacedim=dim>
class MappingCartesian : public Mapping<dim,spacedim>
{
  public:
    virtual
    typename Mapping<dim, spacedim>::InternalDataBase *
    get_data (const UpdateFlags,
	      const Quadrature<dim>& quadrature) const;

    virtual
    typename Mapping<dim, spacedim>::InternalDataBase *
    get_face_data (const UpdateFlags flags,
		   const Quadrature<dim-1>& quadrature) const;

    virtual
    typename Mapping<dim, spacedim>::InternalDataBase *
    get_subface_data (const UpdateFlags flags,
		      const Quadrature<dim-1>& quadrature) const;

    virtual void
    fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
		    const Quadrature<dim>                                     &quadrature,
		    enum CellSimilarity::Similarity                           &cell_similarity,
		    typename Mapping<dim, spacedim>::InternalDataBase         &mapping_data,
		    std::vector<Point<spacedim> >                             &quadrature_points,
		    std::vector<double>                                       &JxW_values,
		    std::vector<Tensor<2,spacedim> >                          &jacobians,
		    std::vector<Tensor<3,spacedim> >                          &jacobian_grads,
		    std::vector<Tensor<2,spacedim> >                          &inverse_jacobians,
	 	    std::vector<Point<spacedim> >   &) const ;


    virtual void
    fill_fe_face_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
			 const unsigned int face_no,
			 const Quadrature<dim-1>& quadrature,
			 typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
			 std::vector<Point<dim> >        &quadrature_points,
			 std::vector<double>             &JxW_values,
			 std::vector<Tensor<1,dim> >        &boundary_form,
			 std::vector<Point<spacedim> >        &normal_vectors) const ;
    virtual void
    fill_fe_subface_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
			    const unsigned int face_no,
			    const unsigned int sub_no,
			    const Quadrature<dim-1>& quadrature,
			    typename Mapping<dim, spacedim>::InternalDataBase &mapping_data,
			    std::vector<Point<dim> >        &quadrature_points,
			    std::vector<double>             &JxW_values,
			    std::vector<Tensor<1,dim> >        &boundary_form,
			    std::vector<Point<spacedim> >        &normal_vectors) const ;

    virtual void
    transform (const VectorSlice<const std::vector<Tensor<1,dim> > > input,
               VectorSlice<std::vector<Tensor<1,spacedim> > > output,
               const typename Mapping<dim,spacedim>::InternalDataBase &internal,
	       const MappingType type) const;

    virtual void
    transform (const VectorSlice<const std::vector<Tensor<2,dim> > > input,
               VectorSlice<std::vector<Tensor<2,spacedim> > > output,
               const typename Mapping<dim,spacedim>::InternalDataBase &internal,
	       const MappingType type) const;    
    
    virtual Point<spacedim>
    transform_unit_to_real_cell (
      const typename Triangulation<dim,spacedim>::cell_iterator &cell,
      const Point<dim>                                 &p) const;

				     /**
				      * Transforms the point @p p on
				      * the real cell to the point
				      * @p p_unit on the unit cell
				      * @p cell and returns @p p_unit.
				      *
				      * Uses Newton iteration and the
				      * @p transform_unit_to_real_cell
				      * function.
				      */
    virtual Point<dim>
    transform_real_to_unit_cell (
      const typename Triangulation<dim,spacedim>::cell_iterator &cell,
      const Point<spacedim>                            &p) const;
    

                                     /**
                                      * Return a pointer to a copy of the
                                      * present object. The caller of this
                                      * copy then assumes ownership of it.
                                      */
    virtual
    Mapping<dim, spacedim> * clone () const;
    
  protected:
				     /** 
				      * Storage for internal data of
				      * the scaling.
				      */
    class InternalData : public Mapping<dim, spacedim>::InternalDataBase
    {
      public:
					 /**
					  * Constructor.
					  */
	InternalData (const Quadrature<dim> &quadrature);

					 /**
					  * Return an estimate (in
					  * bytes) or the memory
					  * consumption of this
					  * object.
					  */
	virtual unsigned int memory_consumption () const;

					 /**
					  * Length of the cell in
					  * different coordinate
					  * directions, <i>h<sub>x</sub></i>,
					  * <i>h<sub>y</sub></i>, <i>h<sub>z</sub></i>.
					  */
	Tensor<1,dim> length;

					 /**
					  * The volume element
					  */
	double volume_element;
	
					 /**
					  * Vector of all quadrature
					  * points. Especially, all
					  * points on all faces.
					  */
	std::vector<Point<dim> > quadrature_points;
    };
    
				     /**
				      * Do the computation for the
				      * <tt>fill_*</tt> functions.
				      */
    void compute_fill (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
		       const unsigned int face_no,
		       const unsigned int sub_no,
		       const enum CellSimilarity::Similarity cell_similarity,
		       InternalData& data,
		       std::vector<Point<dim> > &quadrature_points,
		       std::vector<Point<dim> >& normal_vectors) const;

  private:
    virtual UpdateFlags update_once (const UpdateFlags) const;    
    virtual UpdateFlags update_each (const UpdateFlags) const;
    
				     /**
				      * Value to indicate that a given
				      * face or subface number is
				      * invalid.
				      */
    static const unsigned int invalid_face_number = numbers::invalid_unsigned_int;    
};

/*@}*/

/* -------------- declaration of explicit specializations ------------- */

DEAL_II_NAMESPACE_CLOSE

#endif
