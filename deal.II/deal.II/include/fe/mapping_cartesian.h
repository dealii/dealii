//----------------------------  mapping_Cartesian.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mapping_Cartesian.h  ---------------------------
#ifndef __deal2__mapping_cartesian_h
#define __deal2__mapping_cartesian_h


#include <cmath>
#include <fe/mapping.h>

/**
 * Mapping of an axis-parallel cell.
 *
 * This class maps the unit cell to a grid cell with surfaces parallel
 * to the coordinate lines/planes. It is specifically developed for
 * cartesian meshes. Apply this mapping to a general mesh to get
 * strange results.
 *
 * @author Guido Kanschat, 2001
 */
template <int dim>
class MappingCartesian : public Mapping<dim>
{
  public:
				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual Mapping<dim>::InternalDataBase*
    get_data (const UpdateFlags,
	      const Quadrature<dim>& quadrature) const;

				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual Mapping<dim>::InternalDataBase*
    get_face_data (const UpdateFlags flags,
		   const Quadrature<dim-1>& quadrature) const;

				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual Mapping<dim>::InternalDataBase*
    get_subface_data (const UpdateFlags flags,
		      const Quadrature<dim-1>& quadrature) const;

				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual void
    fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
		    const Quadrature<dim>& quadrature,
		    Mapping<dim>::InternalDataBase &mapping_data,
		    std::vector<Point<dim> >        &quadrature_points,
		    std::vector<double>             &JxW_values) const ;

				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual void
    fill_fe_face_values (const typename DoFHandler<dim>::cell_iterator &cell,
			 const unsigned int face_no,
			 const Quadrature<dim-1>& quadrature,
			 typename Mapping<dim>::InternalDataBase &mapping_data,
			 std::vector<Point<dim> >        &quadrature_points,
			 std::vector<double>             &JxW_values,
			 std::vector<Tensor<1,dim> >        &boundary_form,
			 std::vector<Point<dim> >        &normal_vectors) const ;

				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual void
    fill_fe_subface_values (const typename DoFHandler<dim>::cell_iterator &cell,
			    const unsigned int face_no,
			    const unsigned int sub_no,
			    const Quadrature<dim-1>& quadrature,
			    typename Mapping<dim>::InternalDataBase &mapping_data,
			    std::vector<Point<dim> >        &quadrature_points,
			    std::vector<double>             &JxW_values,
			    std::vector<Tensor<1,dim> >        &boundary_form,
			    std::vector<Point<dim> >        &normal_vectors) const ;


				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual void
    transform_covariant (std::vector<Tensor<1,dim> >       &dst,
			 const std::vector<Tensor<1,dim> > &src,
			 const Mapping<dim>::InternalDataBase &mapping_data,
			 const unsigned int src_offset) const;
    
				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual void
    transform_contravariant (std::vector<Tensor<1,dim> >       &dst,
			     const std::vector<Tensor<1,dim> > &src,
			     const Mapping<dim>::InternalDataBase &mapping_data,
			     const unsigned int src_offset) const;

				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual void
    transform_covariant (std::vector<Point<dim> >       &dst,
			 const std::vector<Point<dim> > &src,
			 const Mapping<dim>::InternalDataBase &mapping_data,
			 const unsigned int src_offset) const;
    
				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual void
    transform_contravariant (std::vector<Point<dim> >       &dst,
			     const std::vector<Point<dim> > &src,
			     const Mapping<dim>::InternalDataBase &mapping_data,
			     const unsigned int src_offset) const;
    
				     /**
				      * Transforms the point @p{p} on
				      * the unit cell to the point
				      * @p{p_real} on the real cell
				      * @p{cell} and returns @p{p_real}.
				      */
//TODO: document meaning of mdata argument    
    virtual Point<dim> transform_unit_to_real_cell (
      const typename Triangulation<dim>::cell_iterator cell,
      const Point<dim> &p,
      const typename Mapping<dim>::InternalDataBase *const mdata=0) const;

				     /**
				      * Transforms the point @p{p} on
				      * the real cell to the point
				      * @p{p_unit} on the unit cell
				      * @p{cell} and returns @p{p_unit}.
				      *
				      * Uses Newton iteration and the
				      * @p{transform_unit_to_real_cell}
				      * function.
				      */
    virtual Point<dim> transform_real_to_unit_cell (
      const typename Triangulation<dim>::cell_iterator cell,
      const Point<dim> &p) const;
    
				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual UpdateFlags update_once (const UpdateFlags) const;
    
				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual UpdateFlags update_each (const UpdateFlags) const;
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidData);

  protected:
    				     /**
				      * Implementation of the
				      * covariant transformation.
				      */
    template <typename tensor_>
    void covariant_transformation (std::vector<tensor_>       &dst,
				   const std::vector<tensor_> &src,
				   const Mapping<dim>::InternalDataBase &mapping_data,
				   const unsigned int src_offset) const;
    
				     /**
				      * Implementation of the
				      * contravariant transformation.
				      */
    template <typename tensor_>
    void contravariant_transformation (std::vector<tensor_>       &dst,
				       const std::vector<tensor_> &src,
				       const Mapping<dim>::InternalDataBase &mapping_data,
				       const unsigned int src_offset) const;

				     /** 
				      * Storage for internal data of
				      * d-linear transformation.
				      */
    class InternalData : public Mapping<dim>::InternalDataBase
    {
      public:
					 /**
					  * Constructor.
					  */
	InternalData (const Quadrature<dim> &quadrature);

					 /**
					  * Length of the cell in
					  * different coordinate
					  * directions.
					  */
	Tensor<1,dim> length;

					 /**
					  * Vector of all quadrature
					  * points. Especially, all
					  * points of all faces.
					  */
	std::vector<Point<dim> > quadrature_points;
	
					 /**
					  * Unit tangential vectors. Used
					  * for the computation of
					  * boundary forms and normal
					  * vectors.
					  *
					  * Filled once.
					  */
        std::vector<std::vector<Tensor<1,dim> > > unit_tangentials;
	
					 /**
					  * Auxiliary vectors for internal use.
					  */
        std::vector<std::vector<Tensor<1,dim> > > aux;
    };
    
				     /**
				      * Do the computation for the
				      * @p{fill_*} functions.
				      */
    void compute_fill (const typename DoFHandler<dim>::cell_iterator
		       &cell,
		       unsigned int face_no,
		       unsigned int sub_no,
		       InternalData& data,
		       std::vector<Point<dim> > &quadrature_points,
		       std::vector<Point<dim> >& normal_vectors) const;

  private:
				     /**
				      * Value to indicate that a given
				      * face or subface number is
				      * invalid.
				      */
    static const unsigned int invalid_face_number = static_cast<unsigned int>(-1);    
};


#endif
