//----------------------------  mapping_q1.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mapping_q1.h  ---------------------------
#ifndef __deal2__mapping_q1_h
#define __deal2__mapping_q1_h


#include <cmath>
#include <fe/mapping.h>

/**
 * Mapping of general quadrilateral/hexahedra by d-linear shape
 * functions.
 *
 * This function maps the unit cell to a general grid cell with
 * straight lines in @p{d} dimensions (remark that in 3D the surfaces
 * may be curved, even if the edges are not). This is the well-known
 * mapping for polyhedral domains.
 *
 * Shape function for this mapping are the same as for the finite
 * element @p{FE_Q} of order 1. Therefore, coupling these two yields
 * an isoparametric element.
 *
 * @author Guido Kanschat, Ralf Hartmann, 2000, 2001
 */
template <int dim>
class MappingQ1 : public Mapping<dim>
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
	InternalData(unsigned int n_shape_functions);

					 /**
					  * Shape function at quadrature
					  * point. Shape functions are
					  * in tensor product order, so
					  * vertices must be reordered
					  * to obtain transformation.
					  */
	double shape (unsigned int qpoint,
		      unsigned int shape_nr) const;
	
					 /**
					  * Shape function at quadrature
					  * point. See above.
					  */
	double &shape (unsigned int qpoint,
		       unsigned int shape_nr);
	
					 /**
					  * Gradient of shape function
					  * in quadrature point. See
					  * above.
					  */
	Tensor<1,dim> derivative (unsigned int qpoint,
				  unsigned int shape_nr) const;

					 /**
					  * Gradient of shape function
					  * in quadrature point. See
					  * above.
					  */
	Tensor<1,dim> &derivative (unsigned int qpoint,
				   unsigned int shape_nr);
	
					 /**
					  * Values of shape
					  * functions. Access by
					  * function @p{shape}.
					  *
					  * Computed once.
					  */
	std::vector<double> shape_values;
	
					 /**
					  * Values of shape function
					  * derivatives. Access by
					  * function @p{derivative}.
					  *
					  * Computed once.
					  */
	std::vector<Tensor<1,dim> > shape_derivatives;
	
					 /**
					  * Tensors of covariant
					  * transformation.
					  *
					  * Computed on each cell.
					  */
	std::vector<Tensor<2,dim> > covariant;
	
					 /**
					  * Tensors of covariant
					  * transformation.
					  *
					  * Computed on each cell.
					  */
	std::vector<Tensor<2,dim> > contravariant;
	
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
					  * Auxuliary vectors for internal use.
					  */
        std::vector<std::vector<Tensor<1,dim> > > aux;

					 /**
					  * Number of shape functions.
					  */
	unsigned int n_shape_functions;

					 /**
					  * Stores the support points of
					  * the mapping shape functions on
					  * the @p{cell_of_current_support_points}.
					  */
	std::vector<Point<dim> > mapping_support_points;
	
					 /**
					  * Stores the cell of which the
					  * @p{mapping_support_points} are
					  * stored.
					  */
	DoFHandler<dim>::cell_iterator cell_of_current_support_points;
	
					 /**
					  * Default value of this flag
					  * is @p{true}. If @p{*this}
					  * is an object of a derived
					  * class, this flag is set to
					  * @p{false}.
					  */
	bool is_mapping_q1_data;
    };
    
				     /**
				      * Do the computations for the
				      * @p{get_face_data}
				      * functions. Here, the data
				      * vectors of @p{InternalData}
				      * are reinitialized to proper
				      * size and shape values and
				      * derivatives are
				      * computed. Furthermore
				      * @p{unit_tangential} vectors of
				      * the face are computed.
				      */
    void compute_face_data (const UpdateFlags flags,
			    const Quadrature<dim> &quadrature,
			    const unsigned int n_orig_q_points,
			    InternalData &data) const;
    
				     /**
				      * Mapping between tensor product
				      * ordering and real ordering of
				      * the vertices.
				      */
    static const unsigned int vertex_mapping[1<<dim];
    
				     /**
				      * Compute shape values and/or
				      * derivatives.
				      *
				      * Calles either the
				      * @p{compute_shapes_virtual} of
				      * this class or that of the
				      * derived class, depending on
				      * whether
				      * @p{data.is_mapping_q1_data}
				      * equals @p{true} or @p{false}.
				      */
    void compute_shapes (const std::vector<Point<dim> > &unit_points,
			 InternalData &data) const;
    
				     /**
				      * Do the computations for the @p{get_data}
				      * functions. Here, the data
				      * vectors of @p{InternalData} are
				      * reinitialized to proper size and
				      * shape values are computed.
				      */
    void compute_data (const UpdateFlags flags,
		       const Quadrature<dim>& quadrature,
		       const unsigned int n_orig_q_points,
		       InternalData& data) const;
    
				     /**
				      * Do the computation for the
				      * @p{fill_*} functions.
				      */
    void compute_fill (const typename DoFHandler<dim>::cell_iterator &cell,
		       const unsigned int   npts,
		       const unsigned int   offset,
		       InternalData        &data,
		       std::vector<Point<dim> > &quadrature_points) const;
    
				     /**
				      * Do the computation for the
				      * @p{fill_*} functions.
				      */
    void compute_fill_face (const typename DoFHandler<dim>::cell_iterator &cell,
			    const unsigned int      face_no,
			    const bool              is_subface,
			    const unsigned int      npts,
			    const unsigned int      offset,
			    const std::vector<double>   &weights,
			    InternalData           &mapping_data,
			    std::vector<Point<dim> >    &quadrature_points,
			    std::vector<double>         &JxW_values,
			    std::vector<Tensor<1,dim> > &boundary_form,
			    std::vector<Point<dim> > &normal_vectors) const;

				     /**
				      * Compute shape values and/or
				      * derivatives.
				      */
    virtual void compute_shapes_virtual (const std::vector<Point<dim> > &unit_points,
					 InternalData &data) const;
  private:

				     /**
				      * Computes the support points of
				      * the mapping. For @p{MappingQ1}
				      * these are the
				      * vertices.
				      */
    virtual void compute_mapping_support_points(
      const typename Triangulation<dim>::cell_iterator &cell,
      std::vector<Point<dim> > &a) const;

				     /**
				      *Number of shape functions
				      */
    static const unsigned int n_shape_functions = 1 << dim;
};


#endif
