//----------------------------  mapping_q.h  ---------------------------
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
//----------------------------  mapping_q.h  ---------------------------
#ifndef __deal2__mapping_q_h
#define __deal2__mapping_q_h


#include <fe/mapping_q1.h>
#include <grid/tria.h>
#include <grid/tria_boundary.h>


template <int dim> class TensorProductPolynomials;
class LagrangeEquidistant;



//TODO: fill_fe_face_values should exist in a version doing all faces
//TODO: to save initialization time.

//TODO: doc: laplace_on_quad_vector, compute_laplace_on_quad, etc all
//TODO: reference to each other, there is no description what they actually
//TODO: do or are good for

/**
 * Mapping class that uses Qp-mappings on boundary AND on inner
 * cells. The mapping shape functions make use of tensor product
 * polynomials with equidistant support points.
 *
 * Make sure elsewhere (e.g. in FEValues) that on inner cells only Q1
 * mappings are used.
 *  
 * @author Ralf Hartmann, Guido Kanschat 2000, 2001
 */
template <int dim>
class MappingQ : public MappingQ1<dim>
{
  public:
				     /**
				      * Constructor.  @p{p} gives the
				      * degree of mapping polynomials
				      * on boundary cells.
				      */
    MappingQ (unsigned int p);

				     /**
				      * Destructor.
				      */
    ~MappingQ ();

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
			 std::vector<Tensor<1,dim> >        &exterior_form,
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
			    std::vector<Tensor<1,dim> >        &exterior_form,
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
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual UpdateFlags update_each (const UpdateFlags) const;    

  private:
				     /** 
				      * Storage for internal data of
				      * Q_degree transformation.
				      */
    class InternalData : public MappingQ1<dim>::InternalData
    {
      public:
					 /**
					  * Constructor.
					  */
	InternalData(unsigned int n_shape_functions);
	
					 /**
					  * Unit normal vectors. Used
					  * for the alternative
					  * computation of the normal
					  * vectors. See doc of the
					  * @p{alternative_normals_computation}
					  * flag.
					  *
					  * Filled (hardcoded) once in
					  * @p{get_face_data}.
					  */
        std::vector<std::vector<Point<dim> > > unit_normals;

					 /**
					  * Flag that is set by the
					  * @p{fill_fe_[[sub]face]_values}
					  * function.
					  *
					  * If this flag is @p{true}
					  * we are on an interior cell
					  * and the
					  * @p{mapping_q1_data} is
					  * used.
					  */
	bool use_mapping_q1_on_current_cell;
	
					 /**
					  * On interior cells
					  * @p{MappingQ1} is used.
					  */
	MappingQ1<dim>::InternalData mapping_q1_data;
    };

				     /**
				      * Calles the
				      * @p{compute_face_data} function
				      * of its base @ref{MappingQ1}
				      * class.
				      *
				      * For the
				      * @p{alternative_normal_computation}
				      * also the @p{unit_normal}
				      * vectors of the face are
				      * computed.
				      */
    void compute_face_data (const UpdateFlags flags,
			    const Quadrature<dim>& quadrature,
			    const unsigned int n_orig_q_points,
			    MappingQ1<dim>::InternalData& data) const;
    
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
			    MappingQ1<dim>::InternalData &mapping_q1_data,
			    std::vector<Point<dim> >    &quadrature_points,
			    std::vector<double>         &JxW_values,
			    std::vector<Tensor<1,dim> > &boundary_form,
			    std::vector<Point<dim> >    &normal_vectors) const;
    
				     /**
				      * Compute shape values and/or
				      * derivatives.
				      */
    virtual void compute_shapes_virtual (const std::vector<Point<dim> > &unit_points,
					 MappingQ1<dim>::InternalData &data) const;

				     /**
				      * This function is needed by the
				      * constructor of @p{MappingQ<dim>}
				      * for @p{dim=} 2 and 3.
				      *
				      * For @p{degree<4} this function
				      * sets the
				      * @p{laplace_on_quad_vector} to
				      * the hardcoded data. For
				      * @p{degree>=4} and MappingQ<2>
				      * this vector is computed.
				      */
    void set_laplace_on_quad_vector(std::vector<std::vector<double> > &loqvs) const;
    
				     /**
				      * This function is needed by the
				      * constructor of @p{MappingQ<3>}.
				      *
				      * For @p{degree==2} this function
				      * sets the
				      * @p{laplace_on_hex_vector} to
				      * the hardcoded data. For
				      * @p{degree>2} this vector is
				      * computed.
				      */
    void set_laplace_on_hex_vector(std::vector<std::vector<double> > &lohvs) const;
    
				     /**
				      * Computes the @p{laplace_on_quad(hex)_vector}.
				      *
				      * Called by the
				      * @p{set_laplace_on_(quad)hex_vector}
				      * functions if the data is not
				      * yet hardcoded.
				      */
    void compute_laplace_vector(std::vector<std::vector<double> > &lvs) const;

				     /**
				      * Takes a
				      * @p{set_laplace_on_hex(quad)_vector}
				      * and applies it to the vector
				      * @p{a} to compute the inner
				      * support points. They are
				      * appended to the vector @p{a}.
				      */
    void apply_laplace_vector(const std::vector<std::vector<double> > &lvs,
			      std::vector<Point<dim> > &a) const;
    
				     /**
				      * Computes the support points of
				      * the mapping.
				      */
    virtual void compute_mapping_support_points(
      const typename Triangulation<dim>::cell_iterator &cell,
      std::vector<Point<dim> > &a) const;

				     /**
				      * Computes all support points of
				      * the mapping shape
				      * functions. The inner support
				      * points (ie. support points in
				      * quads for 2d, in hexes for 3d)
				      * are computed using the
				      * solution of a Laplace equation
				      * with the position of the outer
				      * support points as boundary
				      * values. The outer support
				      * points are all support points
				      * except of the inner ones.
				      */
    void compute_support_points_laplace(
      const typename Triangulation<dim>::cell_iterator &cell,
      std::vector<Point<dim> > &a) const;
    
				     /**
				      * Simple version of the
				      * @p{compute_support_points_laplace}
				      * function. Does not use the
				      * solution to Laplace
				      * equation. Computes the inner
				      * support points by simple
				      * interpolations.*/
    void compute_support_points_simple(
      const typename Triangulation<dim>::cell_iterator &cell,
      std::vector<Point<dim> > &a) const;
    
				     /**
				      * For @p{dim=2,3}. Adds (appends) the
				      * support points of all lines to
				      * the vector a.
				      *
				      * Needed by the
				      * @p{compute_support_points_simple(laplace)}
				      * functions. For @p{dim=1} this
				      * function is empty.
				      */
    void add_line_support_points (const Triangulation<dim>::cell_iterator &cell,
				  std::vector<Point<dim> > &a) const;

				     /**
				      * For @p{dim=3}. Adds (appends) the
				      * support points of all faces (quads in 3d) to
				      * the vector a.
				      *
				      * Needed by the
				      * @p{compute_support_points_laplace}
				      * function. For @p{dim=1} and 2 this
				      * function is empty.
				      */
//TODO: rename function to add_quad_support_points, to unify notation    
    void add_face_support_points(const typename Triangulation<dim>::cell_iterator &cell,
				 std::vector<Point<dim> > &a) const;
    
				     /**
				      * For @p{dim=2} and 3. Simple
				      * version of the
				      * @p{add_face_support_points}
				      * function.
				      *
				      * Needed by the
				      * @p{compute_support_points_simple}
				      */
    void fill_quad_support_points_simple (const Triangulation<dim>::cell_iterator &cell,
					  std::vector<Point<dim> > &a) const;
    
				     /**
				      * Needed by the
				      * @p{laplace_on_quad} function
				      * (for @p{dim==2}). Filled by the
				      * constructor.
				      *
				      * Sizes:
				      * laplace_on_quad_vector.size()=
				      *   number of inner
				      *   unit_support_points
				      * laplace_on_quad_vector[i].size()=
				      *   number of outer
				      *   unit_support_points, i.e.
				      *   unit_support_points on the
				      *   boundary of the quad
				      */
    std::vector<std::vector<double> > laplace_on_quad_vector;
    
				     /**
				      * Needed by the
				      * @p{laplace_on_hex} function
				      * (for @p{dim==3}). Filled by the
				      * constructor.
				      */
    std::vector<std::vector<double> > laplace_on_hex_vector;

				     /**
				      * Exception.
				      */
    DeclException1 (ExcLaplaceVectorNotSet,
		    int,
		    << "laplace_vector not set for degree=" << arg1 << ".");
     
				     /**
				      * Degree @p{p} of the
				      * polynomials used as shape
				      * functions for the Qp mapping
				      * of cells at the boundary.
				      */  
    const unsigned int degree;

				     /**
				      * Number of inner mapping shape
				      * functions.
				      */
    const unsigned int n_inner;

				     /**
				      * Number of mapping shape
				      * functions on the boundary.
				      */
    const unsigned int n_outer;
    
				     /**
				      * Pointer to the
				      * @p{dim}-dimensional tensor
				      * product polynomials used as
				      * shape functions for the Qp
				      * mapping of cells at the
				      * boundary.
				      */
    TensorProductPolynomials<dim> *tensor_pols;
    
    				     /**
				      * Number of the Qp tensor
				      * product shape functions.
				      */
    unsigned int n_shape_functions;

				     /**
				      * Mapping from lexicographic to
				      * to the Qp shape function
				      * numbering. Its size is
				      * @p{dofs_per_cell}.
				      */
    std::vector<unsigned int> renumber;

				     /**
				      * Needed for inner faces.
				      */
//TODO: can we make this variable static?    
    StraightBoundary<dim> straight_boundary;

				     /**
				      * Flag for computing the normal
				      * vectors directly by using a
				      * covariant transformation.
				      * Used to test the covariant
				      * transformation.
				      */
//TODO: why have two ways to compute? if they both work, choose one and remove the other    
    bool alternative_normals_computation;

				     /**
				      * If this flag is set @p{true}
				      * then @p{MappingQ} is used on
				      * all cells, not only on
				      * boundary cells.
				      *
				      * The default value is false.
				      */
//TODO: remove use_mapping_q_on_all_cells as it is set to false in the constructor and never set again    
    bool use_mapping_q_on_all_cells;
};


#endif
