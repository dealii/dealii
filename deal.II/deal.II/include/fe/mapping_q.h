//----------------------------  mapping_q.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal.II authors
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


template <int dim> class TensorProductPolynomials;
class LagrangeEquidistant;




//TODO:[RH] (later) doc: laplace_on_quad_vector, compute_laplace_on_quad, etc all
//      see upcoming paper.


/**
 * Mapping class that uses Qp-mappings on boundary cells. The mapping
 * shape functions make use of tensor product polynomials with
 * equidistant (on the unit cell) support points.
 *
 * For more details about Qp-mappings, see the small `mapping' report
 * in the `Reports' section of `Documentation'.
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
    MappingQ (const unsigned int p);

				     /**
				      * Destructor.
				      */
    virtual ~MappingQ ();
    
				     /**
				      * Transforms the point @p{p} on
				      * the unit cell to the point
				      * @p{p_real} on the real cell
				      * @p{cell} and returns @p{p_real}.
				      */
    virtual Point<dim> transform_unit_to_real_cell (
      const typename Triangulation<dim>::cell_iterator cell,
      const Point<dim> &p) const;
    
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
    virtual void
    transform_covariant (typename std::vector<Tensor<1,dim> >       &dst,
			 const typename std::vector<Tensor<1,dim> > &src,
			 const typename Mapping<dim>::InternalDataBase &mapping_data,
			 const unsigned int src_offset) const;
    
				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual void
    transform_contravariant (typename std::vector<Tensor<1,dim> >       &dst,
			     const typename std::vector<Tensor<1,dim> > &src,
			     const typename Mapping<dim>::InternalDataBase &mapping_data,
			     const unsigned int src_offset) const;

				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual void
    transform_covariant (typename std::vector<Point<dim> >       &dst,
			 const typename std::vector<Point<dim> > &src,
			 const typename Mapping<dim>::InternalDataBase &mapping_data,
			 const unsigned int src_offset) const;
    
				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual void
    transform_contravariant (typename std::vector<Point<dim> >       &dst,
			     const typename std::vector<Point<dim> > &src,
			     const typename Mapping<dim>::InternalDataBase &mapping_data,
			     const unsigned int src_offset) const;

				     /**
				      * Return the degree of the
				      * mapping, i.e. the value which
				      * was passed to the constructor.
				      */
    unsigned int get_degree () const;
    
  protected:
				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual void
    fill_fe_values (const typename DoFHandler<dim>::cell_iterator &cell,
		    const Quadrature<dim>                &quadrature,
		    typename Mapping<dim>::InternalDataBase &mapping_data,
		    typename std::vector<Point<dim> >             &quadrature_points,
		    std::vector<double>                  &JxW_values) const ;

				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual void
    fill_fe_face_values (const typename DoFHandler<dim>::cell_iterator &cell,
			 const unsigned int face_no,
			 const Quadrature<dim-1>& quadrature,
			 typename Mapping<dim>::InternalDataBase &mapping_data,
			 typename std::vector<Point<dim> >        &quadrature_points,
			 std::vector<double>             &JxW_values,
			 typename std::vector<Tensor<1,dim> >        &exterior_form,
			 typename std::vector<Point<dim> >        &normal_vectors) const ;

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
			    typename std::vector<Point<dim> >        &quadrature_points,
			    std::vector<double>             &JxW_values,
			    typename std::vector<Tensor<1,dim> >        &exterior_form,
			    typename std::vector<Point<dim> >        &normal_vectors) const ;

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
	InternalData (const unsigned int n_shape_functions);
	

					 /**
					  * Return an estimate (in
					  * bytes) or the memory
					  * consumption of this
					  * object.
					  */
	virtual unsigned int memory_consumption () const;

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
        typename std::vector<typename std::vector<Point<dim> > > unit_normals;

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
	typename MappingQ1<dim>::InternalData mapping_q1_data;
    };

				     /**
				      * For @p{dim=2,3}. Append the
				      * support points of all shape
				      * functions located on bounding
				      * lines to the vector
				      * @p{a}. Points located on the
				      * line but on vertices are not
				      * included.
				      *
				      * Needed by the
				      * @p{compute_support_points_simple(laplace)}
				      * functions. For @p{dim=1} this
				      * function is empty.
				      *
				      * This function is made virtual
				      * in order to allow derived
				      * classes to choose shape
				      * function support points
				      * differently than the present
				      * class, which chooses the
				      * points as interpolation points
				      * on the boundary.
				      */
    virtual void
    add_line_support_points (const typename Triangulation<dim>::cell_iterator &cell,
			     typename std::vector<Point<dim> > &a) const;

				     /**
				      * For @p{dim=3}. Append the
				      * support points of all shape
				      * functions located on bounding
				      * faces (quads in 3d) to the
				      * vector @p{a}. Points located
				      * on the line but on vertices
				      * are not included.
				      *
				      * Needed by the
				      * @p{compute_support_points_laplace}
				      * function. For @p{dim=1} and 2
				      * this function is empty.
				      *
				      * This function is made virtual
				      * in order to allow derived
				      * classes to choose shape
				      * function support points
				      * differently than the present
				      * class, which chooses the
				      * points as interpolation points
				      * on the boundary.
				      */
    virtual void
    add_quad_support_points(const typename Triangulation<dim>::cell_iterator &cell,
			    typename std::vector<Point<dim> > &a) const;
    
  private:
    
				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual
    typename Mapping<dim>::InternalDataBase *
    get_data (const UpdateFlags,
	      const Quadrature<dim>& quadrature) const;

				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual
    typename Mapping<dim>::InternalDataBase *
    get_face_data (const UpdateFlags flags,
		   const Quadrature<dim-1>& quadrature) const;

				     /**
				      * Implementation of the interface in
				      * @ref{Mapping}.
				      */
    virtual
    typename Mapping<dim>::InternalDataBase *
    get_subface_data (const UpdateFlags flags,
		      const Quadrature<dim-1>& quadrature) const;
    
				     /**
				      * Compute shape values and/or
				      * derivatives.
				      */
    virtual void
    compute_shapes_virtual (const typename std::vector<Point<dim> > &unit_points,
			    typename MappingQ1<dim>::InternalData &data) const;

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
    void
    set_laplace_on_quad_vector(std::vector<std::vector<double> > &loqvs) const;
    
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
				      * @p{laplace_on_hex(quad)_vector}
				      * and applies it to the vector
				      * @p{a} to compute the inner
				      * support points as a linear
				      * combination of the exterior
				      * points.
				      *
				      * The vector @p{a} initially
				      * containts the locations of the
				      * @p{n_outer} points, the
				      * @p{n_inner} computed inner
				      * points are appended.
				      */
    void apply_laplace_vector(const std::vector<std::vector<double> > &lvs,
			      typename std::vector<Point<dim> > &a) const;
    
				     /**
				      * Computes the support points of
				      * the mapping.
				      */
    virtual void compute_mapping_support_points(
      const typename Triangulation<dim>::cell_iterator &cell,
      typename std::vector<Point<dim> > &a) const;

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
				      * values, in order to make the
				      * transformation as smooth as
				      * possible.
				      */
    void compute_support_points_laplace(
      const typename Triangulation<dim>::cell_iterator &cell,
      typename std::vector<Point<dim> > &a) const;
    
				     /**
				      * Simple version of the
				      * @p{compute_support_points_laplace}
				      * function. Does not use the
				      * solution to Laplace
				      * equation. Computes the inner
				      * support points by simple
				      * interpolations.
				      *
				      * This function isn't used in
				      * the code as it was replaced in
				      * @p{compute_mapping_support_points}
				      * by the
				      * @p{compute_support_points_laplace}
				      * function.
				      *
				      * Nethertheless this function is
				      * kept as someone might want to
				      * do some comparative tests.
				      */
    void compute_support_points_simple(
      const typename Triangulation<dim>::cell_iterator &cell,
      typename std::vector<Point<dim> > &a) const;    
    
				     /**
				      * For @p{dim=2} and 3. Simple
				      * version of the
				      * @p{add_face_support_points}
				      * function.
				      *
				      * Needed by the
				      * @p{compute_support_points_simple}
				      */
//TODO:[RH] (later) remove this function altogether?    
    void fill_quad_support_points_simple (const typename Triangulation<dim>::cell_iterator &cell,
					  typename std::vector<Point<dim> > &a) const;
    
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
    const unsigned int n_shape_functions;

				     /**
				      * Mapping from lexicographic to
				      * to the Qp shape function
				      * numbering. Its size is
				      * @p{dofs_per_cell}.
				      */
    std::vector<unsigned int> renumber;

				     /**
				      * If this flag is set @p{true}
				      * then @p{MappingQ} is used on
				      * all cells, not only on
				      * boundary cells.
				      *
				      * The default value is false.
				      *
				      * This flag is kept in the
				      * implementation to allow a fast
				      * switch between the two cases,
				      * as someone might want to do
				      * some comparative tests.
				      */
    static const bool use_mapping_q_on_all_cells = false;
};


/* -------------- declaration of explicit specializations ------------- */


template<> MappingQ<1>::MappingQ (const unsigned int);
template<> MappingQ<1>::~MappingQ ();
template<> void MappingQ<1>::compute_shapes_virtual (
  const std::vector<Point<1> > &unit_points,
  MappingQ1<1>::InternalData   &data) const;
template <> void MappingQ<1>::set_laplace_on_quad_vector(
  std::vector<std::vector<double> > &) const;
template <> void MappingQ<3>::set_laplace_on_hex_vector(
  std::vector<std::vector<double> > &lohvs) const;
template <> void MappingQ<1>::compute_laplace_vector(
  std::vector<std::vector<double> > &) const;
template <> void MappingQ<1>::add_line_support_points (
  const Triangulation<1>::cell_iterator &,
  std::vector<Point<1> > &) const;
template<> void MappingQ<3>::add_quad_support_points(
  const Triangulation<3>::cell_iterator &cell,
  std::vector<Point<3> >                &a) const;
template <> void MappingQ<3>::fill_quad_support_points_simple (
  const Triangulation<3>::cell_iterator &cell,
  std::vector<Point<3> > &a) const;



#endif
