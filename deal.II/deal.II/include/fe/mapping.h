//----------------------------  mapping.h  ---------------------------
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
//----------------------------  mapping.h  ---------------------------
#ifndef __deal2__mapping_h
#define __deal2__mapping_h


#include <cmath>
#include <base/point.h>
#include <base/subscriptor.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <fe/fe_update_flags.h>

template <int dim> class Quadrature;
template <int dim> class FEValuesData;
template <int dim> class FEValues;
template <int dim> class FEFaceValues;
template <int dim> class FESubfaceValues;

/**
 * Abstract basis class for mapping classes.
 *
 * The interface for filling the tables of @ref{FEValues} is provided.
 * Everything else has to happen in derived classes.
 *
 * The following paragraph applies to the implementation of
 * @ref{FEValues}. Usage of the class is as follows: first, call the
 * functionss @p{update_once} and @p{update_each} with the update
 * flags you need. This includes the flags needed by the
 * @ref{FiniteElement}. Then call @p{get_*_data} and with the or'd
 * results.  This will initialize and return some internal data
 * structures.  On the first cell, call @p{fill_fe_*_values} with the
 * result of @p{update_once}. Finally, on each cell, use
 * @p{fill_fe_*_values} with the result of @p{update_each} to compute
 * values for a special cell.
 *
 * A hint to implementators: no function except the two functions
 * @p{update_once} and @p{update_each} may add any flags.
 *
 * @author Guido Kanschat, Ralf Hartmann 2000, 2001
 */
template <int dim>
class Mapping : public Subscriptor
{
  public:
    
				     /**
				      * Virtual destructor.
				      */
    virtual ~Mapping ();
    
				     /**
				      * Transforms the point @p{p} on
				      * the unit cell to the point
				      * @p{p_real} on the real cell
				      * @p{cell} and returns @p{p_real}.
				      */
    virtual Point<dim> transform_unit_to_real_cell (
      const typename Triangulation<dim>::cell_iterator cell,
      const Point<dim> &p) const=0;
    
				     /**
				      * Transforms the point @p{p} on
				      * the real cell to the point
				      * @p{p_unit} on the unit cell
				      * @p{cell} and returns @p{p_unit}.
				      */
    virtual Point<dim> transform_real_to_unit_cell (
      const typename Triangulation<dim>::cell_iterator cell,
      const Point<dim> &p) const=0;
    
				     /**
				      * Class for internal data of finite
				      * element and mapping objects.
				      */
    class InternalDataBase: public Subscriptor
    {
      private:
					 /**
					  * Copy constructor forbidden.
					  */
        InternalDataBase (const InternalDataBase&);

      public:
					 /**
					  * Constructor. Sets
					  * @p{UpdateFlags} to
					  * @p{update_default} and
					  * @p{first_cell} to @p{true}.
					  */
        InternalDataBase ();
	
					 /**
					  * Virtual destructor for
					  * derived classes
					  */
	virtual ~InternalDataBase ();

	                                 /**
					  * Values updated by the constructor or
					  * by reinit.
					  */
	UpdateFlags          update_flags;
	
					 /**
					  * Values computed by
					  * constructor.
					  */
	UpdateFlags          update_once;

					 /**
					  * Values updated on each
					  * cell by reinit.
					  */
	UpdateFlags          update_each;

					 /**
					  * If @p{first_cell==true}
					  * this function returns
					  * @p{update_flags},
					  * i.e. @p{update_once|update_each}.
					  * If @p{first_cell==false}
					  * it returns
					  * @p{update_each}.
					  */
	UpdateFlags  current_update_flags() const;

					 /**
					  * Determine if this is the first
					  * cell visited. Then, we need to
					  * use @p{update_once} to fill
					  * cell-independent fields.
					  */
	bool first_cell;

					 /**
					  * Return an estimate (in
					  * bytes) or the memory
					  * consumption of this
					  * object.
					  */
	virtual unsigned int memory_consumption () const;
    };
    
				     /**
				      * Tranform a field of covariant
				      * vectors.  There must be one
				      * vector for each quadrature
				      * point. Alternatively, for
				      * faces and subfaces, the number
				      * of the first quadrature point
				      * can be given as additional
				      * argument.
				      */
    virtual void transform_covariant (std::vector<Tensor<1,dim> >       &dst,
				      const std::vector<Tensor<1,dim> > &src,
				      const InternalDataBase& internal,
				      const unsigned int src_offset) const = 0;
    
				     /**
				      * Tranform a field of
				      * contravariant vectors.  There
				      * must be one vector for each
				      * quadrature
				      * point. Alternatively, for
				      * faces and subfaces, the number
				      * of the first quadrature point
				      * can be given as additional
				      * argument.
				      */
    virtual void transform_contravariant (std::vector<Tensor<1,dim> >       &dst,
					  const std::vector<Tensor<1,dim> > &src,
					  const InternalDataBase& internal,
					  const unsigned int src_offset) const = 0;
    
				     /**
				      * Tranform a field of covariant vectors.
				      * There must be one vector for each quadrature
				      * point. Alternatively, for faces and subfaces,
				      * the first quadrature point can be
				      * given as additional argument.
				      */
    virtual void transform_covariant (std::vector<Point<dim> >       &dst,
				      const std::vector<Point<dim> > &src,
				      const InternalDataBase& internal,
				      const unsigned int src_offset) const = 0;
    
				     /**
				      * Tranform a field of contravariant vectors.
				      * There must be one vector for each quadrature
				      * point. Alternatively, for faces and subfaces,
				      * the first quadrature point can be
				      * given as additional argument.
				      */
    virtual void transform_contravariant (std::vector<Point<dim> >       &dst,
					  const std::vector<Point<dim> > &src,
					  const InternalDataBase& internal,
					  const unsigned int src_offset) const = 0;

				     /**
				      * Indicate fields to be updated in the
				      * constructor of @ref{FEValues}. Especially,
				      * fields not asked for by @ref{FEValues}, but
				      * computed for efficiency reasons will be
				      * notified here.
				      */
    virtual UpdateFlags update_once (const UpdateFlags) const = 0;
    
				     /**
				      * The same as @p{update_once},
				      * but for the flags to be updated for
				      * each grid cell.
				      */
    virtual UpdateFlags update_each (const UpdateFlags) const = 0;
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidData);


    
  protected:

				     /**
				      * Vector of unit normal
				      * directions. The entry divided by
				      * 2 determines the non-zero
				      * component of the normal vector:
				      * 0 means x, 1 means y and 2 means
				      * z. The entry modulo 2 determines
				      * the orientation of the first
				      * tangential vector in the
				      * cross-product. This has to be
				      * chosen such that the normal
				      * vector points outwards.
				      *
				      * This variable is purely for
				      * internal use and its values are
				      * determined by its usage in the
				      * source code.
				      */
    static const unsigned int normal_directions[2*dim];

  private:
    
				     /**
				      * Prepare internal data
				      * structures and fill in values
				      * independent of the cell.
				      */
    virtual InternalDataBase*
    get_data (const UpdateFlags,
	      const Quadrature<dim>& quadrature) const = 0;

				     /**
				      * Prepare internal data
				      * structure for transformation
				      * of faces and fill in values
				      * independent of the cell.
				      */
    virtual InternalDataBase*
    get_face_data (const UpdateFlags flags,
		   const Quadrature<dim-1>& quadrature) const = 0;
    
				     /**
				      * Prepare internal data
				      * structure for transformation
				      * of children of faces and fill
				      * in values independent of the
				      * cell.
				      */
    virtual InternalDataBase*
    get_subface_data (const UpdateFlags flags,
		      const Quadrature<dim-1>& quadrature) const = 0;


				     /**
				      * Fill the transformation fields
				      * of @p{FEValues}.  Given a grid
				      * cell and the quadrature points
				      * on the unit cell, it computes
				      * all values specified by
				      * @p{flags}. The arrays to be
				      * filled have to have the
				      * correct size.
				      *
				      * Values are split into three
				      * groups: first,
				      * @p{quadrature_points} and
				      * @p{JxW_values} are
				      * filled with the quadrature
				      * rule transformed to the
				      * cell in physical space.
				      *
				      * The second group contains the
				      * matrices needed to transform
				      * vector-valued functions,
				      * namely
				      * @p{covariant_transformation},
				      * @p{contravariant_transformation} and the 
				      * derivatives
				      * @p{covariant_grads}.
				      *
				      */
    virtual void
    fill_fe_values (const typename DoFHandler<dim>::cell_iterator &cell,
		    const Quadrature<dim>& quadrature,
		    InternalDataBase& internal,
		    std::vector<Point<dim> >        &quadrature_points,
		    std::vector<double>             &JxW_values) const = 0;

				     /**
				      * Performs the same as @p{fill_fe_values}
				      * on a face.
				      * Additionally, @p{boundary_form} and
				      * @p{normal_vectors} can be
				      * computed on surfaces. The
				      * boundary form is the vector
				      * product of the image of
				      * coordinate vectors on the
				      * surface of the unit
				      * cell. It is a
				      * vector normal to the surface,
				      * pointing outwards and having
				      * the length of the surface
				      * element.
				      * Therefore, it is more economic
				      * to use the boundary form
				      * instead of the product of the
				      * unit normal and the
				      * transformed quadrature weight.
				      */
    virtual void
    fill_fe_face_values (const typename DoFHandler<dim>::cell_iterator &cell,
			 const unsigned int face_no,
			 const Quadrature<dim-1>& quadrature,
			 InternalDataBase& internal,
			 std::vector<Point<dim> >        &quadrature_points,
			 std::vector<double>             &JxW_values,
			 std::vector<Tensor<1,dim> >        &boundary_form,
			 std::vector<Point<dim> >        &normal_vectors) const = 0;

				     /**
				      * See above.
				      */
    virtual void
    fill_fe_subface_values (const typename DoFHandler<dim>::cell_iterator &cell,
			    const unsigned int face_no,
			    const unsigned int sub_no,
			    const Quadrature<dim-1>& quadrature,
			    InternalDataBase& internal,
			    std::vector<Point<dim> >        &quadrature_points,
			    std::vector<double>             &JxW_values,
			    std::vector<Tensor<1,dim> >        &boundary_form,
			    std::vector<Point<dim> >        &normal_vectors) const = 0;

				     /**
				      * Give class @p{FEValues} access
				      * to the private @p{get_...data}
				      * and @p{fill_fe_...values}
				      * functions.
				      */
  friend class FEValues<dim>;
  friend class FEFaceValues<dim>;
  friend class FESubfaceValues<dim>;
};


/* -------------- declaration of explicit specializations ------------- */


template<> const unsigned int Mapping<1>::normal_directions[2];
template<> const unsigned int Mapping<2>::normal_directions[4];
template<> const unsigned int Mapping<3>::normal_directions[6];


/* ------------------------- inline functions ------------------------- */


template <int dim>
inline
UpdateFlags
Mapping<dim>::InternalDataBase::current_update_flags() const
{
  if (first_cell)
    {
      Assert(update_flags==(update_once|update_each), ExcInternalError());
      return update_flags;
    }
  else
    return update_each;
  
  return update_default;
}



#endif
