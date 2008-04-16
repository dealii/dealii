//---------------------------------------------------------------------------
//    $Id: function_parser.h 14594 2007-03-22 20:17:41Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__fe_function_h
#define __deal2__fe_function_h

#include <base/function.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/mapping_q1.h>
#include <base/function.h>
#include <base/point.h>
#include <base/tensor.h>

#include <lac/vector.h>

DEAL_II_NAMESPACE_OPEN

namespace Functions
{

/**
 * This is an interpolation function for the given dof handler and
 * the given solution vector. The points at which this function can
 * be evaluated MUST be inside the domain of the dof handler, but
 * except from this, no other requirement is given. This function is
 * rather slow, as it needs to construct a quadrature object for the
 * point (or set of points) where you want to evaluate your finite
 * element function. In order to do so, it needs to find out where
 * the points lie.
 *
 * If you know in advance in which cell your points lye, you can
 * accelerate things a bit, by calling set_active_cell before
 * asking for values or gradients of the function. If you don't do
 * this, and your points don't lie in the cell that is currently
 * stored, the function GridTools::find_cell_around_point is called
 * to find out where the point is. You can specify an optional
 * mapping to use when looking for points in the grid. If you don't
 * do so, this function uses a Q1 mapping.
 *
 * Once the FEFieldFunction knows where the points lie, it creates a
 * quadrature formula for those points, and calls
 * FEValues::get_function_values or FEValues::get_function_grads with
 * the given quadrature points.
 *
 * If you only need the quadrature points but not the values of the
 * finite element function (you might want this for the adjoint
 * interpolation), you can also use the function @p
 * compute_point_locations alone.
 *
 * An example of how to use this function is the following:
 *    
 * \code
 * 
 * // Generate two triangulations
 * Triangulation<dim> tria_1;
 * Triangulation<dim> tria_2;
 *
 * // Read the triangulations from files, or build them up, or get
 * // them from some place...  Assume that tria_2 is *entirely*
 * // included in tria_1
 * ...
 *
 * // Associate a dofhandler and a solution to the first
 * // triangulation
 * DoFHandler<dim> dh1(tria_1);
 * Vector<double> solution_1;
 *
 * // Do the same with the second
 * DoFHandler<dim> dh2;
 * Vector<double> solution_2;
 *
 * // Setup the system, assemble matrices, solve problems and get the
 * // nobel prize on the first domain...
 * ...
 *
 * // Now project it to the second domain
 * FEFieldFunction<dim> fe_function_1 (dh_1, solution_1);
 * VectorTools::project(dh_2, constraints_2, quad, fe_function_1, solution_2);
 *
 * // Or interpolate it...
 * Vector<double> solution_3;
 * VectorTools::interpolate(dh_2, fe_function_1, solution_3);
 *
 * \endcode
 *
 * The snippet of code above will work assuming that the second
 * triangulation is entirely included in the first one. 
 *
 * FEFieldFunction is designed to be an easy way to get the results of
 * your computations across different, possibly non matching,
 * grids. No knowledge of the location of the points is assumed in
 * this class, which makes it rely entirely on the
 * GridTools::find_active_cell_around_point utility for its
 * job. However the class can be fed an "educated guess" of where the
 * points that will be computed actually are by using the
 * FEFieldFunction::set_active_cell method, so if you have a smart way to
 * tell where your points are, you will save a lot of computational
 * time by letting this class know.
 *
 * An optimization based on a caching mechanism was used by the
 * author of this class for the implementation of a Finite Element
 * Immersed Boundary Method.
 *
 *  \ingroup functions
 *
 *  \author Luca Heltai, 2006    
 *  
 *  \todo Add hp functionality
 */
  template <int dim, 
	    typename DH=DoFHandler<dim>,
	    typename VECTOR=Vector<double> >
  class FEFieldFunction :  public Function<dim> 
  {
    public:
				       /**
					* Construct a vector
					* function. A smart pointers
					* is stored to the dof
					* handler, so you have to make
					* sure that it make sense for
					* the entire lifetime of this
					* object. The number of
					* components of this functions
					* is equal to the number of
					* components of the finite
					* element object. If a mapping
					* is specified, that is what
					* is used to find out where
					* the points lay. Otherwise
					* the standard Q1 mapping is
					* used.
					*/
      FEFieldFunction (const DH           &dh,
		       const VECTOR       &data_vector, 
		       const Mapping<dim> &mapping = StaticMappingQ1<dim>::mapping);
  
				       /**
					* Set the current cell. If you
					* know in advance where your
					* points lie, you can tell
					* this object by calling this
					* function. This will speed
					* things up a little.
					*/
      void set_active_cell(typename DH::active_cell_iterator &newcell);
  
				       /**
					* Get ONE vector value at the
					* given point. It is
					* inefficient to use single
					* points. If you need more
					* than one at a time, use the
					* vector_value_list()
					* function. For efficiency
					* reasons, it is better if all
					* the points lie on the same
					* cell. This is not mandatory,
					* however it does speed things
					* up.
					*/ 
      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &values) const;

				       /**
					* Return the value of the
					* function at the given
					* point. Unless there is only
					* one component (i.e. the
					* function is scalar), you
					* should state the component
					* you want to have evaluated;
					* it defaults to zero,
					* i.e. the first component.
					* It is inefficient to use
					* single points. If you need
					* more than one at a time, use
					* the vector_value_list
					* function. For efficiency
					* reasons, it is better if all
					* the points lie on the same
					* cell. This is not mandatory,
					* however it does speed things
					* up.
					*/
      virtual double value (const Point< dim > &    p,
			    const unsigned int  component = 0)    const;
  
				       /**
					* Set @p values to the point
					* values of the specified
					* component of the function at
					* the @p points. It is assumed
					* that @p values already has
					* the right size, i.e. the
					* same size as the points
					* array. This is rather
					* efficient if all the points
					* lie on the same cell. If
					* this is not the case, things
					* may slow down a bit.
					*/
      virtual void value_list (const std::vector<Point< dim > > &    points,
			       std::vector< double > &values, 
			       const unsigned int  component = 0)    const;

  
				       /**
					* Set @p values to the point
					* values of the function at
					* the @p points. It is assumed
					* that @p values already has
					* the right size, i.e. the
					* same size as the points
					* array. This is rather
					* efficient if all the points
					* lie on the same cell. If
					* this is not the case, things
					* may slow down a bit.
					*/
      virtual void vector_value_list (const std::vector<Point< dim > > &    points,
				      std::vector< Vector<double> > &values) const;
  
				       /**
					* Return the gradient of all
					* components of the function
					* at the given point.  It is
					* inefficient to use single
					* points. If you need more
					* than one at a time, use the
					* vector_value_list
					* function. For efficiency
					* reasons, it is better if all
					* the points lie on the same
					* cell. This is not mandatory,
					* however it does speed things
					* up.
					*/
      virtual void 
      vector_gradient (const Point< dim > &p, 
		       std::vector< Tensor< 1, dim > > &gradients) const;
  
				       /**
					* Return the gradient of the
					* specified component of the
					* function at the given point.
					* It is inefficient to use
					* single points. If you need
					* more than one at a time, use
					* the vector_value_list
					* function. For efficiency
					* reasons, it is better if all
					* the points lie on the same
					* cell. This is not mandatory,
					* however it does speed things
					* up.
					*/
      virtual Tensor<1,dim> gradient(const Point< dim > &p, 
				     const unsigned int component = 0)const;
  
				       /**
					* Return the gradient of all
					* components of the function
					* at all the given points.
					* This is rather efficient if
					* all the points lie on the
					* same cell. If this is not
					* the case, things may slow
					* down a bit.
					*/
      virtual void 
      vector_gradient_list (const std::vector< Point< dim > > &p, 
			    std::vector< 
			    std::vector< Tensor< 1, dim > > > &gradients) const;

				       /**
					* Return the gradient of the
					* specified component of the
					* function at all the given
					* points.  This is rather
					* efficient if all the points
					* lie on the same cell. If
					* this is not the case, things
					* may slow down a bit.
					*/
      virtual void 
      gradient_list (const std::vector< Point< dim > > &p, 
		     std::vector< Tensor< 1, dim > > &gradients, 
		     const unsigned int component=0) const;
  
				       /**
					* Create quadrature
					* rules. This function groups
					* the points into blocks that
					* live in the same cell, and
					* fills up three vectors: @p
					* cells, @p qpoints and @p
					* maps. The first is a list of
					* the cells that contain the
					* points, the second is a list
					* of quadrature points
					* matching each cell of the
					* first list, and the third
					* contains the index of the
					* given quadrature points,
					* i.e., @p points[maps[3][4]]
					* ends up as the 5th
					* quadrature point in the 4th
					* cell. This is where
					* optimization would
					* help. This function returns
					* the number of cells that
					* contain the given set of
					* points.
					*/
      unsigned int
      compute_point_locations(const std::vector<Point<dim> > &points,
			      std::vector<typename DH::active_cell_iterator > &cells,
			      std::vector<std::vector<Point<dim> > > &qpoints,
			      std::vector<std::vector<unsigned int> > &maps) const;
    
    private:
				       /**
					* Pointer to the dof handler.
					*/
      SmartPointer<const DH> dh;
    
				       /**
					* A reference to the actual
					* data vector.
					*/
      const VECTOR & data_vector;
    
				       /**
					* A reference to the mapping
					* being used.
					*/
      const Mapping<dim> & mapping;
  
				       /**
					* The current cell in which we
					* are evaluating.
					*/
      mutable typename DH::active_cell_iterator cell;
  
				       /**
					* Store the number of
					* components of this function.
					*/
      const unsigned int n_components;
  };
}

DEAL_II_NAMESPACE_CLOSE

#endif
