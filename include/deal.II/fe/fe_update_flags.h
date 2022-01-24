// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_fe_update_flags_h
#define dealii_fe_update_flags_h


#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/point.h>
#include <deal.II/base/table.h>
#include <deal.II/base/tensor.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <int, int>
class FiniteElement;
class UpdateFlags;

namespace internal
{
  constexpr UpdateFlags make_update_flags(int);
}
#endif

/*!@addtogroup feaccess */
/*@{*/

/**
 * The enum type given to the constructors of FEValues, FEFaceValues and
 * FESubfaceValues, telling those objects which data will be needed on each
 * mesh cell.
 *
 * Selecting these flags in a restrictive way is crucial for the efficiency of
 * FEValues::reinit(), FEFaceValues::reinit() and FESubfaceValues::reinit().
 * Therefore, only the flags actually needed should be selected. It is the
 * responsibility of the involved Mapping and FiniteElement to add additional
 * flags according to their own requirements. For instance, most finite
 * elements will add #update_covariant_transformation if #update_gradients is
 * selected.  By default, all flags are off, i.e. no reinitialization will be
 * done.
 *
 * You can select more than one flag by concatenation using the bitwise or
 * operator|(UpdateFlags,UpdateFlags).
 *
 * <h3>Use of these flags flags</h3>
 *
 * More information on the use of this type both in user code as well as
 * internally can be found in the documentation modules on
 * @ref UpdateFlags "The interplay of UpdateFlags, Mapping, and FiniteElement in FEValues"
 * and
 * @ref FE_vs_Mapping_vs_FEValues "How Mapping, FiniteElement, and FEValues work together".
 */
class UpdateFlags
{
	public:

/**
 * Return an object in which all bits are set which are set in the current
 * object or in the input argument @p in.
 */
constexpr UpdateFlags
operator|(const UpdateFlags in) const;

/**
 * Set the bits from the input argument @p in also in the current object.
 */
constexpr UpdateFlags &
operator|=(const UpdateFlags in);

/**
 * Return an object in which all bits are set which are set in the current
 * object as well as in the input argument @p in.
 */
constexpr UpdateFlags
operator&(const UpdateFlags in) const;

/**
 * Clear all the bits in the current object if they are not also set in the
 * input argument @p in.
 */
constexpr UpdateFlags &
operator&=(const UpdateFlags in);

/**
 * Check if the current object and the input argument @p in are equal.
 */ 
    constexpr bool operator==(const UpdateFlags in) const;

/**
 * Check if the current object and the input argument @p in are not equal.
 */
    constexpr bool operator!=(const UpdateFlags in) const;

/**
 *  Check if the current objects contains the flags given in the input argument @p in.
 */
    constexpr bool contains(const UpdateFlags in) const;

  //! No update
  static const UpdateFlags update_default;
  //! Shape function values
  /**
   * Compute the values of the shape functions at the quadrature points on the
   * real space cell. For the usual Lagrange elements, these values are equal
   * to the values of the shape functions at the quadrature points on the unit
   * cell, but they are different for more complicated elements, such as
   * FE_RaviartThomas elements.
   */
  static const UpdateFlags update_values;
  //! Shape function gradients
  /**
   * Compute the gradients of the shape functions in coordinates of the real
   * cell.
   */
  static const UpdateFlags update_gradients;
  //! Second derivatives of shape functions
  /**
   * Compute the second derivatives of the shape functions in coordinates of
   * the real cell.
   */
  static const UpdateFlags update_hessians;
  //! Third derivatives of shape functions
  /**
   * Compute the third derivatives of the shape functions in coordinates of
   * the real cell
   */
  static const UpdateFlags update_3rd_derivatives;
  //! Outer normal vector, not normalized
  /**
   * Vector product of tangential vectors, yielding a normal vector with a
   * length corresponding to the surface element; may be more efficient than
   * computing both.
   */
  static const UpdateFlags update_boundary_forms;
  //! Transformed quadrature points
  /**
   * Compute the quadrature points location in real cell coordinates.
   *
   * FEValues objects take the quadrature point locations on the
   * reference cell as an argument of the constructor (via the
   * Quadrature object). For most finite elements, knowing the
   * location of quadrature points on the reference cell is all that
   * is necessary to evaluate shape functions, evaluate the mapping,
   * and other things. On the other hand, if you want to evaluate a
   * right hand side function $f(\mathbf x_q)$ at quadrature point
   * locations $\mathbf x_q$ on the real cell, you need to pass this
   * flag to the FEValues constructor to make sure you can later
   * access them.
   *
   * In the context of DataPostprocessor,
   * DataPostprocessorInputs::CommonInputs::evaluation_points will be updated.
   */
  static const UpdateFlags update_quadrature_points;
  //! Transformed quadrature weights
  /**
   * Compute the quadrature weights on the real cell, i.e. the weights of the
   * quadrature rule multiplied with the determinant of the Jacobian of the
   * transformation from reference to real cell.
   */
  static const UpdateFlags update_JxW_values;
  //! Normal vectors
  /**
   * Compute the normal vectors, either for a face or for a cell of
   * codimension one. Setting this flag for any other object will raise an
   * error.
   */
  static const UpdateFlags  update_normal_vectors;
  //! Volume element
  /**
   * Compute the Jacobian of the transformation from the reference cell to the
   * real cell.
   */
  static const UpdateFlags  update_jacobians;
  //! Gradient of volume element
  /**
   * Compute the derivatives of the Jacobian of the transformation.
   */
  static const UpdateFlags update_jacobian_grads;
  //! Volume element
  /**
   * Compute the inverse Jacobian of the transformation from the reference
   * cell to the real cell.
   */
  static const UpdateFlags update_inverse_jacobians;
  //! Covariant transformation
  /**
   * Compute all values the Mapping needs to perform a contravariant
   * transformation of vectors. For special mappings like MappingCartesian
   * this may be simpler than #update_inverse_jacobians.
   */
  static const UpdateFlags update_covariant_transformation;
  //! Contravariant transformation
  /**
   * Compute all values the Mapping needs to perform a contravariant
   * transformation of vectors. For special mappings like MappingCartesian
   * this may be simpler than #update_jacobians.
   */
  static const UpdateFlags update_contravariant_transformation;
  //! Shape function values of transformation
  /**
   * Compute the shape function values of the transformation defined by the
   * Mapping.
   */
  static const UpdateFlags update_transformation_values;
  //! Shape function gradients of transformation
  /**
   * Compute the shape function gradients of the transformation defined by the
   * Mapping.
   */
  static const UpdateFlags update_transformation_gradients;
  //! Determinant of the Jacobian
  /**
   * Compute the volume element in each quadrature point.
   */
  static const UpdateFlags update_volume_elements;
  /**
   * Compute the derivatives of the Jacobian of the transformation pushed
   * forward to the real cell coordinates.
   */
  static const UpdateFlags update_jacobian_pushed_forward_grads;
  /**
   * Compute the second derivatives of the Jacobian of the transformation.
   */
  static const UpdateFlags update_jacobian_2nd_derivatives;
  /**
   * Compute the second derivatives of the Jacobian of the transformation
   * pushed forward to the real cell coordinates.
   */
  static const UpdateFlags update_jacobian_pushed_forward_2nd_derivatives;
  /**
   * Compute the third derivatives of the Jacobian of the transformation.
   */
  static const UpdateFlags update_jacobian_3rd_derivatives;
  /**
   * Compute the third derivatives of the Jacobian of the transformation
   * pushed forward to the real cell coordinates.
   */
  static const UpdateFlags update_jacobian_pushed_forward_3rd_derivatives;
  //! Values needed for Piola transform
  /**
   * Combination of the flags needed for Piola transform of Hdiv elements.
   */
  static const UpdateFlags update_piola;
  /**
   * Combination of the flags that require a mapping calculation
   */
  static const UpdateFlags update_mapping;

	private:
  /**
   * An integer representing the flags.
   */
  int field;

  /**
   * Private constructor only to be used internally.
   */
  constexpr UpdateFlags(int i);

  /**
   * Friend declaration so that make_enum can call the private constructor.
   */
  friend constexpr UpdateFlags make_enum(int);
};


/**
 * Output operator which outputs update flags as a set of or'd text values.
 *
 * @ref UpdateFlags
 */
template <class StreamType>
inline StreamType &
operator<<(StreamType &s, const UpdateFlags u)
{
  s << " UpdateFlags|";
  if (u & UpdateFlags::update_values)
    s << "values|";
  if (u & UpdateFlags::update_gradients)
    s << "gradients|";
  if (u & UpdateFlags::update_hessians)
    s << "hessians|";
  if (u & UpdateFlags::update_3rd_derivatives)
    s << "3rd_derivatives|";
  if (u & UpdateFlags::update_quadrature_points)
    s << "quadrature_points|";
  if (u & UpdateFlags::update_JxW_values)
    s << "JxW_values|";
  if (u & UpdateFlags::update_normal_vectors)
    s << "normal_vectors|";
  if (u & UpdateFlags::update_jacobians)
    s << "jacobians|";
  if (u & UpdateFlags::update_inverse_jacobians)
    s << "inverse_jacobians|";
  if (u & UpdateFlags::update_jacobian_grads)
    s << "jacobian_grads|";
  if (u & UpdateFlags::update_covariant_transformation)
    s << "covariant_transformation|";
  if (u & UpdateFlags::update_contravariant_transformation)
    s << "contravariant_transformation|";
  if (u & UpdateFlags::update_transformation_values)
    s << "transformation_values|";
  if (u & UpdateFlags::update_transformation_gradients)
    s << "transformation_gradients|";
  if (u & UpdateFlags::update_jacobian_pushed_forward_grads)
    s << "jacobian_pushed_forward_grads|";
  if (u & UpdateFlags::update_jacobian_2nd_derivatives)
    s << "jacobian_2nd_derivatives|";
  if (u & UpdateFlags::update_jacobian_pushed_forward_2nd_derivatives)
    s << "jacobian_pushed_forward_2nd_derivatives|";
  if (u & UpdateFlags::update_jacobian_3rd_derivatives)
    s << "jacobian_3rd_derivatives|";
  if (u & UpdateFlags::update_jacobian_pushed_forward_3rd_derivatives)
    s << "jacobian_pushed_forward_3rd_derivatives|";

  // TODO: check that 'u' really only has the flags set that are handled above
  return s;
}

// Import all flags into the global namespace.
constexpr UpdateFlags update_default = internal::make_enum(0);
constexpr UpdateFlags UpdateFlags::update_values = internal::make_enum(0x0001);
constexpr UpdateFlags UpdateFlags::update_gradients = internal::make_enum(0x0002);
constexpr UpdateFlags UpdateFlags::update_hessians = internal::make_enum(0x0004);
constexpr UpdateFlags UpdateFlags::update_3rd_derivatives = internal::make_enum(0x0008);
constexpr UpdateFlags UpdateFlags::update_boundary_forms = internal::make_enum(0x0010);
constexpr UpdateFlags UpdateFlags::update_quadrature_points = internal::make_enum(0x0020);
constexpr UpdateFlags UpdateFlags::update_JxW_values = internal::make_enum(0x0040);
constexpr UpdateFlags UpdateFlags::update_normal_vectors = internal::make_enum(0x0080);
constexpr UpdateFlags UpdateFlags::update_jacobians = internal::make_enum(0x01000);
constexpr UpdateFlags UpdateFlags::update_jacobian_grads = internal::make_enum(0x02000:
constexpr UpdateFlags UpdateFlags::update_inverse_jacobians = internal::make_enum(0x0400);
constexpr UpdateFlags UpdateFlags::update_covariant_transformation = internal::make_enum(0x0800);
constexpr UpdateFlags UpdateFlags::update_contravariant_transformation = internal::make_enum(0x1000);
constexpr UpdateFlags UpdateFlags::update_transformation_values = internal::make_enum(0x2000);
constexpr UpdateFlags UpdateFlags::update_transformation_gradients = internal::make_enum(0x4000);
constexpr UpdateFlags UpdateFlags::update_volume_elements = internal::make_enum(0x10000);
constexpr UpdateFlags UpdateFlags::update_jacobian_pushed_forward_grads = internal::make_enum(0x100000);
constexpr UpdateFlags UpdateFlags::update_jacobian_2nd_derivatives = internal::make_enum(0x200000);
constexpr UpdateFlags UpdateFlags::update_jacobian_pushed_forward_2nd_derivatives = internal::make_enum(0x400000);
constexpr UpdateFlags UpdateFlags::update_jacobian_3rd_derivatives = internal::make_enum(0x800000);
constexpr UpdateFlags UpdateFlags::update_jacobian_pushed_forward_3rd_derivatives = internal::make_enum(0x1000000);
constexpr UpdateFlags UpdateFlags::update_piola = update_volume_elements | update_contravariant_transformation,
constexpr UpdateFlags UpdateFlags::update_mapping =
update_quadrature_points | update_JxW_values | update_jacobians |
update_jacobian_grads | update_jacobian_pushed_forward_grads |
update_jacobian_2nd_derivatives |
update_jacobian_pushed_forward_2nd_derivatives |
update_jacobian_3rd_derivatives |
update_jacobian_pushed_forward_3rd_derivatives | update_inverse_jacobians |
update_boundary_forms | update_normal_vectors |
update_covariant_transformation | update_contravariant_transformation |
update_transformation_values | update_transformation_gradients |
update_volume_elements;


// Class member definitions

inline constexpr UpdateFlags
UpdateFlags::operator|(const UpdateFlags in) const
{
  return internal::make_enum(field | in.field);
}



inline constexpr UpdateFlags &
UpdateFlags::operator|=(const UpdateFlags in)
{
  field = field | in.field;
  return *this;
}



inline constexpr UpdateFlags
UpdateFlags::operator&(const UpdateFlags in) const
{
  return internal::make_enum(field & in.field);
}



inline constexpr UpdateFlags &
UpdateFlags::operator&=(const UpdateFlags in)
{
  field = field & in.field;
  return *this;
}



inline constexpr bool
UpdateFlags::operator==(const UpdateFlags in) const
{
  return field == in.field;
}



inline constexpr bool
UpdateFlags::operator!=(const UpdateFlags in) const
{
  return field != in.field;
}


inline constexpr bool
UpdateFlags::contains(const UpdateFlags in)
{
  return (field & in.field)!=0;
}

constexpr UpdateFlags UpdateFlags::update_default = update_default;
constexpr UpdateFlags UpdateFlags::update_values = update_values;
constexpr UpdateFlags UpdateFlags::update_gradients = update_gradients;
constexpr UpdateFlags UpdateFlags::update_hessians = update_hessians;
constexpr UpdateFlags UpdateFlags::update_3rd_derivatives = update_3rd_derivatives;
constexpr UpdateFlags UpdateFlags::update_boundary_forms = update_boundary_forms;
constexpr UpdateFlags UpdateFlags::update_quadrature_points = iupdate_quadrature_points;
constexpr UpdateFlags UpdateFlags::update_JxW_values = update_JxW_values;
constexpr UpdateFlags UpdateFlags::update_normal_vectors = update_normal_vectors;
constexpr UpdateFlags UpdateFlags::update_jacobians = update_jacobians;
constexpr UpdateFlags UpdateFlags::update_jacobian_grads = update_jacobian_grads;
constexpr UpdateFlags UpdateFlags::update_inverse_jacobians = update_inverses_jacobians;
constexpr UpdateFlags UpdateFlags::update_covariant_transformation = update_covariant_transformation;
constexpr UpdateFlags UpdateFlags::update_contravariant_transformation = update_contravariant_transformation;
constexpr UpdateFlags UpdateFlags::update_transformation_values = update_transformation_values;
constexpr UpdateFlags UpdateFlags::update_transformation_gradients = update_transformation_gradients;
constexpr UpdateFlags UpdateFlags::update_volume_elements = update_volume_elements;
constexpr UpdateFlags UpdateFlags::update_jacobian_pushed_forward_grads = update_jacobian_pushed_forward_grad;,
constexpr UpdateFlags UpdateFlags::update_jacobian_2nd_derivatives = update_jacobian_2nd_derivatives;
constexpr UpdateFlags UpdateFlags::update_jacobian_pushed_forward_2nd_derivatives = update_jacobian_pushed_forward_2nd_derivatives;
constexpr UpdateFlags UpdateFlags::update_jacobian_3rd_derivatives = update_jacobian_3rd_derivatives;
constexpr UpdateFlags UpdateFlags::update_jacobian_pushed_forward_3rd_derivatives = update_jacobian_pushed_forward_3rd_derivatives;
constexpr UpdateFlags UpdateFlags::update_piola = update_piola;
constexpr UpdateFlags UpdateFlags::update_mapping = update_mapping;

inline constexpr UpdateFlags::UpdateFlags(int i) : field(in)
{}

/**
 * This enum definition is used for storing similarities of the current cell
 * to the previously visited cell. This information is used for reusing data
 * when calling the method FEValues::reinit() (like derivatives, which do not
 * change if one cell is just a translation of the previous). Currently, this
 * variable does only recognize a translation and an inverted translation (if
 * dim<spacedim). However, this concept makes it easy to add additional states
 * to be detected in FEValues/FEFaceValues for making use of these
 * similarities as well.
 */
namespace CellSimilarity
{
  enum Similarity
  {
    /**
     * The cells differ by something besides a translation or inverted
     * translations.
     */
    none,
    /**
     * The cells differ by a translation.
     */
    translation,
    /**
     * The cells differ by an inverted translation.
     */
    inverted_translation,
    /**
     * The next cell is not valid.
     */
    invalid_next_cell
  };
}


namespace internal
{
  namespace FEValuesImplementation
  {
    /**
     * A class that stores all of the mapping related data used in
     * dealii::FEValues, dealii::FEFaceValues, and dealii::FESubfaceValues
     * objects. Objects of this kind will be given as <i>output</i> argument
     * when dealii::FEValues::reinit() calls Mapping::fill_fe_values() for a
     * given cell, face, or subface.
     *
     * The data herein will then be provided as <i>input</i> argument in the
     * following call to FiniteElement::fill_fe_values().
     *
     * @ingroup feaccess
     */
    template <int dim, int spacedim = dim>
    class MappingRelatedData
    {
    public:
      /**
       * Initialize all vectors to correct size.
       */
      void
      initialize(const unsigned int n_quadrature_points,
                 const UpdateFlags  flags);

      /**
       * Compute and return an estimate for the memory consumption (in bytes)
       * of this object.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Store an array of weights times the Jacobi determinant at the
       * quadrature points. This function is reset each time reinit() is
       * called. The Jacobi determinant is actually the reciprocal value of
       * the Jacobi matrices stored in this class, see the general
       * documentation of this class for more information.
       *
       * However, if this object refers to an FEFaceValues or FESubfaceValues
       * object, then the JxW_values correspond to the Jacobian of the
       * transformation of the face, not the cell, i.e. the dimensionality is
       * that of a surface measure, not of a volume measure. In this case, it
       * is computed from the boundary forms, rather than the Jacobian matrix.
       */
      std::vector<double> JxW_values;

      /**
       * Array of the Jacobian matrices at the quadrature points.
       */
      std::vector<DerivativeForm<1, dim, spacedim>> jacobians;

      /**
       * Array of the derivatives of the Jacobian matrices at the quadrature
       * points.
       */
      std::vector<DerivativeForm<2, dim, spacedim>> jacobian_grads;

      /**
       * Array of the inverse Jacobian matrices at the quadrature points.
       */
      std::vector<DerivativeForm<1, spacedim, dim>> inverse_jacobians;

      /**
       * Array of the derivatives of the Jacobian matrices at the quadrature
       * points, pushed forward to the real cell coordinates.
       */
      std::vector<Tensor<3, spacedim>> jacobian_pushed_forward_grads;

      /**
       * Array of the second derivatives of the Jacobian matrices at the
       * quadrature points.
       */
      std::vector<DerivativeForm<3, dim, spacedim>> jacobian_2nd_derivatives;

      /**
       * Array of the  second derivatives of the Jacobian matrices at the
       * quadrature points, pushed forward to the real cell coordinates.
       */
      std::vector<Tensor<4, spacedim>> jacobian_pushed_forward_2nd_derivatives;

      /**
       * Array of the  third derivatives of the Jacobian matrices at the
       * quadrature points.
       */
      std::vector<DerivativeForm<4, dim, spacedim>> jacobian_3rd_derivatives;

      /**
       * Array of the  third derivatives of the Jacobian matrices at the
       * quadrature points, pushed forward to the real cell coordinates.
       */
      std::vector<Tensor<5, spacedim>> jacobian_pushed_forward_3rd_derivatives;

      /**
       * Array of quadrature points. This array is set up upon calling
       * reinit() and contains the quadrature points on the real element,
       * rather than on the reference element.
       */
      std::vector<Point<spacedim>> quadrature_points;

      /**
       * List of outward normal vectors at the quadrature points.
       */
      std::vector<Tensor<1, spacedim>> normal_vectors;

      /**
       * List of boundary forms at the quadrature points.
       */
      std::vector<Tensor<1, spacedim>> boundary_forms;
    };


    /**
     * A class that stores all of the shape function related data used in
     * dealii::FEValues, dealii::FEFaceValues, and dealii::FESubfaceValues
     * objects. Objects of this kind will be given as <i>output</i> argument
     * when dealii::FEValues::reinit() calls FiniteElement::fill_fe_values().
     *
     * @ingroup feaccess
     */
    template <int dim, int spacedim = dim>
    class FiniteElementRelatedData
    {
    public:
      /**
       * Initialize all vectors to correct size.
       */
      void
      initialize(const unsigned int                  n_quadrature_points,
                 const FiniteElement<dim, spacedim> &fe,
                 const UpdateFlags                   flags);

      /**
       * Compute and return an estimate for the memory consumption (in bytes)
       * of this object.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Storage type for shape values. Each row in the matrix denotes the
       * values of a single shape function at the different points, columns
       * are for a single point with the different shape functions.
       *
       * If a shape function has more than one non-zero component (in deal.II
       * diction: it is non-primitive), then we allocate one row per non-zero
       * component, and shift subsequent rows backward.  Lookup of the correct
       * row for a shape function is thus simple in case the entire finite
       * element is primitive (i.e. all shape functions are primitive), since
       * then the shape function number equals the row number. Otherwise, use
       * the #shape_function_to_row_table array to get at the first row that
       * belongs to this particular shape function, and navigate among all the
       * rows for this shape function using the
       * FiniteElement::get_nonzero_components() function which tells us which
       * components are non-zero and thus have a row in the array presently
       * under discussion.
       */
      using ShapeVector = dealii::Table<2, double>;

      /**
       * Storage type for gradients. The layout of data is the same as for the
       * #ShapeVector data type.
       */
      using GradientVector = dealii::Table<2, Tensor<1, spacedim>>;

      /**
       * Likewise for second order derivatives.
       */
      using HessianVector = dealii::Table<2, Tensor<2, spacedim>>;

      /**
       * And the same also applies to the third order derivatives.
       */
      using ThirdDerivativeVector = dealii::Table<2, Tensor<3, spacedim>>;

      /**
       * Store the values of the shape functions at the quadrature points. See
       * the description of the data type for the layout of the data in this
       * field.
       */
      ShapeVector shape_values;

      /**
       * Store the gradients of the shape functions at the quadrature points.
       * See the description of the data type for the layout of the data in
       * this field.
       */
      GradientVector shape_gradients;

      /**
       * Store the 2nd derivatives of the shape functions at the quadrature
       * points.  See the description of the data type for the layout of the
       * data in this field.
       */
      HessianVector shape_hessians;

      /**
       * Store the 3rd derivatives of the shape functions at the quadrature
       * points.  See the description of the data type for the layout of the
       * data in this field.
       */
      ThirdDerivativeVector shape_3rd_derivatives;

      /**
       * When asked for the value (or gradient, or Hessian) of shape function
       * i's c-th vector component, we need to look it up in the
       * #shape_values, #shape_gradients and #shape_hessians arrays.  The
       * question is where in this array does the data for shape function i,
       * component c reside. This is what this table answers.
       *
       * The format of the table is as follows: - It has dofs_per_cell times
       * n_components entries. - The entry that corresponds to shape function
       * i, component c is <code>i * n_components + c</code>. - The value
       * stored at this position indicates the row in #shape_values and the
       * other tables where the corresponding datum is stored for all the
       * quadrature points.
       *
       * In the general, vector-valued context, the number of components is
       * larger than one, but for a given shape function, not all vector
       * components may be nonzero (e.g., if a shape function is primitive,
       * then exactly one vector component is non-zero, while the others are
       * all zero). For such zero components, #shape_values and friends do not
       * have a row. Consequently, for vector components for which shape
       * function i is zero, the entry in the current table is
       * numbers::invalid_unsigned_int.
       *
       * On the other hand, the table is guaranteed to have at least one valid
       * index for each shape function. In particular, for a primitive finite
       * element, each shape function has exactly one nonzero component and so
       * for each i, there is exactly one valid index within the range
       * <code>[i*n_components, (i+1)*n_components)</code>.
       */
      std::vector<unsigned int> shape_function_to_row_table;
    };
  } // namespace FEValuesImplementation
} // namespace internal


/*@}*/



DEAL_II_NAMESPACE_CLOSE

#endif
