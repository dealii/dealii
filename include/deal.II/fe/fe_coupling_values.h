// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_coupling_values_h
#define dealii_fe_coupling_values_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/general_data_storage.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/point.h>
#include <deal.II/base/std_cxx20/iota_view.h>
#include <deal.II/base/thread_local_storage.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/fe_values_views.h>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <int dim, int spacedim>
class FEValuesBase;
#endif

namespace FEValuesViews
{
  /**
   * A class that stores the renumbering of degrees of freedom and quadrature
   * points for the RenumberedView class. This class is common to all possible
   * renumbered views, and therefore it can be shared between all views that use
   * the same set of renumbering vectors.
   *
   * The renumbering is stored in two vectors, one for the degrees of freedom
   * and one for the quadrature points, and it is used by the
   * FEValuesViews::RenumberedView class.
   */
  struct RenumberingData
  {
    /**
     * Construct a new renumbering data object.
     *
     * The @p dof_renumbering vector is used to renumber the degrees of freedom,
     * while the @p quadrature_renumbering vector is used to renumber the
     * quadrature points. An empty renumbering vector simply means that no
     * renumbering is performed.
     *
     * @note The renumbering vectors are *not* required to match the dimensions
     * of the underlying view, i.e., you could use this class to only run over
     * half of the degrees of freedom, and integrating only over half of the
     * current cell by selecting a subset of quadrature points, if you wish to
     * do so, or to run over some of the degrees of freedom more than once with
     * regard to the underlying view. The important part is that every entry of
     * each renumbering vector is a legal index within the underlying view
     * object. We allow the dof renumbering vector to contain
     * `numbers::invalid_unsigned_int` values, which are ignored, and produce
     * zero values, gradients, hessians, etc. for the corresponding shape
     * function.
     *
     * An example of the renumbering ordering is the following. Assume that you
     * have an FE_Q(1) finite element in dimension two (i.e., four degrees of
     * freedom), and that you are manually taking care of an additional degree
     * of freedom (say, a custom shape function in the middle of the cell),
     * which, for your own convenience, should be numbered locally in the middle
     * of the others. You would like your loops on degrees of freedom to run
     * over all degrees of freedom, and to have the same structure as if the
     * additional degree of freedom was not there. Then the renumbering data
     * object should look like this:
     * @code
     * ...
     * RenumberingData data({{0, 1, numbers::invalid_unsigned_int, 2, 3}});
     * @endcode
     *
     * Using this RenumberingData object, the RenumberedView will return zero
     * when we ask for values, gradients, etc., of the shape function with index
     * two.
     */
    RenumberingData(
      const unsigned int n_inner_dofs = numbers::invalid_unsigned_int,
      const unsigned int n_inner_quadrature_points =
        numbers::invalid_unsigned_int,
      const std::vector<unsigned int> &dof_renumbering        = {},
      const std::vector<unsigned int> &quadrature_renumbering = {});

    /**
     * Make sure we never copy this object, for performance reasons.
     */
    RenumberingData(const RenumberingData &other) = delete;

    /**
     * The number of dofs in the underlying view (before any renumbering).
     */
    const unsigned int n_inner_dofs;

    /**
     * The number of dofs in the renumbered view.
     */
    const unsigned int n_dofs;

    /**
     * The number of quadrature points in the underlying view (before any
     * renumbering).
     */
    const unsigned int n_inner_quadrature_points;

    /**
     * The number of quadrature points in the renumbered view.
     */
    const unsigned int n_quadrature_points;

    /**
     * The renumbering of degrees of freedom.
     */
    const std::vector<unsigned int> dof_renumbering;

    /**
     * The renumbering of quadrature points.
     */
    const std::vector<unsigned int> quadrature_renumbering;

    /**
     * General data storage to store temporary vectors.
     *
     * When the renumbering vectors are non empty, the RenumberedView class may
     * need to construct temporary vectors to store the values of the solution
     * and/or of its gradients with the sizes given by the underlying view
     * object. Unfortunately, we don't know beforehand the Number types with
     * which the vectors will be instantiated, so we cannot use a simple cache
     * internally, and we use a GeneralDataStorage object to avoid allocating
     * memory at each call.
     */
    mutable Threads::ThreadLocalStorage<GeneralDataStorage> data_storage;
  };

  /**
   * A class that provides a renumbered view to a given FEValuesViews object.
   *
   * In general, the order of the degrees of freedom and quadrature points
   * follows the one of the FEValues object, which itself uses the numbering
   * provided by the FiniteElement and Quadrature objects it uses. However, in
   * some cases, it is convenient to group together degrees of freedom and
   * quadrature points in a different order, to select only a subset of the
   * degrees of freedom, or to combine two different sets of degrees of freedom
   * together. This class provides a view to a given FEValuesViews object, where
   * the degrees of freedom and quadrature points are renumbered according to
   * the given RenumberingData object (see there for a documentation of how the
   * renumbering is interpreted).
   *
   * Users will typically not use this class directly, but rather pass an
   * extractor from the FEValuesExtractors namespace to the FECouplingValues
   * class, which returns an object of this type. This is the same mechanism
   * used in the FEValues classes, where passing an extractor returns an
   * FEValuesViews object, and the user rarely instantiates an object of this
   * type.
   */
  template <typename ViewType>
  class RenumberedView
  {
  public:
    /**
     * An alias for the data type of values of the underlying view.
     */
    using value_type = typename ViewType::value_type;

    /**
     * An alias for the data type of gradients of the underlying view.
     */
    using gradient_type = typename ViewType::gradient_type;

    /**
     * An alias for the data type of the product of a @p Number and the values
     * of the underlying view type.
     */
    template <typename Number>
    using solution_value_type =
      typename ViewType::template solution_value_type<Number>;

    /**
     * An alias for the data type of the product of a @p Number and the gradients
     * of the underlying view type.
     */
    template <typename Number>
    using solution_gradient_type =
      typename ViewType::template solution_gradient_type<Number>;

    /**
     * Construct a new RenumberedView object.
     *
     * The renumbering information is taken from the class RenumberingData (see
     * there for a documentation of how the renumbering is interpreted).
     *
     * @note The renumbering information is stored as a reference in this
     * object, so you have to make sure that the RenumberingData object lives
     * longer than this object.
     *
     * @param view The underlying FEValuesViews object.
     * @param data A RenumberingData object, containing renumbering information.
     */
    RenumberedView(const ViewType &view, const RenumberingData &data);

    /**
     * Return the value of the underlying view, for the shape function and
     * quadrature point selected by the arguments.
     *
     * @param shape_function Number of the shape function to be evaluated. Note
     * that this number runs from zero to size of the renumbering vector
     * provided at construction time, or to `dofs_per_cell`, if the renumbering
     * vector is empty.
     *
     * @param q_point Number of the quadrature point at which the function is to
     * be evaluated. Note that this number runs from zero to the size of the
     * renumbering vector provided at construction time, or to
     * `n_quadrature_points`, if the renumbering vector is empty.
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    value_type
    value(const unsigned int shape_function, const unsigned int q_point) const;

    /**
     * Return the gradient of the underlying view, for the shape function and
     * quadrature point selected by the arguments.
     *
     * @param shape_function Number of the shape function to be evaluated. This
     * number runs from zero to size of the renumbering vector provided at
     * construction time, or to `dofs_per_cell`, if the renumbering vector is
     * empty.
     *
     * @param q_point Number of the quadrature point at which the function is to
     * be evaluated. This number runs from zero to the size of the renumbering
     * vector provided at construction time, or to `n_quadrature_points`, if the
     * renumbering vector is empty.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    gradient_type
    gradient(const unsigned int shape_function,
             const unsigned int q_point) const;

    /**
     * Return the values of the underlying view characterized by
     * <tt>fe_function</tt> at the renumbered quadrature points.
     *
     * This function is the equivalent of the FEValuesBase::get_function_values
     * function but it only works on the selected view.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the values of shape functions (i.e., @p value_type) times the
     * type used to store the values of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    template <typename Number>
    void
    get_function_values(const ReadVector<Number>                 &fe_function,
                        std::vector<solution_value_type<Number>> &values) const;

    /**
     * Same as above, but using a vector of renumbered local degree-of-freedom
     * values. In other words, instead of extracting the nodal values of the
     * degrees of freedom located on the current cell from a global vector
     * associated with a DoFHandler object (as the function above does), this
     * function instead takes these local nodal values through its first
     * argument.
     *
     * @param[in] dof_values A vector of local nodal values. This vector must
     *   have a length equal to the size of the dof renumbering vector.
     *
     * @param[out] values A vector of values of the given finite element field,
     *   at the renumbered quadrature points on the current object.
     *
     * @tparam InputVector The @p InputVector type must allow creation of an
     *   ArrayView object from it; this is satisfied by the `std::vector` class,
     *   among others.
     *
     * @dealiiRequiresUpdateFlags{update_values}
     */
    template <class InputVector>
    void
    get_function_values_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_value_type<typename InputVector::value_type>>
        &values) const;

    /**
     * Return the gradients of the underlying view characterized by
     * <tt>fe_function</tt> at the renumbered quadrature points.
     *
     * This function is the equivalent of the
     * FEValuesBase::get_function_gradients function but it only works on the
     * selected view.
     *
     * The data type stored by the output vector must be what you get when you
     * multiply the gradients of shape functions (i.e., @p value_type) times the
     * type used to store the gradients of the unknowns $U_j$ of your finite
     * element vector $U$ (represented by the @p fe_function argument).
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <typename Number>
    void
    get_function_gradients(
      const ReadVector<Number>                    &fe_function,
      std::vector<solution_gradient_type<Number>> &gradients) const;

    /**
     * Same as above, but using a vector of renumbered local degree-of-freedom
     * gradients. In other words, instead of extracting the nodal values of the
     * degrees of freedom located on the current cell from a global vector
     * associated with a DoFHandler object (as the function above does), this
     * function instead takes these local nodal values through its first
     * argument.
     *
     * @param[in] dof_values A vector of local nodal values. This vector must
     *   have a length equal to the size of the dof renumbering vector.
     *
     * @param[out] gradients A vector of gradients of the given finite element
     *   field, at the renumbered quadrature points on the current object.
     *
     * @tparam InputVector The @p InputVector type must allow creation of an
     *   ArrayView object from it; this is satisfied by the `std::vector` class,
     *   among others.
     *
     * @dealiiRequiresUpdateFlags{update_gradients}
     */
    template <class InputVector>
    void
    get_function_gradients_from_local_dof_values(
      const InputVector &dof_values,
      std::vector<solution_gradient_type<typename InputVector::value_type>>
        &gradients) const;

  private:
    /**
     * The data structure that stores the renumbering of the degrees of freedom
     * and of the quadrature points.
     */
    const RenumberingData &data;


    /**
     * Helper function that constructs a unique name for a container, based on a
     * string prefix, on its size, and on the type stored in the container.
     *
     * When the renumbering vectors are non empty, this class may need to
     * construct temporary vectors to store the values of the solution and/or of
     * its gradients with the sizes given by the underlying view object.
     * Unfortunately, we don't know before hand the Number types with which the
     * vectors will be instantiated, so we cannot use a simple cache internally,
     * and we use a GeneralDataStorage object to avoid allocating memory at each
     * call.
     *
     * This function constructs a unique name for temporary containers that will
     * be stored upon the first request in the internal GeneralDataStorage
     * object, to access the
     */
    template <typename Number>
    std::string
    get_unique_container_name(const std::string &prefix,
                              const unsigned int size,
                              const Number      &exemplar_number) const;

    /**
     * Produce an inner vector compatible with the inner view, after copying the
     * values with the correct numbering from the outer vector.
     */
    template <typename InputVector>
    const InputVector &
    outer_to_inner_dofs(const InputVector &outer_vector) const;

    /**
     * Produce an inner vector compatible with the inner view, and zero out its
     * entries if necessary.
     */
    template <typename ValueType>
    std::vector<ValueType> &
    outer_to_inner_values(std::vector<ValueType> &outer_values) const;

    /**
     * Return the outer argument renumbered according to the quadrature
     * renumbering. The values in the inner values are copied to the right
     * position in the outer vector.
     */
    template <typename ValueType>
    void
    inner_to_outer_values(const std::vector<ValueType> &inner_values,
                          std::vector<ValueType>       &outer_values) const;

    /**
     * Store a reference to the underlying view.
     */
    const ViewType &view;
  };
} // namespace FEValuesViews


/**
 * Quadrature coupling options when assembling quadrature formulas for double
 * integrals.
 *
 * When computing the approximation of double integrals of the form
 *
 * \f[
 * \int_{T_1} \int{T_2} K(x_1, x_2) f(x_1) g(x_2) dT_1 dT_2,
 * \f]
 *
 * where $T_1$ and $T_2$ are two arbitrary sets (cells, faces, edges, or any
 * combination thereof), and $K$ is a (possibly singular) coupling kernel, one
 * needs to combine quadrature formulas from two different FEValuesBase objects.
 *
 * This enum class provides a way to specify how the quadrature points and
 * weights should be combined. In general, the two FEValuesBase objects provide
 * different quadrature rules, and these can be interpreted in different ways,
 * depending on the kernel function that is being integrated, and on how the two
 * quadrature rules were constructed.
 *
 * This enum is used in the constructor of FECouplingValues to specify how to
 * interpret and manipulate the quadrature points and weights of the two
 * FEValuesBase objects.
 */
enum class QuadratureCouplingType
{
  /**
   * The FEValuesBase objects provide different quadrature rules, and the
   * resulting quadrature points and weights should be constructed as the tensor
   * product of the two quadrature rules. This is generally used for non-local
   * and non-singular kernels, where the quadrature points of the two
   * FEValuesBase objects are independent of each other.
   */
  tensor_product,

  /**
   * Both FEValuesBase objects provide the same number of (generally different)
   * quadrature points, that should be used as is to implement the double
   * integration. This is equivalent to rewriting double integrals as a single
   * sum over unrolled indices, and it is useful, for example, when performing
   * integration of singular kernels, which require special quadrature rules
   * both on the inside and on the outside integral. An exception is thrown if
   * the two FEValuesBase objects do not provide the same number of quadrature
   * points, but otherwise the two sets of points and weights can be arbitrary.
   */
  unrolled,

  /**
   * The FEValuesBase objects provide the same quadrature rule over the same
   * set, with the same orientation. This is similar to the unrolled case, in
   * the sense that the quadrature formulas are used as is, but both objects are
   * assumed to return the same quadrature points and weights. This is useful
   * when integrating over faces of neighboring cells, and when the reordering
   * of the quadrature points is known to match a priori. An exception is thrown
   * if the quadrature points and weights of the two FEValuesBase objects do not
   * match one-to-one.
   */
  matching,

  /**
   * The FEValuesBase objects provide the same quadrature rule over the same
   * set, but with possibly different orientations. The quadrature points are
   * reordered internally, so that the resulting quadrature points and weights
   * are the same as in the matching case. An exception is thrown if the
   * quadrature points and weights of the two FEValuesBase cannot be reordered
   * to match one-to-one.
   */
  reorder,

  /**
   * The FEValuesBase objects provide partially overlapping quadrature rules
   * over two intersecting sets, with possibly different orientations. The
   * quadrature points are reordered and matched internally, so that the
   * resulting quadrature points and weights are only the ones that match
   * between the two objects. If no overlapping points and weights can be found,
   * an empty set is used to integrate, and the resulting integral is zero. No
   * exception is thrown in this case.
   */
  overlapping,
};

/**
 * DoF coupling options when assembling double integrals.
 *
 * When computing the approximation of double integrals of the form
 *
 * \f[
 * \int_{T_1} \int{T_2} K(x_1, x_2) v_i(x_1) w_j(x_2) dT_1 dT_2,
 * \f]
 *
 * where $T_1$ and $T_2$ are two arbitrary sets (cells, faces, edges, or any
 * combination thereof), and $K$ is a (possibly singular) coupling kernel, one
 * may want to combine degrees from two different FEValuesBase objects (i.e.,
 * basis functions $v_i$ and $w_j$ in the examples above).
 *
 * This enum class provides a way to specify how the degrees of freedom should
 * be combined. There are two cases of interest:
 *
 * 1. the two FEValuesBase objects refer to different DoFHandlers
 * 2. the two FEValuesBase objects refer to the same DoFHandler
 *
 * In the first case, one usually treats the two sets of degrees of freedom as
 * independent of each other, and the resulting matrix is generally rectangular.
 *
 * In the second case, one may choose to treat the two sets of degrees of
 * freedom either as independent or to group them together. A similar approach
 * is used in the FEInterfaceValues class, where the degrees of freedom of the
 * two FEValuesBase objects are grouped together, in a contiguous way, so that
 * the resulting basis functions are interpreted in the following way:
 *
 * \f[
 * \phi_{1,i}(x) = \begin{cases} v_i(x) & \text{ if } i \in [0,n_l) \\
 * 0 & \text{ if } i \in [n_1, n_1+n_2] \end{cases},\quad \phi_{1,i}(x) =
 * \begin{cases} 0(x) & \text{ if } i \in [0,n_1) \\
 * w_{i-n_1}(x) & \text{ if } i \in [n_1, n_1+n_2] \end{cases},
 * \f]
 *
 * where $\phi_{1,i}$ is the first basis function with index $i$ and $n_{1,2}$
 * are the number of local dofs on the first and second FEValuesBase objects.
 *
 * This enum is used in the constructor of FECouplingValues to specify how to
 * interpret and manipulate the local dof indices of the two FEValuesBase
 * objects.
 */
enum class DoFCouplingType
{
  /**
   * The FEValuesBase objects may have different dof indices, possibly indexing
   * different DoFHandler objects, and we are interested in assembling a
   * generally rectangular matrix, where there is no relationship between the
   * two index spaces.
   */
  independent,

  /**
   * The FEValuesBase objects may have different dof indices, but they both
   * index the same DoFHandler objects, and we are interested in assembling them
   * all together in a single matrix. In this coupling type, the DoF indices are
   * grouped together, one after the other, first the first dof indices, and
   * then the second dof indices. This is useful when one wants to assemble four
   * matrices at the same time, corresponding to the four possible combinations
   * of the two FEValuesBase objects, (i.e., first coupled with first, first
   * coupled with second, second coupled with first, and second coupled with
   * second) and then sum them together to obtain the final system. This is
   * similar to what is done in the FEInterfaceValues class.
   */
  contiguous,
};

/**
 * FECouplingValues is a class that facilitates the integration of finite
 * element data between two different finite element objects, possibly living on
 * different grids, and with possibly different topological dimensions (i.e.,
 * cells, faces, edges, and any combination thereof).
 *
 * This class provides a way to simplify the implementation of the following
 * abstract operation:
 *
 * \f[
 * \int_{T_1} \int{T_2} K(x_1, x_2) \phi^1_i(x_1) \phi^2_j(x_2) dT_1 dT_2
 * \f]
 *
 * for three different types of Kernels $K$:
 * - $K(x_1, x_2)$ is a non-singular Kernel function, for example, it is a
 *   function of positive powers $\alpha$ of the distance between the quadrature
 *   points $|x_1-x_2|^\alpha$;
 * - $K(x_1, x_2)$ is a singular Kernel function, for example, it is a function
 *   of negative powers $\alpha$ of the distance between the quadrature points
 *   $|x_1-x_2|^\alpha$;
 * - $K(x_1, x_2)$ is a Dirac delta distribution $\delta(x_1-x_2)$, such that
 *   the integral above is actually a single integral over the intersection of
 *   the two sets $T_1$ and $T_2$.
 *
 * For the first case, one may think that the only natural way to proceed is to
 * compute the double integral by simply nesting two loops:
 * \f[
 * \int_{T_1} \int{T_2} K(x_1, x_2) \phi^1_i(x_1) \phi^2_j(x_2) dT_1 dT_2
 * \approx \sum_{q_1} \sum_{q_2} K(x_1^{q_1}, x_2^{q_2}) \phi^1_i(x_1^{q_1})
 * \phi^2_j(x_2^{q_2}) w_1^{q_1} w_2^{q_2},
 * \f]
 *
 * where $x_1^{q_1}$ and $x_2^{q_2}$ are the quadrature points in $T_1$ and
 * $T_2$ respectively, and $w_1^{q_1}$ and $w_2^{q_2}$ are the corresponding
 * quadrature weights.
 *
 * This, however is not the only way to proceed. In fact, such an integral can
 * be rewritten as a single loop over corresponding elements of two arrays of
 * points with the same length that can be thought of as a single quadrature
 * rule on the set $T_1\times T_2$. For singular kernels, for example, this is
 * often the only way to proceed, since the quadrature formula on $T_1\times
 * T_2$ is usually not written as a tensor product quadrature formula, and one
 * needs to build a custom quadrature formula for this purpose.
 *
 * This class  allows one to treat the three cases above in the same way, and to
 * approximate the integral as follows:
 *
 * \f[
 * \int_{T_1} \int{T_2} K(x_1, x_2) \phi^1_i(x_1) \phi^2_j(x_2) dT_1 dT_2
 * \approx \sum_{i=1}^{N_q} K(x_1^{i}, x_2^{i}) \phi^1_i(x_1^{i})
 * \phi^2_j(x_2^{i}) w_1^{i} w_2^i,
 * \f]
 *
 * Since the triple of objects $(\{q\}, \{w\}, \{\phi\})$ is usually provided by
 * a class derived from the FEValuesBase class, this is the type that the class
 * needs at construction time. $T_1$ and $T_2$ can be two arbitrary cells,
 * faces, or edges belonging to possibly different meshes (or to meshes with
 * different topological dimensions), $\phi^1_i$ and $\phi^2_j$ are basis
 * functions defined on $T_1$ and $T_2$, respectively.
 *
 * The case of the Dirac distribution is when $T_1$ and $T_2$
 * correspond to the common face of two neighboring cells. In this case, this
 * class provides a functionality which is similar to the FEInterfaceValues
 * class, and gives you a way to access values of basis functions on the
 * neighboring cells, as well as their gradients and Hessians, in a unified
 * fashion, on the face.
 *
 * Similarly, this class can be used to couple bulk and surface meshes across
 * the faces of the bulk mesh. In this case, the two FEValuesBase objects will
 * have different topological dimension (i.e., one will be a cell in a
 * co-dimension one triangulation, and the other a face of a bulk grid with
 * co-dimension zero), and the QuadratureCouplingType argument is usually chosen
 * to be QuadratureCouplingType::reorder, since the quadrature points of the two
 * different FEValuesBase objects are not necessarily generated with the same
 * ordering.
 *
 * The type of integral to compute is controlled by the QuadratureCouplingType
 * argument (see the documentation of that enum class for more details), while
 * the type degrees of freedom coupling is controlled by the DoFCouplingType
 * argument (see the documentation of that enum class for more details).
 *
 * As an example usage of this class, consider the a bilinear form of the form:
 * \f[
 * \int_{T_1} \int{T_2} K_1(x_1, x_2) v_i(x_1) u_j(x_2) dT_1 dT_2 +
 * \int_{T_1} \int{T_2} K_2(x_1, x_2) p_i(x_1) q_j(x_2) dT_1 dT_2
 * \f]
 * where the finite dimensional space has two scalar components. We indicate
 * with $v_i$ and $p_i$ the trial functions, and with $u_j$ and
 * $q_j$ the corresponding test functions. $K_1$ and $K_2$ are coupling kernels:
 * such a formulation is used, for example, to write the bilinear forms of
 * Galerkin boundary element methods.
 *
 * The corresponding implementation would look like the following:
 *
 * @code
 * ... // double loop over cells that yields cell_1 and cell_2
 *
 * fe_values_1.reinit(cell_1);
 * fe_values_2.reinit(cell_2);
 *
 * CouplingFEValues<dim, dim, spacedim> cfv(fe_values1, fe_values2,
 *                                          DoFCouplingType::independent,
 *                                          QuadratureCouplingType::tensor_product);
 *
 * FullMatrix<double> local_matrix(fe_values1.dofs_per_cell,
 *                                 fe_values2.dofs_per_cell);
 *
 * // Extractor on first cell
 * const auto v1 = cfv.get_first_extractor(FEValuesExtractor::Scalar(0));
 * const auto p1 = cfv.get_first_extractor(FEValuesExtractor::Scalar(1));
 *
 * // Extractor on second cell
 * const auto u2 = cfv.get_second_extractor(FEValuesExtractor::Scalar(0));
 * const auto q2 = cfv.get_second_extractor(FEValuesExtractor::Scalar(1));
 *
 * ...
 *
 * for (const unsigned int q : cfv.quadrature_point_indices()) {
 *   const auto &[x_q,y_q] = cfv.quadrature_point(q);
 *
 *   for (const unsigned int i : cfv.first_dof_indices())
 *     {
 *       const auto &v_i = cfv[v1].value(i, q);
 *       const auto &p_i = cfv[p1].value(i, q);
 *
 *       for (const unsigned int j : cfv.second_dof_indices())
 *         {
 *           const auto &u_j = cfv[u2].value(j, q);
 *           const auto &q_j = cfv[q2].value(j, q);
 *
 *           local_matrix(i, j) += (K1(x_q, y_q) * v_i * u_j +
 *                                  K2(x_q, y_q) * p_i * q_j) *
 *                                 cfv[0].JxW(q) *
 *                                 cfv[1].JxW(q);
 *         }
 *     }
 * }
 * @endcode
 *
 * In the above loop, the quadrature points of the two FEValuesBase objects are
 * grouped together in a single loop, while the dof indices are kept separate.
 *
 * Internally, this class provides an abstraction to organize coupling terms
 * between two arbitrary FEValuesBase objects, and provides a unified way to
 * access the values of the two basis, and of the two sets of quadrature points
 * and weights. The underlying data is stored in the two FEValuesBase objects,
 * and this class provides access to two RenumberedView objects, which can be
 * obtained by calling the `operator[]` method of this class with the arguments
 * constructed by calling the get_first_extractor() and get_second_extractor()
 * methods of this class with a standard FEValuesExtractors object, i.e.,
 *
 * @code
 * // extract the first scalar component of the basis
 * FEValuesExtractor scalar(0);
 * ...
 *
 * FECouplingValues<dim> fe_coupling(fev1, fev1);
 *
 * // Extractors for the two FEValuesBase objects are returned by the
 * // get_first_extractor() and
 * // get_second_extractor() methods
 *
 * const auto first = fe_coupling.get_first_extractor(scalar);
 * const auto second = fe_coupling.get_second_extractor(scalar);
 *
 * // Now we can use the augmented extractors to access the values of the two
 * // FEValuesBase objects
 * const auto & first_vi = fe_coupling[first].value(i, q);
 * const auto & second_vi = fe_coupling[second].value(i, q);
 * @endcode
 */
template <int dim1, int dim2 = dim1, int spacedim = dim1>
class FECouplingValues
{
public:
  /**
   * Construct the FECouplingValues in an invalid state. Before you can use this
   * object, you must call the reinit() function.
   */
  FECouplingValues();

  /**
   * Construct the FECouplingValues with two arbitrary FEValuesBase objects.
   * This class assumes that the FEValuesBase objects that are given at
   * construction time are initialized and ready to use (i.e., that you have
   * called the reinit() function on them before calling this constructor).
   *
   * Notice that the actual renumbering of the degrees of freedom and quadrature
   * points is done at construction time, or upon calling the reinit() function.
   * If you change the underlying FEValuesBase objects after construction, you
   * must call the reinit() function to update the renumbering. This may or may
   * not be necessary, depending on the type of coupling that you are using.
   *
   * This really depends on the application and on the specific type of
   * coupling. For example, for volume/volume coupling, i.e., for operators with
   * non-local and non-singular kernels of type
   * \f[
   * \int_K \int_T f(\phi_i(x)-\phi_j(y), x-y) dx dy
   * \f]
   * you may initialize FECouplingValues once, and just reinit the underlying
   * FEValuesBase objects on different cells K and T, without the need to
   * recompute the coupling (i.e., the numbering is always the same, and nothing
   * differs from what happened in the first call).
   *
   * For cell/surface coupling, the same cell may couple with different faces,
   * so the renumbering must be really computed from scratch for each pair of
   * FEValuesBase objects, so reinitializing the underlying cells and faces will
   * make the renumbering itself invalid, and FECouplingValues must be
   * reinitialized (o constructed from scratch) after calling
   * `fe_values_1.reinit(K)` and `fe_values_1.reinit(T)`.
   *
   * @param fe_values_1 The first FEValuesBase object.
   * @param fe_values_2 The second FEValuesBase object.
   * @param dof_coupling_type The type of dof coupling to use.
   * @param quadrature_coupling_type The type of quadrature to use for the coupling.
   */
  FECouplingValues(
    const FEValuesBase<dim1, spacedim> &fe_values_1,
    const FEValuesBase<dim2, spacedim> &fe_values_2,
    const DoFCouplingType &dof_coupling_type = DoFCouplingType::independent,
    const QuadratureCouplingType &quadrature_coupling_type =
      QuadratureCouplingType::tensor_product);

  /**
   * Reinitialize the FECouplingValues with two arbitrary FEValuesBase objects.
   * The FEValuesBase objects must be initialized and ready to use, i.e., you
   * must have called the reinit() function on them before calling this method.
   *
   * This method computes the actual renumbering of the degrees of freedom and
   * quadrature points. If you change the underlying FEValuesBase objects after
   * calling this method, you may need to call the reinit() function to update
   * the renumbering. This may or may not be necessary, depending on the type of
   * coupling that you are using.
   *
   * This really depends on the application and on the specific type of
   * coupling. For example, for volume/volume coupling, i.e., for operators with
   * non-local and non-singular kernels of type
   * \f[
   * \int_K \int_T f(\phi_i(x)-\phi_j(y), x-y) dx dy
   * \f]
   * you may initialize FECouplingValues once, and just reinit the underlying
   * FEValuesBase objects on different cells K and T, without the need to
   * recompute the coupling (i.e., the numbering is always the same, and nothing
   * differs from what happened in the first call).
   *
   * For cell/surface coupling, the same cell may couple with different faces,
   * so the renumbering must be really computed from scratch for each pair of
   * FEValuesBase objects, so reinitializing the underlying cells and faces will
   * make the renumbering itself invalid, and FECouplingValues must be
   * reinitialized (o constructed from scratch) after calling
   * `fe_values_1.reinit(K)` and `fe_values_1.reinit(T)`.
   *
   * @param fe_values_1 The first FEValuesBase object.
   * @param fe_values_2 The second FEValuesBase object.
   * @param dof_coupling_type The type of dof coupling to use.
   * @param quadrature_coupling_type The type of quadrature to use for the coupling.
   */
  void
  reinit(
    const FEValuesBase<dim1, spacedim> &fe_values_1,
    const FEValuesBase<dim2, spacedim> &fe_values_2,
    const DoFCouplingType &dof_coupling_type = DoFCouplingType::independent,
    const QuadratureCouplingType &quadrature_coupling_type =
      QuadratureCouplingType::tensor_product);


  /**
   * Return a FirstCoupling object that can be used to extract values from the
   * first FEValuesBase object.
   */
  template <typename Extractor>
  FEValuesExtractors::FirstCoupling<Extractor>
  get_first_extractor(const Extractor &extractor) const;

  /**
   * Return a SecondCoupling object that can be used to extract values from the
   * second FEValuesBase object.
   */
  template <typename Extractor>
  FEValuesExtractors::SecondCoupling<Extractor>
  get_second_extractor(const Extractor &extractor) const;

  /**
   * Return the value of JxW at the given quadrature point.
   *
   * @dealiiRequiresUpdateFlags{update_JxW_values}
   */
  double
  JxW(const unsigned int quadrature_point) const;

  /**
   * Return the two quadrature points in real space at the given quadrature
   * point index, corresponding to a quadrature point in the set $T_1\times
   * T_2$.
   *
   * @dealiiRequiresUpdateFlags{update_quadrature_points}
   */
  std::pair<Point<spacedim>, Point<spacedim>>
  quadrature_point(const unsigned int quadrature_point) const;

  /**
   * Return an object that can be thought of as an array containing all
   * indices from zero to `n_quadrature_points`. This allows to write code
   * using range-based `for` loops.
   *
   * @see CPP11
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  quadrature_point_indices() const;

  /**
   * Return an object that can be thought of as an array containing all
   * indices from zero (inclusive) to `n_first_dofs()` (exclusive).
   * This allows one to write code using range-based `for` loops.
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  first_dof_indices() const;

  /**
   * Return an object that can be thought of as an array containing all
   * indices from zero (inclusive) to `n_second_dofs()` (exclusive).
   * This allows one to write code using range-based `for` loops.
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  second_dof_indices() const;

  /**
   * Return an object that can be thought of as an array containing all indices
   * from zero (inclusive) to `n_coupling_dof_indices()` (exclusive). This
   * allows one to write code using range-based `for` loops.
   *
   * This function only makes sense when the dof coupling type is
   * DoFCouplingType::contiguous, and, in this case, n_coupling_dof_indices() is
   * the  sum of n_first_dofs() and n_second_dofs(). If
   * DoFCouplingType::independent is used, you should call first_dof_indices()
   * or second_dof_indices() instead, and calling this function will throw an
   * exception.
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  coupling_dof_indices() const;

  /**
   * Return the set of joint DoF indices, possibly offset by the given values.
   *
   * For coupled problems, the DoF indices of the two FEValuesBase objects are
   * generally constructed from different DoFHandler objects. According to how
   * the user wants to construct the coupling matrices, the two dof indices may
   * map local dof indices to global dof indices starting from zero (i.e., we
   * are assembling a single rectangular matrix, and the first dof indices run
   * on rows of the matrix, while the second dof indices run on columns of the
   * matrix), or the user may be wanting to map the coupling matrix to a block
   * of a larger block matrix. In this second case the dof indices of either or
   * both blocks may need to be shifted by a fixed amount, to fall within the
   * right range of the target block matrix. This is usually required if the
   * user wants to keep their solution vectors in a single block vector, where
   * each block corresponds to a different DoFHandler object. In this case, all
   * dofs referring to the first FEValuesBase object are offset by the value of
   * dofs_offset_1, while all dofs referring to the second FEValuesBase object
   * are offset by the value of dofs_offset_2, i.e., the following relationship
   * is satisfied:
   * @code
   * fe_values_1.get_dof_indices(dof_indices_1);
   * fe_values_2.get_dof_indices(dof_indices_2);
   * auto coupling_dof_indices = fe_coupling_values.
   *      get_coupling_dof_indices(dof_indices_1, dof_indices_2,
   *                               dofs_offset_1, dofs_offset_2);
   * for(const auto &i : fe_coupling_values.coupling_dof_indices()) {
   *    auto [id1, id2] = fe_coupling_values.coupling_dof_to_dof_indices(i);
   *    if (id1 != numbers::invalid_unsigned_int)
   *      // always true
   *      coupling_dof_indices[i] == dof_indices_1[id1] + dofs_offset_1;
   *    if (id2 != numbers::invalid_unsigned_int)
   *      // always true
   *      coupling_dof_indices[i] == dof_indices_2[id2] + dofs_offset_2;
   * }
   * @endcode
   *
   * This method only makes sense when the dof coupling type is
   * DoFCouplingType::contiguous. If DoFCouplingType::independent is used,
   * calling this function will result in an exception.
   */
  std::vector<types::global_dof_index>
  get_coupling_dof_indices(
    const std::vector<types::global_dof_index> &dof_indices_1,
    const std::vector<types::global_dof_index> &dof_indices_2,
    const types::global_dof_index               dofs_offset_1 = 0,
    const types::global_dof_index               dofs_offset_2 = 0) const;

  /**
   * Convert a coupling dof index into the corresponding local DoF indices of
   * the two FEValuesObjects. If a DoF is only active on one of the
   * FEValuesObjects, the other index will be numbers::invalid_unsigned_int.
   */
  std::pair<unsigned int, unsigned int>
  coupling_dof_to_dof_indices(const unsigned int coupling_dof_index) const;

  /**
   * Convert a quadrature index into the corresponding local quadrature indices
   * of the two FEValuesObjects. Both indices are guaranteed to be valid within
   * the corresponding FEValuesObject.
   */
  std::pair<unsigned int, unsigned int>
  coupling_quadrature_to_quadrature_indices(
    const unsigned int quadrature_point) const;

  /**
   * @name Extractors Methods to extract individual components
   * @{
   */

  /**
   * Create a combined view of the first FECouplingValues object that represents
   * a view of the possibly vector-valued finite element. The concept of views
   * is explained in the documentation of the namespace FEValuesViews.
   */
  template <typename Extractor>
  const FEValuesViews::RenumberedView<
    FEValuesViews::View<dim1, spacedim, Extractor>>
  operator[](
    const FEValuesExtractors::FirstCoupling<Extractor> &extractor) const;

  /**
   * Create a combined view of the second FECouplingValues object that
   * represents a view of the possibly vector-valued finite element. The concept
   * of views is explained in the documentation of the namespace FEValuesViews.
   */
  template <typename Extractor>
  const FEValuesViews::RenumberedView<
    FEValuesViews::View<dim2, spacedim, Extractor>>
  operator[](
    const FEValuesExtractors::SecondCoupling<Extractor> &extractor) const;

  /**
   * @}
   */

  /**
   * Return the number of coupling DoF indices. If DoFCouplingType::independent
   * is used, this function returns numbers::invalid_unsigned_int, otherwise it
   * returns the sum of n_first_dofs() and n_second_dofs().
   */
  unsigned int
  n_coupling_dofs() const;

  /**
   * Return the number of first DoF indices. This generally coincides with
   * n_dofs_per_cell of the first FEValuesBase object.
   */
  unsigned int
  n_first_dofs() const;

  /**
   * Return the number of second DoF indices. This generally coincides with
   * n_dofs_per_cell of the second FEValuesBase object.
   */
  unsigned int
  n_second_dofs() const;

  /**
   * Return the number of quadrature points.
   */
  unsigned int
  n_quadrature_points() const;

private:
  /**
   * The dof coupling type used by this object.
   */
  DoFCouplingType dof_coupling_type;

  /**
   * The quadrature coupling type used by this object.
   */
  QuadratureCouplingType quadrature_coupling_type;

  /**
   * Pointer to first FEValuesBase object.
   */
  ObserverPointer<const FEValuesBase<dim1, spacedim>> first_fe_values;

  /**
   * Pointer to second FEValuesBase object.
   */
  ObserverPointer<const FEValuesBase<dim2, spacedim>> second_fe_values;

  /**
   * Renumbering data for the first FEValuesBase object.
   */
  std::unique_ptr<const FEValuesViews::RenumberingData> first_renumbering_data;

  /**
   * Renumbering data for the second FEValuesBase object.
   */
  std::unique_ptr<const FEValuesViews::RenumberingData> second_renumbering_data;

  /**
   * Number of quadrature points.
   */
  unsigned int n_quadrature_points_;

  /**
   * Number of coupling DoF indices. If DoFCouplingType::independent is used,
   * this is numbers::invalid_unsigned_int, while it is n_first_dofs() +
   * n_second_dofs() in the the case of DoFCouplingType::contiguous.
   */
  unsigned int n_coupling_dofs_;
};


#ifndef DOXYGEN


/*------------------------ Inline functions: namespace FEValuesViews --------*/

namespace FEValuesViews
{
  RenumberingData::RenumberingData(
    const unsigned int               n_inner_dofs,
    const unsigned int               n_inner_quadrature_points,
    const std::vector<unsigned int> &dof_renumbering,
    const std::vector<unsigned int> &quadrature_renumbering)
    : n_inner_dofs(n_inner_dofs)
    , n_dofs(dof_renumbering.empty() ? n_inner_dofs : dof_renumbering.size())
    , n_inner_quadrature_points(n_inner_quadrature_points)
    , n_quadrature_points(quadrature_renumbering.empty() ?
                            n_inner_quadrature_points :
                            quadrature_renumbering.size())
    , dof_renumbering(dof_renumbering)
    , quadrature_renumbering(quadrature_renumbering)
  {
    // Check that the renumbering vectors are valid.
    if constexpr (running_in_debug_mode())
      {
        // While for dofs we admit invalid values, this is not the case for
        // quadrature points.
        for (const auto i : dof_renumbering)
          Assert(i < n_inner_dofs || i == numbers::invalid_unsigned_int,
                 ExcIndexRange(i, 0, n_inner_dofs));

        for (const auto q : quadrature_renumbering)
          AssertIndexRange(q, n_inner_quadrature_points);
      }
  }



  template <typename ViewType>
  RenumberedView<ViewType>::RenumberedView(const ViewType        &view,
                                           const RenumberingData &data)
    : view(view)
    , data(data)
  {}



  template <typename ViewType>
  template <typename Number>
  inline std::string
  RenumberedView<ViewType>::get_unique_container_name(
    const std::string &prefix,
    const unsigned int size,
    const Number      &exemplar_number) const
  {
    return prefix + "_" + Utilities::int_to_string(size) + "_" +
           Utilities::type_to_string(exemplar_number);
  }



  template <typename ViewType>
  typename RenumberedView<ViewType>::value_type
  RenumberedView<ViewType>::value(const unsigned int shape_function,
                                  const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, data.n_dofs);
    AssertIndexRange(q_point, data.n_quadrature_points);

    const auto inner_shape_function = data.dof_renumbering.empty() ?
                                        shape_function :
                                        data.dof_renumbering[shape_function];
    const auto inner_q_point        = data.quadrature_renumbering.empty() ?
                                        q_point :
                                        data.quadrature_renumbering[q_point];
    if (inner_shape_function == numbers::invalid_unsigned_int)
      return value_type(0);
    else
      {
        AssertIndexRange(inner_shape_function, data.n_inner_dofs);
        AssertIndexRange(inner_q_point, data.n_inner_quadrature_points);
        return view.value(inner_shape_function, inner_q_point);
      }
  }



  template <typename ViewType>
  typename RenumberedView<ViewType>::gradient_type
  RenumberedView<ViewType>::gradient(const unsigned int shape_function,
                                     const unsigned int q_point) const
  {
    AssertIndexRange(shape_function, data.n_dofs);
    AssertIndexRange(q_point, data.n_quadrature_points);

    const auto inner_shape_function = data.dof_renumbering.empty() ?
                                        shape_function :
                                        data.dof_renumbering[shape_function];
    const auto inner_q_point        = data.quadrature_renumbering.empty() ?
                                        q_point :
                                        data.quadrature_renumbering[q_point];
    if (inner_shape_function == numbers::invalid_unsigned_int)
      return gradient_type();
    else
      return view.gradient(inner_shape_function, inner_q_point);
  }



  template <typename ViewType>
  template <typename ValueType>
  std::vector<ValueType> &
  RenumberedView<ViewType>::outer_to_inner_values(
    std::vector<ValueType> &outer_values) const
  {
    AssertDimension(outer_values.size(), data.n_quadrature_points);
    if (data.quadrature_renumbering.empty())
      {
        return outer_values;
      }
    else
      {
        const auto name =
          get_unique_container_name("RenumberedView::outer_to_inner_values",
                                    data.n_inner_quadrature_points,
                                    outer_values[0]);
        auto &inner_values =
          data.data_storage.get()
            .template get_or_add_object_with_name<std::vector<ValueType>>(
              name, data.n_inner_quadrature_points);
        return inner_values;
      }
  }



  template <typename ViewType>
  template <typename VectorType>
  const VectorType &
  RenumberedView<ViewType>::outer_to_inner_dofs(
    const VectorType &outer_dofs) const
  {
    AssertDimension(outer_dofs.size(), data.n_dofs);
    if (data.dof_renumbering.empty())
      {
        return outer_dofs;
      }
    else
      {
        const auto name =
          get_unique_container_name("RenumberedView::outer_to_inner_dofs",
                                    data.n_inner_dofs,
                                    outer_dofs[0]);

        auto &inner_dofs = data.data_storage.get()
                             .template get_or_add_object_with_name<VectorType>(
                               name, data.n_inner_dofs);
        for (unsigned int i = 0; i < data.n_dofs; ++i)
          {
            const auto inner_i = data.dof_renumbering[i];
            if (inner_i != numbers::invalid_unsigned_int)
              {
                AssertIndexRange(inner_i, data.n_inner_dofs);
                inner_dofs[inner_i] = outer_dofs[i];
              }
          }
        return inner_dofs;
      }
  }



  template <typename ViewType>
  template <typename ValueType>
  void
  RenumberedView<ViewType>::inner_to_outer_values(
    const std::vector<ValueType> &inner_values,
    std::vector<ValueType>       &outer_values) const
  {
    AssertDimension(outer_values.size(), data.n_quadrature_points);
    AssertDimension(inner_values.size(), data.n_inner_quadrature_points);
    if (data.quadrature_renumbering.empty())
      {
        Assert(&inner_values == &outer_values, ExcInternalError());
        return;
      }
    for (unsigned int i = 0; i < data.quadrature_renumbering.size(); ++i)
      {
        outer_values[i] = inner_values[data.quadrature_renumbering[i]];
      }
  }



  template <typename ViewType>
  template <typename Number>
  void
  RenumberedView<ViewType>::get_function_values(
    const ReadVector<Number>                 &fe_function,
    std::vector<solution_value_type<Number>> &values) const
  {
    auto &inner_values = outer_to_inner_values(values);
    view.get_function_values(fe_function, inner_values);
    inner_to_outer_values(inner_values, values);
  }



  template <typename ViewType>
  template <typename InputVector>
  void
  RenumberedView<ViewType>::get_function_values_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_value_type<typename InputVector::value_type>> &values)
    const
  {
    const auto &inner_dof_values = outer_to_inner_dofs(dof_values);
    auto       &inner_values     = outer_to_inner_values(values);

    view.get_function_values_from_local_dof_values(inner_dof_values,
                                                   inner_values);
    inner_to_outer_values(inner_values, values);
  }


  template <typename ViewType>
  template <typename Number>
  void
  RenumberedView<ViewType>::get_function_gradients(
    const ReadVector<Number>                    &fe_function,
    std::vector<solution_gradient_type<Number>> &gradients) const
  {
    auto &inner_gradients = outer_to_inner_values(gradients);
    view.get_function_gradients(fe_function, inner_gradients);
    inner_to_outer_values(inner_gradients, gradients);
  }



  template <typename ViewType>
  template <typename InputVector>
  void
  RenumberedView<ViewType>::get_function_gradients_from_local_dof_values(
    const InputVector &dof_values,
    std::vector<solution_gradient_type<typename InputVector::value_type>>
      &gradients) const
  {
    const auto &inner_dof_values = outer_to_inner_dofs(dof_values);
    auto       &inner_gradients  = outer_to_inner_values(gradients);

    view.get_function_gradients_from_local_dof_values(inner_dof_values,
                                                      inner_gradients);
    inner_to_outer_values(inner_gradients, gradients);
  }
} // namespace FEValuesViews

/*-------------- Inline functions FECouplingValues ---------------------*/

template <int dim1, int dim2, int spacedim>
std::vector<types::global_dof_index>
FECouplingValues<dim1, dim2, spacedim>::get_coupling_dof_indices(
  const std::vector<types::global_dof_index> &dof_indices_1,
  const std::vector<types::global_dof_index> &dof_indices_2,
  const types::global_dof_index               dofs_offset_1,
  const types::global_dof_index               dofs_offset_2) const
{
  AssertThrow(
    n_coupling_dofs_ != numbers::invalid_unsigned_int,
    ExcMessage(
      "Dofs are independent. You cannot ask for coupling dof indices."));
  AssertDimension(dof_indices_1.size(), first_fe_values->dofs_per_cell);
  AssertDimension(dof_indices_2.size(), second_fe_values->dofs_per_cell);

  std::vector<types::global_dof_index> coupling_dof_indices(
    dof_indices_1.size() + dof_indices_2.size());
  unsigned int idx = 0;
  for (const auto &i : dof_indices_1)
    coupling_dof_indices[idx++] = i + dofs_offset_1;
  for (const auto &i : dof_indices_2)
    coupling_dof_indices[idx++] = i + dofs_offset_2;
  return coupling_dof_indices;
}



template <int dim1, int dim2, int spacedim>
FECouplingValues<dim1, dim2, spacedim>::FECouplingValues()
  : first_fe_values(nullptr)
  , second_fe_values(nullptr)
  , quadrature_coupling_type(QuadratureCouplingType::unrolled)
  , n_quadrature_points_(0)
  , n_coupling_dofs_(numbers::invalid_unsigned_int)
{}



template <int dim1, int dim2, int spacedim>
FECouplingValues<dim1, dim2, spacedim>::FECouplingValues(
  const FEValuesBase<dim1, spacedim> &fe_values_1,
  const FEValuesBase<dim2, spacedim> &fe_values_2,
  const DoFCouplingType              &dof_coupling_type,
  const QuadratureCouplingType       &quadrature_coupling_type)
  : FECouplingValues() // delegate to other constructor
{
  reinit(fe_values_1, fe_values_2, dof_coupling_type, quadrature_coupling_type);
}



template <int dim1, int dim2, int spacedim>
void
FECouplingValues<dim1, dim2, spacedim>::reinit(
  const FEValuesBase<dim1, spacedim> &fe_values_1,
  const FEValuesBase<dim2, spacedim> &fe_values_2,
  const DoFCouplingType              &dof_coupling_type,
  const QuadratureCouplingType       &quadrature_coupling_type)
{
  first_fe_values                = &fe_values_1;
  second_fe_values               = &fe_values_2;
  this->dof_coupling_type        = dof_coupling_type;
  this->quadrature_coupling_type = quadrature_coupling_type;

  // Gather information about the inner objects
  unsigned int first_n_inner_dofs  = fe_values_1.dofs_per_cell;
  unsigned int second_n_inner_dofs = fe_values_2.dofs_per_cell;

  unsigned int first_n_inner_quadrature_points =
    fe_values_1.n_quadrature_points;
  unsigned int second_n_inner_quadrature_points =
    fe_values_2.n_quadrature_points;

  // Initialize counters and renumbering vectors to zero
  std::vector<unsigned int> first_dofs_map;
  std::vector<unsigned int> second_dofs_map;

  std::vector<unsigned int> first_quad_map;
  std::vector<unsigned int> second_quad_map;

  // Local relative tolerance for comparing quadrature points
  const double tol = std::min(fe_values_1.get_cell()->diameter(),
                              fe_values_2.get_cell()->diameter()) *
                     1e-6;

  // Now fill the renumbering vectors
  switch (dof_coupling_type)
    {
      case DoFCouplingType::independent:
        {
          n_coupling_dofs_ = numbers::invalid_unsigned_int;
          first_dofs_map.clear();
          second_dofs_map.clear();
          break;
        }
      case DoFCouplingType::contiguous:
        {
          n_coupling_dofs_ =
            fe_values_1.dofs_per_cell + fe_values_2.dofs_per_cell;

          first_dofs_map.resize(n_coupling_dofs_);
          second_dofs_map.resize(n_coupling_dofs_);

          unsigned int idx = 0;
          for (const unsigned int &i : fe_values_1.dof_indices())
            {
              first_dofs_map[idx]    = i;
              second_dofs_map[idx++] = numbers::invalid_unsigned_int;
            }
          for (const unsigned int &i : fe_values_2.dof_indices())
            {
              first_dofs_map[idx]    = numbers::invalid_unsigned_int;
              second_dofs_map[idx++] = i;
            }
          // Make sure we have the right number of dofs
          AssertDimension(idx, n_coupling_dofs_);
          break;
        }
      default:
        AssertThrow(false, ExcNotImplemented());
    }
  // Compute the quadrature maps
  switch (quadrature_coupling_type)
    {
      case QuadratureCouplingType::tensor_product:
        {
          const auto &quadrature_points_1 = fe_values_1.get_quadrature_points();
          const auto &quadrature_points_2 = fe_values_2.get_quadrature_points();

          n_quadrature_points_ =
            quadrature_points_1.size() * quadrature_points_2.size();

          first_quad_map.resize(n_quadrature_points_);
          second_quad_map.resize(n_quadrature_points_);

          unsigned int idx = 0;
          for (const unsigned int &i : fe_values_1.quadrature_point_indices())
            for (const unsigned int &j : fe_values_2.quadrature_point_indices())
              {
                first_quad_map[idx]  = i;
                second_quad_map[idx] = j;
                ++idx;
              }
          break;
        }
      case QuadratureCouplingType::unrolled:
        {
          Assert(fe_values_1.get_quadrature_points().size() ==
                   fe_values_2.get_quadrature_points().size(),
                 ExcMessage("The two FEValuesBase objects must have the same "
                            "number of quadrature points"));

          n_quadrature_points_ = fe_values_1.get_quadrature_points().size();

          first_quad_map.clear();
          second_quad_map.clear();

          break;
        }
      case QuadratureCouplingType::matching:
        {
          const auto &quadrature_points_1 = fe_values_1.get_quadrature_points();
          const auto &quadrature_points_2 = fe_values_2.get_quadrature_points();

          Assert(quadrature_points_1.size() == quadrature_points_2.size(),
                 ExcMessage("The two FEValuesBase objects must have the same "
                            "number of quadrature points"));

          for (const unsigned int &i : fe_values_1.quadrature_point_indices())
            {
              Assert(quadrature_points_1[i].distance(quadrature_points_2[i]) <
                       tol,
                     ExcMessage(
                       "The two FEValuesBase objects must have the same "
                       "quadrature points"));
            }

          n_quadrature_points_ = fe_values_1.get_quadrature_points().size();

          first_quad_map.clear();
          second_quad_map.clear();

          break;
        }
      case QuadratureCouplingType::reorder:
        {
          const auto &quadrature_points_1 = fe_values_1.get_quadrature_points();
          const auto &quadrature_points_2 = fe_values_2.get_quadrature_points();

          Assert(quadrature_points_1.size() == quadrature_points_2.size(),
                 ExcMessage("The two FEValuesBase objects must have the same "
                            "number of quadrature points"));

          n_quadrature_points_ = fe_values_1.get_quadrature_points().size();

          // The first is the id. The second is renumbered.
          first_quad_map.clear();
          second_quad_map.resize(fe_values_1.get_quadrature_points().size());

          // [TODO]: Avoid quadratic complexity here
          for (const unsigned int &i : fe_values_1.quadrature_point_indices())
            {
              auto id = numbers::invalid_unsigned_int;
              for (const unsigned int &j :
                   fe_values_2.quadrature_point_indices())
                if (quadrature_points_1[i].distance(quadrature_points_2[j]) <
                    tol)
                  {
                    id                 = i;
                    second_quad_map[i] = j;
                    break;
                  }
              Assert(id != numbers::invalid_unsigned_int,
                     ExcMessage(
                       "The two FEValuesBase objects must have the same "
                       "quadrature points, even if not in the same order."));
            }
          break;
        }
      case QuadratureCouplingType::overlapping:
        {
          const auto &quadrature_points_1 = fe_values_1.get_quadrature_points();
          const auto &quadrature_points_2 = fe_values_2.get_quadrature_points();

          // [TODO]: Avoid quadratic complexity here
          for (const unsigned int &i : fe_values_1.quadrature_point_indices())
            {
              for (const unsigned int &j :
                   fe_values_2.quadrature_point_indices())
                if (quadrature_points_1[i].distance(quadrature_points_2[j]) <
                    tol)
                  {
                    first_quad_map.emplace_back(i);
                    second_quad_map.emplace_back(j);
                    break;
                  }
            }
          n_quadrature_points_ = first_quad_map.size();

          break;
        }
      default:
        Assert(false, ExcInternalError());
    }

  first_renumbering_data = std::make_unique<FEValuesViews::RenumberingData>(
    first_n_inner_dofs,
    first_n_inner_quadrature_points,
    first_dofs_map,
    first_quad_map);

  second_renumbering_data = std::make_unique<FEValuesViews::RenumberingData>(
    second_n_inner_dofs,
    second_n_inner_quadrature_points,
    second_dofs_map,
    second_quad_map);
}


template <int dim1, int dim2, int spacedim>
inline double
FECouplingValues<dim1, dim2, spacedim>::JxW(const unsigned int q) const
{
  AssertIndexRange(q, n_quadrature_points_);

  const auto first_q = first_renumbering_data->quadrature_renumbering.empty() ?
                         q :
                         first_renumbering_data->quadrature_renumbering[q];

  if (quadrature_coupling_type == QuadratureCouplingType::tensor_product ||
      quadrature_coupling_type == QuadratureCouplingType::unrolled)
    {
      const auto second_q =
        second_renumbering_data->quadrature_renumbering.empty() ?
          q :
          second_renumbering_data->quadrature_renumbering[q];
      return first_fe_values->JxW(first_q) * second_fe_values->JxW(second_q);
    }
  else
    return first_fe_values->JxW(first_q);
}



template <int dim1, int dim2, int spacedim>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
FECouplingValues<dim1, dim2, spacedim>::quadrature_point_indices() const
{
  return std_cxx20::ranges::iota_view<unsigned int, unsigned int>(
    0U, n_quadrature_points_);
}



template <int dim1, int dim2, int spacedim>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
FECouplingValues<dim1, dim2, spacedim>::coupling_dof_indices() const
{
  AssertThrow(n_coupling_dofs_ != numbers::invalid_unsigned_int,
              ExcMessage(
                "Dofs are independent. You cannot ask for coupling dofs."));
  return std_cxx20::ranges::iota_view<unsigned int, unsigned int>(
    0U, n_coupling_dofs_);
}



template <int dim1, int dim2, int spacedim>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
FECouplingValues<dim1, dim2, spacedim>::first_dof_indices() const
{
  return std_cxx20::ranges::iota_view<unsigned int, unsigned int>(
    0U, n_first_dofs());
}



template <int dim1, int dim2, int spacedim>
inline std_cxx20::ranges::iota_view<unsigned int, unsigned int>
FECouplingValues<dim1, dim2, spacedim>::second_dof_indices() const
{
  return std_cxx20::ranges::iota_view<unsigned int, unsigned int>(
    0U, n_second_dofs());
}



template <int dim1, int dim2, int spacedim>
inline unsigned int
FECouplingValues<dim1, dim2, spacedim>::n_coupling_dofs() const
{
  AssertThrow(n_coupling_dofs_ != numbers::invalid_unsigned_int,
              ExcMessage(
                "Dofs are independent. You cannot ask for coupling dofs."));
  return n_coupling_dofs_;
}



template <int dim1, int dim2, int spacedim>
inline unsigned int
FECouplingValues<dim1, dim2, spacedim>::n_first_dofs() const
{
  return first_renumbering_data->n_dofs;
}



template <int dim1, int dim2, int spacedim>
inline unsigned int
FECouplingValues<dim1, dim2, spacedim>::n_second_dofs() const
{
  return second_renumbering_data->n_dofs;
}



template <int dim1, int dim2, int spacedim>
inline unsigned int
FECouplingValues<dim1, dim2, spacedim>::n_quadrature_points() const
{
  return n_quadrature_points_;
}



template <int dim1, int dim2, int spacedim>
std::pair<Point<spacedim>, Point<spacedim>>
FECouplingValues<dim1, dim2, spacedim>::quadrature_point(
  const unsigned int quadrature_point) const
{
  AssertIndexRange(quadrature_point, n_quadrature_points_);
  const auto first_q = first_fe_values->quadrature_point(
    first_renumbering_data->quadrature_renumbering.empty() ?
      quadrature_point :
      first_renumbering_data->quadrature_renumbering[quadrature_point]);

  const auto second_q = second_fe_values->quadrature_point(
    second_renumbering_data->quadrature_renumbering.empty() ?
      quadrature_point :
      second_renumbering_data->quadrature_renumbering[quadrature_point]);

  return {first_q, second_q};
}



template <int dim1, int dim2, int spacedim>
std::pair<unsigned int, unsigned int>
FECouplingValues<dim1, dim2, spacedim>::
  coupling_quadrature_to_quadrature_indices(
    const unsigned int quadrature_point) const
{
  AssertIndexRange(quadrature_point, n_quadrature_points_);
  const auto first_id =
    first_renumbering_data->quadrature_renumbering.empty() ?
      quadrature_point :
      first_renumbering_data->quadrature_renumbering[quadrature_point];

  const auto second_id =
    second_renumbering_data->quadrature_renumbering.empty() ?
      quadrature_point :
      second_renumbering_data->quadrature_renumbering[quadrature_point];
  return std::make_pair(first_id, second_id);
}



template <int dim1, int dim2, int spacedim>
std::pair<unsigned int, unsigned int>
FECouplingValues<dim1, dim2, spacedim>::coupling_dof_to_dof_indices(
  const unsigned int coupling_dof_index) const
{
  AssertIndexRange(coupling_dof_index, n_coupling_dofs_);
  const auto first_id =
    first_renumbering_data->dof_renumbering.empty() ?
      coupling_dof_index :
      first_renumbering_data->dof_renumbering[coupling_dof_index];

  const auto second_id =
    second_renumbering_data->dof_renumbering.empty() ?
      coupling_dof_index :
      second_renumbering_data->dof_renumbering[coupling_dof_index];
  return std::make_pair(first_id, second_id);
}



template <int dim1, int dim2, int spacedim>
template <typename Extractor>
inline FEValuesExtractors::FirstCoupling<Extractor>
FECouplingValues<dim1, dim2, spacedim>::get_first_extractor(
  const Extractor &extractor) const
{
  return FEValuesExtractors::FirstCoupling<Extractor>(extractor);
}



template <int dim1, int dim2, int spacedim>
template <typename Extractor>
inline FEValuesExtractors::SecondCoupling<Extractor>
FECouplingValues<dim1, dim2, spacedim>::get_second_extractor(
  const Extractor &extractor) const
{
  return FEValuesExtractors::SecondCoupling<Extractor>(extractor);
}



template <int dim1, int dim2, int spacedim>
template <typename Extractor>
inline const FEValuesViews::RenumberedView<
  FEValuesViews::View<dim1, spacedim, Extractor>>
FECouplingValues<dim1, dim2, spacedim>::operator[](
  const FEValuesExtractors::FirstCoupling<Extractor> &extractor) const
{
  return FEValuesViews::RenumberedView<
    typename FEValuesViews::View<dim1, spacedim, Extractor>>(
    (*first_fe_values)[extractor.extractor], *first_renumbering_data);
}



template <int dim1, int dim2, int spacedim>
template <typename Extractor>
inline const FEValuesViews::RenumberedView<
  typename FEValuesViews::View<dim2, spacedim, Extractor>>
FECouplingValues<dim1, dim2, spacedim>::operator[](
  const FEValuesExtractors::SecondCoupling<Extractor> &extractor) const
{
  return FEValuesViews::RenumberedView<
    typename FEValuesViews::View<dim2, spacedim, Extractor>>(
    (*second_fe_values)[extractor.extractor], *second_renumbering_data);
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
