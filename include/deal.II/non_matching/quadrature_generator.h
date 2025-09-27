// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_non_matching_quadrature_generator_h
#define dealii_non_matching_quadrature_generator_h

#include <deal.II/base/config.h>

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_restriction.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/tria.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/non_matching/immersed_surface_quadrature.h>

#include <functional>
#include <optional>

DEAL_II_NAMESPACE_OPEN
namespace NonMatching
{
  namespace internal
  {
    namespace QuadratureGeneratorImplementation
    {
      template <int dim, int spacedim>
      class QGenerator;
    } // namespace QuadratureGeneratorImplementation


    namespace DiscreteQuadratureGeneratorImplementation
    {
      template <int dim>
      class CellWiseFunction;
    } // namespace DiscreteQuadratureGeneratorImplementation
  }   // namespace internal


  /**
   * Struct storing settings for the QuadratureGenerator class.
   */
  struct AdditionalQGeneratorData
  {
    /**
     * Constructor.
     */
    AdditionalQGeneratorData(const unsigned int max_box_splits          = 4,
                             const double lower_bound_implicit_function = 1e-11,
                             const double min_distance_between_roots    = 1e-12,
                             const double limit_to_be_definite          = 1e-11,
                             const double root_finder_tolerance         = 1e-12,
                             const unsigned int max_root_finder_splits  = 2,
                             bool               split_in_half           = true);

    /**
     * The number of times we are allowed to split the incoming box
     * and recurse on each child.
     */
    unsigned int max_box_splits;

    /**
     * For a level set function, $\psi$, the implicit function theorem states
     * that it is possible to write one of the coordinates $x_i$ as a function
     * of the others if
     *
     * $|\frac{\partial \psi}{\partial x_i}| > 0$.
     *
     * In practice, the bound we have for the expression in the left-hand side
     * may be near but not equal to zero due to roundoff errors.
     *
     * This constant is a safety margin, $C$, that states that the implicit
     * function theorem can be used when
     *
     * $|\frac{\partial \psi}{\partial x_i}| > C$
     *
     * Thus this constant must be non-negative.
     */
    double lower_bound_implicit_function;

    /**
     * If two roots are closer to each other than this distance they are
     * merged to one.
     */
    double min_distance_between_roots;

    /**
     * A constant, $C$, controlling when a level set function, $\psi$, is
     * considered positive or negative definite:
     *
     * $\psi(x) >  C \Rightarrow \text{Positive definite}$,
     * $\psi(x) < -C \Rightarrow \text{Negative definite}$.
     */
    double limit_to_be_definite;

    /**
     * Tolerance for convergence of the underlying root finder.
     */
    double root_finder_tolerance;

    /**
     * The number of times the underlying rootfinder is allowed to split
     * an interval, while trying to find multiple roots.
     */
    unsigned int max_root_finder_splits;

    /**
     * This determines how a box is split when this is necessary. If true, the
     * box is split in two, if set to false the box is split into its $2^{dim}$
     * children.
     */
    bool split_in_half;
  };



  /**
   * This class creates immersed quadrature rules over a BoundingBox,
   * $B \subset \mathbb{R}^{dim}$, when the domain is described by a level set
   * function, $\psi$.
   *
   * This class creates quadrature rules for the intersections between the box
   * and the three different regions defined by the level set function. That is,
   * it creates quadrature rules to integrate over the following regions
   * @f[
   * N = \{x \in B : \psi(x) < 0 \}, \\
   * P = \{x \in B : \psi(x) > 0 \}, \\
   * S = \{x \in B : \psi(x) = 0 \}.
   * @f]
   * @image html immersed_quadratures.svg
   *
   * When working with level set functions, the most common is to describe a
   * domain, $\Omega$, as
   * @f[
   * \Omega = \{ x \in \mathbb{R}^{dim} : \psi(x) < 0 \}.
   * @f]
   * Given this, we shall use the name convention that $N$ is the "inside"
   * region (i.e. inside $\Omega$), $P$ is the "outside" region and $S$ is
   * the "surface" region. The "inside" and "outside" quadratures will also be
   * referred to as the "bulk"-quadratures.
   *
   * The underlying algorithm use a 1-dimensional quadrature rule as base for
   * creating the immersed quadrature rules. Gauss-Legendre quadrature
   * (QGauss) is recommended. The constructor takes an hp::QCollection<1>.
   * One can select which 1d-quadrature in the collection should be used
   * through the set_1d_quadrature() function. The number of quadrature points
   * in the constructed quadratures will vary depending on the level set
   * function. More quadrature points will be created if the intersection is
   * "bad", for example, if the zero-contour has a high curvature compared to
   * the size of the box. However, if the number of points in the 1d quadrature
   * is $n$ the number of points will be proportional to $n^{dim}$ in the bulk
   * quadratures and to $n^{dim-1}$ in the surface quadrature. For example,
   * in the 2d-example in the above figure, there are 2 points in the
   * 1d-quadrature. If the 1d-quadrature is a Gauss-Legendre quadrature and the
   * grid has size $h$, the immersed quadratures typically give global errors
   * proportional to $h^{2n}$, both for the bulk and surface integrals. If the
   * 1d-quadrature has positive weights, the weights of the immersed quadratures
   * will also be positive.
   *
   * A detailed description of the underlying algorithm can be found in
   * @cite saye_2015. This implementation has some modifications
   * compared to the algorithm description in the paper. In particular, it
   * builds the three different types of quadratures (inside, outside and
   * surface) simultaneously. Further, the so-called "pruning" step is not yet
   * implemented.
   */
  template <int dim>
  class QuadratureGenerator
  {
  public:
    using AdditionalData = AdditionalQGeneratorData;

    /**
     * Constructor. Each Quadrature<1> in @p quadratures1d can be chosen as base
     * for generating the immersed quadrature rules.
     *
     * @note It is important that each 1d-quadrature rule in the
     * hp::QCollection does not contain the points 0 and 1.
     */
    QuadratureGenerator(
      const hp::QCollection<1> &quadratures1D,
      const AdditionalData     &additional_data = AdditionalData());

    /**
     * Clears the inside, outside and surface quadratures.
     */
    void
    clear_quadratures();

    /**
     * Construct immersed quadratures rules for the incoming level set
     * function over the BoundingBox.
     *
     * To get the constructed quadratures, use the functions
     * get_inside_quadrature(),
     * get_outside_quadrature(),
     * get_surface_quadrature().
     *
     * @note Both value, gradient and hessian need to be implemented on the
     * incoming function.
     */
    void
    generate(const Function<dim> &level_set, const BoundingBox<dim> &box);

    /**
     * Same as above but does not clear quadratures and appends it to the
     * existing quadrature instead.
     */
    void
    generate_append(const Function<dim>    &level_set,
                    const BoundingBox<dim> &box);

    /**
     * Return the quadrature rule for the region
     * $\{x \in B : \psi(x) < 0 \}$
     * created in the previous call to generate().
     * Here, $B$ is BoundingBox passed to generate().
     */
    const Quadrature<dim> &
    get_inside_quadrature() const;

    /**
     * Return the quadrature rule for the region
     * $\{x \in B : \psi(x) > 0 \}$
     * created in the previous call to generate().
     * Here, $B$ is BoundingBox passed to generate().
     */
    const Quadrature<dim> &
    get_outside_quadrature() const;

    /**
     * Return the quadrature rule for the region
     * $\{x \in B : \psi(x) = 0 \}$
     * created in the previous call to generate().
     * Here, $B$ is BoundingBox passed to generate().
     *
     * @note The normal at the quadrature points will be parallel to $\nabla \psi$.
     */
    const ImmersedSurfaceQuadrature<dim> &
    get_surface_quadrature() const;

    /**
     * Set which 1d-quadrature in the collection passed to the constructor
     * should be used to create the immersed quadratures.
     */
    void
    set_1D_quadrature(const unsigned int q_index);

  private:
    /**
     * QuadratureGenerator is mainly used to start up the recursive
     * algorithm. This is the object that actually generates the quadratures.
     */
    internal::QuadratureGeneratorImplementation::QGenerator<dim, dim>
      q_generator;
  };


  /**
   * This class creates immersed quadrature rules over a face, $F$, of a
   * BoundingBox, when the domain is described by a level set function, $\psi$.
   *
   * In the same way as in the QuadratureGenerator class, this class generates
   * quadrature rules to integrate over 3 different regions of the face,
   * $F \subset \mathbb{R}^{dim}$:
   * @f[
   * N = \{x \in F : \psi(x) < 0 \}, \\
   * P = \{x \in F : \psi(x) > 0 \}, \\
   * S = \{x \in F : \psi(x) = 0 \},
   * @f]
   * which are again referred to as the "inside", $N$, "outside", $P$,
   * and "surface" region, $S$. These type of quadrature rules are in general
   * needed in immersed discontinuous Galerkin methods.
   *
   * Under the hood, this class uses the QuadratureGenerator class to build
   * these rules. This is done by restricting the dim-dimensional level set
   * function to the face, thus creating a (dim-1)-dimensional level set
   * function, $\phi$. It then creates the (dim-1)-dimensional quadratures by
   * calling QuadratureGenerator with $\phi$. This means that what holds for the
   * QuadratureGenerator class in general also holds for this class. In
   * particular, if the 1d-quadrature that is used as base contains $n$ points,
   * the number of points will be proportional to $n^{dim-1}$ in the in the
   * inside/outside quadratures and to $n^{dim-2}$ in the surface quadrature.
   */
  template <int dim>
  class FaceQuadratureGenerator
  {
  public:
    using AdditionalData = AdditionalQGeneratorData;

    /**
     * Constructor. Each Quadrature<1> in @p quadratures1d can be chosen as base
     * for generating the immersed quadrature rules.
     *
     * @note It is important that each 1d-quadrature rule in the
     * hp::QCollection does not contain the points 0 and 1.
     */
    FaceQuadratureGenerator(
      const hp::QCollection<1> &quadratures1D,
      const AdditionalData     &additional_data = AdditionalData());

    /**
     * Clears the inside, outside and surface quadratures.
     */
    void
    clear_quadratures();

    /**
     * Construct immersed quadratures rules for the incoming level set
     * function on a given face of the BoundingBox.
     *
     * To get the constructed quadratures, use the functions
     * get_inside_quadrature(),
     * get_outside_quadrature(),
     * get_surface_quadrature().
     *
     * @note Both value, gradient and hessian need to be implemented on the
     * incoming function.
     */
    void
    generate(const Function<dim>    &level_set,
             const BoundingBox<dim> &box,
             const unsigned int      face_index);

    /**
     * Same as above but does not clear quadratures and appends it to the
     * existing quadrature instead.
     */
    void
    generate_append(const Function<dim>    &level_set,
                    const BoundingBox<dim> &box,
                    const unsigned int      face_index);

    /**
     * Return the quadrature rule for the region
     * $\{x \in F : \psi(x) < 0 \}$
     * created in the previous call to generate().
     * Here, $F$ is the face of the BoundingBox passed to generate().
     */
    const Quadrature<dim - 1> &
    get_inside_quadrature() const;

    /**
     * Return the quadrature rule for the region
     * $\{x \in F : \psi(x) > 0 \}$
     * created in the previous call to generate().
     * Here, $F$ is the face of the BoundingBox passed to generate().
     */
    const Quadrature<dim - 1> &
    get_outside_quadrature() const;

    /**
     * Return the quadrature rule for the region
     * $\{x \in F : \psi(x) = 0 \}$
     * created in the previous call to generate().
     * Here, $F$ is the face of the BoundingBox passed to generate().
     *
     * @note The normal at the quadrature points will be parallel to $\nabla \psi$.
     */
    const ImmersedSurfaceQuadrature<dim - 1, dim> &
    get_surface_quadrature() const;

    /**
     * Set which 1d-quadrature in the collection passed to the constructor
     * should be used to create the immersed quadratures.
     */
    void
    set_1D_quadrature(const unsigned int q_index);

  private:
    /**
     * Lower-dimensional quadrature generator used to build the quadratures over
     * the face.
     */
    QuadratureGenerator<dim - 1> quadrature_generator;

    /**
     * The same surface quadrature as created by the quadrature_generator,
     * but having dim-dimensional normals.
     */
    ImmersedSurfaceQuadrature<dim - 1, dim> surface_quadrature;
  };


  /**
   * Specialization of the FaceQuadratureGenerator class for the 1-dimensional
   * case.
   *
   * In 1d, a face is only a point. Thus to generate the immersed
   * quadrature rules we add a single 0-dimensional quadrature point to the
   * inside or outside quadrature rule depending on if the level set function is
   * positive or negative at the face. The added quadrature point will have
   * weight equal to 1. The immersed surface quadrature over a face corresponds
   * to integrating over a dim-1 dimensional curve. Thus, surface quadrature
   * generated by this specialized class is always empty.
   *
   * This class must be specialized in 1d, because the general
   * FaceQuadratureGenerator<dim> class uses the QuadratureGenerator<dim-1>
   * class internally, which does not make sense when dim-1 = 0.
   */
  template <>
  class FaceQuadratureGenerator<1>
  {
  public:
    using AdditionalData = AdditionalQGeneratorData;

    /**
     * Constructor. The incoming hp::QCollection is not used. But this class
     * must have the same signature as the non-specialized class.
     */
    FaceQuadratureGenerator(
      const hp::QCollection<1> &quadratures1D,
      const AdditionalData     &additional_data = AdditionalData());

    /**
     * Does nothing. Exists for compatibility reasons.
     */
    void
    clear_quadratures();

    /**
     * Construct immersed quadratures rules for the incoming level set
     * function on a given face of the BoundingBox.
     *
     * To get the constructed quadratures, use the functions
     * get_inside_quadrature(),
     * get_outside_quadrature(),
     * get_surface_quadrature().
     */
    void
    generate(const Function<1>    &level_set,
             const BoundingBox<1> &box,
             const unsigned int    face_index);

    /**
     * Same as above but does not clear quadratures and appends it to the
     * existing quadrature instead.
     */
    void
    generate_append(const Function<1>    &level_set,
                    const BoundingBox<1> &box,
                    const unsigned int    face_index);

    /**
     * @copydoc FaceQuadratureGenerator<dim>::get_inside_quadrature()
     */
    const Quadrature<0> &
    get_inside_quadrature() const;

    /**
     * @copydoc FaceQuadratureGenerator<dim>::get_outside_quadrature()
     */
    const Quadrature<0> &
    get_outside_quadrature() const;

    /**
     * Return the quadrature rule for the region
     * $\{x \in F : \psi(x) = 0 \}$
     * where, $F$ is the face of the BoundingBox passed to generate().
     *
     * @note In 1d, this quadrature always contains 0 points.
     */
    const ImmersedSurfaceQuadrature<0, 1> &
    get_surface_quadrature() const;

    /**
     * This function does nothing. It only exist to be compatible with
     * FaceQuadratureGenerator<dim>.
     */
    void
    set_1D_quadrature(const unsigned int q_index);

  private:
    /**
     * Quadrature for the region
     * $\{x \in F : \psi(x) < 0 \}$.
     * Created in the last call to generate().
     */
    Quadrature<0> inside_quadrature;

    /**
     * Quadrature for the region
     * $\{x \in F : \psi(x) > 0 \}$.
     * Created in the last call to generate().
     */
    Quadrature<0> outside_quadrature;

    /**
     * Quadrature for the region
     * $\{x \in F : \psi(x) = 0 \}$.
     * This quadrature always contains zero points in 1d.
     */
    const ImmersedSurfaceQuadrature<0, 1> surface_quadrature;
  };


  /**
   * This class generates the same type of immersed quadrature rules as those
   * described in the QuadratureGenerator class. The difference is that this
   * class handles the case when the domain is a discrete level set
   * function, i.e., when the level set function is described as a
   * (DoFHandler, Vector) pair. The generate()-function of this class takes a
   * cell in real space and constructs the immersed quadrature rules in
   * reference space over this cell. These quadrature rules can then be obtained
   * with one of the functions:
   * get_inside_quadrature(),
   * get_outside_quadrature(), and
   * get_surface_quadrature().
   *
   * Internally, the quadrature generation is done by transforming the discrete
   * level set function from real space to reference space and using the same
   * algorithm as in the QuadratureGenerator class.
   */
  template <int dim>
  class DiscreteQuadratureGenerator : public QuadratureGenerator<dim>
  {
  public:
    using AdditionalData = AdditionalQGeneratorData;

    /**
     * Constructor, the discrete level set function is described by the
     * incoming DoFHandler and Vector. Pointers to these are stored
     * internally, so they must have a longer lifetime than the created this
     * class. The hp::QCollection<1> and AdditionalData is passed to the
     * QuadratureGenerator class.
     */
    template <typename Number>
    DiscreteQuadratureGenerator(
      const hp::QCollection<1> &quadratures1D,
      const DoFHandler<dim>    &dof_handler,
      const ReadVector<Number> &level_set,
      const AdditionalData     &additional_data = AdditionalData());

    /**
     * Construct immersed quadratures rules based on the discrete level
     * set vector over the incoming cell.
     *
     * @note The cell needs to belong to the same triangulation as the one
     * associated with the DoFHandler passed to the constructor.
     */
    void
    generate(const typename Triangulation<dim>::active_cell_iterator &cell);

  private:
    /**
     * Construct immersed quadratures for FE_Q_iso_Q1.
     */
    void
    generate_fe_q_iso_q1(const BoundingBox<dim> &unit_box);

    /**
     * Function that describes our level set function in reference space.
     */
    std::unique_ptr<internal::DiscreteQuadratureGeneratorImplementation::
                      CellWiseFunction<dim>>
      reference_space_level_set;
  };

  /**
   * This class generates the same type of immersed quadrature rules as those
   * described in the FaceQuadratureGenerator class. The difference is that this
   * class handles the case when the domain is a discrete level set
   * function, i.e., when the level set function is described as a
   * (DoFHandler, Vector) pair. The generate()-function of this class takes a
   * cell in real space plus the respective face index and constructs the
   * immersed quadrature rules in reference space over this face. These
   * quadrature rules can then be obtained with one of the functions:
   * get_inside_quadrature(),
   * get_outside_quadrature(), and
   * get_surface_quadrature().
   *
   * Internally, the quadrature generation is done by transforming the discrete
   * level set function from real space to reference space and using the same
   * algorithm as in the FaceQuadratureGenerator class.
   */
  template <int dim>
  class DiscreteFaceQuadratureGenerator : public FaceQuadratureGenerator<dim>
  {
  public:
    using AdditionalData = AdditionalQGeneratorData;

    /**
     * Constructor, the discrete level set function is described by the
     * incoming DoFHandler and Vector. Pointers to these are stored
     * internally, so they must have a longer lifetime than the created this
     * class. The hp::QCollection<1> and AdditionalData is passed to the
     * QuadratureGenerator class.
     */
    template <typename Number>
    DiscreteFaceQuadratureGenerator(
      const hp::QCollection<1> &quadratures1D,
      const DoFHandler<dim>    &dof_handler,
      const ReadVector<Number> &level_set,
      const AdditionalData     &additional_data = AdditionalData());

    /**
     * Construct immersed quadratures rules based on the discrete level
     * set vector over the incoming face described by cell and face index.
     *
     * @note The cell needs to belong to the same triangulation as the one
     * associated with the DoFHandler passed to the constructor.
     */
    void
    generate(const typename Triangulation<dim>::active_cell_iterator &cell,
             const unsigned int face_index);

  private:
    /**
     * Construct immersed quadratures for FE_Q_iso_Q1.
     */
    void
    generate_fe_q_iso_q1(const BoundingBox<dim> &unit_box,
                         unsigned int            face_index);

    /**
     * Function that describes our level set function in reference space.
     */
    std::unique_ptr<internal::DiscreteQuadratureGeneratorImplementation::
                      CellWiseFunction<dim>>
      reference_space_level_set;
  };


  namespace internal
  {
    namespace QuadratureGeneratorImplementation
    {
      /**
       * A class that attempts to find multiple distinct roots of a function,
       * $f(x)$, over an interval, $[l, r]$. This is done as follows. If there
       * is a sign change in function value between the interval end points,
       * we solve for the root. If there is no sign change, we attempt to
       * bound the function value away from zero on $[a, b]$, to conclude that
       * no roots exist. If we can't exclude that there are roots, we split
       * the interval in two: $[l, (r + l) / 2]$, $[(r + l) / 2, r]$, and use
       * the same algorithm recursively on each interval. This means that we
       * can typically find 2 distinct roots, but not 3.
       *
       * The bounds on the functions values are estimated using the function
       * taylor_estimate_function_bounds, which approximates the function as a
       * second order Taylor-polynomial around the interval midpoint.
       * When we have a sign change on an interval, this class uses
       * boost::math::tools::toms748_solve for finding roots .
       */
      class RootFinder
      {
      public:
        /**
         * Struct storing settings for the RootFinder class.
         */
        struct AdditionalData
        {
          /**
           * Constructor.
           */
          AdditionalData(const double       tolerance           = 1e-12,
                         const unsigned int max_recursion_depth = 2,
                         const unsigned int max_iterations      = 500);

          /**
           * The tolerance in the stopping criteria for the underlying root
           * finding algorithm boost::math::tools::toms748_solve.
           */
          double tolerance;

          /**
           * The number of times we are allowed to split the interval where we
           * seek roots.
           */
          unsigned int max_recursion_depth;

          /**
           * The maximum number of iterations in
           * boost::math::tools::toms748_solve.
           */
          unsigned int max_iterations;
        };


        /**
         * Constructor.
         */
        RootFinder(const AdditionalData &data = AdditionalData());

        /**
         * For each of the incoming @p functions, attempt to find the roots over
         * the interval defined by @p interval and add these to @p roots.
         * The returned roots will be sorted in ascending order:
         * $x_0 < x_1 <...$ and duplicate roots (with respect to the tolerance
         * in AdditionalData) will be removed.
         */
        void
        find_roots(const std::vector<std::reference_wrapper<const Function<1>>>
                                        &functions,
                   const BoundingBox<1> &interval,
                   std::vector<double>  &roots);

      private:
        /**
         * Attempt to find the roots of the @p function over the interval defined by
         * @p interval and add these to @p roots. @p recursion_depth holds the number
         * of times this function has been called recursively.
         */
        void
        find_roots(const Function<1>    &function,
                   const BoundingBox<1> &interval,
                   const unsigned int    recursion_depth,
                   std::vector<double>  &roots);

        const AdditionalData additional_data;
      };


      /**
       * This is a special Quadrature class with a push_back() method for
       * conveniently adding a point with an associated weight.
       *
       * Since we build the quadrature rules in step-wise fashion,
       * it's easier to use this class than to pass around two vectors:
       * std::vector<Point<dim>>,
       * std::vector<double>.
       * Further, two std::vectors could accidentally end up with different
       * sizes. Using push_back we make sure that the number of points and
       * weights are the same.
       */
      template <int dim>
      class ExtendableQuadrature : public Quadrature<dim>
      {
      public:
        /**
         * Constructor, creates an empty quadrature rule with no
         * points.
         */
        ExtendableQuadrature() = default;

        /**
         * Constructor, copies the incoming Quadrature.
         */
        ExtendableQuadrature(const Quadrature<dim> &quadrature);

        /**
         * Clears weights and points vectors.
         */
        void
        clear();

        /**
         * Add a point with an associated weight to the quadrature.
         */
        void
        push_back(const Point<dim> &point, const double weight);
      };


      /**
       * Type that describes the definiteness of a function over a region.
       */
      enum class Definiteness
      {
        negative,
        positive,
        indefinite
      };


      /**
       * Class that stores quadrature rules to integrate over 4 different
       * regions of a single BoundingBox, $B$. Given multiple level set
       * functions,
       *
       * $\psi_i : \mathbb{R}^{dim} \rightarrow \mathbb{R}$, $i = 0, 1, ...$,
       *
       * the box, $B \subset \mathbb{R}^{dim}$, is partitioned into a
       * "negative", "positive", and "indefinite" region, $B = N \cup P \cup I$,
       * according to the signs of $\psi_i$ over each region:
       *
       * @f[
       * N = \{x \in B : \psi_i(x) < 0, \forall i \}, \\
       * P = \{x \in B : \psi_i(x) > 0, \forall i \}, \\
       * I = B \setminus (\overline{N} \cup \overline{P}).
       * @f]
       *
       * Thus, all $\psi_i$ are positive over $P$ and negative over $N$. Over
       * $I$ the level set functions differ in sign. This class holds quadrature
       * rules for each of these regions. In addition, when there is a single
       * level set function, $\psi$, it holds a surface quadrature for the zero
       * contour of $\psi$:
       *
       * $S = \{x \in B : \psi(x) = 0 \}$.
       *
       * Note that when there is a single level set function, $I$ is empty
       * and $N$ and $P$ are the regions that one typically integrates over in
       * an immersed finite element method.
       */
      template <int dim>
      class QPartitioning
      {
      public:
        /**
         * Return a reference to the "bulk" quadrature with the same name as the
         * member in Definiteness.
         */
        ExtendableQuadrature<dim> &
        quadrature_by_definiteness(const Definiteness definiteness);

        /**
         * Clears all quadratures.
         */
        void
        clear();

        /**
         * Quadrature for the region $\{x \in B : \psi_i(x) < 0 \forall i \}$ of
         * the box, $B$.
         */
        ExtendableQuadrature<dim> negative;

        /**
         * Quadrature for the region $\{x \in B : \psi_i(x) > 0 \forall i \}$ of
         * the box, $B$.
         */
        ExtendableQuadrature<dim> positive;

        /**
         * Quadrature for a region where the level set functions have different
         * sign.
         */
        ExtendableQuadrature<dim> indefinite;

        /**
         * Quadrature for the region $\{x \in B : \psi(x) = 0 \}$ of the
         * box, $B$.
         */
        ImmersedSurfaceQuadrature<dim> surface;
      };


      /**
       * This class is responsible for creating quadrature points for
       * the $dim$-dimensional quadrature partitioning from an
       * $(dim - 1)$-dimensional "indefinite" quadrature (see
       * QPartitioning documentation).
       *
       * To be precise, let $[L, R]$ be the extents of the box in the height
       * function direction and let $I \subset \mathbb{R}^{dim-1}$ be the lower
       * dimensional indefinite region. This class will create quadrature points
       * over $I \times [L, R] \subset \mathbb{R}^{dim}$ and in the case
       * $dim=spacedim$, points for the surface quadrature.
       *
       * For each lower dimensional quadrature point, $(x_I, w_I)$ in the
       * indefinite quadrature, we create several 1d-level set functions by
       * restricting $\psi_j$ to $x_I$. We then partition the interval $[L, R]$
       * into $[y_0, y_1, ..., y_n]$, where $y_0 = L$, $y_n = R$, and the
       * remaining $y_i$ are the roots of the 1d-level set functions in
       * $[L, R]$. Since the level set functions change sign between the
       * roots, each interval belong to different regions in the quadrature
       * partitioning.
       *
       * In each interval, $[y_i, y_{i+1}]$, we distribute points
       * according to the 1d-base quadrature, $(x_q, w_q)$ and take the
       * cartesian product with $(x_I, w_I)$ to create the $dim$-dimensional
       * quadrature points, $(X_q, W_q)$:
       * $X_q = x_I \times (y_i + (y_{i+1} - y_i) x_q)$,
       * $W_q = w_I (y_{i+1} - y_i) w_q$.
       *
       * When $dim=spacedim$, we have a single level set function, $\psi$. Since
       * we have fulfilled the implicit function theorem, there is a single root
       * $y_1 \in [L, R]$. The point, $x_s = x_I \times y_1$, will be added as a
       * point in the surface quadrature. One can show that the correct weight
       * of this point is
       *
       * $w_s = \frac{\|\nabla \psi(x_s)\|}{|\partial_i \psi(x_s)|} w_I$,
       *
       * where $i$ is the height function direction.
       */
      template <int dim, int spacedim>
      class UpThroughDimensionCreator
      {
      public:
        /**
         * Constructor. Takes the same parameters as QuadratureGenerator.
         */
        UpThroughDimensionCreator(
          const hp::QCollection<1>       &q_collection1D,
          const AdditionalQGeneratorData &additional_data);

        /**
         * Create $dim$-dimensional immersed quadratures from the incoming
         * $(dim-1)$-dimensional quadratures and add these to
         * @p q_partitioning.
         */
        void
        generate(const std::vector<std::reference_wrapper<const Function<dim>>>
                                           &level_sets,
                 const BoundingBox<dim>    &box,
                 const Quadrature<dim - 1> &low_dim_quadrature,
                 const unsigned int         height_function_direction,
                 QPartitioning<dim>        &q_partitioning);

        /**
         * Set which 1d-quadrature in the collection passed to the constructor
         * should be used to create the immersed quadratures.
         */
        void
        set_1D_quadrature(const unsigned int q_index);

      private:
        /**
         * Create a surface quadrature point from the lower-dimensional point
         * and add it to surface_quadrature.
         *
         * This function is only called when $dim=spacedim$ and there is a
         * single level set function. At this point there should only be a
         * single root in the interval $[L, R]$
         */
        void
        create_surface_point(
          const Point<dim - 1> &point,
          const double          weight,
          const std::vector<std::reference_wrapper<const Function<dim>>>
                                         &level_sets,
          const BoundingBox<dim>         &box,
          const unsigned int              height_function_direction,
          ImmersedSurfaceQuadrature<dim> &surface_quadrature);

        /**
         * One dimensional quadrature rules used to create the immersed
         * quadratures.
         */
        const ObserverPointer<const hp::QCollection<1>> q_collection1D;

        /**
         * Stores options/settings for the algorithm.
         */
        const AdditionalQGeneratorData additional_data;

        /**
         * Which quadrature rule in the above collection that is used to
         * create the immersed quadrature rules.
         */
        unsigned int q_index;

        /**
         * 1d-functions, that are restrictions of each dim-dimensional level set
         * function passed to generate() to some $(dim-1)$-dimensional point.
         */
        std::vector<Functions::PointRestriction<dim - 1>> point_restrictions;

        /**
         * Class used to find the roots of the above 1d-restrictions.
         */
        RootFinder root_finder;

        /**
         * The roots of the functions in point_restrictions.
         * This will be the values of the height functions, $\{H_i(x_I)\}$ at
         * some lower dimensional quadrature point,
         * $x_I \in \mathbb{R}^{dim-1}$.
         */
        std::vector<double> roots;
      };


      /**
       * Data representing the best choice of height-function direction,
       * which is returned by the function find_best_height_direction.
       *
       * This data consists of a coordinate direction
       *
       * $i \in \{0, ..., dim - 1 \}$,
       *
       * and lower bound on the absolute value of the derivative of some
       * associated function, f, taken in the above coordinate direction. That
       * is, a bound $C$ such that
       *
       * $|\frac{\partial f}{\partial x_i}| > C$,
       *
       * holding over some subset of $\mathbb{R}^{dim}$.
       */
      struct HeightDirectionData
      {
        /**
         * Constructor. Initializes the direction to invalid_unsigned_int and
         * the bound to 0.
         */
        HeightDirectionData();


        /**
         * The height-function direction, described above.
         */
        unsigned int direction;

        /**
         * The lower bound on $|\frac{\partial f}{\partial x_i}|$, described
         * above.
         */
        double min_abs_dfdx;
      };


      /**
       * Base class for the class QGenerator<dim, spacedim> and the
       * one-dimensional specialization QGenerator<1, spacedim>.
       */
      template <int dim, int spacedim>
      class QGeneratorBase
      {
      public:
        QGeneratorBase(const hp::QCollection<1>       &q_collection1D,
                       const AdditionalQGeneratorData &additional_data);

        /**
         * Clear the quadratures created by the previous call to generate().
         */
        void
        clear_quadratures();

        /**
         * Return the created quadratures.
         */
        const QPartitioning<dim> &
        get_quadratures() const;

      protected:
        /**
         * Stores options/settings for the algorithm.
         */
        const AdditionalQGeneratorData additional_data;

        /**
         * Which 1d-quadrature in the collection we should use to generate
         * the immersed quadrature.
         */
        unsigned int q_index;

        /**
         * Index of the quadrature in q_collection1d that should use to
         * generate the immersed quadrature rules.
         */
        const ObserverPointer<const hp::QCollection<1>> q_collection1D;

        /**
         * Quadratures that the derived classes create.
         */
        QPartitioning<dim> q_partitioning;
      };


      /**
       * This class implements the Saye-algorithm cited in the documentation of
       * the QuadratureGenerator class.
       *
       * The generate function takes a number of $dim$-dimensional level set
       * functions, $\psi_i$, and a BoundingBox<dim>, and builds a partitioning
       * of quadratures, as defined in documentation of the QPartitioning class.
       * That is, this class builds an object of type QPartitioning<dim>.
       *
       * If all $\psi_i$ passed to generate can be determined to be positive or
       * negative definite, the QPartitioning will consist of a single
       * quadrature forming a tensor product.
       *
       * If this is not the case, the algorithm uses recursion over the spatial
       * dimension. The spacedim template parameter denotes the dimension we
       * started with and dim denotes on what level we are in the recursion.
       * That is, we first construct a QPartitioning<dim - 1> and then
       * build the higher dimensional quadratures from these. What we in the end
       * actually want is a spacedim-dimensional partitioning of quadratures,
       * for a single level set function, $\psi$.
       *
       * The algorithm is based on the implicit function theorem. Starting with
       * a single level set function, $\psi$, we try to find a direction $i$,
       * such that
       *
       * $|\frac{\partial \psi}{\partial x_i}| > 0$.
       *
       * throughout the whole box. This means that the zero-contour of the
       * level set function can be parameterized by an implicit function
       *
       * $H = H(x_0, ..., x_{i-1}, x_{i+1}, ..., x_{dim-1})$,
       *
       * so that
       *
       * $\psi(..., x_{i-1}, H(..., x_{i-1}, x_{i+1}, ...), x_{i+1}, ...) = 0$,
       *
       * over a subset, $I \subset C \subset \mathbb{R}^{dim-1}$, of the cross
       * section, $C$, of the box (see BoundingBox::cross_section). Here, $I$ is
       * the "indefinite"-region defined in the QPartitioning class. To follow
       * convention in the original paper, we will -refer to $H$ as the
       * "height-function" and to $i$ as the "height-function direction".
       *
       * If a height function direction can be found, we go down in dimension by
       * creating two new level set functions, $\{\psi_0, \psi_1\}$, which are
       * the restriction of $\psi$ to the top and bottom faces of the box (in
       * the height function direction). We then delegate to
       * QGenerator<dim-1, spacedim> to create a QPartitioning<dim-1> over
       * the cross section.
       *
       * When we reach the base case, $dim = 1$, the creation of
       * QPartitioning<1> is simple. See the documentation in specialized
       * class: QGenerator<1, spacedim>.
       *
       * As we go up through the dimensions and create the higher dimensional
       * quadratures, we need to know the function value of the height
       * functions at the lower dimensional quadrature points. Since the
       * functions are implicit, we need to do root-finding on the level set
       * functions to find the function values. For this we use the class
       * UpThroughDimensionCreator, see documentation there.
       *
       * When we have $n$ level set functions (i.e. after having gone
       * down in dimension), we try to find a height function direction,
       * which works for all those $\psi_i$ which are intersected by the zero
       * contour (i.e. those not positive or negative definite).
       * If such a direction exist, we will have a maximum of $n$ associated
       * implicit height functions, $H_j$. Each $H_j$ parametrize the
       * $x_i$-coordinate of the zero-contour over a region, $I_j$. The
       * indefinite region in the lower dimensional partitioning is the union of
       * these $I = \cup_j I_j$.
       *
       * As we try to find a height function direction, we estimate bounds on
       * the gradient components by approximating each component as a 1st-order
       * Taylor-polynomial. If a direction can not be found, the box is split
       * and we recurse on each smaller box. This makes an implicit function
       * more likely to exist since we seek it over a smaller portion of the
       * zero contour. It also makes the estimated bounds tighter since we
       * extrapolate the Taylor-polynomial a shorter distance.
       *
       * Since we can not split a box forever, there is an maximum number of
       * allowed splits on the additional data struct passed to the constructor.
       * If this is reached, the algorithm uses the midpoint method as a last
       * resort.
       */
      template <int dim, int spacedim>
      class QGenerator : public QGeneratorBase<dim, spacedim>
      {
      public:
        /**
         * Constructor. Takes the same parameters QuadratureGenerator.
         */
        QGenerator(const hp::QCollection<1>       &q_collection1D,
                   const AdditionalQGeneratorData &additional_data);

        /**
         * Create immersed quadrature rules over the incoming @p box and add
         * these to the internal QPartitioning<dim> object in the base class.
         * These quadratures can then be obtained using the
         * get_quadratures-function.
         *
         * This function calls itself if the incoming box need to be split.
         * @p n_box_splits counts the number of times this function has called
         * itself.
         */
        void
        generate(const std::vector<std::reference_wrapper<const Function<dim>>>
                                        &level_sets,
                 const BoundingBox<dim> &box,
                 const unsigned int      n_box_splits);

        /**
         * Set which 1d-quadrature in the collection passed to the constructor
         * should be used to create the immersed quadratures.
         */
        void
        set_1D_quadrature(const unsigned int q_index);

      private:
        /**
         * Restricts the incoming level set functions to the top and bottom of
         * the incoming box (w.r.t @p height_function_direction). Then call the
         * lower dimensional QGenerator with the cross section of the box
         * to generate the lower dimensional immersed quadrature rules.
         */
        void
        create_low_dim_quadratures(
          const unsigned int height_function_direction,
          const std::vector<std::reference_wrapper<const Function<dim>>>
                                 &level_sets,
          const BoundingBox<dim> &box,
          const unsigned int      n_box_splits);

        /**
         * Gets the $(dim - 1)$-dimensional quadratures from the lower
         * dimensional algorithm and creates the $dim$-dimensional quadrature
         * rules over the box from the lower dimensional ones.
         */
        void
        create_high_dim_quadratures(
          const unsigned int height_function_direction,
          const std::vector<std::reference_wrapper<const Function<dim>>>
                                 &level_sets,
          const BoundingBox<dim> &box);

        /**
         * Split the incoming box and call generate() recursively with each box.
         * The box is split in 2 or 4 parts depending on the value of
         * AdditionalQGeneratorData::split_in_half.
         */
        void
        split_box_and_recurse(
          const std::vector<std::reference_wrapper<const Function<dim>>>
                                                   &level_sets,
          const BoundingBox<dim>                   &box,
          const std::optional<HeightDirectionData> &direction_data,
          const unsigned int                        n_box_splits);

        /**
         * Uses the midpoint-method to create a quadrature over the box.
         * That is, add a single quadrature point at the center of the box
         * with weight corresponding to the volume of the box.
         *
         * The point is added to the region defined in QPartitioning
         * according to the signs of the level set functions at the center of
         * the box.
         */
        void
        use_midpoint_method(
          const std::vector<std::reference_wrapper<const Function<dim>>>
                                 &level_sets,
          const BoundingBox<dim> &box);

        /**
         * The same algorithm as this, but creating immersed quadratures
         * in one dimension lower.
         */
        QGenerator<dim - 1, spacedim> low_dim_algorithm;

        /**
         * Object responsible for creating the $dim$-dimensional quadratures
         * from
         */
        UpThroughDimensionCreator<dim, spacedim> up_through_dimension_creator;

        /**
         * Stores tensor products of each of the Quadrature<1>'s in
         * q_collection1d.
         */
        hp::QCollection<dim> tensor_products;
      };


      /**
       * The 1d-base case of the recursive algorithm QGenerator<dim, spacedim>.
       *
       * Let $L$ and $R$ be the left and right bounds of the one-dimensional
       * BoundingBox. This interval is partitioned into $[x_0, x_1, ..., x_n]$
       * where $x_0 = L$, $x_n = R$, and the remaining $x_i$ are the roots
       * of the level set functions in the interval $[L, R]$. In each interval,
       * $[x_i, x_{i+1}]$, quadrature points are distributed according to a
       * 1d-quadrature rule. These points are added to one of the regions of
       * QPartitioning determined from the signs of the level set
       * functions on the interval (see documentation of QPartitioning).
       *
       * If spacedim = 1 the points $[x_1, x_{n-1}]$ are also added as surface
       * quadrature points to QPartitioning::surface.
       */
      template <int spacedim>
      class QGenerator<1, spacedim> : public QGeneratorBase<1, spacedim>
      {
      public:
        /**
         * Constructor. Takes the same parameters QuadratureGenerator.
         */
        QGenerator(const hp::QCollection<1>       &quadratures1D,
                   const AdditionalQGeneratorData &additional_data);

        /**
         * Creates quadrature points over the interval defined by the incoming
         * box and adds these quadrature points to the internally stored
         * QPartitioning. These quadratures can then be obtained using
         * the get_quadratures-function.
         */
        void
        generate(const std::vector<std::reference_wrapper<const Function<1>>>
                                      &level_sets,
                 const BoundingBox<1> &box,
                 const unsigned int    n_box_splits);

        /**
         * Set which 1d-quadrature in the collection passed to the constructor
         * should be used to create the immersed quadratures.
         */
        void
        set_1D_quadrature(const unsigned int q_index);

      private:
        /**
         * Adds the point defined by coordinate to the surface quadrature of
         * ImmersedQuadrature with unit weight.
         */
        void
        create_surface_points(
          const std::vector<std::reference_wrapper<const Function<1>>>
            &level_sets);

        /**
         * Class used to find the roots of the functions passed to generate().
         */
        RootFinder root_finder;

        /**
         * Roots of the functions passed to generate().
         */
        std::vector<double> roots;

        /**
         * This would be the height-function direction in higher dimensions,
         * but in 1d there is only one coordinate direction.
         */
        const unsigned int direction = 0;

        /**
         * To reuse the distribute_points_between_roots()-function
         * we need a zero-dimensional quadrature point with unit weight.
         */
        const Point<0> zero_dim_point;
        const double   unit_weight = 1;
      };


      /**
       * Take the tensor product between (point, weight) and @p quadrature1d
       * scaled over [start, end] and add the resulting dim-dimensional
       * quadrature points to @p quadrature.
       *
       * @p component_in_dim specifies which dim-dimensional coordinate
       * quadrature1d should be written to.
       */
      template <int dim>
      void
      tensor_point_with_1D_quadrature(const Point<dim - 1> &point,
                                      const double          weight,
                                      const Quadrature<1>  &quadrature1D,
                                      const double          start,
                                      const double          end,
                                      const unsigned int    component_in_dim,
                                      ExtendableQuadrature<dim> &quadrature);


      /**
       * Checks the sign of the incoming Functions at the incoming point and
       * returns Definiteness::positive/negative if all the functions are
       * positive/negative at the point, otherwise returns
       * Definiteness::indefinite.
       */
      template <int dim>
      Definiteness
      pointwise_definiteness(
        const std::vector<std::reference_wrapper<const Function<dim>>>
                         &functions,
        const Point<dim> &point);


      /**
       * A struct storing the bounds on the function value and bounds
       * on each component of the gradient.
       */
      template <int dim>
      struct FunctionBounds
      {
      public:
        /**
         * Lower and upper bounds on the functions value.
         */
        std::pair<double, double> value;

        /**
         * Lower and upper bounds on each component of the gradient.
         */
        std::array<std::pair<double, double>, dim> gradient;
      };


      /**
       * Returns the max/min bounds on the value, taken over all the entries
       * in the incoming vector of FunctionBounds. That is, given the incoming
       * function bounds, $[L_j, U_j]$, this function returns
       * $[L, U]$,
       * where $L = \min_{j} L_j$ and $U = \max_{j} U_j$.
       */
      template <int dim>
      std::pair<double, double>
      find_extreme_values(
        const std::vector<FunctionBounds<dim>> &all_function_bounds);


      /**
       * Finds the best choice of height function direction, given the
       * FunctionBounds for a number of functions $\{\psi_j\}_{j=0}^{n-1}$.
       * Here, "best" is meant in the sense of the implicit function theorem.
       *
       * Let $J_I$ be the index set of the indefinite functions:
       *
       * $J_I = \{0,..., n - 1\} \setminus \{ j : |\psi_j| > 0 \}$.
       *
       * This function converts the incoming bounds to a lower bound, $L_{ij}$,
       * on the absolute value of each component of the gradient:
       *
       * $|\partial_k \psi_j| > L_{jk}$.
       *
       * and then returns a coordinate direction, $i$, and a lower bound $L$,
       * such that
       *
       * @f[
       * i = \arg \max_{k} \min_{j \in J_I} L_{jk}, \\
       * L =      \max_{k} \min_{j \in J_I} L_{jk}.
       * @f]
       *
       * This means $i$ is a coordinate direction such that all functions
       * intersected by the zero contour (i.e. those belonging to $J_I$) fulfill
       *
       * $|\partial_i \psi_j| > L$.
       *
       * Note that the estimated lower bound, $L$, can be zero or negative. This
       * means that no suitable height function direction exists. If all of the
       * incoming functions are positive or negative definite the returned
       * std::optional is non-set.
       */
      template <int dim>
      std::optional<HeightDirectionData>
      find_best_height_direction(
        const std::vector<FunctionBounds<dim>> &all_function_bounds);

    } // namespace QuadratureGeneratorImplementation


    namespace DiscreteQuadratureGeneratorImplementation
    {
      /**
       * Interface for a scalar Function which has a
       * set_active_cell(..)-function. That is, a function which we in some way
       * need to associate with a given cell in order to evaluate.
       */
      template <int dim>
      class CellWiseFunction : public Function<dim>
      {
      public:
        /**
         * Destructor. Declared to make it virtual.
         */
        virtual ~CellWiseFunction() = default;

        /**
         * Set the cell that the function should be evaluated on.
         */
        virtual void
        set_active_cell(
          const typename Triangulation<dim>::active_cell_iterator &cell) = 0;

        /**
         * Set the dof values and the bounding box of the subcell the
         * function should be evaluated on. Relevant for FE_Q_iso_Q1.
         */
        virtual void
        set_subcell(const std::vector<unsigned int> &mask,
                    const BoundingBox<dim>          &subcell_box) = 0;

        /**
         * Returns flag indicating if the finite element is FE_Q_iso_Q1.
         */
        virtual bool
        is_fe_q_iso_q1() const = 0;

        /**
         * Number of subdivisions of the FE_Q_iso_Q1 element.
         */
        virtual unsigned int
        n_subdivisions() const = 0;
      };

    } // namespace DiscreteQuadratureGeneratorImplementation
  }   // namespace internal

} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif /* dealii_non_matching_quadrature_generator_h */
