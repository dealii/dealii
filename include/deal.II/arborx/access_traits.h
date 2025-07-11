// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_arborx_access_traits_h
#define dealii_arborx_access_traits_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_ARBORX
#  include <deal.II/base/bounding_box.h>

#  include <ArborX.hpp>

#  include <utility>


DEAL_II_NAMESPACE_OPEN

namespace ArborXWrappers
{
#  if ARBORX_VERSION_MAJOR < 2
  /**
   * Base class for Point-based predicates providing basic functionality for
   * derived classes, not supposed to be used on its own.
   */
  class PointPredicate
  {
  protected:
    /**
     * Constructor. @p points is a list of points used by the predicate.
     */
    template <int dim, typename Number>
    PointPredicate(const std::vector<dealii::Point<dim, Number>> &points);

    /**
     * The number of points stored in the structure.
     */
    std::size_t
    size() const;

    /**
     * Return the `i`th Point stored in the object.
     */
    const dealii::Point<3, float> &
    get(unsigned int i) const;

  private:
    std::vector<dealii::Point<3, float>> points;
  };



  /**
   * This class defines a predicate used by ArborXWrappers::BVH to determine
   * for given points which of the bounding boxes used to build the
   * ArborXWrappers::BVH intersect with them.
   * @note The class is not supposed to be used in a polymorphic context.
   */
  class PointIntersectPredicate : private PointPredicate
  {
  public:
    /**
     * Constructor. @p points is a list of points which we are interested in
     * knowing if they intersect ArborXWrappers::BVH bounding boxes.
     */
    template <int dim, typename Number>
    PointIntersectPredicate(
      const std::vector<dealii::Point<dim, Number>> &points);

    // We need these since we inherit privately to avoid polymorphic use.
    using PointPredicate::get;
    using PointPredicate::size;
  };



  /**
   * This class defines a predicate used by ArborXWrappers::BVH to determine
   * for given points which are the nearest bounding boxes/points among the ones
   * used to build the ArborXWrappers::BVH.
   * @note The class is not supposed to be used in a polymorphic context.
   */
  class PointNearestPredicate : private PointPredicate
  {
  public:
    /**
     * Constructor. @p points is a list of points for which we are interested in
     * the @p n_nearest_neighbors in the ArborXWrappers::BVH bounding
     * boxes/points.
     */
    template <int dim, typename Number>
    PointNearestPredicate(const std::vector<dealii::Point<dim, Number>> &points,
                          const unsigned int n_nearest_neighbors);

    /**
     * Return the number of nearest neighbors we are looking for.
     */
    unsigned int
    get_n_nearest_neighbors() const;

    // We need these since we inherit privately to avoid polymorphic use.
    using PointPredicate::get;
    using PointPredicate::size;

  private:
    unsigned int n_nearest_neighbors;
  };



  /**
   * Base class for BoundingBox predicates providing basic functionality for
   * derived classes, not supposed to be used on its own.
   */
  class BoundingBoxPredicate
  {
  protected:
    /**
     * Constructor. @p bounding_boxes is a list of bounding boxes used by the
     * predicate.
     */
    template <int dim, typename Number>
    BoundingBoxPredicate(
      const std::vector<dealii::BoundingBox<dim, Number>> &bounding_boxes);

    /**
     * The number of bounding boxes stored in the structure.
     */
    std::size_t
    size() const;

    /**
     * Return the `i`th BoundingBox stored in the object.
     */
    const dealii::BoundingBox<3, float> &
    get(unsigned int i) const;

  private:
    std::vector<dealii::BoundingBox<3, float>> bounding_boxes;
  };



  /**
   * This class is used by ArborXWrappers::BVH to determine for given bounding
   * boxes which of the bounding boxes used to build the ArborXWrappers::BVH
   * intersect with them.
   * @note The class is not supposed to be used in a polymorphic context.
   */
  class BoundingBoxIntersectPredicate : private BoundingBoxPredicate
  {
  public:
    /**
     * Constructor. @p bounding_boxes is a list of bounding boxes which we are interested in
     * knowing if they intersect ArborXWrappers::BVH bounding boxes.
     */
    template <int dim, typename Number>
    BoundingBoxIntersectPredicate(
      const std::vector<dealii::BoundingBox<dim, Number>> &bounding_boxes);

    // We need these since we inherit privately to avoid polymorphic use.
    using BoundingBoxPredicate::get;
    using BoundingBoxPredicate::size;
  };


  /**
   * This class is used by ArborXWrappers::BVH to determine for given bounding
   * boxes which are the nearest bounding boxes/points among the ones used to
   * build the ArborXWrappers::BVH.
   * @note The class is not supposed to be used in a polymorphic context.
   */
  class BoundingBoxNearestPredicate : private BoundingBoxPredicate
  {
  public:
    /**
     * Constructor. @p bounding_boxes is a list of bounding boxes for which are interested in
     * knowing the @p n_nearest_neighbors nearest bounding boxes used to build the
     * ArborXWrappers::BVH.
     */
    template <int dim, typename Number>
    BoundingBoxNearestPredicate(
      const std::vector<dealii::BoundingBox<dim, Number>> &bounding_boxes,
      const unsigned int                                   n_nearest_neighbors);

    /**
     * Return the number of nearest neighbors we are looking for.
     */
    unsigned int
    get_n_nearest_neighbors() const;

    // We need these since we inherit privately to avoid polymorphic use.
    using BoundingBoxPredicate::get;
    using BoundingBoxPredicate::size;

  private:
    unsigned int n_nearest_neighbors;
  };



  /**
   * Base class for Sphere-based predicates providing basic functionality for
   * derived classes; not supposed to be used on its own.
   */
  class SpherePredicate
  {
  protected:
    /**
     * Constructor. @p spheres is a list of spheres, i.e. a center and radius
     * pair, used by the predicate.
     */
    template <int dim, typename Number>
    SpherePredicate(
      const std::vector<std::pair<dealii::Point<dim, Number>, Number>>
        &spheres);

    /**
     * The number of spheres stored in the structure.
     */
    std::size_t
    size() const;

    /**
     * Return the `i`th sphere stored in the object.
     */
    const std::pair<dealii::Point<3, float>, float> &
    get(unsigned int) const;

  private:
    std::vector<std::pair<dealii::Point<3, float>, float>> spheres;
  };



  /**
   * This class defines a predicate used by ArborXWrappers::BVH to determine
   * for given spheres which of the bounding boxes used to build the
   * ArborXWrappers::BVH intersect with them.
   * @note The class is not supposed to be used in a polymorphic context.
   */
  class SphereIntersectPredicate : private SpherePredicate
  {
  public:
    /**
     * Constructor. @p spheres is a list of spheres which we are interested in
     * knowing if they intersect ArborXWrappers::BVH bounding boxes.
     */
    template <int dim, typename Number>
    SphereIntersectPredicate(
      const std::vector<std::pair<dealii::Point<dim, Number>, Number>>
        &spheres);

    // We need these since we inherit privately to avoid polymorphic use.
    using SpherePredicate::get;
    using SpherePredicate::size;
  };



  /**
   * This class defines a predicate used by ArborXWrappers::BVH to determine,
   * for the given spheres, which are the nearest bounding boxes/points among
   * the ones used to build the ArborXWrappers::BVH.
   * @note The class is not supposed to be used in a polymorphic context.
   */
  class SphereNearestPredicate : private SpherePredicate
  {
  public:
    /**
     * Constructor. @p spheres is a list of spheres for which we are
     * interested in the @p n_nearest_neighbors in the ArborXWrappers::BVH
     * bounding boxes/points.
     */
    template <int dim, typename Number>
    SphereNearestPredicate(
      const std::vector<std::pair<dealii::Point<dim, Number>, Number>> &spheres,
      const unsigned int n_nearest_neighbors);

    /**
     * Return the number of nearest neighbors we are looking for.
     */
    unsigned int
    get_n_nearest_neighbors() const;

    // We need these since we inherit privately to avoid polymorphic use.
    using SpherePredicate::get;
    using SpherePredicate::size;

  private:
    unsigned int n_nearest_neighbors;
  };
#  else
  namespace internal
  {
    template <int dim, typename Number>
    ArborX::Point<dim, Number>
    to_arborx_point(const dealii::Point<dim, Number> &p)
    {
      if constexpr (dim == 1)
        {
          return {p[0]};
        }

      if constexpr (dim == 2)
        {
          return {p[0], p[1]};
        }

      if constexpr (dim == 3)
        {
          return {p[0], p[1], p[2]};
        }
    }

    /**
     * This structure returns the ArborX object associated with the deal.II
     * object stored in a PairValueIndex.
     */
    struct IndexableGetter
    {
      template <int dim, typename Number>
      ArborX::Point<dim, Number>
      operator()(const ArborX::PairValueIndex<dealii::Point<dim, Number>,
                                              unsigned int> &pair) const
      {
        return to_arborx_point(pair.value);
      }



      template <int dim, typename Number>
      ArborX::Box<dim, Number>
      operator()(const ArborX::PairValueIndex<dealii::BoundingBox<dim, Number>,
                                              unsigned int> &pair) const
      {
        const auto boundary_points = pair.value.get_boundary_points();
        const dealii::Point<dim, Number> min_corner = boundary_points.first;
        const dealii::Point<dim, Number> max_corner = boundary_points.second;
        return {to_arborx_point(min_corner), to_arborx_point(max_corner)};
      }



      template <int dim, typename Number>
      ArborX::Sphere<dim, Number>
      operator()(const ArborX::PairValueIndex<
                 std::pair<dealii::Point<dim, Number>, Number>,
                 unsigned int> &pair) const
      {
        return {to_arborx_point(pair.value.first), pair.value.second};
      }
    };



    /**
     * Callback to extract the index of each primitive that satisfies a query.
     */
    struct ExtractIndex
    {
      template <typename Query, typename Value, typename Output>
      KOKKOS_FUNCTION void
      operator()(const Query &, const Value &value, const Output &out) const
      {
        out(value.index);
      }
    };



    /**
     * Callback to extract the index and the rank of each primitive that
     * satisfies a query.
     */
    struct ExtractIndexRank
    {
      unsigned int rank;

      template <typename Predicate, typename Value, typename Output>
      KOKKOS_FUNCTION void
      operator()(const Predicate &,
                 const ArborX::PairValueIndex<Value> &value,
                 const Output                        &out) const
      {
        out({value.index, rank});
      }
    };
  } // namespace internal



  /**
   * This class defines a predicate used by ArborXWrappers::BVH to determine
   * for given points which of the primitives used to build the
   * ArborXWrappers::BVH intersect with them.
   */
  template <int dim, typename Number>
  class PointIntersectPredicate
  {
  public:
    /**
     * Constructor. @p points is a list of points which we are interested in
     * knowing if they intersect ArborXWrappers::BVH primitives.
     */
    PointIntersectPredicate(
      const std::vector<dealii::Point<dim, Number>> &points);

    /**
     * The number of points stored in the structure.
     */
    std::size_t
    size() const;

    /**
     * Return the `i`th Point stored in the object.
     */
    const dealii::Point<dim, Number> &
    get(unsigned int i) const;

    /**
     * A flag that specifies if the predicate is nearest neighbors search.
     */
    static constexpr bool is_nearest = false;

  private:
    std::vector<dealii::Point<dim, Number>> points;
  };



  /**
   * This class defines a predicate used by ArborXWrappers::BVH to determine
   * for given points which are the nearest primitives among the ones
   * used to build the ArborXWrappers::BVH.
   */
  template <int dim, typename Number>
  class PointNearestPredicate
  {
  public:
    /**
     * Constructor. @p points is a list of points for which we are interested in
     * the @p n_nearest_neighbors in the ArborXWrappers::BVH primitives.
     */
    PointNearestPredicate(const std::vector<dealii::Point<dim, Number>> &points,
                          const unsigned int n_nearest_neighbors);

    /**
     * Return the number of nearest neighbors we are looking for.
     */
    unsigned int
    get_n_nearest_neighbors() const;

    /**
     * The number of points stored in the structure.
     */
    std::size_t
    size() const;

    /**
     * Return the `i`th Point stored in the object.
     */
    const dealii::Point<dim, Number> &
    get(unsigned int i) const;


    /**
     * A flag that specifies if the predicate is nearest neighbors search.
     */
    static constexpr bool is_nearest = true;

  private:
    std::vector<dealii::Point<dim, Number>> points;
    unsigned int                            n_nearest_neighbors;
  };



  /**
   * This class is used by ArborXWrappers::BVH to determine for given bounding
   * boxes which of the primitives used to build the ArborXWrappers::BVH
   * intersect with them.
   */
  template <int dim, typename Number>
  class BoundingBoxIntersectPredicate
  {
  public:
    /**
     * Constructor. @p bounding_boxes is a list of bounding boxes which we are interested in
     * knowing if they intersect ArborXWrappers::BVH primitives.
     */
    BoundingBoxIntersectPredicate(
      const std::vector<dealii::BoundingBox<dim, Number>> &bounding_boxes);

    /**
     * The number of points stored in the structure.
     */
    std::size_t
    size() const;

    /**
     * Return the `i`th Point stored in the object.
     */
    const dealii::BoundingBox<dim, Number> &
    get(unsigned int i) const;

    /**
     * A flag that specifies if the predicate is nearest neighbors search.
     */
    static constexpr bool is_nearest = false;

  private:
    std::vector<dealii::BoundingBox<dim, Number>> bounding_boxes;
  };



  /**
   * This class is used by ArborXWrappers::BVH to determine for given bounding
   * boxes which are the nearest primitives among the ones used to
   * build the ArborXWrappers::BVH.
   */
  template <int dim, typename Number>
  class BoundingBoxNearestPredicate
  {
  public:
    /**
     * Constructor. @p bounding_boxes is a list of bounding boxes for which are interested in
     * knowing the @p n_nearest_neighbors nearest primitives used to build the
     * ArborXWrappers::BVH.
     */
    BoundingBoxNearestPredicate(
      const std::vector<dealii::BoundingBox<dim, Number>> &bounding_boxes,
      const unsigned int                                   n_nearest_neighbors);

    /**
     * Return the number of nearest neighbors we are looking for.
     */
    unsigned int
    get_n_nearest_neighbors() const;

    /**
     * The number of points stored in the structure.
     */
    std::size_t
    size() const;

    /**
     * Return the `i`th Point stored in the object.
     */
    const dealii::BoundingBox<dim, Number> &
    get(unsigned int i) const;

    /**
     * A flag that specifies if the predicate is nearest neighbors search.
     */
    static constexpr bool is_nearest = true;

  private:
    std::vector<dealii::BoundingBox<dim, Number>> bounding_boxes;
    unsigned int                                  n_nearest_neighbors;
  };



  /**
   * This class defines a predicate used by ArborXWrappers::BVH to determine
   * for given spheres which of the primitives used to build the
   * ArborXWrappers::BVH intersect with them.
   */
  template <int dim, typename Number>
  class SphereIntersectPredicate
  {
  public:
    /**
     * Constructor. @p spheres is a list of spheres which we are interested in
     * knowing if they intersect ArborXWrappers::BVH primitives.
     */
    SphereIntersectPredicate(
      const std::vector<std::pair<dealii::Point<dim, Number>, Number>>
        &spheres);

    /**
     * The number of points stored in the structure.
     */
    std::size_t
    size() const;

    /**
     * Return the `i`th Point stored in the object.
     */
    const std::pair<dealii::Point<dim, Number>, Number> &
    get(unsigned int) const;

    /**
     * A flag that specifies if the predicate is nearest neighbors search.
     */
    static constexpr bool is_nearest = false;

  private:
    std::vector<std::pair<dealii::Point<dim, Number>, Number>> spheres;
  };



  /**
   * This class defines a predicate used by ArborXWrappers::BVH to determine,
   * for the given spheres, which are the nearest bounding primitives among
   * the ones used to build the ArborXWrappers::BVH.
   */
  template <int dim, typename Number>
  class SphereNearestPredicate
  {
  public:
    /**
     * Constructor. @p spheres is a list of spheres for which we are
     * interested in the @p n_nearest_neighbors in the ArborXWrappers::BVH
     * primitives.
     */
    SphereNearestPredicate(
      const std::vector<std::pair<dealii::Point<dim, Number>, Number>> &spheres,
      const unsigned int n_nearest_neighbors);

    /**
     * Return the number of nearest neighbors we are looking for.
     */
    unsigned int
    get_n_nearest_neighbors() const;

    /**
     * The number of points stored in the structure.
     */
    std::size_t
    size() const;

    /**
     * Return the `i`th Point stored in the object.
     */
    const std::pair<dealii::Point<dim, Number>, Number> &
    get(unsigned int) const;

    /**
     * A flag that specifies if the predicate is nearest neighbors search.
     */
    static constexpr bool is_nearest = true;

  private:
    std::vector<std::pair<dealii::Point<dim, Number>, Number>> spheres;
    unsigned int n_nearest_neighbors;
  };
#  endif
} // namespace ArborXWrappers

DEAL_II_NAMESPACE_CLOSE

/**
 * This namespace contains the implementation of AccessTraits used by ArborX.
 */
namespace ArborX
{
#  if ARBORX_VERSION_MAJOR < 2
  /**
   * This struct allows ArborX to use std::vector<dealii::Point> as
   * primitive.
   */
  template <int dim, typename Number>
  struct AccessTraits<std::vector<dealii::Point<dim, Number>>, PrimitivesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * Return the size of the vector @p v.
     */
    static std::size_t
    size(const std::vector<dealii::Point<dim, Number>> &v);

    /**
     * Return an ArborX::Point from the dealii::Point `v[i]`.
     */
    static Point
    get(const std::vector<dealii::Point<dim, Number>> &v, std::size_t i);
  };



  /**
   * This struct allows ArborX to use std::vector<dealii::BoundingBox> as
   * primitive.
   */
  template <int dim, typename Number>
  struct AccessTraits<std::vector<dealii::BoundingBox<dim, Number>>,
                      PrimitivesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * Return the size of the vector @p v.
     */
    static std::size_t
    size(const std::vector<dealii::BoundingBox<dim, Number>> &v);

    /**
     * Return an ArborX::Box from the dealii::BoundingBox `v[i]`.
     */
    static Box
    get(const std::vector<dealii::BoundingBox<dim, Number>> &v, std::size_t i);
  };



  /**
   * This struct allows ArborX to use
   * std::vector<std::pair<dealii::Point, Number> as primitive.
   */
  template <int dim, typename Number>
  struct AccessTraits<
    std::vector<std::pair<dealii::Point<dim, Number>, Number>>,
    PrimitivesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * Return the size of the vector @p v.
     */
    static std::size_t
    size(const std::vector<std::pair<dealii::Point<dim, Number>, Number>> &v);

    /**
     * Return an ArborX::Sphere from the std::pair<dealii::Point, Number>
     * `v[i]`.
     */
    static Sphere
    get(const std::vector<std::pair<dealii::Point<dim, Number>, Number>> &v,
        std::size_t                                                       i);
  };
#  else
  /**
   * This struct allows ArborX to use std::vector<T> as primitive.
   */
  template <typename T>
  struct AccessTraits<
    std::vector<T>,
    std::enable_if_t<
      std::is_same_v<T, dealii::Point<T::dimension, float>> ||
      std::is_same_v<T, dealii::Point<T::dimension, double>> ||
      std::is_same_v<T, dealii::BoundingBox<T::dimension, float>> ||
      std::is_same_v<T, dealii::BoundingBox<T::dimension, double>> ||
      std::is_same_v<T, std::pair<dealii::Point<T::dimension, float>, float>> ||
      std::is_same_v<T,
                     std::pair<dealii::Point<T::dimension, double>, double>>>>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * Return the size of the vector @p v.
     */
    static std::size_t
    size(const std::vector<T> &v);

    /**
     * Return the ith element of the vector @p v.
     */
    static T
    get(const std::vector<T> &v, std::size_t i);
  };
#  endif



#  if ARBORX_VERSION_MAJOR < 2
  /**
   * This struct allows ArborX to use PointIntersectPredicate as a predicate.
   */
  template <>
  struct AccessTraits<dealii::ArborXWrappers::PointIntersectPredicate,
                      PredicatesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * The number of points stored in @p pt_intersect.
     */
    static std::size_t
    size(const dealii::ArborXWrappers::PointIntersectPredicate &pt_intersect);

    /**
     * Return an ArborX::intersects(ArborX::Point) object constructed from the
     * `i`th dealii::Point stored in @p pt_intersect.
     */
    static auto
    get(const dealii::ArborXWrappers::PointIntersectPredicate &pt_intersect,
        std::size_t                                            i);
  };



  /**
   * This struct allows ArborX to use PointNearestPredicate as a predicate.
   */
  template <>
  struct AccessTraits<dealii::ArborXWrappers::PointNearestPredicate,
                      PredicatesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * The number of points stored in @p pt_nearest.
     */
    static std::size_t
    size(const dealii::ArborXWrappers::PointNearestPredicate &pt_nearest);

    /**
     * Return an ArborX::nearest(ArborX::Point,
     * PointNearestPredicate::get_n_nearest_neighbors) object constructed from
     * the `i`th dealii::Point stored in @p pt_nearest.
     */
    static auto
    get(const dealii::ArborXWrappers::PointNearestPredicate &pt_nearest,
        std::size_t                                          i);
  };



  /**
   * This struct allows ArborX to use BoundingBoxIntersectPredicate as a
   * predicate.
   */
  template <>
  struct AccessTraits<dealii::ArborXWrappers::BoundingBoxIntersectPredicate,
                      PredicatesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * The number of bounding boxes stored in @p bb_intersect.
     */
    static std::size_t
    size(const dealii::ArborXWrappers::BoundingBoxIntersectPredicate
           &bb_intersect);

    /**
     * Return an Arbox::intersects(ArborX::Box) object constructed from the
     * `i`th dealii::BoundingBox stored in @p bb_intersect.
     */
    static auto
    get(
      const dealii::ArborXWrappers::BoundingBoxIntersectPredicate &bb_intersect,
      std::size_t                                                  i);
  };



  /**
   * This struct allows ArborX to use BoundingBoxNearstPredicate as a
   * predicate.
   */
  template <>
  struct AccessTraits<dealii::ArborXWrappers::BoundingBoxNearestPredicate,
                      PredicatesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * The number of bounding boxes stored in @p bb_nearest.
     */
    static std::size_t
    size(const dealii::ArborXWrappers::BoundingBoxNearestPredicate &bb_nearest);

    /**
     * Return an Arbox::nearest(ArborX::Box,
     * BoundingBoxtNearestPredicate::get_n_nearest_neighbors) object constructed
     * from the `i`th dealii::BoundingBox stored in @p bb_nearest.
     */
    static auto
    get(const dealii::ArborXWrappers::BoundingBoxNearestPredicate &bb_nearest,
        std::size_t                                                i);
  };



  /**
   * This struct allows ArborX to use SphereIntersectPredicate as a
   predicate.
   */
  template <>
  struct AccessTraits<dealii::ArborXWrappers::SphereIntersectPredicate,
                      PredicatesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * The number of points stored in @p sph_intersect.
     */
    static std::size_t
    size(const dealii::ArborXWrappers::SphereIntersectPredicate &sph_intersect);

    /**
     * Return an ArborX::intersects(ArborX::Sphere) object constructed from the
     * `i`th sphere stored in @p sph_intersect.
     */
    static auto
    get(const dealii::ArborXWrappers::SphereIntersectPredicate &sph_intersect,
        std::size_t                                             i);
  };



  /**
   * This struct allows ArborX to use SphereNearestPredicate as a predicate.
   */
  template <>
  struct AccessTraits<dealii::ArborXWrappers::SphereNearestPredicate,
                      PredicatesTag>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * The number of spheres stored in @p sph_nearest.
     */
    static std::size_t
    size(const dealii::ArborXWrappers::SphereNearestPredicate &sph_nearest);

    /**
     * Return an ArborX::nearest(ArborX::Sphere,
     * SphereNearestPredicate::get_n_nearest_neighbors) object constructed from
     * the `i`th sphere stored in @p sph_nearest.
     */
    static auto
    get(const dealii::ArborXWrappers::SphereNearestPredicate &sph_nearest,
        std::size_t                                           i);
  };
#  else
  /**
   * This struct allows ArborX to use PointIntersectPredicate as a predicate.
   */
  template <int dim, typename Number>
  struct AccessTraits<
    dealii::ArborXWrappers::PointIntersectPredicate<dim, Number>>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * The number of points stored in @p pt_intersect.
     */
    static std::size_t
    size(const dealii::ArborXWrappers::PointIntersectPredicate<dim, Number>
           &pt_intersect);

    /**
     * Return an ArborX::intersects(ArborX::Point) object constructed from the
     * `i`th dealii::Point stored in @p pt_intersect.
     */
    static auto
    get(const dealii::ArborXWrappers::PointIntersectPredicate<dim, Number>
                   &pt_intersect,
        std::size_t i);
  };



  /**
   * This struct allows ArborX to use PointNearestPredicate as a predicate.
   */
  template <int dim, typename Number>
  struct AccessTraits<
    dealii::ArborXWrappers::PointNearestPredicate<dim, Number>>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * The number of points stored in @p pt_nearest.
     */
    static std::size_t
    size(const dealii::ArborXWrappers::PointNearestPredicate<dim, Number>
           &pt_nearest);

    /**
     * Return an ArborX::nearest(ArborX::Point,
     * PointNearestPredicate::get_n_nearest_neighbors) object constructed from
     * the `i`th dealii::Point stored in @p pt_nearest.
     */
    static auto
    get(const dealii::ArborXWrappers::PointNearestPredicate<dim, Number>
                   &pt_nearest,
        std::size_t i);
  };



  /**
   * This struct allows ArborX to use BoundingBoxIntersectPredicate as a
   * predicate.
   */
  template <int dim, typename Number>
  struct AccessTraits<
    dealii::ArborXWrappers::BoundingBoxIntersectPredicate<dim, Number>>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * The number of bounding boxes stored in @p bb_intersect.
     */
    static std::size_t
    size(
      const dealii::ArborXWrappers::BoundingBoxIntersectPredicate<dim, Number>
        &bb_intersect);

    /**
     * Return an Arbox::intersects(ArborX::Box) object constructed from the
     * `i`th dealii::BoundingBox stored in @p bb_intersect.
     */
    static auto
    get(const dealii::ArborXWrappers::BoundingBoxIntersectPredicate<dim, Number>
                   &bb_intersect,
        std::size_t i);
  };



  /**
   * This struct allows ArborX to use BoundingBoxNearstPredicate as a
   * predicate.
   */
  template <int dim, typename Number>
  struct AccessTraits<
    dealii::ArborXWrappers::BoundingBoxNearestPredicate<dim, Number>>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * The number of bounding boxes stored in @p bb_nearest.
     */
    static std::size_t
    size(const dealii::ArborXWrappers::BoundingBoxNearestPredicate<dim, Number>
           &bb_nearest);

    /**
     * Return an Arbox::nearest(ArborX::Box,
     * BoundingBoxtNearestPredicate::get_n_nearest_neighbors) object constructed
     * from the `i`th dealii::BoundingBox stored in @p bb_nearest.
     */
    static auto
    get(const dealii::ArborXWrappers::BoundingBoxNearestPredicate<dim, Number>
                   &bb_nearest,
        std::size_t i);
  };



  /**
   * This struct allows ArborX to use SphereIntersectPredicate as a
   predicate.
   */
  template <int dim, typename Number>
  struct AccessTraits<
    dealii::ArborXWrappers::SphereIntersectPredicate<dim, Number>>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * The number of points stored in @p sph_intersect.
     */
    static std::size_t
    size(const dealii::ArborXWrappers::SphereIntersectPredicate<dim, Number>
           &sph_intersect);

    /**
     * Return an ArborX::intersects(ArborX::Sphere) object constructed from the
     * `i`th sphere stored in @p sph_intersect.
     */
    static auto
    get(const dealii::ArborXWrappers::SphereIntersectPredicate<dim, Number>
                   &sph_intersect,
        std::size_t i);
  };



  /**
   * This struct allows ArborX to use SphereNearestPredicate as a predicate.
   */
  template <int dim, typename Number>
  struct AccessTraits<
    dealii::ArborXWrappers::SphereNearestPredicate<dim, Number>>
  {
    using memory_space = Kokkos::HostSpace;

    /**
     * The number of spheres stored in @p sph_nearest.
     */
    static std::size_t
    size(const dealii::ArborXWrappers::SphereNearestPredicate<dim, Number>
           &sph_nearest);

    /**
     * Return an ArborX::nearest(ArborX::Sphere,
     * SphereNearestPredicate::get_n_nearest_neighbors) object constructed from
     * the `i`th sphere stored in @p sph_nearest.
     */
    static auto
    get(const dealii::ArborXWrappers::SphereNearestPredicate<dim, Number>
                   &sph_nearest,
        std::size_t i);
  };
#  endif

  // ------------------------------- Inline ----------------------------------//

  // The implementation of AccessTraits<..., PredicatesTag> needs to be in the
  // header file otherwise the return type of auto get() cannot be determined.
  // We use auto because ArborX does not expose the type of intersects

#  if ARBORX_VERSION_MAJOR < 2
  inline std::size_t
  AccessTraits<dealii::ArborXWrappers::PointIntersectPredicate, PredicatesTag>::
    size(const dealii::ArborXWrappers::PointIntersectPredicate &pt_intersect)
  {
    return pt_intersect.size();
  }



  inline auto
  AccessTraits<dealii::ArborXWrappers::PointIntersectPredicate, PredicatesTag>::
    get(const dealii::ArborXWrappers::PointIntersectPredicate &pt_intersect,
        std::size_t                                            i)
  {
    const auto dealii_point = pt_intersect.get(i);
    return intersects(Point{dealii_point[0], dealii_point[1], dealii_point[2]});
  }



  inline std::size_t
  AccessTraits<dealii::ArborXWrappers::PointNearestPredicate, PredicatesTag>::
    size(const dealii::ArborXWrappers::PointNearestPredicate &pt_nearest)
  {
    return pt_nearest.size();
  }



  inline auto
  AccessTraits<dealii::ArborXWrappers::PointNearestPredicate, PredicatesTag>::
    get(const dealii::ArborXWrappers::PointNearestPredicate &pt_nearest,
        std::size_t                                          i)
  {
    const auto dealii_point = pt_nearest.get(i);
    return nearest(Point{dealii_point[0], dealii_point[1], dealii_point[2]},
                   pt_nearest.get_n_nearest_neighbors());
  }



  inline std::size_t
  AccessTraits<dealii::ArborXWrappers::BoundingBoxIntersectPredicate,
               PredicatesTag>::
    size(
      const dealii::ArborXWrappers::BoundingBoxIntersectPredicate &bb_intersect)
  {
    return bb_intersect.size();
  }



  inline auto
  AccessTraits<dealii::ArborXWrappers::BoundingBoxIntersectPredicate,
               PredicatesTag>::
    get(
      const dealii::ArborXWrappers::BoundingBoxIntersectPredicate &bb_intersect,
      std::size_t                                                  i)
  {
    const auto boundary_points = bb_intersect.get(i).get_boundary_points();
    const dealii::Point<3, float> min_corner = boundary_points.first;
    const dealii::Point<3, float> max_corner = boundary_points.second;

    return intersects(Box{{min_corner[0], min_corner[1], min_corner[2]},
                          {max_corner[0], max_corner[1], max_corner[2]}});
  }



  inline std::size_t
  AccessTraits<dealii::ArborXWrappers::BoundingBoxNearestPredicate,
               PredicatesTag>::
    size(const dealii::ArborXWrappers::BoundingBoxNearestPredicate &bb_nearest)
  {
    return bb_nearest.size();
  }



  inline auto
  AccessTraits<dealii::ArborXWrappers::BoundingBoxNearestPredicate,
               PredicatesTag>::
    get(const dealii::ArborXWrappers::BoundingBoxNearestPredicate &bb_nearest,
        std::size_t                                                i)
  {
    const auto boundary_points = bb_nearest.get(i).get_boundary_points();
    const dealii::Point<3, float> min_corner = boundary_points.first;
    const dealii::Point<3, float> max_corner = boundary_points.second;

    return nearest(Box{{min_corner[0], min_corner[1], min_corner[2]},
                       {max_corner[0], max_corner[1], max_corner[2]}},
                   bb_nearest.get_n_nearest_neighbors());
  }



  inline std::size_t
  AccessTraits<dealii::ArborXWrappers::SphereIntersectPredicate,
               PredicatesTag>::
    size(const dealii::ArborXWrappers::SphereIntersectPredicate &sph_intersect)
  {
    return sph_intersect.size();
  }



  inline auto
  AccessTraits<dealii::ArborXWrappers::SphereIntersectPredicate,
               PredicatesTag>::
    get(const dealii::ArborXWrappers::SphereIntersectPredicate &sph_intersect,
        std::size_t                                             i)
  {
    const auto sphere = sph_intersect.get(i);
    return intersects(
      Sphere{{sphere.first[0], sphere.first[1], sphere.first[2]},
             sphere.second});
  }



  inline std::size_t
  AccessTraits<dealii::ArborXWrappers::SphereNearestPredicate, PredicatesTag>::
    size(const dealii::ArborXWrappers::SphereNearestPredicate &sph_nearest)
  {
    return sph_nearest.size();
  }



  inline auto
  AccessTraits<dealii::ArborXWrappers::SphereNearestPredicate, PredicatesTag>::
    get(const dealii::ArborXWrappers::SphereNearestPredicate &sph_nearest,
        std::size_t                                           i)
  {
    const auto sphere = sph_nearest.get(i);
    return nearest(Sphere{{sphere.first[0], sphere.first[1], sphere.first[2]},
                          sphere.second},
                   sph_nearest.get_n_nearest_neighbors());
  }
#  else
  template <int dim, typename Number>
  inline std::size_t
  AccessTraits<dealii::ArborXWrappers::PointIntersectPredicate<dim, Number>>::
    size(const dealii::ArborXWrappers::PointIntersectPredicate<dim, Number>
           &pt_intersect)
  {
    return pt_intersect.size();
  }



  template <int dim, typename Number>
  inline auto
  AccessTraits<dealii::ArborXWrappers::PointIntersectPredicate<dim, Number>>::
    get(const dealii::ArborXWrappers::PointIntersectPredicate<dim, Number>
                   &pt_intersect,
        std::size_t i)
  {
    return intersects(
      dealii::ArborXWrappers::internal::to_arborx_point(pt_intersect.get(i)));
  }



  template <int dim, typename Number>
  inline std::size_t
  AccessTraits<dealii::ArborXWrappers::PointNearestPredicate<dim, Number>>::
    size(const dealii::ArborXWrappers::PointNearestPredicate<dim, Number>
           &pt_nearest)
  {
    return pt_nearest.size();
  }



  template <int dim, typename Number>
  inline auto
  AccessTraits<dealii::ArborXWrappers::PointNearestPredicate<dim, Number>>::get(
    const dealii::ArborXWrappers::PointNearestPredicate<dim, Number>
               &pt_nearest,
    std::size_t i)
  {
    return nearest(dealii::ArborXWrappers::internal::to_arborx_point(
                     pt_nearest.get(i)),
                   pt_nearest.get_n_nearest_neighbors());
  }



  template <int dim, typename Number>
  inline std::size_t
  AccessTraits<
    dealii::ArborXWrappers::BoundingBoxIntersectPredicate<dim, Number>>::
    size(
      const dealii::ArborXWrappers::BoundingBoxIntersectPredicate<dim, Number>
        &bb_intersect)
  {
    return bb_intersect.size();
  }



  template <int dim, typename Number>
  inline auto
  AccessTraits<
    dealii::ArborXWrappers::BoundingBoxIntersectPredicate<dim, Number>>::
    get(const dealii::ArborXWrappers::BoundingBoxIntersectPredicate<dim, Number>
                   &bb_intersect,
        std::size_t i)
  {
    const auto boundary_points = bb_intersect.get(i).get_boundary_points();
    const dealii::Point<dim, Number> min_corner = boundary_points.first;
    const dealii::Point<dim, Number> max_corner = boundary_points.second;

    return intersects(
      Box{dealii::ArborXWrappers::internal::to_arborx_point(min_corner),
          dealii::ArborXWrappers::internal::to_arborx_point(max_corner)});
  }



  template <int dim, typename Number>
  inline std::size_t
  AccessTraits<
    dealii::ArborXWrappers::BoundingBoxNearestPredicate<dim, Number>>::
    size(const dealii::ArborXWrappers::BoundingBoxNearestPredicate<dim, Number>
           &bb_nearest)
  {
    return bb_nearest.size();
  }



  template <int dim, typename Number>
  inline auto
  AccessTraits<
    dealii::ArborXWrappers::BoundingBoxNearestPredicate<dim, Number>>::
    get(const dealii::ArborXWrappers::BoundingBoxNearestPredicate<dim, Number>
                   &bb_nearest,
        std::size_t i)
  {
    const auto boundary_points = bb_nearest.get(i).get_boundary_points();
    const dealii::Point<dim, Number> min_corner = boundary_points.first;
    const dealii::Point<dim, Number> max_corner = boundary_points.second;

    return nearest(
      Box{dealii::ArborXWrappers::internal::to_arborx_point(min_corner),
          dealii::ArborXWrappers::internal::to_arborx_point(max_corner)},
      bb_nearest.get_n_nearest_neighbors());
  }



  template <int dim, typename Number>
  inline std::size_t
  AccessTraits<dealii::ArborXWrappers::SphereIntersectPredicate<dim, Number>>::
    size(const dealii::ArborXWrappers::SphereIntersectPredicate<dim, Number>
           &sph_intersect)
  {
    return sph_intersect.size();
  }



  template <int dim, typename Number>
  inline auto
  AccessTraits<dealii::ArborXWrappers::SphereIntersectPredicate<dim, Number>>::
    get(const dealii::ArborXWrappers::SphereIntersectPredicate<dim, Number>
                   &sph_intersect,
        std::size_t i)
  {
    const auto sphere = sph_intersect.get(i);
    return intersects(
      Sphere{dealii::ArborXWrappers::internal::to_arborx_point(sphere.first),
             sphere.second});
  }



  template <int dim, typename Number>
  inline std::size_t
  AccessTraits<dealii::ArborXWrappers::SphereNearestPredicate<dim, Number>>::
    size(const dealii::ArborXWrappers::SphereNearestPredicate<dim, Number>
           &sph_nearest)
  {
    return sph_nearest.size();
  }



  template <int dim, typename Number>
  inline auto
  AccessTraits<dealii::ArborXWrappers::SphereNearestPredicate<dim, Number>>::
    get(const dealii::ArborXWrappers::SphereNearestPredicate<dim, Number>
                   &sph_nearest,
        std::size_t i)
  {
    const auto sphere = sph_nearest.get(i);
    return nearest(Sphere{dealii::ArborXWrappers::internal::to_arborx_point(
                            sphere.first),
                          sphere.second},
                   sph_nearest.get_n_nearest_neighbors());
  }
#  endif
} // namespace ArborX

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif

#endif
