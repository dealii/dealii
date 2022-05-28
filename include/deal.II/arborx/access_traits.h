// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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
} // namespace ArborXWrappers

DEAL_II_NAMESPACE_CLOSE

/**
 * This namespace contains the implementation of AccessTraits used by ArborX.
 */
namespace ArborX
{
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
     * Return an
     * Arbox::nearest(ArborX::Box,
     * BoundingBoxtNearestPredicate::get_n_nearest_neighbors) object constructed
     * from the
     * `i`th dealii::BoundingBox stored in @p bb_nearest.
     */
    static auto
    get(const dealii::ArborXWrappers::BoundingBoxNearestPredicate &bb_nearest,
        std::size_t                                                i);
  };



  /**
   * This struct allows ArborX to use SphereIntersectPredicate as a predicate.
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
     * SphereNearestPredicate::get_n_nearest_neightbors) object constructed from
     * the `i`th sphere stored in @p sph_nearest.
     */
    static auto
    get(const dealii::ArborXWrappers::SphereNearestPredicate &sph_nearest,
        std::size_t                                           i);
  };

  // ------------------------------- Inline ----------------------------------//

  // The implementation of AccessTraits<..., PredicatesTag> needs to be in the
  // header file otherwise the return type of auto get() cannot be determined.
  // We use auto because ArborX does not expose the type of intersects

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
} // namespace ArborX

#endif

#endif
