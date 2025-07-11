// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_tria_iterator_selector_h
#define dealii_tria_iterator_selector_h


#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class CellAccessor;
template <int, int, int>
class InvalidAccessor;
template <int, int, int>
class TriaAccessor;
template <int dim, int spacedim>
class TriaAccessor<0, dim, spacedim>;
template <typename Accessor>
class TriaRawIterator;
template <typename Accessor>
class TriaIterator;
template <typename Accessor>
class TriaActiveIterator;
#endif

namespace internal
{
  namespace TriangulationImplementation
  {
    template <int dim, int spacedim>
    struct Iterators;

    /**
     * This class implements some types which differ between the dimensions.
     * These are the declarations for the 1d case only. See the
     * @ref Iterators
     * topic for more information.
     *
     * A @p line_iterator is aliased to an iterator operating on the @p
     * lines member variable of a <tt>Triangulation<1></tt> object. An @p
     * active_line_iterator only operates on the active lines. @p
     * raw_line_iterator objects operate on all lines, used or not.
     *
     * Since we are in one dimension, the following identities are declared:
     *  @code
     *    using raw_cell_iterator = raw_line_iterator;
     *    using cell_iterator = line_iterator;
     *    using active_cell_iterator = active_line_iterator;
     *  @endcode
     *
     * To enable the declaration of @p begin_quad and the like in
     * <tt>Triangulation<1></tt>, the @p quad_iterators are declared as
     * iterators over InvalidAccessor. Thus these types exist, but are useless
     * and will certainly make any involuntary use visible. The same holds for
     * hexahedron iterators.
     *
     * The same applies for the @p face_iterator types, since lines have no
     * substructures apart from vertices, which are handled in a different
     * way, however.
     */
    template <int spacedim>
    struct Iterators<1, spacedim>
    {
      using raw_line_iterator =
        TriaRawIterator<dealii::CellAccessor<1, spacedim>>;
      using line_iterator = TriaIterator<dealii::CellAccessor<1, spacedim>>;
      using active_line_iterator =
        TriaActiveIterator<dealii::CellAccessor<1, spacedim>>;

      using raw_quad_iterator =
        TriaRawIterator<dealii::InvalidAccessor<2, 1, spacedim>>;
      using quad_iterator =
        TriaIterator<dealii::InvalidAccessor<2, 1, spacedim>>;
      using active_quad_iterator =
        TriaActiveIterator<dealii::InvalidAccessor<2, 1, spacedim>>;

      using raw_hex_iterator =
        TriaRawIterator<dealii::InvalidAccessor<3, 1, spacedim>>;
      using hex_iterator =
        TriaIterator<dealii::InvalidAccessor<3, 1, spacedim>>;
      using active_hex_iterator =
        TriaActiveIterator<dealii::InvalidAccessor<3, 1, spacedim>>;
    };



    /**
     * This class implements some types which differ between the dimensions.
     * These are the declarations for the 2d case only. See the
     * @ref Iterators
     * topic for more information.
     *
     * A @p line_iterator is aliased to an iterator operating on the @p
     * lines member variable of a <tt>Triangulation<2></tt> object. An @p
     * active_line_iterator only operates on the active lines. @p
     * raw_line_iterator objects operate on all lines, used or not. Using @p
     * active_line_iterators may not be particularly in 2d useful since it
     * only operates on unrefined lines. However, also refined lines may bound
     * unrefined cells if the neighboring cell is refined once more than the
     * present one.
     *
     * Similarly to line iterators, @p quad_iterator, @p raw_quad_iterator and
     * @p active_quad_iterator are declared.
     *
     * To enable the declaration of @p begin_hex and the like in
     * <tt>Triangulation<[12]></tt>, the @p hex_iterators are declared as
     * iterators over InvalidAccessor. Thus these types exist, but are useless
     * and will certainly make any involuntary use visible.
     *
     * Since we are in two dimension, the following identities are declared:
     *  @code
     *    using raw_cell_iterator = raw_quad_iterator;
     *    using cell_iterator = quad_iterator;
     *    using active_cell_iterator = active_quad_iterator;
     *
     *    using raw_face_iterator = raw_line_iterator;
     *    using face_iterator = line_iterator;
     *    using active_face_iterator = active_line_iterator;
     *  @endcode
     */
    template <int spacedim>
    struct Iterators<2, spacedim>
    {
      using raw_line_iterator =
        TriaRawIterator<dealii::TriaAccessor<1, 2, spacedim>>;
      using line_iterator = TriaIterator<dealii::TriaAccessor<1, 2, spacedim>>;
      using active_line_iterator =
        TriaActiveIterator<dealii::TriaAccessor<1, 2, spacedim>>;

      using raw_quad_iterator =
        TriaRawIterator<dealii::CellAccessor<2, spacedim>>;
      using quad_iterator = TriaIterator<dealii::CellAccessor<2, spacedim>>;
      using active_quad_iterator =
        TriaActiveIterator<dealii::CellAccessor<2, spacedim>>;

      using raw_hex_iterator =
        TriaRawIterator<dealii::InvalidAccessor<3, 2, spacedim>>;
      using hex_iterator =
        TriaIterator<dealii::InvalidAccessor<3, 2, spacedim>>;
      using active_hex_iterator =
        TriaActiveIterator<dealii::InvalidAccessor<3, 2, spacedim>>;
    };


    /**
     * This class implements some types which differ between the dimensions.
     * These are the declarations for the 3d case only. See the
     * @ref Iterators
     * topic for more information.
     *
     * For the declarations of the data types, more or less the same holds as
     * for lower dimensions (see <tt>Iterators<[12]></tt>). The dimension
     * specific data types are here, since we are in three dimensions:
     *  @code
     *    using raw_cell_iterator = raw_hex_iterator;
     *    using cell_iterator = hex_iterator;
     *    using active_cell_iterator = active_hex_iterator;
     *
     *    using raw_face_iterator = raw_quad_iterator;
     *    using face_iterator = quad_iterator;
     *    using active_face_iterator = active_quad_iterator;
     *  @endcode
     */
    template <int spacedim>
    struct Iterators<3, spacedim>
    {
      using raw_line_iterator =
        TriaRawIterator<dealii::TriaAccessor<1, 3, spacedim>>;
      using line_iterator = TriaIterator<dealii::TriaAccessor<1, 3, spacedim>>;
      using active_line_iterator =
        TriaActiveIterator<dealii::TriaAccessor<1, 3, spacedim>>;

      using raw_quad_iterator =
        TriaRawIterator<dealii::TriaAccessor<2, 3, spacedim>>;
      using quad_iterator = TriaIterator<dealii::TriaAccessor<2, 3, spacedim>>;
      using active_quad_iterator =
        TriaActiveIterator<dealii::TriaAccessor<2, 3, spacedim>>;

      using raw_hex_iterator =
        TriaRawIterator<dealii::CellAccessor<3, spacedim>>;
      using hex_iterator = TriaIterator<dealii::CellAccessor<3, spacedim>>;
      using active_hex_iterator =
        TriaActiveIterator<dealii::CellAccessor<3, spacedim>>;
    };

  } // namespace TriangulationImplementation

} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_tria_iterator_selector_h
