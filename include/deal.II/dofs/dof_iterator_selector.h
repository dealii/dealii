// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii_dof_iterators_h
#define dealii_dof_iterators_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

template <int, int, int>
class DoFInvalidAccessor;

template <int structdim, typename DoFHandlerType, bool lda>
class DoFAccessor;
template <typename DoFHandlerType, bool lda>
class DoFCellAccessor;

template <typename Accessor>
class TriaRawIterator;
template <typename Accessor>
class TriaIterator;
template <typename Accessor>
class TriaActiveIterator;

namespace internal
{
  namespace DoFHandlerImplementation
  {
    template <typename DoFHandlerType, bool lda = false>
    struct Iterators;

    /**
     * Define some types for DoF handling in one dimension.
     *
     * The types have the same meaning as those declared in
     * internal::TriangulationImplementation::Iterators<1,spacedim>, only the treatment of
     * templates is a little more complicated. See the
     * @ref Iterators
     * module for more information.
     *
     * @author Wolfgang Bangerth, Oliver Kayser-Herold, Guido Kanschat, 1998,
     * 2003, 2008, 2010
     */
    template <template <int, int> class DoFHandlerType, int spacedim, bool lda>
    struct Iterators<DoFHandlerType<1, spacedim>, lda>
    {
      typedef DoFHandlerType<1, spacedim>                   DoFHandler_type;
      typedef dealii::DoFCellAccessor<DoFHandler_type, lda> CellAccessor;
      typedef dealii::DoFAccessor<0, DoFHandler_type, lda>  FaceAccessor;

      typedef TriaRawIterator<CellAccessor>    raw_line_iterator;
      typedef TriaIterator<CellAccessor>       line_iterator;
      typedef TriaActiveIterator<CellAccessor> active_line_iterator;

      typedef TriaRawIterator<DoFInvalidAccessor<2, 1, spacedim>>
                                                               raw_quad_iterator;
      typedef TriaIterator<DoFInvalidAccessor<2, 1, spacedim>> quad_iterator;
      typedef TriaActiveIterator<DoFInvalidAccessor<2, 1, spacedim>>
        active_quad_iterator;

      typedef TriaRawIterator<DoFInvalidAccessor<3, 1, spacedim>>
                                                               raw_hex_iterator;
      typedef TriaIterator<DoFInvalidAccessor<3, 1, spacedim>> hex_iterator;
      typedef TriaActiveIterator<DoFInvalidAccessor<3, 1, spacedim>>
        active_hex_iterator;

      typedef raw_line_iterator    raw_cell_iterator;
      typedef line_iterator        cell_iterator;
      typedef active_line_iterator active_cell_iterator;

      typedef TriaRawIterator<FaceAccessor>    raw_face_iterator;
      typedef TriaIterator<FaceAccessor>       face_iterator;
      typedef TriaActiveIterator<FaceAccessor> active_face_iterator;
    };

    /**
     * Define some types for DoF handling in two dimensions.
     *
     * The types have the same meaning as those declared in
     * internal::TriangulationImplementation::Iterators<2,spacedim>, only the treatment of
     * templates is a little more complicated. See the
     * @ref Iterators
     * module for more information.
     *
     * @author Wolfgang Bangerth, Oliver Kayser-Herold, Guido Kanschat, 1998,
     * 2003, 2008, 2010
     */
    template <template <int, int> class DoFHandlerType, int spacedim, bool lda>
    struct Iterators<DoFHandlerType<2, spacedim>, lda>
    {
      typedef DoFHandlerType<2, spacedim>                   DoFHandler_type;
      typedef dealii::DoFCellAccessor<DoFHandler_type, lda> CellAccessor;
      typedef dealii::DoFAccessor<1, DoFHandler_type, lda>  FaceAccessor;

      typedef TriaRawIterator<FaceAccessor>    raw_line_iterator;
      typedef TriaIterator<FaceAccessor>       line_iterator;
      typedef TriaActiveIterator<FaceAccessor> active_line_iterator;

      typedef TriaRawIterator<CellAccessor>    raw_quad_iterator;
      typedef TriaIterator<CellAccessor>       quad_iterator;
      typedef TriaActiveIterator<CellAccessor> active_quad_iterator;

      typedef TriaRawIterator<DoFInvalidAccessor<3, 2, spacedim>>
                                                               raw_hex_iterator;
      typedef TriaIterator<DoFInvalidAccessor<3, 2, spacedim>> hex_iterator;
      typedef TriaActiveIterator<DoFInvalidAccessor<3, 2, spacedim>>
        active_hex_iterator;

      typedef raw_quad_iterator    raw_cell_iterator;
      typedef quad_iterator        cell_iterator;
      typedef active_quad_iterator active_cell_iterator;

      typedef raw_line_iterator    raw_face_iterator;
      typedef line_iterator        face_iterator;
      typedef active_line_iterator active_face_iterator;
    };

    /**
     * Define some types for DoF handling in three dimensions.
     *
     * The types have the same meaning as those declared in
     * internal::TriangulationImplementation::Iterators<3,spacedim>, only the treatment of
     * templates is a little more complicated. See the
     * @ref Iterators
     * module for more information.
     *
     * @author Wolfgang Bangerth, Oliver Kayser-Herold, Guido Kanschat, 1998,
     * 2003, 2008, 2010
     */
    template <template <int, int> class DoFHandlerType, int spacedim, bool lda>
    struct Iterators<DoFHandlerType<3, spacedim>, lda>
    {
      typedef DoFHandlerType<3, spacedim>                   DoFHandler_type;
      typedef dealii::DoFCellAccessor<DoFHandler_type, lda> CellAccessor;
      typedef dealii::DoFAccessor<2, DoFHandler_type, lda>  FaceAccessor;

      typedef TriaRawIterator<dealii::DoFAccessor<1, DoFHandler_type, lda>>
        raw_line_iterator;
      typedef TriaIterator<dealii::DoFAccessor<1, DoFHandler_type, lda>>
        line_iterator;
      typedef TriaActiveIterator<dealii::DoFAccessor<1, DoFHandler_type, lda>>
        active_line_iterator;

      typedef TriaRawIterator<FaceAccessor>    raw_quad_iterator;
      typedef TriaIterator<FaceAccessor>       quad_iterator;
      typedef TriaActiveIterator<FaceAccessor> active_quad_iterator;

      typedef TriaRawIterator<CellAccessor>    raw_hex_iterator;
      typedef TriaIterator<CellAccessor>       hex_iterator;
      typedef TriaActiveIterator<CellAccessor> active_hex_iterator;

      typedef raw_hex_iterator    raw_cell_iterator;
      typedef hex_iterator        cell_iterator;
      typedef active_hex_iterator active_cell_iterator;

      typedef raw_quad_iterator    raw_face_iterator;
      typedef quad_iterator        face_iterator;
      typedef active_quad_iterator active_face_iterator;
    };
  } // namespace DoFHandlerImplementation
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_dof_iterator_selector_h
