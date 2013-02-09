//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2010, 2011, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__dof_levels_h
#define __deal2__dof_levels_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/dofs/dof_objects.h>
#include <vector>


DEAL_II_NAMESPACE_OPEN

template <int, int> class DoFHandler;
template <int, int> class MGDoFHandler;


namespace internal
{
  namespace DoFHandler
  {


    /**
     * Structure for storing degree of freedom information for cells,
     * organized by levels.
     *
     * We store are cached values for the DoF indices on
     * each cell in#cell_dof_indices_cache, since this is a frequently requested operation. The
     * values are set by DoFCellAccessor::update_cell_dof_indices_cache
     * and are used by DoFCellAccessor::get_dof_indices.
     *
     * Note that vertices are separate from, and in fact have nothing to
     * do with cells. The indices of degrees of freedom located on
     * vertices therefore are not stored here, but rather in member
     * variables of the dealii::DoFHandler class.
     *
     * The indices of degrees of freedom located on lower dimensional
     * objects, i.e. on lines for 2D and on quads and lines for 3D are
     * treated similarly than that on cells. However, theses geometrical
     * objects, which are called faces as a generalisation, are not
     * organised in a hierarchical structure of levels. Therefore, the
     * degrees of freedom located on these objects are stored in separate
     * classes, namely the <tt>DoFFaces</tt> classes.
     *
     * Access to this object is usually
     * through the DoFAccessor::set_dof_index() and
     * DoFAccessor::dof_index() functions or similar functions of derived
     * classes that in turn access the member variables using the
     * DoFHandler::get_dof_index() and corresponding setter functions.
     * Knowledge of the actual data format is therefore
     * encapsulated to the present hierarchy of classes as well as the
     * dealii::DoFHandler class.
     *
     * @author Wolfgang Bangerth, 1998, 2006, Guido Kanschat, 2012
     */
    template <int dim>
    class DoFLevel
    {
    public:
      /**
       * Cache for the DoF indices
       * on cells. The size of this
       * array equals the number of
       * cells on a given level
       * times
       * selected_fe.dofs_per_cell.
       */
      std::vector<types::global_dof_index> cell_dof_indices_cache;

      /**
       * The object containing dof-indices
       * and related access-functions
       */
      DoFObjects<dim> dof_object;


      /**
       * Determine an estimate for the
       * memory consumption (in bytes)
       * of this object.
       */
      std::size_t memory_consumption () const;

      /**
       * Read or write the data of this object to or
       * from a stream for the purpose of serialization
       */
      template <class Archive>
      void serialize(Archive &ar,
                     const unsigned int version);
    };


    template <int dim>
    inline
    std::size_t
    DoFLevel<dim>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (cell_dof_indices_cache) +
              MemoryConsumption::memory_consumption (dof_object));
    }


    template <int dim>
    template <class Archive>
    inline
    void
    DoFLevel<dim>::serialize (Archive &ar,
                              const unsigned int)
    {
      ar &cell_dof_indices_cache;
      ar &dof_object;
    }
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
