// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#ifndef __deal2__persistent_tria_h
#define __deal2__persistent_tria_h


#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/grid/tria.h>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * This class handles the history of a triangulation and can rebuild it after
 * it was deleted some time before. Its main purpose is support for
 * time-dependent problems where one frequently deletes a triangulation due
 * to memory pressure and later wants to rebuild it; this class has all the
 * information to rebuild it exactly as it was before including the mapping
 * of cell numbers to the geometrical cells.
 *
 * Basically, this is a drop-in replacement for the triangulation. Since it
 * is derived from the Triangulation class, it shares all the
 * functionality, but it overrides some virtual functions and adds some
 * functions, too. The main change to the base class is that it overrides
 * the @p execute_coarsening_and_refinement function, where the new version
 * first stores all refinement and coarsening flags and only then calls the
 * respective function of the base class. The stored flags may later be
 * used to restore the grid just as it was before. Some other functions
 * have been extended slightly as well, see their documentation for more
 * information.
 *
 * We note that since the triangulation is created in exactly the same state
 * as it was before, other objects working on it should result in the same
 * state as well. This holds in particular for the DoFHandler object, which
 * will assign the same degrees of freedom to the original cells and the ones
 * after reconstruction of the triangulation. You can therefore safely use data
 * vectors computed on the original grid on the reconstructed grid as well.
 *
 *
 * <h3>Usage</h3>
 * You can use objects of this class almost in the same way as objects of the
 * Triangulation class. One of the few differences is that you can only
 * construct such an object by giving a coarse grid to the constructor. The
 * coarse grid will be used to base the triangulation on, and therefore the
 * lifetime of the coarse grid has to be longer than the lifetime of the
 * object of this class.
 *
 * Basically, usage looks like this:
 * @code
 *   Triangulation<dim> coarse_grid;
 *   ...                     // initialize coarse grid
 *
 *   PersistentTriangulation<dim> grid (coarse_grid);
 *
 *   for (...)
 *     {
 *                           // restore grid from coarse grid
 *                           // and stored refinement flags
 *       grid.restore ();
 *       ...                 // do something with the grid
 *
 *       ...                 // flag some cells for refinement
 *                           // or coarsening
 *       grid.execute_coarsening_and_refinement ();
 *                           // actually refine grid and store
 *                           // the flags
 *
 *       ...                 // so something more with the grid
 *
 *       grid.clear ();      // delete the grid, but keep the
 *                           // refinement flags for later use
 *                           // in grid.restore() above
 *
 *       ...                 // do something where the grid
 *                           // is not needed anymore, e.g.
 *                           // working with another grid
 *     };
 * @endcode
 *
 * Note that initially, the PersistentTriangulation object does not
 * constitute a triangulation; it only becomes one after @p restore is first
 * called. Note also that the @p execute_coarsening_and_refinement stores
 * all necessary flags for later reconstruction using the @p restore function.
 * Triangulation::clear() resets the underlying triangulation to a
 * virgin state, but does not affect the stored refinement flags needed for
 * later reconstruction and does also not touch the coarse grid which is
 * used within restore().
 *
 * @ingroup grid
 * @author Wolfgang Bangerth, 1999
 */
template <int dim, int spacedim=dim>
class PersistentTriangulation : public Triangulation<dim, spacedim>
{
public:
  /**
   * Make the dimension available
   * in function templates.
   */
  static const unsigned int dimension = dim;
  static const unsigned int spacedimension = spacedim;

  /**
   * Build up the triangulation from the
   * coarse grid in future. Copy smoothing
   * flags, etc from that grid as well.
   * Note that the initial state of the
   * triangulation is empty, until
   * @p restore_grid is called for the
   * first time.
   *
   * The coarse grid must persist until
   * the end of this object, since it will
   * be used upon reconstruction of the
   * grid.
   */
  PersistentTriangulation (const Triangulation<dim, spacedim> &coarse_grid);

  /**
   * Copy constructor. This operation
   * is only allowed, if the triangulation
   * underlying the object to be copied
   * is presently empty. Refinement flags
   * as well as the pointer to the
   * coarse grid are copied, however.
   */
  PersistentTriangulation (const PersistentTriangulation<dim, spacedim> &old_tria);

  /**
   * Destructor.
   */
  virtual ~PersistentTriangulation ();

  /**
   * Overloaded version of the same
   * function in the base class which
   * stores the refinement and coarsening
   * flags for later reconstruction of the
   * triangulation and after that calls
   * the respective function of the
   * base class.
   */
  virtual void execute_coarsening_and_refinement ();

  /**
   * Restore the grid according to
   * the saved data. For this, the
   * coarse grid is copied and the
   * grid is stepwise rebuilt using
   * the saved flags.
   *
   * Note that this function will result in
   * an error if the underlying
   * triangulation is not empty, i.e. it
   * will only succeed if this object is
   * newly created or the <tt>clear()</tt>
   * function of the base class was called
   * on it before.
   *
   * Repeatedly calls the
   * <tt>restore(unsigned int)</tt>
   * function in a loop over all
   * refinement steps.
   */
  void restore ();

  /**
   * Differential restore. Performs
   * the @p step_noth local
   * refinement and coarsening step.
   * Step 0 stands for the copying
   * of the coarse grid.
   *
   * This function will only
   * succeed if the triangulation
   * is in just the state it were
   * if restore would have been
   * called from
   * <tt>step=0...step_no-1</tt> before.
   */
  void restore (const unsigned int step_no);

  /**
   * Returns the number of
   * refinement and coarsening
   * steps. This is given by the
   * size of the @p refine_flags
   * vector.
   */
  unsigned int n_refinement_steps () const;

  /**
   * Overload this function to use
   * @p tria as a new coarse grid. The
   * present triangulation and all
   * refinement and coarsening flags
   * storing its history are deleted,
   * and the state of the underlying
   * triangulation is reset to be
   * empty, until @p restore_grid is
   * called the next time.
   *
   * The coarse grid must persist until
   * the end of this object, since it will
   * be used upon reconstruction of the
   * grid.
   */
  virtual void copy_triangulation (const Triangulation<dim, spacedim> &tria);

  /**
   * Throw an error, since this function
   * is not useful in the context of this
   * class.
   */
  virtual void create_triangulation (const std::vector<Point<spacedim> >    &vertices,
                                     const std::vector<CellData<dim> > &cells,
                                     const SubCellData                 &subcelldata);

  /**
   * An overload of the respective function
   * of the base class.
   *
   * Throw an error, since this function
   * is not useful in the context of this
   * class.
   */
  virtual void create_triangulation_compatibility (
    const std::vector<Point<spacedim> >    &vertices,
    const std::vector<CellData<dim> > &cells,
    const SubCellData                 &subcelldata);

  /**
   * Writes all refine and coarsen
   * flags to the ostream @p out.
   */
  virtual void write_flags(std::ostream &out) const;

  /**
   * Reads all refine and coarsen flags
   * that previously were written by
   * <tt>write_flags(...)</tt>. This is especially
   * useful for rebuilding the triangulation
   * after the end or breakdown of a program
   * and its restart.
   */
  virtual void read_flags(std::istream &in);

  /**
   * Clears all flags. Retains the
   * same coarse grid.
   */
  virtual void clear_flags();

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  virtual std::size_t memory_consumption () const;

  /**
   * Exception.
   */
  DeclException0 (ExcTriaNotEmpty);
  /**
   * Exception.
   */
  DeclException0 (ExcFlagsNotCleared);

private:
  /**
   * This grid shall be used as coarse
   * grid.
   */
  SmartPointer<const Triangulation<dim,spacedim>,PersistentTriangulation<dim,spacedim> > coarse_grid;

  /**
   * Vectors holding the refinement and
   * coarsening flags of the different
   * sweeps on this time level. The
   * vectors therefore hold the history
   * of the grid.
   */
  std::vector<std::vector<bool> >   refine_flags;

  /**
   * @ref refine_flags
   */
  std::vector<std::vector<bool> >   coarsen_flags;
};


DEAL_II_NAMESPACE_CLOSE

#endif
