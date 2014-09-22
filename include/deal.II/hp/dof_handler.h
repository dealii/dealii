// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2014 by the deal.II authors
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

#ifndef __deal2__hp_dof_handler_h
#define __deal2__hp_dof_handler_h



#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/iterator_range.h>
#include <deal.II/dofs/function_map.h>
#include <deal.II/dofs/dof_iterator_selector.h>
#include <deal.II/dofs/number_cache.h>
#include <deal.II/hp/fe_collection.h>

#include <vector>
#include <map>
#include <set>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace hp
  {
    class DoFLevel;
    template <int> class DoFIndicesOnFaces;
    template <int> class DoFIndicesOnFacesOrEdges;

    namespace DoFHandler
    {
      struct Implementation;
    }
  }
}

namespace internal
{
  namespace DoFAccessor
  {
    struct Implementation;
  }

  namespace DoFCellAccessor
  {
    struct Implementation;
  }
}



namespace hp
{

  /**
   * Manage the distribution and numbering of the degrees of freedom for
   * hp-FEM algorithms. This class satisfies the requirements outlined in
   * @ref GlossMeshAsAContainer "Meshes as containers".
   *
   * The purpose of this class is to allow for an enumeration of degrees
   * of freedom in the same way as the ::DoFHandler class, but it allows
   * to use a different finite element on every cell. To this end, one
   * assigns an <code>active_fe_index</code> to every cell that indicates which
   * element within a collection of finite elements (represented by an object
   * of type hp::FECollection) is the one that lives on this cell.
   * The class then enumerates the degree of freedom associated with these
   * finite elements on each cell of a triangulation and, if possible,
   * identifies degrees of freedom at the interfaces of cells if they
   * match. If neighboring cells have degrees of freedom along the common
   * interface that do not immediate match (for example, if you have
   * $Q_2$ and $Q_3$ elements meeting at a common face), then one needs
   * to compute constraints to ensure that the resulting finite element
   * space on the mesh remains conforming.
   *
   * The whole process of working with objects of this type is explained in
   * step-27. Many of the algorithms this class implements are described
   * in the @ref hp_paper "hp paper".
   *
   *
   * <h3>Active FE indices and their behavior under mesh refinement</h3>
   *
   * The typical workflow for using this class is to create a mesh,
   * assign an active FE index to every active cell, calls
   * hp::DoFHandler::distribute_dofs(), and then assemble a linear
   * system and solve a problem on this finite element space. However,
   * one can skip assigning active FE indices upon mesh refinement
   * in certain circumstances. In particular, the following rules apply:
   * - Upon mesh refinement, child cells inherit the active FE index
   *   of the parent.
   * - On the other hand, when coarsening cells, the (now active) parent
   *   cell will not have an active FE index set and you will have to
   *   set it explicitly before calling hp::DoFHandler::distribute_dofs().
   *   In particular, to avoid stale information to be used by accident,
   *   this class deletes the active FE index of cells that are
   *   refined after inheriting this index to the children; this implies
   *   that if the children are coarsened away, the old value is no
   *   longer available on the parent cell.
   *
   * @ingroup dofs
   * @ingroup hp
   *
   * @author Wolfgang Bangerth, Oliver Kayser-Herold, 2003, 2004
   */
  template <int dim, int spacedim=dim>
  class DoFHandler : public Subscriptor
  {
    typedef dealii::internal::DoFHandler::Iterators<DoFHandler<dim,spacedim>, false> ActiveSelector;
    typedef dealii::internal::DoFHandler::Iterators<DoFHandler<dim,spacedim>, true> LevelSelector;
  public:
    typedef typename ActiveSelector::CellAccessor         cell_accessor;
    typedef typename ActiveSelector::FaceAccessor         face_accessor;

    typedef typename ActiveSelector::line_iterator        line_iterator;
    typedef typename ActiveSelector::active_line_iterator active_line_iterator;

    typedef typename ActiveSelector::quad_iterator        quad_iterator;
    typedef typename ActiveSelector::active_quad_iterator active_quad_iterator;

    typedef typename ActiveSelector::hex_iterator         hex_iterator;
    typedef typename ActiveSelector::active_hex_iterator  active_hex_iterator;

    /**
     * A typedef that is used to to identify
     * @ref GlossActive "active cell iterators". The
     * concept of iterators is discussed at length in the
     * @ref Iterators "iterators documentation module".
     *
     * The current typedef identifies active cells in a hp::DoFHandler object.
     * While the actual data type of the typedef is hidden behind a few layers
     * of (unfortunately necessary) indirections, it is in essence
     * TriaActiveIterator<DoFCellAccessor>. The TriaActiveIterator
     * class works like a pointer to active objects that when you
     * dereference it yields an object of type DoFCellAccessor.
     * DoFCellAccessor is a class that identifies properties that
     * are specific to cells in a DoFHandler, but it is derived
     * (and consequently inherits) from both DoFAccessor, TriaCellAccessor
     * and TriaAccessor that describe
     * what you can ask of more general objects (lines, faces, as
     * well as cells) in a triangulation and hp::DoFHandler objects.
     *
     * @ingroup Iterators
     */
    typedef typename ActiveSelector::active_cell_iterator active_cell_iterator;

    typedef typename LevelSelector::cell_iterator         level_cell_iterator;

    /**
     * A typedef that is used to to identify cell iterators. The
     * concept of iterators is discussed at length in the
     * @ref Iterators "iterators documentation module".
     *
     * The current typedef identifies cells in a DoFHandler object. Some
     * of these cells may in fact be active (see @ref GlossActive "active cell iterators")
     * in which case they can in fact be asked for the degrees of freedom
     * that live on them. On the other hand, if the cell is not active,
     * any such query will result in an error. Note that this is what distinguishes
     * this typedef from the level_cell_iterator typedef.
     *
     * While the actual data type of the typedef is hidden behind a few layers
     * of (unfortunately necessary) indirections, it is in essence
     * TriaIterator<DoFCellAccessor>. The TriaIterator
     * class works like a pointer to objects that when you
     * dereference it yields an object of type DoFCellAccessor.
     * DoFCellAccessor is a class that identifies properties that
     * are specific to cells in a DoFHandler, but it is derived
     * (and consequently inherits) from both DoFAccessor, TriaCellAccessor
     * and TriaAccessor that describe
     * what you can ask of more general objects (lines, faces, as
     * well as cells) in a triangulation and DoFHandler objects.
     *
     * @ingroup Iterators
     */
    typedef typename ActiveSelector::cell_iterator        cell_iterator;

    typedef typename ActiveSelector::face_iterator        face_iterator;
    typedef typename ActiveSelector::active_face_iterator active_face_iterator;

    typedef typename LevelSelector::CellAccessor          level_cell_accessor;
    typedef typename LevelSelector::FaceAccessor          level_face_accessor;

    typedef typename LevelSelector::face_iterator         level_face_iterator;

    /**
     * Alias the @p FunctionMap type
     * declared elsewhere.
     */
    typedef typename FunctionMap<spacedim>::type FunctionMap;

    /**
     * Make the dimension available
     * in function templates.
     */
    static const unsigned int dimension = dim;

    /**
     * Make the space dimension available
     * in function templates.
     */
    static const unsigned int space_dimension = spacedim;

    /**
     * When the arrays holding the
     * DoF indices are set up, but
     * before they are filled with
     * actual values, they are set to
     * an invalid value, in order to
     * monitor possible
     * problems. This invalid value
     * is the constant defined here.
     *
     * Please note that you should
     * not rely on it having a
     * certain value, but rather take
     * its symbolic name.
     */
    static const types::global_dof_index invalid_dof_index = numbers::invalid_dof_index;

    /**
     * The default index of the
     * finite element to be used on
     * a given cell. For the usual,
     * non-hp dealii::DoFHandler class
     * that only supports the same
     * finite element to be used on
     * all cells, the index of the
     * finite element needs to be
     * the same on all cells
     * anyway, and by convention we
     * pick zero for this
     * value. The situation here is
     * different, since the hp
     * classes support the case
     * where different finite
     * element indices may be used
     * on different cells. The
     * default index consequently
     * corresponds to an invalid
     * value.
     */
    static const unsigned int default_fe_index = numbers::invalid_unsigned_int;


    /**
     * Constructor. Take @p tria as the
     * triangulation to work on.
     */
    DoFHandler (const Triangulation<dim,spacedim> &tria);

    /**
     * Destructor.
     */
    virtual ~DoFHandler ();

    /**
     * Go through the triangulation and "distribute" the degrees of
     * freedoms needed for the given finite element. "Distributing"
     * degrees of freedom involved allocating memory to store the
     * information that describes it (e.g., whether it is located on a
     * vertex, edge, face, etc) and to sequentially enumerate all
     * degrees of freedom. In other words, while the mesh and the finite
     * element object by themselves simply define a finite element space
     * $V_h$, the process of distributing degrees of freedom makes sure
     * that there is a basis for this space and that the shape functions
     * of this basis are enumerated in an indexable, predictable way.
     *
     * The purpose of this function
     * is first discussed in the introduction
     * to the step-2 tutorial program.
     *
     * @note A pointer of the finite element given as argument is
     * stored. Therefore, the lifetime of the finite element object
     * shall be longer than that of this object. If you don't want this
     * behavior, you may want to call the @p clear member function which
     * also releases the lock of this object to the finite element.
     */
    virtual void distribute_dofs (const hp::FECollection<dim,spacedim> &fe);

    /**
     * Go through the triangulation and set
     * the active FE indices of all active
     * cells to the values given in @p
     * active_fe_indices.
     */
    void set_active_fe_indices (const std::vector<unsigned int> &active_fe_indices);

    /**
     * Go through the triangulation and
     * store the active FE indices of all
     * active cells to the vector @p
     * active_fe_indices. This vector is
     * resized, if necessary.
     */
    void get_active_fe_indices (std::vector<unsigned int> &active_fe_indices) const;

    /**
     * Clear all data of this object and
     * especially delete the lock this object
     * has to the finite element used the last
     * time when @p distribute_dofs was called.
     */
    virtual void clear ();

    /**
     * Renumber degrees of freedom based on
     * a list of new dof numbers for all the
     * dofs.
     *
     * @p new_numbers is an array of integers
     * with size equal to the number of dofs
     * on the present grid. It stores the new
     * indices after renumbering in the
     * order of the old indices.
     *
     * This function is called by
     * the functions in
     * DoFRenumbering function
     * after computing the ordering
     * of the degrees of freedom.
     * However, you can call this
     * function yourself, which is
     * necessary if a user wants to
     * implement an ordering scheme
     * herself, for example
     * downwind numbering.
     *
     * The @p new_number array must
     * have a size equal to the
     * number of degrees of
     * freedom. Each entry must
     * state the new global DoF
     * number of the degree of
     * freedom referenced.
     */
    void renumber_dofs (const std::vector<types::global_dof_index> &new_numbers);

    /**
     * Return the maximum number of
     * degrees of freedom a degree of freedom
     * in the given triangulation with the
     * given finite element may couple with.
     * This is the maximum number of entries
     * per line in the system matrix; this
     * information can therefore be used upon
     * construction of the SparsityPattern
     * object.
     *
     * The returned number is not really the
     * maximum number but an estimate based
     * on the finite element and the maximum
     * number of cells meeting at a vertex.
     * The number holds for the constrained
     * matrix also.
     *
     * As for
     * ::DoFHandler::max_couplings_between_dofs(),
     * the result of this function is often
     * not very accurate for 3d and/or high
     * polynomial degrees. The consequences
     * are discussed in the documentation
     * of the module on @ref Sparsity.
     */
    unsigned int max_couplings_between_dofs () const;

    /**
     * Return the number of degrees of freedom
     * located on the boundary another dof on
     * the boundary can couple with.
     *
     * The number is the same as for
     * @p max_coupling_between_dofs in one
     * dimension less.
     *
     * @note The same applies to this function as to max_couplings_per_dofs()
     * as regards the performance of this function. Think about one of the
     * dynamic sparsity pattern classes instead (see @ref Sparsity).
     */
    unsigned int max_couplings_between_boundary_dofs () const;

    /**
    *  @name Cell iterator functions
    */
    /*@{*/
    /**
    * Iterator to the first used
    * cell on level @p level.
    */
    cell_iterator        begin       (const unsigned int level = 0) const;

    /**
     *  Iterator to the first active cell on level @p level. If the
     *  given level does not contain any active cells (i.e., all cells
     *  on this level are further refined, then this function returns
     *  <code>end_active(level)</code> so that loops of the kind
     *  @code
     *    for (cell=dof_handler.begin_active(level); cell!=dof_handler.end_active(level); ++cell)
     *      ...
     *  @endcode
     *  have zero iterations, as may be expected if there are no active
     *  cells on this level.
     */
    active_cell_iterator begin_active(const unsigned int level = 0) const;

    /**
    * Iterator past the end; this
    * iterator serves for
    * comparisons of iterators with
    * past-the-end or
    * before-the-beginning states.
    */
    cell_iterator        end () const;

    /**
    * Return an iterator which is
    * the first iterator not on
    * level. If @p level is the
    * last level, then this returns
    * <tt>end()</tt>.
    */
    cell_iterator        end (const unsigned int level) const;

    /**
     * Return an active iterator which is the first active iterator not
     * on the given level. If @p level is the last level, then this
     * returns <tt>end()</tt>.
     */
    active_cell_iterator end_active (const unsigned int level) const;

    /*@}*/

    /**
     *  @name Cell iterator functions returning ranges of iterators
     */

    /**
     * Return an iterator range that contains all cells (active or not)
     * that make up this DoFHandler. Such a range is useful to
     * initialize range-based for loops as supported by C++11. See the
     * example in the documentation of active_cell_iterators().
     *
     * @return The half open range <code>[this->begin(), this->end())</code>
     *
     * @ingroup CPP11
     */
    IteratorRange<cell_iterator>        cell_iterators () const;

    /**
     * Return an iterator range that contains all active cells
     * that make up this DoFHandler. Such a range is useful to
     * initialize range-based for loops as supported by C++11,
     * see also @ref CPP11 "C++11 standard".
     *
     * Range-based for loops are useful in that they require much less
     * code than traditional loops (see
     * <a href="http://en.wikipedia.org/wiki/C%2B%2B11#Range-based_for_loop">here</a>
     * for a discussion of how they work). An example is that without
     * range-based for loops, one often writes code such as the following:
     * @code
     *   DoFHandler<dim> dof_handler;
     *   ...
     *   typename DoFHandler<dim>::active_cell_iterator
     *     cell = dof_handler.begin_active(),
     *     endc = dof_handler.end();
     *   for (; cell!=endc; ++cell)
     *     {
     *       fe_values.reinit (cell);
     *       ...do the local integration on 'cell'...;
     *     }
     * @endcode
     * Using C++11's range-based for loops, this is now entirely
     * equivalent to the following:
     * @code
     *   DoFHandler<dim> dof_handler;
     *   ...
     *   for (auto cell : dof_handler.active_cell_iterators())
     *     {
     *       fe_values.reinit (cell);
     *       ...do the local integration on 'cell'...;
     *     }
     * @endcode
     * To use this feature, you need a compiler that supports C++11.
     *
     * @return The half open range <code>[this->begin_active(), this->end())</code>
     *
     * @ingroup CPP11
     */
    IteratorRange<active_cell_iterator> active_cell_iterators () const;

    /**
     * Return an iterator range that contains all cells (active or not)
     * that make up the given level of this DoFHandler. Such a range is useful to
     * initialize range-based for loops as supported by C++11. See the
     * example in the documentation of active_cell_iterators().
     *
     * @param[in] level A given level in the refinement hierarchy of this
     *   triangulation.
     * @return The half open range <code>[this->begin(level), this->end(level))</code>
     *
     * @pre level must be less than this->n_levels().
     *
     * @ingroup CPP11
     */
    IteratorRange<cell_iterator>        cell_iterators_on_level (const unsigned int level) const;

    /**
     * Return an iterator range that contains all active cells
     * that make up the given level of this DoFHandler. Such a range is useful to
     * initialize range-based for loops as supported by C++11. See the
     * example in the documentation of active_cell_iterators().
     *
     * @param[in] level A given level in the refinement hierarchy of this
     *   triangulation.
     * @return The half open range <code>[this->begin_active(level), this->end(level))</code>
     *
     * @pre level must be less than this->n_levels().
     *
     * @ingroup CPP11
     */
    IteratorRange<active_cell_iterator> active_cell_iterators_on_level (const unsigned int level) const;

    /*@}*/

    /*---------------------------------------*/


    /**
     * Return the global number of
     * degrees of freedom. If the
     * current object handles all
     * degrees of freedom itself
     * (even if you may intend to
     * solve your linear system in
     * parallel, such as in step-17
     * or step-18), then this number
     * equals the number of locally
     * owned degrees of freedom since
     * this object doesn't know
     * anything about what you want
     * to do with it and believes
     * that it owns every degree of
     * freedom it knows about.
     *
     * On the other hand, if this
     * object operates on a
     * parallel::distributed::Triangulation
     * object, then this function
     * returns the global number of
     * degrees of freedom,
     * accumulated over all
     * processors.
     *
     * In either case, included in
     * the returned number are those
     * DoFs which are constrained by
     * hanging nodes, see @ref constraints.
     */
    types::global_dof_index n_dofs () const;

    /**
    * The number of multilevel
    * dofs on given level. Since
    * hp::DoFHandler does not
    * support multilevel methods
    * yet, this function returns
    * numbers::invalid_unsigned
    * int independent of its argument.
    */
    types::global_dof_index n_dofs(const unsigned int level) const;

    /**
     * Return the number of degrees of freedom
     * located on the boundary.
     */
    types::global_dof_index n_boundary_dofs () const;

    /**
     * Return the number of degrees
     * of freedom located on those
     * parts of the boundary which
     * have a boundary indicator
     * listed in the given set. The
     * reason that a @p map rather
     * than a @p set is used is the
     * same as described in the
     * section on the
     * @p make_boundary_sparsity_pattern
     * function.
     */
    types::global_dof_index
    n_boundary_dofs (const FunctionMap &boundary_indicators) const;

    /**
     * Same function, but with
     * different data type of the
     * argument, which is here simply
     * a list of the boundary
     * indicators under
     * consideration.
     */
    types::global_dof_index
    n_boundary_dofs (const std::set<types::boundary_id> &boundary_indicators) const;

    /**
     * Return the number of
     * degrees of freedom that
     * belong to this
     * process.
     *
     * If this is a sequential job,
     * then the result equals that
     * produced by n_dofs(). On the
     * other hand, if we are
     * operating on a
     * parallel::distributed::Triangulation,
     * then it includes only the
     * degrees of freedom that the
     * current processor owns. Note
     * that in this case this does
     * not include all degrees of
     * freedom that have been
     * distributed on the current
     * processor's image of the mesh:
     * in particular, some of the
     * degrees of freedom on the
     * interface between the cells
     * owned by this processor and
     * cells owned by other
     * processors may be theirs, and
     * degrees of freedom on ghost
     * cells are also not necessarily
     * included.
     */
    types::global_dof_index n_locally_owned_dofs() const;

    /**
     * Return an IndexSet describing
     * the set of locally owned DoFs
     * as a subset of
     * 0..n_dofs(). The number of
     * elements of this set equals
     * n_locally_owned_dofs().
     */
    const IndexSet &locally_owned_dofs() const;


    /**
     * Returns a vector that
     * stores the locally owned
     * DoFs of each processor. If
     * you are only interested in
     * the number of elements
     * each processor owns then
     * n_dofs_per_processor() is
     * a better choice.
     *
     * If this is a sequential job,
     * then the vector has a single
     * element that equals the
     * IndexSet representing the
     * entire range [0,n_dofs()].
     */
    const std::vector<IndexSet> &
    locally_owned_dofs_per_processor () const;

    /**
     * Return a vector that
     * stores the number of
     * degrees of freedom each
     * processor that
     * participates in this
     * triangulation owns
     * locally. The sum of all
     * these numbers equals the
     * number of degrees of
     * freedom that exist
     * globally, i.e. what
     * n_dofs() returns.
     *
     * Each element of the vector
     * returned by this function
     * equals the number of
     * elements of the
     * corresponding sets
     * returned by
     * global_dof_indices().
     *
     * If this is a sequential job,
     * then the vector has a single
     * element equal to n_dofs().
     */
    const std::vector<types::global_dof_index> &
    n_locally_owned_dofs_per_processor () const;

    /**
     * Return a constant reference to
     * the set of finite element
     * objects that are used by this
     * @p DoFHandler.
     */
    const hp::FECollection<dim,spacedim> &get_fe () const;

    /**
     * Return a constant reference to the
     * triangulation underlying this object.
     */
    const Triangulation<dim,spacedim> &get_tria () const;

    /**
     * Determine an estimate for the
     * memory consumption (in bytes)
     * of this object.
     *
     * This function is made virtual,
     * since a dof handler object
     * might be accessed through a
     * pointers to thisr base class,
     * although the actual object
     * might be a derived class.
     */
    virtual std::size_t memory_consumption () const;

    /**
     * Exception
     */
    DeclException0 (ExcInvalidTriangulation);
    /**
     * Exception
     */
    DeclException0 (ExcNoFESelected);
    /**
     * Exception
     */
    DeclException0 (ExcRenumberingIncomplete);
    /**
     * Exception
     */
    DeclException0 (ExcGridsDoNotMatch);
    /**
     * Exception
     */
    DeclException0 (ExcInvalidBoundaryIndicator);
    /**
     * Exception
     */
    DeclException1 (ExcMatrixHasWrongSize,
                    int,
                    << "The matrix has the wrong dimension " << arg1);
    /**
     *  Exception
     */
    DeclException0 (ExcFunctionNotUseful);
    /**
     * Exception
     */
    DeclException1 (ExcNewNumbersNotConsecutive,
                    types::global_dof_index,
                    << "The given list of new dof indices is not consecutive: "
                    << "the index " << arg1 << " does not exist.");
    /**
     * Exception
     */
    DeclException2 (ExcInvalidFEIndex,
                    int, int,
                    << "The mesh contains a cell with an active_fe_index of "
                    << arg1 << ", but the finite element collection only has "
                    << arg2 << " elements");
    /**
     *  Exception
     */
    DeclException1 (ExcInvalidLevel,
                    int,
                    << "The given level " << arg1
                    << " is not in the valid range!");
    /**
     * Exception
     */
    DeclException0 (ExcFacesHaveNoLevel);
    /**
     * The triangulation level you
     * accessed is empty.
     */
    DeclException1 (ExcEmptyLevel,
                    int,
                    << "You tried to do something on level " << arg1
                    << ", but this level is empty.");

  protected:

    /**
     * Address of the triangulation to
     * work on.
     */
    SmartPointer<const Triangulation<dim,spacedim>,DoFHandler<dim,spacedim> > tria;

    /**
     * Store a pointer to the finite
     * element set given latest for
     * the distribution of dofs. In
     * order to avoid destruction of
     * the object before the lifetime
     * of the DoF handler, we
     * subscribe to the finite
     * element object. To unlock the
     * FE before the end of the
     * lifetime of this DoF handler,
     * use the <tt>clear()</tt> function
     * (this clears all data of this
     * object as well, though).
     */
    SmartPointer<const hp::FECollection<dim,spacedim>,hp::DoFHandler<dim,spacedim> > finite_elements;

  private:

    /**
     * Copy constructor. I can see no reason
     * why someone might want to use it, so
     * I don't provide it. Since this class
     * has pointer members, making it private
     * prevents the compiler to provide it's
     * own, incorrect one if anyone chose to
     * copy such an object.
     */
    DoFHandler (const DoFHandler &);

    /**
     * Copy operator. I can see no reason
     * why someone might want to use it, so
     * I don't provide it. Since this class
     * has pointer members, making it private
     * prevents the compiler to provide it's
     * own, incorrect one if anyone chose to
     * copy such an object.
     */
    DoFHandler &operator = (const DoFHandler &);

    class MGVertexDoFs
    {
    public:
      MGVertexDoFs ();
      ~MGVertexDoFs ();
      types::global_dof_index get_index (const unsigned int level, const unsigned int dof_number) const;
      void set_index (const unsigned int level, const unsigned int dof_number, const types::global_dof_index index);
    };

    /**
     * Free all used memory.
     */
    void clear_space ();

    template<int structdim>
    types::global_dof_index get_dof_index (const unsigned int obj_level, const unsigned int obj_index, const unsigned int fe_index, const unsigned int local_index) const;

    template<int structdim>
    void set_dof_index (const unsigned int obj_level, const unsigned int obj_index, const unsigned int fe_index, const unsigned int local_index, const types::global_dof_index global_index) const;

    /**
     *  Create default tables for
     *  the active_fe_indices in
     *  the
     *  dealii::internal::hp::DoFLevel. They
     *  are initialized with a
     *  zero indicator, meaning
     *  that fe[0] is going to be
     *  used by default.  This
     *  method is called before
     *  refinement and before
     *  distribute_dofs is
     *  called. It ensures each
     *  cell has a valid
     *  active_fe_index.
     */
    void create_active_fe_table ();

    /**
     *  Functions that will be triggered
     *  through signals whenever the
     *  triangulation is modified.
     *
     *  Here they are used to
     *  administrate the the
     *  active_fe_fields during the
     *  spatial refinement.
     */
    void pre_refinement_action ();
    void post_refinement_action ();

    /**
     * Compute identities between
     * DoFs located on
     * vertices. Called from
     * distribute_dofs().
     */
    void
    compute_vertex_dof_identities (std::vector<types::global_dof_index> &new_dof_indices) const;

    /**
     * Compute identities between
     * DoFs located on
     * lines. Called from
     * distribute_dofs().
     */
    void
    compute_line_dof_identities (std::vector<types::global_dof_index> &new_dof_indices) const;

    /**
     * Compute identities between
     * DoFs located on
     * quads. Called from
     * distribute_dofs().
     */
    void
    compute_quad_dof_identities (std::vector<types::global_dof_index> &new_dof_indices) const;

    /**
     * Renumber the objects with
     * the given and all lower
     * structural dimensions,
     * i.e. renumber vertices by
     * giving a template argument
     * of zero to the int2type
     * argument, lines and vertices
     * with one, etc.
     *
     * Note that in contrast to the
     * public renumber_dofs()
     * function, these internal
     * functions do not ensure that
     * the new DoFs are
     * contiguously numbered. The
     * function may therefore also
     * be used to assign different
     * DoFs the same number, for
     * example to unify hp DoFs
     * corresponding to different
     * finite elements but
     * co-located on the same
     * entity.
     */
    void renumber_dofs_internal (const std::vector<types::global_dof_index> &new_numbers,
                                 dealii::internal::int2type<0>);

    void renumber_dofs_internal (const std::vector<types::global_dof_index> &new_numbers,
                                 dealii::internal::int2type<1>);

    void renumber_dofs_internal (const std::vector<types::global_dof_index> &new_numbers,
                                 dealii::internal::int2type<2>);

    void renumber_dofs_internal (const std::vector<types::global_dof_index> &new_numbers,
                                 dealii::internal::int2type<3>);

    /**
     * Space to store the DoF
     * numbers for the different
     * levels. Analogous to the
     * <tt>levels[]</tt> tree of
     * the Triangulation objects.
     */
    std::vector<dealii::internal::hp::DoFLevel *> levels;

    /**
     * Space to store the DoF
     * numbers for the faces.
     * Analogous to the
     * <tt>faces</tt> pointer of
     * the Triangulation objects.
     */
    dealii::internal::hp::DoFIndicesOnFaces<dim> *faces;

    /**
     * A structure that contains all
     * sorts of numbers that
     * characterize the degrees of
     * freedom this object works on.
     *
     * For most members of this
     * structure, there is an
     * accessor function in this
     * class that returns its value.
     */
    dealii::internal::DoFHandler::NumberCache number_cache;

    /**
     * Array to store the indices
     * for degrees of freedom
     * located at vertices.
     *
     * The format used here, in the
     * form of a linked list, is
     * the same as used for the
     * arrays used in the
     * internal::hp::DoFLevel
     * hierarchy. Starting indices
     * into this array are provided
     * by the vertex_dofs_offsets
     * field.
     *
     * Access to this field is
     * generally through the
     * DoFAccessor::get_vertex_dof_index() and
     * DoFAccessor::set_vertex_dof_index()
     * functions, encapsulating the
     * actual data format used to
     * the present class.
     */
    std::vector<types::global_dof_index> vertex_dofs;

    /**
     * For each vertex in the
     * triangulation, store the
     * offset within the
     * vertex_dofs array where the
     * dofs for this vertex start.
     *
     * As for that array, the
     * format is the same as
     * described in the
     * documentation of
     * hp::DoFLevel.
     *
     * Access to this field is
     * generally through the
     * Accessor::get_vertex_dof_index() and
     * Accessor::set_vertex_dof_index()
     * functions, encapsulating the
     * actual data format used to
     * the present class.
     */
    std::vector<types::global_dof_index>      vertex_dofs_offsets;

    std::vector<MGVertexDoFs> mg_vertex_dofs;  // we should really remove this field!

    /**
     * Array to store the
     * information if a cell on
     * some level has children or
     * not. It is used by the
     * signal slots as a
     * persistent buffer during the
     * refinement, i.e. from between
     * when pre_refinement_action is
     * called and when post_refinement_action
     * runs.
     */
    std::vector<std::vector<bool> *> has_children;

    /**
     * A list of connections with which this object connects
     * to the triangulation to get information about when the
     * triangulation changes.
     */
    std::vector<boost::signals2::connection> tria_listeners;

    /**
     * Make accessor objects friends.
     */
    template <int, class, bool> friend class dealii::DoFAccessor;
    template <class, bool> friend class dealii::DoFCellAccessor;
    friend struct dealii::internal::DoFAccessor::Implementation;
    friend struct dealii::internal::DoFCellAccessor::Implementation;

    /**
     * Likewise for DoFLevel
     * objects since they need to
     * access the vertex dofs in
     * the functions that set and
     * retrieve vertex dof indices.
     */
    template <int> friend class dealii::internal::hp::DoFIndicesOnFacesOrEdges;
    friend struct dealii::internal::hp::DoFHandler::Implementation;
  };



#ifndef DOXYGEN


  /* ----------------------- Inline functions ---------------------------------- */

  template <int dim, int spacedim>
  inline
  types::global_dof_index
  DoFHandler<dim,spacedim>::n_dofs () const
  {
    return number_cache.n_global_dofs;
  }


  template <int dim, int spacedim>
  inline
  types::global_dof_index
  DoFHandler<dim,spacedim>::n_dofs (const unsigned int) const
  {
    return numbers::invalid_dof_index;
  }


  template <int dim, int spacedim>
  types::global_dof_index
  DoFHandler<dim, spacedim>::n_locally_owned_dofs() const
  {
    return number_cache.n_locally_owned_dofs;
  }


  template <int dim, int spacedim>
  const IndexSet &
  DoFHandler<dim, spacedim>::locally_owned_dofs() const
  {
    return number_cache.locally_owned_dofs;
  }


  template <int dim, int spacedim>
  const std::vector<types::global_dof_index> &
  DoFHandler<dim, spacedim>::n_locally_owned_dofs_per_processor() const
  {
    return number_cache.n_locally_owned_dofs_per_processor;
  }


  template <int dim, int spacedim>
  const std::vector<IndexSet> &
  DoFHandler<dim, spacedim>::locally_owned_dofs_per_processor () const
  {
    return number_cache.locally_owned_dofs_per_processor;
  }



  template<int dim, int spacedim>
  inline
  const hp::FECollection<dim,spacedim> &
  DoFHandler<dim,spacedim>::get_fe () const
  {
    Assert (finite_elements != 0,
            ExcMessage ("No finite element collection is associated with "
                        "this DoFHandler"));
    return *finite_elements;
  }


  template<int dim, int spacedim>
  inline
  const Triangulation<dim,spacedim> &
  DoFHandler<dim,spacedim>::get_tria () const
  {
    return *tria;
  }

  template<int dim, int spacedim>
  inline
  DoFHandler<dim, spacedim>::MGVertexDoFs::MGVertexDoFs()
  {
    Assert (false, ExcNotImplemented ());
  }

  template<int dim, int spacedim>
  inline
  DoFHandler<dim, spacedim>::MGVertexDoFs::~MGVertexDoFs()
  {
    Assert (false, ExcNotImplemented ());
  }

  template<int dim, int spacedim>
  inline
  types::global_dof_index DoFHandler<dim, spacedim>::MGVertexDoFs::get_index (const unsigned int,
      const unsigned int) const
  {
    Assert (false, ExcNotImplemented ());
    return invalid_dof_index;
  }

  template<int dim, int spacedim>
  inline
  void DoFHandler<dim, spacedim>::MGVertexDoFs::set_index (const unsigned int,
                                                           const unsigned int,
                                                           types::global_dof_index)
  {
    Assert (false, ExcNotImplemented ());
  }


#endif

}

DEAL_II_NAMESPACE_CLOSE

#endif
