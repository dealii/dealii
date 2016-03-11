// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2016 by the deal.II authors
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

#ifndef dealii__dof_handler_h
#define dealii__dof_handler_h



#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/iterator_range.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>
#include <deal.II/dofs/block_info.h>
#include <deal.II/dofs/dof_iterator_selector.h>
#include <deal.II/dofs/number_cache.h>
#include <deal.II/dofs/dof_faces.h>
#include <deal.II/dofs/dof_levels.h>
#include <deal.II/dofs/function_map.h>

#include <boost/serialization/split_member.hpp>

#include <vector>
#include <map>
#include <set>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class FiniteElement;
template <int dim, int spacedim> class Triangulation;

namespace internal
{
  namespace DoFHandler
  {
    struct Implementation;

    namespace Policy
    {
      template <int dim, int spacedim> class PolicyBase;
      struct Implementation;
    }
  }

  namespace DoFAccessor
  {
    struct Implementation;
  }

  namespace DoFCellAccessor
  {
    struct Implementation;
  }
}


/**
 * Manage the distribution and numbering of the degrees of freedom for non-
 * multigrid algorithms. This class satisfies the
 * @ref ConceptMeshType "MeshType concept"
 * requirements.
 *
 * It is first used in the step-2 tutorial program.
 *
 * For each vertex, line, quad, etc, this class stores a list of the indices
 * of degrees of freedom living on this object. These indices refer to the
 * unconstrained degrees of freedom, i.e. constrained degrees of freedom are
 * numbered in the same way as unconstrained ones, and are only later
 * eliminated.  This leads to the fact that indices in global vectors and
 * matrices also refer to all degrees of freedom and some kind of condensation
 * is needed to restrict the systems of equations to the unconstrained degrees
 * of freedom only. The actual layout of storage of the indices is described
 * in the dealii::internal::DoFHandler::DoFLevel class documentation.
 *
 * The class offers iterators to traverse all cells, in much the same way as
 * the Triangulation class does. Using the begin() and end() functions (and
 * companions, like begin_active()), one can obtain iterators to walk over
 * cells, and query the degree of freedom structures as well as the
 * triangulation data. These iterators are built on top of those of the
 * Triangulation class, but offer the additional information on degrees of
 * freedom functionality compared to pure triangulation iterators. The order
 * in which dof iterators are presented by the <tt>++</tt> and <tt>--</tt>
 * operators is the same as that for the corresponding iterators traversing
 * the triangulation on which this DoFHandler is constructed.
 *
 * The <tt>spacedim</tt> parameter has to be used if one wants to solve
 * problems on surfaces. If not specified, this parameter takes the default
 * value <tt>=dim</tt> implying that we want to solve problems in a domain
 * whose dimension equals the dimension of the space in which it is embedded.
 *
 *
 * <h3>Distribution of indices for degrees of freedom</h3>
 *
 * The degrees of freedom (`dofs') are distributed on the given triangulation
 * by the function distribute_dofs(). It gets passed a finite element object
 * describing how many degrees of freedom are located on vertices, lines, etc.
 * It traverses the triangulation cell by cell and numbers the dofs of that
 * cell if not yet numbered. For non-multigrid algorithms, only active cells
 * are considered. Active cells are defined to be those cells which have no
 * children, i.e. they are the most refined ones.
 *
 * Since the triangulation is traversed starting with the cells of the
 * coarsest active level and going to more refined levels, the lowest numbers
 * for dofs are given to the largest cells as well as their bounding lines and
 * vertices, with the dofs of more refined cells getting higher numbers.
 *
 * This numbering implies very large bandwidths of the resulting matrices and
 * is thus vastly suboptimal for some solution algorithms. For this reason,
 * the DoFRenumbering class offers several algorithms to reorder the dof
 * numbering according. See there for a discussion of the implemented
 * algorithms.
 *
 *
 * <h3>Interaction with distributed meshes</h3>
 *
 * Upon construction, this class takes a reference to a triangulation object.
 * In most cases, this will be a reference to an object of type Triangulation,
 * i.e. the class that represents triangulations that entirely reside on a
 * single processor. However, it can also be of type
 * parallel::distributed::Triangulation (see, for example, step-32, step-40
 * and in particular the
 * @ref distributed
 * module) in which case the DoFHandler object will proceed to only manage
 * degrees of freedom on locally owned and ghost cells. This process is
 * entirely transparent to the used.
 *
 *
 * <h3>User defined renumbering schemes</h3>
 *
 * The DoFRenumbering class offers a number of renumbering schemes like the
 * Cuthill-McKee scheme. Basically, the function sets up an array in which for
 * each degree of freedom we store the new index this DoF should have after
 * renumbering. Using this array, the renumber_dofs() function of the present
 * class is called, which actually performs the change from old DoF indices to
 * the ones given in the array. In some cases, however, a user may want to
 * compute her own renumbering order; in this case, one can allocate an array
 * with one element per degree of freedom and fill it with the number that the
 * respective degree of freedom shall be assigned. This number may, for
 * example, be obtained by sorting the support points of the degrees of
 * freedom in downwind direction.  Then call the
 * <tt>renumber_dofs(vector<types::global_dof_index>)</tt> function with the
 * array, which converts old into new degree of freedom indices.
 *
 *
 * <h3>Serializing (loading or storing) DoFHandler objects</h3>
 *
 * Like many other classes in deal.II, the DoFHandler class can stream its
 * contents to an archive using BOOST's serialization facilities. The data so
 * stored can later be retrieved again from the archive to restore the
 * contents of this object. This facility is frequently used to save the state
 * of a program to disk for possible later resurrection, often in the context
 * of checkpoint/restart strategies for long running computations or on
 * computers that aren't very reliable (e.g. on very large clusters where
 * individual nodes occasionally fail and then bring down an entire MPI job).
 *
 * The model for doing so is similar for the DoFHandler class as it is for the
 * Triangulation class (see the section in the general documentation of that
 * class). In particular, the load() function does not exactly restore the
 * same state as was stored previously using the save() function. Rather, the
 * function assumes that you load data into a DoFHandler object that is
 * already associated with a triangulation that has a content that matches the
 * one that was used when the data was saved. Likewise, the load() function
 * assumes that the current object is already associated with a finite element
 * object that matches the one that was associated with it when data was
 * saved; the latter can be achieved by calling DoFHandler::distribute_dofs()
 * using the same kind of finite element before re-loading data from the
 * serialization archive.
 *
 * @ingroup dofs
 * @author Wolfgang Bangerth, Markus Buerg, Timo Heister, Guido Kanschat,
 * @date 1998, 1999, 2000, 2012
 */
template <int dim, int spacedim=dim>
class DoFHandler  :  public Subscriptor
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
   * @ref GlossActive "active cell iterators".
   * The concept of iterators is discussed at length in the
   * @ref Iterators "iterators documentation module".
   *
   * The current typedef identifies active cells in a DoFHandler object. While
   * the actual data type of the typedef is hidden behind a few layers of
   * (unfortunately necessary) indirections, it is in essence
   * TriaActiveIterator<DoFCellAccessor>. The TriaActiveIterator class works
   * like a pointer to active objects that when you dereference it yields an
   * object of type DoFCellAccessor. DoFCellAccessor is a class that
   * identifies properties that are specific to cells in a DoFHandler, but it
   * is derived (and consequently inherits) from both DoFAccessor,
   * TriaCellAccessor and TriaAccessor that describe what you can ask of more
   * general objects (lines, faces, as well as cells) in a triangulation and
   * DoFHandler objects.
   *
   * @ingroup Iterators
   */
  typedef typename ActiveSelector::active_cell_iterator active_cell_iterator;

  /**
   * A typedef that is used to to identify cell iterators. The concept of
   * iterators is discussed at length in the
   * @ref Iterators "iterators documentation module".
   *
   * The current typedef identifies cells in a DoFHandler object. Some of
   * these cells may in fact be active (see
   * @ref GlossActive "active cell iterators")
   * in which case they can in fact be asked for the degrees of freedom that
   * live on them. On the other hand, if the cell is not active, any such
   * query will result in an error. Note that this is what distinguishes this
   * typedef from the level_cell_iterator typedef.
   *
   * While the actual data type of the typedef is hidden behind a few layers
   * of (unfortunately necessary) indirections, it is in essence
   * TriaIterator<DoFCellAccessor>. The TriaIterator class works like a
   * pointer to objects that when you dereference it yields an object of type
   * DoFCellAccessor. DoFCellAccessor is a class that identifies properties
   * that are specific to cells in a DoFHandler, but it is derived (and
   * consequently inherits) from both DoFAccessor, TriaCellAccessor and
   * TriaAccessor that describe what you can ask of more general objects
   * (lines, faces, as well as cells) in a triangulation and DoFHandler
   * objects.
   *
   * @ingroup Iterators
   */
  typedef typename ActiveSelector::cell_iterator        cell_iterator;

  typedef typename ActiveSelector::face_iterator        face_iterator;
  typedef typename ActiveSelector::active_face_iterator active_face_iterator;

  typedef typename LevelSelector::CellAccessor          level_cell_accessor;
  typedef typename LevelSelector::FaceAccessor          level_face_accessor;

  typedef typename LevelSelector::cell_iterator         level_cell_iterator;
  typedef typename LevelSelector::face_iterator         level_face_iterator;


  /**
   * Make the dimension available in function templates.
   */
  static const unsigned int dimension = dim;

  /**
   * Make the space dimension available in function templates.
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * When the arrays holding the DoF indices are set up, but before they are
   * filled with actual values, they are set to an invalid value, in order to
   * monitor possible problems. This invalid value is the constant defined
   * here.
   *
   * Please note that you should not rely on it having a certain value, but
   * rather take its symbolic name.
   */
  static const types::global_dof_index invalid_dof_index = numbers::invalid_dof_index;

  /**
   * The default index of the finite element to be used on a given cell. Since
   * the present class only supports the same finite element to be used on all
   * cells, the index of the finite element needs to be the same on all cells
   * anyway, and by convention we pick zero for this value. The situation is
   * different for hp objects (i.e. the hp::DoFHandler class) where different
   * finite element indices may be used on different cells, and the default
   * index there corresponds to an invalid value.
   */
  static const unsigned int default_fe_index = 0;

  /**
   * Standard constructor, not initializing any data. After constructing an
   * object with this constructor, use initialize() to make a valid
   * DoFHandler.
   */
  DoFHandler ();

  /**
   * Constructor. Take @p tria as the triangulation to work on.
   */
  DoFHandler ( const Triangulation<dim,spacedim> &tria);

  /**
   * Destructor.
   */
  virtual ~DoFHandler ();

  /**
   * Assign a Triangulation and a FiniteElement to the DoFHandler and compute
   * the distribution of degrees of freedom over the mesh.
   */
  void initialize(const Triangulation<dim,spacedim> &tria,
                  const FiniteElement<dim,spacedim> &fe);

  /**
   * Go through the triangulation and "distribute" the degrees of freedoms
   * needed for the given finite element. "Distributing" degrees of freedom
   * involved allocating memory to store the information that describes it
   * (e.g., whether it is located on a vertex, edge, face, etc) and to
   * sequentially enumerate all degrees of freedom. In other words, while the
   * mesh and the finite element object by themselves simply define a finite
   * element space $V_h$, the process of distributing degrees of freedom makes
   * sure that there is a basis for this space and that the shape functions of
   * this basis are enumerated in an indexable, predictable way.
   *
   * The purpose of this function is first discussed in the introduction to
   * the step-2 tutorial program.
   *
   * @note A pointer of the finite element given as argument is stored.
   * Therefore, the lifetime of the finite element object shall be longer than
   * that of this object. If you don't want this behavior, you may want to
   * call the @p clear member function which also releases the lock of this
   * object to the finite element.
   */
  virtual void distribute_dofs (const FiniteElement<dim,spacedim> &fe);

  /**
   * Distribute level degrees of freedom on each level for geometric
   * multigrid. The active DoFs need to be distributed using distribute_dofs()
   * before calling this function and the @p fe needs to be identical to the
   * finite element passed to distribute_dofs().
   */
  virtual void distribute_mg_dofs (const FiniteElement<dim, spacedim> &fe);

  /**
   * This function returns whether this DoFHandler has DoFs distributed on
   * each multigrid level or in other words if distribute_mg_dofs() has been
   * called.
   */
  bool has_level_dofs() const;

  /**
   * This function returns whether this DoFHandler has active DoFs. This is
   * equivalent to asking whether (i) distribute_dofs() has been called and
   * (ii) the finite element for which degrees of freedom have been
   * distributed actually has degrees of freedom (which is not the case for
   * FE_Nothing, for example).
   *
   * If this object is based on a parallel::distributed::Triangulation, then
   * the current function returns true if <i>any</i> partition of the parallel
   * DoFHandler object has any degrees of freedom. In other words, the
   * function returns true even if the Triangulation does not own any active
   * cells on the current MPI process, but at least one process owns cells and
   * at least this one process has any degrees of freedom associated with it.
   */
  bool has_active_dofs() const;

  /**
   * After distribute_dofs() with an FESystem element, the block structure of
   * global and level vectors is stored in a BlockInfo object accessible with
   * block_info(). This function initializes the local block structure on each
   * cell in the same object.
   */
  void initialize_local_block_info();

  /**
   * Clear all data of this object and especially delete the lock this object
   * has to the finite element used the last time when @p distribute_dofs was
   * called.
   */
  virtual void clear ();

  /**
   * Renumber degrees of freedom based on a list of new dof numbers for all
   * the dofs.
   *
   * This function is called by the functions in DoFRenumbering function after
   * computing the ordering of the degrees of freedom. This function is
   * called, for example, by the functions in the DoFRenumbering namespace,
   * but it can of course also be called from user code.
   *
   * @arg new_number This array must have a size equal to the number of
   * degrees of freedom owned by the current processor, i.e. the size must be
   * equal to what n_locally_owned_dofs() returns. If only one processor
   * participates in storing the current mesh, then this equals the total
   * number of degrees of freedom, i.e. the result of n_dofs(). The contents
   * of this array are the new global indices for each freedom listed in the
   * IndexSet returned by locally_owned_dofs(). In the case of a sequential
   * mesh this means that the array is a list of new indices for each of the
   * degrees of freedom on the current mesh. In the case that we have a
   * parallel::distributed::Triangulation underlying this DoFHandler object,
   * the array is a list of new indices for all the locally owned degrees of
   * freedom, enumerated in the same order as the currently locally owned
   * DoFs. In other words, assume that degree of freedom <code>i</code> is
   * currently locally owned, then
   * <code>new_numbers[locally_owned_dofs().index_within_set(i)]</code>
   * returns the new global DoF index of <code>i</code>. Since the IndexSet of
   * locally_owned_dofs() is complete in the sequential case, the latter
   * convention for the content of the array reduces to the former in the case
   * that only one processor participates in the mesh.
   */
  void renumber_dofs (const std::vector<types::global_dof_index> &new_numbers);

  /**
   * The same function as above, but renumber the degrees of freedom of a
   * single level of a multigrid hierarchy.
   */
  void renumber_dofs (const unsigned int level,
                      const std::vector<types::global_dof_index> &new_numbers);

  /**
   * Return the maximum number of degrees of freedom a degree of freedom in
   * the given triangulation with the given finite element may couple with.
   * This is the maximum number of entries per line in the system matrix; this
   * information can therefore be used upon construction of the
   * SparsityPattern object.
   *
   * The returned number is not really the maximum number but an estimate
   * based on the finite element and the maximum number of cells meeting at a
   * vertex. The number holds for the constrained matrix as well.
   *
   * The determination of the number of couplings can be done by simple
   * picture drawing. An example can be found in the implementation of this
   * function.
   *
   * @note This function is most often used to determine the maximal row
   * length for sparsity patterns. Unfortunately, while the estimates returned
   * by this function are rather accurate in 1d and 2d, they are often
   * significantly too high in 3d, leading the SparsityPattern class to
   * allocate much too much memory in some cases. Unless someone comes around
   * to improving the present function for 3d, there is not very much one can
   * do about these cases. The typical way to work around this problem is to
   * use an intermediate compressed sparsity pattern that only allocates
   * memory on demand. Refer to the step-2 and step-11 example programs on how
   * to do this. The problem is also discussed in the documentation of the
   * module on
   * @ref Sparsity.
   */
  unsigned int max_couplings_between_dofs () const;

  /**
   * Return the number of degrees of freedom located on the boundary another
   * dof on the boundary can couple with.
   *
   * The number is the same as for max_couplings_between_dofs() in one
   * dimension less.
   *
   * @note The same applies to this function as to max_couplings_per_dofs() as
   * regards the performance of this function. Think about one of the dynamic
   * sparsity pattern classes instead (see
   * @ref Sparsity).
   */
  unsigned int max_couplings_between_boundary_dofs () const;

  /*--------------------------------------*/

  /**
   * @name Cell iterator functions
   */

  /*
   * @{
   */

  /**
   * Iterator to the first used cell on level @p level.
   */
  cell_iterator        begin       (const unsigned int level = 0) const;

  /**
   * Iterator to the first active cell on level @p level. If the given level
   * does not contain any active cells (i.e., all cells on this level are
   * further refined, then this function returns
   * <code>end_active(level)</code> so that loops of the kind
   *  @code
   *    for (cell=dof_handler.begin_active(level); cell!=dof_handler.end_active(level); ++cell)
   *      ...
   *  @endcode
   * have zero iterations, as may be expected if there are no active cells on
   * this level.
   */
  active_cell_iterator begin_active(const unsigned int level = 0) const;

  /**
   * Iterator past the end; this iterator serves for comparisons of iterators
   * with past-the-end or before-the-beginning states.
   */
  cell_iterator        end () const;

  /**
   * Return an iterator which is the first iterator not on the given level. If
   * @p level is the last level, then this returns <tt>end()</tt>.
   */
  cell_iterator end (const unsigned int level) const;

  /**
   * Return an active iterator which is the first active iterator not on the
   * given level. If @p level is the last level, then this returns
   * <tt>end()</tt>.
   */
  active_cell_iterator end_active (const unsigned int level) const;


  /**
   * Iterator to the first used cell on level @p level. This returns a
   * level_cell_iterator that returns level dofs when dof_indices() is called.
   */
  level_cell_iterator begin_mg (const unsigned int level = 0) const;

  /**
   * Iterator past the last cell on level @p level. This returns a
   * level_cell_iterator that returns level dofs when dof_indices() is called.
   */
  level_cell_iterator end_mg (const unsigned int level) const;

  /**
   * Iterator past the end; this iterator serves for comparisons of iterators
   * with past-the-end or before-the-beginning states.
   */
  level_cell_iterator end_mg () const;

  /**
   * @name Cell iterator functions returning ranges of iterators
   */

  /**
   * Return an iterator range that contains all cells (active or not) that
   * make up this DoFHandler. Such a range is useful to initialize range-based
   * for loops as supported by C++11. See the example in the documentation of
   * active_cell_iterators().
   *
   * @return The half open range <code>[this->begin(), this->end())</code>
   *
   * @ingroup CPP11
   */
  IteratorRange<cell_iterator>        cell_iterators () const;

  /**
   * Return an iterator range that contains all active cells that make up this
   * DoFHandler. Such a range is useful to initialize range-based for loops as
   * supported by C++11, see also
   * @ref CPP11 "C++11 standard".
   *
   * Range-based for loops are useful in that they require much less code than
   * traditional loops (see <a href="http://en.wikipedia.org/wiki/C%2B%2B11
   * #Range-based_for_loop">here</a> for a discussion of how they work). An
   * example is that without range-based for loops, one often writes code such
   * as the following:
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
   * Using C++11's range-based for loops, this is now entirely equivalent to
   * the following:
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
   * @return The half open range <code>[this->begin_active(),
   * this->end())</code>
   *
   * @ingroup CPP11
   */
  IteratorRange<active_cell_iterator> active_cell_iterators () const;

  /**
   * Return an iterator range that contains all cells (active or not) that
   * make up this DoFHandler in their level-cell form. Such a range is useful
   * to initialize range-based for loops as supported by C++11. See the
   * example in the documentation of active_cell_iterators().
   *
   * @return The half open range <code>[this->begin_mg(),
   * this->end_mg())</code>
   *
   * @ingroup CPP11
   */
  IteratorRange<level_cell_iterator>  mg_cell_iterators () const;

  /**
   * Return an iterator range that contains all cells (active or not) that
   * make up the given level of this DoFHandler. Such a range is useful to
   * initialize range-based for loops as supported by C++11. See the example
   * in the documentation of active_cell_iterators().
   *
   * @param[in] level A given level in the refinement hierarchy of this
   * triangulation.
   * @return The half open range <code>[this->begin(level),
   * this->end(level))</code>
   *
   * @pre level must be less than this->n_levels().
   *
   * @ingroup CPP11
   */
  IteratorRange<cell_iterator>        cell_iterators_on_level (const unsigned int level) const;

  /**
   * Return an iterator range that contains all active cells that make up the
   * given level of this DoFHandler. Such a range is useful to initialize
   * range-based for loops as supported by C++11. See the example in the
   * documentation of active_cell_iterators().
   *
   * @param[in] level A given level in the refinement hierarchy of this
   * triangulation.
   * @return The half open range <code>[this->begin_active(level),
   * this->end(level))</code>
   *
   * @pre level must be less than this->n_levels().
   *
   * @ingroup CPP11
   */
  IteratorRange<active_cell_iterator> active_cell_iterators_on_level (const unsigned int level) const;

  /**
   * Return an iterator range that contains all cells (active or not) that
   * make up the given level of this DoFHandler in their level-cell form. Such
   * a range is useful to initialize range-based for loops as supported by
   * C++11. See the example in the documentation of active_cell_iterators().
   *
   * @param[in] level A given level in the refinement hierarchy of this
   * triangulation.
   * @return The half open range <code>[this->begin_mg(level),
   * this->end_mg(level))</code>
   *
   * @pre level must be less than this->n_levels().
   *
   * @ingroup CPP11
   *
   */
  IteratorRange<level_cell_iterator> mg_cell_iterators_on_level (const unsigned int level) const;

  /*
   * @}
   */


  /*---------------------------------------*/


  /**
   * Return the global number of degrees of freedom. If the current object
   * handles all degrees of freedom itself (even if you may intend to solve
   * your linear system in parallel, such as in step-17 or step-18), then this
   * number equals the number of locally owned degrees of freedom since this
   * object doesn't know anything about what you want to do with it and
   * believes that it owns every degree of freedom it knows about.
   *
   * On the other hand, if this object operates on a
   * parallel::distributed::Triangulation object, then this function returns
   * the global number of degrees of freedom, accumulated over all processors.
   *
   * In either case, included in the returned number are those DoFs which are
   * constrained by hanging nodes, see
   * @ref constraints.
   */
  types::global_dof_index n_dofs () const;

  /**
   * The (global) number of multilevel degrees of freedom on a given level.
   *
   * If no level degrees of freedom have been assigned to this level, returns
   * numbers::invalid_dof_index. Else returns the number of degrees of freedom
   * on this level.
   */
  types::global_dof_index n_dofs (const unsigned int level) const;

  /**
   * Return the number of degrees of freedom located on the boundary.
   */
  types::global_dof_index n_boundary_dofs () const;

  /**
   * Return the number of degrees of freedom located on those parts of the
   * boundary which have a boundary indicator listed in the given set. The
   * reason that a @p map rather than a @p set is used is the same as
   * described in the section on the @p make_boundary_sparsity_pattern
   * function.
   *
   * The type of boundary_ids equals typename FunctionMap<spacedim,number>::type .
   */
  template<typename number>
  types::global_dof_index
  n_boundary_dofs (const std::map<types::boundary_id, const Function<spacedim,number>*> &boundary_ids) const;

  /**
   * Same function, but with different data type of the argument, which is
   * here simply a list of the boundary indicators under consideration.
   */
  types::global_dof_index
  n_boundary_dofs (const std::set<types::boundary_id> &boundary_ids) const;

  /**
   * Access to an object informing of the block structure of the dof handler.
   *
   * If an FESystem is used in distribute_dofs(), degrees of freedom naturally
   * split into several
   * @ref GlossBlock "blocks".
   * For each base element as many blocks appear as its multiplicity.
   *
   * At the end of distribute_dofs(), the number of degrees of freedom in each
   * block is counted, and stored in a BlockInfo object, which can be accessed
   * here. If you have previously called distribute_mg_dofs(), the same is
   * done on each level of the multigrid hierarchy. Additionally, the block
   * structure on each cell can be generated in this object by calling
   * initialize_local_block_info().
   */
  const BlockInfo &block_info() const;


  /**
   * Return the number of degrees of freedom that belong to this process.
   *
   * If this is a sequential job, then the result equals that produced by
   * n_dofs(). On the other hand, if we are operating on a
   * parallel::distributed::Triangulation, then it includes only the degrees
   * of freedom that the current processor owns. Note that in this case this
   * does not include all degrees of freedom that have been distributed on the
   * current processor's image of the mesh: in particular, some of the degrees
   * of freedom on the interface between the cells owned by this processor and
   * cells owned by other processors may be theirs, and degrees of freedom on
   * ghost cells are also not necessarily included.
   */
  unsigned int n_locally_owned_dofs() const;

  /**
   * Return an IndexSet describing the set of locally owned DoFs as a subset
   * of 0..n_dofs(). The number of elements of this set equals
   * n_locally_owned_dofs().
   */
  const IndexSet &locally_owned_dofs() const;

  /**
   * Returns an IndexSet describing the set of locally owned DoFs used for the
   * given multigrid level as a subset of 0..n_dofs(level).
   */
  const IndexSet &locally_owned_mg_dofs(const unsigned int level) const;


  /**
   * Returns a vector that stores the locally owned DoFs of each processor. If
   * you are only interested in the number of elements each processor owns
   * then n_locally_owned_dofs_per_processor() is a better choice.
   *
   * If this is a sequential job, then the vector has a single element that
   * equals the IndexSet representing the entire range [0,n_dofs()].
   */
  const std::vector<IndexSet> &
  locally_owned_dofs_per_processor () const;

  const std::vector<IndexSet> &
  locally_owned_mg_dofs_per_processor (const unsigned int level) const;

  /**
   * Return a vector that stores the number of degrees of freedom each
   * processor that participates in this triangulation owns locally. The sum
   * of all these numbers equals the number of degrees of freedom that exist
   * globally, i.e. what n_dofs() returns.
   *
   * Each element of the vector returned by this function equals the number of
   * elements of the corresponding sets returned by global_dof_indices().
   *
   * If this is a sequential job, then the vector has a single element equal
   * to n_dofs().
   */
  const std::vector<types::global_dof_index> &
  n_locally_owned_dofs_per_processor () const;

  /**
   * Return a constant reference to the selected finite element object.
   */
  const FiniteElement<dim,spacedim> &get_fe () const;

  /**
   * Return a constant reference to the triangulation underlying this object.
   *
   * @deprecated Use get_triangulation() instead.
   */
  const Triangulation<dim,spacedim> &get_tria () const DEAL_II_DEPRECATED;

  /**
   * Return a constant reference to the triangulation underlying this object.
   */
  const Triangulation<dim,spacedim> &get_triangulation () const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   *
   * This function is made virtual, since a dof handler object might be
   * accessed through a pointers to this base class, although the actual
   * object might be a derived class.
   */
  virtual std::size_t memory_consumption () const;

  /**
   * Write the data of this object to a stream for the purpose of
   * serialization.
   */
  template <class Archive>
  void save (Archive &ar, const unsigned int version) const;

  /**
   * Read the data of this object from a stream for the purpose of
   * serialization.
   */
  template <class Archive>
  void load (Archive &ar, const unsigned int version);

  BOOST_SERIALIZATION_SPLIT_MEMBER()

  /**
   * We are trying to renumber the degrees of freedom, but somehow did not
   * count correctly.
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcRenumberingIncomplete);
  /**
   * Exception
   * @ingroup Exceptions
   */
  DeclException0 (ExcGridsDoNotMatch);
  /**
   * Exception
   * @ingroup Exceptions
   */
  DeclException0 (ExcInvalidBoundaryIndicator);
  /**
   * Exception
   * @ingroup Exceptions
   */
  DeclException1 (ExcNewNumbersNotConsecutive,
                  types::global_dof_index,
                  << "The given list of new dof indices is not consecutive: "
                  << "the index " << arg1 << " does not exist.");
  /**
   * Exception
   * @ingroup Exceptions
   */
  DeclException1 (ExcInvalidLevel,
                  int,
                  << "The given level " << arg1
                  << " is not in the valid range!");
  /**
   * Exception
   * @ingroup Exceptions
   */
  DeclException0 (ExcFacesHaveNoLevel);
  /**
   * The triangulation level you accessed is empty.
   * @ingroup Exceptions
   */
  DeclException1 (ExcEmptyLevel,
                  int,
                  << "You tried to do something on level " << arg1
                  << ", but this level is empty.");


private:
  /**
   * Copy constructor. I can see no reason why someone might want to use it,
   * so I don't provide it. Since this class has pointer members, making it
   * private prevents the compiler to provide it's own, incorrect one if
   * anyone chose to copy such an object.
   */
  DoFHandler (const DoFHandler &);

  /**
   * Copy operator. I can see no reason why someone might want to use it, so I
   * don't provide it. Since this class has pointer members, making it private
   * prevents the compiler to provide it's own, incorrect one if anyone chose
   * to copy such an object.
   */
  DoFHandler &operator = (const DoFHandler &);


  /**
   * An object containing information on the block structure.
   */
  BlockInfo block_info_object;

  /**
   * Address of the triangulation to work on.
   */
  SmartPointer<const Triangulation<dim,spacedim>,DoFHandler<dim,spacedim> >
  tria;

  /**
   * Store a pointer to the finite element given latest for the distribution
   * of dofs. In order to avoid destruction of the object before the lifetime
   * of the DoF handler, we subscribe to the finite element object. To unlock
   * the FE before the end of the lifetime of this DoF handler, use the
   * <tt>clear()</tt> function (this clears all data of this object as well,
   * though).
   */
  SmartPointer<const FiniteElement<dim,spacedim>,DoFHandler<dim,spacedim> >
  selected_fe;

  /**
   * An object that describes how degrees of freedom should be distributed and
   * renumbered.
   */
  std_cxx11::shared_ptr<dealii::internal::DoFHandler::Policy::PolicyBase<dim,spacedim> > policy;

  /**
   * A structure that contains all sorts of numbers that characterize the
   * degrees of freedom this object works on.
   *
   * For most members of this structure, there is an accessor function in this
   * class that returns its value.
   */
  dealii::internal::DoFHandler::NumberCache number_cache;

  /**
   * Data structure like number_cache, but for each multigrid level.
   */
  std::vector<dealii::internal::DoFHandler::NumberCache> mg_number_cache;

  /**
   * A data structure that is used to store the DoF indices associated with a
   * particular vertex. Unlike cells, vertices live on several levels of a
   * multigrid hierarchy; consequently, we need to store DoF indices for each
   * vertex for each of the levels it lives on. This class does this.
   */
  class MGVertexDoFs
  {
  public:
    /**
     * Constructor.
     */
    MGVertexDoFs ();

    /**
     * Destructor.
     */
    ~MGVertexDoFs ();

    /**
     * A function that is called to allocate the necessary amount of memory to
     * store the indices of the DoFs that live on this vertex for the given
     * (inclusive) range of levels.
     */
    void init (const unsigned int coarsest_level,
               const unsigned int finest_level,
               const unsigned int dofs_per_vertex);

    /**
     * Return the coarsest level for which this structure stores data.
     */
    unsigned int get_coarsest_level () const;

    /**
     * Return the finest level for which this structure stores data.
     */
    unsigned int get_finest_level () const;

    /**
     * Return the index of the <code>dof_number</code>th degree of freedom for
     * the given level stored for the current vertex.
     */
    types::global_dof_index
    get_index (const unsigned int level,
               const unsigned int dof_number) const;

    /**
     * Set the index of the <code>dof_number</code>th degree of freedom for
     * the given level stored for the current vertex to <code>index</code>.
     */
    void set_index (const unsigned int level,
                    const unsigned int dof_number,
                    const types::global_dof_index index);

    /**
     * Exception.
     */
    DeclException0 (ExcNoMemory);

  private:
    /**
     * Coarsest level for which this object stores DoF indices.
     */
    unsigned int coarsest_level;

    /**
     * Finest level for which this object stores DoF indices.
     */
    unsigned int finest_level;

    /**
     * A pointer to an array where we store the indices of the DoFs that live
     * on the various levels this vertex exists on.
     */
    types::global_dof_index *indices;

    /**
     * This array stores, for each level starting with coarsest_level, the
     * offset in the <code>indices</code> array where the DoF indices for each
     * level are stored.
     */
    types::global_dof_index *indices_offset;
  };

  void clear_mg_space ();

  /**
   * Free all used memory.
   */
  void clear_space ();

  void reserve_space ();

  template <int structdim>
  types::global_dof_index get_dof_index (const unsigned int obj_level,
                                         const unsigned int obj_index,
                                         const unsigned int fe_index,
                                         const unsigned int local_index) const;

  template<int structdim>
  void set_dof_index (const unsigned int obj_level,
                      const unsigned int obj_index,
                      const unsigned int fe_index,
                      const unsigned int local_index,
                      const types::global_dof_index global_index) const;

  /**
   * Array to store the indices for degrees of freedom located at vertices.
   */
  std::vector<types::global_dof_index> vertex_dofs;

  /**
   * An array to store the indices for level degrees of freedom located at
   * vertices.
   */
  std::vector<MGVertexDoFs> mg_vertex_dofs;

  /**
   * Space to store the DoF numbers for the different levels. Analogous to the
   * <tt>levels[]</tt> tree of the Triangulation objects.
   */
  std::vector<dealii::internal::DoFHandler::DoFLevel<dim>*> levels;

  std::vector<dealii::internal::DoFHandler::DoFLevel<dim>*> mg_levels;

  /**
   * Space to store DoF numbers of faces. They are not stored in
   * <tt>levels</tt> since faces are not organized hierarchically, but in a
   * flat array.
   */
  dealii::internal::DoFHandler::DoFFaces<dim> *faces;

  dealii::internal::DoFHandler::DoFFaces<dim> *mg_faces;

  /**
   * Make accessor objects friends.
   */
  template <int, class, bool> friend class DoFAccessor;
  template <class, bool> friend class DoFCellAccessor;
  friend struct dealii::internal::DoFAccessor::Implementation;
  friend struct dealii::internal::DoFCellAccessor::Implementation;

  friend struct dealii::internal::DoFHandler::Implementation;
  friend struct dealii::internal::DoFHandler::Policy::Implementation;
};




/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <> types::global_dof_index DoFHandler<1>::n_boundary_dofs () const;
template <> types::global_dof_index DoFHandler<1>::n_boundary_dofs (const std::set<types::boundary_id> &) const;

template <> void DoFHandler<1>::renumber_dofs(unsigned int,const std::vector<types::global_dof_index>  &new_numbers);
template <> void DoFHandler<2>::renumber_dofs(unsigned int,const std::vector<types::global_dof_index>  &new_numbers);
template <> void DoFHandler<3>::renumber_dofs(unsigned int,const std::vector<types::global_dof_index>  &new_numbers);


/* ----------------------- Inline functions ---------------------------------- */


template <int dim, int spacedim>
inline
bool
DoFHandler<dim,spacedim>::has_level_dofs() const
{
  return mg_number_cache.size()>0;
}

template <int dim, int spacedim>
inline
bool
DoFHandler<dim,spacedim>::has_active_dofs() const
{
  return number_cache.n_global_dofs>0;
}

template <int dim, int spacedim>
inline
types::global_dof_index
DoFHandler<dim,spacedim>::n_dofs () const
{
  return number_cache.n_global_dofs;
}

template<int dim, int spacedim>
inline
types::global_dof_index DoFHandler<dim, spacedim>::n_dofs (const unsigned int level) const
{
  Assert(has_level_dofs(), ExcMessage("n_dofs(level) can only be called after distribute_mg_dofs()"));
  Assert (level < mg_number_cache.size (), ExcInvalidLevel (level));
  return mg_number_cache[level].n_global_dofs;
}


template <int dim, int spacedim>
unsigned int
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
const IndexSet &
DoFHandler<dim, spacedim>::locally_owned_mg_dofs(const unsigned int level) const
{
  Assert(level < this->get_triangulation().n_global_levels(), ExcMessage("invalid level in locally_owned_mg_dofs"));
  return mg_number_cache[level].locally_owned_dofs;
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

template <int dim, int spacedim>
const std::vector<IndexSet> &
DoFHandler<dim, spacedim>::locally_owned_mg_dofs_per_processor (const unsigned int level) const
{
  Assert(level < this->get_triangulation().n_global_levels(), ExcMessage("invalid level in locally_owned_mg_dofs_per_processor"));
  return mg_number_cache[level].locally_owned_dofs_per_processor;
}


template <int dim, int spacedim>
inline
const FiniteElement<dim,spacedim> &
DoFHandler<dim,spacedim>::get_fe () const
{
  Assert(selected_fe!=0, ExcMessage("You are trying to access the DoFHandler's FiniteElement object before it has been initialized."));
  return *selected_fe;
}



template <int dim, int spacedim>
inline
const Triangulation<dim,spacedim> &
DoFHandler<dim,spacedim>::get_tria () const
{
  Assert(tria != 0, ExcNotInitialized());
  return *tria;
}



template <int dim, int spacedim>
inline
const Triangulation<dim,spacedim> &
DoFHandler<dim,spacedim>::get_triangulation () const
{
  Assert(tria != 0, ExcNotInitialized());
  return *tria;
}



template <int dim, int spacedim>
inline
const BlockInfo &
DoFHandler<dim,spacedim>::block_info () const
{
  return block_info_object;
}


namespace internal
{
  /**
   * returns a string representing the dynamic type of the given argument.
   * This is basically the same what typeid(...).name() does, but it turns out
   * this is broken on Intel 13+.
   *
   * Defined in dof_handler.cc.
   */
  template<int dim, int spacedim>
  std::string policy_to_string(const dealii::internal::DoFHandler::Policy::PolicyBase<dim,spacedim> &policy);

}


template <int dim, int spacedim>
template <class Archive>
void DoFHandler<dim,spacedim>::save (Archive &ar,
                                     const unsigned int) const
{
  ar &block_info_object;
  ar &vertex_dofs;
  ar &number_cache;
  ar &levels;
  ar &faces;

  // write out the number of triangulation cells and later check during
  // loading that this number is indeed correct; same with something that
  // identifies the FE and the policy
  unsigned int n_cells = tria->n_cells();
  std::string  fe_name = selected_fe->get_name();
  std::string  policy_name = internal::policy_to_string(*policy);

  ar &n_cells &fe_name &policy_name;
}


template <int dim, int spacedim>
template <class Archive>
void DoFHandler<dim,spacedim>::load (Archive &ar,
                                     const unsigned int)
{
  ar &block_info_object;
  ar &vertex_dofs;
  ar &number_cache;

  // boost::serialization can restore pointers just fine, but if the
  // pointer object still points to something useful, that object is not
  // destroyed and we end up with a memory leak. consequently, first delete
  // previous content before re-loading stuff
  for (unsigned int i=0; i<levels.size(); ++i)
    delete levels[i];
  levels.resize (0);
  delete faces;
  faces = 0;

  ar &levels;
  ar &faces;

  // these are the checks that correspond to the last block in the save()
  // function
  unsigned int n_cells;
  std::string  fe_name;
  std::string  policy_name;

  ar &n_cells &fe_name &policy_name;

  AssertThrow (n_cells == tria->n_cells(),
               ExcMessage ("The object being loaded into does not match the triangulation "
                           "that has been stored previously."));
  AssertThrow (fe_name == selected_fe->get_name(),
               ExcMessage ("The finite element associated with this DoFHandler does not match "
                           "the one that was associated with the DoFHandler previously stored."));
  AssertThrow (policy_name == internal::policy_to_string(*policy),
               ExcMessage (std::string ("The policy currently associated with this DoFHandler (")
                           + internal::policy_to_string(*policy)
                           +std::string(") does not match the one that was associated with the "
                                        "DoFHandler previously stored (")
                           + policy_name
                           + ")."));
}


template<int dim, int spacedim>
inline
types::global_dof_index DoFHandler<dim, spacedim>::MGVertexDoFs::get_index (
  const unsigned int level,
  const unsigned int dof_number) const
{
  Assert ((level >= coarsest_level) && (level <= finest_level), ExcInvalidLevel (level));
  return indices[indices_offset[level - coarsest_level] + dof_number];
}


template<int dim, int spacedim>
inline
void DoFHandler<dim, spacedim>::MGVertexDoFs::set_index (
  const unsigned int level,
  const unsigned int dof_number,
  const types::global_dof_index index)
{
  Assert ((level >= coarsest_level) && (level <= finest_level), ExcInvalidLevel (level));
  indices[indices_offset[level - coarsest_level] + dof_number] = index;
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
