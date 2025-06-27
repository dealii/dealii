// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_solution_transfer_h
#define dealii_solution_transfer_h


/*----------------------------   solutiontransfer.h     ----------------------*/


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/observer_pointer.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/vector.h>

#include <boost/range/iterator_range_core.hpp>

#include <vector>

DEAL_II_NAMESPACE_OPEN


/**
 * This class implements the transfer of a discrete FE function (e.g. a
 * solution vector) from one mesh to another that is obtained from the first
 * by a single refinement and/or coarsening step. During interpolation the
 * vector is filled with the interpolated values. This class is used in the
 * step-15, step-26, step-31, and step-33 tutorial programs. This class
 * works both for serial and distributed meshes.
 *
 * <h3>Usage</h3>
 *
 * @note It is important to note, that if you use more than one
 * SolutionTransfer object at the same time, that the calls to prepare_*()
 * and interpolate()/deserialize() need to be in the same order.
 *
 * <h3>Transferring a solution</h3> Here VectorType is your favorite
 * vector type, e.g. Vector, LinearAlgebra::distributed::Vector,
 * PETScWrappers::MPI::Vector, TrilinosWrappers::MPI::Vector, or corresponding
 * block vectors.
 * @code
 * SolutionTransfer<dim, VectorType> soltrans(dof_handler);
 * // flag some cells for refinement and coarsening, e.g.
 * GridRefinement::
 *   refine_and_coarsen_fixed_fraction(tria,
 *                                     error_indicators,
 *                                     0.3,
 *                                     0.05);
 * // Use parallel::distributed::GridRefinement for distributed triangulations.
 *
 * // prepare the triangulation,
 * tria.prepare_coarsening_and_refinement();
 *
 * // prepare the SolutionTransfer object for coarsening and refinement
 * // and give the solution vector that we intend to interpolate later,
 * soltrans.prepare_for_coarsening_and_refinement(solution);
 *
 * // actually execute the refinement,
 * tria.execute_coarsening_and_refinement();
 *
 * // redistribute dofs,
 * dof_handler.distribute_dofs(fe);
 *
 * // adjust the vector size so that it can hold all dofs
 * // (for distributed grids, see instructions below),
 * solution.reinit(...);
 *
 * // and interpolate the solution
 * soltrans.interpolate(solution);
 * @endcode
 *
 * <h3>Usage on distributed grids</h3>
 * If the grid is distributed, it is important to note that the old
 * solution(s) must be copied to one that also provides access to the
 * locally relevant DoF values (these values required for the interpolation
 * process):
 * @code
 * // Create initial indexsets pertaining to the grid before refinement
 * const IndexSet &locally_owned_dofs    = dof_handler.locally_owned_dofs();
 * const IndexSet  locally_relevant_dofs =
 *   DoFTools::extract_locally_relevant_dofs(dof_handler);
 *
 * // The solution vector only knows about locally owned DoFs
 * TrilinosWrappers::MPI::Vector solution;
 * solution.reinit(locally_owned_dofs,
 *                 mpi_communicator);
 * ...
 * // Transfer solution to vector that provides access to
 * // locally relevant DoFs
 * TrilinosWrappers::MPI::Vector old_solution;
 * old_solution.reinit(locally_owned_dofs,
 *                     locally_relevant_dofs,
 *                     mpi_communicator);
 * old_solution = solution;
 *
 * // Initialize SolutionTransfer object
 * SolutionTransfer<dim, VectorType> soltrans(dof_handler);
 * soltrans.prepare_for_coarsening_and_refinement(old_solution);
 * ...
 * // Refine grid
 * // Recreate locally_owned_dofs and locally_relevant_dofs index sets
 * ...
 * solution.reinit(locally_owned_dofs, mpi_communicator);
 * soltrans.interpolate(solution);
 * @endcode
 *
 * <h4>Note on ghost elements</h4> In a parallel computation PETSc or
 * Trilinos vector may contain ghost elements or not. For reading in
 * information with prepare_for_coarsening_and_refinement() or
 * prepare_for_serialization() you need to supply vectors with ghost
 * elements, so that all locally_active elements can be read. On the other
 * hand, ghosted vectors are generally not writable, so for calls to
 * interpolate() or deserialize() you need to supply distributed vectors
 * without ghost elements. More precisely, during interpolation the
 * current algorithm writes into all locally active degrees of freedom.
 *
 * Different from PETSc and Trilinos vectors, LinearAlgebra::distributed::Vector
 * allows writing into ghost elements. For a ghosted vector the interpolation
 * step can be accomplished via
 * @code
 * interpolated_solution.zero_out_ghost_values();
 * soltrans.interpolate(interpolated_solution);
 * interpolated_solution.update_ghost_values();
 * @endcode
 *
 * <h3>Use for serialization</h3>
 *
 * This class can be used to serialize and later deserialize a
 * mesh with solution vectors to a file. If you use more than one
 * DoFHandler and therefore more than one SolutionTransfer object, they
 * need to be serialized and deserialized in the same order.
 *
 * If vector has the locally relevant DoFs, serialization works as
 * follows:
 * @code
 * SolutionTransfer<dim, VectorType> sol_trans(dof_handler);
 * sol_trans.prepare_for_serialization(vector);
 *
 * triangulation.save(filename);
 * @endcode
 * For deserialization the vector needs to be a distributed vector
 * (without ghost elements):
 * @code
 * //[create coarse mesh...]
 * triangulation.load(filename);
 *
 * SolutionTransfer<dim, VectorType> sol_trans(dof_handler);
 * sol_trans.deserialize(distributed_vector);
 * @endcode
 *
 *
 * <h3>Note on usage with DoFHandler with hp-capabilities</h3>
 *
 * Since data on DoFHandler objects with hp-capabilities is associated with
 * many different FiniteElement objects, each cell's data has to be
 * processed with its corresponding `future_fe_index`. Further, if
 * refinement is involved, data will be packed on the parent cell with its
 * `future_fe_index` and unpacked later with the same index on its children.
 * If cells get coarsened into one, data will be packed on the children with
 * the least dominant finite element of their common subspace, and unpacked
 * on the parent with this particular finite element (consult
 * hp::FECollection::find_dominated_fe_extended() for more information).
 *
 * Transferring a solution across refinement works exactly like in the
 * non-hp-case. However, when considering serialization, we also have to
 * store the active FE indices in an additional step. A code snippet
 * demonstrating serialization with the SolutionTransfer class with DoFHandler
 * objects with hp-capabilities is provided in the following. Here VectorType is
 * your favorite vector type, e.g. Vector, LinearAlgebra::distributed::Vector,
 * PETScWrappers::MPI::Vector, TrilinosWrappers::MPI::Vector, or corresponding
 * block vectors.
 *
 * If vector has the locally relevant DoFs, serialization works as follows:
 * @code
 * SolutionTransfer<dim, VectorType, DoFHandler<dim,spacedim>>
 *     sol_trans(hp_dof_handler);
 *
 * hp_dof_handler.prepare_for_serialization_of_active_fe_indices();
 * sol_trans.prepare_for_serialization(vector);
 *
 * triangulation.save(filename);
 * @endcode
 *
 * For deserialization the vector needs to be a distributed vector
 * (without ghost elements):
 * @code
 * //[create coarse mesh...]
 * triangulation.load(filename);
 *
 * hp::FECollection<dim,spacedim> fe_collection;
 * //[prepare identical fe_collection...]
 *
 * DoFHandler<dim,spacedim> hp_dof_handler(triangulation);
 * // We need to introduce our dof_handler to the fe_collection
 * // before setting all active FE indices.
 * hp_dof_handler.deserialize_active_fe_indices();
 * hp_dof_handler.distribute_dofs(fe_collection);
 *
 * SolutionTransfer<dim,VectorType,DoFHandler<dim,spacedim>>
 *     sol_trans(hp_dof_handler);
 * sol_trans.deserialize(distributed_vector);
 * @endcode
 *
 *
 * <h3>Interaction with hanging nodes</h3>
 *
 * This class does its best to represent on the new mesh the finite element
 * function that existed on the old mesh, but this may lead to situations
 * where the function on the new mesh is no longer conforming at hanging
 * nodes. To this end, consider a situation of a twice refined mesh that
 * started with a single square cell (i.e., we now have 16 cells). Consider
 * also that we coarsen 4 of the cells back to the first refinement level. In
 * this case, we end up with a mesh that will look as follows if we were to
 * use a $Q_1$ element:
 *
 * @image html hanging_nodes.png ""
 *
 * The process of interpolating from the old to the new mesh would imply that
 * the values of the finite element function will not change on all of the
 * cells that remained as they are (i.e., the fine cells) but that on the
 * coarse cell at the top right, the four values at the vertices are obtained
 * by interpolating down from its former children.  If the original function
 * was not linear, this implies that the marked hanging nodes will retain
 * their old values which, in general, will not lead to a continuous function
 * along the corresponding edges. In other words, the solution vector obtained
 * after SolutionTransfer::interpolate() does not satisfy hanging node
 * constraints: it corresponds to the pointwise interpolation, but not to the
 * interpolation <i>onto the new finite element space that contains
 * constraints from hanging nodes</i>.
 *
 * Whether this is a problem you need to worry about or not depends on your
 * application. The situation is easily corrected, of course, by applying
 * AffineConstraints::distribute() to your solution vector after transfer,
 * using a constraints object computed on the new DoFHandler object (you
 * probably need to create this object anyway if you have hanging nodes). This
 * is also what is done, for example, in step-15.
 *
 * @note This situation can only happen if you do coarsening. If all cells
 * remain as they are or are refined, then SolutionTransfer::interpolate()
 * computes a new vector of nodal values, but the function represented is of
 * course exactly the same because the old finite element space is a subspace
 * of the new one. Thus, if the old function was conforming (i.e., satisfied
 * hanging node constraints), then so does the new one, and it is not
 * necessary to call AffineConstraints::distribute().
 *
 *
 * <h3>Implementation in the context of hp-finite elements</h3>
 *
 * In the case of DoFHandlers with hp-capabilities, nothing defines which of
 * the finite elements that are part of the hp::FECollection associated with
 * the DoFHandler, should be considered on cells that are not active (i.e.,
 * that have children). This is because degrees of freedom are only allocated
 * for active cells and, in fact, it is not allowed to set an active FE index
 * on non-active cells using DoFAccessor::set_active_fe_index().
 *
 * It is, thus, not entirely natural what should happen if, for example, a few
 * cells are coarsened away. This class then implements the following
 * algorithm:
 * - If a cell is refined, then the values of the solution vector(s) are
 *   interpolated before refinement on the to-be-refined cell from the space
 * of the active finite element to the one of the future finite element. These
 *   values are then distributed on the finite element spaces of the children
 *   post-refinement. This may lose information if, for example, the old cell
 *   used a Q2 space and the children use Q1 spaces, or the information may be
 *   prolonged if the parent cell used a Q1 space and the children are Q2s.
 * - If cells are to be coarsened, then the values from the child cells are
 *   interpolated to the parent cell using the largest of the child cell
 * future finite element spaces, which will be identified as the least
 * dominant element following the FiniteElementDomination logic (consult
 *   hp::FECollection::find_dominated_fe_extended() for more information). For
 *   example, if the children of a cell use Q1, Q2 and Q3 spaces, then the
 *   values from the children are interpolated into a Q3 space on the parent
 *   cell. After refinement, this Q3 function on the parent cell is then
 *   interpolated into the space the user has selected for this cell (which
 * may be different from Q3, in this example, if the user has set the active
 * FE index for a different space post-refinement and before calling
 *   DoFHandler::distribute_dofs()).
 *
 * @note In the context of hp-refinement, if cells are coarsened or the
 * polynomial degree is lowered on some cells, then the old finite element
 * space is not a subspace of the new space and you may run into the same
 * situation as discussed above with hanging nodes. You may want to consider
 * calling AffineConstraints::distribute() on the vector obtained by
 * transferring the solution.
 *
 * @ingroup numerics
 */
template <int dim, typename VectorType = Vector<double>, int spacedim = dim>
class SolutionTransfer
{
public:
  /**
   * Constructor.
   *
   * @param[in] dof_handler The DoFHandler on which all operations will
   * happen. At the time when this constructor is called, the DoFHandler
   * still points to the Triangulation before the refinement in question
   *   happens.
   * @param[in] average_values Average the contribututions to the same
   *   DoF coming from different cells. Note: averaging requires an
   * additional communication step, since the valence of the DoF has to be
   * determined.
   */
  SolutionTransfer(const DoFHandler<dim, spacedim> &dof_handler,
                   const bool                       average_values = false);

  /**
   * Destructor.
   */
  ~SolutionTransfer() = default;

  /**
   * Prepare the current object for coarsening and refinement. It
   * stores pointers to the vectors in @p all_in to be interpolated onto the new
   * (refined and/or coarsened) grid, and registers this object for data
   * transfer on the grid.
   */
  void
  prepare_for_coarsening_and_refinement(
    const std::vector<const VectorType *> &all_in);

  /**
   * Same as above but without pointers.
   */
  void
  prepare_for_coarsening_and_refinement(const std::vector<VectorType> &all_in);

  /**
   * Same as the previous function but for only one discrete function to be
   * interpolated.
   */
  void
  prepare_for_coarsening_and_refinement(const VectorType &in);

  /**
   * Interpolate the data previously stored in this object before the mesh
   * was refined or coarsened onto the current set of cells. Do so for
   * each of the vectors provided to
   * prepare_for_coarsening_and_refinement() and write the result into the
   * given set of vectors.
   */
  void
  interpolate(std::vector<VectorType *> &all_out);

  /**
   * Same as above but without pointers.
   */
  void
  interpolate(std::vector<VectorType> &all_out);

  /**
   * Same as the previous function. It interpolates only one function. It
   * assumes the vector has the right size (i.e.
   * <tt>out.size()==n_dofs_refined</tt>)
   *
   * Multiple calling of this function is NOT allowed. Interpolating
   * several functions can be performed in one step by using
   * <tt>interpolate (all_out)</tt>
   */
  void
  interpolate(VectorType &out);

  /**
   * Prepare the serialization of the given vector. The serialization is
   * done by Triangulation::save(). The given vector needs all information
   * on the locally active DoFs (it must be ghosted). See documentation of
   * this class for more information.
   */
  void
  prepare_for_serialization(const VectorType &in);

  /**
   * Same as the function above, only for a list of vectors.
   */
  void
  prepare_for_serialization(const std::vector<const VectorType *> &all_in);

  /**
   * Execute the deserialization of the given vector. This needs to be
   * done after calling Triangulation::load(). The given vector must be a
   * fully distributed vector without ghost elements. See documentation of
   * this class for more information.
   */
  void
  deserialize(VectorType &in);


  /**
   * Same as the function above, only for a list of vectors.
   */
  void
  deserialize(std::vector<VectorType *> &all_in);

  /**
   * Reinit this class to the state that it has directly after calling the
   * constructor.
   */
  DEAL_II_DEPRECATED_EARLY void
  clear();

private:
  /**
   * Pointer to the degree of freedom handler to work with.
   */
  ObserverPointer<const DoFHandler<dim, spacedim>,
                  SolutionTransfer<dim, VectorType, spacedim>>
    dof_handler;

  /**
   * Flag indicating if averaging should be performed.
   */
  const bool average_values;

  /**
   * A vector that stores pointers to all the vectors we are supposed to
   * copy over from the old to the new mesh.
   */
  std::vector<const VectorType *> input_vectors;

  /**
   * The handle that the Triangulation has assigned to this object
   * with which we can access our memory offset and our pack function.
   */
  unsigned int handle;

  /**
   * A callback function used to pack the data on the current mesh into
   * objects that can later be retrieved after refinement, coarsening and
   * repartitioning.
   */
  std::vector<char>
  pack_callback(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellStatus                                            status);

  /**
   * A callback function used to unpack the data on the current mesh that
   * has been packed up previously on the mesh before refinement,
   * coarsening and repartitioning.
   */
  void
  unpack_callback(
    const typename Triangulation<dim, spacedim>::cell_iterator     &cell,
    const CellStatus                                                status,
    const boost::iterator_range<std::vector<char>::const_iterator> &data_range,
    std::vector<VectorType *>                                      &all_out,
    VectorType                                                     &valence);


  /**
   * Registers the pack_callback() function to the Triangulation that has been
   * assigned to the DoFHandler class member and stores the returning handle.
   */
  void
  register_data_attach();
};



namespace Legacy
{
  /**
   * This class implements the transfer of a discrete FE function (e.g. a
   * solution vector) from one mesh to another that is obtained from the first
   * by a single refinement and/or coarsening step. During interpolation the
   * vector is reinitialized to the new size and filled with the interpolated
   * values. This class is used in the step-15, step-26, step-31, and step-33
   * tutorial programs. A version of this class that works on all types of
   * triangulations, including distributed ones, is available as
   * dealii::SolutionTransfer.
   *
   * <h3>Usage</h3>
   *
   * This class implements the algorithms in two different ways:
   * <ul>
   * <li> If the grid will only be refined (i.e. no cells are coarsened) then
   * use @p SolutionTransfer as follows:
   * @code
   * SolutionTransfer<dim, Vector<double> > soltrans(*dof_handler);
   *
   * // flag some cells for refinement, e.g.
   * GridRefinement::refine_and_coarsen_fixed_fraction(*tria,
   *                                                   error_indicators,
   *                                                   0.3,
   *                                                   0);
   * // prepare the triangulation for refinement,
   * tria->prepare_coarsening_and_refinement();
   *
   * // tell the SolutionTransfer object that we intend to do pure refinement,
   * soltrans.prepare_for_pure_refinement();
   *
   * // actually execute the refinement,
   * tria->execute_coarsening_and_refinement();
   *
   * // and redistribute dofs.
   * dof_handler->distribute_dofs (fe);
   * @endcode
   *
   * Then to proceed do
   * @code
   * // take a copy of the solution vector
   * Vector<double> solution_old(solution);
   *
   * // resize solution vector to the correct size, as the @p refine_interpolate
   * // function requires the vectors to be of right sizes
   * solution.reinit(dof_handler->n_dofs());
   *
   * // and finally interpolate
   * soltrans.refine_interpolate(solution_old, solution);
   * @endcode
   *
   * Although the @p refine_interpolate functions are allowed to be called
   * multiple times, e.g. for interpolating several solution vectors, there is
   * the following possibility of interpolating several functions
   * simultaneously.
   * @code
   * std::vector<Vector<double> > solutions_old(n_vectors, Vector<double> (n));
   * ...
   * std::vector<Vector<double> > solutions(n_vectors, Vector<double> (n));
   * soltrans.refine_interpolate(solutions_old, solutions);
   * @endcode
   * This is used in several of the tutorial programs, for example step-31
   * and step-33.
   *
   * <li> If the grid has cells that will be coarsened, then use @p
   * SolutionTransfer as follows:
   * @code
   * SolutionTransfer<dim, Vector<double> > soltrans(*dof_handler);
   *
   * // flag some cells for refinement and coarsening, e.g.
   * GridRefinement::refine_and_coarsen_fixed_fraction(*tria,
   *                                                   error_indicators,
   *                                                   0.3,
   *                                                   0.05);
   *
   * // prepare the triangulation,
   * tria->prepare_coarsening_and_refinement();
   *
   * // prepare the SolutionTransfer object for coarsening and refinement and
   * give
   * // the solution vector that we intend to interpolate later,
   * soltrans.prepare_for_coarsening_and_refinement(solution);
   *
   * // actually execute the refinement,
   * tria->execute_coarsening_and_refinement ();
   *
   * // redistribute dofs,
   * dof_handler->distribute_dofs (fe);
   *
   * // and interpolate the solution
   * Vector<double> interpolate_solution(dof_handler->n_dofs());
   * soltrans.interpolate(solution, interpolated_solution);
   * @endcode
   *
   * If the grid is partitioned across several MPI processes, then it is
   * important to note that the old solution(s) must be copied to one that
   * also provides access to the locally relevant DoF values (these values
   * required for the interpolation process):
   * @code
   * // Create initial indexsets pertaining to the grid before refinement
   * const IndexSet &locally_owned_dofs    = dof_handler.locally_owned_dofs();
   * const IndexSet  locally_relevant_dofs =
   *   DoFTools::extract_locally_relevant_dofs(dof_handler);
   *
   * // The solution vector only knows about locally owned DoFs
   * TrilinosWrappers::MPI::Vector solution;
   * solution.reinit(locally_owned_dofs,
   *                 mpi_communicator);
   * ...
   * // Transfer solution to vector that provides access to locally relevant
   * DoFs TrilinosWrappers::MPI::Vector old_solution;
   * old_solution.reinit(locally_owned_dofs,
   *                     locally_relevant_dofs,
   *                     mpi_communicator);
   * old_solution = solution;
   * ...
   * // Refine grid
   * // Recreate locally_owned_dofs and locally_relevant_dofs index sets
   * ...
   * solution.reinit(locally_owned_dofs, mpi_communicator);
   * soltrans.refine_interpolate(old_solution, solution);
   * @endcode
   *
   * Multiple calls to the function <code>interpolate (const VectorType &in,
   * VectorType &out)</code> are NOT allowed. Interpolating several
   * functions can be performed in one step by using <tt>void interpolate (const
   * vector<VectorType> &all_in, vector<VectorType> &all_out)
   * const</tt>, and using the respective @p
   * prepare_for_coarsening_and_refinement function taking several vectors as
   * input before actually refining and coarsening the triangulation (see
   * there).
   * </ul>
   *
   * For deleting all stored data in @p SolutionTransfer and reinitializing it
   * use the <tt>clear()</tt> function.
   *
   * The template argument @p VectorType denotes the type of data container you
   * want to transfer.
   *
   *
   * <h3>Interpolating in the presence of hanging nodes and boundary values</h3>
   *
   * The interpolation onto the new mesh is a local operation, i.e., it
   * interpolates onto the new mesh only. If that new mesh has hanging nodes,
   * you will therefore get a solution that does not satisfy hanging node
   * constraints. The same is true with boundary values: the interpolated
   * solution will just be the interpolation of the old solution at the
   * boundary, and this may or may not satisfy boundary values at newly
   * introduced boundary nodes.
   *
   * Consequently, you may have to apply hanging node or boundary value
   * constraints after interpolation. step-15 and step-26 have examples of
   * dealing with this.
   *
   *
   * <h3>Implementation</h3>
   *
   * <ul>
   * <li> Solution transfer with only refinement. Assume that we have got a
   * solution vector on the current (original) grid. Each entry of this vector
   * belongs to one of the DoFs of the discretization. If we now refine the grid
   * then the calling of DoFHandler::distribute_dofs() will change at least some
   * of the DoF indices. Hence we need to store the DoF indices of all active
   * cells before the refinement. A pointer for each active cell is used to
   * point to the vector of these DoF indices of that cell. This is done by
   * prepare_for_pure_refinement().
   *
   * In the function <tt>refine_interpolate(in,out)</tt> and on each cell where
   * the pointer is set (i.e. the cells that were active in the original grid)
   * we can now access the local values of the solution vector @p in on that
   * cell by using the stored DoF indices. These local values are interpolated
   * and set into the vector @p out that is at the end the discrete function @p
   * in interpolated on the refined mesh.
   *
   * The <tt>refine_interpolate(in,out)</tt> function can be called multiple
   * times for arbitrary many discrete functions (solution vectors) on the
   * original grid.
   *
   * <li> Solution transfer with coarsening and refinement. After calling
   * Triangulation::prepare_coarsening_and_refinement the coarsen flags of
   * either all or none of the children of a (parent-)cell are set. While
   * coarsening (Triangulation::execute_coarsening_and_refinement) the cells
   * that are not needed any more will be deleted from the Triangulation.
   *
   * For the interpolation from the (to be coarsenend) children to their parent
   * the children cells are needed. Hence this interpolation and the storing of
   * the interpolated values of each of the discrete functions that we want to
   * interpolate needs to take place before these children cells are coarsened
   * (and deleted!!). Again a pointer for each relevant cell is set to point to
   * these values (see below). Additionally the DoF indices of the cells that
   * will not be coarsened need to be stored according to the solution transfer
   * with pure refinement (cf there). All this is performed by
   * <tt>prepare_for_coarsening_and_refinement(all_in)</tt> where the
   * <tt>vector<VectorType> all_in</tt> includes all discrete
   * functions to be interpolated onto the new grid.
   *
   * As we need two different kinds of pointers (<tt>vector<unsigned int> *</tt>
   * for the Dof indices and <tt>vector<VectorType> *</tt> for the
   * interpolated DoF values) we use the @p Pointerstruct that includes both of
   * these pointers and the pointer for each cell points to these @p
   * Pointerstructs. On each cell only one of the two different pointers is used
   * at one time hence we could use a <tt>void * pointer</tt> as
   * <tt>vector<unsigned int> *</tt> at one time and as
   * <tt>vector<VectorType> *</tt> at the other but using this @p
   * Pointerstruct in between makes the use of these pointers more safe and
   * gives better possibility to expand their usage.
   *
   * In <tt>interpolate(all_in, all_out)</tt> the refined cells are treated
   * according to the solution transfer while pure refinement. Additionally, on
   * each cell that is coarsened (hence previously was a parent cell), the
   * values of the discrete functions in @p all_out are set to the stored local
   * interpolated values that are accessible due to the 'vector<VectorType>
   * *' pointer in @p Pointerstruct that is pointed to by the pointer of that
   * cell. It is clear that <tt>interpolate(all_in, all_out)</tt> only can be
   * called with the <tt>vector<VectorType> all_in</tt> that previously was
   * the parameter of the <tt>prepare_for_coarsening_and_refinement(all_in)</tt>
   * function. Hence <tt>interpolate(all_in, all_out)</tt> can (in contrast to
   * <tt>refine_interpolate(in, out)</tt>) only be called once.
   * </ul>
   *
   *
   * <h3>Interaction with hanging nodes</h3>
   *
   * This class does its best to represent on the new mesh the finite element
   * function that existed on the old mesh, but this may lead to situations
   * where the function on the new mesh is no longer conforming at hanging
   * nodes. To this end, consider a situation of a twice refined mesh that
   * started with a single square cell (i.e., we now have 16 cells). Consider
   * also that we coarsen 4 of the cells back to the first refinement level. In
   * this case, we end up with a mesh that will look as follows if we were to
   * use a $Q_1$ element:
   *
   * @image html hanging_nodes.png ""
   *
   * The process of interpolating from the old to the new mesh would imply that
   * the values of the finite element function will not change on all of the
   * cells that remained as they are (i.e., the fine cells) but that on the
   * coarse cell at the top right, the four values at the vertices are obtained
   * by interpolating down from its former children.  If the original function
   * was not linear, this implies that the marked hanging nodes will retain
   * their old values which, in general, will not lead to a continuous function
   * along the corresponding edges. In other words, the solution vector obtained
   * after SolutionTransfer::interpolate() does not satisfy hanging node
   * constraints: it corresponds to the pointwise interpolation, but not to the
   * interpolation <i>onto the new finite element space that contains
   * constraints from hanging nodes</i>.
   *
   * Whether this is a problem you need to worry about or not depends on your
   * application. The situation is easily corrected, of course, by applying
   * AffineConstraints::distribute() to your solution vector after transfer,
   * using a constraints object computed on the new DoFHandler object (you
   * probably need to create this object anyway if you have hanging nodes). This
   * is also what is done, for example, in step-15.
   *
   * @note This situation can only happen if you do coarsening. If all cells
   * remain as they are or are refined, then SolutionTransfer::interpolate()
   * computes a new vector of nodal values, but the function represented is of
   * course exactly the same because the old finite element space is a subspace
   * of the new one. Thus, if the old function was conforming (i.e., satisfied
   * hanging node constraints), then so does the new one, and it is not
   * necessary to call AffineConstraints::distribute().
   *
   *
   * <h3>Implementation in the context of hp-finite elements</h3>
   *
   * In the case of DoFHandlers with hp-capabilities, nothing defines which of
   * the finite elements that are part of the hp::FECollection associated with
   * the DoFHandler, should be considered on cells that are not active (i.e.,
   * that have children). This is because degrees of freedom are only allocated
   * for active cells and, in fact, it is not allowed to set an active FE index
   * on non-active cells using DoFAccessor::set_active_fe_index().
   *
   * It is, thus, not entirely natural what should happen if, for example, a few
   * cells are coarsened away. This class then implements the following
   * algorithm:
   * - If a cell is refined, then the values of the solution vector(s) are
   *   interpolated before refinement on the to-be-refined cell from the space
   * of the active finite element to the one of the future finite element. These
   *   values are then distributed on the finite element spaces of the children
   *   post-refinement. This may lose information if, for example, the old cell
   *   used a Q2 space and the children use Q1 spaces, or the information may be
   *   prolonged if the parent cell used a Q1 space and the children are Q2s.
   * - If cells are to be coarsened, then the values from the child cells are
   *   interpolated to the parent cell using the largest of the child cell
   * future finite element spaces, which will be identified as the least
   * dominant element following the FiniteElementDomination logic (consult
   *   hp::FECollection::find_dominated_fe_extended() for more information). For
   *   example, if the children of a cell use Q1, Q2 and Q3 spaces, then the
   *   values from the children are interpolated into a Q3 space on the parent
   *   cell. After refinement, this Q3 function on the parent cell is then
   *   interpolated into the space the user has selected for this cell (which
   * may be different from Q3, in this example, if the user has set the active
   * FE index for a different space post-refinement and before calling
   *   DoFHandler::distribute_dofs()).
   *
   * @note In the context of hp-refinement, if cells are coarsened or the
   * polynomial degree is lowered on some cells, then the old finite element
   * space is not a subspace of the new space and you may run into the same
   * situation as discussed above with hanging nodes. You may want to consider
   * calling AffineConstraints::distribute() on the vector obtained by
   * transferring the solution.
   *
   * @ingroup numerics
   *
   * @deprecated Use dealii::SolutionTransfer instead.
   */
  template <int dim, typename VectorType = Vector<double>, int spacedim = dim>
  class SolutionTransfer
  {
  public:
    /**
     * Constructor, takes the current DoFHandler as argument.
     */
    SolutionTransfer(const DoFHandler<dim, spacedim> &dof);

    /**
     * Destructor
     */
    ~SolutionTransfer();

    /**
     * Reinit this class to the state that it has directly after calling the
     * Constructor
     */
    void
    clear();

    /**
     * Prepares the @p SolutionTransfer for pure refinement. It stores the dof
     * indices of each cell. After calling this function only calling the @p
     * refine_interpolate functions is allowed.
     */
    void
    prepare_for_pure_refinement();

    /**
     * Prepares the @p SolutionTransfer for coarsening and refinement. It stores
     * the dof indices of each cell and stores the dof values of the vectors in
     * @p all_in in each cell that'll be coarsened. @p all_in includes all
     * vectors that are to be interpolated onto the new (refined and/or
     * coarsenend) grid.
     */
    void
    prepare_for_coarsening_and_refinement(
      const std::vector<VectorType> &all_in);

    /**
     * Same as previous function but for only one discrete function to be
     * interpolated.
     */
    void
    prepare_for_coarsening_and_refinement(const VectorType &in);

    /**
     * This function interpolates the discrete function @p in, which is a vector
     * on the grid before the refinement, to the function @p out which then is a
     * vector on the refined grid. It assumes the vectors having the right sizes
     * (i.e. <tt>in.size()==n_dofs_old</tt>,
     * <tt>out.size()==n_dofs_refined</tt>)
     *
     * Calling this function is allowed only if @p prepare_for_pure_refinement
     * is called and the refinement is executed before. Multiple calling of this
     * function is allowed. e.g. for interpolating several functions.
     */
    void
    refine_interpolate(const VectorType &in, VectorType &out) const;

    /**
     * This function interpolates the discrete functions that are stored in @p
     * all_in onto the refined and/or coarsenend grid. It assumes the vectors in
     * @p all_in denote the same vectors as in @p all_in as parameter of
     * <tt>prepare_for_refinement_and_coarsening(all_in)</tt>. However, there is
     * no way of verifying this internally, so be careful here.
     *
     * Calling this function is allowed only if first
     * Triangulation::prepare_coarsening_and_refinement, second @p
     * SolutionTransfer::prepare_for_coarsening_and_refinement, an then third
     * Triangulation::execute_coarsening_and_refinement are called before.
     * Multiple calling of this function is NOT allowed. Interpolating several
     * functions can be performed in one step.
     *
     * The number of output vectors is assumed to be the same as the number of
     * input vectors. Also, the sizes of the output vectors are assumed to be of
     * the right size (@p n_dofs_refined). Otherwise an assertion will be
     * thrown.
     */
    void
    interpolate(const std::vector<VectorType> &all_in,
                std::vector<VectorType>       &all_out) const;

    /**
     * Same as the previous function. It interpolates only one function. It
     * assumes the vectors having the right sizes (i.e.
     * <tt>in.size()==n_dofs_old</tt>, <tt>out.size()==n_dofs_refined</tt>)
     *
     * Multiple calling of this function is NOT allowed. Interpolating several
     * functions can be performed in one step by using <tt>interpolate (all_in,
     * all_out)</tt>
     */
    void
    interpolate(const VectorType &in, VectorType &out) const;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t
    memory_consumption() const;

    /**
     * Exception
     */
    DeclExceptionMsg(
      ExcNotPrepared,
      "You are attempting an operation for which this object is "
      "not prepared. This may be because you either did not call "
      "one of the prepare_*() functions at all, or because you "
      "called the wrong one for the operation you are currently "
      "attempting.");

    /**
     * Exception
     */
    DeclExceptionMsg(
      ExcAlreadyPrepForRef,
      "You are attempting to call one of the prepare_*() functions "
      "of this object to prepare it for an operation for which it "
      "is already prepared. Specifically, the object was "
      "previously prepared for pure refinement.");

    /**
     * Exception
     */
    DeclExceptionMsg(
      ExcAlreadyPrepForCoarseAndRef,
      "You are attempting to call one of the prepare_*() functions "
      "of this object to prepare it for an operation for which it "
      "is already prepared. Specifically, the object was "
      "previously prepared for both coarsening and refinement.");

  private:
    /**
     * Pointer to the degree of freedom handler to work with.
     */
    ObserverPointer<const DoFHandler<dim, spacedim>,
                    SolutionTransfer<dim, VectorType, spacedim>>
      dof_handler;

    /**
     * Stores the number of DoFs before the refinement and/or coarsening.
     */
    types::global_dof_index n_dofs_old;

    /**
     * Declaration of @p PreparationState that denotes the three possible states
     * of the @p SolutionTransfer: being prepared for 'pure refinement',
     * prepared for 'coarsening and refinement' or not prepared.
     */
    enum PreparationState
    {
      /**
       * The SolutionTransfer is not yet prepared.
       */
      none,
      /**
       * The SolutionTransfer is prepared for purely refinement.
       */
      pure_refinement,
      /**
       * The SolutionTransfer is prepared for coarsening and refinement.
       */
      coarsening_and_refinement
    };

    /**
     * Definition of the respective variable.
     */
    PreparationState prepared_for;


    /**
     * Is used for @p prepare_for_pure_refinement (of course also for
     * @p prepare_for_coarsening_and_refinement) and stores all dof indices of
     * the cells that'll be refined.
     */
    std::vector<std::vector<types::global_dof_index>> indices_on_cell;

    /**
     * All cell data (the dof indices and the dof values) should be accessible
     * from each cell. As each cell has got only one @p user_pointer, multiple
     * pointers to the data need to be packetized in a structure. Note that in
     * our case on each cell either the <tt>vector<unsigned int> indices</tt>
     * (if the cell will be refined) or the <tt>vector<double> dof_values</tt>
     * (if the children of this cell will be deleted) is needed, hence one @p
     * user_pointer should be sufficient, but to allow some error checks and to
     * preserve the user from making user errors the @p user_pointer will be
     * 'multiplied' by this structure.
     */
    struct Pointerstruct
    {
      Pointerstruct()
        : indices_ptr(nullptr)
        , dof_values_ptr(nullptr)
        , active_fe_index(0)
      {}
      Pointerstruct(std::vector<types::global_dof_index> *indices_ptr_in,
                    const unsigned int active_fe_index_in = 0)
        : indices_ptr(indices_ptr_in)
        , dof_values_ptr(nullptr)
        , active_fe_index(active_fe_index_in)
      {}
      Pointerstruct(
        std::vector<Vector<typename VectorType::value_type>> *dof_values_ptr_in,
        const unsigned int active_fe_index_in = 0)
        : indices_ptr(nullptr)
        , dof_values_ptr(dof_values_ptr_in)
        , active_fe_index(active_fe_index_in)
      {}
      std::size_t
      memory_consumption() const;

      std::vector<types::global_dof_index>                 *indices_ptr;
      std::vector<Vector<typename VectorType::value_type>> *dof_values_ptr;
      unsigned int                                          active_fe_index;
    };

    /**
     * Map mapping from level and index of cell to the @p Pointerstructs (cf.
     * there). This map makes it possible to keep all the information needed to
     * transfer the solution inside this object rather than using user pointers
     * of the Triangulation for this purpose.
     */
    std::map<std::pair<unsigned int, unsigned int>, Pointerstruct> cell_map;

    /**
     * Is used for @p prepare_for_coarsening_and_refinement. The interpolated
     * dof values of all cells that'll be coarsened will be stored in this
     * vector.
     */
    std::vector<std::vector<Vector<typename VectorType::value_type>>>
      dof_values_on_cell;
  };
} // namespace Legacy

DEAL_II_NAMESPACE_CLOSE

#endif
