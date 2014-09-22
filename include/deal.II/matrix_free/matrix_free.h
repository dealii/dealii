// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2014 by the deal.II authors
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


#ifndef __deal2__matrix_free_h
#define __deal2__matrix_free_h

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/block_vector_base.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/matrix_free/helper_functions.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/dof_info.h>
#include <deal.II/matrix_free/mapping_info.h>

#ifdef DEAL_II_WITH_THREADS
#include <tbb/task.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#endif

#include <stdlib.h>
#include <memory>
#include <limits>


DEAL_II_NAMESPACE_OPEN



/**
 * This class collects all the data that is stored for the matrix free
 * implementation. The storage scheme is tailored towards several loops
 * performed with the same data, i.e., typically doing many matrix-vector
 * products or residual computations on the same mesh. The class is used in
 * step-37 and step-48.
 *
 * This class does not implement any operations involving finite element basis
 * functions, i.e., regarding the operation performed on the cells. For these
 * operations, the class FEEvaluation is designed to use the data collected in
 * this class.
 *
 * The stored data can be subdivided into three main components:
 *
 * - DoFInfo: It stores how local degrees of freedom relate to global degrees
 *   of freedom. It includes a description of constraints that are evaluated
 *   as going through all local degrees of freedom on a cell.
 *
 * - MappingInfo: It stores the transformations from real to unit cells that
 *   are necessary in order to build derivatives of finite element functions
 *   and find location of quadrature weights in physical space.
 *
 * - ShapeInfo: It contains the shape functions of the finite element,
 *   evaluated on the unit cell.
 *
 * Besides the initialization routines, this class implements only a
 * single operation, namely a loop over all cells (cell_loop()). This
 * loop is scheduled in such a way that cells that share degrees of
 * freedom are not worked on simultaneously, which implies that it is
 * possible to write to vectors (or matrices) in parallel without
 * having to explicitly synchronize access to these vectors and
 * matrices. This class does not implement any shape values, all it
 * does is to cache the respective data. To implement finite element
 * operations, use the class FEEvaluation (or some of the related
 * classes).
 *
 * This class traverses the cells in a different order than the usual
 * Triangulation class in deal.II, in order to be flexible with respect to
 * parallelization in shared memory and vectorization.
 *
 * Vectorization is implemented by merging several topological cells into one
 * so-called macro cell. This enables the application of all cell-related
 * operations for several cells with one CPU instruction and is one of the
 * main features of this framework.
 *
 * @author Katharina Kormann, Martin Kronbichler, 2010, 2011
 */

template <int dim, typename Number=double>
class MatrixFree
{
public:

  /**
   * Collects the options for initialization of the MatrixFree class. The
   * first parameter specifies the MPI communicator to be used, the second the
   * parallelization options in shared memory (task-based parallelism, where
   * one can choose between no parallelism and three schemes that avoid that
   * cells with access to the same vector entries are accessed
   * simultaneously), the third with the block size for task parallel
   * scheduling, the fourth the update flags that should be stored by this
   * class.
   *
   * The fifth parameter specifies the level in the triangulation from which
   * the indices are to be used. If the level is set to
   * numbers::invalid_unsigned_int, the active cells are traversed, and
   * otherwise the cells in the given level. This option has no effect in case
   * a DoFHandler or hp::DoFHandler is given.
   *
   * The parameter @p initialize_plain_indices indicates whether the DoFInfo
   * class should also allow for access to vectors without resolving
   * constraints.
   *
   * The last two parameters allow the user to disable some of the
   * initialization processes. For example, if only the scheduling that avoids
   * touching the same vector/matrix indices simultaneously is to be found,
   * the mapping needs not be initialized. Likewise, if the mapping has
   * changed from one iteration to the next but the topology has not (like
   * when using a deforming mesh with MappingQEulerian), it suffices to
   * initialize the mapping only.
   */
  struct AdditionalData
  {
    /**
     * Collects options for task parallelism.
     */
    enum TasksParallelScheme {none, partition_partition, partition_color, color};

    /**
     * Constructor for AdditionalData.
     */
    AdditionalData (const MPI_Comm            mpi_communicator   = MPI_COMM_SELF,
                    const TasksParallelScheme tasks_parallel_scheme = partition_partition,
                    const unsigned int        tasks_block_size   = 0,
                    const UpdateFlags         mapping_update_flags  = update_gradients | update_JxW_values,
                    const unsigned int level_mg_handler = numbers::invalid_unsigned_int,
                    const bool                store_plain_indices = true,
                    const bool                initialize_indices = true,
                    const bool                initialize_mapping = true)
      :
      mpi_communicator      (mpi_communicator),
      tasks_parallel_scheme (tasks_parallel_scheme),
      tasks_block_size      (tasks_block_size),
      mapping_update_flags  (mapping_update_flags),
      level_mg_handler      (level_mg_handler),
      store_plain_indices   (store_plain_indices),
      initialize_indices    (initialize_indices),
      initialize_mapping    (initialize_mapping)
    {};

    /**
     * Sets the MPI communicator that the parallel layout of the operator
     * should be based upon. Defaults to MPI_COMM_SELF, but should be set to a
     * communicator similar to the one used for a distributed triangulation in
     * order to inform this class over all cells that are present.
     */
    MPI_Comm            mpi_communicator;

    /**
     * Sets the scheme for task parallelism. There are four options
     * available. If set to @p none, the operator application is done in
     * serial without shared memory parallelism. If this class is used
     * together with MPI and MPI is also used for parallelism within the
     * nodes, this flag should be set to @p none. The default value is @p
     * partition_partition, i.e. we actually use multithreading with the first
     * option described below.
     *
     * The first option @p partition_partition is to partition the cells on
     * two levels in onion-skin-like partitions and forming chunks of
     * tasks_block_size after the partitioning. The partitioning finds sets of
     * independent cells that enable working in parallel without accessing the
     * same vector entries at the same time.
     *
     * The second option @p partition_color is to use a partition on the
     * global level and color cells within the partitions (where all chunks
     * within a color are independent). Here, the subdivision into chunks of
     * cells is done before the partitioning, which might give worse
     * partitions but better cache performance if degrees of freedom are not
     * renumbered.
     *
     * The third option @p color is to use a traditional algorithm of coloring
     * on the global level. This scheme is a special case of the second option
     * where only one partition is present. Note that for problems with
     * hanging nodes, there are quite many colors (50 or more in 3D), which
     * might degrade parallel performance (bad cache behavior, many
     * synchronization points).
     */
    TasksParallelScheme tasks_parallel_scheme;

    /**
     * Sets the number of so-called macro cells that should form one
     * partition. If zero size is given, the class tries to find a good size
     * for the blocks based on
     * multithread_info.n_threads() and the number of cells
     * present. Otherwise, the given number is used. If the given number is
     * larger than one third of the number of total cells, this means no
     * parallelism. Note that in the case vectorization is used, a macro cell
     * consists of more than one physical cell.
     */
    unsigned int        tasks_block_size;

    /**
     * This flag is used to determine which quantities should be cached. This
     * class can cache data needed for gradient computations (inverse
     * Jacobians), Jacobian determinants (JxW), quadrature points as well as
     * data for Hessians (derivative of Jacobians). By default, only data for
     * gradients and Jacobian determinants times quadrature weights, JxW, are
     * cached. If quadrature points or second derivatives are needed, they
     * must be specified by this field (even though second derivatives might
     * still be evaluated on Cartesian cells without this option set here,
     * since there the Jacobian describes the mapping completely).
     */
    UpdateFlags         mapping_update_flags;

    /**
     * This option can be used to define whether we work on a certain level of
     * the mesh, and not the active cells. If set to invalid_unsigned_int
     * (which is the default value), the active cells are gone through,
     * otherwise the level given by this parameter. Note that if you specify
     * to work on a level, its dofs must be distributed by using
     * <code>dof_handler.distribute_mg_dofs(fe);</code>.
     */
    unsigned int        level_mg_handler;

    /**
     * Controls whether to allow reading from vectors without resolving
     * constraints, i.e., just read the local values of the vector. By
     * default, this option is disabled, so if you want to use
     * FEEvaluationBase::read_dof_values_plain, this flag needs to be set.
     */
    bool                store_plain_indices;

    /**
     * Option to control whether the indices stored in the DoFHandler should
     * be read and the pattern for task parallelism should be set up in the
     * initialize method of MatrixFree. Defaults to true. Can be disabled in
     * case the mapping should be recomputed (e.g. when using a deforming mesh
     * described through MappingEulerian) but the topology of cells has
     * remained the same.
     */
    bool                initialize_indices;

    /**
     * Option to control whether the mapping information should be computed in
     * the initialize method of MatrixFree. Defaults to true. Can be disabled
     * when only some indices should be set up (e.g. when only a set of
     * independent cells should be computed).
     */
    bool                initialize_mapping;
  };

  /**
   * @name 1: Construction and initialization
   */
  //@{
  /**
   * Default empty constructor. Does nothing.
   */
  MatrixFree ();

  /**
  * Destructor.
  */
  ~MatrixFree();

  /**
  * Extracts the information needed to perform loops over cells. The
  * DoFHandler and ConstraintMatrix describe the layout of degrees of freedom,
  * the DoFHandler and the mapping describe the transformations from unit to
  * real cell, and the finite element underlying the DoFHandler together with
  * the quadrature formula describe the local operations. Note that the finite
  * element underlying the DoFHandler must either be scalar or contain several
  * copies of the same element. Mixing several different elements into one
  * FESystem is not allowed. In that case, use the initialization function
  * with several DoFHandler arguments.
  *
  * The @p IndexSet @p locally_owned_dofs is used to specify the parallel
  * partitioning with MPI. Usually, this needs not be specified, and the other
  * initialization function without and @p IndexSet description can be used,
  * which gets the partitioning information builtin into the DoFHandler.
  */
  template <typename DH, typename Quadrature>
  void reinit (const Mapping<dim>     &mapping,
               const DH               &dof_handler,
               const ConstraintMatrix &constraint,
               const IndexSet         &locally_owned_dofs,
               const Quadrature       &quad,
               const AdditionalData    additional_data = AdditionalData());

  /**
  * Initializes the data structures. Same as above, but with index set stored
  * in the DoFHandler for describing the locally owned degrees of freedom.
  */
  template <typename DH, typename Quadrature>
  void reinit (const Mapping<dim>     &mapping,
               const DH               &dof_handler,
               const ConstraintMatrix &constraint,
               const Quadrature       &quad,
               const AdditionalData    additional_data = AdditionalData());

  /**
  * Initializes the data structures. Same as above, but with mapping @p
  * MappingQ1.
  */
  template <typename DH, typename Quadrature>
  void reinit (const DH               &dof_handler,
               const ConstraintMatrix &constraint,
               const Quadrature       &quad,
               const AdditionalData    additional_data = AdditionalData());

  /**
  * Extracts the information needed to perform loops over cells. The
  * DoFHandler and ConstraintMatrix describe the layout of degrees of freedom,
  * the DoFHandler and the mapping describe the transformations from unit to
  * real cell, and the finite element underlying the DoFHandler together with
  * the quadrature formula describe the local operations. As opposed to the
  * scalar case treated with the other initialization functions, this function
  * allows for problems with two or more different finite elements. The
  * DoFHandlers to each element must be passed as pointers to the
  * initialization function. Note that the finite element underlying an
  * DoFHandler must either be scalar or contain several copies of the same
  * element. Mixing several different elements into one @p FE_System is not
  * allowed.
  *
  * This function also allows for using several quadrature formulas, e.g. when
  * the description contains independent integrations of elements of different
  * degrees. However, the number of different quadrature formulas can be sets
  * independently from the number of DoFHandlers, when several elements are
  * always integrated with the same quadrature formula.
  *
  * The @p IndexSet @p locally_owned_dofs is used to specify the parallel
  * partitioning with MPI. Usually, this needs not be specified, and the other
  * initialization function without and @p IndexSet description can be used,
  * which gets the partitioning information from the DoFHandler. This is the
  * most general initialization function.
  */
  template <typename DH, typename Quadrature>
  void reinit (const Mapping<dim>                         &mapping,
               const std::vector<const DH *>               &dof_handler,
               const std::vector<const ConstraintMatrix *> &constraint,
               const std::vector<IndexSet>                &locally_owned_set,
               const std::vector<Quadrature>              &quad,
               const AdditionalData                        additional_data = AdditionalData());

  /**
  * Initializes the data structures. Same as before, but now the index set
  * description of the locally owned range of degrees of freedom is taken from
  * the DoFHandler.
  */
  template <typename DH, typename Quadrature>
  void reinit (const Mapping<dim>                         &mapping,
               const std::vector<const DH *>               &dof_handler,
               const std::vector<const ConstraintMatrix *> &constraint,
               const std::vector<Quadrature>              &quad,
               const AdditionalData                        additional_data = AdditionalData());

  /**
  * Initializes the data structures. Same as above, but with mapping @p
  * MappingQ1.
  */
  template <typename DH, typename Quadrature>
  void reinit (const std::vector<const DH *>               &dof_handler,
               const std::vector<const ConstraintMatrix *> &constraint,
               const std::vector<Quadrature>              &quad,
               const AdditionalData                        additional_data = AdditionalData());

  /**
  * Initializes the data structures. Same as before, but now the index set
  * description of the locally owned range of degrees of freedom is taken from
  * the DoFHandler. Moreover, only a single quadrature formula is used, as
  * might be necessary when several components in a vector-valued problem are
  * integrated together based on the same quadrature formula.
  */
  template <typename DH, typename Quadrature>
  void reinit (const Mapping<dim>                         &mapping,
               const std::vector<const DH *>               &dof_handler,
               const std::vector<const ConstraintMatrix *> &constraint,
               const Quadrature                           &quad,
               const AdditionalData                        additional_data = AdditionalData());

  /**
  * Initializes the data structures. Same as above, but with mapping @p
  * MappingQ1.
  */
  template <typename DH, typename Quadrature>
  void reinit (const std::vector<const DH *>               &dof_handler,
               const std::vector<const ConstraintMatrix *> &constraint,
               const Quadrature                           &quad,
               const AdditionalData                        additional_data = AdditionalData());

  /**
   * Copy function. Creates a deep copy of all data structures. It is usually
   * enough to keep the data for different operations once, so this function
   * should not be needed very often.
   */
  void copy_from (const MatrixFree<dim,Number> &matrix_free_base);

  /**
   * Clears all data fields and brings the class into a condition similar to
   * after having called the default constructor.
   */
  void clear();

  //@}

  /**
   * @name 2: Loop over cells
   */
  //@{
  /**
   * This method runs the loop over all cells (in parallel) and performs the
   * MPI data exchange on the source vector and destination vector. The first
   * argument indicates a function object that has the following signature:
   * <code>cell_operation (const MatrixFree<dim,Number> &, OutVector &,
   * InVector &, std::pair<unsigned int,unsigned int> &)</code>, where the
   * first argument passes the data of the calling class and the last argument
   * defines the range of cells which should be worked on (typically more than
   * one cell should be worked on in order to reduce overheads).  One can pass
   * a pointer to an object in this place if it has an <code>operator()</code>
   * with the correct set of arguments since such a pointer can be converted
   * to the function object.
   */
  template <typename OutVector, typename InVector>
  void cell_loop (const std_cxx11::function<void (const MatrixFree<dim,Number> &,
                                                  OutVector &,
                                                  const InVector &,
                                                  const std::pair<unsigned int,
                                                  unsigned int> &)> &cell_operation,
                  OutVector      &dst,
                  const InVector &src) const;

  /**
   * This is the second variant to run the loop over all cells, now providing
   * a function pointer to a member function of class @p CLASS with the
   * signature <code>cell_operation (const MatrixFree<dim,Number> &, OutVector
   * &, InVector &, std::pair<unsigned int,unsigned int>&)const</code>. This
   * method obviates the need to call std_cxx11::bind to bind the class into
   * the given function in case the local function needs to access data in the
   * class (i.e., it is a non-static member function).
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void cell_loop (void (CLASS::*function_pointer)(const MatrixFree &,
                                                  OutVector &,
                                                  const InVector &,
                                                  const std::pair<unsigned int,
                                                  unsigned int> &)const,
                  const CLASS    *owning_class,
                  OutVector      &dst,
                  const InVector &src) const;

  /**
   * Same as above, but for class member functions which are non-const.
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void cell_loop (void (CLASS::*function_pointer)(const MatrixFree &,
                                                  OutVector &,
                                                  const InVector &,
                                                  const std::pair<unsigned int,
                                                  unsigned int> &),
                  CLASS          *owning_class,
                  OutVector      &dst,
                  const InVector &src) const;

  /**
   * In the hp adaptive case, a subrange of cells as computed during the cell
   * loop might contain elements of different degrees. Use this function to
   * compute what the subrange for an individual finite element degree is. The
   * finite element degree is associated to the vector component given in the
   * function call.
   */
  std::pair<unsigned int,unsigned int>
  create_cell_subrange_hp (const std::pair<unsigned int,unsigned int> &range,
                           const unsigned int fe_degree,
                           const unsigned int vector_component = 0) const;

  /**
   * In the hp adaptive case, a subrange of cells as computed during the cell
   * loop might contain elements of different degrees. Use this function to
   * compute what the subrange for a given index the hp finite element, as
   * opposed to the finite element degree in the other function.
   */
  std::pair<unsigned int,unsigned int>
  create_cell_subrange_hp_by_index (const std::pair<unsigned int,unsigned int> &range,
                                    const unsigned int fe_index,
                                    const unsigned int vector_component = 0) const;

  //@}

  /**
   * @name 3: Initialization of vectors
   */
  //@{
  /**
   * Initialize function for a general vector. The length of the vector is
   * equal to the total number of degrees in the DoFHandler. If the vector is
   * of class parallel::distributed::Vector@<Number@>, the ghost entries are
   * set accordingly. For vector-valued problems with several DoFHandlers
   * underlying this class, the parameter @p vector_component defines which
   * component is to be used.
   */
  template <typename VectorType>
  void initialize_dof_vector(VectorType &vec,
                             const unsigned int vector_component=0) const;

  /**
   * Initialize function for a distributed vector. The length of the vector is
   * equal to the total number of degrees in the DoFHandler. If the vector is
   * of class parallel::distributed::Vector@<Number@>, the ghost entries are
   * set accordingly. For vector-valued problems with several DoFHandlers
   * underlying this class, the parameter @p vector_component defines which
   * component is to be used.
   */
  template <typename Number2>
  void initialize_dof_vector(parallel::distributed::Vector<Number2> &vec,
                             const unsigned int vector_component=0) const;

  /**
   * Returns the partitioner that represents the locally owned data and the
   * ghost indices where access is needed to for the cell loop. The
   * partitioner is constructed from the locally owned dofs and ghost dofs
   * given by the respective fields. If you want to have specific information
   * about these objects, you can query them with the respective access
   * functions. If you just want to initialize a (parallel) vector, you should
   * usually prefer this data structure as the data exchange information can
   * be reused from one vector to another.
   */
  const std_cxx11::shared_ptr<const Utilities::MPI::Partitioner> &
  get_vector_partitioner (const unsigned int vector_component=0) const;

  /**
   * Returns the set of cells that are oned by the processor.
   */
  const IndexSet &
  get_locally_owned_set (const unsigned int fe_component = 0) const;

  /**
   * Returns the set of ghost cells needed but not owned by the processor.
   */
  const IndexSet &
  get_ghost_set (const unsigned int fe_component = 0) const;

  /**
   * Returns a list of all degrees of freedom that are constrained. The list
   * is returned in MPI-local index space for the locally owned range of the
   * vector, not in global MPI index space that spans all MPI processors. To
   * get numbers in global index space, call
   * <tt>get_vector_partitioner()->local_to_global</tt> on an entry of the
   * vector. In addition, it only returns the indices for degrees of freedom
   * that are owned locally, not for ghosts.
   */
  const std::vector<unsigned int> &
  get_constrained_dofs (const unsigned int fe_component = 0) const;

  /**
   * Calls renumber_dofs function in dof_info which renumbers the degrees
   * of freedom according to the ordering for parallelization.
   */
  void renumber_dofs (std::vector<types::global_dof_index> &renumbering,
                      const unsigned int vector_component = 0);

  //@}

  /**
   * @name 4: General information
   */
  //@{
  /**
   * Returns the number of different DoFHandlers specified at initialization.
   */
  unsigned int n_components () const;

  /**
   * Returns the number of cells this structure is based on. If you are using
   * a usual DoFHandler, it corresponds to the number of (locally owned)
   * active cells. Note that most data structures in this class do not
   * directly act on this number but rather on n_macro_cells() which gives the
   * number of cells as seen when lumping several cells together with
   * vectorization.
   */
  unsigned int n_physical_cells () const;

  /**
   * Returns the number of macro cells that this structure works on, i.e., the
   * number of cell chunks that are worked on after the application of
   * vectorization which in general works on several cells at once. The cell
   * range in @p cell_loop runs from zero to n_macro_cells() (exclusive), so
   * this is the appropriate size if you want to store arrays of data for all
   * cells to be worked on. This number is approximately
   * n_physical_cells()/VectorizedArray::n_array_elements (depending on how
   * many cell chunks that do not get filled up completely).
   */
  unsigned int n_macro_cells () const;

  /**
   * In case this structure was built based on a DoFHandler, this returns the
   * DoFHandler.
   */
  const DoFHandler<dim> &
  get_dof_handler (const unsigned int fe_component = 0) const;

  /**
   * This returns the cell iterator in deal.II speak to a given cell in the
   * renumbering of this structure.
   *
   * Note that the cell iterators in deal.II go through cells differently to
   * what the cell loop of this class does. This is because several cells are
   * worked on together (vectorization), and since cells with neighbors on
   * different MPI processors need to be accessed at a certain time when
   * accessing remote data and overlapping communication with computation.
   */
  typename DoFHandler<dim>::cell_iterator
  get_cell_iterator (const unsigned int macro_cell_number,
                     const unsigned int vector_number,
                     const unsigned int fe_component = 0) const;

  /**
   * This returns the cell iterator in deal.II speak to a given cell in the
   * renumbering of this structure. This function returns an exception in case
   * the structure was not constructed based on an hp::DoFHandler.
   *
   * Note that the cell iterators in deal.II go through cells differently to
   * what the cell loop of this class does. This is because several cells are
   * worked on together (vectorization), and since cells with neighbors on
   * different MPI processors need to be accessed at a certain time when
   * accessing remote data and overlapping communication with computation.
   */
  typename hp::DoFHandler<dim>::active_cell_iterator
  get_hp_cell_iterator (const unsigned int macro_cell_number,
                        const unsigned int vector_number,
                        const unsigned int fe_component = 0) const;

  /**
   * Since this class uses vectorized data types with usually more than one
   * value in the data field, a situation might occur when some components of
   * the vector type do not correspond to an actual cell in the mesh. When
   * using only this class, one usually does not need to bother about that
   * fact since the values are padded with zeros. However, when this class is
   * mixed with deal.II access to cells, care needs to be taken. This function
   * returns @p true if not all @p vectorization_length cells for the given @p
   * macro_cell are real cells. To find out how many cells are actually used,
   * use the function @p n_components_filled.
   */
  bool
  at_irregular_cell (const unsigned int macro_cell_number) const;

  /**
   * Use this function to find out how many cells over the length of
   * vectorization data types correspond to real cells in the mesh. For most
   * given @p macro_cells, this is just @p vectorization_length many, but
   * there might be one or a few meshes (where the numbers do not add up)
   * where there are less such components filled, indicated by the function @p
   * at_irregular_cell.
   */
  unsigned int
  n_components_filled (const unsigned int macro_cell_number) const;

  /**
   * Returns the number of degrees of freedom per cell for a given hp index.
   */
  unsigned int
  get_dofs_per_cell (const unsigned int fe_component = 0,
                     const unsigned int hp_active_fe_index = 0) const;

  /**
   * Returns the number of quadrature points per cell for a given hp index.
   */
  unsigned int
  get_n_q_points (const unsigned int quad_index = 0,
                  const unsigned int hp_active_fe_index = 0) const;

  /**
   * Returns the number of degrees of freedom on each face of the cell for
   * given hp index.
   */
  unsigned int
  get_dofs_per_face (const unsigned int fe_component = 0,
                     const unsigned int hp_active_fe_index = 0) const;

  /**
   * Returns the number of quadrature points on each face of the cell for
   * given hp index.
   */
  unsigned int
  get_n_q_points_face (const unsigned int quad_index = 0,
                       const unsigned int hp_active_fe_index = 0) const;

  /**
   * Returns the quadrature rule for given hp index.
   */
  const Quadrature<dim> &
  get_quadrature (const unsigned int quad_index = 0,
                  const unsigned int hp_active_fe_index = 0) const;

  /**
   * Returns the quadrature rule for given hp index.
   */
  const Quadrature<dim-1> &
  get_face_quadrature (const unsigned int quad_index = 0,
                       const unsigned int hp_active_fe_index = 0) const;

  /**
   * Queries whether or not the indexation has been set.
   */
  bool indices_initialized () const;

  /**
   * Queries whether or not the geometry-related information for the cells has
   * been set.
   */

  bool mapping_initialized () const;

  /**
   * Returns an approximation of the memory consumption of this class in
   * bytes.
   */
  std::size_t memory_consumption() const;

  /**
   * Prints a detailed summary of memory consumption in the different
   * structures of this class to the given output stream.
   */
  template <typename STREAM>
  void print_memory_consumption(STREAM &out) const;

  /**
   * Prints a summary of this class to the given output stream. It is focused
   * on the indices, and does not print all the data stored.
   */
  void print (std::ostream &out) const;

  //@}

  /**
   * @name 5: Access of internal data structure (expert mode)
   */
  //@{
  /**
   * Returns information on task graph.
   */
  const internal::MatrixFreeFunctions::TaskInfo &
  get_task_info () const;

  /**
   * Returns information on system size.
   */
  const internal::MatrixFreeFunctions::SizeInfo &
  get_size_info () const;

  /*
   * Returns geometry-dependent information on the cells.
   */
  const internal::MatrixFreeFunctions::MappingInfo<dim,Number> &
  get_mapping_info () const;

  /**
   * Returns information on indexation degrees of freedom.
   */
  const internal::MatrixFreeFunctions::DoFInfo &
  get_dof_info (const unsigned int fe_component = 0) const;

  /**
   * Returns the number of weights in the constraint pool.
   */
  unsigned int n_constraint_pool_entries() const;

  /**
   * Returns a pointer to the first number in the constraint pool data with
   * index @p pool_index (to be used together with @p constraint_pool_end()).
   */
  const Number *
  constraint_pool_begin (const unsigned int pool_index) const;

  /**
   * Returns a pointer to one past the last number in the constraint pool data
   * with index @p pool_index (to be used together with @p
   * constraint_pool_begin()).
   */
  const Number *
  constraint_pool_end (const unsigned int pool_index) const;

  /**
   * Returns the unit cell information for given hp index.
   */
  const internal::MatrixFreeFunctions::ShapeInfo<Number> &
  get_shape_info (const unsigned int fe_component = 0,
                  const unsigned int quad_index   = 0,
                  const unsigned int hp_active_fe_index = 0,
                  const unsigned int hp_active_quad_index = 0) const;

  //@}

private:

  /**
   * This is the actual reinit function that sets up the indices for the
   * DoFHandler case.
   */
  void internal_reinit (const Mapping<dim>                &mapping,
                        const std::vector<const DoFHandler<dim> *> &dof_handler,
                        const std::vector<const ConstraintMatrix *> &constraint,
                        const std::vector<IndexSet>       &locally_owned_set,
                        const std::vector<hp::QCollection<1> > &quad,
                        const AdditionalData               additional_data);

  /**
   * Same as before but for hp::DoFHandler instead of generic DoFHandler type.
   */
  void internal_reinit (const Mapping<dim>               &mapping,
                        const std::vector<const hp::DoFHandler<dim>*> &dof_handler,
                        const std::vector<const ConstraintMatrix *> &constraint,
                        const std::vector<IndexSet>      &locally_owned_set,
                        const std::vector<hp::QCollection<1> > &quad,
                        const AdditionalData              additional_data);

  /**
   * Initializes the fields in DoFInfo together with the constraint pool that
   * holds all different weights in the constraints (not part of DoFInfo
   * because several DoFInfo classes can have the same weights which
   * consequently only need to be stored once).
   */
  void
  initialize_indices (const std::vector<const ConstraintMatrix *> &constraint,
                      const std::vector<IndexSet> &locally_owned_set);

  /**
   * Initializes the DoFHandlers based on a DoFHandler<dim> argument.
   */
  void initialize_dof_handlers (const std::vector<const DoFHandler<dim>*> &dof_handlers,
                                const unsigned int                         level);

  /**
   * Initializes the DoFHandlers based on a hp::DoFHandler<dim> argument.
   */
  void initialize_dof_handlers (const std::vector<const hp::DoFHandler<dim>*> &dof_handlers,
                                const unsigned int                             level);

  /**
   * This struct defines which DoFHandler has actually been given at
   * construction, in order to define the correct behavior when querying the
   * underlying DoFHandler.
   */
  struct DoFHandlers
  {
    DoFHandlers () : n_dof_handlers (0), level (numbers::invalid_unsigned_int) {};
    std::vector<SmartPointer<const DoFHandler<dim> > >   dof_handler;
    std::vector<SmartPointer<const hp::DoFHandler<dim> > > hp_dof_handler;
    enum ActiveDoFHandler { usual, hp } active_dof_handler;
    unsigned int n_dof_handlers;
    unsigned int level;
  };

  /**
   * Pointers to the DoFHandlers underlying the current problem.
   */
  DoFHandlers dof_handlers;

  /**
   * Contains the information about degrees of freedom on the individual cells
   * and constraints.
   */
  std::vector<internal::MatrixFreeFunctions::DoFInfo> dof_info;

  /**
   * Contains the weights for constraints stored in DoFInfo. Filled into a
   * separate field since several vector components might share similar
   * weights, which reduces memory consumption. Moreover, it obviates template
   * arguments on DoFInfo and keeps it a plain field of indices only.
   */
  std::vector<Number> constraint_pool_data;

  /**
   * Contains an indicator to the start of the ith index in the constraint
   * pool data.
   */
  std::vector<unsigned int> constraint_pool_row_index;

  /**
   * Holds information on transformation of cells from reference cell to real
   * cell that is needed for evaluating integrals.
   */
  internal::MatrixFreeFunctions::MappingInfo<dim,Number> mapping_info;

  /**
   * Contains shape value information on the unit cell.
   */
  Table<4,internal::MatrixFreeFunctions::ShapeInfo<Number> > shape_info;

  /**
   * Describes how the cells are gone through. With the cell level (first
   * index in this field) and the index within the level, one can reconstruct
   * a deal.II cell iterator and use all the traditional things deal.II offers
   * to do with cell iterators.
   */
  std::vector<std::pair<unsigned int,unsigned int> > cell_level_index;

  /**
   * Stores how many cells we have, how many cells that we see after applying
   * vectorization (i.e., the number of macro cells), and MPI-related stuff.
   */
  internal::MatrixFreeFunctions::SizeInfo size_info;

  /**
   * Information regarding the shared memory parallelization.
   */
  internal::MatrixFreeFunctions::TaskInfo task_info;

  /**
   * Stores whether indices have been initialized.
   */
  bool                               indices_are_initialized;

  /**
   * Stores whether indices have been initialized.
   */
  bool                               mapping_is_initialized;
};



/*----------------------- Inline functions ----------------------------------*/

#ifndef DOXYGEN


template <int dim, typename Number>
template <typename VectorType>
inline
void
MatrixFree<dim,Number>::initialize_dof_vector(VectorType &vec,
                                              const unsigned int comp) const
{
  AssertIndexRange (comp, n_components());
  vec.reinit(dof_info[comp].vector_partitioner->size());
}



template <int dim, typename Number>
template <typename Number2>
inline
void
MatrixFree<dim,Number>::initialize_dof_vector(parallel::distributed::Vector<Number2> &vec,
                                              const unsigned int comp) const
{
  AssertIndexRange (comp, n_components());
  vec.reinit(dof_info[comp].vector_partitioner);
}



template <int dim, typename Number>
inline
const std_cxx11::shared_ptr<const Utilities::MPI::Partitioner> &
MatrixFree<dim,Number>::get_vector_partitioner (const unsigned int comp) const
{
  AssertIndexRange (comp, n_components());
  return dof_info[comp].vector_partitioner;
}



template <int dim, typename Number>
inline
const std::vector<unsigned int> &
MatrixFree<dim,Number>::get_constrained_dofs (const unsigned int comp) const
{
  AssertIndexRange (comp, n_components());
  return dof_info[comp].constrained_dofs;
}



template <int dim, typename Number>
inline
unsigned int
MatrixFree<dim,Number>::n_components () const
{
  AssertDimension (dof_handlers.n_dof_handlers, dof_info.size());
  return dof_handlers.n_dof_handlers;
}



template <int dim, typename Number>
inline
const internal::MatrixFreeFunctions::TaskInfo &
MatrixFree<dim,Number>::get_task_info () const
{
  return task_info;
}



template <int dim, typename Number>
inline
const internal::MatrixFreeFunctions::SizeInfo &
MatrixFree<dim,Number>::get_size_info () const
{
  return size_info;
}



template <int dim, typename Number>
inline
unsigned int
MatrixFree<dim,Number>::n_macro_cells () const
{
  return size_info.n_macro_cells;
}



template <int dim, typename Number>
inline
unsigned int
MatrixFree<dim,Number>::n_physical_cells () const
{
  return size_info.n_active_cells;
}



template <int dim, typename Number>
inline
const internal::MatrixFreeFunctions::MappingInfo<dim,Number> &
MatrixFree<dim,Number>::get_mapping_info () const
{
  return mapping_info;
}



template <int dim, typename Number>
inline
const internal::MatrixFreeFunctions::DoFInfo &
MatrixFree<dim,Number>::get_dof_info (unsigned int dof_index) const
{
  AssertIndexRange (dof_index, n_components());
  return dof_info[dof_index];
}



template <int dim, typename Number>
inline
unsigned int
MatrixFree<dim,Number>::n_constraint_pool_entries() const
{
  return constraint_pool_row_index.size()-1;
}



template <int dim, typename Number>
inline
const Number *
MatrixFree<dim,Number>::constraint_pool_begin (const unsigned int row) const
{
  AssertIndexRange (row, constraint_pool_row_index.size()-1);
  return constraint_pool_data.empty() ? 0 :
         &constraint_pool_data[0] + constraint_pool_row_index[row];
}



template <int dim, typename Number>
inline
const Number *
MatrixFree<dim,Number>::constraint_pool_end (const unsigned int row) const
{
  AssertIndexRange (row, constraint_pool_row_index.size()-1);
  return constraint_pool_data.empty() ? 0 :
         &constraint_pool_data[0] + constraint_pool_row_index[row+1];
}



template <int dim, typename Number>
inline
std::pair<unsigned int,unsigned int>
MatrixFree<dim,Number>::create_cell_subrange_hp
(const std::pair<unsigned int,unsigned int> &range,
 const unsigned int degree,
 const unsigned int vector_component) const
{
  AssertIndexRange (vector_component, dof_info.size());
  if (dof_info[vector_component].cell_active_fe_index.empty())
    {
      AssertDimension (dof_info[vector_component].fe_index_conversion.size(),1);
      if (dof_info[vector_component].fe_index_conversion[0].first == degree)
        return range;
      else
        return std::pair<unsigned int,unsigned int> (range.second,range.second);
    }

  const unsigned int fe_index =
    dof_info[vector_component].fe_index_from_degree(degree);
  if (fe_index >= dof_info[vector_component].max_fe_index)
    return std::pair<unsigned int,unsigned int>(range.second, range.second);
  else
    return create_cell_subrange_hp_by_index (range, fe_index, vector_component);
}



template <int dim, typename Number>
inline
std::pair<unsigned int,unsigned int>
MatrixFree<dim,Number>::create_cell_subrange_hp_by_index
(const std::pair<unsigned int,unsigned int> &range,
 const unsigned int fe_index,
 const unsigned int vector_component) const
{
  AssertIndexRange (fe_index, dof_info[vector_component].max_fe_index);
  const std::vector<unsigned int> &fe_indices =
    dof_info[vector_component].cell_active_fe_index;
  if (fe_indices.size() == 0)
    return range;
  else
    {
      // the range over which we are searching must be ordered, otherwise we
      // got a range that spans over too many cells
#ifdef DEBUG
      for (unsigned int i=range.first+1; i<range.second; ++i)
        Assert (fe_indices[i] >= fe_indices[i-1],
                ExcMessage ("Cell range must be over sorted range of fe indices in hp case!"));
      AssertIndexRange(range.first,fe_indices.size()+1);
      AssertIndexRange(range.second,fe_indices.size()+1);
#endif
      std::pair<unsigned int,unsigned int> return_range;
      return_range.first =
        std::lower_bound(&fe_indices[0] + range.first,
                         &fe_indices[0] + range.second, fe_index)
        -&fe_indices[0] ;
      return_range.second =
        std::lower_bound(&fe_indices[0] + return_range.first,
                         &fe_indices[0] + range.second,
                         fe_index + 1)-&fe_indices[0];
      Assert(return_range.first >= range.first &&
             return_range.second <= range.second, ExcInternalError());
      return return_range;
    }
}



template <int dim, typename Number>
inline
void
MatrixFree<dim,Number>::renumber_dofs (std::vector<types::global_dof_index> &renumbering,
                                       const unsigned int vector_component)
{
  AssertIndexRange(vector_component, dof_info.size());
  dof_info[vector_component].renumber_dofs (renumbering);
}



template <int dim, typename Number>
inline
const DoFHandler<dim> &
MatrixFree<dim,Number>::get_dof_handler (const unsigned int dof_index) const
{
  AssertIndexRange (dof_index, n_components());
  if (dof_handlers.active_dof_handler == DoFHandlers::usual)
    {
      AssertDimension (dof_handlers.dof_handler.size(),
                       dof_handlers.n_dof_handlers);
      return *dof_handlers.dof_handler[dof_index];
    }
  else
    {
      Assert (false, ExcNotImplemented());
      // put pseudo return argument to avoid compiler error, but trigger a
      // segfault in case this is only run in optimized mode
      return *dof_handlers.dof_handler[numbers::invalid_unsigned_int];
    }
}



template <int dim, typename Number>
inline
typename DoFHandler<dim>::cell_iterator
MatrixFree<dim,Number>::get_cell_iterator(const unsigned int macro_cell_number,
                                          const unsigned int vector_number,
                                          const unsigned int dof_index) const
{
  const unsigned int vectorization_length=VectorizedArray<Number>::n_array_elements;
#ifdef DEBUG
  AssertIndexRange (dof_index, dof_handlers.n_dof_handlers);
  AssertIndexRange (macro_cell_number, size_info.n_macro_cells);
  AssertIndexRange (vector_number, vectorization_length);
  const unsigned int irreg_filled = dof_info[dof_index].row_starts[macro_cell_number][2];
  if (irreg_filled > 0)
    AssertIndexRange (vector_number, irreg_filled);
#endif

  const DoFHandler<dim> *dofh = 0;
  if (dof_handlers.active_dof_handler == DoFHandlers::usual)
    {
      AssertDimension (dof_handlers.dof_handler.size(),
                       dof_handlers.n_dof_handlers);
      dofh = dof_handlers.dof_handler[dof_index];
    }
  else
    {
      Assert (false, ExcMessage ("Cannot return DoFHandler<dim>::cell_iterator "
                                 "for underlying DoFHandler!"));
    }

  std::pair<unsigned int,unsigned int> index =
    cell_level_index[macro_cell_number*vectorization_length+vector_number];
  return typename DoFHandler<dim>::cell_iterator
         (&dofh->get_tria(), index.first, index.second, dofh);
}



template <int dim, typename Number>
inline
typename hp::DoFHandler<dim>::active_cell_iterator
MatrixFree<dim,Number>::get_hp_cell_iterator(const unsigned int macro_cell_number,
                                             const unsigned int vector_number,
                                             const unsigned int dof_index) const
{
  const unsigned int vectorization_length=VectorizedArray<Number>::n_array_elements;
#ifdef DEBUG
  AssertIndexRange (dof_index, dof_handlers.n_dof_handlers);
  AssertIndexRange (macro_cell_number, size_info.n_macro_cells);
  AssertIndexRange (vector_number, vectorization_length);
  const unsigned int irreg_filled = dof_info[dof_index].row_starts[macro_cell_number][2];
  if (irreg_filled > 0)
    AssertIndexRange (vector_number, irreg_filled);
#endif

  Assert (dof_handlers.active_dof_handler == DoFHandlers::hp,
          ExcNotImplemented());
  const hp::DoFHandler<dim> *dofh = dof_handlers.hp_dof_handler[dof_index];
  std::pair<unsigned int,unsigned int> index =
    cell_level_index[macro_cell_number*vectorization_length+vector_number];
  return typename hp::DoFHandler<dim>::cell_iterator
         (&dofh->get_tria(), index.first, index.second, dofh);
}



template <int dim, typename Number>
inline
bool
MatrixFree<dim,Number>::at_irregular_cell (const unsigned int macro_cell) const
{
  AssertIndexRange (macro_cell, size_info.n_macro_cells);
  return dof_info[0].row_starts[macro_cell][2] > 0;
}



template <int dim, typename Number>
inline
unsigned int
MatrixFree<dim,Number>::n_components_filled (const unsigned int macro_cell) const
{
  AssertIndexRange (macro_cell, size_info.n_macro_cells);
  const unsigned int n_filled = dof_info[0].row_starts[macro_cell][2];
  if (n_filled == 0)
    return VectorizedArray<Number>::n_array_elements;
  else
    return n_filled;
}



template <int dim, typename Number>
inline
unsigned int
MatrixFree<dim,Number>::get_dofs_per_cell(const unsigned int dof_index,
                                          const unsigned int active_fe_index) const
{
  AssertIndexRange (dof_index, dof_info.size());
  return dof_info[dof_index].dofs_per_cell[active_fe_index];
}



template <int dim, typename Number>
inline
unsigned int
MatrixFree<dim,Number>::get_n_q_points(const unsigned int quad_index,
                                       const unsigned int active_fe_index) const
{
  AssertIndexRange (quad_index,
                    mapping_info.mapping_data_gen.size());
  return mapping_info.mapping_data_gen[quad_index].n_q_points[active_fe_index];
}



template <int dim, typename Number>
inline
unsigned int
MatrixFree<dim,Number>::get_dofs_per_face(const unsigned int dof_index,
                                          const unsigned int active_fe_index) const
{
  AssertIndexRange (dof_index, dof_info.size());
  return dof_info[dof_index].dofs_per_face[active_fe_index];
}



template <int dim, typename Number>
inline
unsigned int
MatrixFree<dim,Number>::get_n_q_points_face(const unsigned int quad_index,
                                            const unsigned int active_fe_index) const
{
  AssertIndexRange (quad_index,
                    mapping_info.mapping_data_gen.size());
  return mapping_info.mapping_data_gen[quad_index].n_q_points_face[active_fe_index];
}



template <int dim, typename Number>
inline
const IndexSet &
MatrixFree<dim,Number>::get_locally_owned_set(const unsigned int dof_index) const
{
  AssertIndexRange (dof_index, dof_info.size());
  return dof_info[dof_index].vector_partitioner->locally_owned_range();
}



template <int dim, typename Number>
inline
const IndexSet &
MatrixFree<dim,Number>::get_ghost_set(const unsigned int dof_index) const
{
  AssertIndexRange (dof_index, dof_info.size());
  return dof_info[dof_index].vector_partitioner->ghost_indices();
}



template <int dim, typename Number>
inline
const internal::MatrixFreeFunctions::ShapeInfo<Number> &
MatrixFree<dim,Number>::get_shape_info (const unsigned int index_fe,
                                        const unsigned int index_quad,
                                        const unsigned int active_fe_index,
                                        const unsigned int active_quad_index) const
{
  AssertIndexRange (index_fe, shape_info.size(0));
  AssertIndexRange (index_quad, shape_info.size(1));
  AssertIndexRange (active_fe_index, shape_info.size(2));
  AssertIndexRange (active_quad_index, shape_info.size(3));
  return shape_info(index_fe, index_quad,
                    active_fe_index, active_quad_index);
}



template <int dim, typename Number>
inline
const Quadrature<dim> &
MatrixFree<dim,Number>::get_quadrature (const unsigned int quad_index,
                                        const unsigned int active_fe_index) const
{
  AssertIndexRange (quad_index, mapping_info.mapping_data_gen.size());
  return mapping_info.mapping_data_gen[quad_index].
         quadrature[active_fe_index];
}



template <int dim, typename Number>
inline
const Quadrature<dim-1> &
MatrixFree<dim,Number>::get_face_quadrature (const unsigned int quad_index,
                                             const unsigned int active_fe_index) const
{
  AssertIndexRange (quad_index, mapping_info.mapping_data_gen.size());
  return mapping_info.mapping_data_gen[quad_index].
         face_quadrature[active_fe_index];
}



template <int dim, typename Number>
inline
bool
MatrixFree<dim,Number>::indices_initialized () const
{
  return indices_are_initialized;
}



template <int dim, typename Number>
inline
bool
MatrixFree<dim,Number>::mapping_initialized () const
{
  return mapping_is_initialized;
}



// ------------------------------ reinit functions ---------------------------

namespace internal
{
  namespace MatrixFree
  {
    template <typename DH>
    inline
    std::vector<IndexSet>
    extract_locally_owned_index_sets (const std::vector<const DH *> &dofh,
                                      const unsigned int level)
    {
      std::vector<IndexSet> locally_owned_set;
      locally_owned_set.reserve (dofh.size());
      for (unsigned int j=0; j<dofh.size(); j++)
        if (level == numbers::invalid_unsigned_int)
          locally_owned_set.push_back(dofh[j]->locally_owned_dofs());
        else
          {
            // TODO: not distributed yet
            IndexSet new_set (dofh[j]->n_dofs(level));
            new_set.add_range (0, dofh[j]->n_dofs(level));
            locally_owned_set.push_back(new_set);
          }
      return locally_owned_set;
    }
  }
}



template <int dim, typename Number>
template <typename DH, typename Quad>
void MatrixFree<dim,Number>::
reinit(const DH               &dof_handler,
       const ConstraintMatrix &constraints_in,
       const Quad             &quad,
       const typename MatrixFree<dim,Number>::AdditionalData additional_data)
{
  MappingQ1<dim>                       mapping;
  std::vector<const DH *>               dof_handlers;
  std::vector<const ConstraintMatrix *> constraints;
  std::vector<Quad>          quads;

  dof_handlers.push_back(&dof_handler);
  constraints.push_back (&constraints_in);
  quads.push_back (quad);

  std::vector<IndexSet> locally_owned_sets =
    internal::MatrixFree::extract_locally_owned_index_sets
    (dof_handlers, additional_data.level_mg_handler);
  reinit(mapping, dof_handlers,constraints, locally_owned_sets, quads,
         additional_data);
}



template <int dim, typename Number>
template <typename DH, typename Quad>
void MatrixFree<dim,Number>::
reinit(const Mapping<dim>     &mapping,
       const DH               &dof_handler,
       const ConstraintMatrix &constraints_in,
       const Quad             &quad,
       const typename MatrixFree<dim,Number>::AdditionalData additional_data)
{
  std::vector<const DH *>               dof_handlers;
  std::vector<const ConstraintMatrix *> constraints;
  std::vector<Quad>          quads;

  dof_handlers.push_back(&dof_handler);
  constraints.push_back (&constraints_in);
  quads.push_back (quad);

  std::vector<IndexSet> locally_owned_sets =
    internal::MatrixFree::extract_locally_owned_index_sets
    (dof_handlers, additional_data.level_mg_handler);
  reinit(mapping, dof_handlers,constraints,locally_owned_sets, quads,
         additional_data);
}



template <int dim, typename Number>
template <typename DH, typename Quad>
void MatrixFree<dim,Number>::
reinit(const std::vector<const DH *>               &dof_handler,
       const std::vector<const ConstraintMatrix *> &constraint,
       const std::vector<Quad>                    &quad,
       const typename MatrixFree<dim,Number>::AdditionalData additional_data)
{
  MappingQ1<dim> mapping;
  std::vector<IndexSet> locally_owned_set =
    internal::MatrixFree::extract_locally_owned_index_sets
    (dof_handler, additional_data.level_mg_handler);
  reinit(mapping, dof_handler,constraint,locally_owned_set,
         static_cast<const std::vector<Quadrature<1> >&>(quad),
         additional_data);
}



template <int dim, typename Number>
template <typename DH, typename Quad>
void MatrixFree<dim,Number>::
reinit(const std::vector<const DH *>               &dof_handler,
       const std::vector<const ConstraintMatrix *> &constraint,
       const Quad                                 &quad,
       const typename MatrixFree<dim,Number>::AdditionalData additional_data)
{
  MappingQ1<dim> mapping;
  std::vector<Quad> quads;
  quads.push_back(quad);
  std::vector<IndexSet> locally_owned_set =
    internal::MatrixFree::extract_locally_owned_index_sets
    (dof_handler, additional_data.level_mg_handler);
  reinit(mapping, dof_handler,constraint,locally_owned_set, quads,
         additional_data);
}



template <int dim, typename Number>
template <typename DH, typename Quad>
void MatrixFree<dim,Number>::
reinit(const Mapping<dim>                         &mapping,
       const std::vector<const DH *>               &dof_handler,
       const std::vector<const ConstraintMatrix *> &constraint,
       const Quad                                 &quad,
       const typename MatrixFree<dim,Number>::AdditionalData additional_data)
{
  std::vector<Quad> quads;
  quads.push_back(quad);
  std::vector<IndexSet> locally_owned_set =
    internal::MatrixFree::extract_locally_owned_index_sets
    (dof_handler, additional_data.level_mg_handler);
  reinit(mapping, dof_handler,constraint,locally_owned_set, quads,
         additional_data);
}



template <int dim, typename Number>
template <typename DH, typename Quad>
void MatrixFree<dim,Number>::
reinit(const Mapping<dim>                         &mapping,
       const std::vector<const DH *>  &dof_handler,
       const std::vector<const ConstraintMatrix *> &constraint,
       const std::vector<Quad>              &quad,
       const typename MatrixFree<dim,Number>::AdditionalData additional_data)
{
  std::vector<IndexSet> locally_owned_set =
    internal::MatrixFree::extract_locally_owned_index_sets
    (dof_handler, additional_data.level_mg_handler);
  reinit(mapping, dof_handler,constraint,locally_owned_set,
         quad, additional_data);
}



namespace internal
{
  namespace MatrixFree
  {
    // resolve DoFHandler types

    // MGDoFHandler is deprecated in deal.II but might still be present in
    // user code, so we need to resolve its type (fortunately, it is derived
    // from DoFHandler, so we can static_cast it to a DoFHandler<dim>)
    template <typename DH>
    inline
    std::vector<const dealii::DoFHandler<DH::dimension> *>
    resolve_dof_handler (const std::vector<const DH *> &dof_handler)
    {
      std::vector<const dealii::DoFHandler<DH::dimension> *> conversion(dof_handler.size());
      for (unsigned int i=0; i<dof_handler.size(); ++i)
        conversion[i] = static_cast<const dealii::DoFHandler<DH::dimension> *>(dof_handler[i]);
      return conversion;
    }

    template <int dim>
    inline
    std::vector<const dealii::hp::DoFHandler<dim> *>
    resolve_dof_handler (const std::vector<const dealii::hp::DoFHandler<dim> *> &dof_handler)
    {
      return dof_handler;
    }
  }
}



template <int dim, typename Number>
template <typename DH, typename Quad>
void MatrixFree<dim,Number>::
reinit(const Mapping<dim>                         &mapping,
       const std::vector<const DH *>               &dof_handler,
       const std::vector<const ConstraintMatrix *> &constraint,
       const std::vector<IndexSet>                &locally_owned_set,
       const std::vector<Quad>                    &quad,
       const typename MatrixFree<dim,Number>::AdditionalData additional_data)
{
  // find out whether we use a hp Quadrature or a standard quadrature
  std::vector<hp::QCollection<1> > quad_hp;
  for (unsigned int q=0; q<quad.size(); ++q)
    quad_hp.push_back (hp::QCollection<1>(quad[q]));
  internal_reinit (mapping,
                   internal::MatrixFree::resolve_dof_handler(dof_handler),
                   constraint, locally_owned_set, quad_hp, additional_data);
}



// ------------------------------ implementation of cell_loop ---------------

// internal helper functions that define how to call MPI data exchange
// functions: for generic vectors, do nothing at all. For distributed vectors,
// call update_ghost_values_start function and so on. If we have collections
// of vectors, just do the individual functions of the components. In order to
// keep ghost values consistent (whether we are in read or write mode). the whole situation is a bit complicated by the fact
// that we need to treat block vectors differently, which use some additional
// helper functions to select the blocks and template magic.
namespace internal
{
  template<typename VectorStruct>
  bool update_ghost_values_start_block (const VectorStruct &vec,
                                        const unsigned int channel,
                                        internal::bool2type<true>);
  template<typename VectorStruct>
  void reset_ghost_values_block (const VectorStruct &vec,
                                 const bool          zero_out_ghosts,
                                 internal::bool2type<true>);
  template<typename VectorStruct>
  void update_ghost_values_finish_block (const VectorStruct &vec,
                                         internal::bool2type<true>);
  template<typename VectorStruct>
  void compress_start_block (const VectorStruct &vec,
                             const unsigned int channel,
                             internal::bool2type<true>);
  template<typename VectorStruct>
  void compress_finish_block (const VectorStruct &vec,
                              internal::bool2type<true>);

  template<typename VectorStruct>
  bool update_ghost_values_start_block (const VectorStruct &,
                                        const unsigned int,
                                        internal::bool2type<false>)
  {
    return false;
  }
  template<typename VectorStruct>
  void reset_ghost_values_block (const VectorStruct &,
                                 const bool,
                                 internal::bool2type<false>)
  {}
  template<typename VectorStruct>
  void update_ghost_values_finish_block (const VectorStruct &,
                                         internal::bool2type<false>)
  {}
  template<typename VectorStruct>
  void compress_start_block (const VectorStruct &,
                             const unsigned int,
                             internal::bool2type<false>)
  {}
  template<typename VectorStruct>
  void compress_finish_block (const VectorStruct &,
                              internal::bool2type<false>)
  {}



  // returns true if the vector was in a state without ghost values before,
  // i.e., we need to zero out ghosts in the very end
  template<typename VectorStruct>
  inline
  bool update_ghost_values_start (const VectorStruct &vec,
                                  const unsigned int channel = 0)
  {
    return
      update_ghost_values_start_block(vec, channel,
                                      internal::bool2type<IsBlockVector<VectorStruct>::value>());
  }



  template<typename Number>
  inline
  bool update_ghost_values_start (const parallel::distributed::Vector<Number> &vec,
                                  const unsigned int                  channel = 0)
  {
    bool return_value = !vec.has_ghost_elements();
    vec.update_ghost_values_start(channel);
    return return_value;
  }



  template <typename VectorStruct>
  inline
  bool update_ghost_values_start (const std::vector<VectorStruct> &vec)
  {
    bool return_value = false;
    for (unsigned int comp=0; comp<vec.size(); comp++)
      return_value = update_ghost_values_start(vec[comp], comp);
    return return_value;
  }



  template <typename VectorStruct>
  inline
  bool update_ghost_values_start (const std::vector<VectorStruct *> &vec)
  {
    bool return_value = false;
    for (unsigned int comp=0; comp<vec.size(); comp++)
      return_value = update_ghost_values_start(*vec[comp], comp);
    return return_value;
  }



  template<typename VectorStruct>
  inline
  bool update_ghost_values_start_block (const VectorStruct &vec,
                                        const unsigned int channel,
                                        internal::bool2type<true>)
  {
    bool return_value = false;
    for (unsigned int i=0; i<vec.n_blocks(); ++i)
      return_value = update_ghost_values_start(vec.block(i), channel+509*i);
    return return_value;
  }



  // if the input vector did not have ghosts imported, clear them here again
  // in order to avoid subsequent operations e.g. in linear solvers to work
  // with ghosts all the time
  template<typename VectorStruct>
  inline
  void reset_ghost_values (const VectorStruct &vec,
                           const bool          zero_out_ghosts)
  {
    reset_ghost_values_block(vec, zero_out_ghosts,
                             internal::bool2type<IsBlockVector<VectorStruct>::value>());
  }



  template<typename Number>
  inline
  void reset_ghost_values (const parallel::distributed::Vector<Number> &vec,
                           const bool zero_out_ghosts)
  {
    if (zero_out_ghosts)
      const_cast<parallel::distributed::Vector<Number>&>(vec).zero_out_ghosts();
  }



  template <typename VectorStruct>
  inline
  void reset_ghost_values (const std::vector<VectorStruct> &vec,
                           const bool zero_out_ghosts)
  {
    for (unsigned int comp=0; comp<vec.size(); comp++)
      reset_ghost_values(vec[comp], zero_out_ghosts);
  }



  template <typename VectorStruct>
  inline
  void reset_ghost_values (const std::vector<VectorStruct *> &vec,
                           const bool zero_out_ghosts)
  {
    for (unsigned int comp=0; comp<vec.size(); comp++)
      reset_ghost_values(*vec[comp], zero_out_ghosts);
  }



  template<typename VectorStruct>
  inline
  void reset_ghost_values_block (const VectorStruct &vec,
                                 const bool          zero_out_ghosts,
                                 internal::bool2type<true>)
  {
    for (unsigned int i=0; i<vec.n_blocks(); ++i)
      reset_ghost_values(vec.block(i), zero_out_ghosts);
  }



  template <typename VectorStruct>
  inline
  void update_ghost_values_finish (const VectorStruct &vec)
  {
    update_ghost_values_finish_block(vec,
                                     internal::bool2type<IsBlockVector<VectorStruct>::value>());
  }



  template <typename Number>
  inline
  void update_ghost_values_finish (const parallel::distributed::Vector<Number> &vec)
  {
    vec.update_ghost_values_finish();
  }



  template <typename VectorStruct>
  inline
  void update_ghost_values_finish (const std::vector<VectorStruct> &vec)
  {
    for (unsigned int comp=0; comp<vec.size(); comp++)
      update_ghost_values_finish(vec[comp]);
  }



  template <typename VectorStruct>
  inline
  void update_ghost_values_finish (const std::vector<VectorStruct *> &vec)
  {
    for (unsigned int comp=0; comp<vec.size(); comp++)
      update_ghost_values_finish(*vec[comp]);
  }



  template <typename VectorStruct>
  inline
  void update_ghost_values_finish_block (const VectorStruct &vec,
                                         internal::bool2type<true>)
  {
    for (unsigned int i=0; i<vec.n_blocks(); ++i)
      update_ghost_values_finish(vec.block(i));
  }



  template <typename VectorStruct>
  inline
  void compress_start (VectorStruct &vec,
                       const unsigned int channel = 0)
  {
    compress_start_block (vec, channel,
                          internal::bool2type<IsBlockVector<VectorStruct>::value>());
  }



  template <typename Number>
  inline
  void compress_start (parallel::distributed::Vector<Number> &vec,
                       const unsigned int           channel = 0)
  {
    vec.compress_start(channel);
  }



  template <typename VectorStruct>
  inline
  void compress_start (std::vector<VectorStruct> &vec)
  {
    for (unsigned int comp=0; comp<vec.size(); comp++)
      compress_start (vec[comp], comp);
  }



  template <typename VectorStruct>
  inline
  void compress_start (std::vector<VectorStruct *> &vec)
  {
    for (unsigned int comp=0; comp<vec.size(); comp++)
      compress_start (*vec[comp], comp);
  }



  template <typename VectorStruct>
  inline
  void compress_start_block (VectorStruct      &vec,
                             const unsigned int channel,
                             internal::bool2type<true>)
  {
    for (unsigned int i=0; i<vec.n_blocks(); ++i)
      compress_start(vec.block(i), channel + 500*i);
  }



  template <typename VectorStruct>
  inline
  void compress_finish (VectorStruct &vec)
  {
    compress_finish_block(vec,
                          internal::bool2type<IsBlockVector<VectorStruct>::value>());
  }



  template <typename Number>
  inline
  void compress_finish (parallel::distributed::Vector<Number> &vec)
  {
    vec.compress_finish(::dealii::VectorOperation::add);
  }



  template <typename VectorStruct>
  inline
  void compress_finish (std::vector<VectorStruct> &vec)
  {
    for (unsigned int comp=0; comp<vec.size(); comp++)
      compress_finish(vec[comp]);
  }



  template <typename VectorStruct>
  inline
  void compress_finish (std::vector<VectorStruct *> &vec)
  {
    for (unsigned int comp=0; comp<vec.size(); comp++)
      compress_finish(*vec[comp]);
  }



  template <typename VectorStruct>
  inline
  void compress_finish_block (VectorStruct &vec,
                              internal::bool2type<true>)
  {
    for (unsigned int i=0; i<vec.n_blocks(); ++i)
      compress_finish(vec.block(i));
  }



#ifdef DEAL_II_WITH_THREADS

  // This defines the TBB data structures that are needed to schedule the
  // partition-partition variant

  namespace partition
  {
    template<typename Worker>
    class CellWork : public tbb::task
    {
    public:
      CellWork (const Worker &worker_in,
                const unsigned int partition_in,
                const internal::MatrixFreeFunctions::TaskInfo &task_info_in,
                const bool is_blocked_in)
        :
        worker (worker_in),
        partition (partition_in),
        task_info (task_info_in),
        is_blocked (is_blocked_in)
      {};
      tbb::task *execute ()
      {
        std::pair<unsigned int, unsigned int> cell_range
        (task_info.partition_color_blocks_data[partition],
         task_info.partition_color_blocks_data[partition+1]);
        worker(cell_range);
        if (is_blocked==true)
          dummy->spawn (*dummy);
        return NULL;
      }

      tbb::empty_task *dummy;

    private:
      const Worker      &worker;
      const unsigned int partition;
      const internal::MatrixFreeFunctions::TaskInfo &task_info;
      const bool         is_blocked;
    };



    template<typename Worker>
    class PartitionWork : public tbb::task
    {
    public:
      PartitionWork (const Worker &function_in,
                     const unsigned int partition_in,
                     const internal::MatrixFreeFunctions::TaskInfo &task_info_in,
                     const bool    is_blocked_in = false)
        :
        function (function_in),
        partition (partition_in),
        task_info (task_info_in),
        is_blocked (is_blocked_in)
      {};
      tbb::task *execute ()
      {
        tbb::empty_task *root = new( tbb::task::allocate_root() )
        tbb::empty_task;
        unsigned int evens = task_info.partition_evens[partition];
        unsigned int odds  = task_info.partition_odds[partition];
        unsigned int n_blocked_workers =
          task_info.partition_n_blocked_workers[partition];
        unsigned int n_workers = task_info.partition_n_workers[partition];
        std::vector<CellWork<Worker>*> worker(n_workers);
        std::vector<CellWork<Worker>*> blocked_worker(n_blocked_workers);

        root->set_ref_count(evens+1);
        for (unsigned int j=0; j<evens; j++)
          {
            worker[j] = new(root->allocate_child())
            CellWork<Worker>(function, task_info.
                             partition_color_blocks_row_index[partition]+2*j,
                             task_info, false);
            if (j>0)
              {
                worker[j]->set_ref_count(2);
                blocked_worker[j-1]->dummy = new(worker[j]->allocate_child())
                tbb::empty_task;
                worker[j-1]->spawn(*blocked_worker[j-1]);
              }
            else
              worker[j]->set_ref_count(1);
            if (j<evens-1)
              {
                blocked_worker[j] = new(worker[j]->allocate_child())
                CellWork<Worker>(function, task_info.
                                 partition_color_blocks_row_index
                                 [partition] + 2*j+1, task_info, true);
              }
            else
              {
                if (odds==evens)
                  {
                    worker[evens] = new(worker[j]->allocate_child())
                    CellWork<Worker>(function, task_info.
                                     partition_color_blocks_row_index[partition]+2*j+1,
                                     task_info, false);
                    worker[j]->spawn(*worker[evens]);
                  }
                else
                  {
                    tbb::empty_task *child = new(worker[j]->allocate_child())
                    tbb::empty_task();
                    worker[j]->spawn(*child);
                  }
              }
          }

        root->wait_for_all();
        root->destroy(*root);
        if (is_blocked==true)
          dummy->spawn (*dummy);
        return NULL;
      }

      tbb::empty_task *dummy;

    private:
      const Worker  &function;
      const unsigned int partition;
      const internal::MatrixFreeFunctions::TaskInfo &task_info;
      const bool     is_blocked;
    };

  } // end of namespace partition



  namespace color
  {
    template <typename Worker>
    class CellWork
    {
    public:
      CellWork (const Worker                   &worker_in,
                const internal::MatrixFreeFunctions::TaskInfo &task_info_in)
        :
        worker (worker_in),
        task_info (task_info_in)
      {};
      void operator()(const tbb::blocked_range<unsigned int> &r) const
      {
        for (unsigned int block=r.begin(); block<r.end(); block++)
          {
            std::pair<unsigned int,unsigned int> cell_range;
            if (task_info.position_short_block<block)
              {
                cell_range.first = (block-1)*task_info.block_size+
                                   task_info.block_size_last;
                cell_range.second = cell_range.first + task_info.block_size;
              }
            else
              {
                cell_range.first = block*task_info.block_size;
                cell_range.second = cell_range.first +
                                    ((block == task_info.position_short_block)?
                                     (task_info.block_size_last):(task_info.block_size));
              }
            worker (cell_range);
          }
      }
    private:
      const Worker   &worker;
      const internal::MatrixFreeFunctions::TaskInfo &task_info;
    };


    template<typename Worker>
    class PartitionWork : public tbb::task
    {
    public:
      PartitionWork (const Worker &worker_in,
                     const unsigned int partition_in,
                     const internal::MatrixFreeFunctions::TaskInfo &task_info_in,
                     const bool    is_blocked_in)
        :
        worker (worker_in),
        partition (partition_in),
        task_info (task_info_in),
        is_blocked (is_blocked_in)
      {};
      tbb::task *execute ()
      {
        unsigned int lower = task_info.partition_color_blocks_data[partition],
                     upper = task_info.partition_color_blocks_data[partition+1];
        parallel_for(tbb::blocked_range<unsigned int>(lower,upper,1),
                     CellWork<Worker> (worker,task_info));
        if (is_blocked==true)
          dummy->spawn (*dummy);
        return NULL;
      }

      tbb::empty_task *dummy;

    private:
      const Worker &worker;
      const unsigned int partition;
      const internal::MatrixFreeFunctions::TaskInfo &task_info;
      const bool is_blocked;
    };

  } // end of namespace color


  template<typename VectorStruct>
  class MPIComDistribute : public tbb::task
  {
  public:
    MPIComDistribute (const VectorStruct  &src_in)
      :
      src(src_in)
    {};

    tbb::task *execute ()
    {
      internal::update_ghost_values_finish(src);
      return 0;
    }

  private:
    const VectorStruct &src;
  };



  template<typename VectorStruct>
  class MPIComCompress : public tbb::task
  {
  public:
    MPIComCompress (VectorStruct        &dst_in)
      :
      dst(dst_in)
    {};

    tbb::task *execute ()
    {
      internal::compress_start(dst);
      return 0;
    }

  private:
    VectorStruct &dst;
  };

#endif // DEAL_II_WITH_THREADS

} // end of namespace internal



template <int dim, typename Number>
template <typename OutVector, typename InVector>
inline
void
MatrixFree<dim, Number>::cell_loop
(const std_cxx11::function<void (const MatrixFree<dim,Number> &,
                                 OutVector &,
                                 const InVector &,
                                 const std::pair<unsigned int,
                                 unsigned int> &)> &cell_operation,
 OutVector       &dst,
 const InVector  &src) const
{
  // in any case, need to start the ghost import at the beginning
  bool ghosts_were_not_set = internal::update_ghost_values_start (src);

#ifdef DEAL_II_WITH_THREADS

  // Use multithreading if so requested and if there is enough work to do in
  // parallel (the code might hang if there are less than two chunks!)
  if (task_info.use_multithreading == true && task_info.n_blocks > 3)
    {
      // to simplify the function calls, bind away all arguments except the
      // cell range
      typedef
      std_cxx11::function<void (const std::pair<unsigned int,unsigned int> &range)>
      Worker;

      const Worker func = std_cxx11::bind (std_cxx11::ref(cell_operation),
                                           std_cxx11::cref(*this),
                                           std_cxx11::ref(dst),
                                           std_cxx11::cref(src),
                                           std_cxx11::_1);

      if (task_info.use_partition_partition == true)
        {
          tbb::empty_task *root = new( tbb::task::allocate_root() )
          tbb::empty_task;
          unsigned int evens = task_info.evens;
          unsigned int odds  = task_info.odds;
          root->set_ref_count(evens+1);
          unsigned int n_blocked_workers = task_info.n_blocked_workers;
          unsigned int n_workers = task_info.n_workers;
          std::vector<internal::partition::PartitionWork<Worker>*>
          worker(n_workers);
          std::vector<internal::partition::PartitionWork<Worker>*>
          blocked_worker(n_blocked_workers);
          internal::MPIComCompress<OutVector> *worker_compr =
            new(root->allocate_child())
          internal::MPIComCompress<OutVector>(dst);
          worker_compr->set_ref_count(1);
          for (unsigned int j=0; j<evens; j++)
            {
              if (j>0)
                {
                  worker[j] = new(root->allocate_child())
                  internal::partition::PartitionWork<Worker>
                  (func,2*j,task_info,false);
                  worker[j]->set_ref_count(2);
                  blocked_worker[j-1]->dummy = new(worker[j]->allocate_child())
                  tbb::empty_task;
                  if (j>1)
                    worker[j-1]->spawn(*blocked_worker[j-1]);
                  else
                    worker_compr->spawn(*blocked_worker[j-1]);
                }
              else
                {
                  worker[j] = new(worker_compr->allocate_child())
                  internal::partition::PartitionWork<Worker>
                  (func,2*j,task_info,false);
                  worker[j]->set_ref_count(2);
                  internal::MPIComDistribute<InVector> *worker_dist =
                    new (worker[j]->allocate_child())
                  internal::MPIComDistribute<InVector>(src);
                  worker_dist->spawn(*worker_dist);
                }
              if (j<evens-1)
                {
                  blocked_worker[j] = new(worker[j]->allocate_child())
                  internal::partition::PartitionWork<Worker>
                  (func,2*j+1,task_info,true);
                }
              else
                {
                  if (odds==evens)
                    {
                      worker[evens] = new(worker[j]->allocate_child())
                      internal::partition::PartitionWork<Worker>
                      (func,2*j+1,task_info,false);
                      worker[j]->spawn(*worker[evens]);
                    }
                  else
                    {
                      tbb::empty_task *child = new(worker[j]->allocate_child())
                      tbb::empty_task();
                      worker[j]->spawn(*child);
                    }
                }
            }

          root->wait_for_all();
          root->destroy(*root);
        }
      else // end of partition-partition, start of partition-color
        {
          unsigned int evens = task_info.evens;
          unsigned int odds  = task_info.odds;

          // check whether there is only one partition. if not, build up the
          // tree of partitions
          if (odds > 0)
            {
              tbb::empty_task *root = new( tbb::task::allocate_root() ) tbb::empty_task;
              root->set_ref_count(evens+1);
              unsigned int n_blocked_workers = odds-(odds+evens+1)%2;
              unsigned int n_workers = task_info.partition_color_blocks_data.size()-1-
                                       n_blocked_workers;
              std::vector<internal::color::PartitionWork<Worker>*> worker(n_workers);
              std::vector<internal::color::PartitionWork<Worker>*> blocked_worker(n_blocked_workers);
              unsigned int worker_index = 0, slice_index = 0;
              unsigned int spawn_index =  0, spawn_index_new = 0;
              int spawn_index_child = -2;
              internal::MPIComCompress<OutVector> *worker_compr = new(root->allocate_child())
              internal::MPIComCompress<OutVector>(dst);
              worker_compr->set_ref_count(1);
              for (unsigned int part=0;
                   part<task_info.partition_color_blocks_row_index.size()-1; part++)
                {
                  spawn_index_new = worker_index;
                  if (part == 0)
                    worker[worker_index] = new(worker_compr->allocate_child())
                    internal::color::PartitionWork<Worker>(func,slice_index,task_info,false);
                  else
                    worker[worker_index] = new(root->allocate_child())
                    internal::color::PartitionWork<Worker>(func,slice_index,task_info,false);
                  slice_index++;
                  for (; slice_index<task_info.partition_color_blocks_row_index[part+1];
                       slice_index++)
                    {
                      worker[worker_index]->set_ref_count(1);
                      worker_index++;
                      worker[worker_index] = new (worker[worker_index-1]->allocate_child())
                      internal::color::PartitionWork<Worker>(func,slice_index,task_info,false);
                    }
                  worker[worker_index]->set_ref_count(2);
                  if (part>0)
                    {
                      blocked_worker[(part-1)/2]->dummy =
                        new (worker[worker_index]->allocate_child()) tbb::empty_task;
                      worker_index++;
                      if (spawn_index_child == -1)
                        worker[spawn_index]->spawn(*blocked_worker[(part-1)/2]);
                      else
                        worker[spawn_index]->spawn(*worker[spawn_index_child]);
                      spawn_index = spawn_index_new;
                      spawn_index_child = -2;
                    }
                  else
                    {
                      internal::MPIComDistribute<InVector> *worker_dist =
                        new (worker[worker_index]->allocate_child())
                      internal::MPIComDistribute<InVector>(src);
                      worker_dist->spawn(*worker_dist);
                      worker_index++;
                    }
                  part += 1;
                  if (part<task_info.partition_color_blocks_row_index.size()-1)
                    {
                      if (part<task_info.partition_color_blocks_row_index.size()-2)
                        {
                          blocked_worker[part/2] = new(worker[worker_index-1]->allocate_child())
                          internal::color::PartitionWork<Worker>(func,slice_index,task_info,true);
                          slice_index++;
                          if (slice_index<
                              task_info.partition_color_blocks_row_index[part+1])
                            {
                              blocked_worker[part/2]->set_ref_count(1);
                              worker[worker_index] = new(blocked_worker[part/2]->allocate_child())
                              internal::color::PartitionWork<Worker>(func,slice_index,task_info,false);
                              slice_index++;
                            }
                          else
                            {
                              spawn_index_child = -1;
                              continue;
                            }
                        }
                      for (; slice_index<task_info.partition_color_blocks_row_index[part+1];
                           slice_index++)
                        {
                          if (slice_index>
                              task_info.partition_color_blocks_row_index[part])
                            {
                              worker[worker_index]->set_ref_count(1);
                              worker_index++;
                            }
                          worker[worker_index] = new (worker[worker_index-1]->allocate_child())
                          internal::color::PartitionWork<Worker>(func,slice_index,task_info,false);
                        }
                      spawn_index_child = worker_index;
                      worker_index++;
                    }
                  else
                    {
                      tbb::empty_task *final = new (worker[worker_index-1]->allocate_child())
                      tbb::empty_task;
                      worker[spawn_index]->spawn(*final);
                      spawn_index_child = worker_index-1;
                    }
                }
              if (evens==odds)
                worker[spawn_index]->spawn(*worker[spawn_index_child]);
              root->wait_for_all();
              root->destroy(*root);
            }
          // case when we only have one partition: this is the usual coloring
          // scheme, and we just schedule a parallel for loop for each color
          else
            {
              Assert(evens==1,ExcInternalError());
              internal::update_ghost_values_finish(src);

              for (unsigned int color=0;
                   color < task_info.partition_color_blocks_row_index[1];
                   ++color)
                {
                  unsigned int lower = task_info.partition_color_blocks_data[color],
                               upper = task_info.partition_color_blocks_data[color+1];
                  parallel_for(tbb::blocked_range<unsigned int>(lower,upper,1),
                               internal::color::CellWork<Worker>
                               (func,task_info));
                }

              internal::compress_start(dst);
            }
        }
    }
  else
#endif
    // serial loop
    {
      std::pair<unsigned int,unsigned int> cell_range;

      // First operate on cells where no ghost data is needed (inner cells)
      {
        cell_range.first = 0;
        cell_range.second = size_info.boundary_cells_start;
        cell_operation (*this, dst, src, cell_range);
      }

      // before starting operations on cells that contain ghost nodes (outer
      // cells), wait for the MPI commands to finish
      internal::update_ghost_values_finish(src);

      // For the outer cells, do the same procedure as for inner cells.
      if (size_info.boundary_cells_end > size_info.boundary_cells_start)
        {
          cell_range.first = size_info.boundary_cells_start;
          cell_range.second = size_info.boundary_cells_end;
          cell_operation (*this, dst, src, cell_range);
        }

      internal::compress_start(dst);

      // Finally operate on cells where no ghost data is needed (inner cells)
      if (size_info.n_macro_cells > size_info.boundary_cells_end)
        {
          cell_range.first = size_info.boundary_cells_end;
          cell_range.second = size_info.n_macro_cells;
          cell_operation (*this, dst, src, cell_range);
        }
    }

  // In every case, we need to finish transfers at the very end
  internal::compress_finish(dst);
  internal::reset_ghost_values(src, ghosts_were_not_set);
}



template <int dim, typename Number>
template <typename CLASS, typename OutVector, typename InVector>
inline
void
MatrixFree<dim,Number>::cell_loop
(void (CLASS::*function_pointer)(const MatrixFree<dim,Number> &,
                                 OutVector &,
                                 const InVector &,
                                 const std::pair<unsigned int,
                                 unsigned int> &)const,
 const CLASS    *owning_class,
 OutVector      &dst,
 const InVector &src) const
{
  // here, use std_cxx11::bind to hand a function handler with the appropriate
  // argument to the other loop function
  std_cxx11::function<void (const MatrixFree<dim,Number> &,
                            OutVector &,
                            const InVector &,
                            const std::pair<unsigned int,
                            unsigned int> &)>
  function = std_cxx11::bind<void>(function_pointer,
                                   std_cxx11::cref(*owning_class),
                                   std_cxx11::_1,
                                   std_cxx11::_2,
                                   std_cxx11::_3,
                                   std_cxx11::_4);
  cell_loop (function, dst, src);
}



template <int dim, typename Number>
template <typename CLASS, typename OutVector, typename InVector>
inline
void
MatrixFree<dim,Number>::cell_loop
(void(CLASS::*function_pointer)(const MatrixFree<dim,Number> &,
                                OutVector &,
                                const InVector &,
                                const std::pair<unsigned int,
                                unsigned int> &),
 CLASS          *owning_class,
 OutVector      &dst,
 const InVector &src) const
{
  // here, use std_cxx11::bind to hand a function handler with the appropriate
  // argument to the other loop function
  std_cxx11::function<void (const MatrixFree<dim,Number> &,
                            OutVector &,
                            const InVector &,
                            const std::pair<unsigned int,
                            unsigned int> &)>
  function = std_cxx11::bind<void>(function_pointer,
                                   std_cxx11::ref(*owning_class),
                                   std_cxx11::_1,
                                   std_cxx11::_2,
                                   std_cxx11::_3,
                                   std_cxx11::_4);
  cell_loop (function, dst, src);
}


#endif  // ifndef DOXYGEN



DEAL_II_NAMESPACE_CLOSE

#endif
