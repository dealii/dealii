// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_matrix_free_h
#define dealii_matrix_free_h

#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/thread_local_storage.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_vector_base.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector_operation.h>

#include <deal.II/matrix_free/dof_info.h>
#include <deal.II/matrix_free/mapping_info.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/task_info.h>
#include <deal.II/matrix_free/type_traits.h>
#include <deal.II/matrix_free/vector_data_exchange.h>

#include <cstdlib>
#include <limits>
#include <list>
#include <memory>


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
 * of freedom. It includes a description of constraints that are evaluated as
 * going through all local degrees of freedom on a cell.
 *
 * - MappingInfo: It stores the transformations from real to unit cells that
 * are necessary in order to build derivatives of finite element functions and
 * find location of quadrature weights in physical space.
 *
 * - ShapeInfo: It contains the shape functions of the finite element,
 * evaluated on the unit cell.
 *
 * Besides the initialization routines, this class implements only a single
 * operation, namely a loop over all cells (cell_loop()). This loop is
 * scheduled in such a way that cells that share degrees of freedom are not
 * worked on simultaneously, which implies that it is possible to write to
 * vectors (or matrices) in parallel without having to explicitly synchronize
 * access to these vectors and matrices. This class does not implement any
 * shape values, all it does is to cache the respective data. To implement
 * finite element operations, use the class FEEvaluation (or some of the
 * related classes).
 *
 * This class traverses the cells in a different order than the usual
 * Triangulation class in deal.II, in order to be flexible with respect to
 * parallelization in shared memory and vectorization.
 *
 * Vectorization is implemented by merging several topological cells into one
 * so-called `cell batch`. This enables the application of all cell-related
 * operations for several cells with one CPU instruction and is one of the
 * main features of this framework.
 *
 * For details on usage of this class, see the description of FEEvaluation or
 * the
 * @ref matrixfree "matrix-free topic".
 *
 * @ingroup matrixfree
 */

template <int dim,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
class MatrixFree : public EnableObserverPointer
{
  static_assert(
    std::is_same_v<Number, typename VectorizedArrayType::value_type>,
    "Type of Number and of VectorizedArrayType do not match.");

public:
  /**
   * An alias for the underlying number type specified by the template
   * argument.
   */
  using value_type            = Number;
  using vectorized_value_type = VectorizedArrayType;

  /**
   * The dimension set by the template argument `dim`.
   */
  static constexpr unsigned int dimension = dim;

  /**
   * Collects the options for initialization of the MatrixFree class. The
   * parameter @p tasks_parallel_scheme specifies the
   * parallelization options in shared memory (task-based parallelism, where
   * one can choose between no parallelism and three schemes that avoid that
   * cells with access to the same vector entries are accessed
   * simultaneously), and the parameter @p tasks_block_size the block size for
   * task parallel scheduling. The parameters @p mapping_update_flags,
   * @p mapping_update_flags_boundary_faces, @p mapping_update_flags_inner_faces,
   * and @p mapping_update_flags_faces_by_cells specify the update flags that
   * should be stored by this class.
   *
   * The parameter @p mg_level specifies the level in the triangulation from which
   * the indices are to be used. If the level is set to
   * `numbers::invalid_unsigned_int`, the active cells are traversed, and
   * otherwise the cells in the given level. This option has no effect in case
   * a DoFHandler is given.
   *
   * The parameter @p store_plain_indices indicates whether the DoFInfo
   * class should also allow for access to vectors without resolving
   * constraints.
   *
   * The two parameters @p initialize_indices and @p initialize_mapping allow
   * the user to disable some of the initialization processes. For example, if
   * only the scheduling that avoids touching the same vector/matrix indices
   * simultaneously is to be found, the mapping needs not be
   * initialized. Likewise, if the mapping has changed from one iteration to
   * the next but the topology has not (like when using a deforming mesh with
   * MappingQEulerian), it suffices to initialize the mapping only.
   *
   * The two parameters @p cell_vectorization_categories and
   * @p cell_vectorization_categories_strict control the formation of batches
   * for vectorization over several cells. It is used implicitly when working
   * with hp-adaptivity but can also be useful in other contexts, such as in
   * local time stepping where one would like to control which elements
   * together form a batch of cells. The array @p cell_vectorization_categories
   * is accessed by the number given by cell->active_cell_index() when working
   * on the active cells with @p mg_level set to `numbers::invalid_unsigned_int`
   * and by cell->index() for the level cells. By default, the different
   * categories in @p cell_vectorization_category can be mixed and the algorithm
   * is allowed to merge lower category numbers with the next higher categories
   * if it is necessary inside the algorithm, in order to avoid partially
   * filled SIMD lanes as much as possible. This gives a better utilization of
   * the vectorization but might need special treatment, in particular for
   * face integrals. If set to `true', the algorithm will instead keep
   * different categories separate and not mix them in a single vectorized
   * array.
   *
   * Finally, @p allow_ghosted_vectors_in_loops allows to enable and disable
   * checks and @p communicator_sm gives the MPI communicator to be used
   * if MPI-3.0 shared-memory features should be used.
   */
  struct AdditionalData
  {
    /**
     * Provide the type of the surrounding MatrixFree class.
     */
    using MatrixFreeType = MatrixFree<dim, Number, VectorizedArrayType>;

    /**
     * Collects options for task parallelism. See the documentation of the
     * member variable MatrixFree::AdditionalData::tasks_parallel_scheme for a
     * thorough description.
     */
    enum TasksParallelScheme
    {
      /**
       * Perform application in serial.
       */
      none = internal::MatrixFreeFunctions::TaskInfo::none,
      /**
       * Partition the cells into two levels and afterwards form chunks.
       */
      partition_partition =
        internal::MatrixFreeFunctions::TaskInfo::partition_partition,
      /**
       * Partition on the global level and color cells within the partitions.
       */
      partition_color =
        internal::MatrixFreeFunctions::TaskInfo::partition_color,
      /**
       * Use the traditional coloring algorithm: this is like
       * TasksParallelScheme::partition_color, but only uses one partition.
       */
      color = internal::MatrixFreeFunctions::TaskInfo::color
    };

    /**
     * Constructor for AdditionalData.
     */
    AdditionalData(
      const TasksParallelScheme tasks_parallel_scheme = partition_partition,
      const unsigned int        tasks_block_size      = 0,
      const UpdateFlags         mapping_update_flags  = update_gradients |
                                               update_JxW_values,
      const UpdateFlags  mapping_update_flags_boundary_faces = update_default,
      const UpdateFlags  mapping_update_flags_inner_faces    = update_default,
      const UpdateFlags  mapping_update_flags_faces_by_cells = update_default,
      const unsigned int mg_level            = numbers::invalid_unsigned_int,
      const bool         store_plain_indices = true,
      const bool         initialize_indices  = true,
      const bool         initialize_mapping  = true,
      const bool         overlap_communication_computation    = true,
      const bool         hold_all_faces_to_owned_cells        = false,
      const bool         cell_vectorization_categories_strict = false,
      const bool         allow_ghosted_vectors_in_loops       = true)
      : tasks_parallel_scheme(tasks_parallel_scheme)
      , tasks_block_size(tasks_block_size)
      , mapping_update_flags(mapping_update_flags)
      , mapping_update_flags_boundary_faces(mapping_update_flags_boundary_faces)
      , mapping_update_flags_inner_faces(mapping_update_flags_inner_faces)
      , mapping_update_flags_faces_by_cells(mapping_update_flags_faces_by_cells)
      , mg_level(mg_level)
      , store_plain_indices(store_plain_indices)
      , initialize_indices(initialize_indices)
      , initialize_mapping(initialize_mapping)
      , overlap_communication_computation(overlap_communication_computation)
      , hold_all_faces_to_owned_cells(hold_all_faces_to_owned_cells)
      , cell_vectorization_categories_strict(
          cell_vectorization_categories_strict)
      , allow_ghosted_vectors_in_loops(allow_ghosted_vectors_in_loops)
      , store_ghost_cells(false)
      , communicator_sm(MPI_COMM_SELF)
    {}

    /**
     * Copy constructor.
     */
    AdditionalData(const AdditionalData &other)
      : tasks_parallel_scheme(other.tasks_parallel_scheme)
      , tasks_block_size(other.tasks_block_size)
      , mapping_update_flags(other.mapping_update_flags)
      , mapping_update_flags_boundary_faces(
          other.mapping_update_flags_boundary_faces)
      , mapping_update_flags_inner_faces(other.mapping_update_flags_inner_faces)
      , mapping_update_flags_faces_by_cells(
          other.mapping_update_flags_faces_by_cells)
      , mg_level(other.mg_level)
      , store_plain_indices(other.store_plain_indices)
      , initialize_indices(other.initialize_indices)
      , initialize_mapping(other.initialize_mapping)
      , overlap_communication_computation(
          other.overlap_communication_computation)
      , hold_all_faces_to_owned_cells(other.hold_all_faces_to_owned_cells)
      , cell_vectorization_category(other.cell_vectorization_category)
      , cell_vectorization_categories_strict(
          other.cell_vectorization_categories_strict)
      , allow_ghosted_vectors_in_loops(other.allow_ghosted_vectors_in_loops)
      , store_ghost_cells(other.store_ghost_cells)
      , communicator_sm(other.communicator_sm)
    {}

    /**
     * Copy assignment.
     */
    AdditionalData &
    operator=(const AdditionalData &other) = default;

    /**
     * Set the scheme for task parallelism. There are four options available.
     * If set to @p none, the operator application is done in serial without
     * shared memory parallelism. If this class is used together with MPI and
     * MPI is also used for parallelism within the nodes, this flag should be
     * set to @p none. The default value is @p partition_partition, i.e. we
     * actually use multithreading with the first option described below.
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
     * hanging nodes, there are quite many colors (50 or more in 3d), which
     * might degrade parallel performance (bad cache behavior, many
     * synchronization points).
     *
     * @note Threading support is currently experimental for the case inner
     * face integrals are performed and it is recommended to use MPI
     * parallelism if possible. While the scheme has been verified to work
     * with the `partition_partition` option in case of usual DG elements, no
     * comprehensive tests have been performed for systems of more general
     * elements, like combinations of continuous and discontinuous elements
     * that add face integrals to all terms.
     */
    TasksParallelScheme tasks_parallel_scheme;

    /**
     * Set the number of so-called cell batches that should form one
     * partition. If zero size is given, the class tries to find a good size
     * for the blocks based on MultithreadInfo::n_threads() and the number of
     * cells present. Otherwise, the given number is used. If the given number
     * is larger than one third of the number of total cells, this means no
     * parallelism. Note that in the case vectorization is used, a cell batch
     * consists of more than one physical cell.
     */
    unsigned int tasks_block_size;

    /**
     * This flag determines what data needs to be computed and cached on cells.
     *
     * If your computations require operations like quadrature point locations
     * or Hessians, these need to specified here (update_quadrature_points
     * or update_hessians, respectively). Note that values, gradients, and
     * Jacobian determinants (JxW values) are always computed regardless of the
     * flags specified here.
     *
     * Note that some additional flags might be set automatically (for example
     * second derivatives might be evaluated on Cartesian cells since there
     * the Jacobian describes the mapping completely).
     */
    UpdateFlags mapping_update_flags;

    /**
     * This flag determines the mapping data on boundary faces to be
     * cached. Note that MatrixFree uses a separate loop layout for face
     * integrals in order to effectively vectorize also in the case of hanging
     * nodes (which require different subface settings on the two sides) or
     * some cells in the batch of a VectorizedArray of cells that are adjacent
     * to the boundary and others that are not.
     *
     * If set to a value different from update_general (default), the face
     * information is explicitly built. Currently, MatrixFree supports to
     * cache the following data on faces: inverse Jacobians, Jacobian
     * determinants (JxW), quadrature points, data for Hessians (derivative of
     * Jacobians), and normal vectors.
     *
     * @note In order to be able to perform a `boundary_face_operation` in the
     * MatrixFree::loop(), this field must be set to a value different from
     * UpdateFlags::update_default.
     */
    UpdateFlags mapping_update_flags_boundary_faces;

    /**
     * This flag determines the mapping data on interior faces to be
     * cached. Note that MatrixFree uses a separate loop layout for face
     * integrals in order to effectively vectorize also in the case of hanging
     * nodes (which require different subface settings on the two sides) or
     * some cells in the batch of a VectorizedArray of cells that are adjacent
     * to the boundary and others that are not.
     *
     * If set to a value different from update_general (default), the face
     * information is explicitly built. Currently, MatrixFree supports to
     * cache the following data on faces: inverse Jacobians, Jacobian
     * determinants (JxW), quadrature points, data for Hessians (derivative of
     * Jacobians), and normal vectors.
     *
     * @note In order to be able to perform a `inner_face_operation`
     * in the MatrixFree::loop(), this field must be set to a value different
     * from UpdateFlags::update_default.
     */
    UpdateFlags mapping_update_flags_inner_faces;

    /**
     * This flag determines the mapping data for faces in a different layout
     * with respect to vectorizations. Whereas
     * `mapping_update_flags_inner_faces` and
     * `mapping_update_flags_boundary_faces` trigger building the data in a
     * face-centric way with proper vectorization, the current data field
     * attaches the face information to the cells and their way of
     * vectorization. This is only needed in special situations, as for
     * example for block-Jacobi methods where the full operator to a cell
     * including its faces are evaluated. This data is accessed by
     * <code>FEFaceEvaluation::reinit(cell_batch_index,
     * face_number)</code>. However, currently no coupling terms to neighbors
     * can be computed with this approach because the neighbors are not laid
     * out by the VectorizedArray data layout with an
     * array-of-struct-of-array-type data structures.
     *
     * Note that you should only compute this data field in case you really
     * need it as it more than doubles the memory required by the mapping data
     * on faces.
     *
     * If set to a value different from update_general (default), the face
     * information is explicitly built. Currently, MatrixFree supports to
     * cache the following data on faces: inverse Jacobians, Jacobian
     * determinants (JxW), quadrature points, data for Hessians (derivative of
     * Jacobians), and normal vectors.
     */
    UpdateFlags mapping_update_flags_faces_by_cells;

    /**
     * This option can be used to define whether we work on a certain level of
     * the mesh, and not the active cells. If set to invalid_unsigned_int
     * (which is the default value), the active cells are gone through,
     * otherwise the level given by this parameter. Note that if you specify
     * to work on a level, its dofs must be distributed by using
     * <code>dof_handler.distribute_mg_dofs(fe);</code>.
     */
    unsigned int mg_level;

    /**
     * Controls whether to enable reading from vectors without resolving
     * constraints, i.e., just read the local values of the vector. By
     * default, this option is enabled. In case you want to use
     * FEEvaluationBase::read_dof_values_plain, this flag needs to be set.
     */
    bool store_plain_indices;

    /**
     * Option to control whether the indices stored in the DoFHandler
     * should be read and the pattern for task parallelism should be
     * set up in the initialize method of MatrixFree. The default
     * value is true. Can be disabled in case the mapping should be
     * recomputed (e.g. when using a deforming mesh described through
     * MappingEulerian) but the topology of cells has remained the
     * same.
     */
    bool initialize_indices;

    /**
     * Option to control whether the mapping information should be
     * computed in the initialize method of MatrixFree. The default
     * value is true. Can be disabled when only some indices should be
     * set up (e.g. when only a set of independent cells should be
     * computed).
     */
    bool initialize_mapping;

    /**
     * Option to control whether the loops should overlap communications and
     * computations as far as possible in case the vectors passed to the loops
     * support non-blocking data exchange. In most situations, overlapping is
     * faster in case the amount of data to be sent is more than a few
     * kilobytes. If less data is sent, the communication is latency bound on
     * most clusters (point-to-point latency is around 1 microsecond on good
     * clusters by 2016 standards). Depending on the MPI implementation and
     * the fabric, it may be faster to not overlap and wait for the data to
     * arrive. The default is true, i.e., communication and computation are
     * overlapped.
     */
    bool overlap_communication_computation;

    /**
     * By default, the face part will only hold those faces (and ghost
     * elements behind faces) that are going to be processed locally. In case
     * MatrixFree should have access to all neighbors on locally owned cells,
     * this option enables adding the respective faces at the end of the face
     * range.
     */
    bool hold_all_faces_to_owned_cells;

    /**
     * This data structure allows to assign a fraction of cells to different
     * categories when building the information for vectorization. It is used
     * implicitly when working with hp-adaptivity (with each active index
     * being a category) but can also be useful in other contexts where one
     * would like to control which cells together can form a batch of cells.
     * Such an example is "local time stepping", where cells of different
     * categories progress with different time-step sizes and, as a
     * consequence, can only processed together with cells with the same
     * category.
     *
     * This array is accessed by the number given by cell->active_cell_index()
     * when working on the active cells (with
     * @p mg_level set to numbers::invalid_unsigned_int) and by cell->index()
     * for the level cells.
     *
     * @note This field is empty upon construction of AdditionalData. It is
     * the responsibility of the user to resize this field to
     * `triangulation.n_active_cells()` or `triangulation.n_cells(level)` when
     * filling data.
     */
    std::vector<unsigned int> cell_vectorization_category;

    /**
     * By default, the different categories in @p cell_vectorization_category
     * can be mixed and the algorithm is allowed to merge lower categories with
     * the next higher categories if it is necessary inside the algorithm. This
     * gives a better utilization of the vectorization but might need special
     * treatment, in particular for face integrals. If set to @p true, the
     * algorithm will instead keep different categories separate and not mix
     * them in a single vectorized array.
     */
    bool cell_vectorization_categories_strict;

    /**
     * Assert that vectors passed to the MatrixFree loops are not ghosted.
     * This variable is primarily intended to reveal bugs or performance
     * problems caused by vectors that are involuntarily in ghosted mode,
     * by adding a check that this is not the case. In terms of correctness,
     * the MatrixFree::loop() and MatrixFree::cell_loop() methods support
     * both cases and perform similar operations. In particular, ghost values
     * are always updated on the source vector within the loop, and the
     * difference is only in whether the initial non-ghosted state is restored.
     */
    bool allow_ghosted_vectors_in_loops;

    /**
     * Option to control whether data should be generated on ghost cells.
     * If set to true, the data on ghost cells will be generated.
     * The default value is false.
     */
    bool store_ghost_cells;

    /**
     * Shared-memory MPI communicator. Default: MPI_COMM_SELF.
     */
    MPI_Comm communicator_sm;
  };

  /**
   * @name Construction and initialization
   */
  /** @{ */
  /**
   * Default empty constructor. Does nothing.
   */
  MatrixFree();

  /**
   * Copy constructor, calls copy_from
   */
  MatrixFree(const MatrixFree<dim, Number, VectorizedArrayType> &other);

  /**
   * Destructor.
   */
  ~MatrixFree() override = default;

  /**
   * Extracts the information needed to perform loops over cells. The
   * DoFHandler and AffineConstraints objects describe the layout of degrees
   * of freedom, the DoFHandler and the mapping describe the
   * transformations from unit to real cell, and the finite element
   * underlying the DoFHandler together with the quadrature formula
   * describe the local operations. Note that the finite element underlying
   * the DoFHandler must either be scalar or contain several copies of the
   * same element. Mixing several different elements into one FESystem is
   * not allowed. In that case, use the initialization function with
   * several DoFHandler arguments.
   */
  template <typename QuadratureType, typename number2, typename MappingType>
  void
  reinit(const MappingType                &mapping,
         const DoFHandler<dim>            &dof_handler,
         const AffineConstraints<number2> &constraint,
         const QuadratureType             &quad,
         const AdditionalData             &additional_data = AdditionalData());

  /**
   * Extracts the information needed to perform loops over cells. The
   * DoFHandler and AffineConstraints objects describe the layout of degrees of
   * freedom, the DoFHandler and the mapping describe the transformations from
   * unit to real cell, and the finite element underlying the DoFHandler
   * together with the quadrature formula describe the local operations. As
   * opposed to the scalar case treated with the other initialization
   * functions, this function allows for problems with two or more different
   * finite elements. The DoFHandlers to each element must be passed as
   * pointers to the initialization function. Alternatively, a system of
   * several components may also be represented by a single DoFHandler with an
   * FESystem element. The prerequisite for this case is that each base
   * element of the FESystem must be compatible with the present class, such
   * as the FE_Q or FE_DGQ classes.
   *
   * This function also allows for using several quadrature formulas, e.g.
   * when the description contains independent integrations of elements of
   * different degrees. However, the number of different quadrature formulas
   * can be sets independently from the number of DoFHandlers, when several
   * elements are always integrated with the same quadrature formula.
   */
  template <typename QuadratureType, typename number2, typename MappingType>
  void
  reinit(const MappingType                                     &mapping,
         const std::vector<const DoFHandler<dim> *>            &dof_handler,
         const std::vector<const AffineConstraints<number2> *> &constraint,
         const std::vector<QuadratureType>                     &quad,
         const AdditionalData &additional_data = AdditionalData());

  /**
   * Initializes the data structures. Same as before, but now the index set
   * description of the locally owned range of degrees of freedom is taken
   * from the DoFHandler. Moreover, only a single quadrature formula is used,
   * as might be necessary when several components in a vector-valued problem
   * are integrated together based on the same quadrature formula.
   */
  template <typename QuadratureType, typename number2, typename MappingType>
  void
  reinit(const MappingType                                     &mapping,
         const std::vector<const DoFHandler<dim> *>            &dof_handler,
         const std::vector<const AffineConstraints<number2> *> &constraint,
         const QuadratureType                                  &quad,
         const AdditionalData &additional_data = AdditionalData());

  /**
   * Copy function. Creates a deep copy of all data structures. It is usually
   * enough to keep the data for different operations once, so this function
   * should not be needed very often.
   */
  void
  copy_from(
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free_base);

  /**
   * Refreshes the geometry data stored in the MappingInfo fields when the
   * underlying geometry has changed (e.g. by a mapping that can deform
   * through a change in the spatial configuration like MappingFEField)
   * whereas the topology of the mesh and unknowns have remained the
   * same. Compared to reinit(), this operation only has to re-generate the
   * geometry arrays and can thus be significantly cheaper (depending on the
   * cost to evaluate the geometry).
   */
  void
  update_mapping(const Mapping<dim> &mapping);

  /**
   * Same as above but with hp::MappingCollection.
   */
  void
  update_mapping(const std::shared_ptr<hp::MappingCollection<dim>> &mapping);

  /**
   * Clear all data fields and brings the class into a condition similar to
   * after having called the default constructor.
   */
  void
  clear();

  /** @} */

  /**
   * This class defines the type of data access for face integrals in loop ()
   * that is passed on to the `update_ghost_values` and `compress` functions
   * of the parallel vectors, with the purpose of being able to reduce the
   * amount of data that must be exchanged. The data exchange is a real
   * bottleneck in particular for high-degree DG methods, therefore a more
   * restrictive way of exchange is clearly beneficial. Note that this
   * selection applies to FEFaceEvaluation objects assigned to the exterior
   * side of cells accessing `FaceToCellTopology::exterior_cells` only; all
   * <i>interior</i> objects are available in any case.
   */
  enum class DataAccessOnFaces
  {
    /**
     * The loop does not involve any FEFaceEvaluation access into neighbors,
     * as is the case with only boundary integrals (but no interior face
     * integrals) or when doing mass matrices in a MatrixFree::cell_loop()
     * like setup.
     */
    none,

    /**
     * The loop does only involve FEFaceEvaluation access into neighbors by
     * function values, such as FEFaceEvaluation::gather_evaluate() with
     * argument EvaluationFlags::values, but no access to shape function
     * derivatives (which typically need to access more data). For FiniteElement
     * types where only some of the shape functions have support on a face, such
     * as an FE_DGQ element with Lagrange polynomials with nodes on the element
     * surface, the data exchange is reduced from `(k+1)^dim` to
     * `(k+1)^(dim-1)`.
     */
    values,

    /**
     * Same as above. To be used if data has to be accessed from exterior faces
     * if FEFaceEvaluation was reinitialized by providing the cell batch number
     * and a face number. This configuration is useful in the context of
     * cell-centric loops.
     *
     * @pre AdditionalData::hold_all_faces_to_owned_cells has to enabled.
     */
    values_all_faces,

    /**
     * The loop does involve FEFaceEvaluation access into neighbors by
     * function values and gradients, but no second derivatives, such as
     * FEFaceEvaluation::gather_evaluate() with EvaluationFlags::values and
     * EvaluationFlags::gradients set. For FiniteElement types where only some
     * of the shape functions have non-zero value and first derivative on a
     * face, such as an FE_DGQHermite element, the data exchange is reduced,
     * e.g. from `(k+1)^dim` to `2(k+1)^(dim-1)`. Note that for bases that do
     * not have this special property, the full neighboring data is sent anyway.
     */
    gradients,

    /**
     * Same as above. To be used if data has to be accessed from exterior faces
     * if FEFaceEvaluation was reinitialized by providing the cell batch number
     * and a face number. This configuration is useful in the context of
     * cell-centric loops.
     *
     * @pre AdditionalData::hold_all_faces_to_owned_cells has to enabled.
     */
    gradients_all_faces,

    /**
     * General setup where the user does not want to make a restriction. This
     * is typically more expensive than the other options, but also the most
     * conservative one because the full data of elements behind the faces to
     * be computed locally will be exchanged.
     */
    unspecified
  };

  /**
   * @name Matrix-free loops
   */
  /** @{ */
  /**
   * This method runs the loop over all cells (in parallel) and performs the
   * MPI data exchange on the source vector and destination vector.
   *
   * @param cell_operation `std::function` with the signature <tt>cell_operation
   * (const MatrixFree<dim,Number> &, OutVector &, InVector &,
   * std::pair<unsigned int,unsigned int> &)</tt> where the first argument
   * passes the data of the calling class and the last argument defines the
   * range of cells which should be worked on (typically more than one cell
   * should be worked on in order to reduce overheads).  One can pass a pointer
   * to an object in this place if it has an `operator()` with the correct set
   * of arguments since such a pointer can be converted to the function object.
   *
   * @param dst Destination vector holding the result. If the vector is of
   * type LinearAlgebra::distributed::Vector (or composite objects thereof
   * such as LinearAlgebra::distributed::BlockVector), the loop calls
   * LinearAlgebra::distributed::Vector::compress() at the end of the call
   * internally. For other vectors, including parallel Trilinos or PETSc
   * vectors, no such call is issued. Note that Trilinos/Epetra or PETSc
   * vectors do currently not work in parallel because the present class uses
   * MPI-local index addressing, as opposed to the global addressing implied
   * by those external libraries.
   *
   * @param src Input vector. If the vector is of type
   * LinearAlgebra::distributed::Vector (or composite objects thereof such as
   * LinearAlgebra::distributed::BlockVector), the loop calls
   * LinearAlgebra::distributed::Vector::update_ghost_values() at the start of
   * the call internally to make sure all necessary data is locally
   * available. Note, however, that the vector is reset to its original state
   * at the end of the loop, i.e., if the vector was not ghosted upon entry of
   * the loop, it will not be ghosted upon finishing the loop.
   *
   * @param zero_dst_vector If this flag is set to `true`, the vector `dst`
   * will be set to zero inside the loop. Use this case in case you perform a
   * typical `vmult()` operation on a matrix object, as it will typically be
   * faster than calling `dst = 0;` before the loop separately. This is
   * because the vector entries are set to zero only on subranges of the
   * vector, making sure that the vector entries stay in caches as much as
   * possible.
   */
  template <typename OutVector, typename InVector>
  void
  cell_loop(const std::function<void(
              const MatrixFree<dim, Number, VectorizedArrayType> &,
              OutVector &,
              const InVector &,
              const std::pair<unsigned int, unsigned int> &)> &cell_operation,
            OutVector                                         &dst,
            const InVector                                    &src,
            const bool zero_dst_vector = false) const;

  /**
   * This is the second variant to run the loop over all cells, now providing
   * a function pointer to a member function of class `CLASS`. This method
   * obviates the need to define a lambda function or to call std::bind to bind
   * the class into the given
   * function in case the local function needs to access data in the class
   * (i.e., it is a non-static member function).
   *
   * @param cell_operation Pointer to member function of `CLASS` with the
   * signature <tt>cell_operation (const MatrixFree<dim,Number> &, OutVector &,
   * InVector &, std::pair<unsigned int,unsigned int> &)</tt> where the first
   * argument passes the data of the calling class and the last argument
   * defines the range of cells which should be worked on (typically more than
   * one cell should be worked on in order to reduce overheads).
   *
   * @param owning_class The object which provides the `cell_operation`
   * call. To be compatible with this interface, the class must allow to call
   * `owning_class->cell_operation(...)`.
   *
   * @param dst Destination vector holding the result. If the vector is of
   * type LinearAlgebra::distributed::Vector (or composite objects thereof
   * such as LinearAlgebra::distributed::BlockVector), the loop calls
   * LinearAlgebra::distributed::Vector::compress() at the end of the call
   * internally. For other vectors, including parallel Trilinos or PETSc
   * vectors, no such call is issued. Note that Trilinos/Epetra or PETSc
   * vectors do currently not work in parallel because the present class uses
   * MPI-local index addressing, as opposed to the global addressing implied
   * by those external libraries.
   *
   * @param src Input vector. If the vector is of type
   * LinearAlgebra::distributed::Vector (or composite objects thereof such as
   * LinearAlgebra::distributed::BlockVector), the loop calls
   * LinearAlgebra::distributed::Vector::update_ghost_values() at the start of
   * the call internally to make sure all necessary data is locally
   * available. Note, however, that the vector is reset to its original state
   * at the end of the loop, i.e., if the vector was not ghosted upon entry of
   * the loop, it will not be ghosted upon finishing the loop.
   *
   * @param zero_dst_vector If this flag is set to `true`, the vector `dst`
   * will be set to zero inside the loop. Use this case in case you perform a
   * typical `vmult()` operation on a matrix object, as it will typically be
   * faster than calling `dst = 0;` before the loop separately. This is
   * because the vector entries are set to zero only on subranges of the
   * vector, making sure that the vector entries stay in caches as much as
   * possible.
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  cell_loop(void (CLASS::*cell_operation)(
              const MatrixFree &,
              OutVector &,
              const InVector &,
              const std::pair<unsigned int, unsigned int> &) const,
            const CLASS    *owning_class,
            OutVector      &dst,
            const InVector &src,
            const bool      zero_dst_vector = false) const;

  /**
   * Same as above, but for class member functions which are non-const.
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  cell_loop(void (CLASS::*cell_operation)(
              const MatrixFree &,
              OutVector &,
              const InVector &,
              const std::pair<unsigned int, unsigned int> &),
            CLASS          *owning_class,
            OutVector      &dst,
            const InVector &src,
            const bool      zero_dst_vector = false) const;

  /**
   * This function is similar to the cell_loop with an std::function object to
   * specify to operation to be performed on cells, but adds two additional
   * functors to execute some additional work before and after the cell
   * integrals are computed.
   *
   * The two additional functors work on a range of degrees of freedom,
   * expressed in terms of the degree-of-freedom numbering of the selected
   * DoFHandler `dof_handler_index_pre_post` in MPI-local indices. The
   * arguments to the functors represent a range of degrees of freedom at a
   * granularity of
   * internal::MatrixFreeFunctions::DoFInfo::chunk_size_zero_vector entries
   * (except for the last chunk which is set to the number of locally owned
   * entries) in the form `[first, last)`. The idea of these functors is to
   * bring operations on vectors closer to the point where they accessed in a
   * matrix-free loop, with the goal to increase cache hits by temporal
   * locality. This loop guarantees that the `operation_before_loop` hits all
   * relevant unknowns before they are first touched in the cell_operation
   * (including the MPI data exchange), allowing to execute some vector update
   * that the `src` vector depends upon. The `operation_after_loop` is similar
   * - it starts to execute on a range of DoFs once all DoFs in that range
   * have been touched for the last time by the `cell_operation`
   * (including the MPI data exchange), allowing e.g. to compute some vector
   * operations that depend on the result of the current cell loop in `dst` or
   * want to modify `src`. The efficiency of caching depends on the numbering
   * of the degrees of freedom because of the granularity of the ranges.
   *
   * @param cell_operation Pointer to member function of `CLASS` with the
   * signature <tt>cell_operation (const MatrixFree<dim,Number> &, OutVector &,
   * InVector &, std::pair<unsigned int,unsigned int> &)</tt> where the first
   * argument passes the data of the calling class and the last argument
   * defines the range of cells which should be worked on (typically more than
   * one cell should be worked on in order to reduce overheads).
   *
   * @param owning_class The object which provides the `cell_operation`
   * call. To be compatible with this interface, the class must allow to call
   * `owning_class->cell_operation(...)`.
   *
   * @param dst Destination vector holding the result. If the vector is of
   * type LinearAlgebra::distributed::Vector (or composite objects thereof
   * such as LinearAlgebra::distributed::BlockVector), the loop calls
   * LinearAlgebra::distributed::Vector::compress() at the end of the call
   * internally. For other vectors, including parallel Trilinos or PETSc
   * vectors, no such call is issued. Note that Trilinos/Epetra or PETSc
   * vectors do currently not work in parallel because the present class uses
   * MPI-local index addressing, as opposed to the global addressing implied
   * by those external libraries.
   *
   * @param src Input vector. If the vector is of type
   * LinearAlgebra::distributed::Vector (or composite objects thereof such as
   * LinearAlgebra::distributed::BlockVector), the loop calls
   * LinearAlgebra::distributed::Vector::update_ghost_values() at the start of
   * the call internally to make sure all necessary data is locally
   * available. Note, however, that the vector is reset to its original state
   * at the end of the loop, i.e., if the vector was not ghosted upon entry of
   * the loop, it will not be ghosted upon finishing the loop.
   *
   * @param operation_before_loop This functor can be used to perform an
   * operation on entries of the `src` and `dst` vectors (or other vectors)
   * before the operation on cells first touches a particular DoF according to
   * the general description in the text above. This function is passed a
   * range of the locally owned degrees of freedom on the selected
   * `dof_handler_index_pre_post` (in MPI-local numbering).
   *
   * @param operation_after_loop This functor can be used to perform an
   * operation on entries of the `src` and `dst` vectors (or other vectors)
   * after the operation on cells last touches a particular DoF according to
   * the general description in the text above. This function is passed a
   * range of the locally owned degrees of freedom on the selected
   * `dof_handler_index_pre_post` (in MPI-local numbering).
   *
   * @param dof_handler_index_pre_post Since MatrixFree can be initialized
   * with a vector of DoFHandler objects, each of them will in general have
   * different vector sizes and thus different ranges returned to
   * `operation_before_loop` and `operation_after_loop`. Use this variable to
   * specify which one of the DoFHandler objects the index range should be
   * associated to. Defaults to the `dof_handler_index` 0.
   *
   * @note The close locality of the `operation_before_loop` and
   * `operation_after_loop` is currently only implemented for the MPI-only
   * case. In case threading is enabled, the complete `operation_before_loop`
   * is scheduled before the parallel loop, and `operation_after_loop` is
   * scheduled strictly afterwards, due to the complicated dependencies.
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  cell_loop(void (CLASS::*cell_operation)(
              const MatrixFree &,
              OutVector &,
              const InVector &,
              const std::pair<unsigned int, unsigned int> &) const,
            const CLASS    *owning_class,
            OutVector      &dst,
            const InVector &src,
            const std::function<void(const unsigned int, const unsigned int)>
              &operation_before_loop,
            const std::function<void(const unsigned int, const unsigned int)>
                              &operation_after_loop,
            const unsigned int dof_handler_index_pre_post = 0) const;

  /**
   * Same as above, but for class member functions which are non-const.
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  cell_loop(void (CLASS::*cell_operation)(
              const MatrixFree &,
              OutVector &,
              const InVector &,
              const std::pair<unsigned int, unsigned int> &),
            CLASS          *owning_class,
            OutVector      &dst,
            const InVector &src,
            const std::function<void(const unsigned int, const unsigned int)>
              &operation_before_loop,
            const std::function<void(const unsigned int, const unsigned int)>
                              &operation_after_loop,
            const unsigned int dof_handler_index_pre_post = 0) const;

  /**
   * Same as above, but taking an `std::function` as the `cell_operation`
   * rather than a class member function.
   */
  template <typename OutVector, typename InVector>
  void
  cell_loop(const std::function<void(
              const MatrixFree<dim, Number, VectorizedArrayType> &,
              OutVector &,
              const InVector &,
              const std::pair<unsigned int, unsigned int> &)> &cell_operation,
            OutVector                                         &dst,
            const InVector                                    &src,
            const std::function<void(const unsigned int, const unsigned int)>
              &operation_before_loop,
            const std::function<void(const unsigned int, const unsigned int)>
                              &operation_after_loop,
            const unsigned int dof_handler_index_pre_post = 0) const;

  /**
   * This method runs a loop over all cells (in parallel) and performs the MPI
   * data exchange on the source vector and destination vector. As opposed to
   * the other variants that only runs a function on cells, this method also
   * takes as arguments a function for the interior faces and for the boundary
   * faces, respectively.
   *
   * @param cell_operation `std::function` with the signature <tt>cell_operation
   * (const MatrixFree<dim,Number> &, OutVector &, InVector &,
   * std::pair<unsigned int,unsigned int> &)</tt> where the first argument
   * passes the data of the calling class and the last argument defines the
   * range of cells which should be worked on (typically more than one cell
   * should be worked on in order to reduce overheads). One can pass a pointer
   * to an object in this place if it has an <code>operator()</code> with the
   * correct set of arguments since such a pointer can be converted to the
   * function object.
   *
   * @param inner_face_operation `std::function` with the signature <tt>inner_face_operation
   * (const MatrixFree<dim,Number> &, OutVector &, InVector &,
   * std::pair<unsigned int,unsigned int> &)</tt> in analogy to
   * `cell_operation`, but now the part associated to the work on interior
   * faces. Note that the MatrixFree framework treats periodic faces as interior
   * ones, so they will be assigned their correct neighbor after applying
   * periodicity constraints within the inner_face_operation calls.
   *
   * @param boundary_face_operation `std::function` with the signature
   * <tt>boundary_face_operation (const MatrixFree<dim,Number> &, OutVector &,
   * InVector &, std::pair<unsigned int,unsigned int> &)</tt> in analogy to
   * `cell_operation` and `inner_face_operation`, but now the part
   * associated to the work on boundary faces. Boundary faces are separated by
   * their `boundary_id` and it is possible to query that id using
   * MatrixFree::get_boundary_id(). Note that both interior and faces use the
   * same numbering, and faces in the interior are assigned lower numbers than
   * the boundary faces.
   *
   * @param dst Destination vector holding the result. If the vector is of
   * type LinearAlgebra::distributed::Vector (or composite objects thereof
   * such as LinearAlgebra::distributed::BlockVector), the loop calls
   * LinearAlgebra::distributed::Vector::compress() at the end of the call
   * internally.
   *
   * @param src Input vector. If the vector is of type
   * LinearAlgebra::distributed::Vector (or composite objects thereof such as
   * LinearAlgebra::distributed::BlockVector), the loop calls
   * LinearAlgebra::distributed::Vector::update_ghost_values() at the start of
   * the call internally to make sure all necessary data is locally
   * available. Note, however, that the vector is reset to its original state
   * at the end of the loop, i.e., if the vector was not ghosted upon entry of
   * the loop, it will not be ghosted upon finishing the loop.
   *
   * @param zero_dst_vector If this flag is set to `true`, the vector `dst`
   * will be set to zero inside the loop. Use this case in case you perform a
   * typical `vmult()` operation on a matrix object, as it will typically be
   * faster than calling `dst = 0;` before the loop separately. This is
   * because the vector entries are set to zero only on subranges of the
   * vector, making sure that the vector entries stay in caches as much as
   * possible.
   *
   * @param dst_vector_face_access Set the type of access into the vector
   * `dst` that will happen inside the body of the @p inner_face_operation
   * function. As explained in the description of the DataAccessOnFaces
   * struct, the purpose of this selection is to reduce the amount of data
   * that must be exchanged over the MPI network (or via `memcpy` if within
   * the shared memory region of a node) to gain performance. Note that there
   * is no way to communicate this setting with the FEFaceEvaluation class,
   * therefore this selection must be made at this site in addition to what is
   * implemented inside the `inner_face_operation` function. As a
   * consequence, there is also no way to check that the setting passed to this
   * call is consistent with what is later done by `FEFaceEvaluation`, and it is
   * the user's responsibility to ensure correctness of data.
   *
   * @param src_vector_face_access Set the type of access into the vector
   * `src` that will happen inside the body of the @p inner_face_operation function,
   * in analogy to `dst_vector_face_access`.
   */
  template <typename OutVector, typename InVector>
  void
  loop(
    const std::function<
      void(const MatrixFree<dim, Number, VectorizedArrayType> &,
           OutVector &,
           const InVector &,
           const std::pair<unsigned int, unsigned int> &)> &cell_operation,
    const std::function<void(
      const MatrixFree<dim, Number, VectorizedArrayType> &,
      OutVector &,
      const InVector &,
      const std::pair<unsigned int, unsigned int> &)> &inner_face_operation,
    const std::function<void(
      const MatrixFree<dim, Number, VectorizedArrayType> &,
      OutVector &,
      const InVector &,
      const std::pair<unsigned int, unsigned int> &)> &boundary_face_operation,
    OutVector                                         &dst,
    const InVector                                    &src,
    const bool                                         zero_dst_vector = false,
    const DataAccessOnFaces                            dst_vector_face_access =
      DataAccessOnFaces::unspecified,
    const DataAccessOnFaces src_vector_face_access =
      DataAccessOnFaces::unspecified) const;

  /**
   * This is the second variant to run the loop over all cells, interior
   * faces, and boundary faces, now providing three function pointers to
   * member functions of class @p CLASS with the signature <code>operation
   * (const MatrixFree<dim,Number> &, OutVector &, InVector &,
   * std::pair<unsigned int,unsigned int>&)const</code>. This method obviates
   * the need to define a lambda function or to call std::bind to bind
   * the class into the given
   * function in case the local function needs to access data in the class
   * (i.e., it is a non-static member function).
   *
   * @param cell_operation Pointer to member function of `CLASS` with the
   * signature <tt>cell_operation (const MatrixFree<dim,Number> &, OutVector &,
   * InVector &, std::pair<unsigned int,unsigned int> &)</tt> where the first
   * argument passes the data of the calling class and the last argument
   * defines the range of cells which should be worked on (typically more than
   * one cell should be worked on in order to reduce overheads). Note that the
   * loop will typically split the `cell_range` into smaller pieces and work
   * on `cell_operation`, `inner_face_operation`, and
   * `boundary_face_operation` alternately, in order to increase the potential
   * reuse of vector entries in caches.
   *
   * @param inner_face_operation Pointer to member function of `CLASS` with the
   * signature <tt>inner_face_operation (const MatrixFree<dim,Number> &,
   * OutVector &, InVector &, std::pair<unsigned int,unsigned int> &)</tt> in
   * analogy to `cell_operation`, but now the part associated to the work on
   * interior faces. Note that the MatrixFree framework treats periodic faces as
   * interior ones, so they will be assigned their correct neighbor after
   * applying periodicity constraints within the inner_face_operation
   * calls.
   *
   * @param boundary_face_operation Pointer to member function of `CLASS` with the
   * signature <tt>boundary_face_operation (const MatrixFree<dim,Number> &,
   * OutVector
   * &, InVector &, std::pair<unsigned int,unsigned int> &)</tt> in analogy to
   * `cell_operation` and `inner_face_operation`, but now the part
   * associated to the work on boundary faces. Boundary faces are separated by
   * their `boundary_id` and it is possible to query that id using
   * MatrixFree::get_boundary_id(). Note that both interior and faces use the
   * same numbering, and faces in the interior are assigned lower numbers than
   * the boundary faces.
   *
   * @param owning_class The object which provides the `cell_operation`
   * call. To be compatible with this interface, the class must allow to call
   * `owning_class->cell_operation(...)`,
   * `owning_class->inner_face_operation(...)`, and
   * `owning_class->boundary_face_operation(...)`.
   *
   * @param dst Destination vector holding the result. If the vector is of
   * type LinearAlgebra::distributed::Vector (or composite objects thereof
   * such as LinearAlgebra::distributed::BlockVector), the loop calls
   * LinearAlgebra::distributed::Vector::compress() at the end of the call
   * internally.
   *
   * @param src Input vector. If the vector is of type
   * LinearAlgebra::distributed::Vector (or composite objects thereof such as
   * LinearAlgebra::distributed::BlockVector), the loop calls
   * LinearAlgebra::distributed::Vector::update_ghost_values() at the start of
   * the call internally to make sure all necessary data is locally
   * available. Note, however, that the vector is reset to its original state
   * at the end of the loop, i.e., if the vector was not ghosted upon entry of
   * the loop, it will not be ghosted upon finishing the loop.
   *
   * @param zero_dst_vector If this flag is set to `true`, the vector `dst`
   * will be set to zero inside the loop. Use this case in case you perform a
   * typical `vmult()` operation on a matrix object, as it will typically be
   * faster than calling `dst = 0;` before the loop separately. This is
   * because the vector entries are set to zero only on subranges of the
   * vector, making sure that the vector entries stay in caches as much as
   * possible.
   *
   * @param dst_vector_face_access Set the type of access into the vector
   * `dst` that will happen inside the body of the @p inner_face_operation
   * function. As explained in the description of the DataAccessOnFaces
   * struct, the purpose of this selection is to reduce the amount of data
   * that must be exchanged over the MPI network (or via `memcpy` if within
   * the shared memory region of a node) to gain performance. Note that there
   * is no way to communicate this setting with the FEFaceEvaluation class,
   * therefore this selection must be made at this site in addition to what is
   * implemented inside the `inner_face_operation` function. As a
   * consequence, there is also no way to check that the setting passed to this
   * call is consistent with what is later done by `FEFaceEvaluation`, and it is
   * the user's responsibility to ensure correctness of data.
   *
   * @param src_vector_face_access Set the type of access into the vector
   * `src` that will happen inside the body of the @p inner_face_operation function,
   * in analogy to `dst_vector_face_access`.
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  loop(void (CLASS::*cell_operation)(
         const MatrixFree &,
         OutVector &,
         const InVector &,
         const std::pair<unsigned int, unsigned int> &) const,
       void (CLASS::*inner_face_operation)(
         const MatrixFree &,
         OutVector &,
         const InVector &,
         const std::pair<unsigned int, unsigned int> &) const,
       void (CLASS::*boundary_face_operation)(
         const MatrixFree &,
         OutVector &,
         const InVector &,
         const std::pair<unsigned int, unsigned int> &) const,
       const CLASS            *owning_class,
       OutVector              &dst,
       const InVector         &src,
       const bool              zero_dst_vector = false,
       const DataAccessOnFaces dst_vector_face_access =
         DataAccessOnFaces::unspecified,
       const DataAccessOnFaces src_vector_face_access =
         DataAccessOnFaces::unspecified) const;

  /**
   * Same as above, but for class member functions which are non-const.
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  loop(void (CLASS::*cell_operation)(
         const MatrixFree &,
         OutVector &,
         const InVector &,
         const std::pair<unsigned int, unsigned int> &),
       void (CLASS::*inner_face_operation)(
         const MatrixFree &,
         OutVector &,
         const InVector &,
         const std::pair<unsigned int, unsigned int> &),
       void (CLASS::*boundary_face_operation)(
         const MatrixFree &,
         OutVector &,
         const InVector &,
         const std::pair<unsigned int, unsigned int> &),
       CLASS                  *owning_class,
       OutVector              &dst,
       const InVector         &src,
       const bool              zero_dst_vector = false,
       const DataAccessOnFaces dst_vector_face_access =
         DataAccessOnFaces::unspecified,
       const DataAccessOnFaces src_vector_face_access =
         DataAccessOnFaces::unspecified) const;

  /**
   * This function is similar to the loop method above, but adds two additional
   * functors to execute some additional work before and after the cell, face
   * and boundary integrals are computed.
   *
   * The two additional functors work on a range of degrees of freedom,
   * expressed in terms of the degree-of-freedom numbering of the selected
   * DoFHandler `dof_handler_index_pre_post` in MPI-local indices. The
   * arguments to the functors represent a range of degrees of freedom at a
   * granularity of
   * internal::MatrixFreeFunctions::DoFInfo::chunk_size_zero_vector entries
   * (except for the last chunk which is set to the number of locally owned
   * entries) in the form `[first, last)`. The idea of these functors is to
   * bring operations on vectors closer to the point where they accessed in a
   * matrix-free loop, with the goal to increase cache hits by temporal
   * locality. This loop guarantees that the `operation_before_loop` hits all
   * relevant unknowns before they are first touched by any of the cell, face or
   * boundary operations (including the MPI data exchange), allowing to execute
   * some vector update that the `src` vector depends upon. The
   * `operation_after_loop` is similar - it starts to execute on a range of DoFs
   * once all DoFs in that range have been touched for the last time by the
   * cell, face and boundary operations (including the MPI data exchange),
   * allowing e.g. to compute some vector operations that depend on the result
   * of the current cell loop in `dst` or want to modify `src`. The efficiency
   * of caching depends on the numbering of the degrees of freedom because of
   * the granularity of the ranges.
   *
   * @param cell_operation Pointer to member function of `CLASS` with the
   * signature <tt>cell_operation (const MatrixFree<dim,Number> &, OutVector &,
   * InVector &, std::pair<unsigned int,unsigned int> &)</tt> where the first
   * argument passes the data of the calling class and the last argument
   * defines the range of cells which should be worked on (typically more than
   * one cell should be worked on in order to reduce overheads).
   *
   * @param inner_face_operation Pointer to member function of `CLASS` with the
   * signature <tt>inner_face_operation (const MatrixFree<dim,Number> &,
   * OutVector &, InVector &, std::pair<unsigned int,unsigned int> &)</tt> in
   * analogy to `cell_operation`, but now the part associated to the work on
   * interior faces. Note that the MatrixFree framework treats periodic faces as
   * interior ones, so they will be assigned their correct neighbor after
   * applying periodicity constraints within the inner_face_operation
   * calls.
   *
   * @param boundary_face_operation Pointer to member function of `CLASS` with the
   * signature <tt>boundary_face_operation (const MatrixFree<dim,Number> &,
   * OutVector
   * &, InVector &, std::pair<unsigned int,unsigned int> &)</tt> in analogy to
   * `cell_operation` and `inner_face_operation`, but now the part
   * associated to the work on boundary faces. Boundary faces are separated by
   * their `boundary_id` and it is possible to query that id using
   * MatrixFree::get_boundary_id(). Note that both interior and faces use the
   * same numbering, and faces in the interior are assigned lower numbers than
   * the boundary faces.
   *
   * @param owning_class The object which provides the `cell_operation`
   * call. To be compatible with this interface, the class must allow to call
   * `owning_class->cell_operation(...)`.
   *
   * @param dst Destination vector holding the result. If the vector is of
   * type LinearAlgebra::distributed::Vector (or composite objects thereof
   * such as LinearAlgebra::distributed::BlockVector), the loop calls
   * LinearAlgebra::distributed::Vector::compress() at the end of the call
   * internally. For other vectors, including parallel Trilinos or PETSc
   * vectors, no such call is issued. Note that Trilinos/Epetra or PETSc
   * vectors do currently not work in parallel because the present class uses
   * MPI-local index addressing, as opposed to the global addressing implied
   * by those external libraries.
   *
   * @param src Input vector. If the vector is of type
   * LinearAlgebra::distributed::Vector (or composite objects thereof such as
   * LinearAlgebra::distributed::BlockVector), the loop calls
   * LinearAlgebra::distributed::Vector::update_ghost_values() at the start of
   * the call internally to make sure all necessary data is locally
   * available. Note, however, that the vector is reset to its original state
   * at the end of the loop, i.e., if the vector was not ghosted upon entry of
   * the loop, it will not be ghosted upon finishing the loop.
   *
   * @param operation_before_loop This functor can be used to perform an
   * operation on entries of the `src` and `dst` vectors (or other vectors)
   * before the operation on cells first touches a particular DoF according to
   * the general description in the text above. This function is passed a
   * range of the locally owned degrees of freedom on the selected
   * `dof_handler_index_pre_post` (in MPI-local numbering).
   *
   * @param operation_after_loop This functor can be used to perform an
   * operation on entries of the `src` and `dst` vectors (or other vectors)
   * after the operation on cells last touches a particular DoF according to
   * the general description in the text above. This function is passed a
   * range of the locally owned degrees of freedom on the selected
   * `dof_handler_index_pre_post` (in MPI-local numbering).
   *
   * @param dof_handler_index_pre_post Since MatrixFree can be initialized
   * with a vector of DoFHandler objects, each of them will in general have
   * different vector sizes and thus different ranges returned to
   * `operation_before_loop` and `operation_after_loop`. Use this variable to
   * specify which one of the DoFHandler objects the index range should be
   * associated to. Defaults to the `dof_handler_index` 0.
   *
   * @param dst_vector_face_access Set the type of access into the vector
   * `dst` that will happen inside the body of the @p inner_face_operation
   * function. As explained in the description of the DataAccessOnFaces
   * struct, the purpose of this selection is to reduce the amount of data
   * that must be exchanged over the MPI network (or via `memcpy` if within
   * the shared memory region of a node) to gain performance. Note that there
   * is no way to communicate this setting with the FEFaceEvaluation class,
   * therefore this selection must be made at this site in addition to what is
   * implemented inside the `inner_face_operation` function. As a
   * consequence, there is also no way to check that the setting passed to this
   * call is consistent with what is later done by `FEFaceEvaluation`, and it is
   * the user's responsibility to ensure correctness of data.
   *
   * @param src_vector_face_access Set the type of access into the vector
   * `src` that will happen inside the body of the @p inner_face_operation function,
   * in analogy to `dst_vector_face_access`.
   *
   * @note The close locality of the `operation_before_loop` and
   * `operation_after_loop` is currently only implemented for the MPI-only
   * case. In case threading is enabled, the complete `operation_before_loop`
   * is scheduled before the parallel loop, and `operation_after_loop` is
   * scheduled strictly afterwards, due to the complicated dependencies.
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  loop(void (CLASS::*cell_operation)(
         const MatrixFree &,
         OutVector &,
         const InVector &,
         const std::pair<unsigned int, unsigned int> &) const,
       void (CLASS::*inner_face_operation)(
         const MatrixFree &,
         OutVector &,
         const InVector &,
         const std::pair<unsigned int, unsigned int> &) const,
       void (CLASS::*boundary_face_operation)(
         const MatrixFree &,
         OutVector &,
         const InVector &,
         const std::pair<unsigned int, unsigned int> &) const,
       const CLASS    *owning_class,
       OutVector      &dst,
       const InVector &src,
       const std::function<void(const unsigned int, const unsigned int)>
         &operation_before_loop,
       const std::function<void(const unsigned int, const unsigned int)>
                              &operation_after_loop,
       const unsigned int      dof_handler_index_pre_post = 0,
       const DataAccessOnFaces dst_vector_face_access =
         DataAccessOnFaces::unspecified,
       const DataAccessOnFaces src_vector_face_access =
         DataAccessOnFaces::unspecified) const;

  /**
   * Same as above, but for class member functions which are non-const.
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  loop(void (CLASS::*cell_operation)(
         const MatrixFree &,
         OutVector &,
         const InVector &,
         const std::pair<unsigned int, unsigned int> &),
       void (CLASS::*inner_face_operation)(
         const MatrixFree &,
         OutVector &,
         const InVector &,
         const std::pair<unsigned int, unsigned int> &),
       void (CLASS::*boundary_face_operation)(
         const MatrixFree &,
         OutVector &,
         const InVector &,
         const std::pair<unsigned int, unsigned int> &),
       const CLASS    *owning_class,
       OutVector      &dst,
       const InVector &src,
       const std::function<void(const unsigned int, const unsigned int)>
         &operation_before_loop,
       const std::function<void(const unsigned int, const unsigned int)>
                              &operation_after_loop,
       const unsigned int      dof_handler_index_pre_post = 0,
       const DataAccessOnFaces dst_vector_face_access =
         DataAccessOnFaces::unspecified,
       const DataAccessOnFaces src_vector_face_access =
         DataAccessOnFaces::unspecified) const;

  /**
   * Same as above, but taking an `std::function` as the `cell_operation`,
   * `inner_face_operation` and `boundary_face_operation` rather than a
   * class member function.
   */
  template <typename OutVector, typename InVector>
  void
  loop(
    const std::function<
      void(const MatrixFree<dim, Number, VectorizedArrayType> &,
           OutVector &,
           const InVector &,
           const std::pair<unsigned int, unsigned int> &)> &cell_operation,
    const std::function<void(
      const MatrixFree<dim, Number, VectorizedArrayType> &,
      OutVector &,
      const InVector &,
      const std::pair<unsigned int, unsigned int> &)> &inner_face_operation,
    const std::function<void(
      const MatrixFree<dim, Number, VectorizedArrayType> &,
      OutVector &,
      const InVector &,
      const std::pair<unsigned int, unsigned int> &)> &boundary_face_operation,
    OutVector                                         &dst,
    const InVector                                    &src,
    const std::function<void(const unsigned int, const unsigned int)>
      &operation_before_loop,
    const std::function<void(const unsigned int, const unsigned int)>
                           &operation_after_loop,
    const unsigned int      dof_handler_index_pre_post = 0,
    const DataAccessOnFaces dst_vector_face_access =
      DataAccessOnFaces::unspecified,
    const DataAccessOnFaces src_vector_face_access =
      DataAccessOnFaces::unspecified) const;

  /**
   * This method runs the loop over all cells (in parallel) similarly as
   * cell_loop() does. However, this function is intended to be used
   * for the case if face and boundary integrals should be also
   * evaluated. In contrast to loop(), the user provides only a single function
   * that should contain the cell integral over a cell (or batch of cells when
   * vectorizing) and the face and boundary integrals over all its faces. This
   * is referred to in the literature as `element-centric loop` or `cell-centric
   * loop`.
   *
   * To be able to evaluate all face integrals (with values or gradients
   * from the neighboring cells), all ghost values from neighboring cells are
   * updated. Use
   * FEFaceEvaluation::reinit(cell, face_no) to access quantities on arbitrary
   * faces of a cell and the respective neighbors.
   *
   * @param cell_operation Pointer to member function of `CLASS` with the
   * signature <tt>cell_operation (const MatrixFree<dim,Number> &, OutVector &,
   * InVector &, std::pair<unsigned int,unsigned int> &)</tt> where the first
   * argument passes the data of the calling class and the last argument
   * defines the range of cells which should be worked on (typically more than
   * one cell is passed in from the loop in order to reduce overheads).
   *
   * @param owning_class The object which provides the `cell_operation`
   * call. To be compatible with this interface, the class must allow to call
   * `owning_class->cell_operation(...)`.
   *
   * @param dst Destination vector holding the result. If the vector is of
   * type LinearAlgebra::distributed::Vector (or composite objects thereof
   * such as LinearAlgebra::distributed::BlockVector), the loop calls
   * LinearAlgebra::distributed::Vector::compress() at the end of the call
   * internally.
   *
   * @param src Input vector. If the vector is of type
   * LinearAlgebra::distributed::Vector (or composite objects thereof such as
   * LinearAlgebra::distributed::BlockVector), the loop calls
   * LinearAlgebra::distributed::Vector::update_ghost_values() at the start of
   * the call internally to make sure all necessary data is locally
   * available. Note, however, that the vector is reset to its original state
   * at the end of the loop, i.e., if the vector was not ghosted upon entry of
   * the loop, it will not be ghosted upon finishing the loop.
   *
   * @param zero_dst_vector If this flag is set to `true`, the vector `dst`
   * will be set to zero inside the loop. Use this case in case you perform a
   * typical `vmult()` operation on a matrix object, as it will typically be
   * faster than calling `dst = 0;` before the loop separately. This is
   * because the vector entries are set to zero only on subranges of the
   * vector, making sure that the vector entries stay in caches as much as
   * possible.
   *
   * @param src_vector_face_access Set the type of access into the vector
   * `src` that will happen inside the body of the @p cell_operation function
   * during face integrals.
   * As explained in the description of the DataAccessOnFaces
   * struct, the purpose of this selection is to reduce the amount of data
   * that must be exchanged over the MPI network (or via `memcpy` if within
   * the shared memory region of a node) to gain performance. Note that there
   * is no way to communicate this setting with the FEFaceEvaluation class,
   * therefore this selection must be made at this site in addition to what is
   * implemented inside the `inner_face_operation` function. As a
   * consequence, there is also no way to check that the setting passed to this
   * call is consistent with what is later done by `FEFaceEvaluation`, and it is
   * the user's responsibility to ensure correctness of data.
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  loop_cell_centric(void (CLASS::*cell_operation)(
                      const MatrixFree &,
                      OutVector &,
                      const InVector &,
                      const std::pair<unsigned int, unsigned int> &) const,
                    const CLASS            *owning_class,
                    OutVector              &dst,
                    const InVector         &src,
                    const bool              zero_dst_vector = false,
                    const DataAccessOnFaces src_vector_face_access =
                      DataAccessOnFaces::unspecified) const;

  /**
   * Same as above, but for the class member function which is non-const.
   */
  template <typename CLASS, typename OutVector, typename InVector>
  void
  loop_cell_centric(void (CLASS::*cell_operation)(
                      const MatrixFree &,
                      OutVector &,
                      const InVector &,
                      const std::pair<unsigned int, unsigned int> &),
                    CLASS                  *owning_class,
                    OutVector              &dst,
                    const InVector         &src,
                    const bool              zero_dst_vector = false,
                    const DataAccessOnFaces src_vector_face_access =
                      DataAccessOnFaces::unspecified) const;

  /**
   * Same as above, but with std::function.
   */
  template <typename OutVector, typename InVector>
  void
  loop_cell_centric(
    const std::function<void(const MatrixFree &,
                             OutVector &,
                             const InVector &,
                             const std::pair<unsigned int, unsigned int> &)>
                           &cell_operation,
    OutVector              &dst,
    const InVector         &src,
    const bool              zero_dst_vector = false,
    const DataAccessOnFaces src_vector_face_access =
      DataAccessOnFaces::unspecified) const;

  /**
   * In the hp-adaptive case, a subrange of cells as computed during the cell
   * loop might contain elements of different degrees. Use this function to
   * compute what the subrange for an individual finite element degree is. The
   * finite element degree is associated to the vector component given in the
   * function call.
   */
  std::pair<unsigned int, unsigned int>
  create_cell_subrange_hp(const std::pair<unsigned int, unsigned int> &range,
                          const unsigned int fe_degree,
                          const unsigned int dof_handler_index = 0) const;

  /**
   * In the hp-adaptive case, a subrange of cells as computed during the cell
   * loop might contain elements of different degrees. Use this function to
   * compute what the subrange for a given index the hp-finite element, as
   * opposed to the finite element degree in the other function.
   */
  std::pair<unsigned int, unsigned int>
  create_cell_subrange_hp_by_index(
    const std::pair<unsigned int, unsigned int> &range,
    const unsigned int                           fe_index,
    const unsigned int                           dof_handler_index = 0) const;

  /**
   * In the hp-adaptive case, return number of active FE indices.
   */
  unsigned int
  n_active_fe_indices() const;

  /**
   * In the hp-adaptive case, return the active FE index of a cell range.
   */
  unsigned int
  get_cell_active_fe_index(
    const std::pair<unsigned int, unsigned int> range,
    const unsigned int dof_handler_index = numbers::invalid_unsigned_int) const;

  /**
   * In the hp-adaptive case, return the active FE index of a face range.
   */
  unsigned int
  get_face_active_fe_index(
    const std::pair<unsigned int, unsigned int> range,
    const bool                                  is_interior_face = true,
    const unsigned int dof_handler_index = numbers::invalid_unsigned_int) const;

  /** @} */

  /**
   * @name Initialization of vectors
   */
  /** @{ */
  /**
   * Initialize function for a vector with each entry associated with a cell
   * batch (cell data). For reading and writing the vector use:
   * FEEvaluationBase::read_cell_data() and FEEvaluationBase::set_cell_data().
   */
  template <typename T>
  void
  initialize_cell_data_vector(AlignedVector<T> &vec) const;

  /**
   * Initialize function for a vector with each entry associated with a face
   * batch (face data). For reading and writing the vector use:
   * FEEvaluationBase::read_face_data() and FEEvaluationBase::set_face_data().
   */
  template <typename T>
  void
  initialize_face_data_vector(AlignedVector<T> &vec) const;

  /**
   * Initialize function for a general serial non-block vector.
   * After a call to this function, the length of the vector is equal to the
   * total number of degrees of freedom in the DoFHandler. Vector entries are
   * initialized with zero.
   *
   * If MatrixFree was set up with several DoFHandler objects, the parameter
   * @p dof_handler_index defines which component is to be used.
   *
   * @note Serial vectors also include Trilinos and PETSc vectors; however
   * in these cases, MatrixFree has to be used in a serial context, i.e., the
   * size of the communicator has to be exactly one.
   */
  template <typename VectorType>
  void
  initialize_dof_vector(VectorType        &vec,
                        const unsigned int dof_handler_index = 0) const;

  /**
   * Specialization of the method initialize_dof_vector() for the
   * class LinearAlgebra::distributed::Vector@<Number@>.
   * See the other function with the same name for the general descriptions.
   *
   * @note For the parallel vectors used with MatrixFree and in FEEvaluation, a
   * vector needs to hold all
   * @ref GlossLocallyActiveDof "locally active DoFs"
   * and also some of the
   * @ref GlossLocallyRelevantDof "locally relevant DoFs".
   * The selection of DoFs is such that one can read all degrees of freedom on
   * all locally relevant elements (locally active) plus the degrees of freedom
   * that constraints expand into from the locally owned cells. However, not
   * all locally relevant DoFs are stored because most of them would never be
   * accessed in matrix-vector products and result in too much data sent
   * around which impacts the performance.
   */
  template <typename Number2>
  void
  initialize_dof_vector(LinearAlgebra::distributed::Vector<Number2> &vec,
                        const unsigned int dof_handler_index = 0) const;

  /**
   * Return the partitioner that represents the locally owned data and the
   * ghost indices where access is needed to for the cell loop. The
   * partitioner is constructed from the locally owned dofs and ghost dofs
   * given by the respective fields. If you want to have specific information
   * about these objects, you can query them with the respective access
   * functions. If you just want to initialize a (parallel) vector, you should
   * usually prefer this data structure as the data exchange information can
   * be reused from one vector to another.
   */
  const std::shared_ptr<const Utilities::MPI::Partitioner> &
  get_vector_partitioner(const unsigned int dof_handler_index = 0) const;

  /**
   * Return the set of cells that are owned by the processor.
   */
  const IndexSet &
  get_locally_owned_set(const unsigned int dof_handler_index = 0) const;

  /**
   * Return the set of ghost cells needed but not owned by the processor.
   */
  const IndexSet &
  get_ghost_set(const unsigned int dof_handler_index = 0) const;

  /**
   * Return a list of all degrees of freedom that are constrained. The list
   * is returned in MPI-local index space for the locally owned range of the
   * vector, not in global MPI index space that spans all MPI processors. To
   * get numbers in global index space, call
   * <tt>get_vector_partitioner()->local_to_global</tt> on an entry of the
   * vector. In addition, it only returns the indices for degrees of freedom
   * that are owned locally, not for ghosts.
   */
  const std::vector<unsigned int> &
  get_constrained_dofs(const unsigned int dof_handler_index = 0) const;

  /**
   * Computes a renumbering of degrees of freedom that better fits with the
   * data layout in MatrixFree according to the given layout of data. Note that
   * this function does not re-arrange the information stored in this class,
   * but rather creates a renumbering for consumption of
   * DoFHandler::renumber_dofs. To have any effect a MatrixFree object must be
   * set up again using the renumbered DoFHandler and AffineConstraints. Note
   * that if a DoFHandler calls DoFHandler::renumber_dofs, all information in
   * MatrixFree becomes invalid.
   */
  void
  renumber_dofs(std::vector<types::global_dof_index> &renumbering,
                const unsigned int                    dof_handler_index = 0);

  /** @} */

  /**
   * @name General information
   */
  /** @{ */
  /**
   * Return whether a given FiniteElement @p fe is supported by this class.
   */
  template <int spacedim>
  static bool
  is_supported(const FiniteElement<dim, spacedim> &fe);

  /**
   * Return the number of different DoFHandlers specified at initialization.
   */
  unsigned int
  n_components() const;

  /**
   * For the finite element underlying the DoFHandler specified by @p
   * dof_handler_index, return the number of base elements.
   */
  unsigned int
  n_base_elements(const unsigned int dof_handler_index) const;

  /**
   * Return the number of cells this structure is based on. If you are using a
   * usual DoFHandler, it corresponds to the number of (locally owned) active
   * cells. Note that most data structures in this class do not directly act
   * on this number but rather on n_cell_batches() which gives the number of
   * cells as seen when lumping several cells together with vectorization.
   */
  unsigned int
  n_physical_cells() const;

  /**
   * Return the number of cell batches that this structure works on. The
   * batches are formed by application of vectorization over several cells in
   * general. The cell range in @p cell_loop runs from zero to
   * n_cell_batches() (exclusive), so this is the appropriate size if you want
   * to store arrays of data for all cells to be worked on. This number is
   * approximately `n_physical_cells()/VectorizedArray::%size()`
   * (depending on how many cell batches that do not get filled up completely).
   */
  unsigned int
  n_cell_batches() const;

  /**
   * Return the number of additional cell batches that this structure keeps
   * for face integration. Note that not all cells that are ghosted in the
   * triangulation are kept in this data structure, but only the ones which
   * are necessary for evaluating face integrals from both sides.
   */
  unsigned int
  n_ghost_cell_batches() const;

  /**
   * Return the number of interior face batches that this structure works on.
   * The batches are formed by application of vectorization over several faces
   * in general. The face range in @p loop runs from zero to
   * n_inner_face_batches() (exclusive), so this is the appropriate size if
   * you want to store arrays of data for all interior faces to be worked on.
   * Note that it returns 0 unless mapping_update_flags_inner_faces is set
   * to a value different from  UpdateFlags::update_default.
   */
  unsigned int
  n_inner_face_batches() const;

  /**
   * Return the number of boundary face batches that this structure works on.
   * The batches are formed by application of vectorization over several faces
   * in general. The face range in @p loop runs from n_inner_face_batches() to
   * n_inner_face_batches()+n_boundary_face_batches() (exclusive), so if you
   * need to store arrays that hold data for all boundary faces but not the
   * interior ones, this number gives the appropriate size.
   * Note that it returns 0 unless mapping_update_flags_boundary_faces is set
   * to a value different from UpdateFlags::update_default.
   */
  unsigned int
  n_boundary_face_batches() const;

  /**
   * Return the number of faces that are not processed locally but belong to
   * locally owned faces.
   */
  unsigned int
  n_ghost_inner_face_batches() const;

  /**
   * In order to apply different operators to different parts of the boundary,
   * this method can be used to query the boundary id of a given face in the
   * faces' own sorting by lanes in a VectorizedArray. Only valid for an index
   * indicating a boundary face.
   *
   * @note Alternatively to this function, you can use
   * FEFaceEvaluation::boundary_id() to get the same information if a
   * FEFaceEvaluation object has been set up already.
   */
  types::boundary_id
  get_boundary_id(const unsigned int face_batch_index) const;

  /**
   * Return the boundary ids for the faces within a cell, using the cells'
   * sorting by lanes in the VectorizedArray.
   */
  std::array<types::boundary_id, VectorizedArrayType::size()>
  get_faces_by_cells_boundary_id(const unsigned int cell_batch_index,
                                 const unsigned int face_number) const;

  /**
   * Return the DoFHandler with the index as given to the respective
   * `std::vector` argument in the reinit() function.
   */
  const DoFHandler<dim> &
  get_dof_handler(const unsigned int dof_handler_index = 0) const;

  /**
   * Return the AffineConstraints with the index as given to the
   * respective `std::vector` argument in the reinit() function. Only available
   * if the AffineConstraints objects have the same template parameter Number as
   * MatrixFree. Throws an exception otherwise.
   */
  const AffineConstraints<Number> &
  get_affine_constraints(const unsigned int dof_handler_index = 0) const;

  /**
   * Return the cell iterator in deal.II speak to a given cell batch
   * (populating several lanes in a VectorizedArray) and the lane index within
   * the vectorization across cells in the renumbering of this structure.
   *
   * Note that the cell iterators in deal.II go through cells differently to
   * what the cell loop of this class does. This is because several cells are
   * processed together (vectorization across cells), and since cells with
   * neighbors on different MPI processors need to be accessed at a certain
   * time when accessing remote data and overlapping communication with
   * computation.
   */
  typename DoFHandler<dim>::cell_iterator
  get_cell_iterator(const unsigned int cell_batch_index,
                    const unsigned int lane_index,
                    const unsigned int dof_handler_index = 0) const;

  /**
   * This returns the level and index for the cell that would be returned by
   * get_cell_iterator() for the same arguments `cell_batch_index` and
   * `lane_index`.
   */
  std::pair<int, int>
  get_cell_level_and_index(const unsigned int cell_batch_index,
                           const unsigned int lane_index) const;

  /**
   * Get MatrixFree index associated to a deal.II @p cell. To get
   * the actual cell batch index and lane, do the postprocessing
   * `index / VectorizedArrayType::size()` and `index %
   * VectorizedArrayType::size()`.
   */
  unsigned int
  get_matrix_free_cell_index(
    const typename Triangulation<dim>::cell_iterator &cell) const;

  /**
   * Return the cell iterator in deal.II speak to an interior/exterior cell of
   * a face in a pair of a face batch and lane index. The second element of
   * the pair is the face number so that the face iterator can be accessed:
   * `pair.first->face(pair.second);`
   *
   * Note that the face iterators in deal.II go through cells differently to
   * what the face/boundary loop of this class does. This is because several
   * faces are worked on together (vectorization), and since faces with neighbor
   * cells on different MPI processors need to be accessed at a certain time
   * when accessing remote data and overlapping communication with computation.
   */
  std::pair<typename DoFHandler<dim>::cell_iterator, unsigned int>
  get_face_iterator(const unsigned int face_batch_index,
                    const unsigned int lane_index,
                    const bool         interior     = true,
                    const unsigned int fe_component = 0) const;

  /**
   * Since this class uses vectorized data types with usually more than one
   * value in the data field, a situation might occur when some components of
   * the vector type do not correspond to an actual cell in the mesh. When
   * using only this class, one usually does not need to bother about that
   * fact since the values are padded with zeros. However, when this class is
   * mixed with deal.II access to cells, care needs to be taken. This function
   * returns @p true if not all `n_lanes` cells for the given
   * `cell_batch_index` correspond to actual cells of the mesh and some are
   * merely present for padding reasons. To find out how many cells are
   * actually used, use the function n_active_entries_per_cell_batch().
   */
  bool
  at_irregular_cell(const unsigned int cell_batch_index) const;

  /**
   * This query returns how many cells among the `VectorizedArrayType::size()`
   * many cells within a cell batch to actual cells in the mesh, rather than
   * being present for padding reasons. For most given cell batches in
   * n_cell_batches(), this number is equal to `VectorizedArrayType::size()`,
   * but there might be one or a few cell batches in the mesh (where the
   * numbers do not add up) where only some of the cells within a batch are
   * used, indicated by the function at_irregular_cell().
   */
  unsigned int
  n_active_entries_per_cell_batch(const unsigned int cell_batch_index) const;

  /**
   * Use this function to find out how many faces over the length of
   * vectorization data types correspond to real faces (both interior and
   * boundary faces, as those use the same indexing but with different ranges)
   * in the mesh. For most given indices in n_inner_faces_batches() and
   * n_boundary_face_batches(), this is just @p vectorization_length many, but
   * there might be one or a few meshes (where the numbers do not add up)
   * where there are less such lanes filled.
   */
  unsigned int
  n_active_entries_per_face_batch(const unsigned int face_batch_index) const;

  /**
   * Return the number of degrees of freedom per cell for a given hp-index.
   */
  unsigned int
  get_dofs_per_cell(const unsigned int dof_handler_index  = 0,
                    const unsigned int hp_active_fe_index = 0) const;

  /**
   * Return the number of quadrature points per cell for a given hp-index.
   */
  unsigned int
  get_n_q_points(const unsigned int quad_index         = 0,
                 const unsigned int hp_active_fe_index = 0) const;

  /**
   * Return the number of degrees of freedom on each face of the cell for
   * given hp-index.
   */
  unsigned int
  get_dofs_per_face(const unsigned int dof_handler_index  = 0,
                    const unsigned int hp_active_fe_index = 0) const;

  /**
   * Return the number of quadrature points on each face of the cell for
   * given hp-index.
   */
  unsigned int
  get_n_q_points_face(const unsigned int quad_index         = 0,
                      const unsigned int hp_active_fe_index = 0) const;

  /**
   * Return the quadrature rule for given hp-index.
   */
  const Quadrature<dim> &
  get_quadrature(const unsigned int quad_index         = 0,
                 const unsigned int hp_active_fe_index = 0) const;

  /**
   * Return the quadrature rule for given hp-index.
   */
  const Quadrature<dim - 1> &
  get_face_quadrature(const unsigned int quad_index         = 0,
                      const unsigned int hp_active_fe_index = 0) const;

  /**
   * Return the category the current batch range of cells was assigned to.
   * Categories run between the given values in the field
   * AdditionalData::cell_vectorization_category for the non-hp case
   * and return the active FE index in the hp-adaptive case.
   *
   * @note Following the behaviour of get_cell_category(), we return the
   * maximum category of any cell batch. In the hp case, it is
   * guaranteed that all cells and as a consequence all cell batches in a range
   * have the same category. Otherwise, there may be different categories in
   * different cell batches.
   */
  unsigned int
  get_cell_range_category(
    const std::pair<unsigned int, unsigned int> cell_batch_range,
    const unsigned int dof_handler_index = numbers::invalid_unsigned_int) const;

  /**
   * Return the category of the cells on the two sides of the current batch
   * range of faces.
   */
  std::pair<unsigned int, unsigned int>
  get_face_range_category(
    const std::pair<unsigned int, unsigned int> face_batch_range,
    const unsigned int dof_handler_index = numbers::invalid_unsigned_int) const;

  /**
   * Return the category the current batch of cells was assigned to. Categories
   * run between the given values in the field
   * AdditionalData::cell_vectorization_category for the non-hp case
   * and return the active FE index in the hp-adaptive case.
   *
   * @note In the non-hp case, a category of a cell batch is given
   * as the maximum category of any of its cell. In the hp case or the case that
   * MatrixFree::AdditionalData::cell_vectorization_categories_strict was
   * enabled, it is guaranteed that all cells have the same category.
   */
  unsigned int
  get_cell_category(
    const unsigned int cell_batch_index,
    const unsigned int dof_handler_index = numbers::invalid_unsigned_int) const;

  /**
   * Return the category of the cells on the two sides of the current batch of
   * faces.
   */
  std::pair<unsigned int, unsigned int>
  get_face_category(
    const unsigned int face_batch_index,
    const unsigned int dof_handler_index = numbers::invalid_unsigned_int) const;

  /**
   * Queries whether or not the indexation has been set.
   */
  bool
  indices_initialized() const;

  /**
   * Queries whether or not the geometry-related information for the cells has
   * been set.
   */
  bool
  mapping_initialized() const;

  /**
   * Return the level of the mesh to be worked on. Returns
   * numbers::invalid_unsigned_int if working on active cells.
   */
  unsigned int
  get_mg_level() const;

  /**
   * Return an approximation of the memory consumption of this class in
   * bytes.
   */
  std::size_t
  memory_consumption() const;

  /**
   * Prints a detailed summary of memory consumption in the different
   * structures of this class to the given output stream.
   */
  template <typename StreamType>
  void
  print_memory_consumption(StreamType &out) const;

  /**
   * Prints a summary of this class to the given output stream. It is focused
   * on the indices, and does not print all the data stored.
   */
  void
  print(std::ostream &out) const;

  /** @} */

  /**
   * @name Access of internal data structure
   *
   * Note: Expert mode, interface not stable between releases.
   */
  /** @{ */
  /**
   * Return information on task graph.
   */
  const internal::MatrixFreeFunctions::TaskInfo &
  get_task_info() const;

  /*
   * Return geometry-dependent information on the cells.
   */
  const internal::MatrixFreeFunctions::
    MappingInfo<dim, Number, VectorizedArrayType> &
    get_mapping_info() const;

  /**
   * Return information on indexation degrees of freedom.
   */
  const internal::MatrixFreeFunctions::DoFInfo &
  get_dof_info(const unsigned int dof_handler_index_component = 0) const;

  /**
   * Return the number of weights in the constraint pool.
   */
  unsigned int
  n_constraint_pool_entries() const;

  /**
   * Return a pointer to the first number in the constraint pool data with
   * index @p pool_index (to be used together with @p constraint_pool_end()).
   */
  const Number *
  constraint_pool_begin(const unsigned int pool_index) const;

  /**
   * Return a pointer to one past the last number in the constraint pool data
   * with index @p pool_index (to be used together with @p
   * constraint_pool_begin()).
   */
  const Number *
  constraint_pool_end(const unsigned int pool_index) const;

  /**
   * Return the unit cell information for given hp-index.
   */
  const internal::MatrixFreeFunctions::ShapeInfo<Number> &
  get_shape_info(const unsigned int dof_handler_index_component = 0,
                 const unsigned int quad_index                  = 0,
                 const unsigned int fe_base_element             = 0,
                 const unsigned int hp_active_fe_index          = 0,
                 const unsigned int hp_active_quad_index        = 0) const;

  /**
   * Return the connectivity information of a face.
   */
  const internal::MatrixFreeFunctions::FaceToCellTopology<
    VectorizedArrayType::size()> &
  get_face_info(const unsigned int face_batch_index) const;


  /**
   * Return the table that translates a triple of the cell-batch number,
   * the index of a face within a cell and the index within the cell batch of
   * vectorization into the index within the faces array.
   */
  const Table<3, unsigned int> &
  get_cell_and_face_to_plain_faces() const;

  /**
   * Obtains a scratch data object for internal use. Make sure to release it
   * afterwards by passing the pointer you obtain from this object to the
   * release_scratch_data() function. This interface is used by FEEvaluation
   * objects for storing their data structures.
   *
   * The organization of the internal data structure is a thread-local storage
   * of a list of vectors. Multiple threads will each get a separate storage
   * field and separate vectors, ensuring thread safety. The mechanism to
   * acquire and release objects is similar to the mechanisms used for the
   * local contributions of WorkStream, see
   * @ref workstream_paper "the WorkStream paper".
   */
  AlignedVector<VectorizedArrayType> *
  acquire_scratch_data() const;

  /**
   * Makes the object of the scratchpad available again.
   */
  void
  release_scratch_data(const AlignedVector<VectorizedArrayType> *memory) const;

  /**
   * Obtains a scratch data object for internal use. Make sure to release it
   * afterwards by passing the pointer you obtain from this object to the
   * release_scratch_data_non_threadsafe() function. Note that, as opposed to
   * acquire_scratch_data(), this method can only be called by a single thread
   * at a time, but opposed to the acquire_scratch_data() it is also possible
   * that the thread releasing the scratch data can be different than the one
   * that acquired it.
   */
  AlignedVector<Number> *
  acquire_scratch_data_non_threadsafe() const;

  /**
   * Makes the object of the scratch data available again.
   */
  void
  release_scratch_data_non_threadsafe(
    const AlignedVector<Number> *memory) const;

  /** @} */

private:
  /**
   * This is the actual reinit function that sets up the indices for the
   * DoFHandler case.
   */
  template <typename number2, int q_dim>
  void
  internal_reinit(
    const std::shared_ptr<hp::MappingCollection<dim>>     &mapping,
    const std::vector<const DoFHandler<dim, dim> *>       &dof_handlers,
    const std::vector<const AffineConstraints<number2> *> &constraint,
    const std::vector<IndexSet>                           &locally_owned_set,
    const std::vector<hp::QCollection<q_dim>>             &quad,
    const AdditionalData                                  &additional_data);

  /**
   * Initializes the fields in DoFInfo together with the constraint pool that
   * holds all different weights in the constraints (not part of DoFInfo
   * because several DoFInfo classes can have the same weights which
   * consequently only need to be stored once).
   */
  template <typename number2>
  void
  initialize_indices(
    const std::vector<const AffineConstraints<number2> *> &constraint,
    const std::vector<IndexSet>                           &locally_owned_set,
    const AdditionalData                                  &additional_data);

  /**
   * Initializes the DoFHandlers based on a DoFHandler<dim> argument.
   */
  void
  initialize_dof_handlers(
    const std::vector<const DoFHandler<dim, dim> *> &dof_handlers,
    const AdditionalData                            &additional_data);

  /**
   * Pointers to the DoFHandlers underlying the current problem.
   */
  std::vector<ObserverPointer<const DoFHandler<dim>>> dof_handlers;

  /**
   * Pointers to the AffineConstraints underlying the current problem. Only
   * filled with an AffineConstraints object if objects of the same `Number`
   * template parameter as the `Number` template of MatrixFree is passed to
   * reinit(). Filled with nullptr otherwise.
   */
  std::vector<ObserverPointer<const AffineConstraints<Number>>>
    affine_constraints;

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
  internal::MatrixFreeFunctions::MappingInfo<dim, Number, VectorizedArrayType>
    mapping_info;

  /**
   * Contains shape value information on the unit cell.
   */
  Table<4, internal::MatrixFreeFunctions::ShapeInfo<Number>> shape_info;

  /**
   * Describes how the cells are gone through. With the cell level (first
   * index in this field) and the index within the level, one can reconstruct
   * a deal.II cell iterator and use all the traditional things deal.II offers
   * to do with cell iterators.
   */
  std::vector<std::pair<unsigned int, unsigned int>> cell_level_index;

  /**
   * Conversion from deal.II index (active or level index) to MatrixFree index
   * (inverse of cell_level_index).
   */
  std::vector<unsigned int> mf_cell_indices;

  /**
   * For discontinuous Galerkin, the cell_level_index includes cells that are
   * not on the local processor but that are needed to evaluate the cell
   * integrals. In cell_level_index_end_local, we store the number of local
   * cells.
   */
  unsigned int cell_level_index_end_local;

  /**
   * Stores the basic layout of the cells and faces to be treated, including
   * the task layout for the shared memory parallelization and possible
   * overlaps between communications and computations with MPI.
   */
  internal::MatrixFreeFunctions::TaskInfo task_info;

  /**
   * Vector holding face information. Only initialized if
   * build_face_info=true.
   */
  internal::MatrixFreeFunctions::FaceInfo<VectorizedArrayType::size()>
    face_info;

  /**
   * Stores whether indices have been initialized.
   */
  bool indices_are_initialized;

  /**
   * Stores whether indices have been initialized.
   */
  bool mapping_is_initialized;

  /**
   * Scratchpad memory for use in evaluation. We allow more than one
   * evaluation object to attach to this field (this, the outer
   * std::vector), so we need to keep tracked of whether a certain data
   * field is already used (first part of pair) and keep a list of
   * objects.
   */
  mutable Threads::ThreadLocalStorage<
    std::list<std::pair<bool, AlignedVector<VectorizedArrayType>>>>
    scratch_pad;

  /**
   * Scratchpad memory for use in evaluation and other contexts, non-thread
   * safe variant.
   */
  mutable std::list<std::pair<bool, AlignedVector<Number>>>
    scratch_pad_non_threadsafe;

  /**
   * Stored the level of the mesh to be worked on.
   */
  unsigned int mg_level;

  /**
   * Stores the index of the first DoFHandler that is in hp-mode. If no
   * DoFHandler is in hp-mode, the value is 0.
   */
  unsigned int first_hp_dof_handler_index;
};



/*----------------------- Inline functions ----------------------------------*/

#ifndef DOXYGEN



template <int dim, typename Number, typename VectorizedArrayType>
template <typename T>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::initialize_cell_data_vector(
  AlignedVector<T> &vec) const
{
  vec.resize(this->n_cell_batches() + this->n_ghost_cell_batches());
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename T>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::initialize_face_data_vector(
  AlignedVector<T> &vec) const
{
  vec.resize(this->n_inner_face_batches() + this->n_boundary_face_batches() +
             this->n_ghost_inner_face_batches());
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename VectorType>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::initialize_dof_vector(
  VectorType        &vec,
  const unsigned int comp) const
{
  static_assert(IsBlockVector<VectorType>::value == false,
                "This function is not supported for block vectors.");

  Assert(task_info.n_procs == 1,
         ExcMessage("This function can only be used in serial."));

  AssertIndexRange(comp, n_components());
  vec.reinit(dof_info[comp].vector_partitioner->size());
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename Number2>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::initialize_dof_vector(
  LinearAlgebra::distributed::Vector<Number2> &vec,
  const unsigned int                           comp) const
{
  AssertIndexRange(comp, n_components());
  vec.reinit(dof_info[comp].vector_partitioner, task_info.communicator_sm);
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const std::shared_ptr<const Utilities::MPI::Partitioner> &
MatrixFree<dim, Number, VectorizedArrayType>::get_vector_partitioner(
  const unsigned int comp) const
{
  AssertIndexRange(comp, n_components());
  return dof_info[comp].vector_partitioner;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const std::vector<unsigned int> &
MatrixFree<dim, Number, VectorizedArrayType>::get_constrained_dofs(
  const unsigned int comp) const
{
  AssertIndexRange(comp, n_components());
  return dof_info[comp].constrained_dofs;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_components() const
{
  AssertDimension(dof_handlers.size(), dof_info.size());
  return dof_handlers.size();
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_base_elements(
  const unsigned int dof_no) const
{
  AssertDimension(dof_handlers.size(), dof_info.size());
  AssertIndexRange(dof_no, dof_handlers.size());
  return dof_handlers[dof_no]->get_fe().n_base_elements();
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const internal::MatrixFreeFunctions::TaskInfo &
MatrixFree<dim, Number, VectorizedArrayType>::get_task_info() const
{
  return task_info;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_physical_cells() const
{
  return task_info.n_active_cells;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_cell_batches() const
{
  return *(task_info.cell_partition_data.end() - 2);
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_ghost_cell_batches() const
{
  return *(task_info.cell_partition_data.end() - 1) -
         *(task_info.cell_partition_data.end() - 2);
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_inner_face_batches() const
{
  if (task_info.face_partition_data.empty())
    return 0;
  return task_info.face_partition_data.back();
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_boundary_face_batches() const
{
  if (task_info.face_partition_data.empty())
    return 0;
  return task_info.boundary_partition_data.back() -
         task_info.face_partition_data.back();
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_ghost_inner_face_batches() const
{
  if (task_info.face_partition_data.empty())
    return 0;
  return face_info.faces.size() - task_info.boundary_partition_data.back();
}



template <int dim, typename Number, typename VectorizedArrayType>
inline types::boundary_id
MatrixFree<dim, Number, VectorizedArrayType>::get_boundary_id(
  const unsigned int face_batch_index) const
{
  Assert(face_batch_index >= task_info.boundary_partition_data[0] &&
           face_batch_index < task_info.boundary_partition_data.back(),
         ExcIndexRange(face_batch_index,
                       task_info.boundary_partition_data[0],
                       task_info.boundary_partition_data.back()));
  return types::boundary_id(face_info.faces[face_batch_index].exterior_face_no);
}



template <int dim, typename Number, typename VectorizedArrayType>
inline std::array<types::boundary_id, VectorizedArrayType::size()>
MatrixFree<dim, Number, VectorizedArrayType>::get_faces_by_cells_boundary_id(
  const unsigned int cell_batch_index,
  const unsigned int face_number) const
{
  AssertIndexRange(cell_batch_index, n_cell_batches());
  AssertIndexRange(face_number, ReferenceCells::max_n_faces<dim>());
  Assert(face_info.cell_and_face_boundary_id.size(0) >= n_cell_batches(),
         ExcNotInitialized());
  std::array<types::boundary_id, VectorizedArrayType::size()> result;
  result.fill(numbers::invalid_boundary_id);
  for (unsigned int v = 0;
       v < n_active_entries_per_cell_batch(cell_batch_index);
       ++v)
    result[v] =
      face_info.cell_and_face_boundary_id(cell_batch_index, face_number, v);
  return result;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const internal::MatrixFreeFunctions::
  MappingInfo<dim, Number, VectorizedArrayType> &
  MatrixFree<dim, Number, VectorizedArrayType>::get_mapping_info() const
{
  return mapping_info;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const internal::MatrixFreeFunctions::DoFInfo &
MatrixFree<dim, Number, VectorizedArrayType>::get_dof_info(
  const unsigned int dof_index) const
{
  AssertIndexRange(dof_index, n_components());
  return dof_info[dof_index];
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_constraint_pool_entries() const
{
  return constraint_pool_row_index.size() - 1;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const Number *
MatrixFree<dim, Number, VectorizedArrayType>::constraint_pool_begin(
  const unsigned int row) const
{
  AssertIndexRange(row, constraint_pool_row_index.size() - 1);
  return constraint_pool_data.empty() ?
           nullptr :
           constraint_pool_data.data() + constraint_pool_row_index[row];
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const Number *
MatrixFree<dim, Number, VectorizedArrayType>::constraint_pool_end(
  const unsigned int row) const
{
  AssertIndexRange(row, constraint_pool_row_index.size() - 1);
  return constraint_pool_data.empty() ?
           nullptr :
           constraint_pool_data.data() + constraint_pool_row_index[row + 1];
}



template <int dim, typename Number, typename VectorizedArrayType>
inline std::pair<unsigned int, unsigned int>
MatrixFree<dim, Number, VectorizedArrayType>::create_cell_subrange_hp(
  const std::pair<unsigned int, unsigned int> &range,
  const unsigned int                           degree,
  const unsigned int                           dof_handler_component) const
{
  if (dof_info[dof_handler_component].cell_active_fe_index.empty())
    {
      AssertDimension(
        dof_info[dof_handler_component].fe_index_conversion.size(), 1);
      AssertDimension(
        dof_info[dof_handler_component].fe_index_conversion[0].size(), 1);
      if (dof_info[dof_handler_component].fe_index_conversion[0][0] == degree)
        return range;
      else
        return {range.second, range.second};
    }

  const unsigned int fe_index =
    dof_info[dof_handler_component].fe_index_from_degree(0, degree);
  if (fe_index >= dof_info[dof_handler_component].max_fe_index)
    return {range.second, range.second};
  else
    return create_cell_subrange_hp_by_index(range,
                                            fe_index,
                                            dof_handler_component);
}



template <int dim, typename Number, typename VectorizedArrayType>
inline bool
MatrixFree<dim, Number, VectorizedArrayType>::at_irregular_cell(
  const unsigned int cell_batch_index) const
{
  AssertIndexRange(cell_batch_index, task_info.cell_partition_data.back());
  return VectorizedArrayType::size() > 1 &&
         cell_level_index[(cell_batch_index + 1) * VectorizedArrayType::size() -
                          1] == cell_level_index[(cell_batch_index + 1) *
                                                   VectorizedArrayType::size() -
                                                 2];
}



template <int dim, typename Number, typename VectorizedArrayType>
unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_active_fe_indices() const
{
  return shape_info.size(2);
}


template <int dim, typename Number, typename VectorizedArrayType>
unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_cell_active_fe_index(
  const std::pair<unsigned int, unsigned int> range,
  const unsigned int                          dof_handler_index) const
{
  const unsigned int dof_no =
    dof_handler_index == numbers::invalid_unsigned_int ?
      first_hp_dof_handler_index :
      dof_handler_index;

  const auto &fe_indices = dof_info[dof_no].cell_active_fe_index;

  if (fe_indices.empty() == true ||
      dof_handlers[dof_no]->get_fe_collection().size() == 1)
    return 0;

  const auto index = fe_indices[range.first];

  for (unsigned int i = range.first; i < range.second; ++i)
    AssertDimension(index, fe_indices[i]);

  return index;
}



template <int dim, typename Number, typename VectorizedArrayType>
unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_face_active_fe_index(
  const std::pair<unsigned int, unsigned int> range,
  const bool                                  is_interior_face,
  const unsigned int                          dof_handler_index) const
{
  const unsigned int dof_no =
    dof_handler_index == numbers::invalid_unsigned_int ?
      first_hp_dof_handler_index :
      dof_handler_index;

  const auto &fe_indices = dof_info[dof_no].cell_active_fe_index;

  if (fe_indices.empty() == true)
    return 0;

  if (is_interior_face)
    {
      const unsigned int index =
        fe_indices[face_info.faces[range.first].cells_interior[0] /
                   VectorizedArrayType::size()];

      for (unsigned int i = range.first; i < range.second; ++i)
        AssertDimension(index,
                        fe_indices[face_info.faces[i].cells_interior[0] /
                                   VectorizedArrayType::size()]);

      return index;
    }
  else
    {
      const unsigned int index =
        fe_indices[face_info.faces[range.first].cells_exterior[0] /
                   VectorizedArrayType::size()];

      for (unsigned int i = range.first; i < range.second; ++i)
        AssertDimension(index,
                        fe_indices[face_info.faces[i].cells_exterior[0] /
                                   VectorizedArrayType::size()]);

      return index;
    }
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_active_entries_per_cell_batch(
  const unsigned int cell_batch_index) const
{
  Assert(!dof_info.empty(), ExcNotInitialized());
  AssertIndexRange(cell_batch_index, task_info.cell_partition_data.back());
  const std::vector<unsigned char> &n_lanes_filled =
    dof_info[0].n_vectorization_lanes_filled
      [internal::MatrixFreeFunctions::DoFInfo::dof_access_cell];
  AssertIndexRange(cell_batch_index, n_lanes_filled.size());

  return n_lanes_filled[cell_batch_index];
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::n_active_entries_per_face_batch(
  const unsigned int face_batch_index) const
{
  AssertIndexRange(face_batch_index, face_info.faces.size());
  Assert(!dof_info.empty(), ExcNotInitialized());
  const std::vector<unsigned char> &n_lanes_filled =
    dof_info[0].n_vectorization_lanes_filled
      [internal::MatrixFreeFunctions::DoFInfo::dof_access_face_interior];
  AssertIndexRange(face_batch_index, n_lanes_filled.size());
  return n_lanes_filled[face_batch_index];
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_dofs_per_cell(
  const unsigned int dof_handler_index,
  const unsigned int active_fe_index) const
{
  return dof_info[dof_handler_index].dofs_per_cell[active_fe_index];
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_n_q_points(
  const unsigned int quad_index,
  const unsigned int active_fe_index) const
{
  AssertIndexRange(quad_index, mapping_info.cell_data.size());
  return mapping_info.cell_data[quad_index]
    .descriptor[active_fe_index]
    .n_q_points;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_dofs_per_face(
  const unsigned int dof_handler_index,
  const unsigned int active_fe_index) const
{
  return dof_info[dof_handler_index].dofs_per_face[active_fe_index];
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_n_q_points_face(
  const unsigned int quad_index,
  const unsigned int active_fe_index) const
{
  AssertIndexRange(quad_index, mapping_info.face_data.size());
  return mapping_info.face_data[quad_index]
    .descriptor[active_fe_index]
    .n_q_points;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const IndexSet &
MatrixFree<dim, Number, VectorizedArrayType>::get_locally_owned_set(
  const unsigned int dof_handler_index) const
{
  return dof_info[dof_handler_index].vector_partitioner->locally_owned_range();
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const IndexSet &
MatrixFree<dim, Number, VectorizedArrayType>::get_ghost_set(
  const unsigned int dof_handler_index) const
{
  return dof_info[dof_handler_index].vector_partitioner->ghost_indices();
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const internal::MatrixFreeFunctions::ShapeInfo<Number> &
MatrixFree<dim, Number, VectorizedArrayType>::get_shape_info(
  const unsigned int dof_handler_index,
  const unsigned int index_quad,
  const unsigned int index_fe,
  const unsigned int active_fe_index,
  const unsigned int active_quad_index) const
{
  AssertIndexRange(dof_handler_index, dof_info.size());
  const unsigned int ind =
    dof_info[dof_handler_index].global_base_element_offset + index_fe;
  AssertIndexRange(ind, shape_info.size(0));
  AssertIndexRange(index_quad, shape_info.size(1));
  AssertIndexRange(active_fe_index, shape_info.size(2));
  AssertIndexRange(active_quad_index, shape_info.size(3));
  return shape_info(ind, index_quad, active_fe_index, active_quad_index);
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const internal::MatrixFreeFunctions::FaceToCellTopology<
  VectorizedArrayType::size()> &
MatrixFree<dim, Number, VectorizedArrayType>::get_face_info(
  const unsigned int face_batch_index) const
{
  AssertIndexRange(face_batch_index, face_info.faces.size());
  return face_info.faces[face_batch_index];
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const Table<3, unsigned int> &
MatrixFree<dim, Number, VectorizedArrayType>::get_cell_and_face_to_plain_faces()
  const
{
  return face_info.cell_and_face_to_plain_faces;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const Quadrature<dim> &
MatrixFree<dim, Number, VectorizedArrayType>::get_quadrature(
  const unsigned int quad_index,
  const unsigned int active_fe_index) const
{
  AssertIndexRange(quad_index, mapping_info.cell_data.size());
  return mapping_info.cell_data[quad_index]
    .descriptor[active_fe_index]
    .quadrature;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline const Quadrature<dim - 1> &
MatrixFree<dim, Number, VectorizedArrayType>::get_face_quadrature(
  const unsigned int quad_index,
  const unsigned int active_fe_index) const
{
  AssertIndexRange(quad_index, mapping_info.face_data.size());
  return mapping_info.face_data[quad_index]
    .descriptor[active_fe_index]
    .quadrature;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_cell_range_category(
  const std::pair<unsigned int, unsigned int> range,
  const unsigned int                          dof_handler_index) const
{
  auto result = get_cell_category(range.first, dof_handler_index);

  for (unsigned int i = range.first; i < range.second; ++i)
    result = std::max(result, get_cell_category(i, dof_handler_index));

  return result;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline std::pair<unsigned int, unsigned int>
MatrixFree<dim, Number, VectorizedArrayType>::get_face_range_category(
  const std::pair<unsigned int, unsigned int> range,
  const unsigned int                          dof_handler_index) const
{
  auto result = get_face_category(range.first, dof_handler_index);

  for (unsigned int i = range.first; i < range.second; ++i)
    {
      result.first =
        std::max(result.first, get_face_category(i, dof_handler_index).first);
      result.second =
        std::max(result.second, get_face_category(i, dof_handler_index).second);
    }

  return result;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_cell_category(
  const unsigned int cell_batch_index,
  const unsigned int dof_handler_index) const
{
  AssertIndexRange(0, dof_info.size());

  const unsigned int dof_no =
    dof_handler_index == numbers::invalid_unsigned_int ?
      first_hp_dof_handler_index :
      dof_handler_index;

  AssertIndexRange(cell_batch_index,
                   dof_info[dof_no].cell_active_fe_index.size());
  if (dof_info[dof_no].cell_active_fe_index.empty())
    return 0;
  else
    return dof_info[dof_no].cell_active_fe_index[cell_batch_index];
}



template <int dim, typename Number, typename VectorizedArrayType>
inline std::pair<unsigned int, unsigned int>
MatrixFree<dim, Number, VectorizedArrayType>::get_face_category(
  const unsigned int face_batch_index,
  const unsigned int dof_handler_index) const
{
  const unsigned int dof_no =
    dof_handler_index == numbers::invalid_unsigned_int ?
      first_hp_dof_handler_index :
      dof_handler_index;

  AssertIndexRange(face_batch_index, face_info.faces.size());
  if (dof_info[dof_no].cell_active_fe_index.empty())
    return std::make_pair(0U, 0U);

  std::pair<unsigned int, unsigned int> result = std::make_pair(0U, 0U);
  for (unsigned int v = 0;
       v < VectorizedArrayType::size() &&
       face_info.faces[face_batch_index].cells_interior[v] !=
         numbers::invalid_unsigned_int;
       ++v)
    result.first = std::max(
      result.first,
      dof_info[dof_no].cell_active_fe_index[face_info.faces[face_batch_index]
                                              .cells_interior[v] /
                                            VectorizedArrayType::size()]);
  if (face_info.faces[face_batch_index].cells_exterior[0] !=
      numbers::invalid_unsigned_int)
    for (unsigned int v = 0;
         v < VectorizedArrayType::size() &&
         face_info.faces[face_batch_index].cells_exterior[v] !=
           numbers::invalid_unsigned_int;
         ++v)
      result.second = std::max(
        result.second,
        dof_info[dof_no].cell_active_fe_index[face_info.faces[face_batch_index]
                                                .cells_exterior[v] /
                                              VectorizedArrayType::size()]);
  else
    result.second = numbers::invalid_unsigned_int;
  return result;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline bool
MatrixFree<dim, Number, VectorizedArrayType>::indices_initialized() const
{
  return indices_are_initialized;
}



template <int dim, typename Number, typename VectorizedArrayType>
inline bool
MatrixFree<dim, Number, VectorizedArrayType>::mapping_initialized() const
{
  return mapping_is_initialized;
}


template <int dim, typename Number, typename VectorizedArrayType>
inline unsigned int
MatrixFree<dim, Number, VectorizedArrayType>::get_mg_level() const
{
  return mg_level;
}



template <int dim, typename Number, typename VectorizedArrayType>
AlignedVector<VectorizedArrayType> *
MatrixFree<dim, Number, VectorizedArrayType>::acquire_scratch_data() const
{
  using list_type =
    std::list<std::pair<bool, AlignedVector<VectorizedArrayType>>>;
  list_type &data = scratch_pad.get();
  for (typename list_type::iterator it = data.begin(); it != data.end(); ++it)
    if (it->first == false)
      {
        it->first = true;
        return &it->second;
      }
  data.emplace_front(true, AlignedVector<VectorizedArrayType>());
  return &data.front().second;
}



template <int dim, typename Number, typename VectorizedArrayType>
void
MatrixFree<dim, Number, VectorizedArrayType>::release_scratch_data(
  const AlignedVector<VectorizedArrayType> *scratch) const
{
  using list_type =
    std::list<std::pair<bool, AlignedVector<VectorizedArrayType>>>;
  list_type &data = scratch_pad.get();
  for (typename list_type::iterator it = data.begin(); it != data.end(); ++it)
    if (&it->second == scratch)
      {
        Assert(it->first == true, ExcInternalError());
        it->first = false;
        return;
      }
  AssertThrow(false, ExcMessage("Tried to release invalid scratch pad"));
}



template <int dim, typename Number, typename VectorizedArrayType>
AlignedVector<Number> *
MatrixFree<dim, Number, VectorizedArrayType>::
  acquire_scratch_data_non_threadsafe() const
{
  for (typename std::list<std::pair<bool, AlignedVector<Number>>>::iterator it =
         scratch_pad_non_threadsafe.begin();
       it != scratch_pad_non_threadsafe.end();
       ++it)
    if (it->first == false)
      {
        it->first = true;
        return &it->second;
      }
  scratch_pad_non_threadsafe.push_front(
    std::make_pair(true, AlignedVector<Number>()));
  return &scratch_pad_non_threadsafe.front().second;
}



template <int dim, typename Number, typename VectorizedArrayType>
void
MatrixFree<dim, Number, VectorizedArrayType>::
  release_scratch_data_non_threadsafe(
    const AlignedVector<Number> *scratch) const
{
  for (typename std::list<std::pair<bool, AlignedVector<Number>>>::iterator it =
         scratch_pad_non_threadsafe.begin();
       it != scratch_pad_non_threadsafe.end();
       ++it)
    if (&it->second == scratch)
      {
        Assert(it->first == true, ExcInternalError());
        it->first = false;
        return;
      }
  AssertThrow(false, ExcMessage("Tried to release invalid scratch pad"));
}



// ------------------------------ reinit functions ---------------------------

namespace internal
{
  namespace MatrixFreeImplementation
  {
    template <int dim, int spacedim>
    inline std::vector<IndexSet>
    extract_locally_owned_index_sets(
      const std::vector<const DoFHandler<dim, spacedim> *> &dofh,
      const unsigned int                                    level)
    {
      std::vector<IndexSet> locally_owned_set;
      locally_owned_set.reserve(dofh.size());
      for (unsigned int j = 0; j < dofh.size(); ++j)
        if (level == numbers::invalid_unsigned_int)
          locally_owned_set.push_back(dofh[j]->locally_owned_dofs());
        else
          locally_owned_set.push_back(dofh[j]->locally_owned_mg_dofs(level));
      return locally_owned_set;
    }
  } // namespace MatrixFreeImplementation
} // namespace internal



template <int dim, typename Number, typename VectorizedArrayType>
template <typename QuadratureType, typename number2, typename MappingType>
void
MatrixFree<dim, Number, VectorizedArrayType>::reinit(
  const MappingType                &mapping,
  const DoFHandler<dim>            &dof_handler,
  const AffineConstraints<number2> &constraints_in,
  const QuadratureType             &quad,
  const typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    &additional_data)
{
  std::vector<const DoFHandler<dim, dim> *>       dof_handlers;
  std::vector<const AffineConstraints<number2> *> constraints;

  dof_handlers.push_back(&dof_handler);
  constraints.push_back(&constraints_in);

  std::vector<IndexSet> locally_owned_sets =
    internal::MatrixFreeImplementation::extract_locally_owned_index_sets(
      dof_handlers, additional_data.mg_level);

  std::vector<hp::QCollection<dim>> quad_hp;
  quad_hp.emplace_back(quad);

  internal_reinit(std::make_shared<hp::MappingCollection<dim>>(mapping),
                  dof_handlers,
                  constraints,
                  locally_owned_sets,
                  quad_hp,
                  additional_data);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename QuadratureType, typename number2, typename MappingType>
void
MatrixFree<dim, Number, VectorizedArrayType>::reinit(
  const MappingType                                     &mapping,
  const std::vector<const DoFHandler<dim> *>            &dof_handler,
  const std::vector<const AffineConstraints<number2> *> &constraint,
  const QuadratureType                                  &quad,
  const typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    &additional_data)
{
  std::vector<IndexSet> locally_owned_set =
    internal::MatrixFreeImplementation::extract_locally_owned_index_sets(
      dof_handler, additional_data.mg_level);
  std::vector<hp::QCollection<dim>> quad_hp;
  quad_hp.emplace_back(quad);

  internal_reinit(std::make_shared<hp::MappingCollection<dim>>(mapping),
                  dof_handler,
                  constraint,
                  locally_owned_set,
                  quad_hp,
                  additional_data);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename QuadratureType, typename number2, typename MappingType>
void
MatrixFree<dim, Number, VectorizedArrayType>::reinit(
  const MappingType                                     &mapping,
  const std::vector<const DoFHandler<dim> *>            &dof_handler,
  const std::vector<const AffineConstraints<number2> *> &constraint,
  const std::vector<QuadratureType>                     &quad,
  const typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    &additional_data)
{
  std::vector<IndexSet> locally_owned_set =
    internal::MatrixFreeImplementation::extract_locally_owned_index_sets(
      dof_handler, additional_data.mg_level);
  std::vector<hp::QCollection<dim>> quad_hp;
  for (unsigned int q = 0; q < quad.size(); ++q)
    quad_hp.emplace_back(quad[q]);

  internal_reinit(std::make_shared<hp::MappingCollection<dim>>(mapping),
                  dof_handler,
                  constraint,
                  locally_owned_set,
                  quad_hp,
                  additional_data);
}



// ------------------------------ implementation of loops --------------------

// internal helper functions that define how to call MPI data exchange
// functions: for generic vectors, do nothing at all. For distributed vectors,
// call update_ghost_values_start function and so on. If we have collections
// of vectors, just do the individual functions of the components. In order to
// keep ghost values consistent (whether we are in read or write mode), we
// also reset the values at the end. the whole situation is a bit complicated
// by the fact that we need to treat block vectors differently, which use some
// additional helper functions to select the blocks and template magic.
namespace internal
{
  /**
   * Internal class for exchanging data between vectors.
   */
  template <int dim, typename Number, typename VectorizedArrayType>
  struct VectorDataExchange
  {
    // A shift for the MPI messages to reduce the risk for accidental
    // interaction with other open communications that a user program might
    // set up (parallel vectors support unfinished communication). We let
    // the other vectors use the first 20 assigned numbers and start the
    // matrix-free communication.
    static constexpr unsigned int channel_shift = 20;



    /**
     * Constructor. Takes MF data, flag for face access in DG and
     * number of components.
     */
    VectorDataExchange(
      const dealii::MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
      const typename dealii::MatrixFree<dim, Number, VectorizedArrayType>::
        DataAccessOnFaces vector_face_access,
      const unsigned int  n_components)
      : matrix_free(matrix_free)
      , vector_face_access(
          matrix_free.get_task_info().face_partition_data.empty() ?
            dealii::MatrixFree<dim, Number, VectorizedArrayType>::
              DataAccessOnFaces::unspecified :
            vector_face_access)
      , ghosts_were_set(false)
#  ifdef DEAL_II_WITH_MPI
      , tmp_data(n_components)
      , requests(n_components)
#  endif
    {
      (void)n_components;
      if (this->vector_face_access !=
          dealii::MatrixFree<dim, Number, VectorizedArrayType>::
            DataAccessOnFaces::unspecified)
        for (unsigned int c = 0; c < matrix_free.n_components(); ++c)
          AssertDimension(
            matrix_free.get_dof_info(c).vector_exchanger_face_variants.size(),
            5);
    }



    /**
     * Destructor.
     */
    ~VectorDataExchange() // NOLINT
    {
#  ifdef DEAL_II_WITH_MPI
      for (unsigned int i = 0; i < tmp_data.size(); ++i)
        if (tmp_data[i] != nullptr)
          matrix_free.release_scratch_data_non_threadsafe(tmp_data[i]);
#  endif
    }



    /**
     * Go through all components in MF object and choose the one
     * whose partitioner is compatible with the Partitioner in this component.
     */
    template <typename VectorType>
    unsigned int
    find_vector_in_mf(const VectorType &vec,
                      const bool        check_global_compatibility = true) const
    {
      // case 1: vector was set up with MatrixFree::initialize_dof_vector()
      for (unsigned int c = 0; c < matrix_free.n_components(); ++c)
        if (vec.get_partitioner().get() ==
            matrix_free.get_dof_info(c).vector_partitioner.get())
          return c;

      // case 2: user provided own partitioner (compatibility mode)
      for (unsigned int c = 0; c < matrix_free.n_components(); ++c)
        if (check_global_compatibility ?
              vec.get_partitioner()->is_globally_compatible(
                *matrix_free.get_dof_info(c).vector_partitioner) :
              vec.get_partitioner()->is_compatible(
                *matrix_free.get_dof_info(c).vector_partitioner))
          return c;

      Assert(false,
             ExcNotImplemented("Could not find partitioner that fits vector"));

      return numbers::invalid_unsigned_int;
    }



    /**
     * Get partitioner for the given @p mf_component taking into
     * account vector_face_access set in constructor.
     */
    const internal::MatrixFreeFunctions::VectorDataExchange::Base &
    get_partitioner(const unsigned int mf_component) const
    {
      AssertDimension(matrix_free.get_dof_info(mf_component)
                        .vector_exchanger_face_variants.size(),
                      5);
      if (vector_face_access ==
          dealii::MatrixFree<dim, Number, VectorizedArrayType>::
            DataAccessOnFaces::none)
        return *matrix_free.get_dof_info(mf_component)
                  .vector_exchanger_face_variants[0];
      else if (vector_face_access ==
               dealii::MatrixFree<dim, Number, VectorizedArrayType>::
                 DataAccessOnFaces::values)
        return *matrix_free.get_dof_info(mf_component)
                  .vector_exchanger_face_variants[1];
      else if (vector_face_access ==
               dealii::MatrixFree<dim, Number, VectorizedArrayType>::
                 DataAccessOnFaces::gradients)
        return *matrix_free.get_dof_info(mf_component)
                  .vector_exchanger_face_variants[2];
      else if (vector_face_access ==
               dealii::MatrixFree<dim, Number, VectorizedArrayType>::
                 DataAccessOnFaces::values_all_faces)
        return *matrix_free.get_dof_info(mf_component)
                  .vector_exchanger_face_variants[3];
      else if (vector_face_access ==
               dealii::MatrixFree<dim, Number, VectorizedArrayType>::
                 DataAccessOnFaces::gradients_all_faces)
        return *matrix_free.get_dof_info(mf_component)
                  .vector_exchanger_face_variants[4];
      else
        return *matrix_free.get_dof_info(mf_component).vector_exchanger.get();
    }



    /**
     * Start update_ghost_value for serial vectors
     */
    template <typename VectorType,
              std::enable_if_t<is_not_parallel_vector<VectorType>, VectorType>
                * = nullptr>
    void
    update_ghost_values_start(const unsigned int /*component_in_block_vector*/,
                              const VectorType & /*vec*/)
    {}


    /**
     * Start update_ghost_value for vectors that do not support
     * the split into _start() and finish() stages
     */
    template <typename VectorType,
              std::enable_if_t<!has_update_ghost_values_start<VectorType> &&
                                 !is_not_parallel_vector<VectorType>,
                               VectorType> * = nullptr>
    void
    update_ghost_values_start(const unsigned int component_in_block_vector,
                              const VectorType  &vec)
    {
      (void)component_in_block_vector;
      const bool ghosts_set = vec.has_ghost_elements();

      Assert(matrix_free.get_task_info().allow_ghosted_vectors_in_loops ||
               ghosts_set == false,
             ExcNotImplemented());

      if (ghosts_set)
        {
          ghosts_were_set = true;
          return;
        }

      vec.update_ghost_values();
    }



    /**
     * Start update_ghost_value for vectors that _do_ support
     * the split into _start() and finish() stages, but don't support
     * exchange on a subset of DoFs
     */
    template <typename VectorType,
              std::enable_if_t<has_update_ghost_values_start<VectorType> &&
                                 !has_exchange_on_subset<VectorType>,
                               VectorType> * = nullptr>
    void
    update_ghost_values_start(const unsigned int component_in_block_vector,
                              const VectorType  &vec)
    {
      (void)component_in_block_vector;
      const bool ghosts_set = vec.has_ghost_elements();

      Assert(matrix_free.get_task_info().allow_ghosted_vectors_in_loops ||
               ghosts_set == false,
             ExcNotImplemented());

      if (ghosts_set)
        {
          ghosts_were_set = true;
          return;
        }

      vec.update_ghost_values_start(component_in_block_vector + channel_shift);
    }



    /**
     * Finally, start update_ghost_value for vectors that _do_ support
     * the split into _start() and finish() stages and also support
     * exchange on a subset of DoFs,
     * i.e. LinearAlgebra::distributed::Vector
     */
    template <typename VectorType,
              std::enable_if_t<has_update_ghost_values_start<VectorType> &&
                                 has_exchange_on_subset<VectorType>,
                               VectorType> * = nullptr>
    void
    update_ghost_values_start(const unsigned int component_in_block_vector,
                              const VectorType  &vec)
    {
      static_assert(std::is_same_v<Number, typename VectorType::value_type>,
                    "Type mismatch between VectorType and VectorDataExchange");
      (void)component_in_block_vector;
      const bool ghosts_set = vec.has_ghost_elements();

      Assert(matrix_free.get_task_info().allow_ghosted_vectors_in_loops ||
               ghosts_set == false,
             ExcNotImplemented());

      if (ghosts_set)
        {
          ghosts_were_set = true;
          return;
        }

      if (vec.size() != 0)
        {
#  ifdef DEAL_II_WITH_MPI
          const unsigned int mf_component = find_vector_in_mf(vec);

          const auto &part = get_partitioner(mf_component);

          if (part.n_ghost_indices() == 0 && part.n_import_indices() == 0 &&
              part.n_import_sm_procs() == 0)
            return;

          tmp_data[component_in_block_vector] =
            matrix_free.acquire_scratch_data_non_threadsafe();
          tmp_data[component_in_block_vector]->resize_fast(
            part.n_import_indices());
          AssertDimension(requests.size(), tmp_data.size());

          part.export_to_ghosted_array_start(
            component_in_block_vector * 2 + channel_shift,
            ArrayView<const Number>(vec.begin(), part.locally_owned_size()),
            vec.shared_vector_data(),
            ArrayView<Number>(const_cast<Number *>(vec.begin()) +
                                part.locally_owned_size(),
                              matrix_free.get_dof_info(mf_component)
                                .vector_partitioner->n_ghost_indices()),
            ArrayView<Number>(tmp_data[component_in_block_vector]->begin(),
                              part.n_import_indices()),
            this->requests[component_in_block_vector]);
#  endif
        }
    }



    /**
     * Finish update_ghost_value for vectors that do not support
     * the split into _start() and finish() stages and serial vectors
     */
    template <typename VectorType,
              std::enable_if_t<!has_update_ghost_values_start<VectorType>,
                               VectorType> * = nullptr>
    void
    update_ghost_values_finish(const unsigned int /*component_in_block_vector*/,
                               const VectorType & /*vec*/)
    {}



    /**
     * Finish update_ghost_value for vectors that _do_ support
     * the split into _start() and finish() stages, but don't support
     * exchange on a subset of DoFs
     */
    template <typename VectorType,
              std::enable_if_t<has_update_ghost_values_start<VectorType> &&
                                 !has_exchange_on_subset<VectorType>,
                               VectorType> * = nullptr>
    void
    update_ghost_values_finish(const unsigned int component_in_block_vector,
                               const VectorType  &vec)
    {
      (void)component_in_block_vector;

      if (ghosts_were_set)
        return;

      vec.update_ghost_values_finish();
    }



    /**
     * Finish update_ghost_value for vectors that _do_ support
     * the split into _start() and finish() stages and also support
     * exchange on a subset of DoFs,
     * i.e. LinearAlgebra::distributed::Vector
     */
    template <typename VectorType,
              std::enable_if_t<has_update_ghost_values_start<VectorType> &&
                                 has_exchange_on_subset<VectorType>,
                               VectorType> * = nullptr>
    void
    update_ghost_values_finish(const unsigned int component_in_block_vector,
                               const VectorType  &vec)
    {
      static_assert(std::is_same_v<Number, typename VectorType::value_type>,
                    "Type mismatch between VectorType and VectorDataExchange");
      (void)component_in_block_vector;

      if (ghosts_were_set)
        return;

      if (vec.size() != 0)
        {
#  ifdef DEAL_II_WITH_MPI
          AssertIndexRange(component_in_block_vector, tmp_data.size());
          AssertDimension(requests.size(), tmp_data.size());

          const unsigned int mf_component = find_vector_in_mf(vec);

          const auto &part = get_partitioner(mf_component);

          if (part.n_ghost_indices() != 0 || part.n_import_indices() != 0 ||
              part.n_import_sm_procs() != 0)
            {
              part.export_to_ghosted_array_finish(
                ArrayView<const Number>(vec.begin(), part.locally_owned_size()),
                vec.shared_vector_data(),
                ArrayView<Number>(const_cast<Number *>(vec.begin()) +
                                    part.locally_owned_size(),
                                  matrix_free.get_dof_info(mf_component)
                                    .vector_partitioner->n_ghost_indices()),
                this->requests[component_in_block_vector]);

              matrix_free.release_scratch_data_non_threadsafe(
                tmp_data[component_in_block_vector]);
              tmp_data[component_in_block_vector] = nullptr;
            }
#  endif
        }
      // let vector know that ghosts are being updated and we can read from
      // them
      vec.set_ghost_state(true);
    }



    /**
     * Start compress for serial vectors
     */
    template <typename VectorType,
              std::enable_if_t<is_not_parallel_vector<VectorType>, VectorType>
                * = nullptr>
    void
    compress_start(const unsigned int /*component_in_block_vector*/,
                   VectorType & /*vec*/)
    {}



    /**
     * Start compress for vectors that do not support
     * the split into _start() and finish() stages
     */
    template <typename VectorType,
              std::enable_if_t<!has_compress_start<VectorType> &&
                                 !is_not_parallel_vector<VectorType>,
                               VectorType> * = nullptr>
    void
    compress_start(const unsigned int component_in_block_vector,
                   VectorType        &vec)
    {
      (void)component_in_block_vector;
      Assert(vec.has_ghost_elements() == false, ExcNotImplemented());
      vec.compress(VectorOperation::add);
    }



    /**
     * Start compress for vectors that _do_ support
     * the split into _start() and finish() stages, but don't support
     * exchange on a subset of DoFs
     */
    template <typename VectorType,
              std::enable_if_t<has_compress_start<VectorType> &&
                                 !has_exchange_on_subset<VectorType>,
                               VectorType> * = nullptr>
    void
    compress_start(const unsigned int component_in_block_vector,
                   VectorType        &vec)
    {
      (void)component_in_block_vector;
      Assert(vec.has_ghost_elements() == false, ExcNotImplemented());
      vec.compress_start(component_in_block_vector + channel_shift);
    }



    /**
     * Start compress for vectors that _do_ support
     * the split into _start() and finish() stages and also support
     * exchange on a subset of DoFs,
     * i.e. LinearAlgebra::distributed::Vector
     */
    template <typename VectorType,
              std::enable_if_t<has_compress_start<VectorType> &&
                                 has_exchange_on_subset<VectorType>,
                               VectorType> * = nullptr>
    void
    compress_start(const unsigned int component_in_block_vector,
                   VectorType        &vec)
    {
      static_assert(std::is_same_v<Number, typename VectorType::value_type>,
                    "Type mismatch between VectorType and VectorDataExchange");
      (void)component_in_block_vector;
      Assert(vec.has_ghost_elements() == false, ExcNotImplemented());

      if (vec.size() != 0)
        {
#  ifdef DEAL_II_WITH_MPI
          const unsigned int mf_component = find_vector_in_mf(vec);

          const auto &part = get_partitioner(mf_component);

          if (part.n_ghost_indices() == 0 && part.n_import_indices() == 0 &&
              part.n_import_sm_procs() == 0)
            return;

          tmp_data[component_in_block_vector] =
            matrix_free.acquire_scratch_data_non_threadsafe();
          tmp_data[component_in_block_vector]->resize_fast(
            part.n_import_indices());
          AssertDimension(requests.size(), tmp_data.size());

          part.import_from_ghosted_array_start(
            VectorOperation::add,
            component_in_block_vector * 2 + channel_shift,
            ArrayView<Number>(vec.begin(), part.locally_owned_size()),
            vec.shared_vector_data(),
            ArrayView<Number>(vec.begin() + part.locally_owned_size(),
                              matrix_free.get_dof_info(mf_component)
                                .vector_partitioner->n_ghost_indices()),
            ArrayView<Number>(tmp_data[component_in_block_vector]->begin(),
                              part.n_import_indices()),
            this->requests[component_in_block_vector]);
#  endif
        }
    }



    /**
     * Finish compress for vectors that do not support
     * the split into _start() and finish() stages and serial vectors
     */
    template <
      typename VectorType,
      std::enable_if_t<!has_compress_start<VectorType>, VectorType> * = nullptr>
    void
    compress_finish(const unsigned int /*component_in_block_vector*/,
                    VectorType & /*vec*/)
    {}



    /**
     * Finish compress for vectors that _do_ support
     * the split into _start() and finish() stages, but don't support
     * exchange on a subset of DoFs
     */
    template <typename VectorType,
              std::enable_if_t<has_compress_start<VectorType> &&
                                 !has_exchange_on_subset<VectorType>,
                               VectorType> * = nullptr>
    void
    compress_finish(const unsigned int component_in_block_vector,
                    VectorType        &vec)
    {
      (void)component_in_block_vector;
      vec.compress_finish(VectorOperation::add);
    }



    /**
     * Start compress for vectors that _do_ support
     * the split into _start() and finish() stages and also support
     * exchange on a subset of DoFs,
     * i.e. LinearAlgebra::distributed::Vector
     */
    template <typename VectorType,
              std::enable_if_t<has_compress_start<VectorType> &&
                                 has_exchange_on_subset<VectorType>,
                               VectorType> * = nullptr>
    void
    compress_finish(const unsigned int component_in_block_vector,
                    VectorType        &vec)
    {
      static_assert(std::is_same_v<Number, typename VectorType::value_type>,
                    "Type mismatch between VectorType and VectorDataExchange");
      (void)component_in_block_vector;
      if (vec.size() != 0)
        {
#  ifdef DEAL_II_WITH_MPI
          AssertIndexRange(component_in_block_vector, tmp_data.size());
          AssertDimension(requests.size(), tmp_data.size());

          const unsigned int mf_component = find_vector_in_mf(vec);

          const auto &part = get_partitioner(mf_component);

          if (part.n_ghost_indices() != 0 || part.n_import_indices() != 0 ||
              part.n_import_sm_procs() != 0)
            {
              part.import_from_ghosted_array_finish(
                VectorOperation::add,
                ArrayView<Number>(vec.begin(), part.locally_owned_size()),
                vec.shared_vector_data(),
                ArrayView<Number>(vec.begin() + part.locally_owned_size(),
                                  matrix_free.get_dof_info(mf_component)
                                    .vector_partitioner->n_ghost_indices()),
                ArrayView<const Number>(
                  tmp_data[component_in_block_vector]->begin(),
                  part.n_import_indices()),
                this->requests[component_in_block_vector]);

              matrix_free.release_scratch_data_non_threadsafe(
                tmp_data[component_in_block_vector]);
              tmp_data[component_in_block_vector] = nullptr;
            }

          if (Utilities::MPI::job_supports_mpi())
            {
              const int ierr =
                MPI_Barrier(matrix_free.get_task_info().communicator_sm);
              AssertThrowMPI(ierr);
            }
#  endif
        }
    }



    /**
     * Reset all ghost values for serial vectors
     */
    template <typename VectorType,
              std::enable_if_t<is_not_parallel_vector<VectorType>, VectorType>
                * = nullptr>
    void
    reset_ghost_values(const VectorType & /*vec*/) const
    {}



    /**
     * Reset all ghost values for vector that don't support
     * exchange on a subset of DoFs
     */
    template <typename VectorType,
              std::enable_if_t<!has_exchange_on_subset<VectorType> &&
                                 !is_not_parallel_vector<VectorType>,
                               VectorType> * = nullptr>
    void
    reset_ghost_values(const VectorType &vec) const
    {
      if (ghosts_were_set == true)
        return;

      vec.zero_out_ghost_values();
    }



    /**
     * Reset all ghost values for vector that _do_ support
     * exchange on a subset of DoFs, i.e.
     * LinearAlgebra::distributed::Vector
     */
    template <typename VectorType,
              std::enable_if_t<has_exchange_on_subset<VectorType>, VectorType>
                * = nullptr>
    void
    reset_ghost_values(const VectorType &vec) const
    {
      static_assert(std::is_same_v<Number, typename VectorType::value_type>,
                    "Type mismatch between VectorType and VectorDataExchange");
      if (ghosts_were_set == true)
        return;

      if (vec.size() != 0)
        {
#  ifdef DEAL_II_WITH_MPI
          AssertDimension(requests.size(), tmp_data.size());

          const unsigned int mf_component = find_vector_in_mf(vec);

          const auto &part = get_partitioner(mf_component);

          if (part.n_ghost_indices() > 0)
            {
              part.reset_ghost_values(
                ArrayView<Number>(const_cast<VectorType &>(vec).begin() +
                                    part.locally_owned_size(),
                                  matrix_free.get_dof_info(mf_component)
                                    .vector_partitioner->n_ghost_indices()));
            }

#  endif
        }
      // let vector know that it's not ghosted anymore
      vec.set_ghost_state(false);
    }



    /**
     * Zero out vector region for vector that _do_ support
     * exchange on a subset of DoFs <==> begin() + ind == local_element(ind),
     * i.e. LinearAlgebra::distributed::Vector
     */
    template <typename VectorType,
              std::enable_if_t<has_exchange_on_subset<VectorType>, VectorType>
                * = nullptr>
    void
    zero_vector_region(const unsigned int range_index, VectorType &vec) const
    {
      static_assert(std::is_same_v<Number, typename VectorType::value_type>,
                    "Type mismatch between VectorType and VectorDataExchange");
      if (range_index == numbers::invalid_unsigned_int)
        vec = Number();
      else
        {
          const unsigned int mf_component = find_vector_in_mf(vec, false);
          const internal::MatrixFreeFunctions::DoFInfo &dof_info =
            matrix_free.get_dof_info(mf_component);
          Assert(dof_info.vector_zero_range_list_index.empty() == false,
                 ExcNotInitialized());

          Assert(vec.partitioners_are_compatible(*dof_info.vector_partitioner),
                 ExcInternalError());
          AssertIndexRange(range_index,
                           dof_info.vector_zero_range_list_index.size() - 1);
          for (unsigned int id =
                 dof_info.vector_zero_range_list_index[range_index];
               id != dof_info.vector_zero_range_list_index[range_index + 1];
               ++id)
            std::memset(vec.begin() + dof_info.vector_zero_range_list[id].first,
                        0,
                        (dof_info.vector_zero_range_list[id].second -
                         dof_info.vector_zero_range_list[id].first) *
                          sizeof(Number));
        }
    }



    /**
     * Zero out vector region for vector that do _not_ support exchange on a
     * subset of DoFs <==> begin() + ind == local_element(ind) but are still a
     * vector type
     */
    template <typename VectorType,
              std::enable_if_t<!has_exchange_on_subset<VectorType>, VectorType>
                * = nullptr,
              std::enable_if_t<has_assignment_operator<VectorType>, VectorType>
                * = nullptr>
    void
    zero_vector_region(const unsigned int range_index, VectorType &vec) const
    {
      if (range_index == numbers::invalid_unsigned_int || range_index == 0)
        {
          if constexpr (std::is_same_v<
                          ArrayView<typename VectorType::value_type>,
                          VectorType>)
            {
              for (unsigned int i = 0; i < vec.size(); ++i)
                vec[i] = typename VectorType::value_type();
            }
          else
            vec = typename VectorType::value_type();
        }
    }



    /**
     * Zero out vector region for non-vector types, i.e., classes that do not
     * have operator=(const VectorType::value_type)
     */
    void
    zero_vector_region(const unsigned int, ...) const
    {
      Assert(false,
             ExcNotImplemented(
               "Zeroing is only implemented for vector types "
               "which provide operator=(const VectorType::value_type)"));
    }



    const dealii::MatrixFree<dim, Number, VectorizedArrayType> &matrix_free;
    const typename dealii::MatrixFree<dim, Number, VectorizedArrayType>::
      DataAccessOnFaces vector_face_access;
    bool                ghosts_were_set;
#  ifdef DEAL_II_WITH_MPI
    std::vector<AlignedVector<Number> *>  tmp_data;
    std::vector<std::vector<MPI_Request>> requests;
#  endif
  }; // VectorDataExchange

  template <typename VectorStruct>
  unsigned int
  n_components(const VectorStruct &vec);

  template <typename VectorStruct>
  unsigned int
  n_components_block(const VectorStruct &vec, const std::bool_constant<true>)
  {
    unsigned int components = 0;
    for (unsigned int bl = 0; bl < vec.n_blocks(); ++bl)
      components += n_components(vec.block(bl));
    return components;
  }

  template <typename VectorStruct>
  unsigned int
  n_components_block(const VectorStruct &, const std::bool_constant<false>)
  {
    return 1;
  }

  template <typename VectorStruct>
  unsigned int
  n_components(const VectorStruct &vec)
  {
    return n_components_block(
      vec, std::bool_constant<IsBlockVector<VectorStruct>::value>());
  }

  template <typename VectorStruct>
  inline unsigned int
  n_components(const std::vector<VectorStruct> &vec)
  {
    unsigned int components = 0;
    for (unsigned int comp = 0; comp < vec.size(); ++comp)
      components += n_components_block(
        vec[comp], std::bool_constant<IsBlockVector<VectorStruct>::value>());
    return components;
  }

  template <typename VectorStruct>
  inline unsigned int
  n_components(const std::vector<VectorStruct *> &vec)
  {
    unsigned int components = 0;
    for (unsigned int comp = 0; comp < vec.size(); ++comp)
      components += n_components_block(
        *vec[comp], std::bool_constant<IsBlockVector<VectorStruct>::value>());
    return components;
  }



  // A helper function to identify block vectors with many components where we
  // should not try to overlap computations and communication because there
  // would be too many outstanding communication requests.

  // default value for vectors that do not have communication_block_size
  template <typename VectorStruct,
            std::enable_if_t<!has_communication_block_size<VectorStruct>,
                             VectorStruct> * = nullptr>
  constexpr unsigned int
  get_communication_block_size(const VectorStruct &)
  {
    return numbers::invalid_unsigned_int;
  }



  template <typename VectorStruct,
            std::enable_if_t<has_communication_block_size<VectorStruct>,
                             VectorStruct> * = nullptr>
  constexpr unsigned int
  get_communication_block_size(const VectorStruct &)
  {
    return VectorStruct::communication_block_size;
  }



  template <typename VectorType,
            std::enable_if_t<is_not_parallel_vector<VectorType>, VectorType> * =
              nullptr>
  bool
  has_ghost_elements(const VectorType &vec)
  {
    (void)vec;
    return false;
  }



  template <typename VectorType,
            std::enable_if_t<!is_not_parallel_vector<VectorType>, VectorType>
              * = nullptr>
  bool
  has_ghost_elements(const VectorType &vec)
  {
    return vec.has_ghost_elements();
  }



  // --------------------------------------------------------------------------
  // below we have wrappers to distinguish between block and non-block vectors.
  // --------------------------------------------------------------------------

  //
  // update_ghost_values_start
  //

  // update_ghost_values for block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            std::enable_if_t<IsBlockVector<VectorStruct>::value, VectorStruct>
              * = nullptr>
  void
  update_ghost_values_start(
    const VectorStruct                                   &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger,
    const unsigned int                                    channel = 0)
  {
    if (get_communication_block_size(vec) < vec.n_blocks())
      {
        const bool ghosts_set = vec.has_ghost_elements();

        Assert(exchanger.matrix_free.get_task_info()
                   .allow_ghosted_vectors_in_loops ||
                 ghosts_set == false,
               ExcNotImplemented());

        if (ghosts_set)
          {
            exchanger.ghosts_were_set = true;
            return;
          }

        vec.update_ghost_values();
      }
    else
      {
        for (unsigned int i = 0; i < vec.n_blocks(); ++i)
          update_ghost_values_start(vec.block(i), exchanger, channel + i);
      }
  }



  // update_ghost_values for non-block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            std::enable_if_t<!IsBlockVector<VectorStruct>::value, VectorStruct>
              * = nullptr>
  void
  update_ghost_values_start(
    const VectorStruct                                   &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger,
    const unsigned int                                    channel = 0)
  {
    exchanger.update_ghost_values_start(channel, vec);
  }



  // update_ghost_values_start() for vector of vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  update_ghost_values_start(
    const std::vector<VectorStruct>                      &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    unsigned int component_index = 0;
    for (unsigned int comp = 0; comp < vec.size(); ++comp)
      {
        update_ghost_values_start(vec[comp], exchanger, component_index);
        component_index += n_components(vec[comp]);
      }
  }



  // update_ghost_values_start() for vector of pointers to vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  update_ghost_values_start(
    const std::vector<VectorStruct *>                    &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    unsigned int component_index = 0;
    for (unsigned int comp = 0; comp < vec.size(); ++comp)
      {
        update_ghost_values_start(*vec[comp], exchanger, component_index);
        component_index += n_components(*vec[comp]);
      }
  }



  //
  // update_ghost_values_finish
  //

  // for block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            std::enable_if_t<IsBlockVector<VectorStruct>::value, VectorStruct>
              * = nullptr>
  void
  update_ghost_values_finish(
    const VectorStruct                                   &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger,
    const unsigned int                                    channel = 0)
  {
    if (get_communication_block_size(vec) < vec.n_blocks())
      {
        // do nothing, everything has already been completed in the _start()
        // call
      }
    else
      for (unsigned int i = 0; i < vec.n_blocks(); ++i)
        update_ghost_values_finish(vec.block(i), exchanger, channel + i);
  }



  // for non-block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            std::enable_if_t<!IsBlockVector<VectorStruct>::value, VectorStruct>
              * = nullptr>
  void
  update_ghost_values_finish(
    const VectorStruct                                   &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger,
    const unsigned int                                    channel = 0)
  {
    exchanger.update_ghost_values_finish(channel, vec);
  }



  // for vector of vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  update_ghost_values_finish(
    const std::vector<VectorStruct>                      &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    unsigned int component_index = 0;
    for (unsigned int comp = 0; comp < vec.size(); ++comp)
      {
        update_ghost_values_finish(vec[comp], exchanger, component_index);
        component_index += n_components(vec[comp]);
      }
  }



  // for vector of pointers to vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  update_ghost_values_finish(
    const std::vector<VectorStruct *>                    &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    unsigned int component_index = 0;
    for (unsigned int comp = 0; comp < vec.size(); ++comp)
      {
        update_ghost_values_finish(*vec[comp], exchanger, component_index);
        component_index += n_components(*vec[comp]);
      }
  }



  //
  // compress_start
  //

  // for block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            std::enable_if_t<IsBlockVector<VectorStruct>::value, VectorStruct>
              * = nullptr>
  inline void
  compress_start(
    VectorStruct                                         &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger,
    const unsigned int                                    channel = 0)
  {
    if (get_communication_block_size(vec) < vec.n_blocks())
      vec.compress(VectorOperation::add);
    else
      for (unsigned int i = 0; i < vec.n_blocks(); ++i)
        compress_start(vec.block(i), exchanger, channel + i);
  }



  // for non-block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            std::enable_if_t<!IsBlockVector<VectorStruct>::value, VectorStruct>
              * = nullptr>
  inline void
  compress_start(
    VectorStruct                                         &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger,
    const unsigned int                                    channel = 0)
  {
    exchanger.compress_start(channel, vec);
  }



  // for std::vector of vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  compress_start(
    std::vector<VectorStruct>                            &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    unsigned int component_index = 0;
    for (unsigned int comp = 0; comp < vec.size(); ++comp)
      {
        compress_start(vec[comp], exchanger, component_index);
        component_index += n_components(vec[comp]);
      }
  }



  // for std::vector of pointer to vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  compress_start(
    std::vector<VectorStruct *>                          &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    unsigned int component_index = 0;
    for (unsigned int comp = 0; comp < vec.size(); ++comp)
      {
        compress_start(*vec[comp], exchanger, component_index);
        component_index += n_components(*vec[comp]);
      }
  }



  //
  // compress_finish
  //

  // for block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            std::enable_if_t<IsBlockVector<VectorStruct>::value, VectorStruct>
              * = nullptr>
  inline void
  compress_finish(
    VectorStruct                                         &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger,
    const unsigned int                                    channel = 0)
  {
    if (get_communication_block_size(vec) < vec.n_blocks())
      {
        // do nothing, everything has already been completed in the _start()
        // call
      }
    else
      for (unsigned int i = 0; i < vec.n_blocks(); ++i)
        compress_finish(vec.block(i), exchanger, channel + i);
  }



  // for non-block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            std::enable_if_t<!IsBlockVector<VectorStruct>::value, VectorStruct>
              * = nullptr>
  inline void
  compress_finish(
    VectorStruct                                         &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger,
    const unsigned int                                    channel = 0)
  {
    exchanger.compress_finish(channel, vec);
  }



  // for std::vector of vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  compress_finish(
    std::vector<VectorStruct>                            &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    unsigned int component_index = 0;
    for (unsigned int comp = 0; comp < vec.size(); ++comp)
      {
        compress_finish(vec[comp], exchanger, component_index);
        component_index += n_components(vec[comp]);
      }
  }



  // for std::vector of pointer to vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  compress_finish(
    std::vector<VectorStruct *>                          &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    unsigned int component_index = 0;
    for (unsigned int comp = 0; comp < vec.size(); ++comp)
      {
        compress_finish(*vec[comp], exchanger, component_index);
        component_index += n_components(*vec[comp]);
      }
  }



  //
  // reset_ghost_values:
  //
  // if the input vector did not have ghosts imported, clear them here again
  // in order to avoid subsequent operations e.g. in linear solvers to work
  // with ghosts all the time
  //

  // for block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            std::enable_if_t<IsBlockVector<VectorStruct>::value, VectorStruct>
              * = nullptr>
  inline void
  reset_ghost_values(
    const VectorStruct                                   &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    // return immediately if there is nothing to do.
    if (exchanger.ghosts_were_set == true)
      return;

    for (unsigned int i = 0; i < vec.n_blocks(); ++i)
      reset_ghost_values(vec.block(i), exchanger);
  }



  // for non-block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            std::enable_if_t<!IsBlockVector<VectorStruct>::value, VectorStruct>
              * = nullptr>
  inline void
  reset_ghost_values(
    const VectorStruct                                   &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    // return immediately if there is nothing to do.
    if (exchanger.ghosts_were_set == true)
      return;

    exchanger.reset_ghost_values(vec);
  }



  // for std::vector of vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  reset_ghost_values(
    const std::vector<VectorStruct>                      &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    // return immediately if there is nothing to do.
    if (exchanger.ghosts_were_set == true)
      return;

    for (unsigned int comp = 0; comp < vec.size(); ++comp)
      reset_ghost_values(vec[comp], exchanger);
  }



  // for std::vector of pointer to vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  reset_ghost_values(
    const std::vector<VectorStruct *>                    &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    // return immediately if there is nothing to do.
    if (exchanger.ghosts_were_set == true)
      return;

    for (unsigned int comp = 0; comp < vec.size(); ++comp)
      reset_ghost_values(*vec[comp], exchanger);
  }



  //
  // zero_vector_region
  //

  // for block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            std::enable_if_t<IsBlockVector<VectorStruct>::value, VectorStruct>
              * = nullptr>
  inline void
  zero_vector_region(
    const unsigned int                                    range_index,
    VectorStruct                                         &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    for (unsigned int i = 0; i < vec.n_blocks(); ++i)
      exchanger.zero_vector_region(range_index, vec.block(i));
  }



  // for non-block vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType,
            std::enable_if_t<!IsBlockVector<VectorStruct>::value, VectorStruct>
              * = nullptr>
  inline void
  zero_vector_region(
    const unsigned int                                    range_index,
    VectorStruct                                         &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    exchanger.zero_vector_region(range_index, vec);
  }



  // for std::vector of vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  zero_vector_region(
    const unsigned int                                    range_index,
    std::vector<VectorStruct>                            &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    for (unsigned int comp = 0; comp < vec.size(); ++comp)
      zero_vector_region(range_index, vec[comp], exchanger);
  }



  // for std::vector of pointers to vectors
  template <int dim,
            typename VectorStruct,
            typename Number,
            typename VectorizedArrayType>
  inline void
  zero_vector_region(
    const unsigned int                                    range_index,
    std::vector<VectorStruct *>                          &vec,
    VectorDataExchange<dim, Number, VectorizedArrayType> &exchanger)
  {
    for (unsigned int comp = 0; comp < vec.size(); ++comp)
      zero_vector_region(range_index, *vec[comp], exchanger);
  }



  // Apply a unit matrix operation to constrained DoFs: Default cases where we
  // cannot detect a LinearAlgebra::distributed::Vector, we do not do
  // anything, else we apply the constraints as a unit operation
  template <typename VectorStruct1, typename VectorStruct2>
  inline void
  apply_operation_to_constrained_dofs(const std::vector<unsigned int> &,
                                      const VectorStruct1 &,
                                      VectorStruct2 &)
  {}

  template <typename Number>
  inline void
  apply_operation_to_constrained_dofs(
    const std::vector<unsigned int>                  &constrained_dofs,
    const LinearAlgebra::distributed::Vector<Number> &src,
    LinearAlgebra::distributed::Vector<Number>       &dst)
  {
    for (const unsigned int i : constrained_dofs)
      dst.local_element(i) = src.local_element(i);
  }


  namespace MatrixFreeFunctions
  {
    // struct to select between a const interface and a non-const interface
    // for MFWorker
    template <typename, typename, typename, typename, bool>
    struct InterfaceSelector
    {};

    // Version for constant functions
    template <typename MF,
              typename InVector,
              typename OutVector,
              typename Container>
    struct InterfaceSelector<MF, InVector, OutVector, Container, true>
    {
      using function_type = void (Container::*)(
        const MF &,
        OutVector &,
        const InVector &,
        const std::pair<unsigned int, unsigned int> &) const;
    };

    // Version for non-constant functions
    template <typename MF,
              typename InVector,
              typename OutVector,
              typename Container>
    struct InterfaceSelector<MF, InVector, OutVector, Container, false>
    {
      using function_type =
        void (Container::*)(const MF &,
                            OutVector &,
                            const InVector &,
                            const std::pair<unsigned int, unsigned int> &);
    };
  } // namespace MatrixFreeFunctions



  // A implementation class for the worker object that runs the various
  // operations we want to perform during the matrix-free loop
  template <typename MF,
            typename InVector,
            typename OutVector,
            typename Container,
            bool is_constant>
  class MFWorker : public MFWorkerInterface
  {
  public:
    // An alias to make the arguments further down more readable
    using function_type = typename MatrixFreeFunctions::
      InterfaceSelector<MF, InVector, OutVector, Container, is_constant>::
        function_type;

    // constructor, binds all the arguments to this class
    MFWorker(const MF                            &matrix_free,
             const InVector                      &src,
             OutVector                           &dst,
             const bool                           zero_dst_vector_setting,
             const Container                     &container,
             function_type                        cell_function,
             function_type                        face_function,
             function_type                        boundary_function,
             const typename MF::DataAccessOnFaces src_vector_face_access =
               MF::DataAccessOnFaces::none,
             const typename MF::DataAccessOnFaces dst_vector_face_access =
               MF::DataAccessOnFaces::none,
             const std::function<void(const unsigned int, const unsigned int)>
               &operation_before_loop = {},
             const std::function<void(const unsigned int, const unsigned int)>
                               &operation_after_loop       = {},
             const unsigned int dof_handler_index_pre_post = 0)
      : matrix_free(matrix_free)
      , container(const_cast<Container &>(container))
      , cell_function(cell_function)
      , face_function(face_function)
      , boundary_function(boundary_function)
      , src(src)
      , dst(dst)
      , src_data_exchanger(matrix_free,
                           src_vector_face_access,
                           n_components(src))
      , dst_data_exchanger(matrix_free,
                           dst_vector_face_access,
                           n_components(dst))
      , src_and_dst_are_same(PointerComparison::equal(&src, &dst))
      , zero_dst_vector_setting(zero_dst_vector_setting &&
                                !src_and_dst_are_same)
      , operation_before_loop(operation_before_loop)
      , operation_after_loop(operation_after_loop)
      , dof_handler_index_pre_post(dof_handler_index_pre_post)
    {
      Assert(!has_ghost_elements(dst),
             ExcMessage("The destination vector passed to the matrix-free "
                        "loop is ghosted. This is not allowed."));
    }

    // Runs the cell work. If no function is given, nothing is done
    virtual void
    cell(const std::pair<unsigned int, unsigned int> &cell_range) override
    {
      if (cell_function != nullptr && cell_range.second > cell_range.first)
        for (unsigned int i = 0; i < matrix_free.n_active_fe_indices(); ++i)
          {
            const auto cell_subrange =
              matrix_free.create_cell_subrange_hp_by_index(cell_range, i);

            if (cell_subrange.second <= cell_subrange.first)
              continue;

            (container.*
             cell_function)(matrix_free, this->dst, this->src, cell_subrange);
          }
    }

    virtual void
    cell(const unsigned int range_index) override
    {
      process_range(cell_function,
                    matrix_free.get_task_info().cell_partition_data_hp_ptr,
                    matrix_free.get_task_info().cell_partition_data_hp,
                    range_index);
    }

    virtual void
    face(const unsigned int range_index) override
    {
      process_range(face_function,
                    matrix_free.get_task_info().face_partition_data_hp_ptr,
                    matrix_free.get_task_info().face_partition_data_hp,
                    range_index);
    }

    virtual void
    boundary(const unsigned int range_index) override
    {
      process_range(boundary_function,
                    matrix_free.get_task_info().boundary_partition_data_hp_ptr,
                    matrix_free.get_task_info().boundary_partition_data_hp,
                    range_index);
    }

  private:
    void
    process_range(const function_type             &fu,
                  const std::vector<unsigned int> &ptr,
                  const std::vector<unsigned int> &data,
                  const unsigned int               range_index)
    {
      if (fu == nullptr)
        return;

      AssertIndexRange(range_index + 1, ptr.size());
      for (unsigned int i = ptr[range_index]; i < ptr[range_index + 1]; ++i)
        {
          AssertIndexRange(2 * i + 1, data.size());
          (container.*fu)(matrix_free,
                          this->dst,
                          this->src,
                          std::make_pair(data[2 * i], data[2 * i + 1]));
        }
    }

  public:
    // Starts the communication for the update ghost values operation. We
    // cannot call this update if ghost and destination are the same because
    // that would introduce spurious entries in the destination (there is also
    // the problem that reading from a vector that we also write to is usually
    // not intended in case there is overlap, but this is up to the
    // application code to decide and we cannot catch this case here).
    virtual void
    vector_update_ghosts_start() override
    {
      if (!src_and_dst_are_same)
        internal::update_ghost_values_start(src, src_data_exchanger);
    }

    // Finishes the communication for the update ghost values operation
    virtual void
    vector_update_ghosts_finish() override
    {
      if (!src_and_dst_are_same)
        internal::update_ghost_values_finish(src, src_data_exchanger);
    }

    // Starts the communication for the vector compress operation
    virtual void
    vector_compress_start() override
    {
      internal::compress_start(dst, dst_data_exchanger);
    }

    // Finishes the communication for the vector compress operation
    virtual void
    vector_compress_finish() override
    {
      internal::compress_finish(dst, dst_data_exchanger);
      if (!src_and_dst_are_same)
        internal::reset_ghost_values(src, src_data_exchanger);
    }

    // Zeros the given input vector
    virtual void
    zero_dst_vector_range(const unsigned int range_index) override
    {
      if (zero_dst_vector_setting)
        internal::zero_vector_region(range_index, dst, dst_data_exchanger);
    }

    virtual void
    cell_loop_pre_range(const unsigned int range_index) override
    {
      if (operation_before_loop)
        {
          const internal::MatrixFreeFunctions::DoFInfo &dof_info =
            matrix_free.get_dof_info(dof_handler_index_pre_post);
          if (range_index == numbers::invalid_unsigned_int)
            {
              // Case with threaded loop -> currently no overlap implemented
              dealii::parallel::apply_to_subranges(
                0U,
                dof_info.vector_partitioner->locally_owned_size(),
                operation_before_loop,
                internal::VectorImplementation::minimum_parallel_grain_size);
            }
          else
            {
              AssertIndexRange(range_index,
                               dof_info.cell_loop_pre_list_index.size() - 1);
              for (unsigned int id =
                     dof_info.cell_loop_pre_list_index[range_index];
                   id != dof_info.cell_loop_pre_list_index[range_index + 1];
                   ++id)
                operation_before_loop(dof_info.cell_loop_pre_list[id].first,
                                      dof_info.cell_loop_pre_list[id].second);
            }
        }
    }

    virtual void
    cell_loop_post_range(const unsigned int range_index) override
    {
      if (operation_after_loop)
        {
          // Run unit matrix operation on constrained dofs if we are at the
          // last range
          const std::vector<unsigned int> &partition_row_index =
            matrix_free.get_task_info().partition_row_index;
          if (range_index ==
              partition_row_index[partition_row_index.size() - 2] - 1)
            apply_operation_to_constrained_dofs(
              matrix_free.get_constrained_dofs(dof_handler_index_pre_post),
              src,
              dst);

          const internal::MatrixFreeFunctions::DoFInfo &dof_info =
            matrix_free.get_dof_info(dof_handler_index_pre_post);
          if (range_index == numbers::invalid_unsigned_int)
            {
              // Case with threaded loop -> currently no overlap implemented
              dealii::parallel::apply_to_subranges(
                0U,
                dof_info.vector_partitioner->locally_owned_size(),
                operation_after_loop,
                internal::VectorImplementation::minimum_parallel_grain_size);
            }
          else
            {
              AssertIndexRange(range_index,
                               dof_info.cell_loop_post_list_index.size() - 1);
              for (unsigned int id =
                     dof_info.cell_loop_post_list_index[range_index];
                   id != dof_info.cell_loop_post_list_index[range_index + 1];
                   ++id)
                operation_after_loop(dof_info.cell_loop_post_list[id].first,
                                     dof_info.cell_loop_post_list[id].second);
            }
        }
    }

  private:
    const MF     &matrix_free;
    Container    &container;
    function_type cell_function;
    function_type face_function;
    function_type boundary_function;

    const InVector &src;
    OutVector      &dst;
    VectorDataExchange<MF::dimension,
                       typename MF::value_type,
                       typename MF::vectorized_value_type>
      src_data_exchanger;
    VectorDataExchange<MF::dimension,
                       typename MF::value_type,
                       typename MF::vectorized_value_type>
               dst_data_exchanger;
    const bool src_and_dst_are_same;
    const bool zero_dst_vector_setting;
    const std::function<void(const unsigned int, const unsigned int)>
      operation_before_loop;
    const std::function<void(const unsigned int, const unsigned int)>
                       operation_after_loop;
    const unsigned int dof_handler_index_pre_post;
  };



  /**
   * An internal class to convert three function pointers to the
   * scheme with virtual functions above.
   */
  template <class MF, typename InVector, typename OutVector>
  struct MFClassWrapper
  {
    using function_type =
      std::function<void(const MF &,
                         OutVector &,
                         const InVector &,
                         const std::pair<unsigned int, unsigned int> &)>;

    MFClassWrapper(const function_type cell,
                   const function_type face,
                   const function_type boundary)
      : cell(cell)
      , face(face)
      , boundary(boundary)
    {}

    void
    cell_integrator(const MF                                    &mf,
                    OutVector                                   &dst,
                    const InVector                              &src,
                    const std::pair<unsigned int, unsigned int> &range) const
    {
      if (cell)
        cell(mf, dst, src, range);
    }

    void
    face_integrator(const MF                                    &mf,
                    OutVector                                   &dst,
                    const InVector                              &src,
                    const std::pair<unsigned int, unsigned int> &range) const
    {
      if (face)
        face(mf, dst, src, range);
    }

    void
    boundary_integrator(
      const MF                                    &mf,
      OutVector                                   &dst,
      const InVector                              &src,
      const std::pair<unsigned int, unsigned int> &range) const
    {
      if (boundary)
        boundary(mf, dst, src, range);
    }

    const function_type cell;
    const function_type face;
    const function_type boundary;
  };

} // end of namespace internal



template <int dim, typename Number, typename VectorizedArrayType>
template <typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::cell_loop(
  const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                           OutVector &,
                           const InVector &,
                           const std::pair<unsigned int, unsigned int> &)>
                 &cell_operation,
  OutVector      &dst,
  const InVector &src,
  const bool      zero_dst_vector) const
{
  using Wrapper =
    internal::MFClassWrapper<MatrixFree<dim, Number, VectorizedArrayType>,
                             InVector,
                             OutVector>;
  Wrapper wrap(cell_operation, nullptr, nullptr);
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     Wrapper,
                     true>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           wrap,
           &Wrapper::cell_integrator,
           &Wrapper::face_integrator,
           &Wrapper::boundary_integrator);

  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::cell_loop(
  const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                           OutVector &,
                           const InVector &,
                           const std::pair<unsigned int, unsigned int> &)>
                 &cell_operation,
  OutVector      &dst,
  const InVector &src,
  const std::function<void(const unsigned int, const unsigned int)>
    &operation_before_loop,
  const std::function<void(const unsigned int, const unsigned int)>
                    &operation_after_loop,
  const unsigned int dof_handler_index_pre_post) const
{
  using Wrapper =
    internal::MFClassWrapper<MatrixFree<dim, Number, VectorizedArrayType>,
                             InVector,
                             OutVector>;
  Wrapper wrap(cell_operation, nullptr, nullptr);
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     Wrapper,
                     true>
    worker(*this,
           src,
           dst,
           false,
           wrap,
           &Wrapper::cell_integrator,
           &Wrapper::face_integrator,
           &Wrapper::boundary_integrator,
           DataAccessOnFaces::none,
           DataAccessOnFaces::none,
           operation_before_loop,
           operation_after_loop,
           dof_handler_index_pre_post);

  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::loop(
  const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                           OutVector &,
                           const InVector &,
                           const std::pair<unsigned int, unsigned int> &)>
    &cell_operation,
  const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                           OutVector &,
                           const InVector &,
                           const std::pair<unsigned int, unsigned int> &)>
    &inner_face_operation,
  const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                           OutVector &,
                           const InVector &,
                           const std::pair<unsigned int, unsigned int> &)>
                         &boundary_face_operation,
  OutVector              &dst,
  const InVector         &src,
  const bool              zero_dst_vector,
  const DataAccessOnFaces dst_vector_face_access,
  const DataAccessOnFaces src_vector_face_access) const
{
  using Wrapper =
    internal::MFClassWrapper<MatrixFree<dim, Number, VectorizedArrayType>,
                             InVector,
                             OutVector>;
  Wrapper wrap(cell_operation, inner_face_operation, boundary_face_operation);
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     Wrapper,
                     true>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           wrap,
           &Wrapper::cell_integrator,
           &Wrapper::face_integrator,
           &Wrapper::boundary_integrator,
           src_vector_face_access,
           dst_vector_face_access);

  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::cell_loop(
  void (CLASS::*function_pointer)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &) const,
  const CLASS    *owning_class,
  OutVector      &dst,
  const InVector &src,
  const bool      zero_dst_vector) const
{
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     true>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           *owning_class,
           function_pointer,
           nullptr,
           nullptr);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::cell_loop(
  void (CLASS::*function_pointer)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &) const,
  const CLASS    *owning_class,
  OutVector      &dst,
  const InVector &src,
  const std::function<void(const unsigned int, const unsigned int)>
    &operation_before_loop,
  const std::function<void(const unsigned int, const unsigned int)>
                    &operation_after_loop,
  const unsigned int dof_handler_index_pre_post) const
{
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     true>
    worker(*this,
           src,
           dst,
           false,
           *owning_class,
           function_pointer,
           nullptr,
           nullptr,
           DataAccessOnFaces::none,
           DataAccessOnFaces::none,
           operation_before_loop,
           operation_after_loop,
           dof_handler_index_pre_post);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::loop(
  void (CLASS::*cell_operation)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &) const,
  void (CLASS::*inner_face_operation)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &) const,
  void (CLASS::*boundary_face_operation)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &) const,
  const CLASS            *owning_class,
  OutVector              &dst,
  const InVector         &src,
  const bool              zero_dst_vector,
  const DataAccessOnFaces dst_vector_face_access,
  const DataAccessOnFaces src_vector_face_access) const
{
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     true>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           *owning_class,
           cell_operation,
           inner_face_operation,
           boundary_face_operation,
           src_vector_face_access,
           dst_vector_face_access);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::cell_loop(
  void (CLASS::*function_pointer)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &),
  CLASS          *owning_class,
  OutVector      &dst,
  const InVector &src,
  const bool      zero_dst_vector) const
{
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     false>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           *owning_class,
           function_pointer,
           nullptr,
           nullptr);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::cell_loop(
  void (CLASS::*function_pointer)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &),
  CLASS          *owning_class,
  OutVector      &dst,
  const InVector &src,
  const std::function<void(const unsigned int, const unsigned int)>
    &operation_before_loop,
  const std::function<void(const unsigned int, const unsigned int)>
                    &operation_after_loop,
  const unsigned int dof_handler_index_pre_post) const
{
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     false>
    worker(*this,
           src,
           dst,
           false,
           *owning_class,
           function_pointer,
           nullptr,
           nullptr,
           DataAccessOnFaces::none,
           DataAccessOnFaces::none,
           operation_before_loop,
           operation_after_loop,
           dof_handler_index_pre_post);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::loop(
  void (CLASS::*cell_operation)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &),
  void (CLASS::*inner_face_operation)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &),
  void (CLASS::*boundary_face_operation)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &),
  CLASS                  *owning_class,
  OutVector              &dst,
  const InVector         &src,
  const bool              zero_dst_vector,
  const DataAccessOnFaces dst_vector_face_access,
  const DataAccessOnFaces src_vector_face_access) const
{
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     false>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           *owning_class,
           cell_operation,
           inner_face_operation,
           boundary_face_operation,
           src_vector_face_access,
           dst_vector_face_access);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::loop(
  const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                           OutVector &,
                           const InVector &,
                           const std::pair<unsigned int, unsigned int> &)>
    &cell_operation,
  const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                           OutVector &,
                           const InVector &,
                           const std::pair<unsigned int, unsigned int> &)>
    &inner_face_operation,
  const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                           OutVector &,
                           const InVector &,
                           const std::pair<unsigned int, unsigned int> &)>
                 &boundary_face_operation,
  OutVector      &dst,
  const InVector &src,
  const std::function<void(const unsigned int, const unsigned int)>
    &operation_before_loop,
  const std::function<void(const unsigned int, const unsigned int)>
                         &operation_after_loop,
  const unsigned int      dof_handler_index_pre_post,
  const DataAccessOnFaces dst_vector_face_access,
  const DataAccessOnFaces src_vector_face_access) const
{
  using Wrapper =
    internal::MFClassWrapper<MatrixFree<dim, Number, VectorizedArrayType>,
                             InVector,
                             OutVector>;
  Wrapper wrap(cell_operation, inner_face_operation, boundary_face_operation);
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     Wrapper,
                     true>
    worker(*this,
           src,
           dst,
           false,
           wrap,
           &Wrapper::cell_integrator,
           &Wrapper::face_integrator,
           &Wrapper::boundary_integrator,
           src_vector_face_access,
           dst_vector_face_access,
           operation_before_loop,
           operation_after_loop,
           dof_handler_index_pre_post);

  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::loop(
  void (CLASS::*cell_operation)(const MatrixFree &,
                                OutVector &,
                                const InVector &,
                                const std::pair<unsigned int, unsigned int> &)
    const,
  void (CLASS::*inner_face_operation)(
    const MatrixFree &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &) const,
  void (CLASS::*boundary_face_operation)(
    const MatrixFree &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &) const,
  const CLASS    *owning_class,
  OutVector      &dst,
  const InVector &src,
  const std::function<void(const unsigned int, const unsigned int)>
    &operation_before_loop,
  const std::function<void(const unsigned int, const unsigned int)>
                         &operation_after_loop,
  const unsigned int      dof_handler_index_pre_post,
  const DataAccessOnFaces dst_vector_face_access,
  const DataAccessOnFaces src_vector_face_access) const
{
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     true>
    worker(*this,
           src,
           dst,
           false,
           *owning_class,
           cell_operation,
           inner_face_operation,
           boundary_face_operation,
           src_vector_face_access,
           dst_vector_face_access,
           operation_before_loop,
           operation_after_loop,
           dof_handler_index_pre_post);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::loop(
  void (CLASS::*cell_operation)(const MatrixFree &,
                                OutVector &,
                                const InVector &,
                                const std::pair<unsigned int, unsigned int> &),
  void (CLASS::*inner_face_operation)(
    const MatrixFree &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &),
  void (CLASS::*boundary_face_operation)(
    const MatrixFree &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &),
  const CLASS    *owning_class,
  OutVector      &dst,
  const InVector &src,
  const std::function<void(const unsigned int, const unsigned int)>
    &operation_before_loop,
  const std::function<void(const unsigned int, const unsigned int)>
                         &operation_after_loop,
  const unsigned int      dof_handler_index_pre_post,
  const DataAccessOnFaces dst_vector_face_access,
  const DataAccessOnFaces src_vector_face_access) const
{
  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     false>
    worker(*this,
           src,
           dst,
           false,
           *owning_class,
           cell_operation,
           inner_face_operation,
           boundary_face_operation,
           src_vector_face_access,
           dst_vector_face_access,
           operation_before_loop,
           operation_after_loop,
           dof_handler_index_pre_post);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::loop_cell_centric(
  void (CLASS::*function_pointer)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &) const,
  const CLASS            *owning_class,
  OutVector              &dst,
  const InVector         &src,
  const bool              zero_dst_vector,
  const DataAccessOnFaces src_vector_face_access) const
{
  auto src_vector_face_access_temp = src_vector_face_access;
  if (DataAccessOnFaces::gradients == src_vector_face_access_temp)
    src_vector_face_access_temp = DataAccessOnFaces::gradients_all_faces;
  else if (DataAccessOnFaces::values == src_vector_face_access_temp)
    src_vector_face_access_temp = DataAccessOnFaces::values_all_faces;

  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     true>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           *owning_class,
           function_pointer,
           nullptr,
           nullptr,
           src_vector_face_access_temp,
           DataAccessOnFaces::none);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename CLASS, typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::loop_cell_centric(
  void (CLASS::*function_pointer)(
    const MatrixFree<dim, Number, VectorizedArrayType> &,
    OutVector &,
    const InVector &,
    const std::pair<unsigned int, unsigned int> &),
  CLASS                  *owning_class,
  OutVector              &dst,
  const InVector         &src,
  const bool              zero_dst_vector,
  const DataAccessOnFaces src_vector_face_access) const
{
  auto src_vector_face_access_temp = src_vector_face_access;
  if (DataAccessOnFaces::gradients == src_vector_face_access_temp)
    src_vector_face_access_temp = DataAccessOnFaces::gradients_all_faces;
  else if (DataAccessOnFaces::values == src_vector_face_access_temp)
    src_vector_face_access_temp = DataAccessOnFaces::values_all_faces;

  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     CLASS,
                     false>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           *owning_class,
           function_pointer,
           nullptr,
           nullptr,
           src_vector_face_access_temp,
           DataAccessOnFaces::none);
  task_info.loop(worker);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename OutVector, typename InVector>
inline void
MatrixFree<dim, Number, VectorizedArrayType>::loop_cell_centric(
  const std::function<void(const MatrixFree<dim, Number, VectorizedArrayType> &,
                           OutVector &,
                           const InVector &,
                           const std::pair<unsigned int, unsigned int> &)>
                         &cell_operation,
  OutVector              &dst,
  const InVector         &src,
  const bool              zero_dst_vector,
  const DataAccessOnFaces src_vector_face_access) const
{
  auto src_vector_face_access_temp = src_vector_face_access;
  if (DataAccessOnFaces::gradients == src_vector_face_access_temp)
    src_vector_face_access_temp = DataAccessOnFaces::gradients_all_faces;
  else if (DataAccessOnFaces::values == src_vector_face_access_temp)
    src_vector_face_access_temp = DataAccessOnFaces::values_all_faces;

  using Wrapper =
    internal::MFClassWrapper<MatrixFree<dim, Number, VectorizedArrayType>,
                             InVector,
                             OutVector>;
  Wrapper wrap(cell_operation, nullptr, nullptr);

  internal::MFWorker<MatrixFree<dim, Number, VectorizedArrayType>,
                     InVector,
                     OutVector,
                     Wrapper,
                     true>
    worker(*this,
           src,
           dst,
           zero_dst_vector,
           wrap,
           &Wrapper::cell_integrator,
           &Wrapper::face_integrator,
           &Wrapper::boundary_integrator,
           src_vector_face_access_temp,
           DataAccessOnFaces::none);
  task_info.loop(worker);
}


#endif // ifndef DOXYGEN



DEAL_II_NAMESPACE_CLOSE

#endif
