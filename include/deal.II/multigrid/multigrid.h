// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2017 by the deal.II authors
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

#ifndef dealii_multigrid_h
#define dealii_multigrid_h


#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/multigrid/mg_base.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/distributed/tria.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mg */
/*@{*/

/**
 * Implementation of the multigrid method. The implementation supports both
 * continuous and discontinuous elements and follows the procedure described in
 * the @ref mg_paper "multigrid paper by Janssen and Kanschat".
 *
 * The function which starts a multigrid cycle on the finest level is cycle().
 * Depending on the cycle type chosen with the constructor (see enum Cycle),
 * this function triggers one of the cycles level_v_step() or level_step(),
 * where the latter one can do different types of cycles.
 *
 * Using this class, it is expected that the right hand side has been
 * converted from a vector living on the locally finest level to a multilevel
 * vector. This is a nontrivial operation, usually initiated automatically by
 * the class PreconditionMG and performed by the classes derived from
 * MGTransferBase.
 *
 * @note The interface of this class is still very clumsy. In particular, you
 * will have to set up quite a few auxiliary objects before you can use it.
 * Unfortunately, it seems that this can be avoided only be restricting the
 * flexibility of this class in an unacceptable way.
 *
 * @author Guido Kanschat, 1999 - 2005
 */
template <typename VectorType>
class Multigrid : public Subscriptor
{
public:
  /**
   * List of implemented cycle types.
   */
  enum Cycle
  {
    /// The V-cycle
    v_cycle,
    /// The W-cycle
    w_cycle,
    /// The F-cycle
    f_cycle
  };

  /**
   * A structure that has boost::signal objects for a number of actions that
   * happen for the typical usage of this class.
   *
   * For documentation on signals, see
   * http://www.boost.org/doc/libs/release/libs/signals2 .
   */
  struct Signals
  {
    /**
     * This signal is triggered before (@p start is true) and after (@p start is
     * false) the call to MGTransfer::copy_to_mg which transfers the vector
     * given to it to a multi-level vector.
     */
    boost::signals2::signal<void (const bool start)> transfer_to_mg;

    /**
     * This signal is triggered before (@p start is true) and after (@p start is
     * false) the call to MGTransfer::copy_from_mg which transfers the
     * multi-level vector given to it to a normal vector.
     */
    boost::signals2::signal<void (const bool start)> transfer_to_global;

    /**
     * This signal is triggered before (@p start is true) and after (@p start is
     * false) the call to the coarse solver on @p level.
     */
    boost::signals2::signal<void (const bool start, const unsigned int level)> coarse_solve;

    /**
     * This signal is triggered before (@p start is true) and after (@p start is
     * false) the call to MGTransfer::restrict_and_add() which restricts a
     * vector from @p level to the next coarser one (@p level - 1).
     */
    boost::signals2::signal<void (const bool start, const unsigned int level)> transfer_down_to;

    /**
     * This signal is triggered before (@p start is true) and after (@p start is
     * false) the call to MGTransfer::prolongate() which prolongs a vector to
     * @p level from the next coarser one (@p level - 1).
     */
    boost::signals2::signal<void (const bool start, const unsigned int level)> transfer_up_to;

    /**
     * This signal is triggered before (@p start is true) and after (@p start is
     * false) the call to a pre-smoothing step via MGPreSmoother::apply() on
     * @p level.
     */
    boost::signals2::signal<void (const bool start, const unsigned int level)> pre_smoother_step;

    /**
     * This signal is triggered before (@p start is true) and after (@p start is
     * false) the call to a post-smoothing step via MGPostSmoother::apply() on
     * @p level.
     */
    boost::signals2::signal<void (const bool start, const unsigned int level)> post_smoother_step;
  };

  typedef VectorType vector_type;
  typedef const VectorType const_vector_type;

  /**
   * Constructor. The DoFHandler is used to check whether the provided
   * minlevel and maxlevel are in the range of valid levels.
   * If maxlevel is set to the default value, the highest valid
   * level is used.
   * <tt>transfer</tt> is an object performing prolongation and
   * restriction.
   *
   * This function already initializes the vectors which will be used later in
   * the course of the computations. You should therefore create objects of
   * this type as late as possible.
   *
   * @deprecated Use the other constructor instead.
   * The DoFHandler is actually not needed.
   */
  template <int dim>
#ifndef _MSC_VER
  DEAL_II_DEPRECATED
#endif
  Multigrid(const DoFHandler<dim>              &mg_dof_handler,
            const MGMatrixBase<VectorType>     &matrix,
            const MGCoarseGridBase<VectorType> &coarse,
            const MGTransferBase<VectorType>   &transfer,
            const MGSmootherBase<VectorType>   &pre_smooth,
            const MGSmootherBase<VectorType>   &post_smooth,
            const unsigned int                 minlevel = 0,
            const unsigned int                 maxlevel = numbers::invalid_unsigned_int,
            Cycle                              cycle = v_cycle);

  /**
   * Constructor. <tt>transfer</tt> is an object performing prolongation and
   * restriction. For levels in [minlevel, maxlevel] matrix has to contain
   * valid matrices. By default the maxlevel is set to the maximal valid level.
   *
   * This function already initializes the vectors which will be used later in
   * the course of the computations. You should therefore create objects of
   * this type as late as possible.
   */
  Multigrid(const MGMatrixBase<VectorType>     &matrix,
            const MGCoarseGridBase<VectorType> &coarse,
            const MGTransferBase<VectorType>   &transfer,
            const MGSmootherBase<VectorType>   &pre_smooth,
            const MGSmootherBase<VectorType>   &post_smooth,
            const unsigned int                 minlevel = 0,
            const unsigned int                 maxlevel = numbers::invalid_unsigned_int,
            Cycle                              cycle = v_cycle);

  /**
   * Reinit this class according to #minlevel and #maxlevel.
   */
  void reinit (const unsigned int minlevel,
               const unsigned int maxlevel);

  /**
   * Execute one multigrid cycle. The type of cycle is selected by the
   * constructor argument cycle. See the enum Cycle for available types.
   */
  void cycle ();

  /**
   * Execute one step of the V-cycle algorithm.  This function assumes, that
   * the multilevel vector #defect is filled with the residual of an outer
   * defect correction scheme. This is usually taken care of by
   * PreconditionMG). After vcycle(), the result is in the multilevel vector
   * #solution. See <tt>copy_*_mg</tt> in the MGTools namespace if you want to
   * use these vectors yourself.
   *
   * The actual work for this function is done in level_v_step().
   */
  void vcycle ();

  /**
   * Set additional matrices to correct residual computation at refinement
   * edges. Since we only smoothen in the interior of the refined part of the
   * mesh, the coupling across the refinement edge is missing. This coupling
   * is provided by these two matrices.
   *
   * @note While <tt>edge_out.vmult</tt> is used, for the second argument, we
   * use <tt>edge_in.Tvmult</tt>. Thus, <tt>edge_in</tt> should be assembled
   * in transposed form. This saves a second sparsity pattern for
   * <tt>edge_in</tt>. In particular, for symmetric operators, both arguments
   * can refer to the same matrix, saving assembling of one of them.
   */
  void set_edge_matrices (const MGMatrixBase<VectorType> &edge_out,
                          const MGMatrixBase<VectorType> &edge_in);

  /**
   * Set additional matrices to correct residual computation at refinement
   * edges. These matrices originate from discontinuous Galerkin methods (see
   * FE_DGQ etc.), where they correspond to the edge fluxes at the refinement
   * edge between two levels.
   *
   * @note While <tt>edge_down.vmult</tt> is used, for the second argument, we
   * use <tt>edge_up.Tvmult</tt>. Thus, <tt>edge_up</tt> should be assembled
   * in transposed form. This saves a second sparsity pattern for
   * <tt>edge_up</tt>. In particular, for symmetric operators, both arguments
   * can refer to the same matrix, saving assembling of one of them.
   */
  void set_edge_flux_matrices (const MGMatrixBase<VectorType> &edge_down,
                               const MGMatrixBase<VectorType> &edge_up);

  /**
   * Return the finest level for multigrid.
   */
  unsigned int get_maxlevel() const;

  /**
   * Return the coarsest level for multigrid.
   */
  unsigned int get_minlevel() const;

  /**
   * Set the highest level for which the multilevel method is performed. By
   * default, this is the finest level of the Triangulation. Accepted are
   * values not smaller than the current #minlevel.
   */
  void set_maxlevel (const unsigned int);

  /**
   * Set the coarsest level for which the multilevel method is performed. By
   * default, this is zero. Accepted are non-negative values not larger than
   * the current #maxlevel.
   *
   * If <tt>relative</tt> is <tt>true</tt>, then this function determines the
   * number of levels used, that is, it sets #minlevel to
   * #maxlevel-<tt>level</tt>.
   *
   * @note The mesh on the coarsest level must cover the whole domain. There
   * may not be hanging nodes on #minlevel.
   *
   * @note If #minlevel is set to a nonzero value, do not forget to adjust
   * your coarse grid solver!
   */
  void set_minlevel (const unsigned int level,
                     bool relative = false);

  /**
   * Chance #cycle_type used in cycle().
   */
  void set_cycle(Cycle);

  /**
   * Set the debug level. Higher values will create more debugging output
   * during the multigrid cycles.
   */
  void set_debug (const unsigned int);

  /**
   * Connect a function to Signals::coarse_solve.
   */
  boost::signals2::connection
  connect_coarse_solve(const std::function<void (const bool, const unsigned int)> &slot);

  /**
   * Connect a function to Signals::transfer_down_to.
   */
  boost::signals2::connection
  connect_transfer_down_to(const std::function<void (const bool, const unsigned int)> &slot);

  /**
   * Connect a function to Signals::transfer_up_to.
   */
  boost::signals2::connection
  connect_transfer_up_to(const std::function<void (const bool, const unsigned int)> &slot);

  /**
   * Connect a function to Signals::pre_smoother_step.
   */
  boost::signals2::connection
  connect_pre_smoother_step(const std::function<void (const bool, const unsigned int)> &slot);

  /**
   * Connect a function to Signals::post_smoother_step.
   */
  boost::signals2::connection
  connect_post_smoother_step(const std::function<void (const bool, const unsigned int)> &slot);

private:

  /**
   * Signals for the various actions that the Multigrid algorithm uses.
    */
  Signals signals;

  /**
   * The V-cycle multigrid method. <tt>level</tt> is the level the function
   * starts on. It will usually be called for the highest level from outside,
   * but will then call itself recursively for <tt>level-1</tt>, unless we are
   * on #minlevel where the coarse grid solver solves the problem exactly.
   */
  void level_v_step (const unsigned int level);

  /**
   * The actual W-cycle or F-cycle multigrid method. <tt>level</tt> is the
   * level the function starts on. It will usually be called for the highest
   * level from outside, but will then call itself recursively for
   * <tt>level-1</tt>, unless we are on #minlevel where the coarse grid solver
   * solves the problem exactly.
   */
  void level_step (const unsigned int level, Cycle cycle);

  /**
   * Cycle type performed by the method cycle().
   */
  Cycle cycle_type;

  /**
   * Level for coarse grid solution.
   */
  unsigned int minlevel;

  /**
   * Highest level of cells.
   */
  unsigned int maxlevel;

public:
  /**
   * Input vector for the cycle. Contains the defect of the outer method
   * projected to the multilevel vectors.
   */
  MGLevelObject<VectorType> defect;

  /**
   * The solution update after the multigrid step.
   */
  MGLevelObject<VectorType> solution;

private:
  /**
   * Auxiliary vector.
   */
  MGLevelObject<VectorType> t;

  /**
   * Auxiliary vector for W- and F-cycles. Left uninitialized in V-cycle.
   */
  MGLevelObject<VectorType> defect2;


  /**
   * The matrix for each level.
   */
  SmartPointer<const MGMatrixBase<VectorType>,Multigrid<VectorType> > matrix;

  /**
   * The matrix for each level.
   */
  SmartPointer<const MGCoarseGridBase<VectorType>,Multigrid<VectorType> > coarse;

  /**
   * Object for grid transfer.
   */
  SmartPointer<const MGTransferBase<VectorType>,Multigrid<VectorType> > transfer;

  /**
   * The pre-smoothing object.
   */
  SmartPointer<const MGSmootherBase<VectorType>,Multigrid<VectorType> > pre_smooth;

  /**
   * The post-smoothing object.
   */
  SmartPointer<const MGSmootherBase<VectorType>,Multigrid<VectorType> > post_smooth;

  /**
   * Edge matrix from the interior of the refined part to the refinement edge.
   *
   * @note Only <tt>vmult</tt> is used for these matrices.
   */
  SmartPointer<const MGMatrixBase<VectorType> > edge_out;

  /**
   * Transpose edge matrix from the refinement edge to the interior of the
   * refined part.
   *
   * @note Only <tt>Tvmult</tt> is used for these matrices.
   */
  SmartPointer<const MGMatrixBase<VectorType> > edge_in;

  /**
   * Edge matrix from fine to coarse.
   *
   * @note Only <tt>vmult</tt> is used for these matrices.
   */
  SmartPointer<const MGMatrixBase<VectorType>,Multigrid<VectorType> > edge_down;

  /**
   * Transpose edge matrix from coarse to fine.
   *
   * @note Only <tt>Tvmult</tt> is used for these matrices.
   */
  SmartPointer<const MGMatrixBase<VectorType>,Multigrid<VectorType> > edge_up;

  /**
   * Level for debug output. Defaults to zero and can be set by set_debug().
   */
  unsigned int debug;

  template <int dim, class OtherVectorType, class TRANSFER> friend class PreconditionMG;
};


/**
 * Multi-level preconditioner. Here, we collect all information needed for
 * multi-level preconditioning and provide the standard interface for LAC
 * iterative methods.
 *
 * Furthermore, it needs functions <tt>void copy_to_mg(const VectorType&)</tt>
 * to store @p src in the right hand side of the multi-level method and
 * <tt>void copy_from_mg(VectorType&)</tt> to store the result of the v-cycle
 * in @p dst.
 *
 * If VectorType is in fact a block vector and the TRANSFER object supports
 * use of a separate DoFHandler for each block, this class also allows
 * to be initialized with a separate DoFHandler for each block.
 *
 * @author Guido Kanschat, Daniel Arndt, 1999, 2000, 2001, 2002, 2017
 */
template <int dim, typename VectorType, class TRANSFER>
class PreconditionMG : public Subscriptor
{
public:
  /**
   * Constructor. Arguments are the multigrid object, pre-smoother, post-
   * smoother and coarse grid solver.
   */
  PreconditionMG(const DoFHandler<dim> &dof_handler,
                 Multigrid<VectorType> &mg,
                 const TRANSFER        &transfer);

  /**
   * Same as above in case every component of a block vector
   * uses its own DoFHandler.
   */
  PreconditionMG(const std::vector<const DoFHandler<dim>*> &dof_handler,
                 Multigrid<VectorType> &mg,
                 const TRANSFER        &transfer);

  /**
   * Dummy function needed by other classes.
   */
  bool empty () const;

  /**
   * Preconditioning operator. Calls the @p vcycle function of the @p MG
   * object passed to the constructor.
   *
   * This is the operator used by LAC iterative solvers.
   */
  template <class OtherVectorType>
  void vmult (OtherVectorType       &dst,
              const OtherVectorType &src) const;

  /**
   * Preconditioning operator. Calls the @p vcycle function of the @p MG
   * object passed to the constructor.
   */
  template <class OtherVectorType>
  void vmult_add (OtherVectorType       &dst,
                  const OtherVectorType &src) const;

  /**
   * Transposed preconditioning operator.
   *
   * Not implemented, but the definition may be needed.
   */
  template <class OtherVectorType>
  void Tvmult (OtherVectorType       &dst,
               const OtherVectorType &src) const;

  /**
   * Transposed preconditioning operator.
   *
   * Not implemented, but the definition may be needed.
   */
  template <class OtherVectorType>
  void Tvmult_add (OtherVectorType       &dst,
                   const OtherVectorType &src) const;

  /**
   * Return the partitioning of the range space of this preconditioner, i.e.,
   * the partitioning of the vectors that are result from matrix-vector products.
   * By default, the respective information for the first DoFHandler object
   * are returned.
   */
  IndexSet locally_owned_range_indices(const unsigned int block=0) const;

  /**
   * Return the partitioning of the domain space of this preconditioner, i.e.,
   * the partitioning of the vectors this matrix has to be multiplied with.
   * By default, the respective information for the first DoFHandler object
   * are returned.
   */
  IndexSet locally_owned_domain_indices(const unsigned int block=0) const;

  /**
   * Return the MPI communicator object in use with this preconditioner.
   */
  MPI_Comm get_mpi_communicator() const;

  /**
   * Connect a function to Signals::transfer_to_mg.
   */
  boost::signals2::connection
  connect_transfer_to_mg(const std::function<void (bool)> &slot);

  /**
   * Connect a function to Signals::transfer_to_global.
   */
  boost::signals2::connection
  connect_transfer_to_global(const std::function<void (bool)> &slot);

private:
  /**
   * Associated @p DoFHandler.
   */
  std::vector<SmartPointer<const DoFHandler<dim>,PreconditionMG<dim,VectorType,TRANSFER> > > dof_handler_vector;

  /**
   * Storage for the pointers to the DoFHandler objects
   * without SmartPointer wrapper.
   */
  std::vector<const DoFHandler<dim>*> dof_handler_vector_raw;

  /**
   * The multigrid object.
   */
  SmartPointer<Multigrid<VectorType>,PreconditionMG<dim,VectorType,TRANSFER> > multigrid;

  /**
   * Object for grid transfer.
   */
  SmartPointer<const TRANSFER,PreconditionMG<dim,VectorType,TRANSFER> > transfer;

  /**
   * Flag to indicate if the object is initialized with a single DoFHandler
   * or with one for each block.
   */
  const bool uses_dof_handler_vector;

  /**
  * Signals used by this object
  */
  typename Multigrid<VectorType>::Signals signals;
};

/*@}*/

#ifndef DOXYGEN
/* --------------------------- inline functions --------------------- */


template <typename VectorType>
template <int dim>
Multigrid<VectorType>::Multigrid (const DoFHandler<dim>              &mg_dof_handler,
                                  const MGMatrixBase<VectorType>     &matrix,
                                  const MGCoarseGridBase<VectorType> &coarse,
                                  const MGTransferBase<VectorType>   &transfer,
                                  const MGSmootherBase<VectorType>   &pre_smooth,
                                  const MGSmootherBase<VectorType>   &post_smooth,
                                  const unsigned int                 min_level,
                                  const unsigned int                 max_level,
                                  Cycle                              cycle)
  :
  cycle_type(cycle),
  minlevel(min_level),
  matrix(&matrix, typeid(*this).name()),
  coarse(&coarse, typeid(*this).name()),
  transfer(&transfer, typeid(*this).name()),
  pre_smooth(&pre_smooth, typeid(*this).name()),
  post_smooth(&post_smooth, typeid(*this).name()),
  edge_down(nullptr, typeid(*this).name()),
  edge_up(nullptr, typeid(*this).name()),
  debug(0)
{
  const unsigned int dof_handler_max_level
    = mg_dof_handler.get_triangulation().n_global_levels()-1;
  if (max_level == numbers::invalid_unsigned_int)
    maxlevel = dof_handler_max_level;
  else
    maxlevel = max_level;

  reinit(minlevel, maxlevel);
}



template <typename VectorType>
Multigrid<VectorType>::Multigrid (const MGMatrixBase<VectorType>        &matrix,
                                  const MGCoarseGridBase<VectorType>    &coarse,
                                  const MGTransferBase<VectorType>      &transfer,
                                  const MGSmootherBase<VectorType>      &pre_smooth,
                                  const MGSmootherBase<VectorType>      &post_smooth,
                                  const unsigned int                    min_level,
                                  const unsigned int                    max_level,
                                  Cycle                                 cycle)
  :
  cycle_type(cycle),
  matrix(&matrix, typeid(*this).name()),
  coarse(&coarse, typeid(*this).name()),
  transfer(&transfer, typeid(*this).name()),
  pre_smooth(&pre_smooth, typeid(*this).name()),
  post_smooth(&post_smooth, typeid(*this).name()),
  edge_out(nullptr, typeid(*this).name()),
  edge_in(nullptr, typeid(*this).name()),
  edge_down(nullptr, typeid(*this).name()),
  edge_up(nullptr, typeid(*this).name()),
  debug(0)
{
  if (max_level == numbers::invalid_unsigned_int)
    maxlevel = matrix.get_maxlevel();
  else
    maxlevel = max_level;
  reinit(min_level, maxlevel);
}



template <typename VectorType>
inline
unsigned int
Multigrid<VectorType>::get_maxlevel () const
{
  return maxlevel;
}



template <typename VectorType>
inline
unsigned int
Multigrid<VectorType>::get_minlevel () const
{
  return minlevel;
}


/* --------------------------- inline functions --------------------- */


namespace internal
{
  namespace PreconditionMGImplementation
  {
    template <int dim, typename VectorType, class TRANSFER, typename OtherVectorType>
    typename std::enable_if<TRANSFER::supports_dof_handler_vector>::type
    vmult(const std::vector<const dealii::DoFHandler<dim>*> &dof_handler_vector,
          dealii::Multigrid<VectorType>  &multigrid,
          const TRANSFER &transfer,
          OtherVectorType       &dst,
          const OtherVectorType &src,
          const bool uses_dof_handler_vector,
          const typename dealii::Multigrid<VectorType>::Signals &signals, int)
    {
      signals.transfer_to_mg(true);
      if (uses_dof_handler_vector)
        transfer.copy_to_mg(dof_handler_vector,
                            multigrid.defect,
                            src);
      else
        transfer.copy_to_mg(*dof_handler_vector[0],
                            multigrid.defect,
                            src);
      signals.transfer_to_mg(false);

      multigrid.cycle();

      signals.transfer_to_global(true);
      if (uses_dof_handler_vector)
        transfer.copy_from_mg(dof_handler_vector,
                              dst,
                              multigrid.solution);
      else
        transfer.copy_from_mg(*dof_handler_vector[0],
                              dst,
                              multigrid.solution);
      signals.transfer_to_global(false);
    }

    template <int dim, typename VectorType, class TRANSFER, typename OtherVectorType>
    void
    vmult(const std::vector<const dealii::DoFHandler<dim>*> &dof_handler_vector,
          dealii::Multigrid<VectorType>  &multigrid,
          const TRANSFER &transfer,
          OtherVectorType       &dst,
          const OtherVectorType &src,
          const bool uses_dof_handler_vector,
          const typename dealii::Multigrid<VectorType>::Signals &signals, ...)
    {
      (void) uses_dof_handler_vector;
      Assert (!uses_dof_handler_vector, ExcInternalError());

      signals.transfer_to_mg(true);
      transfer.copy_to_mg(*dof_handler_vector[0],
                          multigrid.defect,
                          src);
      signals.transfer_to_mg(false);

      multigrid.cycle();

      signals.transfer_to_global(true);
      transfer.copy_from_mg(*dof_handler_vector[0],
                            dst,
                            multigrid.solution);
      signals.transfer_to_global(false);
    }

    template <int dim, typename VectorType, class TRANSFER, typename OtherVectorType>
    typename std::enable_if<TRANSFER::supports_dof_handler_vector>::type
    vmult_add(const std::vector<const dealii::DoFHandler<dim>*> &dof_handler_vector,
              dealii::Multigrid<VectorType>  &multigrid,
              const TRANSFER &transfer,
              OtherVectorType       &dst,
              const OtherVectorType &src,
              const bool uses_dof_handler_vector,
              const typename dealii::Multigrid<VectorType>::Signals &signals, int)
    {
      signals.transfer_to_mg(0, true);
      if (uses_dof_handler_vector)
        transfer.copy_to_mg(dof_handler_vector,
                            multigrid.defect,
                            src);
      else
        transfer.copy_to_mg(*dof_handler_vector[0],
                            multigrid.defect,
                            src);
      signals.transfer_to_mg(0, false);

      multigrid.cycle();

      signals.transfer_to_global(0, true);
      if (uses_dof_handler_vector)
        transfer.copy_from_mg_add(dof_handler_vector,
                                  dst,
                                  multigrid.solution);
      else
        transfer.copy_from_mg_add(*dof_handler_vector[0],
                                  dst,
                                  multigrid.solution);
      signals.transfer_to_global(0, false);
    }

    template <int dim, typename VectorType, class TRANSFER, typename OtherVectorType>
    void
    vmult_add(const std::vector<const dealii::DoFHandler<dim>*> &dof_handler_vector,
              dealii::Multigrid<VectorType>  &multigrid,
              const TRANSFER &transfer,
              OtherVectorType       &dst,
              const OtherVectorType &src,
              const bool uses_dof_handler_vector,
              const typename dealii::Multigrid<VectorType>::Signals &signals, ...)
    {
      (void) uses_dof_handler_vector;
      Assert (!uses_dof_handler_vector, ExcInternalError());

      signals.transfer_to_mg(true);
      transfer.copy_to_mg(*dof_handler_vector[0],
                          multigrid.defect,
                          src);
      signals.transfer_to_mg(false);

      multigrid.cycle();

      signals.transfer_to_global(true);
      transfer.copy_from_mg_add(*dof_handler_vector[0],
                                dst,
                                multigrid.solution);
      signals.transfer_to_global(false);
    }
  }
}

template <int dim, typename VectorType, class TRANSFER>
PreconditionMG<dim, VectorType, TRANSFER>
::PreconditionMG(const DoFHandler<dim>  &dof_handler,
                 Multigrid<VectorType>  &mg,
                 const TRANSFER         &transfer)
  :
  dof_handler_vector(1,&dof_handler),
  dof_handler_vector_raw(1,&dof_handler),
  multigrid(&mg),
  transfer(&transfer),
  uses_dof_handler_vector(false)
{}

template <int dim, typename VectorType, class TRANSFER>
PreconditionMG<dim, VectorType, TRANSFER>
::PreconditionMG(const std::vector<const DoFHandler<dim>*>  &dof_handler,
                 Multigrid<VectorType>  &mg,
                 const TRANSFER         &transfer)
  :
  dof_handler_vector(dof_handler.size()),
  dof_handler_vector_raw(dof_handler.size()),
  multigrid(&mg),
  transfer(&transfer),
  uses_dof_handler_vector(true)
{
  for (unsigned int i = 0; i< dof_handler.size() ; ++i)
    {
      dof_handler_vector[i] = dof_handler[i];
      dof_handler_vector_raw[i] = dof_handler[i];
    }
}

template <int dim, typename VectorType, class TRANSFER>
inline bool
PreconditionMG<dim, VectorType, TRANSFER>::empty () const
{
  return false;
}

template <int dim, typename VectorType, class TRANSFER>
template <class OtherVectorType>
void
PreconditionMG<dim, VectorType, TRANSFER>::vmult
(OtherVectorType       &dst,
 const OtherVectorType &src) const
{
  internal::PreconditionMGImplementation::vmult(dof_handler_vector_raw,*multigrid,*transfer,
                                                dst,src,uses_dof_handler_vector,
                                                this->signals, 0);
}


template <int dim, typename VectorType, class TRANSFER>
IndexSet
PreconditionMG<dim, VectorType, TRANSFER>::locally_owned_range_indices(const unsigned int block) const
{
  AssertIndexRange(block,dof_handler_vector.size());
  return dof_handler_vector[block]->locally_owned_dofs();
}


template <int dim, typename VectorType, class TRANSFER>
IndexSet
PreconditionMG<dim, VectorType, TRANSFER>::locally_owned_domain_indices(const unsigned int block) const
{
  AssertIndexRange(block,dof_handler_vector.size());
  return dof_handler_vector[block]->locally_owned_dofs();
}



template <int dim, typename VectorType, class TRANSFER>
MPI_Comm
PreconditionMG<dim, VectorType, TRANSFER>::get_mpi_communicator() const
{
  // currently parallel GMG works with distributed Triangulation only,
  // so it should be a safe bet to use it to query MPI communicator:
  const Triangulation<dim> &tria = dof_handler_vector[0]->get_triangulation();
  const parallel::distributed::Triangulation<dim> *ptria = dynamic_cast<const parallel::distributed::Triangulation<dim> *>(&tria);
  Assert (ptria != nullptr, ExcInternalError());
  return ptria->get_communicator ();
}



template <int dim, typename VectorType, class TRANSFER>
boost::signals2::connection
PreconditionMG<dim, VectorType, TRANSFER>::
connect_transfer_to_mg(const std::function<void (bool)> &slot)
{
  return this->signals.transfer_to_mg.connect(slot);
}



template <int dim, typename VectorType, class TRANSFER>
boost::signals2::connection
PreconditionMG<dim, VectorType, TRANSFER>::
connect_transfer_to_global(
  const std::function<void (bool)> &slot)
{
  return this->signals.transfer_to_global.connect(slot);
}



template <int dim, typename VectorType, class TRANSFER>
template <class OtherVectorType>
void
PreconditionMG<dim, VectorType, TRANSFER>::vmult_add
(OtherVectorType       &dst,
 const OtherVectorType &src) const
{
  internal::PreconditionMGImplementation::vmult_add(dof_handler_vector_raw,*multigrid,*transfer,
                                                    dst,src,uses_dof_handler_vector,
                                                    this->signals, 0);
}


template <int dim, typename VectorType, class TRANSFER>
template <class OtherVectorType>
void
PreconditionMG<dim, VectorType, TRANSFER>::Tvmult
(OtherVectorType &,
 const OtherVectorType &) const
{
  Assert(false, ExcNotImplemented());
}


template <int dim, typename VectorType, class TRANSFER>
template <class OtherVectorType>
void
PreconditionMG<dim, VectorType, TRANSFER>::Tvmult_add
(OtherVectorType &,
 const OtherVectorType &) const
{
  Assert(false, ExcNotImplemented());
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
