// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2016 by the deal.II authors
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

#ifndef dealii__multigrid_h
#define dealii__multigrid_h


#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/multigrid/mg_base.h>
#include <deal.II/base/mg_level_object.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mg */
/*@{*/

/**
 * Implementation of the multigrid method.
 *
 * @warning multigrid on locally refined meshes only works with
 * <b>discontinuous finite elements</b> right now. It is not clear, whether
 * the paradigm of local smoothing we use is applicable to continuous elements
 * with hanging nodes; in fact, most people you meet on conferences seem to
 * deny this.
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

  typedef VectorType vector_type;
  typedef const VectorType const_vector_type;

  /**
   * Constructor. The DoFHandler is used to determine the highest possible
   * level. <tt>transfer</tt> is an object performing prolongation and
   * restriction.
   *
   * This function already initializes the vectors which will be used later in
   * the course of the computations. You should therefore create objects of
   * this type as late as possible.
   */
  template <int dim>
  Multigrid(const DoFHandler<dim>              &mg_dof_handler,
            const MGMatrixBase<VectorType>     &matrix,
            const MGCoarseGridBase<VectorType> &coarse,
            const MGTransferBase<VectorType>   &transfer,
            const MGSmootherBase<VectorType>   &pre_smooth,
            const MGSmootherBase<VectorType>   &post_smooth,
            Cycle                              cycle = v_cycle);

  /**
   * Experimental constructor for cases in which no DoFHandler is available.
   *
   * @warning Not intended for general use.
   */
  Multigrid(const unsigned int                 minlevel,
            const unsigned int                 maxlevel,
            const MGMatrixBase<VectorType>     &matrix,
            const MGCoarseGridBase<VectorType> &coarse,
            const MGTransferBase<VectorType>   &transfer,
            const MGSmootherBase<VectorType>   &pre_smooth,
            const MGSmootherBase<VectorType>   &post_smooth,
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
   * default, this is the finest level of the Triangulation; therefore, this
   * function will only accept arguments smaller than the current #maxlevel
   * and not smaller than the current #minlevel.
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

private:

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
   * Object for grid tranfer.
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

  template<int dim, class OtherVectorType, class TRANSFER> friend class PreconditionMG;
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
 * @author Guido Kanschat, 1999, 2000, 2001, 2002
 */
template<int dim, typename VectorType, class TRANSFER>
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
   * Dummy function needed by other classes.
   */
  bool empty () const;

  /**
   * Preconditioning operator. Calls the @p vcycle function of the @p MG
   * object passed to the constructor.
   *
   * This is the operator used by LAC iterative solvers.
   */
  template<class OtherVectorType>
  void vmult (OtherVectorType       &dst,
              const OtherVectorType &src) const;

  /**
   * Preconditioning operator. Calls the @p vcycle function of the @p MG
   * object passed to the constructor.
   */
  template<class OtherVectorType>
  void vmult_add (OtherVectorType       &dst,
                  const OtherVectorType &src) const;

  /**
   * Tranposed preconditioning operator.
   *
   * Not implemented, but the definition may be needed.
   */
  template<class OtherVectorType>
  void Tvmult (OtherVectorType       &dst,
               const OtherVectorType &src) const;

  /**
   * Tranposed preconditioning operator.
   *
   * Not implemented, but the definition may be needed.
   */
  template<class OtherVectorType>
  void Tvmult_add (OtherVectorType       &dst,
                   const OtherVectorType &src) const;

private:
  /**
   * Associated @p DoFHandler.
   */
  SmartPointer<const DoFHandler<dim>,PreconditionMG<dim,VectorType,TRANSFER> > dof_handler;

  /**
   * The multigrid object.
   */
  SmartPointer<Multigrid<VectorType>,PreconditionMG<dim,VectorType,TRANSFER> > multigrid;

  /**
   * Object for grid tranfer.
   */
  SmartPointer<const TRANSFER,PreconditionMG<dim,VectorType,TRANSFER> > transfer;
};

/*@}*/

#ifndef DOXYGEN
/* --------------------------- inline functions --------------------- */


template <typename VectorType>
template <int dim>
Multigrid<VectorType>::Multigrid (const DoFHandler<dim>          &mg_dof_handler,
                                  const MGMatrixBase<VectorType>     &matrix,
                                  const MGCoarseGridBase<VectorType> &coarse,
                                  const MGTransferBase<VectorType>   &transfer,
                                  const MGSmootherBase<VectorType>   &pre_smooth,
                                  const MGSmootherBase<VectorType>   &post_smooth,
                                  Cycle                              cycle)
  :
  cycle_type(cycle),
  minlevel(0),
  maxlevel(mg_dof_handler.get_triangulation().n_global_levels()-1),
  defect(minlevel,maxlevel),
  solution(minlevel,maxlevel),
  t(minlevel,maxlevel),
  defect2(minlevel,maxlevel),
  matrix(&matrix, typeid(*this).name()),
  coarse(&coarse, typeid(*this).name()),
  transfer(&transfer, typeid(*this).name()),
  pre_smooth(&pre_smooth, typeid(*this).name()),
  post_smooth(&post_smooth, typeid(*this).name()),
  edge_down(0, typeid(*this).name()),
  edge_up(0, typeid(*this).name()),
  debug(0)
{}



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


template<int dim, typename VectorType, class TRANSFER>
PreconditionMG<dim, VectorType, TRANSFER>
::PreconditionMG(const DoFHandler<dim>  &dof_handler,
                 Multigrid<VectorType>  &mg,
                 const TRANSFER         &transfer)
  :
  dof_handler(&dof_handler),
  multigrid(&mg),
  transfer(&transfer)
{}

template<int dim, typename VectorType, class TRANSFER>
inline bool
PreconditionMG<dim, VectorType, TRANSFER>::empty () const
{
  return false;
}

template<int dim, typename VectorType, class TRANSFER>
template<class OtherVectorType>
void
PreconditionMG<dim, VectorType, TRANSFER>::vmult
(OtherVectorType       &dst,
 const OtherVectorType &src) const
{
  transfer->copy_to_mg(*dof_handler,
                       multigrid->defect,
                       src);
  multigrid->cycle();

  transfer->copy_from_mg(*dof_handler,
                         dst,
                         multigrid->solution);
}


template<int dim, typename VectorType, class TRANSFER>
template<class OtherVectorType>
void
PreconditionMG<dim, VectorType, TRANSFER>::vmult_add
(OtherVectorType       &dst,
 const OtherVectorType &src) const
{
  transfer->copy_to_mg(*dof_handler,
                       multigrid->defect,
                       src);
  multigrid->cycle();
  transfer->copy_from_mg_add(*dof_handler,
                             dst,
                             multigrid->solution);
}


template<int dim, typename VectorType, class TRANSFER>
template<class OtherVectorType>
void
PreconditionMG<dim, VectorType, TRANSFER>::Tvmult
(OtherVectorType &,
 const OtherVectorType &) const
{
  Assert(false, ExcNotImplemented());
}


template<int dim, typename VectorType, class TRANSFER>
template<class OtherVectorType>
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
