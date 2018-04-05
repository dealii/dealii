// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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

#ifndef dealii_dof_renumbering_h
#define dealii_dof_renumbering_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/hp/dof_handler.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * Implementation of a number of renumbering algorithms for the degrees of
 * freedom on a triangulation. The functions in this namespace compute
 * new indices for each degree of freedom of a DoFHandler or
 * hp::DoFHandler object, and then call DoFHandler::renumber_dofs() or
 * hp::DoFHandler::renumber_dofs().
 *
 *
 * <h3>Cuthill-McKee like algorithms</h3>
 *
 * Within this class, the Cuthill-McKee algorithm is implemented. It starts at
 * a degree of freedom, searches the other DoFs for those which are coupled
 * with the one we started with and numbers these in a certain way. It then
 * finds the second level of DoFs, namely those that couple with those of the
 * previous level (which were those that coupled with the initial DoF) and
 * numbers these. And so on. For the details of the algorithm, especially the
 * numbering within each level, please see H. R. Schwarz: "Methode der finiten
 * Elemente". The reverse Cuthill-McKee algorithm does the same job, but
 * numbers all elements in the reverse order.
 *
 * These algorithms have one major drawback: they require a good starting
 * point, i.e. the degree of freedom index that will get a new index of zero.
 * The renumbering functions therefore allow the caller to specify such an
 * initial DoF, e.g. by exploiting knowledge of the actual topology of the
 * domain. It is also possible to give several starting indices, which may be
 * used to simulate a simple upstream numbering (by giving the inflow dofs as
 * starting values) or to make preconditioning faster (by letting the
 * Dirichlet boundary indices be starting points).
 *
 * If no starting index is given, one is chosen automatically, namely one with
 * the smallest coordination number (the coordination number is the number of
 * other dofs this dof couples with). This dof is usually located on the
 * boundary of the domain. There is, however, large ambiguity in this when
 * using the hierarchical meshes used in this library, since in most cases the
 * computational domain is not approximated by tilting and deforming elements
 * and by plugging together variable numbers of elements at vertices, but
 * rather by hierarchical refinement. There is therefore a large number of
 * dofs with equal coordination numbers. The renumbering algorithms will
 * therefore not give optimal results.
 *
 * In the book of Schwarz (H.R.Schwarz: Methode der finiten Elemente), it is
 * advised to test many starting points, if possible all with the smallest
 * coordination number and also those with slightly higher numbers. However,
 * this seems only possible for meshes with at most several dozen or a few
 * hundred elements found in small engineering problems of the early 1980s
 * (the second edition was published in 1984), but certainly not with those
 * used in this library, featuring several 10,000 to a few 100,000 elements.
 *
 *
 * <h4>Implementation of renumbering schemes</h4>
 *
 * The renumbering algorithms need quite a lot of memory, since they have to
 * store for each dof with which other dofs it couples. This is done using a
 * SparsityPattern object used to store the sparsity pattern of matrices. It
 * is not useful for the user to do anything between distributing the dofs and
 * renumbering, i.e. the calls to DoFHandler::distribute_dofs and
 * DoFHandler::renumber_dofs should follow each other immediately. If you try
 * to create a sparsity pattern or anything else in between, these will be
 * invalid afterwards.
 *
 * The renumbering may take care of dof-to-dof couplings only induced by
 * eliminating constraints. In addition to the memory consumption mentioned
 * above, this also takes quite some computational time, but it may be
 * switched off upon calling the @p renumber_dofs function. This will then
 * give inferior results, since knots in the graph (representing dofs) are not
 * found to be neighbors even if they would be after condensation.
 *
 * The renumbering algorithms work on a purely algebraic basis, due to the
 * isomorphism between the graph theoretical groundwork underlying the
 * algorithms and binary matrices (matrices of which the entries are binary
 * values) represented by the sparsity patterns. In special, the algorithms do
 * not try to exploit topological knowledge (e.g. corner detection) to find
 * appropriate starting points. This way, however, they work in arbitrary
 * space dimension.
 *
 * If you want to give starting points, you may give a list of dof indices
 * which will form the first step of the renumbering. The dofs of the list
 * will be consecutively numbered starting with zero, i.e. this list is not
 * renumbered according to the coordination number of the nodes. Indices not
 * in the allowed range are deleted. If no index is allowed, the algorithm
 * will search for its own starting point.
 *
 *
 * <h4>Results of renumbering</h4>
 *
 * The renumbering schemes mentioned above do not lead to optimal results.
 * However, after all there is no algorithm that accomplishes this within
 * reasonable time. There are situations where the lack of optimality even
 * leads to worse results than with the original, crude, levelwise numbering
 * scheme; one of these examples is a mesh of four cells of which always those
 * cells are refined which are neighbors to the center (you may call this mesh
 * a `zoom in' mesh). In one such example the bandwidth was increased by about
 * 50 per cent.
 *
 * In most other cases, the bandwidth is reduced significantly. The reduction
 * is the better the less structured the grid is. With one grid where the
 * cells were refined according to a random driven algorithm, the bandwidth
 * was reduced by a factor of six.
 *
 * Using the constraint information usually leads to reductions in bandwidth
 * of 10 or 20 per cent, but may for some very unstructured grids also lead to
 * an increase. You have to weigh the decrease in your case with the time
 * spent to use the constraint information, which usually is several times
 * longer than the `pure' renumbering algorithm.
 *
 * In almost all cases, the renumbering scheme finds a corner to start with.
 * Since there is more than one corner in most grids and since even an
 * interior degree of freedom may be a better starting point, giving the
 * starting point by the user may be a viable way if you have a simple scheme
 * to derive a suitable point (e.g. by successively taking the third child of
 * the cell top left of the coarsest level, taking its third vertex and the
 * dof index thereof, if you want the top left corner vertex). If you do not
 * know beforehand what your grid will look like (e.g. when using adaptive
 * algorithms), searching a best starting point may be difficult, however, and
 * in many cases will not justify the effort.
 *
 *
 * <h3>Component-wise and block-wise numberings</h3>
 *
 * For finite elements composed of several base elements using the FESystem
 * class, or for elements which provide several components themselves, it may
 * be of interest to sort the DoF indices by component. This will then bring
 * out the block matrix structure, since otherwise the degrees of freedom are
 * numbered cell-wise without taking into account that they may belong to
 * different components. For example, one may want to sort degree of freedom
 * for a Stokes discretization so that we first get all velocities and then
 * all the pressures so that the resulting matrix naturally decomposes into a
 * $2\times 2$ system.
 *
 * This kind of numbering may be obtained by calling the component_wise()
 * function of this class. Since it does not touch the order of indices within
 * each component, it may be worthwhile to first renumber using the Cuthill-
 * McKee or a similar algorithm and afterwards renumbering component-wise.
 * This will bring out the matrix structure and additionally have a good
 * numbering within each block.
 *
 * The component_wise() function allows not only to honor enumeration based on
 * vector components, but also allows to group together vector components into
 * "blocks" using a defaulted argument to the various
 * DoFRenumbering::component_wise() functions (see
 * @ref GlossComponent
 * vs
 * @ref GlossBlock
 * for a description of the difference). The blocks designated through this
 * argument may, but do not have to be, equal to the blocks that the finite
 * element reports. For example, a typical Stokes element would be
 * @code
 *   FESystem<dim> stokes_fe (FE_Q<dim>(2), dim,   // dim velocities
 *                            FE_Q<dim>(1), 1);    // one pressure
 * @endcode
 * This element has <code>dim+1</code> vector components and equally many
 * blocks. However, one may want to consider the velocities as one logical
 * block so that all velocity degrees of freedom are enumerated the same way,
 * independent of whether they are $x$- or $y$-velocities. This is done, for
 * example, in step-20 and step-22 as well as several other tutorial programs.
 *
 * On the other hand, if you really want to use block structure reported by
 * the finite element itself (a case that is often the case if you have finite
 * elements that have multiple vector components, e.g. the FE_RaviartThomas or
 * FE_Nedelec elements) then you can use the DoFRenumbering::block_wise
 * instead of the DoFRenumbering::component_wise functions.
 *
 *
 * <h3>Cell-wise numbering</h3>
 *
 * Given an ordered vector of cells, the function cell_wise() sorts the
 * degrees of freedom such that degrees on earlier cells of this vector will
 * occur before degrees on later cells.
 *
 * This rule produces a well-defined ordering for discontinuous Galerkin
 * methods (FE_DGP, FE_DGQ). For continuous methods, we use the additional
 * rule that each degree of freedom is ordered according to the first cell in
 * the ordered vector it belongs to.
 *
 * Applications of this scheme are downstream() and clock_wise_dg(). The first
 * orders the cells according to a downstream direction and then applies
 * cell_wise().
 *
 * @note For DG elements, the internal numbering in each cell remains
 * unaffected. This cannot be guaranteed for continuous elements anymore,
 * since degrees of freedom shared with an earlier cell will be accounted for
 * by the other cell.
 *
 *
 * <h3>Random renumbering</h3>
 *
 * The random() function renumbers degrees of freedom randomly. This function
 * is probably seldom of use, except to check the dependence of solvers
 * (iterative or direct ones) on the numbering of the degrees of freedom.
 *
 *
 * <h3>A comparison of reordering strategies</h3>
 *
 * As a benchmark of comparison, let us consider what the different sparsity
 * patterns produced by the various algorithms when using the $Q_2^d\times
 * Q_1$ element combination typically employed in the discretization of Stokes
 * equations, when used on the mesh obtained in step-22 after one adaptive
 * mesh refinement in 3d. The space dimension together with the coupled finite
 * element leads to a rather dense system matrix with, on average around 180
 * nonzero entries per row. After applying each of the reordering strategies
 * shown below, the degrees of freedom are also sorted using
 * DoFRenumbering::component_wise into velocity and pressure groups; this
 * produces the $2\times 2$ block structure seen below with the large
 * velocity-velocity block at top left, small pressure-pressure block at
 * bottom right, and coupling blocks at top right and bottom left.
 *
 * The goal of reordering strategies is to improve the preconditioner. In
 * step-22 we use a SparseILU to preconditioner for the velocity-velocity
 * block at the top left. The quality of the preconditioner can then be
 * measured by the number of CG iterations required to solve a linear system
 * with this block. For some of the reordering strategies below we record this
 * number for adaptive refinement cycle 3, with 93176 degrees of freedom;
 * because we solve several linear systems with the same matrix in the Schur
 * complement, the average number of iterations is reported. The lower the
 * number the better the preconditioner and consequently the better the
 * renumbering of degrees of freedom is suited for this task. We also state
 * the run-time of the program, in part determined by the number of iterations
 * needed, for the first 4 cycles on one of our machines. Note that the
 * reported times correspond to the run time of the entire program, not just
 * the affected solver; if a program runs twice as fast with one particular
 * ordering than with another one, then this means that the actual solver is
 * actually several times faster.
 *
 * <table> <tr> <td>
 * @image html "reorder_sparsity_step_31_original.png"
 * </td> <td>
 * @image html "reorder_sparsity_step_31_random.png"
 * </td> <td>
 * @image html "reorder_sparsity_step_31_deal_cmk.png"
 * </td> </tr> <tr> <td> Enumeration as produced by deal.II's
 * DoFHandler::distribute_dofs function and no further reordering apart from
 * the component-wise one.
 *
 * With this renumbering, we needed an average of 92.2 iterations for the
 * testcase outlined above, and a runtime of 7min53s. </td> <td> Random
 * enumeration as produced by applying DoFRenumbering::random after calling
 * DoFHandler::distribute_dofs. This enumeration produces nonzero entries in
 * matrices pretty much everywhere, appearing here as an entirely unstructured
 * matrix.
 *
 * With this renumbering, we needed an average of 71 iterations for the
 * testcase outlined above, and a runtime of 10min55s. The longer runtime
 * despite less iterations compared to the default ordering may be due to the
 * fact that computing and applying the ILU requires us to jump back and forth
 * all through memory due to the lack of localization of matrix entries around
 * the diagonal; this then leads to many cache misses and consequently bad
 * timings. </td> <td> Cuthill-McKee enumeration as produced by calling the
 * deal.II implementation of the algorithm provided by
 * DoFRenumbering::Cuthill_McKee after DoFHandler::distribute_dofs.
 *
 * With this renumbering, we needed an average of 57.3 iterations for the
 * testcase outlined above, and a runtime of 6min10s. </td> </td> </tr>
 *
 * <tr> <td>
 * @image html "reorder_sparsity_step_31_boost_cmk.png"
 * </td> <td>
 * @image html "reorder_sparsity_step_31_boost_king.png"
 * </td> <td>
 * @image html "reorder_sparsity_step_31_boost_md.png"
 * </td> </tr> <tr> <td> Cuthill- McKee enumeration as produced by calling the
 * BOOST implementation of the algorithm provided by
 * DoFRenumbering::boost::Cuthill_McKee after DoFHandler::distribute_dofs.
 *
 * With this renumbering, we needed an average of 51.7 iterations for the
 * testcase outlined above, and a runtime of 5min52s. </td> <td> King
 * enumeration as produced by calling the BOOST implementation of the
 * algorithm provided by DoFRenumbering::boost::king_ordering after
 * DoFHandler::distribute_dofs. The sparsity pattern appears denser than with
 * BOOST's Cuthill-McKee algorithm; however, this is only an illusion: the
 * number of nonzero entries is the same, they are simply not as well
 * clustered.
 *
 * With this renumbering, we needed an average of 51.0 iterations for the
 * testcase outlined above, and a runtime of 5min03s. Although the number of
 * iterations is only slightly less than with BOOST's Cuthill-McKee
 * implementation, runtime is significantly less. This, again, may be due to
 * cache effects. As a consequence, this is the algorithm best suited to the
 * testcase, and is in fact used in step-22. </td> <td> Minimum degree
 * enumeration as produced by calling the BOOST implementation of the
 * algorithm provided by DoFRenumbering::boost::minimum_degree after
 * DoFHandler::distribute_dofs. The minimum degree algorithm does not attempt
 * to minimize the bandwidth of a matrix but to minimize the amount of fill-in
 * a LU decomposition would produce, i.e. the number of places in the matrix
 * that would be occupied by elements of an LU decomposition that are not
 * already occupied by elements of the original matrix. The resulting sparsity
 * pattern obviously has an entirely different structure than the ones
 * produced by algorithms trying to minimize the bandwidth.
 *
 * With this renumbering, we needed an average of 58.9 iterations for the
 * testcase outlined above, and a runtime of 6min11s. </td> </tr>
 *
 * <tr> <td>
 * @image html "reorder_sparsity_step_31_downstream.png"
 * </td> <td> </td> <td> </td> </tr> <tr> <td> Downstream enumeration using
 * DoFRenumbering::downstream using a direction that points diagonally through
 * the domain.
 *
 * With this renumbering, we needed an average of 90.5 iterations for the
 * testcase outlined above, and a runtime of 7min05s. </td> <td> </td> <td>
 * </td> </tr> </table>
 *
 *
 * <h3>Multigrid DoF numbering</h3>
 *
 * Most of the algorithms listed above also work on multigrid degree of
 * freedom numberings. Refer to the actual function declarations to get more
 * information on this.
 *
 * @ingroup dofs
 * @author Wolfgang Bangerth, Guido Kanschat, 1998, 1999, 2000, 2004, 2007,
 * 2008
 */
namespace DoFRenumbering
{
  /**
   * Direction based comparator for cell iterators: it returns @p true if the
   * center of the second cell is downstream of the center of the first one
   * with respect to the direction given to the constructor.
   */
  template <class Iterator, int dim>
  struct CompareDownstream
  {
    /**
     * Constructor.
     */
    CompareDownstream (const Tensor<1,dim> &dir)
      :
      dir(dir)
    {}
    /**
     * Return true if c1 less c2.
     */
    bool operator () (const Iterator &c1, const Iterator &c2) const
    {
      const Tensor<1,dim> diff = c2->center() - c1->center();
      return (diff*dir > 0);
    }

  private:
    /**
     * Flow direction.
     */
    const Tensor<1,dim> dir;
  };


  /**
   * Point based comparator for downstream directions: it returns @p true if
   * the second point is downstream of the first one with respect to the
   * direction given to the constructor. If the points are the same with
   * respect to the downstream direction, the point with the lower DoF number
   * is considered smaller.
   */
  template <int dim>
  struct ComparePointwiseDownstream
  {
    /**
     * Constructor.
     */
    ComparePointwiseDownstream (const Tensor<1,dim> &dir)
      :
      dir(dir)
    {}
    /**
     * Return true if c1 less c2.
     */
    bool operator () (const std::pair<Point<dim>,types::global_dof_index> &c1,
                      const std::pair<Point<dim>,types::global_dof_index> &c2) const
    {
      const Tensor<1,dim> diff = c2.first-c1.first;
      return (diff*dir > 0 || (diff*dir==0 && c1.second<c2.second));
    }

  private:
    /**
     * Flow direction.
     */
    const Tensor<1,dim> dir;
  };



  /**
   * A namespace for the implementation of some renumbering algorithms based
   * on algorithms implemented in the Boost Graph Library (BGL) by Jeremy Siek
   * and others.
   *
   * While often slightly slower to compute, the algorithms using BOOST often
   * lead to matrices with smaller bandwidths and sparse ILUs based on this
   * numbering are therefore more efficient.
   *
   * For a comparison of these algorithms with the ones defined in
   * DoFRenumbering, see the comparison section in the documentation of the
   * DoFRenumbering namespace.
   */
  namespace boost
  {
    /**
     * Renumber the degrees of freedom according to the Cuthill-McKee method,
     * eventually using the reverse numbering scheme.
     *
     * See the general documentation of the parent class for details on the
     * different methods.
     *
     * As an example of the results of this algorithm, take a look at the
     * comparison of various algorithms in the documentation of the
     * DoFRenumbering namespace.
     */
    template <typename DoFHandlerType>
    void
    Cuthill_McKee (DoFHandlerType &dof_handler,
                   const bool      reversed_numbering = false,
                   const bool      use_constraints    = false);

    /**
     * Compute the renumbering vector needed by the Cuthill_McKee() function.
     * Does not perform the renumbering on the DoFHandler dofs but returns the
     * renumbering vector.
     */
    template <typename DoFHandlerType>
    void
    compute_Cuthill_McKee (std::vector<types::global_dof_index> &new_dof_indices,
                           const DoFHandlerType &,
                           const bool                            reversed_numbering = false,
                           const bool                            use_constraints    = false);

    /**
     * Renumber the degrees of freedom based on the BOOST implementation of
     * the King algorithm. This often results in slightly larger (by a few
     * percent) bandwidths than the Cuthill-McKee algorithm, but sparse ILUs
     * are often slightly (also by a few percent) better preconditioners.
     *
     * As an example of the results of this algorithm, take a look at the
     * comparison of various algorithms in the documentation of the
     * DoFRenumbering namespace.
     *
     * This algorithm is used in step-22.
     */
    template <typename DoFHandlerType>
    void
    king_ordering (DoFHandlerType &dof_handler,
                   const bool      reversed_numbering = false,
                   const bool      use_constraints    = false);

    /**
     * Compute the renumbering for the King algorithm but do not actually
     * renumber the degrees of freedom in the DoF handler argument.
     */
    template <typename DoFHandlerType>
    void
    compute_king_ordering (std::vector<types::global_dof_index> &new_dof_indices,
                           const DoFHandlerType &,
                           const bool                            reversed_numbering = false,
                           const bool                            use_constraints    = false);

    /**
     * Renumber the degrees of freedom based on the BOOST implementation of
     * the minimum degree algorithm. Unlike the Cuthill-McKee algorithm, this
     * algorithm does not attempt to minimize the bandwidth of a matrix but to
     * minimize the amount of fill-in when doing an LU decomposition. It may
     * sometimes yield better ILUs because of this property.
     *
     * As an example of the results of this algorithm, take a look at the
     * comparison of various algorithms in the documentation of the
     * DoFRenumbering namespace.
     */
    template <typename DoFHandlerType>
    void
    minimum_degree (DoFHandlerType &dof_handler,
                    const bool      reversed_numbering = false,
                    const bool      use_constraints    = false);

    /**
     * Compute the renumbering for the minimum degree algorithm but do not
     * actually renumber the degrees of freedom in the DoF handler argument.
     */
    template <typename DoFHandlerType>
    void
    compute_minimum_degree (std::vector<types::global_dof_index> &new_dof_indices,
                            const DoFHandlerType &,
                            const bool                            reversed_numbering = false,
                            const bool                            use_constraints    = false);
  }

  /**
   * Renumber the degrees of freedom according to the Cuthill-McKee method,
   * possibly using the reverse numbering scheme.
   *
   * See the general documentation of this class for details on the different
   * methods.
   *
   * As an example of the results of this algorithm, take a look at the
   * comparison of various algorithms in the documentation of the
   * DoFRenumbering namespace.
   *
   * @param dof_handler The DoFHandler or hp::DoFHandler object to work on.
   * @param reversed_numbering Whether to use the original Cuthill-McKee
   *   algorithm, or to reverse the ordering.
   * @param use_constraints Whether or not to use hanging node constraints in
   *   determining the reordering of degrees of freedom.
   * @param starting_indices A set of degrees of freedom that form the first
   *   level of renumbered degrees of freedom. If the set is empty, then a
   *   single starting entry is chosen automatically among those that have the
   *   smallest number of others that couple with it.
   *
   * <h4> Operation in parallel </h4>
   *
   * If the given DoFHandler uses a distributed triangulation (i.e., if
   * dof_handler.locally_owned() is not the complete index set), the
   * renumbering is performed on each processor's degrees of freedom
   * individually, without any communication between processors. In other
   * words, the resulting renumbering is an attempt at minimizing the bandwidth
   * of <i>each diagonal block of the matrix corresponding to one processor</i>
   * separately, without making an attempt at minimizing the bandwidth of
   * the global matrix. Furthermore, the renumbering reuses exactly the
   * same set of DoF indices that each processor used before. In other words,
   * if the previous numbering of DoFs on one processor used a contiguous
   * range of DoF indices, then so will the DoFs on that processor after
   * the renumbering, and they will occupy the same range. The same is true
   * if the previous numbering of DoFs on a processor consisted of a number
   * of index ranges or single indices: after renumbering, the locally owned
   * DoFs on that processor will use the exact same indices, just in a
   * different order.
   *
   * In addition, if the DoFHandler is built on a parallel triangulation, then
   * on every processor, the starting indices for renumbering need to be a
   * (possibly empty) subset of the
   * @ref GlossLocallyActiveDof "locally active degrees of freedom". In
   * general, these starting indices will be different on each processor
   * (unless of course you pass an empty list as is the default),
   * and each processor will use them as starting indices for the local
   * renumbering on that processor.
   *
   * The starting indices must be locally active degrees of freedom, but
   * the function will only renumber the locally owned subset of the
   * locally owned DoFs. The function accepts starting indices from the
   * largest set of locally active degrees of freedom because a typical
   * renumbering operation with this function starts with indices that
   * are located on the boundary -- in the case of the current function,
   * that would be the boundary between processor subdomains. Since the
   * degrees of freedom that are located on subdomain interfaces may
   * be owned by either one of the two processors that own the adjacent
   * subdomains, it is not always easy to identify starting indices that
   * are locally owned. On the other hand, all degrees of freedom on subdomain
   * interfaces are locally active, and so the function accepts them as
   * starting indices even though it can only renumber them on a given
   * processor if they are also locally owned.
   */
  template <typename DoFHandlerType>
  void
  Cuthill_McKee (DoFHandlerType                             &dof_handler,
                 const bool                                  reversed_numbering = false,
                 const bool                                  use_constraints    = false,
                 const std::vector<types::global_dof_index> &starting_indices
                 = std::vector<types::global_dof_index>());

  /**
   * Compute the renumbering vector needed by the Cuthill_McKee() function.
   * This function does not perform the renumbering on the DoFHandler DoFs but
   * only returns the renumbering vector.
   *
   * See the Cuthill_McKee() function for an explanation of the arguments.
   */
  template <typename DoFHandlerType>
  void
  compute_Cuthill_McKee (std::vector<types::global_dof_index>       &new_dof_indices,
                         const DoFHandlerType &,
                         const bool                                  reversed_numbering = false,
                         const bool                                  use_constraints    = false,
                         const std::vector<types::global_dof_index> &starting_indices
                         = std::vector<types::global_dof_index>());

  /**
   * Renumber the degrees of freedom according to the Cuthill-McKee method,
   * eventually using the reverse numbering scheme, in this case for a
   * multigrid numbering of degrees of freedom.
   *
   * You can give a triangulation level to which this function is to be
   * applied.  Since with a level-wise numbering there are no hanging nodes,
   * no constraints can be used, so the respective parameter of the previous
   * function is omitted.
   *
   * See the general documentation of this class for details on the different
   * methods.
   */
  template <typename DoFHandlerType>
  void
  Cuthill_McKee (DoFHandlerType                             &dof_handler,
                 const unsigned int                          level,
                 const bool                                  reversed_numbering = false,
                 const std::vector<types::global_dof_index> &starting_indices
                 = std::vector<types::global_dof_index> ());

  /**
   * @name Component-wise numberings
   * @{
   */

  /**
   * Sort the degrees of freedom by vector component. The numbering within
   * each component is not touched, so a degree of freedom with index $i$,
   * belonging to some component, and another degree of freedom with index $j$
   * belonging to the same component will be assigned new indices $n(i)$ and
   * $n(j)$ with $n(i)<n(j)$ if $i<j$ and $n(i)>n(j)$ if $i>j$.
   *
   * You can specify that the components are ordered in a different way than
   * suggested by the FESystem object you use. To this end, set up the vector
   * @p target_component such that the entry at index @p i denotes the number
   * of the target component for dofs with component @p i in the FESystem.
   * Naming the same target component more than once is possible and results
   * in a blocking of several components into one. This is discussed in
   * step-22. If you omit this argument, the same order as given by the finite
   * element is used.
   *
   * If one of the base finite elements from which the global finite element
   * under consideration here, is a non-primitive one, i.e. its shape
   * functions have more than one non-zero component, then it is not possible
   * to associate these degrees of freedom with a single vector component. In
   * this case, they are associated with the first vector component to which
   * they belong.
   *
   * For finite elements with only one component, or a single non-primitive
   * base element, this function is the identity operation.
   */
  template <typename DoFHandlerType>
  void
  component_wise (DoFHandlerType                  &dof_handler,
                  const std::vector<unsigned int> &target_component
                  = std::vector<unsigned int>());


  /**
   * Sort the degrees of freedom by component. It does the same thing as the
   * above function, only that it does this for one single level of a
   * multilevel discretization. The non-multigrid part of the DoFHandler
   * is not touched.
   */
  template <typename DoFHandlerType>
  void
  component_wise (DoFHandlerType                  &dof_handler,
                  const unsigned int               level,
                  const std::vector<unsigned int> &target_component
                  = std::vector<unsigned int>());

  /**
   * Compute the renumbering vector needed by the component_wise() functions.
   * Does not perform the renumbering on the DoFHandler dofs but returns the
   * renumbering vector.
   */
  template <int dim, int spacedim, typename CellIterator>
  types::global_dof_index
  compute_component_wise (std::vector<types::global_dof_index>        &new_dof_indices,
                          const CellIterator                          &start,
                          const typename identity<CellIterator>::type &end,
                          const std::vector<unsigned int>             &target_component,
                          const bool                                   is_level_operation);

  /**
   * @}
   */

  /**
   * @name Block-wise numberings
   * @{
   */

  /**
   * Sort the degrees of freedom by vector block. The numbering within each
   * block is not touched, so a degree of freedom with index $i$, belonging to
   * some block, and another degree of freedom with index $j$ belonging to the
   * same block will be assigned new indices $n(i)$ and $n(j)$ with
   * $n(i)<n(j)$ if $i<j$ and $n(i)>n(j)$ if $i>j$.
   */
  template <int dim, int spacedim>
  void
  block_wise (DoFHandler<dim,spacedim> &dof_handler);

  /**
   * Sort the degrees of freedom by vector block. It does the same thing as
   * the above function, only that it does this for one single level of a
   * multilevel discretization. The non-multigrid part of the DoFHandler
   * is not touched.
   */
  template <int dim, int spacedim>
  void
  block_wise (DoFHandler<dim,spacedim> &dof_handler, const unsigned int level);

  /**
   * Sort the degrees of freedom by block. It does the same thing as the above
   * function.
   *
   * This function only succeeds if each of the elements in the
   * hp::FECollection attached to the hp::DoFHandler argument has exactly the
   * same number of blocks (see
   * @ref GlossBlock "the glossary"
   * for more information). Note that this is not always given: while the
   * hp::FECollection class ensures that all of its elements have the same
   * number of vector components, they need not have the same number of
   * blocks. At the same time, this function here needs to match individual
   * blocks across elements and therefore requires that elements have the same
   * number of blocks and that subsequent blocks in one element have the same
   * meaning as in another element.
   */
  template <int dim, int spacedim>
  void
  block_wise (hp::DoFHandler<dim,spacedim> &dof_handler);

  /**
   * Compute the renumbering vector needed by the block_wise() functions.
   * Does not perform the renumbering on the DoFHandler dofs but returns the
   * renumbering vector.
   */
  template <int dim, int spacedim, class ITERATOR, class ENDITERATOR>
  types::global_dof_index
  compute_block_wise (std::vector<types::global_dof_index> &new_dof_indices,
                      const ITERATOR &start,
                      const ENDITERATOR &end,
                      bool is_level_operation);

  /**
   * @}
   */

  /**
   * @name Various cell-wise numberings
   * @{
   */

  /**
   * Renumber the degrees cell by cell by traversing the cells in
   * @ref GlossZOrder "Z order".
   *
   * There are two reasons to use this function:
   * - It produces a predictable ordering of degrees of freedom
   *   that is independent of how exactly you arrived at a mesh.
   *   In particular, in general the order of cells of a mesh
   *   depends on the order in which cells were marked for
   *   refinement and coarsening during the refinement cycles
   *   the mesh has undergone. On the other hand, the z-order
   *   of cells is independent of the mesh's history, and so yields a
   *   predictable DoF numbering.
   * - For meshes based on parallel::distributed::Triangulation,
   *   the @ref GlossLocallyOwnedCell "locally owned cells" of
   *   each MPI process are contiguous in Z order. That means that
   *   numbering degrees of freedom by visiting cells in Z order yields
   *   @ref GlossLocallyOwnedDof "locally owned DoF indices" that consist
   *   of contiguous ranges for each process. This is also true for the
   *   default ordering of DoFs on such triangulations, but the default
   *   ordering creates an enumeration that also depends on how many
   *   processors participate in the mesh, whereas the one generated
   *   by this function enumerates the degrees of freedom on a particular
   *   cell with indices that will be the same regardless of how many
   *   processes the mesh is split up between.
   *
   * For meshes based on parallel::shared::Triangulation, the situation is
   * more complex. Here, the set of locally owned cells is determined by
   * a partitioning algorithm (selected by passing an object of type
   * parallel::shared::Triangulation::Settings to the constructor of the
   * triangulation), and in general these partitioning algorithms may
   * assign cells to @ref GlossSubdomainId "subdomains" based on
   * decisions that may have nothing to do with the Z order. (Though it
   * is possible to select these flags in a way so that partitioning
   * uses the Z order.) As a consequence, the cells of one subdomain
   * are not contiguous in Z order, and if one renumbered degrees of freedom
   * based on the Z order of cells, one would generally end up with DoF
   * indices that on each processor do not form a contiguous range.
   * This is often inconvenient (for example, because PETSc cannot store
   * vectors and matrices for which the locally owned set of indices
   * is not contiguous), and consequently this function uses the following
   * algorithm for parallel::shared::Triangulation objects:
   * - It determines how many degrees of freedom each processor owns.
   *   This is an invariant under renumbering, and consequently we can
   *   use how many DoFs each processor owns at the beginning of the current
   *   function. Let us call this number $n_P$ for processor $P$.
   * - It determines for each processor a contiguous range of new
   *   DoF indices $[b_P,e_P)$ so that $e_P-b_P=n_P$, $b_0=0$, and
   *   $b_P=e_{P-1}$.
   * - It traverses the <i>locally owned cells</i> in Z order and
   *   renumbers the locally owned degrees of freedom on these cells
   *   so that the new numbers fit within the interval $[b_P,e_P)$.
   * In other words, the <i>locally owned degrees of freedom</i> on each
   * processor are sorted according to the Z order of the locally
   * owned cells they are on, but this property may not hold globally,
   * across cells. This is because the partitioning algorithm may have
   * decided that, for example, processor 0 owns a cell that comes
   * <i>later</i> in Z order than one of the cells assigned to processor 1.
   * On the other hand, the algorithm described above assigns the
   * degrees of freedom on this cell <i>earlier</i> indices than
   * all of the indices owned by processor 1.
   *
   * @note This function generates an ordering that is independent of the previous
   * numbering of degrees of freedom. In other words, any information that may
   * have been produced by a previous call to a renumbering function is
   * ignored.
   */
  template <int dim>
  void
  hierarchical (DoFHandler<dim> &dof_handler);

  /**
   * Renumber degrees of freedom by cell. The function takes a vector of cell
   * iterators (which needs to list <i>all</i> active cells of the DoF handler
   * objects) and will give degrees of freedom new indices based on where in
   * the given list of cells the cell is on which the degree of freedom is
   * located. Degrees of freedom that exist at the interface between two or
   * more cells will be numbered when they are encountered first.
   *
   * Degrees of freedom that are encountered first on the same cell retain
   * their original ordering before the renumbering step.
   *
   * @param[in,out] dof_handler The DoFHandler whose degrees of freedom are to
   * be renumbered.
   * @param[in] cell_order A vector that contains the order of the cells that
   * defines the order in which degrees of freedom should be renumbered.
   *
   * @pre @p cell_order must have size
   * <code>dof_handler.get_triangulation().n_active_cells()</code>. Every
   * active cell iterator of that triangulation needs to be present in @p
   * cell_order exactly once.
   */
  template <typename DoFHandlerType>
  void
  cell_wise (DoFHandlerType &dof_handler,
             const std::vector<typename DoFHandlerType::active_cell_iterator> &cell_order);

  /**
   * Compute a renumbering of degrees of freedom by cell. The function takes a
   * vector of cell iterators (which needs to list <i>all</i> active cells of
   * the DoF handler objects) and will give degrees of freedom new indices
   * based on where in the given list of cells the cell is on which the degree
   * of freedom is located. Degrees of freedom that exist at the interface
   * between two or more cells will be numbered when they are encountered
   * first.
   *
   * Degrees of freedom that are encountered first on the same cell retain
   * their original ordering before the renumbering step.
   *
   * @param[out] renumbering A vector of length
   * <code>dof_handler.n_dofs()</code> that contains for each degree of
   * freedom (in their current numbering) their future DoF index. This vector
   * therefore presents a (very particular) <i>permutation</i> of the current
   * DoF indices.
   * @param[out] inverse_renumbering The reverse of the permutation returned
   * in the previous argument.
   * @param[in] dof_handler The DoFHandler whose degrees of freedom are to be
   * renumbered.
   * @param[in] cell_order A vector that contains the order of the cells that
   * defines the order in which degrees of freedom should be renumbered.
   *
   * @pre @p cell_order must have size
   * <code>dof_handler.get_triangulation().n_active_cells()</code>. Every
   * active cell iterator of that triangulation needs to be present in @p
   * cell_order exactly once. @post For each @p i between zero and
   * <code>dof_handler.n_dofs()</code>, the condition
   * <code>renumbering[inverse_renumbering[i]] == i</code> will hold.
   */
  template <typename DoFHandlerType>
  void
  compute_cell_wise
  (std::vector<types::global_dof_index>                             &renumbering,
   std::vector<types::global_dof_index>                             &inverse_renumbering,
   const DoFHandlerType                                             &dof_handler,
   const std::vector<typename DoFHandlerType::active_cell_iterator> &cell_order);

  /**
   * Like the other cell_wise() function, but for one level of a multilevel
   * enumeration of degrees of freedom.
   */
  template <typename DoFHandlerType>
  void
  cell_wise (DoFHandlerType                                                  &dof_handler,
             const unsigned int                                               level,
             const std::vector<typename DoFHandlerType::level_cell_iterator> &cell_order);

  /**
   * Like the other compute_cell_wise() function, but for one level of a
   * multilevel enumeration of degrees of freedom.
   */
  template <typename DoFHandlerType>
  void
  compute_cell_wise
  (std::vector<types::global_dof_index>                            &renumbering,
   std::vector<types::global_dof_index>                            &inverse_renumbering,
   const DoFHandlerType                                            &dof_handler,
   const unsigned int                                               level,
   const std::vector<typename DoFHandlerType::level_cell_iterator> &cell_order);

  /**
   * @}
   */

  /**
   * @name Directional numberings
   * @{
   */

  /**
   * Downstream numbering with respect to a constant flow direction. If the
   * additional argument @p dof_wise_renumbering is set to @p false, the
   * numbering is performed cell-wise, otherwise it is performed based on the
   * location of the support points.
   *
   * The cells are sorted such that the centers of cells numbered higher are further
   * downstream with respect to the constant vector @p direction than the
   * centers of cells numbered lower. Even if this yields a downstream numbering with
   * respect to the flux on the edges for fairly general grids, this might not
   * be guaranteed for all meshes.
   *
   * If the @p dof_wise_renumbering argument is set to @p false, this function
   * produces a downstream ordering of the mesh cells and calls cell_wise().
   * Therefore, the output only makes sense for Discontinuous Galerkin Finite
   * Elements (all degrees of freedom have to be associated with the interior
   * of the cell in that case) in that case.
   *
   * If @p dof_wise_renumbering is set to @p true, the degrees of freedom are
   * renumbered based on the support point location of the individual degrees
   * of freedom (obviously, the finite element needs to define support points
   * for this to work). The numbering of points with the same position in
   * downstream location (e.g. those parallel to the flow direction, or
   * several dofs within a FESystem) will be unaffected.
   */
  template <typename DoFHandlerType>
  void
  downstream (DoFHandlerType                                  &dof_handler,
              const Tensor<1,DoFHandlerType::space_dimension> &direction,
              const bool                                       dof_wise_renumbering = false);


  /**
   * Cell-wise downstream numbering with respect to a constant flow direction
   * on one level of a multigrid hierarchy. See the other function with the same name.
   */
  template <typename DoFHandlerType>
  void
  downstream (DoFHandlerType                                  &dof_handler,
              const unsigned int                               level,
              const Tensor<1,DoFHandlerType::space_dimension> &direction,
              const bool                                       dof_wise_renumbering = false);

  /**
   * Compute the set of renumbering indices needed by the downstream() function. Does
   * not perform the renumbering on the DoFHandler dofs but returns the
   * renumbering vector.
   */
  template <typename DoFHandlerType>
  void
  compute_downstream (std::vector<types::global_dof_index>            &new_dof_indices,
                      std::vector<types::global_dof_index>            &reverse,
                      const DoFHandlerType                            &dof_handler,
                      const Tensor<1,DoFHandlerType::space_dimension> &direction,
                      const bool                                       dof_wise_renumbering);

  /**
   * Compute the set of renumbering indices needed by the downstream() function. Does
   * not perform the renumbering on the DoFHandler dofs but returns the
   * renumbering vector.
   */
  template <typename DoFHandlerType>
  void
  compute_downstream (std::vector<types::global_dof_index>            &new_dof_indices,
                      std::vector<types::global_dof_index>            &reverse,
                      const DoFHandlerType                            &dof_handler,
                      const unsigned int                               level,
                      const Tensor<1,DoFHandlerType::space_dimension> &direction,
                      const bool                                       dof_wise_renumbering);

  /**
   * Cell-wise clockwise numbering.
   *
   * This function produces a (counter)clockwise ordering of the mesh cells
   * with respect to the hub @p center and calls cell_wise().  Therefore, it
   * only works with Discontinuous Galerkin Finite Elements, i.e. all degrees
   * of freedom have to be associated with the interior of the cell.
   */
  template <typename DoFHandlerType>
  void
  clockwise_dg (DoFHandlerType                               &dof_handler,
                const Point<DoFHandlerType::space_dimension> &center,
                const bool                                    counter = false);

  /**
   * Cell-wise clockwise numbering on one level of a multigrid
   * hierarchy. See the other function with the same name.
   */
  template <typename DoFHandlerType>
  void
  clockwise_dg (DoFHandlerType                               &dof_handler,
                const unsigned int                            level,
                const Point<DoFHandlerType::space_dimension> &center,
                const bool                                    counter = false);

  /**
   * Compute the renumbering vector needed by the clockwise_dg() functions.
   * Does not perform the renumbering on the DoFHandler dofs but returns the
   * renumbering vector.
   */
  template <typename DoFHandlerType>
  void
  compute_clockwise_dg (std::vector<types::global_dof_index>         &new_dof_indices,
                        const DoFHandlerType                         &dof_handler,
                        const Point<DoFHandlerType::space_dimension> &center,
                        const bool                                    counter);

  /**
   * @}
   */

  /**
   * @name Selective and random numberings
   * @{
   */

  /**
   * Sort those degrees of freedom which are tagged with @p true in the @p
   * selected_dofs array to the back of the DoF numbers. The sorting is
   * stable, i.e. the relative order within the tagged degrees of freedom is
   * preserved, as is the relative order within the untagged ones.
   *
   * @pre The @p selected_dofs array must have as many elements as the @p
   * dof_handler has degrees of freedom.
   */
  template <typename DoFHandlerType>
  void
  sort_selected_dofs_back (DoFHandlerType          &dof_handler,
                           const std::vector<bool> &selected_dofs);

  /**
   * Sort those degrees of freedom which are tagged with @p true in the @p
   * selected_dofs array on the level @p level to the back of the DoF numbers.
   * The sorting is stable, i.e. the relative order within the tagged degrees
   * of freedom is preserved, as is the relative order within the untagged
   * ones.
   *
   * @pre The @p selected_dofs array must have as many elements as the @p
   * dof_handler has degrees of freedom on the given level.
   */
  template <typename DoFHandlerType>
  void
  sort_selected_dofs_back (DoFHandlerType          &dof_handler,
                           const std::vector<bool> &selected_dofs,
                           const unsigned int       level);

  /**
   * Compute the renumbering vector needed by the sort_selected_dofs_back()
   * function. Does not perform the renumbering on the DoFHandler dofs but
   * returns the renumbering vector.
   *
   * @pre The @p selected_dofs array must have as many elements as the @p
   * dof_handler has degrees of freedom.
   */
  template <typename DoFHandlerType>
  void
  compute_sort_selected_dofs_back (std::vector<types::global_dof_index> &new_dof_indices,
                                   const DoFHandlerType                 &dof_handler,
                                   const std::vector<bool>              &selected_dofs);

  /**
   * This function computes the renumbering vector on each level needed by the
   * sort_selected_dofs_back() function. Does not perform the renumbering on
   * the DoFHandler dofs but only computes the renumbering and returns the
   * renumbering vector.
   *
   * @pre The @p selected_dofs array must have as many elements as the @p
   * dof_handler has degrees of freedom on the given level.
   */
  template <typename DoFHandlerType>
  void
  compute_sort_selected_dofs_back (std::vector<types::global_dof_index> &new_dof_indices,
                                   const DoFHandlerType                 &dof_handler,
                                   const std::vector<bool>              &selected_dofs,
                                   const unsigned int                    level);

  /**
   * Renumber the degrees of freedom in a random way. The result of this
   * function is repeatable in that two runs of the same program will yield
   * the same result. This is achieved by creating a new random number
   * generator with a fixed seed every time this function is entered. In
   * particular, the function therefore does not rely on an external random
   * number generator for which it would matter how often it has been called
   * before this function (or, for that matter, whether other threads running
   * concurrently to this function also draw random numbers).
   */
  template <typename DoFHandlerType>
  void
  random (DoFHandlerType &dof_handler);

  /**
   * Compute the renumbering vector needed by the random() function. See
   * there for more information on the computed random renumbering.
   *
   * This function does not perform the renumbering on the DoFHandler dofs but
   * returns the renumbering vector.
   */
  template <typename DoFHandlerType>
  void
  compute_random (std::vector<types::global_dof_index> &new_dof_indices,
                  const DoFHandlerType &dof_handler);

  /**
   * @}
   */

  /**
   * @name Numberings based on cell attributes
   * @{
   */

  /**
   * Renumber the degrees of freedom such that they are associated with the
   * subdomain id of the cells they are living on, i.e. first all degrees of
   * freedom that belong to cells with subdomain zero, then all with subdomain
   * one, etc. This is useful when doing parallel computations with a standard
   * Triangulation after assigning subdomain ids using a partitioner (see the
   * GridTools::partition_triangulation function for this). Calling this
   * function is unnecessary when using a parallel::shared::Triangulation or
   * parallel::distributed::Triangulation, as the degrees of freedom are
   * already enumerated according to the MPI process id. Therefore, if the
   * underlying triangulation is of this type then an error will be thrown.
   *
   * Note that degrees of freedom associated with faces, edges, and vertices
   * may be associated with multiple subdomains if they are sitting on
   * partition boundaries. It would therefore be undefined with which
   * subdomain they have to be associated. For this, we use what we get from
   * the DoFTools::get_subdomain_association function.
   *
   * The algorithm is stable, i.e. if two dofs i,j have <tt>i<j</tt> and
   * belong to the same subdomain, then they will be in this order also after
   * reordering.
   */
  template <typename DoFHandlerType>
  void
  subdomain_wise (DoFHandlerType &dof_handler);

  /**
   * Compute the renumbering vector needed by the subdomain_wise() function.
   * Does not perform the renumbering on the @p DoFHandler dofs but returns
   * the renumbering vector.
   */
  template <typename DoFHandlerType>
  void
  compute_subdomain_wise (std::vector<types::global_dof_index> &new_dof_indices,
                          const DoFHandlerType                 &dof_handler);

  /**
   * @}
   */



  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg (ExcDoFHandlerNotInitialized,
                    "The DoFHandler on which this function should work has not "
                    "been initialized, i.e., it doesn't appear that DoF indices "
                    "have been distributed on it.");

  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcInvalidComponentOrder);

  /**
   * The function is only implemented for Discontinuous Galerkin Finite
   * elements.
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcNotDGFEM);
}


DEAL_II_NAMESPACE_CLOSE

#endif
