// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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

#ifndef __deal2__dof_renumbering_h
#define __deal2__dof_renumbering_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/multigrid/mg_dof_handler.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * Implementation of a number of renumbering algorithms for the degrees of
 * freedom on a triangulation.
 *
 * <h3>Cuthill-McKee like algorithms</h3>
 *
 * Within this class, the Cuthill-McKee algorithm is implemented. It
 * starts at a degree of freedom, searches the other DoFs for those
 * which are coupled with the one we started with and numbers these in
 * a certain way. It then finds the second level of DoFs, namely those
 * that couple with those of the previous level (which were those that
 * coupled with the initial DoF) and numbers these. And so on. For the
 * details of the algorithm, especially the numbering within each
 * level, please see H. R. Schwarz:
 * "Methode der finiten Elemente". The reverse Cuthill-McKee algorithm
 * does the same job, but numbers all elements in the reverse order.
 *
 * These algorithms have one major drawback: they require a good starting
 * point, i.e. the degree of freedom index that will get a new index of
 * zero. The renumbering functions therefore allow the caller to specify such
 * an initial DoF, e.g. by exploiting knowledge of the actual topology of the
 * domain. It is also possible to give several starting indices, which may be
 * used to simulate a simple upstream numbering (by giving the inflow dofs as
 * starting values) or to make preconditioning faster (by letting the
 * Dirichlet boundary indices be starting points).
 *
 * If no starting index is given, one is chosen automatically, namely
 * one with the smallest coordination number (the coordination number
 * is the number of other dofs this dof couples with). This dof is
 * usually located on the boundary of the domain. There is, however,
 * large ambiguity in this when using the hierarchical meshes used in
 * this library, since in most cases the computational domain is not
 * approximated by tilting and deforming elements and by plugging
 * together variable numbers of elements at vertices, but rather by
 * hierarchical refinement. There is therefore a large number of dofs
 * with equal coordination numbers. The renumbering algorithms will
 * therefore not give optimal results.
 *
 * In the book of Schwarz (H.R.Schwarz: Methode der finiten Elemente),
 * it is advised to test many starting points, if possible all with
 * the smallest coordination number and also those with slightly
 * higher numbers. However, this seems only possible for meshes with
 * at most several dozen or a few hundred elements found in small
 * engineering problems of the early 1980s (the second edition was
 * published in 1984), but certainly not with those used in this
 * library, featuring several 10,000 to a few 100,000 elements.
 *
 *
 * <h4>Implementation of renumbering schemes</h4>
 *
 * The renumbering algorithms need quite a lot of memory, since they
 * have to store for each dof with which other dofs it couples. This
 * is done using a SparsityPattern object used to store the sparsity
 * pattern of matrices. It is not useful for the user to do anything
 * between distributing the dofs and renumbering, i.e. the calls to
 * DoFHandler::distribute_dofs and DoFHandler::renumber_dofs should
 * follow each other immediately. If you try to create a sparsity
 * pattern or anything else in between, these will be invalid
 * afterwards.
 *
 * The renumbering may take care of dof-to-dof couplings only induced
 * by eliminating constraints. In addition to the memory consumption
 * mentioned above, this also takes quite some computational time, but
 * it may be switched off upon calling the @p renumber_dofs
 * function. This will then give inferior results, since knots in the
 * graph (representing dofs) are not found to be neighbors even if
 * they would be after condensation.
 *
 * The renumbering algorithms work on a purely algebraic basis, due to
 * the isomorphism between the graph theoretical groundwork underlying
 * the algorithms and binary matrices (matrices of which the entries
 * are binary values) represented by the sparsity patterns. In
 * special, the algorithms do not try to exploit topological knowledge
 * (e.g. corner detection) to find appropriate starting points. This
 * way, however, they work in arbitrary space dimension.
 *
 * If you want to give starting points, you may give a list of dof
 * indices which will form the first step of the renumbering. The dofs
 * of the list will be consecutively numbered starting with zero,
 * i.e. this list is not renumbered according to the coordination
 * number of the nodes. Indices not in the allowed range are
 * deleted. If no index is allowed, the algorithm will search for its
 * own starting point.
 *
 *
 * <h4>Results of renumbering</h4>
 *
 * The renumbering schemes mentioned above do not lead to optimal
 * results.  However, after all there is no algorithm that
 * accomplishes this within reasonable time. There are situations
 * where the lack of optimality even leads to worse results than with
 * the original, crude, levelwise numbering scheme; one of these
 * examples is a mesh of four cells of which always those cells are
 * refined which are neighbors to the center (you may call this mesh a
 * `zoom in' mesh). In one such example the bandwidth was increased by
 * about 50 per cent.
 *
 * In most other cases, the bandwidth is reduced significantly. The reduction
 * is the better the less structured the grid is. With one grid where the
 * cells were refined according to a random driven algorithm, the bandwidth
 * was reduced by a factor of six.
 *
 * Using the constraint information usually leads to reductions in bandwidth
 * of 10 or 20 per cent, but may for some very unstructured grids also lead
 * to an increase. You have to weigh the decrease in your case with the time
 * spent to use the constraint information, which usually is several times
 * longer than the `pure' renumbering algorithm.
 *
 * In almost all cases, the renumbering scheme finds a corner to start with.
 * Since there is more than one corner in most grids and since even an
 * interior degree of freedom may be a better starting point, giving the
 * starting point by the user may be a viable way if you have a simple
 * scheme to derive a suitable point (e.g. by successively taking the
 * third child of the cell top left of the coarsest level, taking its
 * third vertex and the dof index thereof, if you want the top left corner
 * vertex). If you do not know beforehand what your grid will look like
 * (e.g. when using adaptive algorithms), searching a best starting point
 * may be difficult, however, and in many cases will not justify the effort.
 *
 *
 * <h3>Component-wise and block-wise numberings</h3>
 *
 * For finite elements composed of several base elements using the FESystem
 * class, or for elements which provide several components themselves, it
 * may be of interest to sort the DoF indices by component. This will then
 * bring out the block matrix structure, since otherwise the degrees of freedom
 * are numbered cell-wise without taking into account that they may belong to
 * different components. For example, one may want to sort degree of freedom for
 * a Stokes discretization so that we first get all velocities and then all
 * the pressures so that the resulting matrix naturally decomposes into a
 * $2\times 2$ system.
 *
 * This kind of numbering may be obtained by calling the
 * component_wise() function of this class. Since it does not touch
 * the order of indices within each component, it may be worthwhile to first
 * renumber using the Cuthill-McKee or a similar algorithm and
 * afterwards renumbering component-wise. This will bring out the
 * matrix structure and additionally have a good numbering within each
 * block.
 *
 * The component_wise() function allows not only to honor enumeration based on
 * vector components, but also allows to group together vector components into
 * "blocks" using a defaulted argument to the various DoFRenumber::component_wise()
 * functions (see @ref GlossComponent vs @ref GlossBlock for a description of
 * the difference). The blocks designated through this argument may, but do not
 * have to be, equal to the blocks that the finite element reports. For example,
 * a typical Stokes element would be
 * @code
 *   FESystem<dim> stokes_fe (FE_Q<dim>(2), dim,   // dim velocities
 *                            FE_Q<dim>(1), 1);    // one pressure
 * @endcode
 * This element has <code>dim+1</code> vector components and equally many
 * blocks. However, one may want to consider the velocities as one logical
 * block so that all velocity degrees of freedom are enumerated the same
 * way, independent of whether they are $x$- or $y$-velocities. This is done,
 * for example, in step-20 and step-22 as well as several other tutorial programs.
 *
 * On the other hand, if you really want to use block structure reported
 * by the finite element itself (a case that is often the case if you have
 * finite elements that have multiple vector components, e.g. the FE_RaviartThomas
 * or FE_Nedelec elements) then you can use the DoFRenumber::block_wise instead
 * of the DoFRenumbering::component_wise functions.
 *
 *
 * <h3>Cell-wise numbering</h3>
 *
 * Given an ordered vector of cells, the function cell_wise()
 * sorts the degrees of freedom such that degrees on earlier cells of
 * this vector will occur before degrees on later cells.
 *
 * This rule produces a well-defined ordering for discontinuous Galerkin
 * methods (FE_DGP, FE_DGQ). For continuous methods, we use the
 * additional rule that each degree of freedom is ordered according to
 * the first cell in the ordered vector it belongs to.
 *
 * Applications of this scheme are downstream() and
 * clock_wise_dg(). The first orders the cells according to a
 * downstream direction and then applies cell_wise().
 *
 * @note For DG elements, the internal numbering in each cell remains
 * unaffected. This cannot be guaranteed for continuous elements
 * anymore, since degrees of freedom shared with an earlier cell will
 * be accounted for by the other cell.
 *
 *
 * <h3>Random renumbering</h3>
 *
 * The random() function renumbers degrees of freedom randomly. This
 * function is probably seldom of use, except to check the dependence of
 * solvers (iterative or direct ones) on the numbering of the degrees
 * of freedom. It uses the @p random_shuffle function from the C++
 * standard library to do its work.
 *
 *
 * <h3>A comparison of reordering strategies</h3>
 *
 * As a benchmark of comparison, let us consider what the different
 * sparsity patterns produced by the various algorithms when using the
 * $Q_2^d\times Q_1$ element combination typically employed in the
 * discretization of Stokes equations, when used on the mesh obtained
 * in step-22 after one adaptive mesh refinement in
 * 3d. The space dimension together with the coupled finite element
 * leads to a rather dense system matrix with, on average around 180
 * nonzero entries per row. After applying each of the reordering
 * strategies shown below, the degrees of freedom are also sorted
 * using DoFRenumbering::component_wise into velocity and pressure
 * groups; this produces the $2\times 2$ block structure seen below
 * with the large velocity-velocity block at top left, small
 * pressure-pressure block at bottom right, and coupling blocks at top
 * right and bottom left.
 *
 * The goal of reordering strategies is to improve the
 * preconditioner. In step-22 we use a SparseILU to
 * preconditioner for the velocity-velocity block at the top left. The
 * quality of the preconditioner can then be measured by the number of
 * CG iterations required to solve a linear system with this
 * block. For some of the reordering strategies below we record this
 * number for adaptive refinement cycle 3, with 93176 degrees of
 * freedom; because we solve several linear systems with the same
 * matrix in the Schur complement, the average number of iterations is
 * reported. The lower the number the better the preconditioner and
 * consequently the better the renumbering of degrees of freedom is
 * suited for this task. We also state the run-time of the program, in
 * part determined by the number of iterations needed, for the first 4
 * cycles on one of our machines. Note that the reported times
 * correspond to the run time of the entire program, not just the
 * affected solver; if a program runs twice as fast with one
 * particular ordering than with another one, then this means that the
 * actual solver is actually several times faster.
 *
 * <table>
 * <tr>
 *   <td>
 *     @image html "reorder_sparsity_step_31_original.png"
 *   </td>
 *   <td>
 *     @image html "reorder_sparsity_step_31_random.png"
 *   </td>
 *   <td>
 *     @image html "reorder_sparsity_step_31_deal_cmk.png"
 *   </td>
 * </tr>
 * <tr>
 *   <td>
 *     Enumeration as produced by deal.II's DoFHandler::distribute_dofs function
 *     and no further reordering apart from the component-wise one.
 *
 *     With this renumbering, we needed an average of 92.2 iterations for the
 *     testcase outlined above, and a runtime of 7min53s.
 *   </td>
 *   <td>
 *     Random enumeration as produced by applying DoFRenumbering::random
 *     after calling DoFHandler::distribute_dofs. This enumeration produces
 *     nonzero entries in matrices pretty much everywhere, appearing here as
 *     an entirely unstructured matrix.
 *
 *     With this renumbering, we needed an average of 71 iterations for the
 *     testcase outlined above, and a runtime of 10min55s. The longer runtime
 *     despite less iterations compared to the default ordering may be due to
 *     the fact that computing and applying the ILU requires us to jump back
 *     and forth all through memory due to the lack of localization of
 *     matrix entries around the diagonal; this then leads to many cache
 *     misses and consequently bad timings.
 *   </td>
 *   <td>
 *     Cuthill-McKee enumeration as produced by calling the deal.II implementation
 *     of the algorithm provided by DoFRenumbering::Cuthill_McKee
 *     after DoFHandler::distribute_dofs.
 *
 *     With this renumbering, we needed an average of 57.3 iterations for the
 *     testcase outlined above, and a runtime of 6min10s.
 *   </td>
 *   </td>
 * </tr>
 *
 * <tr>
 *   <td>
 *     @image html "reorder_sparsity_step_31_boost_cmk.png"
 *   </td>
 *   <td>
 *     @image html "reorder_sparsity_step_31_boost_king.png"
 *   </td>
 *   <td>
 *     @image html "reorder_sparsity_step_31_boost_md.png"
 *   </td>
 * </tr>
 * <tr>
 *   <td>
 *     Cuthill-McKee enumeration as produced by calling the BOOST implementation
 *     of the algorithm provided by DoFRenumbering::boost::Cuthill_McKee
 *     after DoFHandler::distribute_dofs.
 *
 *     With this renumbering, we needed an average of 51.7 iterations for the
 *     testcase outlined above, and a runtime of 5min52s.
 *   </td>
 *   <td>
 *     King enumeration as produced by calling the BOOST implementation
 *     of the algorithm provided by DoFRenumbering::boost::king_ordering
 *     after DoFHandler::distribute_dofs. The sparsity pattern appears
 *     denser than with BOOST's Cuthill-McKee algorithm; however, this is
 *     only an illusion: the number of nonzero entries is the same, they are
 *     simply not as well clustered.
 *
 *     With this renumbering, we needed an average of 51.0 iterations for the
 *     testcase outlined above, and a runtime of 5min03s. Although the number
 *     of iterations is only slightly less than with BOOST's Cuthill-McKee
 *     implementation, runtime is significantly less. This, again, may be due
 *     to cache effects. As a consequence, this is the algorithm best suited
 *     to the testcase, and is in fact used in step-22.
 *   </td>
 *   <td>
 *     Minimum degree enumeration as produced by calling the BOOST implementation
 *     of the algorithm provided by DoFRenumbering::boost::minimum_degree
 *     after DoFHandler::distribute_dofs. The minimum degree algorithm does not
 *     attempt to minimize the bandwidth of a matrix but to minimize the amount
 *     of fill-in a LU decomposition would produce, i.e. the number of places in
 *     the matrix that would be occupied by elements of an LU decomposition that
 *     are not already occupied by elements of the original matrix. The resulting
 *     sparsity pattern obviously has an entirely different structure than the
 *     ones produced by algorithms trying to minimize the bandwidth.
 *
 *     With this renumbering, we needed an average of 58.9 iterations for the
 *     testcase outlined above, and a runtime of 6min11s.
 *   </td>
 * </tr>
 *
 * <tr>
 *   <td>
 *     @image html "reorder_sparsity_step_31_downstream.png"
 *   </td>
 *   <td>
 *   </td>
 *   <td>
 *   </td>
 * </tr>
 * <tr>
 *   <td>
 *     Downstream enumeration using DoFRenumbering::downstream using a
 *     direction that points diagonally through the domain.
 *
 *     With this renumbering, we needed an average of 90.5 iterations for the
 *     testcase outlined above, and a runtime of 7min05s.
 *   </td>
 *   <td>
 *   </td>
 *   <td>
 *   </td>
 * </tr>
 * </table>
 *
 *
 * <h3>Multigrid DoF numbering</h3>
 *
 * Most of the algorithms listed above also work on multigrid degree of freedom
 * numberings. Refer to the actual function declarations to get more
 * information on this.
 *
 * @ingroup dofs
 * @author Wolfgang Bangerth, Guido Kanschat, 1998, 1999, 2000, 2004, 2007, 2008
 */
namespace DoFRenumbering
{
  /**
   * Direction based comparator for
   * cell iterators: it returns @p
   * true if the center of the second
   * cell is downstream of the center
   * of the first one with respect to
   * the direction given to the
   * constructor.
   */
  template <class Iterator, int dim>
  struct CompareDownstream
  {
    /**
     * Constructor.
     */
    CompareDownstream (const Point<dim> &dir)
      :
      dir(dir)
    {}
    /**
     * Return true if c1 less c2.
     */
    bool operator () (const Iterator &c1, const Iterator &c2) const
    {
      const Point<dim> diff = c2->center() - c1->center();
      return (diff*dir > 0);
    }

  private:
    /**
     * Flow direction.
     */
    const Point<dim> dir;
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
    ComparePointwiseDownstream (const Point<dim> &dir)
      :
      dir(dir)
    {}
    /**
     * Return true if c1 less c2.
     */
    bool operator () (const std::pair<Point<dim>,types::global_dof_index> &c1,
                      const std::pair<Point<dim>,types::global_dof_index> &c2) const
    {
      const Point<dim> diff = c2.first-c1.first;
      return (diff*dir > 0 || (diff*dir==0 && c1.second<c2.second));
    }

  private:
    /**
     * Flow direction.
     */
    const Point<dim> dir;
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
    template <class DH>
    void
    Cuthill_McKee (DH                              &dof_handler,
                   const bool                       reversed_numbering = false,
                   const bool                       use_constraints    = false);

    /**
     * Computes the renumbering vector needed by the Cuthill_McKee()
     * function. Does not perform the renumbering on the DoFHandler dofs but
     * returns the renumbering vector.
     */
    template <class DH>
    void
    compute_Cuthill_McKee (std::vector<types::global_dof_index> &new_dof_indices,
                           const DH &,
                           const bool reversed_numbering = false,
                           const bool use_constraints    = false);

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
    template <class DH>
    void
    king_ordering (DH                              &dof_handler,
                   const bool                       reversed_numbering = false,
                   const bool                       use_constraints    = false);

    /**
     * Compute the renumbering for the King algorithm but do not actually
     * renumber the degrees of freedom in the DoF handler argument.
     */
    template <class DH>
    void
    compute_king_ordering (std::vector<types::global_dof_index> &new_dof_indices,
                           const DH &,
                           const bool reversed_numbering = false,
                           const bool use_constraints    = false);

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
    template <class DH>
    void
    minimum_degree (DH                              &dof_handler,
                    const bool                       reversed_numbering = false,
                    const bool                       use_constraints    = false);

    /**
     * Compute the renumbering for the minimum degree algorithm but do not
     * actually renumber the degrees of freedom in the DoF handler argument.
     */
    template <class DH>
    void
    compute_minimum_degree (std::vector<types::global_dof_index> &new_dof_indices,
                            const DH &,
                            const bool reversed_numbering = false,
                            const bool use_constraints    = false);
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
   * If the given DoFHandler uses a distributed triangulation (i.e., if
   * dof_handler.locally_owned() is not the complete index set), the
   * renumbering is performed on each processor's degrees of freedom
   * individually.
   */
  template <class DH>
  void
  Cuthill_McKee (DH                              &dof_handler,
                 const bool                       reversed_numbering = false,
                 const bool                       use_constraints    = false,
                 const std::vector<types::global_dof_index> &starting_indices   = std::vector<types::global_dof_index>());

  /**
   * Computes the renumbering vector needed by the Cuthill_McKee()
   * function. Does not perform the renumbering on the DoFHandler dofs but
   * returns the renumbering vector.
   */
  template <class DH>
  void
  compute_Cuthill_McKee (std::vector<types::global_dof_index> &new_dof_indices,
                         const DH &,
                         const bool reversed_numbering = false,
                         const bool use_constraints    = false,
                         const std::vector<types::global_dof_index> &starting_indices   = std::vector<types::global_dof_index>());

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
  template <class DH>
  void
  Cuthill_McKee (DH &dof_handler,
                 const unsigned int          level,
                 const bool                  reversed_numbering = false,
                 const std::vector<types::global_dof_index> &starting_indices   = std::vector<types::global_dof_index> ());

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
   * of the target component for dofs with component @p i in the
   * FESystem. Naming the same target component more than once is possible and
   * results in a blocking of several components into one. This is discussed
   * in step-22. If you omit this argument, the same order as given by the
   * finite element is used.
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
  *
  * @note A similar function, which renumbered all levels existed for
  * MGDoFHandler. This function was deleted. Thus, you have to call the level
  * function for each level now.
   */
  template <int dim, int spacedim>
  void
  component_wise (DoFHandler<dim,spacedim>        &dof_handler,
                  const std::vector<unsigned int> &target_component
                  = std::vector<unsigned int>());


  /**
   * Sort the degrees of freedom by component. It does the same thing as the
   * above function.
   */
  template <int dim>
  void
  component_wise (hp::DoFHandler<dim>             &dof_handler,
                  const std::vector<unsigned int> &target_component = std::vector<unsigned int> ());

  /**
   * Sort the degrees of freedom by component. It does the same thing as the
   * above function, only that it does this for one single level of a
   * multi-level discretization. The non-multigrid part of the MGDoFHandler is
   * not touched.
   */
  template <class DH>
  void
  component_wise (DH &dof_handler,
                  const unsigned int level,
                  const std::vector<unsigned int> &target_component = std::vector<unsigned int>());


  /**
   * Sort the degrees of freedom by component. It does the same thing as the
   * previous functions, but more: it renumbers not only every level of the
   * multigrid part, but also the global, i.e. non-multigrid components.
   */
  template <int dim>
  void
  component_wise (MGDoFHandler<dim>               &dof_handler,
                  const std::vector<unsigned int> &target_component = std::vector<unsigned int>());

  /**
   * Computes the renumbering vector needed by the component_wise()
   * functions. Does not perform the renumbering on the DoFHandler dofs but
   * returns the renumbering vector.
   */
  template <int dim, int spacedim, class ITERATOR, class ENDITERATOR>
  types::global_dof_index
  compute_component_wise (std::vector<types::global_dof_index> &new_dof_indices,
                          const ITERATOR &start,
                          const ENDITERATOR &end,
                          const std::vector<unsigned int> &target_component,
                          bool is_level_operation);

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
   * Sort the degrees of freedom by block. It does the same thing as the above
   * function.
   *
   * This function only succeeds if each of the elements in the
   * hp::FECollection attached to the hp::DoFHandler argument has exactly the
   * same number of blocks (see @ref GlossBlock "the glossary" for more
   * information). Note that this is not always given: while the
   * hp::FECollection class ensures that all of its elements have the same
   * number of vector components, they need not have the same number of
   * blocks. At the same time, this function here needs to match individual
   * blocks across elements and therefore requires that elements have the same
   * number of blocks and that subsequent blocks in one element have the same
   * meaning as in another element.
   */
  template <int dim>
  void
  block_wise (hp::DoFHandler<dim> &dof_handler);

  /**
   * Sort the degrees of freedom by block. It does the same thing as the above
   * function, only that it does this for one single level of a multi-level
   * discretization. The non-multigrid part of the MGDoFHandler is not
   * touched.
   */
  template <int dim>
  void
  block_wise (MGDoFHandler<dim>  &dof_handler,
              const unsigned int  level);


  /**
   * Sort the degrees of freedom by block. It does the same thing as the
   * previous functions, but more: it renumbers not only every level of the
   * multigrid part, but also the global, i.e. non-multigrid components.
   */
  template <int dim>
  void
  block_wise (MGDoFHandler<dim> &dof_handler);

  /**
   * Computes the renumbering vector needed by the block_wise()
   * functions. Does not perform the renumbering on the DoFHandler dofs but
   * returns the renumbering vector.
   */
  template <int dim, int spacedim, class ITERATOR, class ENDITERATOR>
  types::global_dof_index
  compute_block_wise (std::vector<types::global_dof_index> &new_dof_indices,
                      const ITERATOR &start,
                      const ENDITERATOR &end);

  /**
   * @}
   */

  /**
   * @name Various cell-wise numberings
   * @{
   */

  /**
   * Renumber the degrees cell by cell in hierarchical order (also known as
   * z-order). The main usage is that this guarantees the same ordering
   * independent of the number of processors involved in a parallel
   * distributed computation.
   */
  template <int dim>
  void
  hierarchical (DoFHandler<dim> &dof_handler);

  /**
   * Cell-wise renumbering. This function takes the ordered set of cells in
   * <tt>cell_order</tt>, and makes sure that all degrees of freedom in a cell
   * with higher index are behind all degrees of freedom of a cell with lower
   * index. The order inside a cell block will be the same as before this
   * renumbering.
   */
  template <class DH>
  void
  cell_wise (DH &dof_handler,
             const std::vector<typename DH::active_cell_iterator> &cell_order);

  /**
   * Computes the renumbering vector needed by the cell_wise() function. Does
   * not perform the renumbering on the DoFHandler dofs but returns the
   * renumbering vector.
   */
  template <class DH>
  void
  compute_cell_wise (std::vector<types::global_dof_index> &renumbering,
                     std::vector<types::global_dof_index> &inverse_renumbering,
                     const DH &dof_handler,
                     const std::vector<typename DH::active_cell_iterator> &cell_order);

  /**
   * Cell-wise renumbering on one level. See the other function with the same
   * name.
   */
  template <class DH>
  void
  cell_wise (DH &dof_handler,
             const unsigned int level,
             const std::vector<typename DH::level_cell_iterator> &cell_order);

  /**
   * Computes the renumbering vector needed by the cell_wise() level
   * renumbering function. Does not perform the renumbering on the DoFHandler
   * dofs but returns the renumbering vector.
   */
  template <class DH>
  void
  compute_cell_wise (std::vector<types::global_dof_index> &renumbering,
                     std::vector<types::global_dof_index> &inverse_renumbering,
                     const DH   &dof_handler,
                     const unsigned int         level,
                     const std::vector<typename DH::level_cell_iterator> &cell_order);

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
   * The cells are sorted such that the centers of higher numbers are further
   * downstream with respect to the constant vector @p direction than the
   * centers of lower numbers. Even if this yields a downstream numbering with
   * respect to the flux on the edges for fairly general grids, this might not
   * be guaranteed for all meshes.
   *
   * If the @p dof_wise_renumbering argument is set to @p false, this function
   * produces a downstream ordering of the mesh cells and calls
   * cell_wise(). Therefore, the output only makes sense for Discontinuous
   * Galerkin Finite Elements (all degrees of freedom have to be associated
   * with the interior of the cell in that case) in that case.
   *
   * If @p dof_wise_renumbering is set to @p true, the degrees of freedom are
   * renumbered based on the support point location of the individual degrees
   * of freedom (obviously, the finite element needs to define support points
   * for this to work). The numbering of points with the same position in
   * downstream location (e.g. those parallel to the flow direction, or
   * several dofs within a FESystem) will be unaffected.
   */
  template <class DH>
  void
  downstream (DH               &dof_handler,
              const Point<DH::space_dimension> &direction,
              const bool        dof_wise_renumbering = false);


  /**
   * Cell-wise downstream numbering with respect to a constant flow direction
   * on one level. See the other function with the same name.
   */
  template <class DH>
  void
  downstream (DH &dof_handler,
              const unsigned int level,
              const Point<DH::space_dimension>  &direction,
              const bool         dof_wise_renumbering = false);

  /**
   * @deprecated Use downstream() instead.
   */
  template <class DH>
  void
  downstream_dg (DH &dof,
                 const Point<DH::space_dimension> &direction) DEAL_II_DEPRECATED;

  template <class DH>
  void
  downstream_dg (DH &dof,
                 const Point<DH::space_dimension> &direction)
  {
    downstream(dof, direction);
  }


  /**
   * @deprecated Use downstream() instead.
   */
  template <class DH>
  void
  downstream_dg (DH &dof,
                 unsigned int level,
                 const Point<DH::space_dimension> &direction) DEAL_II_DEPRECATED;

  template <class DH>
  void
  downstream_dg (DH &dof,
                 unsigned int level,
                 const Point<DH::space_dimension> &direction)
  {
    downstream(dof, level, direction);
  }

  /**
   * Computes the renumbering vector needed by the downstream() function. Does
   * not perform the renumbering on the DoFHandler dofs but returns the
   * renumbering vector.
   */
  template <class DH>
  void
  compute_downstream (std::vector<types::global_dof_index> &new_dof_indices,
                      std::vector<types::global_dof_index> &reverse,
                      const DH                  &dof_handler,
                      const Point<DH::space_dimension>          &direction,
                      const bool                 dof_wise_renumbering);

  /**
   * Computes the renumbering vector needed by the downstream() function. Does
   * not perform the renumbering on the DoFHandler dofs but returns the
   * renumbering vector.
   */
  template <class DH>
  void
  compute_downstream (std::vector<types::global_dof_index> &new_dof_indices,
                      std::vector<types::global_dof_index> &reverse,
                      const DH &dof_handler,
                      const unsigned int         level,
                      const Point<DH::space_dimension>          &direction,
                      const bool                 dof_wise_renumbering);

  /**
   * Cell-wise clockwise numbering.
   *
   * This function produces a (counter)clockwise ordering of the mesh cells
   * with respect to the hub @p center and calls cell_wise().  Therefore, it
   * only works with Discontinuous Galerkin Finite Elements, i.e. all degrees
   * of freedom have to be associated with the interior of the cell.
   */
  template <class DH>
  void
  clockwise_dg (DH               &dof_handler,
                const Point<DH::space_dimension> &center,
                const bool        counter = false);

  /**
   * Cell-wise clockwise numbering on one level. See the other function with
   * the same name.
   */
  template <class DH>
  void
  clockwise_dg (DH &dof_handler,
                const unsigned int level,
                const Point<DH::space_dimension> &center,
                const bool counter = false);

  /**
   * Computes the renumbering vector needed by the clockwise_dg()
   * functions. Does not perform the renumbering on the DoFHandler dofs but
   * returns the renumbering vector.
   */
  template <class DH>
  void
  compute_clockwise_dg (std::vector<types::global_dof_index> &new_dof_indices,
                        const DH                  &dof_handler,
                        const Point<DH::space_dimension>          &center,
                        const bool                 counter);

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
  template <class DH>
  void
  sort_selected_dofs_back (DH                      &dof_handler,
                           const std::vector<bool> &selected_dofs);

  /**
   * Sort those degrees of freedom which are tagged with @p true in the @p
   * selected_dofs array on the level @p level to the back of the DoF
   * numbers. The sorting is stable, i.e. the relative order within the tagged
   * degrees of freedom is preserved, as is the relative order within the
   * untagged ones.
   *
   * @pre The @p selected_dofs array must have as many elements as the @p
   * dof_handler has degrees of freedom on the given level.
   */
  template <class DH>
  void
  sort_selected_dofs_back (DH                      &dof_handler,
                           const std::vector<bool> &selected_dofs,
                           const unsigned int       level);

  /**
   * Computes the renumbering vector needed by the sort_selected_dofs_back()
   * function. Does not perform the renumbering on the DoFHandler dofs but
   * returns the renumbering vector.
   *
   * @pre The @p selected_dofs array must have as many elements as the @p
   * dof_handler has degrees of freedom.
   */
  template <class DH>
  void
  compute_sort_selected_dofs_back (std::vector<types::global_dof_index> &new_dof_indices,
                                   const DH                  &dof_handler,
                                   const std::vector<bool>   &selected_dofs);

  /**
   * Computes the renumbering vector on each level needed by the
   * sort_selected_dofs_back() function. Does not perform the renumbering on
   * the MGDoFHandler dofs but returns the renumbering vector.
   *
   * @pre The @p selected_dofs array must have as many elements as the @p
   * dof_handler has degrees of freedom on the given level.
   */
  template <class DH>
  void
  compute_sort_selected_dofs_back (std::vector<types::global_dof_index> &new_dof_indices,
                                   const DH                  &dof_handler,
                                   const std::vector<bool>   &selected_dofs,
                                   const unsigned int         level);

  /**
   * Renumber the degrees of freedom in a random way.
   */
  template <class DH>
  void
  random (DH &dof_handler);

  /**
   * Computes the renumbering vector needed by the random() function. Does not
   * perform the renumbering on the DoFHandler dofs but returns the
   * renumbering vector.
   */
  template <class DH>
  void
  compute_random (std::vector<types::global_dof_index> &new_dof_indices,
                  const DH &dof_handler);

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
   * one, etc. This is useful when doing parallel computations after assigning
   * subdomain ids using a partitioner (see the
   * GridTools::partition_triangulation function for this).
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
  template <class DH>
  void
  subdomain_wise (DH &dof_handler);

  /**
   * Computes the renumbering vector needed by the subdomain_wise()
   * function. Does not perform the renumbering on the @p DoFHandler dofs but
   * returns the renumbering vector.
   */
  template <class DH>
  void
  compute_subdomain_wise (std::vector<types::global_dof_index> &new_dof_indices,
                          const DH                  &dof_handler);

  /**
   * @}
   */

  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcRenumberingIncomplete);
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcInvalidComponentOrder);
  /**
   * The function is only
   * implemented for Discontinuous
   * Galerkin Finite elements.
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcNotDGFEM);
}

/* ------------------------- inline functions -------------- */

#ifndef DOXYGEN
namespace DoFRenumbering
{
  template <class DH>
  void
  inline
  downstream (DH &dof,
              const Point<DH::space_dimension> &direction,
              const bool dof_wise_renumbering)
  {
    std::vector<types::global_dof_index> renumbering(dof.n_dofs());
    std::vector<types::global_dof_index> reverse(dof.n_dofs());
    compute_downstream(renumbering, reverse, dof, direction,
                       dof_wise_renumbering);

    dof.renumber_dofs(renumbering);
  }
}
#endif


DEAL_II_NAMESPACE_CLOSE

#endif
