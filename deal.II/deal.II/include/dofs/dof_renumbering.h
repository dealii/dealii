//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__dof_renumbering_h
#define __deal2__dof_renumbering_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/point.h>
#include <dofs/dof_handler.h>
#include <multigrid/mg_dof_handler.h>

#include <vector>


/**
 * Implementation of a number of renumbering algorithms for the degrees of
 * freedom on a triangulation.
 *
 * <h3>Cuthill-McKee like algorithms</h3>
 *
 * Within this class, the Cuthill-McKee algorithm is implemented. It
 * starts at a degree of freedom, searches the other DoFs for those
 * which are couple with the one we started with and numbers these in
 * a certain way. It then finds the second level of DoFs, namely those
 * that couple with those of the previous level (which were those that
 * coupled with the initial DoF) and numbers these. And so on. For the
 * details of the algorithm, especially the numbering within each
 * level, we refer the reader to the book of Schwarz (H.R.Schwarz:
 * Methode der finiten Elemente). The reverse Cuthill-McKee algorithm
 * does the same job, but numbers all elements in the reverse order.
 *
 * These algorithms have one major drawback: they require a good
 * starting point, i.e. the degree of freedom index afterwards to be
 * numbered zero. This can thus be given by the user, e.g. by
 * exploiting knowledge of the actual topology of the domain. It is
 * also possible to give several starting indices, which may be used
 * to simulate a simple upstream numbering (by giving the inflow dofs
 * as starting values) or to make preconditioning faster (by letting
 * the dirichlet boundary indices be starting points).
 *
 * If no starting index is given, one is chosen by the program, namely
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
 * On the other hand, the need to reduce the bandwidth has decreased
 * since with the mentioned number of cells, only iterative solution
 * methods are able to solve the resulting matrix systems. These,
 * however, are not so demanding with respect to the bandwidth as
 * direct solvers used for smaller problems. Things like upstream
 * numbering become much more important in recent times, so the
 * suboptimality of the renumbering algorithms is not that important
 * any more.
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
 * the original, crude, levelwise numering scheme; one of these
 * examples is a mesh of four cells of which always those cells are
 * refined which are neighbors to the center (you may call this mesh a
 * `zoom in' mesh). In one such example the bandwidth was increased by
 * about 50 per cent.
 *
 * In most other cases, the bandwith is reduced significantly. The reduction
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
 * <h3>Component-wise numbering</h3>
 *
 * For finite elements composed of several base elements using the FESystem
 * class, or for elements which provide several components themselves, it
 * may be of interest to sort the DoF indices by component. This will then
 * bring out the block matrix structure, since otherwise the degrees of freedom
 * are numbered cell-wise without taking into account that they may belong to
 * different components.
 *
 * This kind of numbering may be obtained by calling the
 * component_wise() function of this class. Since it does not touch
 * the order of indices within each, it may be worthwhile to first
 * renumber using the Cuthill-McKee or a similar algorithm and
 * afterwards renumbering component-wise. This will bring out the
 * matrix structure and additionally have a good numbering within each
 * block.
 *
 *
 * <h3>Cell-wise numbering for Discontinuous Galerkin FEM</h3>
 *
 * One advantage of DGFEM is the fact, that it yields invertible
 * blocks on the diagonal of the global matrix. these blocks are in
 * fact the cell matrices and matrix elements outside these blocks are
 * due to fluxes.
 *
 * Still, it may be necessary to apply a downstream numbering of the
 * degrees of freedom. This renumbering may only exchange whole blocks
 * and must not destroy the block structure.
 *
 * Given an ordered vector of cells, the function cell_wise_dg()
 * accomplishes this. Inside the cells, the previous ordering will be
 * preserved, so it may be useful to apply component_wise() first.
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
 * <h3>Multigrid DoF numbering</h3>
 *
 * Most algorithms also work on multigrid degree of freedom
 * numberings. Refer to the actual function declarations to get more
 * information on this.
 *
 * @ingroup dofs
 * @author Wolfgang Bangerth, Guido Kanschat, 1998, 1999, 2000, 2004
 */
class DoFRenumbering 
{
  public:
				     /**
				      * Renumber the degrees of
				      * freedom according to the
				      * Cuthill-McKee method,
				      * eventually using the reverse
				      * numbering scheme.
				      *
				      * See the general documentation
				      * of this class for details on
				      * the different methods.
				      */
    template <int dim>
    static void 
    Cuthill_McKee (DoFHandler<dim>            &dof_handler,
		   const bool                  reversed_numbering = false,
		   const bool                  use_constraints    = false,
		   const std::vector<unsigned int> &starting_indices   = std::vector<unsigned int>());

				     /**
				      * Computes the renumbering
				      * vector needed by the
				      * Cuthill_McKee() function. Does
				      * not perform the renumbering on
				      * the DoFHandler dofs but
				      * returns the renumbering
				      * vector.
				      */    
    template <int dim>
    static void
    compute_Cuthill_McKee (std::vector<unsigned int>& new_dof_indices,
			   const DoFHandler<dim>&,
			   const bool reversed_numbering = false,
			   const bool use_constraints    = false,
			   const std::vector<unsigned int> &starting_indices   = std::vector<unsigned int>());

				     /**
				      * Renumber the degrees of
				      * freedom according to the
				      * Cuthill-McKee method,
				      * eventually using the reverse
				      * numbering scheme, in this case
				      * for a multigrid numbering of
				      * degrees of freedom.
				      *
				      * You can give a triangulation
				      * level to which this function
				      * is to be applied.  Since with
				      * a level-wise numbering there
				      * are no hanging nodes, no
				      * constraints can be used, so
				      * the respective parameter of
				      * the previous function is
				      * ommitted.
				      *
				      * See the general documentation
				      * of this class for details on
				      * the different methods.
				      */
    template <int dim>
    static void 
    Cuthill_McKee (MGDoFHandler<dim>          &dof_handler,
		   const unsigned int          level,
		   const bool                  reversed_numbering = false,
		   const std::vector<unsigned int> &starting_indices   = std::vector<unsigned int> ());

				     /**
				      * Sort the degrees of freedom by
				      * vector component. The
				      * numbering within each
				      * component is not touched, so a
				      * degree of freedom with index
				      * $i$, belonging to some
				      * component, and another degree
				      * of freedom with index $j$
				      * belonging to the same
				      * component will be assigned new
				      * indices $n(i)$ and $n(j)$ with
				      * $n(i)<n(j)$ if $i<j$ and
				      * $n(i)>n(j)$ if $i>j$.
				      *
				      * You can specify that the
				      * components are ordered in a
				      * different way than suggested
				      * by the FESystem object you
				      * use. To this end, set up the
				      * vector @p target_component
				      * such that the entry at index
				      * @p i denotes the number of the
				      * target component for dofs with
				      * component @p i in the
				      * FESystem. Naming the same
				      * component more than once is
				      * possible and results in a
				      * blocking of several components
				      * into one.  If you omit this
				      * argument, the same order as
				      * given by the finite element is
				      * used.
				      *
				      * If one of the base finite
				      * elements from which the global
				      * finite element under
				      * consideration here, is a
				      * non-primitive one, i.e. its
				      * shape functions have more than
				      * one non-zero component, then
				      * it is not possible to
				      * associate these degrees of
				      * freedom with a single vector
				      * component. In this case, they
				      * are associated with the first
				      * vector component to which they
				      * belong.
				      * 
				      * For finite elements with only
				      * one component, or a single
				      * non-primitive base element,
				      * this function is the identity
				      * operation.
				      */
    template <int dim>
    static void
    component_wise (DoFHandler<dim>                 &dof_handler,
		    const std::vector<unsigned int> &target_component
		    = std::vector<unsigned int>());

				     /**
				      * Sort the degrees of freedom by
				      * component. It does the same
				      * thing as the above function,
				      * only that it does this for one
				      * level of a multi-level
				      * discretization.
				      */
    template <int dim>
    static void
    component_wise (MGDoFHandler<dim>&               dof_handler,
		    unsigned int                     level,
		    const std::vector<unsigned int>& target_component = std::vector<unsigned int>());

				     /**
				      * Computes the renumbering
				      * vector needed by the
				      * component_wise()
				      * functions. Does not perform
				      * the renumbering on the
				      * DoFHandler dofs but returns
				      * the renumbering vector.
				      */    
    template <int dim, class ITERATOR, class ENDITERATOR>
    static unsigned int
    compute_component_wise (std::vector<unsigned int>& new_dof_indices,
			    ITERATOR& start,
			    const ENDITERATOR& end,
			    const std::vector<unsigned int> &target_component);
    
				     /**
				      * Cell-wise renumbering for DG
				      * elements.  This function takes
				      * the ordered set of cells in
				      * <tt>cell_order</tt>, and makes
				      * sure that all degrees of
				      * freedom in a cell with higher
				      * index are behind all degrees
				      * of freedom of a cell with
				      * lower index. The order inside
				      * a cell bloock will be the same
				      * as before this renumbering.
				      *
				      * This function only works with
				      * Discontinuous Galerkin Finite
				      * Elements, i.e. all degrees of
				      * freedom have to be associated
				      * with the interior of the cell.
				      */
    template <int dim>
    static void
    cell_wise_dg (DoFHandler<dim> &dof_handler,
		  const std::vector<typename DoFHandler<dim>::cell_iterator> &cell_order);

				     /**
				      * Computes the renumbering
				      * vector needed by the
				      * cell_wise_dg() function. Does
				      * not perform the renumbering on
				      * the DoFHandler dofs but
				      * returns the renumbering
				      * vector.
				      */
    template <int dim>
    static void
    compute_cell_wise_dg (std::vector<unsigned int> &new_dof_indices,
			  const DoFHandler<dim> &dof_handler,
			  const std::vector<typename DoFHandler<dim>::cell_iterator> &cell_order);

				     /**
				      * Cell-wise renumbering on one
				      * level for DG elements. See the
				      * other function with the same
				      * name.
				      */
    template <int dim>
    static void
    cell_wise_dg (MGDoFHandler<dim>   &dof_handler,
		  const unsigned int   level,
		  const std::vector<typename MGDoFHandler<dim>::cell_iterator> &cell_order);
    

				     /**
				      * Cell-wise downstream numbering
				      * with respect to a constant
				      * flow direction.
				      *
				      * The cells are sorted such that
				      * the centers of higher numbers
				      * are further downstream with
				      * respect to the constant vector
				      * @p direction than the centers
				      * of lower numbers. Even if this
				      * yields a downstream numbering
				      * with respect to the flux on
				      * the edges for fairly general
				      * grids, this might not be
				      * guaranteed for all meshes.
				      *
				      * This function produces a
				      * downstream ordering of the
				      * mesh cells and calls
				      * cell_wise_dg().  Therefore, it
				      * only works with Discontinuous
				      * Galerkin Finite Elements,
				      * i.e. all degrees of freedom
				      * have to be associated with the
				      * interior of the cell.
				      */
    template <int dim>
    static void
    downstream_dg (DoFHandler<dim>  &dof_handler,
		   const Point<dim> &direction);

				     /**
				      * Cell-wise downstream numbering
				      * with respect to a constant
				      * flow direction on one
				      * level. See the other function
				      * with the same name.
				      */
    template <int dim>
    static void
    downstream_dg (MGDoFHandler<dim>  &dof_handler,
		   const unsigned int level,
		   const Point<dim> &direction);

    				     /**
				      * Computes the renumbering
				      * vector needed by the
				      * downstream_dg() function. Does
				      * not perform the renumbering on
				      * the DoFHandler dofs but
				      * returns the renumbering
				      * vector.
				      */
    template <int dim>
    static void
    compute_downstream_dg (std::vector<unsigned int> &new_dof_indices,
			   const DoFHandler<dim>  &dof_handler,
			   const Point<dim> &direction);

//TODO:[GK] Documentation!    
				     /**
				      * Cell-wise clockwise numbering.
				      *
				      * This function produces a
				      * (counter)clockwise ordering of
				      * the mesh cells with respect to
				      * the hub @p center and calls
				      * cell_wise_dg().  Therefore, it
				      * only works with Discontinuous
				      * Galerkin Finite Elements,
				      * i.e. all degrees of freedom
				      * have to be associated with the
				      * interior of the cell.
				      */
    template <int dim>
    static void
    clockwise_dg (DoFHandler<dim>  &dof_handler,
		  const Point<dim> &center,
		  const bool counter = false);

				     /**
				      * Cell-wise clockwise numbering
				      * on one level. See the other
				      * function with the same name.
				      */
    template <int dim>
    static void
    clockwise_dg (MGDoFHandler<dim>  &dof_handler,
		  const unsigned int level,
		  const Point<dim> &center,
		  const bool counter = false);

    				     /**
				      * Computes the renumbering
				      * vector needed by the
				      * clockwise_dg() functions. Does
				      * not perform the renumbering on
				      * the DoFHandler dofs but
				      * returns the renumbering
				      * vector.
				      */
    template <int dim>
    static void
    compute_clockwise_dg (std::vector<unsigned int> &new_dof_indices,
			  const DoFHandler<dim>  &dof_handler,
			  const Point<dim> &center,
			  const bool counter);

				     /**
				      * Sort those degrees of freedom
				      * which are tagged with @p true
				      * in the @p selected_dofs array
				      * to the back of the DoF
				      * numbers. The sorting is
				      * stable, i.e. the relative
				      * order within the tagged
				      * degrees of freedom is
				      * preserved, as is the relative
				      * order within the untagged
				      * ones.
				      */
    template <int dim>
    static void
    sort_selected_dofs_back (DoFHandler<dim>         &dof_handler,
			     const std::vector<bool> &selected_dofs);

    				     /**
				      * Computes the renumbering
				      * vector needed by the
				      * sort_selected_dofs_back()
				      * function. Does not perform the
				      * renumbering on the DoFHandler
				      * dofs but returns the
				      * renumbering vector.
				      */
    template <int dim>
    static void
    compute_sort_selected_dofs_back (std::vector<unsigned int> &new_dof_indices,
				     const DoFHandler<dim>&,
				     const std::vector<bool>&selected_dofs);

				     /**
				      * Renumber the degrees of
				      * freedom in a random way.
				      */
    template <int dim>
    static void
    random (DoFHandler<dim> &dof_handler);

    				     /**
				      * Computes the renumbering
				      * vector needed by the random()
				      * function. Does not perform the
				      * renumbering on the DoFHandler
				      * dofs but returns the
				      * renumbering vector.
				      */   
    template <int dim>
    static void
    compute_random (std::vector<unsigned int> &new_dof_indices,
		    const DoFHandler<dim> &dof_handler);

				     /**
				      * Renumber the degrees of
				      * freedom such that they are
				      * associated with the subdomain
				      * id of the cells they are
				      * living on, i.e. first all
				      * degrees of freedom that belong
				      * to cells with subdomain zero,
				      * then all with subdomain one,
				      * etc. This is useful when doing
				      * parallel computations after
				      * assigning subdomain ids using
				      * a partitioner (see the
				      * GridTools::partition_triangulation
				      * function for this).
				      *
				      * Note that degrees of freedom
				      * associated with faces, edges,
				      * and vertices may be associated
				      * with multiple subdomains if
				      * they are sitting on partition
				      * boundaries. It would therefore
				      * be undefined with which
				      * subdomain they have to be
				      * associated. For this, we use
				      * what we get from the
				      * DoFTools::get_subdomain_association
				      * function.
				      *
				      * The algorithm is stable, i.e. if
                                      * two dofs i,j have <tt>i<j</tt> and belong
                                      * to the same subdomain, then they
                                      * will be in this order also after
                                      * reordering.
				      */
    template <int dim>
    static void
    subdomain_wise (DoFHandler<dim> &dof_handler);

    				     /**
				      * Computes the renumbering
				      * vector needed by the
				      * subdomain_wise()
				      * function. Does not perform the
				      * renumbering on the @p
				      * DoFHandler dofs but returns
				      * the renumbering vector.
				      */   
    template <int dim>
    static void
    compute_subdomain_wise (std::vector<unsigned int> &new_dof_indices,
                            const DoFHandler<dim>     &dof_handler);
    
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
};


#endif
