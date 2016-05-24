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

#ifndef dealii__dof_tools_h
#define dealii__dof_tools_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/hp/dof_handler.h>

#include <vector>
#include <set>
#include <map>

DEAL_II_NAMESPACE_OPEN

class BlockMask;
template <int dim, typename Number> class Function;
template <int dim, int spacedim> class FiniteElement;
namespace hp
{
  template <int dim, int spacedim> class MappingCollection;
  template <int dim, int spacedim> class FECollection;
}
template <class MeshType> class InterGridMap;
template <int dim, int spacedim> class Mapping;
class SparsityPattern;
template <int dim, class T> class Table;
template <typename Number> class Vector;

namespace GridTools
{
  template <typename CellIterator> struct PeriodicFacePair;
}

//TODO: map_support_points_to_dofs should generate a multimap, rather than just a map, since several dofs may be located at the same support point

/**
 * This is a collection of functions operating on, and manipulating the
 * numbers of degrees of freedom. The documentation of the member functions
 * will provide more information, but for functions that exist in multiple
 * versions, there are sections in this global documentation stating some
 * commonalities.
 *
 * <h3>Setting up sparsity patterns</h3>
 *
 * When assembling system matrices, the entries are usually of the form
 * $a_{ij} = a(\phi_i, \phi_j)$, where $a$ is a bilinear functional, often an
 * integral. When using sparse matrices, we therefore only need to reserve
 * space for those $a_{ij}$ only, which are nonzero, which is the same as to
 * say that the basis functions $\phi_i$ and $\phi_j$ have a nonempty
 * intersection of their support. Since the support of basis functions is
 * bound only on cells on which they are located or to which they are
 * adjacent, to determine the sparsity pattern it is sufficient to loop over
 * all cells and connect all basis functions on each cell with all other basis
 * functions on that cell.  There may be finite elements for which not all
 * basis functions on a cell connect with each other, but no use of this case
 * is made since no examples where this occurs are known to the author.
 *
 *
 * <h3>DoF numberings on boundaries</h3>
 *
 * When projecting the traces of functions to the boundary or parts thereof,
 * one needs to build matrices and vectors that act only on those degrees of
 * freedom that are located on the boundary, rather than on all degrees of
 * freedom. One could do that by simply building matrices in which the entries
 * for all interior DoFs are zero, but such matrices are always very rank
 * deficient and not very practical to work with.
 *
 * What is needed instead in this case is a numbering of the boundary degrees
 * of freedom, i.e. we should enumerate all the degrees of freedom that are
 * sitting on the boundary, and exclude all other (interior) degrees of
 * freedom. The map_dof_to_boundary_indices() function does exactly this: it
 * provides a vector with as many entries as there are degrees of freedom on
 * the whole domain, with each entry being the number in the numbering of the
 * boundary or DoFHandler::invalid_dof_index if the dof is not on the
 * boundary.
 *
 * With this vector, one can get, for any given degree of freedom, a unique
 * number among those DoFs that sit on the boundary; or, if your DoF was
 * interior to the domain, the result would be DoFHandler::invalid_dof_index.
 * We need this mapping, for example, to build the mass matrix on the boundary
 * (for this, see make_boundary_sparsity_pattern() function, the corresponding
 * section below, as well as the MatrixCreator namespace documentation).
 *
 * Actually, there are two map_dof_to_boundary_indices() functions, one
 * producing a numbering for all boundary degrees of freedom and one producing
 * a numbering for only parts of the boundary, namely those parts for which
 * the boundary indicator is listed in a set of indicators given to the
 * function. The latter case is needed if, for example, we would only want to
 * project the boundary values for the Dirichlet part of the boundary. You
 * then give the function a list of boundary indicators referring to Dirichlet
 * parts on which the projection is to be performed. The parts of the boundary
 * on which you want to project need not be contiguous; however, it is not
 * guaranteed that the indices of each of the boundary parts are continuous,
 * i.e. the indices of degrees of freedom on different parts may be
 * intermixed.
 *
 * Degrees of freedom on the boundary but not on one of the specified boundary
 * parts are given the index DoFHandler::invalid_dof_index, as if they were in
 * the interior. If no boundary indicator was given or if no face of a cell
 * has a boundary indicator contained in the given list, the vector of new
 * indices consists solely of DoFHandler::invalid_dof_index.
 *
 * (As a side note, for corner cases: The question what a degree of freedom on
 * the boundary is, is not so easy.  It should really be a degree of freedom
 * of which the respective basis function has nonzero values on the boundary.
 * At least for Lagrange elements this definition is equal to the statement
 * that the off-point, or what deal.II calls support_point, of the shape
 * function, i.e. the point where the function assumes its nominal value (for
 * Lagrange elements this is the point where it has the function value 1), is
 * located on the boundary. We do not check this directly, the criterion is
 * rather defined through the information the finite element class gives: the
 * FiniteElement class defines the numbers of basis functions per vertex, per
 * line, and so on and the basis functions are numbered after this
 * information; a basis function is to be considered to be on the face of a
 * cell (and thus on the boundary if the cell is at the boundary) according to
 * it belonging to a vertex, line, etc but not to the interior of the cell.
 * The finite element uses the same cell-wise numbering so that we can say
 * that if a degree of freedom was numbered as one of the dofs on lines, we
 * assume that it is located on the line. Where the off-point actually is, is
 * a secret of the finite element (well, you can ask it, but we don't do it
 * here) and not relevant in this context.)
 *
 *
 * <h3>Setting up sparsity patterns for boundary matrices</h3>
 *
 * In some cases, one wants to only work with DoFs that sit on the boundary.
 * One application is, for example, if rather than interpolating non-
 * homogenous boundary values, one would like to project them. For this, we
 * need two things: a way to identify nodes that are located on (parts of) the
 * boundary, and a way to build matrices out of only degrees of freedom that
 * are on the boundary (i.e. much smaller matrices, in which we do not even
 * build the large zero block that stems from the fact that most degrees of
 * freedom have no support on the boundary of the domain). The first of these
 * tasks is done by the map_dof_to_boundary_indices() function (described
 * above).
 *
 * The second part requires us first to build a sparsity pattern for the
 * couplings between boundary nodes, and then to actually build the components
 * of this matrix. While actually computing the entries of these small
 * boundary matrices is discussed in the MatrixCreator namespace, the creation
 * of the sparsity pattern is done by the create_boundary_sparsity_pattern()
 * function. For its work, it needs to have a numbering of all those degrees
 * of freedom that are on those parts of the boundary that we are interested
 * in. You can get this from the map_dof_to_boundary_indices() function. It
 * then builds the sparsity pattern corresponding to integrals like
 * $\int_\Gamma \varphi_{b2d(i)} \varphi_{b2d(j)} dx$, where $i$ and $j$ are
 * indices into the matrix, and $b2d(i)$ is the global DoF number of a degree
 * of freedom sitting on a boundary (i.e., $b2d$ is the inverse of the mapping
 * returned by map_dof_to_boundary_indices() function).
 *
 *
 * @ingroup dofs
 * @author Wolfgang Bangerth, Guido Kanschat and others
 */
namespace DoFTools
{
  /**
   * The flags used in tables by certain <tt>make_*_pattern</tt> functions to
   * describe whether two components of the solution couple in the bilinear
   * forms corresponding to cell or face terms. An example of using these
   * flags is shown in the introduction of step-46.
   *
   * In the descriptions of the individual elements below, remember that these
   * flags are used as elements of tables of size FiniteElement::n_components
   * times FiniteElement::n_components where each element indicates whether
   * two components do or do not couple.
   */
  enum Coupling
  {
    /**
     * Two components do not couple.
     */
    none,
    /**
     * Two components do couple.
     */
    always,
    /**
     * Two components couple only if their shape functions are both nonzero on
     * a given face. This flag is only used when computing integrals over
     * faces of cells, e.g., in DoFTools::make_flux_sparsity_pattern().
     */
    nonzero
  };

  /**
   * @name Functions to support code that generically uses both DoFHandler and hp::DoFHandler
   * @{
   */
  /**
   * Maximal number of degrees of freedom on a cell.
   *
   * @relates DoFHandler
   */
  template <int dim, int spacedim>
  unsigned int
  max_dofs_per_cell (const DoFHandler<dim,spacedim> &dh);

  /**
   * Maximal number of degrees of freedom on a cell.
   *
   * @relates hp::DoFHandler
   */
  template <int dim, int spacedim>
  unsigned int
  max_dofs_per_cell (const hp::DoFHandler<dim,spacedim> &dh);


  /**
   * Maximal number of degrees of freedom on a face.
   *
   * This function exists for both non-hp and hp DoFHandlers, to allow for a
   * uniform interface to query this property.
   *
   * @relates DoFHandler
   */
  template <int dim, int spacedim>
  unsigned int
  max_dofs_per_face (const DoFHandler<dim,spacedim> &dh);

  /**
   * Maximal number of degrees of freedom on a face.
   *
   * This function exists for both non-hp and hp DoFHandlers, to allow for a
   * uniform interface to query this property.
   *
   * @relates hp::DoFHandler
   */
  template <int dim, int spacedim>
  unsigned int
  max_dofs_per_face (const hp::DoFHandler<dim,spacedim> &dh);

  /**
   * Maximal number of degrees of freedom on a vertex.
   *
   * This function exists for both non-hp and hp DoFHandlers, to allow for a
   * uniform interface to query this property.
   *
   * @relates DoFHandler
   */
  template <int dim, int spacedim>
  unsigned int
  max_dofs_per_vertex (const DoFHandler<dim,spacedim> &dh);

  /**
   * Maximal number of degrees of freedom on a vertex.
   *
   * This function exists for both non-hp and hp DoFHandlers, to allow for a
   * uniform interface to query this property.
   *
   * @relates hp::DoFHandler
   */
  template <int dim, int spacedim>
  unsigned int
  max_dofs_per_vertex (const hp::DoFHandler<dim,spacedim> &dh);

  /**
   * Number of vector components in the finite element object used by this
   * DoFHandler.
   *
   * This function exists for both non-hp and hp DoFHandlers, to allow for a
   * uniform interface to query this property.
   *
   * @relates DoFHandler
   */
  template <int dim, int spacedim>
  unsigned int
  n_components (const DoFHandler<dim,spacedim> &dh);

  /**
   * Number of vector components in the finite element object used by this
   * DoFHandler.
   *
   * This function exists for both non-hp and hp DoFHandlers, to allow for a
   * uniform interface to query this property.
   *
   * @relates hp::DoFHandler
   */
  template <int dim, int spacedim>
  unsigned int
  n_components (const hp::DoFHandler<dim,spacedim> &dh);

  /**
   * Find out whether the FiniteElement used by this DoFHandler is primitive
   * or not.
   *
   * This function exists for both non-hp and hp DoFHandlers, to allow for a
   * uniform interface to query this property.
   *
   * @relates DoFHandler
   */
  template <int dim, int spacedim>
  bool
  fe_is_primitive (const DoFHandler<dim,spacedim> &dh);

  /**
   * Find out whether the FiniteElement used by this DoFHandler is primitive
   * or not.
   *
   * This function exists for both non-hp and hp DoFHandlers, to allow for a
   * uniform interface to query this property.
   *
   * @relates hp::DoFHandler
   */
  template <int dim, int spacedim>
  bool
  fe_is_primitive (const hp::DoFHandler<dim,spacedim> &dh);

  /**
   * @}
   */

  /**
   * @name Sparsity pattern generation
   * @{
   */

  /**
   * Compute which entries of a matrix built on the given @p dof_handler may
   * possibly be nonzero, and create a sparsity pattern object that represents
   * these nonzero locations.
   *
   * This function computes the possible positions of non-zero entries in the
   * global system matrix by <i>simulating</i> which entries one would write
   * to during the actual assembly of a matrix. For this, the function assumes
   * that each finite element basis function is non-zero on a cell only if its
   * degree of freedom is associated with the interior, a face, an edge or a
   * vertex of this cell.  As a result, a matrix entry $A_{ij}$ that is
   * computed from two basis functions $\varphi_i$ and $\varphi_j$ with
   * (global) indices $i$ and $j$ (for example, using a bilinear form
   * $A_{ij}=a(\varphi_i,\varphi_j)$) can be non-zero only if these shape
   * functions correspond to degrees of freedom that are defined on at least
   * one common cell. Therefore, this function just loops over all cells,
   * figures out the global indices of all degrees of freedom, and presumes
   * that all matrix entries that couple any of these indices will result in a
   * nonzero matrix entry. These will then be added to the sparsity pattern.
   * As this process of generating the sparsity pattern does not take into
   * account the equation to be solved later on, the resulting sparsity
   * pattern is symmetric.
   *
   * This algorithm makes no distinction between shape functions on each cell,
   * i.e., it simply couples all degrees of freedom on a cell with all other
   * degrees of freedom on a cell. This is often the case, and always a safe
   * assumption. However, if you know something about the structure of your
   * operator and that it does not couple certain shape functions with certain
   * test functions, then you can get a sparser sparsity pattern by calling a
   * variant of the current function described below that allows to specify
   * which vector components couple with which other vector components.
   *
   * The method described above lives on the assumption that coupling between
   * degrees of freedom only happens if shape functions overlap on at least
   * one cell. This is the case with most usual finite element formulations
   * involving conforming elements. However, for formulations such as the
   * Discontinuous Galerkin finite element method, the bilinear form contains
   * terms on interfaces between cells that couple shape functions that live
   * on one cell with shape functions that live on a neighboring cell. The
   * current function would not see these couplings, and would consequently
   * not allocate entries in the sparsity pattern. You would then get into
   * trouble during matrix assembly because you try to write into matrix
   * entries for which no space has been allocated in the sparsity pattern.
   * This can be avoided by calling the DoFTools::make_flux_sparsity_pattern()
   * function instead, which takes into account coupling between degrees of
   * freedom on neighboring cells.
   *
   * There are other situations where bilinear forms contain non-local terms,
   * for example in treating integral equations. These require different
   * methods for building the sparsity patterns that depend on the exact
   * formulation of the problem. You will have to do this yourself then.
   *
   * @param[in] dof_handler The DoFHandler or hp::DoFHandler object that
   * describes which degrees of freedom live on which cells.
   *
   * @param[out] sparsity_pattern The sparsity pattern to be filled with
   * entries.
   *
   * @param[in] constraints The process for generating entries described above
   * is purely local to each cell. Consequently, the sparsity pattern does not
   * provide for matrix entries that will only be written into during the
   * elimination of hanging nodes or other constraints. They have to be taken
   * care of by a subsequent call to ConstraintMatrix::condense().
   * Alternatively, the constraints on degrees of freedom can already be taken
   * into account at the time of creating the sparsity pattern. For this, pass
   * the ConstraintMatrix object as the third argument to the current
   * function. No call to ConstraintMatrix::condense() is then necessary. This
   * process is explained in step-6, step-27, and other tutorial programs.
   *
   * @param[in] keep_constrained_dofs In case the constraints are already
   * taken care of in this function by passing in a ConstraintMatrix object,
   * it is possible to abandon some off-diagonal entries in the sparsity
   * pattern if these entries will also not be written into during the actual
   * assembly of the matrix this sparsity pattern later serves. Specifically,
   * when using an assembly method that uses
   * ConstraintMatrix::distribute_local_to_global(), no entries will ever be
   * written into those matrix rows or columns that correspond to constrained
   * degrees of freedom. In such cases, you can set the argument @p
   * keep_constrained_dofs to @p false to avoid allocating these entries in
   * the sparsity pattern.
   *
   * @param[in] subdomain_id If specified, the sparsity pattern is built only
   * on cells that have a subdomain_id equal to the given argument. This is
   * useful in parallel contexts where the matrix and sparsity pattern (for
   * example a TrilinosWrappers::SparsityPattern) may be distributed and not
   * every MPI process needs to build the entire sparsity pattern; in that
   * case, it is sufficient if every process only builds that part of the
   * sparsity pattern that corresponds to the subdomain_id for which it is
   * responsible. This feature is used in step-32. (This argument is not
   * usually needed for objects of type parallel::distributed::Triangulation
   * because the current function only loops over locally owned cells anyway;
   * thus, this argument typically only makes sense if you want to use the
   * subdomain_id for anything other than indicating which processor owns a
   * cell, for example which geometric component of the domain a cell belongs
   * to.)
   *
   * @note The actual type of the sparsity pattern may be SparsityPattern,
   * DynamicSparsityPattern, BlockSparsityPattern,
   * BlockDynamicSparsityPattern, or any other class that satisfies similar
   * requirements. It is assumed that the size of the sparsity pattern matches
   * the number of degrees of freedom and that enough unused nonzero entries
   * are left to fill the sparsity pattern if the sparsity pattern is of
   * "static" kind (see
   * @ref Sparsity
   * for more information on what this means). The nonzero entries generated
   * by this function are added to possible previous content of the object,
   * i.e., previously added entries are not removed.
   *
   * @note If the sparsity pattern is represented by an object of type
   * SparsityPattern (as opposed to, for example, DynamicSparsityPattern), you
   * need to remember using SparsityPattern::compress() after generating the
   * pattern.
   *
   * @ingroup constraints
   */
  template <typename DoFHandlerType, typename SparsityPatternType>
  void
  make_sparsity_pattern (const DoFHandlerType      &dof_handler,
                         SparsityPatternType       &sparsity_pattern,
                         const ConstraintMatrix    &constraints           = ConstraintMatrix(),
                         const bool                 keep_constrained_dofs = true,
                         const types::subdomain_id  subdomain_id          = numbers::invalid_subdomain_id);

  /**
   * Compute which entries of a matrix built on the given @p dof_handler may
   * possibly be nonzero, and create a sparsity pattern object that represents
   * these nonzero locations.
   *
   * This function is a simple variation on the previous
   * make_sparsity_pattern() function (see there for a description of all of
   * the common arguments), but it provides functionality for vector finite
   * elements that allows to be more specific about which variables couple in
   * which equation.
   *
   * For example, if you wanted to solve the Stokes equations,
   *
   * @f{align*}{
   * -\Delta \mathbf u + \nabla p &= 0,\\ \text{div}\ u &= 0
   * @f}
   *
   * in two space dimensions, using stable Q2/Q1 mixed elements (using the
   * FESystem class), then you don't want all degrees of freedom to couple in
   * each equation. More specifically, in the first equation, only $u_x$ and
   * $p$ appear; in the second equation, only $u_y$ and $p$ appear; and in the
   * third equation, only $u_x$ and $u_y$ appear. (Note that this discussion
   * only talks about vector components of the solution variable and the
   * different equation, and has nothing to do with degrees of freedom, or in
   * fact with any kind of discretization.) We can describe this by the
   * following pattern of "couplings":
   *
   * @f[
   * \left[
   * \begin{array}{ccc}
   *   1 & 0 & 1 \\
   *   0 & 1 & 1 \\
   *   1 & 1 & 0
   * \end{array}
   * \right]
   * @f]
   *
   * where "1" indicates that two variables (i.e., vector components of the
   * FESystem) couple in the respective equation, and a "0" means no coupling.
   * These zeros imply that upon discretization via a standard finite element
   * formulation, we will not write entries into the matrix that, for example,
   * couple pressure test functions with pressure shape functions (and similar
   * for the other zeros above). It is then a waste to allocate memory for
   * these entries in the matrix and the sparsity pattern, and you can avoid
   * this by creating a mask such as the one above that describes this to the
   * (current) function that computes the sparsity pattern. As stated above,
   * the mask shown above refers to components of the composed FESystem,
   * rather than to degrees of freedom or shape functions.
   *
   * This function is designed to accept a coupling pattern, like the one
   * shown above, through the @p couplings parameter, which contains values of
   * type #Coupling. It builds the matrix structure just like the previous
   * function, but does not create matrix elements if not specified by the
   * coupling pattern. If the couplings are symmetric, then so will be the
   * resulting sparsity pattern.
   *
   * There is a complication if some or all of the shape functions of the
   * finite element in use are non-zero in more than one component (in deal.II
   * speak: they are
   * @ref GlossPrimitive "non-primitive finite elements").
   * In this case, the coupling element corresponding to the first non-zero
   * component is taken and additional ones for this component are ignored.
   *
   * @ingroup constraints
   */
  template <typename DoFHandlerType, typename SparsityPatternType>
  void
  make_sparsity_pattern (const DoFHandlerType     &dof_handler,
                         const Table<2, Coupling> &coupling,
                         SparsityPatternType      &sparsity_pattern,
                         const ConstraintMatrix   &constraints           = ConstraintMatrix(),
                         const bool                keep_constrained_dofs = true,
                         const types::subdomain_id subdomain_id         = numbers::invalid_subdomain_id);

  /**
   * Construct a sparsity pattern that allows coupling degrees of freedom on
   * two different but related meshes.
   *
   * The idea is that if the two given DoFHandler objects correspond to two
   * different meshes (and potentially to different finite elements used on
   * these cells), but that if the two triangulations they are based on are
   * derived from the same coarse mesh through hierarchical refinement, then
   * one may set up a problem where one would like to test shape functions
   * from one mesh against the shape functions from another mesh. In
   * particular, this means that shape functions from a cell on the first mesh
   * are tested against those on the second cell that are located on the
   * corresponding cell; this correspondence is something that the
   * IntergridMap class can determine.
   *
   * This function then constructs a sparsity pattern for which the degrees of
   * freedom that represent the rows come from the first given DoFHandler,
   * whereas the ones that correspond to columns come from the second
   * DoFHandler.
   */
  template <typename DoFHandlerType, typename SparsityPatternType>
  void
  make_sparsity_pattern (const DoFHandlerType &dof_row,
                         const DoFHandlerType &dof_col,
                         SparsityPatternType  &sparsity);

  /**
   * Compute which entries of a matrix built on the given @p dof_handler may
   * possibly be nonzero, and create a sparsity pattern object that represents
   * these nonzero locations. This function is a variation of the
   * make_sparsity_pattern() functions above in that it assumes that the
   * bilinear form you want to use to generate the matrix also contains terms
   * that integrate over the <i>faces</i> between cells (i.e., it contains
   * "fluxes" between cells, explaining the name of the function).
   *
   * This function is useful for Discontinuous Galerkin methods where the
   * standard make_sparsity_pattern() function would only create nonzero
   * entries for all degrees of freedom on one cell coupling to all other
   * degrees of freedom on the same cell; however, in DG methods, all or some
   * degrees of freedom on each cell also couple to the degrees of freedom on
   * other cells connected to the current one by a common face. The current
   * function also creates the nonzero entries in the matrix resulting from
   * these additional couplings. In other words, this function computes a
   * strict super-set of nonzero entries compared to the work done by
   * make_sparsity_pattern().
   *
   * @param[in] dof_handler The DoFHandler or hp::DoFHandler object that
   * describes which degrees of freedom live on which cells.
   *
   * @param[out] sparsity_pattern The sparsity pattern to be filled with
   * entries.
   *
   * @note The actual type of the sparsity pattern may be SparsityPattern,
   * DynamicSparsityPattern, BlockSparsityPattern,
   * BlockDynamicSparsityPattern, or any other class that satisfies similar
   * requirements. It is assumed that the size of the sparsity pattern matches
   * the number of degrees of freedom and that enough unused nonzero entries
   * are left to fill the sparsity pattern if the sparsity pattern is of
   * "static" kind (see
   * @ref Sparsity
   * for more information on what this means). The nonzero entries generated
   * by this function are added to possible previous content of the object,
   * i.e., previously added entries are not removed.
   *
   * @note If the sparsity pattern is represented by an object of type
   * SparsityPattern (as opposed to, for example, DynamicSparsityPattern), you
   * need to remember using SparsityPattern::compress() after generating the
   * pattern.
   *
   * @ingroup constraints
   */
  template<typename DoFHandlerType, typename SparsityPatternType>
  void
  make_flux_sparsity_pattern (const DoFHandlerType &dof_handler,
                              SparsityPatternType  &sparsity_pattern);

  /**
   * This function does essentially the same as the other
   * make_flux_sparsity_pattern() function but allows the specification of a
   * number of additional arguments. These carry the same meaning as discussed
   * in the first make_sparsity_pattern() function above.
   *
   * @ingroup constraints
   */
  template<typename DoFHandlerType, typename SparsityPatternType>
  void
  make_flux_sparsity_pattern (const DoFHandlerType      &dof_handler,
                              SparsityPatternType       &sparsity_pattern,
                              const ConstraintMatrix    &constraints,
                              const bool                 keep_constrained_dofs = true,
                              const types::subdomain_id  subdomain_id          = numbers::invalid_unsigned_int);

  /**
   * This function does essentially the same as the other
   * make_flux_sparsity_pattern() function but allows the specification of
   * coupling matrices that state which components of the solution variable
   * couple in each of the equations you are discretizing. This works in
   * complete analogy as discussed in the second make_sparsity_pattern()
   * function above.
   *
   * In fact, this function takes two such masks, one describing which
   * variables couple with each other in the cell integrals that make up your
   * bilinear form, and which variables coupld with each other in the face
   * integrals. If you passed masks consisting of only 1s to both of these,
   * then you would get the same sparsity pattern as if you had called the
   * first of the make_sparsity_pattern() functions above. By setting some of
   * the entries of these masks to zeros, you can get a sparser sparsity
   * pattern.
   *
   * @ingroup constraints
   */
  template <typename DoFHandlerType, typename SparsityPatternType>
  void
  make_flux_sparsity_pattern (const DoFHandlerType    &dof,
                              SparsityPatternType     &sparsity,
                              const Table<2,Coupling> &cell_integrals_mask,
                              const Table<2,Coupling> &face_integrals_mask);

  /**
   * Create the sparsity pattern for boundary matrices. See the general
   * documentation of this class for more information.
   *
   * The function does essentially what the other make_sparsity_pattern()
   * functions do, but assumes that the bilinear form that is used to build
   * the matrix does not consist of domain integrals, but only of integrals
   * over the boundary of the domain.
   */
  template <typename DoFHandlerType, typename SparsityPatternType>
  void
  make_boundary_sparsity_pattern (const DoFHandlerType                       &dof,
                                  const std::vector<types::global_dof_index> &dof_to_boundary_mapping,
                                  SparsityPatternType                        &sparsity_pattern);

  /**
   * This function is a variation of the previous
   * make_boundary_sparsity_pattern() function in which we assume that the
   * boundary integrals that will give rise to the matrix extends only over
   * those parts of the boundary whose boundary indicators are listed in the
   * @p boundary_ids argument to this function.
   *
   * This function could have been written by passing a @p set of boundary_id
   * numbers. However, most of the functions throughout deal.II dealing with
   * boundary indicators take a mapping of boundary indicators and the
   * corresponding boundary function, i.e., a std::map<types::boundary_id, const Function<spacedim,number>*> argument.
   * Correspondingly, this function does the same, though the actual boundary
   * function is ignored here. (Consequently, if you don't have any such
   * boundary functions, just create a map with the boundary indicators you
   * want and set the function pointers to null pointers).
   */
  template <typename DoFHandlerType, typename SparsityPatternType, typename number>
  void
  make_boundary_sparsity_pattern
  (const DoFHandlerType                                              &dof,
   const std::map<types::boundary_id, const Function<DoFHandlerType::space_dimension,number>*> &boundary_ids,
   const std::vector<types::global_dof_index>                        &dof_to_boundary_mapping,
   SparsityPatternType                                               &sparsity);

  /**
   * @}
   */

  /**
   * @name Hanging nodes and other constraints
   * @{
   */

  /**
   * Compute the constraints resulting from the presence of hanging nodes.
   * Hanging nodes are best explained using a small picture:
   *
   * @image html hanging_nodes.png
   *
   * In order to make a finite element function globally continuous, we have
   * to make sure that the dark red nodes have values that are compatible with
   * the adjacent yellow nodes, so that the function has no jump when coming
   * from the small cells to the large one at the top right. We therefore have
   * to add conditions that constrain those "hanging nodes".
   *
   * The object into which these are inserted is later used to condense the
   * global system matrix and right hand side, and to extend the solution
   * vectors from the true degrees of freedom also to the constraint nodes.
   * This function is explained in detail in the
   * @ref step_6 "step-6"
   * tutorial program and is used in almost all following programs as well.
   *
   * This function does not clear the constraint matrix object before use, in
   * order to allow adding constraints from different sources to the same
   * object. You therefore need to make sure it contains only constraints you
   * still want; otherwise call the ConstraintMatrix::clear() function.
   * Likewise, this function does not close the object since you may want to
   * enter other constraints later on yourself.
   *
   * In the hp-case, i.e. when the argument is of type hp::DoFHandler, we
   * consider constraints due to different finite elements used on two sides
   * of a face between cells as hanging nodes as well. In other words, for hp
   * finite elements, this function computes all constraints due to differing
   * mesh sizes (h) or polynomial degrees (p) between adjacent cells.
   *
   * The template argument (and by consequence the type of the first argument
   * to this function) can be either ::DoFHandler or hp::DoFHandler.
   *
   * @ingroup constraints
   */
  template <typename DoFHandlerType>
  void
  make_hanging_node_constraints (const DoFHandlerType &dof_handler,
                                 ConstraintMatrix     &constraints);

  /**
   * This function is used when different variables in a problem are
   * discretized on different grids, where one grid is strictly coarser than
   * the other. An example are optimization problems where the control
   * variable is often discretized on a coarser mesh than the state variable.
   *
   * The function's result can be stated as follows mathematically: Let ${\cal
   * T}_0$ and ${\cal T}_1$ be two meshes where ${\cal T}_1$ results from
   * ${\cal T}_0$ strictly by refining or leaving alone the cells of ${\cal
   * T}_0$. Using the same finite element on both, there are function spaces
   * ${\cal V}_0$ and ${\cal V}_1$ associated with these meshes. Then every
   * function $v_0 \in {\cal V}_0$ can of course also be represented exactly
   * in ${\cal V}_1$ since by construction ${\cal V}_0 \subset {\cal V}_1$.
   * However, not every function in ${\cal V}_1$ can be expressed as a linear
   * combination of the shape functions of ${\cal V}_0$. The functions that
   * can be represented lie in a homogenous subspace of ${\cal V}_1$ (namely,
   * ${\cal V}_0$, of course) and this subspace can be represented by a linear
   * constraint of the form $CV=0$ where $V$ is the vector of nodal values of
   * functions $v\in {\cal V}_1$. In other words, every function $v_h=\sum_j
   * V_j \varphi_j^{(1)} \in {\cal V}_1$ that also satisfies $v_h\in {\cal
   * V}_0$ automatically satisfies $CV=0$. This function computes the matrix
   * $C$ in the form of a ConstraintMatrix object.
   *
   * The construction of these constraints is done as follows: for each of the
   * degrees of freedom (i.e. shape functions) on the coarse grid, we compute
   * its representation on the fine grid, i.e. how the linear combination of
   * shape functions on the fine grid looks like that resembles the shape
   * function on the coarse grid. From this information, we can then compute
   * the constraints which have to hold if a solution of a linear equation on
   * the fine grid shall be representable on the coarse grid. The exact
   * algorithm how these constraints can be computed is rather complicated and
   * is best understood by reading the source code, which contains many
   * comments.
   *
   * The use of this function is as follows: it accepts as parameters two DoF
   * Handlers, the first of which refers to the coarse grid and the second of
   * which is the fine grid. On both, a finite element is represented by the
   * DoF handler objects, which will usually have several vector components,
   * which may belong to different base elements. The second and fourth
   * parameter of this function therefore state which vector component on the
   * coarse grid shall be used to restrict the stated component on the fine
   * grid. The finite element used for the respective components on the two
   * grids needs to be the same. An example may clarify this: consider an
   * optimization problem with controls $q$ discretized on a coarse mesh and a
   * state variable $u$ (and corresponding Lagrange multiplier $\lambda$)
   * discretized on the fine mesh. These are discretized using piecewise
   * constant discontinuous, continuous linear, and continuous linear
   * elements, respectively. Only the parameter $q$ is represented on the
   * coarse grid, thus the DoFHandler object on the coarse grid represents
   * only one variable, discretized using piecewise constant discontinuous
   * elements. Then, the parameter denoting the vector component on the coarse
   * grid would be zero (the only possible choice, since the variable on the
   * coarse grid is scalar). If the ordering of variables in the fine mesh
   * FESystem is $u, q, \lambda$, then the fourth argument of the function
   * corresponding to the vector component would be one (corresponding to the
   * variable $q$; zero would be $u$, two would be $\lambda$).
   *
   * The function also requires an object of type IntergridMap representing
   * how to get from the coarse mesh cells to the corresponding cells on the
   * fine mesh. This could in principle be generated by the function itself
   * from the two DoFHandler objects, but since it is probably available
   * anyway in programs that use different meshes, the function simply takes
   * it as an argument.
   *
   * The computed constraints are entered into a variable of type
   * ConstraintMatrix; previous contents are not deleted.
   */
  template <int dim, int spacedim>
  void
  compute_intergrid_constraints (const DoFHandler<dim,spacedim>                &coarse_grid,
                                 const unsigned int                             coarse_component,
                                 const DoFHandler<dim,spacedim>                &fine_grid,
                                 const unsigned int                             fine_component,
                                 const InterGridMap<DoFHandler<dim,spacedim> > &coarse_to_fine_grid_map,
                                 ConstraintMatrix                              &constraints);


  /**
   * This function generates a matrix such that when a vector of data with as
   * many elements as there are degrees of freedom of this component on the
   * coarse grid is multiplied to this matrix, we obtain a vector with as many
   * elements as there are global degrees of freedom on the fine grid. All the
   * elements of the other vector components of the finite element fields on
   * the fine grid are not touched.
   *
   * Triangulation of the fine grid can be distributed. When called in
   * parallel, each process has to have a copy of the coarse grid. In this
   * case, function returns transfer representation for a set of locally owned
   * cells.
   *
   * The output of this function is a compressed format that can be used to
   * construct corresponding sparse transfer matrix.
   */
  template <int dim, int spacedim>
  void
  compute_intergrid_transfer_representation (const DoFHandler<dim,spacedim>                         &coarse_grid,
                                             const unsigned int                                      coarse_component,
                                             const DoFHandler<dim,spacedim>                         &fine_grid,
                                             const unsigned int                                      fine_component,
                                             const InterGridMap<DoFHandler<dim,spacedim> >          &coarse_to_fine_grid_map,
                                             std::vector<std::map<types::global_dof_index, float> > &transfer_representation);

  /**
   * @}
   */


  /**
   * @name Periodic boundary conditions
   * @{
   */

  /**
   * Insert the (algebraic) constraints due to periodic boundary conditions
   * into a ConstraintMatrix @p constraint_matrix.
   *
   * Given a pair of not necessarily active boundary faces @p face_1 and @p
   * face_2, this functions constrains all DoFs associated with the boundary
   * described by @p face_1 to the respective DoFs of the boundary described
   * by @p face_2. More precisely:
   *
   * If @p face_1 and @p face_2 are both active faces it adds the DoFs of @p
   * face_1 to the list of constrained DoFs in @p constraint_matrix and adds
   * entries to constrain them to the corresponding values of the DoFs on @p
   * face_2. This happens on a purely algebraic level, meaning, the global DoF
   * with (local face) index <tt>i</tt> on @p face_1 gets constraint to the
   * DoF with (local face) index <tt>i</tt> on @p face_2 (possibly corrected
   * for orientation, see below).
   *
   * Otherwise, if @p face_1 and @p face_2 are not active faces, this function
   * loops recursively over the children of @p face_1 and @p face_2. If only
   * one of the two faces is active, then we recursively iterate over the
   * children of the non-active ones and make sure that the solution function
   * on the refined side equals that on the non-refined face in much the same
   * way as we enforce hanging node constraints at places where differently
   * refined cells come together. (However, unlike hanging nodes, we do not
   * enforce the requirement that there be only a difference of one refinement
   * level between the two sides of the domain you would like to be periodic).
   *
   * This routine only constrains DoFs that are not already constrained. If
   * this routine encounters a DoF that already is constrained (for instance
   * by Dirichlet boundary conditions), the old setting of the constraint
   * (dofs the entry is constrained to, inhomogeneities) is kept and nothing
   * happens.
   *
   * The flags in the @p component_mask (see
   * @ref GlossComponentMask)
   * denote which components of the finite element space shall be constrained
   * with periodic boundary conditions. If it is left as specified by the
   * default value all components are constrained. If it is different from the
   * default value, it is assumed that the number of entries equals the number
   * of components of the finite element. This can be used to enforce
   * periodicity in only one variable in a system of equations.
   *
   * @p face_orientation, @p face_flip and @p face_rotation describe an
   * orientation that should be applied to @p face_1 prior to matching and
   * constraining DoFs. This has nothing to do with the actual orientation of
   * the given faces in their respective cells (which for boundary faces is
   * always the default) but instead how you want to see periodicity to be
   * enforced. For example, by using these flags, you can enforce a condition
   * of the kind $u(0,y)=u(1,1-y)$ (i.e., a Moebius band) or in 3d a twisted
   * torus. More precisely, these flags match local face DoF indices in the
   * following manner:
   *
   * In 2d: <tt>face_orientation</tt> must always be <tt>true</tt>,
   * <tt>face_rotation</tt> is always <tt>false</tt>, and face_flip has the
   * meaning of <tt>line_flip</tt>; this implies e.g. for <tt>Q1</tt>:
   *
   * @code
   *
   * face_orientation = true, face_flip = false, face_rotation = false:
   *
   *     face1:           face2:
   *
   *     1                1
   *     |        <-->    |
   *     0                0
   *
   *     Resulting constraints: 0 <-> 0, 1 <-> 1
   *
   *     (Numbers denote local face DoF indices.)
   *
   *
   * face_orientation = true, face_flip = true, face_rotation = false:
   *
   *     face1:           face2:
   *
   *     0                1
   *     |        <-->    |
   *     1                0
   *
   *     Resulting constraints: 1 <-> 0, 0 <-> 1
   * @endcode
   *
   * And similarly for the case of Q1 in 3d:
   *
   * @code
   *
   * face_orientation = true, face_flip = false, face_rotation = false:
   *
   *     face1:           face2:
   *
   *     2 - 3            2 - 3
   *     |   |    <-->    |   |
   *     0 - 1            0 - 1
   *
   *     Resulting constraints: 0 <-> 0, 1 <-> 1, 2 <-> 2, 3 <-> 3
   *
   *     (Numbers denote local face DoF indices.)
   *
   *
   * face_orientation = false, face_flip = false, face_rotation = false:
   *
   *     face1:           face2:
   *
   *     1 - 3            2 - 3
   *     |   |    <-->    |   |
   *     0 - 2            0 - 1
   *
   *     Resulting constraints: 0 <-> 0, 2 <-> 1, 1 <-> 2, 3 <-> 3
   *
   *
   * face_orientation = true, face_flip = true, face_rotation = false:
   *
   *     face1:           face2:
   *
   *     1 - 0            2 - 3
   *     |   |    <-->    |   |
   *     3 - 2            0 - 1
   *
   *     Resulting constraints: 3 <-> 0, 2 <-> 1, 1 <-> 2, 0 <-> 3
   *
   *
   * face_orientation = true, face_flip = false, face_rotation = true
   *
   *     face1:           face2:
   *
   *     0 - 2            2 - 3
   *     |   |    <-->    |   |
   *     1 - 3            0 - 1
   *
   *     Resulting constraints: 1 <-> 0, 3 <-> 1, 0 <-> 2, 2 <-> 3
   *
   * and any combination of that...
   * @endcode
   *
   * Optionally a matrix @p matrix along with an std::vector @p
   * first_vector_components can be specified that describes how DoFs on @p
   * face_1 should be modified prior to constraining to the DoFs of @p face_2.
   * Here, two declarations are possible: If the std::vector @p
   * first_vector_components is non empty the matrix is interpreted as a @p
   * dim $\times$ @p dim rotation matrix that is applied to all vector valued
   * blocks listed in @p first_vector_components of the FESystem. If @p
   * first_vector_components is empty the matrix is interpreted as an
   * interpolation matrix with size no_face_dofs $\times$ no_face_dofs.
   *
   * Detailed information can be found in the see
   * @ref GlossPeriodicConstraints "Glossary entry on periodic boundary conditions".
   *
   * @author Matthias Maier, 2012 - 2015
   */
  template<typename FaceIterator>
  void
  make_periodicity_constraints
  (const FaceIterator                          &face_1,
   const typename identity<FaceIterator>::type &face_2,
   dealii::ConstraintMatrix                    &constraint_matrix,
   const ComponentMask                         &component_mask = ComponentMask(),
   const bool                                  face_orientation = true,
   const bool                                  face_flip = false,
   const bool                                  face_rotation = false,
   const FullMatrix<double>                    &matrix = FullMatrix<double>(),
   const std::vector<unsigned int>             &first_vector_components = std::vector<unsigned int>());



  /**
   * Insert the (algebraic) constraints due to periodic boundary conditions
   * into a ConstraintMatrix @p constraint_matrix.
   *
   * This is the main high level interface for above low level variant of
   * make_periodicity_constraints(). It takes a std::vector @p periodic_faces
   * as argument and applies above make_periodicity_constraints() on each
   * entry. @p periodic_faces can be created by
   * GridTools::collect_periodic_faces.
   *
   * @note For DoFHandler objects that are built on a
   * parallel::distributed::Triangulation object
   * parallel::distributed::Triangulation::add_periodicity has to be called
   * before calling this function..
   *
   * @see
   * @ref GlossPeriodicConstraints "Glossary entry on periodic boundary conditions"
   * and step-45 for further information.
   *
   * @author Daniel Arndt, Matthias Maier, 2013 - 2015
   */
  template<typename DoFHandlerType>
  void
  make_periodicity_constraints
  (const std::vector<GridTools::PeriodicFacePair<typename DoFHandlerType::cell_iterator> >
   &periodic_faces,
   dealii::ConstraintMatrix        &constraint_matrix,
   const ComponentMask             &component_mask = ComponentMask(),
   const std::vector<unsigned int> &first_vector_components = std::vector<unsigned int>());



  /**
   * Insert the (algebraic) constraints due to periodic boundary conditions
   * into a ConstraintMatrix @p constraint_matrix.
   *
   * This function serves as a high level interface for the
   * make_periodicity_constraints() function.
   *
   * Define a 'first' boundary as all boundary faces having boundary_id @p
   * b_id1 and a 'second' boundary consisting of all faces belonging to @p
   * b_id2.
   *
   * This function tries to match all faces belonging to the first boundary
   * with faces belonging to the second boundary with the help of
   * orthogonal_equality().
   *
   * If this matching is successful it constrains all DoFs associated with the
   * 'first' boundary to the respective DoFs of the 'second' boundary
   * respecting the relative orientation of the two faces.
   *
   * @note: This function is a convenience wrapper. It internally calls
   * GridTools::collect_periodic_faces() with the supplied paramaters and
   * feeds the output to above make_periodicity_constraints() variant. If you
   * need more functionality use GridTools::collect_periodic_faces() directly.
   *
   * @see
   * @ref GlossPeriodicConstraints "Glossary entry on periodic boundary conditions"
   * for further information.
   *
   * @author Matthias Maier, 2012
   */
  template<typename DoFHandlerType>
  void
  make_periodicity_constraints
  (const DoFHandlerType     &dof_handler,
   const types::boundary_id  b_id1,
   const types::boundary_id  b_id2,
   const int                 direction,
   dealii::ConstraintMatrix &constraint_matrix,
   const ComponentMask      &component_mask = ComponentMask());



  /**
   * This compatibility version of make_periodicity_constraints only works on
   * grids with cells in
   * @ref GlossFaceOrientation "standard orientation".
   *
   * Instead of defining a 'first' and 'second' boundary with the help of two
   * boundary_ids this function defines a 'left' boundary as all faces with
   * local face index <code>2*dimension</code> and boundary indicator @p b_id
   * and, similarly, a 'right' boundary consisting of all face with local face
   * index <code>2*dimension+1</code> and boundary indicator @p b_id.
   *
   * @note This version of make_periodicity_constraints  will not work on
   * meshes with cells not in
   * @ref GlossFaceOrientation "standard orientation".
   *
   * @note: This function is a convenience wrapper. It internally calls
   * GridTools::collect_periodic_faces() with the supplied paramaters and
   * feeds the output to above make_periodicity_constraints() variant. If you
   * need more functionality use GridTools::collect_periodic_faces() directly.
   *
   * @see
   * @ref GlossPeriodicConstraints "Glossary entry on periodic boundary conditions"
   * for further information.
   */
  template<typename DoFHandlerType>
  void
  make_periodicity_constraints
  (const DoFHandlerType     &dof_handler,
   const types::boundary_id  b_id,
   const int                 direction,
   dealii::ConstraintMatrix &constraint_matrix,
   const ComponentMask      &component_mask = ComponentMask());

  /**
   * Take a vector of values which live on cells (e.g. an error per cell) and
   * distribute it to the dofs in such a way that a finite element field
   * results, which can then be further processed, e.g. for output. You should
   * note that the resulting field will not be continuous at hanging nodes.
   * This can, however, easily be arranged by calling the appropriate @p
   * distribute function of a ConstraintMatrix object created for this
   * DoFHandler object, after the vector has been fully assembled.
   *
   * It is assumed that the number of elements in @p cell_data equals the
   * number of active cells and that the number of elements in @p dof_data
   * equals <tt>dof_handler.n_dofs()</tt>.
   *
   * Note that the input vector may be a vector of any data type as long as it
   * is convertible to @p double.  The output vector, being a data vector on a
   * DoF handler, always consists of elements of type @p double.
   *
   * In case the finite element used by this DoFHandler consists of more than
   * one component, you need to specify which component in the output vector
   * should be used to store the finite element field in; the default is zero
   * (no other value is allowed if the finite element consists only of one
   * component). All other components of the vector remain untouched, i.e.
   * their contents are not changed.
   *
   * This function cannot be used if the finite element in use has shape
   * functions that are non-zero in more than one vector component (in deal.II
   * speak: they are non-primitive).
   */
  template <typename DoFHandlerType, typename Number>
  void
  distribute_cell_to_dof_vector (const DoFHandlerType  &dof_handler,
                                 const Vector<Number>  &cell_data,
                                 Vector<double>        &dof_data,
                                 const unsigned int     component = 0);

  /**
   * @}
   */

  /**
   * @name Identifying subsets of degrees of freedom with particular properties
   * @{
   */

  /**
   * Extract the indices of the degrees of freedom belonging to certain vector
   * components of a vector-valued finite element. The @p component_mask
   * defines which components or blocks of an FESystem are to be extracted
   * from the DoFHandler @p dof. The entries in the output array @p
   * selected_dofs corresponding to degrees of freedom belonging to these
   * components are then flagged @p true, while all others are set to @p
   * false.
   *
   * The size of @p component_mask must be compatible with the number of
   * components in the FiniteElement used by @p dof. The size of @p
   * selected_dofs must equal DoFHandler::n_dofs(). Previous contents of this
   * array are overwritten.
   *
   * If the finite element under consideration is not primitive, i.e., some or
   * all of its shape functions are non-zero in more than one vector component
   * (which holds, for example, for FE_Nedelec or FE_RaviartThomas elements),
   * then shape functions cannot be associated with a single vector component.
   * In this case, if <em>one</em> shape vector component of this element is
   * flagged in @p component_mask (see
   * @ref GlossComponentMask),
   * then this is equivalent to selecting <em>all</em> vector components
   * corresponding to this non-primitive base element.
   */
  template <int dim, int spacedim>
  void
  extract_dofs (const DoFHandler<dim,spacedim> &dof_handler,
                const ComponentMask            &component_mask,
                std::vector<bool>              &selected_dofs);

  /**
   * The same function as above, but for a hp::DoFHandler.
   */
  template <int dim, int spacedim>
  void
  extract_dofs (const hp::DoFHandler<dim,spacedim> &dof_handler,
                const ComponentMask                &component_mask,
                std::vector<bool>                  &selected_dofs);

  /**
   * This function is the equivalent to the DoFTools::extract_dofs() functions
   * above except that the selection of which degrees of freedom to extract is
   * not done based on components (see
   * @ref GlossComponent)
   * but instead based on whether they are part of a particular block (see
   * @ref GlossBlock).
   * Consequently, the second argument is not a ComponentMask but a BlockMask
   * object.
   *
   * @param dof_handler The DoFHandler object from which to extract degrees of
   * freedom
   * @param block_mask The block mask that describes which blocks to consider
   * (see
   * @ref GlossBlockMask)
   * @param selected_dofs A vector of length DoFHandler::n_dofs() in which
   * those entries are true that correspond to the selected blocks.
   */
  template <int dim, int spacedim>
  void
  extract_dofs (const DoFHandler<dim,spacedim> &dof_handler,
                const BlockMask                &block_mask,
                std::vector<bool>              &selected_dofs);

  /**
   * The same function as above, but for a hp::DoFHandler.
   */
  template <int dim, int spacedim>
  void
  extract_dofs (const hp::DoFHandler<dim,spacedim> &dof_handler,
                const BlockMask                    &block_mask,
                std::vector<bool>                  &selected_dofs);

  /**
   * Do the same thing as the corresponding extract_dofs() function for one
   * level of a multi-grid DoF numbering.
   */
  template <typename DoFHandlerType>
  void
  extract_level_dofs (const unsigned int    level,
                      const DoFHandlerType &dof,
                      const ComponentMask  &component_mask,
                      std::vector<bool>    &selected_dofs);

  /**
   * Do the same thing as the corresponding extract_dofs() function for one
   * level of a multi-grid DoF numbering.
   */
  template <typename DoFHandlerType>
  void
  extract_level_dofs (const unsigned int    level,
                      const DoFHandlerType &dof,
                      const BlockMask      &component_mask,
                      std::vector<bool>    &selected_dofs);

  /**
   * Extract all degrees of freedom which are at the boundary and belong to
   * specified components of the solution. The function returns its results in
   * the last non-default-valued parameter which contains @p true if a degree
   * of freedom is at the boundary and belongs to one of the selected
   * components, and @p false otherwise. The function is used in step-15.
   *
   * By specifying the @p boundary_id variable, you can select which boundary
   * indicators the faces have to have on which the degrees of freedom are
   * located that shall be extracted. If it is an empty list, then all
   * boundary indicators are accepted.
   *
   * The size of @p component_mask (see
   * @ref GlossComponentMask)
   * shall equal the number of components in the finite element used by @p
   * dof. The size of @p selected_dofs shall equal
   * <tt>dof_handler.n_dofs()</tt>. Previous contents of this array or
   * overwritten.
   *
   * Using the usual convention, if a shape function is non-zero in more than
   * one component (i.e. it is non-primitive), then the element in the
   * component mask is used that corresponds to the first non-zero components.
   * Elements in the mask corresponding to later components are ignored.
   *
   * @note This function will not work for DoFHandler objects that are built
   * on a parallel::distributed::Triangulation object. The reasons is that the
   * output argument @p selected_dofs has to have a length equal to <i>all</i>
   * global degrees of freedom. Consequently, this does not scale to very
   * large problems. If you need the functionality of this function for
   * parallel triangulations, then you need to use the other
   * DoFTools::extract_boundary_dofs function.
   *
   * @param[in] dof_handler The object that describes which degrees of freedom
   * live on which cell
   * @param[in] component_mask A mask denoting the vector components of the finite
   * element that should be considered (see also
   * @ref GlossComponentMask).
   * @param[out] selected_dofs A vector of booleans that is returned and for which
   * an element will be @p true if the corresponding index is a
   * degree of freedom that is located on the
   * boundary (and correspond to the selected vector components and boundary
   * indicators, depending on the values of the @p component_mask and @p
   * boundary_ids arguments).
   * @param[in] boundary_ids If empty, this function extracts the indices of the
   * degrees of freedom for all parts of the boundary. If it is a non- empty
   * list, then the function only considers boundary faces with the boundary
   * indicators listed in this argument.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <typename DoFHandlerType>
  void
  extract_boundary_dofs (const DoFHandlerType       &dof_handler,
                         const ComponentMask        &component_mask,
                         std::vector<bool>          &selected_dofs,
                         const std::set<types::boundary_id> &boundary_ids = std::set<types::boundary_id>());

  /**
   * This function does the same as the previous one but it returns its result
   * as an IndexSet rather than a std::vector@<bool@>. Thus, it can also be
   * called for DoFHandler objects that are defined on
   * parallel::distributed::Triangulation objects.
   *
   * @note If the DoFHandler object is indeed defined on a
   * parallel::distributed::Triangulation, then the @p selected_dofs index set
   * will contain only those degrees of freedom on the boundary that belong to
   * the locally relevant set (see
   * @ref GlossLocallyRelevantDof "locally relevant DoFs").
   *
   * @param[in] dof_handler The object that describes which degrees of freedom
   * live on which cell
   * @param[in] component_mask A mask denoting the vector components of the finite
   * element that should be considered (see also
   * @ref GlossComponentMask).
   * @param[out] selected_dofs The IndexSet object that is returned and that will
   * contain the indices of degrees of freedom that are located on the
   * boundary (and correspond to the selected vector components and boundary
   * indicators, depending on the values of the @p component_mask and @p
   * boundary_ids arguments).
   * @param[in] boundary_ids If empty, this function extracts the indices of the
   * degrees of freedom for all parts of the boundary. If it is a non- empty
   * list, then the function only considers boundary faces with the boundary
   * indicators listed in this argument.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <typename DoFHandlerType>
  void
  extract_boundary_dofs (const DoFHandlerType       &dof_handler,
                         const ComponentMask        &component_mask,
                         IndexSet                   &selected_dofs,
                         const std::set<types::boundary_id> &boundary_ids = std::set<types::boundary_id>());

  /**
   * This function is similar to the extract_boundary_dofs() function but it
   * extracts those degrees of freedom whose shape functions are nonzero on at
   * least part of the selected boundary. For continuous elements, this is
   * exactly the set of shape functions whose degrees of freedom are defined
   * on boundary faces. On the other hand, if the finite element in used is a
   * discontinuous element, all degrees of freedom are defined in the inside
   * of cells and consequently none would be boundary degrees of freedom.
   * Several of those would have shape functions that are nonzero on the
   * boundary, however. This function therefore extracts all those for which
   * the FiniteElement::has_support_on_face function says that it is nonzero
   * on any face on one of the selected boundary parts.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <typename DoFHandlerType>
  void
  extract_dofs_with_support_on_boundary (const DoFHandlerType   &dof_handler,
                                         const ComponentMask    &component_mask,
                                         std::vector<bool>      &selected_dofs,
                                         const std::set<types::boundary_id> &boundary_ids = std::set<types::boundary_id>());

  /**
   * Extract a vector that represents the constant modes of the DoFHandler for
   * the components chosen by <tt>component_mask</tt> (see
   * @ref GlossComponentMask).
   * The constant modes on a discretization are the null space of a Laplace
   * operator on the selected components with Neumann boundary conditions
   * applied. The null space is a necessary ingredient for obtaining a good
   * AMG preconditioner when using the class
   * TrilinosWrappers::PreconditionAMG.  Since the ML AMG package only works
   * on algebraic properties of the respective matrix, it has no chance to
   * detect whether the matrix comes from a scalar or a vector valued problem.
   * However, a near null space supplies exactly the needed information about
   * the components placement of vector components within the matrix. The null
   * space (or rather, the constant modes) is provided by the finite element
   * underlying the given DoFHandler and for most elements, the null space
   * will consist of as many vectors as there are true arguments in
   * <tt>component_mask</tt> (see
   * @ref GlossComponentMask),
   * each of which will be one in one vector component and zero in all others.
   * However, the representation of the constant function for e.g. FE_DGP is
   * different (the first component on each element one, all other components
   * zero), and some scalar elements may even have two constant modes
   * (FE_Q_DG0). Therefore, we store this object in a vector of vectors, where
   * the outer vector contains the collection of the actual constant modes on
   * the DoFHandler. Each inner vector has as many components as there are
   * (locally owned) degrees of freedom in the selected components. Note that
   * any matrix associated with this null space must have been constructed
   * using the same <tt>component_mask</tt> argument, since the numbering of
   * DoFs is done relative to the selected dofs, not to all dofs.
   *
   * The main reason for this program is the use of the null space with the
   * AMG preconditioner.
   */
  template <typename DoFHandlerType>
  void
  extract_constant_modes (const DoFHandlerType            &dof_handler,
                          const ComponentMask             &component_mask,
                          std::vector<std::vector<bool> > &constant_modes);

  /**
   * @}
   */
  /**
   * @name Hanging nodes
   * @{
   */

  /**
   * Select all dofs that will be constrained by interface constraints, i.e.
   * all hanging nodes.
   *
   * The size of @p selected_dofs shall equal <tt>dof_handler.n_dofs()</tt>.
   * Previous contents of this array or overwritten.
   */
  template <int dim, int spacedim>
  void
  extract_hanging_node_dofs (const DoFHandler<dim,spacedim> &dof_handler,
                             std::vector<bool>              &selected_dofs);
  //@}

  /**
   * @name Parallelization and domain decomposition
   * @{
   */
  /**
   * Flag all those degrees of freedom which are on cells with the given
   * subdomain id. Note that DoFs on faces can belong to cells with differing
   * subdomain ids, so the sets of flagged degrees of freedom are not mutually
   * exclusive for different subdomain ids.
   *
   * If you want to get a unique association of degree of freedom with
   * subdomains, use the @p get_subdomain_association function.
   */
  template <typename DoFHandlerType>
  void
  extract_subdomain_dofs (const DoFHandlerType      &dof_handler,
                          const types::subdomain_id  subdomain_id,
                          std::vector<bool>         &selected_dofs);


  /**
   * Extract the set of global DoF indices that are owned by the current
   * processor. For regular DoFHandler objects, this set is the complete set
   * with all DoF indices. In either case, it equals what
   * DoFHandler::locally_owned_dofs() returns.
   */
  template <typename DoFHandlerType>
  void
  extract_locally_owned_dofs (const DoFHandlerType &dof_handler,
                              IndexSet             &dof_set);


  /**
   * Extract the set of global DoF indices that are active on the current
   * DoFHandler. For regular DoFHandlers, these are all DoF indices, but for
   * DoFHandler objects built on parallel::distributed::Triangulation this set
   * is a superset of DoFHandler::locally_owned_dofs() and contains all DoF
   * indices that live on all locally owned cells (including on the interface
   * to ghost cells). However, it does not contain the DoF indices that are
   * exclusively defined on ghost or artificial cells (see
   * @ref GlossArtificialCell "the glossary").
   *
   * The degrees of freedom identified by this function equal those obtained
   * from the dof_indices_with_subdomain_association() function when called
   * with the locally owned subdomain id.
   */
  template <typename DoFHandlerType>
  void
  extract_locally_active_dofs (const DoFHandlerType &dof_handler,
                               IndexSet             &dof_set);

  /**
   * Extract the set of global DoF indices that are active on the current
   * DoFHandler. For regular DoFHandlers, these are all DoF indices, but for
   * DoFHandler objects built on parallel::distributed::Triangulation this set
   * is the union of DoFHandler::locally_owned_dofs() and the DoF indices on
   * all ghost cells. In essence, it is the DoF indices on all cells that are
   * not artificial (see
   * @ref GlossArtificialCell "the glossary").
   */
  template <typename DoFHandlerType>
  void
  extract_locally_relevant_dofs (const DoFHandlerType &dof_handler,
                                 IndexSet             &dof_set);

  /**
   *
   * For each processor, determine the set of locally owned degrees of freedom
   * as an IndexSet. This function then returns a vector of index sets, where
   * the vector has size equal to the number of MPI processes that participate
   * in the DoF handler object.
   *
   * The function can be used for objects of type dealii::Triangulation or
   * parallel::shared::Triangulation. It will not work for objects of type
   * parallel::distributed::Triangulation since for such triangulations we do
   * not have information about all cells of the triangulation available
   * locally, and consequently can not say anything definitive about the
   * degrees of freedom active on other processors' locally owned cells.
   *
   * @author Denis Davydov, 2015
   */
  template <typename DoFHandlerType>
  std::vector<IndexSet>
  locally_owned_dofs_per_subdomain (const DoFHandlerType &dof_handler);

  /**
   *
   * For each processor, determine the set of locally relevant degrees of
   * freedom as an IndexSet. This function then returns a vector of index
   * sets, where the vector has size equal to the number of MPI processes that
   * participate in the DoF handler object.
   *
   * The function can be used for objects of type dealii::Triangulation or
   * parallel::shared::Triangulation. It will not work for objects of type
   * parallel::distributed::Triangulation since for such triangulations we do
   * not have information about all cells of the triangulation available
   * locally, and consequently can not say anything definitive about the
   * degrees of freedom active on other processors' locally owned cells.
   *
   * @author Jean-Paul Pelteret, 2015
   */
  template <typename DoFHandlerType>
  std::vector<IndexSet>
  locally_relevant_dofs_per_subdomain (const DoFHandlerType &dof_handler);


  /**
   * Same as extract_locally_relevant_dofs() but for multigrid DoFs for the
   * given @p level.
   */
  template <typename DoFHandlerType>
  void
  extract_locally_relevant_level_dofs (const DoFHandlerType &dof_handler,
                                       const unsigned int    level,
                                       IndexSet             &dof_set);


  /**
   * For each degree of freedom, return in the output array to which subdomain
   * (as given by the <tt>cell->subdomain_id()</tt> function) it belongs. The
   * output array is supposed to have the right size already when calling this
   * function.
   *
   * Note that degrees of freedom associated with faces, edges, and vertices
   * may be associated with multiple subdomains if they are sitting on
   * partition boundaries. In these cases, we put them into one of the
   * associated partitions in an undefined way. This may sometimes lead to
   * different numbers of degrees of freedom in partitions, even if the number
   * of cells is perfectly equidistributed. While this is regrettable, it is
   * not a problem in practice since the number of degrees of freedom on
   * partition boundaries is asymptotically vanishing as we refine the mesh as
   * long as the number of partitions is kept constant.
   *
   * This function returns the association of each DoF with one subdomain. If
   * you are looking for the association of each @em cell with a subdomain,
   * either query the <tt>cell->subdomain_id()</tt> function, or use the
   * <tt>GridTools::get_subdomain_association</tt> function.
   *
   * Note that this function is of questionable use for DoFHandler objects
   * built on parallel::distributed::Triangulation since in that case
   * ownership of individual degrees of freedom by MPI processes is controlled
   * by the DoF handler object, not based on some geometric algorithm in
   * conjunction with subdomain id. In particular, the degrees of freedom
   * identified by the functions in this namespace as associated with a
   * subdomain are not the same the DoFHandler class identifies as those it
   * owns.
   */
  template <typename DoFHandlerType>
  void
  get_subdomain_association (const DoFHandlerType             &dof_handler,
                             std::vector<types::subdomain_id> &subdomain);

  /**
   * Count how many degrees of freedom are uniquely associated with the given
   * @p subdomain index.
   *
   * Note that there may be rare cases where cells with the given @p subdomain
   * index exist, but none of its degrees of freedom are actually associated
   * with it. In that case, the returned value will be zero.
   *
   * This function will generate an exception if there are no cells with the
   * given @p subdomain index.
   *
   * This function returns the number of DoFs associated with one subdomain.
   * If you are looking for the association of @em cells with this subdomain,
   * use the <tt>GridTools::count_cells_with_subdomain_association</tt>
   * function.
   *
   * Note that this function is of questionable use for DoFHandler objects
   * built on parallel::distributed::Triangulation since in that case
   * ownership of individual degrees of freedom by MPI processes is controlled
   * by the DoF handler object, not based on some geometric algorithm in
   * conjunction with subdomain id. In particular, the degrees of freedom
   * identified by the functions in this namespace as associated with a
   * subdomain are not the same the DoFHandler class identifies as those it
   * owns.
   */
  template <typename DoFHandlerType>
  unsigned int
  count_dofs_with_subdomain_association (const DoFHandlerType      &dof_handler,
                                         const types::subdomain_id  subdomain);

  /**
   * Count how many degrees of freedom are uniquely associated with the given
   * @p subdomain index.
   *
   * This function does what the previous one does except that it splits the
   * result among the vector components of the finite element in use by the
   * DoFHandler object. The last argument (which must have a length equal to
   * the number of vector components) will therefore store how many degrees of
   * freedom of each vector component are associated with the given subdomain.
   *
   * Note that this function is of questionable use for DoFHandler objects
   * built on parallel::distributed::Triangulation since in that case
   * ownership of individual degrees of freedom by MPI processes is controlled
   * by the DoF handler object, not based on some geometric algorithm in
   * conjunction with subdomain id. In particular, the degrees of freedom
   * identified by the functions in this namespace as associated with a
   * subdomain are not the same the DoFHandler class identifies as those it
   * owns.
   */
  template <typename DoFHandlerType>
  void
  count_dofs_with_subdomain_association (const DoFHandlerType      &dof_handler,
                                         const types::subdomain_id  subdomain,
                                         std::vector<unsigned int> &n_dofs_on_subdomain);

  /**
   * Return a set of indices that denotes the degrees of freedom that live on
   * the given subdomain, i.e. that are on cells owned by the current
   * processor. Note that this includes the ones that this subdomain "owns"
   * (i.e. the ones for which get_subdomain_association() returns a value
   * equal to the subdomain given here and that are selected by the
   * extract_locally_owned_dofs() function) but also all of those that sit on
   * the boundary between the given subdomain and other subdomain. In essence,
   * degrees of freedom that sit on boundaries between subdomain will be in
   * the index sets returned by this function for more than one subdomain.
   *
   * Note that this function is of questionable use for DoFHandler objects
   * built on parallel::distributed::Triangulation since in that case
   * ownership of individual degrees of freedom by MPI processes is controlled
   * by the DoF handler object, not based on some geometric algorithm in
   * conjunction with subdomain id. In particular, the degrees of freedom
   * identified by the functions in this namespace as associated with a
   * subdomain are not the same the DoFHandler class identifies as those it
   * owns.
   */
  template <typename DoFHandlerType>
  IndexSet
  dof_indices_with_subdomain_association (const DoFHandlerType &dof_handler,
                                          const types::subdomain_id subdomain);
  // @}
  /**
   * @name DoF indices on patches of cells
   *
   * Create structures containing a large set of degrees of freedom for small
   * patches of cells. The resulting objects can be used in RelaxationBlockSOR
   * and related classes to implement Schwarz preconditioners and smoothers,
   * where the subdomains consist of small numbers of cells only.
   */
  //@{
  /**
   * Create an incidence matrix that for every cell on a given level of a
   * multilevel DoFHandler flags which degrees of freedom are associated with
   * the corresponding cell. This data structure is a matrix with as many rows
   * as there are cells on a given level, as many columns as there are degrees
   * of freedom on this level, and entries that are either true or false. This
   * data structure is conveniently represented by a SparsityPattern object.
   *
   * @note The ordering of rows (cells) follows the ordering of the standard
   * cell iterators.
   */
  template <typename DoFHandlerType, class SparsityPatternType>
  void make_cell_patches(SparsityPatternType     &block_list,
                         const DoFHandlerType    &dof_handler,
                         const unsigned int       level,
                         const std::vector<bool> &selected_dofs = std::vector<bool>(),
                         types::global_dof_index  offset        = 0);

  /**
   * Create an incidence matrix that for every vertex on a given level of a
   * multilevel DoFHandler flags which degrees of freedom are associated with
   * the adjacent cells. This data structure is a matrix with as many rows as
   * there are vertices on a given level, as many columns as there are degrees
   * of freedom on this level, and entries that are either true or false. This
   * data structure is conveniently represented by a SparsityPattern object.
   * The sparsity pattern may be empty when entering this function and will be
   * reinitialized to the correct size.
   *
   * The function has some boolean arguments (listed below) controlling
   * details of the generated patches. The default settings are those for
   * Arnold-Falk-Winther type smoothers for divergence and curl conforming
   * finite elements with essential boundary conditions. Other applications
   * are possible, in particular changing <tt>boundary_patches</tt> for non-
   * essential boundary conditions.
   *
   * This function returns the <tt>vertex_mapping</tt>,
   * that contains the mapping from the vertex indices to the block indices
   * of the <tt>block_list</tt>. For vertices that do not lead to a vertex patch, the
   * entry in <tt>vertex_mapping</tt> contains the value <tt>invalid_unsigned_int</tt>.
   * If <tt>invert_vertex_mapping</tt> is set to <tt>true</tt>, then the
   * <tt>vertex_mapping</tt> is inverted such that it contains the mapping from
   * the block indices to the corresponding vertex indices.
   *
   * @arg <tt>block_list</tt>: the SparsityPattern into which the patches will
   * be stored.
   *
   * @arg <tt>dof_handler</tt>: the multilevel dof handler providing the
   * topology operated on.
   *
   * @arg <tt>interior_dofs_only</tt>: for each patch of cells around a
   * vertex, collect only the interior degrees of freedom of the patch and
   * disregard those on the boundary of the patch. This is for instance the
   * setting for smoothers of Arnold-Falk-Winther type.
   *
   * @arg <tt>boundary_patches</tt>: include patches around vertices at the
   * boundary of the domain. If not, only patches around interior vertices
   * will be generated.
   *
   * @arg <tt>level_boundary_patches</tt>: same for refinement edges towards
   * coarser cells.
   *
   * @arg <tt>single_cell_patches</tt>: if not true, patches containing a
   * single cell are eliminated.
   *
   * @arg <tt>invert_vertex_mapping</tt>: if true, then the return value
   * contains one vertex index for each block; if false, then the return value
   * contains one block index or <tt>invalid_unsigned_int</tt> for each vertex.
   */
  template <typename DoFHandlerType>
  std::vector<unsigned int>
  make_vertex_patches(SparsityPattern      &block_list,
                      const DoFHandlerType &dof_handler,
                      const unsigned int    level,
                      const bool            interior_dofs_only,
                      const bool            boundary_patches       = false,
                      const bool            level_boundary_patches = false,
                      const bool            single_cell_patches    = false,
                      const bool            invert_vertex_mapping = false);

  /**
   * Create an incidence matrix that for every cell on a given level of a
   * multilevel DoFHandler flags which degrees of freedom are associated with
   * children of this cell. This data structure is conveniently represented by
   * a SparsityPattern object.
   *
   * The function thus creates a sparsity pattern which in each row (with rows
   * corresponding to the cells on this level) lists the degrees of freedom
   * associated to the cells that are the children of this cell. The DoF
   * indices used here are level dof indices of a multilevel hierarchy, i.e.,
   * they may be associated with children that are not themselves active. The
   * sparsity pattern may be empty when entering this function and will be
   * reinitialized to the correct size.
   *
   * The function has some boolean arguments (listed below) controlling
   * details of the generated patches. The default settings are those for
   * Arnold-Falk-Winther type smoothers for divergence and curl conforming
   * finite elements with essential boundary conditions. Other applications
   * are possible, in particular changing <tt>boundary_dofs</tt> for non-
   * essential boundary conditions.
   *
   * @arg <tt>block_list</tt>: the SparsityPattern into which the patches will
   * be stored.
   *
   * @arg <tt>dof_handler</tt>: The multilevel dof handler providing the
   * topology operated on.
   *
   * @arg <tt>interior_dofs_only</tt>: for each patch of cells around a
   * vertex, collect only the interior degrees of freedom of the patch and
   * disregard those on the boundary of the patch. This is for instance the
   * setting for smoothers of Arnold-Falk-Winther type.
   *
   * @arg <tt>boundary_dofs</tt>: include degrees of freedom, which would have
   * excluded by <tt>interior_dofs_only</tt>, but are lying on the boundary of
   * the domain, and thus need smoothing. This parameter has no effect if
   * <tt>interior_dofs_only</tt> is false.
   */
  template <typename DoFHandlerType>
  void make_child_patches(SparsityPattern      &block_list,
                          const DoFHandlerType &dof_handler,
                          const unsigned int    level,
                          const bool            interior_dofs_only,
                          const bool            boundary_dofs = false);

  /**
   * Create a block list with only a single patch, which in turn contains all
   * degrees of freedom on the given level.
   *
   * This function is mostly a closure on level 0 for functions like
   * make_child_patches() and make_vertex_patches(), which may produce an
   * empty patch list.
   *
   * @arg <tt>block_list</tt>: the SparsityPattern into which the patches will
   * be stored.
   *
   * @arg <tt>dof_handler</tt>: The multilevel dof handler providing the
   * topology operated on.
   *
   * @arg <tt>level</tt> The grid level used for building the list.
   *
   * @arg <tt>interior_dofs_only</tt>: if true, exclude degrees of freedom on
   * the boundary of the domain.
   */
  template <typename DoFHandlerType>
  void make_single_patch(SparsityPattern      &block_list,
                         const DoFHandlerType &dof_handler,
                         const unsigned int    level,
                         const bool            interior_dofs_only = false);

  /**
   * @}
   */
  /**
   * @name Counting degrees of freedom and related functions
   * @{
   */

  /**
   * Count how many degrees of freedom out of the total number belong to each
   * component. If the number of components the finite element has is one
   * (i.e. you only have one scalar variable), then the number in this
   * component obviously equals the total number of degrees of freedom.
   * Otherwise, the sum of the DoFs in all the components needs to equal the
   * total number.
   *
   * However, the last statement does not hold true if the finite element is
   * not primitive, i.e. some or all of its shape functions are non-zero in
   * more than one vector component. This applies, for example, to the Nedelec
   * or Raviart-Thomas elements. In this case, a degree of freedom is counted
   * in each component in which it is non-zero, so that the sum mentioned
   * above is greater than the total number of degrees of freedom.
   *
   * This behavior can be switched off by the optional parameter
   * <tt>vector_valued_once</tt>. If this is <tt>true</tt>, the number of
   * components of a nonprimitive vector valued element is collected only in
   * the first component. All other components will have a count of zero.
   *
   * The additional optional argument @p target_component allows for a re-
   * sorting and grouping of components. To this end, it contains for each
   * component the component number it shall be counted as. Having the same
   * number entered several times sums up several components as the same. One
   * of the applications of this argument is when you want to form block
   * matrices and vectors, but want to pack several components into the same
   * block (for example, when you have @p dim velocities and one pressure, to
   * put all velocities into one block, and the pressure into another).
   *
   * The result is returned in @p dofs_per_component. Note that the size of @p
   * dofs_per_component needs to be enough to hold all the indices specified
   * in @p target_component. If this is not the case, an assertion is thrown.
   * The indices not targeted by target_components are left untouched.
   */
  template <typename DoFHandlerType>
  void
  count_dofs_per_component (const DoFHandlerType                 &dof_handler,
                            std::vector<types::global_dof_index> &dofs_per_component,
                            const bool                            vector_valued_once = false,
                            std::vector<unsigned int> target_component
                            = std::vector<unsigned int>());

  /**
   * Count the degrees of freedom in each block. This function is similar to
   * count_dofs_per_component(), with the difference that the counting is done
   * by blocks. See
   * @ref GlossBlock "blocks"
   * in the glossary for details. Again the vectors are assumed to have the
   * correct size before calling this function. If this is not the case, an
   * assertion is thrown.
   *
   * This function is used in the step-22, step-31, and step-32 tutorial
   * programs.
   *
   * @pre The dofs_per_block variable has as many components as the finite
   * element used by the dof_handler argument has blocks, or alternatively as
   * many blocks as are enumerated in the target_blocks argument if given.
   */
  template <typename DoFHandlerType>
  void
  count_dofs_per_block (const DoFHandlerType                 &dof,
                        std::vector<types::global_dof_index> &dofs_per_block,
                        const std::vector<unsigned int>  &target_block
                        = std::vector<unsigned int>());

  /**
   * For each active cell of a DoFHandler or hp::DoFHandler, extract the
   * active finite element index and fill the vector given as second argument.
   * This vector is assumed to have as many entries as there are active cells.
   *
   * For non-hp DoFHandler objects given as first argument, the returned
   * vector will consist of only zeros, indicating that all cells use the same
   * finite element. For a hp::DoFHandler, the values may be different,
   * though.
   */
  template <typename DoFHandlerType>
  void
  get_active_fe_indices (const DoFHandlerType      &dof_handler,
                         std::vector<unsigned int> &active_fe_indices);

  /**
   * Count how many degrees of freedom live on a set of cells (i.e., a patch)
   * described by the argument.
   *
   * Patches are often used in defining error estimators that require the
   * solution of a local problem on the patch surrounding each of the cells of
   * the mesh. You can get a list of cells that form the patch around a given
   * cell using GridTools::get_patch_around_cell(). This function is then
   * useful in setting up the size of the linear system used to solve the
   * local problem on the patch around a cell. The function
   * DoFTools::get_dofs_on_patch() will then help to make the connection
   * between global degrees of freedom and the local ones.
   *
   * @tparam DoFHandlerType A type that is either DoFHandler or
   * hp::DoFHandler. In C++, the compiler can not determine the type of
   * <code>DoFHandlerType</code> from the function call. You need to specify
   * it as an explicit template argument following the function name.
   *
   * @param patch A collection of cells within an object of type
   * DoFHandlerType
   *
   * @return The number of degrees of freedom associated with the cells of
   * this patch.
   *
   * @note In the context of a parallel distributed computation, it only makes
   * sense to call this function on patches around locally owned cells. This
   * is because the neighbors of locally owned cells are either locally owned
   * themselves, or ghost cells. For both, we know that these are in fact the
   * real cells of the complete, parallel triangulation. We can also query the
   * degrees of freedom on these. In other words, this function can only work
   * if all cells in the patch are either locally owned or ghost cells.
   *
   * @author Arezou Ghesmati, Wolfgang Bangerth, 2014
   */
  template <typename DoFHandlerType>
  unsigned int
  count_dofs_on_patch (const std::vector<typename DoFHandlerType::active_cell_iterator> &patch);

  /**
   * Return the set of degrees of freedom that live on a set of cells (i.e., a
   * patch) described by the argument.
   *
   * Patches are often used in defining error estimators that require the
   * solution of a local problem on the patch surrounding each of the cells of
   * the mesh. You can get a list of cells that form the patch around a given
   * cell using GridTools::get_patch_around_cell(). While
   * DoFTools::count_dofs_on_patch() can be used to determine the size of
   * these local problems, so that one can assemble the local system and then
   * solve it, it is still necessary to provide a mapping between the global
   * indices of the degrees of freedom that live on the patch and a local
   * enumeration. This function provides such a local enumeration by returning
   * the set of degrees of freedom that live on the patch.
   *
   * Since this set is returned in the form of a std::vector, one can also
   * think of it as a mapping
   * @code
   *   i -> global_dof_index
   * @endcode
   * where <code>i</code> is an index into the returned vector (i.e., a the
   * <i>local</i> index of a degree of freedom on the patch) and
   * <code>global_dof_index</code> is the global index of a degree of freedom
   * located on the patch. The array returned has size equal to
   * DoFTools::count_dofs_on_patch().
   *
   * @note The array returned is sorted by global DoF index. Consequently, if
   * one considers the index into this array a local DoF index, then the local
   * system that results retains the block structure of the global system.
   *
   * @tparam DoFHandlerType A type that is either DoFHandler or
   * hp::DoFHandler. In C++, the compiler can not determine the type of
   * <code>DoFHandlerType</code> from the function call. You need to specify
   * it as an explicit template argument following the function name.
   *
   * @param patch A collection of cells within an object of type
   * DoFHandlerType
   *
   * @return A list of those global degrees of freedom located on the patch,
   * as defined above.
   *
   * @note In the context of a parallel distributed computation, it only makes
   * sense to call this function on patches around locally owned cells. This
   * is because the neighbors of locally owned cells are either locally owned
   * themselves, or ghost cells. For both, we know that these are in fact the
   * real cells of the complete, parallel triangulation. We can also query the
   * degrees of freedom on these. In other words, this function can only work
   * if all cells in the patch are either locally owned or ghost cells.
   *
   * @author Arezou Ghesmati, Wolfgang Bangerth, 2014
   */
  template <typename DoFHandlerType>
  std::vector<types::global_dof_index>
  get_dofs_on_patch (const std::vector<typename DoFHandlerType::active_cell_iterator> &patch);

  /**
   * @}
   */

  /**
   * Create a mapping from degree of freedom indices to the index of that
   * degree of freedom on the boundary. After this operation,
   * <tt>mapping[dof]</tt> gives the index of the degree of freedom with
   * global number @p dof in the list of degrees of freedom on the boundary.
   * If the degree of freedom requested is not on the boundary, the value of
   * <tt>mapping[dof]</tt> is @p invalid_dof_index. This function is mainly
   * used when setting up matrices and vectors on the boundary from the trial
   * functions, which have global numbers, while the matrices and vectors use
   * numbers of the trial functions local to the boundary.
   *
   * Prior content of @p mapping is deleted.
   */
  template <typename DoFHandlerType>
  void
  map_dof_to_boundary_indices (const DoFHandlerType                 &dof_handler,
                               std::vector<types::global_dof_index> &mapping);

  /**
   * Same as the previous function, except that only those parts of the
   * boundary are considered for which the boundary indicator is listed in the
   * second argument.
   *
   * See the general doc of this class for more information.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <typename DoFHandlerType>
  void
  map_dof_to_boundary_indices (const DoFHandlerType                 &dof_handler,
                               const std::set<types::boundary_id>   &boundary_ids,
                               std::vector<types::global_dof_index> &mapping);

  /**
   * Return a list of support points (see this
   * @ref GlossSupport "glossary entry")
   * for all the degrees of freedom handled by this DoF handler object. This
   * function, of course, only works if the finite element object used by the
   * DoF handler object actually provides support points, i.e. no edge
   * elements or the like. Otherwise, an exception is thrown.
   *
   * @pre The given array must have a length of as many elements as there are
   * degrees of freedom.
   *
   * @note The precondition to this function that the output argument needs to
   * have size equal to the total number of degrees of freedom makes this
   * function unsuitable for the case that the given DoFHandler object derives
   * from a parallel::distributed::Triangulation object.  Consequently, this
   * function will produce an error if called with such a DoFHandler.
   */
  template <int dim, int spacedim>
  void
  map_dofs_to_support_points (const Mapping<dim,spacedim>       &mapping,
                              const DoFHandler<dim,spacedim>    &dof_handler,
                              std::vector<Point<spacedim> >     &support_points);

  /**
   * Same as the previous function but for the hp case.
   */

  template <int dim, int spacedim>
  void
  map_dofs_to_support_points (const dealii::hp::MappingCollection<dim,spacedim>   &mapping,
                              const hp::DoFHandler<dim,spacedim>    &dof_handler,
                              std::vector<Point<spacedim> > &support_points);

  /**
   * This function is a version of the above map_dofs_to_support_points
   * function that doesn't simply return a vector of support points (see this
   * @ref GlossSupport "glossary entry")
   * with one entry for each global degree of freedom, but instead a map that
   * maps from the DoFs index to its location. The point of this function is
   * that it is also usable in cases where the DoFHandler is based on a
   * parallel::distributed::Triangulation object. In such cases, each
   * processor will not be able to determine the support point location of all
   * DoFs, and worse no processor may be able to hold a vector that would
   * contain the locations of all DoFs even if they were known. As a
   * consequence, this function constructs a map from those DoFs for which we
   * can know the locations (namely, those DoFs that are locally relevant (see
   * @ref GlossLocallyRelevantDof "locally relevant DoFs")
   * to their locations.
   *
   * For non-distributed triangulations, the map returned as @p support_points
   * is of course dense, i.e., every DoF is to be found in it.
   *
   * @param mapping The mapping from the reference cell to the real cell on
   * which DoFs are defined.
   * @param dof_handler The object that describes which DoF indices live on
   * which cell of the triangulation.
   * @param support_points A map that for every locally relevant DoF index
   * contains the corresponding location in real space coordinates. Previous
   * content of this object is deleted in this function.
   */
  template <int dim, int spacedim>
  void
  map_dofs_to_support_points (const Mapping<dim,spacedim>       &mapping,
                              const DoFHandler<dim,spacedim>    &dof_handler,
                              std::map<types::global_dof_index, Point<spacedim> >     &support_points);

  /**
   * Same as the previous function but for the hp case.
   */
  template <int dim, int spacedim>
  void
  map_dofs_to_support_points (const dealii::hp::MappingCollection<dim,spacedim>   &mapping,
                              const hp::DoFHandler<dim,spacedim>    &dof_handler,
                              std::map<types::global_dof_index, Point<spacedim> > &support_points);


  /**
   * This is the opposite function to the one above. It generates a map where
   * the keys are the support points of the degrees of freedom, while the
   * values are the DoF indices. For a definition of support points, see this
   * @ref GlossSupport "glossary entry".
   *
   * Since there is no natural order in the space of points (except for the 1d
   * case), you have to provide a map with an explicitly specified comparator
   * object. This function is therefore templatized on the comparator object.
   * Previous content of the map object is deleted in this function.
   *
   * Just as with the function above, it is assumed that the finite element in
   * use here actually supports the notion of support points of all its
   * components.
   */
  template <typename DoFHandlerType, class Comp>
  void
  map_support_points_to_dofs
  (const Mapping<DoFHandlerType::dimension, DoFHandlerType::space_dimension>       &mapping,
   const DoFHandlerType                                                            &dof_handler,
   std::map<Point<DoFHandlerType::space_dimension>, types::global_dof_index, Comp> &point_to_index_map);

  /**
   * Map a coupling table from the user friendly organization by components to
   * the organization by blocks. Specializations of this function for
   * DoFHandler and hp::DoFHandler are required due to the different results
   * of their finite element access.
   *
   * The return vector will be initialized to the correct length inside this
   * function.
   */
  template <int dim, int spacedim>
  void
  convert_couplings_to_blocks (const hp::DoFHandler<dim,spacedim> &dof_handler,
                               const Table<2, Coupling> &table_by_component,
                               std::vector<Table<2,Coupling> > &tables_by_block);

  /**
   * Make a constraint matrix for the constraints that result from zero
   * boundary values on the given boundary indicator.
   *
   * This function constrains all degrees of freedom on the given part of the
   * boundary.
   *
   * A variant of this function with different arguments is used in step-36.
   *
   * @param dof The DoFHandler to work on.
   * @param boundary_id The indicator of that part of the boundary for which
   * constraints should be computed. If this number equals
   * numbers::invalid_boundary_id then all boundaries of the domain will be
   * treated.
   * @param zero_boundary_constraints The constraint object into which the
   * constraints will be written. The new constraints due to zero boundary
   * values will simply be added, preserving any other constraints previously
   * present. However, this will only work if the previous content of that
   * object consists of constraints on degrees of freedom that are not located
   * on the boundary treated here. If there are previously existing
   * constraints for degrees of freedom located on the boundary, then this
   * would constitute a conflict. See the
   * @ref constraints
   * module for handling the case where there are conflicting constraints on
   * individual degrees of freedom.
   * @param component_mask An optional component mask that restricts the
   * functionality of this function to a subset of an FESystem. For non-
   * @ref GlossPrimitive "primitive"
   * shape functions, any degree of freedom is affected that belongs to a
   * shape function where at least one of its nonzero components is affected
   * by the component mask (see
   * @ref GlossComponentMask).
   * If this argument is omitted, all components of the finite element with
   * degrees of freedom at the boundary will be considered.
   *
   * @ingroup constraints
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim, int spacedim, template <int, int> class DoFHandlerType>
  void
  make_zero_boundary_constraints (const DoFHandlerType<dim,spacedim> &dof,
                                  const types::boundary_id            boundary_id,
                                  ConstraintMatrix                   &zero_boundary_constraints,
                                  const ComponentMask                &component_mask = ComponentMask());

  /**
   * Do the same as the previous function, except do it for all parts of the
   * boundary, not just those with a particular boundary indicator. This
   * function is then equivalent to calling the previous one with
   * numbers::invalid_boundary_id as second argument.
   *
   * This function is used in step-36, for example.
   *
   * @ingroup constraints
   */
  template <int dim, int spacedim, template <int, int> class DoFHandlerType>
  void
  make_zero_boundary_constraints (const DoFHandlerType<dim,spacedim> &dof,
                                  ConstraintMatrix                   &zero_boundary_constraints,
                                  const ComponentMask                &component_mask = ComponentMask());


  /**
   * Map a coupling table from the user friendly organization by components to
   * the organization by blocks. Specializations of this function for
   * DoFHandler and hp::DoFHandler are required due to the different results
   * of their finite element access.
   *
   * The return vector will be initialized to the correct length inside this
   * function.
   */
  template <int dim, int spacedim>
  void
  convert_couplings_to_blocks (const DoFHandler<dim,spacedim> &dof_handler,
                               const Table<2, Coupling> &table_by_component,
                               std::vector<Table<2,Coupling> > &tables_by_block);

  /**
   * Given a finite element and a table how the vector components of it couple
   * with each other, compute and return a table that describes how the
   * individual shape functions couple with each other.
   */
  template <int dim, int spacedim>
  Table<2,Coupling>
  dof_couplings_from_component_couplings (const FiniteElement<dim,spacedim> &fe,
                                          const Table<2,Coupling> &component_couplings);

  /**
   * Same function as above for a collection of finite elements, returning a
   * collection of tables.
   *
   * The function currently treats DoFTools::Couplings::nonzero the same as
   * DoFTools::Couplings::always .
   */
  template <int dim, int spacedim>
  std::vector<Table<2,Coupling> >
  dof_couplings_from_component_couplings (const hp::FECollection<dim,spacedim> &fe,
                                          const Table<2,Coupling> &component_couplings);
  /**
   * @todo Write description
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcFiniteElementsDontMatch);
  /**
   * @todo Write description
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcGridNotCoarser);
  /**
   * @todo Write description
   *
   * Exception
   * @ingroup Exceptions
   */
  DeclException0 (ExcGridsDontMatch);
  /**
   * The ::DoFHandler or hp::DoFHandler was not initialized with a finite
   * element. Please call DoFHandler::distribute_dofs() etc. first.
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcNoFESelected);
  /**
   * @todo Write description
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcInvalidBoundaryIndicator);
}



/* ------------------------- inline functions -------------- */

#ifndef DOXYGEN

namespace DoFTools
{
  /**
   * Operator computing the maximum coupling out of two.
   *
   * @relates DoFTools
   */
  inline
  Coupling operator |= (Coupling &c1,
                        const Coupling c2)
  {
    if (c2 == always)
      c1 = always;
    else if (c1 != always && c2 == nonzero)
      return c1 = nonzero;
    return c1;
  }


  /**
   * Operator computing the maximum coupling out of two.
   *
   * @relates DoFTools
   */
  inline
  Coupling operator | (const Coupling c1,
                       const Coupling c2)
  {
    if (c1 == always || c2 == always)
      return always;
    if (c1 == nonzero || c2 == nonzero)
      return nonzero;
    return none;
  }


// ---------------------- inline and template functions --------------------

  template <int dim, int spacedim>
  inline
  unsigned int
  max_dofs_per_cell (const DoFHandler<dim,spacedim> &dh)
  {
    return dh.get_fe().dofs_per_cell;
  }


  template <int dim, int spacedim>
  inline
  unsigned int
  max_dofs_per_face (const DoFHandler<dim,spacedim> &dh)
  {
    return dh.get_fe().dofs_per_face;
  }


  template <int dim, int spacedim>
  inline
  unsigned int
  max_dofs_per_vertex (const DoFHandler<dim,spacedim> &dh)
  {
    return dh.get_fe().dofs_per_vertex;
  }


  template <int dim, int spacedim>
  inline
  unsigned int
  n_components (const DoFHandler<dim,spacedim> &dh)
  {
    return dh.get_fe().n_components();
  }



  template <int dim, int spacedim>
  inline
  bool
  fe_is_primitive (const DoFHandler<dim,spacedim> &dh)
  {
    return dh.get_fe().is_primitive();
  }


  template <int dim, int spacedim>
  inline
  unsigned int
  max_dofs_per_cell (const hp::DoFHandler<dim,spacedim> &dh)
  {
    return dh.get_fe().max_dofs_per_cell ();
  }


  template <int dim, int spacedim>
  inline
  unsigned int
  max_dofs_per_face (const hp::DoFHandler<dim,spacedim> &dh)
  {
    return dh.get_fe().max_dofs_per_face ();
  }


  template <int dim, int spacedim>
  inline
  unsigned int
  max_dofs_per_vertex (const hp::DoFHandler<dim,spacedim> &dh)
  {
    return dh.get_fe().max_dofs_per_vertex ();
  }


  template <int dim, int spacedim>
  inline
  unsigned int
  n_components (const hp::DoFHandler<dim,spacedim> &dh)
  {
    return dh.get_fe()[0].n_components();
  }


  template <int dim, int spacedim>
  inline
  bool
  fe_is_primitive (const hp::DoFHandler<dim,spacedim> &dh)
  {
    return dh.get_fe()[0].is_primitive();
  }


  template <typename DoFHandlerType, class Comp>
  void
  map_support_points_to_dofs
  (
    const Mapping<DoFHandlerType::dimension,DoFHandlerType::space_dimension>        &mapping,
    const DoFHandlerType                                                            &dof_handler,
    std::map<Point<DoFHandlerType::space_dimension>, types::global_dof_index, Comp> &point_to_index_map)
  {
    // let the checking of arguments be
    // done by the function first
    // called
    std::vector<Point<DoFHandlerType::space_dimension> > support_points (dof_handler.n_dofs());
    map_dofs_to_support_points (mapping, dof_handler, support_points);
    // now copy over the results of the
    // previous function into the
    // output arg
    point_to_index_map.clear ();
    for (types::global_dof_index i=0; i<dof_handler.n_dofs(); ++i)
      point_to_index_map[support_points[i]] = i;
  }
}

#endif

DEAL_II_NAMESPACE_CLOSE

#endif
