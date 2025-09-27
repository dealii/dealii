// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_vector_tools_interpolate_h
#define dealii_vector_tools_interpolate_h

#include <deal.II/base/config.h>

#include <deal.II/base/template_constraints.h>
#include <deal.II/base/types.h>

#include <deal.II/fe/component_mask.h>

#include <map>

DEAL_II_NAMESPACE_OPEN

template <typename number>
class AffineConstraints;

template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler;

template <typename number>
class FullMatrix;
template <int dim, typename Number>
class Function;
template <typename MeshType>
DEAL_II_CXX20_REQUIRES(concepts::is_triangulation_or_dof_handler<MeshType>)
class InterGridMap;
template <int dim, int spacedim>
class Mapping;

namespace hp
{
  template <int dim, int spacedim>
  class MappingCollection;
}

namespace VectorTools
{
  /**
   * @name Interpolation and projection
   */
  /** @{ */

  /**
   * Compute the interpolation of @p function at the support points to the
   * finite element space described by the Triangulation and FiniteElement
   * object with which the given DoFHandler argument is initialized. It is
   * assumed that the number of components of @p function matches that of the
   * finite element used by @p dof.
   *
   * Note that you may have to call <tt>hanging_nodes.distribute(vec)</tt>
   * with the hanging nodes from space @p dof afterwards, to make the result
   * continuous again.
   *
   * See the general documentation of this namespace for further information.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate(
    const Mapping<dim, spacedim>                              &mapping,
    const DoFHandler<dim, spacedim>                           &dof,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType                                                &vec,
    const ComponentMask &component_mask = {},
    const unsigned int   level          = numbers::invalid_unsigned_int);

  /**
   * Same as above but in an hp-context.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate(
    const hp::MappingCollection<dim, spacedim>                &mapping,
    const DoFHandler<dim, spacedim>                           &dof,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType                                                &vec,
    const ComponentMask &component_mask = {},
    const unsigned int   level          = numbers::invalid_unsigned_int);


  /**
   * Call the @p interpolate() function above with
   * <tt>mapping=MappingQ@<dim,spacedim@>(1)</tt>.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate(
    const DoFHandler<dim, spacedim>                           &dof,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType                                                &vec,
    const ComponentMask &component_mask = {},
    const unsigned int   level          = numbers::invalid_unsigned_int);

  /**
   * Interpolate different finite element spaces. The interpolation of vector
   * @p data_1 (which is assumed to be ghosted, see
   * @ref GlossGhostedVector)
   * is executed from the FE space represented by @p dof_1
   * to the vector @p data_2 on FE space @p dof_2.
   * The interpolation on each cell is represented by the matrix @p transfer.
   * Curved boundaries are neglected so far.
   *
   * Note that you may have to call <tt>hanging_nodes.distribute(data_2)</tt>
   * with the hanging nodes from space @p dof_2 afterwards, to make the result
   * continuous again.
   *
   * @note Instantiations for this template are provided for some vector types
   * (see the general documentation of the namespace), but only the same
   * vector for InVector and OutVector. Other combinations must be
   * instantiated by hand.
   *
   * @dealiiConceptRequires{concepts::is_dealii_vector_type<InVector>
   *   &&concepts::is_writable_dealii_vector_type<OutVector>}
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  DEAL_II_CXX20_REQUIRES(concepts::is_dealii_vector_type<InVector> &&
                           concepts::is_writable_dealii_vector_type<OutVector>)
  void interpolate(const DoFHandler<dim, spacedim> &dof_1,
                   const DoFHandler<dim, spacedim> &dof_2,
                   const FullMatrix<double>        &transfer,
                   const InVector                  &data_1,
                   OutVector                       &data_2);

  /**
   * This function is a kind of generalization or modification of the very
   * first interpolate() function in the series. It interpolates a set of
   * functions onto the finite element space defined by the DoFHandler argument,
   * where the determination which function to use on each cell is made
   * based on the material id (see
   * @ref GlossMaterialId)
   * of each cell.
   *
   * @param[in] mapping        The mapping to use to determine the location of
   *   support points at which the functions are to be evaluated.
   * @param[in] dof_handler    DoFHandler initialized with Triangulation and
   *   FiniteElement objects and that defines the finite element space.
   * @param[in] function_map   A std::map reflecting the correspondence between
   *   material ids on those cells on which something should be interpolated,
   *   and the functions to be interpolated onto the finite element space.
   * @param[out] dst           The global finie element vector holding the
   *   output of the interpolated values.
   * @param[in] component_mask A mask of components that shall be interpolated.
   *
   * @note If the algorithm encounters a cell whose material id is not listed
   * in the given @p function_map, then @p dst will not be updated in the
   * respective degrees of freedom of the output vector. For example, if
   * @p dst was initialized to zero, then those zeros which correspond to
   * the missed material ids will still remain in @p dst after calling
   * this function.
   *
   * @note Degrees of freedom located on faces between cells of different
   * material ids will get their value by that cell which was called last in
   * the respective loop over cells implemented in this function. Since the
   * order of cells is somewhat arbitrary, you cannot control it. However, if
   * you want to have control over the order in which cells are visited, let us
   * take a
   * look at the following example: Let @p u be a variable of interest which
   * is approximated by some CG finite element. Let @p 0, @p 1 and @p 2 be
   * material ids of cells on the triangulation. Let 0: 0.0, 1: 1.0, 2: 2.0 be
   * the whole @p function_map that you want to pass to this function, where
   * @p key is a material id and @p value is a value of @p u. By using the
   * whole @p function_map you do not really know which values will be
   * assigned to the face DoFs. On the other hand, if you split the whole @p
   * function_map into three smaller independent objects 0: 0.0 and 1: 1.0 and
   * 2: 2.0 and make three distinct calls of this function passing each of
   * these objects separately (the order depends on what you want to get
   * between cells), then each subsequent call will rewrite the intercell @p
   * dofs of the previous one.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate_based_on_material_id(
    const Mapping<dim, spacedim>    &mapping,
    const DoFHandler<dim, spacedim> &dof_handler,
    const std::map<types::material_id,
                   const Function<spacedim, typename VectorType::value_type> *>
                        &function_map,
    VectorType          &dst,
    const ComponentMask &component_mask = {});

  /**
   * Compute the interpolation of a @p dof_handler_1 function @p u1 (i.e.,
   * a function defined on the mesh that underlies @p dof_handler_1, using
   * the finite element associated with that DoFHandler) to
   * a @p dof_handler_2 function @p u2 (i.e., a function defined on the mesh
   * and finite element associated with @p dof_handler_2). This function
   * assumes that @p dof_handler_1 and @p dof_handler_2 are defined on
   * (possibly different) triangulations that have a common coarse grid.
   *
   * dof1 and dof2 need to have the same finite element discretization.
   *
   * Note that for continuous elements on grids with hanging nodes (i.e.
   * locally refined grids) this function does not give the expected output.
   * Indeed, the resulting output vector does not necessarily respect
   * continuity requirements at hanging nodes, due to local cellwise
   * interpolation.
   *
   * For this case (continuous elements on grids with hanging nodes), please
   * use the interpolate_to_different_mesh function with an additional
   * AffineConstraints argument, see below, or make the field conforming
   * yourself by calling the AffineConstraints::distribute() function of your
   * hanging node constraints object.
   *
   * @note This function works with parallel::distributed::Triangulation, but
   * only if the parallel partitioning is the same for both meshes (see the
   * parallel::distributed::Triangulation<dim>::no_automatic_repartitioning
   * flag). In practice, this is rarely the case because two triangulations,
   * partitioned in their own ways, will not typically have corresponding
   * cells owned by the same process, and implementing the interpolation
   * procedure would require transferring data between processes in ways
   * that are difficult to implement efficiently. However, some special
   * cases can more easily be implemented, namely the case where one
   * of the meshes is strictly coarser or finer than the other. For these
   * cases, see the interpolate_to_coarser_mesh() and
   * interpolate_to_finer_mesh().
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate_to_different_mesh(
    const DoFHandler<dim, spacedim> &dof_handler_1,
    const VectorType                &u1,
    const DoFHandler<dim, spacedim> &dof_handler_2,
    VectorType                      &u2);

  /**
   * This is a variation of the previous function that takes an additional
   * constraint object as argument.
   *
   * @p constraints is a hanging node constraints object corresponding to
   * @p dof2. This object is particularly important when interpolating onto
   * continuous elements on grids with hanging nodes (locally refined grids):
   * Without it - due to cellwise interpolation - the resulting output vector
   * does not necessarily respect continuity requirements at hanging nodes.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate_to_different_mesh(
    const DoFHandler<dim, spacedim>                          &dof_handler_1,
    const VectorType                                         &u1,
    const DoFHandler<dim, spacedim>                          &dof_handler_2,
    const AffineConstraints<typename VectorType::value_type> &constraints,
    VectorType                                               &u2);

  /**
   * The same function as above, but takes an InterGridMap object directly as
   * a parameter. This function is useful for interpolating several vectors at
   * the same time.
   *
   * @p intergridmap has to be initialized via InterGridMap::make_mapping
   * pointing from a source DoFHandler to a destination DoFHandler.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate_to_different_mesh(
    const InterGridMap<DoFHandler<dim, spacedim>>            &intergridmap,
    const VectorType                                         &u1,
    const AffineConstraints<typename VectorType::value_type> &constraints,
    VectorType                                               &u2);

  /**
   * This function addresses one of the limitations of the
   * interpolate_to_different_mesh() function, namely that the latter only
   * works on parallel triangulations in very specific cases where both
   * triangulations involved happen to be partitioned in such a way that
   * if a process owns a cell of one mesh, it also needs to own the
   * corresponding parent of child cells of the other mesh. In practice, this is
   * rarely the case.
   *
   * This function does not have this restriction, and consequently also works
   * in parallel, as long as the mesh we interpolate onto can be obtained from
   * the mesh we interpolate off by *coarsening*. In other words, the target
   * mesh is coarser than the source mesh.
   *
   * The function takes an additional constraints argument that is used after
   * interpolation to ensure that the output vector is conforming (that is,
   * that all entries of the output vector conform to their constraints).
   * @p constraints_coarse therefore needs to correspond to
   * @p dof_handler_coarse .
   *
   * The opposite operation, interpolation from a coarser to a finer mesh,
   * is implemented in the interpolate_to_finer_mesh() function.
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate_to_coarser_mesh(
    const DoFHandler<dim, spacedim> &dof_handler_fine,
    const VectorType                &u_fine,
    const DoFHandler<dim, spacedim> &dof_handler_coarse,
    const AffineConstraints<typename VectorType::value_type>
               &constraints_coarse,
    VectorType &u_coarse);

  /**
   * This function addresses one of the limitations of the
   * interpolate_to_different_mesh() function, namely that the latter only
   * works on parallel triangulations in very specific cases where both
   * triangulations involved happen to be partitioned in such a way that
   * if a process owns a cell of one mesh, it also needs to own the
   * corresponding parent of child cells of the other mesh. In practice, this is
   * rarely the case.
   *
   * This function does not have this restriction, and consequently also works
   * in parallel, as long as the mesh we interpolate onto can be obtained from
   * the mesh we interpolate off by *refinement*. In other words, the target
   * mesh is finer than the source mesh.
   *
   * The function takes an additional constraints argument that is used after
   * interpolation to ensure that the output vector is conforming (that is,
   * that all entries of the output vector conform to their constraints).
   * @p constraints_finee therefore needs to correspond to
   * @p dof_handler_fine .
   *
   * The opposite operation, interpolation from a finer to a coarser mesh,
   * is implemented in the interpolate_to_coarser_mesh() function.
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void interpolate_to_finer_mesh(
    const DoFHandler<dim, spacedim> &dof_handler_coarse,
    const VectorType                &u_coarse,
    const DoFHandler<dim, spacedim> &dof_handler_fine,
    const AffineConstraints<typename VectorType::value_type> &constraints_fine,
    VectorType                                               &u_fine);

  /** @} */

  /**
   * Geometrical interpolation
   */
  /** @{ */
  /**
   * Given a DoFHandler containing at least a spacedim vector field, this
   * function interpolates the Triangulation at the support points of a FE_Q()
   * finite element of the same degree as the degree of the required
   * components.
   *
   * Curved manifold are respected, and the resulting VectorType will be
   * geometrically consistent. The resulting map is guaranteed to be
   * interpolatory at the support points of a FE_Q() finite element of the
   * same degree as the degree of the required components.
   *
   * If the underlying finite element is an FE_Q(1)^spacedim, then the
   * resulting @p VectorType is a finite element field representation of the
   * vertices of the Triangulation.
   *
   * The optional ComponentMask argument can be used to specify what
   * components of the FiniteElement to use to describe the
   * geometry. If no mask is specified at construction time, then a
   * default-constructed mask is used, which is then interpreted as
   * saying that the first `spacedim` components of the FiniteElement
   * are assumed to represent the geometry of the problem.
   *
   * This function is only implemented for FiniteElements where the specified
   * components are primitive.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void get_position_vector(const DoFHandler<dim, spacedim> &dh,
                           VectorType                      &vector,
                           const ComponentMask             &mask = {});

  /**
   * Like the above function but also taking @p mapping as argument.
   * This will introduce an additional approximation between the true geometry
   * specified by the manifold if the degree of the mapping is lower than the
   * degree of the finite element in the DoFHandler @p dh, but more
   * importantly it allows to fill location vectors for mappings that do not
   * preserve vertex locations (like Eulerian mappings).
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void get_position_vector(const Mapping<dim, spacedim>    &mapping,
                           const DoFHandler<dim, spacedim> &dh,
                           VectorType                      &vector,
                           const ComponentMask             &mask = {});

  /** @} */
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_interpolate_h
