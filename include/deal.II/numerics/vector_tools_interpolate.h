// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_vector_tools_interpolate_h
#define dealii_vector_tools_interpolate_h

#include <deal.II/base/config.h>

#include <deal.II/fe/component_mask.h>

#include <map>

DEAL_II_NAMESPACE_OPEN

template <typename number>
class AffineConstraints;
template <int dim, int spacedim>
class DoFHandler;
template <typename number>
class FullMatrix;
template <int dim, typename Number>
class Function;
template <class MeshType>
class InterGridMap;
template <int dim, int spacedim>
class Mapping;

namespace VectorTools
{
  /**
   * @name Interpolation and projection
   */
  //@{

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
   * The template argument <code>DoFHandlerType</code> may either be of type
   * DoFHandler or hp::DoFHandler.
   *
   * See the general documentation of this namespace for further information.
   *
   * @todo The @p mapping argument should be replaced by a
   * hp::MappingCollection in case of a hp::DoFHandler.
   */
  template <int dim,
            int spacedim,
            typename VectorType,
            template <int, int> class DoFHandlerType>
  void
  interpolate(
    const Mapping<dim, spacedim> &                             mapping,
    const DoFHandlerType<dim, spacedim> &                      dof,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType &                                               vec,
    const ComponentMask &component_mask = ComponentMask());

  /**
   * Call the @p interpolate() function above with
   * <tt>mapping=MappingQGeneric1@<dim>@()</tt>.
   */
  template <int dim,
            int spacedim,
            typename VectorType,
            template <int, int> class DoFHandlerType>
  void
  interpolate(
    const DoFHandlerType<dim, spacedim> &                      dof,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType &                                               vec,
    const ComponentMask &component_mask = ComponentMask());

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
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void
  interpolate(const DoFHandler<dim, spacedim> &dof_1,
              const DoFHandler<dim, spacedim> &dof_2,
              const FullMatrix<double> &       transfer,
              const InVector &                 data_1,
              OutVector &                      data_2);

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
   * @author Valentin Zingan, 2013
   */
  template <int dim,
            int spacedim,
            typename VectorType,
            template <int, int> class DoFHandlerType>
  void
  interpolate_based_on_material_id(
    const Mapping<dim, spacedim> &       mapping,
    const DoFHandlerType<dim, spacedim> &dof_handler,
    const std::map<types::material_id,
                   const Function<spacedim, typename VectorType::value_type> *>
      &                  function_map,
    VectorType &         dst,
    const ComponentMask &component_mask = ComponentMask());

  /**
   * Compute the interpolation of a @p dof1-function @p u1 to a @p dof2-function
   * @p u2, where @p dof1 and @p dof2 represent different triangulations with
   * a common coarse grid.
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
   * yourself by calling the @p AffineConstraints::distribute function of your
   * hanging node constraints object.
   *
   * @note This function works with parallel::distributed::Triangulation, but
   * only if the parallel partitioning is the same for both meshes (see the
   * parallel::distributed::Triangulation<dim>::no_automatic_repartitioning
   * flag).
   */
  template <int dim,
            int spacedim,
            typename VectorType,
            template <int, int> class DoFHandlerType>
  void
  interpolate_to_different_mesh(const DoFHandlerType<dim, spacedim> &dof1,
                                const VectorType &                   u1,
                                const DoFHandlerType<dim, spacedim> &dof2,
                                VectorType &                         u2);

  /**
   * Compute the interpolation of a @p dof1-function @p u1 to a @p dof2-function
   * @p u2, where @p dof1 and @p dof2 represent different triangulations with
   * a common coarse grid.
   *
   * dof1 and dof2 need to have the same finite element discretization.
   *
   * @p constraints is a hanging node constraints object corresponding to @p
   * dof2. This object is particularly important when interpolating onto
   * continuous elements on grids with hanging nodes (locally refined grids):
   * Without it - due to cellwise interpolation - the resulting output vector
   * does not necessarily respect continuity requirements at hanging nodes.
   */
  template <int dim,
            int spacedim,
            typename VectorType,
            template <int, int> class DoFHandlerType>
  void
  interpolate_to_different_mesh(
    const DoFHandlerType<dim, spacedim> &                     dof1,
    const VectorType &                                        u1,
    const DoFHandlerType<dim, spacedim> &                     dof2,
    const AffineConstraints<typename VectorType::value_type> &constraints,
    VectorType &                                              u2);

  /**
   * The same function as above, but takes an InterGridMap object directly as
   * a parameter. Useful for interpolating several vectors at the same time.
   *
   * @p intergridmap has to be initialized via InterGridMap::make_mapping
   * pointing from a source DoFHandler to a destination DoFHandler.
   */
  template <int dim,
            int spacedim,
            typename VectorType,
            template <int, int> class DoFHandlerType>
  void
  interpolate_to_different_mesh(
    const InterGridMap<DoFHandlerType<dim, spacedim>> &       intergridmap,
    const VectorType &                                        u1,
    const AffineConstraints<typename VectorType::value_type> &constraints,
    VectorType &                                              u2);

  //@}

  /**
   * Geometrical interpolation
   */
  //@{
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
   * components of the FiniteElement to use to describe the geometry. If no
   * mask is specified at construction time, then a default one is used, i.e.,
   * the first spacedim components of the FiniteElement are assumed to
   * represent the geometry of the problem.
   *
   * This function is only implemented for FiniteElements where the specified
   * components are primitive.
   *
   * @author Luca Heltai, 2015
   */
  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename VectorType>
  void
  get_position_vector(const DoFHandlerType<dim, spacedim> &dh,
                      VectorType &                         vector,
                      const ComponentMask &mask = ComponentMask());

  //@}
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_interpolate_h
