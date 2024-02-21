// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_non_matching_fe_immersed_values_h
#define dealii_non_matching_fe_immersed_values_h

#include <deal.II/base/config.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria_iterator.h>

#include <deal.II/non_matching/immersed_surface_quadrature.h>

DEAL_II_NAMESPACE_OPEN
namespace NonMatching
{
  /**
   * Finite element evaluated in the quadrature points of an
   * ImmersedSurfaceQuadrature of a cell.
   *
   * The shape functions values and their derivatives are the same as for an
   * FEValues-object, but the JxW-values are computed with the transformation
   * described in the documentation of ImmersedSurfaceQuadrature. Further, the
   * normal_vector-function returns the normal to the immersed surface.
   *
   * The reinit-function of this class exist mostly to be consistent with the
   * other FEValues-like classes. The immersed quadrature rule will typically
   * vary between each cell of the triangulation. Thus, an
   * FEImmersedSurfaceValues object can, typically, not be reused for different
   * cells.
   *
   * See also documentation in FEValuesBase.
   *
   * @ingroup feaccess
   */
  template <int dim>
  class FEImmersedSurfaceValues : public FEValuesBase<dim, dim>
  {
  public:
    /**
     * Constructor. Gets cell-independent data from mapping and finite element
     * objects, matching the quadrature rule and update flags.
     *
     * @note Currently this class is only implemented for MappingCartesian,
     * MappingQ and MappingFEField.
     */
    FEImmersedSurfaceValues(const Mapping<dim>                   &mapping,
                            const FiniteElement<dim>             &element,
                            const ImmersedSurfaceQuadrature<dim> &quadrature,
                            const UpdateFlags                     update_flags);

    /**
     * Reinitialize quantities (normals, JxW-values, etc) for the given cell of
     * type "iterator into a Triangulation object".
     */
    void
    reinit(const typename Triangulation<dim>::cell_iterator &cell);

    /**
     * Reinitialize quantities (shape function values, gradients, etc) for the
     * given cell of type "iterator into a DoFHandler object", and the finite
     * element associated with this object.
     */
    template <bool level_dof_access>
    void
    reinit(
      const TriaIterator<DoFCellAccessor<dim, dim, level_dof_access>> &cell);

    /**
     * Returns the surface gradient of the shape function with index
     * @p function_no at the quadrature point with index @p quadrature_point.
     *
     * The surface gradient is defined as the projection of the gradient to the
     * tangent plane of the surface:
     * $ \nabla u - (n \cdot \nabla u) n $,
     * where $n$ is the unit normal to the surface.
     *
     * @dealiiRequiresUpdateFlags{update_gradients | update_normal_vectors}
     */
    Tensor<1, dim>
    shape_surface_grad(const unsigned int function_no,
                       const unsigned int quadrature_point) const;

    /**
     * Return one vector component of the surface gradient of the shape function
     * at a quadrature point. See the definition of the surface gradient in the
     * shape_surface_grad function.
     *
     * @p function_no Index of the shape function to be evaluated.
     *
     * @p point_no Index of the quadrature point at which the function is to be
     * evaluated.
     *
     * @p component Vector component to be evaluated.
     *
     * @dealiiRequiresUpdateFlags{update_gradients | update_normal_vectors}
     */
    Tensor<1, dim>
    shape_surface_grad_component(const unsigned int function_no,
                                 const unsigned int quadrature_point,
                                 const unsigned int component) const;

    /**
     * Return a reference to the copy of the quadrature formula stored by this
     * object.
     */
    const NonMatching::ImmersedSurfaceQuadrature<dim> &
    get_quadrature() const;

  protected:
    /**
     * Do work common to the constructors.
     */
    void
    initialize(const UpdateFlags update_flags);

    /**
     * The reinit() functions do only that part of the work that requires
     * knowledge of the type of iterator. After setting present_cell(), they
     * pass on to this function, which does the real work, and which is
     * independent of the actual type of the cell iterator.
     */
    void
    do_reinit();

    /**
     * Copy of the quadrature formula that was passed to the constructor.
     */
    const ImmersedSurfaceQuadrature<dim> quadrature;
  };

} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif /* dealii_non_matching_fe_immersed_values_h */
