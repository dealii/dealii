// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#ifndef dealii_fe_interface_values_h
#define dealii_fe_interface_values_h

#include <deal.II/base/config.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

DEAL_II_NAMESPACE_OPEN

/**
 * FEInterfaceValues is a data structure to access and assemble finite element
 * data on interfaces between two cells of a mesh.
 *
 * It provides a way to access averages, jump terms, and similar operations used
 * in Discontinuous Galerkin methods on a face between two neighboring cells.
 * This allows the computation of typical mesh-dependent linear or bilinear
 * forms in a similar way as FEValues does for cells and FEFaceValues does for
 * faces. In
 * the literature, the faces between neighboring cells are called "inner
 * interfaces" or "facets".
 *
 * Internally, this class provides an abstraction for two FEFaceValues
 * objects (or FESubfaceValues when using adaptive refinement). The class
 * introduces a new "interface dof index" that walks over
 * the union of the dof indices of the two FEFaceValues objects. Helper
 * functions allow translating between the new "interface dof index" and the
 * corresponding "cell index" (0 for the first cell, 1 for the second cell)
 * and "dof index" within that cell.
 *
 * The class is made to be used inside MeshWorker::mesh_loop(). It is intended
 * to be a low level replacement for MeshWorker and LocalIntegrators and a
 * higher level abstraction compared to assembling face terms manually.
 *
 * @author Timo Heister, 2019.
 */
template <int dim, int spacedim = dim>
class FEInterfaceValues
{
public:
  /**
   * Number of quadrature points.
   */
  const unsigned int n_quadrature_points;


  /**
   * Construct the FEInterfaceValues with a single FiniteElement (same on both
   * sides of the facet). The FEFaceValues objects will be initialized with
   * the given @p mapping, @p quadrature, and @p update_flags.
   */
  FEInterfaceValues(const Mapping<dim, spacedim> &      mapping,
                    const FiniteElement<dim, spacedim> &fe,
                    const Quadrature<dim - 1> &         quadrature,
                    const UpdateFlags                   update_flags);

  /**
   * Construct the FEInterfaceValues with a single FiniteElement and
   * a Q1 Mapping.
   *
   * See the constructor above.
   */
  FEInterfaceValues(const FiniteElement<dim, spacedim> &fe,
                    const Quadrature<dim - 1> &         quadrature,
                    const UpdateFlags                   update_flags);

  /**
   * Re-initialize this object to be used on a new interface given by two faces
   * of two neighboring cells. The `cell` and `cell_neighbor` cells will be
   * referred to through `cell_index` zero and one after this call in all places
   * where one needs to identify the two cells adjacent to the interface.
   *
   * Use numbers::invalid_unsigned_int for @p sub_face_no or @p
   * sub_face_no_neighbor to indicate that you want to work on the entire face,
   * not a sub-face.
   *
   * The arguments (including their order) are identical to the @p face_worker arguments
   * in MeshWorker::mesh_loop().
   *
   * @param[in] cell An iterator to the first cell adjacent to the interface.
   * @param[in] face_no An integer identifying which face of the first cell the
   *   interface is on.
   * @param[in] sub_face_no An integer identifying the subface (child) of the
   *   face (identified by the previous two arguments) that the interface
   *   corresponds to. If equal to numbers::invalid_unsigned_int, then the
   *   interface is considered to be the entire face.
   * @param[in] cell_neighbor An iterator to the second cell adjacent to
   *   the interface. The type of this iterator does not have to equal that
   *   of `cell`, but must be convertible to it. This allows using an
   *   active cell iterator for `cell`, and `cell->neighbor(f)` for
   *   `cell_neighbor`, since the return type of `cell->neighbor(f)` is
   *   simply a cell iterator (not necessarily an active cell iterator).
   * @param[in] face_no_neighbor Like `face_no`, just for the neighboring
   *   cell.
   * @param[in] sub_face_no_neighbor Like `sub_face_no`, just for the
   *   neighboring cell.
   */
  template <class CellIteratorType>
  void
  reinit(const CellIteratorType &                         cell,
         const unsigned int                               face_no,
         const unsigned int                               sub_face_no,
         const typename identity<CellIteratorType>::type &cell_neighbor,
         const unsigned int                               face_no_neighbor,
         const unsigned int                               sub_face_no_neighbor);

  /**
   * Re-initialize this object to be used on an interface given by a single face
   * @p face_no of the cell @p cell. This is useful to use FEInterfaceValues
   * on boundaries of the domain.
   *
   * As a consequence, members like jump() will assume a value of zero for the
   * values on the "other" side. Note that no sub_face_number is needed as a
   * boundary face can not neighbor a finer cell.
   *
   * After calling this function at_boundary() will return true.
   */
  template <class CellIteratorType>
  void
  reinit(const CellIteratorType &cell, const unsigned int face_no);

  /**
   * Return a reference to the FEFaceValues or FESubfaceValues object
   * of the specified cell of the interface.
   *
   * The @p cell_index is either 0 or 1 and corresponds to the cell index
   * returned by interface_dof_to_cell_and_dof_index().
   */
  const FEFaceValuesBase<dim, spacedim> &
  get_fe_face_values(const unsigned int cell_index) const;

  /**
   * Return a reference to the quadrature object in use.
   */
  const Quadrature<dim - 1> &
  get_quadrature() const;

  /**
   * Return the update flags set.
   */
  UpdateFlags
  get_update_flags() const;

  /**
   * @name Functions to query information on a given interface
   * @{
   */

  /**
   * Return if the current interface is a boundary face or an internal
   * face with two adjacent cells.
   *
   * See the corresponding reinit() functions for details.
   */
  bool
  at_boundary() const;

  /**
   * Mapped quadrature weight. This value equals the
   * mapped surface element times the weight of the quadrature
   * point.
   *
   * You can think of the quantity returned by this function as the
   * surface element $ds$ in the integral that we implement here by
   * quadrature.
   *
   * @dealiiRequiresUpdateFlags{update_JxW_values}
   */
  double
  JxW(const unsigned int quadrature_point) const;

  /**
   * Return the vector of JxW values for each quadrature point.
   *
   * @dealiiRequiresUpdateFlags{update_JxW_values}
   */
  const std::vector<double> &
  get_JxW_values() const;

  /**
   * Return the normal vector of the interface in each quadrature point.
   *
   * The return value is identical to get_fe_face_values(0).get_normal_vectors()
   * and therefore, are outside normal vectors from the perspective of the
   * first cell of this interface.
   *
   * @dealiiRequiresUpdateFlags{update_normal_vectors}
   */
  const std::vector<Tensor<1, spacedim>> &
  get_normal_vectors() const;

  /**
   * Return a reference to the quadrature points in real space.
   *
   * @dealiiRequiresUpdateFlags{update_quadrature_points}
   */
  const std::vector<Point<spacedim>> &
  get_quadrature_points() const;


  /**
   * Return the number of DoFs (or shape functions) on the current interface.
   *
   * @note This number is only available after a call to reinit() and can change
   * from one call to reinit() to the next. For example, on a boundary interface
   * it is equal to the number of dofs of the single FEFaceValues object, while
   * it is twice that for an interior interface for a DG element. For a
   * continuous element, it is slightly smaller because the two cells on the
   * interface share some of the dofs.
   */
  unsigned
  n_current_interface_dofs() const;

  /**
   * Return the set of joint DoF indices. This includes indices from both cells.
   * If reinit was called with an active cell iterator, the indices are based
   * on the active indices (returned by `DoFCellAccessor::get_dof_indices()` ),
   * in case of level cell (that is, if is_level_cell() return true )
   * the mg dof indices are returned.
   *
   * @note This function is only available after a call to reinit() and can change
   * from one call to reinit() to the next.
   */
  std::vector<types::global_dof_index>
  get_interface_dof_indices() const;

  /**
   * Convert an interface dof index into the corresponding local DoF indices of
   * the two cells. If an interface DoF is only active on one of the
   * cells, the other index will be numbers::invalid_unsigned_int.
   *
   * For discontinuous finite elements each interface dof will correspond to
   * exactly one DoF index.
   *
   * @note This function is only available after a call to reinit() and can change
   * from one call to reinit() to the next.
   */
  std::array<unsigned int, 2>
  interface_dof_to_dof_indices(const unsigned int interface_dof_index) const;

  /**
   * Return the normal in a given quadrature point.
   *
   * The normal points in outwards direction as seen from the first cell of
   * this interface.
   *
   * @dealiiRequiresUpdateFlags{update_normal_vectors}
   */
  Tensor<1, spacedim>
  normal(const unsigned int q_point_index) const;

  /**
   * @}
   */

  /**
   * @name Functions to evaluate data of the shape functions
   * @{
   */

  /**
   * Return component @p component of the value of the shape function
   * with interface dof index @p interface_dof_index in
   * quadrature point @p q_point.
   *
   * The argument @p here_or_there selects between the value on cell 0 (here, @p true)
   * and cell 1 (there, @p false). You can also interpret it as "upstream" (@p true)
   * and "downstream" (@p false) as defined by the direction of the normal
   * vector
   * in this quadrature point. If @p here_or_there is true, the shape
   * functions from the first cell of the interface is used.
   *
   * In other words, this function returns the limit of the value of the shape
   * function in the given quadrature point when approaching it from one of the
   * two cells of the interface.
   *
   * @note This function is typically used to pick the upstream or downstream
   * value based on a direction. This can be achieved by using
   * <code>(direction * normal)>0</code> as the first argument of this
   * function.
   */
  double
  shape_value(const bool         here_or_there,
              const unsigned int interface_dof_index,
              const unsigned int q_point,
              const unsigned int component = 0) const;

  /**
   * Return the jump $\jump{u}=u_{\text{cell0}} - u_{\text{cell1}}$ on the
   * interface
   * for the shape function @p interface_dof_index at the quadrature point
   * @p q_point of component @p component.
   *
   * Note that one can define the jump in
   * different ways (the value "there" minus the value "here", or the other way
   * around; both are used in the finite element literature). The definition
   * here uses "value here minus value there", as seen from the first cell.
   *
   * If this is a boundary face (at_boundary() returns true), then
   * $\jump{u}=u_{\text{cell0}}$.
   */
  double
  jump(const unsigned int interface_dof_index,
       const unsigned int q_point,
       const unsigned int component = 0) const;

  /**
   * Return the average $\average{u}=\frac{1}{2}u_{\text{cell0}} +
   * \frac{1}{2}u_{\text{cell1}}$ on the interface
   * for the shape function @p interface_dof_index at the quadrature point
   * @p q_point of component @p component.
   *
   * If this is a boundary face (at_boundary() returns true), then
   * $\average{u}=u_{\text{cell0}}$.
   */
  double
  average(const unsigned int interface_dof_index,
          const unsigned int q_point,
          const unsigned int component = 0) const;

  /**
   * Return the average of the gradient $\average{\nabla u} = \frac{1}{2}\nabla
   * u_{\text{cell0}} + \frac{1}{2} \nabla u_{\text{cell1}}$ on the interface
   * for the shape function @p interface_dof_index at the quadrature point @p
   * q_point of component @p component.
   *
   * If this is a boundary face (at_boundary() returns true), then
   * $\average{\nabla u}=\nabla u_{\text{cell0}}$.
   */
  Tensor<1, dim>
  average_gradient(const unsigned int interface_dof_index,
                   const unsigned int q_point,
                   const unsigned int component = 0) const;

  /**
   * Return the average of the Hessian $\average{\nabla^2 u} =
   * \frac{1}{2}\nabla^2 u_{\text{cell0}} + \frac{1}{2} \nabla^2
   * u_{\text{cell1}}$ on the interface
   * for the shape function @p interface_dof_index at the quadrature point @p
   * q_point of component @p component.
   *
   * If this is a boundary face (at_boundary() returns true), then
   * $\average{\nabla^2 u}=\nabla^2 u_{\text{cell0}}$.
   */
  Tensor<2, dim>
  average_hessian(const unsigned int interface_dof_index,
                  const unsigned int q_point,
                  const unsigned int component = 0) const;

  /**
   * Return the jump in the gradient $\jump{\nabla u}=\nabla u_{\text{cell0}} -
   * \nabla u_{\text{cell1}}$ on the interface for the shape function @p
   * interface_dof_index at the quadrature point @p q_point of component @p
   * component.
   *
   * If this is a boundary face (at_boundary() returns true), then
   * $\jump{\nabla u}=\nabla u_{\text{cell0}}$.
   */
  Tensor<1, dim>
  jump_gradient(const unsigned int interface_dof_index,
                const unsigned int q_point,
                const unsigned int component = 0) const;

  /**
   * Return the jump in the Hessian $\jump{\nabla^2 u} = \nabla^2
   * u_{\text{cell0}} - \nabla^2 u_{\text{cell1}}$ on the interface for the
   * shape function
   * @p interface_dof_index at the quadrature point @p q_point of component
   * @p component.
   *
   * If this is a boundary face (at_boundary() returns true), then
   * $\jump{\nabla^2 u} = \nabla^2 u_{\text{cell0}}$.
   */
  Tensor<2, dim>
  jump_hessian(const unsigned int interface_dof_index,
               const unsigned int q_point,
               const unsigned int component = 0) const;

  /**
   * Return the jump in the third derivative $\jump{\nabla^3 u} = \nabla^3
   * u_{\text{cell0}} - \nabla^3 u_{\text{cell1}}$ on the interface for the
   * shape function @p interface_dof_index at the quadrature point @p q_point of
   * component @p component.
   *
   * If this is a boundary face (at_boundary() returns true), then
   * $\jump{\nabla^3 u} = \nabla^3 u_{\text{cell0}}$.
   */
  Tensor<3, dim>
  jump_3rd_derivative(const unsigned int interface_dof_index,
                      const unsigned int q_point,
                      const unsigned int component = 0) const;

  /**
   * @}
   */

private:
  /**
   * The list of DoF indices for the current interface, filled in reinit().
   */
  std::vector<types::global_dof_index> interface_dof_indices;

  /**
   * The mapping from interface dof to the two local dof indices of the
   * FeFaceValues objects. If an interface DoF is only active on one of the
   * cells, the other one will have numbers::invalid_unsigned_int.
   */
  std::vector<std::array<unsigned int, 2>> dofmap;

  /**
   * The FEFaceValues object for the current cell.
   */
  FEFaceValues<dim> internal_fe_face_values;

  /**
   * The FEFaceValues object for the current cell if the cell is refined.
   */
  FESubfaceValues<dim> internal_fe_subface_values;

  /**
   * The FEFaceValues object for the neighboring cell.
   */
  FEFaceValues<dim> internal_fe_face_values_neighbor;

  /**
   * The FEFaceValues object for the neighboring cell if the cell is refined.
   */
  FESubfaceValues<dim> internal_fe_subface_values_neighbor;

  /**
   * Pointer to internal_fe_face_values or internal_fe_subface_values,
   * respectively as determined in reinit().
   */
  FEFaceValuesBase<dim> *fe_face_values;

  /**
   * Pointer to internal_fe_face_values_neighbor,
   * internal_fe_subface_values_neighbor, or nullptr, respectively
   * as determined in reinit().
   */
  FEFaceValuesBase<dim> *fe_face_values_neighbor;
};



#ifndef DOXYGEN

/*---------------------- Inline functions ---------------------*/

template <int dim, int spacedim>
FEInterfaceValues<dim, spacedim>::FEInterfaceValues(
  const Mapping<dim, spacedim> &      mapping,
  const FiniteElement<dim, spacedim> &fe,
  const Quadrature<dim - 1> &         quadrature,
  const UpdateFlags                   update_flags)
  : n_quadrature_points(quadrature.size())
  , internal_fe_face_values(mapping, fe, quadrature, update_flags)
  , internal_fe_subface_values(mapping, fe, quadrature, update_flags)
  , internal_fe_face_values_neighbor(mapping, fe, quadrature, update_flags)
  , internal_fe_subface_values_neighbor(mapping, fe, quadrature, update_flags)
  , fe_face_values(nullptr)
  , fe_face_values_neighbor(nullptr)
{}



template <int dim, int spacedim>
FEInterfaceValues<dim, spacedim>::FEInterfaceValues(
  const FiniteElement<dim, spacedim> &fe,
  const Quadrature<dim - 1> &         quadrature,
  const UpdateFlags                   update_flags)
  : n_quadrature_points(quadrature.size())
  , internal_fe_face_values(StaticMappingQ1<dim, spacedim>::mapping,
                            fe,
                            quadrature,
                            update_flags)
  , internal_fe_subface_values(StaticMappingQ1<dim, spacedim>::mapping,
                               fe,
                               quadrature,
                               update_flags)
  , internal_fe_face_values_neighbor(StaticMappingQ1<dim, spacedim>::mapping,
                                     fe,
                                     quadrature,
                                     update_flags)
  , internal_fe_subface_values_neighbor(StaticMappingQ1<dim, spacedim>::mapping,
                                        fe,
                                        quadrature,
                                        update_flags)
  , fe_face_values(nullptr)
  , fe_face_values_neighbor(nullptr)
{}



template <int dim, int spacedim>
template <class CellIteratorType>
void
FEInterfaceValues<dim, spacedim>::reinit(
  const CellIteratorType &                         cell,
  const unsigned int                               face_no,
  const unsigned int                               sub_face_no,
  const typename identity<CellIteratorType>::type &cell_neighbor,
  const unsigned int                               face_no_neighbor,
  const unsigned int                               sub_face_no_neighbor)
{
  if (sub_face_no == numbers::invalid_unsigned_int)
    {
      internal_fe_face_values.reinit(cell, face_no);
      fe_face_values = &internal_fe_face_values;
    }
  else
    {
      internal_fe_subface_values.reinit(cell, face_no, sub_face_no);
      fe_face_values = &internal_fe_subface_values;
    }
  if (sub_face_no_neighbor == numbers::invalid_unsigned_int)
    {
      internal_fe_face_values_neighbor.reinit(cell_neighbor, face_no_neighbor);
      fe_face_values_neighbor = &internal_fe_face_values_neighbor;
    }
  else
    {
      internal_fe_subface_values_neighbor.reinit(cell_neighbor,
                                                 face_no_neighbor,
                                                 sub_face_no_neighbor);
      fe_face_values_neighbor = &internal_fe_subface_values_neighbor;
    }

  // Set up dof mapping and remove duplicates (for continuous elements).
  {
    // Get dof indices first:
    std::vector<types::global_dof_index> v(
      fe_face_values->get_fe().n_dofs_per_cell());
    cell->get_active_or_mg_dof_indices(v);
    std::vector<types::global_dof_index> v2(
      fe_face_values_neighbor->get_fe().n_dofs_per_cell());
    cell_neighbor->get_active_or_mg_dof_indices(v2);



    // Fill a map from the global dof index to the left and right
    // local index.
    std::map<types::global_dof_index, std::pair<unsigned int, unsigned int>>
                                          tempmap;
    std::pair<unsigned int, unsigned int> invalid_entry(
      numbers::invalid_unsigned_int, numbers::invalid_unsigned_int);

    for (unsigned int i = 0; i < v.size(); ++i)
      {
        // If not already existing, add an invalid entry:
        auto result = tempmap.insert(std::make_pair(v[i], invalid_entry));
        result.first->second.first = i;
      }

    for (unsigned int i = 0; i < v2.size(); ++i)
      {
        // If not already existing, add an invalid entry:
        auto result = tempmap.insert(std::make_pair(v2[i], invalid_entry));
        result.first->second.second = i;
      }

    // Transfer from the map to the sorted std::vectors.
    dofmap.resize(tempmap.size());
    interface_dof_indices.resize(tempmap.size());
    unsigned int idx = 0;
    for (auto &x : tempmap)
      {
        interface_dof_indices[idx] = x.first;
        dofmap[idx]                = {{x.second.first, x.second.second}};
        ++idx;
      }
  }
}



template <int dim, int spacedim>
template <class CellIteratorType>
void
FEInterfaceValues<dim, spacedim>::reinit(const CellIteratorType &cell,
                                         const unsigned int      face_no)
{
  internal_fe_face_values.reinit(cell, face_no);
  fe_face_values          = &internal_fe_face_values;
  fe_face_values_neighbor = nullptr;

  interface_dof_indices.resize(fe_face_values->get_fe().n_dofs_per_cell());
  cell->get_active_or_mg_dof_indices(interface_dof_indices);


  dofmap.resize(interface_dof_indices.size());

  for (unsigned int i = 0; i < interface_dof_indices.size(); ++i)
    {
      dofmap[i] = {{i, numbers::invalid_unsigned_int}};
    }
}



template <int dim, int spacedim>
inline double
FEInterfaceValues<dim, spacedim>::JxW(const unsigned int q) const
{
  Assert(fe_face_values != nullptr,
         ExcMessage("This call requires a call to reinit() first."));
  return fe_face_values->JxW(q);
}



template <int dim, int spacedim>
const std::vector<double> &
FEInterfaceValues<dim, spacedim>::get_JxW_values() const
{
  Assert(fe_face_values != nullptr,
         ExcMessage("This call requires a call to reinit() first."));
  return fe_face_values->get_JxW_values();
}



template <int dim, int spacedim>
const std::vector<Tensor<1, spacedim>> &
FEInterfaceValues<dim, spacedim>::get_normal_vectors() const
{
  Assert(fe_face_values != nullptr,
         ExcMessage("This call requires a call to reinit() first."));
  return fe_face_values->get_normal_vectors();
}



template <int dim, int spacedim>
const Quadrature<dim - 1> &
FEInterfaceValues<dim, spacedim>::get_quadrature() const
{
  return internal_fe_face_values.get_quadrature();
}



template <int dim, int spacedim>
const std::vector<Point<spacedim>> &
FEInterfaceValues<dim, spacedim>::get_quadrature_points() const
{
  Assert(fe_face_values != nullptr,
         ExcMessage("This call requires a call to reinit() first."));
  return fe_face_values->get_quadrature_points();
}



template <int dim, int spacedim>
UpdateFlags
FEInterfaceValues<dim, spacedim>::get_update_flags() const
{
  return internal_fe_face_values.get_update_flags();
}



template <int dim, int spacedim>
unsigned
FEInterfaceValues<dim, spacedim>::n_current_interface_dofs() const
{
  Assert(
    interface_dof_indices.size() > 0,
    ExcMessage(
      "n_current_interface_dofs() is only available after a call to reinit()."));
  return interface_dof_indices.size();
}



template <int dim, int spacedim>
bool
FEInterfaceValues<dim, spacedim>::at_boundary() const
{
  return fe_face_values_neighbor == nullptr;
}



template <int dim, int spacedim>
std::vector<types::global_dof_index>
FEInterfaceValues<dim, spacedim>::get_interface_dof_indices() const
{
  return interface_dof_indices;
}



template <int dim, int spacedim>
std::array<unsigned int, 2>
FEInterfaceValues<dim, spacedim>::interface_dof_to_dof_indices(
  const unsigned int interface_dof_index) const
{
  AssertIndexRange(interface_dof_index, n_current_interface_dofs());
  return dofmap[interface_dof_index];
}



template <int dim, int spacedim>
const FEFaceValuesBase<dim, spacedim> &
FEInterfaceValues<dim, spacedim>::get_fe_face_values(
  const unsigned int cell_index) const
{
  AssertIndexRange(cell_index, 2);
  Assert(
    cell_index == 0 || !at_boundary(),
    ExcMessage(
      "You are on a boundary, so you can only ask for the first FEFaceValues object."));

  return (cell_index == 0) ? *fe_face_values : *fe_face_values_neighbor;
}



template <int dim, int spacedim>
Tensor<1, spacedim>
FEInterfaceValues<dim, spacedim>::normal(const unsigned int q_point_index) const
{
  return fe_face_values->normal_vector(q_point_index);
}



template <int dim, int spacedim>
double
FEInterfaceValues<dim, spacedim>::shape_value(
  const bool         here_or_there,
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (here_or_there && dof_pair[0] != numbers::invalid_unsigned_int)
    return get_fe_face_values(0).shape_value_component(dof_pair[0],
                                                       q_point,
                                                       component);
  if (!here_or_there && dof_pair[1] != numbers::invalid_unsigned_int)
    return get_fe_face_values(1).shape_value_component(dof_pair[1],
                                                       q_point,
                                                       component);

  return 0.0;
}



template <int dim, int spacedim>
double
FEInterfaceValues<dim, spacedim>::jump(const unsigned int interface_dof_index,
                                       const unsigned int q_point,
                                       const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  double value = 0.0;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += get_fe_face_values(0).shape_value_component(dof_pair[0],
                                                         q_point,
                                                         component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value -= get_fe_face_values(1).shape_value_component(dof_pair[1],
                                                         q_point,
                                                         component);
  return value;
}



template <int dim, int spacedim>
double
FEInterfaceValues<dim, spacedim>::average(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return 1.0 * get_fe_face_values(0).shape_value_component(dof_pair[0],
                                                             q_point,
                                                             component);

  double value = 0.0;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(0).shape_value_component(dof_pair[0],
                                                               q_point,
                                                               component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(1).shape_value_component(dof_pair[1],
                                                               q_point,
                                                               component);

  return value;
}



template <int dim, int spacedim>
Tensor<1, dim>
FEInterfaceValues<dim, spacedim>::average_gradient(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_grad_component(dof_pair[0],
                                                      q_point,
                                                      component);

  Tensor<1, dim> value;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(0).shape_grad_component(dof_pair[0],
                                                              q_point,
                                                              component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(1).shape_grad_component(dof_pair[1],
                                                              q_point,
                                                              component);

  return value;
}



template <int dim, int spacedim>
Tensor<2, dim>
FEInterfaceValues<dim, spacedim>::average_hessian(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_hessian_component(dof_pair[0],
                                                         q_point,
                                                         component);

  Tensor<2, dim> value;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(0).shape_hessian_component(dof_pair[0],
                                                                 q_point,
                                                                 component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value += 0.5 * get_fe_face_values(1).shape_hessian_component(dof_pair[1],
                                                                 q_point,
                                                                 component);

  return value;
}



template <int dim, int spacedim>
Tensor<1, dim>
FEInterfaceValues<dim, spacedim>::jump_gradient(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_grad_component(dof_pair[0],
                                                      q_point,
                                                      component);

  Tensor<1, dim> value;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += 1.0 * get_fe_face_values(0).shape_grad_component(dof_pair[0],
                                                              q_point,
                                                              component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value += -1.0 * get_fe_face_values(1).shape_grad_component(dof_pair[1],
                                                               q_point,
                                                               component);

  return value;
}



template <int dim, int spacedim>
Tensor<2, dim>
FEInterfaceValues<dim, spacedim>::jump_hessian(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_hessian_component(dof_pair[0],
                                                         q_point,
                                                         component);

  Tensor<2, dim> value;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value += 1.0 * get_fe_face_values(0).shape_hessian_component(dof_pair[0],
                                                                 q_point,
                                                                 component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value += -1.0 * get_fe_face_values(1).shape_hessian_component(dof_pair[1],
                                                                  q_point,
                                                                  component);

  return value;
}


template <int dim, int spacedim>
Tensor<3, dim>
FEInterfaceValues<dim, spacedim>::jump_3rd_derivative(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto dof_pair = dofmap[interface_dof_index];

  if (at_boundary())
    return get_fe_face_values(0).shape_3rd_derivative_component(dof_pair[0],
                                                                q_point,
                                                                component);

  Tensor<3, dim> value;

  if (dof_pair[0] != numbers::invalid_unsigned_int)
    value +=
      1.0 * get_fe_face_values(0).shape_3rd_derivative_component(dof_pair[0],
                                                                 q_point,
                                                                 component);
  if (dof_pair[1] != numbers::invalid_unsigned_int)
    value +=
      -1.0 * get_fe_face_values(1).shape_3rd_derivative_component(dof_pair[1],
                                                                  q_point,
                                                                  component);

  return value;
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
