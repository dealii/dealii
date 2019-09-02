// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

DEAL_II_NAMESPACE_OPEN

/**
 * FEInterfaceValues is a data structure to access and assemble Finite Element
 * data on interfaces between two cells of a mesh.
 *
 * It provides a way to access average and jump terms used in
 * Discontinuous Galerkin methods on faces between two neighboring cells. In
 * the literature these are called "inner interfaces" or "facets".
 *
 * Internally, this class provides an abstraction for two FEFaceValues
 * objects. The class introduces a new "interface dof index" that walks over
 * the union of the dof indices of the two FEFaceValues objects. Helper
 * functions allow translating between the new "interface dof index" and the
 * corresponding "fe index" and "dof index" of the FE.
 *
 * The class is made to be used inside MeshWorker::mesh_loop. It is intended
 * to be a low level replacement for MeshWorker and LocalIntegrators and a
 * higher level abstraction compared to assembling face terms manually.
 *
 *
 * @author Timo Heister, 2019.
 */
template <int dim, int spacedim = dim>
class FEInterfaceValues
{
public:
  /**
   * Construct an FEInterfaceValues from existing FEFaceValues objects.
   */
  FEInterfaceValues(const FEFaceValues<dim> &   fe,
                    const FESubfaceValues<dim> &fe_sub,
                    const FEFaceValues<dim> &   fe_neighbor,
                    const FESubfaceValues<dim> &fe_sub_neighbor);

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
   * of two neighboring cells.
   *
   * The arguments (including their order) is identical to the @p face_worker arguments
   * in MeshWorker::mesh_loop.
   *
   * Use numbers::invalid_unsigned_int for @p sub_face_no or @p
   * sub_face_no_neighbor to indicate that no subface is in use.
   */
  template <class Iterator>
  void
  reinit(const Iterator &    cell,
         const unsigned int &face_no,
         const unsigned int &sub_face_no,
         const Iterator &    cell_neighbor,
         const unsigned int &face_no_neighbor,
         const unsigned int &sub_face_no_neighbor);

  /**
   * Re-initialize this object to be used on a interface given by a single face
   * of a cell. This is useful to use FEInterfaceValues to work on boundaries of
   * the domain.
   *
   * As a consequence, queries like jump() will assume a value of zero for the
   * values on the "other" side.
   */
  template <class Iterator>
  void
  reinit(const Iterator &cell, const unsigned int &face_no);

  /**
   * Return the vector of JxW values for each quadrature point.
   */
  const std::vector<double> &
  get_JxW_values() const;

  /**
   * Return the normal vector of the interface in each quadrature point.
   */
  const std::vector<Tensor<1, spacedim>> &
  get_normal_vectors() const;

  /**
   * Return a reference to the quadrature object in use.
   */
  const Quadrature<dim - 1> &
  get_quadrature() const;

  /**
   * Return a reference to the quadrature points.
   */
  const std::vector<Point<spacedim>> &
  get_quadrature_points() const;

  /**
   * Return the update flags set.
   */
  const UpdateFlags
  get_update_flags() const;

  /**
   * Return the number of DoFs on the current interface. This is the sum of the
   * number of DoFs on the two adjacent cells.
   */
  unsigned
  n_interface_dofs() const;

  /**
   * Return if the current interface is an internal
   * face with two adjacent cells or a boundary face.
   */
  bool
  at_boundary() const;

  /**
   * Return the set of joint DoF indices. This includes indices from both cells.
   */
  std::vector<types::global_dof_index>
  get_interface_dof_indices() const;

  /**
   * Translate a local "interface dof index" to a local "dof index" of the
   * underlying FEFaceValues.
   */
  unsigned int
  interface_dof_index_to_fe_dof_index(unsigned int interface_dof_index) const;

  /**
   * For a given "interface dof index" return whether it belongs to the cell (0)
   * or to the neighbor (1).
   */
  unsigned int
  interface_dof_index_to_fe_index(const unsigned int interface_dof_index) const;

  /**
   * Convert a FiniteElement index @p fe_index (0=cell, 1=neighbor) and dof
   * index @p face_dof_index into an interface DoF index.
   *
   * The inverse of this operation is done with
   * interface_dof_index_to_fe_index() and
   * interface_dof_index_to_fe_dof_index().
   */
  unsigned int
  interface_dof_index(const unsigned int fe_index,
                      const unsigned int face_dof_index) const;

  /**
   * Return the FEFaceValues object of the cell or the neighboring cell
   *
   * If @p cell_or_neighbor is 0, return the cell FEFaceValues object , if it is 1, return the
   * neighboring cell FEFaceValues object.
   *
   * @note The argument @p cell_or_neighbor is returned by
   * interface_dof_index_to_fe_index().
   */
  const FEFaceValuesBase<dim, spacedim> &
  get_fe_values(const unsigned int cell_or_neighbor) const;

  /**
   * Return the FEFaceValue object of the current cell
   */
  const FEFaceValuesBase<dim, spacedim> &
  get_fe_values() const;

  /**
   * Return the FEFaceValue object of the current neighboring cell
   */
  const FEFaceValuesBase<dim, spacedim> &
  get_fe_values_neighbor() const;

  /**
   * Return the normal in a given quadrature point.
   */
  Tensor<1, spacedim>
  normal(const unsigned int q_point_index) const;

  /**
   * Return the current value (if @p current_cell is true) or the neighboring value
   * (if @p current_cell is false) on the interface for the shape function @p
   * interface_dof_index in the quadrature point @p q_point of component @p
   * component.
   *
   * @note This function is typically used to pick the upstream or downstream
   * value based on a direction. This can be achieved by using
   * <code>(direction * normal)>0</code> as the first argument of this
   * function.
   */
  double
  shape_value(const bool         current_cell,
              const unsigned int interface_dof_index,
              const unsigned int q_point,
              const unsigned int component = 0) const;

  /**
   * Return the jump $[u]=u_{\text{cell}} - u_{\text{neighbor}}$ on the
   * interface
   * for the shape function @p interface_dof_index in the quadrature point
   * @p q_point of component @p component.
   */
  double
  jump(const unsigned int interface_dof_index,
       const unsigned int q_point,
       const unsigned int component = 0) const;

  /**
   * Return the average $\{u\}=frac{1}{2}u_{\text{cell}} +
   * \frac{1}{2}u_{\text{neighbor}}$ on the interface
   * for the shape function @p interface_dof_index in the quadrature point
   * @p q_point of component @p component.
   */
  double
  average(const unsigned int interface_dof_index,
          const unsigned int q_point,
          const unsigned int component = 0) const;

private:
  /**
   * The list of DoF indices for the current interface, filled in reinit().
   */
  std::vector<types::global_dof_index> interface_dof_indices;

  /**
   * The number of DoFs belonging to the current cell
   */
  unsigned int n_dofs_fe;

  /**
   * The number of DoFs belonging to the neighboring cell
   */
  unsigned int n_dofs_fe_neighbor;

  /**
   * Pointer to internal_fe_face_values or internal_fe_subface_values,
   * respectively.
   */
  FEFaceValuesBase<dim> *fe_face_values = nullptr;

  /**
   * Pointer to internal_fe_face_values_neighbor,
   * internal_fe_subface_values_neighbor, or nullptr, respectively.
   */
  FEFaceValuesBase<dim> *fe_face_values_neighbor = nullptr;

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
};



#ifndef DOXYGEN

/*---------------------- Inline functions ---------------------*/



template <int dim, int spacedim>
FEInterfaceValues<dim, spacedim>::FEInterfaceValues(
  const FEFaceValues<dim> &   fe,
  const FESubfaceValues<dim> &fe_sub,
  const FEFaceValues<dim> &   fe_neighbor,
  const FESubfaceValues<dim> &fe_sub_neighbor)
  : internal_fe_face_values(fe)
  , internal_fe_subface_values(fe_sub)
  , internal_fe_face_values_neighbor(fe_neighbor)
  , internal_fe_subface_values_neighbor(fe_sub_neighbor)
  , fe_face_values(nullptr)
  , fe_face_values_neighbor(nullptr)
{
  n_dofs_fe          = fe.get_fe().n_dofs_per_cell();
  n_dofs_fe_neighbor = fe_neighbor.get_fe().n_dofs_per_cell();
}



template <int dim, int spacedim>
FEInterfaceValues<dim, spacedim>::FEInterfaceValues(
  const Mapping<dim, spacedim> &      mapping,
  const FiniteElement<dim, spacedim> &fe,
  const Quadrature<dim - 1> &         quadrature,
  const UpdateFlags                   update_flags)
  : internal_fe_face_values(mapping, fe, quadrature, update_flags)
  , internal_fe_subface_values(mapping, fe, quadrature, update_flags)
  , internal_fe_face_values_neighbor(mapping, fe, quadrature, update_flags)
  , internal_fe_subface_values_neighbor(mapping, fe, quadrature, update_flags)
  , fe_face_values(nullptr)
  , fe_face_values_neighbor(nullptr)
{
  n_dofs_fe          = fe.n_dofs_per_cell();
  n_dofs_fe_neighbor = fe.n_dofs_per_cell();
}



template <int dim, int spacedim>
FEInterfaceValues<dim, spacedim>::FEInterfaceValues(
  const FiniteElement<dim, spacedim> &fe,
  const Quadrature<dim - 1> &         quadrature,
  const UpdateFlags                   update_flags)
  : internal_fe_face_values(StaticMappingQ1<dim, spacedim>::mapping,
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
{
  n_dofs_fe          = fe.n_dofs_per_cell();
  n_dofs_fe_neighbor = fe.n_dofs_per_cell();
}



template <int dim, int spacedim>
template <class Iterator>
void
FEInterfaceValues<dim, spacedim>::reinit(
  const Iterator &    cell,
  const unsigned int &face_no,
  const unsigned int &sub_face_no,
  const Iterator &    cell_neighbor,
  const unsigned int &face_no_neighbor,
  const unsigned int &sub_face_no_neighbor)
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

  // Fill the DoF indices:
  std::vector<types::global_dof_index> v(n_dofs_fe);
  cell->get_dof_indices(v);
  std::vector<types::global_dof_index> v2(n_dofs_fe_neighbor);
  cell_neighbor->get_dof_indices(v2);

  interface_dof_indices.clear();
  interface_dof_indices.insert(interface_dof_indices.end(), v.begin(), v.end());
  interface_dof_indices.insert(interface_dof_indices.end(),
                               v2.begin(),
                               v2.end());
}



template <int dim, int spacedim>
template <class Iterator>
void
FEInterfaceValues<dim, spacedim>::reinit(const Iterator &    cell,
                                         const unsigned int &face_no)
{
  internal_fe_face_values.reinit(cell, face_no);
  fe_face_values          = &internal_fe_face_values;
  fe_face_values_neighbor = nullptr;

  interface_dof_indices.resize(n_dofs_fe);
  cell->get_dof_indices(interface_dof_indices);
}



template <int dim, int spacedim>
const std::vector<double> &
FEInterfaceValues<dim, spacedim>::get_JxW_values() const
{
  return internal_fe_face_values.get_JxW_values();
}



template <int dim, int spacedim>
const std::vector<Tensor<1, spacedim>> &
FEInterfaceValues<dim, spacedim>::get_normal_vectors() const
{
  return internal_fe_face_values.get_normal_vectors();
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
  return internal_fe_face_values.get_quadrature_points();
}



template <int dim, int spacedim>
const UpdateFlags
FEInterfaceValues<dim, spacedim>::get_update_flags() const
{
  return internal_fe_face_values.get_update_flags();
}



template <int dim, int spacedim>
unsigned
FEInterfaceValues<dim, spacedim>::n_interface_dofs() const
{
  return n_dofs_fe + (at_boundary() ? 0 : n_dofs_fe_neighbor);
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
  Assert(interface_dof_indices.size() == n_interface_dofs(),
         ExcInternalError());
  return interface_dof_indices;
}



template <int dim, int spacedim>
unsigned int
FEInterfaceValues<dim, spacedim>::interface_dof_index_to_fe_dof_index(
  unsigned int interface_dof_index) const
{
  AssertIndexRange(interface_dof_index, n_interface_dofs());
  return (interface_dof_index >= n_dofs_fe) ?
           (interface_dof_index - n_dofs_fe) :
           interface_dof_index;
}



template <int dim, int spacedim>
unsigned int
FEInterfaceValues<dim, spacedim>::interface_dof_index_to_fe_index(
  const unsigned int interface_dof_index) const
{
  AssertIndexRange(interface_dof_index, n_interface_dofs());
  return (interface_dof_index < n_dofs_fe) ? 0 : 1;
}



template <int dim, int spacedim>
unsigned int
FEInterfaceValues<dim, spacedim>::interface_dof_index(
  const unsigned int fe_index,
  const unsigned int face_dof_index) const
{
  Assert(fe_index <= 1, ExcMessage("fe_index should be 0 or 1"));
  if (at_boundary())
    Assert(fe_index == 0, ExcMessage("A boundary facet only has FE index 0."));

  if (fe_index == 0)
    {
      Assert(face_dof_index < n_dofs_fe, ExcMessage("invalid face_dof_idx"));
      return face_dof_index;
    }
  else if (fe_index == 1)
    {
      Assert(face_dof_index < n_dofs_fe_neighbor,
             ExcMessage("invalid face_dof_idx"));
      return face_dof_index + n_dofs_fe;
    }

  // return something invalid:
  return numbers::invalid_unsigned_int;
}



template <int dim, int spacedim>
const FEFaceValuesBase<dim, spacedim> &
FEInterfaceValues<dim, spacedim>::get_fe_values(
  const unsigned int cell_or_neighbor) const
{
  AssertIndexRange(cell_or_neighbor, 2);
  return (cell_or_neighbor == 0) ? get_fe_values() : get_fe_values_neighbor();
}



template <int dim, int spacedim>
const FEFaceValuesBase<dim, spacedim> &
FEInterfaceValues<dim, spacedim>::get_fe_values() const
{
  return *fe_face_values;
}



template <int dim, int spacedim>
const FEFaceValuesBase<dim, spacedim> &
FEInterfaceValues<dim, spacedim>::get_fe_values_neighbor() const
{
  Assert(!at_boundary(), ExcMessage("Not possible for boundary Facet."));
  return *fe_face_values_neighbor;
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
  const bool         current_cell,
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const unsigned int shape_fct =
    interface_dof_index_to_fe_dof_index(interface_dof_index);
  const unsigned int fe_idx =
    interface_dof_index_to_fe_index(interface_dof_index);

  if (current_cell && fe_idx == 0)
    return get_fe_values().shape_value_component(shape_fct, q_point, component);
  if (!current_cell && fe_idx == 1)
    return get_fe_values_neighbor().shape_value_component(shape_fct,
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
  const unsigned int shape_fct =
    interface_dof_index_to_fe_dof_index(interface_dof_index);
  const unsigned int fe_idx =
    interface_dof_index_to_fe_index(interface_dof_index);

  if (fe_idx == 0)
    return get_fe_values().shape_value_component(shape_fct, q_point, component);
  else
    {
      if (at_boundary())
        return 0.0;
      else
        return -get_fe_values_neighbor().shape_value_component(shape_fct,
                                                               q_point,
                                                               component);
    }
}



template <int dim, int spacedim>
double
FEInterfaceValues<dim, spacedim>::average(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const unsigned int shape_fct =
    interface_dof_index_to_fe_dof_index(interface_dof_index);
  const unsigned int fe_idx =
    interface_dof_index_to_fe_index(interface_dof_index);

  if (fe_idx == 0)
    return 0.5 *
           get_fe_values().shape_value_component(shape_fct, q_point, component);
  else
    {
      if (at_boundary())
        return 0.0;
      else
        return 0.5 * get_fe_values_neighbor().shape_value_component(shape_fct,
                                                                    q_point,
                                                                    component);
    }
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
