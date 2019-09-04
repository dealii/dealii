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
 * The class is made to be used inside MeshWorker::mesh_loop. It is intended
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
   * Construct an FEInterfaceValues from existing FEFaceValues objects.
   *
   * The arguments to this functions are used to construct the internal
   * FEFaceValues and FESubfaceValues objects for the two cells of the
   * interface. As the arguments are copied, there are no requirements
   * on the lifetime of the arguments.
   */
  FEInterfaceValues(const FEFaceValues<dim> &   fe_face_values,
                    const FESubfaceValues<dim> &fe_subface_values);

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
   * sub_face_no_neighbor to indicate that you want to work on the entire face,
   * not a sub-face.
   */
  template <class CellIteratorType>
  void
  reinit(const CellIteratorType &cell,
         const unsigned int &    face_no,
         const unsigned int &    sub_face_no,
         const CellIteratorType &cell_neighbor,
         const unsigned int &    face_no_neighbor,
         const unsigned int &    sub_face_no_neighbor);

  /**
   * Re-initialize this object to be used on a interface given by a single face
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
  reinit(const CellIteratorType &cell, const unsigned int &face_no);

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
   * Return the vector of JxW values for each quadrature point.
   */
  const std::vector<double> &
  get_JxW_values() const;

  /**
   * Return the normal vector of the interface in each quadrature point.
   *
   * The return value is identical to get_fe_face_values(0).get_normal_vectors()
   * and therefore, are outside normal vectors from the perspective of the
   * first cell of this interface.
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
   * Return if the current interface is a boundary face or an internal
   * face with two adjacent cells.
   *
   * See the corresponding reinit() functions for details.
   */
  bool
  at_boundary() const;

  /**
   * Return the set of joint DoF indices. This includes indices from both cells.
   */
  std::vector<types::global_dof_index>
  get_interface_dof_indices() const;

  /**
   * Convert a cell index @p cell_index (0 or 1) and dof
   * index @p fe_dof_index into an interface dof index.
   *
   * This is the inverse operation to interface_dof_to_cell_and_dof_index().
   */
  unsigned int
  cell_and_dof_to_interface_dof_index(const unsigned int cell_index,
                                      const unsigned int fe_dof_index) const;

  /**
   * Convert an interface dof index into the corresponding cell index
   * and dof index of that cell returned as a pair.
   *
   * This is the inverse operation to cell_and_dof_to_interface_dof_index().
   */
  std::pair<unsigned int, unsigned int>
  interface_dof_to_cell_and_dof_index(
    const unsigned int interface_dof_index) const;


  /**
   * Return the normal in a given quadrature point.
   *
   * The normal points in outwards direction as seen from the first cell of
   * this interface.
   */
  Tensor<1, spacedim>
  normal(const unsigned int q_point_index) const;

  /**
   * Return component @p component of the value of the shape function
   * with interface dof index @p interface_dof_index in
   * quadrature point @p q_point.
   *
   * The argument @p use_upstream_value selects between the upstream value and
   * the downstream value as defined by the direction of the normal vector
   * in this quadrature point. If @p use_upstream_value is true, the shape
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
  shape_value(const bool         use_upstream_value,
              const unsigned int interface_dof_index,
              const unsigned int q_point,
              const unsigned int component = 0) const;

  /**
   * Return the jump $[u]=u_{\text{cell1}} - u_{\text{cell2}}$ on the
   * interface
   * for the shape function @p interface_dof_index at the quadrature point
   * @p q_point of component @p component.
   */
  double
  jump(const unsigned int interface_dof_index,
       const unsigned int q_point,
       const unsigned int component = 0) const;

  /**
   * Return the average $\{u\}=frac{1}{2}u_{\text{cell}} +
   * \frac{1}{2}u_{\text{neighbor}}$ on the interface
   * for the shape function @p interface_dof_index at the quadrature point
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
  FEFaceValuesBase<dim> *fe_face_values = nullptr;

  /**
   * Pointer to internal_fe_face_values_neighbor,
   * internal_fe_subface_values_neighbor, or nullptr, respectively
   * as determined in reinit().
   */
  FEFaceValuesBase<dim> *fe_face_values_neighbor = nullptr;
};



#ifndef DOXYGEN

/*---------------------- Inline functions ---------------------*/



template <int dim, int spacedim>
FEInterfaceValues<dim, spacedim>::FEInterfaceValues(
  const FEFaceValues<dim> &   fe_face_values,
  const FESubfaceValues<dim> &fe_subface_values)
  : internal_fe_face_values(fe_face_values)
  , internal_fe_subface_values(fe_subface_values)
  , internal_fe_face_values_neighbor(fe_face_values)
  , internal_fe_subface_values_neighbor(fe_subface_values)
  , fe_face_values(nullptr)
  , fe_face_values_neighbor(nullptr)
{}



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
{}



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
{}



template <int dim, int spacedim>
template <class CellIteratorType>
void
FEInterfaceValues<dim, spacedim>::reinit(
  const CellIteratorType &cell,
  const unsigned int &    face_no,
  const unsigned int &    sub_face_no,
  const CellIteratorType &cell_neighbor,
  const unsigned int &    face_no_neighbor,
  const unsigned int &    sub_face_no_neighbor)
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
  std::vector<types::global_dof_index> v(
    fe_face_values->get_fe().n_dofs_per_cell);
  cell->get_dof_indices(v);
  std::vector<types::global_dof_index> v2(
    fe_face_values_neighbor->get_fe().n_dofs_per_cell);
  cell_neighbor->get_dof_indices(v2);

  interface_dof_indices = std::move(v);
  interface_dof_indices.insert(interface_dof_indices.end(),
                               v2.begin(),
                               v2.end());
}



template <int dim, int spacedim>
template <class CellIteratorType>
void
FEInterfaceValues<dim, spacedim>::reinit(const CellIteratorType &cell,
                                         const unsigned int &    face_no)
{
  internal_fe_face_values.reinit(cell, face_no);
  fe_face_values          = &internal_fe_face_values;
  fe_face_values_neighbor = nullptr;

  interface_dof_indices.resize(fe_face_values->get_fe().n_dofs_per_cell);
  cell->get_dof_indices(interface_dof_indices);
}



template <int dim, int spacedim>
const std::vector<double> &
FEInterfaceValues<dim, spacedim>::get_JxW_values() const
{
  return fe_face_values->get_JxW_values();
}



template <int dim, int spacedim>
const std::vector<Tensor<1, spacedim>> &
FEInterfaceValues<dim, spacedim>::get_normal_vectors() const
{
  return fe_face_values->get_normal_vectors();
}



template <int dim, int spacedim>
const Quadrature<dim - 1> &
FEInterfaceValues<dim, spacedim>::get_quadrature() const
{
  return fe_face_values->get_quadrature();
}



template <int dim, int spacedim>
const std::vector<Point<spacedim>> &
FEInterfaceValues<dim, spacedim>::get_quadrature_points() const
{
  return fe_face_values->get_quadrature_points();
}



template <int dim, int spacedim>
const UpdateFlags
FEInterfaceValues<dim, spacedim>::get_update_flags() const
{
  return fe_face_values->get_update_flags();
}



template <int dim, int spacedim>
unsigned
FEInterfaceValues<dim, spacedim>::n_interface_dofs() const
{
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
std::pair<unsigned int, unsigned int>
FEInterfaceValues<dim, spacedim>::interface_dof_to_cell_and_dof_index(
  const unsigned int interface_dof_index) const
{
  AssertIndexRange(interface_dof_index, n_interface_dofs());
  const unsigned int n_dofs_fe  = fe_face_values->get_fe().n_dofs_per_cell;
  const unsigned int cell_index = (interface_dof_index < n_dofs_fe) ? 0 : 1;
  const unsigned int dof_index =
    (cell_index == 1) ? (interface_dof_index - n_dofs_fe) : interface_dof_index;
  return std::make_pair(cell_index, dof_index);
}



template <int dim, int spacedim>
unsigned int
FEInterfaceValues<dim, spacedim>::cell_and_dof_to_interface_dof_index(
  const unsigned int cell_index,
  const unsigned int fe_dof_index) const
{
  Assert(cell_index <= 1, ExcMessage("cell_index should be 0 or 1"));
  if (at_boundary())
    Assert(cell_index == 0,
           ExcMessage("A boundary facet only has cell index 0."));

  const unsigned int n_dofs_fe = fe_face_values->get_fe().n_dofs_per_cell;
  if (cell_index == 0)
    {
      Assert(fe_dof_index < n_dofs_fe, ExcMessage("invalid fe_dof_index"));
      return fe_dof_index;
    }
  else if (cell_index == 1)
    {
      Assert(fe_dof_index < fe_face_values_neighbor->get_fe().n_dofs_per_cell,
             ExcMessage("invalid fe_dof_index"));
      return fe_dof_index + n_dofs_fe;
    }

  // return something invalid:
  return numbers::invalid_unsigned_int;
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
  const bool         current_cell,
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto cell_and_dof_idx =
    interface_dof_to_cell_and_dof_index(interface_dof_index);

  if (current_cell && cell_and_dof_idx.second == 0)
    return get_fe_face_values(0).shape_value_component(cell_and_dof_idx.first,
                                                       q_point,
                                                       component);
  if (!current_cell && cell_and_dof_idx.second == 1)
    return get_fe_face_values(1).shape_value_component(cell_and_dof_idx.first,
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
  const auto cell_and_dof_idx =
    interface_dof_to_cell_and_dof_index(interface_dof_index);

  if (cell_and_dof_idx.first == 0)
    return get_fe_face_values(0).shape_value_component(cell_and_dof_idx.second,
                                                       q_point,
                                                       component);
  else
    {
      if (at_boundary())
        return 0.0;
      else
        return -get_fe_face_values(1).shape_value_component(
          cell_and_dof_idx.second, q_point, component);
    }
}



template <int dim, int spacedim>
double
FEInterfaceValues<dim, spacedim>::average(
  const unsigned int interface_dof_index,
  const unsigned int q_point,
  const unsigned int component) const
{
  const auto cell_and_dof_idx =
    interface_dof_to_cell_and_dof_index(interface_dof_index);

  if (cell_and_dof_idx.first == 0)
    return 0.5 * get_fe_face_values(0).shape_value_component(
                   cell_and_dof_idx.second, q_point, component);
  else
    {
      if (at_boundary())
        return 0.0;
      else
        return 0.5 * get_fe_face_values(1).shape_value_component(
                       cell_and_dof_idx.second, q_point, component);
    }
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
