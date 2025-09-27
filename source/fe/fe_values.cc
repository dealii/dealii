// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/array_view.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/thread_management.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <boost/container/small_vector.hpp>

#include <iomanip>
#include <memory>
#include <type_traits>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  template <int dim, int spacedim>
  inline std::vector<unsigned int>
  make_shape_function_to_row_table(const FiniteElement<dim, spacedim> &fe)
  {
    std::vector<unsigned int> shape_function_to_row_table(
      fe.n_dofs_per_cell() * fe.n_components(), numbers::invalid_unsigned_int);
    unsigned int row = 0;
    for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
      {
        // loop over all components that are nonzero for this particular
        // shape function. if a component is zero then we leave the
        // value in the table unchanged (at the invalid value)
        // otherwise it is mapped to the next free entry
        unsigned int nth_nonzero_component = 0;
        for (unsigned int c = 0; c < fe.n_components(); ++c)
          if (fe.get_nonzero_components(i)[c] == true)
            {
              shape_function_to_row_table[i * fe.n_components() + c] =
                row + nth_nonzero_component;
              ++nth_nonzero_component;
            }
        row += fe.n_nonzero_components(i);
      }

    return shape_function_to_row_table;
  }
} // namespace internal



namespace internal
{
  namespace FEValuesImplementation
  {
    template <int dim, int spacedim>
    void
    FiniteElementRelatedData<dim, spacedim>::initialize(
      const unsigned int                  n_quadrature_points,
      const FiniteElement<dim, spacedim> &fe,
      const UpdateFlags                   flags)
    {
      // initialize the table mapping from shape function number to
      // the rows in the tables storing the data by shape function and
      // nonzero component
      this->shape_function_to_row_table =
        dealii::internal::make_shape_function_to_row_table(fe);

      // count the total number of non-zero components accumulated
      // over all shape functions
      unsigned int n_nonzero_shape_components = 0;
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        n_nonzero_shape_components += fe.n_nonzero_components(i);
      Assert(n_nonzero_shape_components >= fe.n_dofs_per_cell(),
             ExcInternalError());

      // with the number of rows now known, initialize those fields
      // that we will need to their correct size
      if (flags & update_values)
        {
          this->shape_values.reinit(n_nonzero_shape_components,
                                    n_quadrature_points);
          this->shape_values.fill(numbers::signaling_nan<double>());
        }

      if (flags & update_gradients)
        {
          this->shape_gradients.reinit(n_nonzero_shape_components,
                                       n_quadrature_points);
          this->shape_gradients.fill(
            numbers::signaling_nan<Tensor<1, spacedim>>());
        }

      if (flags & update_hessians)
        {
          this->shape_hessians.reinit(n_nonzero_shape_components,
                                      n_quadrature_points);
          this->shape_hessians.fill(
            numbers::signaling_nan<Tensor<2, spacedim>>());
        }

      if (flags & update_3rd_derivatives)
        {
          this->shape_3rd_derivatives.reinit(n_nonzero_shape_components,
                                             n_quadrature_points);
          this->shape_3rd_derivatives.fill(
            numbers::signaling_nan<Tensor<3, spacedim>>());
        }
    }



    template <int dim, int spacedim>
    std::size_t
    FiniteElementRelatedData<dim, spacedim>::memory_consumption() const
    {
      return (
        MemoryConsumption::memory_consumption(shape_values) +
        MemoryConsumption::memory_consumption(shape_gradients) +
        MemoryConsumption::memory_consumption(shape_hessians) +
        MemoryConsumption::memory_consumption(shape_3rd_derivatives) +
        MemoryConsumption::memory_consumption(shape_function_to_row_table));
    }
  } // namespace FEValuesImplementation
} // namespace internal

/*------------------------------- FEValues -------------------------------*/
#ifndef DOXYGEN

template <int dim, int spacedim>
const unsigned int FEValues<dim, spacedim>::integral_dimension;



template <int dim, int spacedim>
FEValues<dim, spacedim>::FEValues(const Mapping<dim, spacedim>       &mapping,
                                  const FiniteElement<dim, spacedim> &fe,
                                  const Quadrature<dim>              &q,
                                  const UpdateFlags update_flags)
  : FEValuesBase<dim, spacedim>(q.size(),
                                fe.n_dofs_per_cell(),
                                update_default,
                                mapping,
                                fe)
  , quadrature(q)
{
  initialize(update_flags);
}



template <int dim, int spacedim>
FEValues<dim, spacedim>::FEValues(const Mapping<dim, spacedim>       &mapping,
                                  const FiniteElement<dim, spacedim> &fe,
                                  const hp::QCollection<dim>         &q,
                                  const UpdateFlags update_flags)
  : FEValues(mapping, fe, q[0], update_flags)
{
  AssertDimension(q.size(), 1);
}



template <int dim, int spacedim>
FEValues<dim, spacedim>::FEValues(const FiniteElement<dim, spacedim> &fe,
                                  const Quadrature<dim>              &q,
                                  const UpdateFlags update_flags)
  : FEValuesBase<dim, spacedim>(
      q.size(),
      fe.n_dofs_per_cell(),
      update_default,
      fe.reference_cell().template get_default_linear_mapping<dim, spacedim>(),
      fe)
  , quadrature(q)
{
  initialize(update_flags);
}



template <int dim, int spacedim>
FEValues<dim, spacedim>::FEValues(const FiniteElement<dim, spacedim> &fe,
                                  const hp::QCollection<dim>         &q,
                                  const UpdateFlags update_flags)
  : FEValues(fe, q[0], update_flags)
{
  AssertDimension(q.size(), 1);
}



template <int dim, int spacedim>
void
FEValues<dim, spacedim>::initialize(const UpdateFlags update_flags)
{
  // You can compute normal vectors to the cells only in the
  // codimension one case.
  if (dim != spacedim - 1)
    Assert((update_flags & update_normal_vectors) == false,
           ExcMessage("You can only pass the 'update_normal_vectors' "
                      "flag to FEFaceValues or FESubfaceValues objects, "
                      "but not to an FEValues object unless the "
                      "triangulation it refers to is embedded in a higher "
                      "dimensional space."));

  const UpdateFlags flags = this->compute_update_flags(update_flags);

  // initialize the base classes
  if (flags & update_mapping)
    this->mapping_output.initialize(this->max_n_quadrature_points, flags);
  this->finite_element_output.initialize(this->max_n_quadrature_points,
                                         *this->fe,
                                         flags);

  // then get objects into which the FE and the Mapping can store
  // intermediate data used across calls to reinit. we can do this in parallel
  Threads::Task<
    std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>>
    fe_get_data = Threads::new_task([&]() {
      return this->fe->get_data(flags,
                                *this->mapping,
                                quadrature,
                                this->finite_element_output);
    });

  Threads::Task<
    std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>>
    mapping_get_data;
  if (flags & update_mapping)
    mapping_get_data = Threads::new_task(
      [&]() { return this->mapping->get_data(flags, quadrature); });

  this->update_flags = flags;

  // then collect answers from the two task above
  this->fe_data = std::move(fe_get_data.return_value());
  if (flags & update_mapping)
    this->mapping_data = std::move(mapping_get_data.return_value());
  else
    this->mapping_data =
      std::make_unique<typename Mapping<dim, spacedim>::InternalDataBase>();
}



template <int dim, int spacedim>
void
FEValues<dim, spacedim>::reinit(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell)
{
  // Check that mapping and reference cell type are compatible:
  Assert(this->get_mapping().is_compatible_with(cell->reference_cell()),
         ExcMessage(
           "You are trying to call FEValues::reinit() with a cell of type " +
           cell->reference_cell().to_string() +
           " with a Mapping that is not compatible with it."));

  // no FE in this cell, so no assertion
  // necessary here
  this->maybe_invalidate_previous_present_cell(cell);
  this->check_cell_similarity(cell);

  this->present_cell = {cell};

  // this was the part of the work that is dependent on the actual
  // data type of the iterator. now pass on to the function doing
  // the real work.
  do_reinit();
}



template <int dim, int spacedim>
template <bool lda>
void
FEValues<dim, spacedim>::reinit(
  const TriaIterator<DoFCellAccessor<dim, spacedim, lda>> &cell)
{
  // assert that the finite elements passed to the constructor and
  // used by the DoFHandler used by this cell, are the same
  Assert(static_cast<const FiniteElementData<dim> &>(*this->fe) ==
           static_cast<const FiniteElementData<dim> &>(cell->get_fe()),
         (typename FEValuesBase<dim, spacedim>::ExcFEDontMatch()));

  // Check that mapping and reference cell type are compatible:
  Assert(this->get_mapping().is_compatible_with(cell->reference_cell()),
         ExcMessage(
           "You are trying to call FEValues::reinit() with a cell of type " +
           cell->reference_cell().to_string() +
           " with a Mapping that is not compatible with it."));

  this->maybe_invalidate_previous_present_cell(cell);
  this->check_cell_similarity(cell);

  this->present_cell = {cell};

  // this was the part of the work that is dependent on the actual
  // data type of the iterator. now pass on to the function doing
  // the real work.
  do_reinit();
}



template <int dim, int spacedim>
void
FEValues<dim, spacedim>::do_reinit()
{
  // first call the mapping and let it generate the data
  // specific to the mapping. also let it inspect the
  // cell similarity flag and, if necessary, update
  // it
  if (this->update_flags & update_mapping)
    {
      this->cell_similarity =
        this->get_mapping().fill_fe_values(this->present_cell,
                                           this->cell_similarity,
                                           quadrature,
                                           *this->mapping_data,
                                           this->mapping_output);
    }

  // then call the finite element and, with the data
  // already filled by the mapping, let it compute the
  // data for the mapped shape function values, gradients,
  // etc.
  this->get_fe().fill_fe_values(this->present_cell,
                                this->cell_similarity,
                                this->quadrature,
                                this->get_mapping(),
                                *this->mapping_data,
                                this->mapping_output,
                                *this->fe_data,
                                this->finite_element_output);
}



template <int dim, int spacedim>
std::size_t
FEValues<dim, spacedim>::memory_consumption() const
{
  return (FEValuesBase<dim, spacedim>::memory_consumption() +
          MemoryConsumption::memory_consumption(quadrature));
}

#endif
/*------------------------------- FEFaceValuesBase --------------------------*/
#ifndef DOXYGEN

template <int dim, int spacedim>
FEFaceValuesBase<dim, spacedim>::FEFaceValuesBase(
  const unsigned int                  dofs_per_cell,
  const UpdateFlags                   flags,
  const Mapping<dim, spacedim>       &mapping,
  const FiniteElement<dim, spacedim> &fe,
  const Quadrature<dim - 1>          &quadrature)
  : FEFaceValuesBase<dim, spacedim>(dofs_per_cell,
                                    flags,
                                    mapping,
                                    fe,
                                    hp::QCollection<dim - 1>(quadrature))
{}



template <int dim, int spacedim>
FEFaceValuesBase<dim, spacedim>::FEFaceValuesBase(
  const unsigned int dofs_per_cell,
  const UpdateFlags,
  const Mapping<dim, spacedim>       &mapping,
  const FiniteElement<dim, spacedim> &fe,
  const hp::QCollection<dim - 1>     &quadrature)
  : FEValuesBase<dim, spacedim>(quadrature.max_n_quadrature_points(),
                                dofs_per_cell,
                                update_default,
                                mapping,
                                fe)
  , present_face_index(numbers::invalid_unsigned_int)
  , quadrature(quadrature)
{
  Assert(quadrature.size() == 1 ||
           quadrature.size() == fe.reference_cell().n_faces(),
         ExcInternalError());
}



template <int dim, int spacedim>
const std::vector<Tensor<1, spacedim>> &
FEFaceValuesBase<dim, spacedim>::get_boundary_forms() const
{
  Assert(this->update_flags & update_boundary_forms,
         (typename FEValuesBase<dim, spacedim>::ExcAccessToUninitializedField(
           "update_boundary_forms")));
  return this->mapping_output.boundary_forms;
}



template <int dim, int spacedim>
std::size_t
FEFaceValuesBase<dim, spacedim>::memory_consumption() const
{
  return (FEValuesBase<dim, spacedim>::memory_consumption() +
          MemoryConsumption::memory_consumption(quadrature));
}
#endif

/*------------------------------- FEFaceValues -------------------------------*/

#ifndef DOXYGEN

template <int dim, int spacedim>
const unsigned int FEFaceValues<dim, spacedim>::dimension;



template <int dim, int spacedim>
const unsigned int FEFaceValues<dim, spacedim>::integral_dimension;



template <int dim, int spacedim>
FEFaceValues<dim, spacedim>::FEFaceValues(
  const Mapping<dim, spacedim>       &mapping,
  const FiniteElement<dim, spacedim> &fe,
  const Quadrature<dim - 1>          &quadrature,
  const UpdateFlags                   update_flags)
  : FEFaceValues<dim, spacedim>(mapping,
                                fe,
                                hp::QCollection<dim - 1>(quadrature),
                                update_flags)
{}



template <int dim, int spacedim>
FEFaceValues<dim, spacedim>::FEFaceValues(
  const Mapping<dim, spacedim>       &mapping,
  const FiniteElement<dim, spacedim> &fe,
  const hp::QCollection<dim - 1>     &quadrature,
  const UpdateFlags                   update_flags)
  : FEFaceValuesBase<dim, spacedim>(fe.n_dofs_per_cell(),
                                    update_flags,
                                    mapping,
                                    fe,
                                    quadrature)
{
  initialize(update_flags);
}



template <int dim, int spacedim>
FEFaceValues<dim, spacedim>::FEFaceValues(
  const FiniteElement<dim, spacedim> &fe,
  const Quadrature<dim - 1>          &quadrature,
  const UpdateFlags                   update_flags)
  : FEFaceValues<dim, spacedim>(fe,
                                hp::QCollection<dim - 1>(quadrature),
                                update_flags)
{}



template <int dim, int spacedim>
FEFaceValues<dim, spacedim>::FEFaceValues(
  const FiniteElement<dim, spacedim> &fe,
  const hp::QCollection<dim - 1>     &quadrature,
  const UpdateFlags                   update_flags)
  : FEFaceValuesBase<dim, spacedim>(
      fe.n_dofs_per_cell(),
      update_flags,
      fe.reference_cell().template get_default_linear_mapping<dim, spacedim>(),
      fe,
      quadrature)
{
  initialize(update_flags);
}



template <int dim, int spacedim>
void
FEFaceValues<dim, spacedim>::initialize(const UpdateFlags update_flags)
{
  const UpdateFlags flags = this->compute_update_flags(update_flags);

  // initialize the base classes
  if (flags & update_mapping)
    this->mapping_output.initialize(this->max_n_quadrature_points, flags);
  this->finite_element_output.initialize(this->max_n_quadrature_points,
                                         *this->fe,
                                         flags);

  // then get objects into which the FE and the Mapping can store
  // intermediate data used across calls to reinit. this can be done in parallel

  std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase> (
    FiniteElement<dim, spacedim>::*finite_element_get_face_data)(
    const UpdateFlags,
    const Mapping<dim, spacedim> &,
    const hp::QCollection<dim - 1> &,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &) const = &FiniteElement<dim, spacedim>::get_face_data;

  std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase> (
    Mapping<dim, spacedim>::*mapping_get_face_data)(
    const UpdateFlags, const hp::QCollection<dim - 1> &) const =
    &Mapping<dim, spacedim>::get_face_data;


  Threads::Task<
    std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>>
    fe_get_data = Threads::new_task(finite_element_get_face_data,
                                    *this->fe,
                                    flags,
                                    *this->mapping,
                                    this->quadrature,
                                    this->finite_element_output);
  Threads::Task<
    std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>>
    mapping_get_data;
  if (flags & update_mapping)
    mapping_get_data = Threads::new_task(mapping_get_face_data,
                                         *this->mapping,
                                         flags,
                                         this->quadrature);

  this->update_flags = flags;

  // then collect answers from the two task above
  this->fe_data = std::move(fe_get_data.return_value());
  if (flags & update_mapping)
    this->mapping_data = std::move(mapping_get_data.return_value());
  else
    this->mapping_data =
      std::make_unique<typename Mapping<dim, spacedim>::InternalDataBase>();
}



template <int dim, int spacedim>
template <bool lda>
void
FEFaceValues<dim, spacedim>::reinit(
  const TriaIterator<DoFCellAccessor<dim, spacedim, lda>> &cell,
  const unsigned int                                       face_no)
{
  // assert that the finite elements passed to the constructor and
  // used by the DoFHandler used by this cell, are the same
  Assert(static_cast<const FiniteElementData<dim> &>(*this->fe) ==
           static_cast<const FiniteElementData<dim> &>(
             cell->get_dof_handler().get_fe(cell->active_fe_index())),
         (typename FEValuesBase<dim, spacedim>::ExcFEDontMatch()));

  AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);

  this->maybe_invalidate_previous_present_cell(cell);
  this->present_cell = {cell};

  // this was the part of the work that is dependent on the actual
  // data type of the iterator. now pass on to the function doing
  // the real work.
  do_reinit(face_no);
}



template <int dim, int spacedim>
template <bool lda>
void
FEFaceValues<dim, spacedim>::reinit(
  const TriaIterator<DoFCellAccessor<dim, spacedim, lda>>    &cell,
  const typename Triangulation<dim, spacedim>::face_iterator &face)
{
  const auto face_n = cell->face_iterator_to_index(face);
  reinit(cell, face_n);
}



template <int dim, int spacedim>
void
FEFaceValues<dim, spacedim>::reinit(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no)
{
  AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);

  this->maybe_invalidate_previous_present_cell(cell);
  this->present_cell = {cell};

  // this was the part of the work that is dependent on the actual
  // data type of the iterator. now pass on to the function doing
  // the real work.
  do_reinit(face_no);
}



template <int dim, int spacedim>
void
FEFaceValues<dim, spacedim>::reinit(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const typename Triangulation<dim, spacedim>::face_iterator &face)
{
  const auto face_n = cell->face_iterator_to_index(face);
  reinit(cell, face_n);
}



template <int dim, int spacedim>
void
FEFaceValues<dim, spacedim>::do_reinit(const unsigned int face_no)
{
  this->present_face_no = face_no;

  // first of all, set the present_face_index (if available)
  const typename Triangulation<dim, spacedim>::cell_iterator cell =
    this->present_cell;
  this->present_face_index = cell->face_index(face_no);

  if (this->update_flags & update_mapping)
    {
      this->get_mapping().fill_fe_face_values(this->present_cell,
                                              face_no,
                                              this->quadrature,
                                              *this->mapping_data,
                                              this->mapping_output);
    }

  this->get_fe().fill_fe_face_values(this->present_cell,
                                     face_no,
                                     this->quadrature,
                                     this->get_mapping(),
                                     *this->mapping_data,
                                     this->mapping_output,
                                     *this->fe_data,
                                     this->finite_element_output);

  const_cast<unsigned int &>(this->n_quadrature_points) =
    this->quadrature[this->quadrature.size() == 1 ? 0 : face_no].size();
}


/* ---------------------------- FESubFaceValues ---------------------------- */


template <int dim, int spacedim>
const unsigned int FESubfaceValues<dim, spacedim>::dimension;



template <int dim, int spacedim>
const unsigned int FESubfaceValues<dim, spacedim>::integral_dimension;



template <int dim, int spacedim>
FESubfaceValues<dim, spacedim>::FESubfaceValues(
  const Mapping<dim, spacedim>       &mapping,
  const FiniteElement<dim, spacedim> &fe,
  const Quadrature<dim - 1>          &quadrature,
  const UpdateFlags                   update_flags)
  : FEFaceValuesBase<dim, spacedim>(fe.n_dofs_per_cell(),
                                    update_flags,
                                    mapping,
                                    fe,
                                    quadrature)
{
  initialize(update_flags);
}



template <int dim, int spacedim>
FESubfaceValues<dim, spacedim>::FESubfaceValues(
  const Mapping<dim, spacedim>       &mapping,
  const FiniteElement<dim, spacedim> &fe,
  const hp::QCollection<dim - 1>     &quadrature,
  const UpdateFlags                   update_flags)
  : FESubfaceValues(mapping, fe, quadrature[0], update_flags)
{
  AssertDimension(quadrature.size(), 1);
}



template <int dim, int spacedim>
FESubfaceValues<dim, spacedim>::FESubfaceValues(
  const FiniteElement<dim, spacedim> &fe,
  const Quadrature<dim - 1>          &quadrature,
  const UpdateFlags                   update_flags)
  : FEFaceValuesBase<dim, spacedim>(
      fe.n_dofs_per_cell(),
      update_flags,
      fe.reference_cell().template get_default_linear_mapping<dim, spacedim>(),
      fe,
      quadrature)
{
  initialize(update_flags);
}



template <int dim, int spacedim>
FESubfaceValues<dim, spacedim>::FESubfaceValues(
  const FiniteElement<dim, spacedim> &fe,
  const hp::QCollection<dim - 1>     &quadrature,
  const UpdateFlags                   update_flags)
  : FESubfaceValues(fe, quadrature[0], update_flags)
{
  AssertDimension(quadrature.size(), 1);
}



template <int dim, int spacedim>
void
FESubfaceValues<dim, spacedim>::initialize(const UpdateFlags update_flags)
{
  const UpdateFlags flags = this->compute_update_flags(update_flags);

  // initialize the base classes
  if (flags & update_mapping)
    this->mapping_output.initialize(this->max_n_quadrature_points, flags);
  this->finite_element_output.initialize(this->max_n_quadrature_points,
                                         *this->fe,
                                         flags);

  // then get objects into which the FE and the Mapping can store
  // intermediate data used across calls to reinit. this can be done
  // in parallel
  Threads::Task<
    std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>>
    fe_get_data =
      Threads::new_task(&FiniteElement<dim, spacedim>::get_subface_data,
                        *this->fe,
                        flags,
                        *this->mapping,
                        this->quadrature[0],
                        this->finite_element_output);
  Threads::Task<
    std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>>
    mapping_get_data;
  if (flags & update_mapping)
    mapping_get_data =
      Threads::new_task(&Mapping<dim, spacedim>::get_subface_data,
                        *this->mapping,
                        flags,
                        this->quadrature[0]);

  this->update_flags = flags;

  // then collect answers from the two task above
  this->fe_data = std::move(fe_get_data.return_value());
  if (flags & update_mapping)
    this->mapping_data = std::move(mapping_get_data.return_value());
  else
    this->mapping_data =
      std::make_unique<typename Mapping<dim, spacedim>::InternalDataBase>();
}



template <int dim, int spacedim>
template <bool lda>
void
FESubfaceValues<dim, spacedim>::reinit(
  const TriaIterator<DoFCellAccessor<dim, spacedim, lda>> &cell,
  const unsigned int                                       face_no,
  const unsigned int                                       subface_no)
{
  // assert that the finite elements passed to the constructor and
  // used by the DoFHandler used by this cell, are the same
  Assert(static_cast<const FiniteElementData<dim> &>(*this->fe) ==
           static_cast<const FiniteElementData<dim> &>(
             cell->get_dof_handler().get_fe(cell->active_fe_index())),
         (typename FEValuesBase<dim, spacedim>::ExcFEDontMatch()));
  AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);
  // We would like to check for subface_no < cell->face(face_no)->n_children(),
  // but unfortunately the current function is also called for
  // faces without children (see tests/fe/mapping.cc). Therefore,
  // we must use following workaround of two separate assertions
  Assert(cell->face(face_no)->has_children() ||
           subface_no < GeometryInfo<dim>::max_children_per_face,
         ExcIndexRange(subface_no,
                       0,
                       GeometryInfo<dim>::max_children_per_face));
  Assert(!cell->face(face_no)->has_children() ||
           subface_no < cell->face(face_no)->n_active_descendants(),
         ExcIndexRange(subface_no,
                       0,
                       cell->face(face_no)->n_active_descendants()));
  Assert(cell->has_children() == false,
         ExcMessage("You can't use subface data for cells that are "
                    "already refined. Iterate over their children "
                    "instead in these cases."));

  this->maybe_invalidate_previous_present_cell(cell);
  this->present_cell = {cell};

  // this was the part of the work that is dependent on the actual
  // data type of the iterator. now pass on to the function doing
  // the real work.
  do_reinit(face_no, subface_no);
}



template <int dim, int spacedim>
template <bool lda>
void
FESubfaceValues<dim, spacedim>::reinit(
  const TriaIterator<DoFCellAccessor<dim, spacedim, lda>>    &cell,
  const typename Triangulation<dim, spacedim>::face_iterator &face,
  const typename Triangulation<dim, spacedim>::face_iterator &subface)
{
  reinit(cell,
         cell->face_iterator_to_index(face),
         face->child_iterator_to_index(subface));
}



template <int dim, int spacedim>
void
FESubfaceValues<dim, spacedim>::reinit(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const unsigned int                                          subface_no)
{
  AssertIndexRange(face_no, GeometryInfo<dim>::faces_per_cell);
  // We would like to check for subface_no < cell->face(face_no)->n_children(),
  // but unfortunately the current function is also called for
  // faces without children for periodic faces, which have hanging nodes on
  // the other side (see include/deal.II/matrix_free/mapping_info.templates.h).
  AssertIndexRange(subface_no,
                   (cell->has_periodic_neighbor(face_no) ?
                      cell->periodic_neighbor(face_no)
                        ->face(cell->periodic_neighbor_face_no(face_no))
                        ->n_children() :
                      cell->face(face_no)->n_children()));

  this->maybe_invalidate_previous_present_cell(cell);
  this->present_cell = {cell};

  // this was the part of the work that is dependent on the actual
  // data type of the iterator. now pass on to the function doing
  // the real work.
  do_reinit(face_no, subface_no);
}



template <int dim, int spacedim>
void
FESubfaceValues<dim, spacedim>::reinit(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const typename Triangulation<dim, spacedim>::face_iterator &face,
  const typename Triangulation<dim, spacedim>::face_iterator &subface)
{
  reinit(cell,
         cell->face_iterator_to_index(face),
         face->child_iterator_to_index(subface));
}



template <int dim, int spacedim>
void
FESubfaceValues<dim, spacedim>::do_reinit(const unsigned int face_no,
                                          const unsigned int subface_no)
{
  this->present_face_no = face_no;

  // first of all, set the present_face_index (if available)
  const typename Triangulation<dim, spacedim>::cell_iterator cell =
    this->present_cell;

  if (!cell->face(face_no)->has_children())
    // no subfaces at all, so set present_face_index to this face rather
    // than any subface
    this->present_face_index = cell->face_index(face_no);
  else if (dim != 3)
    this->present_face_index = cell->face(face_no)->child_index(subface_no);
  else
    {
      // this is the same logic we use in cell->neighbor_child_on_subface(). See
      // there for an explanation of the different cases
      unsigned int subface_index = numbers::invalid_unsigned_int;
      switch (cell->subface_case(face_no))
        {
          case internal::SubfaceCase<3>::case_x:
          case internal::SubfaceCase<3>::case_y:
          case internal::SubfaceCase<3>::case_xy:
            subface_index = cell->face(face_no)->child_index(subface_no);
            break;
          case internal::SubfaceCase<3>::case_x1y2y:
          case internal::SubfaceCase<3>::case_y1x2x:
            subface_index = cell->face(face_no)
                              ->child(subface_no / 2)
                              ->child_index(subface_no % 2);
            break;
          case internal::SubfaceCase<3>::case_x1y:
          case internal::SubfaceCase<3>::case_y1x:
            switch (subface_no)
              {
                case 0:
                case 1:
                  subface_index =
                    cell->face(face_no)->child(0)->child_index(subface_no);
                  break;
                case 2:
                  subface_index = cell->face(face_no)->child_index(1);
                  break;
                default:
                  DEAL_II_ASSERT_UNREACHABLE();
              }
            break;
          case internal::SubfaceCase<3>::case_x2y:
          case internal::SubfaceCase<3>::case_y2x:
            switch (subface_no)
              {
                case 0:
                  subface_index = cell->face(face_no)->child_index(0);
                  break;
                case 1:
                case 2:
                  subface_index =
                    cell->face(face_no)->child(1)->child_index(subface_no - 1);
                  break;
                default:
                  DEAL_II_ASSERT_UNREACHABLE();
              }
            break;
          default:
            DEAL_II_ASSERT_UNREACHABLE();
            break;
        }
      Assert(subface_index != numbers::invalid_unsigned_int,
             ExcInternalError());
      this->present_face_index = subface_index;
    }

  // now ask the mapping and the finite element to do the actual work
  if (this->update_flags & update_mapping)
    {
      this->get_mapping().fill_fe_subface_values(this->present_cell,
                                                 face_no,
                                                 subface_no,
                                                 this->quadrature[0],
                                                 *this->mapping_data,
                                                 this->mapping_output);
    }

  this->get_fe().fill_fe_subface_values(this->present_cell,
                                        face_no,
                                        subface_no,
                                        this->quadrature[0],
                                        this->get_mapping(),
                                        *this->mapping_data,
                                        this->mapping_output,
                                        *this->fe_data,
                                        this->finite_element_output);
}
#endif

/*------------------------- Explicit Instantiations --------------------------*/

#include "fe/fe_values.inst"

DEAL_II_NAMESPACE_CLOSE
