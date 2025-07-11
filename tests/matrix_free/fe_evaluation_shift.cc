// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include "../tests.h"


// Test FEEvaluationShift which accesses the data of the neighbor of a cell for
// FE_Q
//
// @note This test is only run if vectorization is enabled.

namespace dealii
{
  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number,
            typename VectorizedArrayType>
  class FEEvaluationShift : public FEEvaluation<dim,
                                                fe_degree,
                                                n_q_points_1d,
                                                n_components_,
                                                Number,
                                                VectorizedArrayType>
  {
  public:
    using BaseClass = FEEvaluation<dim,
                                   fe_degree,
                                   n_q_points_1d,
                                   n_components_,
                                   Number,
                                   VectorizedArrayType>;

    explicit FEEvaluationShift(
      const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
      const unsigned int                                  dof_no  = 0,
      const unsigned int                                  quad_no = 0,
      const unsigned int first_selected_component                 = 0,
      const unsigned int active_fe_index   = numbers::invalid_unsigned_int,
      const unsigned int active_quad_index = numbers::invalid_unsigned_int)
      : BaseClass(matrix_free,
                  dof_no,
                  quad_no,
                  first_selected_component,
                  active_fe_index,
                  active_quad_index)
      , shape_info_base(this->data){

        };

    const internal::MatrixFreeFunctions::ShapeInfo<Number> *shape_info_base;

    void
    reinit(unsigned int cell_batch_index,
           const internal::MatrixFreeFunctions::ShapeInfo<Number> *shape_info =
             nullptr);

    void
    reinit(unsigned int cell_batch_index,
           unsigned int face_number,
           const internal::MatrixFreeFunctions::ShapeInfo<Number> *shape_info =
             nullptr);
  };

  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number,
            typename VectorizedArrayType>
  inline void
  FEEvaluationShift<dim,
                    fe_degree,
                    n_q_points_1d,
                    n_components_,
                    Number,
                    VectorizedArrayType>::
    reinit(const unsigned int cell_batch_index,
           const internal::MatrixFreeFunctions::ShapeInfo<Number> *shape_info)
  {
    Assert(this->mapped_geometry == nullptr,
           ExcMessage(
             "FEEvaluation was initialized without a matrix-free object."
             " Integer indexing is not possible"));
    if (this->mapped_geometry != nullptr)
      return;

    Assert(this->dof_info != nullptr, ExcNotInitialized());
    Assert(this->mapping_data != nullptr, ExcNotInitialized());
    this->cell = cell_batch_index;
    this->cell_type =
      this->matrix_free->get_mapping_info().get_cell_type(cell_batch_index);

    const unsigned int offsets =
      this->mapping_data->data_index_offsets[cell_batch_index];
    this->jacobian = &this->mapping_data->jacobians[0][offsets];
    this->J_value  = &this->mapping_data->JxW_values[offsets];

    for (unsigned int i = 0; i < VectorizedArrayType::size(); ++i)
      this->cell_ids[i] = cell_batch_index * VectorizedArrayType::size() + i;

    this->interior_face = true;

    if (shape_info != nullptr)
      {
        AssertDimension(shape_info->n_q_points, this->data->n_q_points);
        AssertDimension(shape_info->dofs_per_component_on_cell,
                        this->data->dofs_per_component_on_cell);
        this->data = shape_info;
      }
    else
      this->data = this->shape_info_base;

    if constexpr (running_in_debug_mode())
      {
        this->is_reinitialized           = true;
        this->dof_values_initialized     = false;
        this->values_quad_initialized    = false;
        this->gradients_quad_initialized = false;
        this->hessians_quad_initialized  = false;
      }
  }

  template <int dim,
            int fe_degree,
            int n_q_points_1d,
            int n_components_,
            typename Number,
            typename VectorizedArrayType>
  inline void
  FEEvaluationShift<dim,
                    fe_degree,
                    n_q_points_1d,
                    n_components_,
                    Number,
                    VectorizedArrayType>::
    reinit(const unsigned int cell_batch_index,
           const unsigned int face_number,
           const internal::MatrixFreeFunctions::ShapeInfo<Number> *shape_info)
  {
    Assert(
      this->quad_no <
        this->matrix_free->get_mapping_info().face_data_by_cells.size(),
      ExcMessage(
        "You must set MatrixFree::AdditionalData::mapping_update_flags_faces_by_cells to use the present reinit method."));
    AssertIndexRange(face_number, GeometryInfo<dim>::faces_per_cell);
    AssertIndexRange(cell_batch_index,
                     this->matrix_free->get_mapping_info().cell_type.size());
    Assert(this->mapped_geometry == nullptr,
           ExcMessage(
             "FEEvaluation was initialized without a matrix-free object."
             " Integer indexing is not possible"));
    if (this->mapped_geometry != nullptr)
      return;
    Assert(this->matrix_free != nullptr, ExcNotInitialized());

    this->cell_type =
      this->matrix_free->get_mapping_info().cell_type[cell_batch_index];
    this->cell                 = numbers::invalid_unsigned_int;
    this->face_orientations[0] = 0;
    this->subface_index        = GeometryInfo<dim>::max_children_per_cell;
    this->face_numbers[0]      = face_number;
    this->dof_access_index =
      internal::MatrixFreeFunctions::DoFInfo::dof_access_cell;

    this->interior_face = false;

    const unsigned int offsets =
      this->mapping_data->data_index_offsets[cell_batch_index];
    this->jacobian = &this->mapping_data->jacobians[0][offsets];
    this->J_value  = &this->mapping_data->JxW_values[offsets];

    unsigned int n_lanes = VectorizedArrayType::size();
    for (unsigned int i = 0;
         i <
         this->matrix_free->n_active_entries_per_cell_batch(cell_batch_index);
         ++i)
      {
        // compute actual (non vectorized) cell ID
        const unsigned int cell_this = cell_batch_index * n_lanes + i;
        // compute face ID
        unsigned int face_index =
          this->matrix_free->get_cell_and_face_to_plain_faces()(
            cell_batch_index, face_number, i);

        const auto &faces =
          this->matrix_free->get_face_info(face_index / n_lanes);
        // get cell ID on both sides of face
        auto cell_m = faces.cells_interior[face_index % n_lanes];
        auto cell_p = faces.cells_exterior[face_index % n_lanes];

        const bool is_interior_face = cell_m != cell_this;

        Assert(cell_m == cell_this || cell_p == cell_this, ExcInternalError());

        // compare the IDs with the given cell ID
        if (is_interior_face)
          {
            this->cell_ids[i]     = cell_m; // neighbor has the other ID
            this->face_numbers[i] = faces.interior_face_no;
          }
        else
          {
            this->cell_ids[i]     = cell_p;
            this->face_numbers[i] = faces.exterior_face_no;
          }
      }

    if (shape_info != nullptr)
      {
        AssertDimension(shape_info->n_q_points, this->data->n_q_points);
        AssertDimension(shape_info->dofs_per_component_on_cell,
                        this->data->dofs_per_component_on_cell);
        this->data = shape_info;
      }
    else
      this->data = this->shape_info_base;

    if constexpr (running_in_debug_mode())
      {
        this->is_reinitialized           = true;
        this->dof_values_initialized     = false;
        this->values_quad_initialized    = false;
        this->gradients_quad_initialized = false;
        this->hessians_quad_initialized  = false;
      }
  }
} // namespace dealii

template <typename Number>
Quadrature<1>
shift_1d_quadrature(const Quadrature<1> &quadrature_1D, const Number shift)
{
  const auto &points  = quadrature_1D.get_points();
  const auto &weights = quadrature_1D.get_weights();

  std::vector<Point<1>> points_shift;

  for (const auto &i : points)
    points_shift.push_back(i + Point<1>(shift));

  return {points_shift, weights};
}


template <int dim,
          int fe_degree,
          int n_points                 = fe_degree + 1,
          typename Number              = double,
          typename VectorizedArrayType = VectorizedArray<Number>>
void
test(const dealii::FE_Poly<dim> &fe)
{
  if (VectorizedArrayType::size() == 1)
    return;

  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  unsigned int x_repetitions = 4;

  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_rectangle(tria,
                                            {x_repetitions, 1},
                                            Point<dim>(0.0, 0.0),
                                            Point<dim>(1. * x_repetitions, 1.0),
                                            true);

  std::vector<dealii::GridTools::PeriodicFacePair<
    typename dealii::Triangulation<dim>::cell_iterator>>
    periodic_faces;

  if (dim >= 1)
    dealii::GridTools::collect_periodic_faces(tria, 0, 1, 0, periodic_faces);

  if (dim >= 2)
    dealii::GridTools::collect_periodic_faces(tria, 2, 3, 1, periodic_faces);

  tria.add_periodicity(periodic_faces);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  MappingQ<dim> mapping(1);
  QGauss<1>     quad(n_points);

  AffineConstraints<Number> constraint;

  std::vector<dealii::GridTools::PeriodicFacePair<
    typename dealii::DoFHandler<dim>::cell_iterator>>
    periodic_faces_dof;

  if (dim >= 1)
    dealii::GridTools::collect_periodic_faces(
      dof_handler, 0, 1, 0, periodic_faces_dof);

  if (dim >= 2)
    dealii::GridTools::collect_periodic_faces(
      dof_handler, 2, 3, 1, periodic_faces_dof);

  DoFTools::make_periodicity_constraints<dim, dim, Number>(periodic_faces_dof,
                                                           constraint);
  constraint.close();

  using MF = MatrixFree<dim, Number, VectorizedArrayType>;

  typename MF::AdditionalData additional_data;
  additional_data.mapping_update_flags                = update_values;
  additional_data.mapping_update_flags_inner_faces    = update_values;
  additional_data.mapping_update_flags_boundary_faces = update_values;
  additional_data.mapping_update_flags_faces_by_cells = update_values;
  additional_data.hold_all_faces_to_owned_cells       = true;

  MF matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraint, quad, additional_data);

  VectorType dst;

  matrix_free.initialize_dof_vector(dst);

  std::array<internal::MatrixFreeFunctions::ShapeInfo<Number>, 2 * dim>
    shape_info_shift;

  const QGauss<1> quadrature_1D(n_points);

  const Quadrature<1> quadrature_1D_shift_minus =
    shift_1d_quadrature(quadrature_1D, -1.);
  const Quadrature<1> quadrature_1D_shift_plus =
    shift_1d_quadrature(quadrature_1D, +1.);

  dealii::internal::MatrixFreeFunctions::ShapeInfo<Number> shape_info_base(
    quadrature_1D, fe);
  dealii::internal::MatrixFreeFunctions::ShapeInfo<Number>
    shape_info_shift_plus(quadrature_1D_shift_plus, fe);
  dealii::internal::MatrixFreeFunctions::ShapeInfo<Number>
    shape_info_shift_minus(quadrature_1D_shift_minus, fe);

  for (unsigned int f = 0; f < 2 * dim; ++f)
    {
      shape_info_shift[f].data.assign(dim, shape_info_base.data[0]);

      shape_info_shift[f].n_q_points = shape_info_base.n_q_points;
      shape_info_shift[f].dofs_per_component_on_cell =
        shape_info_base.dofs_per_component_on_cell;

      switch (f)
        {
          case 0:
            shape_info_shift[f].data[0] = shape_info_shift_minus.data[0];
            break;
          case 1:
            shape_info_shift[f].data[0] = shape_info_shift_plus.data[0];
            break;
          case 2:
            shape_info_shift[f].data[1] = shape_info_shift_minus.data[0];
            break;
          case 3:
            shape_info_shift[f].data[1] = shape_info_shift_plus.data[0];
            break;
          case 4:
            shape_info_shift[f].data[2] = shape_info_shift_minus.data[0];
            break;
          case 5:
            shape_info_shift[f].data[2] = shape_info_shift_plus.data[0];
            break;
          default:
            AssertThrow(false, ExcNotImplemented());
        }
    }

  FEEvaluationShift<dim, fe_degree, n_points, 1, Number, VectorizedArrayType>
    phi_m(matrix_free);
  FEEvaluationShift<dim, fe_degree, n_points, 1, Number, VectorizedArrayType>
    phi_p(matrix_free);

  // fill dst vector
  matrix_free.template cell_loop<VectorType, VectorType>(
    [&](const auto &, auto &dst, const auto &src, const auto range) {
      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          phi_m.reinit(cell);
          for (unsigned int v = 0;
               v < matrix_free.n_active_entries_per_cell_batch(cell);
               ++v)
            {
              const auto cell_no = cell * VectorizedArrayType::size() + v;
              for (unsigned int i = 0; i < phi_m.static_dofs_per_component; ++i)
                phi_m.begin_dof_values()[i][v] = cell_no;
            }
          phi_m.evaluate(EvaluationFlags::values);
          for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
            phi_m.submit_value(phi_m.get_value(q), q);
          phi_m.integrate_scatter(EvaluationFlags::values, dst);
        }
    },
    dst,
    dst);

  // read out dst vector and print cell solutions
  matrix_free.template loop_cell_centric<VectorType, VectorType>(
    [&](const auto &, auto &dst, const auto &src, const auto range) {
      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          phi_m.reinit(cell);
          phi_m.read_dof_values(dst);

          for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
            {
              deallog << "lane " << v << std::endl;
              for (unsigned int i = 0; i < phi_m.static_dofs_per_component; ++i)
                deallog << phi_m.begin_dof_values()[i][v] << " ";
              deallog << std::endl;
            }

          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
               face++)
            {
              deallog << "face " << face << std::endl;
              phi_p.reinit(cell, face);
              phi_p.read_dof_values(dst);
              phi_p.evaluate(EvaluationFlags::values);
              for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
                {
                  deallog << "lane " << v << std::endl;
                  for (unsigned int i = 0; i < phi_p.static_dofs_per_component;
                       ++i)
                    deallog << phi_p.begin_dof_values()[i][v] << " ";
                  deallog << std::endl;

                  for (unsigned int i = 0; i < phi_p.static_n_q_points; ++i)
                    deallog << phi_p.begin_values()[i][v] << " ";
                  deallog << std::endl;
                }

              deallog << "shifted  shape info" << std::endl;
              phi_m.reinit(cell, &shape_info_shift[face]);
              phi_m.read_dof_values(dst);

              for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
                {
                  deallog << "lane " << v << std::endl;
                  for (unsigned int i = 0; i < phi_m.static_dofs_per_component;
                       ++i)
                    deallog << phi_m.begin_dof_values()[i][v] << " ";
                  deallog << std::endl;
                }

              phi_p.reinit(cell, face, &shape_info_shift[face]);
              phi_p.read_dof_values(dst);
              phi_p.evaluate(EvaluationFlags::values);
              for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
                {
                  deallog << "lane " << v << std::endl;
                  for (unsigned int i = 0; i < phi_p.static_dofs_per_component;
                       ++i)
                    deallog << phi_p.begin_dof_values()[i][v] << " ";
                  deallog << std::endl;

                  for (unsigned int i = 0; i < phi_p.static_n_q_points; ++i)
                    deallog << phi_p.begin_values()[i][v] << " ";
                  deallog << std::endl;
                }
              deallog << std::endl;
            }
        }
    },
    dst,
    dst,
    false,
    MF::DataAccessOnFaces::gradients);
}

int
main()
{
  initlog();

  constexpr int dim        = 2;
  constexpr int degree     = 1;
  constexpr int n_q_points = degree + 1;

  test<dim,
       degree,
       n_q_points,
       double,
       VectorizedArray<double,
                       VectorizedArray<double>::size() < 2 ?
                         VectorizedArray<double>::size() :
                         2>>(FE_Q<2>(1));

  test<dim,
       degree,
       n_q_points,
       double,
       VectorizedArray<double,
                       VectorizedArray<double>::size() < 2 ?
                         VectorizedArray<double>::size() :
                         2>>(FE_DGQ<2>(1));
}
