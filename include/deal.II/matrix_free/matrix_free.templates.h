// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2020 by the deal.II authors
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

#ifndef dealii_matrix_free_templates_h
#define dealii_matrix_free_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_consensus_algorithms.h>
#include <deal.II/base/polynomials_piecewise.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_poly.h>
#include <deal.II/fe/fe_q_dg0.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/matrix_free/face_info.h>
#include <deal.II/matrix_free/face_setup_internal.h>
#include <deal.II/matrix_free/matrix_free.h>

#ifdef DEAL_II_WITH_TBB
#  include <deal.II/base/parallel.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  include <tbb/concurrent_unordered_map.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
#endif

#include <fstream>


DEAL_II_NAMESPACE_OPEN



// --------------------- MatrixFree -----------------------------------

template <int dim, typename Number, typename VectorizedArrayType>
MatrixFree<dim, Number, VectorizedArrayType>::MatrixFree()
  : Subscriptor()
  , indices_are_initialized(false)
  , mapping_is_initialized(false)
  , mg_level(numbers::invalid_unsigned_int)
{}



template <int dim, typename Number, typename VectorizedArrayType>
MatrixFree<dim, Number, VectorizedArrayType>::MatrixFree(
  const MatrixFree<dim, Number, VectorizedArrayType> &other)
  : Subscriptor()
{
  copy_from(other);
}



template <int dim, typename Number, typename VectorizedArrayType>
std::pair<unsigned int, unsigned int>
MatrixFree<dim, Number, VectorizedArrayType>::create_cell_subrange_hp_by_index(
  const std::pair<unsigned int, unsigned int> &range,
  const unsigned int                           fe_index,
  const unsigned int                           dof_handler_index) const
{
  if (dof_info[dof_handler_index].max_fe_index == 0)
    return range;

  AssertIndexRange(fe_index, dof_info[dof_handler_index].max_fe_index);
  const std::vector<unsigned int> &fe_indices =
    dof_info[dof_handler_index].cell_active_fe_index;

  if (fe_indices.empty() == true ||
      dof_handlers[dof_handler_index]->get_fe_collection().size() == 1)
    return range;
  else
    {
      // the range over which we are searching must be ordered, otherwise we
      // got a range that spans over too many cells
#ifdef DEBUG
      for (unsigned int i = range.first + 1; i < range.second; ++i)
        Assert(
          fe_indices[i] >= fe_indices[i - 1],
          ExcMessage(
            "Cell range must be over sorted range of FE indices in hp-case!"));
      AssertIndexRange(range.first, fe_indices.size() + 1);
      AssertIndexRange(range.second, fe_indices.size() + 1);
#endif
      std::pair<unsigned int, unsigned int> return_range;
      return_range.first = std::lower_bound(fe_indices.begin() + range.first,
                                            fe_indices.begin() + range.second,
                                            fe_index) -
                           fe_indices.begin();
      return_range.second =
        std::lower_bound(fe_indices.begin() + return_range.first,
                         fe_indices.begin() + range.second,
                         fe_index + 1) -
        fe_indices.begin();
      Assert(return_range.first >= range.first &&
               return_range.second <= range.second,
             ExcInternalError());
      return return_range;
    }
}



namespace
{
  class FaceRangeCompartor
  {
  public:
    FaceRangeCompartor(const std::vector<unsigned int> &fe_indices,
                       const bool                       include,
                       const bool                       only_face_type)
      : fe_indices(fe_indices)
      , include(include)
      , only_face_type(only_face_type)
    {}

    template <int vectorization_width>
    bool
    operator()(const internal::MatrixFreeFunctions::FaceToCellTopology<
                 vectorization_width> &           face,
               const std::array<unsigned int, 2> &fe_index)
    {
      const unsigned int face_type = face.face_type;

      if (only_face_type)
        return include ? (face_type <= fe_index[0]) : (face_type < fe_index[0]);

      const std::array<unsigned int, 2> fe_index_face = {
        {face_type, fe_indices[face.cells_interior[0] / vectorization_width]}};

      return include ? (fe_index_face <= fe_index) : (fe_index_face < fe_index);
    }

    template <int vectorization_width>
    bool
    operator()(const internal::MatrixFreeFunctions::FaceToCellTopology<
                 vectorization_width> &           face,
               const std::array<unsigned int, 3> &fe_index)
    {
      const unsigned int face_type = face.face_type;

      if (only_face_type)
        return include ? (face_type <= fe_index[0]) : (face_type < fe_index[0]);

      const std::array<unsigned int, 3> fe_index_face = {
        {face_type,
         fe_indices[face.cells_interior[0] / vectorization_width],
         fe_indices[face.cells_exterior[0] / vectorization_width]}};

      return include ? (fe_index_face <= fe_index) : (fe_index_face < fe_index);
    }

  private:
    const std::vector<unsigned int> &fe_indices;
    const bool                       include;
    const bool                       only_face_type;
  };
} // namespace



template <int dim, typename Number, typename VectorizedArrayType>
void
MatrixFree<dim, Number, VectorizedArrayType>::renumber_dofs(
  std::vector<types::global_dof_index> &renumbering,
  const unsigned int                    dof_handler_index)
{
  AssertIndexRange(dof_handler_index, dof_info.size());
  dof_info[dof_handler_index].compute_dof_renumbering(renumbering);
}



template <int dim, typename Number, typename VectorizedArrayType>
const DoFHandler<dim> &
MatrixFree<dim, Number, VectorizedArrayType>::get_dof_handler(
  const unsigned int dof_handler_index) const
{
  AssertIndexRange(dof_handler_index, n_components());

  return *(dof_handlers[dof_handler_index]);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename DoFHandlerType>
const DoFHandlerType &
MatrixFree<dim, Number, VectorizedArrayType>::get_dof_handler(
  const unsigned int dof_handler_index) const
{
  AssertIndexRange(dof_handler_index, n_components());

  auto dh =
    dynamic_cast<const DoFHandlerType *>(&*dof_handlers[dof_handler_index]);
  Assert(dh != nullptr, ExcNotInitialized());

  return *dh;
}



template <int dim, typename Number, typename VectorizedArrayType>
typename DoFHandler<dim>::cell_iterator
MatrixFree<dim, Number, VectorizedArrayType>::get_cell_iterator(
  const unsigned int cell_batch_index,
  const unsigned int lane_index,
  const unsigned int dof_handler_index) const
{
  AssertIndexRange(dof_handler_index, dof_handlers.size());
  AssertIndexRange(cell_batch_index, task_info.cell_partition_data.back());
  AssertIndexRange(lane_index, n_components_filled(cell_batch_index));

  std::pair<unsigned int, unsigned int> index =
    cell_level_index[cell_batch_index * VectorizedArrayType::size() +
                     lane_index];
  return typename DoFHandler<dim>::cell_iterator(
    &dof_handlers[dof_handler_index]->get_triangulation(),
    index.first,
    index.second,
    &*dof_handlers[dof_handler_index]);
}



template <int dim, typename Number, typename VectorizedArrayType>
std::pair<int, int>
MatrixFree<dim, Number, VectorizedArrayType>::get_cell_level_and_index(
  const unsigned int cell_batch_index,
  const unsigned int lane_index) const
{
  AssertIndexRange(cell_batch_index, task_info.cell_partition_data.back());
  AssertIndexRange(lane_index, n_components_filled(cell_batch_index));

  std::pair<int, int> level_index_pair =
    cell_level_index[cell_batch_index * VectorizedArrayType::size() +
                     lane_index];
  return level_index_pair;
}



template <int dim, typename Number, typename VectorizedArrayType>
std::pair<typename DoFHandler<dim>::cell_iterator, unsigned int>
MatrixFree<dim, Number, VectorizedArrayType>::get_face_iterator(
  const unsigned int face_batch_index,
  const unsigned int lane_index,
  const bool         interior,
  const unsigned int fe_component) const
{
  AssertIndexRange(fe_component, dof_handlers.size());
  AssertIndexRange(face_batch_index,
                   n_inner_face_batches() +
                     (interior ? n_boundary_face_batches() : 0));

  AssertIndexRange(lane_index,
                   n_active_entries_per_face_batch(face_batch_index));

  const internal::MatrixFreeFunctions::FaceToCellTopology<
    VectorizedArrayType::size()>
    face2cell_info = get_face_info(face_batch_index);

  const unsigned int cell_index = interior ?
                                    face2cell_info.cells_interior[lane_index] :
                                    face2cell_info.cells_exterior[lane_index];

  std::pair<unsigned int, unsigned int> index = cell_level_index[cell_index];

  return {typename DoFHandler<dim>::cell_iterator(
            &dof_handlers[fe_component]->get_triangulation(),
            index.first,
            index.second,
            &*dof_handlers[fe_component]),
          interior ? face2cell_info.interior_face_no :
                     face2cell_info.exterior_face_no};
}



template <int dim, typename Number, typename VectorizedArrayType>
typename DoFHandler<dim>::active_cell_iterator
MatrixFree<dim, Number, VectorizedArrayType>::get_hp_cell_iterator(
  const unsigned int cell_batch_index,
  const unsigned int lane_index,
  const unsigned int dof_handler_index) const
{
  AssertIndexRange(dof_handler_index, dof_handlers.size());
  AssertIndexRange(cell_batch_index, task_info.cell_partition_data.back());
  AssertIndexRange(lane_index, n_components_filled(cell_batch_index));

  std::pair<unsigned int, unsigned int> index =
    cell_level_index[cell_batch_index * VectorizedArrayType::size() +
                     lane_index];
  return typename DoFHandler<dim>::cell_iterator(
    &dof_handlers[dof_handler_index]->get_triangulation(),
    index.first,
    index.second,
    &*dof_handlers[dof_handler_index]);
}



template <int dim, typename Number, typename VectorizedArrayType>
void
MatrixFree<dim, Number, VectorizedArrayType>::copy_from(
  const MatrixFree<dim, Number, VectorizedArrayType> &v)
{
  clear();
  dof_handlers               = v.dof_handlers;
  dof_info                   = v.dof_info;
  constraint_pool_data       = v.constraint_pool_data;
  constraint_pool_row_index  = v.constraint_pool_row_index;
  mapping_info               = v.mapping_info;
  shape_info                 = v.shape_info;
  cell_level_index           = v.cell_level_index;
  cell_level_index_end_local = v.cell_level_index_end_local;
  task_info                  = v.task_info;
  face_info                  = v.face_info;
  indices_are_initialized    = v.indices_are_initialized;
  mapping_is_initialized     = v.mapping_is_initialized;
  mg_level                   = v.mg_level;
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename number2, int q_dim>
void
MatrixFree<dim, Number, VectorizedArrayType>::internal_reinit(
  const std::shared_ptr<hp::MappingCollection<dim>> &    mapping,
  const std::vector<const DoFHandler<dim, dim> *> &      dof_handler,
  const std::vector<const AffineConstraints<number2> *> &constraint,
  const std::vector<IndexSet> &                          locally_owned_dofs,
  const std::vector<hp::QCollection<q_dim>> &            quad,
  const typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    &additional_data)
{
  // Store the level of the mesh to be worked on.
  this->mg_level = additional_data.mg_level;

  // Reads out the FE information and stores the shape function values,
  // gradients and Hessians for quadrature points.
  {
    unsigned int n_components = 0;
    for (unsigned int no = 0; no < dof_handler.size(); ++no)
      n_components += dof_handler[no]->get_fe(0).n_base_elements();
    const unsigned int n_quad             = quad.size();
    unsigned int       n_fe_in_collection = 0;
    for (unsigned int no = 0; no < dof_handler.size(); ++no)
      n_fe_in_collection =
        std::max(n_fe_in_collection,
                 dof_handler[no]->get_fe_collection().size());
    unsigned int n_quad_in_collection = 0;
    for (unsigned int q = 0; q < n_quad; ++q)
      n_quad_in_collection = std::max(n_quad_in_collection, quad[q].size());
    shape_info.reinit(TableIndices<4>(
      n_components, n_quad, n_fe_in_collection, n_quad_in_collection));
    for (unsigned int no = 0, c = 0; no < dof_handler.size(); no++)
      for (unsigned int b = 0; b < dof_handler[no]->get_fe(0).n_base_elements();
           ++b, ++c)
        for (unsigned int fe_no = 0;
             fe_no < dof_handler[no]->get_fe_collection().size();
             ++fe_no)
          for (unsigned int nq = 0; nq < n_quad; nq++)
            for (unsigned int q_no = 0; q_no < quad[nq].size(); ++q_no)
              shape_info(c, nq, fe_no, q_no)
                .reinit(quad[nq][q_no], dof_handler[no]->get_fe(fe_no), b);
  }

  if (additional_data.initialize_indices == true)
    {
      clear();
      Assert(dof_handler.size() > 0, ExcMessage("No DoFHandler is given."));
      AssertDimension(dof_handler.size(), constraint.size());
      AssertDimension(dof_handler.size(), locally_owned_dofs.size());

      // set variables that are independent of FE
      if (Utilities::MPI::job_supports_mpi() == true)
        {
          const parallel::TriangulationBase<dim> *dist_tria =
            dynamic_cast<const parallel::TriangulationBase<dim> *>(
              &(dof_handler[0]->get_triangulation()));
          task_info.communicator = dist_tria != nullptr ?
                                     dist_tria->get_communicator() :
                                     MPI_COMM_SELF;
          task_info.my_pid =
            Utilities::MPI::this_mpi_process(task_info.communicator);
          task_info.n_procs =
            Utilities::MPI::n_mpi_processes(task_info.communicator);

          task_info.communicator_sm = additional_data.communicator_sm;
        }
      else
        {
          task_info.communicator    = MPI_COMM_SELF;
          task_info.communicator_sm = MPI_COMM_SELF;
          task_info.my_pid          = 0;
          task_info.n_procs         = 1;
        }

      initialize_dof_handlers(dof_handler, additional_data);
      for (unsigned int no = 0; no < dof_handler.size(); ++no)
        {
          dof_info[no].store_plain_indices =
            additional_data.store_plain_indices;
          dof_info[no].global_base_element_offset =
            no > 0 ? dof_info[no - 1].global_base_element_offset +
                       dof_handler[no - 1]->get_fe(0).n_base_elements() :
                     0;
        }

        // initialize the basic multithreading information that needs to be
        // passed to the DoFInfo structure
#ifdef DEAL_II_WITH_TBB
      if (additional_data.tasks_parallel_scheme != AdditionalData::none &&
          MultithreadInfo::n_threads() > 1)
        {
          task_info.scheme =
            internal::MatrixFreeFunctions::TaskInfo::TasksParallelScheme(
              static_cast<int>(additional_data.tasks_parallel_scheme));
          task_info.block_size = additional_data.tasks_block_size;
        }
      else
#endif
        task_info.scheme = internal::MatrixFreeFunctions::TaskInfo::none;

      // set dof_indices together with constraint_indicator and
      // constraint_pool_data. It also reorders the way cells are gone through
      // (to separate cells with overlap to other processors from others
      // without).
      initialize_indices(constraint, locally_owned_dofs, additional_data);
    }

  // initialize bare structures
  else if (dof_info.size() != dof_handler.size())
    {
      initialize_dof_handlers(dof_handler, additional_data);
      std::vector<unsigned int>  dummy;
      std::vector<unsigned char> dummy2;
      task_info.vectorization_length = VectorizedArrayType::size();
      task_info.n_active_cells       = cell_level_index.size();
      task_info.create_blocks_serial(
        dummy, 1, false, dummy, false, dummy, dummy, dummy2);

      for (unsigned int i = 0; i < dof_info.size(); ++i)
        {
          Assert(dof_handler[i]->get_fe_collection().size() == 1,
                 ExcNotImplemented());
          dof_info[i].dimension = dim;
          dof_info[i].n_base_elements =
            dof_handler[i]->get_fe(0).n_base_elements();
          dof_info[i].n_components.resize(dof_info[i].n_base_elements);
          dof_info[i].start_components.resize(dof_info[i].n_base_elements + 1);
          for (unsigned int c = 0; c < dof_info[i].n_base_elements; ++c)
            {
              dof_info[i].n_components[c] =
                dof_handler[i]->get_fe(0).element_multiplicity(c);
              for (unsigned int l = 0; l < dof_info[i].n_components[c]; ++l)
                dof_info[i].component_to_base_index.push_back(c);
              dof_info[i].start_components[c + 1] =
                dof_info[i].start_components[c] + dof_info[i].n_components[c];
            }
          dof_info[i].dofs_per_cell.push_back(
            dof_handler[i]->get_fe(0).n_dofs_per_cell());

          // if indices are not initialized, the cell_level_index might not be
          // divisible by the vectorization length. But it must be for
          // mapping_info...
          while (cell_level_index.size() % VectorizedArrayType::size() != 0)
            cell_level_index.push_back(cell_level_index.back());
        }
    }


  // subdivide cell, face and boundary face partitioner data, s.t., all
  // ranges have the same active_fe_indices
  if (task_info.scheme != internal::MatrixFreeFunctions::TaskInfo::
                            TasksParallelScheme::partition_color)
    {
      // ... for cells
      {
        auto &ptr  = task_info.cell_partition_data_hp_ptr;
        auto &data = task_info.cell_partition_data_hp;

        ptr  = {0};
        data = {};
        data.reserve(2 * task_info.cell_partition_data.size());

        for (unsigned int i = 0; i < task_info.cell_partition_data.size() - 1;
             ++i)
          {
            if (task_info.cell_partition_data[i + 1] >
                task_info.cell_partition_data[i])
              {
                std::pair<unsigned int, unsigned int> range{
                  task_info.cell_partition_data[i],
                  task_info.cell_partition_data[i + 1]};

                if (range.second > range.first)
                  for (unsigned int i = 0; i < this->n_active_fe_indices(); ++i)
                    {
                      const auto cell_subrange =
                        this->create_cell_subrange_hp_by_index(range, i);

                      if (cell_subrange.second <= cell_subrange.first)
                        continue;

                      data.push_back(cell_subrange.first);
                      data.push_back(cell_subrange.second);
                    }
              }

            ptr.push_back(data.size() / 2);
          }
      }

      // ... for inner faces
      if (task_info.face_partition_data.empty() == false)
        {
          auto &ptr  = task_info.face_partition_data_hp_ptr;
          auto &data = task_info.face_partition_data_hp;

          ptr  = {0};
          data = {};
          data.reserve(2 * task_info.face_partition_data.size());

          const auto create_inner_face_subrange_hp_by_index =
            [&](const std::pair<unsigned int, unsigned int> &range,
                const unsigned int                           face_type,
                const unsigned int                           fe_index_interior,
                const unsigned int                           fe_index_exterior,
                const unsigned int                           dof_handler_index =
                  0) -> std::pair<unsigned int, unsigned int> {
            const unsigned int n_face_types =
              std::max<unsigned int>(dim - 1, 1);

            if ((dof_info[dof_handler_index].max_fe_index == 0 ||
                 dof_handlers[dof_handler_index]->get_fe_collection().size() ==
                   1) &&
                n_face_types == 1)
              return range;

            AssertIndexRange(fe_index_interior,
                             dof_info[dof_handler_index].max_fe_index);
            AssertIndexRange(fe_index_exterior,
                             dof_info[dof_handler_index].max_fe_index);
            const std::vector<unsigned int> &fe_indices =
              dof_info[dof_handler_index].cell_active_fe_index;
            if (fe_indices.empty() == true && n_face_types == 1)
              return range;
            else
              {
                const bool only_face_type =
                  dof_info[dof_handler_index].max_fe_index == 0 ||
                  dof_handlers[dof_handler_index]->get_fe_collection().size() ==
                    1 ||
                  fe_indices.empty() == true;

                std::pair<unsigned int, unsigned int> return_range;
                return_range.first =
                  std::lower_bound(
                    face_info.faces.begin() + range.first,
                    face_info.faces.begin() + range.second,
                    std::array<unsigned int, 3>{
                      {face_type, fe_index_interior, fe_index_exterior}},
                    FaceRangeCompartor(fe_indices, false, only_face_type)) -
                  face_info.faces.begin();
                return_range.second =
                  std::lower_bound(
                    face_info.faces.begin() + return_range.first,
                    face_info.faces.begin() + range.second,
                    std::array<unsigned int, 3>{
                      {face_type, fe_index_interior, fe_index_exterior}},
                    FaceRangeCompartor(fe_indices, true, only_face_type)) -
                  face_info.faces.begin();
                Assert(return_range.first >= range.first &&
                         return_range.second <= range.second,
                       ExcInternalError());
                return return_range;
              }
          };

          for (unsigned int i = 0; i < task_info.face_partition_data.size() - 1;
               ++i)
            {
              if (task_info.face_partition_data[i + 1] >
                  task_info.face_partition_data[i])
                {
                  std::pair<unsigned int, unsigned int> range{
                    task_info.face_partition_data[i],
                    task_info.face_partition_data[i + 1]};

                  if (range.second > range.first)
                    for (unsigned int t = 0;
                         t < std::max<unsigned int>(dim - 1, 1);
                         ++t)
                      for (unsigned int i = 0; i < this->n_active_fe_indices();
                           ++i)
                        for (unsigned int j = 0;
                             j < this->n_active_fe_indices();
                             ++j)
                          {
                            const auto subrange =
                              create_inner_face_subrange_hp_by_index(range,
                                                                     t,
                                                                     i,
                                                                     j);

                            if (subrange.second <= subrange.first)
                              continue;

                            data.push_back(subrange.first);
                            data.push_back(subrange.second);
                          }
                }

              ptr.push_back(data.size() / 2);
            }
        }

      // ... for boundary faces
      if (task_info.boundary_partition_data.empty() == false)
        {
          auto &ptr  = task_info.boundary_partition_data_hp_ptr;
          auto &data = task_info.boundary_partition_data_hp;

          ptr  = {0};
          data = {};
          data.reserve(2 * task_info.boundary_partition_data.size());

          const auto create_boundary_face_subrange_hp_by_index =
            [&](const std::pair<unsigned int, unsigned int> &range,
                const unsigned int                           face_type,
                const unsigned int                           fe_index,
                const unsigned int                           dof_handler_index =
                  0) -> std::pair<unsigned int, unsigned int> {
            const unsigned int n_face_types =
              std::max<unsigned int>(dim - 1, 1);

            if ((dof_info[dof_handler_index].max_fe_index == 0 ||
                 dof_handlers[dof_handler_index]->get_fe_collection().size() ==
                   1) &&
                n_face_types == 1)
              return range;

            AssertIndexRange(fe_index,
                             dof_info[dof_handler_index].max_fe_index);
            const std::vector<unsigned int> &fe_indices =
              dof_info[dof_handler_index].cell_active_fe_index;
            if (fe_indices.empty() == true && n_face_types == 1)
              return range;
            else
              {
                const bool only_face_type =
                  dof_info[dof_handler_index].max_fe_index == 0 ||
                  dof_handlers[dof_handler_index]->get_fe_collection().size() ==
                    1 ||
                  fe_indices.empty() == true;

                std::pair<unsigned int, unsigned int> return_range;
                return_range.first =
                  std::lower_bound(
                    face_info.faces.begin() + range.first,
                    face_info.faces.begin() + range.second,
                    std::array<unsigned int, 2>{{face_type, fe_index}},
                    FaceRangeCompartor(fe_indices, false, only_face_type)) -
                  face_info.faces.begin();
                return_range.second =
                  std::lower_bound(
                    face_info.faces.begin() + return_range.first,
                    face_info.faces.begin() + range.second,
                    std::array<unsigned int, 2>{{face_type, fe_index}},
                    FaceRangeCompartor(fe_indices, true, only_face_type)) -
                  face_info.faces.begin();
                Assert(return_range.first >= range.first &&
                         return_range.second <= range.second,
                       ExcInternalError());
                return return_range;
              }
          };

          for (unsigned int i = 0;
               i < task_info.boundary_partition_data.size() - 1;
               ++i)
            {
              if (task_info.boundary_partition_data[i + 1] >
                  task_info.boundary_partition_data[i])
                {
                  std::pair<unsigned int, unsigned int> range{
                    task_info.boundary_partition_data[i],
                    task_info.boundary_partition_data[i + 1]};

                  if (range.second > range.first)
                    for (unsigned int t = 0;
                         t < std::max<unsigned int>(dim - 1, 1);
                         ++t)
                      for (unsigned int i = 0; i < this->n_active_fe_indices();
                           ++i)
                        {
                          const auto cell_subrange =
                            create_boundary_face_subrange_hp_by_index(range,
                                                                      t,
                                                                      i);

                          if (cell_subrange.second <= cell_subrange.first)
                            continue;

                          data.push_back(cell_subrange.first);
                          data.push_back(cell_subrange.second);
                        }
                }

              ptr.push_back(data.size() / 2);
            }
        }
    }

  // Evaluates transformations from unit to real cell, Jacobian determinants,
  // quadrature points in real space, based on the ordering of the cells
  // determined in @p extract_local_to_global_indices.
  if (additional_data.initialize_mapping == true)
    {
      if (dof_handler.size() > 1)
        {
          // check if all DoHandlers are in the same hp-mode; and if hp-
          // capabilities are enabled: check if active_fe_indices of all
          // DoFHandler are the same.
          for (unsigned int i = 1; i < dof_handler.size(); ++i)
            {
              Assert(dof_handler[0]->has_hp_capabilities() ==
                       dof_handler[i]->has_hp_capabilities(),
                     ExcNotImplemented());

              if (dof_handler[0]->has_hp_capabilities())
                {
                  Assert(dof_info[0].cell_active_fe_index ==
                           dof_info[i].cell_active_fe_index,
                         ExcNotImplemented());
                }
            }
        }

      mapping_info.initialize(
        dof_handler[0]->get_triangulation(),
        cell_level_index,
        face_info,
        dof_handler[0]->has_hp_capabilities() ?
          dof_info[0].cell_active_fe_index :
          std::vector<unsigned int>(),
        mapping,
        quad,
        additional_data.mapping_update_flags,
        additional_data.mapping_update_flags_boundary_faces,
        additional_data.mapping_update_flags_inner_faces,
        additional_data.mapping_update_flags_faces_by_cells);

      mapping_is_initialized = true;
    }
}



template <int dim, typename Number, typename VectorizedArrayType>
void
MatrixFree<dim, Number, VectorizedArrayType>::update_mapping(
  const Mapping<dim> &mapping)
{
  update_mapping(std::make_shared<hp::MappingCollection<dim>>(mapping));
}



template <int dim, typename Number, typename VectorizedArrayType>
void
MatrixFree<dim, Number, VectorizedArrayType>::update_mapping(
  const std::shared_ptr<hp::MappingCollection<dim>> &mapping)
{
  AssertDimension(shape_info.size(1), mapping_info.cell_data.size());
  mapping_info.update_mapping(dof_handlers[0]->get_triangulation(),
                              cell_level_index,
                              face_info,
                              dof_info[0].cell_active_fe_index,
                              mapping);
}



template <int dim, typename Number, typename VectorizedArrayType>
template <int spacedim>
bool
MatrixFree<dim, Number, VectorizedArrayType>::is_supported(
  const FiniteElement<dim, spacedim> &fe)
{
  return internal::MatrixFreeFunctions::ShapeInfo<double>::is_supported(fe);
}



namespace internal
{
  namespace MatrixFreeFunctions
  {
    // steps through all children and adds the active cells recursively
    template <typename InIterator>
    void
    resolve_cell(const InIterator &                                  cell,
                 std::vector<std::pair<unsigned int, unsigned int>> &cell_its,
                 const unsigned int subdomain_id)
    {
      if (cell->has_children())
        for (unsigned int child = 0; child < cell->n_children(); ++child)
          resolve_cell(cell->child(child), cell_its, subdomain_id);
      else if (subdomain_id == numbers::invalid_subdomain_id ||
               cell->subdomain_id() == subdomain_id)
        {
          Assert(cell->is_active(), ExcInternalError());
          cell_its.emplace_back(cell->level(), cell->index());
        }
    }
  } // namespace MatrixFreeFunctions
} // namespace internal



template <int dim, typename Number, typename VectorizedArrayType>
void
MatrixFree<dim, Number, VectorizedArrayType>::initialize_dof_handlers(
  const std::vector<const DoFHandler<dim, dim> *> &dof_handler_in,
  const AdditionalData &                           additional_data)
{
  cell_level_index.clear();
  dof_handlers.resize(dof_handler_in.size());
  for (unsigned int no = 0; no < dof_handler_in.size(); ++no)
    dof_handlers[no] = dof_handler_in[no];

  dof_info.resize(dof_handlers.size());
  for (unsigned int no = 0; no < dof_handlers.size(); ++no)
    dof_info[no].vectorization_length = VectorizedArrayType::size();

  const unsigned int n_mpi_procs = task_info.n_procs;
  const unsigned int my_pid      = task_info.my_pid;

  const Triangulation<dim> &tria  = dof_handlers[0]->get_triangulation();
  const unsigned int        level = additional_data.mg_level;
  if (level == numbers::invalid_unsigned_int)
    {
      if (n_mpi_procs == 1)
        cell_level_index.reserve(tria.n_active_cells());
      // For serial Triangulations always take all cells
      const unsigned int subdomain_id =
        (dynamic_cast<const parallel::TriangulationBase<dim> *>(
           &dof_handlers[0]->get_triangulation()) != nullptr) ?
          my_pid :
          numbers::invalid_subdomain_id;

      // Go through cells on zeroth level and then successively step down into
      // children. This gives a z-ordering of the cells, which is beneficial
      // when setting up neighboring relations between cells for thread
      // parallelization
      for (const auto &cell : tria.cell_iterators_on_level(0))
        internal::MatrixFreeFunctions::resolve_cell(cell,
                                                    cell_level_index,
                                                    subdomain_id);

      Assert(n_mpi_procs > 1 ||
               cell_level_index.size() == tria.n_active_cells(),
             ExcInternalError());
    }
  else
    {
      AssertIndexRange(level, tria.n_global_levels());
      if (level < tria.n_levels())
        {
          cell_level_index.reserve(tria.n_cells(level));
          for (const auto &cell : tria.cell_iterators_on_level(level))
            if (cell->level_subdomain_id() == my_pid)
              cell_level_index.emplace_back(cell->level(), cell->index());
        }
    }

  // All these are cells local to this processor. Therefore, set
  // cell_level_index_end_local to the size of cell_level_index.
  cell_level_index_end_local = cell_level_index.size();
}



namespace internal
{
#ifdef DEAL_II_WITH_TBB

  inline void
  fill_index_subrange(
    const unsigned int                                        begin,
    const unsigned int                                        end,
    const std::vector<std::pair<unsigned int, unsigned int>> &cell_level_index,
    tbb::concurrent_unordered_map<std::pair<unsigned int, unsigned int>,
                                  unsigned int> &             map)
  {
    if (cell_level_index.empty())
      return;
    unsigned int cell = begin;
    if (cell == 0)
      map.insert(std::make_pair(cell_level_index[cell++], 0U));
    for (; cell < end; ++cell)
      if (cell_level_index[cell] != cell_level_index[cell - 1])
        map.insert(std::make_pair(cell_level_index[cell], cell));
  }

  template <int dim>
  inline void
  fill_connectivity_subrange(
    const unsigned int                                        begin,
    const unsigned int                                        end,
    const dealii::Triangulation<dim> &                        tria,
    const std::vector<std::pair<unsigned int, unsigned int>> &cell_level_index,
    const tbb::concurrent_unordered_map<std::pair<unsigned int, unsigned int>,
                                        unsigned int> &       map,
    DynamicSparsityPattern &connectivity_direct)
  {
    std::vector<types::global_dof_index> new_indices;
    for (unsigned int cell = begin; cell < end; ++cell)
      {
        new_indices.clear();
        typename dealii::Triangulation<dim>::cell_iterator dcell(
          &tria, cell_level_index[cell].first, cell_level_index[cell].second);
        for (auto f : dcell->face_indices())
          {
            // Only inner faces couple different cells
            if (dcell->at_boundary(f) == false &&
                dcell->neighbor_or_periodic_neighbor(f)->level_subdomain_id() ==
                  dcell->level_subdomain_id())
              {
                std::pair<unsigned int, unsigned int> level_index(
                  dcell->neighbor_or_periodic_neighbor(f)->level(),
                  dcell->neighbor_or_periodic_neighbor(f)->index());
                auto it = map.find(level_index);
                if (it != map.end())
                  {
                    const unsigned int neighbor_cell = it->second;
                    if (neighbor_cell != cell)
                      new_indices.push_back(neighbor_cell);
                  }
              }
          }
        std::sort(new_indices.begin(), new_indices.end());
        connectivity_direct.add_entries(cell,
                                        new_indices.begin(),
                                        std::unique(new_indices.begin(),
                                                    new_indices.end()));
      }
  }

  inline void
  fill_connectivity_indirect_subrange(
    const unsigned int            begin,
    const unsigned int            end,
    const DynamicSparsityPattern &connectivity_direct,
    DynamicSparsityPattern &      connectivity)
  {
    std::vector<types::global_dof_index> new_indices;
    for (unsigned int block = begin; block < end; ++block)
      {
        new_indices.clear();
        for (DynamicSparsityPattern::iterator it =
               connectivity_direct.begin(block);
             it != connectivity_direct.end(block);
             ++it)
          {
            new_indices.push_back(it->column());
            for (DynamicSparsityPattern::iterator it_neigh =
                   connectivity_direct.begin(it->column());
                 it_neigh != connectivity_direct.end(it->column());
                 ++it_neigh)
              if (it_neigh->column() != block)
                new_indices.push_back(it_neigh->column());
          }
        std::sort(new_indices.begin(), new_indices.end());
        connectivity.add_entries(block,
                                 new_indices.begin(),
                                 std::unique(new_indices.begin(),
                                             new_indices.end()));
      }
  }

#endif

  template <int dim, typename number>
  std::vector<bool>
  compute_dof_info(
    const std::vector<const dealii::AffineConstraints<number> *> &constraint,
    const std::vector<IndexSet> &locally_owned_dofs,
    const std::vector<SmartPointer<const dealii::DoFHandler<dim>>> &dof_handler,
    const Table<2, MatrixFreeFunctions::ShapeInfo<double>> &        shape_infos,
    const unsigned int               cell_level_index_end_local,
    const unsigned int               mg_level,
    const bool                       hold_all_faces_to_owned_cells,
    const std::vector<unsigned int> &cell_vectorization_category,
    const bool                       cell_vectorization_categories_strict,
    const bool                       do_face_integrals,
    const bool                       overlap_communication_computation,
    MatrixFreeFunctions::TaskInfo &  task_info,
    std::vector<std::pair<unsigned int, unsigned int>> &cell_level_index,
    std::vector<MatrixFreeFunctions::DoFInfo> &         dof_info,
    MatrixFreeFunctions::FaceSetup<dim> &               face_setup,
    MatrixFreeFunctions::ConstraintValues<double> &     constraint_values,
    const bool use_vector_data_exchanger_full)
  {
    if (do_face_integrals)
      face_setup.initialize(dof_handler[0]->get_triangulation(),
                            mg_level,
                            hold_all_faces_to_owned_cells,
                            cell_level_index);

    const unsigned int n_dof_handlers = dof_handler.size();
    const unsigned int n_active_cells = cell_level_index.size();

    const Triangulation<dim> &tria = dof_handler[0]->get_triangulation();

    AssertDimension(n_dof_handlers, locally_owned_dofs.size());
    AssertDimension(n_dof_handlers, constraint.size());

    std::vector<types::global_dof_index>                local_dof_indices;
    std::vector<std::vector<std::vector<unsigned int>>> lexicographic(
      n_dof_handlers);

    std::vector<bool> is_fe_dg(n_dof_handlers, false);

    bool cell_categorization_enabled = !cell_vectorization_category.empty();

    for (unsigned int no = 0; no < n_dof_handlers; ++no)
      {
        const dealii::hp::FECollection<dim> &fes =
          dof_handler[no]->get_fe_collection();

        if (fes.size() > 1)
          {
            Assert(cell_vectorization_category.empty(), ExcNotImplemented());
            dof_info[no].cell_active_fe_index.resize(
              n_active_cells, numbers::invalid_unsigned_int);
          }
        else if (cell_categorization_enabled == true)
          dof_info[no].cell_active_fe_index.resize(
            n_active_cells, numbers::invalid_unsigned_int);

        is_fe_dg[no] = fes[0].n_dofs_per_vertex() == 0;

        lexicographic[no].resize(fes.size());

        dof_info[no].fe_index_conversion.resize(fes.size());
        dof_info[no].max_fe_index = fes.size();

        dof_info[no].component_dof_indices_offset.clear();
        dof_info[no].component_dof_indices_offset.resize(fes.size());
        for (unsigned int fe_index = 0; fe_index < fes.size(); ++fe_index)
          {
            const FiniteElement<dim> &fe = fes[fe_index];
            // cache number of finite elements and dofs_per_cell
            dof_info[no].dofs_per_cell.push_back(fe.n_dofs_per_cell());
            dof_info[no].dofs_per_face.push_back(fe.n_dofs_per_face(
              0)); // we assume that all faces have the same number of dofs
            dof_info[no].dimension       = dim;
            dof_info[no].n_base_elements = fe.n_base_elements();
            dof_info[no].n_components.resize(dof_info[no].n_base_elements);
            dof_info[no].start_components.resize(dof_info[no].n_base_elements +
                                                 1);
            dof_info[no].component_to_base_index.clear();
            dof_info[no].component_dof_indices_offset[fe_index].push_back(0);
            dof_info[no].fe_index_conversion[fe_index].clear();
            for (unsigned int c = 0; c < dof_info[no].n_base_elements; ++c)
              {
                dof_info[no].n_components[c] = fe.element_multiplicity(c);
                for (unsigned int l = 0; l < dof_info[no].n_components[c]; ++l)
                  {
                    dof_info[no].component_to_base_index.push_back(c);
                    dof_info[no]
                      .component_dof_indices_offset[fe_index]
                      .push_back(dof_info[no]
                                   .component_dof_indices_offset[fe_index]
                                   .back() +
                                 fe.base_element(c).n_dofs_per_cell());
                    dof_info[no].fe_index_conversion[fe_index].push_back(
                      fe.base_element(c).degree);
                  }
                dof_info[no].start_components[c + 1] =
                  dof_info[no].start_components[c] +
                  dof_info[no].n_components[c];
                const auto &lex =
                  shape_infos(dof_info[no].global_base_element_offset + c,
                              fe_index)
                    .lexicographic_numbering;
                lexicographic[no][fe_index].insert(
                  lexicographic[no][fe_index].end(), lex.begin(), lex.end());
              }

            AssertDimension(lexicographic[no][fe_index].size(),
                            dof_info[no].dofs_per_cell[fe_index]);
            AssertDimension(
              dof_info[no].component_dof_indices_offset[fe_index].size() - 1,
              dof_info[no].start_components.back());
            AssertDimension(
              dof_info[no].component_dof_indices_offset[fe_index].back(),
              dof_info[no].dofs_per_cell[fe_index]);
          }

        // set locally owned range for each component
        Assert(locally_owned_dofs[no].is_contiguous(), ExcNotImplemented());
        dof_info[no].vector_partitioner =
          std::make_shared<Utilities::MPI::Partitioner>(locally_owned_dofs[no],
                                                        task_info.communicator);

        if (use_vector_data_exchanger_full == false)
          dof_info[no].vector_exchanger =
            std::make_shared<internal::MatrixFreeFunctions::VectorDataExchange::
                               PartitionerWrapper>(
              dof_info[no].vector_partitioner);
        else
          dof_info[no].vector_exchanger = std::make_shared<
            internal::MatrixFreeFunctions::VectorDataExchange::Full>(
            dof_info[no].vector_partitioner, task_info.communicator_sm);

        // initialize the arrays for indices
        const unsigned int n_components_total =
          dof_info[no].start_components.back();
        dof_info[no].row_starts.resize(n_active_cells * n_components_total + 1);
        dof_info[no].row_starts[0].first  = 0;
        dof_info[no].row_starts[0].second = 0;
        dof_info[no].dof_indices.reserve(
          (n_active_cells * dof_info[no].dofs_per_cell[0] * 3) / 2);

        // cache the constrained indices for use in matrix-vector products and
        // the like
        const types::global_dof_index
          start_index = dof_info[no].vector_partitioner->local_range().first,
          end_index   = dof_info[no].vector_partitioner->local_range().second;
        for (types::global_dof_index i = start_index; i < end_index; ++i)
          if (constraint[no]->is_constrained(i) == true)
            dof_info[no].constrained_dofs.push_back(
              static_cast<unsigned int>(i - start_index));
      }

    // extract all the global indices associated with the computation, and form
    // the ghost indices
    std::vector<unsigned int> subdomain_boundary_cells;
    for (unsigned int counter = 0; counter < n_active_cells; ++counter)
      {
        bool cell_at_subdomain_boundary =
          (face_setup.at_processor_boundary.size() > counter &&
           face_setup.at_processor_boundary[counter]) ||
          (overlap_communication_computation == false && task_info.n_procs > 1);

        for (unsigned int no = 0; no < n_dof_handlers; ++no)
          {
            // read indices from active cells
            if (mg_level == numbers::invalid_unsigned_int)
              {
                const DoFHandler<dim> *dofh = &*dof_handler[no];
                typename DoFHandler<dim>::active_cell_iterator cell_it(
                  &tria,
                  cell_level_index[counter].first,
                  cell_level_index[counter].second,
                  dofh);
                const unsigned int fe_index =
                  dofh->get_fe_collection().size() > 1 ?
                    cell_it->active_fe_index() :
                    0;
                if (dofh->get_fe_collection().size() > 1)
                  dof_info[no].cell_active_fe_index[counter] = fe_index;
                local_dof_indices.resize(dof_info[no].dofs_per_cell[fe_index]);
                cell_it->get_dof_indices(local_dof_indices);
                dof_info[no].read_dof_indices(local_dof_indices,
                                              lexicographic[no][fe_index],
                                              *constraint[no],
                                              counter,
                                              constraint_values,
                                              cell_at_subdomain_boundary);
                if (dofh->get_fe_collection().size() == 1 &&
                    cell_categorization_enabled)
                  {
                    AssertIndexRange(cell_it->active_cell_index(),
                                     cell_vectorization_category.size());
                    dof_info[no].cell_active_fe_index[counter] =
                      cell_vectorization_category[cell_it->active_cell_index()];
                  }
              }
            // we are requested to use a multigrid level
            else
              {
                const DoFHandler<dim> *dofh = dof_handler[no];
                AssertIndexRange(mg_level, tria.n_levels());
                typename DoFHandler<dim>::cell_iterator cell_it(
                  &tria,
                  cell_level_index[counter].first,
                  cell_level_index[counter].second,
                  dofh);
                local_dof_indices.resize(dof_info[no].dofs_per_cell[0]);
                cell_it->get_mg_dof_indices(local_dof_indices);
                dof_info[no].read_dof_indices(local_dof_indices,
                                              lexicographic[no][0],
                                              *constraint[no],
                                              counter,
                                              constraint_values,
                                              cell_at_subdomain_boundary);
                if (cell_categorization_enabled)
                  {
                    AssertIndexRange(cell_it->index(),
                                     cell_vectorization_category.size());
                    dof_info[no].cell_active_fe_index[counter] =
                      cell_vectorization_category[cell_level_index[counter]
                                                    .second];
                  }
              }
          }

        // if we found dofs on some FE component that belong to other
        // processors, the cell is added to the boundary cells.
        if (cell_at_subdomain_boundary == true &&
            counter < cell_level_index_end_local)
          subdomain_boundary_cells.push_back(counter);
      }

    task_info.n_active_cells = cell_level_index_end_local;
    task_info.n_ghost_cells  = n_active_cells - cell_level_index_end_local;

    // Finalize the creation of the ghost indices
    {
      std::vector<unsigned int> cells_with_ghosts(subdomain_boundary_cells);
      for (unsigned int c = cell_level_index_end_local; c < n_active_cells; ++c)
        cells_with_ghosts.push_back(c);
      for (unsigned int no = 0; no < n_dof_handlers; ++no)
        {
          if (do_face_integrals && mg_level != numbers::invalid_unsigned_int)
            {
              // in case of adaptivity, go through the cells on the next finer
              // level and check whether we need to get read access to some of
              // those entries for the mg flux matrices
              std::vector<types::global_dof_index> dof_indices;
              if (mg_level + 1 < tria.n_global_levels())
                for (const auto &cell :
                     dof_handler[no]->cell_iterators_on_level(mg_level + 1))
                  if (cell->level_subdomain_id() == task_info.my_pid)
                    for (const unsigned int f : cell->face_indices())
                      if ((cell->at_boundary(f) == false ||
                           cell->has_periodic_neighbor(f) == true) &&
                          cell->level() >
                            cell->neighbor_or_periodic_neighbor(f)->level() &&
                          cell->neighbor_or_periodic_neighbor(f)
                              ->level_subdomain_id() != task_info.my_pid)
                        {
                          dof_indices.resize(
                            cell->neighbor_or_periodic_neighbor(f)
                              ->get_fe()
                              .n_dofs_per_cell());
                          cell->neighbor_or_periodic_neighbor(f)
                            ->get_mg_dof_indices(dof_indices);
                          for (const auto dof_index : dof_indices)
                            dof_info[no].ghost_dofs.push_back(dof_index);
                        }
            }
          dof_info[no].assign_ghosts(cells_with_ghosts,
                                     task_info.communicator_sm,
                                     use_vector_data_exchanger_full);
        }
    }

    bool hp_functionality_enabled = false;
    for (const auto &dh : dof_handler)
      if (dh->get_fe_collection().size() > 1)
        hp_functionality_enabled = true;
    const unsigned int         n_lanes = task_info.vectorization_length;
    std::vector<unsigned int>  renumbering;
    std::vector<unsigned char> irregular_cells;

    Assert(
      task_info.scheme == internal::MatrixFreeFunctions::TaskInfo::none ||
        cell_vectorization_category.empty(),
      ExcMessage(
        "You explicitly requested re-categorization of cells; however, this "
        "feature is not available if threading is enabled. Please disable "
        "threading in MatrixFree by setting "
        "MatrixFree::Additional_data.tasks_parallel_scheme = MatrixFree<dim, double>::AdditionalData::none."));

    if (task_info.scheme == internal::MatrixFreeFunctions::TaskInfo::none)
      {
        const bool strict_categories =
          cell_vectorization_categories_strict || hp_functionality_enabled;
        unsigned int max_dofs_per_cell = 0;
        for (const auto &info : dof_info)
          for (const auto &dofs : info.dofs_per_cell)
            max_dofs_per_cell = std::max(max_dofs_per_cell, dofs);

        // Detect cells with the same parent to make sure they get scheduled
        // together in the loop, which increases data locality.
        std::vector<unsigned int> parent_relation(
          task_info.n_active_cells + task_info.n_ghost_cells,
          numbers::invalid_unsigned_int);
        std::map<std::pair<int, int>, std::vector<unsigned int>> cell_parents;
        for (unsigned int c = 0; c < cell_level_index_end_local; ++c)
          if (cell_level_index[c].first > 0)
            {
              typename Triangulation<dim>::cell_iterator cell(
                &tria, cell_level_index[c].first, cell_level_index[c].second);
              Assert(cell->level() > 0, ExcInternalError());
              cell_parents[std::make_pair(cell->parent()->level(),
                                          cell->parent()->index())]
                .push_back(c);
            }
        unsigned int position = 0;
        for (const auto &it : cell_parents)
          if (it.second.size() == GeometryInfo<dim>::max_children_per_cell)
            {
              for (auto i : it.second)
                parent_relation[i] = position;
              ++position;
            }
        task_info.create_blocks_serial(subdomain_boundary_cells,
                                       max_dofs_per_cell,
                                       hp_functionality_enabled,
                                       dof_info[0].cell_active_fe_index,
                                       strict_categories,
                                       parent_relation,
                                       renumbering,
                                       irregular_cells);
      }
    else
      {
        task_info.make_boundary_cells_divisible(subdomain_boundary_cells);

        // For strategy with blocking before partitioning: reorganize the
        // indices in order to overlap communication in MPI with computations:
        // Place all cells with ghost indices into one chunk. Also reorder cells
        // so that we can parallelize by threads
        task_info.initial_setup_blocks_tasks(subdomain_boundary_cells,
                                             renumbering,
                                             irregular_cells);
        task_info.guess_block_size(dof_info[0].dofs_per_cell[0]);

        unsigned int n_cell_batches_before =
          *(task_info.cell_partition_data.end() - 2);
        unsigned int n_ghost_slots =
          *(task_info.cell_partition_data.end() - 1) - n_cell_batches_before;

        unsigned int start_nonboundary = numbers::invalid_unsigned_int;
        if (task_info.scheme ==
              internal::MatrixFreeFunctions::TaskInfo::partition_color ||
            task_info.scheme == internal::MatrixFreeFunctions::TaskInfo::color)
          {
            // set up partitions. if we just use coloring without partitions, do
            // nothing here, assume all cells to belong to the zero partition
            // (that we otherwise use for MPI boundary cells)
            if (task_info.scheme ==
                internal::MatrixFreeFunctions::TaskInfo::color)
              {
                start_nonboundary =
                  task_info.n_procs > 1 ?
                    std::min(((task_info.cell_partition_data[2] -
                               task_info.cell_partition_data[1] +
                               task_info.block_size - 1) /
                              task_info.block_size) *
                               task_info.block_size,
                             task_info.cell_partition_data[3]) :
                    0;
              }
            else
              {
                if (task_info.n_procs > 1)
                  {
                    task_info.cell_partition_data[1] = 0;
                    task_info.cell_partition_data[2] =
                      task_info.cell_partition_data[3];
                  }
                start_nonboundary = task_info.cell_partition_data.back();
              }

            if (hp_functionality_enabled)
              {
                irregular_cells.resize(0);
                irregular_cells.resize(task_info.cell_partition_data.back() +
                                       2 * dof_info[0].max_fe_index);
                std::vector<std::vector<unsigned int>> renumbering_fe_index;
                renumbering_fe_index.resize(dof_info[0].max_fe_index);
                unsigned int counter;
                n_cell_batches_before = 0;
                for (counter = 0;
                     counter < std::min(start_nonboundary * n_lanes,
                                        task_info.n_active_cells);
                     counter++)
                  {
                    AssertIndexRange(counter, renumbering.size());
                    AssertIndexRange(renumbering[counter],
                                     dof_info[0].cell_active_fe_index.size());
                    renumbering_fe_index
                      [dof_info[0].cell_active_fe_index[renumbering[counter]]]
                        .push_back(renumbering[counter]);
                  }
                counter = 0;
                for (unsigned int j = 0; j < dof_info[0].max_fe_index; j++)
                  {
                    for (const auto jj : renumbering_fe_index[j])
                      renumbering[counter++] = jj;
                    irregular_cells[renumbering_fe_index[j].size() / n_lanes +
                                    n_cell_batches_before] =
                      renumbering_fe_index[j].size() % n_lanes;
                    n_cell_batches_before +=
                      (renumbering_fe_index[j].size() + n_lanes - 1) / n_lanes;
                    renumbering_fe_index[j].resize(0);
                  }

                for (counter = start_nonboundary * n_lanes;
                     counter < task_info.n_active_cells;
                     counter++)
                  {
                    renumbering_fe_index
                      [dof_info[0].cell_active_fe_index.empty() ?
                         0 :
                         dof_info[0].cell_active_fe_index[renumbering[counter]]]
                        .push_back(renumbering[counter]);
                  }
                counter = start_nonboundary * n_lanes;
                for (unsigned int j = 0; j < dof_info[0].max_fe_index; j++)
                  {
                    for (const auto jj : renumbering_fe_index[j])
                      renumbering[counter++] = jj;
                    irregular_cells[renumbering_fe_index[j].size() / n_lanes +
                                    n_cell_batches_before] =
                      renumbering_fe_index[j].size() % n_lanes;
                    n_cell_batches_before +=
                      (renumbering_fe_index[j].size() + n_lanes - 1) / n_lanes;
                  }
                AssertIndexRange(n_cell_batches_before,
                                 task_info.cell_partition_data.back() +
                                   2 * dof_info[0].max_fe_index + 1);
                irregular_cells.resize(n_cell_batches_before + n_ghost_slots);
                *(task_info.cell_partition_data.end() - 2) =
                  n_cell_batches_before;
                *(task_info.cell_partition_data.end() - 1) =
                  n_cell_batches_before + n_ghost_slots;
              }
          }

        task_info.n_blocks = (*(task_info.cell_partition_data.end() - 2) +
                              task_info.block_size - 1) /
                             task_info.block_size;

        DynamicSparsityPattern connectivity;
        connectivity.reinit(task_info.n_active_cells, task_info.n_active_cells);
        if (do_face_integrals)
          {
#ifdef DEAL_II_WITH_TBB
            // step 1: build map between the index in the matrix-free context
            // and the one in the triangulation
            tbb::concurrent_unordered_map<std::pair<unsigned int, unsigned int>,
                                          unsigned int>
              map;
            dealii::parallel::apply_to_subranges(
              0,
              cell_level_index.size(),
              [&cell_level_index, &map](const unsigned int begin,
                                        const unsigned int end) {
                fill_index_subrange(begin, end, cell_level_index, map);
              },
              50);

            // step 2: Make a list for all blocks with other blocks that write
            // to the cell (due to the faces that are associated to it)
            DynamicSparsityPattern connectivity_direct(connectivity.n_rows(),
                                                       connectivity.n_cols());
            dealii::parallel::apply_to_subranges(
              0,
              task_info.n_active_cells,
              [&cell_level_index, &tria, &map, &connectivity_direct](
                const unsigned int begin, const unsigned int end) {
                fill_connectivity_subrange<dim>(
                  begin, end, tria, cell_level_index, map, connectivity_direct);
              },
              20);
            connectivity_direct.symmetrize();

            // step 3: Include also interaction between neighbors one layer away
            // because faces might be assigned to cells differently
            dealii::parallel::apply_to_subranges(
              0,
              task_info.n_active_cells,
              [&connectivity_direct, &connectivity](const unsigned int begin,
                                                    const unsigned int end) {
                fill_connectivity_indirect_subrange(begin,
                                                    end,
                                                    connectivity_direct,
                                                    connectivity);
              },
              20);
#endif
          }
        if (task_info.n_active_cells > 0)
          dof_info[0].make_connectivity_graph(task_info,
                                              renumbering,
                                              connectivity);

        task_info.make_thread_graph(dof_info[0].cell_active_fe_index,
                                    connectivity,
                                    renumbering,
                                    irregular_cells,
                                    hp_functionality_enabled);

        Assert(irregular_cells.size() >= task_info.cell_partition_data.back(),
               ExcInternalError());

        irregular_cells.resize(task_info.cell_partition_data.back() +
                               n_ghost_slots);
        if (n_ghost_slots > 0)
          {
            for (unsigned int i = task_info.cell_partition_data.back();
                 i < task_info.cell_partition_data.back() + n_ghost_slots - 1;
                 ++i)
              irregular_cells[i] = 0;
            irregular_cells.back() = task_info.n_ghost_cells % n_lanes;
          }

        {
          unsigned int n_cells = 0;
          for (unsigned int i = 0; i < task_info.cell_partition_data.back();
               ++i)
            n_cells += irregular_cells[i] > 0 ? irregular_cells[i] : n_lanes;
          AssertDimension(n_cells, task_info.n_active_cells);
          n_cells = 0;
          for (unsigned int i = task_info.cell_partition_data.back();
               i < n_ghost_slots + task_info.cell_partition_data.back();
               ++i)
            n_cells += irregular_cells[i] > 0 ? irregular_cells[i] : n_lanes;
          AssertDimension(n_cells, task_info.n_ghost_cells);
        }

        task_info.cell_partition_data.push_back(
          task_info.cell_partition_data.back() + n_ghost_slots);
      }

      // Finally perform the renumbering. We also want to group several cells
      // together to a batch of cells for SIMD (vectorized) execution (where the
      // arithmetic operations of several cells will then be done
      // simultaneously).
#ifdef DEBUG
    {
      AssertDimension(renumbering.size(),
                      task_info.n_active_cells + task_info.n_ghost_cells);
      std::vector<unsigned int> sorted_renumbering(renumbering);
      std::sort(sorted_renumbering.begin(), sorted_renumbering.end());
      for (unsigned int i = 0; i < sorted_renumbering.size(); ++i)
        Assert(sorted_renumbering[i] == i, ExcInternalError());
    }
#endif
    {
      std::vector<std::pair<unsigned int, unsigned int>> cell_level_index_old;
      cell_level_index.swap(cell_level_index_old);
      cell_level_index.reserve(task_info.cell_partition_data.back() * n_lanes);
      unsigned int position_cell = 0;
      for (unsigned int i = 0; i < task_info.cell_partition_data.back(); ++i)
        {
          unsigned int n_comp =
            (irregular_cells[i] > 0) ? irregular_cells[i] : n_lanes;
          for (unsigned int j = 0; j < n_comp; ++j)
            cell_level_index.push_back(
              cell_level_index_old[renumbering[position_cell + j]]);

          // generate a cell and level index also when we have not filled up
          // vectorization_length cells. This is needed for MappingInfo when the
          // transformation data is initialized. We just set the value to the
          // last valid cell in that case.
          for (unsigned int j = n_comp; j < n_lanes; ++j)
            cell_level_index.push_back(
              cell_level_index_old[renumbering[position_cell + n_comp - 1]]);
          position_cell += n_comp;
        }
      AssertDimension(position_cell,
                      task_info.n_active_cells + task_info.n_ghost_cells);
      AssertDimension(cell_level_index.size(),
                      task_info.cell_partition_data.back() * n_lanes);
    }

    std::vector<unsigned int> constraint_pool_row_index;
    constraint_pool_row_index.resize(1, 0);
    for (const auto &it : constraint_values.constraints)
      constraint_pool_row_index.push_back(constraint_pool_row_index.back() +
                                          it.first.size());

    for (unsigned int no = 0; no < n_dof_handlers; ++no)
      dof_info[no].reorder_cells(task_info,
                                 renumbering,
                                 constraint_pool_row_index,
                                 irregular_cells);

    return is_fe_dg;
  }
} // namespace internal



template <int dim, typename Number, typename VectorizedArrayType>
template <typename number2>
void
MatrixFree<dim, Number, VectorizedArrayType>::initialize_indices(
  const std::vector<const AffineConstraints<number2> *> &constraint,
  const std::vector<IndexSet> &                          locally_owned_dofs,
  const AdditionalData &                                 additional_data)
{
  // insert possible ghost cells and construct face topology
  const bool do_face_integrals =
    (additional_data.mapping_update_flags_inner_faces |
     additional_data.mapping_update_flags_boundary_faces) != update_default;
  internal::MatrixFreeFunctions::FaceSetup<dim> face_setup;

  // create a vector with the dummy information about dofs in ShapeInfo
  // without the template of VectorizedArrayType
  Table<2, internal::MatrixFreeFunctions::ShapeInfo<double>> shape_info_dummy(
    shape_info.size(0), shape_info.size(2));
  {
    Quadrature<dim> quad(QGauss<dim>(1));
    Quadrature<dim> quad_simplex(QGaussSimplex<dim>(1));
    for (unsigned int no = 0, c = 0; no < dof_handlers.size(); no++)
      for (unsigned int b = 0;
           b < dof_handlers[no]->get_fe(0).n_base_elements();
           ++b, ++c)
        for (unsigned int fe_no = 0;
             fe_no < dof_handlers[no]->get_fe_collection().size();
             ++fe_no)
          shape_info_dummy(c, fe_no).reinit(
            dof_handlers[no]->get_fe(fe_no).reference_cell() ==
                ReferenceCells::get_hypercube<dim>() ?
              quad :
              quad_simplex,
            dof_handlers[no]->get_fe(fe_no),
            b);
  }

  const unsigned int n_lanes     = VectorizedArrayType::size();
  task_info.vectorization_length = n_lanes;
  internal::MatrixFreeFunctions::ConstraintValues<double> constraint_values;
  const std::vector<bool> is_fe_dg = internal::compute_dof_info(
    constraint,
    locally_owned_dofs,
    dof_handlers,
    shape_info_dummy,
    cell_level_index_end_local,
    additional_data.mg_level,
    additional_data.hold_all_faces_to_owned_cells,
    additional_data.cell_vectorization_category,
    additional_data.cell_vectorization_categories_strict,
    do_face_integrals,
    additional_data.overlap_communication_computation,
    task_info,
    cell_level_index,
    dof_info,
    face_setup,
    constraint_values,
    additional_data.communicator_sm != MPI_COMM_SELF);

  // set constraint pool from the std::map and reorder the indices
  std::vector<const std::vector<double> *> constraints(
    constraint_values.constraints.size());
  unsigned int length = 0;
  for (const auto &it : constraint_values.constraints)
    {
      AssertIndexRange(it.second, constraints.size());
      constraints[it.second] = &it.first;
      length += it.first.size();
    }
  constraint_pool_data.clear();
  constraint_pool_data.reserve(length);
  constraint_pool_row_index.reserve(constraint_values.constraints.size() + 1);
  constraint_pool_row_index.resize(1, 0);
  for (const auto &constraint : constraints)
    {
      Assert(constraint != nullptr, ExcInternalError());
      constraint_pool_data.insert(constraint_pool_data.end(),
                                  constraint->begin(),
                                  constraint->end());
      constraint_pool_row_index.push_back(constraint_pool_data.size());
    }

  AssertDimension(constraint_pool_data.size(), length);

  // Finally resort the faces and collect several faces for vectorization
  if ((additional_data.mapping_update_flags_inner_faces |
       additional_data.mapping_update_flags_boundary_faces) != update_default)
    {
      face_setup.generate_faces(dof_handlers[0]->get_triangulation(),
                                cell_level_index,
                                task_info);
      if (additional_data.mapping_update_flags_inner_faces != update_default)
        Assert(face_setup.refinement_edge_faces.empty(),
               ExcNotImplemented("Setting up data structures on MG levels with "
                                 "hanging nodes is currently not supported."));
      face_info.faces.clear();

      std::vector<bool> hard_vectorization_boundary(
        task_info.face_partition_data.size(), false);
      if (task_info.scheme == internal::MatrixFreeFunctions::TaskInfo::none &&
          task_info.partition_row_index[2] <
            task_info.face_partition_data.size())
        hard_vectorization_boundary[task_info.partition_row_index[2]] = true;
      else
        std::fill(hard_vectorization_boundary.begin(),
                  hard_vectorization_boundary.end(),
                  true);

      internal::MatrixFreeFunctions::collect_faces_vectorization(
        face_setup.inner_faces,
        hard_vectorization_boundary,
        task_info.face_partition_data,
        face_info.faces,
        dof_info[0].cell_active_fe_index);

      // on boundary faces, we must also respect the vectorization boundary of
      // the inner faces because we might have dependencies on ghosts of
      // remote vector entries for continuous elements
      internal::MatrixFreeFunctions::collect_faces_vectorization(
        face_setup.boundary_faces,
        hard_vectorization_boundary,
        task_info.boundary_partition_data,
        face_info.faces,
        dof_info[0].cell_active_fe_index);

      // for the other ghosted faces, there are no scheduling restrictions
      hard_vectorization_boundary.clear();
      hard_vectorization_boundary.resize(
        task_info.ghost_face_partition_data.size(), false);
      internal::MatrixFreeFunctions::collect_faces_vectorization(
        face_setup.inner_ghost_faces,
        hard_vectorization_boundary,
        task_info.ghost_face_partition_data,
        face_info.faces,
        dof_info[0].cell_active_fe_index);
      hard_vectorization_boundary.clear();
      hard_vectorization_boundary.resize(
        task_info.refinement_edge_face_partition_data.size(), false);
      internal::MatrixFreeFunctions::collect_faces_vectorization(
        face_setup.refinement_edge_faces,
        hard_vectorization_boundary,
        task_info.refinement_edge_face_partition_data,
        face_info.faces,
        dof_info[0].cell_active_fe_index);

      cell_level_index.resize(
        cell_level_index.size() +
        VectorizedArrayType::size() *
          (task_info.refinement_edge_face_partition_data[1] -
           task_info.refinement_edge_face_partition_data[0]));

      for (auto &di : dof_info)
        di.compute_face_index_compression(face_info.faces);

      // build the inverse map back from the faces array to
      // cell_and_face_to_plain_faces
      face_info.cell_and_face_to_plain_faces.reinit(
        TableIndices<3>(task_info.cell_partition_data.back(),
                        GeometryInfo<dim>::faces_per_cell,
                        VectorizedArrayType::size()),
        true);
      face_info.cell_and_face_to_plain_faces.fill(
        numbers::invalid_unsigned_int);
      face_info.cell_and_face_boundary_id.reinit(
        TableIndices<3>(task_info.cell_partition_data.back(),
                        GeometryInfo<dim>::faces_per_cell,
                        VectorizedArrayType::size()),
        true);
      face_info.cell_and_face_boundary_id.fill(numbers::invalid_boundary_id);

      for (unsigned int f = 0; f < task_info.ghost_face_partition_data.back();
           ++f)
        for (unsigned int v = 0; v < VectorizedArrayType::size() &&
                                 face_info.faces[f].cells_interior[v] !=
                                   numbers::invalid_unsigned_int;
             ++v)
          {
            TableIndices<3> index(face_info.faces[f].cells_interior[v] /
                                    VectorizedArrayType::size(),
                                  face_info.faces[f].interior_face_no,
                                  face_info.faces[f].cells_interior[v] %
                                    VectorizedArrayType::size());

            // Assert(cell_and_face_to_plain_faces(index) ==
            // numbers::invalid_unsigned_int,
            //       ExcInternalError("Should only visit each face once"));
            face_info.cell_and_face_to_plain_faces(index) =
              f * VectorizedArrayType::size() + v;
            if (face_info.faces[f].cells_exterior[v] !=
                numbers::invalid_unsigned_int)
              {
                TableIndices<3> index(face_info.faces[f].cells_exterior[v] /
                                        VectorizedArrayType::size(),
                                      face_info.faces[f].exterior_face_no,
                                      face_info.faces[f].cells_exterior[v] %
                                        VectorizedArrayType::size());
                // Assert(cell_and_face_to_plain_faces(index) ==
                // numbers::invalid_unsigned_int,
                //       ExcInternalError("Should only visit each face once"));
                face_info.cell_and_face_to_plain_faces(index) =
                  f * VectorizedArrayType::size() + v;
              }
            else
              face_info.cell_and_face_boundary_id(index) =
                types::boundary_id(face_info.faces[f].exterior_face_no);
          }

      // compute tighter index sets for various sets of face integrals
      unsigned int count = 0;
      for (auto &di : dof_info)
        di.compute_tight_partitioners(
          shape_info_dummy,
          *(task_info.cell_partition_data.end() - 2) *
            VectorizedArrayType::size(),
          VectorizedArrayType::size(),
          face_setup.inner_faces,
          face_setup.inner_ghost_faces,
          is_fe_dg[count++] && additional_data.hold_all_faces_to_owned_cells,
          task_info.communicator_sm,
          task_info.communicator_sm != MPI_COMM_SELF);
    }

  for (auto &di : dof_info)
    di.compute_vector_zero_access_pattern(task_info, face_info.faces);

#ifdef DEAL_II_WITH_MPI
  {
    // non-buffering mode is only supported if the indices of all cells are
    // contiguous for all dof_info objects.
    bool is_non_buffering_sm_supported = true;
    for (const auto &di : dof_info)
      {
        is_non_buffering_sm_supported &= di.dofs_per_cell.size() == 1;
        for (const auto &v : di.index_storage_variants[2])
          is_non_buffering_sm_supported &=
            (v == internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
                    contiguous);
      }

    is_non_buffering_sm_supported =
      Utilities::MPI::min(static_cast<unsigned int>(
                            is_non_buffering_sm_supported),
                          task_info.communicator);

    const MPI_Comm communicator_sm = this->task_info.communicator_sm;

    if (is_non_buffering_sm_supported)
      {
        // gather the ranks of the shared-memory domain
        const auto sm_procs = Utilities::MPI::mpi_processes_within_communicator(
          task_info.communicator, communicator_sm);

        // get my rank within the shared-memory domain
        const unsigned int my_sm_pid = std::distance(
          sm_procs.begin(),
          std::find(sm_procs.begin(), sm_procs.end(), task_info.my_pid));

        // process any dof_info
        const auto  di_index = 0;
        const auto &di       = dof_info[di_index];

        std::array<std::vector<std::pair<unsigned int, unsigned int>>, 3>
          cell_indices_contiguous_sm;

        // collect for each cell: sm rank + offset
        cell_indices_contiguous_sm[2] = [&]() {
          // result
          std::vector<std::pair<unsigned int, unsigned int>> cells(
            di.dof_indices_contiguous[2].size(),
            {numbers::invalid_unsigned_int, numbers::invalid_unsigned_int});

          // locally-owned cells
          std::vector<std::pair<unsigned int, types::global_cell_index>>
            cells_locally_owned;
          cells_locally_owned.reserve(n_cell_batches() * n_lanes);

          // shared-ghost cells (categorized according to sm-rank)
          std::vector<
            std::vector<std::pair<unsigned int, types::global_cell_index>>>
            cells_shared_ghosts(sm_procs.size());

          for (unsigned int cell = 0;
               cell < n_cell_batches() + n_ghost_cell_batches();
               ++cell)
            for (unsigned int v = 0; v < this->n_components_filled(cell); ++v)
              {
                const unsigned int index = cell * n_lanes + v;

                const auto cell_iterator = get_cell_iterator(cell, v, di_index);

                // determine global cell index
                const unsigned int local_cell_index =
                  di.dof_indices_contiguous[2][index] /
                  this->get_dofs_per_cell(di_index);

                const types::global_cell_index global_cell_index =
                  (additional_data.mg_level == numbers::invalid_unsigned_int) ?
                    cell_iterator->global_active_cell_index() :
                    cell_iterator->global_level_cell_index();

                // determine sm rank
                unsigned int sm_rank = my_sm_pid;

                if (cell < n_cell_batches())
                  {
                    // locally-owned cell
                    cells_locally_owned.emplace_back(local_cell_index,
                                                     global_cell_index);
                  }
                else
                  {
                    // ghost cell
                    const auto ptr =
                      std::find(sm_procs.begin(),
                                sm_procs.end(),
                                additional_data.mg_level ==
                                    numbers::invalid_unsigned_int ?
                                  cell_iterator->subdomain_id() :
                                  cell_iterator->level_subdomain_id());

                    if (ptr != sm_procs.end())
                      {
                        // shared ghost cell:
                        sm_rank = std::distance(sm_procs.begin(), ptr);
                        cells_shared_ghosts[sm_rank].emplace_back(
                          index, global_cell_index);
                      }
                    else
                      {
                        // remote ghost cell: nothing to do since the values are
                        // communicated (i.e. locally available)
                      }
                  }

                // write back result
                cells[cell * n_lanes + v] = {sm_rank, local_cell_index};
              }

          std::sort(cells_locally_owned.begin(),
                    cells_locally_owned.end(),
                    [](const auto &a, const auto &b) {
                      return a.second < b.second;
                    });

          // get offsets of shared cells
          Utilities::MPI::ConsensusAlgorithms::
            AnonymousProcess<dealii::types::global_dof_index, unsigned int>
              process(
                [&]() {
                  std::vector<unsigned int> targets;
                  for (unsigned int i = 0; i < cells_shared_ghosts.size(); ++i)
                    if (cells_shared_ghosts[i].size() > 0)
                      targets.push_back(i);
                  return targets;
                },
                [&](const auto other_rank, auto &send_buffer) {
                  auto &source = cells_shared_ghosts[other_rank];
                  std::sort(source.begin(),
                            source.end(),
                            [](const auto &a, const auto &b) {
                              return a.second < b.second;
                            });

                  send_buffer.reserve(source.size());
                  for (const auto &i : source)
                    send_buffer.push_back(i.second);
                },
                [&](const auto &other_rank,
                    const auto &buffer_recv,
                    auto &      request_buffer) {
                  (void)other_rank;

                  request_buffer.reserve(buffer_recv.size());

                  unsigned int j = 0;

                  for (unsigned int i = 0; i < buffer_recv.size(); ++i)
                    {
                      for (; j < cells_locally_owned.size(); ++j)
                        if (cells_locally_owned[j].second == buffer_recv[i])
                          {
                            request_buffer.push_back(
                              cells_locally_owned[j].first);
                            break;
                          }
                    }

                  AssertDimension(request_buffer.size(), buffer_recv.size());
                },
                [&](const auto other_rank, auto &recv_buffer) {
                  recv_buffer.resize(cells_shared_ghosts[other_rank].size());
                },
                [&](const auto other_rank, const auto &recv_buffer) {
                  for (unsigned int i = 0; i < recv_buffer.size(); ++i)
                    {
                      cells[cells_shared_ghosts[other_rank][i].first] = {
                        other_rank, recv_buffer[i]};
                    }
                });

          Utilities::MPI::ConsensusAlgorithms::Selector<
            dealii::types::global_dof_index,
            unsigned int>(process, communicator_sm)
            .run();

          return cells;
        }();

        // process faces (0: interior; 1: exterior)
        for (unsigned int i = 0; i < 2; ++i)
          {
            cell_indices_contiguous_sm[i].assign(
              di.dof_indices_contiguous[i].size(),
              {numbers::invalid_unsigned_int, numbers::invalid_unsigned_int});

            for (unsigned int face_index = 0;
                 face_index < cell_indices_contiguous_sm[i].size();
                 ++face_index)
              {
                const unsigned int cell_index =
                  (i == 0) ? get_face_info(face_index / n_lanes)
                               .cells_interior[face_index % n_lanes] :
                             get_face_info(face_index / n_lanes)
                               .cells_exterior[face_index % n_lanes];

                if (cell_index == numbers::invalid_unsigned_int)
                  continue;

                cell_indices_contiguous_sm[i][face_index] =
                  cell_indices_contiguous_sm[2][cell_index];
              }
          }

        for (auto &di : dof_info)
          di.compute_shared_memory_contiguous_indices(
            cell_indices_contiguous_sm);
      }
  }
#endif

  indices_are_initialized = true;
}



template <int dim, typename Number, typename VectorizedArrayType>
void
MatrixFree<dim, Number, VectorizedArrayType>::clear()
{
  dof_info.clear();
  mapping_info.clear();
  cell_level_index.clear();
  task_info.clear();
  dof_handlers.clear();
  face_info.clear();
  indices_are_initialized = false;
  mapping_is_initialized  = false;
}



template <int dim, typename Number, typename VectorizedArrayType>
std::size_t
MatrixFree<dim, Number, VectorizedArrayType>::memory_consumption() const
{
  std::size_t memory = MemoryConsumption::memory_consumption(dof_info);
  memory += MemoryConsumption::memory_consumption(cell_level_index);
  memory += MemoryConsumption::memory_consumption(face_info);
  memory += MemoryConsumption::memory_consumption(shape_info);
  memory += MemoryConsumption::memory_consumption(constraint_pool_data);
  memory += MemoryConsumption::memory_consumption(constraint_pool_row_index);
  memory += MemoryConsumption::memory_consumption(task_info);
  memory += sizeof(*this);
  memory += mapping_info.memory_consumption();
  return memory;
}



template <int dim, typename Number, typename VectorizedArrayType>
template <typename StreamType>
void
MatrixFree<dim, Number, VectorizedArrayType>::print_memory_consumption(
  StreamType &out) const
{
  out << "  Memory matrix-free data total: --> ";
  task_info.print_memory_statistics(out, memory_consumption());
  out << "   Memory cell index:                ";
  task_info.print_memory_statistics(
    out, MemoryConsumption::memory_consumption(cell_level_index));
  if (Utilities::MPI::sum(face_info.faces.size(), task_info.communicator) > 0)
    {
      out << "   Memory face indicators:           ";
      task_info.print_memory_statistics(
        out, MemoryConsumption::memory_consumption(face_info.faces));
    }
  for (unsigned int j = 0; j < dof_info.size(); ++j)
    {
      out << "   Memory DoFInfo component " << j << std::endl;
      dof_info[j].print_memory_consumption(out, task_info);
    }

  out << "   Memory mapping info" << std::endl;
  mapping_info.print_memory_consumption(out, task_info);

  out << "   Memory unit cell shape data:      ";
  task_info.print_memory_statistics(
    out, MemoryConsumption::memory_consumption(shape_info));
  if (task_info.scheme != internal::MatrixFreeFunctions::TaskInfo::none)
    {
      out << "   Memory task partitioning info:    ";
      task_info.print_memory_statistics(
        out, MemoryConsumption::memory_consumption(task_info));
    }
}



template <int dim, typename Number, typename VectorizedArrayType>
void
MatrixFree<dim, Number, VectorizedArrayType>::print(std::ostream &out) const
{
  // print indices local to global
  for (unsigned int no = 0; no < dof_info.size(); ++no)
    {
      out << "\n-- Index data for component " << no << " --" << std::endl;
      dof_info[no].print(constraint_pool_data, constraint_pool_row_index, out);
      out << std::endl;
    }
}



DEAL_II_NAMESPACE_CLOSE

#endif
