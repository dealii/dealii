// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Select different components from three vectors on the same DoFHandler.

#include <deal.II/base/function.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_faces.h>
#include <deal.II/numerics/data_out_resample.h>
#include <deal.II/numerics/data_out_rotation.h>
#include <deal.II/numerics/vector_tools.h>

#include <iomanip>

#include "../tests.h"

using namespace dealii::DataComponentInterpretation;
using Interpretation =
  dealii::DataComponentInterpretation::DataComponentInterpretation;

// Expose the metadata as well as the patches for independent checks.
template <typename Base>
struct Output : Base
{
  using Base::Base;
  using Base::get_dataset_names;
  using Base::get_nonscalar_data_ranges;
  using Base::get_patches;
};


// Record actual metadata and sampled values, not just a success marker.
// Values are constant per component; assertions below check every sample.
template <typename OutputType>
void
log_output(const OutputType &output)
{
  const auto names = output.get_dataset_names();
  deallog << "Fields: " << names.size() << std::endl;
  for (unsigned int c = 0; c < names.size(); ++c)
    deallog << names[c] << " = " << output.get_patches().front().data(c, 0)
            << std::endl;
  const auto ranges = output.get_nonscalar_data_ranges();
  deallog << "Nonscalar ranges: " << ranges.size() << std::endl;
  for (const auto &range : ranges)
    deallog << std::get<0>(range) << ' ' << std::get<1>(range) << ' '
            << std::get<2>(range) << std::endl;
}


template <typename Number>
Number
value(const unsigned int component, const unsigned int vector)
{
  if constexpr (numbers::NumberTraits<Number>::is_complex)
    return Number(10 * vector + component + 1., -10. * vector - component - 1.);
  else
    return 10 * vector + component + 1;
}


template <int dim, typename Number, typename OutputType>
void
check_output(OutputType                        &output,
             const DoFHandler<dim>             &dh,
             const std::vector<Vector<Number>> &vectors)
{
  std::vector<Interpretation> velocity(dim + 1, component_is_part_of_vector),
    pressure(dim + 1, component_is_ignored);
  velocity.back() = component_is_ignored;
  pressure.back() = component_is_scalar;
  std::vector<std::string> names(dim + 1, "velocity_a");
  names.back() = "ignored_pressure";
  output.add_data_vector(dh, vectors[0], names, velocity);
  std::fill(names.begin(), names.end(), "velocity_b");
  names.back() = "ignored_pressure";
  output.add_data_vector(dh, vectors[1], names, velocity);
  std::fill(names.begin(), names.end(), "ignored_velocity");
  names.back() = "pressure";
  output.add_data_vector(dh, vectors[2], names, pressure);

  constexpr bool           complex = numbers::NumberTraits<Number>::is_complex;
  std::vector<double>      expected;
  std::vector<std::string> expected_names;
  for (unsigned int v = 0; v < 3; ++v)
    for (unsigned int part = 0; part < (complex ? 2 : 1); ++part)
      for (unsigned int c = (v == 2 ? dim : 0); c < (v == 2 ? dim + 1 : dim);
           ++c)
        {
          expected.push_back((part == 0 ? 1 : -1) * (10 * v + c + 1.));
          expected_names.push_back(
            (v == 0 ? "velocity_a" :
             v == 1 ? "velocity_b" :
                      "pressure") +
            std::string(complex ? (part == 0 ? "_re" : "_im") : ""));
        }

  output.build_patches(2);
  AssertThrow(output.get_dataset_names() == expected_names, ExcInternalError());
  const auto ranges = output.get_nonscalar_data_ranges();
  AssertDimension(ranges.size(), complex ? 4 : 2);
  for (unsigned int r = 0; r < ranges.size(); ++r)
    {
      AssertDimension(std::get<0>(ranges[r]), r * dim);
      AssertDimension(std::get<1>(ranges[r]), (r + 1) * dim - 1);
      AssertThrow(std::get<2>(ranges[r]) == expected_names[r * dim],
                  ExcInternalError());
    }
  AssertThrow(!output.get_patches().empty(), ExcInternalError());
  for (const auto &patch : output.get_patches())
    {
      AssertDimension(patch.data.n_rows(),
                      expected.size() + (patch.points_are_available ? dim : 0));
      for (unsigned int row = 0; row < expected.size(); ++row)
        for (unsigned int q = 0; q < patch.data.n_cols(); ++q)
          AssertThrow(std::abs(patch.data(row, q) - expected[row]) < 1.e-12,
                      ExcInternalError());
    }
  if constexpr (!complex)
    log_output(output);
  std::ostringstream vtk;
  output.write_vtk(vtk);
  AssertThrow(vtk.str().find("ignored_") == std::string::npos,
              ExcInternalError());

  // Rebuilding must not filter already compacted patches a second time.
  output.build_patches(2);
  std::ostringstream rebuilt;
  output.write_vtk(rebuilt);
  AssertThrow(vtk.str() == rebuilt.str(), ExcInternalError());

  output.clear();
  output.add_data_vector(dh,
                         vectors[0],
                         std::vector<std::string>(dim + 1, "ignored_all"),
                         std::vector<Interpretation>(dim + 1,
                                                     component_is_ignored));
  output.build_patches(2);
  AssertThrow(output.get_dataset_names().empty(), ExcInternalError());
  AssertThrow(output.get_nonscalar_data_ranges().empty(), ExcInternalError());
  if constexpr (!complex)
    {
      deallog << "All ignored:" << std::endl;
      log_output(output);
    }
  std::ostringstream mesh;
  output.write_vtk(mesh);
  AssertThrow(mesh.str().find("ignored_all") == std::string::npos,
              ExcInternalError());
}


template <int dim, typename Number>
void
check()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 1., 2.);
  FESystem<dim>   fe(FE_Q<dim>(1), dim, FE_Q<dim>(1), 1);
  DoFHandler<dim> dh(tria);
  dh.distribute_dofs(fe);
  std::vector<Vector<Number>>          vectors(3, Vector<Number>(dh.n_dofs()));
  std::vector<types::global_dof_index> indices(fe.n_dofs_per_cell());
  for (const auto &cell : dh.active_cell_iterators())
    {
      cell->get_dof_indices(indices);
      for (unsigned int i = 0; i < indices.size(); ++i)
        for (unsigned int v = 0; v < vectors.size(); ++v)
          vectors[v][indices[i]] =
            value<Number>(fe.system_to_component_index(i).first, v);
    }

  Output<DataOut<dim>> output;
  if constexpr (!numbers::NumberTraits<Number>::is_complex)
    deallog << "DataOut " << dim << "d:" << std::endl;
  check_output(output, dh, vectors);
  if constexpr (!numbers::NumberTraits<Number>::is_complex)
    {
      Output<DataOutFaces<dim>> faces;
      deallog << "DataOutFaces " << dim << "d:" << std::endl;
      check_output(faces, dh, vectors);
    }
}


class Postprocessor : public DataPostprocessor<2>
{
public:
  std::vector<std::string>
  get_names() const override
  {
    return {"ignored", "processed", "ignored"};
  }

  std::vector<Interpretation>
  get_data_component_interpretation() const override
  {
    return {component_is_ignored, component_is_scalar, component_is_ignored};
  }

  UpdateFlags
  get_needed_update_flags() const override
  {
    return update_values;
  }

  void
  evaluate_vector_field(const DataPostprocessorInputs::Vector<2> &inputs,
                        std::vector<Vector<double>> &computed) const override
  {
    for (unsigned int q = 0; q < computed.size(); ++q)
      for (unsigned int c = 0; c < 3; ++c)
        computed[q][c] = 2 * inputs.solution_values[q][c];
  }
};


// Exercise ignored fields before and after a tensor, cell data, and the
// rotation output path (where patch_spacedim differs from dim).
void
check_other_fields()
{
  Triangulation<2> tria;
  GridGenerator::hyper_cube(tria, 1., 2.);
  FESystem<2>   fe(FE_Q<2>(1), 6);
  DoFHandler<2> dh(tria);
  dh.distribute_dofs(fe);
  Vector<double>                       values(dh.n_dofs());
  std::vector<types::global_dof_index> indices(fe.n_dofs_per_cell());
  dh.begin_active()->get_dof_indices(indices);
  for (unsigned int i = 0; i < indices.size(); ++i)
    values[indices[i]] = fe.system_to_component_index(i).first + 1;
  std::vector<std::string> names(6, "tensor");
  names.front() = names.back() = "ignored";
  std::vector<Interpretation> interpretation(6, component_is_part_of_tensor);
  interpretation.front() = interpretation.back() = component_is_ignored;

  Output<DataOut<2>> output;
  output.add_data_vector(dh, values, names, interpretation);
  Vector<double> cell_values(1);
  cell_values[0] = 42;
  output.add_data_vector(cell_values,
                         "ignored_cell",
                         DataOut<2>::type_cell_data,
                         {component_is_ignored});
  output.add_data_vector(cell_values, "cell", DataOut<2>::type_cell_data);
  output.build_patches(2);
  AssertDimension(output.get_dataset_names().size(), 5);
  const auto ranges = output.get_nonscalar_data_ranges();
  AssertDimension(ranges.size(), 1);
  AssertDimension(std::get<0>(ranges[0]), 0);
  AssertDimension(std::get<1>(ranges[0]), 3);
  AssertThrow(std::get<3>(ranges[0]) == component_is_part_of_tensor,
              ExcInternalError());
  for (const auto &patch : output.get_patches())
    for (unsigned int q = 0; q < patch.data.n_cols(); ++q)
      {
        for (unsigned int c = 0; c < 4; ++c)
          AssertThrow(std::abs(patch.data(c, q) - (c + 2.)) < 1.e-12,
                      ExcInternalError());
        AssertThrow(patch.data(4, q) == 42., ExcInternalError());
      }

  deallog << "Tensor and cell data:" << std::endl;
  log_output(output);

  Postprocessor      postprocessor;
  Output<DataOut<2>> processed;
  processed.add_data_vector(dh, values, postprocessor);
  processed.build_patches(2);
  AssertThrow(processed.get_dataset_names() ==
                std::vector<std::string>{"processed"},
              ExcInternalError());
  for (const auto &patch : processed.get_patches())
    {
      AssertDimension(patch.data.n_rows(), 1);
      for (unsigned int q = 0; q < patch.data.n_cols(); ++q)
        AssertThrow(std::abs(patch.data(0, q) - 4.) < 1.e-12,
                    ExcInternalError());
    }

  deallog << "Postprocessor:" << std::endl;
  log_output(processed);

  for (unsigned int c = 1; c < 5; ++c)
    {
      interpretation[c] = component_is_scalar;
      names[c]          = "scalar_" + std::to_string(c);
    }
#ifdef DEAL_II_WITH_MPI
  MappingQ1<2>                     mapping;
  Output<DataOutResample<2, 2, 2>> resampled(tria, mapping);
  resampled.add_data_vector(dh, values, names, interpretation);
  resampled.build_patches(mapping, 2);
  AssertDimension(resampled.get_dataset_names().size(), 4);
  for (const auto &patch : resampled.get_patches())
    for (unsigned int c = 0; c < 4; ++c)
      for (unsigned int q = 0; q < patch.data.n_cols(); ++q)
        AssertThrow(std::abs(patch.data(c, q) - (c + 2.)) < 1.e-12,
                    ExcInternalError());

#endif

  Output<DataOutRotation<2>> rotation;
  rotation.add_data_vector(dh, values, names, interpretation);
  rotation.build_patches(4, 2);
  AssertDimension(rotation.get_dataset_names().size(), 4);
  for (const auto &patch : rotation.get_patches())
    {
      AssertDimension(patch.data.n_rows(), 4);
      for (unsigned int c = 0; c < 4; ++c)
        for (unsigned int q = 0; q < patch.data.n_cols(); ++q)
          AssertThrow(std::abs(patch.data(c, q) - (c + 2.)) < 1.e-12,
                      ExcInternalError());
    }
  deallog << "DataOutRotation:" << std::endl;
  log_output(rotation);
}


#ifdef DEAL_II_WITH_MPI
// Use affine, component-dependent data to check values in the files against
// their written point coordinates, independently of FE and DoF numbering.
class ParallelFields : public Function<2>
{
public:
  explicit ParallelFields(const unsigned int source)
    : Function<2>(3)
    , source(source)
  {}

  double
  value(const Point<2> &point, const unsigned int component) const override
  {
    return 100 * source + 10 * (component + 1) + point[0] +
           (component + 1) * point[1];
  }

private:
  const unsigned int source;
};


std::string
read_output_file(const std::string &filename)
{
  std::ifstream input(filename);
  AssertThrow(input, ExcIO());
  return {std::istreambuf_iterator<char>(input),
          std::istreambuf_iterator<char>()};
}


std::vector<double>
read_vtu_array(const std::string &xml, const std::string &key)
{
  const auto position = xml.find(key);
  AssertThrow(position != std::string::npos, ExcMessage("Missing " + key));
  const auto begin = xml.find('>', position);
  const auto end   = xml.find("</DataArray>", begin);
  AssertThrow(begin != std::string::npos && end != std::string::npos, ExcIO());
  // For Points, skip the enclosing element's opening tag as well.
  const auto data_begin = key == "<Points>" ? xml.find('>', begin + 1) : begin;
  std::istringstream  data(xml.substr(data_begin + 1, end - data_begin - 1));
  std::vector<double> values;
  double              value;
  while (data >> value)
    values.push_back(value);
  AssertThrow(data.eof() && !values.empty(), ExcIO());
  return values;
}


void
check_parallel_output()
{
  const MPI_Comm     comm = MPI_COMM_WORLD;
  const unsigned int rank = Utilities::MPI::this_mpi_process(comm);
  AssertThrow(Utilities::MPI::n_mpi_processes(comm) == 2,
              ExcMessage("Run this test with exactly two MPI processes."));

  // Shared triangulations still distribute ownership of cells and DoFs, but
  // do not require p4est. Z-order gives a reproducible partition on two ranks.
  parallel::shared::Triangulation<2> tria(
    comm,
    Triangulation<2>::none,
    false,
    parallel::shared::Triangulation<2>::partition_zorder);
  // Use four coarse cells so sibling grouping cannot move all cells to rank 0.
  GridGenerator::subdivided_hyper_cube(tria, 2);
  FESystem<2>   fe(FE_Q<2>(1), 2, FE_Q<2>(1), 1);
  DoFHandler<2> dh(tria);
  dh.distribute_dofs(fe);
  const auto owned    = dh.locally_owned_dofs();
  const auto relevant = DoFTools::extract_locally_relevant_dofs(dh);
  std::vector<LinearAlgebra::distributed::Vector<double>> vectors(3);
  for (unsigned int v = 0; v < vectors.size(); ++v)
    {
      LinearAlgebra::distributed::Vector<double> values(owned, relevant, comm);
      VectorTools::interpolate(dh, ParallelFields(v), values);
      vectors[v].reinit(owned, relevant, comm);
      vectors[v] = values;
      vectors[v].update_ghost_values();
      // Ownership of interface DoFs may leave one rank without ghosts;
      // require the distributed field as a whole to exercise ghost exchange.
      AssertThrow(Utilities::MPI::sum(
                    vectors[v].get_partitioner()->n_ghost_indices(), comm) > 0,
                  ExcInternalError());
    }

  Output<DataOut<2>>                output;
  const std::vector<Interpretation> velocity = {component_is_part_of_vector,
                                                component_is_part_of_vector,
                                                component_is_ignored};
  const std::vector<Interpretation> pressure = {component_is_ignored,
                                                component_is_ignored,
                                                component_is_scalar};
  output.add_data_vector(dh,
                         vectors[0],
                         {"velocity_a", "velocity_a", "ignored_pressure"},
                         velocity);
  output.add_data_vector(dh,
                         vectors[1],
                         {"velocity_b", "velocity_b", "ignored_pressure"},
                         velocity);
  output.add_data_vector(dh,
                         vectors[2],
                         {"ignored_velocity", "ignored_velocity", "pressure"},
                         pressure);
  output.build_patches();
  AssertThrow(output.get_patches().size() == 2, ExcInternalError());

  DataOutBase::VtkFlags flags;
  flags.print_date_and_time = false;
  flags.compression_level   = DataOutBase::CompressionLevel::plain_text;
  output.set_flags(flags);
  const auto pvtu_name =
    output.write_vtu_with_pvtu_record("", "selected", 0, comm, 1);
  // The writer does not synchronize the completion of all piece files.
  MPI_Barrier(comm);

  if (rank == 0)
    {
      const auto pvtu = read_output_file(pvtu_name);
      AssertThrow(pvtu.find("ignored_") == std::string::npos,
                  ExcInternalError());
      // Ignore the timestamp comment, but record the complete PVTU schema.
      deallog << pvtu.substr(pvtu.find("<VTKFile")) << std::endl;
      for (unsigned int piece = 0; piece < 2; ++piece)
        {
          const std::string filename =
            "selected_0." + std::to_string(piece) + ".vtu";
          AssertThrow(pvtu.find("Source=\"" + filename + "\"") !=
                        std::string::npos,
                      ExcInternalError());
          const auto xml = read_output_file(filename);
          AssertThrow(xml.find("ignored_") == std::string::npos,
                      ExcInternalError());
          const auto points = read_vtu_array(xml, "<Points>");
          const auto a      = read_vtu_array(xml, "Name=\"velocity_a\"");
          const auto b      = read_vtu_array(xml, "Name=\"velocity_b\"");
          const auto p      = read_vtu_array(xml, "Name=\"pressure\"");
          AssertThrow(points.size() == 24 && a.size() == points.size() &&
                        b.size() == points.size() &&
                        p.size() * 3 == points.size(),
                      ExcInternalError());
          for (unsigned int q = 0; q < p.size(); ++q)
            {
              const Point<2> point(points[3 * q], points[3 * q + 1]);
              for (unsigned int c = 0; c < 2; ++c)
                {
                  AssertThrow(std::abs(a[3 * q + c] -
                                       ParallelFields(0).value(point, c)) <
                                1.e-5,
                              ExcInternalError());
                  AssertThrow(std::abs(b[3 * q + c] -
                                       ParallelFields(1).value(point, c)) <
                                1.e-5,
                              ExcInternalError());
                }
              AssertThrow(a[3 * q + 2] == 0. && b[3 * q + 2] == 0. &&
                            points[3 * q + 2] == 0.,
                          ExcInternalError());
              AssertThrow(std::abs(p[q] - ParallelFields(2).value(point, 2)) <
                            1.e-5,
                          ExcInternalError());
            }
          deallog << filename << ": " << p.size()
                  << " points, all field values verified" << std::endl;
          deallog << "First values: " << a[0] << ' ' << a[1] << ' ' << b[0]
                  << ' ' << b[1] << ' ' << p[0] << std::endl;
        }
    }
}
#endif


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);
#ifdef DEAL_II_WITH_MPI
  if (Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) > 1)
    {
      mpi_initlog();
      check_parallel_output();
      return 0;
    }
#endif
  initlog();
  deallog << std::fixed << std::setprecision(0);
  check<2, double>();
  check<3, double>();
  check_other_fields();
#ifdef DEAL_II_WITH_COMPLEX_VALUES
  check<2, std::complex<double>>();
  check<3, std::complex<double>>();
#endif
  deallog << "OK" << std::endl;
}
