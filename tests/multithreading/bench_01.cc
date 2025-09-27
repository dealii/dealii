// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test several WorkStream implementations (TBB, taskflow, sequential) and
// allow timing benchmarks (run this manually to see std output and see below
// for settings to change).

#include <deal.II/base/multithread_info.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

#ifdef DEAL_II_WITH_TASKFLOW
#  include "taskflow.h"
#endif

// number of global refinements:
const unsigned int n_refinements = 3;
// how many times to run each benchmark before averaging
const unsigned int n_runs = 1;
// maximum number of threads to test (to speed up testing)
const unsigned int n_max_threads = 4;


template <int dim>
class AdvectionField : public TensorFunction<1, dim>
{
public:
  virtual Tensor<1, dim>
  value(const Point<dim> &p) const override;
};

template <int dim>
Tensor<1, dim>
AdvectionField<dim>::value(const Point<dim> &p) const
{
  Point<dim> value;
  value[0] = 2;
  for (unsigned int i = 1; i < dim; ++i)
    value[i] = 1 + 0.8 * std::sin(8. * numbers::PI * p[0]);

  return value;
}
template <int dim>
class RightHandSide : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;

private:
  static const Point<dim> center_point;
};


template <>
const Point<1> RightHandSide<1>::center_point = Point<1>(-0.75);

template <>
const Point<2> RightHandSide<2>::center_point = Point<2>(-0.75, -0.75);

template <>
const Point<3> RightHandSide<3>::center_point = Point<3>(-0.75, -0.75, -0.75);



// The only new thing here is that we check for the value of the
// <code>component</code> parameter. As this is a scalar function, it is
// obvious that it only makes sense if the desired component has the index
// zero, so we assert that this is indeed the
// case. <code>ExcIndexRange</code> is a global predefined exception
// (probably the one most often used, we therefore made it global instead of
// local to some class), that takes three parameters: the index that is
// outside the allowed range, the first element of the valid range and the
// one past the last (i.e. again the half-open interval so often used in the
// C++ standard library):
template <int dim>
double
RightHandSide<dim>::value(const Point<dim>  &p,
                          const unsigned int component) const
{
  (void)component;
  Assert(component == 0, ExcIndexRange(component, 0, 1));
  const double diameter = 0.1;
  return ((p - center_point).norm_square() < diameter * diameter ?
            0.1 / std::pow(diameter, dim) :
            0.0);
}
template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;
};


template <int dim>
double
BoundaryValues<dim>::value(const Point<dim>  &p,
                           const unsigned int component) const
{
  (void)component;
  Assert(component == 0, ExcIndexRange(component, 0, 1));

  const double sine_term = std::sin(16. * numbers::PI * p.norm_square());
  const double weight    = std::exp(5. * (1. - p.norm_square()));
  return weight * sine_term;
}

template <int dim>
struct AssemblyScratchData
{
  AssemblyScratchData(const FiniteElement<dim> &fe);
  AssemblyScratchData(const AssemblyScratchData &scratch_data);

  // FEValues and FEFaceValues are expensive objects to set up, so we
  // include them in the scratch object so that as much data is reused
  // between cells as possible.
  FEValues<dim>     fe_values;
  FEFaceValues<dim> fe_face_values;

  // We also store a few vectors that we will populate with values on each
  // cell. Setting these objects up is, in the usual case, cheap; however,
  // they require memory allocations, which can be expensive in
  // multithreaded applications. Hence we keep them here so that
  // computations on a cell do not require new allocations.
  std::vector<double>         rhs_values;
  std::vector<Tensor<1, dim>> advection_directions;
  std::vector<double>         face_boundary_values;
  std::vector<Tensor<1, dim>> face_advection_directions;

  // Finally, we need objects that describe the problem's data:
  AdvectionField<dim> advection_field;
  RightHandSide<dim>  right_hand_side;
  BoundaryValues<dim> boundary_values;
};

template <int dim>
AssemblyScratchData<dim>::AssemblyScratchData(const FiniteElement<dim> &fe)
  : fe_values(fe,
              QGauss<dim>(fe.degree + 1),
              update_values | update_gradients | update_quadrature_points |
                update_JxW_values)
  , fe_face_values(fe,
                   QGauss<dim - 1>(fe.degree + 1),
                   update_values | update_quadrature_points |
                     update_JxW_values | update_normal_vectors)
  , rhs_values(fe_values.get_quadrature().size())
  , advection_directions(fe_values.get_quadrature().size())
  , face_boundary_values(fe_face_values.get_quadrature().size())
  , face_advection_directions(fe_face_values.get_quadrature().size())
{}



template <int dim>
AssemblyScratchData<dim>::AssemblyScratchData(
  const AssemblyScratchData &scratch_data)
  : fe_values(scratch_data.fe_values.get_fe(),
              scratch_data.fe_values.get_quadrature(),
              update_values | update_gradients | update_quadrature_points |
                update_JxW_values)
  , fe_face_values(scratch_data.fe_face_values.get_fe(),
                   scratch_data.fe_face_values.get_quadrature(),
                   update_values | update_quadrature_points |
                     update_JxW_values | update_normal_vectors)
  , rhs_values(scratch_data.rhs_values.size())
  , advection_directions(scratch_data.advection_directions.size())
  , face_boundary_values(scratch_data.face_boundary_values.size())
  , face_advection_directions(scratch_data.face_advection_directions.size())
{}

struct AssemblyCopyData
{
  FullMatrix<double>                   cell_matrix;
  Vector<double>                       cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
};

template <int dim>
void
local_assemble_system(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  AssemblyScratchData<dim>                             &scratch_data,
  AssemblyCopyData                                     &copy_data)
{
  // We define some abbreviations to avoid unnecessarily long lines:
  const unsigned int dofs_per_cell =
    scratch_data.fe_values.get_fe().dofs_per_cell;
  const unsigned int n_q_points =
    scratch_data.fe_values.get_quadrature().size();
  const unsigned int n_face_q_points =
    scratch_data.fe_face_values.get_quadrature().size();

  // We declare cell matrix and cell right hand side...
  copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
  copy_data.cell_rhs.reinit(dofs_per_cell);

  // ... an array to hold the global indices of the degrees of freedom of
  // the cell on which we are presently working...
  copy_data.local_dof_indices.resize(dofs_per_cell);

  // ... then initialize the <code>FEValues</code> object...
  scratch_data.fe_values.reinit(cell);

  // ... obtain the values of right hand side and advection directions
  // at the quadrature points...
  scratch_data.advection_field.value_list(
    scratch_data.fe_values.get_quadrature_points(),
    scratch_data.advection_directions);
  scratch_data.right_hand_side.value_list(
    scratch_data.fe_values.get_quadrature_points(), scratch_data.rhs_values);

  // ... set the value of the streamline diffusion parameter as
  // described in the introduction...
  const double delta = 0.1 * cell->diameter();

  // ... and assemble the local contributions to the system matrix and
  // right hand side as also discussed above:
  for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        // Alias the AssemblyScratchData object to keep the lines from
        // getting too long:
        const auto &sd = scratch_data;
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          copy_data.cell_matrix(i, j) +=
            ((sd.fe_values.shape_value(i, q_point) +           // (phi_i +
              delta * (sd.advection_directions[q_point] *      // delta beta
                       sd.fe_values.shape_grad(i, q_point))) * // grad phi_i)
             sd.advection_directions[q_point] *                // beta
             sd.fe_values.shape_grad(j, q_point)) *            // grad phi_j
            sd.fe_values.JxW(q_point);                         // dx

        copy_data.cell_rhs(i) +=
          (sd.fe_values.shape_value(i, q_point) +           // (phi_i +
           delta * (sd.advection_directions[q_point] *      // delta beta
                    sd.fe_values.shape_grad(i, q_point))) * // grad phi_i)
          sd.rhs_values[q_point] *                          // f
          sd.fe_values.JxW(q_point);                        // dx
      }

  // Besides the cell terms which we have built up now, the bilinear
  // form of the present problem also contains terms on the boundary of
  // the domain. Therefore, we have to check whether any of the faces of
  // this cell are on the boundary of the domain, and if so assemble the
  // contributions of this face as well. Of course, the bilinear form
  // only contains contributions from the <code>inflow</code> part of
  // the boundary, but to find out whether a certain part of a face of
  // the present cell is part of the inflow boundary, we have to have
  // information on the exact location of the quadrature points and on
  // the direction of flow at this point; we obtain this information
  // using the FEFaceValues object and only decide within the main loop
  // whether a quadrature point is on the inflow boundary.
  for (const auto &face : cell->face_iterators())
    if (face->at_boundary())
      {
        // Ok, this face of the present cell is on the boundary of the
        // domain. Just as for the usual FEValues object which we have
        // used in previous examples and also above, we have to
        // reinitialize the FEFaceValues object for the present face:
        scratch_data.fe_face_values.reinit(cell, face);

        // For the quadrature points at hand, we ask for the values of
        // the inflow function and for the direction of flow:
        scratch_data.boundary_values.value_list(
          scratch_data.fe_face_values.get_quadrature_points(),
          scratch_data.face_boundary_values);
        scratch_data.advection_field.value_list(
          scratch_data.fe_face_values.get_quadrature_points(),
          scratch_data.face_advection_directions);

        // Now loop over all quadrature points and see whether this face is on
        // the inflow or outflow part of the boundary. The normal
        // vector points out of the cell: since the face is at
        // the boundary, the normal vector points out of the domain,
        // so if the advection direction points into the domain, its
        // scalar product with the normal vector must be negative (to see why
        // this is true, consider the scalar product definition that uses a
        // cosine):
        for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point)
          if (scratch_data.fe_face_values.normal_vector(q_point) *
                scratch_data.face_advection_directions[q_point] <
              0.)
            // If the face is part of the inflow boundary, then compute the
            // contributions of this face to the global matrix and right
            // hand side, using the values obtained from the
            // FEFaceValues object and the formulae discussed in the
            // introduction:
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  copy_data.cell_matrix(i, j) -=
                    (scratch_data.face_advection_directions[q_point] *
                     scratch_data.fe_face_values.normal_vector(q_point) *
                     scratch_data.fe_face_values.shape_value(i, q_point) *
                     scratch_data.fe_face_values.shape_value(j, q_point) *
                     scratch_data.fe_face_values.JxW(q_point));

                copy_data.cell_rhs(i) -=
                  (scratch_data.face_advection_directions[q_point] *
                   scratch_data.fe_face_values.normal_vector(q_point) *
                   scratch_data.face_boundary_values[q_point] *
                   scratch_data.fe_face_values.shape_value(i, q_point) *
                   scratch_data.fe_face_values.JxW(q_point));
              }
      }

  // The final piece of information the copy routine needs is the global
  // indices of the degrees of freedom on this cell, so we end by writing
  // them to the local array:
  cell->get_dof_indices(copy_data.local_dof_indices);
}



template <int dim>
void
assemble()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(n_refinements);
  std::cout << "   Number of active cells:              "
            << triangulation.n_active_cells() << std::endl;

  DoFHandler<dim> dof_handler(triangulation);
  FE_Q<dim>       fe(5);


  AffineConstraints<double> hanging_node_constraints;
  hanging_node_constraints.close();
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;
  Vector<double>       system_rhs;
  {
    dof_handler.distribute_dofs(fe);

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    hanging_node_constraints,
                                    /*keep_constrained_dofs =*/false);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
    system_rhs.reinit(dof_handler.n_dofs());
    std::cout << "   Number of degrees of freedom:        "
              << dof_handler.n_dofs() << std::endl;
  }

  const unsigned int n_phys_cores =
    std::min(n_max_threads, MultithreadInfo::n_cores());
  std::cout << "MultithreadInfo::n_cores()=" << n_phys_cores
            << " (limited from " << MultithreadInfo::n_cores() << ')'
            << std::endl;
  std::cout << "MultithreadInfo::n_threads()=" << MultithreadInfo::n_threads()
            << std::endl;



  Timer timer;

  auto worker = [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
                    AssemblyScratchData<dim> &scratch_data,
                    AssemblyCopyData         &copy_data) {
    local_assemble_system(cell, scratch_data, copy_data);
  };

  auto copier = [&](const AssemblyCopyData &copy_data) {
    hanging_node_constraints.distribute_local_to_global(
      copy_data.cell_matrix,
      copy_data.cell_rhs,
      copy_data.local_dof_indices,
      system_matrix,
      system_rhs);
  };


  double reference_l2 = -1.;
  {
    // reference sequential implementation

    double avg = 0.;

    for (unsigned int c = 0; c < n_runs; ++c)
      {
        system_rhs = 0.;
        timer.reset();
        timer.start();
        WorkStream::internal::sequential::run(dof_handler.begin_active(),
                                              dof_handler.end(),
                                              worker,
                                              copier,
                                              AssemblyScratchData<dim>(fe),
                                              AssemblyCopyData());
        timer.stop();
        const double time = timer.last_wall_time();
        avg += time;
        std::cout << time << ' ' << std::flush;
        reference_l2 = system_rhs.l2_norm();
      }
    avg /= n_runs;
    std::cout << " avg: " << avg << std::endl;
  }

  using Iterator = typename DoFHandler<dim>::active_cell_iterator;
  std::vector<std::vector<Iterator>> graph;

  {
    // make graph coloring
    timer.reset();
    timer.start();

    graph = GraphColoring::make_graph_coloring(
      dof_handler.begin_active(),
      dof_handler.end(),
      [&](const Iterator &cell) -> std::vector<types::global_dof_index> {
        std::vector<types::global_dof_index> local_dof_indices(
          fe.dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);
        return local_dof_indices;
      });
    timer.stop();
    const double time = timer.last_wall_time();
    std::cout << "build graph coloring time: " << time << std::endl;
  }


#ifdef DEAL_II_WITH_TASKFLOW
  {
    std::cout << "** TASKFLOW v1 **" << std::endl;

    for (unsigned int n_cores = n_phys_cores; n_cores > 0; n_cores /= 2)
      {
        MultithreadInfo::set_thread_limit(n_cores);

        std::cout << "n_cores " << n_cores;
        std::cout << ' ' << std::flush;
        double avg = 0.;

        for (unsigned int c = 0; c < n_runs; ++c)
          {
            system_rhs = 0.;
            timer.reset();
            timer.start();

            taskflow_v1::run(dof_handler.begin_active(),
                             dof_handler.end(),
                             worker,
                             copier,
                             AssemblyScratchData<dim>(fe),
                             AssemblyCopyData());

            timer.stop();
            const double time = timer.last_wall_time();
            avg += time;
            std::cout << time << ' ' << std::flush;
            Assert(abs(reference_l2 - system_rhs.l2_norm()) < 1e-10,
                   ExcInternalError());
          }
        avg /= n_runs;
        std::cout << " avg: " << avg << std::endl;
      }
  }
#endif

#ifdef DEAL_II_WITH_TBB
  {
    std::cout << "** TBB **" << std::endl;

    for (unsigned int n_cores = n_phys_cores; n_cores > 0; n_cores /= 2)
      {
        MultithreadInfo::set_thread_limit(n_cores);

        std::cout << "n_cores " << n_cores;
        std::cout << ' ' << std::flush;
        double avg = 0.;

        for (unsigned int c = 0; c < n_runs; ++c)
          {
            system_rhs = 0.;
            timer.reset();
            timer.start();
            WorkStream::internal::tbb_no_coloring::run(
              dof_handler.begin_active(),
              dof_handler.end(),
              worker,
              copier,
              AssemblyScratchData<dim>(fe),
              AssemblyCopyData(),
              2 * MultithreadInfo::n_threads(),
              8);


            timer.stop();
            const double time = timer.last_wall_time();
            avg += time;
            std::cout << time << ' ' << std::flush;
            Assert(abs(reference_l2 - system_rhs.l2_norm()) < 1e-10,
                   ExcInternalError());
          }
        avg /= n_runs;
        std::cout << " avg: " << avg << std::endl;
      }
  }

  {
    std::cout << "** TBB colored **" << std::endl;
    for (unsigned int n_cores = n_phys_cores; n_cores > 0; n_cores /= 2)
      {
        MultithreadInfo::set_thread_limit(n_cores);

        std::cout << "n_cores " << n_cores;
        std::cout << ' ' << std::flush;
        double avg = 0.;

        for (unsigned int c = 0; c < n_runs; ++c)
          {
            system_rhs = 0.;
            timer.reset();
            timer.start();
            WorkStream::internal::tbb_colored::run(graph,
                                                   worker,
                                                   copier,
                                                   AssemblyScratchData<dim>(fe),
                                                   AssemblyCopyData(),
                                                   8);

            timer.stop();
            const double time = timer.last_wall_time();
            avg += time;
            std::cout << time << ' ' << std::flush;
            Assert(abs(reference_l2 - system_rhs.l2_norm()) < 1e-10,
                   ExcInternalError());
          }
        avg /= n_runs;
        std::cout << " avg: " << avg << std::endl;
      }
  }

#endif

  {
    std::cout << "** WorkStream **" << std::endl;
    for (unsigned int n_cores = n_phys_cores; n_cores > 0; n_cores /= 2)
      {
        MultithreadInfo::set_thread_limit(n_cores);

        std::cout << "n_cores " << n_cores;
        std::cout << ' ' << std::flush;
        double avg = 0.;

        for (unsigned int c = 0; c < n_runs; ++c)
          {
            system_rhs = 0.;
            timer.reset();
            timer.start();
            WorkStream::run(dof_handler.begin_active(),
                            dof_handler.end(),
                            worker,
                            copier,
                            AssemblyScratchData<dim>(fe),
                            AssemblyCopyData());

            timer.stop();
            const double time = timer.last_wall_time();
            avg += time;
            std::cout << time << ' ' << std::flush;
            Assert(abs(reference_l2 - system_rhs.l2_norm()) < 1e-10,
                   ExcInternalError());
          }
        avg /= n_runs;
        std::cout << " avg: " << avg << std::endl;
      }
  }

  {
    std::cout << "** WorkStream colored **" << std::endl;
    for (unsigned int n_cores = n_phys_cores; n_cores > 0; n_cores /= 2)
      {
        MultithreadInfo::set_thread_limit(n_cores);

        std::cout << "n_cores " << n_cores;
        std::cout << ' ' << std::flush;
        double avg = 0.;

        for (unsigned int c = 0; c < n_runs; ++c)
          {
            system_rhs = 0.;
            timer.reset();
            timer.start();
            WorkStream::run(graph,
                            worker,
                            copier,
                            AssemblyScratchData<dim>(fe),
                            AssemblyCopyData());

            timer.stop();
            const double time = timer.last_wall_time();
            avg += time;
            std::cout << time << ' ' << std::flush;
            Assert(abs(reference_l2 - system_rhs.l2_norm()) < 1e-10,
                   ExcInternalError());
          }
        avg /= n_runs;
        std::cout << " avg: " << avg << std::endl;
      }
  }
}



int
main()
{
  initlog();
  // remove the limit normally applied to all tests:
  MultithreadInfo::set_thread_limit();

  assemble<2>();
  deallog << "ok" << std::endl;
}
