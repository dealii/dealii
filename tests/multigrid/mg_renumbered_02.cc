// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/meshworker/local_results.h>
#include <deal.II/meshworker/loop.h>
#include <deal.II/meshworker/output.h>

#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_component.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

#include <fstream>
#include <sstream>
#include <deal.II/base/std_cxx11/bind.h>

using namespace dealii;


template <int dim, typename number, int spacedim>
void
reinit_vector (const dealii::MGDoFHandler<dim,spacedim> &mg_dof,
               MGLevelObject<dealii::Vector<number> > &v)
{
  for (unsigned int level=v.min_level();
       level<=v.max_level(); ++level)
    {
      unsigned int n = mg_dof.n_dofs (level);
      v[level].reinit(n);
    }
}

template <int dim>
void initialize (const MGDoFHandler<dim> &dof,
                 Vector<double> &u)
{
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  for (typename MGDoFHandler<dim>::active_cell_iterator
       cell = dof.begin_active();
       cell != dof.end(); ++cell)
    {
      cell->get_dof_indices(dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          const unsigned int comp = dof.get_fe().system_to_component_index(i).first;
          u(dof_indices[i]) = comp+1;
        }
    }
}


template <int dim>
void initialize (const MGDoFHandler<dim> &dof,
                 MGLevelObject<Vector<double> > &u)
{
  unsigned int counter=0;
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  typename MGDoFHandler<dim>::cell_iterator
  cell = dof.begin(0);
  cell->get_mg_dof_indices(dof_indices);
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    u[0](dof_indices[i]) = ++counter;
}


template <int dim>
void diff (Vector<double> &diff, const MGDoFHandler<dim> &dof,
           const Vector<double> &v, const unsigned int level)
{
  diff.reinit (v);
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<unsigned int> mg_dof_indices(dofs_per_cell);
  const unsigned int n_comp = dof.get_fe().n_components();
  for (typename MGDoFHandler<dim>::cell_iterator
       cell= dof.begin(level);
       cell != dof.end(level); ++cell)
    {
      cell->get_mg_dof_indices(mg_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell/n_comp; ++i)
        {
          const unsigned int ni = n_comp* i;
          const unsigned int i0 = mg_dof_indices[ni];
          const unsigned int i1 = mg_dof_indices[ni+1];
          const unsigned int i2 = mg_dof_indices[ni+2];
          diff(i0) = 2*v(i0) - v(i1);
          diff(i1) = 3*v(i1) - 2*v(i2);
          diff(i2) = v(i2) - 3*v(i0);
        }
    }
}


template <int dim>
class OutputCreator : public Subscriptor
{
public:
  void cell(MeshWorker::DoFInfo<dim> &dinfo,
            MeshWorker::IntegrationInfo<dim> &info);

};

template <int dim>
void OutputCreator<dim>::cell(
  MeshWorker::DoFInfo<dim> &dinfo,
  MeshWorker::IntegrationInfo<dim> &info)
{
  const FEValuesBase<dim> &fe = info.fe_values();
  const std::vector<std::vector<double> > &uh = info.values[0];

  const unsigned int square_root = static_cast<unsigned int>(std::pow(fe.n_quadrature_points, 1./dim)+.5);
  for (unsigned int k1=0; k1<square_root; ++k1)
    {
      for (unsigned int k2=0; k2<square_root; ++k2)
        {
          for (unsigned int d=0; d<dim; ++d)
            dinfo.quadrature_value(k1*square_root+k2, d) = fe.quadrature_point(k1*square_root+k2)[d];
          for (unsigned int i=0; i<uh.size(); ++i)
            dinfo.quadrature_value(k1*square_root+k2, i+dim) = uh[i][k1*square_root+k2];
        }
    }
}


template <int dim>
class LaplaceProblem
{
public:
  typedef MeshWorker::IntegrationInfo<dim> CellInfo;

  LaplaceProblem (const unsigned int deg);
  void run ();

private:
  void setup_system ();
  void test_boundary ();
  void output_gpl(const MGDoFHandler<dim> &dof,
                  MGLevelObject<Vector<double> > &v);
  void refine_local ();

  Triangulation<dim>   triangulation;
  const MappingQ1<dim> mapping;
  FESystem<dim>        fe;
  MGDoFHandler<dim>    mg_dof_handler_renumbered;

  const unsigned int degree;
  std::vector<std::set<types::global_dof_index> >
  boundary_indices_renumbered;

};


template <int dim>
LaplaceProblem<dim>::LaplaceProblem (const unsigned int deg) :
  triangulation (Triangulation<dim>::limit_level_difference_at_vertices),
  fe (FE_Q<dim> (deg),3),
  mg_dof_handler_renumbered (triangulation),
  degree(deg)
{}


template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  mg_dof_handler_renumbered.distribute_dofs (fe);

  const unsigned int nlevels = triangulation.n_levels();

  DoFHandler<dim> &dof = mg_dof_handler_renumbered;
  DoFRenumbering::component_wise (dof);
  for (unsigned int level=0; level<nlevels; ++level)
    DoFRenumbering::component_wise (mg_dof_handler_renumbered, level);

  boundary_indices_renumbered.resize(nlevels);
  for (unsigned int l=0; l<nlevels; ++l)
    boundary_indices_renumbered[l].clear();

  deallog << "Number of degrees of freedom: "
          << mg_dof_handler_renumbered.n_dofs();

  for (unsigned int l=0; l<triangulation.n_levels(); ++l)
    deallog << "   " << 'L' << l << ": "
            << mg_dof_handler_renumbered.n_dofs(l);
  deallog  << std::endl;
}

template <int dim>
void
LaplaceProblem<dim>::output_gpl(const MGDoFHandler<dim> &dof,
                                MGLevelObject<Vector<double> > &v)
{
  MeshWorker::IntegrationInfoBox<dim> info_box;
  const unsigned int n_gauss_points = dof.get_fe().tensor_degree();
  QTrapez<1> trapez;
  QIterated<dim> quadrature(trapez, n_gauss_points);
  info_box.cell_quadrature = quadrature;
  info_box.initialize_update_flags();
  UpdateFlags update_flags = update_quadrature_points | update_values | update_gradients;
  info_box.add_update_flags(update_flags, true, true, true, true);

  NamedData<MGLevelObject<Vector<double> >* > data;
  data.add(&v, "mg_vector");
  info_box.cell_selector.add("mg_vector");
  info_box.initialize(fe, mapping, data);

  MeshWorker::DoFInfo<dim> dof_info(dof);

  MeshWorker::Assembler::GnuplotPatch assembler;
  assembler.initialize(n_gauss_points+1, dof.get_fe().n_components());

  for (unsigned int l=0; l<triangulation.n_levels(); ++l)
    {
      OutputCreator<dim> matrix_integrator;
      MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, MeshWorker::IntegrationInfoBox<dim> > (
        dof.begin(l), dof.end(l),
        dof_info, info_box,
        std_cxx11::bind(&OutputCreator<dim>::cell, &matrix_integrator, std_cxx11::_1, std_cxx11::_2),
        0,
        0,
        assembler);
    }
}

template <int dim>
void LaplaceProblem<dim>::test_boundary ()
{
  typename FunctionMap<dim>::type      dirichlet_boundary;
  ZeroFunction<dim>                    dirichlet_bc(fe.n_components());
  dirichlet_boundary[0] =             &dirichlet_bc;
  MGTools::make_boundary_list (mg_dof_handler_renumbered, dirichlet_boundary,
                               boundary_indices_renumbered);

  MGLevelObject<Vector<double> > u(0, triangulation.n_levels()-1);
  MGLevelObject<Vector<double> > d(0, triangulation.n_levels()-1);

  for (unsigned int l=0; l<triangulation.n_levels(); ++l)
    {
      u[l].reinit(mg_dof_handler_renumbered.n_dofs(l));
      u[l] = 1.;

      for (std::set<types::global_dof_index>::const_iterator i = boundary_indices_renumbered[l].begin();
           i!=boundary_indices_renumbered[l].end(); ++i)
        {
          u[l](*i) = 0;
        }
    }
  output_gpl(mg_dof_handler_renumbered, u);
}



template <int dim>
void LaplaceProblem<dim>::refine_local ()
{
  bool cell_refined = false;
  for (typename Triangulation<dim>::active_cell_iterator
       cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    {
      for (unsigned int vertex=0;
           vertex < GeometryInfo<dim>::vertices_per_cell;
           ++vertex)
        {
          const Point<dim> p = cell->vertex(vertex);
          const Point<dim> origin = (dim == 2 ?
                                     Point<dim>(0,0) :
                                     Point<dim>(0,0,0));
          const double dist = p.distance(origin);
          if (dist<0.25/numbers::PI)
            {
              cell->set_refine_flag ();
              cell_refined = true;
              break;
            }
        }
    }
  //Wenn nichts verfeinert wurde bisher, global verfeinern!
  if (!cell_refined)
    for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      cell->set_refine_flag();


  triangulation.execute_coarsening_and_refinement ();
}

template <int dim>
void LaplaceProblem<dim>::run ()
{
  for (unsigned int cycle=0; cycle<6; ++cycle)
    {
      deallog << "Cycle " << cycle << std::endl;

      if (cycle == 0)
        {
          GridGenerator::hyper_cube(triangulation, -1, 1);
          triangulation.refine_global (1);
        }
      //triangulation.refine_global (1);
      refine_local ();
      setup_system ();
      test_boundary();
    };
}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  LaplaceProblem<2> laplace_problem_2d(1);
  laplace_problem_2d.run ();
}
