// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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


// Test that we can use MeshWorker also in 1d. test by Scott Miller

#include "../tests.h"
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/base/std_cxx11/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iomanip>

using namespace dealii;


//! Solve the advection equation:  \dot{u} + \div{\mathbf{c} u} = 0.0,
//! with c={1}, {1,0}, {1,0,0} in d=1,2,3

//! Domain x \in [0,1], y,z \in [0,0.01]
//! Initial condition:  u = 0.0
//! Boundary condition at x=0:  u=1

//! Use a DG formulation with upwind fluxes


namespace Advection
{
  using namespace dealii;

  /********************************************
   * ADVECTION PROBLEM
   ********************************************/
  template <int dim>
  class AdvectionProblem
  {
  public:
    AdvectionProblem ();
    ~AdvectionProblem ();
    void run ();

  private:

    const MappingQ1<dim> mapping;

    void setup_system ();

    void integrate_cell_term (MeshWorker::DoFInfo<dim> &dinfo,
                              MeshWorker::IntegrationInfo<dim> &info);

    void integrate_boundary_term (MeshWorker::DoFInfo<dim> &dinfo,
                                  MeshWorker::IntegrationInfo<dim> &info);

    void integrate_face_term (MeshWorker::DoFInfo<dim> &dinfo1,
                              MeshWorker::DoFInfo<dim> &dinfo2,
                              MeshWorker::IntegrationInfo<dim> &info1,
                              MeshWorker::IntegrationInfo<dim> &info2);

    void output_results (int timestep) const;

    void create_grid ();

    void assemble_rhs (Vector<double> &solution,
                       Vector<double> &residual);

    // For problems with non-diagonal mass matrices
//    void assemble_mass_matrix_and_multiply (Vector<double>& solution,
//                                          Vector<double>& residual);

    const double wavespeed;


    // DATA:
    Triangulation<dim>   triangulation;
    DoFHandler<dim>      dof_handler;

    FESystem<dim>        fe;

    Vector<double>       solution, stage;

    const FEValuesExtractors::Scalar upos;
  };


  template <int dim>
  AdvectionProblem<dim>::AdvectionProblem ()
    :
    mapping(),
    wavespeed(1.0),
    dof_handler (triangulation),
    fe (FE_DGQ<dim>(0), 1),// p=0, and solving for a scalar
    upos(0)
  {}

  template <int dim>
  AdvectionProblem<dim>::~AdvectionProblem ()
  {
    dof_handler.clear ();
  }


  template < >
  void AdvectionProblem<1>::create_grid()
  {
    double ll_x=0.;
    double ur_x=1.;

    int n_cells_x = 10;

    const Point<1> LowerLeft (ll_x),
          UpperRight (ur_x);

    // Define the subdivisions in the x1 and x2 coordinates.
    std::vector<unsigned int> subdivisions(1);
    subdivisions[0] =   n_cells_x;

    GridGenerator::subdivided_hyper_rectangle(triangulation,
                                              subdivisions,
                                              LowerLeft,
                                              UpperRight,
                                              true);

  }//create_grid()

  template < >
  void AdvectionProblem<2>::create_grid()
  {
    // dim==1 implemented:
    double ll_x=0., ll_y=0.;
    double ur_x=1., ur_y=0.01;

    int n_cells_x = 100;
    int n_cells_y = 1;

    const Point<2> LowerLeft (ll_x, ll_y),
          UpperRight (ur_x, ur_y );

    // Define the subdivisions in the x1 and x2 coordinates.
    std::vector<unsigned int> subdivisions(2);
    subdivisions[0] =   n_cells_x;
    subdivisions[1] =   n_cells_y;

    GridGenerator::subdivided_hyper_rectangle(triangulation,
                                              subdivisions,
                                              LowerLeft,
                                              UpperRight,
                                              true);

  }//create_grid()

  template <>
  void AdvectionProblem<3>::create_grid()
  {
    double ll_x=0., ll_y=0., ll_z=0.;
    double ur_x=1., ur_y=0.01, ur_z=0.01;

    int n_cells_x = 100;
    int n_cells_y = 1;
    int n_cells_z = 1;

    const Point<3> LowerLeft (ll_x, ll_y, ll_z),
          UpperRight (ur_x, ur_y, ur_z);

    // Define the subdivisions in the x1 and x2 coordinates.
    std::vector<unsigned int> subdivisions(3);
    subdivisions[0] =   n_cells_x;
    subdivisions[1] =   n_cells_y;
    subdivisions[2] =   n_cells_z;

    GridGenerator::subdivided_hyper_rectangle(triangulation,
                                              subdivisions,
                                              LowerLeft,
                                              UpperRight,
                                              true);

  }//create_grid()

  template <int dim>
  void AdvectionProblem<dim>::setup_system ()
  {
    dof_handler.distribute_dofs (fe);
    solution.reinit(dof_handler.n_dofs());
    stage.reinit(dof_handler.n_dofs());

    solution = 0.0;
    stage = 0.0;
  }//setup_system


  template <int dim>
  void AdvectionProblem<dim>::assemble_rhs (Vector<double> &solution,
                                            Vector<double> &residual)
  {
    const unsigned int n_gauss_points = std::ceil(((2.0*fe.degree) +1)/2);

    MeshWorker::IntegrationInfoBox<dim> info_box;

    info_box.initialize_gauss_quadrature(n_gauss_points,
                                         n_gauss_points,
                                         n_gauss_points);

    info_box.initialize_update_flags();
    UpdateFlags update_flags = update_quadrature_points |
                               update_values |
                               update_gradients;

    info_box.add_update_flags(update_flags, true, true, true, true);

    NamedData<Vector<double>* > solution_data;

    Vector<double> *u = &solution;

    solution_data.add(u, "solution");
    info_box.cell_selector.add("solution", true, true, false);
    info_box.boundary_selector.add("solution", true, false, false);
    info_box.face_selector.add("solution", true, false, false);

    info_box.initialize(fe, mapping, solution_data);

//deallog<<"\nWe are now going to attend construction of  MeshWorker::DoFInfo..."<<std::endl;
    MeshWorker::DoFInfo<dim> dof_info(dof_handler);
//deallog<<"\nApparently it DoFInfo was constructed fine!"<<std::endl;

    MeshWorker::Assembler::ResidualSimple<Vector<double> > assembler;
    NamedData<Vector<double>* > data;
    Vector<double> *rhs = &residual;
    data.add(rhs, "Residual");
    assembler.initialize(data);

    MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>, MeshWorker::IntegrationInfoBox<dim> >
    (dof_handler.begin_active(), dof_handler.end(),
     dof_info, info_box,
     std_cxx11::bind(&AdvectionProblem<dim>::integrate_cell_term,
                     this, std_cxx11::_1, std_cxx11::_2),
     std_cxx11::bind(&AdvectionProblem<dim>::integrate_boundary_term,
                     this, std_cxx11::_1, std_cxx11::_2),
     std_cxx11::bind(&AdvectionProblem<dim>::integrate_face_term,
                     this, std_cxx11::_1, std_cxx11::_2, std_cxx11::_3, std_cxx11::_4),
     assembler, true);

  }//assemble_system

  template <int dim>
  void AdvectionProblem<dim>::integrate_cell_term (MeshWorker::DoFInfo<dim> &dinfo,
                                                   MeshWorker::IntegrationInfo<dim> &info)
  {
    const FEValuesBase<dim> &fe_v = info.fe_values();

    Vector<double> &cell_rhs = dinfo.vector(0).block(0);

    const unsigned int   dofs_per_cell = fe_v.dofs_per_cell;
    const unsigned int   n_q_points    = fe_v.n_quadrature_points;

    FullMatrix<double> cell_matrix(dofs_per_cell,dofs_per_cell);
    cell_matrix = 0.0;

    const std::vector<std::vector<double> > &values = info.values[0];

    std::vector<double> u(values[0]);

    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {

            cell_rhs(i) -= wavespeed*(u[q_point]*fe_v[upos].gradient(i,q_point)[0])*fe_v.JxW(q_point);

            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                cell_matrix(i,j) += fe_v[upos].value(i,q_point) *
                                    fe_v[upos].value(j,q_point) *
                                    fe_v.JxW(q_point);
              }
          }//i
      }//q_point

  }//integrate_cell_term

  template <int dim>
  void AdvectionProblem<dim>::integrate_boundary_term (MeshWorker::DoFInfo<dim> &dinfo,
                                                       MeshWorker::IntegrationInfo<dim> &info)
  {
    const unsigned int boundary_id = dinfo.face->boundary_indicator();

    // We only have a non-zero boundary contribution at the
    // x=0 boundary
    if (boundary_id != 0)
      return;

    const FEValuesBase<dim> &fe_v = info.fe_values();

    Vector<double> &cell_rhs = dinfo.vector(0).block(0);

    const unsigned int   dofs_per_cell = fe_v.dofs_per_cell;
    const unsigned int   n_q_points    = fe_v.n_quadrature_points;

    double boundary_flux = -1.0;

    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {

            cell_rhs(i) += wavespeed*boundary_flux*
                           fe_v[upos].value(i,q_point)*fe_v.JxW(q_point);

          }//i
      }//q_point

  }//integrate_boundary_term

  template <int dim>
  void AdvectionProblem<dim>::integrate_face_term (MeshWorker::DoFInfo<dim> &dinfo1,
                                                   MeshWorker::DoFInfo<dim> &dinfo2,
                                                   MeshWorker::IntegrationInfo<dim> &info1,
                                                   MeshWorker::IntegrationInfo<dim> &info2)
  {
    const FEValuesBase<dim> &fe_v_1 = info1.fe_values();
    const FEValuesBase<dim> &fe_v_2 = info2.fe_values();

    Vector<double> &cell_vector_1 = dinfo1.vector(0).block(0);
    Vector<double> &cell_vector_2 = dinfo2.vector(0).block(0);

    const unsigned int   dofs_per_cell = fe_v_1.dofs_per_cell;
    const unsigned int   n_q_points    = fe_v_1.n_quadrature_points;

    std::vector<double> &u_1 = info1.values[0][0];
    std::vector<double> &u_2 = info1.values[0][0];

    double flux;

    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {

        if (fe_v_1.normal_vector(q_point)[0]>0)
          flux = u_1[q_point]*fe_v_1.normal_vector(q_point)[0];
        else
          flux = u_2[q_point]*fe_v_1.normal_vector(q_point)[0];

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {

            cell_vector_1(i) += wavespeed*flux*fe_v_1[upos].value(i,q_point)*fe_v_1.JxW(q_point);

            cell_vector_2(i) -= wavespeed*flux*fe_v_2[upos].value(i,q_point)*fe_v_1.JxW(q_point);

          }//i
      }//q_point

  }//integrate_face_term


  template <int dim>
  void AdvectionProblem<dim>::output_results (int timestep) const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);

    std::vector<std::string> solution_names(1, "u");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation (1, DataComponentInterpretation::component_is_scalar);

    data_out.add_data_vector (solution,
                              solution_names,
                              DataOut<dim,DoFHandler<dim> >::type_dof_data,
                              interpretation);

    data_out.build_patches (fe.degree);
    data_out.write_gnuplot (deallog.get_file_stream());

  }//output_results


  template <int dim>
  void AdvectionProblem<dim>::run ()
  {
    // Make the mesh
    create_grid();

    // Setup the system
    setup_system();

    deallog << "\tNumber of active cells:       "
            << triangulation.n_active_cells() << std::endl;

    deallog << "\tNumber of degrees of freedom: "
            << dof_handler.n_dofs() << std::endl;

    const double delta_t = 0.005;
    const double n_dt = 10;

    double inv_cell_vol = 1.0/std::pow(0.01, dim);

    for (unsigned int dt=0; dt<n_dt; ++dt)
      {

        assemble_rhs(solution, stage);

        stage *= delta_t;
        stage *= inv_cell_vol;
        solution -= stage;
        stage = 0.0;
      }

    output_results (n_dt);

  }//AdvectionProblem::run()

}//namespace


int main ()
{
  const std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console (0);

  Advection::AdvectionProblem<1> advection_problem;
  advection_problem.run ();
}
