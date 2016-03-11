// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2015 by the deal.II authors
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



// like _01, but with the real testcase by Krishna Garikipati

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/hp/q_collection.h>

#include <fstream>
#include <vector>


#define DIMS 3
#define totalDOF 4

template <int dim>
class diffusionMechanics
{
public:
  diffusionMechanics (const unsigned int mech_degree, const unsigned int diff_degree);
  ~diffusionMechanics();
  void run ();

  //Input lengths of rectangular prism
  double alen, blen, clen;
private:
  enum
  {
    omega1_domain_id,
    omega2_domain_id
  };

  static bool
  cell_is_in_omega1_domain (const typename hp::DoFHandler<dim>::cell_iterator &cell);
  static bool
  cell_is_in_omega2_domain (const typename hp::DoFHandler<dim>::cell_iterator &cell);
  void set_active_fe_indices ();
  void setup_dofs ();
  void setup_system();

  const unsigned int   mech_degree;
  const unsigned int   diff_degree;

  Triangulation<dim>   triangulation;
  FESystem<dim>        omega1_fe;
  FESystem<dim>        omega2_fe;
  hp::FECollection<dim> fe_collection;
  hp::DoFHandler<dim>      dof_handler;
  QGauss<dim>          quadrature_formula;
  QGauss<dim-1>        face_quadrature_formula;
  SparsityPattern      sparsity_pattern;
  std::map<types::global_dof_index, double>    boundary_values;

  Vector<double>       system_rhs, U, Un, dU, U0;

  //solution variables
  std::vector<std::string> nodal_solution_names;
};

// Constructor
template <int dim>
diffusionMechanics<dim>::diffusionMechanics (const unsigned int mech_degree, const unsigned int diff_degree):
  mech_degree (mech_degree),
  diff_degree (diff_degree),
  omega1_fe (FE_Q<dim>(mech_degree), dim,
             FE_Q<dim>(diff_degree),   1),
  omega2_fe (FE_Q<dim>(mech_degree), dim,
             FE_Nothing<dim>(),   1),
  dof_handler (triangulation), quadrature_formula(3), face_quadrature_formula(2)
{

  //Nodal Solution names
  for (unsigned int i=0; i<dim; ++i)
    {
      nodal_solution_names.push_back("u");
    }
  nodal_solution_names.push_back("c");

  //FE object
  fe_collection.push_back(omega1_fe);
  fe_collection.push_back(omega2_fe);
  alen = 10.0, blen = 10.0, clen = 4.0;
}
template <int dim>
diffusionMechanics<dim>::~diffusionMechanics ()
{
  dof_handler.clear ();
}

//Initial conditions - use only if FE_Nothing not used. Else an error occurs with VectorTools::interpolate
template <int dim>
class InitialConditions: public Function<dim>
{
public:
  InitialConditions (): Function<dim>(totalDOF) {}
  void vector_value (const Point<dim>   &p, Vector<double>   &values) const
  {
    Assert (values.size() == totalDOF, ExcDimensionMismatch (values.size(), totalDOF));
    values(totalDOF-4)=0; // u=0
    values(totalDOF-3)=0;
    values(totalDOF-2)=0;
    values(totalDOF-1)=1.0;
  };
};
// Boolean functions for specifying sub-domains
template <int dim>
bool
diffusionMechanics<dim>::
cell_is_in_omega1_domain (const typename hp::DoFHandler<dim>::cell_iterator &cell)
{
  return (cell->material_id() == omega1_domain_id);
}
template <int dim>
bool
diffusionMechanics<dim>::
cell_is_in_omega2_domain (const typename hp::DoFHandler<dim>::cell_iterator &cell)
{
  return (cell->material_id() == omega2_domain_id);
}

// Set active fe indices in each sub-domain
template <int dim>
void
diffusionMechanics<dim>::set_active_fe_indices()
{
  for (typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    {
      if (cell_is_in_omega1_domain(cell))
        cell->set_active_fe_index(0);
      else if (cell_is_in_omega2_domain(cell))
        cell->set_active_fe_index(1);
      else
        Assert(false, ExcNotImplemented());
    }
}

// Call function to set active indices and distribute dofs
template <int dim>
void
diffusionMechanics<dim>::setup_dofs()
{
  set_active_fe_indices ();
  dof_handler.distribute_dofs(fe_collection);
}

//====================done====================================

//Setup
template <int dim>
void diffusionMechanics<dim>::setup_system()
{
  dof_handler.distribute_dofs (fe_collection);
  sparsity_pattern.reinit (dof_handler.n_dofs(), dof_handler.n_dofs(), dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();
  U.reinit (dof_handler.n_dofs());
  Un.reinit (dof_handler.n_dofs());
  dU.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
  U0.reinit (dof_handler.n_dofs());
  deallog << "   Number of active cells:       " << triangulation.n_active_cells() << std::endl;
  deallog << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;
}


//Run
template <int dim>
void diffusionMechanics<dim>::run ()
{

  //subdivided_hyper
  std::vector<unsigned int> repetitions(3);
  repetitions[0] = 2;
  repetitions[1] = 2;
  repetitions[2] = 4;
  GridGenerator::subdivided_hyper_rectangle (triangulation, repetitions,
                                             Point<3>(0.0,0.0,0.0),
                                             Point<3>(alen,blen,clen));

  //Mark cells by sub-domain
  for (typename Triangulation<dim>::active_cell_iterator cell = dof_handler.begin_active();
       cell != dof_handler.end(); ++cell)
    if (cell->center()[dim-1] >= 0.75*clen)
      cell->set_material_id(omega1_domain_id);
    else
      cell->set_material_id(omega2_domain_id);
  setup_dofs();
  setup_system();

  //Initial conditions
  VectorTools::interpolate(dof_handler, InitialConditions<dim>(), Un);
  //  VectorTools::interpolate(dof_handler, ZeroFunction<dim>(fe_collection.n_components()), Un);
  U=Un;
  U0=Un;

  deallog << Un.l2_norm() << std::endl;
}


int main ()
{
  std::ofstream logfile("output");
  logfile.precision (3);
  deallog << std::setprecision(3);

  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  diffusionMechanics<DIMS> problem(1,1);
  problem.run ();
}

