// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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


/*  Test code to help fixing

DoFTools::extract_boundary_dofs
and

VectorTools::interpolate_boundary_values



The underlying FE is a FESystem with one Nedelec
component and one scalar component.

We would like to use
DoFTools::extract_boundary_dofs
with a component_mask disabling the
Nedelec and the scalar component respectively.

Furthermore, we would like to use
VectorTools::interpolate_boundary_values
with a component_mask enabling only
the scalar component.

author: Anna Schneebeli, February 2003

*/




#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>

#include <deal.II/grid/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>


#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>
#include <iomanip>


template <int dim>
class ImposeBC
{
public:
  ImposeBC();
  ~ImposeBC();
  void run ();

private:

  void get_ready ();
  void test_extract_boundary_DoFs ();
  void test_interpolate_BC ();

  Triangulation<dim>   triangulation;

  // We use a FE-System with 2 components:
  // a vector-valued one for the
  // field-variable u and a scalar one for
  // the pressure-variable p
  FESystem<dim>        fe;
  DoFHandler<dim>      dof_handler;

  unsigned int    n_u_dofs;
  unsigned int    n_p_dofs;

};


// some boundary function for the scalar component
template <int dim>
class BoundaryFunction : public Function<dim>
{
public:
  BoundaryFunction ();

  virtual void vector_value (const Point<dim>   &p,
                             Vector<double> &values) const;

};


template <int dim>
BoundaryFunction<dim>::BoundaryFunction () :
  Function<dim> (dim+1) {}



template <int dim>
inline
void BoundaryFunction<dim>::vector_value (const Point<dim> &,
                                          Vector<double> &values) const
{

  Assert (values.size() == dim+1,
          ExcDimensionMismatch (values.size(), dim+1));

  values = 0;
  values(dim) = 1.;
}





// Construct FE with first component: Nedelec-Element,
// second component: Q1_Element
template <int dim>
ImposeBC<dim>::ImposeBC() :
  fe (FE_Nedelec<dim>(0), 1, FE_Q<dim>(1), 1),
  dof_handler (triangulation)
{}



template <int dim>
ImposeBC<dim>::~ImposeBC()
{
  dof_handler.clear ();
}


template <int dim>
void ImposeBC<dim>::get_ready ()
{
  dof_handler.distribute_dofs (fe);
  std::vector<types::global_dof_index>  dofs_per_comp (fe.n_components());
  DoFTools::count_dofs_per_component(dof_handler, dofs_per_comp);

  // For an FESystem with Nedelec-elements as
  // first component and bilinear elements as
  // component we have:
  // dofs_per_comp[0] = dofs_per_comp[1] = # Ned-DoFs
  // dofs_per_comp[2] = # Q1-DoFs
  n_u_dofs = dofs_per_comp[0];
  n_p_dofs = dofs_per_comp[2];

}


template <int dim>
void ImposeBC<dim>::test_extract_boundary_DoFs ()
{

  std::map<types::global_dof_index,double> boundary_values;
  std::vector<bool>   bc_component_select(dim + 1);

  // extract boundary DoFs for the Nedelec-component
  // and impose zero boundary condition
  bc_component_select[0] = true;
  bc_component_select[1] = true;
  bc_component_select[2] = false;

  std::vector<bool> ned_boundary_dofs (dof_handler.n_dofs());
  std::set<types::boundary_id> boundary_indicators;
  boundary_indicators.insert (0);
  DoFTools::extract_boundary_dofs (dof_handler,
                                   bc_component_select,
                                   ned_boundary_dofs,
                                   boundary_indicators);


  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    if (ned_boundary_dofs[i] == true)
      boundary_values[i] = 0.;
}



template <int dim>
void ImposeBC<dim>::test_interpolate_BC ()
{

  std::map<types::global_dof_index,double> boundary_values;
  std::vector<bool>   bc_component_select(dim + 1, false);


  // impose inhomogeneous boundary condition
  // on the scalar variable
  bc_component_select.back() = true;

  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            BoundaryFunction<dim>(),
                                            boundary_values,
                                            bc_component_select);



  // check
  // (the pressure is assumed to be set to 1
  // on the boundary)
  std::vector<bool> p_boundary_dofs (dof_handler.n_dofs());
  std::set<types::boundary_id> boundary_indicators;
  boundary_indicators.insert (0);
  DoFTools::extract_boundary_dofs (dof_handler,
                                   bc_component_select,
                                   p_boundary_dofs,
                                   boundary_indicators);
  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
    {
      // error: pressure boundary DoF
      // i has not been set to the
      // correct value
      //
      // or:
      //
      // nedelec boundary DoF i has
      // wrongly been set to some
      // value
      Assert ((p_boundary_dofs[i] && boundary_values[i] == 1.)
              ||
              (!(p_boundary_dofs[i])  && boundary_values[i] != 1.),
              ExcInternalError());

      deallog << boundary_values[i] << ' ';
    }
  deallog << std::endl;
}




template <int dim>
void ImposeBC<dim>::run ()
{
  GridGenerator::hyper_cube(triangulation, -1,1);
  triangulation.refine_global (1);
  triangulation.begin_active()->set_refine_flag ();
  triangulation.execute_coarsening_and_refinement ();

  deallog << "   Number of active cells:       "
          << triangulation.n_active_cells()
          << std::endl;

  get_ready ();

  deallog << "  Total number of degrees of freedom: "
          << dof_handler.n_dofs()
          << std::endl
          << "   Number of degrees of freedom for the field variable U: "
          << n_u_dofs
          << std::endl
          << "   Number of degrees of freedom for the pressure variable p: "
          << n_p_dofs
          << std::endl;

  test_extract_boundary_DoFs ();
  test_interpolate_BC ();
}


int main ()
{
  try
    {
      std::ofstream logfile("output");
      deallog << std::setprecision (2);
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      ImposeBC<2>().run ();
      ImposeBC<3>().run ();
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };

  return 0;
}
