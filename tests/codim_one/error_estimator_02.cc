// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2015 by the deal.II authors
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


// check that the kelly error estimator uses correct normals on surfaces
// with kinks by interpolating a function that is inside the FE space
// and should produce zero errors.



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>


template<int dim>
class MyFunction : public Function<dim>
{
public:
  MyFunction () : Function<dim>(1) {}

  virtual double value (const Point<dim>   &p,
                        const unsigned int  component) const
  {
    return p(0)*p(0)+p(0)*p(1);
  }

  virtual void   vector_value (const Point<dim>   &p,
                               Vector<double>     &values) const
  {
    values(0) = value(p,0);
    values(1) = value(p,1);
  }
};



template <int dim>
Quadrature<dim-1> &
get_q_face ()
{
  static QGauss<dim-1> q(4);
  return q;
}

template<>
Quadrature<0> &
get_q_face <1>()
{
  Quadrature<0> *q = 0;
  return *q;
}

template <int dim, int spacedim>
void make_mesh (Triangulation<dim,spacedim> &tria)
{
  // two faces of a hyper_cube
  Triangulation<spacedim,spacedim> volume_mesh;
  GridGenerator::hyper_cube(volume_mesh, 0, 1, true);
  std::set<types::boundary_id> boundary_ids;
  boundary_ids.insert (1);
  boundary_ids.insert (2);
  GridGenerator::extract_boundary_mesh (volume_mesh, tria,
                                            boundary_ids);
  tria.refine_global (2);
}




template <int dim, int spacedim>
void
check ()
{
  MyFunction<spacedim> function;

  Triangulation<dim,spacedim> tria;
  make_mesh (tria);

  FE_Q<dim,spacedim> element(2);
  DoFHandler<dim,spacedim> dof(tria);
  dof.distribute_dofs(element);

  MappingQ<dim,spacedim> mapping(3);
  Quadrature<dim-1> &q_face = get_q_face<dim>();

  std::map<types::boundary_id,const Function<spacedim>*> neumann_bc;
  neumann_bc[0] = &function;

  Vector<double> v (dof.n_dofs());
  VectorTools::interpolate (mapping, dof, function, v);

  Vector<float> error (tria.n_active_cells());

  KellyErrorEstimator<dim,spacedim>::estimate (mapping, dof, q_face,
					       typename FunctionMap<spacedim>::type(),
                                               v, error);
  deallog << "Estimated error indicators:" << std::endl;
  for (unsigned int i=0; i<error.size(); ++i)
    deallog << error(i) << std::endl;

  {
    DataOut<dim,DoFHandler<dim,spacedim> > data_out;
    data_out.attach_dof_handler (dof);
    data_out.add_data_vector (v,
                              "solution",
                              DataOut<dim,DoFHandler<dim,spacedim> >::type_dof_data);
    data_out.add_data_vector (error,
                              "error");
    data_out.build_patches ();
    std::string filename = spacedim == 2 ?
          "solution-2d-" :
          "solution-3d-";
    filename += Utilities::int_to_string(0,2) + ".vtk";
    std::ofstream output (filename.c_str());
    data_out.write_vtk (output);
    
  }
  
  deallog << "OK" << std::endl;
}


int main ()
{
  initlog();
  
  deallog.push ("1d");
  check<1,2> ();
  deallog.pop();
  deallog.push ("2d");
  check<2,3> ();
  deallog.pop();
}
