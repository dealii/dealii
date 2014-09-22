// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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


// crash with fe nedelec reported on mailing list (01/14/2014)
/*
#0  0x00007fffeba50425 in __GI_raise (sig=<optimized out>) at ../nptl/sysdeps/unix/sysv/linux/raise.c:64
#1  0x00007fffeba53b8b in __GI_abort () at abort.c:91
#2  0x00007ffff4f95895 in dealii::deal_II_exceptions::internals::abort (exc=..., nothrow=false)
    at /ssd/deal-trunk/deal.II/source/base/exceptions.cc:323
#3  0x00007ffff18ef959 in dealii::deal_II_exceptions::internals::issue_error<dealii::StandardExceptions::ExcIndexRange> (
    handling=dealii::deal_II_exceptions::internals::abort_on_exception, file=0x7ffff6035000 "/ssd/deal-trunk/deal.II/source/fe/fe.cc", 
    line=557, 
    function=0x7ffff6039240 "unsigned int dealii::FiniteElement<dim, spacedim>::face_to_cell_index(unsigned int, unsigned int, bool, bool, bool) const [with int dim = 3, int spacedim = 3]", cond=0x7ffff6035db0 "face_index < this->dofs_per_face", 
    exc_name=0x7ffff6035d78 "ExcIndexRange(face_index, 0, this->dofs_per_face)", e=...)
    at /ssd/deal-trunk/deal.II/include/deal.II/base/exceptions.h:272
#4  0x00007ffff31708c1 in dealii::FiniteElement<3, 3>::face_to_cell_index (this=0xa16f70, face_index=4, face=0, face_orientation=true, 
    face_flip=false, face_rotation=false) at /ssd/deal-trunk/deal.II/source/fe/fe.cc:556
#5  0x00007ffff2d7fe1c in dealii::VectorTools::internals::compute_edge_projection<dealii::TriaActiveIterator<dealii::DoFCellAccessor<dealii::DoFHandler<3, 3>, false> > > (cell=..., face=0, line=3, hp_fe_values=..., boundary_function=..., first_vector_component=0, 
    dof_values=std::vector of length 8, capacity 8 = {...}, dofs_processed=std::vector<bool> of length 8, capacity 64 = {...})
    at /ssd/deal-trunk/deal.II/include/deal.II/numerics/vector_tools.templates.h:2944
#6  0x00007ffff2d77846 in dealii::VectorTools::project_boundary_values_curl_conforming<3> (dof_handler=..., first_vector_component=0, 
    boundary_function=..., boundary_component=0 '\000', constraints=..., mapping=...)
    at /ssd/deal-trunk/deal.II/include/deal.II/numerics/vector_tools.templates.h:3543
#7  0x000000000042b479 in test<3> (fe=...) at /ssd/deal-trunk/tests/fe/crash_01.cc:74
#8  0x0000000000429a4e in main () at /ssd/deal-trunk/tests/fe/crash_01.cc:91




*/


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_nedelec.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/function.h>

#include <vector>
#include <fstream>
#include <string>


template <int dim>
void test(FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, -1.0, 1.0);
  tria.refine_global(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs (fe);
  ConstraintMatrix constraints;
  ZeroFunction<dim> boundary_values(fe.n_components());
  VectorTools::project_boundary_values_curl_conforming (dof_handler, 0, boundary_values, 0, constraints);
}



int
main()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog << std::setprecision(7);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  FE_Nedelec<3> fe1(0);  // works
  test<3>(fe1);
  FESystem<3> fe2(FE_Nedelec<3> (0), 2); // crash
  test<3>(fe2);

  return 0;
}
