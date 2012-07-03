//-----------------------------------------------------------------------------
//    $Id: mapping_real_to_unit_q1.cc 25655 2012-06-27 21:40:00Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2006, 2007, 2010, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

/*

  reduced failing test. Original exception in transform_real_to_unit_cell:
  
#0  0x00007fffeec931a0 in __cxa_throw () from /usr/lib64/libstdc++.so.6
#1  0x00007ffff5215d60 in dealii::deal_II_exceptions::internals::issue_error_throw<dealii::Mapping<2, 2>::ExcTransformationFailed> (
    file=0x7ffff6b47ce8 "/w/kwang/kwang/deal.II/source/fe/mapping_q1.cc", line=1803, 
    function=0x7ffff6b54e80 "dealii::Point<dim, double> dealii::MappingQ1<dim, spacedim>::transform_real_to_unit_cell_internal(const typename dealii::Triangulation<dim, spacedim>::cell_iterator&, const dealii::Point<spacedim>&, c"..., cond=0x7ffff6b47b3d "false", 
    exc_name=0x7ffff6b47fd8 "(typename Mapping<dim,spacedim>::ExcTransformationFailed())", e=...)
    at /w/kwang/kwang/deal.II/include/deal.II/base/exceptions.h:304
#2  0x00007ffff5201d5f in dealii::MappingQ1<2, 2>::transform_real_to_unit_cell_internal (this=0x7ffff7dd8160, cell={
  triangulation = 0x7fffffffd2b0,
  level = 1,
  index = 4
}, p={-0.27999999999999992, -0.62999999999999989}, initial_p_unit={1.4871340295728577, -0.21318360089489402}, mdata=...)
    at /w/kwang/kwang/deal.II/source/fe/mapping_q1.cc:1803
#3  0x00007ffff51ff24a in dealii::MappingQ1<2, 2>::transform_real_to_unit_cell (this=0x7ffff7dd8160, cell={
  triangulation = 0x7fffffffd2b0,
  level = 1,
  index = 4
}, p={-0.27999999999999992, -0.62999999999999989}) at /w/kwang/kwang/deal.II/source/fe/mapping_q1.cc:1683
#4  0x00007ffff4037ea2 in dealii::Functions::FEFieldFunction<2, dealii::DoFHandler<2, 2>, dealii::Vector<double> >::compute_point_locations (this=0x7fffffffcc80, points=..., cells=..., qpoints=..., maps=...)
    at /w/kwang/kwang/deal.II/include/deal.II/numerics/fe_field_function.templates.h:443
#5  0x00007ffff4035808 in dealii::Functions::FEFieldFunction<2, dealii::DoFHandler<2, 2>, dealii::Vector<double> >::vector_value_list (this=0x7fffffffcc80, points=..., values=...) at /w/kwang/kwang/deal.II/include/deal.II/numerics/fe_field_function.templates.h:208
#6  0x00007ffff40355e5 in dealii::Functions::FEFieldFunction<2, dealii::DoFHandler<2, 2>, dealii::Vector<double> >::value_list (
    this=0x7fffffffcc80, points=..., values=..., component=0)
    at /w/kwang/kwang/deal.II/include/deal.II/numerics/fe_field_function.templates.h:246
#7  0x0000000000450fd5 in Step6<2>::compute_measure (this=0x7fffffffd2b0) at step-6.cc:627
#8  0x000000000045874c in Step6<2>::run (this=0x7fffffffd2b0) at step-6.cc:705
#9  0x000000000045127f in main () at step-6.cc:728


 Note that test() only looks at the one cell that causes the failure.
 
test2() is the original code using FEFieldFunction. Note that it depends on the contents of the points in the value_list (see comments).
 
 */

#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/fe_field_function.h>


template<int dim>
void test()
{
  const HyperBallBoundary<dim> boundary_description;

  Triangulation<dim>   triangulation;
  GridGenerator::hyper_ball (triangulation);
  triangulation.set_boundary (0, boundary_description);
  triangulation.refine_global (1);
  MappingQ1< dim > mapping;


  Point<dim> p(-0.27999999999999992, -0.62999999999999989);
  typename Triangulation<dim>::active_cell_iterator cell;
  
  for (cell = triangulation.begin_active();
       cell != triangulation.end(); cell++)
    {
      if (cell->level()==1 && cell->index()==4)
	mapping.transform_real_to_unit_cell(cell, p);
      

    }
  
  triangulation.set_boundary (0);
  
  deallog << "OK" << std::endl;
}


template<int dim>
void test2()
{
  const HyperBallBoundary<dim> boundary_description;

  Triangulation<dim>   triangulation;
  GridGenerator::hyper_ball (triangulation);
  triangulation.set_boundary (0, boundary_description);
  triangulation.refine_global (1);
  MappingQ1< dim > mapping;


  Point<dim> p(-0.27999999999999992, -0.62999999999999989);

  FE_Q<dim> fe(2);
  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  Vector<double> solution(dof_handler.n_dofs());
  Functions::FEFieldFunction<2> fe_function (dof_handler, solution);
  fe_function.value (p); //this works <<<<<<<<<<<

  std::vector<Point<dim> > points(19*19);
  std::vector<double> m (19*19);
  
  if (1) //works if changed to "if (0)"   <<<<<<<<<
    for (unsigned int i = 0; i < 19; i++)
    for (unsigned int j = 0; j < 19; j++)
      { /// all points are inside
        points[19*i+j] (0) = -0.7 + (i + 1) * .07;
        points[19*i+j] (1) = -0.7 + (j + 1) * .07;
      }
  points[95]=p;
  fe_function.value_list (points, m); // <<<< this fails at point[95] but only if the other points are filled in?!
  
  triangulation.set_boundary (0);
  
  deallog << "OK" << std::endl;
}


int
main()
{
  std::ofstream logfile ("mapping_real_to_unit_02/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2>();
  test2<2>();

  return 0;
}



