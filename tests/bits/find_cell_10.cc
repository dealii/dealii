// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



//strange hang/crash in find_active_cell_around_point depending on which
//mapping is used

// this should also not hang, but throw an Exception. Call stack:
/*
#2  0x00007fffed53f20a in __backtrace_symbols (array=0x7fffffffc610, size=9)
    at ../sysdeps/generic/elf/backtracesyms.c:52
#3  0x00007ffff585cd9d in dealii::ExceptionBase::set_fields (this=0x7fffffffc870,
    f=0x7ffff668dec0 "/scratch/deal-trunk/deal.II/source/fe/mapping_q1.cc", l=1803,
    func=0x7ffff6699be0 "dealii::Point<dim> dealii::MappingQ1<dim, spacedim>::transform_real_to_unit_cell_internal(const typename dealii::Triangulation<dim, spacedim>::cell_iterator&, const dealii::Point<spacedim>&, const dea"..., c=0x7ffff668dea5 "false", e=0x7ffff668e1c0 "(typename Mapping<dim,spacedim>::ExcTransformationFailed())")
    at /scratch/deal-trunk/deal.II/source/base/exceptions.cc:124
#4  0x00007ffff4663508 in dealii::deal_II_exceptions::internals::issue_error_throw<dealii::Mapping<2, 2>::ExcTransformationFailed> (file=0x7ffff668dec0 "/scratch/deal-trunk/deal.II/source/fe/mapping_q1.cc", line=1803,
    function=0x7ffff6699be0 "dealii::Point<dim> dealii::MappingQ1<dim, spacedim>::transform_real_to_unit_cell_internal(const typename dealii::Triangulation<dim, spacedim>::cell_iterator&, const dealii::Point<spacedim>&, const dea"..., cond=0x7ffff668dea5 "false",
    exc_name=0x7ffff668e1c0 "(typename Mapping<dim,spacedim>::ExcTransformationFailed())", e=...)
    at /scratch/deal-trunk/deal.II/include/deal.II/base/exceptions.h:245
#5  0x00007ffff46505cb in dealii::MappingQ1<2, 2>::transform_real_to_unit_cell_internal (this=0x7ffff7dd9560,
    cell=..., p=..., initial_p_unit=..., mdata=...) at /scratch/deal-trunk/deal.II/source/fe/mapping_q1.cc:1803
#6  0x00007ffff464dd6b in dealii::MappingQ1<2, 2>::transform_real_to_unit_cell (this=0x7ffff7dd9560, cell=...,
    p=...) at /scratch/deal-trunk/deal.II/source/fe/mapping_q1.cc:1683
#7  0x00007ffff5a1d467 in dealii::GridTools::find_active_cell_around_point<2, dealii::DoFHandler, 2> (
    mapping=..., container=..., p=...) at /scratch/deal-trunk/deal.II/source/grid/grid_tools.cc:892
#8  0x000000000040e49d in test () at find_cell_10.cc:95
#9  0x000000000040e5ec in main (argc=1, argv=0x7fffffffdf48) at find_cell_10.cc:103



 */

#include "../tests.h"

#include <stdio.h>
#include <cstdlib>

#include <base/quadrature_lib.h>
#include <fe/mapping_q.h>
#include <base/function.h>
#include <base/logstream.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_tools.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <grid/grid_in.h>


#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <sstream>
#include <time.h>

using namespace dealii;

void test()
{
  Triangulation<2>     triangulation;

  GridIn<2> grid_in1;
  grid_in1.attach_triangulation (triangulation);
  std::ifstream input_file1(SOURCE_DIR "/grids/mesh.msh");
  grid_in1.read_msh(input_file1);

  Point< 2 > ePos;
  ePos(0) = 0.0653630060373507487669897386695;
  ePos(1) = 1125.59175030825804242340382189;

  MappingQ<2> mapping(1);
  MappingQ1<2> &mapping2 = StaticMappingQ1< 2 >::mapping;
  deallog << "1:" << std::endl;
  GridTools::find_active_cell_around_point (mapping, triangulation, ePos);
  deallog << "2:" << std::endl;
  //this second call seems to hang/crash:
  GridTools::find_active_cell_around_point (mapping2, triangulation, ePos);
  deallog << "done" << std::endl;
}

int main (int argc, char **argv)
{
  initlog();

  test();

  return 0;
}
