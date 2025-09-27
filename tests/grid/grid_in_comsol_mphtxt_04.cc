// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// COMSOL's .mphtxt files are often written on Windows, with Windows
// line endings, and then read on Unix systems. We used to trip over
// this -- check that we can read these kinds of files correctly
// today.
//
// To address the fact that version control systems sometimes silently
// change line endings to match the current system, this test just
// keeps an entire input file as part of the program, using \r\n line
// endings.

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"


const std::string mphtxt_file =
  "# Created by COMSOL Multiphysics.\r\n"
  "\r\n"
  "# Major & minor version\r\n"
  "0 1 \r\n"
  "1 # number of tags\r\n"
  "# Tags\r\n"
  "5 mesh1 \r\n"
  "1 # number of types\r\n"
  "# Types\r\n"
  "3 obj \r\n"
  "\r\n"
  "# --------- Object 0 ----------\r\n"
  "\r\n"
  "0 0 1 \r\n"
  "4 Mesh # class\r\n"
  "4 # version\r\n"
  "2 # sdim\r\n"
  "45 # number of mesh vertices\r\n"
  "0 # lowest mesh vertex index\r\n"
  "\r\n"
  "# Mesh vertex coordinates\r\n"
  "-0.92387960315740103 -0.38268326181023815 \r\n"
  "-0.70710705326435586 -0.70710650910863437 \r\n"
  "-1 0 \r\n"
  "-0.71406619972883534 -0.14205755756489793 \r\n"
  "-0.92387930866193591 0.3826839727847815 \r\n"
  "-0.59085324079676915 -0.37787505516896647 \r\n"
  "-0.42256406953621362 -0.60025528016845497 \r\n"
  "-0.17325295128021034 -0.65914392922937293 \r\n"
  "-0.3826839727847815 -0.92387930866193591 \r\n"
  "-0.6708725011093043 0.13101574920094311 \r\n"
  "-0.60122732614083518 0.40458527745291445 \r\n"
  "-0.45124977636911656 -0.098559526304571682 \r\n"
  "-0.70710650910863437 0.70710705326435586 \r\n"
  "-0.25493018820753449 -0.35128404104151728 \r\n"
  "-0.023164464492933783 -0.46258210857632936 \r\n"
  "0.11663081135853207 -0.66734668931627006 \r\n"
  "0 -1 \r\n"
  "-0.37495215999746051 0.57437821075560735 \r\n"
  "-0.33481239575510269 0.23944232883743233 \r\n"
  "-0.15519782964058473 -0.037943022808229039 \r\n"
  "-0.13680767722345469 0.7138742226585264 \r\n"
  "-0.38268326181023815 0.92387960315740103 \r\n"
  "0.031771014663086801 -0.2226293650346963 \r\n"
  "0.16400335633700858 -0.42598880075648382 \r\n"
  "0.38756623838007437 -0.56017610270489793 \r\n"
  "0.38268326181023815 -0.92387960315740103 \r\n"
  "-0.087075789946497451 0.45738077916301395 \r\n"
  "-0.040879464685735556 0.20216795620370126 \r\n"
  "0.14632396550773411 0.025672007644215806 \r\n"
  "0.13856986744657646 0.67220390032116684 \r\n"
  "0 1 \r\n"
  "0.34720721242110097 -0.23358638149716446 \r\n"
  "0.60780574566186374 -0.40314665931644073 \r\n"
  "0.70710650910863437 -0.70710705326435586 \r\n"
  "0.22369877004942476 0.34334510990034345 \r\n"
  "0.44582115830928287 0.086376558607743625 \r\n"
  "0.40440866092777095 0.60478154208149548 \r\n"
  "0.3826839727847815 0.92387930866193591 \r\n"
  "0.67327603747372811 -0.13217954817745833 \r\n"
  "0.92387930866193591 -0.3826839727847815 \r\n"
  "0.5699836410937158 0.37737288241539707 \r\n"
  "0.71494279978818498 0.13996316592589308 \r\n"
  "0.70710705326435586 0.70710650910863437 \r\n"
  "1 0 \r\n"
  "0.92387960315740103 0.38268326181023815 \r\n"
  "\r\n"
  "3 # number of element types\r\n"
  "\r\n"
  "# Type #0\r\n"
  "\r\n"
  "3 vtx # type name\r\n"
  "\r\n"
  "\r\n"
  "1 # number of vertices per element\r\n"
  "4 # number of elements\r\n"
  "# Elements\r\n"
  "2 \r\n"
  "16 \r\n"
  "30 \r\n"
  "43 \r\n"
  "\r\n"
  "4 # number of geometric entity indices\r\n"
  "# Geometric entity indices\r\n"
  "0 \r\n"
  "1 \r\n"
  "2 \r\n"
  "3 \r\n"
  "\r\n"
  "# Type #1\r\n"
  "\r\n"
  "3 edg # type name\r\n"
  "\r\n"
  "\r\n"
  "2 # number of vertices per element\r\n"
  "16 # number of elements\r\n"
  "# Elements\r\n"
  "1 0 \r\n"
  "0 2 \r\n"
  "2 4 \r\n"
  "8 1 \r\n"
  "4 12 \r\n"
  "16 8 \r\n"
  "12 21 \r\n"
  "25 16 \r\n"
  "21 30 \r\n"
  "33 25 \r\n"
  "30 37 \r\n"
  "39 33 \r\n"
  "37 42 \r\n"
  "43 39 \r\n"
  "42 44 \r\n"
  "44 43 \r\n"
  "\r\n"
  "16 # number of geometric entity indices\r\n"
  "# Geometric entity indices\r\n"
  "0 \r\n"
  "0 \r\n"
  "1 \r\n"
  "0 \r\n"
  "1 \r\n"
  "0 \r\n"
  "1 \r\n"
  "2 \r\n"
  "1 \r\n"
  "2 \r\n"
  "3 \r\n"
  "2 \r\n"
  "3 \r\n"
  "2 \r\n"
  "3 \r\n"
  "3 \r\n"
  "\r\n"
  "# Type #2\r\n"
  "\r\n"
  "3 tri # type name\r\n"
  "\r\n"
  "\r\n"
  "3 # number of vertices per element\r\n"
  "72 # number of elements\r\n"
  "# Elements\r\n"
  "3 2 0 \r\n"
  "5 0 1 \r\n"
  "3 0 5 \r\n"
  "6 5 1 \r\n"
  "6 1 8 \r\n"
  "6 8 7 \r\n"
  "9 4 2 \r\n"
  "3 9 2 \r\n"
  "10 4 9 \r\n"
  "3 11 9 \r\n"
  "3 5 11 \r\n"
  "10 12 4 \r\n"
  "6 13 5 \r\n"
  "6 7 13 \r\n"
  "11 5 13 \r\n"
  "14 13 7 \r\n"
  "14 7 15 \r\n"
  "7 8 16 \r\n"
  "15 7 16 \r\n"
  "10 17 12 \r\n"
  "10 18 17 \r\n"
  "10 9 18 \r\n"
  "9 11 18 \r\n"
  "19 18 11 \r\n"
  "11 13 19 \r\n"
  "17 21 12 \r\n"
  "20 21 17 \r\n"
  "14 22 13 \r\n"
  "22 19 13 \r\n"
  "23 14 15 \r\n"
  "23 22 14 \r\n"
  "23 15 24 \r\n"
  "15 16 25 \r\n"
  "24 15 25 \r\n"
  "20 17 26 \r\n"
  "26 17 18 \r\n"
  "27 26 18 \r\n"
  "27 18 19 \r\n"
  "22 28 19 \r\n"
  "27 19 28 \r\n"
  "20 26 29 \r\n"
  "20 30 21 \r\n"
  "20 29 30 \r\n"
  "23 24 31 \r\n"
  "23 31 22 \r\n"
  "22 31 28 \r\n"
  "32 31 24 \r\n"
  "24 25 33 \r\n"
  "32 24 33 \r\n"
  "26 34 29 \r\n"
  "27 34 26 \r\n"
  "27 28 34 \r\n"
  "28 31 35 \r\n"
  "28 35 34 \r\n"
  "36 29 34 \r\n"
  "29 37 30 \r\n"
  "36 37 29 \r\n"
  "32 38 31 \r\n"
  "35 31 38 \r\n"
  "32 33 39 \r\n"
  "32 39 38 \r\n"
  "36 34 40 \r\n"
  "35 40 34 \r\n"
  "41 35 38 \r\n"
  "41 40 35 \r\n"
  "36 42 37 \r\n"
  "36 40 42 \r\n"
  "38 39 43 \r\n"
  "41 38 43 \r\n"
  "41 43 44 \r\n"
  "40 44 42 \r\n"
  "41 44 40 \r\n"
  "\r\n"
  "72 # number of geometric entity indices\r\n"
  "# Geometric entity indices\r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n"
  "1 \r\n";



template <int dim, int spacedim>
void
comsol_grid()
{
  Triangulation<dim, spacedim> tria;
  GridIn<dim, spacedim>        grid_in;
  grid_in.attach_triangulation(tria);
  std::istringstream input_file(mphtxt_file);
  grid_in.read_comsol_mphtxt(input_file);

  deallog << "  " << tria.n_active_cells() << " active cells" << std::endl;

  int hash  = 0;
  int index = 0;
  for (typename Triangulation<dim, spacedim>::active_cell_iterator c =
         tria.begin_active();
       c != tria.end();
       ++c, ++index)
    for (const unsigned int i : c->vertex_indices())
      hash += (index * i * c->vertex_index(i)) % (tria.n_active_cells() + 1);
  deallog << "  hash=" << hash << std::endl;

#if 0
  {
    static int    output_file = 0;
    GridOut       go;
    std::ofstream o("x-" + std::to_string(output_file) + ".gnuplot");
    go.write_gnuplot(tria, o);

    ++output_file;
  }

#endif
}


int
main()
{
  initlog();

  try
    {
      comsol_grid<2, 2>();
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
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
      deallog << std::endl
              << std::endl
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
