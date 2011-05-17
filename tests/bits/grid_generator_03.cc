//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// Test GridGenerator::subdivided_hyper_rectangle with vector of step
// sizes and the table argument for material id

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/table.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>
#include <iomanip>


template <int dim> Table<dim,unsigned char> material_ids();

template <> Table<1,unsigned char> material_ids<1>()
{
  Table<1,unsigned char> t(2);
  for (unsigned int i=0; i<2; ++i)
    t[i] = 1;
  return t;
}


template <> Table<2,unsigned char> material_ids<2>()
{
  Table<2,unsigned char> t(2,3);
  for (unsigned int i=0; i<2; ++i)
    for (unsigned int j=0; j<3; ++j)
      t[i][j] = 1;
				   // produce a hole in the middle
  t[1][1] = (unsigned char)(-1);
  return t;
}


template <> Table<3,unsigned char> material_ids<3>()
{
  Table<3,unsigned char> t(2,3,4);
  for (unsigned int i=0; i<2; ++i)
    for (unsigned int j=0; j<3; ++j)
      for (unsigned int k=0; k<4; ++k)
	t[i][j][k] = 1;
				   // produce a hole in the middle
  t[1][1][1] = (unsigned char)(-1);
  t[1][1][2] = (unsigned char)(-1);
  return t;
}



template<int dim>
void test(std::ostream& out)
{
  Point<dim> p1;
  p1[0] = 2.;
  if (dim>1) p1[1] = -1.;
  if (dim>2) p1[2] = 0.;
  
  Point<dim> p2;
  p2[0] = 3.;
  if (dim>1) p2[1] = 2.;
  if (dim>2) p2[2] = 4.;

  Point<dim> p3;
  p3[0] = 2.;
  if (dim>1) p3[1] = 1.;
  if (dim>2) p3[2] = 4.;
  
  GridOut go;

                                   // uniformly subdivided mesh
  if (true)
    {
      deallog << "subdivided_hyper_rectangle" << std::endl;
      Triangulation<dim> tr;
      std::vector<std::vector<double> > sub(dim);
      for (unsigned int i=0; i<dim; ++i)
        sub[i] = std::vector<double> (i+2, (p2[i]-p1[i])/(i+2));

      GridGenerator::subdivided_hyper_rectangle(tr, sub, p1, material_ids<dim>(),
						dim!=1);
      if (tr.n_cells() > 0)
	go.write_gnuplot(tr, out);
    }


                                   // non-uniformly subdivided mesh
  if (true)
    {
      deallog << "subdivided_hyper_rectangle" << std::endl;
      Triangulation<dim> tr;
      std::vector<std::vector<double> > sub(dim);
      for (unsigned int i=0; i<dim; ++i)
        {
          sub[i] = std::vector<double> (i+2, (p2[i]-p1[i])/(i+2));
          sub[i][0] /= 2;
          sub[i].back() *= 1.5;
        }

      GridGenerator::subdivided_hyper_rectangle(tr, sub, p1, material_ids<dim>(),
						dim!=1);
      if (tr.n_cells() > 0)
	go.write_gnuplot(tr, out);
    }
}


int main()
{
  std::ofstream logfile("grid_generator_03/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);
  
  deallog.push("1d");
  test<1>(logfile);
  deallog.pop();
  deallog.push("2d");
  test<2>(logfile);
  deallog.pop();
  deallog.push("3d");
  test<3>(logfile);
  deallog.pop();
}
