//---------------------------------------------------------------------------
//    $Id: coarse_grid_01.cc 16363 2008-06-16 20:42:57Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include "../tests.h"
#include <base/logstream.h>
#include <base/utilities.h>
#include <distributed/tria.h>
#include <grid/tria.h>
#include <grid/grid_out.h>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <unistd.h>

void cat_file(const char * filename)
{
  std::ifstream in(filename);
  while (in)
    {
      std::string s;
      std::getline(in, s);
      deallog.get_file_stream() << s << "\n";
    }
}


template <int dim>
void write_vtk (const parallel::distributed::Triangulation<dim> &tria,
		const char *                               dirname,
		const char *                               filename)
{
  deallog << "Checksum: "
	  << tria.get_checksum ()
	  << std::endl;

  chdir (dirname);
  tria.write_mesh_vtk (filename);

  if (false)
    {

      int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);


				       // copy the .pvtu and .vtu files
				       // into the logstream
      if (myid==0)
	{
	  std::ifstream in((std::string(filename) + ".pvtu").c_str());
	  while (in)
	    {
	      std::string s;
	      std::getline(in, s);
	      deallog.get_file_stream() << s << "\n";
	    }
	}

      if (myid==0)
	{
	  std::ifstream in((std::string(filename) + "_0000.vtu").c_str());
	  while (in)
	    {
	      std::string s;
	      std::getline(in, s);
	      deallog.get_file_stream() << s << "\n";
	    }
	}
    }


  chdir ("..");
}


template <int dim>
void assert_tria_equal(const char* testdir, const Triangulation<dim> & a, const Triangulation<dim> & b)
{
  Assert (a.n_active_cells() == b.n_active_cells(),
	      ExcInternalError());

  std::string file1=std::string(testdir)+"/tmp_grid1";
  std::string file2=std::string(testdir)+"/tmp_grid2";

  std::ofstream out1(file1.c_str());
  GridOut().write_gnuplot (a, out1);
  std::ofstream out2(file2.c_str());
  GridOut().write_gnuplot (b, out2);

  out1.close();
  out2.close();

				   //compare the two files
  std::string cmd = std::string("diff -q ")+file1+std::string(" ")+file2;
  Assert (system(cmd.c_str()) == 0, ExcInternalError());

				   //and delete them
  std::remove (file1.c_str());
  std::remove (file2.c_str());
}


template <typename T>
LogStream &
operator << (LogStream &out,
	     const std::vector<T> &v)
{
  for (unsigned int i=0; i<v.size(); ++i)
    out << v[i] << (i == v.size()-1 ? "" : " ");
  return out;
}
