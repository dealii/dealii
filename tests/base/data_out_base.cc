//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

#include <base/data_out_base.h>
#include <base/logstream.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#define WRITE(type) DataOutBase::write_ ## type (patches, names, type ## flags, out)

template <int dim, int spacedim>
void
write_patches(const std::vector<DataOutBase::Patch<dim,spacedim> >& patches,
	      std::ostream& out)
{
  std::vector<std::string> names(2);
  names[0] = std::string("first name");
  names[1] = std::string("last name");

  DataOutBase::DXFlags dxflags;
//    DataOutBase::EpsFlags epsflags;
//    DataOutBase::GnuplotFlags gnuplotflags;
//    DataOutBase::GmvFlags gmvflags;
//    DataOutBase::PovrayFlags povrayflags;
//    DataOutBase::UcdFlags ucdflags;
//    DataOutBase::VtkFlags vtkflags;
  
  WRITE(dx);
  
}

template<int dim>
struct PatchInfo
{
    static double vertices[GeometryInfo<dim>::vertices_per_cell][3];
    static double offsets[GeometryInfo<dim>::vertices_per_cell][3];
    static int neighbors[GeometryInfo<dim>::vertices_per_cell][GeometryInfo<dim>::faces_per_cell];
};


double PatchInfo<1>::vertices[2][3]
= {{0, 0, 0}
  ,{3, 6, 9}
};

double PatchInfo<2>::vertices[4][3]
= {{0, 0, 0}
  ,{3, 0, 3}
  ,{3, 3, 6}
  ,{3, 0, 3}
};


double PatchInfo<3>::vertices[8][3]
= {{0, 0, 0}
  ,{3, 0, 0}
  ,{3, 3, 0}
  ,{0, 3, 0}
  ,{0, 0, 3}
  ,{3, 0, 3}
  ,{3, 3, 3}
  ,{0, 3, 3}
};


int PatchInfo<1>::neighbors[2][2]
= {{-1,  1}
  ,{ 0, -1}
};



int PatchInfo<2>::neighbors[4][4]
= {{-1,  1,  2, -1}
  ,{-1, -1,  3,  0}
  ,{ 0,  3, -1, -1}
  ,{ 1, -1, -1,  2}  
};



int PatchInfo<3>::neighbors[8][6]
= {{-1, -1, -1, -1, -1, -1}
  ,{-1, -1, -1, -1, -1, -1}
  ,{-1, -1, -1, -1, -1, -1}
  ,{-1, -1, -1, -1, -1, -1}
  ,{-1, -1, -1, -1, -1, -1}
  ,{-1, -1, -1, -1, -1, -1}
  ,{-1, -1, -1, -1, -1, -1}
  ,{-1, -1, -1, -1, -1, -1}  
};


double PatchInfo<1>::offsets[2][3]
= {{ 0, 0, 0}
  ,{ 3, 0, 0}
};


double PatchInfo<2>::offsets[4][3]
= {{ 0, 0, 0}
  ,{ 3, 0, 0}
  ,{ 0, 3, 0}
  ,{ 3, 3, 0}
};


double PatchInfo<3>::offsets[8][3]
= {{ 0, 0, 0}
  ,{ 3, 0, 0}
  ,{ 0, 3, 0}
  ,{ 3, 3, 0}
  ,{ 0, 0, 3}
  ,{ 3, 0, 3}
  ,{ 0, 3, 3}
  ,{ 3, 3, 3}
};



template <int dim, int spacedim>
void
create_patches(std::vector<DataOutBase::Patch<dim,spacedim> >& patches)
{
  const unsigned int ncells = GeometryInfo<dim>::vertices_per_cell;
  const unsigned int nsub = 3;
  
  patches.resize(ncells);

  for (unsigned int c = 0; c < ncells;++c)
    {
      DataOutBase::Patch<dim, spacedim>& p = patches[c];
      p.patch_index = c;
      p.n_subdivisions = nsub;

      for (unsigned int i=0;i<ncells;++i)
	for (unsigned int j=0;j<spacedim;++j)
	  p.vertices[i](j) = PatchInfo<dim>::vertices[i][j]
			     + PatchInfo<dim>::offsets[c][j];

      for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
	p.neighbors[i] = (unsigned int) PatchInfo<dim>::neighbors[c][i];

      unsigned int ndata = 1;
      for (unsigned int i=0;i<dim;++i)
	ndata *= nsub+1;

      p.data.reinit(2,ndata);
      for (unsigned int i=0;i<ndata;++i)
	{
	  p.data(0,i) = i;
	  p.data(1,i) = -(float) i;
	}
    }
}


int main()
{
  std::ofstream logfile("data_out_base.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  std::vector<DataOutBase::Patch<1,1> > patch11;
  std::vector<DataOutBase::Patch<1,2> > patch12;
  std::vector<DataOutBase::Patch<1,3> > patch13;
  std::vector<DataOutBase::Patch<2,2> > patch22;
  std::vector<DataOutBase::Patch<2,3> > patch23;
  std::vector<DataOutBase::Patch<3,3> > patch33;

  create_patches(patch33);
  write_patches(patch33, logfile);
}
