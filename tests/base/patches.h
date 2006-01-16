//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// Function for dreating a set of patches with different subgrids.

// To be included after all headers


double cell_coordinates [3][8] =
{
      { 0, 1, 0, 1, 0, 1, 0, 1 },
      { 0, 0, 1, 1, 0, 0, 1, 1 },
      { 0, 0, 0, 0, 1, 1, 1, 1 }
};

      
template <int dim, int spacedim>
void
create_patches(std::vector<DataOutBase::Patch<dim, spacedim> > & patches)
{
  for (unsigned int p=0;p<patches.size();++p)
    {
      DataOutBase::Patch<dim, spacedim>& patch = patches[p];
      
      const unsigned int nsub = p+1;
      const unsigned int nsubp = nsub+1;
      
      patch.n_subdivisions = nsub;
      for (unsigned int v=0;v<GeometryInfo<dim>::vertices_per_cell;++v)
	for (unsigned int d=0;d<spacedim;++d)
	  patch.vertices[v](d) = p+cell_coordinates[d][v]
				 + ((d>=dim) ? v : 0);
      
      unsigned int n1 = (dim>0) ? nsubp : 1;
      unsigned int n2 = (dim>1) ? nsubp : 1;
      unsigned int n3 = (dim>2) ? nsubp : 1;
      unsigned int n4 = (dim>3) ? nsubp : 1;
      patch.data.reinit(5, n1*n2*n3*n4);
      
      for (unsigned int i4=0;i4<n4;++i4)
	for (unsigned int i3=0;i3<n3;++i3)
	  for (unsigned int i2=0;i2<n2;++i2)
	    for (unsigned int i1=0;i1<n1;++i1)
	      {
		const unsigned int i = i1+nsubp*(i2+nsubp*(i3+nsubp*i4));
		const float x1 = 2.*i1/n1-1.;
		const float x2 = 2.*i2/n2-1.;
		const float x3 = 2.*i3/n3-1.;
		const float x4 = 2.*i4/n4-1.;
		
		patch.data(0,i) = p+x1;
		patch.data(1,i) = p+x2;
		patch.data(2,i) = p+x3;
		patch.data(3,i) = p+x4;
		patch.data(4,i) = i;
	      }
    }
}

