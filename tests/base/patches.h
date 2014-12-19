// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2010 by the deal.II authors
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

// Function for creating a set of patches with different subgrids.

// To be included after all headers


double cell_coordinates [3][8] =
{
  { 0, 1, 0, 1, 0, 1, 0, 1 },
  { 0, 0, 1, 1, 0, 0, 1, 1 },
  { 0, 0, 0, 0, 1, 1, 1, 1 }
};


template <int dim, int spacedim>
void
create_patches(std::vector<DataOutBase::Patch<dim, spacedim> > &patches)
{
  for (unsigned int p=0; p<patches.size(); ++p)
    {
      DataOutBase::Patch<dim, spacedim> &patch = patches[p];

      const unsigned int nsub = p+1;
      const unsigned int nsubp = nsub+1;

      patch.n_subdivisions = nsub;
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
        for (unsigned int d=0; d<spacedim; ++d)
          patch.vertices[v](d) = p+cell_coordinates[d][v]
                                 + ((d>=dim) ? v : 0);

      unsigned int n1 = (dim>0) ? nsubp : 1;
      unsigned int n2 = (dim>1) ? nsubp : 1;
      unsigned int n3 = (dim>2) ? nsubp : 1;
      unsigned int n4 = (dim>3) ? nsubp : 1;
      patch.data.reinit(5, n1*n2*n3*n4);

      for (unsigned int i4=0; i4<n4; ++i4)
        for (unsigned int i3=0; i3<n3; ++i3)
          for (unsigned int i2=0; i2<n2; ++i2)
            for (unsigned int i1=0; i1<n1; ++i1)
              {
                const unsigned int i = i1+nsubp*(i2+nsubp*(i3+nsubp*i4));
                const float x1 = 1.*i1/nsub;
                const float x2 = 1.*i2/nsub;
                const float x3 = 1.*i3/nsub;
                const float x4 = 1.*i4/nsub;

                patch.data(0,i) = p+x1;
                patch.data(1,i) = p+x2;
                patch.data(2,i) = p+x3;
                patch.data(3,i) = p+x4;
                patch.data(4,i) = i;
              }
      patch.patch_index = p;
    }
}

// Do this only if the necessary headers were included

#if defined(__deal2__quadrature_lib_h) && defined(__deal2__function_lib_h)

template <int dim>
void
create_continuous_patches(
  std::vector<DataOutBase::Patch<dim, dim> > &patches,
  unsigned int n_cells,
  unsigned int n_sub)
{
  unsigned int n1 = (dim>=1) ? n_cells : 1;
  unsigned int n2 = (dim>=2) ? n_cells : 1;
  unsigned int n3 = (dim>=3) ? n_cells : 1;

  QTrapez<dim> trapez;
  QTrapez<1> trapez1d;
  QIterated<dim> trapezsub(trapez1d, n_sub);

  Point<dim> midpoint;
  for (unsigned int d=0; d<dim; ++d)
    midpoint(d) = n_cells/2.;

  Functions::CutOffFunctionCinfty<dim> function(2., midpoint);

  for (unsigned int i3=0; i3<n3; ++i3)
    for (unsigned int i2=0; i2<n2; ++i2)
      for (unsigned int i1=0; i1<n1; ++i1)
        {
          DataOutBase::Patch<dim, dim> patch;
          patch.n_subdivisions = n_sub;
          for (unsigned int k=0; k<trapez.size(); ++k)
            {
              Point<dim> p = trapez.point(k);
              if (dim>=1) p(0) += i1;
              if (dim>=2) p(1) += i2;
              if (dim>=3) p(2) += i3;
              patch.vertices[k] = p;
            }
          std::vector<Point<dim> > points = trapezsub.get_points();
          for (unsigned int k=0; k<points.size(); ++k)
            {
              if (dim>=1) points[k](0) += i1;
              if (dim>=2) points[k](1) += i2;
              if (dim>=3) points[k](2) += i3;
            }
          std::vector<double> values(points.size());
          function.value_list(points, values);
          patch.data.reinit(1, points.size());
          for (unsigned int k=0; k<points.size(); ++k)
            patch.data(0,k) = values[k];
          patches.push_back(patch);
        }
}

#endif
