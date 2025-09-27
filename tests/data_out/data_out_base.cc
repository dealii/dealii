// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2003 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/data_out_base.h>

#include <string>
#include <vector>

#include "../tests.h"

// TODO: Several functions are commented out since implementations are missing

#define WRITE(type)                                                       \
  {                                                                       \
    out << "Writing " #type " dimensions " << dim << ',' << spacedim      \
        << std::endl;                                                     \
    DataOutBase::write_##type(patches, names, vectors, type##flags, out); \
  }

template <int dim, int spacedim>
void
write_patches(const std::vector<DataOutBase::Patch<dim, spacedim>> &patches,
              std::ostream                                         &out)
{
  std::vector<std::string> names(2);
  names[0] = std::string("first name");
  names[1] = std::string("last name");

  DataOutBase::DXFlags                   dxflags;
  DataOutBase::EpsFlags                  epsflags;
  DataOutBase::GnuplotFlags              gnuplotflags;
  DataOutBase::GmvFlags                  gmvflags;
  DataOutBase::PovrayFlags               povrayflags;
  DataOutBase::UcdFlags                  ucdflags(true);
  DataOutBase::VtkFlags                  vtkflags;
  DataOutBase::Deal_II_IntermediateFlags deal_II_intermediateflags;

  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    vectors;

  WRITE(dx);
  if (dim == 2)
    WRITE(eps);
  WRITE(gnuplot);
  WRITE(gmv);
  if (dim == 2 && spacedim == 2)
    WRITE(povray);
  WRITE(ucd);
  WRITE(vtk);
  WRITE(deal_II_intermediate);
  //  if (dim==2)
  //    WRITE(eps);
}

template <int dim>
struct PatchInfo
{
  static double vertices[GeometryInfo<dim>::vertices_per_cell][3];
  static double offsets[GeometryInfo<dim>::vertices_per_cell][3];
  static int    neighbors[GeometryInfo<dim>::vertices_per_cell]
                      [GeometryInfo<dim>::faces_per_cell];
};


template <>
double PatchInfo<1>::vertices[2][3] = {{0, 0, 0}, {3, 6, 9}};

template <>
double PatchInfo<2>::vertices[4][3] = {{0, 0, 0},
                                       {3, 0, 3},
                                       {3, 0, 3},
                                       {3, 3, 6}};


template <>
double PatchInfo<3>::vertices[8][3] = {{0, 0, 0},
                                       {3, 0, 0},
                                       {0, 3, 0},
                                       {3, 3, 0},
                                       {0, 0, 3},
                                       {3, 0, 3},
                                       {0, 3, 3},
                                       {3, 3, 3}};


template <>
int PatchInfo<1>::neighbors[2][2] = {{-1, 1}, {0, -1}};



template <>
int PatchInfo<2>::neighbors[4][4] = {{-1, 1, 2, -1},
                                     {-1, -1, 3, 0},
                                     {0, 3, -1, -1},
                                     {1, -1, -1, 2}};



template <>
int PatchInfo<3>::neighbors[8][6] = {{-1, -1, -1, -1, -1, -1},
                                     {-1, -1, -1, -1, -1, -1},
                                     {-1, -1, -1, -1, -1, -1},
                                     {-1, -1, -1, -1, -1, -1},
                                     {-1, -1, -1, -1, -1, -1},
                                     {-1, -1, -1, -1, -1, -1},
                                     {-1, -1, -1, -1, -1, -1},
                                     {-1, -1, -1, -1, -1, -1}};


template <>
double PatchInfo<1>::offsets[2][3] = {{0, 0, 0}, {3, 0, 0}};


template <>
double PatchInfo<2>::offsets[4][3] = {{0, 0, 0},
                                      {3, 0, 0},
                                      {0, 3, 0},
                                      {3, 3, 0}};


template <>
double PatchInfo<3>::offsets[8][3] = {{0, 0, 0},
                                      {3, 0, 0},
                                      {0, 3, 0},
                                      {3, 3, 0},
                                      {0, 0, 3},
                                      {3, 0, 3},
                                      {0, 3, 3},
                                      {3, 3, 3}};



template <int dim, int spacedim>
void
create_patches(std::vector<DataOutBase::Patch<dim, spacedim>> &patches)
{
  const unsigned int ncells = GeometryInfo<dim>::vertices_per_cell;
  const unsigned int nsub   = 3;

  patches.resize(ncells);

  for (unsigned int c = 0; c < ncells; ++c)
    {
      DataOutBase::Patch<dim, spacedim> &p = patches[c];
      p.patch_index                        = c;
      p.n_subdivisions                     = nsub;
      p.reference_cell = ReferenceCells::get_hypercube<dim>();

      for (unsigned int i = 0; i < ncells; ++i)
        for (unsigned int j = 0; j < spacedim; ++j)
          p.vertices[i][j] =
            PatchInfo<dim>::vertices[i][j] + PatchInfo<dim>::offsets[c][j];

      for (const unsigned int i : GeometryInfo<dim>::face_indices())
        p.neighbors[i] = (unsigned int)PatchInfo<dim>::neighbors[c][i];

      unsigned int ndata = 1;
      for (unsigned int i = 0; i < dim; ++i)
        ndata *= nsub + 1;

      p.data.reinit(2, ndata);
      for (unsigned int i = 0; i < ndata; ++i)
        {
          p.data(0, i) = i;
          p.data(1, i) = -(float)i;
        }
    }
}


template <int dim, int spacedim>
void
test(std::ostream &out)
{
  std::vector<DataOutBase::Patch<dim, spacedim>> patches;
  create_patches(patches);
  write_patches(patches, out);
}


int
main()
{
  initlog();
  auto &logfile = deallog.get_file_stream();

  test<1, 2>(logfile);
  test<2, 2>(logfile);
  test<2, 3>(logfile);
  test<3, 3>(logfile);
}
