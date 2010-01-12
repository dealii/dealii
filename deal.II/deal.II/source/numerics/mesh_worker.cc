//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <multigrid/mg_tools.h>
#include <fe/fe.h>
#include <fe/fe_tools.h>
#include <numerics/mesh_worker.h>
#include <base/quadrature_lib.h>

using namespace dealii;
using namespace MeshWorker;

template <int dim>
IntegrationWorker<dim>::IntegrationWorker ()
{
  cell_flags = update_JxW_values;
  bdry_flags = UpdateFlags(update_JxW_values | update_normal_vectors);
  face_flags = bdry_flags;
  ngbr_flags = update_default;
}


template<int dim>
void
IntegrationWorker<dim>::initialize_update_flags ()
{
  if (cell_selector.has_values() != 0) cell_flags |= update_values;
  if (cell_selector.has_gradients() != 0) cell_flags |= update_gradients;
  if (cell_selector.has_hessians() != 0) cell_flags |= update_hessians;
  
  if (bdry_selector.has_values() != 0) bdry_flags |= update_values;
  if (bdry_selector.has_gradients() != 0) bdry_flags |= update_gradients;
  if (bdry_selector.has_hessians() != 0) bdry_flags |= update_hessians;
  
  if (face_selector.has_values() != 0) face_flags |= update_values;
  if (face_selector.has_gradients() != 0) face_flags |= update_gradients;
  if (face_selector.has_hessians() != 0) face_flags |= update_hessians;
  
  if (face_selector.has_values() != 0) ngbr_flags |= update_values;
  if (face_selector.has_gradients() != 0) ngbr_flags |= update_gradients;
  if (face_selector.has_hessians() != 0) ngbr_flags |= update_hessians;  
}


template<int dim>
void
IntegrationWorker<dim>::add_update_flags(
  const UpdateFlags flags, bool cell, bool bdry, bool face, bool ngbr)
{
  if (cell) cell_flags |= flags;
  if (bdry) bdry_flags |= flags;
  if (face) face_flags |= flags;
  if (ngbr) ngbr_flags |= flags;  
}

  
template<int dim>
void
IntegrationWorker<dim>::initialize_gauss_quadrature(
  unsigned int cp,
  unsigned int bp,
  unsigned int fp)
{
  cell_quadrature = QGauss<dim>(cp);
  bdry_quadrature = QGauss<dim-1>(bp);
  face_quadrature = QGauss<dim-1>(fp);

}

#if deal_II_dimension > 1
template class IntegrationWorker<deal_II_dimension>;
#endif

