//----------------------------  fe_lib.dg.constant.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_lib.dg.constant.cc  ---------------------------


#include <fe/fe_lib.dg.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <grid/geometry_info.h>
#include <algorithm>

// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
using namespace std;
#endif




template <int dim>
FEDG_Q0<dim>::FEDG_Q0 () :
		FEQ1Mapping<dim> (0, 
				  (dim==1 ? 1 : 0),
				  (dim==2 ? 1 : 0),
				  (dim==3 ? 1 : 0),
				  1,
				  std::vector<bool> (1, true))
{
  for (unsigned int i=0; i<GeometryInfo<dim>::children_per_cell; ++i)
    { 
      restriction[i](0,0) = 1./GeometryInfo<dim>::children_per_cell;
      prolongation[i](0,0) = 1.0;
    }
};


#if deal_II_dimension == 1


template <>
void
FEDG_Q0<1>::get_face_support_points (const DoFHandler<1>::face_iterator &,
				     std::vector<Point<1> >  &) const
{
  Assert (false, ExcInternalError());
};

#endif



template <int dim>
inline
double
FEDG_Q0<dim>::shape_value (const unsigned int i,
			   const Point<dim>&) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));
  return 1.;
};



template <int dim>
inline
Tensor<1,dim>
FEDG_Q0<dim>::shape_grad (const unsigned int i,
			  const Point<dim>&) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));
  return Tensor<1,dim> ();
};



template <int dim>
inline
Tensor<2,dim>
FEDG_Q0<dim>::shape_grad_grad (const unsigned int i,
			       const Point<dim> &) const
{
  Assert((i<dofs_per_cell), ExcIndexRange(i, 0, dofs_per_cell));

  return Tensor<2,dim>();
};



template <int dim>
void FEDG_Q0<dim>::get_local_mass_matrix (const typename DoFHandler<dim>::cell_iterator &cell,
					  FullMatrix<double> &local_mass_matrix) const
{
  Assert (local_mass_matrix.n() == dofs_per_cell,
	  FiniteElementBase<dim>::ExcWrongFieldDimension(local_mass_matrix.n(),
							 dofs_per_cell));
  Assert (local_mass_matrix.m() == dofs_per_cell,
	  FiniteElementBase<dim>::ExcWrongFieldDimension(local_mass_matrix.m(),
							 dofs_per_cell));

  local_mass_matrix(0,0) = cell->measure();
};



template <int dim>
void
FEDG_Q0<dim>::get_unit_support_points (typename std::vector<Point<dim> > &unit_points) const
{
  Assert (unit_points.size() == dofs_per_cell,
	  FiniteElementBase<dim>::ExcWrongFieldDimension (unit_points.size(), dofs_per_cell));
  for (unsigned int d=0; d<dim; ++d)
    unit_points[0](d) = 0.5;
};



template <int dim>
void
FEDG_Q0<dim>::get_support_points (const typename DoFHandler<dim>::cell_iterator &cell,
				  typename std::vector<Point<dim> >  &support_points) const
{
  Assert (support_points.size() == dofs_per_cell,
	  FiniteElementBase<dim>::ExcWrongFieldDimension (support_points.size(), dofs_per_cell));
  
  support_points[0] = cell->center();
};



template <int dim>
void
FEDG_Q0<dim>::get_face_support_points (const typename DoFHandler<dim>::face_iterator &,
				       typename std::vector<Point<dim> >  &support_points) const
{
  Assert ((support_points.size() == 0),
	  FiniteElementBase<dim>::ExcWrongFieldDimension (support_points.size(),0));
};


// explicit instantiations

template class FEDG_Q0<deal_II_dimension>;
