//----------------------------  equation.cc  ---------------------------
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
//----------------------------  equation.cc  ---------------------------


#include <numerics/assembler.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>
#include <fe/fe.h>
#include <lac/full_matrix.h>
#include <lac/vector.h>
#include <lac/sparse_matrix.h>

template <int dim>
Equation<dim>::Equation (const unsigned int n_equations) :
		n_eq(n_equations) {};


template <int dim>
void Equation<dim>::assemble (FullMatrix<double>          &,
			      Vector<double>           &,
			      const FEValues<dim> &,
			      const typename DoFHandler<dim>::cell_iterator &) const
{
  Assert (false, ExcPureVirtualFunctionCalled());
};


template <int dim>
void Equation<dim>::assemble (FullMatrix<double>          &,
			      const FEValues<dim> &,
			      const typename DoFHandler<dim>::cell_iterator &) const
{
  Assert (false, ExcPureVirtualFunctionCalled());
};


template <int dim>
void Equation<dim>::assemble (Vector<double>           &,
			      const FEValues<dim> &,
			      const typename DoFHandler<dim>::cell_iterator &) const
{
  Assert (false, ExcPureVirtualFunctionCalled());
};


template class Equation<deal_II_dimension>;
