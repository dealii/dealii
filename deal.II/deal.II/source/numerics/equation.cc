// $Id$

#include <numerics/assembler.h>
#include <grid/tria_iterator.h>
#include <grid/tria_iterator.templates.h>
#include <fe/fe.h>
#include <lac/fullmatrix.h>
#include <lac/vector.h>
#include <lac/sparsematrix.h>

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
