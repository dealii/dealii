//----------------------------  mg_smoother.templates.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mg_smoother.templates.h  ---------------------------
#ifndef __deal2__mg_dof_handler_templates_h
#define __deal2__mg_dof_handler_templates_h



template<class MATRIX, class VECTOR>
template<int dim>
MGSmootherRelaxation<MATRIX, VECTOR>
::MGSmootherRelaxation (const MGDoFHandler<dim>&     mg_dof,
			const MGLevelObject<MATRIX>& matrix,
			const function_ptr           relaxation,
			const unsigned int           steps,
			const double                 omega)
		:
		MGSmootherContinuous(mg_dof, steps),
		matrix(&matrix),
		relaxation(relaxation),
		omega(omega)
{};


template<class MATRIX, class VECTOR>
void
MGSmootherRelaxation<MATRIX, VECTOR>
::smooth (const unsigned int level,
	  VECTOR&            u,
	  const VECTOR&      rhs) const
{
  h1.reinit(u);
  h2.reinit(u);
  for(unsigned i=0;i<get_steps();++i)
    {
      (*matrix)[level].residual(h1, u, rhs);
      set_zero_interior_boundary(level,h1);
      ((*matrix)[level].*relaxation)(h2,h1,omega);
      set_zero_interior_boundary(level,h2);
      u.add(h2);
    }
}


#endif
