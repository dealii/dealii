// $Id$


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

