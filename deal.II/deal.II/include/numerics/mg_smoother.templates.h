// $Id$

template<typename number>
template<int dim>
MGSmootherRelaxation<number>::
MGSmootherRelaxation(const MGDoFHandler<dim> &mg_dof,
		     const MGMatrix<SparseMatrix<number> >& matrix,
		     function_ptr relaxation,
		     unsigned int steps,
		     double omega)
		:
		MGSmoother(mg_dof, steps),
		matrix(&matrix),
		relaxation(relaxation),
		omega(omega)
{}

template<typename number>
void
MGSmootherRelaxation<number>::smooth (const unsigned int   level,
		       Vector<double>       &u,
		       const Vector<double> &rhs) const
{
  h1.reinit(u);
  h2.reinit(u);
  for(unsigned i=0;i<get_steps();++i)
    {
      (*matrix)[level].residual(h1, u, rhs);
      ((*matrix)[level].*relaxation)(h2,h1,omega);
      set_zero_interior_boundary(level,h2);
      u.add(h2);
    }
}

