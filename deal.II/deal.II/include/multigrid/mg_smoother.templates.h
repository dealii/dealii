// $Id$


template<typename number>
template<int dim>
MGSmootherRelaxation<number>::MGSmootherRelaxation (const MGDoFHandler<dim>                       &mg_dof,
		      const MGLevelObject<SparseMatrix<number> > &matrix,
		      const function_ptr                             relaxation,
		      const unsigned int                             steps,
		      const double                                   omega)
		:
		MGSmoother(mg_dof, steps),
		matrix(&matrix),
		relaxation(relaxation),
		omega(omega)
{};


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
      set_zero_interior_boundary(level,h1);
      ((*matrix)[level].*relaxation)(h2,h1,omega);
      set_zero_interior_boundary(level,h2);
      u.add(h2);
    }
}

