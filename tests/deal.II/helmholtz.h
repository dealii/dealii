// helmholtz.h,v 1.2 2000/02/09 20:01:34 wolf Exp

#include "../tests.h"
#include <multigrid/mg_base.h>

/**
 * Generator for system matrix and right hand side.
 */
class Helmholtz
{
  public:
    template<int dim, class MATRIX, class VECTOR>
    void build_all(MATRIX& A,
		   VECTOR& f,
		   const DoFHandler<dim>& dof,
		   const Quadrature<dim>& quadrature,
		   const Function<dim>& rhs);

    template<int dim, class MATRIX>
    void build_mgmatrix(MGLevelObject<MATRIX>& A,
		   const MGDoFHandler<dim>& dof,
		   const Quadrature<dim>& quadrature);
};

#include "helmholtz1.th"
#include "helmholtz1mg.th"

