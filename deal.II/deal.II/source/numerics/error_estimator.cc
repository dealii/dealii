/*      $Id$                 */


#include <numerics/error-estimator.h>
#include <grid/dof.h>
#include <lac/dvector.h>



void KellyErrorEstimator<1>::estimate_error (const DoFHandler<1> &,
					     const dVector &,
					     dVector &) const {
  Assert(false, ExcNotImplemented());
};



template <int dim>
void KellyErrorEstimator<dim>::estimate_error (const DoFHandler<dim> &dof,
					       const dVector         &solution,
					       dVector               &error) const {
  ;
};

