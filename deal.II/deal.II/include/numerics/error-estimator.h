/*----------------------------   error-estimator.h     ---------------------------*/
/*      $Id$                 */
#ifndef __error-estimator_H
#define __error-estimator_H
/*----------------------------   error-estimator.h     ---------------------------*/


template <int dim>
class KellyErrorEstimator {
  public:
    void estimate_error (const DoFHandler<dim> &dof,
			 const dVector         &solution,
			 const dVector         &error) const;
};



/*----------------------------   error-estimator.h     ---------------------------*/
/* end of #ifndef __error-estimator_H */
#endif
/*----------------------------   error-estimator.h     ---------------------------*/
