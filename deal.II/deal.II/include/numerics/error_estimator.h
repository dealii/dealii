/*----------------------------   error_estimator.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __error_estimator_H
#define __error_estimator_H
/*----------------------------   error_estimator.h     ---------------------------*/


#include <base/exceptions.h>
#include <base/function.h>
#include <grid/forward_declarations.h>
#include <dofs/dof_handler.h>
#include <map>

				 // if multithreaded set number
				 // of threads to use.
				 // The default number of threads can
				 // be set during the compiling with
				 // the flag N_THREADS.
#ifdef DEAL_II_USE_MT
 #ifndef N_THREADS
  #define N_THREADS 4
 #endif
#else
 #ifdef N_THREADS
  #undef N_THREADS
 #endif
 #define N_THREADS 1
#endif


/**
 *  Implementation of the error estimator by Kelly, Gago, Zienkiewicz
 *  and Babuska.  This error estimator tries to approximate the error
 *  per cell by integration of the jump of the gradient of the
 *  solution along the faces of each cell.  It can be understood as a
 *  gradient recovery estimator; see the survey of Ainsworth for a
 *  complete discussion.
 *
 *  It seem as if this error estimator should only be valid for linear trial
 *  spaces, and there are indications that for higher order trial spaces the
 *  integrals computed here show superconvergence properties, i.e. they tend
 *  to zero faster than the error itself, thus ruling out the values as error
 *  indicators.
 *
 *  The error estimator really only estimates the error for the generalized
 *  Poisson equation $-\nabla\cdot a(x) \nabla u = f$ with either Dirichlet
 *  boundary conditions or generalized Neumann boundary conditions involving
 *  the conormal derivative $a\frac{du}{dn} = g$.
 *
 *  The error estimator returns a vector of estimated errors per cell which
 *  can be used to feed the #Triangulation<dim>::refine_*# functions. This
 *  vector contains elements of data type #float#, rather than #double#,
 *  since accuracy is not so important here, and since this can save rather
 *  a lot of memory, when using many cells.
 *
 *  
 *  \subsection{Implementation}
 *
 *  In principle, the implementation of the error estimation is simple: let
 *  $$ \eta_K^2 =
 *  \frac h{24} \int_{\partial K} \left[a \frac{\partial u_h}{\partial n}\right]^2 do
 *  $$
 *  be the error estimator for cell $K$. $[\cdot]$ denotes the jump of the
 *  argument at the face. In the paper of Ainsworth, $h$ is divided by $24$,
 *  but this factor is a bit esoteric, stemming from interpolation estimates
 *  and stability constants which may hold for the Poisson problem, but may
 *  not hold for more general situations. In the implementation, this factor
 *  is considered, but may lead to wrong results. You may scale the vector
 *  appropriately afterwards.
 *
 *  To perform the integration, use is made of the #FEFaceValues# and
 *  #FESubfaceValues# classes. The integration is performed by looping
 *  over all cells and integrating over faces that are not yet treated.
 *  This way we avoid integration on faces twice, once for each time we
 *  visit one of the adjacent cells. In a second loop over all cells, we
 *  sum up the contributions of the faces (which are the integrated
 *  square of the jumps) of each cell and take the square root.
 *
 *  The integration is done using a quadrature formula on the face.
 *  For linear trial functions (#FELinear#), the #Gauss2# or even the
 *  #Midpoint# rule will suffice. For higher order elements, it is
 *  necessary to utilize higher order quadrature formulae as well.
 *
 *  We store the contribution of each face in a #map#, as provided by the
 *  C++ standard library, with the iterator pointing to that face being the
 *  key into the map. In fact, we do not store the indicator per face, but
 *  only the integral listed above. When looping the second time over all
 *  cells, we have to sum up the contributions of the faces, multiply them
 *  with $\frac h{24}$ and take the square root. By doing the multiplication
 *  with $h$ in the second loop, we avoid problems to decide with which $h$
 *  to multiply, that of the cell on the one or that of the cell on the other
 *  side of the face.
 *
 *  $h$ is taken to be the greatest length of the diagonals of the cell. For
 *  more or less uniform cells without deformed angles, this coincides with
 *  the diameter of the cell.
 *  
 *
 *  \subsection{Vector-valued functions}
 *
 *  If the finite element field for which the error is to be estimated
 *  is vector-valued, i.e. the finite element has more than one
 *  component, then all the above can be applied to all or only some
 *  components at the same time. The main function of this class takes
 *  a list of flags denoting the components for which components the
 *  error estimator is to be applied; by default, it is a list of only
 *  #true#s, meaning that all variables shall be treated.
 *
 *  In case the different components of a field have different
 *  physical meaning (for example velocity and pressure in the Stokes
 *  equations), it would be senseless to use the same coefficient for
 *  all components. In that case, you can pass a function with as many
 *  components as there are components in the finite element field and
 *  each component of the error estimator will then be weighted by the
 *  respective component in this coefficient function. In the other
 *  case, when all components have the same meaning (for example the
 *  displacements in Lame's equations of elasticity), you can specify
 *  a scalar coefficient which will then be used for all components.
 *
 *
 *  \subsection{Boundary values}
 *  
 *  If the face is at the boundary, i.e. there is no neighboring cell to which
 *  the jump in the gradiend could be computed, there are two possibilities:
 *  \begin{itemize}
 *  \item The face belongs to a Dirichlet boundary. Then the face is not
 *    considered, which can be justified looking at a dual problem technique and
 *    should hold exactly if the boundary can be approximated exactly by the
 *    finite element used (i.e. it is a linear boundary for linear finite elements,
 *    quadratic for isoparametric quadratic elements, etc). For boundaries which
 *    can not be exactly approximated, one should consider the difference
 *    $z-z_h$ on the face, $z$ being a dual problem's solution which is zero at
 *    the true boundary and $z_h$ being an approximation, which in most cases
 *    will be zero on the numerical boundary. Since on the numerical boundary
 *    $z$ will not be zero in general, we would get another term here, but this
 *    one is neglected for practical reasons, in the hope that the error made
 *    here will tend to zero faster than the energy error we wish to estimate.
 *
 *    Though no integration is necessary, in the list of face contributions we
 *    store a zero for this face, which makes summing up the contributions of
 *    the different faces to the cells easier.
 *
 *  \item The face belongs to a Neumann boundary.  In this case, the
 *    contribution of the face $F\in\partial K$ looks like
 *    $$ \int_F \left|g-a\frac{\partial u_h}{\partial n}\right|^2 ds $$
 *    where $g$ is the Neumann boundary function. If the finite element is
 *    vector-valued, then obviously the function denoting the Neumann boundary
 *    conditions needs to be vector-valued as well.
 *
 *  \item No other boundary conditions are considered.
 *  \end{itemize}
 *  The object describing the boundary conditions is obtained from the
 *  triangulation.
 *
 *  Thanks go to Franz-Theo Suttmeier for clarifications about boundary
 *  conditions.
 *
 *  
 *  \subsection{Handling of hanging nodes}
 *  
 *  The integration along faces with hanging nodes is quite tricky, since one
 *  of the elements has to be shifted one level up or down. See the
 *  documentation for the #FESubfaceValues# class for more information about
 *  technical issues regarding this topic.
 *
 *  In praxi, since we integrate over each face only once, we do this when we
 *  are on the coarser one of the two cells adjacent to a subface (a subface
 *  is defined to be the child of a face; seen from the coarse cell, it is a
 *  subface, while seen from the refined cell it is one of its faces). The
 *  reason is that finding neighborship information is a bit easier then, but
 *  that's all practical reasoning, nothing fundamental.
 *
 *  Since we integrate from the coarse side of the face, we have the mother
 *  face readily at hand and store the result of the integration over that
 *  mother face (being the sum of the integrals along the subfaces) in the
 *  abovementioned map of integrals as well. This consumes some memory more
 *  than needed, but makes the summing up of the face contributions to the
 *  cells easier, since then we have the information from all faces of all
 *  cells at hand and need not think about explicitely determining whether
 *  a face was refined or not. The same applies for boundary faces, see
 *  above.
 *  
 *  @author Wolfgang Bangerth, 1998, 1999; parallelization by Thomas Richter, 2000
 */
template <int dim>
class KellyErrorEstimator
{
  public:

				     /**
				      *	Declare a data type which denotes a
				      *	mapping between a boundary indicator
				      *	and the function denoting the boundary
				      *	values on this part of the boundary.
				      *	Only one boundary function may be given
				      *	for each boundary indicator, which is
				      *	guaranteed by the #map# data type.
				      */    
    typedef map<unsigned char,const Function<dim>*> FunctionMap;

    
				     /**
				      * Implementation of the error
				      * estimator described above. You
				      * may give a coefficient, but
				      * there is a default value which
				      * denotes the constant
				      * coefficient with value
				      * one. The coefficient function
				      * may either be a scalar one, in
				      * which case it is used for all
				      * components of the finite
				      * element, or a vector-valued
				      * one with as many components as
				      * there are in the finite
				      * element; in the latter case,
				      * each component is weighted by
				      * the respective component in
				      * the coefficient.
				      *
				      * You might give a list of
				      * components you want to
				      * evaluate, in case the finite
				      * element used by the
				      * #DoFHandler# object is
				      * vector-valued. You then have
				      * to set those entries to true
				      * in the bit-vector
				      * #component_mask# for which the
				      * respective component is to be
				      * used in the error
				      * estimator. The default is to
				      * use all components, which is
				      * done by either providing a
				      * bit-vector with all-set
				      * entries, or an empty
				      * bit-vector.
				      *
				      * The estimator supports
				      * multithreading and splits the
				      * cells to #N_THREADS# (default)
				      * threads. The number of threads
				      * to be used in multithreaded
				      * mode can be set with the last
				      * parameter of the error
				      * estimator.  Multithreading is
				      * only implemented in two and
				      * three dimensions.
				      */
    static void estimate (const DoFHandler<dim>   &dof,
			  const Quadrature<dim-1> &quadrature,
			  const FunctionMap       &neumann_bc,
			  const Vector<double>    &solution,
			  Vector<float>           &error,
			  const vector<bool>      &component_mask_ = vector<bool>(),
			  const Function<dim>     *coefficients   = 0,
			  unsigned int            n_threads=N_THREADS);
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidBoundaryIndicator);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidComponentMask);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidCoefficient);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidBoundaryFunction);
    


  private:


    				     /**
				      * Declare a data type to represent the
				      * mapping between faces and integrated
				      * jumps of gradients. See the general
				      * doc of this class for more information.
				      */
    typedef map<typename DoFHandler<dim>::face_iterator,double> FaceIntegrals;


				     /**
				      * Redeclare an active cell iterator.
				      * This is simply for convenience.
				      */
    typedef DoFHandler<dim>::active_cell_iterator active_cell_iterator;


    				     /**
				      * All data needed by the several functions
				      * of the error estimator is gathered in
				      * this struct. It is passed as a reference
				      * to the seperate functions in this class.
				      *
				      * The reason for invention of
				      * this object is two-fold:
				      * first, class member data is
				      * not possible because no real
				      * object is created (all
				      * functions are #static#), which
				      * is a historical
				      * reason. Second, if we don't
				      * collect the data the various
				      * functions need somewhere at a
				      * central place, that would mean
				      * that the functions would have
				      * to allocate them upon
				      * need. However, then some
				      * variables would be allocated
				      * over and over again, which can
				      * take a significant amount of
				      * time (10-20 per cent) and most
				      * importantly, memory allocation
				      * requires synchronisation in
				      * multithreaded mode. While that
				      * is done by the C++ library and
				      * has not to be handcoded, it
				      * nevertheless seriously damages
				      * tha ability to efficiently run
				      * the functions of this class in
				      * parallel, since they are quite
				      * often blocked by these
				      * synchronisation points.
				      */
    struct Data
    {
	const DoFHandler<dim>   &dof;
	const Quadrature<dim-1> &quadrature;
	const FunctionMap       &neumann_bc;
	const Vector<double>    &solution;
	vector<bool>            component_mask;
	const Function<dim>     *coefficients;
	unsigned int            n_threads;
    
	DoFHandler<dim>::active_cell_iterator endc;
	unsigned int            n_components;
	unsigned int            n_q_points;
	FaceIntegrals           face_integrals;  

					 /**
					  * A vector to store the jump of the
					  * normal vectors in the quadrature
					  * points.
					  * There is one vector for every
					  * thread used in the estimator.
					  * The allocation of memory has to
					  * be global to enable fast use
					  * of multithreading	  
					  */ 
	vector< vector<vector<double> > >         phi;

					 /**
					  * A vector for the gradients of
					  * the finite element function
					  * on one cell
					  *
					  * Let psi be a short name for
					  * #a grad u_h#, where the second
					  * index be the component of the
					  * finite element, and the first
					  * index the number of the
					  * quadrature point.
					  */
	vector< vector<vector<Tensor<1,dim> > > > psi;

					 /**
					  * The same vector for a neighbor cell
					  */
	vector< vector<vector<Tensor<1,dim> > > > neighbor_psi;

					 /**
					  * The normal vectors of the finite
					  * element function on one face
					  */
	vector< vector<Point<dim> > > normal_vectors;

					 /**
					  * Two arrays needed for the
					  * values of coefficients in
					  * the jumps, if they are
					  * given.
					  */
	vector< vector<double> >          coefficient_values1;
	vector< vector<Vector<double> > > coefficient_values;

					 /**
					  * Array for the products of
					  * Jacobian determinants and
					  * weights of quadraturs
					  * points.
					  */
	vector< vector<double> >          JxW_values;

					 /**
					  * A constructor of the
					  * class Data. All variables are
					  * passed as references.
					  */
	Data(const DoFHandler<dim>   &dof,
	     const Quadrature<dim-1> &quadrature,
	     const FunctionMap       &neumann_bc,
	     const Vector<double>    &solution,
	     vector<bool>             component_mask_,
	     const Function<dim>     *coefficients,
	     unsigned int             n_threads);
    };

    
				     /**
				      * Computates the error on all cells
				      * of the domain with the number n,
				      * satisfying
				      * #n=this_thread (mod n_threads)#
				      * This enumeration is chosen to
				      * generate a random distribution
				      * of all cells.
				      *
				      * This function is only needed
				      * in two or three dimensions.
				      * The error estimator in one
				      * dimension is implemented
				      * seperatly.
				      */
    static void estimate_some (Data &data,
			       const unsigned int this_thread);
    				
				     /**
				      * Actually do the computation on a face
				      * which has no hanging nodes (it is
				      * regular),  i.e.
				      * either on the other side there is
				      * nirvana (face is at boundary), or
				      * the other side's refinement level
				      * is the same as that of this side,
				      * then handle the integration of
				      * these both cases together.
				      *
				      * The meaning of the parameters becomes
				      * clear when looking at the source
				      * code. This function is only
				      * externalized from #estimate_error#
				      * to avoid ending up with a function
				      * of 500 lines of code.
				      */


    static void integrate_over_regular_face (Data                       &data,
					     const unsigned int          this_thread,
					     const active_cell_iterator &cell,
					     const unsigned int          face_no,
					     FEFaceValues<dim>          &fe_face_values_cell,
					     FEFaceValues<dim>          &fe_face_values_neighbor);
    
    
				     /**
				      * The same applies as for the function
				      * above, except that integration is
				      * over face #face_no# of #cell#, where
				      * the respective neighbor is refined,
				      * so that the integration is a bit more
				      * complex.
				      */
    static void integrate_over_irregular_face (Data                       &data,
					       const unsigned int          this_thread,
					       const active_cell_iterator &cell,
					       const unsigned int          face_no,
					       FEFaceValues<dim>          &fe_face_values,
					       FESubfaceValues<dim>       &fe_subface_values);
};





/*----------------------------   error_estimator.h     ---------------------------*/
/* end of #ifndef __error_estimator_H */
#endif
/*----------------------------   error_estimator.h     ---------------------------*/


