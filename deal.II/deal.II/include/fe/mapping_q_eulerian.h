//---------------------------------------------------------------------------
//    $Id: mapping_q.h 15711 2008-02-05 20:43:14Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__mapping_q_eulerian_h
#define __deal2__mapping_q_eulerian_h

#include <base/smartpointer.h>
#include <base/thread_management.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <fe/mapping_q.h> 


DEAL_II_NAMESPACE_OPEN


/*!@addtogroup mapping */
/*@{*/

/**
 * This class is an extension of the MappingQ1Eulerian
 * class to higher order Qp mappings.  It is useful
 * when one wants to calculate shape function information on
 * a domain that is deforming as the computation proceeds.
 *
 * <h3>Usage</h3> 
 *
 * The constructor of this class takes three arguments: the polynomial
 * degree of the desire Qp mapping, a reference to
 * the vector that defines the mapping from the initial
 * configuration to the current configuration, and a reference to the
 * DoFHandler. The most common case is to use the solution 
 * vector for the problem under consideration as the shift vector.  
 * The key reqirement is that the number of components 
 * of the given vector field be equal to (or possibly greater than) the 
 * number of space dimensions. If there are more components than space 
 * dimensions (for example, if one is working with a coupled problem
 * where there are additional solution variables), the
 * first <tt>dim</tt> components are assumed to represent the displacement 
 * field, and the remaining components are ignored.  If this assumption
 * does not hold one may need to set up a separate DoFHandler on
 * the triangulation and associate the desired shift vector to it.
 *
 * Typically, the DoFHandler operates on a finite element that
 * is constructed as a system element (FESystem) from continuous FE_Q()
 * objects. An example is shown below:
 * @verbatim
 *    FESystem<dim> fe(FE_Q<dim>(2), dim, FE_Q<dim>(1), 1);
 *    DoFHandler<dim> dof_handler(triangulation);
 *    dof_handler.distribute_dofs(fe);
 *    Vector<double> soln_vector(dof_handler.n_dofs());
 *    MappingQEulerian<dim> q2_mapping(2,soln_vector,dof_handler);
 * @endverbatim
 *
 * In this example, our element consists of <tt>(dim+1)</tt> components.
 * Only the first <tt>dim</tt> components will be used, however, to define
 * the Q2 mapping.  The remaining components are ignored.
 *
 * Note that it is essential to call the distribute_dofs(...) function
 * before constructing a mapping object.
 *
 * Also note that since the vector of shift values and the dof handler are
 * only associated to this object at construction time, you have to
 * make sure that whenever you use this object, the given objects
 * still represent valid data.
 *
 * To enable the use of the MappingQ1Eulerian class also in the context
 * of parallel codes using the PETSc wrapper classes, the type of
 * the vector can be specified as template parameter <tt>EulerVectorType</tt>
 * Not specifying this template argument in applications using the PETSc
 * vector classes leads to the construction of a copy of the vector
 * which is not acccessible afterwards!
 *
 * @author Joshua White, 2008
 */
template <int dim, class EulerVectorType = Vector<double>, int spacedim=dim >
class MappingQEulerian : public MappingQ<dim, spacedim>
{
  public:
                                     /**
                                      * Constructor. The first argument is
                                      * the polynomical degree of the desired
                                      * Qp mapping.  It then takes a
                                      * <tt>Vector<double> &</tt> to specify the
                                      * transformation of the domain
                                      * from the reference to
                                      * the current configuration.
                                      * The organization of the
                                      * elements in the @p Vector
                                      * must follow the concept how
                                      * deal.II stores solutions that
                                      * are associated to a
                                      * triangulation.  This is
                                      * automatically the case if the
                                      * @p Vector represents the
                                      * solution of the previous step
                                      * of a nonlinear problem.
                                      * Alternatively, the @p Vector
                                      * can be initialized by
                                      * <tt>DoFAccessor::set_dof_values()</tt>.
                                      */

    MappingQEulerian (const unsigned int     degree,
                      const EulerVectorType  &euler_vector,
                      const DoFHandler<dim>  &euler_dof_handler);

                                     /**
                                      * Return a pointer to a copy of the
                                      * present object. The caller of this
                                      * copy then assumes ownership of it.
                                      */
    virtual
    Mapping<dim,spacedim> * clone () const;

                                     /**
                                      * Exception
                                      */

    DeclException0 (ExcWrongNoOfComponents);

				     /**
                                      * Exception
                                      */

    DeclException0 (ExcInactiveCell);
    
				     /**
                                      * Exception
                                      */

    DeclException2 (ExcWrongVectorSize, int, int,
		    << "Vector has wrong size " << arg1
		    << "-- expected size " << arg2);

  protected:
				     /**
				      * Implementation of the interface in
				      * MappingQ. Overrides the function in
				      * the base class, since we cannot use
				      * any cell similarity for this class.
				      */
    virtual void
    fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
		    const Quadrature<dim>                                     &quadrature,
		    typename Mapping<dim,spacedim>::InternalDataBase          &mapping_data,
		    typename std::vector<Point<spacedim> >                    &quadrature_points,
		    std::vector<double>                                       &JxW_values,
		    std::vector<Tensor<2,spacedim> >                          &jacobians,
		    std::vector<Tensor<3,spacedim> >                          &jacobian_grads,
		    std::vector<Tensor<2,spacedim> >                          &inverse_jacobians,
		    std::vector<Point<spacedim> >                             &cell_normal_vectors,
		    enum CellSimilarity::Similarity                           &cell_similarity) const;

                                     /**
                                      * Reference to the vector of
                                      * shifts.
                                      */

    const EulerVectorType &euler_vector;
    
                                     /**
                                      * Pointer to the DoFHandler to
                                      * which the mapping vector is
                                      * associated.
                                      */

    const SmartPointer<const DoFHandler<dim> > euler_dof_handler;


  private:   
 
                                     /**
                                      * Special quadrature rule used
                                      * to define the support points
                                      * in the reference configuration.
                                      */

    class SupportQuadrature : public Quadrature<dim>
    { 
      public:
					 /**
					  * Constructor, with an argument
					  * defining the desired polynomial
					  * degree.
					  */

        SupportQuadrature (const unsigned int map_degree); 

    };

				     /**
				      * A member variable holding the
				      * quadrature points in the right
				      * order.
				      */
    const SupportQuadrature support_quadrature;

				     /**
                                      * FEValues object used to query the
                                      * the given finite element field
                                      * at the support points in the
                                      * reference configuration.
				      *
				      * The variable is marked as
				      * mutable since we have to call
				      * FEValues::reinit from
				      * compute_mapping_support_points,
				      * a function that is 'const'.
                                      */
    mutable FEValues<dim> fe_values;

				     /**
				      * A variable to guard access to
				      * the fe_values variable.
				      */
    mutable Threads::ThreadMutex fe_values_mutex;
    
                                     /**
                                      * Compute the positions of the
                                      * support points in the current
                                      * configuration
                                      */
    virtual void compute_mapping_support_points(
      const typename Triangulation<dim>::cell_iterator &cell,
      std::vector<Point<dim> > &a) const;
    
};

/*@}*/


DEAL_II_NAMESPACE_CLOSE


#endif // __deal2__mapping_q_eulerian_h

