//----------------------------  mapping_q1_eulerian.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001 by the deal.II authors and Michael Stadler
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mapping_q1_eulerian.h  ---------------------------
#ifndef __deal2__mapping_q1_eulerian_h
#define __deal2__mapping_q1_eulerian_h

#include <base/config.h>
#include <base/smartpointer.h>
#include <fe/mapping_q1.h> 


/**
 * Eulerian mapping of general unit cells by d-linear shape
 * functions. Each cell is thus shifted in space by values given to
 * the mapping through a finite element field.
 *
 * @sect3{Usage} 
 *
 * The constructor of this class takes two arguments: a reference to
 * the vector that defines the mapping from the reference
 * configuration to the current configuration and a reference to the
 * @ref{DoFHandler}. The vector should then represent a (flattened out
 * version of a) vector valued field defined at nodes defined by the 
 * the @ref{DoFHandler}, where the number of components of the vector
 * field equals the number of space dimensions. Thus, the 
 * @ref{DoFHandler} shall operate on a finite element that has as many 
 * components as space dimensions. As an additional requirement, we 
 * impose that it have as many degree of freedom per vertex as there
 * are space dimensions; since this object only evaluates the finite 
 * element field at the vertices, the values
 * of all other degrees of freedom (not associated to vertices) are
 * ignored. These requirements are met if the finite element which the
 * given @ref{DoFHandler} operates on is constructed as a system
 * element (@ref{FESystem}) from @p{dim} continuous @ref{FE_Q}
 * objects.
 *
 * In many cases, the shift vector will also be the solution vector of
 * the problem under investigation. If this is not the case (i.e. the
 * number of components of the solution variable is not equal to the
 * space dimension, e.g. for scalar problems in @p{dim>1} where the
 * Eulerian coordinates only give a background field) or for coupled
 * problems where more variables are computed than just the flow
 * field), then a different @ref{DoFHandler} has to be set up on the
 * given triangulation, and the shift vector has then to be associated
 * to it.
 *
 * An example is shown below:
 * @begin{verbatim}
 *    FESystem<dim> fe(FE_Q<dim>(1), dim);
 *    DoFHandler<dim> flowfield_dof_handler(triangulation);
 *    flowfield_dof_handler.distribute_dofs(fe);
 *    Vector<double> map_points(flowfield_dof_handler.n_dofs());
 *    MappingQ1Eulerian<dim> mymapping(map_points, flowfield_dof_handler);
 * @end{verbatim}
 *
 * Note that since the vector of shift values and the dof handler are
 * only associated to this object at construction time, you have to
 * make sure that whenever you use this object, the given objects
 * still represent valid data.
 *
 * @author Michael Stadler, 2001
 */
template <int dim>
class MappingQ1Eulerian : public MappingQ1<dim>
{
  public:

				     /**
				      * Constructor. It takes a
				      * @p{Vector<double> &} as its
				      * first argument to specify the
				      * transformation of the whole
				      * problem from the reference to
				      * the current configuration.
				      * The organization of the
				      * elements in the @p{Vector}
				      * must follow the concept how
				      * deal.II stores solutions that
				      * are associated to a
				      * triangulation.  This is
				      * automatically the case if the
				      * @p{Vector} represents the
				      * solution of the previous step
				      * of a nonlinear problem.
				      * Alternatively, the @p{Vector}
				      * can be initialized by
				      * @p{DoFObjectAccessor::set_dof_values()}.
				      */
    MappingQ1Eulerian ( const Vector<double>  &euler_transform_vectors,
		        const DoFHandler<dim> &shiftmap_dof_handler);


				     /**
				      * Exception
				      */
    DeclException0 (ExcWrongNoOfComponents);
    
				     /**
				      * Exception.
				      */
    DeclException2 (ExcWrongVectorSize,
		    int, int,
		    << "Vector has wrong size " << arg1
		    << ", expected size " << arg2);

				     /**
				      * Exception.
				      */
    DeclException0 (ExcInactiveCell);



  protected:

				     /**
				      * Reference to the vector of
				      * shifts.
				      */
    const Vector<double> &euler_transform_vectors;
    
                                     /**
				      * Pointer to the DoFHandler to
				      * which the mapping vector is
				      * associated.
				      */
    const SmartPointer<const DoFHandler<dim> > shiftmap_dof_handler;
    

  private:    
				     /**
				      * Computes the support points of
				      * the mapping. For
				      * @p{MappingQ1Eulerian} these
				      * are the vertices.
				      */
    virtual void compute_mapping_support_points(
      const typename Triangulation<dim>::cell_iterator &cell,
      typename std::vector<Point<dim> > &a) const;
    
};


#endif
