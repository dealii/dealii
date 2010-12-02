//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__data_postprocessor_h
#define __deal2__data_postprocessor_h



#include <base/subscriptor.h>
#include <base/tensor.h>
#include <lac/vector.h>
#include <fe/fe_update_flags.h>
#include <numerics/data_component_interpretation.h>

#include <vector>
#include <string>

DEAL_II_NAMESPACE_OPEN


/**
 * For the (graphical) output of a FE solution one frequently wants to include
 * derived quantities, which are calculated from the values of the solution
 * and possibly the first and second derivates of the solution. Examples are
 * the calculation Mach numbers from velocity and density in supersonic flow
 * computations, or the computation of the magnitude of a complex-valued
 * solution as demonstrated in step-29. Other uses are shown in
 * step-33. This class offers the interface to perform such
 * postprocessing. Given the values and derivatives of the solution at those
 * points where we want to generated output, the functions of this class can
 * be overloaded to compute new quantities.
 *
 * A data vector and an object of a derived class can be given to the
 * <tt>DataOut::add_data_vector</tt> function, which will write the derived
 * quantities instead of the provided data to the output file. Note, that the
 * DataPostprocessor has to live until @p build_patches has been
 * called. DataOutFaces and DataOutRotation can be used as well.
 *
 * In order not to perform needless calculations, DataPostprocessor
 * has to provide the information, which input data is needed for the
 * calculation of the derived quantities, i.e. whether it needs the
 * values, the first derivative and/or the second derivative of the
 * provided data. DataPostprocessor objects which are used in
 * combination with a DataOutFaces object can also ask for the normal
 * vectors at each point. The information, which data is needed has to
 * be provided via the UpdateFlags returned by the virtual function @p
 * get_needed_update_flags. It is your responsibility to use only
 * those values which were updated in the calculation of derived
 * quantities. The DataOut object will provide references to the
 * requested data in the call to compute_derived_quantities_scalar()
 * or compute_derived_quantities_vector() (DataOut decides which of
 * the two functions to call depending on whether the finite element
 * in use has only a single, or multiple vector components).
 *
 * Furthermore, derived classes have to implement the @p get_names and @p
 * n_output_variables functions, where the number of output variables returned
 * by the latter function has to match the size of the vector returned by the
 * former. Furthermore, this number has to match the number of computed
 * quantities, of course.
 *
 * @ingroup output
 * @author Tobias Leicht, 2007
 */
template <int dim>
class DataPostprocessor: public Subscriptor
{
  public:
				     /**
				      * Virtual desctructor for safety. Does not
				      * do anything so far.
				      */
    virtual ~DataPostprocessor ();

				     /**
				      * @deprecated
				      *
				      * This function only exists for backward
				      * compatibility as this is the interface
				      * provided by previous versions of the
				      * library. The default implementation of
				      * the other function of same name below
				      * calls this function by simply dropping
				      * the argument that denotes the
				      * evaluation points. Since this function
				      * might at one point go away, you should
				      * overload the other function instead.
				      */
    virtual
    void
    compute_derived_quantities_scalar (const std::vector<double>         &uh,
				       const std::vector<Tensor<1,dim> > &duh,
				       const std::vector<Tensor<2,dim> > &dduh,
				       const std::vector<Point<dim> >    &normals,
				       std::vector<Vector<double> >      &computed_quantities) const;
				       
				     /**
				      * This is the main function which actually
				      * performs the postprocessing. The last
				      * argument is a reference to the
				      * postprocessed data which has correct
				      * size already and must be filled by this
				      * function. @p uh is a reference to a
				      * vector of data values at all points, @p
				      * duh the same for gradients, @p dduh for
				      * second derivatives and @p normals is a
				      * reference to the normal vectors. Note,
				      * that the last four references will only
				      * contain valid data, if the respective
				      * flags are returned by @p
				      * get_needed_update_flags, otherwise those
				      * vectors will be in an unspecified
				      * state. @p normals will always be an
				      * empty vector when working on cells, not
				      * on faces.
				      *
				      * This function is called when
				      * the original data vector
				      * represents scalar data,
				      * i.e. the finite element in use
				      * has only a single vector
				      * component.
				      */
    virtual
    void
    compute_derived_quantities_scalar (const std::vector<double>         &uh,
				       const std::vector<Tensor<1,dim> > &duh,
				       const std::vector<Tensor<2,dim> > &dduh,
				       const std::vector<Point<dim> >    &normals,
				       const std::vector<Point<dim> >    &evaluation_points,
				       std::vector<Vector<double> >      &computed_quantities) const;
    
				     /**
				      * @deprecated
				      *
				      * This function only exists for backward
				      * compatibility as this is the interface
				      * provided by previous versions of the
				      * library. The default implementation of
				      * the other function of same name below
				      * calls this function by simply dropping
				      * the argument that denotes the
				      * evaluation points. Since this function
				      * might at one point go away, you should
				      * overload the other function instead.
				      */
    virtual
    void
    compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
				       const std::vector<std::vector<Tensor<1,dim> > > &duh,
				       const std::vector<std::vector<Tensor<2,dim> > > &dduh,
				       const std::vector<Point<dim> >                  &normals,
				       std::vector<Vector<double> >                    &computed_quantities) const;
				       
				     /**
				      * Same as the
				      * compute_derived_quantities_scalar()
				      * function, but this function is called
				      * when the original data vector
				      * represents vector data, i.e. the
				      * finite element in use has multiple
				      * vector components.
				      */
    virtual
    void
    compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
				       const std::vector<std::vector<Tensor<1,dim> > > &duh,
				       const std::vector<std::vector<Tensor<2,dim> > > &dduh,
				       const std::vector<Point<dim> >                  &normals,
				       const std::vector<Point<dim> >                  &evaluation_points,
				       std::vector<Vector<double> >                    &computed_quantities) const;

				     /**
				      * Return the vector of strings describing
				      * the names of the computed quantities.
				      */
    virtual std::vector<std::string> get_names () const=0;
    
				     /**
				      * This functions returns
				      * information about how the
				      * individual components of
				      * output files that consist of
				      * more than one data set are to
				      * be interpreted.
				      *
				      * For example, if one has a finite
				      * element for the Stokes equations in
				      * 2d, representing components (u,v,p),
				      * one would like to indicate that the
				      * first two, u and v, represent a
				      * logical vector so that later on when
				      * we generate graphical output we can
				      * hand them off to a visualization
				      * program that will automatically know
				      * to render them as a vector field,
				      * rather than as two separate and
				      * independent scalar fields.
				      *
				      * The default implementation of this
				      * function returns a vector of values
				      * DataComponentInterpretation::component_is_scalar,
				      * indicating that all output components
				      * are independent scalar
				      * fields. However, if a derived class
				      * produces data that represents vectors,
				      * it may return a vector that contains
				      * values
				      * DataComponentInterpretation::component_is_part_of_vector. In
				      * the example above, one would return a
				      * vector with components
				      * (DataComponentInterpretation::component_is_part_of_vector,
				      * DataComponentInterpretation::component_is_part_of_vector,
				      * DataComponentInterpretation::component_is_scalar)
				      * for (u,v,p).
				      */
    virtual
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    get_data_component_interpretation () const;
    
				     /**
				      * Return, which data has to be provided to
				      * compute the derived quantities. This has
				      * to be a combination of @p update_values,
				      * @p update_gradients and @p
				      * update_hessians. If the
				      * DataPostprocessor is to be used in
				      * combination with DataOutFaces, you may
				      * also ask for a update of normals via the
				      * @p update_normal_vectors flag.
				      */
    virtual UpdateFlags get_needed_update_flags () const=0;

				     /**
				      * Number of postprocessed
				      * variables. Has to match the
				      * number of entries filled by
				      * compute_derived_quantities_scalar()
				      * or
				      * compute_derived_quantities_vector()
				      * as well as the size of the
				      * vector of names returned by
				      * get_names().
				      */
    virtual unsigned int n_output_variables() const=0;
    
};

DEAL_II_NAMESPACE_CLOSE

#endif
