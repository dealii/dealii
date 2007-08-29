//---------------------------------------------------------------------------
//    $Id: function_parser.h 14594 2007-03-22 20:17:41Z bangerth $
//    Version: $Name$
//
//    Copyright (C) 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/utilities.h>
#include <base/logstream.h>
#include <grid/grid_tools.h>
#include <fe/fe_values.h>
#include <numerics/fe_field_function.h>


DEAL_II_NAMESPACE_OPEN

namespace Functions
{

  template <int dim, typename DH, typename VECTOR>
  FEFieldFunction<dim, DH, VECTOR>::FEFieldFunction (const DH &mydh, 
						     const VECTOR &myv,
						     const Mapping<dim> &mymapping)
		  : 
		  Function<dim>(mydh.get_fe().n_components()),
		  dh(&mydh, "FEFieldFunction"),
		  data_vector(myv),
		  mapping(mymapping),
		  n_components(mydh.get_fe().n_components())
  {
    cell = dh->begin_active();
  }


  
  template <int dim, typename DH, typename VECTOR>
  void
  FEFieldFunction<dim, DH, VECTOR>::
  set_active_cell(typename DH::active_cell_iterator &newcell)
  {
    cell = newcell;
  }


  
  template <int dim, typename DH, typename VECTOR>
  void FEFieldFunction<dim, DH, VECTOR>::vector_value (const Point<dim> &p,
						       Vector<double>   &values) const 
  { 
    Assert (values.size() == n_components, 
	    ExcDimensionMismatch(values.size(), n_components));
    Point<dim> qp = mapping.transform_real_to_unit_cell(cell, p);
  
				     // Check if we already have all we need
    if (!GeometryInfo<dim>::is_inside_unit_cell(qp))
      {
	const std::pair<typename DH::active_cell_iterator, Point<dim> > my_pair 
	  = GridTools::find_active_cell_around_point (mapping, *dh, p); 
	cell = my_pair.first;
	qp = my_pair.second;
      }				    
  
				     // Now we can find out about the point
    Quadrature<dim> quad(qp);
    FEValues<dim> fe_v(mapping, dh->get_fe(), quad, 
		       update_values);
    fe_v.reinit(cell);
    std::vector< Vector<double> > vvalues (1, values);
    fe_v.get_function_values(data_vector, vvalues);
    values = vvalues[0];
  }


  
  template <int dim, typename DH, typename VECTOR>
  double
  FEFieldFunction<dim, DH, VECTOR>::value (const Point<dim>   &p,
					   const unsigned int comp) const
  { 
    Vector<double> values(n_components);
    vector_value(p, values);
    return values(comp);
  }


  template <int dim, typename DH, typename VECTOR>
  void
  FEFieldFunction<dim, DH, VECTOR>::vector_gradient 
  (const Point<dim> &p,
   std::vector<Tensor<1,dim> > &gradients) const 
  { 
    Assert (gradients.size() == n_components, 
	    ExcDimensionMismatch(gradients.size(), n_components));
    Point<dim> qp = mapping.transform_real_to_unit_cell(cell, p);
  
				     // Check if we already have all we need
    if (!GeometryInfo<dim>::is_inside_unit_cell(qp))
      {
	std::pair<typename DH::active_cell_iterator, Point<dim> > my_pair  
	  = GridTools::find_active_cell_around_point (mapping, *dh, p); 
	cell = my_pair.first;
	qp = my_pair.second;
      }				    
  
				     // Now we can find out about the point
    Quadrature<dim> quad(qp);
    FEValues<dim> fe_v(mapping, dh->get_fe(), quad, 
		       update_gradients);
    fe_v.reinit(cell);
    std::vector< std::vector<Tensor<1,dim> > > vgrads 
      (1,  std::vector<Tensor<1,dim> >(n_components) );
    fe_v.get_function_grads(data_vector, vgrads);
    gradients = vgrads[0];
  }


  
  template <int dim, typename DH, typename VECTOR>
  Tensor<1,dim> FEFieldFunction<dim, DH, VECTOR>::gradient 
  (const Point<dim>   &p, unsigned int comp) const
  { 
    std::vector<Tensor<1,dim> > grads(n_components);
    vector_gradient(p, grads);
    return grads[comp];
  }

				   // Now the list versions
				   // ==============================

  template <int dim, typename DH, typename VECTOR>
  void
  FEFieldFunction<dim, DH, VECTOR>::
  vector_value_list (const std::vector<Point< dim > > &    points,
		     std::vector< Vector<double> > &values) const
  { 
    Assert(points.size() == values.size(),
	   ExcDimensionMismatch(points.size(), values.size()));

    std::vector<typename DH::active_cell_iterator > cells;
    std::vector<std::vector<Point<dim> > > qpoints;
    std::vector<std::vector<unsigned int> > maps;
  
    unsigned int ncells = compute_point_locations(points, cells, qpoints, maps);
  
				     // Now gather all the informations we need
    for (unsigned int i=0; i<ncells; ++i)
      {
					 // Number of quadrature points on this cell
	unsigned int nq = qpoints[i].size();
    
					 // Construct a quadrature formula
	std::vector< double > ww(nq, 1./((double) nq));
	Quadrature<dim> quad(qpoints[i], ww);
    
					 // Get a function value object
	FEValues<dim> fe_v(mapping, dh->get_fe(), quad, 
			   update_values);
	fe_v.reinit(cells[i]);
	std::vector< Vector<double> > vvalues (nq, Vector<double>(n_components));
	fe_v.get_function_values(data_vector, vvalues);
	for (unsigned int q=0; q<nq; ++q)
	  values[maps[i][q]] = vvalues[q];
      }
  }

  template <int dim, typename DH, typename VECTOR>
  void
  FEFieldFunction<dim, DH, VECTOR>::
  value_list (const std::vector<Point< dim > > &points,
	      std::vector< double > &values, 
	      const unsigned int  component) const
  { 
    Assert(points.size() == values.size(),
	   ExcDimensionMismatch(points.size(), values.size()));
    std::vector< Vector<double> > vvalues(points.size(), Vector<double>(n_components));
    vector_value_list(points, vvalues);
    for (unsigned int q=0; q<points.size(); ++q)
      values[q] = vvalues[q](component);
  }


  
  template <int dim, typename DH, typename VECTOR>
  void
  FEFieldFunction<dim, DH, VECTOR>::
  vector_gradient_list (const std::vector<Point< dim > > &    points,
			std::vector< 
			std::vector< Tensor<1,dim> > > &values) const
  { 
    Assert(points.size() == values.size(),
	   ExcDimensionMismatch(points.size(), values.size()));

    std::vector<typename DH::active_cell_iterator > cells;
    std::vector<std::vector<Point<dim> > > qpoints;
    std::vector<std::vector<unsigned int> > maps;
  
    unsigned int ncells = compute_point_locations(points, cells, qpoints, maps);
  
				     // Now gather all the informations we need
    for (unsigned int i=0; i<ncells; ++i)
      {
					 // Number of quadrature points on this cell
	unsigned int nq = qpoints[i].size();
    
					 // Construct a quadrature formula
	std::vector< double > ww(nq, 1./((double) nq));
	Quadrature<dim> quad(qpoints[i], ww);
    
					 // Get a function value object
	FEValues<dim> fe_v(mapping, dh->get_fe(), quad, 
			   update_gradients);
	fe_v.reinit(cells[i]);
	std::vector< std::vector<Tensor<1,dim> > >
	  vgrads (nq, std::vector<Tensor<1,dim> >(n_components));
	fe_v.get_function_grads(data_vector, vgrads);
	for (unsigned int q=0; q<nq; ++q)
	  values[maps[i][q]] = vgrads[q];
      }
  }

  template <int dim, typename DH, typename VECTOR>
  void
  FEFieldFunction<dim, DH, VECTOR>::
  gradient_list (const std::vector<Point< dim > > &points,
		 std::vector< Tensor<1,dim> > &values, 
		 const unsigned int  component) const
  { 
    Assert(points.size() == values.size(),
	   ExcDimensionMismatch(points.size(), values.size()));
    std::vector< std::vector<Tensor<1,dim> > >
      vvalues(points.size(), std::vector<Tensor<1,dim> >(n_components));
    vector_gradient_list(points, vvalues);
    for (unsigned int q=0; q<points.size(); ++q)
      values[q] = vvalues[q][component];
  }

  

  template <int dim, typename DH, typename VECTOR>
  unsigned int FEFieldFunction<dim, DH, VECTOR>::
  compute_point_locations(const std::vector<Point<dim> > &points,
			  std::vector<typename DH::active_cell_iterator > &cells,
			  std::vector<std::vector<Point<dim> > > &qpoints,
			  std::vector<std::vector<unsigned int> > &maps) const
  {
				     // How many points are here?
    const unsigned int np = points.size();
  
				     // Reset output maps.
    cells.clear();
    qpoints.clear();
    maps.clear();
    
				     // Now the easy case.
    if (np==0) return 0;
    
				     // Keep track of the points we
				     // found
    std::vector<bool> point_flags(np, false);
  
				     // Set this to true untill all
				     // points have been classified
    bool left_over = true;
  
				     // Current quadrature point
    Point<dim> qp = mapping.transform_real_to_unit_cell(cell, points[0]);
  
				     // Check if we already have a
				     // valid cell for the first point
    if (!GeometryInfo<dim>::is_inside_unit_cell(qp))
      {
	const std::pair<typename DH::active_cell_iterator, Point<dim> >
	  my_pair  = GridTools::find_active_cell_around_point
	  (mapping, *dh, points[0]); 
	cell = my_pair.first;
	qp = my_pair.second;
	point_flags[0] = true;
      }
    
				     // Put in the first point.
    cells.push_back(cell);
    qpoints.push_back(std::vector<Point<dim> >(1, qp));
    maps.push_back(std::vector<unsigned int> (1, 0));
  
				     // Check if we need to do anything else
    if (points.size() > 1)
      left_over = true;
    else
      left_over = false;

  
				     // This is the first index of a non processed point
    unsigned int first_outside = 1;
  
				     // And this is the index of the current cell
    unsigned int c = 0;
  
    while (left_over == true)
      {
					 // Assume this is the last one
	left_over = false;
	Assert(first_outside < np,
	       ExcIndexRange(first_outside, 0, np));
    
					 // If we found one in this cell, keep looking in the same cell
	for (unsigned int p=first_outside; p<np; ++p) 
	  if (point_flags[p] == false) {
	    Point<dim> qpoint =  mapping.transform_real_to_unit_cell(cell, points[p]);
	    if (GeometryInfo<dim>::is_inside_unit_cell(qpoint))
	      {
		point_flags[p] = true;
		qpoints[c].push_back(qpoint);
		maps[c].push_back(p);
	      }
	    else
	      {
						 // Set things up for next round 
		if (left_over == false)
		  first_outside = p;
		left_over = true;
	      }
	  }
					 // If we got here and there is
					 // no left over, we are
					 // done. Else we need to find
					 // the next cell
	if (left_over == true)
	  {
	    const std::pair<typename DH::active_cell_iterator, Point<dim> > my_pair  
	      = GridTools::find_active_cell_around_point (mapping, *dh, points[first_outside]); 
	    cells.push_back(my_pair.first);
	    qpoints.push_back(std::vector<Point<dim> >(1, my_pair.second));
	    maps.push_back(std::vector<unsigned int>(1, first_outside));
	    c++;
	    point_flags[first_outside] = true;
					     // And check if we can exit the loop now
	    if (first_outside == np-1)
	      left_over = false;
	  }			
      }
  
				     // Augment of one the number of cells
    ++c;
				     // Debug Checking
    Assert(c == cells.size(), ExcInternalError());
  
    Assert(c == maps.size(),
	   ExcDimensionMismatch(c, maps.size()));
  
    Assert(c == qpoints.size(),
	   ExcDimensionMismatch(c, qpoints.size()));
  
#ifdef DEBUG
    unsigned int qps = 0;
				     // The number of points in all
				     // the cells must be the same as
				     // the number of points we
				     // started off from.
    for (unsigned int n=0; n<c; ++n)
      {
	Assert(qpoints[n].size() == maps[n].size(),
	       ExcDimensionMismatch(qpoints[n].size(), maps[n].size()));
	qps += qpoints[n].size();
      }
    Assert(qps == np,
	   ExcDimensionMismatch(qps, np));
#endif
    
    return c;
  }
}

DEAL_II_NAMESPACE_CLOSE
