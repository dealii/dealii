// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/base/utilities.h>
#include <deal.II/base/logstream.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/vector_tools.h>


DEAL_II_NAMESPACE_OPEN

namespace Functions
{

  template <int dim, typename DoFHandlerType, typename VectorType>
  FEFieldFunction<dim, DoFHandlerType, VectorType>::FEFieldFunction
  (const DoFHandlerType &mydh,
   const VectorType     &myv,
   const Mapping<dim>   &mymapping)
    :
    Function<dim,typename VectorType::value_type>(mydh.get_fe().n_components()),
    dh(&mydh, "FEFieldFunction"),
    data_vector(myv),
    mapping(mymapping),
    cell_hint(dh->end()),
    n_components(mydh.get_fe().n_components())
  {
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::
  set_active_cell(const typename DoFHandlerType::active_cell_iterator &newcell)
  {
    cell_hint.get() = newcell;
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  void FEFieldFunction<dim, DoFHandlerType, VectorType>::vector_value (const Point<dim> &p,
      Vector<typename VectorType::value_type>   &values) const
  {
    Assert (values.size() == n_components,
            ExcDimensionMismatch(values.size(), n_components));
    typename DoFHandlerType::active_cell_iterator cell = cell_hint.get();
    if (cell == dh->end())
      cell = dh->begin_active();

    boost::optional<Point<dim> >
    qp = get_reference_coordinates (cell, p);
    if (!qp)
      {
        const std::pair<typename dealii::internal::ActiveCellIterator<dim, dim, DoFHandlerType>::type, Point<dim> > my_pair
          = GridTools::find_active_cell_around_point (mapping, *dh, p);
        AssertThrow (my_pair.first->is_locally_owned(),
                     VectorTools::ExcPointNotAvailableHere());

        cell = my_pair.first;
        qp = my_pair.second;
      }

    cell_hint.get() = cell;

    // Now we can find out about the point
    Quadrature<dim> quad(qp.get());
    FEValues<dim> fe_v(mapping, cell->get_fe(), quad,
                       update_values);
    fe_v.reinit(cell);
    std::vector< Vector<typename VectorType::value_type> >
    vvalues (1, Vector<typename VectorType::value_type>(values.size()));
    fe_v.get_function_values(data_vector, vvalues);
    values = vvalues[0];
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  typename VectorType::value_type
  FEFieldFunction<dim, DoFHandlerType, VectorType>::value (const Point<dim>   &p,
                                                           const unsigned int comp) const
  {
    Vector<typename VectorType::value_type> values(n_components);
    vector_value(p, values);
    return values(comp);
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::
  vector_gradient (const Point<dim>            &p,
                   std::vector<Tensor<1,dim,typename VectorType::value_type> > &gradients) const
  {
    typedef typename VectorType::value_type number;
    Assert (gradients.size() == n_components,
            ExcDimensionMismatch(gradients.size(), n_components));
    typename DoFHandlerType::active_cell_iterator cell = cell_hint.get();
    if (cell == dh->end())
      cell = dh->begin_active();

    boost::optional<Point<dim> >
    qp = get_reference_coordinates (cell, p);
    if (!qp)
      {
        const std::pair<typename dealii::internal::ActiveCellIterator<dim, dim, DoFHandlerType>::type, Point<dim> > my_pair
          = GridTools::find_active_cell_around_point (mapping, *dh, p);
        AssertThrow (my_pair.first->is_locally_owned(),
                     VectorTools::ExcPointNotAvailableHere());

        cell = my_pair.first;
        qp = my_pair.second;
      }

    cell_hint.get() = cell;

    // Now we can find out about the point
    Quadrature<dim> quad(qp.get());
    FEValues<dim> fe_v(mapping, cell->get_fe(), quad,
                       update_gradients);
    fe_v.reinit(cell);

    if (n_components == 1)
      {
        // the size of the @p gradients coincidentally coincides
        // with the number of quadrature points we evaluate the function at.
        fe_v.get_function_gradients(data_vector, gradients);
      }
    else
      {
        // Unfortunately we still need a temporary argument as we want to
        // evaluate a gradient of a (generally) multicomponent function at
        // a single quadrature point. Note that the first std::vector<> is related
        // to the number of quadrature points (always one here), whereas the second
        // to the number of components.
        std::vector< std::vector<Tensor<1,dim,number> > > vgrads
        (1,  std::vector<Tensor<1,dim,number> >(n_components) );
        fe_v.get_function_gradients(data_vector, vgrads);
        gradients = vgrads[0];
      }
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  Tensor<1,dim,typename VectorType::value_type>
  FEFieldFunction<dim, DoFHandlerType, VectorType>::
  gradient (const Point<dim>   &p,
            const unsigned int comp) const
  {
    std::vector<Tensor<1,dim,typename VectorType::value_type> > grads(n_components);
    vector_gradient(p, grads);
    return grads[comp];
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::
  vector_laplacian (const Point<dim> &p,
                    Vector<typename VectorType::value_type>   &values) const
  {
    Assert (values.size() == n_components,
            ExcDimensionMismatch(values.size(), n_components));
    typename DoFHandlerType::active_cell_iterator cell = cell_hint.get();
    if (cell == dh->end())
      cell = dh->begin_active();

    boost::optional<Point<dim> >
    qp = get_reference_coordinates (cell, p);
    if (!qp)
      {
        const std::pair<typename dealii::internal::ActiveCellIterator<dim, dim, DoFHandlerType>::type, Point<dim> > my_pair
          = GridTools::find_active_cell_around_point (mapping, *dh, p);
        AssertThrow (my_pair.first->is_locally_owned(),
                     VectorTools::ExcPointNotAvailableHere());

        cell = my_pair.first;
        qp = my_pair.second;
      }

    cell_hint.get() = cell;

    // Now we can find out about the point
    Quadrature<dim> quad(qp.get());
    FEValues<dim> fe_v(mapping, cell->get_fe(), quad,
                       update_hessians);
    fe_v.reinit(cell);
    std::vector< Vector<typename VectorType::value_type> >
    vvalues (1, Vector<typename VectorType::value_type>(values.size()));
    fe_v.get_function_laplacians(data_vector, vvalues);
    values = vvalues[0];
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  typename VectorType::value_type FEFieldFunction<dim, DoFHandlerType, VectorType>::laplacian
  (const Point<dim>   &p,
   const unsigned int  comp) const
  {
    Vector<typename VectorType::value_type> lap(n_components);
    vector_laplacian(p, lap);
    return lap[comp];
  }


  // Now the list versions
  // ==============================

  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::
  vector_value_list (const std::vector<Point< dim > > &points,
                     std::vector< Vector<typename VectorType::value_type> >    &values) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    std::vector<typename DoFHandlerType::active_cell_iterator > cells;
    std::vector<std::vector<Point<dim> > > qpoints;
    std::vector<std::vector<unsigned int> > maps;

    unsigned int ncells = compute_point_locations(points, cells, qpoints, maps);
    hp::MappingCollection<dim> mapping_collection (mapping);
    hp::FECollection<dim> fe_collection (dh->get_fe ());
    hp::QCollection<dim> quadrature_collection;
    // Create quadrature collection
    for (unsigned int i=0; i<ncells; ++i)
      {
        // Number of quadrature points on this cell
        unsigned int nq = qpoints[i].size();
        // Construct a quadrature formula
        std::vector< double > ww(nq, 1./((double) nq));

        quadrature_collection.push_back (Quadrature<dim> (qpoints[i], ww));
      }
    // Get a function value object
    hp::FEValues<dim> fe_v(mapping_collection, fe_collection, quadrature_collection,
                           update_values);
    // Now gather all the informations we need
    for (unsigned int i=0; i<ncells; ++i)
      {
        fe_v.reinit(cells[i], i, 0);
        const unsigned int nq = qpoints[i].size();
        std::vector< Vector<typename VectorType::value_type> > vvalues (nq, Vector<typename VectorType::value_type>(n_components));
        fe_v.get_present_fe_values ().get_function_values(data_vector, vvalues);
        for (unsigned int q=0; q<nq; ++q)
          values[maps[i][q]] = vvalues[q];
      }
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::
  value_list (const std::vector<Point< dim > > &points,
              std::vector< typename VectorType::value_type >            &values,
              const unsigned int                component) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));
    std::vector< Vector<typename VectorType::value_type> > vvalues(points.size(), Vector<typename VectorType::value_type>(n_components));
    vector_value_list(points, vvalues);
    for (unsigned int q=0; q<points.size(); ++q)
      values[q] = vvalues[q](component);
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::
  vector_gradient_list (const std::vector<Point< dim > >           &points,
                        std::vector<std::vector< Tensor<1,dim, typename VectorType::value_type> > > &values) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    std::vector<typename DoFHandlerType::active_cell_iterator > cells;
    std::vector<std::vector<Point<dim> > > qpoints;
    std::vector<std::vector<unsigned int> > maps;

    unsigned int ncells = compute_point_locations(points, cells, qpoints, maps);
    hp::MappingCollection<dim> mapping_collection (mapping);
    hp::FECollection<dim> fe_collection (dh->get_fe ());
    hp::QCollection<dim> quadrature_collection;
    // Create quadrature collection
    for (unsigned int i=0; i<ncells; ++i)
      {
        // Number of quadrature points on this cell
        unsigned int nq = qpoints[i].size();
        // Construct a quadrature formula
        std::vector< double > ww(nq, 1./((double) nq));

        quadrature_collection.push_back (Quadrature<dim> (qpoints[i], ww));
      }
    // Get a function value object
    hp::FEValues<dim> fe_v(mapping_collection, fe_collection, quadrature_collection,
                           update_gradients);
    // Now gather all the informations we need
    for (unsigned int i=0; i<ncells; ++i)
      {
        fe_v.reinit(cells[i], i, 0);
        const unsigned int nq = qpoints[i].size();
        std::vector< std::vector<Tensor<1,dim,typename VectorType::value_type> > >
        vgrads (nq, std::vector<Tensor<1,dim,typename VectorType::value_type> >(n_components));
        fe_v.get_present_fe_values ().get_function_gradients(data_vector, vgrads);
        for (unsigned int q=0; q<nq; ++q)
          {
            const unsigned int s = vgrads[q].size();
            values[maps[i][q]].resize(s);
            for (unsigned int l=0; l<s; l++)
              values[maps[i][q]][l] = vgrads[q][l];
          }
      }
  }

  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::
  gradient_list (const std::vector<Point< dim > > &points,
                 std::vector< Tensor<1,dim, typename VectorType::value_type> >     &values,
                 const unsigned int                component) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));
    std::vector< std::vector<Tensor<1,dim, typename VectorType::value_type> > >
    vvalues(points.size(), std::vector<Tensor<1,dim,typename VectorType::value_type> >(n_components));
    vector_gradient_list(points, vvalues);
    for (unsigned int q=0; q<points.size(); ++q)
      values[q] = vvalues[q][component];
  }


  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::
  vector_laplacian_list (const std::vector<Point< dim > > &points,
                         std::vector< Vector<typename VectorType::value_type> >    &values) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));

    std::vector<typename DoFHandlerType::active_cell_iterator> cells;
    std::vector<std::vector<Point<dim> > > qpoints;
    std::vector<std::vector<unsigned int> > maps;

    unsigned int ncells = compute_point_locations(points, cells, qpoints, maps);
    hp::MappingCollection<dim> mapping_collection (mapping);
    hp::FECollection<dim> fe_collection (dh->get_fe ());
    hp::QCollection<dim> quadrature_collection;
    // Create quadrature collection
    for (unsigned int i=0; i<ncells; ++i)
      {
        // Number of quadrature points on this cell
        unsigned int nq = qpoints[i].size();
        // Construct a quadrature formula
        std::vector< double > ww(nq, 1./((double) nq));

        quadrature_collection.push_back (Quadrature<dim> (qpoints[i], ww));
      }
    // Get a function value object
    hp::FEValues<dim> fe_v(mapping_collection, fe_collection, quadrature_collection,
                           update_hessians);
    // Now gather all the informations we need
    for (unsigned int i=0; i<ncells; ++i)
      {
        fe_v.reinit(cells[i], i, 0);
        const unsigned int nq = qpoints[i].size();
        std::vector< Vector<typename VectorType::value_type> > vvalues (nq, Vector<typename VectorType::value_type>(n_components));
        fe_v.get_present_fe_values ().get_function_laplacians(data_vector, vvalues);
        for (unsigned int q=0; q<nq; ++q)
          values[maps[i][q]] = vvalues[q];
      }
  }

  template <int dim, typename DoFHandlerType, typename VectorType>
  void
  FEFieldFunction<dim, DoFHandlerType, VectorType>::
  laplacian_list (const std::vector<Point<dim> > &points,
                  std::vector<typename VectorType::value_type>            &values,
                  const unsigned int              component) const
  {
    Assert(points.size() == values.size(),
           ExcDimensionMismatch(points.size(), values.size()));
    std::vector< Vector<typename VectorType::value_type> > vvalues(points.size(), Vector<typename VectorType::value_type>(n_components));
    vector_laplacian_list(points, vvalues);
    for (unsigned int q=0; q<points.size(); ++q)
      values[q] = vvalues[q](component);
  }



  template <int dim, typename DoFHandlerType, typename VectorType>
  unsigned int FEFieldFunction<dim, DoFHandlerType, VectorType>::
  compute_point_locations
  (const std::vector<Point<dim> >                              &points,
   std::vector<typename DoFHandlerType::active_cell_iterator > &cells,
   std::vector<std::vector<Point<dim> > >                      &qpoints,
   std::vector<std::vector<unsigned int> >                     &maps) const
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

    // Set this to true until all
    // points have been classified
    bool left_over = true;

    // Current quadrature point
    typename DoFHandlerType::active_cell_iterator cell = cell_hint.get();
    if (cell == dh->end())
      cell = dh->begin_active();

    {
      // see if the point is
      // inside the
      // cell. there are two
      // ways that
      // transform_real_to_unit_cell
      // can indicate that a
      // point is outside: by
      // returning
      // coordinates that lie
      // outside the
      // reference cell, or
      // by throwing an
      // exception. handle
      // both
      boost::optional<Point<dim> >
      qp = get_reference_coordinates (cell, points[0]);
      if (!qp)
        {
          const std::pair<typename dealii::internal::ActiveCellIterator<dim, dim, DoFHandlerType>::type, Point<dim> >
          my_pair  = GridTools::find_active_cell_around_point
                     (mapping, *dh, points[0]);
          AssertThrow (my_pair.first->is_locally_owned(),
                       VectorTools::ExcPointNotAvailableHere());

          cell = my_pair.first;
          qp = my_pair.second;
          point_flags[0] = true;
        }

      // Put in the first point.
      cells.push_back(cell);
      qpoints.push_back(std::vector<Point<dim> >(1, qp.get()));
      maps.push_back(std::vector<unsigned int> (1, 0));
    }


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
          if (point_flags[p] == false)
            {
              // same logic as above
              const boost::optional<Point<dim> >
              qp = get_reference_coordinates (cells[c], points[p]);
              if (qp)
                {
                  point_flags[p] = true;
                  qpoints[c].push_back(qp.get());
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
            const std::pair<typename dealii::internal::ActiveCellIterator<dim, dim, DoFHandlerType>::type, Point<dim> > my_pair
              = GridTools::find_active_cell_around_point (mapping, *dh, points[first_outside]);
            AssertThrow (my_pair.first->is_locally_owned(),
                         VectorTools::ExcPointNotAvailableHere());

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


  template <int dim, typename DoFHandlerType, typename VectorType>
  boost::optional<Point<dim> >
  FEFieldFunction<dim, DoFHandlerType, VectorType>::
  get_reference_coordinates (const typename DoFHandlerType::active_cell_iterator &cell,
                             const Point<dim>                                    &point) const
  {
    try
      {
        Point<dim> qp =  mapping.transform_real_to_unit_cell(cell, point);
        if (GeometryInfo<dim>::is_inside_unit_cell(qp))
          return qp;
        else
          return boost::optional<Point<dim> >();
      }
    catch (const typename Mapping<dim>::ExcTransformationFailed &)
      {
        // transformation failed, so
        // assume the point is
        // outside
        return boost::optional<Point<dim> >();
      }
  }

}

DEAL_II_NAMESPACE_CLOSE
