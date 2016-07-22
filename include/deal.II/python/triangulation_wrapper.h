// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

#ifndef dealii__triangulation_wrapper_h
#define dealii__triangulation_wrapper_h

#include <boost/python.hpp>
#include <deal.II/python/point_wrapper.h>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  class CellAccessorWrapper;

  class TriangulationWrapper
  {
  public:
    typedef std::vector<CellAccessorWrapper>::iterator iterator;

    /**
     * Constructor. Takes a string @p dim with one of the following values
     * "2D", "2d", "3D", or "3d".
     */
    TriangulationWrapper(const std::string &dim);

    /**
     * Constructor. Takes a string @p dim with one of the following values
     * "2D", "2d", "3D", or "3d" and a string @p spacedim with one of the
     * following values "2D", "2d", "3D", or "3d". The dimension of @p spacedim
     * must be larger than the dimension of @p dim
     */
    TriangulationWrapper(const std::string &dim, const std::string &spacedim);

    /**
     * Destructor.
     */
    ~TriangulationWrapper();

    /**
     * Return the number of active cells.
     */
    unsigned int n_active_cells() const;

    /**
     * Generate a hyper cube (square in 2D and cube in 3D) with exactly one
     * cell.
     */
    void generate_hyper_cube(const double left = 0.,
                             const double right = 1.,
                             const bool   colorize = false);

    /**
     * Generate a simplex with (dim+1) vertices and mesh cells.
     */
    void generate_simplex(boost::python::list &vertices);

    /**
     * Same as hyper_cube but not only one cell is created but each coordinate
     * direction is subdivided in @p repetitions cells.
     */
    void generate_subdivided_hyper_cube(const unsigned int repetitions,
                                        const double       left = 0.,
                                        const double       right = 1.);

    /**
     * Generate a coordinate-parallel brick from the two diagonally opposite
     * corners points @p1 and @p2.
     */
    void generate_hyper_rectangle(PointWrapper &p1,
                                  PointWrapper &p2,
                                  const bool    colorize = false);

    /**
     * Generate a coordinate-parallel brick from the two diagonally opposite
     * corners points @p1 and @p2. In direction i, repetitions[i] cells are
     * created.
     */
    void generate_subdivided_hyper_rectangle(boost::python::list &repetitions,
                                             PointWrapper        &p1,
                                             PointWrapper        &p2,
                                             const bool           colorize = false);

    /**
     * Generate a hyperball, i.e. a circle or a ball around @p center with
     * given @p radius.
     */
    void generate_hyper_ball(PointWrapper &center,
                             const double  radius = 1.);

    /**
     * Shift each vertex of the Triangulation by the given @p shift_list.
     */
    void shift(boost::python::list &shift_list);

    /**
     * Given two triangulations, create the triangulation that contains the
     * cells of both triangulations.
     */
    void merge_triangulations(TriangulationWrapper &triangulation_1,
                              TriangulationWrapper &triangulation_2);

    /**
     * Refine all the cells @p n times.
     */
    void refine_global(const unsigned int n);

    /**
     * Execute both refinement and coarsening of the Triangulation.
     */
    void execute_coarsening_and_refinement();

    /**
     * Return the list of active cell accessors associated to the underlying
     * Triangulation.
     */
    boost::python::list active_cells();

    /**
     * Write mesh to the output file @filename according to the given data format.
     */
    void write(const std::string &filename, const std::string format) const;

    /**
     * Write the Triangulation in file.
     */
    void save(const std::string &filename) const;

    /**
     * Load the Triangulation from a file.
     */
    void load(const std::string &filename);

    /**
     * Return the dimension of the underlying Triangulation object.
     */
    int get_dim() const;

    /**
     * Return the space dimension of the underlying Triangulation object.
     */
    int get_spacedim() const;

    /**
     * Return a pointer that can be casted to the underlying Triangulation
     * object.
     */
    void *get_triangulation();

  private:
    /**
     * Helper function for the contructors.
     */
    void setup(const std::string &dimension, const std::string &spacedimension);

    /**
     * Flag is true if the vector cells has been updated.
     */
    bool updated;

    /**
     * Dimension of the underlying Triangulation object.
     */
    int dim;

    /**
     * Space dimension of the underlying Triangulation object.
     */
    int spacedim;

    /**
     * Pointer that can be casted to the underlying Triangulation object.
     */
    void *triangulation;
  };


//-------------------- inline functions -----------------------//



  inline
  int TriangulationWrapper::get_dim() const
  {
    return dim;
  }



  inline
  int TriangulationWrapper::get_spacedim() const
  {
    return spacedim;
  }



  inline
  void *TriangulationWrapper::get_triangulation()
  {
    return triangulation;
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
