// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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


#ifndef __deal2__mesh_worker_output_h
#define __deal2__mesh_worker_output_h

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/base/named_data.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/base/mg_level_object.h>


DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  namespace Assembler
  {

    /**
     * A class that, instead of assembling into a matrix or vector,
     * outputs the results on a cell to a gnuplot patch.
     *
     * This assembler expects that LocalResults contains quadrature values
     * set with LocalResults::quadrature_value(). When it is initialized
     * with the number of quadrature points in a single (!) space
     * direction and the number of data fields to be displayed, it
     * initializes LocalResults automatically. The number of data fields
     * in local results will be increased by dim in order to accommodate
     * for the coordinates of the data points.
     *
     * While data slots for the space coordinates are allocated
     * automatically, these coordinates are not entered. It is up to the
     * user to enter the coordinates in the first dim data entries at
     * every point. This adds the flexibility to output transformed
     * coordinates or even something completely different.
     *
     * @note In the current implementation, only cell data can be written.
     *
     * @author Guido Kanschat
     * @date 2011, 2012
     */
    class GnuplotPatch
    {
    public:
      /**
       * Constructor.
       */
      GnuplotPatch();

      /**
       * Initialize for writing
       * <i>n</i> data vectors. The
       * number of points is the
       * number of quadrature
       * points in a single
       * direction in a tensor
       * product formula. It must
       * match the number in the
       * actual Quadrature used to
       * create the patches. The
       * total number of data
       * vectors produced is <tt>n+dim</tt>
       * and the first dim should be
       * the space coordinates of
       * the points. Nevertheless,
       * it is up to the user to
       * set these values to
       * whatever is desired.
       */
      void initialize (const unsigned int n_points,
                       const unsigned int n_vectors);

      /**
       * Set the stream #os to
       * which data is written. If
       * no stream is selected with
       * this function, data goes
       * to @p deallog.
       */
      void initialize_stream (std::ostream &stream);

      /**
       * Initialize the local data
       * in the
       * DoFInfo
       * object used later for
       * assembling.
       *
       * The info object refers to
       * a cell if
       * <code>!face</code>, or
       * else to an interior or
       * boundary face.
       */
      template <int dim>
      void initialize_info(DoFInfo<dim> &info, bool face);

      /**
       * Write the patch to the
       * output stream.
       */
      template<int dim>
      void assemble(const DoFInfo<dim> &info);

      /**
       * @warning Not implemented yet
       */
      template<int dim>
      void assemble(const DoFInfo<dim> &info1,
                    const DoFInfo<dim> &info2);

    private:
      /**
       * Write the object T either
       * to the stream #os, if initialize_stream()
       * has been called, or to
       * @p deallog if no pointer has
       * been set.
       */
      template<typename T>
      void write(const T &t) const;

      /**
       * Write an end-of-line marker either
       * to the stream #os, if initialize_stream
       * has been called, or to
       * @p deallog if no pointer has
       * been set.
       */
      void write_endl () const;

      /**
       * The number of output
       * components in each point.
       */
      unsigned int n_vectors;
      /**
       * The number of points in
       * one direction.
       */
      unsigned int n_points;

      /**
       * Stream to which output is to be written. Set by initialize_stream().
       */
      std::ostream *os;
    };

//----------------------------------------------------------------------//

    template <typename T>
    inline void
    GnuplotPatch::write(const T &d) const
    {
      if (os == 0)
        deallog << d;
      else
        (*os) << d;
    }


    inline void
    GnuplotPatch::write_endl() const
    {
      if (os == 0)
        deallog << std::endl;
      else
        (*os) << std::endl;
    }


    inline
    GnuplotPatch::GnuplotPatch()
      :
      os(0)
    {}


    inline void
    GnuplotPatch::initialize (const unsigned int np,
                              const unsigned int nv)
    {
      n_vectors = nv;
      n_points = np;
    }


    inline void
    GnuplotPatch::initialize_stream (std::ostream &stream)
    {
      os = &stream;
    }


    template <int dim>
    inline void
    GnuplotPatch::initialize_info(DoFInfo<dim> &info, bool face)
    {
      if (face)
        info.initialize_quadrature(Utilities::fixed_power<dim-1>(n_points), n_vectors+dim);
      else
        info.initialize_quadrature(Utilities::fixed_power<dim>(n_points), n_vectors+dim);
    }


    template <int dim>
    inline void
    GnuplotPatch::assemble(const DoFInfo<dim> &info)
    {
      const unsigned int np = info.n_quadrature_points();
      const unsigned int nv = info.n_quadrature_values();
      const unsigned int patch_dim = (info.face_number == numbers::invalid_unsigned_int)
                                     ? dim : (dim-1);
      const unsigned int row_length = n_points;
      // If patches are 1D, end the
      // patch after a row, else end
      // it after a square
      const unsigned int row_length2 = (patch_dim==1) ? row_length : (row_length*row_length);

//      AssertDimension(np, Utilities::fixed_power<dim>(n_points));
      AssertDimension(nv, n_vectors+dim);


      for (unsigned int k=0; k<np; ++k)
        {
          if (k % row_length == 0)
            write_endl();
          if (k % row_length2 == 0)
            write_endl();

          for (unsigned int i=0; i<nv; ++i)
            {
              write(info.quadrature_value(k,i));
              write('\t');
            }
          write_endl();
        }
    }


    template <int dim>
    inline void
    GnuplotPatch::assemble(const DoFInfo<dim> &info1, const DoFInfo<dim> &info2)
    {
      assemble(info1);
      assemble(info2);
    }
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
