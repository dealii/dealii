//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__mesh_worker_output_h
#define __deal2__mesh_worker_output_h

#include <numerics/mesh_worker_info.h>
#include <base/named_data.h>
#include <base/smartpointer.h>
#include <lac/block_vector.h>
#include <multigrid/mg_level_object.h>


DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  namespace Assembler
  {
    
/**
 * A class that, instead of assembling into a matrix or vector,
 * outputs the results on a cell to a gnuplot patch. This works only
 * for elements with support points. The first dim data vectors will
 * be the coordinates, the following are the data.
 *
 * @note In the current implementation, only cell data can be written.
 */
    class GnuplotPatch
    {
					 /**
					  * Constructor.
					  */
	GnuplotPatch();

					 /**
					  * Initialize for writing
					  * <i>n</i> data vectors. The
					  * total number of data
					  * vectors produced is n+dim
					  * and the first dim will be
					  * the space coordinates of
					  * the points.
					  */
	void initialize(unsigned int dim, unsigned int n_vectors, unsigned int n_points);
	
					 /**
					  * Set the stream #os to
					  * which data is written. If
					  * no stream is selected with
					  * this function, data goes
					  * to #deallog.
					  */
	void initialize_stream(std::ostream& stream);
	
					 /**
					  * Initialize the local data
					  * in the
					  * DoFInfo
					  * object used later for
					  * assembling.
					  *
					  * The second parameter is
					  * used to distinguish
					  * between the data used on
					  * cells and boundary faces
					  * on the one hand and
					  * interior faces on the
					  * other. Interior faces may
					  * require additional data
					  * being initialized.
					  */
	template <int dim>
	void initialize_info(DoFInfo<dim>& info,
			     bool interior_face);
	
					 /**
					  * Assemble the local values
					  * into the global vectors.
					  */
	template<int dim>
	void assemble(const DoFInfo<dim>& info);
	
					 /**
					  * Assemble both local values
					  * into the global vectors.
					  */
	template<int dim>
	void assemble(const DoFInfo<dim>& info1,
		      const DoFInfo<dim>& info2);

					 /**
					  * The value of the ith entry
					  * in #results.
					  */
	number operator() (unsigned int i) const;
      private:
					 /**
					  * Write the object T either
					  * to the stream #os, if the
					  * pointer is nonzero, or to
					  * #deallog if no pointer has
					  * been set.
					  */
	template<typename T>
	void write(const T& t) const;
	
	unsigned int patch_dimension;
	unsigned int n_vectors;
	unsigned int n_points;
	
	std::ostream* os;
    };

//----------------------------------------------------------------------//
    
    template<typename T>
    inline void
    GnuplotPatch::write(const T& d) const
    {
      if (os == 0)
	deallog << d;
      else
	(*os) << d;
    }
    
    
    inline void
    GnuplotPatch::GnuplotPatch()
		    :
		    os(0)
    {}
    
    
    inline void
    GnuplotPatch::initialize(unsigned int dim, unsigned int nv, unsigned int np)
    {
      patch_dimension = dim;
      n_vectors = n;
      n_points = np;
    }

    
    inline void
    GnuplotPatch::initialize_stream(std::ostream& stream)
    {
      os = &stream;
    }


    template <int dim>
    inline void
    GnuplotPatch::initialize_info(DoFInfo<dim>& info,
				  bool interior_face)
    {
      AssertDimension(patch_dimension,dim);
      info.initialize_quadrature(n_points, n_vectors);
    }

    
    template <int dim>
    inline void
    GnuplotPatch::assemble(const DoFInfo<dim>& info)
    {
      
    }
    
    
    template <int dim>
    inline void
    GnuplotPatch::assemble(const DoFInfo<dim>& info1, const DoFInfo<dim>& info2)
    {
      
    }
  }
}

#endif
