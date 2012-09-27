//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2006, 2007, 2008, 2009, 2010, 2011, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__mesh_worker_local_integrator_h
#define __deal2__mesh_worker_local_integrator_h

#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/std_cxx1x/function.h>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  template <int dim, int spacedim, typename number> class DoFInfo;
  template <int dim, int spacedim> class IntegrationInfo;
  
/**
 * A local integrator object, which can be used to simplify the call
 * of loop(). Instead of providing the three local integration
 * functions separately, we bundle them as virtual functions in this
 * class.
 *
 * Additionally, since we cannot have a virtual null function, we
 * provide flags, which allow us to indicate, whether we want to
 * integrate on boundary and interior faces. Thes flags are true by
 * default, but can be modified by applications to speed up the loop.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat
 * @date 2012
 */
  template <int dim, int spacedim=dim, typename number=double>
    class LocalIntegrator : public Subscriptor
  {
    public:
				       /**
					* The constructor setting
					* default values.
					*/
      LocalIntegrator();

				       /**
					* The empty virtual destructor.
					*/
      ~LocalIntegrator();
      
				       /**
					* Virtual function for
					* integrating on cells.
					*/
      virtual void cell(DoFInfo<dim, spacedim, number>& dinfo,
			IntegrationInfo<dim, spacedim>& info) const = 0;
				       /**
					* Virtual function for
					* integrating on boundary faces.
					*/
      virtual void boundary(DoFInfo<dim, spacedim, number>& dinfo,
			    IntegrationInfo<dim, spacedim>& info) const = 0;
				       /**
					* Virtual function for
					* integrating on interior faces.
					*/
      virtual void face(DoFInfo<dim, spacedim, number>& dinfo1,
			DoFInfo<dim, spacedim, number>& dinfo2,
			IntegrationInfo<dim, spacedim>& info1,
			IntegrationInfo<dim, spacedim>& info2) const = 0;

				       /**
					* The flag indicating whether
					* the cell integrator cell()
					* is to be used in the
					* loop. Defaults to
					* <tt>true</tt>.
					*/
      bool use_cell;

				       /**
					* The flag indicating whether
					* the boundary integrator
					* boundary() is to be
					* used in the loop. Defaults
					* to <tt>true</tt>.
					*/
      bool use_boundary;

				       /**
					* The flag indicating whether
					* the interior face integrator
					* face() is to be used in the
					* loop. Defaults to
					* <tt>true</tt>.
					*/
      bool use_face;
      
  };
}



DEAL_II_NAMESPACE_CLOSE

#endif
