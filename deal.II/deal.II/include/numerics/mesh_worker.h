//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__mesh_worker_h
#define __deal2__mesh_worker_h

#include <base/config.h>
#include <base/std_cxx1x/function.h>
#include <base/geometry_info.h>
#include <base/named_data.h>
#include <lac/block_indices.h>
#include <base/mg_level_object.h>
#include <numerics/mesh_worker_vector_selector.h>
#include <numerics/mesh_worker_info.h>
#include <numerics/mesh_worker_assembler.h>


DEAL_II_NAMESPACE_OPEN

template<int,int> class DoFHandler;
template<int,int> class MGDoFHandler;

/**
 * A collection of functions and classes for the mesh loops that are
 * an ubiquitous part of each finite element program.
 *
 * The workhorse of this namespace is the loop() function, which implements a
 * completely generic loop over all mesh cells.
 *
 * The loop() depends on certain objects handed to it as
 * arguments. These objects are of two types, info objects like
 * DoFInfo and IntegrationInfo and worker objects like LocalWorker and
 * IntegrationWorker.
 *
 * Worker objects usually do two different jobs: first, they compute
 * the local contribution of a cell or face to the global
 * operation. Second, they assemble this local contribution into the
 * global result, whether a functional, a form or a bilinear
 * form. While the first job is particular to the problem being
 * solved, the second is generic and only depends on the data
 * structures. Therefore, base classes for workers assembling into
 * global data are provided in the namespace Assembler.
 *
 * <h3>Simplified interfaces</h3>
 *
 * Since the loop() is fairly general, a specialization
 * integration_loop() is available, which is a wrapper around loop()
 * with a simplified interface.
 *
 * The integration_loop() function loop takes most of the information
 * that it needs to pass to loop() from an IntegrationInfoBox
 * object. Its use is explained in step-12, but in
 * short it requires functions that do the local integration on a
 * cell, interior or boundary face, and it needs an object (called
 * "assembler") that copies these local contributions into the global
 * matrix and right hand side objects.
 *
 * Before we can run the integration loop, we have to initialize
 * several data structures in our IntegrationWorker and assembler
 * objects. For instance, we have to decide on the quadrature rule or
 * we may need more than the default update flags.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
namespace MeshWorker
{
/**
 * Template for a class for the objects doing the actual work on cells
 * and faces.
 *
 * This class can serve as a base class for the actual worker class,
 * since, while we do not use virtual functions, we provide the
 * necessary interface for the mesh loops and the DoFInfo
 * class here. Thus, they do not have to be reprogrammed.
 *
 * In particular, the mesh loops will require data elements
 * #interior_fluxes and #boundary_fluxes to determine whether the
 * loop over faces will be started at all.
 *
 * @author Guido Kanschat, 2009
 */
  template <int dim>
  class LocalWorker
  {
    public:
				       /**
					* Constructor, setting
					* #interior_fluxes and
					* #boundary_fluxes to @p
					* true.
					*/
      LocalWorker ();

				       /**
					* Do the work on a cell.
					*/
      void cell(DoFInfo<dim>& cell);

				       /**
					* Do the work on a boundary face.
					*/
      void boundary(DoFInfo<dim>& face);

				       /**
					* Do the work on an interior face.
					*/
      void face(DoFInfo<dim>& face1, DoFInfo<dim>& face2);

				       /**
					* Computations on interior
					* faces are necessary.
					*/
      bool interior_fluxes;

				       /**
					* Computations on interior
					* faces are necessary.
					*/
      bool boundary_fluxes;
  };

/**
 * Worker object for integration of functionals, residuals or matrices.
 *
 * In order to allow for sufficient generality, a few steps have to be
 * undertaken to use this class. The constructor will only fill some
 * default values.
 *
 * First, you should consider if you need values from any vectors in a
 * NamedData object. If so, fill the VectorSelector objects
 * #cell_selector, #boundary_selector and #face_selector with their names
 * and the data type (value, gradient, Hessian) to be extracted.
 *
 * Afterwards, you will need to consider UpdateFlags for FEValues
 * objects. A good start is initialize_update_flags(), which looks at
 * the selectors filled before and adds all the flags needed to get
 * the selection. Additional flags can be set with add_update_flags().
 *
 * Finally, we need to choose quadrature formulas. If you choose to
 * use Gauss formulas only, use initialize_gauss_quadrature() with
 * appropriate values. Otherwise, you can fill the variables
 * #cell_quadrature, #boundary_quadrature and #face_quadrature directly.
 *
 * In order to save time, you can set the variables #boundary_fluxes
 * and #interior_fluxes of the base class to false, thus telling the
 * Meshworker::loop() not to loop over those faces.
 *
 * All the information in here is used to set up IntegrationInfo
 * objects correctly, typically in an IntegrationInfoBox.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
  template <int dim>
  class IntegrationWorker : public LocalWorker<dim>
  {
    public:
				       /**
					* The info type expected by a
					* cell integrator.
					*/
      typedef IntegrationInfo<dim, dim> CellInfo;

				       /**
					* The info type expected by a
					* face integrator.
					*/
      typedef IntegrationInfo<dim, dim> FaceInfo;

				       /**
					* Initialize default values.
					*/
      IntegrationWorker();
				       /**
					* Initialize the
					* VectorSelector objects
					* #cell_selector,
					* #boundary_selector and
					* #face_selector in order to
					* save computational
					* eeffort. If no selectors
					* are used, then values for
					* all named vectors in
					* DoFInfo::global_data will be
					* computed in all quadrature
					* points.
					*
					* This function will also
					* add UpdateFlags to the
					* flags stored in this class.
					*/
      void initialize_update_flags();

				       /**
					* Add additional values for update.
					*/
      void add_update_flags(const UpdateFlags flags, bool cell = true,
			    bool boundary = true, bool face = true,
			    bool neighbor = true);

				       /** Assign n-point Gauss
					* quadratures to each of the
					* quadrature rules. Here, a
					* size of zero points means
					* that no loop over these grid
					* entities should be
					* performed.
					*/
      void initialize_gauss_quadrature(unsigned int n_cell_points,
				       unsigned int n_boundary_points,
				       unsigned int n_face_points);

				       /**
					* Select the vectors from
					* DoFInfo::global_data
					* that should be computed in
					* the quadrature points on cells.
					*/
      MeshWorker::VectorSelector cell_selector;

				       /**
					* Select the vectors from
					* DoFInfo::global_data
					* that should be computed in
					* the quadrature points on
					* boundary faces.
					*/
      MeshWorker::VectorSelector boundary_selector;

				       /**
					* Select the vectors from
					* DoFInfo::global_data
					* that should be computed in
					* the quadrature points on
					* interior faces.
					*/
      MeshWorker::VectorSelector face_selector;

				       /**
					* The set of update flags
					* for boundary cell integration.
					*
					* Defaults to
					* #update_JxW_values.
					*/
      UpdateFlags cell_flags;
				       /**
					* The set of update flags
					* for boundary face integration.
					*
					* Defaults to
					* #update_JxW_values and
					* #update_normal_vectors.
					*/
      UpdateFlags boundary_flags;

				       /**
					* The set of update flags
					* for interior face integration.
					*
					* Defaults to
					* #update_JxW_values and
					* #update_normal_vectors.
					*/
      UpdateFlags face_flags;

				       /**
					* The set of update flags
					* for interior face integration.
					*
					* Defaults to
					* #update_default, since
					* quadrature weights are
					* taken from the other cell.
					*/
      UpdateFlags neighbor_flags;

				       /**
					* The quadrature rule used
					* on cells.
					*/
      Quadrature<dim> cell_quadrature;

				       /**
					* The quadrature rule used
					* on boundary faces.
					*/
      Quadrature<dim-1> boundary_quadrature;

				       /**
					* The quadrature rule used
					* on interior faces.
					*/
      Quadrature<dim-1> face_quadrature;
  };


//----------------------------------------------------------------------//
  template <int dim>
  inline
  LocalWorker<dim>::LocalWorker()
		  :
		  interior_fluxes(true), boundary_fluxes(true)
  {}
}


DEAL_II_NAMESPACE_CLOSE

#endif
