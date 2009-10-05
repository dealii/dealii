//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2006, 2007, 2008, 2009 by Guido Kanschat
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__mesh_worker_h
#define __deal2__mesh_worker_h

#include <base/geometry_info.h>
#include <base/named_data.h>
#include <lac/block_indices.h>
#include <multigrid/mg_level_object.h>
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
 * This integration_loop() is complemented by a class
 * AssemblingIntegrator, which separates local integration from
 * assembling. In order to use it, follow this recipe: fist, create a
 * class responsible for the local integration:
 *
 * @code
 * template <int dim>
 * class LocalIntegrator : public Subscriptor
 * {
 *   void cell(typename MeshWorker::IntegrationWorker<dim>::CellInfo& info) const;
 *   void bdry(typename MeshWorker::IntegrationWorker<dim>::FaceInfo& info) const;
 *   void face(typename MeshWorker::IntegrationWorker<dim>::FaceInfo& info1,
 *             typename MeshWorker::IntegrationWorker<dim>::FaceInfo& info2) const;
 * };
 * @endcode
 * The three functions in there must have exactly this signature and
 * should do the integration on cells, boundary and interior faces,
 * respectively.
 *
 * Then, create the AssemblingIntegrator object, where the second
 * template argument decides on what kind of global data is
 * assembled. In the following, we decide to assemble a simple
 * SparseMatrix with no block structure:
 *
 * @code
 * LocalIntegrator<dim> loc;
 * MeshWorker::AssemblingIntegrator<dim, MeshWorker::Assembler::MatrixSimple<SparseMatrix<double> >, LocalIntegrator<dim> >
 *   integrator(loc);
 * @endcode
 *
 * Before we can run the integration loop, we have to initialize
 * several data structures in our AssemblingIntegrator. For instance,
 * we have to decide on the quadrature rule or we may need more than
 * the default update flags.
 *
 * @code
 * integrator.initialize_gauss_quadrature(2,2,2);
 * integrator.add_update_flags(update_values | update_gradients, true, true, true, true);
 * integrator.initialize(system_matrix);
 * @endcode
 *
 * Finally, we set up the structures needed by the integration_loop()
 * and run it.
 *
 * @code
 * MeshWorker::IntegrationInfoBox<dim> info_box(dof_handler);
 * info_box.initialize(integrator, fe, mapping);
 *
 * MeshWorker::integration_loop(dof_handler.begin_active(), dof_handler.end(),
 *                              info_box, integrator);
 * @endcode
 *
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
      void bdry(DoFInfo<dim>& face);
      
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
 */
  template <int dim>
  class IntegrationWorker : public LocalWorker<dim>
  {
    public:
				       /**
					* The info type expected by a
					* cell integrator.
					*/
      typedef IntegrationInfo<dim, FEValuesBase<dim, dim>, dim> CellInfo;
      
				       /**
					* The info type expected by a
					* face integrator.
					*/
      typedef IntegrationInfo<dim, FEFaceValuesBase<dim, dim>, dim> FaceInfo;
      
				       /**
					* Initialize default values.
					*/
      IntegrationWorker();
				       /**
					* Initialize the
					* VectorSelector objects
					* #cell_selector,
					* #bdry_selector and
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
      void initialize_selectors(const VectorSelector& cell_selector,
				const VectorSelector& bdry_selector,
				const VectorSelector& face_selector);

				       /**
					* Add a vector to some or
					* all selectors.
					*/
      void add_selector(const std::string& name, bool values, bool gradients, bool hessians,
			bool cell, bool bdry, bool face);

				       /**
					* Add additional values for update.
					*/
      void add_update_flags(const UpdateFlags flags, bool cell = true,
			    bool bdry = true, bool face = true,
			    bool ngbr = true);
	
				       /** Assign n-point Gauss
					* quadratures to each of the
					* quadrature rules. Here, a
					* size of zero points means
					* that no loop over these grid
					* entities should be
					* performed.
					*/
      void initialize_gauss_quadrature(unsigned int n_cell_points,
				       unsigned int n_bdry_points,
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
      MeshWorker::VectorSelector bdry_selector;
	
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
      UpdateFlags bdry_flags;
	
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
      UpdateFlags ngbr_flags;
	
				       /**
					* The quadrature rule used
					* on cells.
					*/
      Quadrature<dim> cell_quadrature;
      
				       /**
					* The quadrature rule used
					* on boundary faces.
					*/
      Quadrature<dim-1> bdry_quadrature;
      
				       /**
					* The quadrature rule used
					* on interior faces.
					*/
      Quadrature<dim-1> face_quadrature;
  };
      
  
/**
 * A worker class for integration_loop(), separating assembling and
 * local integration. Objects of this class rely on the base class
 * provided by the template argument ASSEMBLER for assembling local
 * data into global. The integration of local data is delegated to the
 * other template parameter class INTEGRATOR.
 *
 * In order to use this class, it will be necessary to create an
 * INTEGRATOR class following this template:
 *
 * @code
 * template <int dim>
 * class MyLocalOperator : public Subscriptor
 * {
 *   public:
 *     void cell(typename IntegrationWorker<dim>::CellInfo& info) const;
 *     void bdry(typename IntegrationWorker<dim>::FaceInfo& info) const;
 *     void face(typename IntegrationWorker<dim>::FaceInfo& info1,
 *               typename IntegrationWorker<dim>::FaceInfo& info2) const;
 * };
 * @endcode
 *
 * This class will do whatever your problem requires locally on each
 * cell and/or face. Once this class is defined, you choose a suitable
 * assembler for your problem from the Assembler namespace and set up
 * objects:
 *
 * @code
 * MyLocalOperator<dim> myop;
 *
 * AssemblingIntegrator<dim, Assembler::MyAssembler, MyLocalOperator<dim> >
 *   integrator(myop);
 * @endcode
 *
 * You do the necessary initializations of this @p integrator and then
 * you have a worker object suitable for integration_loop().
 *
 * @author Guido Kanschat, 2009
 */
  template <int dim, class ASSEMBLER, class INTEGRATOR>
  class AssemblingIntegrator :
      public IntegrationWorker<dim>,
      public ASSEMBLER
  {
    public:
				       /**
					* Constructor, initializing
					* the local data.
					*/
      AssemblingIntegrator(const INTEGRATOR& local);
				       /**
					* The cell operator called by
					* integration_loop().
					*/
      void cell(typename IntegrationWorker<dim>::CellInfo& info);
      
				       /**
					* The local boundary operator
					* called by
					* integration_loop().
					*/
      void bdry(typename IntegrationWorker<dim>::FaceInfo& info);
      
				       /**
					* The interior face operator
					* called by
					* integration_loop().
					*/
      void face(typename IntegrationWorker<dim>::FaceInfo& info1,
		typename IntegrationWorker<dim>::FaceInfo& info2);
      
    private:
				       /**
					* Pointer to the object doing
					* local integration.
					*/
      SmartPointer<const INTEGRATOR> local;
  };
  
//----------------------------------------------------------------------//
  template <int dim>
  inline
  LocalWorker<dim>::LocalWorker()
		  :
		  interior_fluxes(true), boundary_fluxes(true)
  {}

//----------------------------------------------------------------------//

  template <int dim, class ASSEMBLER, class INTEGRATOR>
  inline
  AssemblingIntegrator<dim,ASSEMBLER,INTEGRATOR>::AssemblingIntegrator(const INTEGRATOR& local)
		  :
		  local(&local, typeid(*this).name())
  {}

  
  template <int dim, class ASSEMBLER, class INTEGRATOR>
  inline void
  AssemblingIntegrator<dim,ASSEMBLER,INTEGRATOR>::cell(
    typename IntegrationWorker<dim>::CellInfo& info)
  {
    local->cell(info);
    ASSEMBLER::assemble(info);
  }

  
  template <int dim, class ASSEMBLER, class INTEGRATOR>
  inline void
  AssemblingIntegrator<dim,ASSEMBLER,INTEGRATOR>::bdry(
    typename IntegrationWorker<dim>::FaceInfo& info)
  {
    local->bdry(info);
    ASSEMBLER::assemble(info);
  }

  
  template <int dim, class ASSEMBLER, class INTEGRATOR>
  inline void
  AssemblingIntegrator<dim,ASSEMBLER,INTEGRATOR>::face(
    typename IntegrationWorker<dim>::FaceInfo& info1,
    typename IntegrationWorker<dim>::FaceInfo& info2)
  {
    local->face(info1, info2);
    ASSEMBLER::assemble(info1, info2);
  }  
}

DEAL_II_NAMESPACE_CLOSE

#endif
