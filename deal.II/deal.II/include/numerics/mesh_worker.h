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
 * <h3>Template argument types</h3>
 *
 * The functions loop() and cell_action() take some arguments which
 * are template parameters. Let us list the minimum requirements for
 * these classes here and describe their properties.
 *
 * <h4>ITERATOR</h4>
 *
 * Any object that has an <tt>operator++()</tt> and points to a
 * TriaObjectAccessor.
 *
 * <h4>DOFINFO</h4>
 *
 * For an example implementation, refer to the class template DoFInfo.
 * In order to work with cell_action() and loop(), DOFINFO needs to
 * follow the following interface.
 * @code
 * class DOFINFO
 * {
 *   private:
 *     DOFINFO();
 *     DOFINFO(const DOFINFO&);
 *     DOFINFO& operator=(const DOFINFO&);
 *
 *   public:
 *     template <class CellIt>
 *     void reinit(const CellIt& c);
 *
 *     template <class CellIt, class FaceIt>
 *     void reinit(const CellIt& c, const FaceIt& f, unsigned int n);
 *
 *     template <class CellIt, class FaceIt>
 *     void reinit(const CellIt& c, const FaceIt& f, unsigned int n,
 *     unsigned int s);
 *
 *   friend template class DoFInfoBox<int dim, DOFINFO>;
 * };
 * @endcode
 *
 * The three private functions are called by DoFInfoBox and should not
 * be needed elsewhere. Obviously, they can be made public and then
 * the friend declaration at the end may be missing.
 *
 * Additionally, you will need at least one public constructor. Furthermore
 * DOFINFO is pretty useless yet: functions to interface with
 * INTEGRATIONINFO and ASSEMBLER are needed.
 *
 * DOFINFO objects are gathered in a DoFInfoBox. In those objects, we
 * store the results of local operations on each cel and its
 * faces. Once all this information has been gathered, an ASSEMBLER is
 * used to assemble it into golbal data.
 *
 * <h4>INFOBOX</h4>
 *
 * This type is exemplified in IntegrationInfoBox. It collects the
 * input data for actions on cells and faces in INFO objects (see
 * below). It provides the following interface to loop() and
 * cell_action():
 *
 * @code
 * class INFOBOX
 * {
 *   public:
 *     template <int dim, class DOFINFO>
 *     void post_cell(const DoFInfoBox<dim, DOFINFO>&);
 *
 *     template <int dim, class DOFINFO>
 *     void post_faces(const DoFInfoBox<dim, DOFINFO>&);
 *
 *     INFO cell;
 *     INFO boundary;
 *     INFO face;
 *     INFO subface;
 *     INFO neighbor;
 * };
 * @endcode
 *
 * The main purpose of this class is gathering the five INFO objects,
 * which contain the temporary data used on each cell or face. The
 * requirements on these objects are listed below. Here, we only note
 * that there need to be these 5 objects with the names listed above.
 *
 * The two function templates are call back functions called in
 * cell_action(). The first is called before the faces are worked on,
 * the second after the faces.
 *
 * <h4>INFO</h4>
 *
 * See IntegrationInfo for an example of these objects. They contain
 * the temporary data needed on each cell or face to compute the
 * result. The MeshWorker only uses the interface
 *
 * @code
 * class INFO
 * {
 *   public:
 *     void reinit(const DOFINFO& i);
 * };
 * @endcode
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
 * In order to save time, you can set the variables boundary_fluxes
 * and interior_fluxes of the base class to false, thus telling the
 * Meshworker::loop() not to loop over those faces.
 *
 * All the information in here is used to set up IntegrationInfo
 * objects correctly, typically in an IntegrationInfoBox.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
  template <int dim>
  class IntegrationWorker
  {
    public:
				       /**
					* The info type expected by a
					* cell integrator.
					*/
      typedef IntegrationInfo<dim, dim> CellInfo;

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
}


DEAL_II_NAMESPACE_CLOSE

#endif
