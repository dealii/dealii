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

#ifndef __deal2__mesh_worker_vector_info_h
#define __deal2__mesh_worker_vector_info_h

#include <deal.II/base/named_data.h>
#include <deal.II/lac/block_indices.h>
#include <deal.II/numerics/mesh_worker_info.h>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  /**
   * This class is a simple object that can be used in
   * MeshWorker::loop(). For all vectors provided in the argument to
   * initialize_data(), this class generates the local data vector on a
   * cell specified by reinit().
   *
   * This objects of class as part of VectorInfoBox conform to the
   * interface for INFOBOX laid out in the documentation of MeshWorker.
   *
   * See IntegrationInfo for an alternative.
   *
   * @author Guido Kanschat
   * @date 2012
   */
  class VectorInfo
  {
  public:
    /**
    * Initialize the data
    * vector and cache the
    * selector.
    */
    void initialize_data(const NamedData<BlockVector<double>*> &data);

    /**
    * Reinitialize the local data
    * vectors to represent the
    * information on the given cell.
    */
    template<int dim, int spacedim>
    void reinit(const MeshWorker::DoFInfo<dim, spacedim> &cell);

    /**
    * Return local data vector
    * <i>i</i>, which selects the
    * local data from vector
    * <i>i</i> of the
    * <code>data</code> object
    * used in initialize_data().
    */
    BlockVector<double> &operator() (unsigned int i);

  private:
    std::vector<BlockVector<double> > local_data;
    /**
    * The global data vector
    * used to compute function
    * values in quadrature
    * points.
    */
    SmartPointer<const NamedData<BlockVector<double>*> > global_data;
  };


  /**
   * A class conforming to the INFOBOX interface in the documentation of
   * MeshWorker. This class provides local operators in
   * MeshWorker::loop() with the local values of the degrees of freedom
   * of global vectors, for instance for calculations based on cell
   * matrices.
   *
   * For an alternative providing the values in quadrature points see
   * IntegrationInfoBox.
   *
   * @author Guido Kanschat
   * @date 2012
   */
  class VectorInfoBox
  {
  public:
    typedef VectorInfo CellInfo;

    void initialize_data(const NamedData<BlockVector<double>*> &data);

    template <int dim, class DOFINFO>
    void post_cell(const MeshWorker::DoFInfoBox<dim, DOFINFO> &)
    {}

    template <int dim, class DOFINFO>
    void post_faces(const MeshWorker::DoFInfoBox<dim, DOFINFO> &)
    {}

    VectorInfo cell;
    VectorInfo boundary;
    VectorInfo face;
    VectorInfo subface;
    VectorInfo neighbor;
  };


  inline
  void
  VectorInfo::initialize_data(const NamedData<BlockVector<double>*> &data)
  {
    global_data = &data;
    local_data.resize(global_data->size());
  }


  template<int dim, int spacedim>
  inline
  void
  VectorInfo::reinit(const MeshWorker::DoFInfo<dim, spacedim> &i)
  {
    const NamedData<BlockVector<double>*> &gd = *global_data;

    for (unsigned int k=0; k<local_data.size(); ++k)
      {
        const BlockVector<double> &v = *gd(k);

        local_data[k].reinit(i.block_info->local());
        for (unsigned int j=0; j<local_data[k].size(); ++j)
          local_data[k](j) = v(i.indices[j]);
      }
  }


  inline
  BlockVector<double> &
  VectorInfo::operator() (unsigned int i)
  {
    AssertIndexRange(i, local_data.size());
    return local_data[i];
  }


  inline
  void
  VectorInfoBox::initialize_data(const NamedData<BlockVector<double>*> &data)
  {
    cell.initialize_data(data);
    boundary.initialize_data(data);
    face.initialize_data(data);
    subface.initialize_data(data);
    neighbor.initialize_data(data);
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif
