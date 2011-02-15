//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2003, 2004, 2007, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

// Provide scratch data objects with no contents for the MeshWorker::loop()

#include <numerics/mesh_worker_info.h>

class EmptyInfo
{
  public:
    template <class DOFINFO>
    void reinit(const DOFINFO&)
      {}
};


class EmptyInfoBox
{
   public:
    typedef EmptyInfo CellInfo;
    template <int dim, class DOFINFO>
    void post_cell(const MeshWorker::DoFInfoBox<dim, DOFINFO>&)
      {}
    
    template <int dim, class DOFINFO>
    void post_faces(const MeshWorker::DoFInfoBox<dim, DOFINFO>&)
      {}
    
    EmptyInfo cell;
    EmptyInfo boundary;
    EmptyInfo face;
    EmptyInfo subface;
    EmptyInfo neighbor;
};


