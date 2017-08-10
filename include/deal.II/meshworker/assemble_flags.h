// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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

#ifndef dealii__mesh_worker_assemble_flags_h
#define dealii__mesh_worker_assemble_flags_h


#include <deal.II/base/config.h>
#include <vector>


DEAL_II_NAMESPACE_OPEN

/*!@addtogroup MeshWorker */
/*@{*/

namespace MeshWorker
{

  /**
   * The enum type given to the mesh_loop() function, telling that function which
   * elements need to be assembled.
   *
   * You can select more than one flag by concatenation using the bitwise or
   * `operator|(AssembleFlags,AssembleFlags)`.
   *
   * @author Luca Heltai, 2017.
   */
  enum AssembleFlags
  {
    //! No update
    assemble_default = 0,
    //! Own cells
    /**
     * Assemble on cells.
     */
    assemble_own_cells = 0x0001,
    //! Ghost cells
    /**
     * Assemble on ghost cells.
     */
    assemble_ghost_cells = 0x0002,
    //! Own faces once
    /**
     * Assemble on own interior faces, visiting each face only once
     */
    assemble_own_interior_faces_once = 0x0004,
    //! Own faces both
    /**
     * Assemble on own interior faces, visiting each interior face twice
     */
    assemble_own_interior_faces_both = 0x0008,
    //! Ghost faces once
    /**
     * Assemble on faces between a locally owned cell and a ghost cell, making
     * sure that only one of the processes will assemble these faces (from the
     * finer side or the process with the lower mpi rank)
     */
    assemble_ghost_faces_once = 0x0010,
    //! Ghost faces both
    /**
     * Assemble on faces between a locally owned cell and a ghost cell, both
     * processes will assemble these faces. Note that these faces are never
     * assembled from both sides on a single process.
     */
    assemble_ghost_faces_both = 0x0020,
    //! Assemble cells before faces
    /**
     * Assemble cell integrals before face integrals.
     */
    assemble_cells_first = 0x0040,
    //! Assemble boundary faces
    /**
     * Assemble on boundary faces
     */
    assemble_boundary_faces = 0x0080,
    //! Assemble own interior aces
    /**
     * Assemble own interior faces, either interior ones or on the boundary.
     */
    assemble_own_interior_faces = assemble_own_interior_faces_both | assemble_own_interior_faces_once,
    //! Assemble own faces
    /**
     * Assemble own faces, either interior ones or on the boundary.
     */
    assemble_own_faces = assemble_own_interior_faces | assemble_boundary_faces,
    //! Assemble ghost faces
    /**
     * Assemble ghost faces
     */
    assemble_ghost_faces = assemble_ghost_faces_both | assemble_ghost_faces_once,
  };


  /**
   * Output operator which outputs update flags as a set of or'd text values.
   *
   * @ref AssembleFlags
   */
  template <class StreamType>
  inline
  StreamType &operator << (StreamType &s, AssembleFlags u)
  {
    s << " AssembleFlags";
    if (u & assemble_own_cells        ) s << "|own_cells"        ;
    if (u & assemble_own_interior_faces_once   ) s << "|own_faces_once"   ;
    if (u & assemble_own_interior_faces_both   ) s << "|own_faces_both"   ;
    if (u & assemble_ghost_cells      ) s << "|ghost_cells"      ;
    if (u & assemble_ghost_faces_once ) s << "|ghost_faces_once" ;
    if (u & assemble_ghost_faces_both ) s << "|ghost_faces_both" ;
    if (u & assemble_boundary_faces   ) s << "|boundary_faces"   ;
    return s;
  }


  /**
   * Global operator which returns an object in which all bits are set which are
   * either set in the first or the second argument. This operator exists since
   * if it did not then the result of the bit-or <tt>operator |</tt> would be an
   * integer which would in turn trigger a compiler warning when we tried to
   * assign it to an object of type AssembleFlags.
   *
   * @ref AssembleFlags
   */
  inline
  AssembleFlags
  operator | (AssembleFlags f1, AssembleFlags f2)
  {
    return static_cast<AssembleFlags> (
             static_cast<unsigned int> (f1) |
             static_cast<unsigned int> (f2));
  }




  /**
   * Global operator which sets the bits from the second argument also in the
   * first one.
   *
   * @ref AssembleFlags
   */
  inline
  AssembleFlags &
  operator |= (AssembleFlags &f1, AssembleFlags f2)
  {
    f1 = f1 | f2;
    return f1;
  }


  /**
   * Global operator which returns an object in which all bits are set which are
   * set in the first as well as the second argument. This operator exists since
   * if it did not then the result of the bit-and <tt>operator &</tt> would be
   * an integer which would in turn trigger a compiler warning when we tried to
   * assign it to an object of type AssembleFlags.
   *
   * @ref AssembleFlags
   */
  inline
  AssembleFlags
  operator & (AssembleFlags f1, AssembleFlags f2)
  {
    return static_cast<AssembleFlags> (
             static_cast<unsigned int> (f1) &
             static_cast<unsigned int> (f2));
  }


  /**
   * Global operator which clears all the bits in the first argument if they are
   * not also set in the second argument.
   *
   * @ref AssembleFlags
   */
  inline
  AssembleFlags &
  operator &= (AssembleFlags &f1, AssembleFlags f2)
  {
    f1 = f1 & f2;
    return f1;
  }
}

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif
