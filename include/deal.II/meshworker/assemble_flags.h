// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mesh_worker_assemble_flags_h
#define dealii_mesh_worker_assemble_flags_h


#include <deal.II/base/config.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup MeshWorker
 * @{
 */

namespace MeshWorker
{
  /**
   * The enum type given to the mesh_loop() function, telling that function
   * which elements need to be assembled.
   *
   * You can select more than one flag by concatenation using the bitwise or
   * <code>operator|(AssembleFlags,AssembleFlags)</code>.
   */
  enum AssembleFlags
  {
    /**
     * Do Nothing.
     */
    assemble_nothing = 0,
    /**
     * Assemble on locally owned cells.
     */
    assemble_own_cells = 0x0001,
    /**
     * Assemble on ghost cells.
     */
    assemble_ghost_cells = 0x0002,
    /**
     * Assemble on interior faces between two locally owned cells,
     * visiting each face only once.
     */
    assemble_own_interior_faces_once = 0x0004,
    /**
     * Assemble on interior faces between two locally owned cells,
     * visiting each interior face twice, once from each of the two
     * adjacent cells.
     */
    assemble_own_interior_faces_both = 0x0008,
    /**
     * Assemble on faces between a locally owned cell and a ghost cell, making
     * sure that only one of the processes will assemble these faces (from the
     * finer side or the process with the lower mpi rank).
     */
    assemble_ghost_faces_once = 0x0010,
    /**
     * Assemble on faces between a locally owned cell and a ghost cell. Both
     * processes will assemble these faces. Note that they are never
     * assembled from both sides on a single process.
     */
    assemble_ghost_faces_both = 0x0020,
    /**
     * Assemble on boundary faces of the locally owned cells.
     */
    assemble_boundary_faces = 0x0040,

    /**
     * By default we assemble cell integrals before face integrals. If this
     * flag is specified, cells will be assembled after faces and boundaries.
     */
    cells_after_faces = 0x0080,

    /**
     * Combination of flags to determine if any work on cells is done.
     */
    work_on_cells = assemble_own_cells | assemble_ghost_cells,

    /**
     * Combination of flags to determine if any work is done on faces.
     */
    work_on_faces = assemble_own_interior_faces_once |
                    assemble_own_interior_faces_both |
                    assemble_ghost_faces_once | assemble_ghost_faces_both,

    /**
     * Combination of flags to determine if any work is done on the boundary
     * faces.
     */
    work_on_boundary = assemble_boundary_faces,
  };


  /**
   * Output operator which outputs assemble flags as a set of or'd text values.
   *
   * @ref AssembleFlags
   */
  template <typename StreamType>
  inline StreamType &
  operator<<(StreamType &s, AssembleFlags u)
  {
    s << " AssembleFlags";
    if (u & assemble_own_cells)
      s << "|own_cells";
    if (u & assemble_own_interior_faces_once)
      s << "|own_faces_once";
    if (u & assemble_own_interior_faces_both)
      s << "|own_faces_both";
    if (u & assemble_ghost_cells)
      s << "|ghost_cells";
    if (u & assemble_ghost_faces_once)
      s << "|ghost_faces_once";
    if (u & assemble_ghost_faces_both)
      s << "|ghost_faces_both";
    if (u & assemble_boundary_faces)
      s << "|boundary_faces";
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
  inline AssembleFlags
  operator|(AssembleFlags f1, AssembleFlags f2)
  {
    return static_cast<AssembleFlags>(static_cast<unsigned int>(f1) |
                                      static_cast<unsigned int>(f2));
  }



  /**
   * Global operator which sets the bits from the second argument also in the
   * first one.
   *
   * @ref AssembleFlags
   */
  inline AssembleFlags &
  operator|=(AssembleFlags &f1, AssembleFlags f2)
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
  inline AssembleFlags
  operator&(AssembleFlags f1, AssembleFlags f2)
  {
    return static_cast<AssembleFlags>(static_cast<unsigned int>(f1) &
                                      static_cast<unsigned int>(f2));
  }


  /**
   * Global operator which clears all the bits in the first argument if they are
   * not also set in the second argument.
   *
   * @ref AssembleFlags
   */
  inline AssembleFlags &
  operator&=(AssembleFlags &f1, AssembleFlags f2)
  {
    f1 = f1 & f2;
    return f1;
  }
} // namespace MeshWorker

/** @} */

DEAL_II_NAMESPACE_CLOSE

#endif
