//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

/**
 * @page DEALSoftware Third party software deal.II interfaces with
 *
 * @note Names of software products are trademarks of their respective owners.
 *
 * This is a list of software packages that deal.II interacts with in some
 * way, either by supporting a data format of another program, or by actively
 * calling functions of other packages. Configuration of interfaces, if
 * necessary at all, is described in the deal.II ReadMe file.
 *
 *
 * @section SodtwarePrePost Pre- and Postprocessing
 *
 * @subsection SoftwareVis Visualization tools
 *
 * This is a list of visualization formats and software that deal.II
 * supports. Data in these formats is written by the DataOutBase and in parts
 * by the GridOut classes (see @ref output).
 *
 *
 * @subsubsection SoftwareAVS AVS Express
 *
 * deal.II reads and writes the UCD format specified in the AVS
 * Express user manual. See http://www.avs.com/ for more details on
 * this software.
 *
 * @subsubsection SoftwareGMV GMV
 *
 * Output for the General Mesh Viewer can be produced. See
 * http://www-xdiv.lanl.gov/XCM/gmv/GMVHome.html for details on GMV
 *
 * @subsubsection SoftwareGnuplot gnuplot
 *
 * Two-dimensional data can be written in a format suitable for
 * gnuplot, even on locally refined and unstructured meshes, See
 * http://www.gnuplot.info/
 *
 *
 * @subsubsection SoftwareOpenDX OpenDX
 *
 * The former IBM Visual Data Explorer, now an OpenSource project at
 * http://www.opendx.org/.
 *
 *
 * @subsubsection SoftwarePovray Povray
 *
 * While it is not actually taylored to scientific visualization, you
 * may be able to produce impressive pictures of three-dimensional
 * deformed bodies with http://www.povray.org/
 *
 *
 * @subsubsection SoftwareTecplot Tecplot
 *
 * deal.II writes textual and binary files for Tecplot. See
 * http://www.tecplot.com for more details on this software.
 *
 *
 * @subsubsection SoftwareXFig XFig
 *
 * Though not a visualization tool at all, you can nevertheless write grids in
 * XFig format 3.2. Those can be postprocessed manually within XFig and
 * written in many different graphics formats. See http://www.xfig.org/
 *
 *
 * @subsubsection SoftwareVTK VTK
 *
 * The Visualisation Toolkit VTK is an open format supported by a
 * number of visualization projects such as ParaView, VisIt, or
 * MayaVi. See http://public.kitware.com/VTK/ . The file format is described
 * at http://vtk.org/pdf/file-formats.pdf .
 *
 *
 *
 * @section SoftwareLibs Libraries used inside deal.II
 *
 * @subsection SoftwarePETSc PETSc
 *
 * PETSc is a library that supports, among other things, a large number of
 * linear algebra data structures and algorithms, in much the same way as we
 * do in the linear algebra classes of deal.II (see @ref LAC). However, PETSc
 * goes beyond what we have to offer in that it has more algorithms (for
 * example algebraic multigrid) and most importantly it works in parallel on
 * distributed memory clusters, using MPI.
 *
 * In order to support parallel computations in deal.II, we have
 * written interfaces to many PETSc functions and data structures in
 * the PETScWrapper namespace, that allow the use of PETSc in much the
 * same way as deal.II's own linear algebra classes are used. The use
 * of these wrappers is explained in the step-17 and step-18 example
 * programs, as well as in the @ref PETScWrappers module. The <a
 * href="../../readme-petsc-trilinos.html">PETSc and Trilinos
 * ReadMe</a> file explains how to configure deal.II to use PETSc.
 *
 * PETSc can be obtained from http://www.mcs.anl.gov/petsc/.
 *
 *
 * @subsection SoftwareTrilinos Trilinos
 *
 * Trilinos, like @ref SoftwarePETSc, is a library that supports, among other
 * things of interest to numerical computations, a large number of linear
 * algebra data structures and algorithms. It, too, can work on parallel
 * clusters.
 *
 * Interfaces to Trilinos exist in the TrilinosWrappers namespace,
 * making matrices, vectors, and solvers look like the corresponding
 * deal.II classes. Their use is explained in the @ref step_31
 * "step-31", step-32, and step-33
 * tutorial programs.  The <a
 * href="../../readme-petsc-trilinos.html">PETSc and Trilinos
 * ReadMe</a> file explains how to configure deal.II to use this
 * feature.
 *
 * Trilinos can be obtained from http://trilinos.sandia.gov.
 *
 *
 * @subsection SoftwareMETIS METIS
 *
 * METIS is a tool that allows to partition a graph into chunks of roughly
 * equal size. We use it to subdivide a domain into blocks that have about the
 * same number of cells, when distributing work for parallel programs. METIS
 * can be obtained from http://www-users.cs.umn.edu/~karypis/metis/index.html
 *
 *
 * @subsection SoftwareUMFPACK UMFPACK
 *
 * UMFPACK is a sparse direct solver and is included by permission with
 * deal.II distributions. To configure its available, read the README file.
 *
 * UMFPACK can be obtained from
 * http://www.cise.ufl.edu/research/sparse/umfpack/
 */
