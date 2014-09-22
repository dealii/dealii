// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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


/**
@defgroup geodynamics The geodynamics demonstration suite

deal.II's @ref Tutorial "tutorial" contains a set of programs that together
form the geodynamics demonstration suite. The idea of these programs is to
demonstrate techniques for advanced finite element software using
applications from geodynamics, i.e. the investigation of processes in the
solid earth. By doing so, these programs are supposed to provide a basis
for more specialized, dedicated programs that can solve actual geodynamics
problems, for example as part of the work of graduate students or
postdocs. A more thorough discussion of the motivation for these programs
follows below.

Currently, the geodynamics testsuite contains the following
programs:

- step-8: Elasticity
- step-16: A %parallel elasticity solver
- step-20: Porous media flow
- step-21: Multiphase flow through porous media
- step-22: Stokes flow
- step-31: Thermal convection (Boussinesq flow)
- step-32: A %parallel Boussinesq solver for mantle convection

Some of these programs were developed under contract from the California
Institute of Technology with support by the National Science Foundation
under Award No. EAR-0426271, the grant that funded the <a target="_top"
href="http://www.geodynamics.org">Computational Infrastructure in
Geodynamics</a> initiative. The recipient, Wolfgang Bangerth, gratefully
acknowledges this source of support.


<h3>Rationale</h3>

Adaptive mesh refinement (AMR) has long been identified as a key technology
that would aid in the accurate and efficient numerical solution of a number of
geodynamics applications. It has been discussed in the geodynamics community
for several years and has been a continuous topic on the task list of CIG
since its inception. Yet, relatively little has happened in this direction so
far. Only recently have there been attempts to use AMR in geodynamics: CIG
sponsored a workshop on AMR technique in Boulder in October 2007, and a
collaboration between George Biros, Omar Ghattas, Mike Gurnis, and Shijie
Zhong's groups is currently developing a %parallel adaptive mantle convection
solver.

One of the reasons for the slow adoption of AMR techniques in geodynamics is
the relatively steep initial hurdle: codes have to provide the data structures
and algorithms to deal with adaptive meshes, finite elements have to be able
to deal with hanging nodes, etc. To do so efficiently and in sufficient
generality adds several 10,000 lines of code to finite element programs, too
much for the average student to do within the time frame of a dissertation. On
the other hand, there are libraries that provide the infrastructure code on
which applications supporting AMR can rapidly be built. deal.II
of course provides exactly this infrastructure.

The goal of the geodynamics testsuite is to write programs for a variety of
topics relevant to geodynamics. Continuing in the style of the existing tutorial
programs -- an extensive introduction explaining the background and
formulation of an application as well as the concepts of the numerical scheme
used in its solution; detailed comments throughout the code explaining
implementation details; and a section showing numerical results -- we intend to
provide the resulting programs as well-documented applications solving model
problems. In particular, they are aimed at the following goals:
<ul>
<li> <i>Starting points:</i> The existing tutorial of deal.II has proven to
  be an excellent starting point for graduate students and researchers to
  jump-start developing their own applications. By providing programs that are
  already close to the targeted application, first results can often be
  obtained very quickly, both maintaining the initial enthusiasm during
  development as well as allowing to spend research time on implementing
  application specific behavior rather than using months of work on basic
  infrastructure code supporting AMR.

  Supporting this point is the fact that although there are currently at least
  170 publications presenting results obtained with deal.II, we are aware of
  only a handful of applications that have been built with deal.II from
  scratch; all others have started as modifications of one of the tutorial
  programs.

<li> <i>Training:</i> The tutorial programs we propose to write will
  provide students and researchers with a reference implementation of current
  numerical technology such as AMR, higher order elements, sophisticated
  linear and nonlinear solvers, stabilization techniques, etc. Providing these
  as starting points for further development by others will also serve the
  goal of training a new generation of geodynamicists in modern numerical
  algorithms.

<li> <i>Extending equations and formulations:</i> In deal.II, it is fairly
  simple to extend a set of equations by another equation, for example an
  additional advected quantity that enters the existing equations as a right
  hand side or in one of the coefficients. Since applications typically use
  blocked matrices rather than the one-big-matrix-for-everything approach, it
  is also not complicated to find suitable linear solvers for augmented
  equations. Consequently, deal.II is a good tool for trying out more complex
  formulations of problems, or more complete models and their effects on the
  accuracy of solutions.

<li> <i>Rapid prototyping and benchmarking:</i> deal.II provides many
  interchangeable components that allow rapid prototyping of finite element
  kinds and orders, stabilization techniques, or linear solvers. For example,
  typically only a few lines of code have to be changed to replace low-order
  by high-order elements. Through this, it becomes relatively simple to try
  out higher order elements, a different block elimination solver, or a
  different stabilization technique. In turn, this may help in benchmarking
  applications both regarding computing times to solve as well as concerning
  the accuracy of numerical solutions.

  The applications in this module will already have been benchmarked for
  correctness. Existing tutorial programs typically employ simpler rather than
  more complicated solver schemes for exposition but frequently suggest more
  complicated schemes including hints on how they might be implemented in an
  appendix. 

<li> <i>Try algorithms:</i> The rapid prototyping abilities of deal.II may
  also help in determining best algorithms on the scale of programs to which
  deal.II is applicable, and then to implement this particular algorithm
  (without the ability to change it easily) in a dedicated program that can
  run on larger scale machines. For example, a small mantle convection code
  built on deal.II may be used to determine whether second order elements are
  useful for this purpose (see, for example, the results shown in
  step-31). If so, then one may use this result to implement
  second, rather than first, order elements in dedicated, large-scale mantle
  convection codes such as that which
  Ghattas and Zhong are building and that may run on 10,000s of processors, a
  range currently unattainable by deal.II.
</ul>


*/
