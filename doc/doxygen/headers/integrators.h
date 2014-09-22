// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
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
 * @defgroup Integrators Integrators
 *
 * A collection of namespaces and functions which simplify the coding
 * of forms and bilinear forms on finite element spaces. Functions for
 * two different purposes are gathered here: the abstract integration
 * on finite element meshes in MeshWorker and actual implementations
 * of the integration of cell and face terms for concrete problems in
 * LocalIntegrators.
 *
 * @note Documentation on coding conventions, relations between
 * classes, and details of the implementation is found in the
 * documentation of namespaces in this module.
 *
 * <h3>Integration on finite element meshes</h3>
 *
 * When we integrate a function or a functional on a finite element
 * space, the structure of the integration loop is always the same. We
 * have between 3 and 5 nested loops, from outside to inside:
 * <ol>
 * <li> Loop over all cells
 * <li> Optionally, loop over all faces to compute fluxes
 * <li> Loop over all quadrature points of the cell/face
 * <li> Optionally, loop over all test functions to compute forms
 * <li> Optionally, loop over all trial functions to compute bilinear
 * forms
 * </ol>
 *
 * These loops naturally fall into two classes, namely the computation
 * of cell and face contributions (loops 3 to 5), and the outer loops
 * over the mesh objects, often referred to as <em>assembling</em>.
 *
 * Support for the outer loop in deal.II can be found in the namespace
 * MeshWorker (see the documentation there). In order to support the
 * cell and face contributions (referred to as local contributions
 * from now on), deal.II offers FEValuesBase and its derived
 * classes. While the outer loop is generic (with exception of the
 * data types), the computation of local contributions is problem
 * dependent. Therefore, no generic algorithm is possible
 * here. Nevertheless, we can define a generic interface for functions
 * for this purpose and provide a library of local integrators for use
 * in applications. These are collected in the namespace
 * LocalIntegrators
 */
