// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2006 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/**
 * @defgroup memory Memory handling
 *
 * This group has some basic classes and namespaces for memory
 * handling. The EnableObserverPointer and ObserverPointer classes are used for
 * counted memory handling, i.e. whenever a ObserverPointer is set to
 * point to an object, it increases a counter in that object; when the
 * pointer is set to point elsewhere, it decreases it again. This way,
 * one always knows how many users of an object there still are. While
 * this is rarely useful in itself, it is used to generate an
 * exception if an object is destroyed while a pointer somewhere is
 * still pointing to it, as any access through that pointer at a later
 * time would otherwise lead to access of invalid memory regions.
 *
 * In contrast to this, the MemoryConsumption namespace provides
 * functions that can be used to determine the memory consumption of
 * objects. For some simple classes, like the standard library
 * containers, it directly determines how much memory they need (or at
 * least gives an estimate). For deal.II classes, it uses the
 * <code>memory_consumption</code> member function that most classes
 * have.
 *
 * @ingroup utilities
 */
