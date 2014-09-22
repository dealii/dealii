// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2013 by the deal.II authors
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
 * @defgroup memory Memory handling
 *
 * This group has some basic classes and namespaces for memory
 * handling. The Subscriptor and SmartPointer classes are used for
 * counted memory handling, i.e. whenever a SmartPointer is set to
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
 * objects. For some simple classes, like the STL containers, it
 * directly determines how much memory they need (or at least gives an
 * estimate). For deal.II classes, it uses the
 * <code>memory_consumption</code> member function that most classes
 * have.
 *
 * @ingroup utilities
 */
