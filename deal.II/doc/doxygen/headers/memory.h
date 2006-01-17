//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

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
 */
