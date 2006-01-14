//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

/**
 * @defgroup Exceptions Exceptions and assertions
 *
 * This module contains classes that are used in the exception mechanism of
 * deal.II. Exceptions are used in two different ways:
 *
 * <ul>
 * 
 *   <li> Static assertions: These are checks that are only enabled in debug
 *   mode, not in optimized (or production) mode. They are meant to check that
 *   parameters to functions satisfy certain properties and similar
 *   assertions. For example, static assertions are used to make sure that two
 *   vectors that are added together have the same number of components --
 *   everything else would not make any sense anyway.
 *
 *   Such checks are performed by the Assert macro in several thousand places
 *   within the library. Also, several tutorial programs starting with step-5
 *   show how to do this.
 *
 *   If a static assertion is violated, the exception mechanism generates an
 *   exception of a type that indicates what exactly goes wrong, displays
 *   appropriate information, and then aborts the program -- if you try to add
 *   two vectors of different length, there is nothing that can be done within
 *   the program to cope with the situation, you have to go fix the program
 *   code instead. The exceptions of this module are used to indicate the
 *   reason for the failure.
 *
 *
 *   <li> Dynamic assertions: These are used to check dynamic features, such
 *   as whether an output file can be written to. These are things that can't
 *   be checked statically, i.e. they may change from program run to program
 *   run. It is therefore insufficient to only check these situations in debug
 *   mode.
 *
 *   Rather, one has to check them every time during execution of a
 *   program. Within deal.II, this is done using the AssertThrow macro
 *   introduced in step-9, step-13, and following tutorial programs. The macro
 *   checks a condition, and if violated throws an exception of one of the
 *   types declared in this module, using the C++ <code>throw</code>
 *   mechanism. Since these are run-time exceptions, this gives the program
 *   the chance to catch the exception and, for example, write the output to a
 *   writable file instead.
 * </ul>
 *   
 * @author Wolfgang Bangerth, 1998-2006
 */
