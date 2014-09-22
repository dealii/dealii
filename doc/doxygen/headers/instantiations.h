// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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
 * @page Instantiations Template instantiations
 *
 * Instantiation of complex class and function templates is expensive both in
 * terms of compile time and disk space. Therefore, we try to separate
 * declaration and implementation of templates as far as possible, and make
 * sure that implementations are read by the compiler only when 
 * necessary.
 *
 * Template classes in <tt>deal.II</tt> can be grouped into three categories,
 * depending on the number of probable different instantiations. These three
 * groups are discussed in the following.
 * 
 *
 * @section Inst1 Known and fixed number of instantiations
 *
 * These are the classes having template parameters with very
 * predictable values. The typical prototype is
 * @code
 * template <int dim> class Function;
 * @endcode
 *
 * Here, we have a small number of instantiations (<code>dim = 1,2,3</code>)
 * known at the time of design of the library. Therefore, member functions of
 * this class are defined in a <tt>.cc</tt> file in the source directory and
 * we instantiate the template for these known values explicitly in the source
 * file.
 *
 * From an application viewpoint, all you actually get to see then is the
 * declaration of the template. Actual instantiations of member functions
 * happens inside the library and is done when you compile the library, not
 * when you compile your application code.
 *
 * For these classes, adding instantiations for new parameters involves
 * changing the library. However, this is rarely needed, of course, unless you
 * are not content with computing only in 1d, 2d, or 3d.
 * 
 *
 * @subsection Inst1a Available instances
 *
 * If the template parameter is <tt>dim</tt>, the available instances
 * are for <tt>dim=1,2,3</tt>, if there is no other information.
 *
 * There are other cases of classes (not depending on the spatial
 * dimension) for which only a certain, small number of template
 * arguments is supported and explicit instantiations are provided in
 * the library. In particular, this includes all the linear algebra
 * classes that are templatized on the type of the scalar underlying
 * stored values: we only support <code>double</code>,
 * <code>float</code>, and in some cases <code>long double</code>,
 * <code>std::complex@<double@></code>,
 * <code>std::complex@<float@></code>, and <code>std::complex@<long
 * double@></code>.
 * 
 *
 * @section Inst2 A few instantiations, most of which are known
 *
 * These are class templates usually having a small number of instantiations,
 * but additional instantiations may be necessary. Therefore, a set of
 * instantiations for the most likely parameters is provided precompiled in
 * the libraries, but the implementation of the templates are provided in a
 * special header file so that it is accessible in case someone wants to
 * instantiate it for an unforeseen argument.
 *
 * Typical examples for this would be some of the linear algebra classes that
 * take a vector type as template argument. They would be instantiated within
 * the library for <code>Vector&lt;double&gt;</code>,
 * <code>Vector&lt;float&gt;</code>, <code>BlockVector&lt;double&gt;</code>,
 * and <code>BlockVector&lt;float&gt;</code>, for example. However, they may
 * also be used with other vector types as long as they satisfy certain
 * interfaces, including vector types that are not part of the library but
 * possibly defined in an application program. In such a case, applications
 * can instantiate these templates by hand as described in the next section.
 * 
 *
 * @subsection Inst2c Creating new instances
 *
 * Choose one of your source files to provide the required
 * instantiations. Say that you want the class template <tt>XXXX</tt>,
 * defined in the header file <tt>xxxx.h</tt>, instantiated with the
 * template parameter <tt>Lager</tt>. Then, your file should contain
 * the lines
 * @code
 *                   // Include class template declaration
 * #include <xxxx.h>
 *                   // Include member definitions
 * #include <xxxx.templates.h>
 *
 * ...
 *
 * template class XXXX<Lager>;
 * @endcode
 *
 * 
 * @subsection Inst2p Provided instances
 *
 * Like with the classes in section @ref Inst1, the instances provided in the
 * library are often listed in the documentation of that class in a form
 * similar to this:
 @verbatim
 Template Instantiations: some  (<p1>a,b,c<p2>)
 @endverbatim
 *
 *
 * @section Inst3 Many unknown instantiations
 *
 * These are the classes, where no reasonable predetermined set of instances
 * exists. Therefore, all member definitions are included in the header file
 * and are instantiated wherever needed.  An example would be the SmartPointer
 * class template that can be used with virtually any template argument.
 */
