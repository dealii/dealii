/**
 * @page Instantiations Template instantiations
 *
 * The instantiation of complex function templates costs a lot of
 * compile time and blows up the size of the resulting object
 * files. Therefore, declaration and implementation of these functions
 * should be split and the implementation should be read by the
 * compiler only when really necessary. We achieve compile time
 * efficiency by implementing the following concept:
 *
 * Template classes in <tt>deal.II</tt> are split into three
 * categories, depending on the number of probable different
 * instantiations. Accordingly, the function definitions are in
 * different places.
 *
 * @section Inst1 Few instantiations
 *
 * These are the classes having template parameters with very
 * predictable values. The typical prototype is
 * @code
 * template <int dim> class Function;
 * @endcode
 *
 * Here, we have a small number of instantiations (dim = 1,2,3)
 * known at the time of design of the library. Therefore, member
 * functions of this class are defined in a <tt>.cc</tt> file in the
 * source directory and all instantiations are in the source file.
 *
 * This means that it is in fact not the template being part of the
 * library, but its instantiation; in the case of triangulation
 * objects of the <tt>deal.II</tt> library, implementations may be
 * completely different indeed. There, we even put instantiations for
 * different template parameters into different libraries.
 *
 * For these classes, adding instantiations for new parameters
 * involves changing the library.
 *
 * @subsection Inst1a Available instances
 *
 * If the template parameter is <tt>dim</tt>, the available instances
 * are for <tt>dim=1,2,3</tt>, if there is no other information.
 *
 * For other classes, the list of available instances is often listed
 * in the line
 @verbatim
 Template Instantiations: few  (<p1> <p2>)
 @endverbatim
 *
 * @section Inst2 Some instantiations
 *
 * These are class templates usually having a small number of
 * instantiations, but additional instantiations may be
 * necessary. Therefore, a set of instantiations for standard
 * parameters is provided precompiled in the libraries.
 *
 * Generic example of such a class is
 * @code
 * template<typename number> class Polynomials::Polynomial<number>;
 * @endcode
 *
 * Here, the instantiations provided are for <tt>number</tt> types
 * <tt>float</tt>, <tt>double</tt> and <tt>long double</tt>. Still,
 * there are more number types and somebody might want to use complex
 * polynomials. Therefore, complementing the declaration of
 * Polynomials::Polynomial in <tt>polynomial.h</tt>, there is a file
 * <tt>polynomial.templates.h</tt>, providing the implementation
 * source code.
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
 * @subsection Inst2p Provided instances
 *
 * Like with the classes in section Inst1, the instances provided in
 * the library are often listed in the line pointing to this page,
 * that is for example:
 @verbatim
 Template Instantiations: some  (<p1> <p2>)
 @endverbatim
 *
 * @section Inst3 Many instantiations
 *
 * These are the classes, where no reasonable predetermined set of
 * instances exists. Therefore, all member definitions are included
 * into the header file and are instantiated wherever needed.
 * For an example, see
 * @code
 * template<typename T> class SmartPointer;
 * @endcode
 */
