//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__conditional_ostream_h
#define __deal2__conditional_ostream_h



#ifdef HAVE_STD_OSTREAM_HEADER
#  include <ostream>
#else
#  include <iostream>
#endif


/**
 * A class that allows printing to an output stream, e.g. @p std::cout,
 * depending on the ConditionalOStream object being active (default)
 * or not. The condition of this object can be changed by
 * set_condition() and in the constructor.
 *
 * This class is mostly useful in parallel computations. Ordinarily, you would
 * use @p std::cout to print messages like what the program is presently
 * doing, or the number of degrees of freedom in each step. However, in
 * parallel programs, this means that each of the MPI processes write to the
 * screen, which yields many repetitions of the same text. To avoid it, one
 * would have to have a designated process, say the one with MPI process
 * number zero, do the output, and guard each write statement with an
 * if-condition. This becomes cumbersome and clutters up the code. Rather than
 * doing so, the present class can be used: objects of its type act just like
 * a standard output stream, but they only print something based on a
 * condition that can be set to, for example, <tt>mpi_process==0</tt>, so that
 * only one process has a true condition and in all other processes writes to
 * this object just disappear in nirvana.
 *
 * The usual usage of this class is as follows:
 *
 * @code
 * ConditionalOStream pout(std::cout, this_mpi_process==0);
 *
 *                                  // all processes print following
 *                                  // information to standard output
 * std::cout << "Reading parameter file on process "
 *           << this_mpi_process << std::endl;
 *
 *                                  // following is printed by
 *                                  // process 0 only
 * pout << "Solving ..." << std::endl;
 * solve();
 * pout << "done" << std::endl;
 * @endcode
 *
 * Here, `Reading parameter file on process xy' is printed by each
 * process separately. In contrast to that, `Solving ...' and `done'
 * is printed to standard output only once, namely by process 0.
 *
 * This class is not derived from ostream. Therefore
 * @code
 * system_matrix.print_formatted(pout);
 * @endcode
 * is <em>not</em> possible. Instead use the is_active() funtion for a
 * work-around:
 *
 * @code
 * if (pout.is_active())
 *   system_matrix.print_formatted(cout);
 * @endcode
 *
 * @author Ralf Hartmann, Wolfgang Bangerth, 2004
 */
class ConditionalOStream
{
  public:
				     /**
				      * Constructor. Set the stream to which
				      * we want to write, and the condition
				      * based on which writes are actually
				      * forwarded. Per default the condition
				      * of an object is active.
				      */
    ConditionalOStream (std::ostream &stream,
                        const bool    active = true);

				     /**
				      * Depending on the
				      * <tt>active</tt> flag set the
				      * condition of this stream to
				      * active (true) or non-active
				      * (false). An object of this
				      * class prints to <tt>cout</tt>
				      * if and only if its condition
				      * is active.
				      */
    void set_condition (const bool active);

				     /**
				      * Return the condition of the object.
				      */
    bool is_active() const;

    				     /**
				      * Output a constant something through
				      * this stream. This function must be @p
				      * const so that member objects of this
				      * type can also be used from @p const
				      * member functions of the surrounding
				      * class.
				      */
    template <typename T>
    const ConditionalOStream &
    operator << (const T &t) const;

				     /**
				      * Treat ostream manipulators. This
				      * function must be @p const so that
				      * member objects of this type can also
				      * be used from @p const member functions
				      * of the surrounding class.
				      *
				      * Note that compilers want to see this
				      * treated differently from the general
				      * template above since functions like @p
				      * std::endl are actually overloaded and
				      * can't be bound directly to a template
				      * type.
				      */
    const ConditionalOStream &
    operator<< (std::ostream& (*p) (std::ostream&)) const;

  private:
				     /**
				      * Pointer to <tt>cout</tt>. This
				      * class could easily be extended
				      * to treat streams different
				      * to the standard output.
				      *
				      * This variable must be @p mutable so
				      * that we can write to it in above @p
				      * const @p operator<< functions. For the
				      * reason why they, in turn, need to be
				      * @p const, see there.
				      */
    mutable std::ostream  &output_stream;

				     /**
				      * Stores the actual condition
				      * the object is in.
				      */
    bool active_flag;
};


// --------------------------- inline and template functions -----------

template <class T>
inline
const ConditionalOStream &
ConditionalOStream::operator<< (const T& t) const
{
  if (active_flag == true)
    output_stream << t;

  return *this;
}


inline
const ConditionalOStream &
ConditionalOStream::operator<< (std::ostream& (*p) (std::ostream&)) const
{
  if (active_flag == true)
    output_stream << p;

  return *this;
}


#endif
