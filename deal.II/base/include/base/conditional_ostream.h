//-----------------------  conditional_ostream.h  ----------------------
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
//-----------------------  conditional_ostream.h  ----------------------
#ifndef __deal2__conditional_ostream_h
#define __deal2__conditional_ostream_h



#ifdef HAVE_STD_OSTREAM_HEADER
#  include <ostream>
#else
#  include <iostream>
#endif


/**
 * A class that allows printing to standard output, i.e. cout,
 * depending on the ConditionalOStream object being active (default)
 * or not. The condition of this object can be changed by
 * set_condition().
 *
 * The usual usage of this class is through the pregenerated object
 * <tt>pout</tt> which can be used within parallel computations as
 * follows:
 *
 * @code
 * pout.set_condition(this_mpi_process==0);
 *
 *                                  // all processes print following
 *                                  // information to standard output
 * cout << "Reading parameter file on process " << this_mpi_process << endl;
 *
 *                                  // following is printed by
 *                                  // process 0 only
 * pout << "Solving ..." << endl;
 * solve();
 * pout << "done" << endl;
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
 * @author Ralf Hartmann, 2004
 */
class ConditionalOStream
{
  public:
				     /**
				      * Constructor. Per default the
				      * condition of an object is
				      * active.
				      */
    ConditionalOStream();

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
    void set_condition(bool active);

				     /**
				      * Give read access to the
				      * condition of the object.
				      */
    bool is_active() const;

    				     /**
				      * Output a constant something
				      * through this stream.
				      */
    template <typename T>
    ConditionalOStream & operator << (const T &t);

				     /**
				      * Treat ostream manipulators.
				      */
    ConditionalOStream & operator<< (std::ostream& (*p) (std::ostream&));

  private:
				     /**
				      * Pointer to <tt>cout</tt>. This
				      * class could easily be extended
				      * to treat streams different
				      * to the standard output.
				      */
    std::ostream  *std_out;

				     /**
				      * Stores the actual condition
				      * the object is in.
				      */
    bool active_flag;
};


template <class T>
inline
ConditionalOStream &
ConditionalOStream::operator<< (const T& t)
{
  if (active_flag)
    *std_out << t;

  return *this;
}


inline
ConditionalOStream &
ConditionalOStream::operator<< (std::ostream& (*p) (std::ostream&))
{
  if (active_flag)
    *std_out << p;

  return *this;
}


/**
 * Pregenerated object for conditional output to <tt>cout</tt>. Prints
 * to standard output depending on <tt>pout</tt> being active (default)
 * or not. The output can be disabled by set_condition(false).
 *
 * This object is particularly useful in the context of parallel
 * computations in order to avoid multiple but identical output by
 * several processes.
 *
 * @code
 * pout.set_condition(this_mpi_process==0);
 *
 *                                  // all processes print following
 *                                  // information to standard output
 * cout << "Reading parameter file on process " << this_mpi_process << endl;
 *
 *                                  // following is printed by
 *                                  // process 0 only
 * pout << "Solving ..." << endl;
 * solve();
 * pout << "done" << endl;
 * @endcode
 *
 * Here, `Reading parameter file on process xy' is printed by each
 * process separately. In contrast to that, `Solving ...' and `done'
 * is printed to standard output only once, namely by process 0.
 *
 * @author Ralf Hartmann, 2004
 */
extern ConditionalOStream pout;

#endif
