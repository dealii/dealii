/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
#ifndef MESQUITE_TIMER_HPP
#define MESQUITE_TIMER_HPP

#ifdef WIN32
#pragma warning ( 4 : 4786)
#endif

#include "Mesquite.hpp"
#include "MsqDebug.hpp"

#include <vector>
#include <utility>
#include <string>
#include <iosfwd>

namespace MESQUITE_NS
{
  class MESQUITE_EXPORT Timer
  {
  public:
    Timer();

    void reset();//resets the timer as if it were just created
    
    double since_last_check(); //- return time in seconds since last
                               //- call to since_last_check().  Note that
                               //- calling since_birth() doesn't count
                               //- as a check.  The first time this function
                               //- is called, it returns the time since birth.
    
    double since_birth() const;//- return time in seconds since
                               //- object was created.
    
  private:
    double atBirth;      //- Time at birth
    double atLastCheck;        //- Time at last call to since_last_check()
  };
  
  
  class MESQUITE_EXPORT StopWatch
  {
  public:
      // Creates the stopwatch.  The stopwatch is stopped
      // until start() is called.
    StopWatch() 
        :isRunning(false), totalTime(0.0), numStarts(0)
      {}
    
      // Starts the stopwatch.  If it was already running,
      // this function does nothing.
    void start();
    
      // Stops the stopwatch.  If it was not already running,
      // this function does nothing
    void stop();
    
      // Stops the stopwatch and resets the total_time() to zero.
    void reset();
    
      // Returns the total accumulated time.  If the stopwatch
      // is currently running, the time between the last start()
      // and the current time IS included in total_time().
    double total_time() const;    
      
      /*! \brief Returns the number of times this StopWatch has
        been started.*/
    int number_of_starts() const{
        return numStarts;
      }

    
  private:
    bool isRunning;
    double timeAtLastStart;
    double totalTime;
    int numStarts;
  };
  
  class StopWatchCollection
  {
  public:
    typedef size_t Key;
    
      // Create a new collection
    MESQUITE_EXPORT StopWatchCollection()
      {}
    
      // Add a stopwatch to the collection.  Returns a non-zero
      // StopWatchCollection::Key if it succeeds, zero if it fails.
      // If a StopWatch with the given name already exists in the
      // collection, the Key of the existing StopWatch is returned
      // if 'fail_if_exists' is false, or zero is returned if
      // 'fail_if_exists' is true.
    MESQUITE_EXPORT Key add(const std::string &name, bool fail_if_exists = true);
    
      // Gets the Key for an existing stopwatch.  If a stopwatch
      // with the given name does not exist, function returns zero.
    MESQUITE_EXPORT Key get_key(const std::string &name) const;

      //!Gets the string associated with a key
    MESQUITE_EXPORT std::string get_string(const Key key){
        return mEntries[key-1].first;}
      //!Gets the string associated with a key      
    MESQUITE_EXPORT void get_string(const Key key, std::string &new_string){
      new_string=mEntries[key-1].first;}
    
      // Remove a specific stopwatch.
    MESQUITE_EXPORT void remove(const Key key);
    MESQUITE_EXPORT void remove(const std::string &name)
      { remove(get_key(name)); }
    
      // start a specific stopwatch
    MESQUITE_EXPORT void start(const Key key);
    MESQUITE_EXPORT void start(const std::string &name)
      { start(get_key(name)); }
    
      // stop a specific stopwatch
    MESQUITE_EXPORT void stop(const Key key);
    MESQUITE_EXPORT void stop(const std::string &name)
      { stop(get_key(name)); }
    
      // reset a specific stopwatch
    MESQUITE_EXPORT void reset(const Key key);
    MESQUITE_EXPORT void reset(const std::string &name)
      { reset(get_key(name)); }
    
      // Get the total time for a specific stopwatch, zero if
      // the stopwatch doesn't exist.
    MESQUITE_EXPORT double total_time(const Key key) const;
    MESQUITE_EXPORT double total_time(const std::string &name) const
      { return total_time(get_key(name)); }
      // Get the number of times a StopWatch was started.
    MESQUITE_EXPORT int number_of_starts(const Key key) const;
    MESQUITE_EXPORT int number_of_starts(const std::string &name) const
      { return number_of_starts(get_key(name));}
    
      //Gets the number of stop watches in the collection
    MESQUITE_EXPORT int number_of_stop_watches(){
      return (int) mEntries.size();}

    MESQUITE_EXPORT void get_keys_sorted_by_time(std::vector<Key> &sorted_keys);
    
  private:
    std::vector< std::pair<std::string, StopWatch> > mEntries;
  };
  
  MESQUITE_EXPORT std::ostream& operator<<( std::ostream&, StopWatchCollection& coll );
  
    // A stopWatchCollection available anywhere
  extern Mesquite::StopWatchCollection GlobalStopWatches;

  MESQUITE_EXPORT inline void print_timing_diagnostics( int debugflag )
    { MSQ_DBGOUT(debugflag) << GlobalStopWatches; }

  MESQUITE_EXPORT inline void print_timing_diagnostics( std::ostream& stream )
    { stream << GlobalStopWatches; }


class FunctionTimer
{
  public:
    inline FunctionTimer( StopWatchCollection::Key key ) : mKey( key ) {}
    inline void start()      { GlobalStopWatches.start( mKey ); }
    inline ~FunctionTimer()  { GlobalStopWatches.stop( mKey ); }
  private:
    StopWatchCollection::Key mKey;
    // Don't allow any of this stuff (make them private)
    void* operator new(size_t size);
    FunctionTimer( const FunctionTimer& );
    FunctionTimer& operator=( const FunctionTimer& );
};

#ifdef MSQ_USE_FUNCTION_TIMERS
  #define MSQ_FUNCTION_TIMER( NAME ) \
    static Mesquite::StopWatchCollection::Key _mesquite_timer_key = \
      Mesquite::GlobalStopWatches.add( NAME, false );               \
    FunctionTimer _mesquite_timer( _mesquite_timer_key );           \
    _mesquite_timer.start()
#else
  #define MSQ_FUNCTION_TIMER( NAME )
#endif 

}  // namespace Mesquite

#endif
