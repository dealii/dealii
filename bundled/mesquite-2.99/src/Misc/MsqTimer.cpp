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
#include "MsqTimer.hpp"

#include <iostream>
#include <iomanip>

// Create the global collection of stop watches
Mesquite::StopWatchCollection Mesquite::GlobalStopWatches;

#ifdef HAVE_TIMES
#  include <sys/times.h>
#  include <unistd.h>
#  include <limits.h>
#  ifndef CLK_TCK
#    ifdef _SC_CLK_TCK
#      define CLK_TCK sysconf(_SC_CLK_TCK)
#    else
#      include <sys/param.h>
#      ifdef HZ
#        define CLK_TCK HZ
#      else
#        error times(3) w/out CLK_TCK.  Please report this.
#        undef HAVE_TIMES
#      endif
#    endif
#  endif
#endif

#ifdef HAVE_TIMES
   static inline double now()
   {
     tms t;
     times( &t );
     return (double)(t.tms_utime + t.tms_stime) / CLK_TCK;
   }
#elif defined(HAVE_CLOCK)
#  include <ctime>
   static inline double now()
   {
     return (double)std::clock() / CLOCKS_PER_SEC; 
   }
#endif



Mesquite::Timer::Timer() 
    : atBirth(now())
{
  atLastCheck = atBirth;
}

void Mesquite::Timer::reset()
{
  atBirth=now();
  atLastCheck = atBirth;
}

double Mesquite::Timer::since_last_check()
{
  double right_now = now();
  double rv = right_now - atLastCheck;
  atLastCheck = right_now;
  return rv;
}

double Mesquite::Timer::since_birth() const
{
  return now() - atBirth;
}

void Mesquite::StopWatch::start()
{
  if (!isRunning)
  {
    isRunning = true;
    timeAtLastStart=now();
    ++numStarts;
  }
}

void Mesquite::StopWatch::stop()
{
  if (isRunning)
  {
    isRunning = false;
    totalTime += now() - timeAtLastStart;
  }
}

void Mesquite::StopWatch::reset()
{
  isRunning=false;
  totalTime=0;
  numStarts=0;
}

double Mesquite::StopWatch::total_time() const
{
  double rv = totalTime;
  if (isRunning)
    rv += now() - timeAtLastStart;
  return rv;
}

Mesquite::StopWatchCollection::Key Mesquite::StopWatchCollection::add(
  const std::string &name,
  bool fail_if_exists)
{
    // Don't allow empty name
  if (name == "")
    return 0;
  
  Key key = get_key(name);
  
    // If the named stopwatch doesn't exist...
  if (!key)
  {
      // See if there is an unused existing stopwatch
    size_t i;
    for (i = 0; i < mEntries.size(); i++)
    {
      if (mEntries[i].first == "")
      {
        mEntries[i].first = name;
        mEntries[i].second.reset();
        break;
      }
    }
      // If not, create a new one
    if (i == mEntries.size())
    {
      mEntries.push_back(std::pair<std::string, StopWatch>(name, StopWatch()));
    }
    key = i+1;
  }
    // If it already existed...
  else if (fail_if_exists)
    key = 0;
  
  return key;
}


Mesquite::StopWatchCollection::Key Mesquite::StopWatchCollection::get_key(
  const std::string &name) const
{
  Key key = 0;
  
  for (size_t i = 0; i < mEntries.size(); i++)
  {
    if (mEntries[i].first == name)
    {
      key = i + 1;
      break;
    }
  }

  return key;
}

void Mesquite::StopWatchCollection::remove(
  const Mesquite::StopWatchCollection::Key key)
{
    // Get rid of anything at the end of the list
  if (key == mEntries.size())
  {
    mEntries.pop_back();
    while (!mEntries.empty() && mEntries.back().first == "")
    {
      mEntries.pop_back();
    }
  }
  
  else if (key > 0 && key < mEntries.size())
  {
      // If in the middle of the list, set its name to ""
    mEntries[key-1].first = "";
  }
}


void Mesquite::StopWatchCollection::start(
  const Mesquite::StopWatchCollection::Key key)
{
  if (key > 0 &&
      key <= mEntries.size() &&
      mEntries[key-1].first != "")
    mEntries[key-1].second.start();
}

void Mesquite::StopWatchCollection::stop(
  const Mesquite::StopWatchCollection::Key key)
{
  if (key > 0 &&
      key <= mEntries.size() &&
      mEntries[key-1].first != "")
    mEntries[key-1].second.stop();
}

void Mesquite::StopWatchCollection::reset(
  const Mesquite::StopWatchCollection::Key key)
{
  if (key > 0 &&
      key <= mEntries.size())
    mEntries[key-1].second.reset();
}


double Mesquite::StopWatchCollection::total_time(
  const Mesquite::StopWatchCollection::Key key) const
{
  if (key > 0 &&
      key <= mEntries.size() &&
      mEntries[key-1].first != "")
    return mEntries[key-1].second.total_time();
  else
    return 0.0;
}

int Mesquite::StopWatchCollection::number_of_starts(
  const Mesquite::StopWatchCollection::Key key) const
{
  if (key > 0 &&
      key <= mEntries.size() &&
      mEntries[key-1].first != "")
    return mEntries[key-1].second.number_of_starts();
  else
    return 0;
}
/*! Fills a vector of StopWatchCollection::Key in which the Keys are ordered
  by the associated StopWatch's total_time.  The key associated with the
  largest total_time StopWatch is in the first position of the vector.  The
  key associated with the smallest total_time StopWatch is in the last
  position of the vector.*/
void Mesquite::StopWatchCollection::get_keys_sorted_by_time(
  std::vector<Key> &sorted_keys)
{
  int num_watches=mEntries.size();
  int *sorted_indices=new int[num_watches];
  int i=0;
  int counter=0;
  for(i=0;i<num_watches;++i){
    sorted_indices[i]=0;
  }
  double current_max;
  int index_to_max;
    //While we haven't added all of the Keys to the vector
  while(counter<num_watches){
    current_max=-1;
    index_to_max=-1;
      //loop over the times and find the largest remaining
    for(i=0;i<num_watches;++i){
      if(mEntries[i].second.total_time()>current_max && sorted_indices[i]==0){
        current_max=mEntries[i].second.total_time();
        index_to_max=i;
      }
    }
      //Add the key associated with index_to_max and any subsequent
      //keys which are associated with a StopWatch that has a total
      //time equal to current_max;
    for(i=index_to_max;i<num_watches;++i){
      if(mEntries[i].second.total_time()>=current_max && sorted_indices[i]==0)
      {
        counter++;
        sorted_indices[i]=counter;
        sorted_keys.push_back(i+1);
      }
    }
  }
    //clean up
  delete[] sorted_indices;
}

// Originally in MsqMessage.cpp
// Moved here and converted to an ostream operator 
// by J.Kraftcheck, 2004-10-18
std::ostream& Mesquite::operator<<( std::ostream& str,
                                          Mesquite::StopWatchCollection& coll )
{
  std::vector<Mesquite::StopWatchCollection::Key> sorted_keys;
  Mesquite::GlobalStopWatches.get_keys_sorted_by_time(sorted_keys);
  int number_of_keys=sorted_keys.size();
  int i =0;
  str<<"\nTIME        | NUM. STARTS | TIMER NAME ("<<number_of_keys<<" timers)\n";
  for(i=0;i<number_of_keys;++i){
	  str<<std::setiosflags(std::ios::left)
             <<std::setw(13)
             <<Mesquite::GlobalStopWatches.total_time(sorted_keys[i])
             <<" "
             <<std::setw(13)
             <<Mesquite::GlobalStopWatches.number_of_starts(sorted_keys[i])
             <<" "
             <<Mesquite::GlobalStopWatches.get_string(sorted_keys[i])
             <<std::endl;
  }
  return str;
}


