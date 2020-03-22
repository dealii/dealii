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
   
  ***************************************************************** */

#include "MsqDebug.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

namespace MESQUITE_NS
{

/* IMPORTANT: Be careful to initialize the InitializeFlags object
              *LAST* as it will access the other static data. */

std::vector<std::ostream*> MsqDebug::streams;
std::vector<bool> MsqDebug::flags;

#ifdef MSQ_ENABLE_DEBUG
static void parse_flag_string( const char* str, bool flag_val )
{
  char* s = strdup(str);
  const char delim[] = ";,";
  for (const char* p = strtok( s, delim ); p; p = strtok( 0, delim ) ) {
    int beg, end;
    switch (sscanf( p, "%u - %u", &beg, &end )) {
      case 1:
        end = beg;
        // fall through 
      case 2:
        for ( ; beg <= end; ++beg)
          MsqDebug::set( beg, flag_val );
        break;
    }
  }
  free( s );
}

MsqDebug::InitializeFlags::InitializeFlags( )
{
  const unsigned flag_array[] = { MSQ_ENABLE_DEBUG, 0 };
  size_t length = sizeof(flag_array) / sizeof(unsigned) - 1;
  while (length > 0)
    MsqDebug::set( flag_array[--length], true );
  
  const char* envstr = getenv( "MESQUITE_DEBUG" );
  if (envstr) 
    parse_flag_string( envstr, true );
  envstr = getenv( "MESQUITE_NO_DEBUG" );
  if (envstr) 
    parse_flag_string( envstr, false );
}
#else
MsqDebug::InitializeFlags::InitializeFlags( ) {}
#endif

MsqDebug::InitializeFlags MsqDebug::init;


bool MsqDebug::get( unsigned flag )
{
  return flag < flags.size() && flags[flag];
}

void MsqDebug::set( unsigned flag, bool state )
{
  if (state)
  {
    if (flag >= flags.size())
    {
      flags.resize(flag+1);
    }
    flags[flag] = true;
  }
  else
  {
    if (flag < flags.size())
      flags[flag] = false;
  }
}

void MsqDebug::disable_all()
{
  flags.clear();
}

std::ostream& MsqDebug::get_stream( unsigned flag )
{
  if (flag < streams.size())
    return *streams[flag];
  else
    return std::cout;
}

void MsqDebug::set_stream( unsigned flag, std::ostream& stream )
{
  if (flag >= streams.size())
  {
    size_t old_size = streams.size();
    streams.resize( flag );
    for (unsigned i = old_size; i < flag; ++i)
      streams[i] = &std::cout;
  }
  streams[flag] = &stream;
}


void MsqDebug::FormatPrinter::print( const char* format, ... ) const
{
  if (!MsqDebug::get( myFlag ))
    return;
  
  char buffer[512];
  
#if defined(HAVE_VSNPRINTF)
  va_list args;
  va_start( args, format );
  vsnprintf( buffer, sizeof(buffer), format, args );
  va_end( args );
#elif defined(HAVE__VSNPRINTF)
  va_list args;
  va_start( args, format );
  _vsnprintf( buffer, sizeof(buffer), format, args );
  va_end( args );
#elif defined(HAVE_VSPRINTF) 
  va_list args;
  va_start( args, format );
  vsprintf( buffer, format, args );
  va_end( args );
#else
  strncpy( buffer, format, sizeof(buffer) );
  buffer[sizeof(buffer)-1] = '\0';
#endif

  MsqDebug::get_stream( myFlag ) << buffer;
}

}

