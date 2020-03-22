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

#include "MsqInterrupt.hpp"
#include <signal.h>

namespace MESQUITE_NS {

unsigned MsqInterrupt::instanceCount = 0;
bool MsqInterrupt::sawInterrupt = false;
MsqInterrupt::InterruptMode MsqInterrupt::interruptMode = MsqInterrupt::AUTO;

extern "C" { typedef void (*msq_sig_handler_t)(int); }
msq_sig_handler_t oldHandler = SIG_ERR;


extern "C" void msq_sigint_handler( int ) 
{  
  MsqInterrupt::set_interrupt(); 
  if (oldHandler != SIG_DFL && oldHandler != SIG_IGN)
    oldHandler( SIGINT );
  MsqInterrupt::set_handler();
}
   

void MsqInterrupt::set_handler()
{
  oldHandler = signal( SIGINT, &msq_sigint_handler );
  if (MsqInterrupt::interruptMode == MsqInterrupt::AUTO &&
      (oldHandler == SIG_DFL || oldHandler == SIG_IGN))
  {
    signal( SIGINT, oldHandler );
    oldHandler = SIG_ERR;
  }  
}

void MsqInterrupt::disable( MsqError& /*err*/ )
{
  interruptMode = IGNORE;
  if (instanceCount && SIG_ERR != oldHandler)
  {
    signal( SIGINT, oldHandler );
    oldHandler = SIG_ERR;
  }
  sawInterrupt = false;
}

void MsqInterrupt::enable( MsqError& /*err*/ )
{
  sawInterrupt = false;
  interruptMode = CATCH;
  if (instanceCount && SIG_ERR == oldHandler)
    set_handler();
}

void MsqInterrupt::allow( MsqError& /*err*/ )
{
  sawInterrupt = false;
  interruptMode = AUTO;
  if (instanceCount && SIG_ERR == oldHandler)
    set_handler();
}

MsqInterrupt::MsqInterrupt()
{
  if (!instanceCount)
  {
    if (IGNORE != interruptMode)
      set_handler();
    sawInterrupt = false;
  }  
  ++instanceCount;
}

MsqInterrupt::~MsqInterrupt()
{
  if (!--instanceCount && SIG_ERR != oldHandler)
  {
    signal( SIGINT, oldHandler );
    oldHandler = SIG_ERR;
  }
  sawInterrupt = false;
}

} // namespace Mesquite


