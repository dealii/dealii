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

#include <base/conditional_ostream.h>


ConditionalOStream::ConditionalOStream(std::ostream &stream,
                                       const bool    active)
                :
		output_stream (stream),
		active_flag(active)
{}


void ConditionalOStream::set_condition(bool flag)
{
  active_flag = flag;
}


bool ConditionalOStream::is_active() const
{
  return active_flag;
}
