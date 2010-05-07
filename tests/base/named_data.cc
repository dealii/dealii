//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------


#include "../tests.h"
#include <base/named_data.h>
#include <base/logstream.h>

#include <fstream>
#include <boost/shared_ptr.hpp>

template<typename DATA>
void
test_const(const NamedData<DATA>& data)
{
  NamedData<DATA> newdata;
  newdata.merge(data);
  data.print(deallog);
  
  deallog << "Const should be true: " << (newdata.is_const() ? "true" : "false") << std::endl;
}

template<typename DATA>
void
test_selector(const NamedData<DATA>& data)
{
  NamedSelection select;
  select.add("P1");
  select.add("P6");
  select.add("P5");  
  select.add("P3");
  select.initialize(data);

  for (unsigned int i=0;i<select.size();++i)
    {
      deallog << "Selection " << i << ':';
      if (select(i) != deal_II_numbers::invalid_unsigned_int)
	deallog << data.name(select(i)) << ' ' << *data(select(i));
      deallog << std::endl;
    }
}

void
test_shared_pointer()
{
  NamedData<std_cxx1x::shared_ptr<double> > data;
  deallog << "const" << data.is_const() << std::endl;
  std_cxx1x::shared_ptr<double> p = std_cxx1x::shared_ptr<double>(new double);
  data.add(p, "P1");
  p = std_cxx1x::shared_ptr<double>(new double);
  data.add(p, "P2");
  p = std_cxx1x::shared_ptr<double>(new double);
  data.add(p, "P3");
  p = std_cxx1x::shared_ptr<double>(new double);
  data.add(p, "P4");
  p = std_cxx1x::shared_ptr<double>(new double);
  data.add(p, "P5");

  for (unsigned int i=0;i<data.size();++i)
    *data(i) = i+0.5;

  test_const(data);
  test_selector(data);
}


int main ()
{
  std::ofstream logfile("named_data/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
//  deallog.log_cerr();

  test_shared_pointer();
}
