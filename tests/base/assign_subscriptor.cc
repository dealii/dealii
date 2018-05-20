// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

// check that Subscriptor objects need to be empty before assigning.

#include "../tests.h"
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>
#include <iostream>
#include <vector>

class Test : public Subscriptor
{};

int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();

  // should work
  {
    Subscriptor subscriptor_1;
    Subscriptor subscriptor_2;

    SmartPointer<Subscriptor> smart(&subscriptor_1);

    subscriptor_2 = subscriptor_1;
  }

  try
    {
      Subscriptor subscriptor_1;
      Subscriptor subscriptor_2;

      SmartPointer<Subscriptor> smart(&subscriptor_2);

      subscriptor_2 = subscriptor_1;
    }
  catch(ExceptionBase& e)
    {
      deallog << e.get_exc_name() << std::endl;
    }

  try
    {
      Subscriptor subscriptor_1;
      Subscriptor subscriptor_2;

      SmartPointer<Subscriptor> smart(&subscriptor_1);

      subscriptor_2 = std::move(subscriptor_1);
    }
  catch(ExceptionBase& e)
    {
      deallog << e.get_exc_name() << std::endl;
    }

  try
    {
      Subscriptor subscriptor_1;
      Subscriptor subscriptor_2;

      SmartPointer<Subscriptor> smart(&subscriptor_2);

      subscriptor_2 = std::move(subscriptor_1);
    }
  catch(ExceptionBase& e)
    {
      deallog << e.get_exc_name() << std::endl;
    }
}
