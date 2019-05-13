// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



// Check the behavior of the SmartPointer-Subscriptor pair
// for copy and move semantics.


#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>

#include <iostream>
#include <vector>

#include "../tests.h"

class Test : public Subscriptor
{};

int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();

  {
    deallog << "Checking copy assignment" << std::endl;

    Subscriptor               subscriptor_1;
    Subscriptor               subscriptor_2;
    SmartPointer<Subscriptor> smart_pointer_1(&subscriptor_1);
    SmartPointer<Subscriptor> smart_pointer_2(&subscriptor_2);

    subscriptor_2 = subscriptor_1;

    deallog << "Checking smart_pointer_1" << std::endl;
    try
      {
        const auto dummy_1 = *smart_pointer_1;
        (void)dummy_1;
      }
    catch (ExceptionBase &e)
      {
        deallog << e.get_exc_name() << std::endl;
      }

    deallog << "Checking smart_pointer_2" << std::endl;
    try
      {
        const auto dummy_2 = *smart_pointer_2;
        (void)dummy_2;
      }
    catch (ExceptionBase &e)
      {
        deallog << e.get_exc_name() << std::endl;
      }
    deallog << std::endl;
  }

  {
    deallog << "Checking copy construction" << std::endl;

    Subscriptor               subscriptor_1;
    SmartPointer<Subscriptor> smart_pointer_1(&subscriptor_1);

    Subscriptor subscriptor_2(subscriptor_1);

    deallog << "Checking smart_pointer_1" << std::endl;
    try
      {
        const auto dummy_1 = *smart_pointer_1;
        (void)dummy_1;
      }
    catch (ExceptionBase &e)
      {
        deallog << e.get_exc_name() << std::endl;
      }
    deallog << std::endl;
  }

  {
    deallog << "Checking move assignment" << std::endl;

    Subscriptor               subscriptor_1;
    Subscriptor               subscriptor_2;
    SmartPointer<Subscriptor> smart_pointer_1(&subscriptor_1);
    SmartPointer<Subscriptor> smart_pointer_2(&subscriptor_2);

    subscriptor_2 = std::move(subscriptor_1);

    deallog << "Checking smart_pointer_1" << std::endl;
    try
      {
        const auto dummy_1 = *smart_pointer_1;
        (void)dummy_1;
      }
    catch (ExceptionBase &e)
      {
        deallog << e.get_exc_name() << std::endl;
      }

    deallog << "Checking smart_pointer_2" << std::endl;
    try
      {
        const auto dummy_2 = *smart_pointer_2;
        (void)dummy_2;
      }
    catch (ExceptionBase &e)
      {
        deallog << e.get_exc_name() << std::endl;
      }
    deallog << std::endl;
  }

  {
    deallog << "Checking move construction" << std::endl;

    Subscriptor               subscriptor_1;
    SmartPointer<Subscriptor> smart_pointer_1(&subscriptor_1);

    Subscriptor subscriptor_2(std::move(subscriptor_1));

    deallog << "Checking smart_pointer_1" << std::endl;
    try
      {
        const auto dummy_1 = *smart_pointer_1;
        (void)dummy_1;
      }
    catch (ExceptionBase &e)
      {
        deallog << e.get_exc_name() << std::endl;
      }
  }
}
