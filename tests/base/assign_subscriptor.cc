// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check the behavior of the ObserverPointer-EnableObserverPointer
// pair for copy and move semantics.


#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/observer_pointer.h>

#include <iostream>
#include <vector>

#include "../tests.h"

class Test : public EnableObserverPointer
{};

int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();

  {
    deallog << "Checking copy assignment" << std::endl;

    EnableObserverPointer                  subscriptor_1;
    EnableObserverPointer                  subscriptor_2;
    ObserverPointer<EnableObserverPointer> smart_pointer_1(&subscriptor_1);
    ObserverPointer<EnableObserverPointer> smart_pointer_2(&subscriptor_2);

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

    EnableObserverPointer                  subscriptor_1;
    ObserverPointer<EnableObserverPointer> smart_pointer_1(&subscriptor_1);

    EnableObserverPointer subscriptor_2(subscriptor_1);

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

    EnableObserverPointer                  subscriptor_1;
    EnableObserverPointer                  subscriptor_2;
    ObserverPointer<EnableObserverPointer> smart_pointer_1(&subscriptor_1);
    ObserverPointer<EnableObserverPointer> smart_pointer_2(&subscriptor_2);

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

    EnableObserverPointer                  subscriptor_1;
    ObserverPointer<EnableObserverPointer> smart_pointer_1(&subscriptor_1);

    EnableObserverPointer subscriptor_2(std::move(subscriptor_1));

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
