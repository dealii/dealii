// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test the core functionality of the GeneralDataStorage class


#include <deal.II/algorithms/general_data_storage.h>

#include <utility>

#include "../tests.h"

int
main()
{
  initlog();

  GeneralDataStorage data;

  deallog << "Add by copy" << std::endl;
  {
    data.reset();

    // Create new data instance
    const double val_1 = 1.0;
    data.add_unique_copy("value", val_1);
    const double &val_2 = data.get_object_with_name<double>("value");
    Assert(data.stores_object_with_name("value"), ExcInternalError());
    Assert(val_2 == val_1, ExcInternalError());
    Assert(&val_2 != &val_1, ExcInternalError());

    // Allowed overwrite of existing data
    const double val_3 = 2.0;
    data.add_or_overwrite_copy("value", val_3);
    const double &val_4 = data.get_object_with_name<double>("value");
    Assert(data.stores_object_with_name("value"), ExcInternalError());
    Assert(val_4 == val_3, ExcInternalError());
    Assert(&val_4 != &val_3, ExcInternalError());

    // Create new data instance using alternative function
    const double val_5 = 3.0;
    data.add_or_overwrite_copy("value_2", val_5);
    const double &val_6 = data.get_object_with_name<double>("value_2");
    Assert(data.stores_object_with_name("value_2"), ExcInternalError());
    Assert(val_6 == val_5, ExcInternalError());
    Assert(&val_6 != &val_5, ExcInternalError());

    deallog << "Size: " << data.size() << std::endl;
  }

  deallog << "Add by reference" << std::endl;
  {
    data.reset();

    // Create new data instance
    const double val_1 = 1.0;
    data.add_unique_reference("value", val_1);
    const double &val_2 = data.get_object_with_name<const double>("value");
    Assert(data.stores_object_with_name("value"), ExcInternalError());
    Assert(val_2 == val_1, ExcInternalError());
    Assert(&val_2 == &val_1, ExcInternalError());

    // Allowed overwrite of existing data
    const double val_3 = 2.0;
    data.add_or_overwrite_reference("value", val_3);
    const double &val_4 = data.get_object_with_name<const double>("value");
    Assert(data.stores_object_with_name("value"), ExcInternalError());
    Assert(val_4 == val_3, ExcInternalError());
    Assert(&val_4 == &val_3, ExcInternalError());

    // Create new data instance using alternative function
    const double val_5 = 3.0;
    data.add_or_overwrite_reference("value_2", val_5);
    const double &val_6 = data.get_object_with_name<const double>("value_2");
    Assert(data.stores_object_with_name("value_2"), ExcInternalError());
    Assert(val_6 == val_5, ExcInternalError());
    Assert(&val_6 == &val_5, ExcInternalError());

    deallog << "Size: " << data.size() << std::endl;
  }

  deallog << "Add or construct" << std::endl;
  {
    data.reset();

    using Type = std::pair<double, double>;

    // Create new data instance
    const Type &val_1 =
      data.get_or_add_object_with_name<Type>("value", 1.0, 2.0);
    Assert(data.stores_object_with_name("value"), ExcInternalError());
    Assert(val_1 == Type({1.0, 2.0}), ExcInternalError());

    // Should not overwrite existing data
    const Type &val_2 =
      data.get_or_add_object_with_name<Type>("value", Type(3.0, 4.0));
    Assert(data.stores_object_with_name("value"), ExcInternalError());
    Assert(val_2 == Type({1.0, 2.0}), ExcInternalError());

    deallog << "Size: " << data.size() << std::endl;
  }

  deallog << "Merge" << std::endl;
  {
    data.reset();

    const double val_1 = 1.0;
    data.add_unique_copy("value", val_1);

    GeneralDataStorage data_2;
    data_2.add_unique_copy("value", 2.0); // Duplicate
    data_2.add_unique_copy("value_2", 3.0);

    deallog << "Data pre-merge:" << std::endl;
    data.print_info(deallog);
    deallog << "Size: " << data.size() << std::endl;
    deallog << "Data 2 pre-merge:" << std::endl;
    data_2.print_info(deallog);
    deallog << "Size: " << data_2.size() << std::endl;

    data.merge(data_2);

    deallog << "Data post-merge:" << std::endl;
    data.print_info(deallog);
    deallog << "Size: " << data.size() << std::endl;
  }

  deallog << "Ambiguous construction" << std::endl;
  {
    // Pass Arguments by lvalue reference
    {
      const double  val_1 = 1.0;
      const double &val_2 =
        data.get_or_add_object_with_name<double>("value", val_1);
      const std::string &str =
        data.get_or_add_object_with_name<std::string>("empty string");
    }

    // Pass Arguments by rvalue reference
    {
      double        val_1 = 1.0;
      const double &val_2 =
        data.get_or_add_object_with_name<double>("value", std::move(val_1));
      const std::string &str =
        data.get_or_add_object_with_name<std::string>("empty string");
    }

    // Pass Arguments ambiguously
    {
      const double &val_2 =
        data.get_or_add_object_with_name<double>("value", 1.0);
      const std::string &str =
        data.get_or_add_object_with_name<std::string>("empty string");
    }
  }


  deal_II_exceptions::disable_abort_on_exception();

  deallog << "Try to overwrite existing entry: Copy" << std::endl;
  {
    data.reset();

    const double val_1 = 1.0;
    data.add_unique_copy("value", val_1);

    try
      {
        const double val_2 = 1.0;
        data.add_unique_copy("value", val_2);
      }
    catch (const GeneralDataStorage::ExcNameHasBeenFound &exc)
      {
        deallog << exc.what() << std::endl;
      }
  }

  deallog << "Try to overwrite existing entry: Reference" << std::endl;
  {
    data.reset();

    const double val_1 = 1.0;
    data.add_unique_reference("value", val_1);

    try
      {
        const double val_2 = 2.0;
        data.add_unique_reference("value", val_2);
      }
    catch (const GeneralDataStorage::ExcNameHasBeenFound &exc)
      {
        deallog << exc.what() << std::endl;
      }
  }

  deallog << "Fetch non-existing entry" << std::endl;
  {
    data.reset();

    try
      {
        data.get_object_with_name<double>("value");
      }
    catch (const GeneralDataStorage::ExcNameNotFound &exc)
      {
        deallog << exc.what() << std::endl;
      }
  }

  deallog << "Access removed entry (reference)" << std::endl;
  {
    data.reset();

    const double val_1 = 1.0;
    data.add_unique_reference("value", val_1);
    data.remove_object_with_name("value");

    try
      {
        data.get_object_with_name<double>("value");
      }
    catch (const GeneralDataStorage::ExcNameNotFound &exc)
      {
        deallog << exc.what() << std::endl;
      }
  }

  deallog << "Access removed entry (copy)" << std::endl;
  {
    data.reset();

    data.add_unique_copy("value", 1.0);
    data.remove_object_with_name("value");

    try
      {
        data.get_object_with_name<double>("value");
      }
    catch (const GeneralDataStorage::ExcNameNotFound &exc)
      {
        deallog << exc.what() << std::endl;
      }
  }
}
