// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_algorithms_general_data_storage_h
#define dealii_algorithms_general_data_storage_h

#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>

#include <boost/core/demangle.hpp>

#include <algorithm>
#include <any>
#include <map>
#include <ostream>
#include <string>
#include <typeinfo>

DEAL_II_NAMESPACE_OPEN


/**
 * This class facilitates the storage of any general data.
 * It offers a mechanism to store any amount of data, of any type,
 * which is then made accessible by an identifier string.
 *
 * When using this class, please cite @cite SartoriGiulianiBardelloni-2018-a.
 */
class GeneralDataStorage : public EnableObserverPointer
{
public:
  /**
   * Default constructor.
   */
  GeneralDataStorage() = default;

  /**
   * Copy constructor.
   */
  GeneralDataStorage(const GeneralDataStorage &) = default;

  /**
   * Move constructor.
   */
  GeneralDataStorage(GeneralDataStorage &&) = default;

  /**
   * Number of objects stored by this class instance.
   */
  std::size_t
  size() const;

  /**
   * Merge the contents of @p other_data with this object.
   */
  void
  merge(const GeneralDataStorage &other_data);

  /**
   * Print the contents of the internal cache to the @p stream.
   *
   * Each key and value pair in the @p any_data map are printed on an
   * individual line, with the <tt>std::string</tt> key listed first
   * followed by the demangled <tt>type_id</tt> of the associated mapped
   * type.
   */
  template <class Stream>
  void
  print_info(Stream &stream);

  /**
   * Clear all data stored in this class instance.
   *
   * When you call this function, it destroys all objects you asked to be stored
   * as copies, and it forgets about the references to data you asked to store
   * by reference. As a consequence, you are now free to destroy the objects to
   * which references were stored at whatever time you want -- before or after
   * the current `GeneralDataStorage` object is destroyed.
   *
   * To clarify this point, consider the following small example:
   *
   * @code
   *   GeneralDataStorage data;
   *
   *   {
   *     const double some_number = ...;
   *     data.add_unique_reference("value", some_number);
   *
   *     // Adding either of these next two lines could fix the
   *     // issue, by removing the association of some_number with data:
   *     // data.remove_object_with_name("value");
   *     // data.reset();
   *   } // some_number goes out of scope here
   *
   *   const double some_other_number
   *     = data.get_object_with_name<double>("value"); // Invalid call
   * @endcode
   *
   * In the code above, the @p data  object has a longer scope than
   * <tt>some_number</tt>. By the time we fetch the <tt>"value"</tt> from
   * @p data , the reference to @p some_number is no longer valid.
   *
   * Similarly, for data copied into a GeneralDataStorage object one should
   * consider the scope under which it remains valid:
   *
   * @code
   *   double* ptr_to_some_number = null_ptr;
   *
   *   {
   *     GeneralDataStorage data;
   *     const double some_number = ...;
   *     data.add_unique_copy("value", some_number);
   *
   *     ptr_to_some_number = &(data.get_object_with_name<double>("value"));
   *   } // The copy to some_number goes out of scope here
   *
   *   const double some_other_number
   *     = *ptr_to_some_number; // Invalid call
   * @endcode
   *
   * Similar to the first example, we must be conscious of the fact that the
   * copies of any @p Type stored by @p data only remains valid while the
   * GeneralDataStorage instance in which it is stored is alive.
   *
   * Furthermore, as elucidated in the last example, the copy of the
   * class instance (owned by GeneralDataStorage) that is being pointed to
   * is no longer alive when the reset() function is called, or when it is
   * removed via a call to remove_object_with_name().
   *
   * @code
   *   GeneralDataStorage data;
   *   double* ptr_to_some_number = null_ptr;
   *
   *   {
   *     const double some_number = ...;
   *     data.add_unique_copy("value", some_number);
   *
   *     ptr_to_some_number = &(data.get_object_with_name<double>("value"));
   *
   *     // The copy to some_number would go out of scope when either of
   *     // following two calls are made:
   *     data.remove_object_with_name("value");
   *     data.reset();
   *   }
   *
   *   const double some_other_number
   *     = *ptr_to_some_number; // Invalid call
   * @endcode
   */
  void
  reset();

  /**
   * @name Data storage and access
   * @{
   */

  /**
   * Store internally a copy of the given object. The copied object is
   * owned by this class, and is accessible via reference through the
   * get_object_with_name() method.
   *
   * This function ensures that no @p entry with the given @p name is
   * already stored by this class instance.
   */
  template <typename Type>
  void
  add_unique_copy(const std::string &name, const Type &entry);

  /**
   * Store internally a copy of the given object. The copied object is
   * owned by this class, and is accessible via reference through the
   * get_object_with_name() method.
   *
   * This function does not perform any checks to ensure that the @p entry
   * with the given @p name is already stored by this class instance. If the
   * @p name does in fact point to existing data, then this is overwritten.
   */
  template <typename Type>
  void
  add_or_overwrite_copy(const std::string &name, const Type &entry);

  /**
   * Add a reference to an already existing object. The object is not
   * owned by this class, and the user has to guarantee that the
   * referenced object lives longer than this class instance. The stored
   * reference is accessible through the get_object_with_name() method.
   *
   * This function ensures that no @p entry with the given @p name is
   * already stored by this class instance.
   */
  template <typename Type>
  void
  add_unique_reference(const std::string &name, Type &entry);

  /**
   * Add a reference to an already existing object. The object is not
   * owned by this class, and the user has to guarantee that the
   * referenced object lives longer than this class instance. The stored
   * reference is accessible through the get_object_with_name() method.
   *
   * This function does not perform any checks to ensure that the @p entry
   * with the given @p name is already stored by this class instance. If the
   * @p name does in fact point to existing data, then this is overwritten.
   */
  template <typename Type>
  void
  add_or_overwrite_reference(const std::string &name, Type &entry);

  /**
   * Return a reference to the object with given name. If the object does
   * not exist, then the input @p arguments will be used to construct an object
   * of the given @p Type and a reference to this new object then be returned.
   *
   * A copy of an object of type @p Type , which is owned by this class
   * instance, is generated by calling its constructor with the given set of
   * arguments. For this function, the @p arguments are passed as
   * <tt>lvalue</tt> references.
   */
  template <typename Type, typename Arg, typename... Args>
  Type &
  get_or_add_object_with_name(const std::string &name,
                              Arg               &argument,
                              Args &...arguments);

  /**
   * Return a reference to the object with given name. If the object does
   * not exist, then the input @p arguments will be used to construct an object
   * of the given @p Type and a reference to this new object then be returned.
   *
   * Same as above for a single argument.
   */
  template <typename Type, typename Arg>
  Type &
  get_or_add_object_with_name(const std::string &name, Arg &argument);

  /**
   * Return a reference to the object with given name. If the object does
   * not exist, then the input @p arguments will be used to construct an object
   * of the given @p Type and a reference to this new object then be returned.
   *
   * A copy of an object of type @p Type , which is owned by this class
   * instance, is generated by calling its constructor with the given set of
   * arguments. In contrast to the previous function of the same name, for
   * this function the @p arguments are passed as <tt>rvalue</tt> references.
   */
  template <typename Type, typename Arg, typename... Args>
  Type &
  get_or_add_object_with_name(const std::string &name,
                              Arg              &&argument,
                              Args &&...arguments);

  /**
   * Return a reference to the object with given name. If the object does
   * not exist, then the input @p arguments will be used to construct an object
   * of the given @p Type and a reference to this new object then be returned.
   *
   * Same as above for a single argument.
   */
  template <typename Type, typename Arg>
  Type &
  get_or_add_object_with_name(const std::string &name, Arg &&argument);

  /**
   * Same as above for default constructors.
   */
  template <typename Type>
  Type &
  get_or_add_object_with_name(const std::string &name);

  /**
   * Return a reference to the object with given name.
   *
   * This function throws an exception if either an object with the given name
   * is not stored in this class, or if the object with the given name is
   * neither of the exact specified @p Type nor a pointer to the @p Type.
   */
  template <typename Type>
  Type &
  get_object_with_name(const std::string &name);

  /**
   * Return a constant reference to the object with the given name.
   *
   * This function throws an exception if either an object with the given name
   * is not stored in this class, or if the object with the given name is
   * neither of the exact specified @p Type nor a pointer to the @p Type.
   */
  template <typename Type>
  const Type &
  get_object_with_name(const std::string &name) const;

  /**
   * Find out if we store an object with given name.
   */
  bool
  stores_object_with_name(const std::string &name) const;

  /**
   * Remove the object with given name.
   */
  void
  remove_object_with_name(const std::string &name);

  /** @} */

  /**
   * An entry with this name does not exist in the internal std::any map.
   */
  DeclException1(ExcNameNotFound,
                 std::string,
                 << "No entry with the name " << arg1 << " exists.");

  /**
   * An entry with this name does not exist in the internal std::any map.
   */
  DeclException1(ExcNameHasBeenFound,
                 std::string,
                 << "An entry with the name " << arg1 << " already exists.");

  /**
   * The requested type and the stored type are different.
   */
  DeclException3(ExcTypeMismatch,
                 std::string,
                 const char *,
                 const char *,
                 << "The stored type for entry with name \"" << arg1 << "\" is "
                 << arg2 << " but you requested type " << arg3 << '.');

private:
  /**
   * Arbitrary user data, identified by a string.
   */
  std::map<std::string, std::any> any_data;
};


/*----------------- Inline and template methods -----------------*/


#ifndef DOXYGEN


template <class Stream>
void
GeneralDataStorage::print_info(Stream &os)
{
  for (const auto &it : any_data)
    {
      os << it.first << '\t' << '\t'
         << boost::core::demangle(it.second.type().name()) << std::endl;
    }
}


template <typename Type>
void
GeneralDataStorage::add_unique_copy(const std::string &name, const Type &entry)
{
  AssertThrow(!stores_object_with_name(name), ExcNameHasBeenFound(name));
  add_or_overwrite_copy(name, entry);
}


template <typename Type>
void
GeneralDataStorage::add_or_overwrite_copy(const std::string &name,
                                          const Type        &entry)
{
  any_data[name] = entry;
}


template <typename Type>
void
GeneralDataStorage::add_unique_reference(const std::string &name, Type &entry)
{
  AssertThrow(!stores_object_with_name(name), ExcNameHasBeenFound(name));
  add_or_overwrite_reference(name, entry);
}


template <typename Type>
void
GeneralDataStorage::add_or_overwrite_reference(const std::string &name,
                                               Type              &entry)
{
  Type *ptr      = &entry;
  any_data[name] = ptr;
}


template <typename Type>
Type &
GeneralDataStorage::get_object_with_name(const std::string &name)
{
  AssertThrow(stores_object_with_name(name), ExcNameNotFound(name));

  Type *p = nullptr;

  if (any_data[name].type() == typeid(Type *))
    {
      p = std::any_cast<Type *>(any_data[name]);
    }
  else if (any_data[name].type() == typeid(Type))
    {
      p = std::any_cast<Type>(&any_data[name]);
    }
  else
    {
      AssertThrow(false,
                  ExcTypeMismatch(name,
                                  any_data[name].type().name(),
                                  typeid(Type).name()));
    }

  return *p;
}


template <typename Type>
const Type &
GeneralDataStorage::get_object_with_name(const std::string &name) const
{
  AssertThrow(stores_object_with_name(name), ExcNameNotFound(name));

  const auto it = any_data.find(name);

  if (it->second.type() == typeid(Type *))
    {
      const Type *p = std::any_cast<Type *>(it->second);
      return *p;
    }
  else if (it->second.type() == typeid(Type))
    {
      const Type *p = std::any_cast<Type>(&it->second);
      return *p;
    }
  else
    {
      AssertThrow(false,
                  ExcTypeMismatch(name,
                                  it->second.type().name(),
                                  typeid(Type).name()));
      const Type *p = nullptr;
      return *p;
    }
}



template <typename Type, typename Arg>
Type &
GeneralDataStorage::get_or_add_object_with_name(const std::string &name,
                                                Arg               &argument)
{
  if (!stores_object_with_name(name))
    add_unique_copy(name, Type(argument));

  return get_object_with_name<Type>(name);
}



template <typename Type, typename Arg, typename... Args>
Type &
GeneralDataStorage::get_or_add_object_with_name(const std::string &name,
                                                Arg               &argument,
                                                Args &...arguments)
{
  if (!stores_object_with_name(name))
    add_unique_copy(name, Type(argument, arguments...));

  return get_object_with_name<Type>(name);
}



template <typename Type, typename Arg>
Type &
GeneralDataStorage::get_or_add_object_with_name(const std::string &name,
                                                Arg              &&argument)
{
  if (!stores_object_with_name(name))
    add_unique_copy(name, Type(std::forward<Arg>(argument)));

  return get_object_with_name<Type>(name);
}



template <typename Type, typename Arg, typename... Args>
Type &
GeneralDataStorage::get_or_add_object_with_name(const std::string &name,
                                                Arg              &&argument,
                                                Args &&...arguments)
{
  if (!stores_object_with_name(name))
    add_unique_copy(name,
                    Type(std::forward<Arg>(argument),
                         std::forward<Args>(arguments)...));

  return get_object_with_name<Type>(name);
}


template <typename Type>
Type &
GeneralDataStorage::get_or_add_object_with_name(const std::string &name)
{
  if (!stores_object_with_name(name))
    add_unique_copy(name, Type());

  return get_object_with_name<Type>(name);
}


#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif // dealii_algorithms_general_data_storage_h
