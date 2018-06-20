#!/bin/bash
## ---------------------------------------------------------------------
##
## Copyright (C) 2018 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# This file contains a number of common functions used all indent scripts
#

#
# This function checks that we are in the root directory and that
# clang-format is available. It is called by both indent and indent-all
# to ensure that the rest of the indentation pipeline makes sense.
#

checks() {
  if test ! -d source -o ! -d include -o ! -d examples ; then
    echo "*** This script must be run from the top-level directory of deal.II."
    exit 1
  fi

  # Add the location 'download_clang_format' or 'compile_clang_format'
  # installs clang-format to to the local PATH.
  CLANG_FORMAT_PATH="$(cd "$(dirname "$0")" && pwd)/programs/clang-6/bin"
  export PATH="${CLANG_FORMAT_PATH}:${PATH}"

  if ! [ -x "$(command -v clang-format)" ]; then
    echo "***   No clang-format program found."
    echo "***"
    echo "***   You can run the './contrib/utilities/download_clang_format'"
    echo "***   script, or the './contrib/utilities/compile_clang_format' script "
    echo "***   to install a compatible binary into './contrib/utilities/programs'."
    exit 1
  fi

  # Make sure to have the right version. We know that clang-6.0.0
  # and clang-6.0.1 work. Hence, test for clang-6.0.
  CLANG_FORMAT_VERSION="$(clang-format --version)"
  CLANG_FORMAT_MAJOR_VERSION=$(echo "${CLANG_FORMAT_VERSION}" | sed 's/^[^0-9]*\([0-9]*\).*$/\1/g')
  CLANG_FORMAT_MINOR_VERSION=$(echo "${CLANG_FORMAT_VERSION}" | sed 's/^[^0-9]*[0-9]*\.\([0-9]*\).*$/\1/g')

  if [ "${CLANG_FORMAT_MAJOR_VERSION}" -ne 6 ] || [ "${CLANG_FORMAT_MINOR_VERSION}" -ne 0 ]; then
    echo "***   This indent script requires clang-format version 6.0,"
    echo "***   but version ${CLANG_FORMAT_MAJOR_VERSION}.${CLANG_FORMAT_MINOR_VERSION} was found instead."
    echo "***"
    echo "***   You can run the './contrib/utilities/download_clang_format'"
    echo "***   script, or the './contrib/utilities/compile_clang_format' script "
    echo "***   to install a compatible binary into './contrib/utilities/programs'."
    exit 1
  fi
}

#
# Mac OSX's mktemp doesn't know the --tmpdir option without argument. So,
# let's do all of that by hand:
#
export TMPDIR="${TMPDIR:-/tmp}"

export REPORT_ONLY="${REPORT_ONLY:-false}"

#
# If REPORT_ONLY is set to "true", this function reports a formatting issue
# if file "${1}" and tmpfile "${2}" don't match (using the error message
# "${3}"), or, if set to "false" silently replaces file "${1}" with "${2}".
#

fix_or_report()
{
  file="${1}"
  tmpfile="${2}"
  message="${3}"

  if ! diff -q "${file}" "${tmpfile}" >/dev/null; then
    if ${REPORT_ONLY}; then
      echo "    ${file}  -  ${message}"
    else
      mv "${tmpfile}" "${file}"
    fi
  fi
}
export -f fix_or_report

#
# In order to format .cc and .h files we have to make sure that we override
# the source/header file only if the actual contents changed.
# Unfortunately, clang-format isn't exactly helpful there. Thus, use a
# temporary file and diff as a workaround.
#

format_file()
{
  file="${1}"
  tmpfile="$(mktemp "${TMPDIR}/$(basename "$1").tmp.XXXXXXXX")"

  clang-format "${file}" > "${tmpfile}"
  fix_or_report "${file}" "${tmpfile}" "file indented incorrectly"
  rm -f "${tmpfile}"
}
export -f format_file

#
# In order to format .inst.in files, we need to replace \{ and \} by a
# sentinel because clang-format happily strips away the backslash. Further,
# there is a minor caveat: clang-format automatically renames a comment
# after a closing bracket for a namespace to "} // namespace FooBar". Thus
# use "namespace" as keyword:
#

format_inst()
{
  file="${1}"
  tmpfile="$(mktemp "${TMPDIR}/$(basename "$1").tmp.XXXXXXXX")"

  cp "${file}" "${tmpfile}"
  sed -i -e 's#\\{#{ // namespace#g' "${tmpfile}"
  sed -i -e 's#\\}#} // namespace#g' "${tmpfile}"

  #
  # Yes, we have to call clang-format in this weird way to ensure that it
  # picks up and uses the .clang-format style-sheet in the current
  # directory.
  #
  clang-format < "${tmpfile}" > "${tmpfile}new"

  sed -i -e 's#{ // namespace#\\{#g' "${tmpfile}new"
  sed -i -e 's#}[ ]*// namespace.*#\\}#g' "${tmpfile}new"

  fix_or_report "${file}" "${tmpfile}new" "file indented incorrectly"
  rm -f "${tmpfile}" "${tmpfile}new"
}
export -f format_inst

#
# Convert DOS formatted files to unix file format by stripping out
# carriage returns (15=0x0D):
#

dos_to_unix()
{
  file="${1}"
  tmpfile="$(mktemp "${TMPDIR}/$(basename "$1").tmp.XXXXXXXX")"

  tr -d '\015' <"${file}" >"${tmpfile}"

  fix_or_report "${file}" "${tmpfile}" "file has non-unix line-ending '\\r\\n'"
  rm -f "${tmpfile}" "${tmpfile}"
}
export -f dos_to_unix

#
# Fix permissions
#

fix_permissions()
{
  file="${1}"

  case "${OSTYPE}" in
    darwin*)
      PERMISSIONS="$(stat -f '%a' ${file})"
      ;;
    *)
      PERMISSIONS="$(stat -c '%a' ${file})"
      ;;
  esac

  if [ "${PERMISSIONS}" != "644" ]; then
    if ${REPORT_ONLY}; then
      echo "    ${file}  -  file has incorrect permissions"
    else
      chmod 644 "${file}"
    fi
  fi
}
export -f fix_permissions

#
# Collect all files found in a list of directories "${1}$" matching a
# regular expression "${2}$", and process them with a command "${3}" on 10
# threads in parallel.
#
# The command line is a bit complicated, so let's discuss the more
# complicated arguments:
# - For 'find', -print0 makes sure that file names are separated by \0
#   characters, as opposed to the default \n. That's because, surprisingly,
#   \n is a valid character in a file name, whereas \0 is not -- so it
#   serves as a good candidate to separate individual file names.
# - For 'xargs', -0 does the opposite: it separates filenames that are
#   delimited by \0
# - the options "-n 1 -P 10" make sure that the following script with be
#   called exactly with one file name as argument at a time, but we allow
#   execution for up to 10 times in parallel
#

process()
{
  case "${OSTYPE}" in
    darwin*)
      find -E ${1} -regex "${2}" -print0 |
        xargs -0 -n 1 -P 10 -I {} bash -c "${3} {}"
      ;;
    *)
      find ${1} -regextype egrep -regex "${2}" -print0 |
        xargs -0 -n 1 -P 10 -I {} bash -c "${3} {}"
      ;;
  esac
}

#
# Variant of above function that only processes files that have changed
# since the last merge commit to master. For this, we collect all files
# that
#  - are new
#  - have changed since the last merge commit to master
#

process_changed()
{
  LAST_MERGE_COMMIT="$(git log --format="%H" --merges --max-count=1 master)"

  case "${OSTYPE}" in
    darwin*)
      XARGS="xargs -E"
      ;;
    *)
      XARGS="xargs --no-run-if-empty -d"
      ;;
  esac

  ( git ls-files --others --exclude-standard -- ${1};
    git diff --name-only $LAST_MERGE_COMMIT -- ${1} )|
      sort -u |
      grep -E "^${2}$" |
      ${XARGS} '\n' -n 1 -P 10 -I {} bash -c "${3} {}"
}
