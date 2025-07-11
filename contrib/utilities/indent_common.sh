#!/bin/bash
## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2018 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# This file contains a number of common functions used all indent scripts
#

#
# This function checks that we are in the root directory and that
# clang-format is available. It is called by both indent and indent-all
# to ensure that the rest of the indentation pipeline makes sense.
#

#
# DEAL_II_CLANG_FORMAT can be used to change the default version of
# clang-format
#

export DEAL_II_CLANG_FORMAT="${DEAL_II_CLANG_FORMAT:-clang-format}"

checks() {
  if test ! -d source -o ! -d include -o ! -d examples ; then
    echo "*** This script must be run from the top-level directory of deal.II."
    exit 1
  fi

  # Add the location 'download_clang_format' or 'compile_clang_format'
  # installs clang-format to the local PATH.
  CLANG_FORMAT_PATH="$(cd "$(dirname "$0")" && pwd)/programs/clang-16/bin"
  export PATH="${CLANG_FORMAT_PATH}:${PATH}"

  if ! [ -x "$(command -v "${DEAL_II_CLANG_FORMAT}")" ]; then
    echo "***   No clang-format program found, or found with the wrong version."
    echo "***"
    echo "***   You can run the './contrib/utilities/download_clang_format'"
    echo "***   script, or the './contrib/utilities/compile_clang_format' script "
    echo "***   to install a compatible binary into './contrib/utilities/programs'."
    exit 1
  fi

  # Make sure to have the right version.
  CLANG_FORMAT_VERSION="$(${DEAL_II_CLANG_FORMAT} --version)"
  CLANG_FORMAT_MAJOR_VERSION=$(echo "${CLANG_FORMAT_VERSION}" | sed 's/^[^0-9]*\([0-9]*\).*$/\1/g')
  CLANG_FORMAT_MINOR_VERSION=$(echo "${CLANG_FORMAT_VERSION}" | sed 's/^[^0-9]*[0-9]*\.\([0-9]*\).*$/\1/g')

  if [ "${CLANG_FORMAT_MAJOR_VERSION}" -ne 16 ] || [ "${CLANG_FORMAT_MINOR_VERSION}" -ne 0 ]; then
    echo "***   This indent script requires clang-format version 16.0,"
    echo "***   but version ${CLANG_FORMAT_MAJOR_VERSION}.${CLANG_FORMAT_MINOR_VERSION} was found instead."
    echo "***"
    echo "***   You can run the './contrib/utilities/download_clang_format'"
    echo "***   script, or the './contrib/utilities/compile_clang_format' script "
    echo "***   to install a compatible binary into './contrib/utilities/programs'."
    exit 1
  fi


  # check formatting of usernames and email addresses, examples that will be detected:
  # not-using-a-name <a@b.com>
  # John Doe <doe@macbook.local>
  # Jane Doe <a@nodomain>
  #
  # For commits already in the history, please see .mailmap in the root directory.
  #
  # Note that we currently allow email addresses of the form
  # Luca Heltai <luca-heltai@users.noreply.github.com>
  # as these are generated when using the website to commit.
  #
  # Finally, to stay sane, just go back until the beginning of 2019 for now.
  #
  # first user names:
  git log --since "2019-01-01" --format="%aN" --no-merges | sort -u | while read name ; do
      words=($name)
      if [ "${#words[@]}" -lt "2" -a "$name" != "dependabot[bot]" ]; then
	  echo "Invalid author '$name' without firstname and lastname"
	  echo ""
	  echo "For possible solutions, consult the following page:"
	  echo "      https://github.com/dealii/dealii/wiki/Commit-authorship"
	  exit 2
      fi
  done || exit 2

  # now emails:
  git log --since "2019-01-01" --format="%aE" --no-merges | sort -u | while read email ; do
      words=($name)
      if ! echo "$email" | grep -q "\."; then
	  echo "invalid email '$email'"
          echo ""
          echo "hint: for possible solutions, consult the webpage:"
          echo "      https://github.com/dealii/dealii/wiki/Commit-authorship"
	  exit 3
      fi
      if ! echo "$email" | grep -q -v -e "\.local$"; then
	  echo "invalid email '$email'"
          echo ""
          echo "hint: for possible solutions, consult the webpage:"
          echo "      https://github.com/dealii/dealii/wiki/Commit-authorship"
	  exit 3
      fi
  done || exit 3

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
# In order to format .cc and .h files we have to make sure that we overwrite
# the source/header file only if the actual contents changed.
# Unfortunately, clang-format isn't exactly helpful there. Thus, use a
# temporary file and diff as a workaround.
#

format_file()
{
  file="${1}"
  tmpfile="$(mktemp "${TMPDIR}/$(basename "$1").tmp.XXXXXXXX")"

  "${DEAL_II_CLANG_FORMAT}" "${file}" > "${tmpfile}"
  fix_or_report "${file}" "${tmpfile}" "file indented incorrectly"
  rm -f "${tmpfile}"
}
export -f format_file


#
# Also format Python files, assuming that we have the 'black'
# autoindenter available. Black is better than clang-format at not
# changing the time stamp of a file if the file wasn't actually
# changed. But, for simplicity, just use the same scheme above
# so we can use the same reporting mechanism.
#
# If 'black' isn't available, just don't do anything.
#
# (There is one file, contrib/utilities/dotgdbinit.py, that contains
# some Python code, but starts with some GDB instruction that confuses
# the indenter. Exclude that file.)
#
format_python_file()
{
  if which black > /dev/null ; then
    file="${1}"
    if test "$file" = "contrib/utilities/dotgdbinit.py" ; then
      return ;
    fi
      
    tmpfile="$(mktemp "${TMPDIR}/$(basename "$1").tmp.XXXXXXXX")"
    cp "${file}" "${tmpfile}"
    black -q "${tmpfile}"
    fix_or_report "${file}" "${tmpfile}" "file indented incorrectly"
    rm -f "${tmpfile}"
  fi
}
export -f format_python_file


#
# Remove trailing whitespace.
#

remove_trailing_whitespace()
{
  file="${1}"
  tmpfile="$(mktemp "${TMPDIR}/$(basename "$1").tmp.XXXXXXXX")"

  #
  # Mac OS uses BSD sed (other than GNU sed in Linux),
  # so it doesn't recognize \s as 'spaces' or + as 'one or more'.
  #
  sed 's/[[:space:]]\{1,\}$//g' "${file}" > "${tmpfile}"
  if ! diff -q "${file}" "${tmpfile}" >/dev/null; then
    mv "${tmpfile}" "${file}"
  fi
  rm -f "${tmpfile}"
}
export -f remove_trailing_whitespace

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
  "${DEAL_II_CLANG_FORMAT}" < "${tmpfile}" > "${tmpfile}new"

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
  rm -f "${tmpfile}"
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
# Collect all files found in a list of directories "${1}" matching a
# regular expression "${2}", and process them with a command "${3}" on 10
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
# - the option "-P 10" starts up to 10 processes in parallel. -0 implies '-L 1'
#   (one argument to each command) so each launch of clang-format corresponds
#   to exactly one file.
#

process()
{
  directories=$1
  case "${OSTYPE}" in
    darwin*)
      find -E ${directories} -regex "${2}" -print0 |
        xargs -0 -P 10 -I {} bash -c "${3} {}"
      ;;
    *)
      find ${directories} -regextype egrep -regex "${2}" -print0 |
        xargs -0 -P 10 -I {} bash -c "${3} {}"
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
  COMMON_ANCESTOR_WITH_MASTER="$(git merge-base "${LAST_MERGE_COMMIT}" HEAD)"

  case "${OSTYPE}" in
    darwin*)
      XARGS="xargs -E"
      ;;
    *)
      XARGS="xargs --no-run-if-empty -d"
      ;;
  esac

  ( git ls-files --others --exclude-standard -- ${1};
    git diff --name-only $COMMON_ANCESTOR_WITH_MASTER -- ${1} )|
      sort -u |
      xargs -n 1 ls -d 2>/dev/null |
      grep -E "^${2}$" |
      ${XARGS} '\n' -P 10 -I {} bash -c "${3} {}"
}

#
# Ensure only a single newline at end of files
#
ensure_single_trailing_newline()
{
  file="${1}"

  if test $(tail -c1 "$file") ; then
    echo '' >> "$file"
  fi
}
export -f ensure_single_trailing_newline
