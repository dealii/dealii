#!/bin/bash
## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2018 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

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
  if ${DEBUG_OUTPUT}; then
    local TIMEFORMAT="DEBUG_OUTPUT: checks: %Rs"
    time DEBUG_OUTPUT=false checks "$@"
    return
  fi

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
}

#
# Print the commit number at which the current feature branch forked off
# from the branch "master" branch (or any other branch given by
# DEAL_II_BASE_BRANCH). We use this commit in process_changed() to identify
# all changed files on the feature branch.
#

feature_branch_base()
{
  if ${DEBUG_OUTPUT}; then
    local TIMEFORMAT="DEBUG_OUTPUT: feature_branch_base: %Rs"
    time DEBUG_OUTPUT=false feature_branch_base "$@"
    return
  fi

  #
  # We need to use a heuristic here:
  #
  # We want to identify a reasonable base commit, but this is not as easy
  # as it seems. We, collect *every* local and remote branch called
  # "master" (or DEAL_II_BASE_BRANCH) and compute the "merge base" (see
  # $ man git-merge-base) between our current branch and the candidate
  # "master" branch. We then use the "most recent" one as determined by
  # $ git merge-base --is-ancestor.
  #
  BASE=""
  for ref in $(git for-each-ref --format="%(refname)" \
                 "refs/heads/${DEAL_II_BASE_BRANCH:-master}" \
                 "refs/remotes/*/${DEAL_II_BASE_BRANCH:-master}"); do

    # Determine the most recent merge commit:
    LAST_MERGE_COMMIT="$(git log --format="%H" --merges --max-count=1 "${ref}")"
    [ -n "${LAST_MERGE_COMMIT}" ] || continue

    # Compute merge base with HEAD:
    CANDIDATE="$(git merge-base "${LAST_MERGE_COMMIT}" HEAD)" || continue

    # All candidates are ancestors of HEAD. Keep the most recent one.
    if [ -z "${BASE}" ] || git merge-base --is-ancestor "${BASE}" "${CANDIDATE}"; then
      BASE="${CANDIDATE}"
    fi
  done
  echo "${BASE}"
}

#
# Print a regular expression that matches a string if (and only if) the
# whole string matches one of the POSIX extended regular expressions
# listed in file "${1}", one per line. Empty lines and lines starting with
# "#" are ignored. Return 1 if one of the expressions is invalid.
#

construct_exception_regex()
{
  regex=""
  while read -r entry; do
    case "${entry}" in
      ""|"#"*) continue ;;
    esac
    # Test the entry: bash returns 2 for an invalid regular expression
    { [[ "" =~ ${entry} ]]; } 2>/dev/null
    if [ "$?" -eq 2 ]; then
      echo "*** Invalid regular expression '${entry}' in ${1}" >&2
      return 1
    fi
    regex="${regex:+${regex}|}(${entry})"
  done < "${1}"
  echo "^(${regex:-})$"
}

#
# Check that commit authors have a proper name (consisting of at least two
# words) and a proper email address for all commits given by the git log
# arguments "${@}". Examples that will be detected:
#   not-using-a-name <a@b.com>
#   John Doe <doe@macbook.local>
#   Jane Doe <a@nodomain>
#   John Smith <12345678+JohnSmith@users.noreply.github.com>
#
# Anonymized github addresses carry no contact information and are thus
# rejected as well.
#
# Commit authors in the history can be fixed by an entry in .mailmap in
# the root directory. Authors that cannot be fixed this way are listed in
# contrib/utilities/indent_author_exceptions.
#

check_commit_authors()
{
  if ${DEBUG_OUTPUT}; then
    local TIMEFORMAT="DEBUG_OUTPUT: check_commit_authors [${1}]: %Rs"
    time DEBUG_OUTPUT=false check_commit_authors "$@"
    return
  fi

  exceptions="$(construct_exception_regex contrib/utilities/indent_author_exceptions)" || return 1

  failed="$(
    git log --no-merges --format="%aN <%aE>" "${@}" | sort -u | while read -r identity; do
      # Skip if we match an exception from "indent_author_exceptions":
      if [[ "${identity}" =~ ${exceptions} ]]; then
        continue
      fi

      name="${identity% <*}"
      email="${identity##*<}"
      email="${email%>}"

      read -r -a words <<< "${name}"
      if [ "${#words[@]}" -lt "2" ]; then
        echo "    invalid author '${name}' without firstname and lastname"
      fi

      case "${email}" in
        *.*) ;;
        *) echo "    invalid email '${email}'" ;;
      esac

      case "${email}" in
        *.local)
          echo "    invalid email '${email}'" ;;
        *users.noreply.github.com)
          echo "    invalid email '${email}' (anonymized github address)" ;;
      esac
    done
  )"

  if [ -n "${failed}" ]; then
    echo "*** The following commit authors are invalid:"
    echo ""
    echo "${failed}"
    echo ""
    echo "hint: anonymized github addresses are not accepted. Turn off"
    echo "      'Keep my email address private' in your github settings,"
    echo "      or set a real address with 'git config user.email', and"
    echo "      then amend or rebase the commits of this pull request."
    echo ""
    echo "hint: for possible solutions, consult the webpage:"
    echo "      https://github.com/dealii/dealii/wiki/Commit-authorship"
    return 1
  fi
}

#
# Check that file "${1}" starts with our license header, using the comment
# prefix "${2}" ("//" or "##"). If "${3}" is "true" the header may be
# preceded by a "#!" line, optionally followed by a single empty line. The
# header has to match exactly, except that additional "Copyright (C) ..."
# lines may follow the one of the deal.II authors. The function prints the
# first line that does not match.
#
# This is written in plain awk without regex intervals so that it works
# with the awk variants found on Linux and Mac OS.
#

check_license_header_of_file()
{
  awk -v p="${2}" -v allow_shebang="${3}" -v file="${1}" '
    BEGIN {
      dashes = "-----------------------------------------------------------------------------"
      n = 0
      want[++n] = p " " dashes
      want[++n] = p
      want[++n] = p " SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later"
      want[++n] = p " Copyright (C) YYYY - YYYY by the deal.II authors"
      copyright_line = n
      want[++n] = p
      want[++n] = p " This file is part of the deal.II library."
      want[++n] = p
      want[++n] = p " Detailed license information governing the source code and contributions"
      want[++n] = p " can be found in LICENSE.md and CONTRIBUTING.md at the top level directory."
      want[++n] = p
      want[++n] = p " " dashes

      year = "[0-9][0-9][0-9][0-9]"
      dealii_copyright = "^" p " Copyright [(]C[)] " year "( - " year ")? by the deal[.]II authors$"
      other_copyright = "^" p " Copyright [(]C[)] ."

      i = 1
      done = 0
      line = 0
      shebang = 0
    }
    {
      sub(/\r$/, "")
    }
    NR == 1 && allow_shebang == "true" && /^#!/ {
      shebang = 1
      next
    }
    NR == 2 && shebang && $0 == "" {
      next
    }
    {
      if (i == copyright_line) {
        if ($0 !~ dealii_copyright) {
          line = NR
          exit
        }
      } else if (i == copyright_line + 1 && $0 ~ other_copyright) {
        next
      } else if ($0 != want[i]) {
        line = NR
        exit
      }
      if (++i > n) {
        done = 1
        exit
      }
    }
    END {
      if (!done) {
        # If we ran out of lines, the expected line is the one after the last:
        if (line == 0)
          line = NR + 1
        printf "    %s:%d: expected \"%s\"\n", file, line, want[i]
        exit 1
      }
    }' "${1}"
}
export -f check_license_header_of_file


#
# Check that file "${1}" starts with our license header (see
# check_license_header_of_file() above). The comment prefix is determined
# from the file name. Files that match one of the regular expressions in
# contrib/utilities/indent_license_exceptions are skipped. The function
# expects these regular expressions to be combined in the variable
# LICENSE_HEADER_EXCEPTIONS (see setup_license_header_exceptions() below).
#

check_license_header()
{
  file="${1}"

  if [ -z "${LICENSE_HEADER_EXCEPTIONS:-}" ]; then
    echo "*** check_license_header: LICENSE_HEADER_EXCEPTIONS is not set" >&2
    return 1
  fi

  # Intentionally unquoted so that ${LICENSE_HEADER_EXCEPTIONS} is used as
  # regular expression:
  if [[ "${file}" =~ ${LICENSE_HEADER_EXCEPTIONS} ]]; then
    return 0
  fi

  case "${file}" in
    *.h|*.cc|*.ccm|*.cu|*.cuh|*.inst.in)
      check_license_header_of_file "${file}" "//" "false" ;;
    *.cmake|CMakeLists.txt|*/CMakeLists.txt)
      check_license_header_of_file "${file}" "##" "false" ;;
    *.py|*.sh|*.pl)
      check_license_header_of_file "${file}" "##" "true" ;;
  esac
}
export -f check_license_header

#
# Read the regular expressions in contrib/utilities/indent_license_exceptions
# once and export them combined into a single regular expression for
# check_license_header(), which runs in separate processes:
#

setup_license_header_exceptions()
{
  LICENSE_HEADER_EXCEPTIONS="$(construct_exception_regex contrib/utilities/indent_license_exceptions)" || return 1
  export LICENSE_HEADER_EXCEPTIONS
}

#
# Print the license header that check_license_header() expects:
#

print_license_header_template()
{
  echo ""
  echo "*** The files listed above do not start with a valid license header."
  echo "*** All files have to start with the following header (with \"##\""
  echo "*** instead of \"//\" for CMake, Python, Perl and shell scripts, where a"
  echo "*** \"#!\" line, optionally followed by an empty line, may come first),"
  echo "*** where YYYY are the years of creation and last modification, or only"
  echo "*** a single year:"
  echo ""
  echo "// -----------------------------------------------------------------------------"
  echo "//"
  echo "// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later"
  echo "// Copyright (C) YYYY - YYYY by the deal.II authors"
  echo "//"
  echo "// This file is part of the deal.II library."
  echo "//"
  echo "// Detailed license information governing the source code and contributions"
  echo "// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory."
  echo "//"
  echo "// -----------------------------------------------------------------------------"
  echo ""
  echo "*** Additional \"// Copyright (C) ...\" lines may follow the line for the"
  echo "*** deal.II authors."
}

#
# Mac OSX's mktemp doesn't know the --tmpdir option without argument. So,
# let's do all of that by hand:
#
export TMPDIR="${TMPDIR:-/tmp}"

export REPORT_ONLY="${REPORT_ONLY:-false}"

export DEBUG_OUTPUT="${DEBUG_OUTPUT:-false}"

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
  if which black > /dev/null 2>&1 ; then
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
  if ${DEBUG_OUTPUT}; then
    local TIMEFORMAT="DEBUG_OUTPUT: process ${3} [${1}] [${2}]: %Rs"
    time DEBUG_OUTPUT=false process "$@"
    return
  fi

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
# in the current feature branch (see feature_branch_base()). For this, we
# collect all files that
#  - are new
#  - have changed since the feature branch forked off
#

process_changed()
{
  if ${DEBUG_OUTPUT}; then
    local TIMEFORMAT="DEBUG_OUTPUT: process_changed ${3} [${1}] [${2}]: %Rs"
    time DEBUG_OUTPUT=false process_changed "$@"
    return
  fi

  BASE="$(feature_branch_base)"

  case "${OSTYPE}" in
    darwin*)
      XARGS="xargs -E"
      ;;
    *)
      XARGS="xargs --no-run-if-empty -d"
      ;;
  esac

  ( git ls-files --others --exclude-standard -- ${1};
    git diff --name-only ${BASE} -- ${1} )|
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
