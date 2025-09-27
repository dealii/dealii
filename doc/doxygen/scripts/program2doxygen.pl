## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2006 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------


# Skip header lines at the top of the file, such as copyright notices and
# license information, if the file is a step-xx.cc tutorial. Here, xx
# is either a number like step-32 or a number plus a letter like step-12b.
# Don't skip for other files such as code-gallery files
use Getopt::Long;
GetOptions("prefix:s" => \$prefix);

if ($ARGV[0] =~ /step-\d+[a-z]?.cc/)
{
  $_ = <>;
  while ( m!^/\*!  ||  m!\s*\*! || m/^$/ ) {
      $_ = <>;
  }
}


# have four states, in which the program can be:
# comment-mode, program-mode, skip-mode, and code-fragment-mode
$comment_mode          = 0;
$program_mode          = 1;
$skip_mode             = 2;
$code_fragment_mode    = 3;
$state =  $comment_mode;

print " * \n";

do {

    # substitute tabs
    s/\t/        /g;

    # Remove "//" if it is present at the end of a line. These comments exist
    # to force a specific formatting of the line for clang-format.
    s!//$!!g;

    # We are entering code fragments
    if (($state != $code_fragment_mode) && m!^\s*//\s*\@code!)
    {
        # Cleanly close @code segments in program mode:
        if ($state == $program_mode)
        {
            print " * \@endcode\n";
            print " * \n";
        }
        $state = $code_fragment_mode;
        # Make sure we distinguish from the other code style
        print " * <div class=CodeFragmentInTutorialComment>\n";
    }
    # We shall skip something...
    elsif (($state != $skip_mode) && m!^\s*//\s*\@cond SKIP!)
    {
        # Cleanly close @code segments in program mode:
        if ($state == $program_mode)
        {
            print " * \@endcode\n";
            print " * \n";
        }
        $state = $skip_mode;
    }
    elsif (($state == $program_mode) && m!^\s*//!)
    {
        print " * \@endcode\n";
        print " * \n";
        $state = $comment_mode;
    }
    # if in comment mode and no comment line: toggle state.
    # don't do so, if only a blank line
    elsif (($state == $comment_mode) && !m!^\s*//! && !m!^\s*$!)
    {
        print " * \n";
        print " * \@code\n";
        $state = $program_mode;
    }

    if ($state == $comment_mode)
    {
        # in comment mode: first skip leading whitespace and
        # comment // signs
        s!\s*//\s*(.*)\n!$1!;

        # second, replace section headers, and generate addressable
        # anchor
        if ( /\@sect/ ) {
           s!\@sect(\d)\{(.*)\}\s*$!<h$1>$2</h$1>!g;
           $sect_name = $2;

           # for the anchor, use the name of the section but discard
           # everything except for letters, numbers, and underscores
           $sect_name =~ s/[^a-zA-Z0-9_]//g;

           $_ = "\n * <a name=\"$prefix-$sect_name\"></a> \n * $_";
        }

        # finally print this line
        print " * $_\n";

        # if empty line, introduce paragraph break
        print " * \n" if  $_ =~ m!^\s*$!;
    }
    elsif ($state == $program_mode)
    {
        # in program mode, output the program line. the only thing we need
        # to do is to avoid $ signs because that confuses doxygen. since
        # we don't want formulas rendered in the program text anyway,
        # simply replace them by spaces (it would be nice to suppress their
        # meaning somehow, but I don't know how...)
        s/\$//g;

        # Then print the line. doxygen has the annoying habit of eating
        # the maximal number of spaces at the front of each code block,
        # leading to visually wrong indentation if one, for example, has
        # ```
        #   if (cond)
        #   {
        # ```
        # in one block, and then
        # ```
        #     some_function();
        # ```
        # in the next block -- the call to some_function() is not shown any
        # further to the right than the if(cond) before. Work around this by
        # prefixing all code lines with a non-printing Unicode space 0x00A0
        # that doxygen interprets as the first non-space character for
        # determining indentation, but that does not actually print as anything
        # other than a space (and that compilers appear to successfully ignore
        # when one copy-pastes code snippets from the generated doxygen pages
        # into an editor).
        print " *   $_";   # Note the (invisible) Unicode space after the '* '
    }
    elsif ($state == $skip_mode)
    {
        # This is the end of a @cond - @endcond block, so back to
        # comment_mode:
        if (m!^\s*//\s*\@endcond!)
        {
            $state = $comment_mode;
        }
    }
    elsif ($state == $code_fragment_mode)
    {
        # If this is the end of a @code - @endcode block, go back to
        # comment_mode:
        if (m!^\s*//\s*\@endcode!)
        {
            $state = $comment_mode;
        }
        # in code fragment mode: only skip leading whitespace and
        # comment // signs
        s!^\s*//(.*)\n!$1!;

        # finally print this line
        print " *$_\n";

        if (m!^ *\s*\@endcode!)
        {
            print " * </div>\n";
        }
   }
 } while (<>);

if ($state == $program_mode) {
   print " * \@endcode\n";
}
