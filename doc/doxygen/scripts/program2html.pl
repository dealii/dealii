## ---------------------------------------------------------------------
##
## Copyright (C) 2006 - 2021 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

print "<a name=\"CommProg\"></a>\n";
print "<h1> The commented program</h1>\n";

# ignore comments at the start of the program. this includes subversion
# tags and copyright notices.
$_ = <>;
while ( m!^/[/\*]!  ||  m!\s*\*! || m/^$/ ) {
    $_ = <>;
}

# have two states, in which the program can be:
# comment-mode and program-mode
$comment_mode = 0;
$program_mode = 1;
$state =  $comment_mode;

print "<p>\n";

do {
    # substitute special characters
    s/&/&amp;/g;
    s/</&lt;/g;
    s/>/&gt;/g;
    s/\t/        /g;

    if (($state == $program_mode) && m!^\s*//!)
    {
	$state = $comment_mode;
	print "</code></pre>\n";
	print "\n";
	print "<p>\n";
    }
    # if in comment mode and no comment line: toggle state.
    # don't do so, if only a blank line
    elsif (($state == $comment_mode) && !m!^\s*//! && !m!^\s*$!)
    {
	$state = $program_mode;
	print "</p>\n";
	print "\n";
	print "<pre><code>\n";
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

	   $_ = "\n<a name=\"$sect_name\"></a>" . $_;
	}

	# replace quotation marks by the appropriate HTML quotation marks
	s!``!&#8220;!g;
	s!''!&#8221;!g;

        # replace double dashes in comments by &mdash;
	s!--!&mdash;!g;

	# finally print this line
	print $_, "\n";

	# if empty line, introduce paragraph break
	print "</p>\n\n<p>" if  $_ =~ m!^\s*$!;
    }
    else
    {
	print "        $_";
    }
} while (<>);

if ($state == $program_mode) {
   print "</code></pre>\n";
}
