## ---------------------------------------------------------------------
##
## Copyright (C) 2007 - 2014 by the deal.II authors
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


# If we find a heading in a .dox file, create an HTML anchor for it.

while (<>) {
    if ( /<h.>(.*)<\/h.>\s*/ ) {
	$reftext = $1;

	# for the anchor, use the name of the section but discard
	# everything except for letters, numbers, and underscores
	$reftext =~ s/[^a-zA-Z0-9_]//g;

	print "<a name=\"$reftext\"></a>$_\n";
    } else {
        print;
    }
}
