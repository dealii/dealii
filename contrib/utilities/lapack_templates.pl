## ---------------------------------------------------------------------
##
## Copyright (C) 2005 - 2013 by the deal.II authors
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
# This perl script translates lapack_templates.h.in to lapack_templates.h
#
# In the *.in file, every BLAS/LAPACK function which is defined for
# double precision, i.e. having a name like 'dfoo_', is expanded to
# itself plus the same function for single precision, namely
# 'sfoo_'. Additionally, a C++ function 'foo' without the prefix
# letter and the trailing underscore is generated, such that the
# fortran functions can easily be called from templates. The
# implementation of this function is modified due to the configure
# variables 'HAVE_DFOO_' and 'HAVE_DFOO_': if these are set, then the
# lapack functions 'dfoo_' and 'sfoo_' will be called, if not, an
# exception will be thrown.
#
# Therefore, in order to be able to call a LAPACK function, the
# functions have to be tested by configure. Search for the section
# "Check for LAPACK..." in deal.II/configure.in and add the functions
# 'dfoo_' and 'sfoo_' to the tests at the end of that section.
#


my $templates;
my $double;


print << 'EOT'
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __LAPACK_TEMPLATES_H
#define __LAPACK_TEMPLATES_H

#include <deal.II/base/config.h>
#include <deal.II/lac/lapack_support.h>

extern "C"
{
EOT
    ;

while(<>)
{
    # Write comment lines literally
    if (m'^\s*//')
    {
	print;
	next;
    }
    # Lines of the form 'typename functionname (...'
    # where functionname is of the form d..._,
    # that is a double precision LAPACK function
    if (m'\s*(\w+)\s+d(\w+)_\s*\(')
    {
	$double = $_;
	my $type = $1;
	my $name = $2;
	my $capname = $name;
	$capname =~ tr/[a-z]/[A-Z]/;
	while (<>)
	{
	    $double .= $_;
	    last if (m';');
	}
	my $single = $double;
	$single =~ s/d$name/s$name/;
	$single =~ s/double/float/g;
	print $double,$single;

	$double =~ m/\(([^\)]*)/;
	my $args = $1;
	# The arglist for the C++ function
	$args =~ s/\s+/ /g;
	# The arglist handed down to the FORTRAN function
	$args2 = $args;
	# Fortunately, all arguments are pointers, so we can use the *
	# to separate data type and argument name
	$args2 =~ s/\w+\*//g;
	$args2 =~ s/const//g;
	$args2 =~ s/\s//g;
	# The arglist of the empty C++ function
	$args0 = $args;
	$args0 =~ s/\*[^,]*,/\*,/g;
	$args0 =~ s/\*[^,]*$/\*/g;
	
	# First, do the general template None of these functions is
	# implemented, but they allow us to link for instance with
	# long double lapack support
	my $numbers = 1;
	my $argst = $args0;
	my $typet = $type;
	while ($argst =~ s/double/number$numbers/)
	{
	    $numbers++;
	}
	$typet =~ s/double/number1/g;

	$templates .= "\n\n/// Template wrapper for LAPACK functions d$name and s$name\n";
	$templates .= "template<typename number1";
	for (my $i=2;$i<$numbers;++$i)
	{
	    $templates .= ", typename number$i";
	}
	$templates .= ">\ninline $typet\n$name ($argst)\n";
	$templates .= "{\n  Assert (false, ExcNotImplemented());\n}";

	# The specialization for double. Note that we always have two
	# versions, one implemented and calling LAPACK, the other not
	# implemented
	$templates .= "\n\n#ifdef HAVE_D$capname\_";
	$templates .= "\ninline $type\n$name ($args)\n{\n  d$name\_ ($args2);\n}\n";
	$templates .= "#else\ninline $type\n$name ($args0)\n";
	$templates .= "{\n  Assert (false, LAPACKSupport::ExcMissing(\"d$name\"));\n}\n#endif\n";

	$args =~ s/double/float/g;
	$args0 =~ s/double/float/g;
	$type =~ s/double/float/g;
	$templates .= "\n\n#ifdef HAVE_S$capname\_";
	$templates .= "\ninline $type\n$name ($args)\n{\n  s$name\_ ($args2);\n}\n";
	$templates .= "#else\ninline $type\n$name ($args0)\n";
	$templates .= "{\n  Assert (false, LAPACKSupport::ExcMissing(\"s$name\"));\n}\n#endif\n";
    }
}

print "\n}\n\nDEAL_II_NAMESPACE_OPEN\n";

print "$templates\n\nDEAL_II_NAMESPACE_CLOSE\n\n#endif\n";
