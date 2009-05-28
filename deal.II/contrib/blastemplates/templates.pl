#---------------------------------------------------------------------------
#    $Id$
#    Version: $Name$
#
#    Copyright (C) 2005, 2006, 2007, 2008, 2009 by the deal authors
#
#    This file is subject to QPL and may not be  distributed
#    without copyright and license information. Please refer
#    to the file deal.II/doc/license.html for the  text  and
#    further information on this license.
#
#---------------------------------------------------------------------------

my $templates;
my $double;


print << 'EOT'
//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    This file was automatically generated from lapack_templates.h.in
//    See blastemplates in the deal.II contrib directory
//
//    Copyright (C) 2005, 2006, 2007, 2008, 2009 by the deal authors
//
//    This file is subject to QPL and may not be distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __LAPACK_TEMPLATES_H
#define __LAPACK_TEMPLATES_H

#include <lac/lapack_support.h>

using namespace dealii;

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

print "\n}\n";

print "\n$templates\n\n#endif\n";
