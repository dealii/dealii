#---------------------------------------------------------------------------
#    $Id$
#    Version: $Name$
#
#    Copyright (C) 2005 by the deal authors
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
//    This file was automatically generated from blas.h.in
//
//    Copyright (C) 2005 by the deal authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __BLAS_TEMPLATES_H
#define __BLAS_TEMPLATES_H

extern "C"
{
EOT
    ;

while(<>)
{
    if (m'^\s*//')
    {
	print;
	next;
    }
    if (m'\s*(\w+)\s+d(\w+)_\s*\(')
    {
	$double = $_;
	my $type = $1;
	my $name = $2;
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
	$args =~ s/\s+/ /g;
	$args2 = $args;
	$args2 =~ s/\w+\*//g;
	
	$templates .= "\n\ninline $type\n$name ($args)\n{\n  d$name\_ ($args2);\n}\n";
	$args =~ s/double/float/g;
	$type =~ s/double/float/g;
	$templates .= "\n\ninline $type\n$name ($args)\n{\n  s$name\_ ($args2);\n}\n";
    }
}

print "\n}\n";

print "\n$templates\n\n#endif\n";
