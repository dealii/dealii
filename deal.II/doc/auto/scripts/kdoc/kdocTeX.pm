#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
package kdocTeX;

require Ast;

## LaTeX Output.
# Peeks in main: quiet, version, classSource

BEGIN {
	$modname = "kdocTeX";
}

sub dumpDoc {
	my( $name, $nodem $outputdir ) = @_;

	print "Generating LaTeX documentation.\n" unless $main::quiet;

	open(TEX,"> $outputdir/$name.tex") 
		|| die "Couldn't create $outputdir/".$name.".tex.\n";

	print TEX<<EOF;
\\documentclass[12pt]\{article\}
\\usepackage\{a4\}
\\usepackage\{makeidx\}
\\title\{$name Library Index\}
\\date\{\\today\}
\\def\\mtl{{\\tt\\~\\relax}}
\\author{Generated with KDOC $main::version}
\\makeindex
\\begin\{document\}
\\maketitle
\\pagestyle\{headings\}
\\tableofcontents
\\newpage
\\sloppy
EOF

	$node->Visit($modname);

	foreach $class ( @{$classes} )  {

		$class->Visit($modname);

		$inPre = 0;
		$doc = "";

		foreach $line ( split /\n/, $Description ) {
			if( $line =~ /<(pre|PRE)>/ ) {
				$inPre = 1;
				$doc .= "\\begin{verbatim}\n";
				next;
			}
			if( $line =~ m#</(pre|PRE)># ) {
				$inPre = 0;
				$doc .= "\\end{verbatim}\n";
				next;
			}
			if( !$inPre ) {
				$line =~ s/\@ref(\s*(#|::))?//g;
				$line =~ s/([~#\$\\\&%#_\{\}])/\\$1/g;
				$line =~ s/<(code|pre)>/\\texttt\{/g;
				$line =~ s/<\/(code|pre)>/}/g;
				$line =~ s/([<>])/\\texttt\{$1\} /g;
				$line =~ s/~/\\texttt{~}\\relax /g;
			}

			$doc .= $line."\n";
		}

		$ClassName = escape( $astNodeName );
		$Header =  escape( $Header );
		$ClassSee = escape( $ClassSee );

		print TEX<<EOF;
\\section[$ClassName]\{\\emph\{$ClassName\} class reference \\small\\sf ($Header)\}
\\vspace{-0.6cm}
\\hrulefill
\\index\{$ClassName\}
EOF

		if( $class->{ "TmplArgs" } ne "" ) {
			print TEX "\\begin{description}\n", 
				"\\item{template form:}",
				"\\texttt{template <", escape( $class->{ "TmplArgs" } ),
				"> ",
				 $ClassName, "}\n";

			print TEX "\\end{description}\n";
		}

		if( defined $class->{ "Ancestors" } ){
			print TEX "\\begin{description}\n",
			"\\item{inherits:}";

			foreach $ancestor ( @{$Ancestors} ) {
				$ancName = $ancestor->{"astNodeName"};

				$str = " ".$ancName;

				if( defined $main::classSource{ $ancName } ) {
					$str .="(".$main::classSource{$ancName}.
							")";
				}

				print TEX escape( $str ),"\n";
			}

			print TEX "\\end{description}\n";
		}


		if( $Author ne "" || $Version ne "" || $ClassShort ne "" 
				|| $ClassSee ne "" ) {
			print TEX "\\begin{description}\n";

			print TEX "\\item[Description:] ", escape( $ClassShort ),"\n"
				if $ClassShort ne "";
			print TEX "\\item[Version:] ", escape( $Version ),"\n"
				if $Version ne "";

			print TEX "\\item[Author:] ", escape( $Author ),"\n"
				if $Author ne "";		

			print TEX "\\item[See Also:] ", escape( $ClassSee ),"\n"
				if $ClassSee ne "";		

			print TEX "\\end{description}\n";
		}

		print TEX "$doc\n\n" if $Description ne "";

		dumpMembers( "public members", $public ) 
			if defined $class->{"public"};
		dumpMembers( "public slots", $public_slots )
			if defined $class->{"public_slots"};
		dumpMembers( "protected members", $protected )
			if defined $class->{"protected"};
		dumpMembers( "protected slots", $protected_slots )
			if defined $class->{"protected_slots"};
		dumpMembers( "signals", $signals )
			if defined $class->{"signals"};

		if( $main::dontDoPrivate == 0 ) {
			dumpMembers( "private members", $private )
				if defined $class->{"private"};
			dumpMembers( "private slots", $private_slots )
				if defined $class->{"private_slots"};
		}

		Ast::UnVisit();
	}

	Ast::UnVisit();

	print TEX "\\printindex\n\\end{document}\n";

}

sub dumpMembers
{
	my( $access, $nodes ) = @_;

	print TEX "\\subsection{$access}\n";

	foreach $member ( @{$nodes} ) {

		$member->Visit($modname);
		
		$inPre = 0;
		$doc = "";

		foreach $line ( split /\n/, $Description ) {
			if( $line =~ /<(pre|PRE)>/ ) {
				$inPre = 1;
				$doc .= "\\begin{verbatim}\n";
				next;
			}
			if( $line =~ /<\/(pre|PRE)>/ ) {
				$inPre = 0;
				$doc .= "\\end{verbatim}\n";
				next;
			}
			if( !$inPre ) {
				$line =~ s/\@ref(\s*(#|::))?//g;
				$line =~ s/([~#\$\\\&%#_\{\}])/\\$1/g;
				$line =~ s/<(code|pre)>/\\texttt\{/g;
				$line =~ s/<\/(code|pre)>/}/g;
				$line =~ s/([<>])/\\texttt\{$1\} /g;
				$line =~ s/~/\\texttt{~}\\relax /g;
			}

			$doc .= $line."\n";
		}

		$astNodeName = escape( $astNodeName );


		print TEX<<EOF;
\\subsubsection*{$ClassName\:\:$astNodeName}
\\vspace{-0.65cm}
\\hrulefill
\\index{$ClassName!$astNodeName}
\\index{$astNodeName!$ClassName}
\\begin{flushleft}
EOF

		if( $Keyword eq "method" ) {
			$ReturnType = escape( $ReturnType );
			$Parameters = escape( $Parameters );
			$Parameters =~ s/\n+/ /g;

			print TEX "\\texttt{$ReturnType$astNodeName(",
				"$Parameters)$Const;}\n";
		}
		elsif ( $Keyword eq "property" ) {
			$Type = escape( $Type );
			print TEX "\\texttt{$Type$astNodeName;}\n";
		}
		elsif ( $Keyword eq "typedef" ) {
			$Type = escape( $Type );
			print TEX "\\texttt{$Type$astNodeName;}\n";
		}
		elsif ( $Keyword eq "enum" ) {
			$Constants = escape( $Constants );
			$Constants =~ s/\n+/ /g;

			print TEX<<EOF
\\texttt{enum $astNodeName\\{
$Constants
\\};}

EOF
		}

		print TEX<<EOF;
\\end{flushleft}

$doc

EOF
		if ( $Keyword eq "method" && 
			($Returns ne "" || $Exceptions ne "") ) {

			if ( $Returns ne "" ) {
				$Returns = escape ( $Returns );
				print TEX<<EOF;
\\begin{description}
\\item[Returns:] $Returns
\\end{description}
EOF
			}

			if ( $Exceptions ne "" ) {
				$Exceptions = escape ( $Exceptions );
				print TEX<<EOF;
\\begin{description}
\\item[Throws:] $Exceptions
\\end{description}
EOF
			}

		}
		Ast::UnVisit();
	}
}

sub escape
{
	my( $str ) = @_;

	$str =~ s/\@ref\s+//g;
	$str =~ s/([#\$\\\&%#_\{\}])/\\$1/g;
	$str =~ s/([<>]+)/\\texttt\{$1\}/g;
	$str =~ s/~/\\mtl /g;

	return $str;
}

1;
