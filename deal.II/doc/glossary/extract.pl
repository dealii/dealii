$/ = "\n\n";

print << 'EOT'
<HTML>
<HEAD>
<TITLE>deal.II glossary</TITLE>
<link rel=stylesheet type="text/css" href="glossary.css">
</HEAD>
<BODY>
<H1><ACRONYM>deal.II</ACRONYM> Glossary</H1>

<DL>
EOT
    ;

while(<>)
{
    s/\s+/ /g;
    while ( /<<([^>]+)>>/g ) {
	$link_name = $1;
	$link_name =~ s/\s/_/g;
	s/<<([^>]+)>>/<A HREF="#$link_name">$1<\/A>/g;
    }
    while ( /([^:]+):(.*)/g ) {
	$link_name = $1;
	$link_name =~ s/\s/_/g;
	s/([^:]+):(.*)/<DT><A NAME="$link_name" CLASS="name">$1<\/A><\/DT>\n\t<DD>$2<\/DD>/;
    }
    print "\n$_\n";
}

print << 'EOT'
</DL>
</BODY>
EOT
    ;
