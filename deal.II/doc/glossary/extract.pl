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
    s/<<([^>]+)>>/<A HREF="#$1">$1<\/A>/g;
    s/([^:]+):(.*)/<DT><A CLASS="name" NAME="$1">$1<\/A><\/DT>\n\t<DD>$2<\/DD>/;
    print "\n$_\n";
}

print << 'EOT'
</DL>
</BODY>
EOT
    ;
