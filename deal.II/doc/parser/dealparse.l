
%option stack

%{
#define YYSTYPE string
#include <string>
#include <iostream>

#include "dealparse.tab.h"

extern int yydebug;

%}

%x COMMENT
%x LINECOMMENT
%x EXC_IGNORE
%x INSTRING

NAME  [a-zA-Z][_a-zA-Z0-9]*
INT   [+-]?[0-9]+
WHITE [ \t\n\f\v]
CLASS (class)|(struct)|(union)
ACCESS (public)|(protected)|(private)
OP \+|\-|\/|\||\^|(==)|(!=)|(<<)|(>>)

%%


"//<"               { return COMTEMP; }

\n?"/*"             { yy_push_state(COMMENT); }
<COMMENT>.
<COMMENT>\n
<COMMENT>"*/"       { yy_pop_state(); }

DeclException[0-9]  { return(DECLEXCEPTION); }
<EXC_IGNORE>.
<EXC_IGNORE>\n
<EXC_IGNORE>");" { yy_pop_state(); }

"//"                { yy_push_state(LINECOMMENT); }
<LINECOMMENT>.
<LINECOMMENT>\n     { yy_pop_state(); }

 /* Ignore preprocessor commands */

^#.*\n

 /* Ignore friend declarations so far */

friend[^;]+; {}

 /*** keywords ***/

const    { return CONST; }
virtual  { return VIRTUAL; }
operator { return OPERATOR; }
unsigned { return UNSIGNED; }
template { return TEMPLATE; }
enum     { return ENUM; }
friend   { return FRIEND; }

{CLASS}  { yylval = string(yytext); return CLASS; }
{ACCESS} { yylval = string(yytext); return ACCESS; }

{WHITE} {}
{OP}     { yylval = string(yytext); if (yydebug)
             cerr << "[" << yylval << "]"; return OP; }

 /*** literals ***/

 /* Parse strings: does not account for embedded \" */

\"[^\"]*\" { yylval = string(yytext); return STRING; }

{INT}    { yylval = string(yytext); return INT; }
0[xX][0-9a-fA-F]+ { yylval = string(yytext); return INT; }
{INT}\.({INT})?(E[+-]?{INT})?  { yylval = string(yytext); return FLOAT; }
\.{INT}(E[+-]?{INT})?  { yylval = string(yytext); return FLOAT; }
{NAME}   { yylval = string(yytext); if (yydebug)
             cerr << "[" << yylval << "]"; return IDENTIFIER; }
. { yylval = string(yytext); return yytext[0]; }

%%

