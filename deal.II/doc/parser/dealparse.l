
%option stack

%{
#define YYSTYPE string
#include <string>
#include <iostream>

#include "dealparse.tab.h"

extern int yydebug;
extern string doctext;
extern string cvs_id;

%}

%x COMMENT
%x LINECOMMENT
%x EXC_IGNORE
%x INSTRING
%x INFUNCTION
%x BASEINIT
%x DOCCOMMENT

NAME  [_a-zA-Z][_a-zA-Z0-9]*
INT   [+-]?[0-9]+
WHITE [ \t\n\f\v]
CLASS (class)|(struct)|(union)
ACCESS (public)|(protected)|(private)
OP \+|\-|\/|\||\^|(==)|(!=)|(<<)|(>>)|(->)|(\+\+)|(\-\-)

%%

 /* Commentary template arguments for the parser */

"//<"               { return COMTEMP; }

"/**"		{ yy_push_state(DOCCOMMENT); }
<DOCCOMMENT>\*\/       { yy_pop_state(); }
<DOCCOMMENT>.|\n	{ doctext += string(yytext); }

"/*"             { yy_push_state(COMMENT); }
<COMMENT>.
<COMMENT>\n
<COMMENT>"*/"       { yy_pop_state(); }

DeclException[0-9]  { yylval = string("@") + string(yytext); return(DECLEXCEPTION); }
<EXC_IGNORE>.
<EXC_IGNORE>\n
<EXC_IGNORE>");" { yy_pop_state(); }

"//"                { yy_push_state(LINECOMMENT); }
<LINECOMMENT>.
<LINECOMMENT>\n     { yy_pop_state(); }

 /* scan CVS id in file preamble */

"$Id$" { cvs_id = string(yytext); }
<LINECOMMENT>"$Id$" { cvs_id = string(yytext); }
<COMMENT>"$Id$" { cvs_id = string(yytext); }


 /* Ignore function bodies */

<INFUNCTION>"{" { yy_push_state(INFUNCTION); }
<INFUNCTION>"}" { yy_pop_state(); }
<INFUNCTION>.
<INFUNCTION>\n

 /* Ignore base initializers */
<BASEINIT>"{" { yy_pop_state(); yy_push_state(INFUNCTION); }
<BASEINIT>.
<BASEINIT>\n

 /* Ignore preprocessor commands */

^#.*(\\\n.*)*\n

 /* Ignore friend declarations so far */

friend[^;]+; { return FRIEND; }

 /*** keywords ***/

const     { return CONST; }
virtual   { return VIRTUAL; }
operator  { return OPERATOR; }
unsigned  { return UNSIGNED; }
static    { return STATIC; }
extern    { return EXTERN; }
mutable   { return MUTABLE; }
template  { return TEMPLATE; }
enum      { return ENUM; }
friend    { return FRIEND; }
typedef   { return TYPEDEF; }
namespace { return NAMESPACE; }
using     { return USING; }
typename {}

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

"..."    { yylval = string(yytext); return THREEDOTS; }
%%

void enterfunction()
{
  yy_push_state(INFUNCTION);
}

void enterbaseinitializers()
{
  yy_push_state(BASEINIT);
}
