%{
#define YYSTYPE string
#define OUTW 25

#include <string>
#include <iomanip>
#include <stack>

string SPC(" ");
string class_name; // Contains the last classname parsed
stack<string> class_stack;
stack<string> access_stack;

int yylex();
void yy_push_fargs();

void yyerror(const char* s)
  {
    cerr << s << endl;
  }
%}

%token CLASS
%token ACCESS
%token IDENTIFIER
%token FARG
%token CONST
%token VIRTUAL
%token ENUM
%token FRIEND
%token OPERATOR
%token OP
%token TEMPLATE
%token UNSIGNED
%token INT
%token FLOAT
%token CHAR
%token STRING
%token COMTEMP
%token DECLEXCEPTION
%token INLINE
%%

all:  declaration_list
  | declaration_list INLINE { return 0; }
;

declaration_list: declaration
  | declaration_list { cout << "@@@" << endl; } declaration
  | declaration_list ';'
;

declaration: class_declaration ';'
  | function_declaration ';'
  | function_declaration '{' '}' ';'
  | function_declaration '{' '}'
  | variable_declaration ';'
  | enum_declaration ';'
  | template_declaration class_declaration ';'
    { cout << setw(OUTW) << "Template-Description:" << $1 << endl; }
  | template_declaration function_declaration ';'
    { cout << setw(OUTW) << "Template-Description:" << $1 << endl; }
  | deal_exception_declaration ';'
;

variable_declaration:
  vartype identifier
    { cout << setw(OUTW) << "Variable-Declaration:" << $2 << endl
	   << setw(OUTW) << "Type:" << $1 << endl
	   << setw(OUTW) << "Access:" << access_stack.top() << endl; }
  | vartype identifier array_dimensions
    { cout << setw(OUTW) << "Variable-Declaration:" << $2 << endl
	   << setw(OUTW) << "Type:" << $1 << endl
	   << setw(OUTW) << "Access:" << access_stack.top() << endl
	   << setw(OUTW) << "Array-Dimension:" << $3 << endl; }
;

array_dimensions:
  '[' ']' { $$ = string("[]"); }
  | '[' identifier ']' { $$ = string("[") + $2 + string("]"); }
  | '[' INT ']' { $$ = string("[") + $2 + string("]"); }
  | array_dimensions array_dimensions { $$ = $1 + $2; }
;

function_declaration:
  function_name argument_declaration
    {
      cout << setw(OUTW) << "Function-Definition:" << $1 << endl
	   << setw(OUTW) << "Function-Parameters:" << $2 << endl
	   << setw(OUTW) << "In-Class:" << class_stack.top() << endl
	   << setw(OUTW) << "Access:" << access_stack.top() << endl;
    }
  | function_name argument_declaration CONST
    {
      cout << setw(OUTW) << "Function-Definition:" << $1 << endl
	   << setw(OUTW) << "Function-Parameters:" << $2 << endl
	   << setw(OUTW) << "Const:" << "const" << endl
	   << setw(OUTW) << "In-Class:" << class_stack.top() << endl
	   << setw(OUTW) << "Access:" << access_stack.top() << endl;
    }
  | vartype function_declaration
    {
      cout << setw(OUTW) << "Return-Type:" << $1 << endl;
    }
;

function_name:
  identifier
  | '~' identifier { $$ = string("~") + $2; }
  | OPERATOR operator { $$ = string("operator") + $2; }
;

operator: OP
  | '='
  | '&'
  | '*'
  | '<'
  | '>'
  | '(' ')' { $$ = string("()"); }
  | '[' ']' { $$ = string("[]"); }
;

argument_declaration:
  '(' ')' { $$ = string(""); }
  | '(' arguments ')'
    { $$ = $2; }
;

arguments:
  argument
  | arguments ',' argument { $$ = $1 + string("@") + $3; }
;

argument:
  vartype
  | vartype identifier { $$ = $1 + SPC + $2; }
  | vartype identifier '=' default_arg
    { $$ = $1 + SPC + $2 + string(" = ") + $4; }
;

default_arg:
    literal
    | identifier
;

enum_declaration:
    ENUM IDENTIFIER '{' enum_list '}'
;

enum_list: /* empty */
    |  enum_list ',' enumerator
;

enumerator: IDENTIFIER
  | IDENTIFIER '=' literal
;

class_declaration: class_head
  | class_head
    {
      class_name = class_stack.top() + string("::") + class_name;
      class_stack.push(class_name);
      access_stack.push(string("private"));
    } '{' class_body '}'
    {
      class_stack.pop();
      access_stack.pop();
      cout << "@@@" << endl << setw(OUTW) << "Class-Definition:" << $1 << endl
	   << setw(OUTW) << "In-Class:" << class_stack.top() << endl;
    }
;

class_head:
  CLASS identifier { $$ = $1 + SPC + $2; class_name = $2;}
  | CLASS identifier inheritance
    { $$ = $1 + SPC + $2 + $3; class_name = $2;}
;

inheritance:
  ':' inheritance_list { $$ = $2; }
;

inheritance_list:
  ancestor
  | inheritance_list ',' ancestor { $$ = $1 + $3; }
;

ancestor:
  virtualopt ACCESS virtualopt identifier
  { $$ = string("@") + $1 + SPC + $2 + SPC + $3 + SPC + $4; }
;

class_body: /* empty */
  | member_declaration_list
;

member_declaration_list: member_declaration
  | member_declaration_list { cout << "@@@" << endl; } member_declaration 
  | member_declaration_list ';'
;

member_declaration:
  declaration
  | VIRTUAL declaration
    {
      cout << setw(OUTW) << "Virtual:" << "virtual" << endl;
    }
  | VIRTUAL function_declaration '=' INT
    {
      cout << setw(OUTW) << "Virtual:" << "pure" << endl;
    }
  | ACCESS ':' { access_stack.top() = $1; }
;

vartype:
  identifier { $$ = $1; }
  | UNSIGNED IDENTIFIER { $$ = string("unsigned ") + $2; }
  | CONST vartype { $$ = string("const ") + $2; }
/*  | vartype CONST { $$ = $1 +  string(" const"); } */
  | vartype '&' { $$ = $1 + string("&"); }
  | vartype '*' { $$ = $1 + string("*"); }
;

virtualopt: VIRTUAL { $$ = string("virtual"); }
  | /* empty */ { $$ = string(""); }
;

identifier:
  IDENTIFIER
  | identifier ':' ':' identifier { $$ = $1 + string("::") + $4; }
  | identifier template_args { $$ = $1 + $2; }
  | IDENTIFIER COMTEMP template_arglist '>'
    { $$ = $1 + string("<") + $3 + string(">"); }
;

expression:
  literal
  | identifier
  | OP expression { $$ = $1 + $2; }
  | expression OP expression { $$ = $1 + $2 + $3; }
  | expression INT { $$ = $1 + $2; }
;

literal:
  INT
  | FLOAT
  | CHAR
  | STRING
;

template_args:
  '<' template_arglist '>' { $$ = string("<") + $2 + string(">"); }
;

template_arglist:
  template_arg
  | template_arglist ',' template_arg { $$ = $1 + string(",") + $3; }
;

template_arg:
  vartype
  | expression
;

template_declaration:
  TEMPLATE '<' arguments '>' { $$ = $3; }
;

deal_exception_declaration:
  DECLEXCEPTION '(' IDENTIFIER deal_exception_arglist ')'
;

deal_exception_arglist: /* empty */
  | deal_exception_arglist ',' deal_exception_arg
;

deal_exception_arg: identifier
  | deal_exception_output_declaration
  | literal
;

deal_exception_output_declaration: OP
  | deal_exception_output_declaration IDENTIFIER
  | deal_exception_output_declaration literal
  | deal_exception_output_declaration OP
;

%%

main(int argc, const char** argv)
{
  string toplevel("");
  class_stack.push(toplevel);
  access_stack.push(toplevel);
  if (argc>1)
    yydebug=atoi(argv[1]);
  else
    yydebug=0;
  yyparse();
}
