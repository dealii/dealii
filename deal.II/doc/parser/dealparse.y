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
void enterfunction();
void enterbaseinitializers();  

string doctext;
string cvs_id;
  
void printdoc(ostream& s)
  {
    s << doctext << endl;
    doctext = string("");
  }
  
      
void yyerror(const char* s)
  {
    cerr << s << endl;
  }
%}

%token NAMESPACE
%token USING
%token CLASS
%token TYPEDEF
%token ACCESS
%token IDENTIFIER
%token FARG
%token CONST
%token STATIC
%token EXTERN
%token MUTABLE
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
%token THREEDOTS
%token ENDOFDECL
%%

all:  declaration_list
  | declaration_list INLINE { return 0; }
  | declaration_list ENDOFDECL { return 0; }
;

declaration_list: declaration
  | declaration_list { cout << "@@@" << endl; } declaration
  | declaration_list NAMESPACE identifier '{' declaration_list '}' ';'
  | declaration_list ';'
;

declaration:
  class_declaration ';'
  | function_declaration '{' { enterfunction(); }
  | function_declaration ':' { enterbaseinitializers(); }
  | function_declaration ';'
  | variable_declaration ';'
  | enum_declaration ';'
  | typedef_declaration ';'
  | template_declaration
    { cout << setw(OUTW) << "@Template-Description:" << $1 << endl; } class_declaration ';'
  | template_declaration function_declaration ';'
    { cout << setw(OUTW) << "@Template-Description:" << $1 << endl; }
  | template_declaration function_declaration '{' { enterfunction(); 
      cout << setw(OUTW) << "@Template-Description:" << $1 << endl; }
  | template_declaration function_declaration ':' { enterbaseinitializers(); 
      cout << setw(OUTW) << "@Template-Description:" << $1 << endl; }
  | using_declaration ';';
  | deal_exception_declaration ';'
  | ';'
;

using_declaration:
  USING NAMESPACE identifier
;

variable_declaration:
  vartype identifier
    { cout << setw(OUTW) << "@Variable-Declaration:" << $2 << endl
	   << cvs_id << endl
	   << setw(OUTW) << "@Type:" << $1 << endl
	   << setw(OUTW) << "@In-Class:" << class_stack.top() << endl
	   << setw(OUTW) << "@Access:" << access_stack.top() << endl;
      printdoc(cout); }
  | vartype identifier array_dimensions
    { cout << setw(OUTW) << "@Variable-Declaration:" << $2 << endl
	   << cvs_id << endl
	   << setw(OUTW) << "@Type:" << $1 << endl
	   << setw(OUTW) << "@In-Class:" << class_stack.top() << endl
	   << setw(OUTW) << "@Access:" << access_stack.top() << endl
	   << setw(OUTW) << "@Array-Dimension:" << $3 << endl;
      printdoc(cout); }
  | variable_declaration '=' expression
/*  | enum_declaration identifier*/
;

typedef_declaration:
    TYPEDEF vartype identifier
      { cout << setw(OUTW) << "@Typedef:" << $3 << endl
	   << cvs_id << endl
	     << setw(OUTW) << "@Type:" << $2 << endl
	     << setw(OUTW) << "@Access:" << access_stack.top() << endl;
	printdoc(cout); }
    | TYPEDEF  vartype '(' '*' identifier ')' argument_declaration
      { cout << setw(OUTW) << "@Typedef-Functionpointer:" << $4 << endl
	     << setw(OUTW) << "@Type:" << $2<< endl
	     << setw(OUTW) << "@Access:" << access_stack.top() << endl;
	printdoc(cout); }
    | TYPEDEF  vartype '(' identifier ':' ':' '*' identifier ')' argument_declaration
      { cout << setw(OUTW) << "@Typedef-Functionpointer:" << $8 << endl
	     << setw(OUTW) << "@Type:" << $2<< endl
	     << setw(OUTW) << "@Access:" << access_stack.top() << endl;
	printdoc(cout); }
    | TYPEDEF  vartype '(' '*' identifier ')' argument_declaration CONST
      { cout << setw(OUTW) << "@Typedef-Functionpointer:" << $4 << endl
	     << setw(OUTW) << "@Type:" << $2<< endl
	     << setw(OUTW) << "@Access:" << access_stack.top() << endl;
	printdoc(cout); }
    | TYPEDEF  vartype '(' identifier ':' ':' '*' identifier ')' argument_declaration CONST
      { cout << setw(OUTW) << "@Typedef-Functionpointer:" << $8 << endl
	     << setw(OUTW) << "@Type:" << $2<< endl
	     << setw(OUTW) << "@Access:" << access_stack.top() << endl;
	printdoc(cout); }
    | TYPEDEF vartype identifier array_dimensions
      { cout << setw(OUTW) << "@Typedef:" << $3 << endl
	     << cvs_id << endl
	     << setw(OUTW) << "@Type:" << $2 << endl
	     << setw(OUTW) << "@Access:" << access_stack.top() << endl
	     << setw(OUTW) << "@Array-Dimension:" << $4 << endl;
	printdoc(cout); }
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
      cout << setw(OUTW) << "@Function-Definition:" << $1 << endl
	   << cvs_id << endl
	   << setw(OUTW) << "@Function-Parameters:" << $2 << endl
	   << setw(OUTW) << "@In-Class:" << class_stack.top() << endl
	   << setw(OUTW) << "@Access:" << access_stack.top() << endl;
      printdoc(cout);
    }
  | function_name argument_declaration CONST
    {
      cout << setw(OUTW) << "@Function-Definition:" << $1 << endl
	   << cvs_id << endl
	   << setw(OUTW) << "@Function-Parameters:" << $2 << endl
	   << setw(OUTW) << "@Const:" << "const" << endl
	   << setw(OUTW) << "@In-Class:" << class_stack.top() << endl
	   << setw(OUTW) << "@Access:" << access_stack.top() << endl;
      printdoc(cout);
    }
  | vartype function_declaration
    {
      cout << setw(OUTW) << "@Return-Type:" << $1 << endl;
    }
;

function_name:
  identifier
  | '~' identifier { $$ = string("~") + $2; }
  | OPERATOR operator { $$ = string("operator") + $2; }
  | OPERATOR vartype { $$ = string("operator") + $2; }
  | identifier ':' ':' function_name { $$ = $1 + string("::") + $4; }
;

operator: OP
  | operator '=' { $$ = $1 + string("="); }
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
  | '(' THREEDOTS ')' { $$ = $2; }
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
  | vartype  '=' default_arg
    { $$ = $1 + SPC + string(" = ") + $3; }
  | CLASS IDENTIFIER { $$ = $1 + SPC + $2; }
  | vartype '(' identifier ')' argument_declaration
    { $$ = $1 + $2 + $3 + $4 + string(" (") + $5 + string(")"); }
;

default_arg:
  expression
;

enum_declaration:
    ENUM IDENTIFIER '{' enum_list '}'
;

enum_list: /* empty */
    | enumerator
    |  enum_list ',' enumerator { $$ = $1 + string("@") + $3; }
;

enumerator: IDENTIFIER
  | IDENTIFIER '=' expression
;

class_declaration: class_head
  | class_head
    {
      cout << setw(OUTW) << "@Class-Definition:" << $1 << endl
	   << cvs_id << endl
	   << setw(OUTW) << "@In-Class:" << class_stack.top() << endl;
      printdoc(cout);
      class_name = class_stack.top() + string("::") + class_name;
      class_stack.push(class_name);
      access_stack.push(string("private"));
    } '{' class_body '}'
    {
      class_stack.pop();
      access_stack.pop();
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
      cout << setw(OUTW) << "@Virtual:" << "virtual" << endl;
    }
  | VIRTUAL function_declaration '=' INT
    {
      cout << setw(OUTW) << "@Virtual:" << "pure" << endl;
    }
  | ACCESS ':' { access_stack.top() = $1; }
  | friend_declaration
;

friend_declaration:
  FRIEND /* declaration  { cout << setw(OUTW) << "Friend:" << endl; } */
  | template_declaration FRIEND /* declaration { cout << setw(OUTW) << "Friend:" << endl; } */
;

vartype:
  identifier { $$ = $1; }
  | UNSIGNED IDENTIFIER { $$ = string("unsigned ") + $2; }
  | CONST vartype { $$ = string("const ") + $2; }
  | vartype '*' CONST { $$ = $1 +  string("* const"); }
  | vartype '&' { $$ = $1 + string("&"); }
  | vartype '*' { $$ = $1 + string("*"); }
  | STATIC vartype { $$ = string("static ") + $2; }
  | MUTABLE vartype { $$ = string("mutable ") + $2; }
  | EXTERN vartype { $$ = $2; }

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
  | expression '*' expression { $$ = $1 + $2 + $3; }
  | expression INT { $$ = $1 + $2; }
  | identifier argument_declaration { $$ = $1 + $2; }
  | '(' expression ')' { $$ = $1 + $2 + $3; }
;

literal:
  INT
  | FLOAT
  | CHAR
  | STRING
;

template_args:
  '<' '>' { $$ = string("<>"); }
  | '<' template_arglist '>' { $$ = string("<") + $2 + string(">"); }
;

template_arglist:
  template_arg
  | template_arglist ',' template_arg { $$ = $1 + string(",") + $3; }
;

template_arg:
  vartype
  | expression
  | CLASS identifier { $$ = $1 + SPC + $2; }
;

template_statement:
  TEMPLATE '<' arguments '>' { $$ = $3; }
  | TEMPLATE '<'  '>' { $$ = SPC; }
;

template_declaration:
  template_statement
  | template_declaration template_statement { $$ = $1 + $2; }
;

deal_exception_declaration:
  DECLEXCEPTION '(' IDENTIFIER deal_exception_arglist ')'
  { cout << setw(OUTW) << $1 << ':' << $3 << endl
	 << setw(OUTW) << "@In-Class:" << class_stack.top() << endl
	 << setw(OUTW) << "@Access:" << access_stack.top() << endl;
    printdoc(cout); }
;

deal_exception_arglist: /* empty */
  | deal_exception_arglist ',' deal_exception_arg
;

deal_exception_arg:
  vartype
  | deal_exception_output_declaration
  | literal
;

deal_exception_output_declaration: OP
  | deal_exception_output_declaration vartype
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
