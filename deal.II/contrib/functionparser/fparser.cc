//==============================
// Function parser v2.7 by Warp
//==============================

// Comment out the following line if your compiler supports the (non-standard)
// asinh, acosh and atanh functions and you want them to be supported. If
// you are not sure, just leave it (those function will then not be supported).
#define NO_ASINH


// Uncomment the following line to disable the eval() function if it could
// be too dangerous in the target application:
//#define DISABLE_EVAL


// Comment this line out if you are not going to use the optimizer and want
// a slightly smaller library. The Optimize() method can still be called,
// but it will not do anything.
// If you are unsure, just leave it. It won't slow down the other parts of
// the library.
#define SUPPORT_OPTIMIZER


//============================================================================

#include "fparser.h"

#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>

using namespace std;

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

namespace
{
// The functions must be in alphabetical order:
    enum OPCODE
    {
        cAbs, cAcos,
#ifndef NO_ASINH
        cAcosh,
#endif
        cAsin,
#ifndef NO_ASINH
        cAsinh,
#endif
        cAtan,
        cAtan2,
#ifndef NO_ASINH
        cAtanh,
#endif
        cCeil, cCos, cCosh, cCot, cCsc,
#ifndef DISABLE_EVAL
        cEval,
#endif
        cExp, cFloor, cIf, cInt, cLog, cLog10, cMax, cMin,
        cSec, cSin, cSinh, cSqrt, cTan, cTanh,

// These do not need any ordering:
        cImmed, cJump,
        cNeg, cAdd, cSub, cMul, cDiv, cMod, cPow,
        cEqual, cLess, cGreater, cAnd, cOr,

        cDeg, cRad,

        cFCall, cPCall,

#ifdef SUPPORT_OPTIMIZER
        cVar, cDup, cInv,
#endif

        VarBegin
    };

    struct FuncDefinition
    {
        const char* name;
        unsigned nameLength;
        unsigned opcode;
        unsigned params;

        // This is basically strcmp(), but taking 'nameLength' as string
        // length (not ending '\0'):
        bool operator<(const FuncDefinition& rhs) const
        {
            for(unsigned i = 0; i < nameLength; ++i)
            {
                if(i == rhs.nameLength) return false;
                const char c1 = name[i], c2 = rhs.name[i];
                if(c1 < c2) return true;
                if(c2 < c1) return false;
            }
            return nameLength < rhs.nameLength;
        }
    };


// This list must be in alphabetical order:
    const FuncDefinition Functions[]=
    {
        { "abs", 3, cAbs, 1 },
        { "acos", 4, cAcos, 1 },
#ifndef NO_ASINH
        { "acosh", 5, cAcosh, 1 },
#endif
        { "asin", 4, cAsin, 1 },
#ifndef NO_ASINH
        { "asinh", 5, cAsinh, 1 },
#endif
        { "atan", 4, cAtan, 1 },
        { "atan2", 5, cAtan2, 2 },
#ifndef NO_ASINH
        { "atanh", 5, cAtanh, 1 },
#endif
        { "ceil", 4, cCeil, 1 },
        { "cos", 3, cCos, 1 },
        { "cosh", 4, cCosh, 1 },
        { "cot", 3, cCot, 1 },
        { "csc", 3, cCsc, 1 },
#ifndef DISABLE_EVAL
        { "eval", 4, cEval, 0 },
#endif
        { "exp", 3, cExp, 1 },
        { "floor", 5, cFloor, 1 },
        { "if", 2, cIf, 0 },
        { "int", 3, cInt, 1 },
        { "log", 3, cLog, 1 },
        { "log10", 5, cLog10, 1 },
        { "max", 3, cMax, 2 },
        { "min", 3, cMin, 2 },
        { "sec", 3, cSec, 1 },
        { "sin", 3, cSin, 1 },
        { "sinh", 4, cSinh, 1 },
        { "sqrt", 4, cSqrt, 1 },
        { "tan", 3, cTan, 1 },
        { "tanh", 4, cTanh, 1 }
    };

    const unsigned FUNC_AMOUNT = sizeof(Functions)/sizeof(Functions[0]);


    // BCB4 does not implement the standard lower_bound function.
    // This is used instead:
    const FuncDefinition* fp_lower_bound(const FuncDefinition* first,
                                         const FuncDefinition* last,
                                         const FuncDefinition& value)
    {
        while(first < last)
        {
            const FuncDefinition* middle = first+(last-first)/2;
            if(*middle < value) first = middle+1;
            else last = middle;
        }
        return last;
    }


    // Returns a pointer to the FuncDefinition instance which 'name' is
    // the same as the one given by 'F'. If no such function name exists,
    // returns 0.
    inline const FuncDefinition* FindFunction(const char* F)
    {
        FuncDefinition func = { F, 0, 0, 0 };
        while(isalnum(F[func.nameLength])) ++func.nameLength;
        if(func.nameLength)
        {
            const FuncDefinition* found =
                fp_lower_bound(Functions, Functions+FUNC_AMOUNT, func);
            if(found == Functions+FUNC_AMOUNT || func < *found)
                return 0;
            return found;
        }
        return 0;
    }
}


//---------------------------------------------------------------------------
// Copy-on-write method
//---------------------------------------------------------------------------
inline void FunctionParser::copyOnWrite()
{
    if(data->referenceCounter > 1)
    {
        Data* oldData = data;
        data = new Data(*oldData);
        --(oldData->referenceCounter);
        data->referenceCounter = 1;
    }
}


//---------------------------------------------------------------------------
// Constructors and destructors
//---------------------------------------------------------------------------
//===========================================================================
FunctionParser::FunctionParser():
    parseErrorType(FP_NO_ERROR), evalErrorType(0),
    data(new Data)
{
    data->referenceCounter = 1;
}

FunctionParser::~FunctionParser()
{
    if(--(data->referenceCounter) == 0)
    {
        delete data;
    }
}

FunctionParser::FunctionParser(const FunctionParser& cpy):
    parseErrorType(cpy.parseErrorType),
    evalErrorType(cpy.evalErrorType),
    data(cpy.data)
{
    ++(data->referenceCounter);
}

FunctionParser& FunctionParser::operator=(const FunctionParser& cpy)
{
    if(data != cpy.data)
    {
        if(--(data->referenceCounter) == 0) delete data;

        parseErrorType = cpy.parseErrorType;
        evalErrorType = cpy.evalErrorType;
        data = cpy.data;

        ++(data->referenceCounter);
    }

    return *this;
}


FunctionParser::Data::Data():
    useDegreeConversion(false),
    ByteCode(0), ByteCodeSize(0),
    Immed(0), ImmedSize(0),
    Stack(0), StackSize(0)
{}

FunctionParser::Data::~Data()
{
    if(ByteCode) { delete[] ByteCode; ByteCode=0; }
    if(Immed) { delete[] Immed; Immed=0; }
    if(Stack) { delete[] Stack; Stack=0; }
}

// Makes a deep-copy of Data:
FunctionParser::Data::Data(const Data& cpy):
    varAmount(cpy.varAmount), useDegreeConversion(cpy.useDegreeConversion),
    Variables(cpy.Variables), Constants(cpy.Constants),
    FuncPtrNames(cpy.FuncPtrNames), FuncPtrs(cpy.FuncPtrs),
    FuncParserNames(cpy.FuncParserNames), FuncParsers(cpy.FuncParsers),
    ByteCode(0), ByteCodeSize(cpy.ByteCodeSize),
    Immed(0), ImmedSize(cpy.ImmedSize),
    Stack(0), StackSize(cpy.StackSize)
{
    if(ByteCodeSize) ByteCode = new unsigned[ByteCodeSize];
    if(ImmedSize) Immed = new double[ImmedSize];
    if(StackSize) Stack = new double[StackSize];

    for(unsigned i=0; i<ByteCodeSize; ++i) ByteCode[i] = cpy.ByteCode[i];
    for(unsigned i=0; i<ImmedSize; ++i) Immed[i] = cpy.Immed[i];

    // No need to copy the stack contents because it's obsolete outside Eval()
}


//---------------------------------------------------------------------------
// Function parsing
//---------------------------------------------------------------------------
//===========================================================================
namespace
{
    // Error messages returned by ErrorMsg():
    const char* ParseErrorMessage[]=
    {
        "Syntax error",                             // 0
        "Mismatched parenthesis",                   // 1
        "Missing ')'",                              // 2
        "Empty parentheses",                        // 3
        "Syntax error: Operator expected",          // 4
        "Not enough memory",                        // 5
        "An unexpected error ocurred. Please make a full bug report "
        "to warp@iki.fi",                           // 6
        "Syntax error in parameter 'Vars' given to "
        "FunctionParser::Parse()",                  // 7
        "Illegal number of parameters to function", // 8
        "Syntax error: Premature end of string",    // 9
        "Syntax error: Expecting ( after function", // 10
        ""
    };


    // Parse variables
    bool ParseVars(const string& Vars, map<string, unsigned>& dest)
    {
        unsigned varNumber = VarBegin;
        unsigned ind1 = 0, ind2;

        while(ind1 < Vars.size())
        {
            if(!isalpha(Vars[ind1]) && Vars[ind1]!='_') return false;
            for(ind2=ind1+1; ind2<Vars.size() && Vars[ind2]!=','; ++ind2)
                if(!isalnum(Vars[ind2]) && Vars[ind2]!='_') return false;
            const string varName = Vars.substr(ind1, ind2-ind1);

            if(dest.insert(make_pair(varName, varNumber++)).second == false)
                return false;

            ind1 = ind2+1;
        }
        return true;
    }
}

bool FunctionParser::isValidName(const std::string& name) const
{
    if(name.empty() || (!isalpha(name[0]) && name[0] != '_')) return false;
    for(unsigned i=0; i<name.size(); ++i)
        if(!isalnum(name[i]) && name[i] != '_') return false;

    if(FindFunction(name.c_str())) return false;

    return true;
}


// Constants:
bool FunctionParser::AddConstant(const string& name, double value)
{
    if(isValidName(name))
    {
        const char* n = name.c_str();
        if(FindVariable(n, data->FuncParserNames) !=
           data->FuncParserNames.end() ||
           FindVariable(n, data->FuncPtrNames) !=
           data->FuncPtrNames.end())
            return false;

        copyOnWrite();

        data->Constants[name] = value;
        return true;
    }
    return false;
}

// Function pointers
bool FunctionParser::AddFunction(const std::string& name,
                                 FunctionPtr func, unsigned paramsAmount)
{
    if(paramsAmount == 0) return false; // Currently must be at least one

    if(isValidName(name))
    {
        const char* n = name.c_str();
        if(FindVariable(n, data->FuncParserNames) !=
           data->FuncParserNames.end() ||
           FindConstant(n) != data->Constants.end())
            return false;

        copyOnWrite();

        data->FuncPtrNames[name] = data->FuncPtrs.size();
        data->FuncPtrs.push_back(Data::FuncPtrData(func, paramsAmount));
        return true;
    }
    return false;
}

bool FunctionParser::checkRecursiveLinking(const FunctionParser* fp) const
{
    if(fp == this) return true;
    for(unsigned i=0; i<fp->data->FuncParsers.size(); ++i)
        if(checkRecursiveLinking(fp->data->FuncParsers[i])) return true;
    return false;
}

bool FunctionParser::AddFunction(const std::string& name,
                                 FunctionParser& parser)
{
    if(parser.data->varAmount == 0) // Currently must be at least one
        return false;

    if(isValidName(name))
    {
        const char* n = name.c_str();
        if(FindVariable(n, data->FuncPtrNames) != data->FuncPtrNames.end() ||
           FindConstant(n) != data->Constants.end())
            return false;

        if(checkRecursiveLinking(&parser)) return false;

        copyOnWrite();

        data->FuncParserNames[name] = data->FuncParsers.size();
        data->FuncParsers.push_back(&parser);
        return true;
    }
    return false;
}



// Main parsing function
// ---------------------
int FunctionParser::Parse(const std::string& Function,
                          const std::string& Vars,
                          bool useDegrees)
{
    copyOnWrite();

    data->Variables.clear();

    if(!ParseVars(Vars, data->Variables))
    {
        parseErrorType = INVALID_VARS;
        return Function.size();
    }
    data->varAmount = data->Variables.size(); // this is for Eval()

    const char* Func = Function.c_str();

    parseErrorType = FP_NO_ERROR;

    int Result = CheckSyntax(Func);
    if(Result>=0) return Result;

    data->useDegreeConversion = useDegrees;
    if(!Compile(Func)) return Function.size();

    data->Variables.clear();

    parseErrorType = FP_NO_ERROR;
    return -1;
}

namespace
{
    // Is given char an operator?
    inline bool IsOperator(int c)
    {
        return strchr("+-*/%^=<>&|",c)!=NULL;
    }

    // skip whitespace
    inline void sws(const char* F, int& Ind)
    {
        while(F[Ind] && isspace(F[Ind])) ++Ind;
    }
}

// Returns an iterator to the variable with the same name as 'F', or to
// Variables.end() if no such variable exists:
inline FunctionParser::Data::VarMap_t::const_iterator
FunctionParser::FindVariable(const char* F, const Data::VarMap_t& vars) const
{
    if(vars.size())
    {
        unsigned ind = 0;
        while(isalnum(F[ind]) || F[ind] == '_') ++ind;
        if(ind)
        {
            string name(F, ind);
            return vars.find(name);
        }
    }
    return vars.end();
}

inline FunctionParser::Data::ConstMap_t::const_iterator
FunctionParser::FindConstant(const char* F) const
{
    if(data->Constants.size())
    {
        unsigned ind = 0;
        while(isalnum(F[ind]) || F[ind] == '_') ++ind;
        if(ind)
        {
            string name(F, ind);
            return data->Constants.find(name);
        }
    }
    return data->Constants.end();
}

//---------------------------------------------------------------------------
// Check function string syntax
// ----------------------------
int FunctionParser::CheckSyntax(const char* Function)
{
    const Data::VarMap_t& Variables = data->Variables;
    const Data::ConstMap_t& Constants = data->Constants;
    const Data::VarMap_t& FuncPtrNames = data->FuncPtrNames;
    const Data::VarMap_t& FuncParserNames = data->FuncParserNames;

    vector<int> functionParenthDepth;

    int Ind=0, ParenthCnt=0, c;
    char* Ptr;

    while(true)
    {
        sws(Function, Ind);
        c=Function[Ind];

// Check for valid operand (must appear)

        // Check for leading -
        if(c=='-') { sws(Function, ++Ind); c=Function[Ind]; }
        if(c==0) { parseErrorType=PREMATURE_EOS; return Ind; }

        // Check for math function
        bool foundFunc = false;
        const FuncDefinition* fptr = FindFunction(&Function[Ind]);
        if(fptr)
        {
            Ind += fptr->nameLength;
            foundFunc = true;
        }
        else
        {
            // Check for user-defined function
            Data::VarMap_t::const_iterator fIter =
                FindVariable(&Function[Ind], FuncPtrNames);
            if(fIter != FuncPtrNames.end())
            {
                Ind += fIter->first.size();
                foundFunc = true;
            }
            else
            {
                Data::VarMap_t::const_iterator pIter =
                    FindVariable(&Function[Ind], FuncParserNames);
                if(pIter != FuncParserNames.end())
                {
                    Ind += pIter->first.size();
                    foundFunc = true;
                }
            }
        }

        if(foundFunc)
        {
            sws(Function, Ind);
            c = Function[Ind];
            if(c!='(') { parseErrorType=EXPECT_PARENTH_FUNC; return Ind; }
            functionParenthDepth.push_back(ParenthCnt+1);
        }

        // Check for opening parenthesis
        if(c=='(')
        {
            ++ParenthCnt;
            sws(Function, ++Ind);
            if(Function[Ind]==')') { parseErrorType=EMPTY_PARENTH; return Ind;}
            continue;
        }

        // Check for number
        if(isdigit(c) || (c=='.' && isdigit(Function[Ind+1])))
        {
            strtod(&Function[Ind], &Ptr);
            Ind += int(Ptr-&Function[Ind]);
            sws(Function, Ind);
            c = Function[Ind];
        }
        else
        { // Check for variable
            Data::VarMap_t::const_iterator vIter =
                FindVariable(&Function[Ind], Variables);
            if(vIter != Variables.end())
                Ind += vIter->first.size();
            else
            {
                // Check for constant
                Data::ConstMap_t::const_iterator cIter =
                    FindConstant(&Function[Ind]);
                if(cIter != Constants.end())
                    Ind += cIter->first.size();
                else
                { parseErrorType=SYNTAX_ERROR; return Ind; }
            }
            sws(Function, Ind);
            c = Function[Ind];
        }

        // Check for closing parenthesis
        while(c==')')
        {
            if(functionParenthDepth.size() &&
               functionParenthDepth.back() == ParenthCnt)
                functionParenthDepth.pop_back();
            if((--ParenthCnt)<0) { parseErrorType=MISM_PARENTH; return Ind; }
            sws(Function, ++Ind);
            c=Function[Ind];
        }

// If we get here, we have a legal operand and now a legal operator or
// end of string must follow

        // Check for EOS
        if(c==0) break; // The only way to end the checking loop without error
        // Check for operator
        if(!IsOperator(c) &&
           (c != ',' || functionParenthDepth.empty() ||
            functionParenthDepth.back() != ParenthCnt))
        { parseErrorType=EXPECT_OPERATOR; return Ind; }

// If we get here, we have an operand and an operator; the next loop will
// check for another operand (must appear)
        ++Ind;
    } // while

    // Check that all opened parentheses are also closed
    if(ParenthCnt>0) { parseErrorType=MISSING_PARENTH; return Ind; }

// The string is ok
    parseErrorType=FP_NO_ERROR;
    return -1;
}


// Compile function string to bytecode
// -----------------------------------
bool FunctionParser::Compile(const char* Function)
{
    if(data->ByteCode) { delete[] data->ByteCode; data->ByteCode=0; }
    if(data->Immed) { delete[] data->Immed; data->Immed=0; }
    if(data->Stack) { delete[] data->Stack; data->Stack=0; }

    vector<unsigned> byteCode; byteCode.reserve(1024);
    tempByteCode = &byteCode;

    vector<double> immed; immed.reserve(1024);
    tempImmed = &immed;

    data->StackSize = StackPtr = 0;

    CompileExpression(Function, 0);
    if(parseErrorType != FP_NO_ERROR) return false;

    data->ByteCodeSize = byteCode.size();
    data->ImmedSize = immed.size();

    if(data->ByteCodeSize)
    {
        data->ByteCode = new unsigned[data->ByteCodeSize];
        memcpy(data->ByteCode, &byteCode[0],
               sizeof(unsigned)*data->ByteCodeSize);
    }
    if(data->ImmedSize)
    {
        data->Immed = new double[data->ImmedSize];
        memcpy(data->Immed, &immed[0],
               sizeof(double)*data->ImmedSize);
    }
    if(data->StackSize)
        data->Stack = new double[data->StackSize];

    return true;
}


inline void FunctionParser::AddCompiledByte(unsigned c)
{
    tempByteCode->push_back(c);
}

inline void FunctionParser::AddImmediate(double i)
{
    tempImmed->push_back(i);
}

inline void FunctionParser::AddFunctionOpcode(unsigned opcode)
{
    if(data->useDegreeConversion)
        switch(opcode)
        {
          case cCos:
          case cCosh:
          case cCot:
          case cCsc:
          case cSec:
          case cSin:
          case cSinh:
          case cTan:
          case cTanh:
              AddCompiledByte(cRad);
        }

    AddCompiledByte(opcode);

    if(data->useDegreeConversion)
        switch(opcode)
        {
          case cAcos:
#ifndef NO_ASINH
          case cAcosh:
          case cAsinh:
          case cAtanh:
#endif
          case cAsin:
          case cAtan:
          case cAtan2:
              AddCompiledByte(cDeg);
        }
}

inline void FunctionParser::incStackPtr()
{
    if(++StackPtr > data->StackSize) ++(data->StackSize);
}


// Compile if()
int FunctionParser::CompileIf(const char* F, int ind)
{
    int ind2 = CompileExpression(F, ind, true); // condition
    sws(F, ind2);
    if(F[ind2] != ',') { parseErrorType=ILL_PARAMS_AMOUNT; return ind2; }
    AddCompiledByte(cIf);
    unsigned curByteCodeSize = tempByteCode->size();
    AddCompiledByte(0); // Jump index; to be set later
    AddCompiledByte(0); // Immed jump index; to be set later

    --StackPtr;

    ind2 = CompileExpression(F, ind2+1, true); // then
    sws(F, ind2);
    if(F[ind2] != ',') { parseErrorType=ILL_PARAMS_AMOUNT; return ind2; }
    AddCompiledByte(cJump);
    unsigned curByteCodeSize2 = tempByteCode->size();
    unsigned curImmedSize2 = tempImmed->size();
    AddCompiledByte(0); // Jump index; to be set later
    AddCompiledByte(0); // Immed jump index; to be set later

    --StackPtr;

    ind2 = CompileExpression(F, ind2+1, true); // else
    sws(F, ind2);
    if(F[ind2] != ')') { parseErrorType=ILL_PARAMS_AMOUNT; return ind2; }

    // Set jump indices
    (*tempByteCode)[curByteCodeSize] = curByteCodeSize2+1;
    (*tempByteCode)[curByteCodeSize+1] = curImmedSize2;
    (*tempByteCode)[curByteCodeSize2] = tempByteCode->size()-1;
    (*tempByteCode)[curByteCodeSize2+1] = tempImmed->size();

    return ind2+1;
}

int FunctionParser::CompileFunctionParams(const char* F, int ind,
                                          unsigned requiredParams)
{
    unsigned curStackPtr = StackPtr;
    int ind2 = CompileExpression(F, ind);

    if(StackPtr != curStackPtr+requiredParams)
    { parseErrorType=ILL_PARAMS_AMOUNT; return ind; }

    StackPtr -= requiredParams - 1;
    sws(F, ind2);
    return ind2+1; // F[ind2] is ')'
}

// Compiles element
int FunctionParser::CompileElement(const char* F, int ind)
{
    sws(F, ind);
    char c = F[ind];

    if(c == '(')
    {
        ind = CompileExpression(F, ind+1);
        sws(F, ind);
        return ind+1; // F[ind] is ')'
    }

    if(isdigit(c) || c=='.' /*|| c=='-'*/) // Number
    {
        const char* startPtr = &F[ind];
        char* endPtr;
        double val = strtod(startPtr, &endPtr);
        AddImmediate(val);
        AddCompiledByte(cImmed);
        incStackPtr();
        return ind+(endPtr-startPtr);
    }

    if(isalpha(c) || c == '_') // Function, variable or constant
    {
        const FuncDefinition* func = FindFunction(F+ind);
        if(func) // is function
        {
            int ind2 = ind + func->nameLength;
            sws(F, ind2); // F[ind2] is '('
            if(strcmp(func->name, "if") == 0) // "if" is a special case
            {
                return CompileIf(F, ind2+1);
            }

#ifndef DISABLE_EVAL
            unsigned requiredParams =
                strcmp(func->name, "eval") == 0 ?
                data->Variables.size() : func->params;
#else
            unsigned requiredParams = func->params;
#endif
            ind2 = CompileFunctionParams(F, ind2+1, requiredParams);
            AddFunctionOpcode(func->opcode);
            return ind2; // F[ind2-1] is ')'
        }

        Data::VarMap_t::const_iterator vIter =
            FindVariable(F+ind, data->Variables);
        if(vIter != data->Variables.end()) // is variable
        {
            AddCompiledByte(vIter->second);
            incStackPtr();
            return ind + vIter->first.size();
        }

        Data::ConstMap_t::const_iterator cIter = FindConstant(F+ind);
        if(cIter != data->Constants.end()) // is constant
        {
            AddImmediate(cIter->second);
            AddCompiledByte(cImmed);
            incStackPtr();
            return ind + cIter->first.size();
        }

        Data::VarMap_t::const_iterator fIter =
            FindVariable(F+ind, data->FuncPtrNames);
        if(fIter != data->FuncPtrNames.end()) // is user-defined func pointer
        {
            unsigned index = fIter->second;

            int ind2 = ind + fIter->first.length();
            sws(F, ind2); // F[ind2] is '('

            ind2 = CompileFunctionParams(F, ind2+1,
                                         data->FuncPtrs[index].params);

            AddCompiledByte(cFCall);
            AddCompiledByte(index);
            return ind2;
        }

        Data::VarMap_t::const_iterator pIter =
            FindVariable(F+ind, data->FuncParserNames);
        if(pIter != data->FuncParserNames.end()) // is user-defined func parser
        {
            unsigned index = pIter->second;

            int ind2 = ind + pIter->first.length();
            sws(F, ind2); // F[ind2] is '('

            ind2 = CompileFunctionParams
                (F, ind2+1, data->FuncParsers[index]->data->varAmount);

            AddCompiledByte(cPCall);
            AddCompiledByte(index);
            return ind2;
        }
    }

    parseErrorType = UNEXPECTED_ERROR;
    return ind;
}

// Compiles '^'
int FunctionParser::CompilePow(const char* F, int ind)
{
    int ind2 = CompileElement(F, ind);
    sws(F, ind2);

    while(F[ind2] == '^')
    {
        ind2 = CompileUnaryMinus(F, ind2+1);
        sws(F, ind2);
        AddCompiledByte(cPow);
        --StackPtr;
    }

    return ind2;
}

// Compiles unary '-'
int FunctionParser::CompileUnaryMinus(const char* F, int ind)
{
    sws(F, ind);
    if(F[ind] == '-')
    {
        int ind2 = ind+1;
        sws(F, ind2);
        ind2 = CompilePow(F, ind2);
        sws(F, ind2);

        // if we are negating a constant, negate the constant itself:
        if(tempByteCode->back() == cImmed)
            tempImmed->back() = -tempImmed->back();

        // if we are negating a negation, we can remove both:
        else if(tempByteCode->back() == cNeg)
            tempByteCode->pop_back();

        else
            AddCompiledByte(cNeg);

        return ind2;
    }

    int ind2 = CompilePow(F, ind);
    sws(F, ind2);
    return ind2;
}

// Compiles '*', '/' and '%'
int FunctionParser::CompileMult(const char* F, int ind)
{
    int ind2 = CompileUnaryMinus(F, ind);
    sws(F, ind2);
    char op;

    while((op = F[ind2]) == '*' || op == '/' || op == '%')
    {
        ind2 = CompileUnaryMinus(F, ind2+1);
        sws(F, ind2);
        switch(op)
        {
          case '*': AddCompiledByte(cMul); break;
          case '/': AddCompiledByte(cDiv); break;
          case '%': AddCompiledByte(cMod); break;
        }
        --StackPtr;
    }

    return ind2;
}

// Compiles '+' and '-'
int FunctionParser::CompileAddition(const char* F, int ind)
{
    int ind2 = CompileMult(F, ind);
    sws(F, ind2);
    char op;

    while((op = F[ind2]) == '+' || op == '-')
    {
        ind2 = CompileMult(F, ind2+1);
        sws(F, ind2);
        AddCompiledByte(op=='+' ? cAdd : cSub);
        --StackPtr;
    }

    return ind2;
}

// Compiles '=', '<' and '>'
int FunctionParser::CompileComparison(const char* F, int ind)
{
    int ind2 = CompileAddition(F, ind);
    sws(F, ind2);
    char op;

    while((op = F[ind2]) == '=' || op == '<' || op == '>')
    {
        ind2 = CompileAddition(F, ind2+1);
        sws(F, ind2);
        switch(op)
        {
          case '=': AddCompiledByte(cEqual); break;
          case '<': AddCompiledByte(cLess); break;
          case '>': AddCompiledByte(cGreater); break;
        }
        --StackPtr;
    }

    return ind2;
}

// Compiles '&'
int FunctionParser::CompileAnd(const char* F, int ind)
{
    int ind2 = CompileComparison(F, ind);
    sws(F, ind2);

    while(F[ind2] == '&')
    {
        ind2 = CompileComparison(F, ind2+1);
        sws(F, ind2);
        AddCompiledByte(cAnd);
        --StackPtr;
    }

    return ind2;
}

// Compiles '|'
int FunctionParser::CompileOr(const char* F, int ind)
{
    int ind2 = CompileAnd(F, ind);
    sws(F, ind2);

    while(F[ind2] == '|')
    {
        ind2 = CompileAnd(F, ind2+1);
        sws(F, ind2);
        AddCompiledByte(cOr);
        --StackPtr;
    }

    return ind2;
}

// Compiles ','
int FunctionParser::CompileExpression(const char* F, int ind, bool stopAtComma)
{
    int ind2 = CompileOr(F, ind);
    sws(F, ind2);

    if(stopAtComma) return ind2;

    while(F[ind2] == ',')
    {
        ind2 = CompileOr(F, ind2+1);
        sws(F, ind2);
    }

    return ind2;
}


// Return parse error message
// --------------------------
const char* FunctionParser::ErrorMsg() const
{
    if(parseErrorType != FP_NO_ERROR) return ParseErrorMessage[parseErrorType];
    return 0;
}

//---------------------------------------------------------------------------
// Function evaluation
//---------------------------------------------------------------------------
//===========================================================================
namespace
{
    inline int doubleToInt(double d)
    {
        return d<0 ? -int((-d)+.5) : int(d+.5);
    }

    inline double Min(double d1, double d2)
    {
        return d1<d2 ? d1 : d2;
    }
    inline double Max(double d1, double d2)
    {
        return d1>d2 ? d1 : d2;
    }


    inline double DegreesToRadians(double degrees)
    {
        return degrees*(M_PI/180.0);
    }
    inline double RadiansToDegrees(double radians)
    {
        return radians*(180.0/M_PI);
    }
}

double FunctionParser::Eval(const double* Vars)
{
    const unsigned* const ByteCode = data->ByteCode;
    const double* const Immed = data->Immed;
    double* const Stack = data->Stack;
    const unsigned ByteCodeSize = data->ByteCodeSize;
    unsigned IP, DP=0;
    int SP=-1;

    for(IP=0; IP<ByteCodeSize; ++IP)
    {
        switch(ByteCode[IP])
        {
// Functions:
          case   cAbs: Stack[SP] = fabs(Stack[SP]); break;
          case  cAcos: if(Stack[SP] < -1 || Stack[SP] > 1)
                       { evalErrorType=4; return 0; }
                       Stack[SP] = acos(Stack[SP]); break;
#ifndef NO_ASINH
          case cAcosh: Stack[SP] = acosh(Stack[SP]); break;
#endif
          case  cAsin: if(Stack[SP] < -1 || Stack[SP] > 1)
                       { evalErrorType=4; return 0; }
                       Stack[SP] = asin(Stack[SP]); break;
#ifndef NO_ASINH
          case cAsinh: Stack[SP] = asinh(Stack[SP]); break;
#endif
          case  cAtan: Stack[SP] = atan(Stack[SP]); break;
          case cAtan2: Stack[SP-1] = atan2(Stack[SP-1], Stack[SP]);
                       --SP; break;
#ifndef NO_ASINH
          case cAtanh: Stack[SP] = atanh(Stack[SP]); break;
#endif
          case  cCeil: Stack[SP] = ceil(Stack[SP]); break;
          case   cCos: Stack[SP] = cos(Stack[SP]); break;
          case  cCosh: Stack[SP] = cosh(Stack[SP]); break;

          case   cCot:
              {
                  double t = tan(Stack[SP]);
                  if(t == 0) { evalErrorType=1; return 0; }
                  Stack[SP] = 1/t; break;
              }
          case   cCsc:
              {
                  double s = sin(Stack[SP]);
                  if(s == 0) { evalErrorType=1; return 0; }
                  Stack[SP] = 1/s; break;
              }


#ifndef DISABLE_EVAL
          case  cEval:
              {
                  data->Stack = new double[data->StackSize];
                  double retVal = Eval(&Stack[SP-data->varAmount+1]);
                  delete[] data->Stack;
                  data->Stack = Stack;
                  SP -= data->varAmount-1;
                  Stack[SP] = retVal;
                  break;
              }
#endif

          case   cExp: Stack[SP] = exp(Stack[SP]); break;
          case cFloor: Stack[SP] = floor(Stack[SP]); break;

          case    cIf:
              {
                  unsigned jumpAddr = ByteCode[++IP];
                  unsigned immedAddr = ByteCode[++IP];
                  if(doubleToInt(Stack[SP]) == 0)
                  {
                      IP = jumpAddr;
                      DP = immedAddr;
                  }
                  --SP; break;
              }

          case   cInt: Stack[SP] = floor(Stack[SP]+.5); break;
          case   cLog: if(Stack[SP] <= 0) { evalErrorType=3; return 0; }
                       Stack[SP] = log(Stack[SP]); break;
          case cLog10: if(Stack[SP] <= 0) { evalErrorType=3; return 0; }
                       Stack[SP] = log10(Stack[SP]); break;
          case   cMax: Stack[SP-1] = Max(Stack[SP-1], Stack[SP]);
                       --SP; break;
          case   cMin: Stack[SP-1] = Min(Stack[SP-1], Stack[SP]);
                       --SP; break;
          case   cSec:
              {
                  double c = cos(Stack[SP]);
                  if(c == 0) { evalErrorType=1; return 0; }
                  Stack[SP] = 1/c; break;
              }
          case   cSin: Stack[SP] = sin(Stack[SP]); break;
          case  cSinh: Stack[SP] = sinh(Stack[SP]); break;
          case  cSqrt: if(Stack[SP] < 0) { evalErrorType=2; return 0; }
                       Stack[SP] = sqrt(Stack[SP]); break;
          case   cTan: Stack[SP] = tan(Stack[SP]); break;
          case  cTanh: Stack[SP] = tanh(Stack[SP]); break;


// Misc:
          case cImmed: Stack[++SP] = Immed[DP++]; break;
          case  cJump: DP = ByteCode[IP+2];
                       IP = ByteCode[IP+1];
                       break;

// Operators:
          case   cNeg: Stack[SP] = -Stack[SP]; break;
          case   cAdd: Stack[SP-1] += Stack[SP]; --SP; break;
          case   cSub: Stack[SP-1] -= Stack[SP]; --SP; break;
          case   cMul: Stack[SP-1] *= Stack[SP]; --SP; break;
          case   cDiv: if(Stack[SP] == 0) { evalErrorType=1; return 0; }
                       Stack[SP-1] /= Stack[SP]; --SP; break;
          case   cMod: if(Stack[SP] == 0) { evalErrorType=1; return 0; }
                       Stack[SP-1] = fmod(Stack[SP-1], Stack[SP]);
                       --SP; break;
          case   cPow: Stack[SP-1] = pow(Stack[SP-1], Stack[SP]);
                       --SP; break;

          case cEqual: Stack[SP-1] = (Stack[SP-1] == Stack[SP]);
                       --SP; break;
          case  cLess: Stack[SP-1] = (Stack[SP-1] < Stack[SP]);
                       --SP; break;
          case cGreater: Stack[SP-1] = (Stack[SP-1] > Stack[SP]);
                         --SP; break;
          case   cAnd: Stack[SP-1] =
                           (doubleToInt(Stack[SP-1]) &&
                            doubleToInt(Stack[SP]));
                       --SP; break;
          case    cOr: Stack[SP-1] =
                           (doubleToInt(Stack[SP-1]) ||
                            doubleToInt(Stack[SP]));
                       --SP; break;

// Degrees-radians conversion:
          case   cDeg: Stack[SP] = RadiansToDegrees(Stack[SP]); break;
          case   cRad: Stack[SP] = DegreesToRadians(Stack[SP]); break;

// User-defined function calls:
          case cFCall:
              {
                  unsigned index = ByteCode[++IP];
                  unsigned params = data->FuncPtrs[index].params;
                  double retVal =
                      data->FuncPtrs[index].ptr(&Stack[SP-params+1]);
                  SP -= params-1;
                  Stack[SP] = retVal;
                  break;
              }

          case cPCall:
              {
                  unsigned index = ByteCode[++IP];
                  unsigned params = data->FuncParsers[index]->data->varAmount;
                  double retVal =
                      data->FuncParsers[index]->Eval(&Stack[SP-params+1]);
                  SP -= params-1;
                  Stack[SP] = retVal;
                  break;
              }


#ifdef SUPPORT_OPTIMIZER
          case   cVar: break; // Paranoia. These should never exist
          case   cDup: Stack[SP+1] = Stack[SP]; ++SP; break;
          case   cInv:
              if(Stack[SP] == 0.0) { evalErrorType=1; return 0; }
              Stack[SP] = 1.0/Stack[SP];
              break;
#endif

// Variables:
          default:
              Stack[++SP] = Vars[ByteCode[IP]-VarBegin];
        }
    }

    evalErrorType=0;
    return Stack[SP];
}


#ifdef FUNCTIONPARSER_SUPPORT_DEBUG_OUTPUT
namespace
{
    inline void printHex(std::ostream& dest, unsigned n)
    {
        dest.width(8); dest.fill('0'); hex(dest); //uppercase(dest);
        dest << n;
    }
}

void FunctionParser::PrintByteCode(std::ostream& dest) const
{
    const unsigned* const ByteCode = data->ByteCode;
    const double* const Immed = data->Immed;

    for(unsigned IP=0, DP=0; IP<data->ByteCodeSize; ++IP)
    {
        printHex(dest, IP);
        dest << ": ";

        unsigned opcode = ByteCode[IP];

        switch(opcode)
        {
          case cIf:
              dest << "jz\t";
              printHex(dest, ByteCode[IP+1]+1);
              dest << endl;
              IP += 2;
              break;

          case cJump:
              dest << "jump\t";
              printHex(dest, ByteCode[IP+1]+1);
              dest << endl;
              IP += 2;
              break;
          case cImmed:
              dest.precision(10);
              dest << "push\t" << Immed[DP++] << endl;
              break;

          case cFCall:
              {
                  unsigned index = ByteCode[++IP];
                  Data::VarMap_t::const_iterator iter =
                      data->FuncPtrNames.begin();
                  while(iter->second != index) ++iter;
                  dest << "call\t" << iter->first << endl;
                  break;
              }

          case cPCall:
              {
                  unsigned index = ByteCode[++IP];
                  Data::VarMap_t::const_iterator iter =
                      data->FuncParserNames.begin();
                  while(iter->second != index) ++iter;
                  dest << "call\t" << iter->first << endl;
                  break;
              }

          default:
              if(opcode < VarBegin)
              {
                  string n;
                  switch(opcode)
                  {
                    case cNeg: n = "neg"; break;
                    case cAdd: n = "add"; break;
                    case cSub: n = "sub"; break;
                    case cMul: n = "mul"; break;
                    case cDiv: n = "div"; break;
                    case cMod: n = "mod"; break;
                    case cPow: n = "pow"; break;
                    case cEqual: n = "eq"; break;
                    case cLess: n = "lt"; break;
                    case cGreater: n = "gt"; break;
                    case cAnd: n = "and"; break;
                    case cOr: n = "or"; break;
                    case cDeg: n = "deg"; break;
                    case cRad: n = "rad"; break;

#ifndef DISABLE_EVAL
                    case cEval: n = "call\t0"; break;
#endif

#ifdef SUPPORT_OPTIMIZER
                    case cVar: n = "(var)"; break;
                    case cDup: n = "dup"; break;
                    case cInv: n = "inv"; break;
#endif

                    default: n = Functions[opcode-cAbs].name;
                  }
                  dest << n << endl;
              }
              else
              {
                  dest << "push\tVar" << opcode-VarBegin << endl;
              }
        }
    }
}
#endif


//========================================================================
// Optimization code was contributed by Bisqwit (http://iki.fi/bisqwit/)
//========================================================================
#ifdef SUPPORT_OPTIMIZER

#include <list>
#include <utility>

#define CONSTANT_E     2.71828182845904509080  // exp(1)
#define CONSTANT_PI    M_PI                    // atan2(0,-1)
#define CONSTANT_L10   2.30258509299404590109  // log(10)
#define CONSTANT_L10I  0.43429448190325176116  // 1/log(10)
#define CONSTANT_L10E  CONSTANT_L10I           // log10(e)
#define CONSTANT_L10EI CONSTANT_L10            // 1/log10(e)
#define CONSTANT_DR    (180.0 / M_PI)          // 180/pi
#define CONSTANT_RD    (M_PI / 180.0)          // pi/180

namespace {
class compres
{
    // states: 0=false, 1=true, 2=unknown
public:
    compres(bool b) : state(b) {}
    compres(char v) : state(v) {}
    // is it?
    operator bool() const { return state != 0; }
    // is it not?
    bool operator! () const { return state != 1; }
    bool operator==(bool b) const { return state != !b; }
    bool operator!=(bool b) const { return state != b; }
private:
    char state;
};

const compres maybe = (char)2;

struct CodeTree;

class SubTree
{
    CodeTree *tree;
    bool sign;  // Only possible when parent is cAdd or cMul

    inline void flipsign() { sign = !sign; }
public:
    SubTree();
    SubTree(double value);
    SubTree(const SubTree &b);
    SubTree(const CodeTree &b);

    ~SubTree();
    const SubTree &operator= (const SubTree &b);
    const SubTree &operator= (const CodeTree &b);

    bool getsign() const { return sign; }

    const CodeTree* operator-> () const { return tree; }
    const CodeTree& operator* () const { return *tree; }
    struct CodeTree* operator-> () { return tree; }
    struct CodeTree& operator* () { return *tree; }

    bool operator< (const SubTree& b) const;
    bool operator== (const SubTree& b) const;
    void Negate(); // Note: Parent must be cAdd
    void Invert(); // Note: Parent must be cMul

    void CheckConstNeg();
    void CheckConstInv();
};

bool IsNegate(const SubTree &p1, const SubTree &p2);
bool IsInverse(const SubTree &p1, const SubTree &p2);

typedef list<SubTree> paramlist;

struct CodeTreeData
{
    paramlist args;

private:
    unsigned op;       // Operation
    double value;      // In case of cImmed
    unsigned var;      // In case of cVar
    unsigned funcno;   // In case of cFCall, cPCall

public:
    CodeTreeData() : op(cAdd) {}
    ~CodeTreeData() {}

    void SetOp(unsigned newop)     { op=newop; }
    void SetFuncNo(unsigned newno) { funcno=newno; }
    unsigned GetFuncNo() const { return funcno; }

    bool IsFunc() const  { return op == cFCall || op == cPCall; }
    bool IsImmed() const { return op == cImmed; }
    bool IsVar() const   { return op == cVar; }
    inline unsigned GetOp() const { return op; }
    inline double GetImmed() const
    {
        return value;
    }
    inline unsigned GetVar() const
    {
        return var;
    }

    void AddParam(const SubTree &p)
    {
        args.push_back(p);
    }
    void SetVar(unsigned v)
    {
        args.clear();
        op  = cVar;
        var = v;
    }
    void SetImmed(double v)
    {
        args.clear();
        op       = cImmed;
        value    = orig = v;
        inverted = negated = false;
    }
    void NegateImmed()
    {
        negated = !negated;
        UpdateValue();
    }
    void InvertImmed()
    {
        inverted = !inverted;
        UpdateValue();
    }

    bool IsOriginal() const { return !(IsInverted() || IsNegated()); }
    bool IsInverted() const { return inverted; }
    bool IsNegated() const { return negated; }
    bool IsInvertedOriginal() const { return IsInverted() && !IsNegated(); }
    bool IsNegatedOriginal() const { return !IsInverted() && IsNegated(); }

private:
    void UpdateValue()
    {
        value = orig;
        if(IsInverted()) { value = 1.0 / value;
                           // FIXME: potential divide by zero.
                         }
        if(IsNegated()) value = -value;
    }

    double orig;
    bool inverted;
    bool negated;
protected:
    // Ensure we don't accidentally copy this
    void operator=(const CodeTreeData &b);
};


class CodeTreeDataPtr
{
    typedef pair<CodeTreeData, unsigned> p_t;
    typedef p_t* pp;
    mutable pp p;

    void Alloc()   const { ++p->second; }
    void Dealloc() const { if(!--p->second) delete p; p = 0; }

    void PrepareForWrite()
    {
        // We're ready if we're the only owner.
        if(p->second == 1) return;

        // Then make a clone.
        p_t *newtree = new p_t(p->first, 1);
        // Forget the old
        Dealloc();
        // Keep the new
        p = newtree;
    }

public:
    CodeTreeDataPtr() : p(new p_t) { p->second = 1; }
    CodeTreeDataPtr(const CodeTreeDataPtr &b): p(b.p) { Alloc(); }
    ~CodeTreeDataPtr() { Dealloc(); }
    const CodeTreeDataPtr &operator= (const CodeTreeDataPtr &b)
    {
        b.Alloc();
        Dealloc();
        p = b.p;
        return *this;
    }
    const CodeTreeData *operator-> () const { return &p->first; }
    const CodeTreeData &operator*  () const { return p->first; }
    CodeTreeData *operator-> () { PrepareForWrite(); return &p->first; }
    CodeTreeData &operator*  () { PrepareForWrite(); return p->first; }

    void Shock();
};


#define CHECKCONSTNEG(item, op) \
    ((op)==cMul) \
       ? (item).CheckConstInv() \
       : (item).CheckConstNeg()

struct CodeTree
{
    CodeTreeDataPtr data;

private:
    typedef paramlist::iterator pit;
    typedef paramlist::const_iterator pcit;

    /*
    template<unsigned v> inline void chk() const
    {
    }
    */

public:
    const pcit GetBegin() const { return data->args.begin(); }
    const pcit GetEnd()   const { return data->args.end(); }
    const pit GetBegin() { return data->args.begin(); }
    const pit GetEnd()   { return data->args.end(); }
    const SubTree& getp0() const { /*chk<1>();*/pcit tmp=GetBegin();               return *tmp; }
    const SubTree& getp1() const { /*chk<2>();*/pcit tmp=GetBegin(); ++tmp;        return *tmp; }
    const SubTree& getp2() const { /*chk<3>();*/pcit tmp=GetBegin(); ++tmp; ++tmp; return *tmp; }
    unsigned GetArgCount() const { return data->args.size(); }
    void Erase(const pit p)      { data->args.erase(p); }

    SubTree& getp0() { /*chk<1>();*/pit tmp=GetBegin();               return *tmp; }
    SubTree& getp1() { /*chk<2>();*/pit tmp=GetBegin(); ++tmp;        return *tmp; }
    SubTree& getp2() { /*chk<3>();*/pit tmp=GetBegin(); ++tmp; ++tmp; return *tmp; }

    // set
    void SetImmed(double v) { data->SetImmed(v); }
    void SetOp(unsigned op) { data->SetOp(op); }
    void SetVar(unsigned v) { data->SetVar(v); }
    // get
    double GetImmed() const { return data->GetImmed(); }
    unsigned GetVar() const { return data->GetVar(); }
    unsigned GetOp() const  { return data->GetOp(); }
    // test
    bool IsImmed() const { return data->IsImmed(); }
    bool IsVar()   const { return data->IsVar(); }
    // act
    void AddParam(const SubTree &p) { data->AddParam(p); }
    void NegateImmed() { data->NegateImmed(); } // don't use when op!=cImmed
    void InvertImmed() { data->InvertImmed(); } // don't use when op!=cImmed

    compres NonZero() const { if(!IsImmed()) return maybe;
                              return GetImmed() != 0.0; }
    compres IsPositive() const { if(!IsImmed()) return maybe;
                                 return GetImmed() > 0.0; }

private:
    struct ConstList
    {
        double voidvalue;
        list<pit> cp;
        double value;
        unsigned size() const { return cp.size(); }
    };
    struct ConstList BuildConstList();
    void KillConst(const ConstList &cl)
    {
        for(list<pit>::const_iterator i=cl.cp.begin(); i!=cl.cp.end(); ++i)
            Erase(*i);
    }
    void FinishConst(const ConstList &cl)
    {
        if(cl.value != cl.voidvalue && cl.size() > 1) AddParam(cl.value);
        if(cl.value == cl.voidvalue || cl.size() > 1) KillConst(cl);
    }

public:
    CodeTree() {}
    CodeTree(double v) { SetImmed(v); }

    CodeTree(unsigned op, const SubTree &p)
    {
        SetOp(op);
        AddParam(p);
    }
    CodeTree(unsigned op, const SubTree &p1, const SubTree &p2)
    {
        SetOp(op);
        AddParam(p1);
        AddParam(p2);
    }

    bool operator== (const CodeTree& b) const;
    bool operator< (const CodeTree& b) const;

private:
    bool IsSortable() const
    {
        switch(GetOp())
        {
            case cAdd:  case cMul:
            case cEqual:
            case cAnd: case cOr:
            case cMax: case cMin:
                return true;
            default:
                return false;
        }
    }
    void SortIfPossible()
    {
        if(IsSortable())
        {
            data->args.sort();
        }
    }

    void ReplaceWithConst(double value)
    {
        SetImmed(value);

        /* REMEMBER TO CALL CheckConstInv / CheckConstNeg
         * FOR PARENT SubTree, OR MAYHEM HAPPENS
         */
    }

    void ReplaceWith(const CodeTree &b)
    {
        // If b is child of *this, mayhem
        // happens. So we first make a clone
        // and then proceed with copy.
        CodeTreeDataPtr tmp = b.data;
        tmp.Shock();
        data = tmp;
    }

    void ReplaceWith(unsigned op, const SubTree &p)
    {
        ReplaceWith(CodeTree(op, p));
    }

    void ReplaceWith(unsigned op, const SubTree &p1, const SubTree &p2)
    {
        ReplaceWith(CodeTree(op, p1, p2));
    }

    void OptimizeConflict()
    {
        // This optimization does this: x-x = 0, x/x = 1, a+b-a = b.

        if(GetOp() == cAdd || GetOp() == cMul)
        {
        Redo:
            pit a, b;
            for(a=GetBegin(); a!=GetEnd(); ++a)
            {
                for(b=GetBegin(); ++b != GetEnd(); )
                {
                    const SubTree &p1 = *a;
                    const SubTree &p2 = *b;

                    if(GetOp() == cMul ? IsInverse(p1,p2)
                                       : IsNegate(p1,p2))
                    {
                        // These parameters complement each others out
                        Erase(b);
                        Erase(a);
                        goto Redo;
                    }
                }
            }
        }
        OptimizeRedundant();
    }

    void OptimizeRedundant()
    {
        // This optimization does this: min()=0, max()=0, add()=0, mul()=1

        if(!GetArgCount())
        {
            if(GetOp() == cAdd || GetOp() == cMin || GetOp() == cMax)
                ReplaceWithConst(0);
            else if(GetOp() == cMul)
                ReplaceWithConst(1);
            return;
        }

        // And this: mul(x) = x, min(x) = x, max(x) = x, add(x) = x

        if(GetArgCount() == 1)
        {
            if(GetOp() == cMul || GetOp() == cAdd || GetOp() == cMin || GetOp() == cMax)
                if(!getp0().getsign())
                    ReplaceWith(*getp0());
        }

        OptimizeDoubleNegations();
    }

    void OptimizeDoubleNegations()
    {
        if(GetOp() == cAdd)
        {
            // Eschew double negations

            // If any of the elements is cMul
            // and has a numeric constant, negate
            // the constant and negate sign.

            for(pit a=GetBegin(); a!=GetEnd(); ++a)
            {
                SubTree &pa = *a;
                if(pa.getsign()
                && pa->GetOp() == cMul)
                {
                    CodeTree &p = *pa;
                    for(pit b=p.GetBegin();
                            b!=p.GetEnd(); ++b)
                    {
                        SubTree &pb = *b;
                        if(pb->IsImmed())
                        {
                            pb.Negate();
                            pa.Negate();
                            break;
                        }
                    }
                }
            }
        }

        if(GetOp() == cMul)
        {
            // If any of the elements is cPow
            // and has a numeric exponent, negate
            // the exponent and negate sign.

            for(pit a=GetBegin(); a!=GetEnd(); ++a)
            {
                SubTree &pa = *a;
                if(pa.getsign() && pa->GetOp() == cPow)
                {
                    CodeTree &p = *pa;
                    if(p.getp1()->IsImmed())
                    {
                        // negate ok for pow when op=cImmed
                        p.getp1().Negate();
                        pa.Negate();
                    }
                }
            }
        }
    }

    void OptimizeConstantMath1()
    {
        // This optimization does three things:
        //      - For adding groups:
        //          Constants are added together.
        //      - For multiplying groups:
        //          Constants are multiplied together.
        //      - For function calls:
        //          If all parameters are constants,
        //          the call is replaced with constant value.

        // First, do this:
        OptimizeAddMulFlat();

        switch(GetOp())
        {
            case cAdd:
            {
                ConstList cl = BuildConstList();
                FinishConst(cl);
                break;
            }
            case cMul:
            {
                ConstList cl = BuildConstList();

                if(cl.value == 0.0) ReplaceWithConst(0.0);
                else FinishConst(cl);

                break;
            }
            #define ConstantUnaryFun(token, fun) \
                case token: { const SubTree &p0 = getp0(); \
                    if(p0->IsImmed()) ReplaceWithConst(fun(p0->GetImmed())); \
                    break; }
            #define ConstantBinaryFun(token, fun) \
                case token: { const SubTree &p0 = getp0(); \
                              const SubTree &p1 = getp1(); \
                    if(p0->IsImmed() && \
                       p1->IsImmed()) ReplaceWithConst(fun(p0->GetImmed(), p1->GetImmed())); \
                    break; }

            // FIXME: potential invalid parameters for functions
            //        can cause exceptions here

            ConstantUnaryFun(cAbs,   fabs);
            ConstantUnaryFun(cAcos,  acos);
            ConstantUnaryFun(cAsin,  asin);
            ConstantUnaryFun(cAtan,  atan);
            ConstantUnaryFun(cCeil,  ceil);
            ConstantUnaryFun(cCos,   cos);
            ConstantUnaryFun(cCosh,  cosh);
            ConstantUnaryFun(cFloor, floor);
            ConstantUnaryFun(cLog,   log);
            ConstantUnaryFun(cSin,   sin);
            ConstantUnaryFun(cSinh,  sinh);
            ConstantUnaryFun(cTan,   tan);
            ConstantUnaryFun(cTanh,  tanh);
            ConstantBinaryFun(cAtan2, atan2);
            ConstantBinaryFun(cMax,   Max);
            ConstantBinaryFun(cMin,   Min);
            ConstantBinaryFun(cMod,   fmod); // not a func, but belongs here too
            ConstantBinaryFun(cPow,   pow);

            case cNeg:
            case cSub:
            case cDiv:
                /* Unreached (nonexistent operator)
                 * TODO: internal error here?
                 */
                break;

            case cCot:
            case cCsc:
            case cSec:
            case cDeg:
            case cRad:
            case cLog10:
            case cSqrt:
            case cExp:
                /* Unreached (nonexistent function)
                 * TODO: internal error here?
                 */
                 break;
        }

        OptimizeConflict();
    }

    void OptimizeAddMulFlat()
    {
        // This optimization flattens the topography of the tree.
        //   Examples:
        //       x + (y+z) = x+y+z
        //       x * (y/z) = x*y/z
        //       x / (y/z) = x/y*z

        if(GetOp() == cAdd || GetOp() == cMul)
        {
            // If children are same type as parent add them here
            for(pit b, a=GetBegin(); a!=GetEnd(); a=b)
            {
                const SubTree &pa = *a;  b=a; ++b;
                if(pa->GetOp() != GetOp()) continue;

                // Child is same type
                for(pcit c=pa->GetBegin();
                         c!=pa->GetEnd();
                         ++c)
                {
                    const SubTree &pb = *c;
                    if(pa.getsign())
                    {
                        // +a -(+b +c)
                        // means b and c will be negated

                        SubTree tmp = pb;
                        if(GetOp() == cMul)
                            tmp.Invert();
                        else
                            tmp.Negate();
                        AddParam(tmp);
                    }
                    else
                        AddParam(pb);
                }
                Erase(a);

                // Note: OptimizeConstantMath1() would be a good thing to call next.
            }
        }
    }

    void OptimizeLinearCombine()
    {
        // This optimization does the following:
        //
        //   x*x*x*x -> x^4
        //   x+x+x+x -> x*4
        //   x*x     -> x^2
        //   x/z/z   ->
        //

        // Remove conflicts first, so we don't have to worry about signs.
        OptimizeConflict();

        bool didchanges = false;
        if(GetOp() == cAdd || GetOp() == cMul)
        {
        Redo:
            for(pit a=GetBegin(); a!=GetEnd(); ++a)
            {
                const SubTree &pa = *a;

                list<pit> poslist;

                for(pit b=a; ++b!=GetEnd(); )
                {
                    const SubTree &pb = *b;
                    if(*pa == *pb)
                        poslist.push_back(b);
                }

                unsigned min = 2;
                if(poslist.size() >= min)
                {
                    SubTree arvo = pa;
                    bool negate = arvo.getsign();

                    double factor = poslist.size() + 1;

                    if(negate)
                    {
                        arvo.Negate();
                        factor = -factor;
                    }

                    CodeTree tmp(GetOp()==cAdd ? cMul : cPow,
                                 arvo,
                                 factor);

                    list<pit>::const_iterator j;
                    for(j=poslist.begin(); j!=poslist.end(); ++j)
                        Erase(*j);
                    poslist.clear();

                    *a = tmp;
                    didchanges = true;
                    goto Redo;
                }
            }
        }
        if(didchanges)
        {
            // As a result, there might be need for this:
            OptimizeAddMulFlat();
            // And this:
            OptimizeRedundant();
        }
    }

    void OptimizeLogarithm()
    {
        /*
            This is basic logarithm math:
              pow(X,Y)/log(Y) = X
              log(X)/log(Y) = logY(X)
              log(X^Y)      = log(X)*Y
              log(X*Y)      = log(X)+log(Y)
              exp(log(X)*Y) = X^Y

            This function does these optimizations:
               pow(const_E, log(x))   = x
               pow(const_E, log(x)*y) = x^y
               pow(10,      log(x)*const_L10I*y) = x^y
               pow(z,       log(x)/log(z)*y)     = x^y

            And this:
               log(x^z)             = z * log(x)
            Which automatically causes these too:
               log(pow(const_E, x))         = x
               log(pow(y,       x))         = x * log(y)
               log(pow(pow(const_E, y), x)) = x*y

            And it does this too:
               log(x) + log(y) + log(z) = log(x * y * z)
               log(x * exp(y)) = log(x) + y

        */

        // Must be already in exponential form.

        // Optimize exponents before doing something.
        OptimizeExponents();

        if(GetOp() == cLog)
        {
            // We should have one parameter for log() function.
            // If we don't, we're screwed.

            const SubTree &p = getp0();

            if(p->GetOp() == cPow)
            {
                // Found log(x^y)
                SubTree p0 = p->getp0(); // x
                SubTree p1 = p->getp1(); // y

                // Build the new logarithm.
                CodeTree tmp(GetOp(), p0);  // log(x)

                // Become log(x) * y
                ReplaceWith(cMul, tmp, p1);
            }
            else if(p->GetOp() == cMul)
            {
                // Redefine &p nonconst
                SubTree &p = getp0();

                p->OptimizeAddMulFlat();
                p->OptimizeExponents();
                CHECKCONSTNEG(p, p->GetOp());

                list<SubTree> adds;

                for(pit b, a = p->GetBegin();
                           a != p->GetEnd(); a=b)
                {
                    SubTree &pa = *a;  b=a; ++b;
                    if(pa->GetOp() == cPow
                    && pa->getp0()->IsImmed()
                    && pa->getp0()->GetImmed() == CONSTANT_E)
                    {
                        adds.push_back(pa->getp1());
                        p->Erase(a);
                        continue;
                    }
                }
                if(adds.size())
                {
                    CodeTree tmp(cAdd, *this);

                    list<SubTree>::const_iterator i;
                    for(i=adds.begin(); i!=adds.end(); ++i)
                        tmp.AddParam(*i);

                    ReplaceWith(tmp);
                }
            }
        }
        if(GetOp() == cAdd)
        {
            // Check which ones are logs.
            list<pit> poslist;

            for(pit a=GetBegin(); a!=GetEnd(); ++a)
            {
                const SubTree &pa = *a;
                if(pa->GetOp() == cLog)
                    poslist.push_back(a);
            }

            if(poslist.size() >= 2)
            {
                CodeTree tmp(cMul, 1.0); // eek

                list<pit>::const_iterator j;
                for(j=poslist.begin(); j!=poslist.end(); ++j)
                {
                    const SubTree &pb = **j;
                    // Take all of its children
                    for(pcit b=pb->GetBegin();
                             b!=pb->GetEnd();
                             ++b)
                    {
                        SubTree tmp2 = *b;
                        if(pb.getsign()) tmp2.Negate();
                        tmp.AddParam(tmp2);
                    }
                    Erase(*j);
                }
                poslist.clear();

                AddParam(CodeTree(cLog, tmp));
            }
            // Done, hopefully
        }
        if(GetOp() == cPow)
        {
            const SubTree &p0 = getp0();
            SubTree &p1 = getp1();

            if(p0->IsImmed() && p0->GetImmed() == CONSTANT_E
            && p1->GetOp() == cLog)
            {
                // pow(const_E, log(x)) = x
                ReplaceWith(*(p1->getp0()));
            }
            else if(p1->GetOp() == cMul)
            {
                //bool didsomething = true;

                pit poslogpos; bool foundposlog = false;
                pit neglogpos; bool foundneglog = false;

                ConstList cl = p1->BuildConstList();

                for(pit a=p1->GetBegin(); a!=p1->GetEnd(); ++a)
                {
                    const SubTree &pa = *a;
                    if(pa->GetOp() == cLog)
                    {
                        if(!pa.getsign())
                        {
                            foundposlog = true;
                            poslogpos   = a;
                        }
                        else if(*p0 == *(pa->getp0()))
                        {
                            foundneglog = true;
                            neglogpos   = a;
                        }
                    }
                }

                if(p0->IsImmed()
                && p0->GetImmed() == 10.0
                && cl.value == CONSTANT_L10I
                && foundposlog)
                {
                    SubTree base = (*poslogpos)->getp0();
                    p1->KillConst(cl);
                    p1->Erase(poslogpos);
                    p1->OptimizeRedundant();
                    SubTree mul = p1;

                    ReplaceWith(cPow, base, mul);

                    // FIXME: what optimizations should be done now?
                    return;
                }

                // Put back the constant
                FinishConst(cl);

                if(p0->IsImmed()
                && p0->GetImmed() == CONSTANT_E
                && foundposlog)
                {
                    SubTree base = (*poslogpos)->getp0();
                    p1->Erase(poslogpos);

                    p1->OptimizeRedundant();
                    SubTree mul = p1;

                    ReplaceWith(cPow, base, mul);

                    // FIXME: what optimizations should be done now?
                    return;
                }

                if(foundposlog
                && foundneglog
                && *((*neglogpos)->getp0()) == *p0)
                {
                    SubTree base = (*poslogpos)->getp0();
                    p1->Erase(poslogpos);
                    p1->Erase(neglogpos);

                    p1->OptimizeRedundant();
                    SubTree mul = p1;

                    ReplaceWith(cPow, base, mul);

                    // FIXME: what optimizations should be done now?
                    return;
                }
            }
        }
    }

    void OptimizeFunctionCalls()
    {
        /* Goals: sin(asin(x)) = x
         *        cos(acos(x)) = x
         *        tan(atan(x)) = x
         * NOTE:
         *   Do NOT do these:
         *     asin(sin(x))
         *     acos(cos(x))
         *     atan(tan(x))
         *   Because someone might want to wrap the angle.
         */
        // FIXME: TODO
    }

    void OptimizePowMulAdd()
    {
        // x^3 * x -> x^4
        // x*3 + x -> x*4
        // FIXME: Do those

        // x^1 -> x
        if(GetOp() == cPow)
        {
            const SubTree &base     = getp0();
            const SubTree &exponent = getp1();

            if(exponent->IsImmed())
            {
                if(exponent->GetImmed() == 1.0)
                    ReplaceWith(*base);
                else if(exponent->GetImmed() == 0.0
                     && base->NonZero())
                    ReplaceWithConst(1.0);
            }
        }
    }

    void OptimizeExponents()
    {
        /* Goals:
         *     (x^y)^z   -> x^(y*z)
         *     x^y * x^z -> x^(y+z)
         */
        // First move to exponential form.
        OptimizeLinearCombine();

        bool didchanges = false;

    Redo:
        if(GetOp() == cPow)
        {
            // (x^y)^z   -> x^(y*z)

            const SubTree &p0 = getp0();
            const SubTree &p1 = getp1();
            if(p0->GetOp() == cPow)
            {
                CodeTree tmp(cMul, p0->getp1(), p1);
                tmp.Optimize();

                ReplaceWith(cPow, p0->getp0(), tmp);

                didchanges = true;
                goto Redo;
            }
        }
        if(GetOp() == cMul)
        {
            // x^y * x^z -> x^(y+z)

            for(pit a=GetBegin(); a!=GetEnd(); ++a)
            {
                const SubTree &pa = *a;

                if(pa->GetOp() != cPow) continue;

                list<pit> poslist;

                for(pit b=a; ++b != GetEnd(); )
                {
                    const SubTree &pb = *b;
                    if(pb->GetOp() == cPow
                    && *(pa->getp0())
                    == *(pb->getp0()))
                    {
                        poslist.push_back(b);
                    }
                }

                if(poslist.size() >= 1)
                {
                    poslist.push_back(a);

                    CodeTree base = *(pa->getp0());

                    CodeTree exponent(cAdd, 0.0); //eek

                    // Collect all exponents to cAdd
                    list<pit>::const_iterator i;
                    for(i=poslist.begin(); i!=poslist.end(); ++i)
                    {
                        const SubTree &pb = **i;

                        SubTree tmp2 = pb->getp1();
                        if(pb.getsign()) tmp2.Invert();

                        exponent.AddParam(tmp2);
                    }

                    exponent.Optimize();

                    CodeTree result(cPow, base, exponent);

                    for(i=poslist.begin(); i!=poslist.end(); ++i)
                        Erase(*i);
                    poslist.clear();

                    AddParam(result); // We're cMul, remember

                    didchanges = true;
                    goto Redo;
                }
            }
        }

        OptimizePowMulAdd();

        if(didchanges)
        {
            // As a result, there might be need for this:
            OptimizeConflict();
        }
    }

    void OptimizeLinearExplode()
    {
        // x^2 -> x*x
        // But only if x is just a simple thing

        // Won't work on anything else.
        if(GetOp() != cPow) return;

        // TODO TODO TODO
    }

    void OptimizePascal()
    {
#if 0    // Too big, too specific, etc

        // Won't work on anything else.
        if(GetOp() != cAdd) return;

        // Must be done after OptimizeLinearCombine();

        // Don't need pascal triangle
        // Coefficient for x^a * y^b * z^c = 3! / (a! * b! * c!)

        // We are greedy and want other than just binomials
        // FIXME

        // note: partial ones are also nice
        //     x*x + x*y + y*y
        //   = (x+y)^2 - x*y
        //
        //     x x * x y * + y y * +
        // ->  x y + dup * x y * -
#endif
    }

public:

    void Optimize();

    void Assemble(vector<unsigned> &byteCode,
                  vector<double>   &immed) const;

    void FinalOptimize()
    {
        // First optimize each parameter.
        for(pit a=GetBegin(); a!=GetEnd(); ++a)
            (*a)->FinalOptimize();

        /* These things are to be done:
         *
         * x * CONSTANT_DR        -> cDeg(x)
         * x * CONSTANT_RD        -> cRad(x)
         * pow(x, 0.5)            -> sqrt(x)
         * log(x) * CONSTANT_L10I -> log10(x)
         * pow(CONSTANT_E, x)     -> exp(x)
         * inv(sin(x))            -> csc(x)
         * inv(cos(x))            -> sec(x)
         * inv(tan(x))            -> cot(x)
         */


        if(GetOp() == cPow)
        {
            const SubTree &p0 = getp0();
            const SubTree &p1 = getp1();
            if(p0->GetOp()    == cImmed
            && p0->GetImmed() == CONSTANT_E)
            {
                ReplaceWith(cExp, p1);
            }
            else if(p1->GetOp()    == cImmed
                 && p1->GetImmed() == 0.5)
            {
                ReplaceWith(cSqrt, p0);
            }
        }
        if(GetOp() == cMul)
        {
            if(GetArgCount() == 1 && getp0().getsign())
            {
                /***/if(getp0()->GetOp() == cSin)ReplaceWith(cCsc, getp0()->getp0());
                else if(getp0()->GetOp() == cCos)ReplaceWith(cSec, getp0()->getp0());
                else if(getp0()->GetOp() == cTan)ReplaceWith(cCot, getp0()->getp0());
            }
        }
        // Separate "if", because op may have just changed
        if(GetOp() == cMul)
        {
            CodeTree *found_log = 0;

            ConstList cl = BuildConstList();

            for(pit a=GetBegin(); a!=GetEnd(); ++a)
            {
                SubTree &pa = *a;
                if(pa->GetOp() == cLog && !pa.getsign())
                    found_log = &*pa;
            }
            if(cl.value == CONSTANT_L10I && found_log)
            {
                // Change the log() to log10()
                found_log->SetOp(cLog10);
                // And forget the constant
                KillConst(cl);
            }
            else if(cl.value == CONSTANT_DR)
            {
                OptimizeRedundant();
                ReplaceWith(cDeg, *this);
            }
            else if(cl.value == CONSTANT_RD)
            {
                OptimizeRedundant();
                ReplaceWith(cRad, *this);
            }
            else FinishConst(cl);
        }

        SortIfPossible();
    }
};

void CodeTreeDataPtr::Shock()
{
 /*
    PrepareForWrite();
    paramlist &p2 = (*this)->args;
    for(paramlist::iterator i=p2.begin(); i!=p2.end(); ++i)
    {
        (*i)->data.Shock();
    }
 */
}

CodeTree::ConstList CodeTree::BuildConstList()
{
    ConstList result;
    result.value     =
    result.voidvalue = GetOp()==cMul ? 1.0 : 0.0;

    list<pit> &cp = result.cp;
    for(pit b, a=GetBegin(); a!=GetEnd(); a=b)
    {
        SubTree &pa = *a;  b=a; ++b;
        if(!pa->IsImmed()) continue;

        double thisvalue = pa->GetImmed();
        if(thisvalue == result.voidvalue)
        {
            // This value is no good, forget it
            Erase(a);
            continue;
        }
        if(GetOp() == cMul)
            result.value *= thisvalue;
        else
            result.value += thisvalue;
        cp.push_back(a);
    }
    if(GetOp() == cMul)
    {
        /*
          Jos joku niistä arvoista on -1 eikä se ole ainoa arvo,
          niin joku muu niistä arvoista negatoidaan.
        */
        for(bool done=false; cp.size() > 1 && !done; )
        {
            done = true;
            for(list<pit>::iterator b,a=cp.begin(); a!=cp.end(); a=b)
            {
                b=a; ++b;
                if((**a)->GetImmed() == -1.0)
                {
                    Erase(*a);
                    cp.erase(a);

                    // take randomly something
                    (**cp.begin())->data->NegateImmed();
                    if(cp.size() < 2)break;
                    done = false;
                }
            }
        }
    }
    return result;
}

void CodeTree::Assemble
   (vector<unsigned> &byteCode,
    vector<double>   &immed) const
{
    #define AddCmd(op) byteCode.push_back((op))
    #define AddConst(v) do { \
        byteCode.push_back(cImmed); \
        immed.push_back((v)); \
    } while(0)

    if(IsVar())
    {
        AddCmd(GetVar());
        return;
    }
    if(IsImmed())
    {
        AddConst(GetImmed());
        return;
    }

    switch(GetOp())
    {
        case cAdd:
        case cMul:
        {
            unsigned opcount = 0;
            for(pcit a=GetBegin(); a!=GetEnd(); ++a)
            {
                const SubTree &pa = *a;

                if(opcount < 2) ++opcount;

                bool pnega = pa.getsign();

                bool done = false;
                if(pa->IsImmed())
                {
                    if(GetOp() == cMul
                    && pa->data->IsInverted()
                    && (pnega || opcount==2)
                      )
                    {
                        CodeTree tmp = *pa;
                        tmp.data->InvertImmed();
                        tmp.Assemble(byteCode, immed);
                        pnega = !pnega;
                        done = true;
                    }
                    else if(GetOp() == cAdd
                    && (pa->data->IsNegatedOriginal()
                //     || pa->GetImmed() < 0
                       )
                    && (pnega || opcount==2)
                           )
                    {
                        CodeTree tmp = *pa;
                        tmp.data->NegateImmed();
                        tmp.Assemble(byteCode, immed);
                        pnega = !pnega;
                        done = true;
                    }
                }
                if(!done)
                    pa->Assemble(byteCode, immed);

                if(opcount == 2)
                {
                    unsigned tmpop = GetOp();
                    if(pnega) // negate
                    {
                        tmpop = (tmpop == cMul) ? cDiv : cSub;
                    }
                    AddCmd(tmpop);
                }
                else if(pnega)
                {
                    if(GetOp() == cMul) AddCmd(cInv);
                    else AddCmd(cNeg);
                }
            }
            break;
        }
        case cIf:
        {
            // If the parameter amount is != 3, we're screwed.
            getp0()->Assemble(byteCode, immed);

            unsigned ofs = byteCode.size();
            AddCmd(cIf);
            AddCmd(0); // code index
            AddCmd(0); // immed index

            getp1()->Assemble(byteCode, immed);

            byteCode[ofs+1] = byteCode.size()+2;
            byteCode[ofs+2] = immed.size();

            ofs = byteCode.size();
            AddCmd(cJump);
            AddCmd(0); // code index
            AddCmd(0); // immed index

            getp2()->Assemble(byteCode, immed);

            byteCode[ofs+1] = byteCode.size()-1;
            byteCode[ofs+2] = immed.size();

            break;
        }
        case cFCall:
        {
            // If the parameter count is invalid, we're screwed.
            for(pcit a=GetBegin(); a!=GetEnd(); ++a)
            {
                const SubTree &pa = *a;
                pa->Assemble(byteCode, immed);
            }
            AddCmd(GetOp());
            AddCmd(data->GetFuncNo());
            break;
        }
        case cPCall:
        {
            // If the parameter count is invalid, we're screwed.
            for(pcit a=GetBegin(); a!=GetEnd(); ++a)
            {
                const SubTree &pa = *a;
                pa->Assemble(byteCode, immed);
            }
            AddCmd(GetOp());
            AddCmd(data->GetFuncNo());
            break;
        }
        default:
        {
            // If the parameter count is invalid, we're screwed.
            for(pcit a=GetBegin(); a!=GetEnd(); ++a)
            {
                const SubTree &pa = *a;
                pa->Assemble(byteCode, immed);
            }
            AddCmd(GetOp());
            break;
        }
    }
}

void CodeTree::Optimize()
{
    // Phase:
    //   Phase 0: Do local optimizations.
    //   Phase 1: Optimize each.
    //   Phase 2: Do local optimizations again.

    for(unsigned phase=0; phase<=2; ++phase)
    {
        if(phase == 1)
        {
            // Optimize each parameter.
            for(pit a=GetBegin(); a!=GetEnd(); ++a)
            {
                (*a)->Optimize();
                CHECKCONSTNEG(*a, GetOp());
            }
            continue;
        }
        if(phase == 0 || phase == 2)
        {
            // Do local optimizations.

            OptimizeConstantMath1();
            OptimizeLogarithm();
            OptimizeFunctionCalls();
            OptimizeExponents();
            OptimizeLinearExplode();
            OptimizePascal();

            /* Optimization paths:

               doublenegations=
               redundant= * doublenegations
               conflict= * redundant
               addmulflat=
               constantmath1= addmulflat * conflict
               linearcombine= conflict * addmulflat¹ redundant¹
               powmuladd=
               exponents= linearcombine * powmuladd conflict¹
               logarithm= exponents *
               functioncalls= IDLE
               linearexplode= IDLE
               pascal= IDLE

               * = actions here
               ¹ = only if made changes
            */
        }
    }
}


bool CodeTree::operator== (const CodeTree& b) const
{
    if(GetOp() != b.GetOp()) return false;
    if(IsImmed()) if(GetImmed()  != b.GetImmed())  return false;
    if(IsVar())   if(GetVar()    != b.GetVar())    return false;
    if(data->IsFunc())
        if(data->GetFuncNo() != b.data->GetFuncNo()) return false;
    return data->args == b.data->args;
}

bool CodeTree::operator< (const CodeTree& b) const
{
    if(GetArgCount() != b.GetArgCount())
        return GetArgCount() > b.GetArgCount();

    if(GetOp() != b.GetOp())
    {
        // sort immeds last
        if(IsImmed() != b.IsImmed()) return IsImmed() < b.IsImmed();

        return GetOp() < b.GetOp();
    }

    if(IsImmed())
    {
        if(GetImmed() != b.GetImmed()) return GetImmed() < b.GetImmed();
    }
    if(IsVar() && GetVar() != b.GetVar())
    {
        return GetVar() < b.GetVar();
    }
    if(data->IsFunc() && data->GetFuncNo() != b.data->GetFuncNo())
    {
        return data->GetFuncNo() < b.data->GetFuncNo();
    }

    pcit i = GetBegin(), j = b.GetBegin();
    for(; i != GetEnd(); ++i, ++j)
    {
        const SubTree &pa = *i, &pb = *j;

        if(!(pa == pb))
            return pa < pb;
    }
    return false;
}


bool IsNegate(const SubTree &p1, const SubTree &p2) /*const */
{
    if(p1->IsImmed() && p2->IsImmed())
    {
        return p1->GetImmed() == -p2->GetImmed();
    }
    if(p1.getsign() == p2.getsign()) return false;
    return *p1 == *p2;
}
bool IsInverse(const SubTree &p1, const SubTree &p2) /*const*/
{
    if(p1->IsImmed() && p2->IsImmed())
    {
        // FIXME: potential divide by zero.
        return p1->GetImmed() == 1.0 / p2->GetImmed();
    }
    if(p1.getsign() == p2.getsign()) return false;
    return *p1 == *p2;
}

SubTree::SubTree() : tree(new CodeTree), sign(false)
{
}

SubTree::SubTree(const SubTree &b) : tree(new CodeTree(*b.tree)), sign(b.sign)
{
}

#define SubTreeDecl(p1, p2) \
    SubTree::SubTree p1 : tree(new CodeTree p2), sign(false) { }

SubTreeDecl( (const CodeTree &b), (b) )
SubTreeDecl( (double value),             (value) )

#undef SubTreeDecl

SubTree::~SubTree()
{
    delete tree; tree=0;
}

const SubTree &SubTree::operator= (const SubTree &b)
{
    sign = b.sign;
    CodeTree *oldtree = tree;
    tree = new CodeTree(*b.tree);
    delete oldtree;
    return *this;
}
const SubTree &SubTree::operator= (const CodeTree &b)
{
    sign = false;
    CodeTree *oldtree = tree;
    tree = new CodeTree(b);
    delete oldtree;
    return *this;
}

bool SubTree::operator< (const SubTree& b) const
{
    if(getsign() != b.getsign()) return getsign() < b.getsign();
    return *tree < *b.tree;
}
bool SubTree::operator== (const SubTree& b) const
{
    return sign == b.sign && *tree == *b.tree;
}
void SubTree::Negate() // Note: Parent must be cAdd
{
    flipsign();
    CheckConstNeg();
}
void SubTree::CheckConstNeg()
{
    if(tree->IsImmed() && getsign())
    {
        tree->NegateImmed();
        sign = false;
    }
}
void SubTree::Invert() // Note: Parent must be cMul
{
    flipsign();
    CheckConstInv();
}
void SubTree::CheckConstInv()
{
    if(tree->IsImmed() && getsign())
    {
        tree->InvertImmed();
        sign = false;
    }
}

}//namespace

void FunctionParser::MakeTree(void *r) const
{
    // Dirty hack. Should be fixed.
    CodeTree* result = static_cast<CodeTree*>(r);

    vector<CodeTree> stack(1);

    #define GROW(n) do { \
        stacktop += n; \
        if(stack.size() <= stacktop) stack.resize(stacktop+1); \
    } while(0)

    #define EAT(n, opcode) do { \
        unsigned newstacktop = stacktop-n; \
        stack[stacktop].SetOp((opcode)); \
        for(unsigned a=0, b=(n); a<b; ++a) \
            stack[stacktop].AddParam(stack[newstacktop+a]); \
        stack.erase(stack.begin() + newstacktop, \
                    stack.begin() + stacktop); \
        stacktop = newstacktop; GROW(1); \
    } while(0)

    #define ADDCONST(n) do { \
        stack[stacktop].SetImmed((n)); \
        GROW(1); \
    } while(0)

    unsigned stacktop=0;

    list<unsigned> labels;

    const unsigned* const ByteCode = data->ByteCode;
    const unsigned ByteCodeSize = data->ByteCodeSize;
    const double* const Immed = data->Immed;

    for(unsigned IP=0, DP=0; ; ++IP)
    {
        while(labels.size() > 0
        && *labels.begin() == IP)
        {
            // The "else" of an "if" ends here
            EAT(3, cIf);
            labels.erase(labels.begin());
        }

        if(IP >= ByteCodeSize)
        {
            break;
        }

        unsigned opcode = ByteCode[IP];

        if(opcode == cIf)
        {
            IP += 2;
        }
        else if(opcode == cJump)
        {
            labels.push_front(ByteCode[IP+1]+1);
            IP += 2;
        }
        else if(opcode == cImmed)
        {
            ADDCONST(Immed[DP++]);
        }
        else if(opcode < VarBegin)
        {
            switch(opcode)
            {
                // Unary operators
                case cNeg:
                {
                    EAT(1, cAdd); // Unary minus is negative adding.
                    stack[stacktop-1].getp0().Negate();
                    break;
                }
                // Binary operators
                case cSub:
                {
                    EAT(2, cAdd); // Minus is negative adding
                    stack[stacktop-1].getp1().Negate();
                    break;
                }
                case cDiv:
                {
                    EAT(2, cMul); // Divide is inverse multiply
                    stack[stacktop-1].getp1().Invert();
                    break;
                }

                // ADD ALL TWO PARAMETER NON-FUNCTIONS HERE
                case cAdd: case cMul:
                case cMod: case cPow:
                case cEqual: case cLess: case cGreater:
                case cAnd: case cOr:
                    EAT(2, opcode);
                    break;

                case cFCall:
                {
                    unsigned index = ByteCode[++IP];
                    unsigned params = data->FuncPtrs[index].params;
                    EAT(params, opcode);
                    stack[stacktop-1].data->SetFuncNo(index);
                    break;
                }
                case cPCall:
                {
                    unsigned index = ByteCode[++IP];
                    unsigned params =
                        data->FuncParsers[index]->data->varAmount;
                    EAT(params, opcode);
                    stack[stacktop-1].data->SetFuncNo(index);
                    break;
                }

                // Converted to cMul on fly
                case cDeg:
                    ADDCONST(CONSTANT_DR);
                    EAT(2, cMul);
                    break;

                // Converted to cMul on fly
                case cRad:
                    ADDCONST(CONSTANT_RD);
                    EAT(2, cMul);
                    break;

                // Functions
                default:
                {
                    const FuncDefinition& func = Functions[opcode-cAbs];

                    unsigned paramcount = func.params;
#ifndef DISABLE_EVAL
                    if(opcode == cEval) paramcount = data->varAmount;
#endif
                    if(opcode == cSqrt)
                    {
                        // Converted on fly: sqrt(x) = x^0.5
                        opcode = cPow;
                        paramcount = 2;
                        ADDCONST(0.5);
                    }
                    if(opcode == cExp)
                    {
                        // Converted on fly: exp(x) = CONSTANT_E^x

                        opcode = cPow;
                        paramcount = 2;
                        // reverse the parameters... kludgey
                        stack[stacktop] = stack[stacktop-1];
                        stack[stacktop-1].SetImmed(CONSTANT_E);
                        GROW(1);
                    }
                    bool do_inv = false;
                    if(opcode == cCot) { do_inv = true; opcode = cTan; }
                    if(opcode == cCsc) { do_inv = true; opcode = cSin; }
                    if(opcode == cSec) { do_inv = true; opcode = cCos; }

                    bool do_log10 = false;
                    if(opcode == cLog10)
                    {
                        // Converted on fly: log10(x) = log(x) * CONSTANT_L10I
                        opcode = cLog;
                        do_log10 = true;
                    }
                    EAT(paramcount, opcode);
                    if(do_log10)
                    {
                        ADDCONST(CONSTANT_L10I);
                        EAT(2, cMul);
                    }
                    if(do_inv)
                    {
                        // Unary cMul, inverted. No need for "1.0"
                        EAT(1, cMul);
                        stack[stacktop-1].getp0().Invert();
                    }
                    break;
                }
            }
        }
        else
        {
            stack[stacktop].SetVar(opcode);
            GROW(1);
        }
    }

    if(!stacktop)
    {
        // ERROR: Stack does not have any values!
        return;
    }

    --stacktop; // Ignore the last element, it is always nop (cAdd).

    if(stacktop > 0)
    {
        // ERROR: Stack has too many values!
        return;
    }

    // Okay, the tree is now stack[0]
    *result = stack[0];
}

void FunctionParser::Optimize()
{
    copyOnWrite();

    CodeTree tree;
    MakeTree(&tree);

    // Do all sorts of optimizations
    tree.Optimize();
    // Last changes before assembly
    tree.FinalOptimize();

    // Now rebuild from the tree.

    vector<unsigned> byteCode;
    vector<double> immed;

#if 0
    byteCode.resize(Comp.ByteCodeSize);
    for(unsigned a=0; a<Comp.ByteCodeSize; ++a)byteCode[a] = Comp.ByteCode[a];

    immed.resize(Comp.ImmedSize);
    for(unsigned a=0; a<Comp.ImmedSize; ++a)immed[a] = Comp.Immed[a];
#else
    byteCode.clear(); immed.clear();
    tree.Assemble(byteCode, immed);
#endif

    delete[] data->ByteCode; data->ByteCode = 0;
    if((data->ByteCodeSize = byteCode.size()) > 0)
    {
        data->ByteCode = new unsigned[data->ByteCodeSize];
        for(unsigned a=0; a<byteCode.size(); ++a)
            data->ByteCode[a] = byteCode[a];
    }

    delete[] data->Immed; data->Immed = 0;
    if((data->ImmedSize = immed.size()) > 0)
    {
        data->Immed = new double[data->ImmedSize];
        for(unsigned a=0; a<immed.size(); ++a)
            data->Immed[a] = immed[a];
    }
}


#else /* !SUPPORT_OPTIMIZER */

/* keep the linker happy */
void FunctionParser::MakeTree(CodeTree *) const {}
void FunctionParser::Optimize()
{
    // Do nothing if no optimizations are supported.
}
#endif
