/* NOTE:
  Do not include this file in your project. The fparser.cc file #includes
this file internally and thus you don't need to do anything (other than keep
this file in the same directory as fparser.cc).

  Part of this file is generated code (by using the make_function_name_parser
utility, found in the development version of this library). It's not intended
to be modified by hand.
*/

        unsigned nameLength = 0;
        const unsigned maximumNameLength = 0x80000000U-8;
        /*
        Due to the manner the identifier lengths are returned from
        the readOpcode() function, the maximum supported length for
        identifiers is 0x7FFFFFFF bytes. We minus 8 here to add some
        buffer, because of the multibyteness of UTF-8.
        Function names are limited to 0xFFFF bytes instead, but because
        function names that long just are not defined, the point is moot.
        */
        const unsigned char* const uptr = (const unsigned char*) input;
        typedef signed char schar;
        while(likely(nameLength < maximumNameLength))
        {
            unsigned char byte = uptr[nameLength+0];
            /* Handle the common case of A-Za-z first */
            if(byte >= 0x40)
            {
                if(byte < 0x80) // 0x40..0x7F - most common case
                {
                    // Valid characters in 40..7F: A-Za-z_
                    // Valid bitmask for 40..5F: 01111111111111111111111111100001
                    // Valid bitmask for 60..7F: 01111111111111111111111111100000
                    if(sizeof(unsigned long) == 8)
                    {
                        const unsigned n = sizeof(unsigned long)*8-32;
                        // ^ avoids compiler warning when not 64-bit
                        unsigned long masklow6bits = 1UL << (byte & 0x3F);
                        if(masklow6bits & ~((1UL << 0) | (0x0FUL << (0x1B  ))
                                          | (1UL << n) | (0x1FUL << (0x1B+n))))
                            { ++nameLength; continue; }
                    }
                    else
                    {
                        unsigned masklow5bits = 1 << (byte & 0x1F);
                        if((masklow5bits & ~(1 | (0x1F << 0x1B))) || byte == '_')
                            { ++nameLength; continue; }
                    }
                    break;
                }
                if(byte < 0xF0)
                {
                    if(byte < 0xE0)
                    {
                        if(byte < 0xC2) break; // 0x80..0xC1
                        if(byte == 0xC2 && uptr[nameLength+1]==0xA0) break; // skip nbsp
                        // C2-DF - next common case when >= 0x40
                        // Valid sequence: C2-DF 80-BF
                        if(schar(uptr[nameLength+1]) > schar(0xBF)) break;
                        nameLength += 2;
                        continue;
                    }
                    if(byte == 0xE0) // E0
                    {
                        // Valid sequence: E0 A0-BF 80-BF
                        if((unsigned char)(uptr[nameLength+1] - 0xA0) > (0xBF-0xA0)) break;
                    }
                    else
                    {
                        if(byte == 0xED) break; // ED is invalid
                        // Valid sequence: E1-EC 80-BF 80-BF
                        //            And: EE-EF 80-BF 80-BF
                        if(byte == 0xE2)
                        {
                            // break on various space characters
                            if(uptr[nameLength+1] == 0x80
                            && (schar(uptr[nameLength+2]) <= schar(0x8B)
                            || (uptr[nameLength+2] == 0xAF))) break;
                            if(uptr[nameLength+1] == 0x81
                            && uptr[nameLength+2] == 0x9F) break;
                        } else
                        if(byte == 0xE3 && uptr[nameLength+1] == 0x80
                        && uptr[nameLength+2] == 0x80) break; // this too

                        if(schar(uptr[nameLength+1]) > schar(0xBF)) break;
                    }
                    if(schar(uptr[nameLength+2]) > schar(0xBF)) break;
                    nameLength += 3;
                    continue;
                }
                if(byte == 0xF0) // F0
                {
                    // Valid sequence: F0 90-BF 80-BF 80-BF
                    if((unsigned char)(uptr[nameLength+1] - 0x90) > (0xBF-0x90)) break;
                }
                else
                {
                    if(byte > 0xF4) break; // F5-FF are invalid
                    if(byte == 0xF4) // F4
                    {
                        // Valid sequence: F4 80-8F
                        if(schar(uptr[nameLength+1]) > schar(0x8F)) break;
                    }
                    else
                    {
                        // F1-F3
                        // Valid sequence: F1-F3 80-BF 80-BF 80-BF
                        if(schar(uptr[nameLength+1]) > schar(0xBF)) break;
                    }
                }
                if(schar(uptr[nameLength+2]) > schar(0xBF)) break;
                if(schar(uptr[nameLength+3]) > schar(0xBF)) break;
                nameLength += 4;
                continue;
            }
            if(nameLength > 0)
            {
                if(sizeof(unsigned long) == 8)
                {
                    // Valid bitmask for 00..1F: 00000000000000000000000000000000
                    // Valid bitmask for 20..3F: 00000000000000001111111111000000
                    const unsigned n = sizeof(unsigned long)*8-32;
                    // ^ avoids compiler warning when not 64-bit
                    unsigned long masklow6bits = 1UL << byte;
                    if(masklow6bits & (((1UL << 10)-1UL) << (16+n)))
                        { ++nameLength; continue; }
                }
                else
                {
                    if(byte >= '0' && byte <= '9')
                        { ++nameLength; continue; }
                }
            }
            break;
        }

        /* This function generated with make_function_name_parser.cc */
#define lO l3 lH
#define lN switch(
#define lM l4 lH
#define lL if('i' l5
#define lK 'n' l5
#define lJ 0x80000003U;
#define lI l1 3]={
#define lH case
#define lG 0x80000005U;
#define lF )==0)l0(
#define lE l8 3;}lH
#define lD std::memcmp(uptr+
#define lC l2 3 lF
#define lB lA 1]){lH
#define lA :lN uptr[
#define l9 'a' lB
#define l8 default:l0
#define l7 lG l0 5;}lH
#define l6 <<16)|
#define l5 ==uptr[
#define l4 lJ l0 3;
#define l3 0x80000004U;l0 4;
#define l2 lD 1,tmp,
#define l1 static const char tmp[
#define l0 return
lN
nameLength){lH
2:lL
0]&&'f' l5
1])l0(cIf
l6
0x80000002U;l0
2;lH
3
lA
0]){lH
l9'b':if('s' l5
2])l0(cAbs
l6
lM'r':if('g' l5
2])l0(cArg
l6
l4
lE'c' lB'o' lA
2]){lH's':l0(cCos
l6
lJ
lH't':l0(cCot
l6
lJ
lE's':if('c' l5
2])l0(cCsc
l6
l4
lE'e':if('x' l5
1]&&'p' l5
2])l0(cExp
l6
lM'i':if(lK
1]&&'t' l5
2])l0(cInt
l6
lM'l':if('o' l5
1]&&'g' l5
2])l0(cLog
l6
lM'm' lB'a':if('x' l5
2])l0(cMax
l6
lM'i':if(lK
2])l0(cMin
l6
l4
lE'p':if('o' l5
1]&&'w' l5
2])l0(cPow
l6
lM's' lB'e':if('c' l5
2])l0(cSec
l6
lM'i':if(lK
2])l0(cSin
l6
l4
lE't':if('a' l5
1]&&lK
2])l0(cTan
l6
l4
lE
4
lA
0]){lH
l9'c':if('o' l5
2]&&'s' l5
3])l0(cAcos
l6
lO's':lL
2]&&lK
3])l0(cAsin
l6
lO't':if('a' l5
2]&&lK
3])l0(cAtan
l6
l3
l8
4;}
lH'c' lB'b':if('r' l5
2]&&'t' l5
3])l0(cCbrt
l6
lO'e':lL
2]&&'l' l5
3])l0(cCeil
l6
lO'o' lA
2]){lH'n':if('j' l5
3])l0(cConj
l6
lO's':if('h' l5
3])l0(cCosh
l6
l3
l8
4;}
l8
4;}
lH'e':{lI'x','p','2'}
;if(lC
cExp2
l6
l3}
lH'i':{lI'm','a','g'}
;if(lC
cImag
l6
l3}
lH'l':{lI'o','g','2'}
;if(lC
cLog2
l6
l3}
lH'r':{lI'e','a','l'}
;if(lC
cReal
l6
l3}
lH's' lB'i':if(lK
2]&&'h' l5
3])l0(cSinh
l6
lO'q':if('r' l5
2]&&'t' l5
3])l0(cSqrt
l6
l3
l8
4;}
lH't':{lI'a','n','h'}
;if(lC
cTanh
l6
l3}
l8
4;}
lH
5
lA
0]){lH
l9'c':{lI'o','s','h'}
;if(lD
2,tmp,3
lF
cAcosh
l6
l7's':{lI'i','n','h'}
;if(lD
2,tmp,3
lF
cAsinh
l6
l7't':if('a' l5
2]){if(lK
3]){lN
uptr[4]){lH'2':l0(cAtan2
l6
lG
lH'h':l0(cAtanh
l6
lG
l8
5;}
}
l0
5;}
l0
5;l8
5;}
lH'f':{l1
4]={'l','o','o','r'}
;if(l2
4
lF
cFloor
l6
l7'h':{l1
4]={'y','p','o','t'}
;if(l2
4
lF
cHypot
l6
l7'l':{l1
4]={'o','g','1','0'}
;if(l2
4
lF
cLog10
l6
l7'p':{l1
4]={'o','l','a','r'}
;if(l2
4
lF
cPolar
l6
l7't':{l1
4]={'r','u','n','c'}
;if(l2
4
lF
cTrunc
l6
lG
l0
5;}
l8
5;}
default:break;}
l0
nameLength;
