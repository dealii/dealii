/***************************************************************************\
|* Function Parser for C++ v4.5                                            *|
|*-------------------------------------------------------------------------*|
|* Function optimizer                                                      *|
|*-------------------------------------------------------------------------*|
|* Copyright: Joel Yliluoma                                                *|
|*                                                                         *|
|* This library is distributed under the terms of the                      *|
|* GNU Lesser General Public License version 3.                            *|
|* (See lgpl.txt and gpl.txt for the license text.)                        *|
\***************************************************************************/

/* NOTE:
 This file contains generated code (from the optimizer sources) and is
 not intended to be modified by hand. If you want to modify the optimizer,
 download the development version of the library.
*/

#include "fpconfig.hh"
#ifdef FP_SUPPORT_OPTIMIZER
#include "fparser.hh"
#include "extrasrc/fptypes.hh"
#include "extrasrc/fpaux.hh"
#define l64 y91 a),
#define l54 if(n5 iK3
#define l44 ,{ReplaceParams,
#define l34 cTan iO
#define l24 "Found "
#define l14 ,cPow,
#define l04 stackpos
#define iZ3 sim.nY 1,
#define iY3 ,tree,info
#define iX3 "dup(%u) "
#define iW3 "%d, cost "
#define iV3 "PUSH "i13(
#define iU3 "immed "<<
#define iT3 mFuncParsers
#define iS3 e62{assert
#define iR3 stderr
#define iQ3 sep2=" "
#define iP3 FPHASH_CONST
#define iO3 .SetParamsMove(
#define iN3 cache_needed[
#define iM3 fprintf
#define iL3 ::cout<<"Applying "
#define iK3 HANDLE_UNARY_CONST_FUNC
#define iJ3 1)i61){
#define iI3 c_count
#define iH3 s_count
#define iG3 2)lS 2*
#define iF3 tmp.yH1
#define iE3 tmp2.nJ
#define iD3 ,tZ 2,
#define iC3 (p0 n51&&
#define iB3 else nM
#define iA3 max.val
#define i93 eK2 if(
#define i83 tree c9
#define i73 sim.x61
#define i63 ].swap(
#define i53 codes[b
#define i43 whydump
#define i33 );}case
#define i23 :if(eQ3
#define i13 ;DumpTree
#define i03 :{lX1 r
#define tZ3 nparams
#define tY3 cLog iO
#define tX3 l4 16,1,
#define tW3 l4 0,1,
#define tV3 0x12 nH
#define tU3 ,0,0x16},{
#define tT3 nQ 0,
#define tS3 cAbs nQ
#define tR3 false;}
#define tQ3 nE1 cF1
#define tP3 xF1++b)
#define tO3 =false;
#define tN3 {data->
#define tM3 info.
#define tL3 b.Value)
#define tK3 b.Opcode
#define tJ3 ,tB info
#define tI3 tree xD
#define tH3 xD cMul);
#define tG3 ParamHolder
#define tF3 size()
#define tE3 .second
#define tD3 ]tE3
#define tC3 ].first
#define tB3 Ne_Mask
#define tA3 Gt_Mask
#define t93 Lt_Mask
#define t83 opcode,
#define t73 resize(
#define t63 t33 nC1;
#define t53 xD cond nC
#define t43 );}else
#define t33 lE2 a)
#define t23 AddFrom(
#define t13 max.known
#define t03 =fp_pow(
#define eZ3 (tree)!=
#define eY3 );if(
#define eX3 public:
#define eW3 pclone
#define eV3 Others
#define eU3 --cN1.
#define eT3 cOr,l6
#define eS3 newpow
#define eR3 if(p0 yN
#define eQ3 &*(*x5 n71){
#define eP3 ,lT 1,0,
#define eO3 lE2 1)
#define eN3 (op1==
#define eM3 change
#define eL3 (count
#define eK3 133,2,
#define eJ3 Needs
#define eI3 byteCode
#define eH3 (p1 nC1
#define eG3 lS1 nC==
#define eF3 cLog2by
#define eE3 factor_t
#define eD3 (*x5)[0].info
#define eC3 tree nC1
#define eB3 value1
#define eA3 Finite
#define e93 fp_mod(
#define e83 )const
#define e73 (e83{
#define e63 else{if(
#define e53 xK());nD
#define e43 TreeCountItem
#define e33 if(op==
#define e23 yN val)<
#define e13 p yN val
#define e03 c9 ifp2
#define cZ3 cAbsNot
#define cY3 stackptr
#define cX3 cLog);xL
#define cW3 switch(tM
#define cV3 p1 cK p1);
#define cU3 );p1 c9 ifp1
#define cT3 .empty()
#define cS3 opcodes
#define cR3 did_muli
#define cQ3 xD data.
#define cP3 &Value){
#define cO3 yH const
#define cN3 used[b]
#define cM3 sizeof(
#define cL3 cAbsIf,
#define cK3 cNotNot,
#define cJ3 18,2 eN1
#define cI3 c22 iE,
#define cH3 450998,
#define cG3 cAtan2,
#define cF3 cExp2 nQ
#define cE3 lJ 2},0,
#define cD3 middle2
#define cC3 ::string
#define cB3 fp_equal(
#define cA3 ;if(cB3
#define c93 (iE2 l8 1)xC1
#define c83 FP_GetOpcodeName(
#define c73 default:
#define c63 n31;case
#define c53 tG));nV
#define c43 range nJ3
#define c33 range<yN2
#define c23 range xK
#define c13 Ge0Lt1
#define c03 Gt0Le1
#define yZ3 cAdd lS2
#define yY3 ==cOr tE2
#define yX3 IsLogicalValue y91
#define yW3 cAbsIf)
#define yV3 iterator
#define yU3 begin();
#define yT3 TreeSet
#define yS3 parent
#define yR3 insert(i
#define yQ3 newrel
#define yP3 for iJ nS
#define yO3 ;tree tK2
#define yN3 ;nV l5::
#define yM3 ;for iA1 a=
#define yL3 param.
#define yK3 &param=*
#define yJ3 break;l73
#define yI3 xE lT 2,
#define yH3 ;}static yU1
#define yG3 ::ByteCodeSynth xK
#define yF3 break;}
#define yE3 synth lO2
#define yD3 synth.xM
#define yC3 b_needed
#define yB3 cachepos
#define yA3 half=
#define y93 ,iB,1,l51+1);
#define y83 131,4,1,
#define y73 131,8,1,
#define y63 4,1,2,1,
#define y53 )n31 lH
#define y43 ,1,562 nP
#define y33 ){case
#define y23 (yL3
#define y13 ){switch(
#define y03 (cond yV
#define xZ3 {i31 lE
#define xY3 template lL
#define xX3 iY2 tF3
#define xW3 template lX
#define xV3 ::vector
#define xU3 FindPos(
#define xT3 src_pos
#define xS3 reserve(
#define xR3 nJ1 void
#define xQ3 treeptr
#define xP3 tA1 void
#define xO3 ImmedTag
#define xN3 a,const
#define xM3 RefCount
#define xL3 Birth();
#define xK3 typename
#define xJ3 7168,
#define xI3 cost_t
#define xH3 fpdata
#define xG3 middle
#define xF3 };enum
#define xE3 tG1 1))){
#define xD3 sqrt_cost
#define xC3 const int
#define xB3 maxValue1
#define xA3 minValue1
#define x93 maxValue0
#define x83 minValue0
#define x73 ValueType
#define x63 yN n3 0),
#define x53 =true;lI2
#define x43 yK known
#define x33 yK n3 0),
#define x23 abs_mul
#define x13 l8 a));
#define x03 pos_set
#define nZ3 ContainsOtherCandidates
#define nY3 goto cY
#define nX3 l14 xV1
#define nW3 t2=!t2;}
#define nV3 bool t2 tO3
#define nU3 (yK val)
#define nT3 result
#define nS3 .tL1 n]
#define nR3 1)tG1 1));
#define nQ3 nT3 yN known
#define nP3 nT3 n51
#define nO3 ),child);
#define nN3 mulgroup)
#define nM3 nP3 tO3 if(
#define nL3 nT3 yN val
#define nK3 nT3 tF2
#define nJ3 xK nT3
#define nI3 const std::eT
#define nH3 const char*
#define nG3 nE3 IsNever nV2
#define nF3 nE3 IsAlways;if(
#define nE3 ))return
#define nD3 y3 cAdd);
#define nC3 subtree
#define nB3 invtree
#define nA3 MakeHash(
#define n93 yR false;
#define n83 rulenumit
#define n73 cAnd l3
#define n63 282870 nP
#define n53 ,tree nX
#define n43 ,1,538 nP
#define n33 cAnd,l6
#define n23 MakeEqual
#define n13 MakeTrue,
#define n03 newbase
#define lZ3 branch1op
#define lY3 branch2op
#define lX3 overlap
#define lW3 nE1 c9 y6
#define lV3 truth_b
#define lU3 truth_a
#define lT3 if y91 0)
#define lS3 found_dup
#define lR3 (mulgroup
#define lQ3 cH1 eJ&
#define lP3 rangeutil
#define lO3 Plan_Has(
#define lN3 StackMax)
#define lM3 ;if(half
#define lL3 ;}void
#define lK3 cGreater,
#define lJ3 ,cGreater
#define lI3 )lL3
#define lH3 nT3 xQ
#define lG3 const xA2
#define lF3 namespace
#define lE3 ByteCode[
#define lD3 inverted
#define lC3 IsNever:
#define lB3 .known&&
#define lA3 ==cNot||
#define l93 iftree
#define l83 if(i01==
#define l73 }switch(
#define l63 depcodes
#define l53 explicit
#define l43 l2 16 eI
#define l33 l91 2,1,
#define l23 cCosh nQ
#define l13 n91 1,
#define l03 VarBegin
#define iZ2 Params[a].
#define iY2 Params.
#define iX2 ].data);
#define iW2 PlusInf
#define iV2 Value_t
#define iU2 ;iV2
#define iT2 (nT3
#define iS2 synth);
#define iR2 lJ2 apos==
#define iQ2 default_function_handling
#define iP2 xD i01);
#define iO2 tE3)
#define iN2 begin(),
#define iM2 cond_add
#define iL2 cond_mul
#define iK2 cond_and
#define iJ2 ){eJ r;r xD
#define iI2 func lL1
#define iH2 const eH
#define iG2 bool eA1
#define iF2 unsigned
#define iE2 leaf1
#define iD2 costree
#define iC2 sintree
#define iB2 leaf_count
#define iA2 sub_params
#define i92 tG1-1)))xC
#define i82 printf(
#define i72 swap(tmp);
#define i62 cbrt_count
#define i52 sqrt_count
#define i42 ,std::cout
#define i32 pcall_tree
#define i22 after_powi
#define i12 GetHash().
#define i02 *x5 n71,info
#define tZ2 eQ1=0;a eX2
#define tY2 eQ1;if t51
#define tX2 grammar
#define tW2 cCos iO
#define tV2 l14 l2 0,2,
#define tU2 cEqual eW1
#define tT2 0x12},{{3,
#define tS2 cNeg,lT 1,
#define tR2 ,bool abs)
#define tQ2 sim.Eat(
#define tP2 .what nU1
#define tO2 tJ2 nT3(
#define tN2 ,n42 l2 c22
#define tM2 if y91 1)i61&&
#define tL2 ;range.xE2
#define tK2 .DelParam(
#define tJ2 ){iV2
#define tI2 nR tJ2 tmp=
#define tH2 );nD lC
#define tG2 iC3 p0 tF2>=x22
#define tF2 .min.val
#define tE2 )?0:1))l7
#define tD2 .n_int_sqrt
#define tC2 ;}else{x5=new
#define tB2 (half&63)-1;
#define tA2 ),0},{
#define t92 ),Value
#define t82 data;data.
#define t72 MakeNEqual
#define t62 (exponent
#define t52 Dump(std::
#define t42 isInteger(
#define t32 Comparison
#define t22 needs_flip
#define t12 value]
#define t02 ~size_t(0)
#define eZ2 xJ1 xS+1);
#define eY2 TopLevel)
#define eX2 <xV;++a)
#define eW2 }}return
#define eV2 vector<bool>
#define eU2 *)&*start_at;
#define eT2 Rule&rule,
#define eS2 >::res,b8<
#define eR2 continue;}
#define eQ2 c9 tree);
#define eP2 mul_item
#define eO2 innersub
#define eN2 cbrt_cost
#define eM2 best_cost
#define eL2 fp_min(yL
#define eK2 tree t31);}
#define eJ2 tree))cZ
#define eI2 condition
#define eH2 per_item
#define eG2 item_type
#define eF2 first2
#define eE2 l4 18,1,
#define eD2 lJ 1},0,
#define eC2 Decision
#define eB2 not_tree
#define eA2 group_by
#define e92 exponent=
#define e82 ->second
#define e72 :xL iV2(
#define e62 eJ&tree)
#define e52 targetpos
#define e42 ParamSpec
#define e32 rhs.hash2;}
#define e22 rhs.hash1
#define e12 struct
#define e02 Forget()
#define cZ2 &&cond e8))
#define cY2 source_tree
#define cX2 nC==cLog2&&
#define cW2 ,eL1 0x1 nH
#define cV2 )lS 3*3*
#define cU2 (tK==eA3&&
#define cT2 );yF3
#define cS2 ;iB.Remember(
#define cR2 <tW,xI3>
#define cQ2 );sim.nY 2,
#define cP2 )continue
#define cO2 CodeTree lW
#define cN2 p1_evenness
#define cM2 isNegative(
#define cL2 yN known&&(
#define cK2 neg_set
#define cJ2 cNop,cNop}}
#define cI2 cTanh,cNop,
#define cH2 NewHash
#define cG2 >e12 cA<
#define cF2 matches
#define cE2 .match_tree
#define cD2 cH1 void*)&
#define cC2 cSin iO
#define cB2 cTan nQ
#define cA2 ,cLog nQ
#define c92 cCos nQ
#define c82 long value
#define c72 (std::move(
#define c62 (xX l8 a)xF
#define c52 negated
#define c42 Specializer
#define c32 params
#define c22 18,2,
#define c12 coshtree
#define c02 sinhtree
#define yZ2 best_score
#define yY2 mulvalue
#define yX2 pow_item
#define yW2 subgroup
#define yV2 nC==cPow&&tO
#define yU2 PowiResult
#define yT2 ))y53
#define yS2 0));c23
#define yR2 maxValue
#define yQ2 minValue
#define yP2 div_tree
#define yO2 pow_tree
#define yN2 iV2 nT
#define yM2 preserve
#define yL2 PullResult()
#define yK2 dup_or_fetch
#define yJ2 nominator]
#define yI2 test_order
#define yH2 parampair
#define yG2 yH2,
#define yF2 .param_count
#define yE2 minimum_need
#define yD2 shift(index)
#define yC2 {std::cout<<
#define yB2 rulenumber
#define yA2 l2 16,2,
#define y92 cTanh nQ
#define y82 l14 l2 2,2,
#define y72 ,cPow iO
#define y62 cSinh nQ
#define y52 cIf,n91 3,
#define y42 cInv,lT 1,
#define y32 constraints=
#define y22 GetDepth()
#define y12 factor_immed
#define y02 changes
#define xZ2 c9 cond l8
#define xY2 Become y91
#define xX2 for(xK3
#define xW2 exp_diff
#define xV2 ExponentInfo
#define xU2 lower_bound(
#define xT2 factor
#define xS2 is_logical
#define xR2 newrel_and
#define xQ2 )){data x7 lN
#define xP2 Suboptimal
#define xO2 eK[c eA
#define xN2 res_stackpos
#define xM2 half_pos
#define xL2 ifdata.ofs
#define xK2 i8 push_back(
#define xJ2 n31 true;}
#define xI2 >>1)):(
#define xH2 CodeTreeData
#define xG2 exponent)
#define xF2 yN known)
#define xE2 multiply(
#define xD2 )y3 cPow)
#define xC2 var_trees
#define xB2 nB OPCODE
#define xA2 CodeTree&
#define x92 parent_opcode
#define x82 GetParam(a eS
#define x72 {cW start_at;
#define x62 changed=true;
#define x52 log2_exponent
#define x42 iV2(2)));
#define x32 lR,tH)tH2
#define x22 iV2(0.0))
#define x12 dup_fetch_pos
#define x02 IsNever cL lC
#define nZ2 cSin nQ
#define nY2 Value_EvenInt
#define nX2 AddCollection
#define nW2 ConditionType
#define nV2 n31 Unknown;}
#define nU2 1 y7 i8 size(
#define nT2 (iF2
#define nS2 iA|nT2)
#define nR2 SpecialOpcode
#define nQ2 =i e82.
#define nP2 yN known&&p
#define nO2 assimilated
#define nN2 denominator
#define nM2 fraction
#define nL2 DUP_BOTH();
#define nK2 0x80000000u
#define nJ2 IsDescendantOf
#define nI2 TreeCounts
#define nH2 SetOpcode(
#define nG2 found_log2
#define nF2 div_params
#define nE2 immed_sum
#define nD2 lE3++IP]
#define nC2 sim n92
#define nB2 ;sim.Push(
#define nA2 ,cMul l3
#define n92 .Eat(1,
#define n82 OPCODE(opcode)
#define n72 FactorStack xK
#define n62 Rehash(false);
#define n52 IsAlways cL lC
#define n42 cEqual,
#define n32 cNotNot nQ
#define n22 cNot nQ
#define n12 DumpHashesFrom
#define n02 replacing_slot
#define lZ2 RefParams
#define lY2 if_always[
#define lX2 WhatDoWhenCase
#define lW2 exponent_immed
#define lV2 new_base_immed
#define lU2 base_immed
#define lT2 );std::cout<<
#define lS2 ||op1==
#define lR2 l8 0));
#define lQ2 data[a tD3
#define lP2 if(newrel_or==
#define lO2 .AddOperation(
#define lN2 DUP_ONE(apos);
#define lM2 flipped
#define lL2 .UseGetNeeded(
#define lK2 ,fp_max(yL);
#define lJ2 else if(
#define lI2 lJ2!nT3
#define lH2 synth.PushImmed(
#define lG2 (*x5)[a].info
#define lF2 ,PowiCache&iB,
#define lE2 tree l8
#define lD2 <<tree.i12
#define lC2 const eJ&
#define lB2 ,l1 0x0},{{3,
#define lA2 ,0,0x4},{{1,
#define l92 {case IsAlways:
#define l82 e3 2,131,
#define l72 (IfData&ifdata
#define l62 [xS-1-offset].
#define l52 lE3 a
#define l42 Immed.tF3);
#define l32 OptimizedUsing
#define l22 Var_or_Funcno
#define l12 l22;
#define l02 GetParams(
#define iZ1 crc32_t
#define iY1 signed_chain
#define iX1 MinusInf
#define iW1 n_immeds
#define iV1 stack.tF3
#define iU1 FindClone(xI
#define iT1 lE3 IP]
#define iS1 GetOpcode())
#define iR1 needs_rehash
#define iQ1 AnyWhere_Rec
#define iP1 ~iF2(0)
#define iO1 41,42,43,44,
#define iN1 p1_logical_b
#define iM1 p0_logical_b
#define iL1 p1_logical_a
#define iK1 p0_logical_a
#define iJ1 *const func)
#define iI1 synth.DoDup(
#define iH1 cache_needed
#define iG1 e3 2,1,e3 2,
#define iF1 treelist
#define iE1 has_bad_balance
#define iD1 {if(GetOpcode()
#define iC1 if(remaining[a])
#define iB1 (*x5 n71=r.specs;if(r.found){
#define iA1 (size_t
#define i91 }yF3 case
#define i81 for iA1 b=0;b<
#define i71 for iA1 a y0
#define i61 .IsImmed()
#define i51 i61){if(
#define i41 eE3 xT2
#define i31 template<
#define i21 divgroup
#define i11 2*2*2)lS 3
#define i01 tree nC
#define tZ1 xW3 void
#define tY1 yD3 1
#define tX1 ,iB,eU,iS2
#define tW1 fp_abs(iA3))
#define tV1 fp_abs(min.val)
#define tU1 tJ 2},0,0x0},{{
#define tT1 Oneness_NotOne|
#define tS1 Value_IsInteger
#define tR1 Constness_Const
#define tQ1 ;xI iO3
#define tP1 l32(
#define tO1 reltype
#define tN1 SequenceOpcodes
#define tM1 synth.Find(
#define tL1 sep_list[
#define tK1 TreeCountType xK
#define tJ1 nI2.erase(cs_it);
#define tI1 nE3 true
#define tH1 yK set(fp_floor);
#define tG1 ,iV2(
#define tF1 );exponent
#define tE1 back().thenbranch
#define tD1 grammar_rules[*r]
#define tC1 goto fail;}
#define tB1 l1 0x4 nH
#define tA1 lZ2);
#define t91 IsDefined()){
#define t81 eJ tmp;tmp xD
#define t71 >(iV2(1),
#define t61 x22{xO
#define t51 (&*start_at){x5=(
#define t41 DelParam(a);}
#define t31 tK2 a
#define t21 a)i61)
#define t11 n82);
#define t01 0.5)xD2;lC
#define eZ1 std::cout<<"POP "
#define eY1 stack[iV1-
#define eX1 stack.push_back(
#define eW1 ,l0 2,
#define eV1 CollectionSet xK
#define eU1 xK())y3 cMul);lC
#define eT1 .Rehash();
#define eS1 ParsePowiMuli
#define eR1 MaxChildDepth
#define eQ1 iF2 a
#define eP1 std::pair<It,It>
#define eO1 cLess,l2
#define eN1 ,n63 281856,
#define eM1 cNEqual,
#define eL1 cPow,lA
#define eK1 Sign_Negative
#define eJ1 Value_Logical
#define eI1 yD MakeFalse,{l5
#define eH1 new_factor_immed
#define eG1 occurance_pos
#define eF1 exponent_hash
#define eE1 exponent_list
#define eD1 CollectMulGroup(
#define eC1 source_set
#define eB1 exponent,yT3
#define eA1 operator
#define e91 FindAndDup(tree);
#define e81 ParamSpec_Extract
#define e71 retry_anyparams_3
#define e61 retry_anyparams_2
#define e51 needlist_cached_t
#define e41 lJ 2}lA2
#define e31 lJ 1}lA2
#define e21 CodeTreeImmed xK(
#define e11 by_float_exponent
#define e01 fp_equal t62
#define cZ1 new_exp
#define cY1 end()&&i->first==
#define cX1 return BecomeZero;
#define cW1 return BecomeOne;
#define cV1 if(lQ.tF3<=n1)
#define cU1 addgroup
#define cT1 found_log2by
#define cS1 nC==cZ3)
#define cR1 if(keep_powi
#define cQ1 l22)
#define cP1 branch1_backup
#define cO1 branch2_backup
#define cN1 NeedList
#define cM1 tree.SetParam(
#define cL1 exponent_map
#define cK1 plain_set
#define cJ1 rangehalf
#define cI1 ;synth.StackTopIs(
#define cH1 (const
#define cG1 cH1 iV2&
#define cF1 eT1 tI3 op1);tree.DelParams()
#define cE1 .GetParamCount()==
#define cD1 LightWeight(
#define cC1 if(value
#define cB1 xW3 yX
#define cA1 xW3 static
#define c91 m yN set(fp_ceil);tV
#define c81 m yN val
#define c71 m yN known
#define c61 )yH2 tE3;
#define c51 ::MakeTrue
#define c41 should_regenerate=true;
#define c31 should_regenerate,
#define c21 Collection
#define c11 RelationshipResult
#define c01 Subdivide_Combine(
#define yZ1 e83 yR
#define yY1 rhs yZ1 hash1
#define yX1 PositionalParams,0
#define yW1 best_sep_factor
#define yV1 needlist_cached
#define yU1 inline iF2
#define yT1 252421 nP 24830
#define yS1 t83 bool pad
#define yR1 MakesInteger(
#define yQ1 mulgroup.
#define yP1 .AddParamMove(
#define yO1 =comp.AddItem(atree
#define yN1 tree i61 cL
#define yM1 ByteCodeSynth xK&synth)
#define yL1 lF3 FPoptimizer_Optimize
#define yK1 (long double)
#define yJ1 const iV2&value
#define yI1 best_sep_cost
#define yH1 AddParamMove(tree);
#define yG1 MultiplicationRange
#define yF1 pihalf_limits
#define yE1 n_stacked
#define yD1 cH2.hash1
#define yC1 AnyParams_Rec
#define yB1 ):e2(),std xV3<
#define yA1 i13 cJ
#define y91 (lE2
#define y81 Become(value l8 0))
#define y71 always_sincostan
#define y61 Recheck_RefCount_Div
#define y51 Recheck_RefCount_Mul
#define y41 mulgroup;mulgroup xD
#define y31 MultiplyAndMakeLong(
#define y21 iV2(0)
#define y11 covers_plus1
#define y01 if(synth.FindAndDup(
#define xZ1 SynthesizeParam(
#define xY1 grammar_func
#define xX1 cOr l3 16,1,
#define xW1 252180 nP 281854
#define xV1 l2 0,2,165888 nP
#define xU1 l1 tV3
#define xT1 Modulo_Radians},
#define xS1 PositionType
#define xR1 CollectionResult
#define xQ1 ;m yN template set_if<
#define xP1 ByteCode,size_t&IP,size_t limit,size_t y1
#define xO1 iT2))break;nT3*=
#define xN1 xW3 bool
#define xM1 const_offset
#define xL1 inline TriTruthValue
#define xK1 stacktop_desired
#define xJ1 SetStackTop(
#define xI1 }inline
#define xH1 FPoptimizer_ByteCode
#define xG1 1)?(poly^(
#define xF1 GetParamCount();
#define xE1 public e2,public std xV3<
#define xD1 y21)
#define xC1 xF leaf2 l8
#define xB1 c22 529654 nP
#define xA1 cond_type
#define x91 fphash_value_t
#define x81 Recheck_RefCount_RDiv
#define x71 cMul);tmp.nJ 0));tmp.
#define x61 SwapLastTwoInStack();
#define x51 fPExponentIsTooLarge(
#define x41 CollectMulGroup_Item(
#define x31 pair<iV2,yT3>
#define x21 nL xJ1 xS-1);
#define x11 covers_full_cycle
#define x01 AssembleSequence(
#define nZ1 <<std::dec<<")";}
#define nY1 yD MakeNotP1,l5::
#define nX1 yD MakeNotP0,l5::
#define nW1 {DataP slot_holder(xY[
#define nV1 :return p tF2
#define nU1 !=xB)if(TestCase(
#define nT1 &&IsLogicalValue(
#define nS1 std::pair<T1,T2>&
#define nR1 i31 xK3
#define nQ1 has_good_balance_found
#define nP1 n_occurrences
#define nO1 found_log2_on_exponent
#define nN1 covers_minus1
#define nM1 needs_resynth
#define nL1 immed_product
#define nK1 ,2,1 xN if(found[data.
#define nJ1 ;xW3
#define nI1 yJ3 bitmask&
#define nH1 Sign_Positive
#define nG1 SetParamMove(
#define nF1 CodeTreeImmed(iV2(
#define nE1 changed_if
#define nD1 0,xG2;DelParam(1);
#define nC1 .GetImmed()
#define nB1 for iA1 a=0;a<c8++a)
#define nA1 iF2 index
#define n91 lA 0x4},{{
#define n81 iI1 found[data.
#define n71 )[a].start_at
#define n61 ,cIf,l0 3,
#define n51 .min.known
#define n41 .xF1 a-->0;)if(
#define n31 ;return
#define n21 n31 xP2;}
#define n11 >::Optimize(){}
#define n01 n_as_tanh_param
#define lZ1 opposite=
#define lY1 x91(
#define lX1 MatchResultType
#define lW1 needs_sincos
#define lV1 resulting_exponent
#define lU1 val):Value(Value::
#define lT1 Unknown:c73;}
#define lS1 GetParam(a)
#define lR1 inverse_nominator]
#define lQ1 cMul l3 0,1,
#define lP1 tree yP1
#define lO1 AddFunctionOpcode(
#define lN1 SetParams(l02));
#define lM1 o<<"("<<std::hex<<data.
#define lL1 (val);else*this=model;}
#define lK1 IfBalanceGood(
#define lJ1 n_as_tan_param
#define lI1 changed_exponent
#define lH1 &&e13<iV2(
#define lG1 inverse_denominator
#define lF1 ;cH2.hash2+=
#define lE1 xK(rule.repl_param_list,
#define lD1 retry_positionalparams_2
#define lC1 situation_flags&
#define lB1 518 nP 400412,
#define lA1 7168 nP 401798
#define l91 }},{ProduceNewTree,
#define l81 data.subfunc_opcode
#define l71 i61 tJ2
#define l61 CopyOnWrite();
#define l51 recursioncount
#define l41 PlanNtimesCache(
#define l31 >){int mStackPtr=0;
#define l21 FPoptimizer_Grammar
#define l11 AddOperation(cInv,1,1 xN}
#define l01 GetPositivityInfo eZ3
#define iZ ParamSpec_SubFunctionData
#define iY iA1 a=c8 a-->0;)
#define iX PositionalParams_Rec
#define iW yD MakeNotNotP1,l5::
#define iV ,cMul x4
#define iU eT1 tI3 iE2 nC);tree.
#define iT yD MakeNotNotP0,l5::
#define iS lE3 xL2+
#define iR DumpTreeWithIndent(*this);
#define iQ switch(type y33 cond_or:
#define iP ;if(eY2 tM3 SaveMatchedParamIndex(
#define iO ,l4 2,1,
#define iN CalculateResultBoundaries(
#define iM i31 iF2 Compare>
#define iL cE3 0x0},{{1,
#define iK edited_powgroup
#define iJ iA1 a=xF1 a
#define iI has_unknown_max
#define iH has_unknown_min
#define iG static const c23
#define iF synthed_tree
#define iE 408964 nP 24963
#define iD SelectedParams,0},0,0x0},{{
#define iC collections
#define iB cache
#define iA ;xK2 nK2
#define i9 )iA);
#define i8 ByteCode.
#define i7 ;pow xD cLog);tree tH3
#define i6 goto ReplaceTreeWithOne;case
#define i5 ]);yE3
#define i4 !=xB)return lY2
#define i3 e11.data
#define i2 l53 xH2(
#define i1 needs_sinhcosh
#define i0 cAdd x4
#define tZ cAdd l3 0,
#define tY tA2 iV2(
#define tX xW3 nA
#define tW int_exponent_t
#define tV return m;}case
#define tU MakeFalse,l5::
#define tT lT2 std::endl;DumpHashes(
#define tS ,ByteCode,IP,limit,y1,stack);
#define tR matched_params
#define tQ [n1 tC3=true;lQ[n1 tD3
#define tP l21::Grammar*
#define tO powgroup l8
#define tN );p2 cK p2);tI3 l93 nC);cZ}
#define tM GetLogicalValue y91
#define tL t02&&found[data.
#define tK valueType
#define tJ xE AnyParams,
#define tI =iN lE2
#define tH eO3 nC1
#define tG (tH)
#define tF nF1(
#define tE has_mulgroups_remaining
#define tD by_exponent
#define tC const iZ
#define tB MatchInfo xK&
#define tA Rehash();iA2.push_back(
#define t9 best_factor
#define t8 RootPowerTable xK::RootPowers[
#define t7 MatchPositionSpec_AnyParams xK
#define t6 cLessOrEq,l2
#define t5 lF3 FPoptimizer_CodeTree
#define t4 n_as_sinh_param
#define t3 n_as_cosh_param
#define t2 is_signed
#define t1 result_positivity
#define t0 yN known tO3
#define eZ biggest_minimum
#define eY 142455 nP 141449,
#define eX cond_tree
#define eW else_tree
#define eV then_tree
#define eU sequencing
#define eT string c83
#define eS );bool needs_cow=GetRefCount()>1;
#define eR nM iV2(-nR3
#define eQ {AdoptChildrenWithSameOpcode(tree);
#define eP lJ3,l2 c22
#define eO ;AddParamMove(
#define eN yN known&&p0 yN val<=fp_const_negativezero xK())
#define eM tree.GetParamCount()
#define eL :goto ReplaceTreeWithZero;case
#define eK relationships
#define eJ CodeTree xK
#define eI ,2,122999 nP 139399,
#define eH std xV3<eJ>
#define eG if_stack
#define eF (l02));yQ1 Rehash();
#define eE n_as_sin_param
#define eD n_as_cos_param
#define eC PowiResolver::
#define eB cIf,tW3
#define eA ].relationship
#define e9 PACKED_GRAMMAR_ATTRIBUTE;
#define e8 .BalanceGood
#define e7 AddParamMove(yW2
#define e6 back().endif_location
#define e5 x91 key
#define e4 AddParamMove(mul);
#define e3 130,1,
#define e2 MatchPositionSpecBase
#define e1 l53 CodeTree(
#define e0 smallest_maximum
#define cZ goto redo;
#define cY ReplaceTreeWithParam0;
#define cX factor_needs_rehashing
#define cW MatchPositionSpecBaseP
#define cV xK3 tK1::yV3
#define cU e81 xK(nN.param_list,
#define cT 243,244,245,246,249,250,251,253,255,256,257,258,259}};}
#define cS ];};extern"C"{
#define cR 79,122,123,160,161,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,
#define cQ 27,28,29,30,31,32,33,35,36,
#define cP const ParamSpec_SubFunction
#define cO const ParamSpec_ParamHolder
#define cN otherhalf
#define cM StackState
#define cL )return false;
#define cK eT1 lP1
#define cJ (tree lT2"\n";
#define cI const SequenceOpCode xK
#define cH paramholder_matches[nZ]
#define cG MatchPositionSpec_PositionalParams xK
#define cF lC2 tree,std::ostream&o
#define cE paramholder_matches.
#define cD iV2(1.5)*fp_const_pi xK()
#define cC xB,l5::Never yD xB,l5::Never}
#define cB CalculatePowiFactorCost(
#define cA ImmedHashGenerator
#define c9 .AddParam(
#define c8 eM;
#define c7 ::map<fphash_t,std::set<std cC3> >
#define c6 ComparisonSetBase::
#define c5 AddParamMove(comp.cK1[a].value);
#define c4 ,yA2 528503 nP 24713,
#define c3 T1,xK3 T2>inline iG2()(
#define c2 has_nonlogical_values
#define c1 from_logical_context)
#define c0 AnyParams,0 l91
#define yZ for iA1 a=xX.xF1 a-->0;)
#define yY POWI_CACHE_SIZE
#define yX static inline eJ
#define yW ++IP;eR2 if(iT1==cS3.
#define yV .FoundChild
#define yU BalanceResultType
#define yT xM3(0),Opcode(
#define yS );void lO1 iF2 t83 c42<
#define yR {return
#define yQ const yR data->
#define yP +=fp_const_twopi xK();
#define yO fp_const_twopi xK()eY3
#define yN .max.
#define yM eJ tmp,tmp2;tmp2 xD
#define yL fp_sin(min),fp_sin(max))
#define yK m.min.
#define yJ lO2 GetOpcode(),
#define yI for iA1 a=0;a<xF1++a){if(
#define yH static void nA3 nB fphash_t&cH2,
#define yG MatchPositionSpec_AnyWhere
#define yF if y23 data.match_type==
#define yE void OutFloatHex(std::ostream&o,
#define yD },{l5::
#define yC AddParam(CodeTreeImmed(
#define yB cGreaterOrEq,
#define yA ,xK3 eJ::
#define y9 xK model=cJ1 xK()){if(known
#define y8 AssembleSequence_Subdivide(
#define y7 ]=nK2|iF2(
#define y6 branch2
#define y5 (n83 r=range.first;r!=range tE3;++r){
#define y4 !=t02){n81
#define y3 ;tQ2 2,
#define y2 iF2 c;iF2 short l[
#define y1 factor_stack_base
#define y0 =0;a<yS3.xF1++a)if(
#define xZ ,eL1 0x4 nH
#define xY data->Params
#define xX branch1
#define xW );eR2 if(list.first nC1==iV2(
#define xV nN yF2
#define xU =fp_cosh nU3;c81=fp_cosh(c81);
#define xT {nI2.erase(i);eR2
#define xS StackTop
#define xR FPOPT_autoptr
#define xQ +=nT3 n31 nT3;}xW3 inline iV2
#define xP int_exponent
#define xO tree.ReplaceWithImmed(
#define xN )cI1*this)n31;}
#define xM GetStackTop()-
#define xL sim.AddConst(
#define xK <iV2>
#define xJ .SetParam(0,l93 lR2 eJ p1;p1 xD
#define xI newnode
#define xH has_highlevel_opcodes
#define xG eD2 0x0},{{
#define xF .IsIdenticalTo(
#define xE ,cAdd,
#define xD .nH2
#define xC {if(needs_cow){l61 goto
#define xB Unchanged
#define xA cF=std::cout
#define x9 best_selected_sep
#define x8 i31>void FunctionParserBase<
#define x7 ->Recalculate_Hash_NoRecursion();}
#define x6 xF1++a)if(ApplyGrammar(tX2,t33,
#define x5 position
#define x4 l3 2,1,
#define x3 ;iE3 0));tmp xD cInv);tmp yP1 tmp2)n31
#define x2 nB1{c23
#define x1 std xV3<CodeTree>
#define x0 TestImmedConstraints y23 constraints,tree)cL
#define nZ paramholder_index
#define nY x61 tQ2
#define nX )){tree.FixIncompleteHashes();}
#define nW );cR1){nC2 cInv cT2 xL-1 xD2;lC
#define nV return true;case
#define nU occurance_counts
#define nT >p tI a)eY3 p.
#define nS -->0;){lC2 powgroup=lS1;if(powgroup
#define nR lT3 i61
#define nQ ,l0 1,
#define nP ,{2,
#define nO const FPoptimizer_CodeTree::eJ&tree
#define nN model_tree
#define nM return c23(
#define nL ){using lF3 FUNCTIONPARSERTYPES;
#define nK eH&lZ2
#define nJ AddParam y91
#define nI ConstantFolding_LogicCommon(tree,c6
#define nH },{{2,
#define nG nR1 Ref>inline void xR<Ref>::
#define nF AnyParams,1},0,0x0},{{
#define nE ):data(new xH2 xK(
#define nD goto do_return;}
#define nC .GetOpcode()
#define nB FUNCTIONPARSERTYPES::
#define nA xH2 xK::xH2(
#define n9 b;}};i31>e12 Comp<nB
#define n8 l22(),Params(),Hash(),Depth(1),tP1 0){}
#define n7 SynthesizeByteCode(iS2
#define n6 while(ApplyGrammar(cD2
#define n5 GetIntegerInfo y91 0))==IsAlways)nY3
#define n4 ;lP1 nE1)xJ2
#define n3 template set_if<cGreater>(iV2(
#define n2 DumpParams xK y23 data.param_list,yL3 data yF2,o);
#define n1 restholder_index
#define n0 eJ exponent;exponent xD cMul tF1 c9
#define lZ lR eY3 fp_nequal(tmp,xD1){xO iV2(1)/tmp);nD}lC
#define lY :if(ParamComparer xK()(Params[1],Params[0])){std::swap(Params[0],Params[1]);Opcode=
#define lX <xK3 iV2>
#define lW xK tmp;tmp xD cPow);tmp.nJ 0));tmp.yC iV2(
#define lV tR1,0x0},
#define lU AddParamMove(pow l8 1));pow tK2 1);pow eT1 tree.nG1 0,pow);goto NowWeAreMulGroup;}
#define lT GroupFunction,0},lV{{
#define lS tG1 1)/iV2(
#define lR lE2 0)nC1
#define lQ restholder_matches
#define lP yD1|=key;x91 crc=(key>>10)|(key<<(64-10))lF1((~lY1 crc))*3)^1234567;}};
#define lO nE1;nE1 iP2 nE1.AddParamMove y91 0));nE1 c9 xX l8
#define lN xW3 eJ::CodeTree(
#define lM cM1 0,lE2 0)lR2 cM1 1,CodeTreeImmed(
#define lL lX void ByteCodeSynth xK::lO1 iF2 t83 c42<
#define lK cMul,lT 2,
#define lJ cMul,AnyParams,
#define lI y91 0)i61&&eO3 i61){xO
#define lH iN tmp i33
#define lG :eM3=comp.AddRelationship(atree l8 0),atree l8 1),c6
#define lF cPow,l0 2
#define lE xK3 iV2>inline iG2()cG1 xN3 iV2&b)yR a
#define lD {c23 m tI 0));
#define lC break;case
#define lB tZ1 eJ::
#define lA yX1},0,
#define l9 l1 0x0 nH
#define l8 .GetParam(
#define l7 ;eJ nE1;nE1 iP2 nE1 iO3 tree.l02));nE1 eT1 tI3
#define l6 SelectedParams,0},0,0x0 nH
#define l5 RangeComparisonData
#define l4 yX1 l91
#define l3 ,AnyParams,0}}l44
#define l2 yX1}}l44
#define l1 cMul,SelectedParams,0},0,
#define l0 lA 0x0},{{
#ifdef _MSC_VER
typedef
iF2
int
iZ1;
#else
#include <stdint.h>
typedef
uint_least32_t
iZ1;
#endif
lF3
crc32{enum{startvalue=0xFFFFFFFFUL,poly=0xEDB88320UL}
;i31
iZ1
crc>e12
b8{enum{b1=(crc&xG1
crc
xI2
crc>>1),b2=(b1&xG1
b1
xI2
b1>>1),b3=(b2&xG1
b2
xI2
b2>>1),b4=(b3&xG1
b3
xI2
b3>>1),b5=(b4&xG1
b4
xI2
b4>>1),b6=(b5&xG1
b5
xI2
b5>>1),b7=(b6&xG1
b6
xI2
b6>>1),res=(b7&xG1
b7
xI2
b7>>1)}
;}
;inline
iZ1
update(iZ1
crc,iF2
b){
#define B4(n) b8<n eS2 n+1 eS2 n+2 eS2 n+3>::res
#define R(n) B4(n),B4(n+4),B4(n+8),B4(n+12)
static
const
iZ1
table[256]={R(0x00),R(0x10),R(0x20),R(0x30),R(0x40),R(0x50),R(0x60),R(0x70),R(0x80),R(0x90),R(0xA0),R(0xB0),R(0xC0),R(0xD0),R(0xE0),R(0xF0)}
;
#undef R
#undef B4
return((crc>>8))^table[(crc^b)&0xFF];xI1
iZ1
calc_upd(iZ1
c,const
iF2
char*buf,size_t
size){iZ1
value=c;for
iA1
p=0;p<size;++p)value=update(value,buf[p])n31
value;xI1
iZ1
calc
cH1
iF2
char*buf,size_t
size)yR
calc_upd(startvalue,buf,size);}
}
#ifndef FPOptimizerAutoPtrHH
#define FPOptimizerAutoPtrHH
nR1
Ref>class
xR{eX3
xR():p(0){}
xR(Ref*b):p(b){xL3}
xR
cH1
xR&b):p(b.p){xL3
xI1
Ref&eA1*(yZ1*p;xI1
Ref*eA1->(yZ1
p;}
xR&eA1=(Ref*b){Set(b)n31*this;}
xR&eA1=cH1
xR&b){Set(b.p)n31*this;}
#ifdef __GXX_EXPERIMENTAL_CXX0X__
xR(xR&&b):p(b.p){b.p=0;}
xR&eA1=(xR&&b){if(p!=b.p){e02;p=b.p;b.p=0;}
return*this;}
#endif
~xR(){e02
lL3
UnsafeSetP(Ref*newp){p=newp
lL3
swap(xR<Ref>&b){Ref*tmp=p;p=b.p;b.p=tmp;}
private:inline
static
void
Have(Ref*p2);inline
void
e02;inline
void
xL3
inline
void
Set(Ref*p2);private:Ref*p;}
;nG
e02{if(!p)return;p->xM3-=1;if(!p->xM3)delete
p;}
nG
Have(Ref*p2){if(p2)++(p2->xM3);}
nG
Birth(){Have(p);}
nG
Set(Ref*p2){Have(p2);e02;p=p2;}
#endif
#include <utility>
e12
Compare2ndRev{nR1
T>inline
iG2()cH1
T&xN3
T&b
yZ1
a
tE3>b
tE3;}
}
;e12
Compare1st{nR1
c3
const
nS1
xN3
nS1
b
yZ1
a.first<b.first;}
nR1
c3
const
nS1
a,T1
b
yZ1
a.first<b;}
nR1
c3
T1
xN3
nS1
b
yZ1
a<b.first;}
}
;
#ifndef FPoptimizerHashHH
#define FPoptimizerHashHH
#ifdef _MSC_VER
typedef
iF2
long
long
x91;
#define FPHASH_CONST(x) x##ULL
#else
#include <stdint.h>
typedef
uint_fast64_t
x91;
#define FPHASH_CONST(x) x##ULL
#endif
lF3
FUNCTIONPARSERTYPES{e12
fphash_t{x91
hash1,hash2;fphash_t():hash1(0),hash2(0){}
fphash_t
cH1
x91&xN3
x91&b):hash1(a),hash2(b){}
iG2==cH1
fphash_t&yY1==e22&&hash2==e32
iG2!=cH1
fphash_t&yY1!=e22||hash2!=e32
iG2<cH1
fphash_t&yY1!=e22?hash1<e22:hash2<e32}
;}
#endif
#ifndef FPOptimizer_CodeTreeHH
#define FPOptimizer_CodeTreeHH
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
lF3
l21{e12
Grammar;}
lF3
xH1{xW3
class
ByteCodeSynth;}
t5{xW3
class
CodeTree
nJ1
e12
xH2
nJ1
class
CodeTree{typedef
xR<xH2
xK>DataP;DataP
data;eX3
CodeTree();~CodeTree();e12
OpcodeTag{}
;e1
xB2
o,OpcodeTag);e12
FuncOpcodeTag{}
;e1
xB2
o,iF2
f,FuncOpcodeTag);e12
xO3{}
;e1
const
iV2&v,xO3);
#ifdef __GXX_EXPERIMENTAL_CXX0X__
e1
iV2&&v,xO3);
#endif
e12
VarTag{}
;e1
iF2
varno,VarTag);e12
CloneTag{}
;e1
lG3
b,CloneTag);void
GenerateFrom
cH1
xK3
FunctionParserBase
xK::Data&data,bool
keep_powi=false);void
GenerateFrom
cH1
xK3
FunctionParserBase
xK::Data&data,const
x1&xC2,bool
keep_powi=false);void
SynthesizeByteCode(std
xV3<iF2>&eI3,std
xV3
xK&immed,size_t&stacktop_max);void
SynthesizeByteCode(xH1
yG3&synth,bool
MustPopTemps=true
e83;size_t
SynthCommonSubExpressions(xH1::yM1
const;void
SetParams
cH1
x1&xP3
SetParamsMove(x1&tA1
CodeTree
GetUniqueRef();
#ifdef __GXX_EXPERIMENTAL_CXX0X__
void
SetParams(x1&&tA1
#endif
void
SetParam
iA1
which,lG3
b);void
nG1
size_t
which,xA2
b);void
AddParam
cH1
xA2
param);void
AddParamMove(xA2
param);void
AddParams
cH1
x1&xP3
AddParamsMove(x1&xP3
AddParamsMove(x1&lZ2,size_t
n02);void
DelParam
iA1
index);void
DelParams();void
Become
cH1
xA2
b);inline
size_t
GetParamCount(yZ1
l02).tF3;xI1
xA2
GetParam
iA1
n)yR
l02)[n];xI1
lG3
GetParam
iA1
n
yZ1
l02)[n];xI1
void
nH2
xB2
o)tN3
Opcode=o;xI1
xB2
GetOpcode()yQ
Opcode;xI1
nB
fphash_t
GetHash()yQ
Hash;xI1
const
x1&l02
yZ1
xY;xI1
x1&l02)yR
xY;xI1
size_t
y22
yQ
Depth;xI1
const
iV2&GetImmed()yQ
Value;xI1
iF2
GetVar()yQ
l12
xI1
iF2
GetFuncNo()yQ
l12
xI1
bool
IsDefined(yZ1
GetOpcode()!=nB
cNop;xI1
bool
IsImmed(yZ1
GetOpcode()==nB
cImmed;xI1
bool
IsVar(yZ1
GetOpcode()==nB
l03;xI1
iF2
GetRefCount()yQ
xM3
lL3
ReplaceWithImmed
cG1
i);void
Rehash(bool
constantfolding=true);void
Sort();inline
void
Mark_Incompletely_Hashed()tN3
Depth=0;xI1
bool
Is_Incompletely_Hashed()yQ
Depth==0;xI1
const
tP
GetOptimizedUsing()yQ
l32;xI1
void
SetOptimizedUsing
cH1
tP
g)tN3
l32=g;}
bool
RecreateInversionsAndNegations(bool
prefer_base2=false);void
FixIncompleteHashes();void
swap(xA2
b){data.swap(b.data);}
bool
IsIdenticalTo
cH1
xA2
b
e83;void
l61}
nJ1
e12
xH2{int
xM3;xB2
Opcode
iU2
Value;iF2
l12
eH
Params;nB
fphash_t
Hash;size_t
Depth;const
tP
l32;xH2();xH2
cH1
xH2&b);i2
xB2
o);i2
xB2
o,iF2
f);i2
const
iV2&i);
#ifdef __GXX_EXPERIMENTAL_CXX0X__
i2
iV2&&i);xH2(xH2&&b);
#endif
bool
IsIdenticalTo
cH1
xH2&b
e83;void
Sort();void
Recalculate_Hash_NoRecursion();private:void
eA1=cH1
xH2&b);}
nJ1
yX
CodeTreeImmed
cG1
i)yR
eJ(i
yA
xO3());}
#ifdef __GXX_EXPERIMENTAL_CXX0X__
cB1
CodeTreeImmed(iV2&&i)yR
eJ
c72
i)yA
xO3());}
#endif
cB1
CodeTreeOp(xB2
opcode)yR
eJ(opcode
yA
OpcodeTag());}
cB1
CodeTreeFuncOp(xB2
t83
iF2
f)yR
eJ(t83
f
yA
FuncOpcodeTag());}
cB1
CodeTreeVar
nT2
varno)yR
eJ(varno
yA
VarTag());}
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
tZ1
DumpHashes(xA)xR3
DumpTree(xA)xR3
DumpTreeWithIndent(xA,const
std
cC3&indent="\\"
);
#endif
}
#endif
#endif
#ifndef FPOptimizer_GrammarHH
#define FPOptimizer_GrammarHH
#include <iostream>
t5{xW3
class
CodeTree;}
lF3
l21{enum
ImmedConstraint_Value{ValueMask=0x07,Value_AnyNum=0x0,nY2=0x1,Value_OddInt=0x2,tS1=0x3,Value_NonInteger=0x4,eJ1=0x5
xF3
ImmedConstraint_Sign{SignMask=0x18,Sign_AnySign=0x00,nH1=0x08,eK1=0x10,Sign_NoIdea=0x18
xF3
ImmedConstraint_Oneness{OnenessMask=0x60,Oneness_Any=0x00,Oneness_One=0x20,Oneness_NotOne=0x40
xF3
ImmedConstraint_Constness{ConstnessMask=0x180,Constness_Any=0x00,tR1=0x80,Constness_NotConst=0x100
xF3
Modulo_Mode{Modulo_None=0,Modulo_Radians=1
xF3
Situation_Flags{LogicalContextOnly=0x01,NotForIntegers=0x02,OnlyForIntegers=0x04,OnlyForComplex=0x08,NotForComplex=0x10
xF3
nR2{NumConstant,tG3,SubFunction
xF3
ParamMatchingType{PositionalParams,SelectedParams,AnyParams,GroupFunction
xF3
RuleType{ProduceNewTree,ReplaceParams}
;
#ifdef __GNUC__
# define PACKED_GRAMMAR_ATTRIBUTE __attribute__((packed))
#else
# define PACKED_GRAMMAR_ATTRIBUTE
#endif
typedef
std::pair<nR2,const
void*>e42
nJ1
e42
e81
nT2
paramlist,nA1)nJ1
bool
ParamSpec_Compare
cH1
void*xN3
void*b,nR2
type);iF2
ParamSpec_GetDepCode
cH1
e42&b);e12
ParamSpec_ParamHolder{nA1:8;iF2
constraints:9;iF2
depcode:15;}
e9
xW3
e12
ParamSpec_NumConstant{iV2
constvalue;iF2
modulo;}
;e12
iZ{iF2
param_count:2;iF2
param_list:30;xB2
subfunc_opcode:8;ParamMatchingType
match_type:3;iF2
n1:5;}
e9
e12
ParamSpec_SubFunction{iZ
data;iF2
constraints:9;iF2
depcode:7;}
e9
e12
Rule{RuleType
ruletype:2;iF2
situation_flags:5;iF2
repl_param_count:2+9;iF2
repl_param_list:30;iZ
match_tree;}
e9
e12
Grammar{iF2
rule_count;iF2
short
rule_list[999
cS
extern
const
Rule
grammar_rules[];}
tZ1
DumpParam
cH1
e42&p,std::ostream&o=std::cout)xR3
DumpParams
nT2
paramlist,iF2
count,std::ostream&o=std::cout);}
#endif
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#define CONSTANT_POS_INF HUGE_VAL
#define CONSTANT_NEG_INF (-HUGE_VAL)
lF3
FUNCTIONPARSERTYPES{xW3
inline
iV2
fp_const_pihalf()yR
fp_const_pi
xK()*iV2(0.5);}
xW3
inline
iV2
fp_const_twopi(tO2
fp_const_pi
xK());lH3
fp_const_twoe(tO2
fp_const_e
xK());lH3
fp_const_twoeinv(tO2
fp_const_einv
xK());lH3
fp_const_negativezero()yR-fp_epsilon
xK();}
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
#include <iostream>
yL1{using
lF3
l21;using
t5;using
lF3
FUNCTIONPARSERTYPES
nJ1
class
MatchInfo{eX3
std
xV3<std::pair<bool,eH> >lQ;eH
paramholder_matches;std
xV3<iF2>tR;eX3
MatchInfo():lQ(),paramholder_matches(),tR(){}
eX3
bool
SaveOrTestRestHolder
nT2
n1,iH2&iF1){cV1{lQ.t73
n1+1);lQ
tQ=iF1
xJ2
if(lQ[n1
tC3==false){lQ
tQ=iF1
xJ2
iH2&found=lQ[n1
tD3;if(iF1.tF3!=found.tF3
cL
for
iA1
a=0;a<iF1.tF3;++a)if(!iF1[a]xF
found[a])cL
return
true
lL3
SaveRestHolder
nT2
n1,eH&iF1){cV1
lQ.t73
n1+1);lQ
tQ.swap(iF1);}
bool
SaveOrTestParamHolder
nT2
nZ,lC2
xQ3){if(cE
tF3<=nZ){cE
xS3
nZ+1);cE
t73
nZ);cE
push_back(xQ3)xJ2
if(!cH.t91
cH=xQ3
xJ2
return
xQ3
xF
cH
lI3
SaveMatchedParamIndex(nA1){tR.push_back(index);}
lC2
GetParamHolderValueIfFound
nT2
nZ
e83{static
const
eJ
dummytree;if(cE
tF3<=nZ)return
dummytree
n31
cH;}
lC2
GetParamHolderValue
nT2
nZ
yZ1
cH;}
bool
HasRestHolder
nT2
n1
yZ1
lQ.tF3>n1&&lQ[n1
tC3==true;}
iH2&GetRestHolderValues
nT2
n1
e83{static
iH2
empty_result;cV1
return
empty_result
n31
lQ[n1
tD3;}
const
std
xV3<iF2>&GetMatchedParamIndexes(yZ1
tR
lL3
swap(tB
b){lQ.swap(b.lQ);cE
swap(b.paramholder_matches);tR.swap(b.tR);}
tB
eA1=cH1
tB
b){lQ=b.lQ;paramholder_matches=b.paramholder_matches;tR=b.tR
n31*this;}
}
;class
e2;typedef
xR<e2>cW;class
e2{eX3
int
xM3;eX3
e2():xM3(0){}
virtual~e2(){}
}
;e12
lX1{bool
found;cW
specs;lX1(bool
f):found(f),specs(){}
lX1(bool
f,const
cW&s):found(f),specs(s){}
}
xR3
SynthesizeRule
cH1
eT2
eJ&tree
tJ3)nJ1
lX1
TestParam
cH1
e42&yG2
lC2
tree,const
cW&start_at
tJ3)nJ1
lX1
TestParams(tC&nN,lC2
tree,const
cW&start_at
tJ3,bool
eY2
nJ1
bool
ApplyGrammar
cH1
Grammar&tX2,FPoptimizer_CodeTree::eJ&tree,bool
from_logical_context=false)xR3
ApplyGrammars(FPoptimizer_CodeTree::e62
nJ1
bool
IsLogisticallyPlausibleParamsMatch(tC&c32,lC2
tree);}
lF3
l21{tZ1
DumpMatch
cH1
eT2
nO,const
FPoptimizer_Optimize::tB
info,bool
DidMatch,std::ostream&o=std::cout)xR3
DumpMatch
cH1
eT2
nO,const
FPoptimizer_Optimize::tB
info,nH3
i43,std::ostream&o=std::cout);}
#endif
#include <string>
nI3
l21::nR2
yS1=false);nI3
xB2
yS1=false);
#include <string>
#include <sstream>
#include <assert.h>
#include <iostream>
using
lF3
l21;using
lF3
FUNCTIONPARSERTYPES;nI3
l21::nR2
yS1){
#if 1
nH3
p=0;switch(opcode
y33
NumConstant:p="NumConstant"
;lC
tG3:p="ParamHolder"
;lC
SubFunction:p="SubFunction"
;yF3
std::ostringstream
tmp;assert(p);tmp<<p;if(pad)while(tmp.str().tF3<12)tmp<<' 'n31
tmp.str();
#else
std::ostringstream
tmp;tmp<<opcode;if(pad)while(tmp.str().tF3<5)tmp<<' 'n31
tmp.str();
#endif
}
nI3
xB2
yS1){
#if 1
nH3
p=0;switch(opcode
y33
cAbs:p="cAbs"
;lC
cAcos:p="cAcos"
;lC
cAcosh:p="cAcosh"
;lC
cArg:p="cArg"
;lC
cAsin:p="cAsin"
;lC
cAsinh:p="cAsinh"
;lC
cAtan:p="cAtan"
;lC
cAtan2:p="cAtan2"
;lC
cAtanh:p="cAtanh"
;lC
cCbrt:p="cCbrt"
;lC
cCeil:p="cCeil"
;lC
cConj:p="cConj"
;lC
cCos:p="cCos"
;lC
cCosh:p="cCosh"
;lC
cCot:p="cCot"
;lC
cCsc:p="cCsc"
;lC
cExp:p="cExp"
;lC
cExp2:p="cExp2"
;lC
cFloor:p="cFloor"
;lC
cHypot:p="cHypot"
;lC
cIf:p="cIf"
;lC
cImag:p="cImag"
;lC
cInt:p="cInt"
;lC
cLog:p="cLog"
;lC
cLog2:p="cLog2"
;lC
cLog10:p="cLog10"
;lC
cMax:p="cMax"
;lC
cMin:p="cMin"
;lC
cPolar:p="cPolar"
;lC
cPow:p="cPow"
;lC
cReal:p="cReal"
;lC
cSec:p="cSec"
;lC
cSin:p="cSin"
;lC
cSinh:p="cSinh"
;lC
cSqrt:p="cSqrt"
;lC
cTan:p="cTan"
;lC
cTanh:p="cTanh"
;lC
cTrunc:p="cTrunc"
;lC
cImmed:p="cImmed"
;lC
cJump:p="cJump"
;lC
cNeg:p="cNeg"
;lC
cAdd:p="cAdd"
;lC
cSub:p="cSub"
;lC
cMul:p="cMul"
;lC
cDiv:p="cDiv"
;lC
cMod:p="cMod"
;lC
cEqual:p="cEqual"
;lC
cNEqual:p="cNEqual"
;lC
cLess:p="cLess"
;lC
cLessOrEq:p="cLessOrEq"
;lC
cGreater:p="cGreater"
;lC
cGreaterOrEq:p="cGreaterOrEq"
;lC
cNot:p="cNot"
;lC
cAnd:p="cAnd"
;lC
cOr:p="cOr"
;lC
cDeg:p="cDeg"
;lC
cRad:p="cRad"
;lC
cFCall:p="cFCall"
;lC
cPCall:p="cPCall"
;break;
#ifdef FP_SUPPORT_OPTIMIZER
case
cFetch:p="cFetch"
;lC
cPopNMov:p="cPopNMov"
;lC
eF3:p="cLog2by"
;lC
cNop:p="cNop"
;break;
#endif
case
cSinCos:p="cSinCos"
;lC
cSinhCosh:p="cSinhCosh"
;lC
cZ3:p="cAbsNot"
;lC
cAbsNotNot:p="cAbsNotNot"
;lC
cAbsAnd:p="cAbsAnd"
;lC
cAbsOr:p="cAbsOr"
;lC
cAbsIf:p="cAbsIf"
;lC
cDup:p="cDup"
;lC
cInv:p="cInv"
;lC
cSqr:p="cSqr"
;lC
cRDiv:p="cRDiv"
;lC
cRSub:p="cRSub"
;lC
cNotNot:p="cNotNot"
;lC
cRSqrt:p="cRSqrt"
;lC
l03:p="VarBegin"
;yF3
std::ostringstream
tmp;assert(p);tmp<<p;if(pad)while(tmp.str().tF3<12)tmp<<' 'n31
tmp.str();
#else
std::ostringstream
tmp;tmp<<opcode;if(pad)while(tmp.str().tF3<5)tmp<<' 'n31
tmp.str();
#endif
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
#ifndef FP_GENERATING_POWI_TABLE
enum{MAX_POWI_BYTECODE_LENGTH=20}
;
#else
enum{MAX_POWI_BYTECODE_LENGTH=999}
;
#endif
enum{MAX_MULI_BYTECODE_LENGTH=3}
;lF3
xH1{xW3
class
ByteCodeSynth{eX3
ByteCodeSynth():ByteCode(),Immed(),cM(),xS(0),StackMax(0){i8
xS3
64);Immed.xS3
8);cM.xS3
16
lI3
Pull(std
xV3<iF2>&bc,std
xV3
xK&imm,size_t&StackTop_max){for(eQ1=0;a<i8
tF3;++a){l52]&=~nK2;}
i8
swap(bc);Immed.swap(imm);StackTop_max=StackMax;}
size_t
GetByteCodeSize(yZ1
i8
tF3;}
size_t
GetStackTop(yZ1
xS
lL3
PushVar
nT2
varno){xK2
varno);eZ2}
void
PushImmed(iV2
immed
nL
xK2
cImmed);Immed.push_back(immed);eZ2}
void
StackTopIs(nO,int
offset=0){if((int)xS>offset){cM
l62
first=true;cM
l62
second=tree;}
}
bool
IsStackTop(nO,int
offset=0
yZ1(int)xS>offset&&cM
l62
first&&cM
l62
second
xF
tree);xI1
void
EatNParams
nT2
eat_count){xS-=eat_count
lL3
ProducedNParams
nT2
produce_count){xJ1
xS+produce_count
lI3
DoPopNMov
iA1
e52,size_t
srcpos
nL
xK2
cPopNMov)nS2
e52)nS2
srcpos);xJ1
srcpos+1);cM[e52]=cM[srcpos];xJ1
e52+1
lI3
DoDup
iA1
xT3
nL
if(xT3==xS-1){xK2
cDup
t43{xK2
cFetch)nS2
xT3);}
eZ2
cM[xS-1]=cM[xT3];}
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
i31
int>void
Dump(){std::ostream&o=std::cout;o<<"Stack state now("
<<xS<<"):\n"
yM3
0;a<xS;++a){o<<a<<": "
;if(cM[a
tC3){nO=cM[a
tD3;o<<'['<<std::hex<<(void*)(&tree.l02))<<std::dec<<','<<tree.GetRefCount()<<']'i13(tree,o
t43
o<<"?"
;o<<"\n"
;}
o<<std::flush;}
#endif
size_t
xU3
nO
e83{for
iA1
a=xS;a-->0;)if(cM[a
tC3&&cM[a
tD3
xF
tree
nE3
a
n31
t02;}
bool
Find(nO
yZ1
xU3
tree)!=t02;}
bool
FindAndDup(nO){size_t
pos=xU3
tree
eY3
pos!=t02){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<l24"duplicate at ["
<<pos<<"]: "
i13(tree
lT2" -- issuing cDup or cFetch\n"
;
#endif
DoDup(pos)xJ2
return
tR3
e12
IfData{size_t
ofs;}
;void
SynthIfStep1
l72,xB2
op
x21
xL2=i8
tF3;xK2
op
i9
xK2
nK2
lI3
SynthIfStep2
l72
x21
iS
nU2)+2);iS
2
y7
l42
xL2=i8
tF3;xK2
cJump
i9
xK2
nK2
lI3
SynthIfStep3
l72
x21
i8
back()|=nK2;iS
nU2)-1);iS
2
y7
l42
eZ2
for
iA1
a=0;a<xL2;++a){if(l52]==cJump&&l52+1]==(nK2|(xL2-1))){l52+nU2)-1);l52+2
y7
l42
l73
l52]y33
cAbsIf:case
cIf:case
cJump:case
cPopNMov:a+=2;lC
cFCall:case
cPCall:case
cFetch:a+=1;break;c73
yF3}
}
protected:void
xJ1
size_t
value){xS=value;if(xS>lN3{StackMax=xS;cM.t73
lN3;}
}
protected:std
xV3<iF2>ByteCode;std
xV3
xK
Immed;std
xV3<std::pair<bool,FPoptimizer_CodeTree::eJ> >cM;size_t
xS;size_t
StackMax;private:void
incStackPtr(){if(xS+2>lN3
cM.t73
StackMax=xS+2);}
i31
bool
IsIntType,bool
IsComplexType>e12
c42{}
;eX3
void
AddOperation
nT2
t83
iF2
eat_count,iF2
produce_count=1){EatNParams(eat_count);lO1
opcode);ProducedNParams(produce_count
lI3
lO1
iF2
t83
c42<false,false>yS
false,true>yS
true,false>yS
true,true>);inline
void
lO1
iF2
opcode){lO1
t83
c42<bool(nB
IsIntType
xK::nT3),bool(nB
IsComplexType
xK::nT3)>());}
}
nJ1
e12
SequenceOpCode
nJ1
e12
tN1{static
cI
AddSequence;static
cI
MulSequence;}
xR3
x01
long
count,cI&eU,yM1;}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
lF3
FUNCTIONPARSERTYPES;lF3
xH1{xW3
e12
SequenceOpCode{iV2
basevalue;iF2
op_flip;iF2
op_normal,op_normal_flip;iF2
op_inverse,op_inverse_flip;}
nJ1
cI
tN1
xK::AddSequence={y21,cNeg
xE
cAdd,cSub,cRSub}
nJ1
cI
tN1
xK::MulSequence={iV2(1),cInv,cMul,cMul,cDiv,cRDiv}
;
#define findName(a,b,c) "var"
#define TryCompilePowi(o) false
#define mData this
#define mByteCode ByteCode
#define mImmed Immed
xY3
false,false
l31
# define FP_FLOAT_VERSION 1
# define FP_COMPLEX_VERSION 0
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
xY3
true,false
l31
# define FP_FLOAT_VERSION 0
# define FP_COMPLEX_VERSION 0
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
#ifdef FP_SUPPORT_COMPLEX_NUMBERS
xY3
false,true
l31
# define FP_FLOAT_VERSION 1
# define FP_COMPLEX_VERSION 1
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
xY3
true,true
l31
# define FP_FLOAT_VERSION 0
# define FP_COMPLEX_VERSION 1
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
#endif
#undef findName
#undef mImmed
#undef mByteCode
#undef mData
#undef TryCompilePowi
}
using
lF3
xH1;
#define POWI_TABLE_SIZE 256
#define POWI_WINDOW_SIZE 3
lF3
xH1{
#ifndef FP_GENERATING_POWI_TABLE
extern
const
iF2
char
powi_table[POWI_TABLE_SIZE];const
#endif
iF2
char
powi_table[POWI_TABLE_SIZE]={0,1,1,1,2,1,2,1,y63
4,1,2,y73
2,1,y63
8,eK3
y83
15,1,16,1,2,1,4,1,2,y73
2,1,4,eK3
1,16,1,25,y83
27,5,8,3,2,1,30,1,31,3,32,1,2,1,y63
8,1,2,y83
39,1,16,137,2,1,4,eK3
y73
45,135,4,31,2,5,32,1,2,131,50,1,51,1,8,3,2,1,54,1,55,3,16,1,57,133,4,137,2,135,60,1,61,3,62,133,63,1,iG1
131,iG1
139,l82
e3
30,1,130,137,2,31,l82
e3
e3
130,eK3
1,e3
e3
2,1,130,133,iG1
61,130,133,62,139,130,137,e3
l82
e3
e3
iG1
131,e3
e3
130,131,2,133,l82
130,141,e3
130,eK3
1,e3
5,135,e3
l82
e3
l82
130,133,130,141,130,131,e3
e3
2,131}
;}
static
xC3
yY=256;
#define FPO(x)
lF3{class
PowiCache{private:int
iB[yY];int
iH1[yY];eX3
PowiCache():iB(),iH1(){iB[1]=1;}
bool
Plan_Add(c82,int
count){cC1>=yY
cL
iH1[t12+=count
n31
iB[t12!=0
lL3
lO3
c82){cC1<yY)iB[t12=1
lL3
Start
iA1
value1_pos){for(int
n=2;n<yY;++n)iB[n]=-1;Remember(1,value1_pos);DumpContents();}
int
Find(c82
e83{cC1<yY){if(iB[t12>=0){FPO(iM3(iR3,"* I found %ld from cache (%u,%d)\n",value,(unsigned)cache[value],iN3 value]))n31
iB[t12;eW2-1
lL3
Remember(c82,size_t
l04){cC1>=yY)return;FPO(iM3(iR3,"* Remembering that %ld can be found at %u (%d uses remain)\n",value,(unsigned)l04,iN3 value]));iB[t12=(int)l04
lL3
DumpContents
e73
FPO(for(int a=1;a<POWI_CACHE_SIZE;++a)if(cache[a]>=0||iN3 a]>0){iM3(iR3,"== cache: sp=%d, val=%d, needs=%d\n",cache[a],a,iN3 a]);})}
int
UseGetNeeded(c82){cC1>=0&&value<yY)return--iH1[t12
n31
0;}
}
nJ1
size_t
y8
long
count
lF2
cI&eU,yM1
xR3
c01
size_t
apos,long
aval,size_t
bpos,long
bval
lF2
iF2
cumulation_opcode,iF2
cimulation_opcode_flip,yM1;void
l41
c82
lF2
int
need_count,int
l51=0){cC1<1)return;
#ifdef FP_GENERATING_POWI_TABLE
if(l51>32)throw
false;
#endif
if(iB.Plan_Add(value,need_count
nE3;long
yA3
1;cC1<POWI_TABLE_SIZE){yA3
powi_table[t12
lM3&128){half&=127
lM3&64)yA3-tB2
FPO(iM3(iR3,"value=%ld, half=%ld, otherhalf=%ld\n",value,half,value/half));l41
half
y93
iB.lO3
half)n31;}
lJ2
half&64){yA3-tB2}
}
else
cC1&1)yA3
value&((1<<POWI_WINDOW_SIZE)-1);else
yA3
value/2;long
cN=value-half
lM3>cN||half<0)std::swap(half,cN);FPO(iM3(iR3,"value=%ld, half=%ld, otherhalf=%ld\n",value,half,otherhalf))lM3==cN){l41
half,iB,2,l51+1
t43{l41
half
y93
l41
cN>0?cN:-cN
y93}
iB.lO3
value);}
xW3
size_t
y8
c82
lF2
cI&eU,yM1{int
yB3=iB.Find(value
eY3
yB3>=0)yR
yB3;}
long
yA3
1;cC1<POWI_TABLE_SIZE){yA3
powi_table[t12
lM3&128){half&=127
lM3&64)yA3-tB2
FPO(iM3(iR3,"* I want %ld, my plan is %ld * %ld\n",value,half,value/half));size_t
xM2=y8
half
tX1
if(iB
lL2
half)>0||xM2!=tY1){iI1
xM2)cS2
half,tY1);}
x01
value/half,eU,iS2
size_t
l04=tY1
cS2
value,l04);iB.DumpContents()n31
l04;}
lJ2
half&64){yA3-tB2}
}
else
cC1&1)yA3
value&((1<<POWI_WINDOW_SIZE)-1);else
yA3
value/2;long
cN=value-half
lM3>cN||half<0)std::swap(half,cN);FPO(iM3(iR3,"* I want %ld, my plan is %ld + %ld\n",value,half,value-half))lM3==cN){size_t
xM2=y8
half
tX1
c01
xM2,half,xM2,half,iB,eU.op_normal,eU.op_normal_flip,iS2}
else{long
part1=half;long
part2=cN>0?cN:-cN;size_t
part1_pos=y8
part1
tX1
size_t
part2_pos=y8
part2
tX1
FPO(iM3(iR3,"Subdivide(%ld: %ld, %ld)\n",value,half,otherhalf));c01
part1_pos,part1,part2_pos,part2,iB,cN>0?eU.op_normal:eU.op_inverse,cN>0?eU.op_normal_flip:eU.op_inverse_flip,iS2}
size_t
l04=tY1
cS2
value,l04);iB.DumpContents()n31
l04;}
tZ1
c01
size_t
apos,long
aval,size_t
bpos,long
bval
lF2
iF2
cumulation_opcode,iF2
cumulation_opcode_flip,yM1{int
a_needed=iB
lL2
aval);int
yC3=iB
lL2
bval);bool
lM2
tO3
#define DUP_BOTH() do{if(apos<bpos){size_t tmp=apos;apos=bpos;bpos=tmp;lM2=!lM2;}FPO(iM3(iR3,"-> "iX3 iX3"op\n",(unsigned)apos,(unsigned)bpos));iI1 apos);iI1 apos==bpos?tY1:bpos);}while(0)
#define DUP_ONE(p) do{FPO(iM3(iR3,"-> "iX3"op\n",(unsigned)p));iI1 p);}while(0)
if(a_needed>0){if(yC3>0){nL2}
e63
bpos!=tY1)nL2
else{lN2
lM2=!lM2;}
}
}
lJ2
yC3>0){if(apos!=tY1)nL2
else
DUP_ONE(bpos);}
e63
apos==bpos&&apos==tY1)lN2
iR2
tY1&&bpos==yD3
2){FPO(iM3(iR3,"-> op\n"));lM2=!lM2;}
iR2
yD3
2&&bpos==tY1)FPO(iM3(iR3,"-> op\n"));iR2
tY1)DUP_ONE(bpos);lJ2
bpos==tY1){lN2
lM2=!lM2;}
else
nL2}
yE3
lM2?cumulation_opcode_flip:cumulation_opcode,2);}
tZ1
cD1
long
count,cI&eU,yM1{while
eL3<256){int
yA3
xH1::powi_table[count]lM3&128){half&=127;cD1
half,eU,iS2
count/=half;}
else
yF3
if
eL3==1)return;if(!eL3&1)){yE3
cSqr,1);cD1
count/2,eU,iS2}
else{iI1
tY1);cD1
count-1,eU,iS2
yE3
cMul,2);}
}
}
lF3
xH1{tZ1
x01
long
count,cI&eU,yM1{if
eL3==0)lH2
eU.basevalue);else{bool
t22
tO3
if
eL3<0){t22=true;count=-count;}
if(false)cD1
count,eU,iS2
lJ2
count>1){PowiCache
iB;l41
count,iB,1);size_t
xK1=synth.GetStackTop();iB.Start(tY1);FPO(iM3(iR3,"Calculating result for %ld...\n",count));size_t
xN2=y8
count
tX1
size_t
n_excess=yD3
xK1;if(n_excess>0||xN2!=xK1-1){synth.DoPopNMov(xK1-1,xN2);}
}
if(t22)yE3
eU.op_flip,1);}
}
}
#endif
#ifndef FPOptimizer_ValueRangeHH
#define FPOptimizer_ValueRangeHH
t5{lF3
lP3{iM
e12
Comp{}
;i31>e12
Comp<nB
cLess>xZ3<n9
cLessOrEq>xZ3<=n9
cGreater>xZ3>n9
cGreaterOrEq>xZ3>=n9
cEqual>xZ3==n9
cNEqual>xZ3!=b;}
}
;}
xW3
e12
cJ1{iV2
val;bool
known;cJ1():val(),known(false){}
cJ1
cG1
v):val(v),known(true){xI1
void
set
cG1
v){known=true;val=v
lL3
set(iV2(iJ1(iV2),cJ1
y9)val=iI2
void
set(iV2(iJ1
cG1),cJ1
y9)val=iI2
iM
void
set_if(iV2
v
tG1
iJ1(iV2),cJ1
y9&&lP3::Comp<Compare>()(val,v))val=iI2
iM
void
set_if
cG1
v
tG1
iJ1
cG1),cJ1
y9&&lP3::Comp<Compare>()(val,v))val=iI2}
nJ1
e12
range{cJ1
xK
min,max;range():min(),max(){}
range(iV2
mi,iV2
ma):min(mi),max(ma){}
range(bool,iV2
ma):min(),max(ma){}
range(iV2
mi,bool):min(mi),max(){}
void
set_abs();void
set_neg();}
nJ1
bool
IsLogicalTrueValue
cH1
c23&p
tR2
nJ1
bool
IsLogicalFalseValue
cH1
c23&p
tR2;}
#endif
#ifndef FPOptimizer_RangeEstimationHH
#define FPOptimizer_RangeEstimationHH
t5{enum
TriTruthValue{IsAlways,IsNever,Unknown}
nJ1
c23
iN
lC2
tree)nJ1
bool
IsLogicalValue
cH1
e62
nJ1
TriTruthValue
GetIntegerInfo
cH1
e62
nJ1
xL1
GetEvennessInfo
cH1
e62{if(!tree
i61)return
Unknown;yJ1=eC3;if(isEvenInteger(value
nF3
isOddInteger(value
nG3
xW3
xL1
GetPositivityInfo
cH1
e62{c23
p=iN
tree
eY3
p
n51&&p
tF2>=iV2(nF3
p
yN
known
lH1
nG3
xW3
xL1
GetLogicalValue
lQ3
tree
tR2{c23
p=iN
tree
eY3
IsLogicalTrueValue(p,abs
nF3
IsLogicalFalseValue(p,abs
nG3}
#endif
#ifndef FPOptimizer_ConstantFoldingHH
#define FPOptimizer_ConstantFoldingHH
t5{tZ1
ConstantFolding(e62;}
#endif
lF3{using
lF3
FUNCTIONPARSERTYPES;using
t5;e12
ComparisonSetBase{enum{t93=0x1,Eq_Mask=0x2,Le_Mask=0x3,tA3=0x4,tB3=0x5,Ge_Mask=0x6}
;static
int
Swap_Mask(int
m)yR(m&Eq_Mask)|((m&t93)?tA3:0)|((m&tA3)?t93:0);}
enum
c11{Ok,BecomeZero,BecomeOne,xP2
xF3
nW2{cond_or,iK2,iL2,iM2}
;}
nJ1
e12
ComparisonSet:public
ComparisonSetBase{e12
t32{eJ
a;eJ
b;int
relationship;t32():a(),b(),relationship(){}
}
;std
xV3<t32>eK;e12
Item{eJ
value;bool
c52;Item():value(),c52(false){}
}
;std
xV3<Item>cK1;int
xM1;ComparisonSet():eK(),cK1(),xM1(0){}
c11
AddItem
lQ3
a,bool
c52,nW2
type){for
iA1
c=0;c<cK1.tF3;++c)if(cK1[c].value
xF
a)){if(c52!=cK1[c].c52){iQ
cW1
case
iM2:cK1.erase(cK1.begin()+c);xM1+=1
n31
xP2;case
iK2:case
iL2:cX1
eW2
xP2;}
Item
pole;pole.value=a;pole.c52=c52;cK1.push_back(pole)n31
Ok;}
c11
AddRelationship(eJ
a,eJ
b,int
tO1,nW2
type){iQ
if(tO1==7)cW1
lC
iM2:if(tO1==7){xM1+=1
n21
lC
iK2:case
iL2:if(tO1==0)cX1
yF3
if(!(a.GetHash()<b.GetHash())){a.swap(b);tO1=Swap_Mask(tO1);}
for
iA1
c=0;c<eK.tF3;++c){if(eK[c].a
xF
a)&&eK[c].b
xF
b)){iQ{int
yQ3=xO2|tO1;if(yQ3==7)cW1
xO2=yQ3;yF3
case
iK2:case
iL2:{int
yQ3=xO2&tO1;if(yQ3==0)cX1
xO2=yQ3;yF3
case
iM2:{int
newrel_or=xO2|tO1;int
xR2=xO2&tO1;lP2
5&&xR2==0){xO2=tB3
n21
lP2
7&&xR2==0){xM1+=1;eK.erase(eK.begin()+c)n21
lP2
7&&xR2==Eq_Mask){xO2=Eq_Mask;xM1+=1
n21
eR2}
return
xP2;}
}
t32
comp;comp.a=a;comp.b=b;comp.relationship=tO1;eK.push_back(comp)n31
Ok;}
}
;nR1
iV2,xK3
CondType>bool
ConstantFolding_LogicCommon(eJ&tree,CondType
xA1,bool
xS2){bool
should_regenerate
tO3
ComparisonSet
xK
comp;nB1{xK3
c6
c11
eM3=c6
Ok;lC2
atree=t33;switch(atree
nC
y33
cEqual
lG
Eq_Mask,xA1);lC
cNEqual
lG
tB3,xA1);lC
cLess
lG
t93,xA1);lC
cLessOrEq
lG
Le_Mask,xA1);lC
cGreater
lG
tA3,xA1);lC
cGreaterOrEq
lG
Ge_Mask,xA1);lC
cNot:eM3
yO1
l8
0),true,xA1);lC
cNotNot:eM3
yO1
l8
0),false,xA1);break;c73
if(xS2||IsLogicalValue(atree))eM3
yO1,false,xA1);l73
eM3){ReplaceTreeWithZero:xO
0)n31
true;ReplaceTreeWithOne:xO
1);nV
c6
Ok:lC
c6
BecomeZero
eL
c6
BecomeOne:i6
c6
xP2:c41
yF3}
if(should_regenerate){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_LogicCommon: "
yA1
#endif
if(xS2){tree.DelParams(t43{for
iY{lC2
atree=lE2
a
eY3
IsLogicalValue(atree))eK2}
for
iA1
a=0;a<comp.cK1.tF3;++a){if(comp.cK1[a].c52
iJ2
cNot);r.c5
r
cK
r);}
lJ2!xS2
iJ2
cNotNot);r.c5
r
cK
r
t43
tree.c5}
for
iA1
a=0;a<comp.eK.tF3;++a
iJ2
cNop);switch(comp.eK[a
eA
y33
c6
t93:r
xD
cLess);lC
c6
Eq_Mask:r
xD
cEqual);lC
c6
tA3:r
xD
cGreater);lC
c6
Le_Mask:r
xD
cLessOrEq);lC
c6
tB3:r
xD
cNEqual);lC
c6
Ge_Mask:r
xD
cGreaterOrEq
cT2
r
yP1
comp.eK[a].a);r
yP1
comp.eK[a].b);r
cK
r);}
if(comp.xM1!=0)tree.yC
iV2(comp.xM1)));
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_LogicCommon: "
yA1
#endif
return
true;}
return
tR3
xN1
ConstantFolding_AndLogic(iS3(tree.GetOpcode()==cAnd||tree.GetOpcode()==cAbsAnd)n31
nI
iK2,true);}
xN1
ConstantFolding_OrLogic(iS3(tree.GetOpcode()==cOr||tree.GetOpcode()==cAbsOr)n31
nI
cond_or,true);}
xN1
ConstantFolding_AddLogicItems(iS3(tree.GetOpcode()==cAdd)n31
nI
iM2,false);}
xN1
ConstantFolding_MulLogicItems(iS3(tree.GetOpcode()==cMul)n31
nI
iL2,false);}
}
#include <vector>
#include <map>
#include <algorithm>
lF3{using
lF3
FUNCTIONPARSERTYPES;using
t5;e12
CollectionSetBase{enum
xR1{Ok,xP2}
;}
nJ1
e12
CollectionSet:public
CollectionSetBase{e12
c21{eJ
value;eJ
xT2;bool
cX;c21():value(),xT2(),cX(false){}
c21
lQ3
v,lC2
f):value(v),xT2(f),cX(false){}
}
;std::multimap<fphash_t,c21>iC;typedef
xK3
std::multimap<fphash_t,c21>::yV3
xS1;CollectionSet():iC(){}
xS1
FindIdenticalValueTo
lQ3
value){fphash_t
hash=value.GetHash();for(xS1
i=iC.xU2
hash);i!=iC.cY1
hash;++i){cC1
xF
i
e82.value
nE3
i;}
return
iC.end();}
bool
Found
cH1
xS1&b)yR
b!=iC.end();}
xR1
AddCollectionTo
lQ3
xT2,const
xS1&into_which){c21&c=into_which
e82;if(c.cX)c.xT2
c9
xT2);else{eJ
add;add
xD
cAdd);add
yP1
c.xT2);add
c9
xT2);c.xT2.swap(add);c.cX=true;}
return
xP2;}
xR1
nX2
lQ3
value,lC2
xT2){const
fphash_t
hash=value.GetHash();xS1
i=iC.xU2
hash);for(;i!=iC.cY1
hash;++i){if(i
e82.value
xF
value
nE3
AddCollectionTo(xT2,i);}
iC.yR3,std::make_pair(hash,c21(value,xT2)))n31
Ok;}
xR1
nX2
lQ3
a)yR
nX2(a,nF1
1)));}
}
nJ1
e12
ConstantExponentCollection{typedef
eH
yT3;typedef
std::x31
xV2;std
xV3<xV2>data;ConstantExponentCollection():data(){}
void
MoveToSet_Unique
cG1
eB1&eC1){data.push_back(std::x31(eB1()));data.back()tE3.swap(eC1
lI3
MoveToSet_NonUnique
cG1
eB1&eC1){xK3
std
xV3<xV2>::yV3
i=std::xU2
data.iN2
data.end(),exponent,Compare1st()eY3
i!=data.cY1
xG2{i
e82.yR3
e82.end(),eC1.iN2
eC1.end()t43{data.yR3,std::x31
t62,eC1));}
}
bool
Optimize(){bool
changed
tO3
std::sort(data.iN2
data.end(),Compare1st());redo:for
iA1
a=0;a<data.tF3;++a
tJ2
exp_a=data[a
tC3
cA3
exp_a
tG1
1))cP2;for
iA1
b=a+1;b<data.tF3;++b
tJ2
exp_b=data[b
tC3
iU2
xW2=exp_b-exp_a;if(xW2>=fp_abs(exp_a))break
iU2
exp_diff_still_probable_integer=xW2*iV2(16
eY3
t42
exp_diff_still_probable_integer)&&!(t42
exp_b)&&!t42
xW2))){yT3&a_set=lQ2;yT3&b_set=data[b
tD3;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantExponentCollection iteration:\n"
;t52
cout);
#endif
if(isEvenInteger(exp_b)&&!isEvenInteger(xW2+exp_a)){eJ
tmp2;tmp2
tH3
tmp2
iO3
b_set);tmp2
eT1
t81
cAbs);tmp
yP1
tmp2);tmp
eT1
b_set.t73
1);b_set[0].i72}
a_set.insert(a_set.end(),b_set.iN2
b_set.end());yT3
b_copy=b_set;data.erase(data.begin()+b);MoveToSet_NonUnique(xW2,b_copy);x62
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantExponentCollection iteration:\n"
;t52
cout);
#endif
cZ}
eW2
changed;}
#ifdef DEBUG_SUBSTITUTIONS
void
t52
ostream&out){for
iA1
a=0;a<data.tF3;++a){out.precision(12);out<<data[a
tC3<<": "
;i81
lQ2.tF3;++b){if(b>0)out<<'*'i13(lQ2[b],out);}
out<<std::endl;}
}
#endif
}
nJ1
static
eJ
x41
eJ&value,bool&xH
y13
value
nC
y33
cPow:{eJ
e92
value
l8
1);value.y81
n31
exponent;}
case
cRSqrt:value.y81;xH=true
n31
nF1-0.5));case
cInv:value.y81;xH=true
n31
nF1-1));c73
yF3
return
nF1
1));}
cA1
void
eD1
eV1&mul,lC2
tree,lC2
xT2,bool&c31
bool&xH){nB1{eJ
value
y91
a));eJ
exponent(x41
value,xH)eY3!xT2
i61||xT2
nC1!=iV2(1.0)){eJ
cZ1;cZ1
tH3
cZ1
c9
xG2;cZ1
c9
xT2);cZ1
eT1
exponent.swap(cZ1);}
#if 0 /* FIXME: This does not work */
cC1
nC==cMul){if(1){bool
exponent_is_even=exponent
i61&&isEvenInteger
t62
nC1);i81
value.tP3{bool
tmp
tO3
eJ
val(value
l8
b));eJ
exp(x41
val,tmp)eY3
exponent_is_even||(exp
i61&&isEvenInteger(exp
nC1))){eJ
cZ1;cZ1
tH3
cZ1
c9
xG2;cZ1
yP1
exp);cZ1.ConstantFolding(eY3!cZ1
i61||!isEvenInteger(cZ1
nC1)){goto
cannot_adopt_mul;}
}
}
}
eD1
mul,value,exponent,c31
xH
t43
cannot_adopt_mul:
#endif
{if(mul.nX2(value,xG2==CollectionSetBase::xP2)c41}
}
}
xN1
ConstantFolding_MulGrouping(e62{bool
xH
tO3
bool
should_regenerate
tO3
eV1
mul;eD1
mul,tree,nF1
1)),c31
xH);typedef
std::pair<eJ,eH>eE1;typedef
std::multimap<fphash_t,eE1>cL1;cL1
tD;xX2
eV1::xS1
j=mul.iC.yU3
j!=mul.iC.end();++j){eJ&value=j
e82.value;eJ&e92
j
e82.xT2;if(j
e82.cX)exponent
eT1
const
fphash_t
eF1=exponent.GetHash();xK3
cL1::yV3
i=tD.xU2
eF1);for(;i!=tD.cY1
eF1;++i)if(i
e82.first
xF
xG2){if(!exponent
i61||!e01
nC1
tG1
1)))c41
i
e82
tE3.push_back(value);goto
skip_b;}
tD.yR3,std::make_pair(eF1,std::make_pair
t62,eH
iA1(1),value))));skip_b:;}
#ifdef FP_MUL_COMBINE_EXPONENTS
ConstantExponentCollection
xK
e11;xX2
cL1::yV3
j,i=tD.yU3
i!=tD.end();i=j){j=i;++j;eE1&list=i
e82;if(list.first
l71
e92
list.first
nC1;if(!t62==xD1)e11.MoveToSet_Unique
t62,list
iO2;tD.erase(i);}
}
if(e11.Optimize())c41
#endif
if(should_regenerate){eJ
before=tree;before.l61
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_MulGrouping: "
i13(before
lT2"\n"
;
#endif
tree.DelParams();xX2
cL1::yV3
i=tD.yU3
i!=tD.end();++i){eE1&list=i
e82;
#ifndef FP_MUL_COMBINE_EXPONENTS
if(list.first
l71
e92
list.first
nC1;if
t62==xD1
continue;if(e01
xE3
tree.AddParamsMove(list
iO2;eR2}
#endif
eJ
mul;mul
tH3
mul
iO3
list
iO2;mul
eT1
if(xH&&list.first
i51
list.first
nC1==iV2(1)/iV2(3)){eJ
cbrt;cbrt
xD
cCbrt);cbrt.e4
cbrt
cK
cbrt
xW
0.5)){eJ
sqrt;sqrt
xD
cSqrt);sqrt.e4
sqrt
cK
sqrt
xW-0.5)){eJ
rsqrt;rsqrt
xD
cRSqrt);rsqrt.e4
rsqrt
cK
rsqrt
xW-1)){eJ
inv;inv
xD
cInv);inv.e4
inv
cK
inv);eR2}
eJ
pow;pow
xD
cPow);pow.e4
pow
yP1
list.first);pow
cK
pow);}
#ifdef FP_MUL_COMBINE_EXPONENTS
tD.clear()yM3
0;a<i3.tF3;++a
tJ2
e92
i3[a
tC3;if(e01
xE3
tree.AddParamsMove(i3[a]iO2;eR2
eJ
mul;mul
tH3
mul
iO3
i3[a]iO2;mul
eT1
eJ
pow;pow
xD
cPow);pow.e4
pow.yC
xG2);pow
cK
pow);}
#endif
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_MulGrouping: "
yA1
#endif
return!tree
xF
before);}
return
tR3
xN1
ConstantFolding_AddGrouping(e62{bool
should_regenerate
tO3
eV1
add;nB1{if
y91
a)nC==cMul
cP2;if(add.nX2
y91
a))==CollectionSetBase::xP2)c41}
std::eV2
remaining(eM);size_t
tE=0;nB1{lC2
mulgroup=t33;if
lR3
nC==cMul){i81
yQ1
tP3{if
lR3
l8
b)i61
cP2;xK3
eV1::xS1
c=add.FindIdenticalValueTo
lR3
l8
b)eY3
add.Found(c)){eJ
tmp
lR3
yA
CloneTag());tmp
tK2
b);tmp
eT1
add.AddCollectionTo(tmp,c);c41
goto
done_a;}
}
remaining[a]=true;tE+=1;done_a:;}
}
if(tE>0){if(tE>1){std
xV3<std::pair<eJ,size_t> >nU;std::multimap<fphash_t,size_t>eG1;bool
lS3
tO3
nB1
iC1{i81
t33.tP3{lC2
p=t33
l8
b);const
fphash_t
p_hash=p.GetHash();for(std::multimap<fphash_t,size_t>::const_iterator
i=eG1.xU2
p_hash);i!=eG1.cY1
p_hash;++i){if(nU[i
e82
tC3
xF
p)){nU[i
e82
tD3+=1;lS3=true;goto
found_mulgroup_item_dup;}
}
nU.push_back(std::make_pair(p,size_t(1)));eG1.insert(std::make_pair(p_hash,nU.tF3-1));found_mulgroup_item_dup:;}
}
if(lS3){eJ
eA2;{size_t
max=0;for
iA1
p=0;p<nU.tF3;++p)if(nU[p
tD3<=1)nU[p
tD3=0;else{nU[p
tD3*=nU[p
tC3.y22;if(nU[p
tD3>max){eA2=nU[p
tC3;max=nU[p
tD3;}
}
}
eJ
group_add;group_add
xD
cAdd);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Duplicate across some trees: "
i13(eA2
lT2" in "
yA1
#endif
nB1
iC1
i81
t33.tP3
if(eA2
xF
t33
l8
b))){eJ
tmp
y91
a)yA
CloneTag());tmp
tK2
b);tmp
eT1
group_add
yP1
tmp);remaining[a]tO3
yF3
group_add
eT1
eJ
group;group
tH3
group
yP1
eA2);group
yP1
group_add);group
eT1
add.nX2(group);c41}
}
nB1
iC1{if(add.nX2
y91
a))==CollectionSetBase::xP2)c41}
}
if(should_regenerate){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_AddGrouping: "
yA1
#endif
tree.DelParams();xX2
eV1::xS1
j=add.iC.yU3
j!=add.iC.end();++j){eJ&value=j
e82.value;eJ&coeff=j
e82.xT2;if(j
e82.cX)coeff
eT1
if(coeff
i51
cB3
coeff
nC1,xD1
cP2
cA3
coeff
nC1
xE3
lP1
value);eR2}
eJ
mul;mul
tH3
mul
yP1
value);mul
yP1
coeff);mul
cK
mul);}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_AddGrouping: "
yA1
#endif
return
true;}
return
tR3}
lF3{using
lF3
FUNCTIONPARSERTYPES;using
t5
nJ1
bool
ConstantFolding_IfOperations(iS3(tree.GetOpcode()==cIf||tree.GetOpcode()==cAbsIf);for(;;){lT3
nC==cNot){tI3
cIf);lE2
0).xY2
0)lR2
eO3.swap
y91
2)t43
lT3
cS1{tI3
yW3;lE2
0).xY2
0)lR2
eO3.swap
y91
2)t43
yJ3
tM
0),i01==yW3)l92
tree.xY2
1));nV
lC3
tree.xY2
2));nV
lT1
lT3
nC==cIf||lE2
0)nC==yW3{eJ
cond=lE2
0);eJ
lU3;lU3
t53==cIf?cNotNot:cAbsNotNot);lU3
xZ2
1));ConstantFolding(lU3);eJ
lV3;lV3
t53==cIf?cNotNot:cAbsNotNot);lV3
xZ2
2));ConstantFolding(lV3
eY3
lU3
i61||lV3
i61){eJ
eV;eV
t53);eV
xZ2
1));eV.nJ
1));eV.nJ
2));eV
eT1
eJ
eW;eW
t53);eW
xZ2
2));eW.nJ
1));eW.nJ
2));eW
eT1
tree
t53);cM1
0,cond
lR2
tree.nG1
1,eV);tree.nG1
2,eW)xJ2}
if
y91
1)nC==lE2
2)nC&&y91
1)nC==cIf||eO3
nC==yW3){eJ&iE2=eO3;eJ&leaf2=lE2
2
eY3
iE2
l8
0)xC1
0))&&c93
1))||iE2
l8
2)xC1
2)))){eJ
eV;eV
iP2
eV.nJ
0));eV
c9
iE2
l8
1));eV
c9
leaf2
l8
1));eV
eT1
eJ
eW;eW
iP2
eW.nJ
0));eW
c9
iE2
l8
2));eW
c9
leaf2
l8
2));eW
iU
SetParam(0,iE2
lR2
tree.nG1
1,eV);tree.nG1
2,eW)xJ2
if
c93
1))&&iE2
l8
2)xC1
2))){eJ
eX;eX
iP2
eX.AddParamMove
y91
0));eX
c9
iE2
lR2
eX
c9
leaf2
lR2
eX
iU
nG1
0,eX);cM1
2,iE2
l8
2));cM1
1,iE2
l8
1))xJ2
if
c93
2))&&iE2
l8
2)xC1
1))){eJ
eB2;eB2
xD
leaf2
nC==cIf?cNot:cZ3);eB2
c9
leaf2
lR2
eB2
eT1
eJ
eX;eX
iP2
eX.AddParamMove
y91
0));eX
c9
iE2
lR2
eX
yP1
eB2);eX
iU
nG1
0,eX);cM1
2,iE2
l8
2));cM1
1,iE2
l8
1))xJ2}
eJ&xX=eO3;eJ&y6=lE2
2
eY3
xX
xF
y6)){tree.xY2
1))xJ2
const
OPCODE
op1=xX
nC;const
OPCODE
op2=y6
nC;if
eN3
op2){if(xX
cE1
1){eJ
lO
0));lW3
lR2
tQ3
n4
if(xX
cE1
2&&y6
cE1
2){if(xX
l8
0)xF
y6
l8
0))){eJ
param0=xX
l8
0);eJ
lO
1));lW3
l8
1));tQ3;lP1
param0)n4
if(xX
l8
1)xF
y6
l8
1))){eJ
param1=xX
l8
1);eJ
lO
0));lW3
lR2
tQ3;lP1
nE1);lP1
param1)xJ2}
if
eN3
yZ3
cMul
lS2
cAnd
lS2
cOr
lS2
cAbsAnd
lS2
cAbsOr
lS2
cMin
lS2
cMax){eH
lX3;yZ{for
iA1
b=y6.xF1
b-->0;){if
c62
y6
l8
b))){if(lX3
cT3){xX.l61
y6.l61}
lX3.push_back(xX
x13
y6
tK2
b);xX
t31
cT2}
}
if(!lX3
cT3){xX
eT1
y6.Rehash()l7
op1);tree
iO3
lX3)n4}
}
if
eN3
yZ3
cMul||eN3
cAnd
nT1
y6))||eN3
cOr
nT1
y6))){yZ
if
c62
y6)){xX.l61
xX
t31);xX
eT1
eJ
cO1=y6;y6=tF
op1==yZ3
cOr
tE2
op1);lP1
cO1)n4}
if(eN3
cAnd
lS2
cOr)&&op2==cNotNot){eJ&lY3=y6
l8
0);yZ
if
c62
lY3)){xX.l61
xX
t31);xX
eT1
eJ
cO1=lY3;y6=tF
op1
yY3
op1);lP1
cO1)n4}
if(op2==cAdd||op2==cMul||(op2==cAnd
nT1
xX))||(op2==cOr
nT1
xX))){for
iA1
a=y6
n41
y6
l8
a)xF
xX)){y6.l61
y6
t31);y6
eT1
eJ
cP1=xX;xX=tF
op2==cAdd||op2
yY3
op2);lP1
cP1)n4}
if((op2==cAnd||op2==cOr)&&op1==cNotNot){eJ&lZ3=xX
l8
0)yM3
y6
n41
y6
l8
a)xF
lZ3)){y6.l61
y6
t31);y6
eT1
eJ
cP1=lZ3;xX=tF
op2
yY3
op2);lP1
cP1)n4}
return
tR3}
#include <limits>
lF3{using
lF3
FUNCTIONPARSERTYPES;using
t5
nJ1
int
maxFPExponent()yR
std::numeric_limits
xK::max_exponent;}
xN1
x51
iV2
base,iV2
xG2{if(base<xD1
return
true
cA3
base,xD1||cB3
base
tG1
1))cL
return
exponent>=iV2(maxFPExponent
xK())/fp_log2(base);}
xN1
ConstantFolding_PowOperations(iS3(tree.GetOpcode()==cPow);nR&&eO3
l71
const_value
t03
lR,tH);xO
const_value)n31
tR3
tM2
cB3
tH
xE3
tree.xY2
0))xJ2
nR&&cB3
lR
xE3
xO
1)n31
tR3
nR&&eO3
nC==cMul){bool
y02=false
iU2
lU2=lR;eJ
mulgroup=eO3
yM3
mulgroup
n41
mulgroup
l8
a)l71
imm=mulgroup
l8
a)nC1;{if(x51
lU2,imm))break
iU2
lV2
t03
lU2,imm)cA3
lV2,xD1)break;if(!y02){y02=true;yQ1
l61}
lU2=lV2;yQ1
DelParam(a
cT2}
if(y02){yQ1
Rehash();
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before pow-mul change: "
yA1
#endif
lE2
0).Become(e21
lU2));eO3.Become
lR3);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After pow-mul change: "
yA1
#endif
}
}
tM2
lE2
0)nC==cMul
tJ2
lW2=tH
iU2
y12=1.0;bool
y02
tO3
eJ&mulgroup=lE2
0)yM3
mulgroup
n41
mulgroup
l8
a)l71
imm=mulgroup
l8
a)nC1;{if(x51
imm,lW2))break
iU2
eH1
t03
imm,lW2)cA3
eH1,xD1)break;if(!y02){y02=true;yQ1
l61}
y12*=eH1;yQ1
DelParam(a
cT2}
if(y02){yQ1
Rehash();eJ
eS3;eS3
xD
cPow);eS3
iO3
tree.l02));eS3.n62
tree
tH3
lP1
eS3);i83
e21
y12))xJ2}
lT3
nC==cPow&&eO3
i61&&lE2
0)l8
1)l71
a=lE2
0)l8
1)nC1
iU2
b=tH
iU2
c=a*b;if(isEvenInteger(a)&&!isEvenInteger(c)){eJ
n03;n03
xD
cAbs);n03.nJ
0)lR2
n03
eT1
tree.nG1
0,n03
t43
cM1
0,lE2
0)lR2
cM1
1,e21
c));}
return
tR3}
lF3{using
lF3
FUNCTIONPARSERTYPES;using
t5;e12
l5{enum
eC2{MakeFalse=0,MakeTrue=1,t72=2,n23=3,MakeNotNotP0=4,MakeNotNotP1=5,MakeNotP0=6,MakeNotP1=7,xB=8
xF3
lX2{Never=0,Eq0=1,Eq1=2,c03=3,c13=4}
;eC2
if_identical;eC2
lY2
4];e12{eC2
what:4;lX2
when:4;}
iK1,iL1,iM1,iN1
nJ1
eC2
Analyze
lQ3
a,lC2
b
e83{if(a
xF
b
nE3
if_identical;c23
p0=iN
a);c23
p1=iN
b);eR3
known&&p1
n51){eR3
val<p1
tF2&&lY2
0]i4
0];eR3
val<=p1
tF2&&lY2
1]i4
1];}
if
iC3
p1
xF2{if(p0
tF2>p1
yN
val&&lY2
2]i4
2];if(p0
tF2>=p1
yN
val&&lY2
3]i4
3];}
if(IsLogicalValue(a)){if(iK1
tP2
iK1.when,p1
nE3
iK1.what;if(iM1
tP2
iM1.when,p1
nE3
iM1.what;}
if(IsLogicalValue(b)){if(iL1
tP2
iL1.when,p0
nE3
iL1.what;if(iN1
tP2
iN1.when,p0
nE3
iN1.what;}
return
xB;}
cA1
bool
TestCase(lX2
when,const
c23&p){if(!p
n51||!p
yN
known
cL
switch(when
y33
Eq0
nV1==iV2(0.0)&&e13==p
tF2;case
Eq1
nV1==iV2(1.0)&&e13==e13;case
c03
nV1>y21&&e13<=iV2(1);case
c13
nV1>=y21
lH1
1);c73;}
return
tR3}
;lF3
RangeComparisonsData{static
const
l5
Data[6]={{l5
c51,{l5::tU
xB,l5::tU
xB
iT
Eq1
iW
Eq1
nX1
Eq0
nY1
Eq0}
eI1
c51,l5::xB,l5
c51,l5::xB
iT
Eq0
iW
Eq0
nX1
Eq1
nY1
Eq1}
eI1
c51,l5::t72,l5::tU
MakeFalse
nX1
c03
iW
c13
yD
cC
yD
n13{l5::xB,l5
c51,l5::tU
n23
nX1
c13
iW
c03
yD
cC
eI1::tU
tU
n13
l5::t72
iT
c13
nY1
c03
yD
cC
yD
n13{l5::tU
n23,l5::xB,l5
c51
iT
c03
nY1
c13
yD
cC}
}
;}
xN1
ConstantFolding_Comparison(e62{using
lF3
RangeComparisonsData;assert(tree.GetOpcode()>=cEqual&&tree.GetOpcode()<=cGreaterOrEq);switch(Data[i01-cEqual].Analyze
y91
0),eO3)y33
l5::MakeFalse:xO
0);nV
l5
c51:xO
1)yN3
n23:tI3
cEqual)yN3
t72:tI3
cNEqual)yN3
MakeNotNotP0:tI3
cNotNot)yO3
1)yN3
MakeNotNotP1:tI3
cNotNot)yO3
0)yN3
MakeNotP0:tI3
cNot)yO3
1)yN3
MakeNotP1:tI3
cNot)yO3
0)yN3
xB:;}
if
y91
1)i61)switch
y91
0)nC
y33
cAsin:lM
fp_sin
c53
cAcos:lM
fp_cos
tG));tI3
i01==cLess?cGreater:i01==cLessOrEq?cGreaterOrEq:i01==cGreater?cLess:i01==cGreaterOrEq?cLessOrEq:i01);nV
cAtan:lM
fp_tan
c53
cLog:lM
fp_exp
c53
cSinh:lM
fp_asinh
c53
cTanh:if(fp_less(fp_abs
tG
xE3
lM
fp_atanh
tG))xJ2
break;c73
yF3
return
tR3}
#include <list>
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
using
lF3
FUNCTIONPARSERTYPES;lF3{
#ifdef DEBUG_SUBSTITUTIONS
yE
double
d){union{double
d;uint_least64_t
h;}
t82
d=d;lM1
h
nZ1
#ifdef FP_SUPPORT_FLOAT_TYPE
yE
float
f){union{float
f;uint_least32_t
h;}
t82
f=f;lM1
h
nZ1
#endif
#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
yE
long
double
ld){union{long
double
ld;e12{uint_least64_t
a;iF2
short
b;}
s;}
t82
ld=ld;lM1
s.b<<data.s.a
nZ1
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
yE
long
ld){o<<"("
<<std::hex<<ld
nZ1
#endif
#endif
}
t5{lN
nE)){}
lN
const
iV2&i
yA
xO3
nE
i)){data
x7
#ifdef __GXX_EXPERIMENTAL_CXX0X__
lN
iV2&&i
yA
xO3
nE
std::move(i))){data
x7
#endif
lN
iF2
v
yA
VarTag
nE
l03,v
xQ2
xB2
o
yA
OpcodeTag
nE
o
xQ2
xB2
o,iF2
f
yA
FuncOpcodeTag
nE
o,f
xQ2
lC2
b
yA
CloneTag
nE*b.data)){}
xW3
eJ::~CodeTree(){}
lB
ReplaceWithImmed
cG1
i){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Replacing "
i13(*this
eY3
IsImmed())OutFloatHex(std::cout,GetImmed()lT2" with const value "
<<i;OutFloatHex(std::cout,i
lT2"\n"
;
#endif
data=new
xH2
xK(i);}
xW3
e12
ParamComparer{iG2()lQ3
a,lC2
b
e83{if(a.y22!=b.y22)return
a.y22<b.y22
n31
a.GetHash()<b.GetHash();}
}
xR3
xH2
xK::Sort(y13
Opcode
y33
cAdd:case
cMul:case
cMin:case
cMax:case
cAnd:case
cAbsAnd:case
cOr:case
cAbsOr:case
cHypot:case
cEqual:case
cNEqual:std::sort(iY2
iN2
iY2
end(),ParamComparer
xK());lC
cLess
lY
cGreater;}
lC
cLessOrEq
lY
cGreaterOrEq;}
lC
cGreater
lY
cLess;}
lC
cGreaterOrEq
lY
cLessOrEq;}
break;c73
yF3}
lB
AddParam
lQ3
param){xY.push_back(param);}
lB
AddParamMove(eJ&param){xY.push_back(eJ());xY.back().swap(param);}
lB
SetParam
iA1
which,lC2
b)nW1
which
iX2
xY[which]=b;}
lB
nG1
size_t
which,eJ&b)nW1
which
iX2
xY[which
i63
b);}
lB
AddParams
cH1
nK){xY.insert(xY.end(),lZ2.iN2
lZ2.end());}
lB
AddParamsMove(nK){size_t
endpos=xY.tF3,added=lZ2.tF3;xY.t73
endpos+added,eJ());for
iA1
p=0;p<added;++p)xY[endpos+p
i63
lZ2[p]);}
lB
AddParamsMove(nK,size_t
n02)nW1
n02
iX2
DelParam(n02);AddParamsMove(tA1}
lB
SetParams
cH1
nK){eH
tmp(tA1
xY.i72}
lB
SetParamsMove(nK){xY.swap(tA1
lZ2.clear();}
#ifdef __GXX_EXPERIMENTAL_CXX0X__
lB
SetParams(eH&&lZ2){SetParamsMove(tA1}
#endif
lB
DelParam
iA1
index){eH&Params=xY;
#ifdef __GXX_EXPERIMENTAL_CXX0X__
iY2
erase(iY2
begin()+index);
#else
Params[index].data=0;for
iA1
p=index;p+1<xX3;++p)Params[p].data.UnsafeSetP(&*Params[p+1
iX2
Params[xX3-1].data.UnsafeSetP(0);iY2
t73
xX3-1);
#endif
}
lB
DelParams(){xY.clear();}
xN1
eJ::IsIdenticalTo
lQ3
b
e83{if(&*data==&*b.data)return
true
n31
data->IsIdenticalTo(*b.data);}
xN1
xH2
xK::IsIdenticalTo
cH1
xH2
xK&b
e83{if(Hash!=b.Hash
cL
if(Opcode!=tK3
cL
switch(Opcode
y33
cImmed:return
cB3
Value,tL3;case
l03:return
l22==b.l12
case
cFCall:case
cPCall:if(l22!=b.l22
cL
break;c73
yF3
if(xX3!=b.xX3
cL
for
iA1
a=0;a<xX3;++a){if(!Params[a]xF
b.Params[a])cL}
return
true;}
lB
Become
lQ3
b){if(&b!=this&&&*data!=&*b.data){DataP
tmp=b.data;l61
data.i72}
}
lB
CopyOnWrite(){if(GetRefCount()>1)data=new
xH2
xK(*data);}
xW3
eJ
eJ::GetUniqueRef(){if(GetRefCount()>1)return
eJ(*this,CloneTag())n31*this;}
tX):yT
cNop
t92(),n8
tX
const
xH2&b):yT
tK3
t92(tL3,l22(b.cQ1,Params(b.Params),Hash(b.Hash),Depth(b.Depth),tP1
b.l32){}
tX
const
iV2&i):yT
cImmed
t92(i),n8
#ifdef __GXX_EXPERIMENTAL_CXX0X__
tX
xH2
xK&&b):yT
tK3
t92
c72
tL3),l22(b.cQ1,Params
c72
b.Params)),Hash(b.Hash),Depth(b.Depth),tP1
b.l32){}
tX
iV2&&i):yT
cImmed
t92
c72
i)),n8
#endif
tX
xB2
o):yT
o
t92(),n8
tX
xB2
o,iF2
f):yT
o
t92(),l22(f),Params(),Hash(),Depth(1),tP1
0){}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <iostream>
using
lF3
FUNCTIONPARSERTYPES;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
lF3{tZ1
n12(nO,std
c7&done,std::ostream&o){nB1
n12
l64
done,o);std::ostringstream
buf
i13(tree,buf);done[tree.GetHash()].insert(buf.str());}
}
#endif
t5{
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
tZ1
DumpHashes(cF){std
c7
done;n12(tree,done,o);for(std
c7::const_iterator
i=done.yU3
i!=done.end();++i){const
std::set<std
cC3>&flist=i
e82;if(flist.tF3!=1)o<<"ERROR - HASH COLLISION?\n"
;for(std::set<std
cC3>::const_iterator
j=flist.yU3
j!=flist.end();++j){o<<'['<<std::hex<<i->first.hash1<<','<<i->first.hash2<<']'<<std::dec;o<<": "
<<*j<<"\n"
;}
}
}
tZ1
DumpTree(cF){nH3
iQ3;switch(i01
y33
cImmed:o<<eC3
c63
l03:o<<"Var"
<<(tree.GetVar()-l03)c63
cAdd:iQ3"+"
;lC
cMul:iQ3"*"
;lC
cAnd:iQ3"&"
;lC
cOr:iQ3"|"
;lC
cPow:iQ3"^"
;break;c73
iQ3;o<<c83
i01);l83
cFCall||i01==cPCall)o<<':'<<tree.GetFuncNo();}
o<<'(';if(eM<=1&&sep2[1])o<<(sep2+1)<<' ';nB1{if(a>0)o<<' 'i13
l64
o
eY3
a+1<eM)o<<sep2;}
o<<')';}
tZ1
DumpTreeWithIndent(cF,const
std
cC3&indent){o<<'['<<std::hex<<(void*)(&tree.l02))<<std::dec<<','<<tree.GetRefCount()<<']';o<<indent<<'_';switch(i01
y33
cImmed:o<<"cImmed "
<<eC3;o<<'\n'c63
l03:o<<"VarBegin "
<<(tree.GetVar()-l03);o<<'\n'n31;c73
o<<c83
i01);l83
cFCall||i01==cPCall)o<<':'<<tree.GetFuncNo();o<<'\n';}
nB1{std
cC3
ind=indent;for
iA1
p=0;p<ind.tF3;p+=2)if(ind[p]=='\\')ind[p]=' ';ind+=(a+1<eM)?" |"
:" \\"
;DumpTreeWithIndent
l64
o,ind);}
o<<std::flush;}
#endif
}
#endif
using
lF3
l21;using
lF3
FUNCTIONPARSERTYPES;
#include <cctype>
lF3
l21{iF2
ParamSpec_GetDepCode
cH1
e42&b
y13
b.first
y33
tG3:{cO*s=(cO*)b
tE3
n31
s->depcode;}
case
SubFunction:{cP*s=(cP*)b
tE3
n31
s->depcode;}
c73
yF3
return
0;}
tZ1
DumpParam
cH1
e42&yG2
std::ostream&o){static
const
char
ParamHolderNames[][2]={"%"
,"&"
,"x"
,"y"
,"z"
,"a"
,"b"
,"c"
}
;iF2
y32
0;switch(yH2.first
y33
NumConstant:{const
ParamSpec_NumConstant
xK
yK3
cH1
ParamSpec_NumConstant
xK*c61
using
lF3
FUNCTIONPARSERTYPES;o.precision(12);o<<yL3
constvalue;yF3
case
tG3:{cO
yK3(cO*c61
o<<ParamHolderNames[yL3
index];y32
yL3
constraints;yF3
case
SubFunction:{cP
yK3(cP*c61
y32
yL3
constraints;yF
GroupFunction){if
y23
l81==cNeg){o<<"-"
;n2}
lJ2
yL3
l81==cInv){o<<"/"
;n2}
else{std
cC3
opcode=c83(xB2)yL3
l81).substr(1)yM3
0;a<opcode.tF3;++a)opcode[a]=(char)std::toupper(opcode[a]);o<<opcode<<"( "
;n2
o<<" )"
;}
}
else{o<<'('<<c83(xB2)yL3
l81)<<' ';yF
PositionalParams)o<<'[';yF
SelectedParams)o<<'{';n2
if
y23
data.n1!=0)o<<" <"
<<yL3
data.n1<<'>';yF
PositionalParams)o<<"]"
;yF
SelectedParams)o<<"}"
;o<<')';}
yF3
l73
ImmedConstraint_Value(constraints&ValueMask)y33
ValueMask:lC
Value_AnyNum:lC
nY2:o<<"@E"
;lC
Value_OddInt:o<<"@O"
;lC
tS1:o<<"@I"
;lC
Value_NonInteger:o<<"@F"
;lC
eJ1:o<<"@L"
;yJ3
ImmedConstraint_Sign(constraints&SignMask)y33
SignMask:lC
Sign_AnySign:lC
nH1:o<<"@P"
;lC
eK1:o<<"@N"
;yJ3
ImmedConstraint_Oneness(constraints&OnenessMask)y33
OnenessMask:lC
Oneness_Any:lC
Oneness_One:o<<"@1"
;lC
Oneness_NotOne:o<<"@M"
;yJ3
ImmedConstraint_Constness(constraints&ConstnessMask)y33
ConstnessMask:lC
tR1:if(yH2.first==tG3){cO
yK3(cO*c61
if
y23
index<2)yF3
o<<"@C"
;lC
Constness_NotConst:o<<"@V"
;lC
Oneness_Any:yF3}
tZ1
DumpParams
nT2
paramlist,iF2
count,std::ostream&o){for(eQ1=0;a<count;++a){if(a>0)o<<' ';const
e42&param=e81
xK(paramlist,a);DumpParam
xK(param,o);iF2
depcode=ParamSpec_GetDepCode(param
eY3
depcode!=0)o<<"@D"
<<depcode;}
}
}
#include <algorithm>
using
lF3
l21;using
lF3
FUNCTIONPARSERTYPES;lF3{cO
plist_p[37]={{2,0,0x0}
nP
0,0x4}
nP
nH1,0x0}
nP
eK1|Constness_NotConst,0x0}
nP
Sign_NoIdea,0x0}
nP
eJ1,0x0}
,{3,Sign_NoIdea,0x0}
,{3,0,0x0}
,{3,eJ1,0x0}
,{3,0,0x8}
,{3,Value_OddInt,0x0}
,{3,Value_NonInteger,0x0}
,{3,nY2,0x0}
,{3,nH1,0x0}
,{0,eK1|lV{0,lV{0,nH1|lV{0,nY2|lV{0,tR1,0x1}
,{0,tS1|nH1|lV{0,tT1
tR1,0x1}
,{0,tT1
lV{0,Oneness_One|lV{0,eJ1|lV{1,lV{1,nY2|lV{1,tT1
lV{1,tS1|lV{1,nH1|lV{1,eK1|lV{6,0,0x0}
,{4,0,0x0}
,{4,tS1,0x0}
,{4,lV{4
tU3
5,0,0x0}
,{5,lV}
nJ1
e12
plist_n_container{static
const
ParamSpec_NumConstant
xK
plist_n[20];}
nJ1
const
ParamSpec_NumConstant
xK
plist_n_container
xK::plist_n[20]={{iV2(-2
tY-1
tY-0.5
tY-0.25
tY
0
tA2
fp_const_deg_to_rad
xK(tA2
fp_const_einv
xK(tA2
fp_const_log10inv
xK(tY
0.5
tA2
fp_const_log2
xK(tY
1
tA2
fp_const_log2inv
xK(tY
2
tA2
fp_const_log10
xK(tA2
fp_const_e
xK(tA2
fp_const_rad_to_deg
xK(tA2-fp_const_pihalf
xK(),xT1{y21,xT1{fp_const_pihalf
xK(),xT1{fp_const_pi
xK(),xT1}
;cP
plist_s[517]={{{1,15,tS2
398,tS2
477,tS2
15,cNeg,GroupFunction,0}
,tR1,0x1}
,{{1,15,y42
24,y42
465,y42
466,y42
498,cInv,lT
2,327995
xE
l0
2,48276
xE
l6
260151
xE
l6
470171
xE
l6
169126
xE
l6
48418
xE
l6
1328
xE
l6
283962
xE
l6
169275
xE
l6
39202
xE
l6
283964
xE
l6
283973
xE
l6
476619
xE
l6
296998
xE
l6
47
xE
SelectedParams,0}
,0,0x4
nH
161839
xE
l6
25036
xE
l6
35847
xE
l6
60440
xE
l6
30751
xE
l6
270599
xE
l6
60431
xE
l6
259119
xE
l6
183474
xE
l6
332066
xE
l6
7168
xE
l6
197632
xE
l6
291840
xE
l6
283648
xE
l6
238866
xE
l6
239902
xE
l6
31751
xE
l6
244743
xE
l6
384022
xE
SelectedParams,0}
,0,0x4
nH
385262
xE
l6
386086
xE
l6
393254
xE
SelectedParams,0}
,0,0x5
nH
393254
xE
l6
386095
xE
l6
387312
xE
l6
18662
xE
l6
61670
xE
l6
387397
xE
l6
247855
xE
SelectedParams,0}
,0,0x1
nH
342063
xE
l6
297007
xE
l6
15820
xE
l6
393263
xE
l6
393263
xE
SelectedParams,0}
,0,0x5
nH
161847
xE
l6
258103
xE
l6
249073
xE
l6
249076
xE
iD
0,0
xE
nF
0,0
tU1
1,45
xE
nF
1,53
xE
nF
1,54
xE
nF
1,55
xE
nF
1,56
xE
nF
1,26
xE
nF
1,259
tJ
1}
tU3{1,253
xE
nF
1,272
tU1
1,323
tJ
1}
tU3{1,0
xE
nF
1,21
xE
nF
1,447
tJ
1}
lA2
449
tJ
1}
lA2
0
tJ
1}
lA2
0
tJ
2}
lA2
15
xE
nF
1,24
tJ
2}
,0,0x0
nH
58392
tU1
0,0
tJ
1}
,nH1,0x0
nH
24591
yI3
33807
yI3
48143
yI3
285720
yI3
290840
yI3
305152,l9
312400,l9
39202,l9
121894,l9
421926,l9
429094,l9
443430,l9
317834,l9
329098,l9
7633,l9
7706,l9
7730,l9
38,l9
50587,l9
406528,l9
24583,l9
31751,l9
405511,l9
321551,xU1
327713,l9
322596,l9
88361,l9
335174,l9
327050,l9
493606,l9
496678,l9
503846,l9
516134,l9
7217,l9
333875,l9
336896,l9
524326,l9
509952,l9
286727,l9
90127,l9
131087,l9
296976,tB1
324623,l1
0x14
nH
332815,l1
0x10}
,{{3,7340056,tB1
289092,l9
92176,xU1
337935
lB2
7340060,l1
tT2
7340176,l9
338959
lB2
7340061,xU1
7206,l9
xJ3
l9
357414,l9
368678,l9
370745,l1
0x7}
,{{3,7340177,l9
39277,tB1
426398,l1
tT2
40272286,xU1
490910,l1
tT2
40336798,xU1
50600,l9
426462,xU1
490974,xU1
370726,l1
0x6
nH
371750,l1
0x6
nH
428070
lB2
40336862,xU1
38378,l9
50671
lB2
47662080,l9
477184,l9
568320,l9
371727,l1
0x7}
,{{3,15779306,l9
370703,l1
0x7
nH
39277,l9
39279,l1
0x4}
,{{3,15779238,l9
39338,tB1
436262,l9
508966,l9
39409,tB1
296998,tB1
35847,l9
15,tB1
377894,l9
386063,l1
0x1
nH
15,l9
7192,l9
122904,l9
121880,l9
30751,l9
57,l9
7456,l9
15674
lB2
67579935,l9
39237,l9
58768,l9
62924,l9
121856,l9
15760
lB2
64009216,l1
0x0}
,{{0,0,xG
0,0,iL
2,e31
2,e41
3,e31
3,e41
38,xG
1,38,iL
14,xG
1,57,xG
1,16,eD2
0x0
nH
471103,eD2
0x1}
,{{1,303,xG
1,323,cE3
0x0
nH
471363,eD2
0x16}
,{{1,293,e31
294,e41
295,xG
1,296,iL
400,xG
1,0,xG
1,460,xG
1,465,xG
1,16,eD2
0x1}
,{{1,57,cE3
0x1}
,{{1,0,iL
21,xG
1,15,eD2
0x0
nH
24591,xG
1,24,iL
517,cE3
0x0
nH
46095,lK
46104,lK
15397,lK
287789,lK
66584,lK
404763,lK
62504,lK
15409,lK
39951,lK
24591,lK
33807,lK
50200,lK
62509,lK
50176,lF,178176,eL1
tV3
283648,lF,19456,lF,27648,lF,89088,lF,86016,lF,488448,lF,14342,lF,58375,lF,46147
xZ
46151,lF,284679,lF,7183,lF,46159
xZ
38993
xZ
50265,lF,50249,lF,283808,lF,284835,lF,24822,lF,10240,lF,11264,lF,7170,lF,xJ3
lF,17408,lF,164864,lF,237568,lF,242688,eL1
0x14
nH
476160,lF,25607,lF,121871,lF,50252,lF,39374,lF,50183,lF,7192,lF,121887,lF,252979,lF,46155,lF,38919,lF,50267,lF,50268,lF,50253,lF,46190,lF,50295,lF,7563,eL1
0x10
nH
416811,lF,416819,lF,40046,lF,46191
xZ
415795,lF,40047
xZ
415787,lF,39015,eL1
0x5
nH
39326
xZ
39326,lF,39332,eL1
0x5
nH
39333
cW2
50590
xZ
50590,lF,39338
xZ
39338,lF,39335,eL1
0x5
nH
15786
xZ
146858,lF,39372,lF,39379,lF,39380,lF,39390
xZ
50654
xZ
50654,lF,24,eL1
0x6
nH
62,lF,24,lF,62,eL1
0x6
nH
43,lF,43
xZ
51,lF,51
xZ
50269,lF,50176
xZ
50270,lF,39159,lF,39183
xZ
7168
xZ
31744,lF,99328,lF,31746,lF,100376,lF,39409
xZ
39411
xZ
39411,lF,39420,lF,39420
xZ
15,lF,39025,eL1
0x5
nH
39422,lF,16384,lF,62853,lF,15360,lF,15
cW2
16,lF,7183
cW2
7172
l14
yX1}
,nH1,0x0
nH
24591
l14
lT
2,50200
l14
lT
2,63521
l14
lT
2,62500
l14
lT
2,50453
l14
lT
2,62488,cPow
eP3
tS3
7,tS3
194,tS3
0,cAcos
tT3
cAcosh
tT3
cAsin
tT3
cAsinh
nQ
119,cAsinh
tT3
cAtan
eW1
306176,cAtan2
eW1
xJ3
cAtan2
tT3
cAtanh
nQ
246,cCeil
tT3
cCeil,l13
0,c92
0,cCos,l13
7,c92
91,c92
92,c92
119,c92
236,c92
255,c92
214,l23
236,l23
464,l23
0,cCosh,l13
0,l23
0,cExp
nQ
7,cExp
nQ
91,cExp
tT3
cF3
7,cF3
91,cF3
246,cFloor
tT3
cFloor,lA
0x4
nH
309540,cHypot
eW1
316708,cHypot
eW1
316724,cHypot,l0
3,32513024,y52
34627584
n61
31493120,y52
89213952
n61
149042176
n61
246647808
n61
301234176
n61
494360576
n61
498558976
n61
62933520
n61
62933520,y52
62933526
n61
62933526,y52
24670208
n61
579378176
n61
573578240
n61
32513024
n61
566254592
n61
7900160
n61
588822528,cIf
nQ
119,cInt
nQ
246
cA2
0
cA2
7
cA2
31
cA2
194
cA2
363
cA2
15,cLog,lT
1,24,cLog
eP3
cLog10
tT3
cLog2
eW1
xJ3
cMax
eW1
35847,cMax
eW1
30751,cMax
tT3
cMax,AnyParams,1}
,0,0x4
nH
xJ3
cMin
eW1
35847,cMin
eW1
30751,cMin
tT3
cMin,AnyParams,1}
,0,0x4
nH
24591,cMin
eP3
nZ2
7,nZ2
91,nZ2
92,nZ2
119,nZ2
149,nZ2
231,cSin,lA
0x5}
,{{1,246,nZ2
255,nZ2
254,nZ2
0,cSin,l13
273,cSin,lA
0x1}
,{{1,214,y62
231,cSinh,lA
0x5}
,{{1,246,y62
254,y62
255,y62
464,y62
0,cSinh,l13
0,y62
15,cSqrt
eP3
cB2
0,cTan,l13
115,cTan,l13
116,cB2
231,cB2
246,cB2
273,cB2
254,cB2
255,cB2
0,y92
0,cTanh,l13
213,y92
231,y92
246,y92
254,y92
255,y92
0,cTrunc
eW1
15384,cSub,lT
2,15384,cDiv,lT
2,476626,cDiv,lT
2,121913,tU2
xJ3
n42
lA
tV3
xJ3
tU2
31744,n42
lA
0x20
nH
31751,n42
lA
0x24
nH
31751,tU2
121913,eM1
l0
2,xJ3
cLess,lA
tV3
41984,cLess,lA
0x4
nH
41984,cLess
eW1
7,cLess
eW1
xJ3
cLessOrEq
eW1
296182,cLessOrEq
eW1
7168
lJ3,lA
tV3
41984
lJ3,lA
0x4
nH
41984
lJ3
eW1
7
lJ3
eW1
xJ3
yB
l0
2,296182,cGreaterOrEq
tT3
n22
245,n22
7,n22
550,n22
553,n22
554,n22
556,n22
31,n22
559,n22
15,n22
560,cNot
eW1
7706,n33
xJ3
n33
35847,n33
30751,n33
463903,n33
466975,cAnd,iD
0,0,cAnd,nF
2,xJ3
eT3
7706,eT3
35847,eT3
463903,eT3
466975,eT3
30751,cOr,iD
1,0,n32
91,n32
131,n32
245,n32
215,n32
246,cDeg
nQ
246,cRad
eW1
xJ3
cAbsAnd,l6
xJ3
cAbsOr,iD
1,0,cZ3
tT3
cAbsNotNot,l0
3,32513024,cL3
lA
0x0}
,}
;}
lF3
l21{const
Rule
grammar_rules[262]={{ProduceNewTree,17,1,0,{1,0,cAbs,eE2
409,{1,146,cAtan,eE2
403
nP
1324,cG3
eE2
405
nP
307201,cG3
l2
c22
253174
nP
255224,cG3
l2
c22
259324
nP
257274,cG3
eE2
152,{1,252,cCeil
iO
486,{1,68,tW2
482,{1,122,tW2
483,{1,124,tW2
151,{1,125,tW2
419,{1,123,tW2
0,{1,403,cCos,l2
2,1,246,{1,252,cCos,l2
18,1,0,{1,400,tW2
301,{1,406,cCosh,l2
2,1,246,{1,252,cCosh,l2
18,1,0,{1,400,cCosh
iO
458,{1,121,cFloor,eE2
150,{1,252,cFloor,tW3
156,{3,7382016,eB
549,{3,8430592,eB
556,{3,8436736,eB
157,{3,42998784,eB
550,{3,42999808,eB
562,{3,43039744,eB
557,{3,49291264,eB
538,{3,49325056,eB
469,{3,1058318,eB
473,{3,1058324,eB
473,{3,9438734,eB
469,{3,9438740,cIf,l2
0,3,32542225,{3,36732434,cIf,l2
0,3,32542231,{3,36732440,cIf,tX3
573,{3,32513026,cIf,tX3
515,{3,455505423,cIf,tX3
515,{3,433506837,cIf
iO
78,{1,256,tY3
69,{1,258,tY3
404,{1,72,tY3
159,{1,147,cLog,l2
0,1,0
nP
487425,cMax
l3
16,1,445
nP
cH3
cMax
l3
0,1,0
nP
483329,cMin
l3
16,1,446
nP
cH3
cMin,c0
0,1,153
nP
24832
l14
tW3
153
nP
25854
l14
tW3
154
nP
129039
nX3
32055
nX3
32056
nX3
32057
tV2
166288
nP
32137
nX3
33082
tV2
7168
nP
12688
tV2
7434
nP
12553
y72
435
nP
46146
y72
436
nP
46154
y72
437
nP
46150
y72
169
nP
83983
y72
168
nP
130082
y72
175
nP
133154
y82
476160
nP
471055
y82
274432
nP
273423
y82
251904
nP
266274
y82
251904
nP
263186
y72
171,{1,252,cC2
421,{1,68,cC2
151,{1,122,cC2
419,{1,124,cC2
170,{1,125,cC2
482,{1,123,cC2
0,{1,405,cC2
172,{1,252,cSinh
iO
328,{1,404,cSinh
iO
173,{1,252,l34
0,{1,408,l34
176,{1,410,l34
177,{1,252,cTanh,l2
0,1,442
nP
449551,tZ
1,441
nP
cH3
tZ
1,167
nP
268549,tZ
1,181
nP
276749,tZ
1,180
nP
276500
iD3
190770
nP
189622
iD3
194748
nP
193723
iD3
202943
nP
196795
iD3
59699
nP
298148
iD3
59714
nP
325815
iD3
59724
nP
343224
xE
c0
2,1,337,{1,333
tJ
1
l33
336,{1,338
tJ
1}
}
l44
2,1,340
nP
1363,i0
342
nP
1365,i0
463
nP
472524,i0
47
nP
356711,i0
349
nP
200751,i0
360
nP
199727,i0
480
nP
207053,i0
481
nP
208077,i0
417
nP
211144,i0
209
nP
211145,i0
418
nP
215240,i0
212
nP
212329,i0
204
nP
373097,i0
211
nP
372944,i0
217
nP
201944,i0
221
nP
223448,i0
367
nP
508329,i0
219
nP
508126,i0
224
nP
225705,i0
223
nP
225776,i0
365
nP
230825,i0
426
nP
377057,i0
497
nP
377054,i0
497
nP
204201,i0
426
nP
375280,i0
224
nP
375006,cAdd
l3
2,2,407781
nP
233698,cAdd
l3
2,2,59763
nP
233842,tZ
1,372
nP
1397,lQ1
95
nP
24705,lQ1
96
nP
24708,lQ1
444
nP
449551,lQ1
443
nP
cH3
lQ1
100
nP
101750,lQ1
108
nP
106821,lQ1
105
nP
103749
nA2
0,2,110607
nP
108869
nA2
0,2,107535
nP
109893,lJ
0
l33
112
nP
111634,cMul,SelectedParams,0
l33
567,{1,52,lJ
1
l33
568,{1,42,lJ
1}
}
l44
2,1,467
nP
45516
iV
356
nP
51555
iV
468
nP
49612
iV
357
nP
47459
iV
429
nP
438699
iV
432
nP
441774
iV
486
nP
498726
iV
494
nP
504870
iV
382
nP
435579
iV
497
nP
435709
iV
426
nP
508287
iV
414
nP
500092
iV
499
nP
352744
iV
345
nP
367092
iV
381
nP
425318
iV
478
nP
425460
iV
47
nP
512501
iV
505
nP
355817
iV
47
nP
516598
iV
507
nP
518182
iV
508
nP
358896
iV
351
nP
388605
iV
511
nP
360939
iV
503
nP
354788
iV
514
nP
525350
iV
510
nP
394342
iV
386
nP
351346
nA2
2,2,363004
nP
361968
nA2
16,1,117
nP
1157
nA2
16,1,118
nP
1158
nA2
16,1,402
nP
411024
nA2
16,2,58768
nP
1472
nA2
16,2,15760
nP
1474
nA2
17,1,0,{1,400
nA2
17,1,57,{1,14,lJ
0
l91
4
n43
41,n42
l4
4,1,0
nP
5167,n42
yA2
41984
nP
409641,cEqual
c4
n42
yA2
eY
n42
l43
n42
l2
xB1
24849
tN2
iE
tN2
n63
281873
tN2
lA1
tN2
lB1
n42
l4
4
y43
41,eM1
l4
4
n43
5167,eM1
yA2
41984
nP
409641,cNEqual
c4
eM1
yA2
eY
eM1
l43
eM1
l2
xB1
24849,eM1
l2
cI3
eM1
l2
c22
n63
281873,eM1
l2
c22
lA1,eM1
l2
c22
lB1
cNEqual
c4
eO1
16,2,eY
eO1
16
eI
cLess,eE2
571
nP
46080,eO1
xB1
24832,eO1
c22
yT1,eO1
cI3
eO1
cJ3
eO1
c22
xW1,eO1
c22
lA1,eO1
c22
lB1
cLess,l4
20
y43
409641,cLess
c4
t6
16,2,eY
t6
16
eI
cLessOrEq,eE2
565
nP
409615,t6
xB1
24832,t6
c22
yT1,t6
cI3
t6
cJ3
t6
c22
xW1,t6
c22
lA1,t6
c22
lB1
cLessOrEq,l4
20
y43
409647,cLessOrEq
c4
lK3
yA2
eY
lK3
l43
lK3
eE2
539
nP
409615
eP
529654
nP
24832
eP
yT1
eP
iE
eP
n63
281856
eP
xW1
eP
lA1
eP
lB1
lK3
l4
20
n43
409647
lJ3
c4
yB
yA2
eY
yB
l43
yB
eE2
572
nP
46080,yB
l2
xB1
24832,yB
l2
c22
yT1,yB
l2
cI3
yB
l2
cJ3
yB
l2
c22
xW1,yB
l2
c22
lA1,yB
l2
c22
lB1
yB
l4
20
n43
409641,yB
l4
4,1,519,{1,137,cNot,tX3
571,{1,2,cNot,l2
0,1,452
nP
cH3
n73
0,2,537097,{3,547892744,cAnd,c0
16,1,566,{1,5,cAnd,AnyParams,1}
}
l44
16,1,569
nP
13314,n73
16,1,544
nP
553498,n73
16,1,546
nP
462369,n73
16,1,548
nP
466465,n73
0,1,457
nP
cH3
xX1
570
nP
13314,xX1
563
nP
8197,xX1
541
nP
553498,xX1
542
nP
462369,xX1
543
nP
466465,xX1
564
nP
143365,cOr,c0
4,1,525,{1,137,cK3
tX3
572,{1,2,cK3
l4
17,1,0,{1,0,cK3
eE2
537,{1,256,cAbsNotNot,c0
18,1,531,{1,254,cAbsNotNot,c0
0,1,572,{3,43039744,cL3
tW3
571,{3,49325056,cL3
tX3
454,{3,32513586,cL3
l2
16,3,32542225,{3,36732434,cL3
yX1}
}
,}
;e12
grammar_optimize_abslogical_type{y2
9
cS
grammar_optimize_abslogical_type
grammar_optimize_abslogical={9,{34,192,228,238,242,247,254,260,261}
}
;}
e12
grammar_optimize_ignore_if_sideeffects_type{y2
59
cS
grammar_optimize_ignore_if_sideeffects_type
grammar_optimize_ignore_if_sideeffects={59,{0,20,21,22,23,24,25,26,cQ
iO1
78,cR
cT
e12
grammar_optimize_nonshortcut_logical_evaluation_type{y2
56
cS
grammar_optimize_nonshortcut_logical_evaluation_type
grammar_optimize_nonshortcut_logical_evaluation={56,{0,25,cQ
iO1
78,cR
241,243,244,245,246,248,249,250,251,252,253,255,256,257,258,259}
}
;}
e12
grammar_optimize_recreate_type{y2
22
cS
grammar_optimize_recreate_type
grammar_optimize_recreate={22,{18,55,56,57,80,81,82,83,84,85,117,118,120,121,130,131,132,133,134,135,136,137}
}
;}
e12
grammar_optimize_round1_type{y2
125
cS
grammar_optimize_round1_type
grammar_optimize_round1={125,{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,19,25,cQ
37,38,iO1
45,46,47,48,49,50,51,52,53,54,58,59,60,61,62,63,64,65,66,67,68,69,70,71,78,79,80,81,82,83,84,85,86,87,88,93,94,95,96,97,98,99,100,101,117,118,119,120,121,122,123,124,125,126,127,128,129,138,160,161,162,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,cT
e12
grammar_optimize_round2_type{y2
103
cS
grammar_optimize_round2_type
grammar_optimize_round2={103,{0,15,16,17,25,cQ
39,40,iO1
45,46,47,48,49,50,51,52,53,54,59,60,72,73,78,79,86,87,88,89,90,91,92,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,119,122,123,124,125,126,127,128,139,159,160,161,162,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,cT
e12
grammar_optimize_round3_type{y2
79
cS
grammar_optimize_round3_type
grammar_optimize_round3={79,{74,75,76,77,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,170,171,172,173,174,175,176,177,181,182,183,184,185,186,187,188,189,190,191,193,194,195,196,197,198,199,201,202,203,205,206,207,208,209,210,211,213,214,215,217,218,219,220,221,222,223,225,226,227,229,230,231,232,233,234,235}
}
;}
e12
grammar_optimize_round4_type{y2
12
cS
grammar_optimize_round4_type
grammar_optimize_round4={12,{18,55,56,57,130,131,132,133,134,135,136,137}
}
;}
e12
grammar_optimize_shortcut_logical_evaluation_type{y2
53
cS
grammar_optimize_shortcut_logical_evaluation_type
grammar_optimize_shortcut_logical_evaluation={53,{0,25,cQ
iO1
78,cR
cT}
lF3
l21{xW3
e42
e81
nT2
paramlist,nA1){index=(paramlist>>(index*10))&1023;if(index>=57)return
e42(SubFunction,cD2
plist_s[index-57]eY3
index>=37)return
e42(NumConstant,cD2
plist_n_container
xK::plist_n[index-37])n31
e42(tG3,cD2
plist_p[index]);}
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <stdio.h>
#include <algorithm>
#include <map>
#include <sstream>
using
lF3
FUNCTIONPARSERTYPES;using
lF3
l21;using
t5;using
yL1;lF3{nR1
It,xK3
T,xK3
Comp>eP1
MyEqualRange(It
first,It
last,const
T&val,Comp
comp){size_t
len=last-first;while(len>0){size_t
yA3
len/2;It
xG3(first);xG3+=half;if(comp(*xG3,val)){first=xG3;++first;len=len-half-1;}
lJ2
comp(val,*xG3)){len=half;}
else{It
left(first);{It&eF2=left;It
last2(xG3);size_t
len2=last2-eF2;while(len2>0){size_t
half2=len2/2;It
cD3(eF2);cD3+=half2;if(comp(*cD3,val)){eF2=cD3;++eF2;len2=len2-half2-1;}
else
len2=half2;}
}
first+=len;It
right(++xG3);{It&eF2=right;It&last2=first;size_t
len2=last2-eF2;while(len2>0){size_t
half2=len2/2;It
cD3(eF2);cD3+=half2;if(comp(val,*cD3))len2=half2;else{eF2=cD3;++eF2;len2=len2-half2-1;}
eW2
eP1(left,right);eW2
eP1(first,first);}
xW3
e12
OpcodeRuleCompare{iG2()lQ3
tree,iF2
yB2
e83{const
Rule&rule=grammar_rules[yB2]n31
i01<rule
cE2.subfunc_opcode;}
iG2()nT2
yB2,lC2
tree
e83{const
Rule&rule=grammar_rules[yB2]n31
rule
cE2.subfunc_opcode<i01;}
}
nJ1
bool
TestRuleAndApplyIfMatch
cH1
eT2
eJ&tree,bool
c1{MatchInfo
xK
info;lX1
found(false,cW()eY3(rule.lC1
LogicalContextOnly)&&!c1{tC1
if(nB
IsIntType
xK::nT3){if(rule.lC1
NotForIntegers)tC1
e63
rule.lC1
OnlyForIntegers)tC1
if(nB
IsComplexType
xK::nT3){if(rule.lC1
NotForComplex)tC1
e63
rule.lC1
OnlyForComplex)tC1
for(;;){
#ifdef DEBUG_SUBSTITUTIONS
#endif
found=TestParams(rule
cE2,tree,found.specs,info,true
eY3
found.found)break;if(!&*found.specs){fail:;
#ifdef DEBUG_SUBSTITUTIONS
DumpMatch(rule
iY3,false);
#endif
return
tR3}
#ifdef DEBUG_SUBSTITUTIONS
DumpMatch(rule
iY3,true);
#endif
SynthesizeRule(rule
iY3)xJ2}
yL1{xN1
ApplyGrammar
cH1
Grammar&tX2,eJ&tree,bool
c1{if(tree.GetOptimizedUsing()==&tX2){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Already optimized:  "
i13(tree
lT2"\n"
<<std::flush;
#endif
return
tR3
if(true){bool
changed
tO3
switch(i01
y33
cNot:case
cNotNot:case
cAnd:case
cOr:for
iA1
a=0;a<tree.x6
true))x62
lC
cIf:case
cAbsIf:if(ApplyGrammar(tX2,lE2
0),i01==cIf))x62
for
iA1
a=1;a<tree.x6
c1)x62
break;c73
for
iA1
a=0;a<tree.x6
false))x62}
if(changed){tree.Mark_Incompletely_Hashed()xJ2}
typedef
const
iF2
short*n83;std::pair<n83,n83>range=MyEqualRange(tX2.rule_list,tX2.rule_list+tX2.rule_count,tree,OpcodeRuleCompare
xK());std
xV3<iF2
short>rules;rules.xS3
range
tE3-range.first);for
y5
if(IsLogisticallyPlausibleParamsMatch(tD1
cE2,tree))rules.push_back(*r);}
range.first=!rules
cT3?&rules[0]:0;range
tE3=!rules
cT3?&rules[rules.tF3-1]+1:0;if(range.first!=range
iO2{
#ifdef DEBUG_SUBSTITUTIONS
if(range.first!=range
iO2
yC2"Input ("
<<c83
i01)<<")["
<<eM<<"]"
;if(c1
std::cout<<"(Logical)"
;iF2
first=iP1,prev=iP1;nH3
sep=", rules "
;for
y5
if(first==iP1)first=prev=*r;lJ2*r==prev+1)prev=*r;else
yC2
sep<<first;sep=","
;if(prev!=first)std::cout<<'-'<<prev;first=prev=*r;}
}
if(first!=iP1)yC2
sep<<first;if(prev!=first)std::cout<<'-'<<prev;}
std::cout<<": "
i13(tree
lT2"\n"
<<std::flush;}
#endif
bool
changed
tO3
for
y5
#ifndef DEBUG_SUBSTITUTIONS
if(!IsLogisticallyPlausibleParamsMatch(tD1
cE2,tree)cP2;
#endif
if(TestRuleAndApplyIfMatch(tD1,tree,c1){x62
yF3}
if(changed){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Changed."
<<std::endl;std::cout<<"Output: "
i13(tree
lT2"\n"
<<std::flush;
#endif
tree.Mark_Incompletely_Hashed()xJ2}
tree.SetOptimizedUsing(&tX2)n31
tR3
xN1
ApplyGrammar
cH1
void*p,FPoptimizer_CodeTree::e62
yR
ApplyGrammar(*cH1
Grammar*)p,tree);}
tZ1
ApplyGrammars(FPoptimizer_CodeTree::e62{
#ifdef DEBUG_SUBSTITUTIONS
std
iL3"grammar_optimize_round1\n"
;
#endif
n6
grammar_optimize_round1
n53
#ifdef DEBUG_SUBSTITUTIONS
std
iL3"grammar_optimize_round2\n"
;
#endif
n6
grammar_optimize_round2
n53
#ifdef DEBUG_SUBSTITUTIONS
std
iL3"grammar_optimize_round3\n"
;
#endif
n6
grammar_optimize_round3
n53
#ifndef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
#ifdef DEBUG_SUBSTITUTIONS
std
iL3"grammar_optimize_nonshortcut_logical_evaluation\n"
;
#endif
n6
grammar_optimize_nonshortcut_logical_evaluation
n53
#endif
#ifdef DEBUG_SUBSTITUTIONS
std
iL3"grammar_optimize_round4\n"
;
#endif
n6
grammar_optimize_round4
n53
#ifdef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
#ifdef DEBUG_SUBSTITUTIONS
std
iL3"grammar_optimize_shortcut_logical_evaluation\n"
;
#endif
n6
grammar_optimize_shortcut_logical_evaluation
n53
#endif
#ifdef FP_ENABLE_IGNORE_IF_SIDEEFFECTS
#ifdef DEBUG_SUBSTITUTIONS
std
iL3"grammar_optimize_ignore_if_sideeffects\n"
;
#endif
n6
grammar_optimize_ignore_if_sideeffects
n53
#endif
#ifdef DEBUG_SUBSTITUTIONS
std
iL3"grammar_optimize_abslogical\n"
;
#endif
n6
grammar_optimize_abslogical
n53
#undef C
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
#include <algorithm>
#include <assert.h>
#include <cstring>
#include <cmath>
#include <memory> /* for auto_ptr */
using
lF3
FUNCTIONPARSERTYPES;using
lF3
l21;using
t5;using
yL1;lF3{xN1
TestImmedConstraints
nT2
bitmask,lC2
tree
y13
bitmask&ValueMask
y33
Value_AnyNum:case
ValueMask:lC
nY2:if(GetEvennessInfo
eZ3
n52
Value_OddInt:if(GetEvennessInfo
eZ3
x02
tS1:if(GetIntegerInfo
eZ3
n52
Value_NonInteger:if(GetIntegerInfo
eZ3
x02
eJ1:if(!IsLogicalValue(tree)cL
nI1
SignMask
y33
Sign_AnySign:lC
nH1:if(l01
n52
eK1:if(l01
x02
Sign_NoIdea:if(l01
Unknown
cL
nI1
OnenessMask
y33
Oneness_Any:case
OnenessMask:lC
Oneness_One:if(!yN1
if(!cB3
fp_abs(eC3)tG1
1))cL
lC
Oneness_NotOne:if(!yN1
if(cB3
fp_abs(eC3)tG1
1))cL
nI1
ConstnessMask
y33
Constness_Any:lC
tR1:if(!yN1
lC
Constness_NotConst:if(yN1
yF3
return
true;}
i31
iF2
extent,iF2
nbits,xK3
eG2=iF2
int>e12
nbitmap{private:static
const
iF2
bits_in_char=8;static
const
iF2
eH2=(cM3
eG2)*bits_in_char)/nbits;eG2
data[(extent+eH2-1)/eH2];eX3
void
inc(nA1,int
by=1){data[pos(index)]+=by*eG2(1<<yD2);xI1
void
dec(nA1){inc(index,-1);}
int
get(nA1
yZ1(data[pos(index)]>>yD2)&mask()yH3
pos(nA1)yR
index/eH2
yH3
shift(nA1)yR
nbits*(index%eH2)yH3
mask()yR(1<<nbits)-1
yH3
mask(nA1)yR
mask()<<yD2;}
}
;e12
eJ3{int
SubTrees:8;int
eV3:8;int
yE2:8;int
Immeds:8;nbitmap<l03,2>SubTreesDetail;eJ3(){std::memset(this,0,cM3*this));}
eJ3
cH1
eJ3&b){std::memcpy(this,&b,cM3
b));}
eJ3&eA1=cH1
eJ3&b){std::memcpy(this,&b,cM3
b))n31*this;}
}
nJ1
eJ3
CreateNeedList_uncached(tC&c32){eJ3
cN1;for(eQ1=0;a<c32
yF2;++a){const
e42&yH2=e81
xK(c32.param_list,a);switch(yH2.first
y33
SubFunction:{cP
yK3(cP*c61
yF
GroupFunction)++cN1.Immeds;else{++cN1.SubTrees;assert(param.data.subfunc_opcode<VarBegin);cN1.SubTreesDetail.inc
y23
l81);}
++cN1.yE2;yF3
case
NumConstant:case
tG3:++cN1.eV3;++cN1.yE2;break;eW2
cN1;}
xW3
eJ3&CreateNeedList(tC&c32){typedef
std::map<tC*,eJ3>e51;static
e51
yV1;e51::yV3
i=yV1.xU2&c32
eY3
i!=yV1.cY1&c32)return
i
e82
n31
yV1.yR3,std::make_pair(&c32,CreateNeedList_uncached
xK(c32)))e82;}
xW3
eJ
CalculateGroupFunction
cH1
e42&yG2
const
tB
info
y13
yH2.first
y33
NumConstant:{const
ParamSpec_NumConstant
xK
yK3
cH1
ParamSpec_NumConstant
xK*)yH2
tE3
n31
CodeTreeImmed
y23
constvalue
i33
tG3:{cO
yK3(cO*)yH2
tE3
n31
tM3
GetParamHolderValueIfFound
y23
index
i33
SubFunction:{cP
yK3(cP*c61
eJ
nT3;nT3
xD
yL3
l81);nT3.l02).xS3
yL3
data
yF2);for(eQ1=0;a<yL3
data
yF2;++a){eJ
tmp(CalculateGroupFunction(e81
xK
y23
data.param_list,a),info));nT3
yP1
tmp);}
nT3.Rehash()n31
nT3;eW2
eJ();}
}
yL1{xN1
IsLogisticallyPlausibleParamsMatch(tC&c32,lC2
tree){eJ3
cN1(CreateNeedList
xK(c32));size_t
tZ3=c8
if(tZ3<size_t(cN1.yE2))n93}
for
iA1
a=0;a<tZ3;++a){iF2
opcode=t33
nC;switch(opcode
y33
cImmed:if(cN1.Immeds>0)eU3
Immeds;else
eU3
eV3;lC
l03:case
cFCall:case
cPCall:eU3
eV3;break;c73
assert(opcode<VarBegin);if(cN1.SubTrees>0&&cN1.SubTreesDetail.get(opcode)>0){eU3
SubTrees;cN1.SubTreesDetail.dec(opcode
t43
eU3
eV3;}
}
if(cN1.Immeds>0||cN1.SubTrees>0||cN1.eV3>0)n93}
if(c32.match_type!=AnyParams){if(0||cN1.SubTrees<0||cN1.eV3<0)n93
eW2
true;}
xW3
lX1
TestParam
cH1
e42&yG2
lC2
tree,const
cW&start_at
tJ3
y13
yH2.first
y33
NumConstant:{const
ParamSpec_NumConstant
xK
yK3
cH1
ParamSpec_NumConstant
xK*c61
if(!yN1
iV2
imm=eC3;switch
y23
modulo
y33
Modulo_None:lC
Modulo_Radians:imm=e93
imm,yO
imm<xD1
imm
yP
if(imm>fp_const_pi
xK())imm-=fp_const_twopi
xK(cT2
return
cB3
imm,yL3
constvalue
i33
tG3:{cO
yK3(cO*c61
if(!x0
return
tM3
SaveOrTestParamHolder
y23
index,tree
i33
SubFunction:{cP
yK3(cP*c61
yF
GroupFunction){if(!x0
eJ
xY1=CalculateGroupFunction(yG2
info);
#ifdef DEBUG_SUBSTITUTIONS
DumpHashes(xY1
lT2*cH1
void**)&xY1
nC1;std::cout<<"\n"
;std::cout<<*cH1
void**)&eC3;std::cout<<"\n"
;DumpHashes(tree
lT2"Comparing "
i13(xY1
lT2" and "
i13(tree
lT2": "
;std::cout<<(xY1
xF
tree)?"true"
:"false"
lT2"\n"
;
#endif
return
xY1
xF
tree);}
e63!&*start_at){if(!x0
if(i01!=yL3
l81
cL}
return
TestParams
y23
data,tree,start_at,info,false);}
eW2
tR3
xW3
e12
iX
x72
MatchInfo
xK
info;iX():start_at(),info(){}
}
nJ1
class
MatchPositionSpec_PositionalParams:xE1
iX
xK>{eX3
l53
MatchPositionSpec_PositionalParams
iA1
n
yB1
iX
xK>(n){}
}
;e12
iQ1
x72
iQ1():start_at(){}
}
;class
yG:xE1
iQ1>{eX3
iF2
trypos;l53
yG
iA1
n
yB1
iQ1>(n),trypos(0){}
}
nJ1
lX1
TestParam_AnyWhere
cH1
e42&yG2
lC2
tree,const
cW&start_at
tJ3,std::eV2&used,bool
eY2{xR<yG>x5;tY2
yG
eU2
a=x5->trypos;goto
retry_anywhere_2
tC2
yG(eM);a=0;}
for(;a<c8++a){if(used[a]cP2;retry_anywhere
i03=TestParam(yG2
t33,(i02);iB1
used[a]=true
iP
a);x5->trypos=a
n31
lX1(true,&*x5);}
}
retry_anywhere_2
i23
goto
retry_anywhere;eW2
tR3
xW3
e12
yC1
x72
MatchInfo
xK
info;std::eV2
used;l53
yC1
iA1
tZ3):start_at(),info(),used(tZ3){}
}
nJ1
class
MatchPositionSpec_AnyParams:xE1
yC1
xK>{eX3
l53
MatchPositionSpec_AnyParams
iA1
n,size_t
m
yB1
yC1
xK>(n,yC1
xK(m)){}
}
nJ1
lX1
TestParams(tC&nN,lC2
tree,const
cW&start_at
tJ3,bool
eY2{if(nN.match_type!=AnyParams){if(xV!=eM
cL}
if(!IsLogisticallyPlausibleParamsMatch(nN,tree))n93
l73
nN.match_type
y33
PositionalParams:{xR<cG>x5;tY2
cG
eU2
a=xV-1;goto
lD1
tC2
cG(xV);a=0;}
for(;a
eX2{lG2=info;retry_positionalparams
i03=TestParam(cU
a),t33,(i02);iB1
eR2}
lD1
i23
info=lG2;goto
retry_positionalparams;}
if(a>0){--a;goto
lD1;}
info=eD3
n31
tR3
if(eY2
for(tZ2
tM3
SaveMatchedParamIndex(a)n31
lX1(true,&*x5
i33
SelectedParams:case
AnyParams:{xR<t7>x5;std::eV2
used(eM);std
xV3<iF2>l63(xV);std
xV3<iF2>yI2(xV);for(tZ2{const
e42
yH2=cU
a);l63[a]=ParamSpec_GetDepCode(yH2);}
{iF2
b=0;for(tZ2
if(l63[a]!=0)yI2[b++]=a;for(tZ2
if(l63[a]==0)yI2[b++]=a;}
tY2
t7
eU2
if(xV==0){a=0;goto
retry_anyparams_4;}
a=xV-1;goto
e61
tC2
t7(xV,eM);a=0;if(xV!=0){eD3=info;(*x5)[0].used=used;}
}
for(;a
eX2{if(a>0){lG2=info;(*x5)[a].used=used;}
retry_anyparams
i03=TestParam_AnyWhere
xK(cU
yI2[a]),tree,(i02,used,eY2;iB1
eR2}
e61
i23
info=lG2;used=(*x5)[a].used;goto
retry_anyparams;}
e71:if(a>0){--a;goto
e61;}
info=eD3
n31
tR3
retry_anyparams_4:if(nN.n1!=0){if(!TopLevel||!tM3
HasRestHolder(nN.n1)){eH
cF2;cF2.xS3
eM);for
nT2
b=0;b<c8++b){if(cN3
cP2;cF2.push_back
y91
b));cN3=true
iP
b);}
if(!tM3
SaveOrTestRestHolder(nN.n1,cF2)){goto
e71;}
}
else{iH2&cF2=tM3
GetRestHolderValues(nN.n1)yM3
0;a<cF2.tF3;++a){bool
found
tO3
for
nT2
b=0;b<c8++b){if(cN3
cP2;if(cF2[a]xF
lE2
b))){cN3=true
iP
b);found=true;yF3}
if(!found){goto
e71;}
}
eW2
lX1(true,xV?&*x5:0
i33
GroupFunction:yF3
return
tR3}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
#include <algorithm>
#include <assert.h>
using
t5;using
yL1;lF3{xW3
eJ
xZ1
const
e42&yG2
tB
info,bool
inner=true
y13
yH2.first
y33
NumConstant:{const
ParamSpec_NumConstant
xK
yK3
cH1
ParamSpec_NumConstant
xK*)yH2
tE3
n31
CodeTreeImmed
y23
constvalue
i33
tG3:{cO
yK3(cO*)yH2
tE3
n31
tM3
GetParamHolderValue
y23
index
i33
SubFunction:{cP
yK3(cP*c61
eJ
tree;tI3
yL3
l81);for(eQ1=0;a<yL3
data
yF2;++a){eJ
nparam=xZ1
e81
xK
y23
data.param_list,a),info,true);lP1
nparam);}
if
y23
data.n1!=0){eH
trees(tM3
GetRestHolderValues
y23
data.n1));tree.AddParamsMove(trees
eY3
eM==1){assert(tree.GetOpcode()==cAdd||tree.GetOpcode()==cMul||tree.GetOpcode()==cMin||tree.GetOpcode()==cMax||tree.GetOpcode()==cAnd||tree.GetOpcode()==cOr||tree.GetOpcode()==cAbsAnd||tree.GetOpcode()==cAbsOr);tree.xY2
0));}
lJ2
eM==0
y13
i01
y33
cAdd:case
cOr:tree=nF1
0));lC
cMul:case
cAnd:tree=nF1
1));c73
yF3}
}
if(inner)tree.Rehash()n31
tree;eW2
eJ();}
}
yL1{tZ1
SynthesizeRule
cH1
eT2
eJ&tree
tJ3
y13
rule.ruletype
y33
ProduceNewTree:{tree.Become(xZ1
e81
lE1
0),info,false)cT2
case
ReplaceParams:c73{std
xV3<iF2>list=tM3
GetMatchedParamIndexes();std::sort(list.iN2
list.end())yM3
list.tF3;a-->0;)tree
tK2
list[a]);for(eQ1=0;a<rule.repl_param_count;++a){eJ
nparam=xZ1
e81
lE1
a),info,true);lP1
nparam);}
yF3}
}
}
#endif
#ifdef DEBUG_SUBSTITUTIONS
#include <sstream>
#include <cstring>
using
lF3
FUNCTIONPARSERTYPES;using
lF3
l21;using
t5;using
yL1;lF3
l21{tZ1
DumpMatch
cH1
eT2
lC2
tree,const
tB
info,bool
DidMatch,std::ostream&o){DumpMatch(rule
iY3,DidMatch?l24"match"
:l24"mismatch"
,o);}
tZ1
DumpMatch
cH1
eT2
lC2
tree,const
tB
info,nH3
i43,std::ostream&o){static
const
char
ParamHolderNames[][2]={"%"
,"&"
,"x"
,"y"
,"z"
,"a"
,"b"
,"c"
}
;o<<i43<<" (rule "
<<(&rule-grammar_rules)<<")"
<<":\n  Pattern    : "
;{e42
tmp;tmp.first=SubFunction;ParamSpec_SubFunction
tmp2;tmp2.data=rule
cE2;tmp
tE3=cD2
tmp2;DumpParam
xK(tmp,o);}
o<<"\n  Replacement: "
;DumpParams
lE1
rule.repl_param_count,o);o<<"\n"
;o<<"  Tree       : "
i13(tree,o);o<<"\n"
;if(!std::strcmp(i43,l24"match"
))DumpHashes(tree,o)yM3
0;a<tM3
cE
tF3;++a){if(!tM3
paramholder_matches[a].IsDefined()cP2;o<<"           "
<<ParamHolderNames[a]<<" = "
i13(tM3
paramholder_matches[a],o);o<<"\n"
;}
i81
tM3
lQ.tF3;++b){if(!tM3
lQ[b
tC3
cP2
yM3
0;a<tM3
lQ[b
tD3.tF3;++a){o<<"         <"
<<b<<"> = "
i13(tM3
lQ[b
tD3[a],o);o<<std::endl;}
}
o<<std::flush;}
}
#endif
#include <list>
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
using
lF3
FUNCTIONPARSERTYPES;lF3{xN1
MarkIncompletes(FPoptimizer_CodeTree::e62{if(tree.Is_Incompletely_Hashed(tI1;bool
iR1
tO3
nB1
iR1|=MarkIncompletes
y91
a)eY3
iR1)tree.Mark_Incompletely_Hashed()n31
iR1;}
tZ1
FixIncompletes(FPoptimizer_CodeTree::e62{if(tree.Is_Incompletely_Hashed()){nB1
FixIncompletes
y91
a));tree
eT1}
}
}
t5{lB
Sort()tN3
Sort();}
lB
Rehash(bool
constantfolding){if(constantfolding)ConstantFolding(*this);else
Sort();data
x7
xW3
e12
cA{cO3
iV2
cP3
yD1=0;
#if 0
long
double
value=Value;e5=crc32::calc(cH1
iF2
char*)&value,cM3
value));key^=(key<<24);
#elif 0
union{e12{iF2
char
filler1[16]iU2
v;iF2
char
filler2[16];}
buf2;e12{iF2
char
filler3[cM3
iV2)+16-c
M3
x91)];e5;}
buf1;}
data;memset(&data,0,cM3
data));data.buf2.v=Value;e5=data.buf1.key;
#else
int
exponent
iU2
nM2=std::frexp(Value,&xG2;e5=nT2
t62+0x8000)&0xFFFF
eY3
nM2<0){nM2=-nM2;key=key^0xFFFF;}
else
key+=0x10000;nM2-=iV2(0.5);key<<=39;key|=lY1(nM2+nM2)*iV2(1u<<31))<<8;
#endif
lP
#ifdef FP_SUPPORT_COMPLEX_NUMBERS
nR1
T
cG2
std::complex<T> >{cO3
std::complex<T>cP3
cA<T>::nA3
cH2,Value.real());nB
fphash_t
temp;cA<T>::nA3
temp,Value.imag());yD1^=temp.hash2;cH2.hash2^=temp.hash1;}
}
;
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
i31
cG2
long>{yH
long
Value){e5=Value;lP
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
i31
cG2
GmpInt>{cO3
GmpInt
cP3
e5=Value.toInt();lP
#endif
tZ1
xH2
xK::Recalculate_Hash_NoRecursion(){fphash_t
cH2(lY1
Opcode)<<56,Opcode*iP3(0x1131462E270012B));Depth=1;switch(Opcode
y33
cImmed:{cA
xK::nA3
cH2,Value
cT2
case
l03:{yD1|=lY1
cQ1<<48
lF1((lY1
cQ1)*11)^iP3(0x3A83A83A83A83A0);yF3
case
cFCall:case
cPCall:{yD1|=lY1
cQ1<<48
lF1((~lY1
cQ1)*7)^3456789;}
c73{size_t
eR1=0
yM3
0;a<xX3;++a){if(iZ2
y22>eR1)eR1=iZ2
y22;yD1+=((iZ2
i12
hash1*(a+1))>>12)lF1
iZ2
i12
hash1
lF1(3)*iP3(0x9ABCD801357);cH2.hash2*=iP3(0xECADB912345)lF1(~iZ2
i12
hash2)^4567890;}
Depth+=eR1;}
}
if(Hash!=cH2){Hash=cH2;l32=0;}
}
lB
FixIncompleteHashes(){MarkIncompletes(*this);FixIncompletes(*this);}
}
#endif
#include <cmath>
#include <list>
#include <cassert>
#ifdef FP_SUPPORT_OPTIMIZER
using
lF3
FUNCTIONPARSERTYPES;lF3{using
t5
nJ1
bool
x01
lC2
tree,long
count,const
xH1::SequenceOpCode
xK&eU,xH1
yG3&synth,size_t
max_bytecode_grow_length);static
const
e12
SinCosTanDataType{OPCODE
whichopcode;OPCODE
inverse_opcode;enum{nominator,nN2,inverse_nominator,lG1}
;OPCODE
codes[4];}
SinCosTanData[12]={{cTan,cCot,{cSin,cCos,cCsc,cSec}
}
,{cCot,cCot,{cCos,cSin,cSec,cCsc}
}
,{cCos,cSec,{cSin,cTan,cCsc,cCot}
}
,{cSec,cCos,{cTan,cSin,cCot,cCsc}
}
,{cSin,cCsc,{cCos,cCot,cSec,cTan}
}
,{cCsc,cSin,{cCot,cCos,cTan,cSec}
}
,{cI2{cSinh,cCosh,cJ2,{cSinh,cNop,{cI2
cNop,cCosh}
}
,{cCosh,cNop,{cSinh,cI2
cNop}
}
,{cNop,cTanh,{cCosh,cSinh,cJ2,{cNop,cSinh,{cNop,cTanh,cCosh,cNop}
}
,{cNop,cCosh,{cTanh,cSinh,cJ2}
;}
t5{lB
SynthesizeByteCode(std
xV3<iF2>&ByteCode,std
xV3
xK&Immed,size_t&stacktop_max){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Making bytecode for:\n"
;iR
#endif
while(RecreateInversionsAndNegations()){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"One change issued, produced:\n"
;iR
#endif
FixIncompleteHashes();using
yL1;using
lF3
l21;const
void*g=cD2
grammar_optimize_recreate;while(ApplyGrammar(*cH1
Grammar*)g,*this)){FixIncompleteHashes();}
}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Actually synthesizing, after recreating inv/neg:\n"
;iR
#endif
xH1
yG3
synth;SynthesizeByteCode(synth,false);synth.Pull(ByteCode,Immed,stacktop_max);}
lB
SynthesizeByteCode(xH1
yG3&synth,bool
MustPopTemps
e83{y01*this))yR;}
for
iA1
a=0;a<12;++a){const
SinCosTanDataType&data=SinCosTanData[a];if(data.whichopcode!=cNop)iD1!=data.whichopcode
cP2;CodeTree
nB3;nB3.lN1
nB3
cQ3
inverse_opcode);nB3.n62
y01
nB3)){synth.l11
else
iD1!=cInv
cP2;if(GetParam(0)nC!=data.inverse_opcode
cP2;y01
GetParam(0))){synth.l11
size_t
found[4];i81
4;++b){CodeTree
tree;if(data.i53]==cNop){tI3
cInv);CodeTree
nC3;nC3.lN1
nC3
cQ3
i53^2]);nC3.n62
lP1
nC3
t43{tree.lN1
tree
cQ3
i53]);}
tree.n62
found[b]=synth.xU3
tree);}
if(found[data.yJ2!=tL
nN2]y4
yJ2);n81
nN2
i5
cDiv
nK1
yJ2!=tL
lG1]y4
yJ2);n81
lG1
i5
cMul
nK1
lR1!=tL
lG1]y4
lR1);n81
lG1
i5
cRDiv
nK1
lR1!=tL
nN2]y4
lR1);n81
nN2
i5
cMul,2,1);synth.l11
size_t
n_subexpressions_synthesized=SynthCommonSubExpressions(iS2
switch(iS1{case
l03:synth.PushVar(GetVar());lC
cImmed:lH2
GetImmed());lC
cAdd:case
cMul:case
cMin:case
cMax:case
cAnd:case
cOr:case
cAbsAnd:case
cAbsOr:iD1==cMul){bool
cR3
tO3
yI
lS1
i61&&isLongInteger(lS1
nC1)){c82=makeLongInteger(lS1
nC1);CodeTree
tmp(*this,xK3
CodeTree::CloneTag());tmp
t31);tmp
eT1
if(x01
tmp,value,xH1::tN1
xK::AddSequence,synth,MAX_MULI_BYTECODE_LENGTH)){cR3=true;yF3}
}
if(cR3)yF3
int
yE1=0;std::eV2
done(GetParamCount(),false);CodeTree
iF;iF
xD
iS1;for(;;){bool
found
tO3
yI
done[a]cP2;if(synth.IsStackTop(lS1)){found=true;done[a]=true;lS1.n7
iF
c9
lS1
eY3++yE1>1){synth
yJ
2);iF.Rehash(false)cI1
iF);yE1=yE1-2+1;}
}
}
if(!found)yF3
yI
done[a]cP2;lS1.n7
iF
c9
lS1
eY3++yE1>1){synth
yJ
2);iF.Rehash(false)cI1
iF);yE1=yE1-2+1;}
}
if(yE1==0
y13
iS1{case
cAdd:case
cOr:case
cAbsOr:lH2
0);lC
cMul:case
cAnd:case
cAbsAnd:lH2
1);lC
cMin:case
cMax:lH2
0);break;c73
yF3++yE1;}
assert(n_stacked==1);yF3
case
cPow:{lG3
p0=GetParam(0);lG3
p1=GetParam(1
eY3!p1
i61||!isLongInteger
eH3)||!x01
p0,makeLongInteger
eH3),xH1::tN1
xK::MulSequence,synth,MAX_POWI_BYTECODE_LENGTH)){p0.n7
p1.n7
synth
yJ
2);i91
cIf:case
cAbsIf:{xK3
xH1
yG3::IfData
ifdata;GetParam(0).n7
synth.SynthIfStep1(ifdata,iS1;GetParam(1).n7
synth.SynthIfStep2(ifdata);GetParam(2).n7
synth.SynthIfStep3(ifdata
cT2
case
cFCall:case
cPCall:{for
iA1
a=0;a<xF1++a)lS1.n7
synth
yJ
nT2)GetParamCount());yE3
nK2|GetFuncNo(),0,0
cT2
c73{for
iA1
a=0;a<xF1++a)lS1.n7
synth
yJ
nT2)GetParamCount()cT2}
synth.StackTopIs(*this
eY3
MustPopTemps&&n_subexpressions_synthesized>0){size_t
top=synth.GetStackTop();synth.DoPopNMov(top-1-n_subexpressions_synthesized,top-1);}
}
}
lF3{xN1
x01
lC2
tree,long
count,const
xH1::SequenceOpCode
xK&eU,xH1
yG3&synth,size_t
max_bytecode_grow_length){if
eL3!=0){xH1
yG3
backup=synth;tree.n7
size_t
bytecodesize_backup=synth.GetByteCodeSize();xH1::x01
count,eU,iS2
size_t
bytecode_grow_amount=synth.GetByteCodeSize()-bytecodesize_backup;if(bytecode_grow_amount>max_bytecode_grow_length){synth=backup
n31
tR3
return
true;}
else{xH1::x01
count,eU,synth)xJ2}
}
#endif
#include <cmath>
#include <cassert>
#ifdef FP_SUPPORT_OPTIMIZER
using
lF3
FUNCTIONPARSERTYPES;lF3{using
t5;
#define FactorStack std xV3
const
e12
PowiMuliType{iF2
opcode_square;iF2
opcode_cumulate;iF2
opcode_invert;iF2
opcode_half;iF2
opcode_invhalf;}
iseq_powi={cSqr,cMul,cInv,cSqrt,cRSqrt}
,iseq_muli={iP1
xE
cNeg,iP1,iP1}
nJ1
iV2
eS1
cH1
PowiMuliType&cS3,const
std
xV3<iF2>&xP1,n72&stack
tO2
1);while(IP<limit){if(iT1==cS3.opcode_square){if(!isInteger
xO1
2;yW
opcode_invert){nT3=-nT3;yW
opcode_half){if
iT2>y21&&isEvenInteger
xO1
iV2(0.5);yW
opcode_invhalf){if
iT2>y21&&isEvenInteger
xO1
iV2(-0.5);++IP;eR2
size_t
x12=IP
iU2
lhs(1
eY3
iT1==cFetch){nA1=nD2;if(index<y1||size_t(index-y1)>=iV1){IP=x12;yF3
lhs=stack[index-y1];goto
yK2;}
if(iT1==cDup){lhs=nT3;goto
yK2;yK2:eX1
nT3);++IP
iU2
subexponent=eS1(cS3
tS
if(IP>=limit||iT1!=cS3.opcode_cumulate){IP=x12;yF3++IP;stack.pop_back();nT3+=lhs*subexponent;eR2
yF3
return
nT3;}
xW3
iV2
ParsePowiSequence
cH1
std
xV3<iF2>&xP1){n72
stack;eX1
iV2(1))n31
eS1(iseq_powi
tS}
xW3
iV2
ParseMuliSequence
cH1
std
xV3<iF2>&xP1){n72
stack;eX1
iV2(1))n31
eS1(iseq_muli
tS}
xW3
class
CodeTreeParserData{eX3
l53
CodeTreeParserData(bool
k_powi):stack(),clones(),keep_powi(k_powi){}
void
Eat
iA1
tZ3,OPCODE
opcode){eJ
xI;xI
xD
opcode);eH
c32=Pop(tZ3)tQ1
c32
eY3!keep_powi)switch(opcode
y33
cTanh:{eJ
sinh,cosh;sinh
xD
cSinh);sinh
c9
xI
lR2
sinh
eT1
cosh
xD
cCosh);cosh
yP1
xI
lR2
cosh
eT1
eJ
pow;pow
xD
cPow);pow
yP1
cosh);pow.yC
iV2(-1)));pow
eT1
xI
tH3
xI.nG1
0,sinh);xI
yP1
pow
cT2
case
cTan:{eJ
sin,cos;sin
xD
cSin);sin
c9
xI
lR2
sin
eT1
cos
xD
cCos);cos
yP1
xI
lR2
cos
eT1
eJ
pow;pow
xD
cPow);pow
yP1
cos);pow.yC
iV2(-1)));pow
eT1
xI
tH3
xI.nG1
0,sin);xI
yP1
pow
cT2
case
cPow:{lC2
p0=xI
l8
0);lC2
p1=xI
l8
1
eY3
p1
nC==cAdd){eH
mulgroup(p1.GetParamCount())yM3
0;a<p1.xF1++a){eJ
pow;pow
xD
cPow);pow
c9
p0);pow
c9
p1
x13
pow
eT1
mulgroup[a
i63
pow);}
xI
xD
cMul)tQ1
nN3;}
yF3
c73
yF3
xI.Rehash(!keep_powi);iU1,false);
#ifdef DEBUG_SUBSTITUTIONS
eZ1<<tZ3<<", "
<<c83
opcode)<<"->"
<<c83
xI
nC)<<": "
iV3
xI
tT
xI);
#endif
eX1
xI
lI3
EatFunc
iA1
tZ3,OPCODE
t83
iF2
funcno){eJ
xI=CodeTreeFuncOp
xK(t83
funcno);eH
c32=Pop(tZ3)tQ1
c32);xI.n62
#ifdef DEBUG_SUBSTITUTIONS
eZ1<<tZ3<<", "
iV3
xI
tT
xI);
#endif
iU1);eX1
xI
lI3
AddConst(yJ1){eJ
xI=CodeTreeImmed(value);iU1);Push(xI
lI3
AddVar
nT2
varno){eJ
xI=CodeTreeVar
xK(varno);iU1);Push(xI
lI3
SwapLastTwoInStack(){eY1
1
i63
eY1
2]lI3
Dup(){Fetch(iV1-1
lI3
Fetch
iA1
which){Push(stack[which]);}
nR1
T>void
Push(T
tree){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<iV3
tree
tT
tree);
#endif
eX1
tree
lI3
PopNMov
iA1
target,size_t
source){stack[target]=stack[source];stack.t73
target+1);}
eJ
yL2{clones.clear();eJ
nT3(stack.back());stack.t73
iV1-1)n31
nT3;}
eH
Pop
iA1
n_pop){eH
nT3(n_pop);for
nT2
n=0;n<n_pop;++n)nT3[n
i63
eY1
n_pop+n]);
#ifdef DEBUG_SUBSTITUTIONS
for
iA1
n=n_pop;n-->0;){eZ1
i13
iT2[n]tT
nT3[n]);}
#endif
stack.t73
iV1-n_pop)n31
nT3;}
size_t
GetStackTop(yZ1
iV1;}
private:void
FindClone(eJ&,bool=true)yR;}
private:eH
stack;std::multimap<fphash_t,eJ>clones;bool
keep_powi;private:CodeTreeParserData
cH1
CodeTreeParserData&);CodeTreeParserData&eA1=cH1
CodeTreeParserData&);}
nJ1
e12
IfInfo{eJ
eI2;eJ
thenbranch;size_t
endif_location;IfInfo():eI2(),thenbranch(),endif_location(){}
}
;}
t5{lB
GenerateFrom
cH1
xK3
FunctionParserBase
xK::Data&xH3,bool
keep_powi){eH
xC2;xC2.xS3
xH3.mVariablesAmount);for
nT2
n=0;n<xH3.mVariablesAmount;++n){xC2.push_back(CodeTreeVar
xK(n+l03));}
GenerateFrom(xH3,xC2,keep_powi);}
lB
GenerateFrom
cH1
xK3
FunctionParserBase
xK::Data&xH3,const
x1&xC2,bool
keep_powi){const
std
xV3<iF2>&ByteCode=xH3.mByteCode;const
std
xV3
xK&Immed=xH3.mImmed;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"ENTERS GenerateFrom()\n"
;
#endif
CodeTreeParserData
xK
sim(keep_powi);std
xV3<IfInfo
xK>eG;for
iA1
IP=0,DP=0;;++IP){i22:while(!eG
cT3&&(eG.e6==IP||(IP<i8
tF3&&iT1==cJump&&eG.tE1.IsDefined()))){CodeTree
elsebranch=sim.yL2
nB2
eG.back().eI2)nB2
eG.tE1)nB2
elsebranch);tQ2
3,cIf);eG.pop_back();}
if(IP>=i8
tF3)break;iF2
opcode=iT1;if((opcode==cSqr||opcode==cDup||(opcode==cInv&&!IsIntType
xK::nT3)||opcode==cNeg||opcode==cSqrt||opcode==cRSqrt||opcode==cFetch)){size_t
was_ip=IP
iU2
e92
ParsePowiSequence
xK(ByteCode,IP,eG
cT3?i8
tF3:eG.e6,sim.xM
1);if
t62!=iV2(1.0)){xL
exponent
xD2;goto
i22;}
if(opcode==cDup||opcode==cFetch||opcode==cNeg
tJ2
xT2=ParseMuliSequence
xK(ByteCode,IP,eG
cT3?i8
tF3:eG.e6,sim.xM
1
eY3
xT2!=iV2(1.0)){xL
xT2)y3
cMul);goto
i22;}
}
IP=was_ip;}
if(n82>=l03){nA1=opcode-l03
nB2
xC2[index]t43{switch(n82
y33
cIf:case
cAbsIf:{eG.t73
eG.tF3+1);CodeTree
res(sim.yL2);eG.back().eI2.swap(res);eG.e6=i8
tF3;IP+=2;eR2
case
cJump:{CodeTree
res(sim.yL2);eG.tE1.swap(res);eG.e6=lE3
IP+1]+1;IP+=2;eR2
case
cImmed:xL
Immed[DP++]);lC
cDup:sim.Dup();lC
cNop:lC
cFCall:{iF2
funcno=nD2;assert(funcno<fpdata.mFuncPtrs.size());iF2
c32=xH3.mFuncPtrs[funcno].mParams;sim.EatFunc(c32,n82,funcno
cT2
case
cPCall:{iF2
funcno=nD2;assert(funcno<fpdata.iT3.size());const
FunctionParserBase
xK&p=*xH3.iT3[funcno].mParserPtr;iF2
c32=xH3.iT3[funcno].mParams;x1
paramlist=sim.Pop(c32);CodeTree
i32;i32.GenerateFrom(*p.mData,paramlist)nB2
i32
cT2
case
cInv:xL
1
cQ2
cDiv);lC
cNeg:nC2
cNeg);break;xL
0
cQ2
cSub);lC
cSqr:xL
2
xD2;lC
cSqrt
e72
t01
cRSqrt
e72-t01
cCbrt
e72
1)/iV2(3)xD2;lC
cDeg:xL
fp_const_rad_to_deg
eU1
cRad:xL
fp_const_deg_to_rad
eU1
cExp:cR1)goto
iQ2;xL
fp_const_e
xK()cQ2
cPow);lC
cExp2:cR1)goto
iQ2;xL
2.0
cQ2
cPow);lC
cCot:nC2
cTan
nW
cCsc:nC2
cSin
nW
cSec:nC2
cCos
nW
cInt:
#ifndef __x86_64
cR1){nC2
cInt
cT2
#endif
xL
iV2(0.5))nD3
nC2
cFloor);lC
cLog10:nC2
cX3
fp_const_log10inv
eU1
cLog2:nC2
cX3
fp_const_log2inv
eU1
eF3:iZ3
cX3
fp_const_log2inv
xK());tQ2
3,cMul);lC
cHypot:xL
2
xD2;i73
xL
2
xD2
nD3
xL
iV2(t01
cSinCos:sim.Dup();nC2
cSin);iZ3
cCos);lC
cSinhCosh:sim.Dup();nC2
cSinh);iZ3
cCosh);lC
cRSub:i73
case
cSub:cR1){tQ2
2,cSub
cT2
xL-1)y3
cMul)nD3
lC
cRDiv:i73
case
cDiv:cR1||IsIntType
xK::nT3){tQ2
2,cDiv
cT2
xL-1
xD2
y3
cMul);lC
cAdd:case
cMul:case
cMod:case
cPow:case
cEqual:case
cLess:case
cGreater:case
cNEqual:case
cLessOrEq:case
cGreaterOrEq:case
cAnd:case
cOr:case
cAbsAnd:case
cAbsOr:tQ2
2,t11
lC
cNot:case
cNotNot:case
cZ3:case
cAbsNotNot:nC2
t11
lC
cFetch:sim.Fetch(nD2);lC
cPopNMov:{iF2
stackOffs_target=nD2;iF2
stackOffs_source=nD2;sim.PopNMov(stackOffs_target,stackOffs_source
cT2
c73
iQ2:;iF2
funcno=opcode-cAbs;assert(funcno<FUNC_AMOUNT);const
FuncDefinition&func=Functions[funcno];tQ2
func.c32,t11
yF3}
}
Become(sim.yL2);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Produced tree:\n"
;iR
#endif
}
}
#endif
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
#include <assert.h>
#define FP_MUL_COMBINE_EXPONENTS
lF3{using
lF3
FUNCTIONPARSERTYPES;using
t5
nJ1
static
void
AdoptChildrenWithSameOpcode(e62{
#ifdef DEBUG_SUBSTITUTIONS
bool
nO2
tO3
#endif
for
iY
if
y91
a)nC==i01){
#ifdef DEBUG_SUBSTITUTIONS
if(!nO2)yC2"Before assimilation: "
yA1
nO2=true;}
#endif
tree.AddParamsMove
y91
a).GetUniqueRef().l02),a);}
#ifdef DEBUG_SUBSTITUTIONS
if(nO2)yC2"After assimilation:   "
yA1}
#endif
}
}
t5{tZ1
ConstantFolding(e62{tree.Sort();
#ifdef DEBUG_SUBSTITUTIONS
void*cY3=0;std::cout<<"["
<<(&cY3)<<"]Runs ConstantFolding for: "
yA1
DumpHashes(tree
lT2
std::flush;
#endif
if(false){redo:;tree.Sort();
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"["
<<(&cY3)<<"]Re-runs ConstantFolding: "
yA1
DumpHashes(tree);
#endif
}
if(i01!=cImmed){c23
p=iN
tree
eY3
p
n51&&p
nP2
tF2==e13){xO
p
tF2);nD}
if(false){ReplaceTreeWithOne:xO
iV2(1));goto
do_return;ReplaceTreeWithZero:xO
xD1;goto
do_return;ReplaceTreeWithParam0:
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before replace: "
;std::cout<<std::hex<<'['lD2
hash1<<','lD2
hash2<<']'<<std::dec
yA1
#endif
tree.xY2
0));
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After replace: "
;std::cout<<std::hex<<'['lD2
hash1<<','lD2
hash2<<']'<<std::dec
yA1
#endif
cZ
l73
i01
y33
cImmed:lC
l03:lC
cAnd:case
cAbsAnd:eQ
bool
c2
tO3
for
iY{if(!yX3
a)))c2=true;cW3
a),i01==cAbsAnd)y33
IsNever
eL
IsAlways:tree
t31);lC
lT1
l73
eM
y33
0:i6
1:tI3
i01==cAnd?cNotNot:cAbsNotNot);cZ
c73
l83
cAnd||!c2)if(ConstantFolding_AndLogic(eJ2
i91
cOr:case
cAbsOr:eQ
bool
c2
tO3
for
iY{if(!yX3
a)))c2=true;cW3
a),i01==cAbsOr))l92
i6
lC3
tree
t31);lC
lT1
l73
eM
y33
0
eL
1:tI3
i01==cOr?cNotNot:cAbsNotNot);cZ
c73
l83
cOr||!c2)if(ConstantFolding_OrLogic(eJ2
i91
cNot:case
cZ3:{iF2
lZ1
0;switch
y91
0)nC
y33
cEqual:lZ1
cNEqual;lC
cNEqual:lZ1
cEqual;lC
cLess:lZ1
cGreaterOrEq;lC
cGreater:lZ1
cLessOrEq;lC
cLessOrEq:lZ1
cGreater;lC
cGreaterOrEq:lZ1
cLess;lC
cNotNot:lZ1
cNot;lC
cNot:lZ1
cNotNot;lC
cZ3:lZ1
cAbsNotNot;lC
cAbsNotNot:lZ1
cZ3;break;c73
yF3
if(opposite){tI3
OPCODE(opposite));tree.SetParamsMove
y91
0).GetUniqueRef().l02));cZ
l73
tM
0),tree
cS1
y33
IsAlways
eL
lC3
i6
lT1
l83
cNot&&GetPositivityInfo
y91
0))==IsAlways)tI3
cZ3);lT3
nC==cIf||lE2
0)nC==yW3{eJ
l93=lE2
0);lC2
ifp1=l93
l8
1);lC2
ifp2=l93
l8
2
eY3
ifp1
nC
lA3
ifp1
cS1{tree
xJ
ifp1
nC==cNot?cNotNot:cAbsNotNot
cU3
lR2
cV3
eJ
p2;p2
iP2
p2
e03
tN
if(ifp2
nC
lA3
ifp2
cS1{tree
xJ
i01
cU3);cV3
eJ
p2;p2
xD
ifp2
nC==cNot?cNotNot:cAbsNotNot);p2
e03
l8
0)tN
i91
cNotNot:case
cAbsNotNot:{if(yX3
0)))nY3
cW3
0),i01==cAbsNotNot)y33
IsNever
eL
IsAlways:i6
lT1
l83
cNotNot&&GetPositivityInfo
y91
0))==IsAlways)tI3
cAbsNotNot);lT3
nC==cIf||lE2
0)nC==yW3{eJ
l93=lE2
0);lC2
ifp1=l93
l8
1);lC2
ifp2=l93
l8
2
eY3
ifp1
nC
lA3
ifp1
cS1{cM1
0,l93
lR2
i83
ifp1);eJ
p2;p2
iP2
p2
e03
tN
if(ifp2
nC
lA3
ifp2
cS1{tree
xJ
i01
cU3);cV3
tree
e03);tI3
l93
nC);cZ}
i91
cIf:case
cAbsIf:{if(ConstantFolding_IfOperations(eJ2
yF3
case
cMul:{NowWeAreMulGroup:;AdoptChildrenWithSameOpcode(tree)iU2
nL1=iV2(1);size_t
iW1=0;bool
nM1
tO3
nB1{if(!lE2
t21
continue
iU2
immed=t63
if(immed==xD1
goto
ReplaceTreeWithZero;nL1*=immed;++iW1;}
if(iW1>1||(iW1==1&&cB3
nL1
tG1
1))))nM1=true;if(nM1){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"cMul: Will add new "
iU3
nL1<<"\n"
;
#endif
for
iY
if
y91
t21{
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<" - For that, deleting "
iU3
t63
std::cout<<"\n"
;
#endif
i93!cB3
nL1
tG1
1)))i83
e21
nL1));l73
eM
y33
0:i6
1:nY3
c73
if(ConstantFolding_MulGrouping(eJ2
if(ConstantFolding_MulLogicItems(eJ2
i91
cAdd:eQ
iV2
nE2=0.0;size_t
iW1=0;bool
nM1
tO3
nB1{if(!lE2
t21
continue
iU2
immed=t63
nE2+=immed;++iW1;}
if(iW1>1||(iW1==1&&nE2==xD1)nM1=true;if(nM1){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"cAdd: Will add new "
iU3
nE2<<"\n"
;std::cout<<"In: "
yA1
#endif
for
iY
if
y91
t21{
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<" - For that, deleting "
iU3
t63
std::cout<<"\n"
;
#endif
i93!(nE2==x22)i83
e21
nE2));l73
eM
y33
0
eL
1:nY3
c73
if(ConstantFolding_AddGrouping(eJ2
if(ConstantFolding_AddLogicItems(eJ2
i91
cMin:eQ
size_t
yM2=0;c23
e0;nB1{while(a+1<eM&&t33
xF
lE2
a+1)))tree
t31+1);c33
max
lB3(!e0
yN
known||(e13)<e0
yN
val)){e0
yN
val=e13;e0
yN
known=true;yM2=a;}
}
if(e0
xF2
for
iY{c33
min
lB3
a!=yM2&&p
tF2>=e0
yN
val)i93
eM==1){nY3
i91
cMax:eQ
size_t
yM2=0;c23
eZ;nB1{while(a+1<eM&&t33
xF
lE2
a+1)))tree
t31+1);c33
min
lB3(!eZ
n51||p
tF2>eZ
tF2)){eZ
tF2=p
tF2;eZ
n51=true;yM2=a;}
}
if(eZ
n51){for
iY{c33
max
lB3
a!=yM2&&(e13)<eZ
tF2){eK2}
}
if(eM==1){nY3
i91
cEqual:case
cNEqual:case
cLess:case
cGreater:case
cLessOrEq:case
cGreaterOrEq:if(ConstantFolding_Comparison(eJ2
lC
cAbs:{c23
p0
tI
0));if
tG2
nY3
if(p0
eN{tree
tH3
tree.yC
iV2(1)));goto
NowWeAreMulGroup;}
lT3
nC==cMul){lC2
p=lE2
0);eH
x03;eH
cK2
yM3
0;a<p.xF1++a){p0=iN
p
x13
if
tG2{x03.push_back(p
x13}
if(p0
eN{cK2.push_back(p
x13}
}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Abs: mul group has "
<<x03.tF3<<" pos, "
<<cK2.tF3<<"neg\n"
;
#endif
if(!x03
cT3||!cK2
cT3){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"AbsReplace-Before: "
i13(tree
lT2"\n"
<<std::flush;DumpHashes(tree
i42);
#endif
eJ
eW3;eW3
xD
cMul)yM3
0;a<p.xF1++a){p0=iN
p
x13
if(tG2||(p0
eN){}
else
eW3
c9
p
x13}
eW3
eT1
eJ
x23;x23
xD
cAbs);x23
yP1
eW3);x23
eT1
eJ
y41
cMul);mulgroup
yP1
x23);yQ1
AddParamsMove(x03
eY3!cK2
cT3){if(cK2.tF3%2)yQ1
yC
iV2(-1)));yQ1
AddParamsMove(cK2);}
tree.Become
lR3);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"AbsReplace-After: "
i13(tree
i42
lT2"\n"
<<std::flush;DumpHashes(tree
i42);
#endif
goto
NowWeAreMulGroup;}
}
yF3
#define HANDLE_UNARY_CONST_FUNC(funcname) nR){xO funcname(lR));nD
case
cLog:iK3(fp_log);lT3
nC==cPow){eJ
pow=lE2
0
eY3
GetPositivityInfo(pow
l8
0))==IsAlways){pow.CopyOnWrite()i7
tree.lU
if(GetEvennessInfo(pow
l8
1))==IsAlways){pow.l61
eJ
abs;abs
xD
cAbs);abs
yP1
pow
lR2
abs.Rehash()i7
pow.nG1
0,abs);tree.lU}
else
lT3
nC==cAbs){eJ
pow=lE2
0)l8
0
eY3
pow
nC==cPow){pow.l61
eJ
abs;abs
xD
cAbs);abs
yP1
pow
lR2
abs.Rehash()i7
pow.nG1
0,abs);tree.lU}
lC
cAcosh:iK3(fp_acosh);lC
cAsinh:iK3(fp_asinh);lC
cAtanh:iK3(fp_atanh);lC
cAcos:iK3(fp_acos);lC
cAsin:iK3(fp_asin);lC
cAtan:iK3(fp_atan);lC
cCosh:iK3(fp_cosh);lC
cSinh:iK3(fp_sinh);lC
cTanh:iK3(fp_tanh);lC
cSin:iK3(fp_sin);lC
cCos:iK3(fp_cos);lC
cTan:iK3(fp_tan);lC
cCeil:l54(fp_ceil);lC
cTrunc:l54(fp_trunc);lC
cFloor:l54(fp_floor);lC
cInt:l54(fp_int);lC
cCbrt:iK3(fp_cbrt);lC
cSqrt:iK3(fp_sqrt);lC
cExp:iK3(fp_exp);lC
cLog2:iK3(fp_log2);lC
cLog10:iK3(fp_log10);lC
eF3:if
lI
fp_log2(lR)*tH
tH2
cArg:iK3(fp_arg);lC
cConj:iK3(fp_conj);lC
cImag:iK3(fp_imag);lC
cReal:iK3(fp_real);lC
cPolar:if
lI
fp_polar(x32
cMod:if
lI
e93
x32
cAtan2:{c23
p0
tI
yS2
p1
tI
1));nR&&cB3
lR,xD1){if(p1
cL2
p1
e23
xD1{xO
fp_const_pi
e53
if(p1
n51&&p1
tF2>=t61
xD1;nD}
tM2
cB3
tH,xD1){if(p0
cL2
p0
e23
xD1{xO-fp_const_pihalf
e53
if
iC3
p0
tF2>xD1{xO
fp_const_pihalf
e53}
if
lI
fp_atan2(lR,tH));nD
if((p1
n51&&p1
tF2>xD1||(p1
cL2
p1
e23
fp_const_negativezero
xK())){eJ
yO2;yO2
xD
cPow);yO2.AddParamMove
y91
1));yO2.yC
iV2(-1)));yO2
eT1
eJ
yP2;yP2
tH3
yP2.AddParamMove
y91
0));yP2
yP1
yO2);yP2
eT1
tI3
cAtan);tree.nG1
0,yP2)yO3
1);i91
cPow:{if(ConstantFolding_PowOperations(eJ2
yF3
case
cDiv:nR&&eO3
i61&&tH!=t61
lR/tH
tH2
cInv:nR&&lR!=t61
iV2(1)/lR
tH2
cSub:if
lI
lR-tH
tH2
cNeg:nR){xO-lR
tH2
cRad:nR){xO
RadiansToDegrees(lR)tH2
cDeg:nR){xO
DegreesToRadians(lR)tH2
cSqr:nR){xO
lR*lR
tH2
cExp2:iK3(fp_exp2);lC
cRSqrt:nR){xO
iV2(1)/fp_sqrt(lR)tH2
cCot:tI2
fp_tan(lZ
cSec:tI2
fp_cos(lZ
cCsc:tI2
fp_sin(lZ
cHypot:if
lI
fp_hypot(x32
cRDiv:case
cRSub:case
cDup:case
cFetch:case
cPopNMov:case
cSinCos:case
cSinhCosh:case
cNop:case
cJump:lC
cPCall:case
cFCall:yF3
do_return:;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"["
<<(&cY3)<<"]Done ConstantFolding, result: "
yA1
DumpHashes(tree);
#endif
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
t5{tZ1
c23::set_abs(nL
bool
has_negative=!min.known||min.val<iV2();bool
has_positive=!t13||iA3>iV2();bool
crosses_axis=has_negative&&has_positive;cJ1
xK
newmax;if(min
lB3
t13)newmax.set(fp_max(tV1,tW1
eY3
crosses_axis)min.set(iV2());e63
min
lB3
t13)min.set(fp_min(tV1,tW1);lJ2
min.known)min.set(tV1);else
min.set(tW1;}
max=newmax;}
tZ1
c23::set_neg(){std::swap(min,max);min.val=-min.val;iA3=-iA3;}
xN1
IsLogicalTrueValue
cH1
c23&p
tR2{if(nB
IsIntType
xK::nT3){if(p
n51&&p
tF2>=iV2(1
tI1;if(!abs&&p
nP2
yN
val<=iV2(-1
tI1;}
e63
p
n51&&p
tF2>=iV2(0.5
tI1;if(!abs&&p
nP2
yN
val<=iV2(-0.5
tI1;}
return
tR3
xN1
IsLogicalFalseValue
cH1
c23&p
tR2{if(nB
IsIntType
xK::nT3){if(abs)return
p
yN
known
lH1
1);else
return
p
n51&&p
nP2
tF2>iV2(-1)lH1
1);}
e63
abs)return
p
yN
known
lH1
0.5);else
return
p
n51&&p
nP2
tF2>iV2(-0.5)lH1
0.5);}
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
lF3
FUNCTIONPARSERTYPES;using
t5;t5{xW3
c23
iN
lC2
tree)
#ifdef DEBUG_SUBSTITUTIONS_extra_verbose
{using
lF3
FUNCTIONPARSERTYPES;c23
tmp=CalculateResultBoundaries_do(tree
lT2"Estimated boundaries: "
;if(tmp
n51)std::cout<<tmp
tF2;else
std::cout<<"-inf"
;std::cout<<" .. "
;if(tmp
xF2
std::cout<<tmp
yN
val;else
std::cout<<"+inf"
;std::cout<<": "
i13(tree
lT2
std::endl
n31
tmp;}
xW3
c23
eJ::CalculateResultBoundaries_do
cH1
e62
#endif
{iG
yF1(-fp_const_pihalf
xK(),fp_const_pihalf
xK());iG
pi_limits(-fp_const_pi
xK(),fp_const_pi
xK());iG
abs_pi_limits(y21,fp_const_pi
xK());iG
plusminus1_limits(iV2(-nR3
using
lF3
std;switch(i01
y33
cImmed:nM
eC3,eC3);case
cAnd:case
cAbsAnd:case
cOr:case
cAbsOr:case
cNot:case
cZ3:case
cNotNot:case
cAbsNotNot:case
cEqual:case
cNEqual:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:{nM
y21
tG1
1)i33
cAbs:lD
m.set_abs();tV
cLog:lD
x33
fp_log);m
x63
fp_log);tV
cLog2:lD
x33
fp_log2);m
x63
fp_log2);tV
cLog10:lD
x33
fp_log10);m
x63
fp_log10);tV
cAcosh:lD
yK
template
set_if<cGreaterOrEq
t71
fp_acosh)xQ1
cGreaterOrEq
t71
fp_acosh);tV
cAsinh:lD
yK
set(fp_asinh);m
yN
set(fp_asinh);tV
cAtanh:lD
yK
n3-1),fp_atanh)xQ1
cLess
t71
fp_atanh);tV
cAcos:lD
nM(c71&&(c81)<iV2(1))?fp_acos(c81):y21,(x43&&nU3>=iV2(-1))?fp_acos
nU3:fp_const_pi
xK()i33
cAsin:lD
yK
n3-1),fp_asin,yF1
tF2)xQ1
cLess
t71
fp_asin,yF1
yN
val);tV
cAtan:lD
yK
set(fp_atan,yF1
tF2);m
yN
set(fp_atan,yF1
yN
val);tV
cAtan2:{nR&&cB3
lR,xD1)yR
abs_pi_limits;}
tM2
cB3
tH,xD1)yR
yF1;}
return
pi_limits;}
case
cSin:lD
bool
x11=!x43||!c71||(c81-yK
val)>=(yO
x11)eR
iV2
min=e93
yK
val,yO
min<xD1
min
yP
iV2
max=e93
c81,yO
max<xD1
max
yP
if(max<min)max
yP
bool
y11=(min<=fp_const_pihalf
xK()&&max>=fp_const_pihalf
xK());bool
nN1=(min<=cD&&max>=cD
eY3
y11&&nN1)eR
if(nN1)nM
iV2(-1)lK2
if(y11)nM
eL2
tG1
1));nM
eL2
lK2}
case
cCos:lD
if(x43)yK
val+=fp_const_pihalf
xK(eY3
c71)c81+=fp_const_pihalf
xK();bool
x11=!x43||!c71||(c81-yK
val)>=(yO
x11)eR
iV2
min=e93
yK
val,yO
min<xD1
min
yP
iV2
max=e93
c81,yO
max<xD1
max
yP
if(max<min)max
yP
bool
y11=(min<=fp_const_pihalf
xK()&&max>=fp_const_pihalf
xK());bool
nN1=(min<=cD&&max>=cD
eY3
y11&&nN1)eR
if(nN1)nM
iV2(-1)lK2
if(y11)nM
eL2
tG1
1));nM
eL2
lK2}
case
cTan:{nM
i33
cCeil:lD
c91
cFloor:lD
tH1
tV
cTrunc:lD
tH1
c91
cInt:lD
tH1
c91
cSinh:lD
yK
set(fp_sinh);m
yN
set(fp_sinh);tV
cTanh:lD
yK
set(fp_tanh,plusminus1_limits.min);m
yN
set(fp_tanh,plusminus1_limits.max);tV
cCosh:lD
if(x43){if(c71){if(yK
val>=y21&&c81>=xD1{yK
val
xU}
lJ2
nU3<y21&&c81>=xD1{iV2
tmp
xU
if(tmp>c81)c81=tmp;yK
val=iV2(1
t43{yK
val
xU
std::swap(yK
val,c81);}
}
e63
yK
val>=xD1{m
t0
yK
val=fp_cosh
nU3;}
else{m
t0
yK
val=iV2(1);}
}
}
else{x43=true;yK
val=iV2(1
eY3
c71){yK
val=fp_cosh(c81);m
t0}
else
m
t0}
tV
cIf:case
cAbsIf:{c23
res1
tI
1));c23
res2
tI
2)eY3!res2
n51)res1
n51
tO3
lJ2
res1
n51&&(res2
tF2)<res1
tF2)res1
tF2=res2
tF2;if(!res2
xF2
res1
t0
lJ2
res1
cL2
res2
yN
val)>res1
yN
val)res1
yN
val=res2
yN
val
n31
res1;}
case
cMin:{bool
iH
tO3
bool
iI
tO3
c43;x2
m
tI
a)eY3!x43)iH
x53
n51||nU3<nK3)nK3=yK
val;if(!c71)iI
x53
yN
known||(c81)<nL3)nL3=c81;}
if(iH)nM3
iI)nT3
t0
return
nT3;}
case
cMax:{bool
iH
tO3
bool
iI
tO3
c43;x2
m
tI
a)eY3!x43)iH
x53
n51||yK
val>nK3)nK3=yK
val;if(!c71)iI
x53
yN
known||c81>nL3)nL3=c81;}
if(iH)nM3
iI)nT3
t0
return
nT3;}
case
cAdd:{c43(y21,xD1;x2
item
tI
a)eY3
item
n51)nK3+=item
tF2;else
nM3
item
xF2
nL3+=item
yN
val;else
nT3
t0
if(!nP3&&!nT3
xF2
yF3
if
iT2
n51&&nQ3&&nK3>nL3)std::swap
iT2
tF2,nL3)n31
nT3;}
case
cMul:{e12
Value{enum
x73{eA3,iX1,iW2}
;x73
tK
iU2
value;Value(x73
t):tK(t),value(0){}
Value(iV2
v):tK(eA3),value(v){}
bool
cM2
yZ1
tK==iX1||cU2
value<xD1
lL3
eA1*=cH1
Value&rhs){if
cU2
rhs.tK==eA3)value*=rhs.value;else
tK=(cM2)!=rhs.cM2)?iX1:iW2);}
iG2<cH1
Value&rhs
yZ1(tK==iX1&&rhs.tK!=iX1)||cU2(rhs.tK==iW2||(rhs.tK==eA3&&value<rhs.value)));}
}
;e12
yG1{Value
yQ2,yR2;yG1():yQ2(Value::iW2),yR2(Value::iX1){}
void
xE2
Value
eB3,const
Value&value2){eB3*=value2;if(eB3<yQ2)yQ2=eB3;if(yR2<eB3)yR2=eB3;}
}
;c43(iV2(nR3
x2
item
tI
a)eY3!item
n51&&!item
xF2
nM);Value
x83=nP3?Value
iT2.min.lU1
iX1);Value
x93=nQ3?Value
iT2
yN
lU1
iW2);Value
xA3=item
n51?Value(item.min.lU1
iX1);Value
xB3=item
yN
known?Value(item
yN
lU1
iW2);yG1
range
tL2
x83,xA3)tL2
x83,xB3)tL2
x93,xA3)tL2
x93,xB3
eY3
range.yQ2.tK==Value::eA3)nK3=range.yQ2.value;else
nM3
range.yR2.tK==Value::eA3)nL3=range.yR2.value;else
nT3
t0
if(!nP3&&!nT3
xF2
yF3
if
iT2
n51&&nQ3&&nK3>nL3)std::swap
iT2
tF2,nL3)n31
nT3;}
case
cMod:{c23
x
tI
yS2
y
tI
1)eY3
y
xF2{if(y
yN
val>=xD1{if(!x
n51||(x
tF2)<xD1
nM-y
yN
val,y
yN
val);iB3
y21,y
yN
val);}
e63!x
yN
known||(x
yN
val)>=xD1
nM
y
yN
val,-y
yN
val);iB3
y
yN
val,fp_const_negativezero
xK());}
}
iB3
i33
cPow:{tM2
tH==xD1{nM
iV2(nR3}
nR&&lR==xD1{nM
y21,xD1;}
nR&&cB3
lR
xE3
nM
iV2(nR3}
tM2
tH>y21&&GetEvennessInfo
y91
1))==IsAlways
tJ2
e92
tH;c23
tmp
tI
yS2
nT3;nP3=true;nK3=0;if(tmp
n51&&tmp
tF2>=xD1
nK3
t03
tmp
tF2,xG2;lJ2
tmp
yN
known&&tmp
yN
val<=xD1
nK3
t03
tmp
yN
val,xG2;nT3
t0
if(tmp
n51&&tmp
xF2{nQ3=true;nL3=fp_max(fp_abs(tmp
tF2),fp_abs(tmp
yN
val));nL3=fp_pow
iT2
yN
val,xG2;}
return
nT3;}
c23
p0
tI
yS2
p1
tI
1));TriTruthValue
p0_positivity=iC3(p0
tF2)>=xD1?IsAlways:(p0
cL2
p0
e23
y21?lC3
Unknown);TriTruthValue
cN2=GetEvennessInfo
y91
1));TriTruthValue
t1=Unknown;switch(p0_positivity)l92
t1=IsAlways;lC
lC3{t1=cN2;yF3
c73
switch(cN2)l92
t1=IsAlways;lC
lC3
lC
Unknown:{tM2!isInteger
tG&&tH>=xD1{t1=IsAlways;}
yF3}
l73
t1)l92{iV2
min=y21;if
iC3
p1
n51){min
t03
p0
tF2,p1
tF2
eY3
p0
tF2<y21&&(!p1
yN
known||p1
yN
val>=xD1&&min>=xD1
min=y21;}
if
iC3
p0
tF2>=y21&&p0
yN
known&&p1
xF2{iV2
max
t03
p0
yN
val,p1
yN
val
eY3
min>max)std::swap(min,max);nM
min,max);}
nM
min,false
i33
lC3{nM
false,fp_const_negativezero
xK());}
c73{yF3
i91
cNeg:lD
m.set_neg();tV
cSub:{yM
cNeg);iE3
1));tmp
xD
cAdd);tmp.nJ
0));tmp
yP1
tmp2
y53
cInv:{cO2-1
yT2
cDiv:{yM
cInv);iE3
1));tmp
xD
x71
AddParamMove(tmp2
y53
cRad:{t81
x71
yC
fp_const_rad_to_deg
xK(yT2
cDeg:{t81
x71
yC
fp_const_deg_to_rad
xK(yT2
cSqr:{cO2
2
yT2
cExp:{t81
cPow);tmp.yC
fp_const_e
xK()));tmp.nJ
0)y53
cExp2:{t81
cPow);tmp.yC
x42
tmp.nJ
0)y53
cCbrt:lD
yK
set(fp_cbrt);m
yN
set(fp_cbrt);tV
cSqrt:lD
if(x43)yK
val=nU3<y21?0:fp_sqrt
nU3;if(c71)c81=(c81)<y21?0:fp_sqrt(c81);tV
cRSqrt:{cO2-0.5
yT2
cHypot:{eJ
xsqr,ysqr,add,sqrt;xsqr.nJ
0));xsqr.yC
x42
ysqr.nJ
1));ysqr.yC
x42
xsqr
xD
cPow);ysqr
xD
cPow);add
yP1
xsqr);add
yP1
ysqr);add
xD
cAdd);sqrt
yP1
add);sqrt
xD
cSqrt)n31
iN
sqrt
i33
eF3:{yM
cLog2);iE3
0));tmp
tH3
tmp
yP1
tmp2);tmp.nJ
1)y53
cCot:{yM
cTan)x3
lH
cSec:{yM
cCos)x3
lH
cCsc:{yM
cSin)x3
iN
tmp);}
lC
cRDiv:case
cRSub:case
cDup:case
cFetch:case
cPopNMov:case
cSinCos:case
cSinhCosh:case
cNop:case
cJump:case
l03:lC
cArg:case
cConj:case
cImag:case
cReal:case
cPolar:lC
cPCall:lC
cFCall:yF3
nM);}
xW3
TriTruthValue
GetIntegerInfo
cH1
e62{switch(i01
y33
cImmed:return
t42
eC3)?IsAlways:IsNever;case
cFloor:case
cCeil:case
cTrunc:case
cInt:return
IsAlways;case
cAnd:case
cOr:case
cNot:case
cNotNot:case
cEqual:case
cNEqual:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:return
IsAlways;case
cIf:{TriTruthValue
a=GetIntegerInfo
y91
1));TriTruthValue
b=GetIntegerInfo
y91
2)eY3
a==b)return
a
nV2
case
cAdd:case
cMul:{for
iY
if(GetIntegerInfo
y91
a))!=IsAlways)return
Unknown
n31
IsAlways;}
c73
yF3
return
Unknown;}
xN1
IsLogicalValue
cH1
e62{switch(i01
y33
cImmed:return
cB3
eC3,xD1||cB3
eC3
tG1
1));case
cAnd:case
cOr:case
cNot:case
cNotNot:case
cAbsAnd:case
cAbsOr:case
cZ3:case
cAbsNotNot:case
cEqual:case
cNEqual:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:nV
cMul:{for
iY
if(!yX3
a))cL
return
true;}
case
cIf:case
cAbsIf:yR
yX3
1))nT1
lE2
2));}
c73
yF3
return
tR3}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
lF3
FUNCTIONPARSERTYPES;
#if defined(__x86_64) || __cplusplus < 201100
# define CBRT_IS_SLOW
#endif
#if defined(DEBUG_POWI) || defined(DEBUG_SUBSTITUTIONS)
#include <cstdio>
#endif
lF3
xH1{extern
const
iF2
char
powi_table[256];}
lF3{using
t5
nJ1
bool
IsOptimizableUsingPowi(long
immed,long
penalty=0){xH1
yG3
synth;synth.PushVar(l03);size_t
bytecodesize_backup=synth.GetByteCodeSize();xH1::x01
immed,xH1::tN1
xK::MulSequence,iS2
size_t
bytecode_grow_amount=synth.GetByteCodeSize()-bytecodesize_backup
n31
bytecode_grow_amount<size_t(MAX_POWI_BYTECODE_LENGTH-penalty);}
tZ1
ChangeIntoRootChain(eJ&tree,bool
lD3,long
i52,long
i62){while(i62>0){t81
cCbrt);iF3
tmp
eT1
tree.i72--i62;}
while(i52>0){t81
cSqrt
eY3
lD3){tmp
xD
cRSqrt);lD3
tO3}
iF3
tmp
eT1
tree.i72--i52;}
if(lD3){t81
cInv);iF3
tree.i72}
}
xW3
e12
RootPowerTable{static
const
iV2
RootPowers[(1+4)*(1+3)];}
nJ1
const
iV2
t8(1+4)*(1+3)]={iV2(1)lS
iG3
iG3
2*iG3
i11)lS
3*2)lS
3*2*2)lS
3*i11*2*i11*3
cV2
2
cV2
2*2
cV2
i11*3*2*i11*3*3
cV2
3*2
cV2
3*2*2
cV2
3*i11*3*3*2*2*2*2)}
;e12
PowiResolver{static
const
iF2
MaxSep=4;static
xC3
MaxOp=5;typedef
int
eE3;typedef
long
xI3;typedef
long
tW;e12
yU2{yU2():n_int_sqrt(0),n_int_cbrt(0),sep_list(),lV1(0){}
int
n_int_sqrt;int
n_int_cbrt;int
tL1
MaxSep];tW
lV1;}
nJ1
static
yU2
CreatePowiResult(iV2
xG2{yU2
nT3;eE3
t9=FindIntegerFactor(xG2;if(t9==0){
#ifdef DEBUG_POWI
i82"no factor found for %Lg\n"
,yK1
xG2;
#endif
return
nT3;}
nT3.lV1=y31
exponent,t9);xI3
eM2=EvaluateFactorCost(t9,0,0,0)+cB
nT3.lV1);int
iH3=0;int
iI3=0;int
mul_count=0;
#ifdef DEBUG_POWI
i82"orig = %Lg\n"
,yK1
xG2;i82"plain factor = "
iW3"%ld\n"
,(int)t9,(long)eM2);
#endif
for
nT2
n_s=0;n_s<MaxSep;++n_s){int
x9=0;xI3
yI1=eM2;eE3
yW1=t9;for(int
s=1;s<MaxOp*4;++s){
#ifdef CBRT_IS_SLOW
if(s>=MaxOp)break;
#endif
int
n_sqrt=s%MaxOp;int
n_cbrt=s/MaxOp;if(n_sqrt+n_cbrt>4
cP2
iU2
lI1=exponent;lI1-=t8
s];i41=FindIntegerFactor(lI1
eY3
xT2!=0){tW
xP=y31
lI1,xT2);xI3
cost=EvaluateFactorCost(xT2,iH3+n_sqrt,iI3+n_cbrt,mul_count+1)+cB
xP);
#ifdef DEBUG_POWI
i82"Candidate sep %u (%d*sqrt %d*cbrt)factor = "
iW3"%ld (for %Lg to %ld)\n"
,s,n_sqrt,n_cbrt,xT2,(long)cost,yK1
lI1,(long)xP);
#endif
if(cost<yI1){x9=s;yW1=xT2;yI1=cost;}
}
}
if(!x9)break;
#ifdef DEBUG_POWI
i82"CHOSEN sep %u (%d*sqrt %d*cbrt)factor = "
iW3"%ld, exponent %Lg->%Lg\n"
,x9,x9%MaxOp,x9/MaxOp,yW1,yI1,yK1(xG2,yK1
t62-t8
x9]));
#endif
nT3.tL1
n_s]=x9;exponent-=t8
x9];iH3+=x9%MaxOp;iI3+=x9/MaxOp;eM2=yI1;t9=yW1;mul_count+=1;}
nT3.lV1=y31
exponent,t9);
#ifdef DEBUG_POWI
i82"resulting exponent is %ld (from exponent=%Lg, best_factor=%Lg)\n"
,nT3.lV1,yK1
exponent,yK1
t9);
#endif
while(t9%2==0){++nT3
tD2;t9/=2;}
while(t9%3==0){++nT3.n_int_cbrt;t9/=3;}
return
nT3;}
private:static
xI3
cB
tW
xP){static
std::map
cR2
iB;if(xP<0){xI3
cost=22
n31
cost+cB-xP);}
std::map
cR2::yV3
i=iB.xU2
xP
eY3
i!=iB.cY1
xP)return
i
e82;std::pair
cR2
nT3(xP,0.0);xI3&cost=nT3
tE3;while(xP>1){int
xT2=0;if(xP<256){xT2=xH1::powi_table[xP];if(xT2&128)xT2&=127;else
xT2=0;if(xT2&64)xT2=-(xT2&63)-1;}
if(xT2){cost+=cB
xT2);xP/=xT2;eR2
if(!(xP&1)){xP/=2;cost+=6;}
else{cost+=7;xP-=1;}
}
iB.yR3,nT3)n31
cost;}
cA1
tW
y31
yJ1,i41)yR
makeLongInteger(value*iV2(xT2));}
cA1
bool
yR1
yJ1,i41
tJ2
v=value*iV2(xT2)n31
isLongInteger(v);}
cA1
eE3
FindIntegerFactor(yJ1){i41=(2*2*2*2);
#ifdef CBRT_IS_SLOW
#else
xT2*=(3*3*3);
#endif
eE3
nT3=0;if(yR1
value,xT2)){nT3=xT2;while((xT2%2)==0&&yR1
value,xT2/2))nT3=xT2/=2;while((xT2%3)==0&&yR1
value,xT2/3))nT3=xT2/=3;}
#ifdef CBRT_IS_SLOW
if
iT2==0){if(yR1
value,3
nE3
3;}
#endif
return
nT3;}
static
int
EvaluateFactorCost(int
xT2,int
s,int
c,int
nmuls){xC3
xD3=6;
#ifdef CBRT_IS_SLOW
xC3
eN2=25;
#else
xC3
eN2=8;
#endif
int
nT3=s*xD3+c*eN2;while(xT2%2==0){xT2/=2;nT3+=xD3;}
while(xT2%3==0){xT2/=3;nT3+=eN2;}
nT3+=nmuls
n31
nT3;}
}
;}
t5{xN1
eJ::RecreateInversionsAndNegations(bool
prefer_base2){bool
changed=false
yM3
0;a<xF1++a)if(lS1.RecreateInversionsAndNegations(prefer_base2))x62
if(changed){exit_changed:Mark_Incompletely_Hashed()xJ2
switch(iS1{case
cMul:{eH
nF2;eJ
nG2,cT1;if(true){bool
nO1=false
iU2
x52=0;yP3
yV2
0)cX2
tO
iJ3
nO1=true;x52=tO
1)nC1;yF3}
if(nO1
tJ2
immeds=1.0;yP3
i61){immeds*=powgroup
nC1;t41}
for
iJ-->0;){eJ&powgroup=lS1;if(powgroup
yV2
0)cX2
tO
iJ3
eJ&log2=tO
0);log2.l61
log2
xD
eF3);log2.yC
fp_pow(immeds
tG1
1)/x52)));log2
eT1
yF3}
}
}
yP3
yV2
iJ3
lC2
exp_param=tO
1)iU2
e92
exp_param
nC1;if(e01
tG1-1))){l61
nF2.push_back(lS1
lR2
t41
lJ2
exponent<y21&&t42
xG2){eJ
iK;iK
xD
cPow);iK
c9
tO
0));iK.yC-xG2);iK
eT1
nF2.push_back(iK);l61
t41}
lJ2
powgroup
cX2!nG2.t91
nG2=tO
0);l61
t41
lJ2
powgroup
nC==eF3&&!cT1.t91
cT1=powgroup;l61
t41}
if(!nF2
cT3){x62
eJ
i21;i21
tH3
i21
iO3
nF2);i21
eT1
eJ
y41
cMul);yQ1
SetParamsMove
eF
if(yQ1
IsImmed()&&fp_equal
lR3
nC1
xE3
nH2
cInv)eO
i21);}
e63
yQ1
y22>=i21.y22){nH2
cDiv)eO
nN3
eO
i21
t43{nH2
cRDiv)eO
i21)eO
nN3;}
}
}
if(nG2.t91
eJ
y41
iS1;yQ1
SetParamsMove
eF
while(yQ1
RecreateInversionsAndNegations(prefer_base2))yQ1
FixIncompleteHashes();nH2
eF3)eO
nG2)eO
nN3;x62}
if(cT1.t91
eJ
y41
cMul);mulgroup
yP1
cT1
l8
1));yQ1
AddParamsMove
eF
while(yQ1
RecreateInversionsAndNegations(prefer_base2))yQ1
FixIncompleteHashes();DelParams();nH2
eF3)eO
cT1
l8
0))eO
nN3;x62
i91
cAdd:{eH
iA2;for
iJ-->0;)if(eG3
cMul){nV3
y51:;eJ&mulgroup=x82
for
iA1
b=yQ1
xF1
b-->0;){if
lR3
l8
b)l71
xT2=mulgroup
l8
b)nC1
cA3
xT2
i92
y51;}
yQ1
l61
yQ1
DelParam(b);nW3
lJ2
cB3
xT2
tG1-2)))xC
y51;}
yQ1
l61
yQ1
DelParam(b);yQ1
yC
x42
nW3}
}
if(t2){yQ1
tA
nN3;t41}
lJ2
eG3
cDiv&&!IsIntType
xK::nT3){nV3
y61:;eJ&i21=x82
if(i21
l8
0)i51
cB3
i21
l8
0)nC1
i92
y61;}
i21.l61
i21
tK2
0);i21
xD
cInv);nW3}
if(t2)xC
y61;}
i21.tA
i21);t41}
lJ2
eG3
cRDiv&&!IsIntType
xK::nT3){nV3
x81:;eJ&i21=x82
if(i21
l8
1)i51
cB3
i21
l8
1)nC1
i92
x81;}
i21.l61
i21
tK2
1);i21
xD
cInv);nW3}
if(t2)xC
x81;}
i21.tA
i21);t41}
if(!iA2
cT3){
#ifdef DEBUG_SUBSTITUTIONS
i82"Will make a Sub conversion in:\n"
);fflush(stdout);iR
#endif
eJ
yW2;yW2
xD
cAdd);yW2
iO3
iA2);yW2
eT1
eJ
cU1;cU1
xD
cAdd);cU1
iO3
l02));cU1
eT1
if(cU1
i61&&cB3
cU1
nC1,xD1){nH2
cNeg);e7);}
e63
cU1.y22==1){nH2
cRSub);e7)eO
cU1);}
lJ2
yW2
nC==cAdd){nH2
cSub)eO
cU1);e7
lR2
for
iA1
a=1;a<yW2.xF1++a){eJ
eO2;eO2
xD
cSub);eO2
iO3
l02));eO2.Rehash(false)eO
eO2);e7
x13}
}
else{nH2
cSub)eO
cU1);e7);}
}
#ifdef DEBUG_SUBSTITUTIONS
i82"After Sub conversion:\n"
);fflush(stdout);iR
#endif
i91
cPow:{lC2
p0=GetParam(0);lC2
p1=GetParam(1
eY3
p1
i51
p1
nC1!=y21&&!t42
p1
nC1)){eC
yU2
r=eC
CreatePowiResult(fp_abs
eH3)eY3
r.lV1!=0){bool
iY1
tO3
if
eH3<y21&&r.tL1
0]==0&&r
tD2>0){iY1=true;}
#ifdef DEBUG_POWI
i82"Will resolve powi %Lg as powi(chain(%d,%d),%ld)"
,yK1
fp_abs
eH3),r
tD2,r.n_int_cbrt,r.lV1);for
nT2
n=0;n<eC
MaxSep;++n){if(r
nS3==0)break;int
n_sqrt=r
nS3%eC
MaxOp;int
n_cbrt=r
nS3/eC
MaxOp;i82"*chain(%d,%d)"
,n_sqrt,n_cbrt);}
i82"\n"
);
#endif
eJ
cY2=GetParam(0);eJ
yX2=cY2;yX2.l61
ChangeIntoRootChain(yX2,iY1,r
tD2,r.n_int_cbrt);yX2
eT1
eJ
pow;if(r.lV1!=1){pow
xD
cPow);pow
yP1
yX2);pow.yC
iV2(r.lV1))t43
pow.swap(yX2);eJ
mul;mul
tH3
mul
yP1
pow);for
nT2
n=0;n<eC
MaxSep;++n){if(r
nS3==0)break;int
n_sqrt=r
nS3%eC
MaxOp;int
n_cbrt=r
nS3/eC
MaxOp;eJ
eP2=cY2;eP2.l61
ChangeIntoRootChain(eP2,false,n_sqrt,n_cbrt);eP2
eT1
mul
yP1
eP2);}
if
eH3<y21&&!iY1){mul
eT1
nH2
cInv);nG1
0,mul);DelParam(1
t43{nH2
cMul);SetParamsMove(mul.l02));}
#ifdef DEBUG_POWI
iR
#endif
x62
yF3}
}
if(GetOpcode()==cPow&&(!p1
i61||!isLongInteger
eH3)||!IsOptimizableUsingPowi
xK(makeLongInteger
eH3)))){if(p0
i61&&p0
nC1>x22{if(prefer_base2
tJ2
yY2=fp_log2(p0
nC1)cA3
yY2
xE3
DelParam(0
t43{n0
e21
yY2)tF1
c9
p1
tF1
eT1
nG1
nD1}
nH2
cExp2);x62}
else{iV2
yY2=fp_log(p0
nC1)cA3
yY2
xE3
DelParam(0
t43{n0
e21
yY2)tF1
c9
p1
tF1
eT1
nG1
nD1}
nH2
cExp);x62}
}
lJ2
GetPositivityInfo(p0)==IsAlways){if(prefer_base2){eJ
log;log
xD
cLog2);log
c9
p0);log
eT1
n0
p1
tF1
yP1
log
tF1
eT1
nH2
cExp2);nG1
nD1
x62}
else{eJ
log;log
xD
cLog);log
c9
p0);log
eT1
n0
p1
tF1
yP1
log
tF1
eT1
nH2
cExp);nG1
nD1
x62}
}
i91
cDiv:{if(GetParam(0)i61&&cB3
GetParam(0)nC1
xE3
nH2
cInv);DelParam(0);}
yF3
c73
yF3
if(changed)goto
exit_changed
n31
changed;}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
lF3
FUNCTIONPARSERTYPES;lF3{using
t5;class
e43{size_t
nP1;size_t
eD;size_t
eE;size_t
lJ1;size_t
t3;size_t
t4;size_t
n01;eX3
e43():nP1(0),eD(0),eE(0),lJ1(0),t3(0),t4(0),n01(0){}
void
t23
OPCODE
op){nP1+=1;e33
cCos)++eD;e33
cSin)++eE;e33
cSec)++eD;e33
cCsc)++eE;e33
cTan)++lJ1;e33
cCot)++lJ1;e33
cSinh)++t4;e33
cCosh)++t3;e33
cTanh)++n01;}
size_t
GetCSEscore
e73
size_t
nT3=nP1
n31
nT3;}
int
NeedsSinCos
e73
bool
y71=(nP1==(eD+eE+lJ1)eY3(lJ1&&(eE||eD))||(eE&&eD)){if(y71)return
1
n31
2;}
return
0;}
int
NeedsSinhCosh
e73
bool
y71=(nP1==(t3+t4+n01)eY3(n01&&(t4||t3))||(t4&&t3)){if(y71)return
1
n31
2;}
return
0;}
size_t
MinimumDepth
e73
size_t
n_sincos=std::min(eD,eE);size_t
n_sinhcosh=std::min(t3,t4
eY3
n_sincos==0&&n_sinhcosh==0)return
2
n31
1;}
}
nJ1
class
TreeCountType:public
std::multimap<fphash_t,std::pair<e43,eJ> >{}
xR3
FindTreeCounts(tK1&nI2,lC2
tree,OPCODE
x92,bool
skip_root=false){cV
i=nI2.xU2
tree.GetHash()eY3!skip_root){bool
found
tO3
for(;i!=nI2.cY1
tree.GetHash();++i){if(tree
xF
i
e82
iO2){i
e82.first.t23
x92);found=true;yF3}
if(!found){e43
count;count.t23
x92);nI2.yR3,std::make_pair(tree.GetHash(),std::make_pair
eL3,tree)));}
}
nB1
FindTreeCounts(nI2,t33,i01);}
e12
yU{bool
BalanceGood;bool
FoundChild;}
nJ1
yU
lK1
lC2
root,lC2
child){if(root
xF
child)){yU
nT3={true,true}
n31
nT3;}
yU
nT3={true,false}
;if(root
nC==cIf||root
nC==yW3{yU
cond=lK1
root
l8
0
nO3
yU
xX=lK1
root
l8
1
nO3
yU
y6=lK1
root
l8
2
nO3
if
y03||xX
yV||y6
yV){nT3
yV=true;}
nT3
e8=((xX
yV==y6
yV)||y03
cZ2&&(cond
e8||(xX
yV&&y6
yV))&&(xX
e8||y03
cZ2&&(y6
e8||y03
cZ2;}
else{bool
iE1
tO3
bool
nQ1
tO3
for
iA1
b=root.GetParamCount(),a=0;a<b;++a){yU
tmp=lK1
root
l8
a
nO3
if(tmp
yV)nT3
yV=true;if(tmp
e8==false)iE1=true;lJ2
tmp
yV)nQ1=true;}
if(iE1&&!nQ1)nT3
e8
tO3}
return
nT3;}
xN1
nZ3
lQ3
within,lC2
tree,const
xH1
yG3&synth,const
tK1&nI2){for
iA1
b=eM,a=0;a<b;++a){lC2
leaf=t33;cV
synth_it;xX2
tK1::const_iterator
i=nI2.yU3
i!=nI2.end();++i){if(i->first!=leaf.GetHash()cP2;const
e43&occ
nQ2
first;size_t
score=occ.GetCSEscore();lC2
candidate
nQ2
second;if(tM1
candidate)cP2;if(leaf.y22<occ.MinimumDepth()cP2;if(score<2
cP2;if(lK1
within,leaf)e8==false
cP2
xJ2
if(nZ3(within,leaf,synth,nI2
tI1;}
return
tR3
xN1
nJ2
lQ3
yS3,lC2
expr){i71
yS3
l8
a)xF
expr
tI1;i71
nJ2(yS3
l8
a),expr
tI1
n31
tR3
xN1
GoodMomentForCSE
lQ3
yS3,lC2
expr){if(yS3
nC==cIf)return
true;i71
yS3
l8
a)xF
expr
tI1;size_t
iB2=0;i71
nJ2(yS3
l8
a),expr))++iB2
n31
iB2!=1;}
}
t5{xW3
size_t
eJ::SynthCommonSubExpressions(xH1::yM1
const{if(GetParamCount()==0)return
0;size_t
stacktop_before=synth.GetStackTop();tK1
nI2;FindTreeCounts(nI2,*this,GetOpcode(),true);for(;;){size_t
yZ2=0;
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"Finding a CSE candidate, root is:"
<<std::endl;DumpHashes(*this);
#endif
cV
cs_it(nI2.end());for(cV
j=nI2.yU3
j!=nI2.end();){cV
i(j++);const
e43&occ
nQ2
first;size_t
score=occ.GetCSEscore();lC2
tree
nQ2
second;
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"Score "
<<score<<":\n"
<<std::flush;DumpTreeWithIndent(tree);
#endif
if(tM1
tree))xT
if(tree.y22<occ.MinimumDepth())xT
if(score<2)xT
if(lK1*this,tree)e8==false)xT
if(nZ3(*this,tree,synth,nI2)){eR2
if(!GoodMomentForCSE(*this,tree))xT
score*=tree.y22;if(score>yZ2){yZ2=score;cs_it=i;}
}
if(yZ2<=0){
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"No more CSE candidates.\n"
<<std::flush;
#endif
yF3
lC2
tree=cs_it
e82
tE3;
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<l24"Common Subexpression:"
i13
xK(tree
lT2
std::endl;
#endif
#if 0
int
lW1=occ.NeedsSinCos();int
i1=occ.NeedsSinhCosh();eJ
iC2,iD2,c02,c12;if(lW1){iC2
eQ2
iC2
xD
cSin);iC2
eT1
iD2
eQ2
iD2
xD
cCos);iD2
eT1
if(tM1
iC2)||tM1
iD2)){if(lW1==2){tJ1
eR2
lW1=0;}
}
if(i1){c02
eQ2
c02
xD
cSinh);c02
eT1
c12
eQ2
c12
xD
cCosh);c12
eT1
if(tM1
c02)||tM1
c12)){if(i1==2){tJ1
eR2
i1=0;}
}
#endif
tree.SynthesizeByteCode(synth,false);tJ1
#ifdef DEBUG_SUBSTITUTIONS_CSE
synth.template
Dump<0>(lT2"Done with Common Subexpression:"
i13
xK(tree
lT2
std::endl;
#endif
#if 0
if(lW1){if(lW1==2||i1){synth.e91}
yE3
cSinCos,1,2)cI1
iC2,1)cI1
iD2,0);}
if(i1){if(lW1)synth.e91
if(i1==2){synth.e91}
yE3
cSinhCosh,1,2)cI1
c02,1)cI1
c12,0);}
#endif
}
return
yD3
stacktop_before;}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
tZ1
FunctionParserBase
xK::Optimize(){using
t5;l61
eJ
tree;tree.GenerateFrom(*mData);FPoptimizer_Optimize::ApplyGrammars(tree);std
xV3<iF2>eI3;std
xV3
xK
immed;size_t
stacktop_max=0;tree.SynthesizeByteCode(eI3,immed,stacktop_max
eY3
mData->mStackSize!=stacktop_max){mData->mStackSize=iF2(stacktop_max);
#if !defined(FP_USE_THREAD_SAFE_EVAL) && \
    !defined(FP_USE_THREAD_SAFE_EVAL_WITH_ALLOCA)
mData->mStack.t73
stacktop_max);
#endif
}
mData->mByteCode.swap(eI3);mData->mImmed.swap(immed);}
#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
x8
MpfrFloat
n11
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
x8
GmpInt
n11
#endif
#ifdef FP_SUPPORT_COMPLEX_DOUBLE_TYPE
x8
std::complex<double>n11
#endif
#ifdef FP_SUPPORT_COMPLEX_FLOAT_TYPE
x8
std::complex<float>n11
#endif
#ifdef FP_SUPPORT_COMPLEX_LONG_DOUBLE_TYPE
x8
std::complex<long
double>n11
#endif
FUNCTIONPARSER_INSTANTIATE_TYPES
#endif

#endif
