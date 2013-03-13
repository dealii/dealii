/***************************************************************************\
|* Function Parser for C++ v4.5.1                                          *|
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
#define l14 {Value_t
#define l04 ),info,
#define iZ3 {data->
#define iY3 ==cNot||
#define iX3 "Found "
#define iW3 info.lQ[b
#define iV3 ;o<<"\n";
#define iU3 )yR 2,cPow
#define iT3 ;}static yU1
#define iS3 (tree)!=
#define iR3 ),Value(
#define iQ3 stackpos
#define iP3 minValue0
#define iO3 "dup(%u) "
#define iN3 =true;iO1
#define iM3 eR{assert
#define iL3 "%d, cost "
#define iK3 "PUSH " iL2
#define iJ3 "immed "<<
#define iI3 mFuncParsers
#define iH3 param.data
#define iG3 stderr
#define iF3 sep2=" "
#define iE3 FPHASH_CONST
#define iD3 cache_needed[
#define iC3 fprintf
#define iB3 "Applying "
#define iA3 FUNCTIONPARSER_INSTANTIATE_OPTIMIZE
#define i93 FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE
#define i83 HANDLE_UNARY_CONST_FUNC
#define i73 {if(n11
#define i63 second;
#define i53 within,
#define i43 AddFrom(
#define i33 tH3==
#define i23 c_count
#define i13 MaxOp
#define i03 else nM
#define tZ3 =tZ cL3
#define tY3 .nJ 0));
#define tX3 .nJ 1));
#define tW3 lD m.
#define tV3 ;a<tree.xB
#define tU3 iZ);}if
#define tT3 cS ifp2
#define tS3 sim.xB1
#define tR3 );tO tL1
#define tQ3 ].swap(
#define tP3 codes[b
#define tO3 whydump
#define tN3 for lS1
#define tM3 for(;a<
#define tL3 nparams
#define tK3 281856,
#define tJ3 cTan i7
#define tI3 l2 2,2,
#define tH3 ;if(op
#define tG3 l4 0,1,
#define tF3 0x12 nH
#define tE3 FixIncompleteHashes()
#define tD3 cSinh,
#define tC3 ,cTan nS
#define tB3 ,iX 1,
#define tA3 nS 0,
#define t93 cAbs nS
#define t83 ;case
#define t73 i1 t83
#define t63 Params[
#define t53 Params(
#define t43 &&p nL3
#define t33 fp_pow(
#define t23 false;}
#define t13 lD1 y4.
#define t03 cAbsIf)
#define eZ3 l91++b)
#define eY3 lX1 nC==
#define eX3 ;exponent
#define eW3 =false;
#define eV3 tN cS1);
#define eU3 .empty()
#define eT3 nS3))n92
#define eS3 (opcode==
#define eR3 <data eJ3;++
#define eQ3 l14 nS3(
#define eP3 ByteCodeSynth xF
#define eO3 std xM3<
#define eN3 middle2
#define eM3 std::string
#define eL3 SetOpcode(
#define eK3 ;for iZ1
#define eJ3 .size()
#define eI3 ].second
#define eH3 ].first
#define eG3 return p
#define eF3 TreeCountItem
#define eE3 }return
#define eD3 Ne_Mask
#define eC3 Gt_Mask
#define eB3 Lt_Mask
#define eA3 const eK
#define e93 .n7 synth.
#define e83 range nE3
#define e73 eK3 a=
#define e63 ,eQ,c92
#define e53 eO|lT1)
#define e43 public:
#define e33 Rehash(
#define e23 switch(
#define e13 pclone
#define e03 cOr,l6
#define cZ3 newpow
#define cY3 change
#define cX3 (Value_t
#define cW3 (count
#define cV3 133,2,
#define cU3 ,tree
#define cT3 byteCode
#define cS3 child)
#define cR3 (p1 x41)
#define cQ3 cLog2by
#define cP3 long iD1
#define cO3 factor_t
#define cN3 (p0 yI&&
#define cM3 value1
#define cL3 a));if(!
#define cK3 fp_mod(
#define cJ3 else{if(
#define cI3 )min.set(
#define cH3 xF());nD
#define cG3 max.known
#define cF3 known&&
#define cE3 tX p2;p2
#define cD3 {tree.x4
#define cC3 cAbsNot
#define cB3 ,x53 l8
#define cA3 e23 lY1
#define c93 xW 0)nC
#define c83 IsLogicalValue(xW
#define c73 }switch
#define c63 stackptr
#define c53 cLog);xG
#define c43 l8 0));
#define c33 opcodes
#define c23 did_muli
#define c13 &Value){
#define c03 yD const
#define yZ3 (param.
#define yY3 ){e23
#define yX3 :if(&*lE1){
#define yW3 :{n31 r=
#define yV3 xY3 e23
#define yU3 param=*(
#define yT3 cAbsIf,
#define yS3 cNotNot,
#define yR3 l4 16,1,
#define yQ3 yO3 info
#define yP3 lE1=r.specs;if(r.found){
#define yO3 *x8)[a].
#define yN3 cLess,cH
#define yM3 default_function_handling
#define yL3 l4 20,1,
#define yK3 l4 4,1,
#define yJ3 450998,
#define yI3 cExp2 nS
#define yH3 lJ 2},0,
#define yG3 )lS 3*3*2*2
#define yF3 default:
#define yE3 range<nV
#define yD3 range xF
#define yC3 Ge0Lt1
#define yB3 Gt0Le1
#define yA3 cAdd lV2
#define y93 if(op1==
#define y83 iterator
#define y73 begin();
#define y63 TreeSet
#define y53 parent
#define y43 insert(i
#define y33 newrel
#define y23 y03 xX3 xO2
#define y13 y03 iF2;if(
#define y03 ))return
#define xZ3 ;if(half
#define xY3 break;}
#define xX3 IsNever
#define xW3 e62 eO1
#define xV3 synth.xH
#define xU3 b_needed
#define xT3 cachepos
#define xS3 ,lE1,info
#define xR3 half=
#define xQ3 131,4,1,
#define xP3 131,8,1,
#define xO3 ,iM,1,lC1+1);
#define xN3 4,1,2,1,
#define xM3 ::vector
#define xL3 1 y5 n41
#define xK3 FindPos(
#define xJ3 src_pos
#define xI3 xM1 xQ+
#define xH3 reserve(
#define xG3 treeptr
#define xF3 .resize(
#define xE3 tO1 void
#define xD3 ImmedTag
#define xC3 a,const
#define xB3 RefCount
#define xA3 Birth();
#define x93 typename
#define x83 unsigned
#define x73 template
#define x63 7168,
#define x53 leaf1
#define x43 cost_t
#define x33 fpdata
#define x23 middle
#define x13 sqrt_cost
#define x03 const int
#define nZ3 mul_count
#define nY3 );nZ l5::
#define nX3 e62 2)));
#define nW3 maxValue1
#define nV3 minValue1
#define nU3 maxValue0
#define nT3 ValueType
#define nS3 result
#define nR3 nS3;}case
#define nQ3 nS3.tL
#define nP3 ::res,b8<
#define nO3 eG3 yM
#define nN3 .what nS1
#define nM3 e62 1),e62 1));
#define nL3 .max.val
#define nK3 nS3 yM2
#define nJ3 nS3 e61
#define nI3 nS3 yI
#define nH3 nS3 nL3
#define nG3 nS3 tG1
#define nF3 nS3 yM
#define nE3 xF nS3
#define nD3 tO n3 0),
#define nC3 i8);nD lC
#define nB3 abs_mul
#define nA3 xM3<x83>&e71
#define n93 l8 a));
#define n83 pos_set
#define n73 goto e1
#define n63 p1.lT2 p1
#define n53 [funcno].
#define n43 eN1[++IP]
#define n33 sim.x2 1,
#define n23 {sim.Eat(
#define n13 eN1[IP]==
#define n03 subtree
#define lZ3 invtree
#define lY3 MakeHash(
#define lX3 parampair
#define lW3 rulenumit
#define lV3 cAnd l3
#define lU3 ,cMul l3
#define lT3 cAnd,l6
#define lS3 x9 lT 2,
#define lR3 },{{1,
#define lQ3 cEqual,
#define lP3 lQ3 lA
#define lO3 t01},{{3,
#define lN3 MakeEqual
#define lM3 nC1,l5::
#define lL3 nC1,{l5::
#define lK3 newbase
#define lJ3 fp_equal(
#define lI3 branch1op
#define lH3 branch2op
#define lG3 if(lJ3
#define lF3 l8 a)xI
#define lE3 overlap
#define lD3 truth_b
#define lC3 truth_a
#define lB3 found_dup
#define lA3 nQ r;r iH
#define l93 void set(
#define l83 {tP1 lE
#define l73 rangeutil
#define l63 Plan_Has(
#define l53 StackMax)
#define l43 i1 true;}
#define l33 namespace
#define l23 (cond yX
#define l13 inverted
#define l03 xX3:
#define iZ2 iftree
#define iY2 depcodes
#define iX2 explicit
#define iW2 cCosh nS
#define iV2 t01 nH
#define iU2 VarBegin
#define iT2 t63 a]
#define iS2 iR2 size()
#define iR2 Params.
#define iQ2 ].data);
#define iP2 i8)));nZ
#define iO2 yQ1.SubTrees
#define iN2 yQ1.Others
#define iM2 );synth
#define iL2 ;DumpTree(
#define iK2 ;Value_t
#define iJ2 begin(),
#define iI2 cond_add
#define iH2 cond_mul
#define iG2 cond_and
#define iF2 IsAlways
#define iE2 func lQ1
#define iD2 bool eM1
#define iC2 Forget()
#define iB2 .second);
#define iA2 Optimize()
#define i92 costree
#define i82 sintree
#define i72 leaf_count
#define i62 &&cond eD)
#define i52 .tT1 n]
#define i42 =GetParam(
#define i32 sub_params
#define i22 nC==cLog2&&
#define i12 nC==cPow&&
#define i02 printf(
#define tZ2 cbrt_count
#define tY2 sqrt_count
#define tX2 cPow i7
#define tW2 ,cPow,
#define tV2 ,cGreater
#define tU2 exponent);
#define tT2 Finite
#define tS2 min.n3 0),
#define tR2 p1 cS ifp1
#define tQ2 yR 2,cAdd)
#define tP2 pcall_tree
#define tO2 after_powi
#define tN2 GetHash().
#define tM2 yP t23
#define tL2 params)
#define tK2 grammar
#define tJ2 cEqual t11
#define tI2 cLog nS
#define tH2 cNeg,lT 1,
#define tG2 ),0},{
#define tF2 std::move(
#define tE2 iH cond nC
#define tD2 tree iH
#define tC2 ){eL3
#define tB2 tree))cJ
#define tA2 );t0=!t0;}
#define t92 tmp c91 tree
#define t82 nQ tmp;tmp iH
#define t72 tree nC
#define t62 MakeNEqual
#define t52 )?0:1))l7
#define t42 isInteger(
#define t32 Comparison
#define t22 needs_flip
#define t12 (half&63)-1;
#define t02 value]
#define eZ2 lT1 opcode
#define eY2 )lS 3*3*3*2
#define eX2 cS tree);
#define eW2 mul_item
#define eV2 innersub
#define eU2 cbrt_cost
#define eT2 best_cost
#define eS2 condition
#define eR2 nominator
#define eQ2 per_item
#define eP2 item_type
#define eO2 first2
#define eN2 l4 18,1,
#define eM2 cIf,iX 3,
#define eL2 lJ 1},0,
#define eK2 tJ 1},0,
#define eJ2 Decision
#define eI2 not_tree
#define eH2 (mulgroup
#define eG2 (lR));nD lC
#define eF2 Become(xW
#define eE2 group_by
#define eD2 exponent=
#define eC2 ->second
#define eB2 targetpos
#define eA2 ParamSpec
#define e92 rhs.hash2;}
#define e82 rhs.hash1
#define e72 struct
#define e62 Value_t(
#define e52 .n_int_sqrt
#define e42 const std::eP
#define e32 const char*
#define e22 nT 409641,
#define e12 ,xF1);lC
#define e02 );xY3
#define cZ2 if(t72==
#define cY2 eO3 bool>
#define cX2 ,(long double)
#define cW2 ContainsOtherCandidates
#define cV2 std::cout
#define cU2 source_tree
#define cT2 GetParam eS
#define cS2 <tP,x43>
#define cR2 p1_evenness
#define cQ2 isNegative(
#define cP2 neg_set
#define cO2 }else{x8=new
#define cN2 );}void
#define cM2 cNop,cNop}}
#define cL2 cTanh,cNop,
#define cK2 NewHash
#define cJ2 >e72 cB<
#define cI2 matches
#define cH2 {cV2<<
#define cG2 iS1 void*)&
#define cF2 cGreater,cH
#define cE2 cSin i7
#define cD2 cCos nS
#define cC2 ,t21 0x1 nH
#define cB2 +=1 i1 n91;
#define cA2 negated
#define c92 synth);
#define c82 Specializer
#define c72 ifdata.ofs
#define c62 (IfData&ifdata
#define c52 .push_back(
#define c42 ;}data;data.
#define c32 );sim.x2 2,
#define c22 nP)l14 tmp=
#define c12 (*x8)[0].info
#define c02 CodeTree
#define yZ2 c02 xF
#define yY2 coshtree
#define yX2 sinhtree
#define yW2 best_score
#define yV2 mulvalue
#define yU2 pow_item
#define yT2 subgroup
#define yS2 PowiResult
#define yR2 .match_tree
#define yQ2 )l43
#define yP2 0));yD3
#define yO2 maxValue
#define yN2 minValue
#define yM2 yI eW3 if(
#define yL2 fp_min(yL,
#define yK2 div_tree
#define yJ2 pow_tree
#define yI2 preserve
#define yH2 ,cCos i7
#define yG2 (rule cU3,info
#define yF2 e62 0.5)
#define yE2 PullResult()
#define yD2 dup_or_fetch
#define yC2 e33 false
#define yB2 test_order
#define yA2 lX3,
#define y92 .param_count
#define y82 shift(index)
#define y72 rulenumber
#define y62 cTanh nS
#define y52 cSinh nS
#define y42 cInv,lT 1,
#define y32 constraints=
#define y22 factor_immed
#define y12 changes
#define y02 n81 cS y4 l8
#define xZ2 cS leaf2 l8
#define xY2 cS x53 l8
#define xX2 cS cond l8
#define xW2 exp_diff
#define xV2 ExponentInfo
#define xU2 lower_bound(
#define xT2 factor
#define xS2 is_logical
#define xR2 newrel_and
#define xQ2 tH[c eE
#define xP2 ;iM.Remember(
#define xO2 i1 Unknown;}
#define xN2 res_stackpos
#define xM2 half_pos
#define xL2 >>1)):(
#define xK2 CodeTreeData
#define xJ2 multiply(
#define xI2 tO known)
#define xH2 var_trees
#define xG2 erase(cs_it);
#define xF2 parent_opcode
#define xE2 log2_exponent
#define xD2 yB swap(tmp);
#define xC2 Value(Value::
#define xB2 dup_fetch_pos
#define xA2 a;if(eK1){x8=
#define x92 {cZ start_at;
#define x82 xX3 cQ lC
#define x72 cSin nS
#define x62 Value_EvenInt
#define x52 )){data xC
#define x42 MakeFalse,{l5
#define x32 if(list.first
#define x22 AddCollection
#define x12 ConditionType
#define x02 DUP_ONE(apos)
#define nZ2 SpecialOpcode
#define nY2 =i eC2.
#define nX2 IsDefined()
#define nW2 fp_max(yL);
#define nV2 (tree,cV2
#define nU2 e62-1)
#define nT2 assimilated
#define nS2 denominator
#define nR2 fraction
#define nQ2 l2 18,2,
#define nP2 .GetDepth()
#define nO2 iH t72)
#define nN2 xI leaf2 l8
#define nM2 DUP_BOTH();
#define nL2 x73 lL
#define nK2 -1-offset].
#define nJ2 tree.GetHash()
#define nI2 TreeCounts
#define nH2 ,e62 1))){
#define nG2 bool t0 eW3
#define nF2 found_log2
#define nE2 div_params
#define nD2 immed_sum
#define nC2 :sim.Eat(1,
#define nB2 OPCODE(opcode)
#define nA2 ;sim.Push(
#define n92 break;nS3*=
#define n82 FactorStack xF
#define n72 iF2 cQ lC
#define n62 cLessOrEq,
#define n52 282870 nT
#define n42 cNotNot nS
#define n32 cNot nS
#define n22 replacing_slot
#define n12 RefParams
#define n02 if_always[
#define lZ2 WhatDoWhenCase
#define lY2 exponent_immed
#define lX2 new_base_immed
#define lW2 base_immed
#define lV2 ||op1==
#define lU2 remaining
#define lT2 Rehash t8
#define lS2 data[a eI3
#define lR2 lT2 r);}
#define lQ2 if(newrel_or==
#define lP2 .UseGetNeeded(
#define lO2 e7 2,131,
#define lN2 Immed eJ3);
#define lM2 OptimizedUsing
#define lL2 Var_or_Funcno
#define lK2 lL2;
#define lJ2 GetParams(
#define lI2 crc32_t
#define lH2 signed_chain
#define lG2 MinusInf
#define lF2 synth.Find(
#define lE2 );cV2<<
#define lD2 return true;
#define lC2 n_immeds
#define lB2 stack eJ3
#define lA2 FindClone(xM
#define l92 GetOpcode())
#define l82 needs_rehash
#define l72 AnyWhere_Rec
#define l62 minimum_need
#define l52 ~x83(0)
#define l42 41,42,43,44,
#define l32 p1_logical_b
#define l22 p0_logical_b
#define l12 p1_logical_a
#define l02 p0_logical_a
#define iZ1 (size_t
#define iY1 )e73
#define iX1 TopLevel)
#define iW1 .e33)
#define iV1 {pow.CopyOnWrite
#define iU1 nB OPCODE
#define iT1 const yZ2
#define iS1 (const
#define iR1 iS1 yZ2&
#define iQ1 ,PowiCache&iM,
#define iP1 else if(
#define iO1 iP1!nS3
#define iN1 synth.DoDup(
#define iM1 cache_needed
#define iL1 e7 2,1,e7 2,
#define iK1 [c72+
#define iJ1 treelist
#define iI1 IsDescendantOf(
#define iH1 has_bad_balance
#define iG1 )tN mulgroup)
#define iF1 .SetParamsMove(
#define iE1 cO3 xT2
#define iD1 double)exponent
#define iC1 {if(GetOpcode()
#define iB1 cV2<<std::y11
#define iA1 cV2<<"POP "
#define i91 DelParam(
#define i81 set(fp_ceil);tK
#define i71 fp_abs(max.val))
#define i61 fp_abs(min.val)
#define i51 cNEqual
#define i41 },0,0x0},{{
#define i31 tJ 2 i41
#define i21 Oneness_NotOne|
#define i11 Value_IsInteger
#define i01 Constness_Const
#define tZ1 DumpHashesFrom(
#define tY1 reltype
#define tX1 const Value_t&i)
#define tW1 const c02&
#define tV1 const Value_t&v
#define tU1 SequenceOpcodes
#define tT1 sep_list[
#define tS1 ;else_tree
#define tR1 goto fail;}
#define tQ1 l1 0x4 nH
#define tP1 x73<
#define tO1 n12);
#define tN1 ){case iF2:
#define tM1 grammar_rules[*r]
#define tL1 x73 set_if<
#define tK1 x73 lX
#define tJ1 TreeCountType xF
#define tI1 GetParamCount(nX
#define tH1 >(e62 1),
#define tG1 nL3)
#define tF1 )lX3.second
#define tE1 e62 0.0)){nU
#define tD1 .IsImmed()
#define tC1 a)tD1)
#define tB1 nB2);
#define tA1 stack[lB2-
#define t91 stack c52
#define t81 )eX3 iW1;
#define t71 synth.PushImmed(
#define t61 MaxChildDepth
#define t51 repl_param_list,
#define t41 const Rule&rule,
#define t31 std::pair<It,It>
#define t21 cPow,lA
#define t11 ,l0 2,
#define t01 ,l1 0x12
#define eZ1 Sign_Negative
#define eY1 Value_Logical
#define eX1 yZ2&b
#define eW1 new_factor_immed
#define eV1 occurance_pos
#define eU1 exponent_hash
#define eT1 exponent_list
#define eS1 CollectionSet xF
#define eR1 CollectMulGroup(
#define eQ1 source_set
#define eP1 exponent,y63
#define eO1 *const func)(
#define eN1 ByteCode
#define eM1 operator
#define eL1 FindAndDup(tree);
#define eK1 &*start_at
#define eJ1 <<nJ2.
#define eI1 xW 1)tD1&&
#define eH1 retry_anyparams_3
#define eG1 retry_anyparams_2
#define eF1 e6(),eO3
#define eE1 needlist_cached_t
#define eD1 yH3 0x4 lR3
#define eC1 eL2 0x4 lR3
#define eB1 ),lM2(
#define eA1 CodeTreeImmed xF(
#define e91 GetParamCount()==
#define e81 .back().thenbranch
#define e71 eN1,size_t&IP,size_t limit,size_t y2
#define e61 .cG3
#define e51 (lR,xW i8));nD
#define e41 ;flipped=!flipped;}
#define e31 ,l1 0x0},{{3,
#define e21 }xY3 case
#define e11 for iZ1 b=0;b<
#define e01 by_float_exponent
#define cZ1 lJ3 exponent
#define cY1 new_exp
#define cX1 end()&&i->first==
#define cW1 return BecomeZero;
#define cV1 =comp.AddItem(atree
#define cU1 return BecomeOne;
#define cT1 if(lQ eJ3<=n1)
#define cS1 addgroup
#define cR1 found_log2by
#define cQ1 nC==cC3)
#define cP1 ParsePowiMuli(
#define cO1 lL2)
#define cN1 eL 529654 nT
#define cM1 branch1_backup
#define cL1 branch2_backup
#define cK1 exponent_map
#define cJ1 plain_set
#define cI1 rangehalf
#define cH1 LightWeight(
#define cG1 xV3 1
#define cF1 divgroup
#define cE1 ,iM e63
#define cD1 if(value
#define cC1 tK1 c0
#define cB1 tK1 static
#define cA1 mulgroup.
#define c91 .AddParamMove(
#define c81 yB AddParamMove(
#define c71 ;n81 cY op1 yB DelParams
#define c61 should_regenerate=true;
#define c51 should_regenerate,
#define c41 Collection
#define c31 RelationshipResult
#define c21 Subdivide_Combine(
#define c11 long value
#define c01 )const yP
#define yZ1 rhs c01 hash1
#define yY1 for iZ1 a xY
#define yX1 best_sep_factor
#define yW1 SynthesizeParam
#define yV1 needlist_cached
#define yU1 inline x83
#define yT1 opcode,bool pad
#define yS1 changed=true;
#define yR1 );xM iF1
#define yQ1 NeedList
#define yP1 tK1 bool
#define yO1 ;tK1
#define yN1 i91 a);}
#define yM1 MakesInteger(
#define yL1 const Value_t&value
#define yK1 best_sep_cost
#define yJ1 MultiplicationRange
#define yI1 .min.set(fp_floor);
#define yH1 pihalf_limits
#define yG1 yR 2,cMul);lC
#define yF1 n_stacked
#define yE1 cK2.hash1
#define yD1 AnyParams_Rec
#define yC1 ;synth.StackTopIs(
#define yB1 synth.AddOperation(
#define yA1 continue;
#define y91 Become(value l8 0))
#define y81 }inline
#define y71 PositionalParams,0}
#define y61 always_sincostan
#define y51 Recheck_RefCount_Div
#define y41 Recheck_RefCount_Mul
#define y31 MultiplyAndMakeLong(
#define y21 covers_plus1
#define y11 endl;DumpHashes(
#define y01 if(synth.FindAndDup(
#define xZ1 grammar_func
#define xY1 tree tD1 cQ
#define xX1 cOr l3 16,1,
#define xW1 252421 nT 24830,
#define xV1 l2 0,2,165888 nT
#define xU1 Modulo_Radians},
#define xT1 yB SetParam(
#define xS1 PositionType
#define xR1 CollectionResult
#define xQ1 const_offset
#define xP1 inline TriTruthValue
#define xO1 stacktop_desired
#define xN1 iU3);lC
#define xM1 SetStackTop(
#define xL1 tK1 void
#define xK1 FPoptimizer_ByteCode
#define xJ1 1)?(poly^(
#define xI1 eP3&synth)
#define xH1 e62 0)
#define xG1 xH1)
#define xF1 cond_type
#define xE1 fphash_value_t
#define xD1 Recheck_RefCount_RDiv
#define xC1 cMul);tmp tY3 tmp.
#define xB1 SwapLastTwoInStack();
#define xA1 SetParams(lJ2)
#define x91 fPExponentIsTooLarge(
#define x81 CollectMulGroup_Item(
#define x71 pair<Value_t,y63>
#define x61 (x53 l8 1)nN2
#define x51 for iZ1 a=0;a<yC++a)
#define x41 .GetImmed()
#define x31 {int mStackPtr=0;
#define x21 nL xM1 xQ-1);
#define x11 covers_full_cycle
#define x01 AssembleSequence(
#define nZ1 inverse_nominator
#define nY1 252180 nT 281854,
#define nX1 <<std::dec<<")";}
#define nW1 {DataP slot_holder(y1[
#define nV1 },{l5::MakeNotP1,l5::
#define nU1 },{l5::MakeNotP0,l5::
#define nT1 },{l5::MakeNotNotP1,l5::
#define nS1 !=xJ)if(TestCase(
#define nR1 AddFunctionOpcode
#define nQ1 public e6,public eO3
#define nP1 yZ2&tree,
#define nO1 std::pair<T1,T2>&
#define nN1 tP1 x93
#define nM1 has_good_balance_found
#define nL1 n_occurrences
#define nK1 found_log2_on_exponent
#define nJ1 covers_minus1
#define nI1 needs_resynth
#define nH1 immed_product
#define nG1 l33 FPoptimizer_Optimize
#define nF1 (ParamSpec_Extract xF(
#define nE1 yV3 bitmask&
#define nD1 Sign_Positive
#define nC1 ::MakeTrue
#define nB1 SetParamMove(
#define nA1 CodeTreeImmed(e62
#define n91 Suboptimal
#define n81 changed_if
#define n71 n_as_tanh_param
#define n61 opposite=
#define n51 xE1(
#define n41 eN1 eJ3
#define n31 MatchResultType
#define n21 yP yZ2(
#define n11 needs_sincos
#define n01 resulting_exponent
#define lZ1 Unknown:yF3;}
#define lY1 GetLogicalValue(xW
#define lX1 GetParam(a)
#define lW1 cMul l3 0,1,
#define lV1 IsImmed())l14
#define lU1 void FunctionParserBase
#define lT1 (x83
#define lS1 lT1 a=0;a<xT;++a)
#define lR1 o<<"("<<std::hex<<data.
#define lQ1 (val);else*this=model;}
#define lP1 IfBalanceGood(
#define lO1 n_as_tan_param
#define lN1 changed_exponent
#define lM1 inverse_denominator
#define lL1 ;cK2.hash2+=
#define lK1 retry_positionalparams_2
#define lJ1 x83 index
#define lI1 situation_flags&
#define lH1 518 nT 400412,
#define lG1 data.subfunc_opcode
#define lF1 },{l5::MakeNotNotP0,l5::
#define lE1 (yO3 start_at
#define lD1 CopyOnWrite();
#define lC1 recursioncount
#define lB1 PlanNtimesCache(
#define lA1 FPoptimizer_Grammar
#define l91 GetParamCount();
#define l81 GetPositivityInfo iS3
#define l71 ParamSpec_SubFunctionData
#define l61 t43<e62
#define l51 (tree.GetParamCount()
#define l41 iZ1 a=yC a-->0;)
#define l31 PositionalParams_Rec
#define l21 DumpTreeWithIndent(*this);
#define l11 e23 type){case cond_or:
#define l01 tP1 x83 Compare>
#define iZ tree.i91 a
#define iY .l91 a-->0;)if(
#define iX lA 0x4},{{
#define iW lJ 2 i41 1,
#define iV edited_powgroup
#define iU has_unknown_max
#define iT has_unknown_min
#define iS static const yD3
#define iR if(keep_powi
#define iQ synthed_tree
#define iP 7168 nT 401798,
#define iO SelectedParams,0 i41
#define iN collections
#define iM cache
#define iL ,cIf,l0 3,
#define iK ,eN1,IP,limit,y2,stack);
#define iJ by_exponent
#define iI );p2.lT2 p2 yB eL3 iZ2 nC);cJ}
#define iH .eL3
#define iG mulgroup;mulgroup iH
#define iF (p0 e61&&p0 nL3 i2
#define iE cN3 p0 yM>=e62 0.0))
#define iD goto ReplaceTreeWithOne t83
#define iC !=xJ)return n02
#define iB e01.data
#define iA iX2 xK2(
#define i9 needs_sinhcosh
#define i8 1)x41
#define i7 ,l4 2,1,
#define i6 cAdd l3 0,
#define i5 tG2 e62
#define i4 tK1 n9
#define i3 MakeFalse,l5::
#define i2 <=fp_const_negativezero xF())
#define i1 ;return
#define i0 )i1 lH
#define tZ CalculateResultBoundaries(xW
#define tY p0=CalculateResultBoundaries(
#define tX ;yZ2
#define tW ,2,1 iM2.xR if(found[data.
#define tV 408964 nT 24963,
#define tU 528503 nT 24713,
#define tT matched_params
#define tS [n1 eH3=true;lQ[n1 eI3
#define tR lA1::Grammar*
#define tQ AddOperation(cInv,1,1 iM2.xR}
#define tP int_exponent_t
#define tO m.max.
#define tN ;AddParamMove(
#define tM nM nU2,e62 1));
#define tL cG3 eW3
#define tK return m;}case
#define tJ x9 AnyParams,
#define tI powgroup l8
#define tH relationships
#define tG ]!=~size_t(0)){synth.yT
#define tF }},{ProduceNewTree,2,1,
#define tE nA1(
#define tD has_mulgroups_remaining
#define tC MatchInfo xF&
#define tB e33);i32 c52
#define tA best_factor
#define t9 RootPowerTable xF::RootPowers[
#define t8 (c81
#define t7 :goto ReplaceTreeWithZero t83
#define t6 MatchPositionSpec_AnyParams xF
#define t5 l33 FPoptimizer_CodeTree
#define t4 n_as_sinh_param
#define t3 n_as_cosh_param
#define t2 ();pow iH cLog yB eL3 cMul
#define t1 0,tU2 i91 1);
#define t0 is_signed
#define eZ result_positivity
#define eY biggest_minimum
#define eX const l71
#define eW 122999 nT 139399,
#define eV 142455 nT 141449,
#define eU cond_tree
#define eT then_tree
#define eS (a);bool needs_cow=GetRefCount()>1;
#define eR yZ2&tree)
#define eQ sequencing
#define eP string FP_GetOpcodeName(
#define eO );eN1 c52 0x80000000u
#define eN {if(needs_cow){lD1 goto
#define eM (lJ2));cA1 e33);
#define eL ,nQ2
#define eK eO3 yZ2>
#define eJ if_stack
#define eI n_as_sin_param
#define eH n_as_cos_param
#define eG PowiResolver::
#define eF cIf,tG3
#define eE ].relationship
#define eD .BalanceGood
#define eC AddParamMove(yT2
#define eB valueType
#define eA back().endif_location
#define e9 xE1 key
#define e8 AddParamMove(mul);
#define e7 130,1,
#define e6 MatchPositionSpecBase
#define e5 iX2 c02(
#define e4 smallest_maximum
#define e3 ]!=~size_t(0)&&found[data.
#define e2 }PACKED_GRAMMAR_ATTRIBUTE;
#define e1 ReplaceTreeWithParam0;
#define e0 factor_needs_rehashing
#define cZ MatchPositionSpecBaseP
#define cY .e33 yB eL3
#define cX x93 tJ1::y83
#define cW ParamSpec_Extract xF(nN.param_list,
#define cV }x32 x41==e62
#define cU 243,244,245,246,249,250,251,253,255,256,257,258,259}};}
#define cT ];};extern"C"{
#define cS .AddParam(
#define cR iL2 tree lE2"\n";
#define cQ )return false;
#define cP 79,122,123,160,161,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,
#define cO 27,28,29,30,31,32,33,35,36,
#define cN const ParamSpec_SubFunction
#define cM const ParamSpec_ParamHolder
#define cL otherhalf
#define cK :{AdoptChildrenWithSameOpcode(tree);
#define cJ goto redo;
#define cI StackState
#define cH l2 16,2,
#define cG const SequenceOpCode xF
#define cF MatchPositionSpec_PositionalParams xF
#define cE const nP1 std::ostream&o
#define cD e62 1.5)*fp_const_pi xF()
#define cC CalculatePowiFactorCost(
#define cB ImmedHashGenerator
#define cA paramholder_matches
#define c9 ::map<fphash_t,std::set<eM3> >
#define c8 ComparisonSetBase::
#define c7 AddParamMove(comp.cJ1[a].value);
#define c6 T1,x93 T2>inline iD2()(
#define c5 has_nonlogical_values
#define c4 from_logical_context)
#define c3 AnyParams,0}},{ProduceNewTree,
#define c2 for iZ1 a=xV.l91 a-->0;)
#define c1 POWI_CACHE_SIZE
#define c0 static inline yZ2
#define yZ ++IP;yA1}if(n13 c33.
#define yY },{l5::xJ,l5::Never},{l5::xJ,l5::Never}}
#define yX .FoundChild
#define yW BalanceResultType
#define yV {yZ2 tmp;tmp iH
#define yU yB1 GetOpcode(),
#define yT DoDup(found[data.
#define yS xB3(0),Opcode(
#define yR ;sim.Eat(
#define yQ );void nR1 eZ2,c82<
#define yP {return
#define yO const yP data->
#define yN +=fp_const_twopi xF();
#define yM .min.val
#define yL fp_sin(min),fp_sin(max))
#define yK fp_const_twopi xF());if(
#define yJ {yZ2 tmp,tmp2;tmp2 iH
#define yI .min.known
#define yH for iZ1 a=0;a<l91++a){if(
#define yG MatchPositionSpec_AnyWhere
#define yF if yZ3 data.match_type==
#define yE void OutFloatHex(std::ostream&o,
#define yD {static void lY3 nB fphash_t&cK2,
#define yC tree.l91
#define yB );tree.
#define yA AddParam(CodeTreeImmed(
#define y9 cGreaterOrEq,
#define y8 ,x93 yZ2::
#define y7 xF model=cI1 xF()){if(known
#define y6 AssembleSequence_Subdivide(
#define y5 ]=0x80000000u|x83(
#define y4 branch2
#define y3 x83 c;x83 short l[
#define y2 factor_stack_base
#define y1 data->Params
#define y0 (lW3 r=range.first;r!=range.i63++r){
#define xZ {nI2.erase(i);yA1}
#define xY =0;a<y53.l91++a)if(
#define xX ,t21 0x4 nH
#define xW tree l8
#define xV branch1
#define xU ,nU2))eN
#define xT nN y92
#define xS =fp_cosh(m yM);tO val=fp_cosh(tO val);
#define xR StackTopIs(*this)i1;}
#define xQ StackTop
#define xP FPOPT_autoptr
#define xO +=nS3 i1 nS3;}tK1 inline Value_t
#define xN int_exponent
#define xM newnode
#define xL lJ 1 i41
#define xK has_highlevel_opcodes
#define xJ Unchanged
#define xI .IsIdenticalTo(
#define xH GetStackTop()-
#define xG sim.AddConst(
#define xF <Value_t>
#define xE cE=cV2
#define xD best_selected_sep
#define xC ->Recalculate_Hash_NoRecursion();}
#define xB l91++a)if(ApplyGrammar(tK2,xW a),
#define xA lU3 2,1,
#define x9 ,cAdd,
#define x8 position
#define x7 x51{yD3
#define x6 ;tmp2 tY3 tmp iH cInv);tmp c91 tmp2)i1
#define x5 eO3 c02>
#define x4 SetParam(0,iZ2 l8 0))tX p1;p1 iH
#define x3 TestImmedConstraints yZ3 constraints cU3)cQ
#define x2 SwapLastTwoInStack()yR
#define x1 n23 1,cInv e02 xG-1 xN1
#define x0 paramholder_index
#define nZ lD2 case
#define nY occurance_counts
#define nX );a-->0;){iT1&powgroup=lX1;if(powgroup
#define nW )){tree.tE3;}
#define nV Value_t>p=tZ a));if(p.
#define nU tree.ReplaceWithImmed(
#define nT ,{2,
#define nS ,l0 1,
#define nR ,cAdd l3 2,1,
#define nQ ){yZ2
#define nP if(xW 0)tD1
#define nO const FPoptimizer_CodeTree::yZ2&tree
#define nN model_tree
#define nM return yD3(
#define nL ){using l33 FUNCTIONPARSERTYPES;
#define nK eK&n12
#define nJ AddParam(xW
#define nI ConstantFolding_LogicCommon(tree,c8
#define nH },{{2,
#define nG nN1 Ref>inline void xP<Ref>::
#define nF AnyParams,1 i41
#define nE ):data(new xK2 xF(
#define nD goto do_return;}
#define nC .GetOpcode()
#define nB FUNCTIONPARSERTYPES::
#define nA b;}};tP1>e72 Comp<nB
#define n9 xK2 xF::xK2(
#define n8 lL2(),t53),Hash(),Depth(1 eB1 0){}
#define n7 SynthesizeByteCode(c92
#define n6 while(ApplyGrammar(cG2
#define n5 GetIntegerInfo(xW 0))==iF2)n73
#define n4 ;tree c91 n81 yQ2
#define n3 tL1 cGreater>(e62
#define n2 DumpParams xF yZ3 data.param_list,iH3 y92,o);
#define n1 restholder_index
#define n0 yZ2 exponent eX3 iH cMul)eX3 cS
#define lZ lR);if(fp_nequal(tmp,xG1){nU e62 1)/tmp);nD}lC
#define lY :if(ParamComparer xF()(t63 1],t63 0])){std::swap(t63 0],t63 1]);Opcode=
#define lX <x93 Value_t>
#define lW xF tmp;tmp iH cPow);tmp tY3 tmp.yA e62
#define lV i01,0x0},
#define lU AddParamMove(pow l8 1));pow.i91 1);pow.e33 yB nB1 0,pow);goto NowWeAreMulGroup;}
#define lT GroupFunction,0},lV{{
#define lS ,e62 1)/e62
#define lR xW 0)x41
#define lQ restholder_matches
#define lP yE1|=key;xE1 crc=(key>>10)|(key<<(64-10))lL1((~n51 crc))*3)^1234567;}};
#define lO n81;n81 nO2;n81 c91 xW 0));n81 cS xV l8
#define lN tK1 yZ2::c02(
#define lM tree.SetParam(0,xW 0)l8 0)xT1 1,CodeTreeImmed(
#define lL lX void eP3::nR1 eZ2,c82<
#define lK cMul,lT 2,
#define lJ cMul,AnyParams,
#define lI (xW 0)tD1&&xW 1)tD1){nU
#define lH CalculateResultBoundaries(tmp);}case
#define lG :cY3=comp.AddRelationship(atree l8 0),atree l8 1),c8
#define lF cPow,l0 2
#define lE x93 Value_t>inline iD2()iS1 Value_t&xC3 Value_t&b)yP a
#define lD {yD3 m=tZ 0));
#define lC break t83
#define lB xL1 yZ2::
#define lA y71,0,
#define l9 l1 0x0 nH
#define l8 .GetParam(
#define l7 tX n81;n81 nO2;n81 iF1 tree.lJ2));n81 cY
#define l6 SelectedParams,0},0,0x0 nH
#define l5 RangeComparisonData
#define l4 y71},{ProduceNewTree,
#define l3 ,AnyParams,0}},{ReplaceParams,
#define l2 y71},{ReplaceParams,
#define l1 cMul,SelectedParams,0},0,
#define l0 lA 0x0},{{
#ifdef _MSC_VER
typedef
x83
int
lI2;
#else
#include <stdint.h>
typedef
uint_least32_t
lI2;
#endif
l33
crc32{enum{startvalue=0xFFFFFFFFUL,poly=0xEDB88320UL}
;tP1
lI2
crc>e72
b8{enum{b1=(crc&xJ1
crc
xL2
crc>>1),b2=(b1&xJ1
b1
xL2
b1>>1),b3=(b2&xJ1
b2
xL2
b2>>1),b4=(b3&xJ1
b3
xL2
b3>>1),b5=(b4&xJ1
b4
xL2
b4>>1),b6=(b5&xJ1
b5
xL2
b5>>1),b7=(b6&xJ1
b6
xL2
b6>>1),res=(b7&xJ1
b7
xL2
b7>>1)}
;}
;inline
lI2
update(lI2
crc,x83
b){
#define B4(n) b8<n>nP3 n+1>nP3 n+2>nP3 n+3>::res
#define R(n) B4(n),B4(n+4),B4(n+8),B4(n+12)
static
const
lI2
table[256]={R(0x00),R(0x10),R(0x20),R(0x30),R(0x40),R(0x50),R(0x60),R(0x70),R(0x80),R(0x90),R(0xA0),R(0xB0),R(0xC0),R(0xD0),R(0xE0),R(0xF0)}
;
#undef R
#undef B4
return((crc>>8))^table[(crc^b)&0xFF];y81
lI2
calc_upd(lI2
c,const
x83
char*buf,size_t
size){lI2
value=c
eK3
p=0;p<size;++p)value=update(value,buf[p])i1
value;y81
lI2
calc
iS1
x83
char*buf,size_t
size)yP
calc_upd(startvalue,buf,size);}
}
#ifndef FPOptimizerAutoPtrHH
#define FPOptimizerAutoPtrHH
nN1
Ref>class
xP{e43
xP():p(0){}
xP(Ref*b):p(b){xA3}
xP
iS1
xP&b):p(b.p){xA3
y81
Ref&eM1*(c01*p;y81
Ref*eM1->(c01
p;}
xP&eM1=(Ref*b){Set(b)i1*this;}
xP&eM1=iS1
xP&b){Set(b.p)i1*this;}
#ifdef __GXX_EXPERIMENTAL_CXX0X__
xP(xP&&b):p(b.p){b.p=0;}
xP&eM1=(xP&&b){if(p!=b.p){iC2;p=b.p;b.p=0;eE3*this;}
#endif
~xP(){Forget(cN2
UnsafeSetP(Ref*newp){p=newp;}
void
swap(xP<Ref>&b){Ref*tmp=p;p=b.p;b.p=tmp;}
private:inline
static
void
Have(Ref*p2);inline
void
iC2;inline
void
xA3
inline
void
Set(Ref*p2);private:Ref*p;}
;nG
iC2{if(!p)return;p->xB3-=1;if(!p->xB3)delete
p;}
nG
Have(Ref*p2){if(p2)++(p2->xB3);}
nG
Birth(){Have(p);}
nG
Set(Ref*p2){Have(p2);iC2;p=p2;}
#endif
#include <utility>
e72
Compare2ndRev{nN1
T>inline
iD2()iS1
T&xC3
T&b
c01
a.second>b.i63}
}
;e72
Compare1st{nN1
c6
const
nO1
xC3
nO1
b
c01
a.first<b.first;}
nN1
c6
const
nO1
a,T1
b
c01
a.first<b;}
nN1
c6
T1
xC3
nO1
b
c01
a<b.first;}
}
;
#ifndef FPoptimizerHashHH
#define FPoptimizerHashHH
#ifdef _MSC_VER
typedef
x83
long
long
xE1;
#define FPHASH_CONST(x) x##ULL
#else
#include <stdint.h>
typedef
uint_fast64_t
xE1;
#define FPHASH_CONST(x) x##ULL
#endif
l33
FUNCTIONPARSERTYPES{e72
fphash_t{xE1
hash1,hash2;fphash_t():hash1(0),hash2(0){}
fphash_t
iS1
xE1&xC3
xE1&b):hash1(a),hash2(b){}
iD2==iS1
fphash_t&yZ1==e82&&hash2==e92
iD2!=iS1
fphash_t&yZ1!=e82||hash2!=e92
iD2<iS1
fphash_t&yZ1!=e82?hash1<e82:hash2<e92}
;}
#endif
#ifndef FPOptimizer_CodeTreeHH
#define FPOptimizer_CodeTreeHH
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
l33
lA1{e72
Grammar;}
l33
xK1{tK1
class
ByteCodeSynth;}
t5{tK1
class
c02
yO1
e72
xK2
yO1
class
c02{typedef
xP<xK2
xF>DataP;DataP
data;e43
c02();~c02();e72
OpcodeTag{}
;e5
iU1
o,OpcodeTag);e72
FuncOpcodeTag{}
;e5
iU1
o,x83
f,FuncOpcodeTag);e72
xD3{}
;e5
tV1,xD3);
#ifdef __GXX_EXPERIMENTAL_CXX0X__
e5
Value_t&&v,xD3);
#endif
e72
VarTag{}
;e5
x83
varno,VarTag);e72
CloneTag{}
;e5
tW1
b,CloneTag);void
GenerateFrom
iS1
x93
FunctionParserBase
xF::Data&data,bool
keep_powi=false);void
GenerateFrom
iS1
x93
FunctionParserBase
xF::Data&data,const
x5&xH2,bool
keep_powi=false);void
SynthesizeByteCode(eO3
x83>&cT3,std
xM3
xF&immed,size_t&stacktop_max);void
SynthesizeByteCode(xK1::eP3&synth,bool
MustPopTemps=true)const;size_t
SynthCommonSubExpressions(xK1::xI1
const;void
SetParams
iS1
x5&xE3
SetParamsMove(x5&tO1
c02
GetUniqueRef();
#ifdef __GXX_EXPERIMENTAL_CXX0X__
void
SetParams(x5&&tO1
#endif
void
SetParam
iZ1
which,tW1
b);void
nB1
size_t
which,c02&b);void
AddParam(tW1
param);void
AddParamMove(c02&param);void
AddParams
iS1
x5&xE3
AddParamsMove(x5&xE3
AddParamsMove(x5&n12,size_t
n22);void
i91
size_t
index);void
DelParams();void
Become(tW1
b);inline
size_t
GetParamCount(c01
lJ2)eJ3;y81
c02&GetParam
iZ1
n)yP
lJ2)[n];y81
tW1
GetParam
iZ1
n
c01
lJ2)[n];y81
void
eL3
iU1
o)iZ3
Opcode=o;y81
iU1
GetOpcode()yO
Opcode;y81
nB
fphash_t
GetHash()yO
Hash;y81
const
x5&lJ2
c01
y1;y81
x5&lJ2)yP
y1;y81
size_t
GetDepth()yO
Depth;y81
const
Value_t&GetImmed()yO
Value;y81
x83
GetVar()yO
lK2
y81
x83
GetFuncNo()yO
lK2
y81
bool
IsDefined(c01
GetOpcode()!=nB
cNop;y81
bool
IsImmed(c01
GetOpcode()==nB
cImmed;y81
bool
IsVar(c01
GetOpcode()==nB
iU2;y81
x83
GetRefCount()yO
xB3;}
void
ReplaceWithImmed(tX1;void
e33
bool
constantfolding=true);void
Sort();inline
void
Mark_Incompletely_Hashed()iZ3
Depth=0;y81
bool
Is_Incompletely_Hashed()yO
Depth==0;y81
const
tR
GetOptimizedUsing()yO
lM2;y81
void
SetOptimizedUsing
iS1
tR
g)iZ3
lM2=g;}
bool
RecreateInversionsAndNegations(bool
prefer_base2=false);void
tE3;void
swap(c02&b){data.swap(b.data);}
bool
IsIdenticalTo(tW1
b)const;void
lD1}
yO1
e72
xK2{int
xB3;iU1
Opcode
iK2
Value;x83
lK2
eK
Params;nB
fphash_t
Hash;size_t
Depth;const
tR
lM2;xK2();xK2
iS1
xK2&b);iA
iU1
o);iA
iU1
o,x83
f);iA
tX1;
#ifdef __GXX_EXPERIMENTAL_CXX0X__
iA
Value_t&&i);xK2(xK2&&b);
#endif
bool
IsIdenticalTo
iS1
xK2&b)const;void
Sort();void
Recalculate_Hash_NoRecursion();private:void
eM1=iS1
xK2&b);}
yO1
c0
CodeTreeImmed(tX1
n21
i
y8
xD3());}
#ifdef __GXX_EXPERIMENTAL_CXX0X__
cC1
CodeTreeImmed
cX3&&i)n21
tF2
i)y8
xD3());}
#endif
cC1
CodeTreeOp(iU1
opcode)n21
opcode
y8
OpcodeTag());}
cC1
CodeTreeFuncOp(iU1
opcode,x83
f)n21
opcode,f
y8
FuncOpcodeTag());}
cC1
CodeTreeVar
lT1
varno)n21
varno
y8
VarTag());}
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
xL1
DumpHashes(xE);xL1
DumpTree(xE);xL1
DumpTreeWithIndent(xE,const
eM3&indent="\\"
);
#endif
}
#endif
#endif
#ifndef FPOptimizer_GrammarHH
#define FPOptimizer_GrammarHH
#include <iostream>
t5{tK1
class
c02;}
l33
lA1{enum
ImmedConstraint_Value{ValueMask=0x07,Value_AnyNum=0x0,x62=0x1,Value_OddInt=0x2,i11=0x3,Value_NonInteger=0x4,eY1=0x5}
;enum
ImmedConstraint_Sign{SignMask=0x18,Sign_AnySign=0x00,nD1=0x08,eZ1=0x10,Sign_NoIdea=0x18}
;enum
ImmedConstraint_Oneness{OnenessMask=0x60,Oneness_Any=0x00,Oneness_One=0x20,Oneness_NotOne=0x40}
;enum
ImmedConstraint_Constness{ConstnessMask=0x180,Constness_Any=0x00,i01=0x80,Constness_NotConst=0x100}
;enum
Modulo_Mode{Modulo_None=0,Modulo_Radians=1}
;enum
Situation_Flags{LogicalContextOnly=0x01,NotForIntegers=0x02,OnlyForIntegers=0x04,OnlyForComplex=0x08,NotForComplex=0x10}
;enum
nZ2{NumConstant,ParamHolder,SubFunction}
;enum
ParamMatchingType{PositionalParams,SelectedParams,AnyParams,GroupFunction}
;enum
RuleType{ProduceNewTree,ReplaceParams}
;
#ifdef __GNUC__
# define PACKED_GRAMMAR_ATTRIBUTE __attribute__((packed))
#else
# define PACKED_GRAMMAR_ATTRIBUTE
#endif
typedef
std::pair<nZ2,const
void*>eA2
yO1
eA2
ParamSpec_Extract
lT1
paramlist,lJ1)yO1
bool
ParamSpec_Compare
iS1
void*xC3
void*b,nZ2
type);x83
ParamSpec_GetDepCode
iS1
eA2&b);e72
ParamSpec_ParamHolder{lJ1:8;x83
constraints:9;x83
depcode:15;e2
tK1
e72
ParamSpec_NumConstant
l14
constvalue;x83
modulo;}
;e72
l71{x83
param_count:2;x83
param_list:30;iU1
subfunc_opcode:8;ParamMatchingType
match_type:3;x83
n1:5;e2
e72
ParamSpec_SubFunction{l71
data;x83
constraints:9;x83
depcode:7;e2
e72
Rule{RuleType
ruletype:2;x83
situation_flags:5;x83
repl_param_count:2+9;x83
repl_param_list:30;l71
match_tree;e2
e72
Grammar{x83
rule_count;x83
short
rule_list[999
cT
extern
const
Rule
grammar_rules[];}
xL1
DumpParam
iS1
eA2&p,std::ostream&o=cV2);xL1
DumpParams
lT1
paramlist,x83
count,std::ostream&o=cV2);}
#endif
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#define CONSTANT_POS_INF HUGE_VAL
#define CONSTANT_NEG_INF (-HUGE_VAL)
l33
FUNCTIONPARSERTYPES{tK1
inline
Value_t
fp_const_pihalf()yP
fp_const_pi
xF()*yF2;}
tK1
inline
Value_t
fp_const_twopi()eQ3
fp_const_pi
xF());nS3
xO
fp_const_twoe()eQ3
fp_const_e
xF());nS3
xO
fp_const_twoeinv()eQ3
fp_const_einv
xF());nS3
xO
fp_const_negativezero()yP-Epsilon
xF::value;}
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
#include <iostream>
nG1{using
l33
lA1;using
t5;using
l33
FUNCTIONPARSERTYPES
yO1
class
MatchInfo{e43
eO3
std::pair<bool,eK> >lQ;eK
cA;eO3
x83>tT;e43
MatchInfo():lQ(),cA(),tT(){}
e43
bool
SaveOrTestRestHolder
lT1
n1,eA3&iJ1){cT1{lQ
xF3
n1+1);lQ
tS=iJ1
l43
if(lQ[n1
eH3==false){lQ
tS=iJ1
l43
eA3&found=lQ[n1
eI3;if(iJ1
eJ3!=found
eJ3
cQ
for
iZ1
a=0;a<iJ1
eJ3;++a)if(!iJ1[a]xI
found[a])cQ
lD2}
void
SaveRestHolder
lT1
n1,eK&iJ1){cT1
lQ
xF3
n1+1);lQ
tS.swap(iJ1);}
bool
SaveOrTestParamHolder
lT1
x0,iT1&xG3){if(cA
eJ3<=x0){cA.xH3
x0+1);cA
xF3
x0);cA
c52
xG3
yQ2
if(!cA[x0].nX2){cA[x0]=xG3
l43
return
xG3
xI
cA[x0]cN2
SaveMatchedParamIndex(lJ1){tT
c52
index);}
iT1&GetParamHolderValueIfFound
lT1
x0)const{static
iT1
dummytree;if(cA
eJ3<=x0)return
dummytree
i1
cA[x0];}
iT1&GetParamHolderValue
lT1
x0
c01
cA[x0];}
bool
HasRestHolder
lT1
n1
c01
lQ
eJ3>n1&&lQ[n1
eH3==true;}
eA3&GetRestHolderValues
lT1
n1)const{static
eA3
empty_result;cT1
return
empty_result
i1
lQ[n1
eI3;}
const
eO3
x83>&GetMatchedParamIndexes(c01
tT;}
void
swap(tC
b){lQ.swap(b.lQ);cA.swap(b.cA);tT.swap(b.tT);}
tC
eM1=iS1
tC
b){lQ=b.lQ;cA=b.cA;tT=b.tT
i1*this;}
}
;class
e6;typedef
xP<e6>cZ;class
e6{e43
int
xB3;e43
e6():xB3(0){}
virtual~e6(){}
}
;e72
n31{bool
found;cZ
specs;n31(bool
f):found(f),specs(){}
n31(bool
f,const
cZ&s):found(f),specs(s){}
}
;xL1
SynthesizeRule(t41
nP1
tC
info)yO1
n31
TestParam
iS1
eA2&yA2
const
nP1
const
cZ&start_at,tC
info)yO1
n31
TestParams(eX&nN,const
nP1
const
cZ&start_at,tC
info,bool
iX1
yO1
bool
ApplyGrammar
iS1
Grammar&tK2,FPoptimizer_CodeTree::nP1
bool
from_logical_context=false);xL1
ApplyGrammars(FPoptimizer_CodeTree::eR
yO1
bool
IsLogisticallyPlausibleParamsMatch(eX&params,const
eR;}
l33
lA1{xL1
DumpMatch(t41
nO,const
FPoptimizer_Optimize::tC
info,bool
DidMatch,std::ostream&o=cV2);xL1
DumpMatch(t41
nO,const
FPoptimizer_Optimize::tC
info,e32
tO3,std::ostream&o=cV2);}
#endif
#include <string>
e42
lA1::nZ2
yT1=false);e42
iU1
yT1=false);
#include <string>
#include <sstream>
#include <assert.h>
#include <iostream>
using
l33
lA1;using
l33
FUNCTIONPARSERTYPES;e42
lA1::nZ2
yT1){
#if 1
e32
p=0;e23
opcode){case
NumConstant:p="NumConstant"
;lC
ParamHolder:p="ParamHolder"
;lC
SubFunction:p="SubFunction"
;xY3
std::ostringstream
tmp;assert(p);tmp<<p;if(pad)while(tmp.str()eJ3<12)tmp<<' '
i1
tmp.str();
#else
std::ostringstream
tmp;tmp<<opcode;if(pad)while(tmp.str()eJ3<5)tmp<<' '
i1
tmp.str();
#endif
}
e42
iU1
yT1){
#if 1
e32
p=0;e23
opcode){case
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
i51:p="cNEqual"
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
cQ3:p="cLog2by"
;lC
cNop:p="cNop"
;break;
#endif
case
cSinCos:p="cSinCos"
;lC
cSinhCosh:p="cSinhCosh"
;lC
cC3:p="cAbsNot"
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
iU2:p="VarBegin"
;xY3
std::ostringstream
tmp;assert(p);tmp<<p;if(pad)while(tmp.str()eJ3<12)tmp<<' '
i1
tmp.str();
#else
std::ostringstream
tmp;tmp<<opcode;if(pad)while(tmp.str()eJ3<5)tmp<<' '
i1
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
;l33
xK1{tK1
class
ByteCodeSynth{e43
ByteCodeSynth():eN1(),Immed(),cI(),xQ(0),StackMax(0){eN1.xH3
64);Immed.xH3
8);cI.xH3
16
cN2
Pull(eO3
x83>&bc,std
xM3
xF&imm,size_t&StackTop_max){for
lT1
a=0;a<n41;++a){eN1[a]&=~0x80000000u;}
eN1.swap(bc);Immed.swap(imm);StackTop_max=StackMax;}
size_t
GetByteCodeSize(c01
n41;}
size_t
GetStackTop(c01
xQ;}
void
PushVar
lT1
varno){eN1
c52
varno);xI3
1
cN2
PushImmed
cX3
immed
nL
eN1
c52
cImmed);Immed
c52
immed);xI3
1
cN2
StackTopIs(nO,int
offset=0){if((int)xQ>offset){cI[xQ
nK2
first=true;cI[xQ
nK2
second=tree;}
}
bool
IsStackTop(nO,int
offset=0
c01(int)xQ>offset&&cI[xQ
nK2
first&&cI[xQ
nK2
second
xI
tree);y81
void
EatNParams
lT1
eat_count){xQ-=eat_count;}
void
ProducedNParams
lT1
produce_count){xI3
produce_count
cN2
DoPopNMov
iZ1
eB2,size_t
srcpos
nL
eN1
c52
cPopNMov
e53
eB2
e53
srcpos);xM1
srcpos+1);cI[eB2]=cI[srcpos];xM1
eB2+1
cN2
DoDup
iZ1
xJ3
nL
if(xJ3==xQ-1){eN1
c52
cDup);}
else{eN1
c52
cFetch
e53
xJ3);}
xI3
1);cI[xQ-1]=cI[xJ3];}
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
tP1
int>void
Dump(){std::ostream&o=cV2;o<<"Stack state now("
<<xQ<<"):\n"
e73
0;a<xQ;++a){o<<a<<": "
;if(cI[a
eH3){nO=cI[a
eI3;o<<'['<<std::hex<<(void*)(&tree.lJ2))<<std::dec<<','<<tree.GetRefCount()<<']'
iL2
tree,o);}
else
o<<"?"
iV3}
o<<std::flush;}
#endif
size_t
xK3
nO)const{for
iZ1
a=xQ;a-->0;)if(cI[a
eH3&&cI[a
eI3
xI
tree
y03
a
i1~size_t(0);}
bool
Find(nO
c01
xK3
tree)!=~size_t(0);}
bool
FindAndDup(nO){size_t
pos=xK3
tree);if(pos!=~size_t(0)){
#ifdef DEBUG_SUBSTITUTIONS
cV2<<iX3"duplicate at ["
<<pos<<"]: "
iL2
tree
lE2" -- issuing cDup or cFetch\n"
;
#endif
DoDup(pos
yQ2
return
t23
e72
IfData{size_t
ofs;}
;void
SynthIfStep1
c62,iU1
op
x21
c72=n41;eN1
c52
op
eO
eO
cN2
SynthIfStep2
c62
x21
eN1
iK1
xL3+2);eN1
iK1
2
y5
lN2
c72=n41;eN1
c52
cJump
eO
eO
cN2
SynthIfStep3
c62
x21
eN1.back()|=0x80000000u;eN1
iK1
xL3-1);eN1
iK1
2
y5
lN2
xI3
1
iY1
0;a<c72;++a){if(eN1[a]==cJump&&eN1[a+1]==(0x80000000u|(c72-1))){eN1[a+xL3-1);eN1[a+2
y5
lN2
c73(eN1[a]){case
cAbsIf:case
cIf:case
cJump:case
cPopNMov:a+=2;lC
cFCall:case
cPCall:case
cFetch:a+=1;break;yF3
xY3}
}
protected:void
xM1
size_t
value){xQ=value;if(xQ>l53{StackMax=xQ;cI
xF3
l53;}
}
protected:eO3
x83>eN1;std
xM3
xF
Immed;eO3
std::pair<bool,FPoptimizer_CodeTree::yZ2> >cI;size_t
xQ;size_t
StackMax;private:void
incStackPtr(){if(xQ+2>l53
cI
xF3
StackMax=xQ+2);}
tP1
bool
IsIntType,bool
IsComplexType>e72
c82{}
;e43
void
AddOperation
eZ2,x83
eat_count,x83
produce_count=1){EatNParams(eat_count);nR1(opcode);ProducedNParams(produce_count
cN2
nR1
eZ2,c82<false,false>yQ
false,true>yQ
true,false>yQ
true,true>);inline
void
nR1
eZ2){nR1(opcode,c82<bool(nB
IsIntType
xF::nS3),bool(nB
IsComplexType
xF::nS3)>());}
}
yO1
e72
SequenceOpCode
yO1
e72
tU1{static
cG
AddSequence;static
cG
MulSequence;}
;xL1
x01
long
count,cG&eQ,xI1;}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
l33
FUNCTIONPARSERTYPES;l33
xK1{tK1
e72
SequenceOpCode
l14
basevalue;x83
op_flip;x83
op_normal,op_normal_flip;x83
op_inverse,op_inverse_flip;}
yO1
cG
tU1
xF::AddSequence={xH1,cNeg
x9
cAdd,cSub,cRSub}
yO1
cG
tU1
xF::MulSequence={e62
1),cInv,cMul,cMul,cDiv,cRDiv}
;
#define findName(a,b,c) "var"
#define TryCompilePowi(o) false
#define mData this
#define mByteCode eN1
#define mImmed Immed
nL2
false,false>)x31
# define FP_FLOAT_VERSION 1
# define FP_COMPLEX_VERSION 0
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
nL2
true,false>)x31
# define FP_FLOAT_VERSION 0
# define FP_COMPLEX_VERSION 0
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
#ifdef FP_SUPPORT_COMPLEX_NUMBERS
nL2
false,true>)x31
# define FP_FLOAT_VERSION 1
# define FP_COMPLEX_VERSION 1
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
nL2
true,true>)x31
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
l33
xK1;
#define POWI_TABLE_SIZE 256
#define POWI_WINDOW_SIZE 3
l33
xK1{
#ifndef FP_GENERATING_POWI_TABLE
extern
const
x83
char
powi_table[POWI_TABLE_SIZE];const
#endif
x83
char
powi_table[POWI_TABLE_SIZE]={0,1,1,1,2,1,2,1,xN3
4,1,2,xP3
2,1,xN3
8,cV3
xQ3
15,1,16,1,2,1,4,1,2,xP3
2,1,4,cV3
1,16,1,25,xQ3
27,5,8,3,2,1,30,1,31,3,32,1,2,1,xN3
8,1,2,xQ3
39,1,16,137,2,1,4,cV3
xP3
45,135,4,31,2,5,32,1,2,131,50,1,51,1,8,3,2,1,54,1,55,3,16,1,57,133,4,137,2,135,60,1,61,3,62,133,63,1,iL1
131,iL1
139,lO2
e7
30,1,130,137,2,31,lO2
e7
e7
130,cV3
1,e7
e7
2,1,130,133,iL1
61,130,133,62,139,130,137,e7
lO2
e7
e7
iL1
131,e7
e7
130,131,2,133,lO2
130,141,e7
130,cV3
1,e7
5,135,e7
lO2
e7
lO2
130,133,130,141,130,131,e7
e7
2,131}
;}
static
x03
c1=256;
#define FPO(x)
l33{class
PowiCache{private:int
iM[c1];int
iM1[c1];e43
PowiCache():iM(),iM1(){iM[1]=1;}
bool
Plan_Add(c11,int
count){cD1>=c1
cQ
iM1[t02+=count
i1
iM[t02!=0;}
void
l63
c11){cD1<c1)iM[t02=1;}
void
Start
iZ1
value1_pos){for(int
n=2;n<c1;++n)iM[n]=-1;Remember(1,value1_pos);DumpContents();}
int
Find(c11)const{cD1<c1){if(iM[t02>=0){FPO(iC3(iG3,"* I found %ld from cache (%u,%d)\n",value,(unsigned)cache[value],iD3 value]))i1
iM[t02;}
eE3-1;}
void
Remember(c11,size_t
iQ3){cD1>=c1)return;FPO(iC3(iG3,"* Remembering that %ld can be found at %u (%d uses remain)\n",value,(unsigned)iQ3,iD3 value]));iM[t02=(int)iQ3;}
void
DumpContents()const{FPO(for(int a=1;a<POWI_CACHE_SIZE;++a)if(cache[a]>=0||iD3 a]>0){iC3(iG3,"== cache: sp=%d, val=%d, needs=%d\n",cache[a],a,iD3 a]);})}
int
UseGetNeeded(c11){cD1>=0&&value<c1)return--iM1[t02
i1
0;}
}
yO1
size_t
y6
long
count
iQ1
cG&eQ,xI1;xL1
c21
size_t
apos,long
aval,size_t
bpos,long
bval
iQ1
x83
cumulation_opcode,x83
cimulation_opcode_flip,xI1;void
lB1
c11
iQ1
int
need_count,int
lC1=0){cD1<1)return;
#ifdef FP_GENERATING_POWI_TABLE
if(lC1>32)throw
false;
#endif
if(iM.Plan_Add(value,need_count
y03;long
xR3
1;cD1<POWI_TABLE_SIZE){xR3
powi_table[t02
xZ3&128){half&=127
xZ3&64)xR3-t12
FPO(iC3(iG3,"value=%ld, half=%ld, otherhalf=%ld\n",value,half,value/half));lB1
half
xO3
iM.l63
half)i1;}
iP1
half&64){xR3-t12}
}
else
cD1&1)xR3
value&((1<<POWI_WINDOW_SIZE)-1);else
xR3
value/2;long
cL=value-half
xZ3>cL||half<0)std::swap(half,cL);FPO(iC3(iG3,"value=%ld, half=%ld, otherhalf=%ld\n",value,half,otherhalf))xZ3==cL){lB1
half,iM,2,lC1+1);}
else{lB1
half
xO3
lB1
cL>0?cL:-cL
xO3}
iM.l63
value);}
tK1
size_t
y6
c11
iQ1
cG&eQ,xI1{int
xT3=iM.Find(value);if(xT3>=0)yP
xT3;}
long
xR3
1;cD1<POWI_TABLE_SIZE){xR3
powi_table[t02
xZ3&128){half&=127
xZ3&64)xR3-t12
FPO(iC3(iG3,"* I want %ld, my plan is %ld * %ld\n",value,half,value/half));size_t
xM2=y6
half
cE1
if(iM
lP2
half)>0||xM2!=cG1){iN1
xM2)xP2
half,cG1);}
x01
value/half
e63
size_t
iQ3=cG1
xP2
value,iQ3);iM.DumpContents()i1
iQ3;}
iP1
half&64){xR3-t12}
}
else
cD1&1)xR3
value&((1<<POWI_WINDOW_SIZE)-1);else
xR3
value/2;long
cL=value-half
xZ3>cL||half<0)std::swap(half,cL);FPO(iC3(iG3,"* I want %ld, my plan is %ld + %ld\n",value,half,value-half))xZ3==cL){size_t
xM2=y6
half
cE1
c21
xM2,half,xM2,half,iM,eQ.op_normal,eQ.op_normal_flip,c92}
else{long
part1=half;long
part2=cL>0?cL:-cL;size_t
part1_pos=y6
part1
cE1
size_t
part2_pos=y6
part2
cE1
FPO(iC3(iG3,"Subdivide(%ld: %ld, %ld)\n",value,half,otherhalf));c21
part1_pos,part1,part2_pos,part2,iM,cL>0?eQ.op_normal:eQ.op_inverse,cL>0?eQ.op_normal_flip:eQ.op_inverse_flip,c92}
size_t
iQ3=cG1
xP2
value,iQ3);iM.DumpContents()i1
iQ3;}
xL1
c21
size_t
apos,long
aval,size_t
bpos,long
bval
iQ1
x83
cumulation_opcode,x83
cumulation_opcode_flip,xI1{int
a_needed=iM
lP2
aval);int
xU3=iM
lP2
bval);bool
flipped
eW3
#define DUP_BOTH() do{if(apos<bpos){size_t tmp=apos;apos=bpos;bpos=tmp e41 FPO(iC3(iG3,"-> " iO3 iO3"op\n",(unsigned)apos,(unsigned)bpos));iN1 apos);iN1 apos==bpos?cG1:bpos);}while(0)
#define DUP_ONE(p) do{FPO(iC3(iG3,"-> " iO3"op\n",(unsigned)p));iN1 p);}while(0)
if(a_needed>0){if(xU3>0){nM2}
cJ3
bpos!=cG1)nM2
else{x02
e41}
}
iP1
xU3>0){if(apos!=cG1)nM2
else
DUP_ONE(bpos);}
cJ3
apos==bpos&&apos==cG1)x02;iP1
apos==cG1&&bpos==xV3
2){FPO(iC3(iG3,"-> op\n"))e41
iP1
apos==xV3
2&&bpos==cG1)FPO(iC3(iG3,"-> op\n"));iP1
apos==cG1)DUP_ONE(bpos);iP1
bpos==cG1){x02
e41
else
nM2}
yB1
flipped?cumulation_opcode_flip:cumulation_opcode,2);}
xL1
cH1
long
count,cG&eQ,xI1{while
cW3<256){int
xR3
xK1::powi_table[count]xZ3&128){half&=127;cH1
half
e63
count/=half;}
else
xY3
if
cW3==1)return;if(!cW3&1)){yB1
cSqr,1);cH1
count/2
e63}
else{iN1
cG1);cH1
count-1
e63
yB1
cMul,2);}
}
}
l33
xK1{xL1
x01
long
count,cG&eQ,xI1{if
cW3==0)t71
eQ.basevalue);else{bool
t22
eW3
if
cW3<0){t22=true;count=-count;}
if(false)cH1
count
e63
iP1
count>1){PowiCache
iM;lB1
count,iM,1);size_t
xO1=synth.GetStackTop();iM.Start(cG1);FPO(iC3(iG3,"Calculating result for %ld...\n",count));size_t
xN2=y6
count
cE1
size_t
n_excess=xV3
xO1;if(n_excess>0||xN2!=xO1-1){synth.DoPopNMov(xO1-1,xN2);}
}
if(t22)yB1
eQ.op_flip,1);}
}
}
#endif
#ifndef FPOptimizer_ValueRangeHH
#define FPOptimizer_ValueRangeHH
t5{l33
l73{l01
e72
Comp{}
;tP1>e72
Comp<nB
cLess>l83<nA
cLessOrEq>l83<=nA
cGreater>l83>nA
cGreaterOrEq>l83>=nA
cEqual>l83==nA
i51>l83!=b;}
}
;}
tK1
e72
cI1
l14
val;bool
known;cI1():val(),known(false){}
cI1(tV1):val(v),known(true){y81
l93
tV1){known=true;val=v;}
l93
xW3
Value_t),cI1
y7)val=iE2
l93
xW3
const
Value_t&),cI1
y7)val=iE2
l01
void
set_if
cX3
v,xW3
Value_t),cI1
y7&&l73::Comp<Compare>()(val,v))val=iE2
l01
void
set_if(tV1,xW3
const
Value_t&),cI1
y7&&l73::Comp<Compare>()(val,v))val=iE2}
yO1
e72
range{cI1
xF
min,max;range():min(),max(){}
range
cX3
mi,Value_t
ma):min(mi),max(ma){}
range(bool,Value_t
ma):min(),max(ma){}
range
cX3
mi,bool):min(mi),max(){}
void
set_abs();void
set_neg();}
yO1
bool
IsLogicalTrueValue
iS1
yD3&p,bool
abs)yO1
bool
IsLogicalFalseValue
iS1
yD3&p,bool
abs);}
#endif
#ifndef FPOptimizer_RangeEstimationHH
#define FPOptimizer_RangeEstimationHH
t5{enum
TriTruthValue{iF2,xX3,Unknown}
yO1
yD3
CalculateResultBoundaries
iS1
eR
yO1
bool
IsLogicalValue
iS1
eR
yO1
TriTruthValue
GetIntegerInfo
iS1
eR
yO1
xP1
GetEvennessInfo
iS1
eR{if(!tree
tD1)return
Unknown;yL1=tree
x41;if(nB
isEvenInteger(value
y13
nB
isOddInteger(value
y23
tK1
xP1
GetPositivityInfo
iS1
eR{yD3
p=CalculateResultBoundaries(tree);if(p
yI&&p
yM>=e62
y13
p
e61
l61
y23
tK1
xP1
GetLogicalValue
iS1
nP1
bool
abs){yD3
p=CalculateResultBoundaries(tree);if(IsLogicalTrueValue(p,abs
y13
IsLogicalFalseValue(p,abs
y23}
#endif
#ifndef FPOptimizer_ConstantFoldingHH
#define FPOptimizer_ConstantFoldingHH
t5{xL1
ConstantFolding(eR;}
#endif
l33{using
l33
FUNCTIONPARSERTYPES;using
t5;e72
ComparisonSetBase{enum{eB3=0x1,Eq_Mask=0x2,Le_Mask=0x3,eC3=0x4,eD3=0x5,Ge_Mask=0x6}
;static
int
Swap_Mask(int
m)yP(m&Eq_Mask)|((m&eB3)?eC3:0)|((m&eC3)?eB3:0);}
enum
c31{Ok,BecomeZero,BecomeOne,n91}
;enum
x12{cond_or,iG2,iH2,iI2}
;}
yO1
e72
ComparisonSet:public
ComparisonSetBase{e72
t32{yZ2
a
tX
b;int
relationship;t32():a(),b(),relationship(){}
}
;eO3
t32>tH;e72
Item{yZ2
value;bool
cA2;Item():value(),cA2(false){}
}
;eO3
Item>cJ1;int
xQ1;ComparisonSet():tH(),cJ1(),xQ1(0){}
c31
AddItem
iR1
a,bool
cA2,x12
type){for
iZ1
c=0;c<cJ1
eJ3;++c)if(cJ1[c].value
xI
a)){if(cA2!=cJ1[c].cA2){l11
cU1
case
iI2:cJ1.erase(cJ1.begin()+c);xQ1
cB2
case
iG2:case
iH2:cW1}
eE3
n91;}
Item
pole;pole.value=a;pole.cA2=cA2;cJ1
c52
pole)i1
Ok;}
c31
AddRelationship(yZ2
a,yZ2
b,int
tY1,x12
type){l11
if(tY1==7)cU1
lC
iI2:if(tY1==7){xQ1
cB2}
lC
iG2:case
iH2:if(tY1==0)cW1
xY3
if(!(a.GetHash()<b.GetHash())){a.swap(b);tY1=Swap_Mask(tY1);}
for
iZ1
c=0;c<tH
eJ3;++c){if(tH[c].a
xI
a)&&tH[c].b
xI
b)){l11{int
y33=xQ2|tY1;if(y33==7)cU1
xQ2=y33;xY3
case
iG2:case
iH2:{int
y33=xQ2&tY1;if(y33==0)cW1
xQ2=y33;xY3
case
iI2:{int
newrel_or=xQ2|tY1;int
xR2=xQ2&tY1;lQ2
5&&xR2==0){xQ2=eD3
i1
n91;}
lQ2
7&&xR2==0){xQ1+=1;tH.erase(tH.begin()+c)i1
n91;}
lQ2
7&&xR2==Eq_Mask){xQ2=Eq_Mask;xQ1
cB2}
yA1}
eE3
n91;}
}
t32
comp;comp.a=a;comp.b=b;comp.relationship=tY1;tH
c52
comp)i1
Ok;}
}
;nN1
Value_t,x93
CondType>bool
ConstantFolding_LogicCommon(nP1
CondType
xF1,bool
xS2){bool
should_regenerate
eW3
ComparisonSet
xF
comp;x51{x93
c8
c31
cY3=c8
Ok;iT1&atree=xW
a);e23
atree
nC){case
cEqual
lG
Eq_Mask
e12
i51
lG
eD3
e12
cLess
lG
eB3
e12
cLessOrEq
lG
Le_Mask
e12
cGreater
lG
eC3
e12
cGreaterOrEq
lG
Ge_Mask
e12
cNot:cY3
cV1
l8
0),true
e12
cNotNot:cY3
cV1
l8
0),false,xF1);break;yF3
if(xS2||IsLogicalValue(atree))cY3
cV1,false,xF1);c73(cY3){ReplaceTreeWithZero:nU
0)i1
true;ReplaceTreeWithOne:nU
1);nZ
c8
Ok:lC
c8
BecomeZero
t7
c8
BecomeOne:iD
c8
n91:c61
xY3}
if(should_regenerate){
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"Before ConstantFolding_LogicCommon: "
cR
#endif
if(xS2){tree.DelParams();}
else{for
l41{iT1&atree=xW
a);if(IsLogicalValue(atree))iZ);}
}
for
iZ1
a=0;a<comp.cJ1
eJ3;++a){if(comp.cJ1[a].cA2
lA3
cNot);r.c7
r.lR2
iP1!xS2
lA3
cNotNot);r.c7
r.lR2
else
tree.c7}
for
iZ1
a=0;a<comp.tH
eJ3;++a
lA3
cNop);e23
comp.tH[a
eE){case
c8
eB3:r
iH
cLess);lC
c8
Eq_Mask:r
iH
cEqual);lC
c8
eC3:r
iH
cGreater);lC
c8
Le_Mask:r
iH
cLessOrEq);lC
c8
eD3:r
iH
i51);lC
c8
Ge_Mask:r
iH
cGreaterOrEq
e02
r
c91
comp.tH[a].a);r
c91
comp.tH[a].b);r.lR2
if(comp.xQ1!=0)tree.yA
e62
comp.xQ1)));
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"After ConstantFolding_LogicCommon: "
cR
#endif
lD2
eE3
t23
yP1
ConstantFolding_AndLogic(iM3(tree.GetOpcode()==cAnd||tree.GetOpcode()==cAbsAnd)i1
nI
iG2,true);}
yP1
ConstantFolding_OrLogic(iM3(tree.GetOpcode()==cOr||tree.GetOpcode()==cAbsOr)i1
nI
cond_or,true);}
yP1
ConstantFolding_AddLogicItems(iM3(tree.GetOpcode()==cAdd)i1
nI
iI2,false);}
yP1
ConstantFolding_MulLogicItems(iM3(tree.GetOpcode()==cMul)i1
nI
iH2,false);}
}
#include <vector>
#include <map>
#include <algorithm>
l33{using
l33
FUNCTIONPARSERTYPES;using
t5;e72
CollectionSetBase{enum
xR1{Ok,n91}
;}
yO1
e72
CollectionSet:public
CollectionSetBase{e72
c41{yZ2
value
tX
xT2;bool
e0;c41():value(),xT2(),e0(false){}
c41
iR1
v,iT1&f):value(v),xT2(f),e0(false){}
}
;std::multimap<fphash_t,c41>iN;typedef
x93
std::multimap<fphash_t,c41>::y83
xS1;CollectionSet():iN(){}
xS1
FindIdenticalValueTo
iR1
value){fphash_t
hash=value.GetHash();for(xS1
i=iN.xU2
hash);i!=iN.cX1
hash;++i){cD1
xI
i
eC2.value
y03
i;eE3
iN.end();}
bool
Found
iS1
xS1&b)yP
b!=iN.end();}
xR1
AddCollectionTo
iR1
xT2,const
xS1&into_which){c41&c=into_which
eC2;if(c.e0)c.xT2
cS
xT2);else{yZ2
add;add
iH
cAdd);add
c91
c.xT2);add
cS
xT2);c.xT2.swap(add);c.e0=true;eE3
n91;}
xR1
x22
iR1
value,iT1&xT2){const
fphash_t
hash=value.GetHash();xS1
i=iN.xU2
hash);for(;i!=iN.cX1
hash;++i){if(i
eC2.value
xI
value
y03
AddCollectionTo(xT2,i);}
iN.y43,std::make_pair(hash,c41(value,xT2)))i1
Ok;}
xR1
x22
iR1
a)yP
x22(a,nA1
1)));}
}
yO1
e72
ConstantExponentCollection{typedef
eK
y63;typedef
std::x71
xV2;eO3
xV2>data;ConstantExponentCollection():data(){}
void
MoveToSet_Unique
iS1
Value_t&eP1&eQ1){data
c52
std::x71(eP1()));data.back().second.swap(eQ1
cN2
MoveToSet_NonUnique
iS1
Value_t&eP1&eQ1){x93
eO3
xV2>::y83
i=std::xU2
data.iJ2
data.end(),exponent,Compare1st());if(i!=data.cX1
exponent){i
eC2.y43
eC2.end(),eQ1.iJ2
eQ1.end());}
else{data.y43,std::x71(exponent,eQ1));}
}
bool
iA2{bool
changed
eW3
std::sort(data.iJ2
data.end(),Compare1st());redo:for
iZ1
a=0;a
eR3
a)l14
exp_a=data[a
eH3;lG3
exp_a,e62
1)))yA1
for
iZ1
b=a+1;b
eR3
b)l14
exp_b=data[b
eH3
iK2
xW2=exp_b-exp_a;if(xW2>=fp_abs(exp_a))break
iK2
exp_diff_still_probable_integer=xW2*e62
16);if(t42
exp_diff_still_probable_integer)&&!(t42
exp_b)&&!t42
xW2))){y63&a_set=lS2;y63&b_set=data[b
eI3;
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"Before ConstantExponentCollection iteration:\n"
;Dump(cV2);
#endif
if(isEvenInteger(exp_b)&&!isEvenInteger(xW2+exp_a)nQ
tmp2;tmp2
iH
cMul);tmp2
iF1
b_set);tmp2
iW1
tX
tmp;tmp
iH
cAbs);tmp
c91
tmp2);tmp
iW1;b_set
xF3
1);b_set[0
tQ3
tmp);}
a_set.insert(a_set.end(),b_set.iJ2
b_set.end());y63
b_copy=b_set;data.erase(data.begin()+b);MoveToSet_NonUnique(xW2,b_copy);yS1
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"After ConstantExponentCollection iteration:\n"
;Dump(cV2);
#endif
cJ}
}
eE3
changed;}
#ifdef DEBUG_SUBSTITUTIONS
void
Dump(std::ostream&out){for
iZ1
a=0;a
eR3
a){out.precision(12);out<<data[a
eH3<<": "
;e11
lS2
eJ3;++b){if(b>0)out<<'*'
iL2
lS2[b],out);}
out<<std::endl;}
}
#endif
}
yO1
static
yZ2
x81
yZ2&value,bool&xK
yY3
value
nC){case
cPow:{yZ2
eD2
value
l8
1);value.y91
i1
exponent;}
case
cRSqrt:value.y91;xK=true
i1
nA1-0.5))t83
cInv:value.y91;xK=true
i1
nA1-1));yF3
xY3
return
nA1
1));}
cB1
void
eR1
eS1&mul,const
nP1
iT1&xT2,bool&c51
bool&xK){for
iZ1
a=0;a<yC++a
nQ
value(xW
a))tX
exponent(x81
value,xK));if(!xT2
tD1||xT2
x41!=e62
1.0)nQ
cY1;cY1
iH
cMul);cY1
cS
tU2
cY1
cS
xT2);cY1
iW1
eX3.swap(cY1);}
#if 0 /* FIXME: This does not work */
cD1
nC==cMul){if(1){bool
exponent_is_even=exponent
tD1&&isEvenInteger(exponent
x41);e11
value.eZ3{bool
tmp=false
tX
val(value
l8
b))tX
exp(x81
val,tmp));if(exponent_is_even||(exp
tD1&&isEvenInteger(exp
x41))nQ
cY1;cY1
iH
cMul);cY1
cS
tU2
cY1
c91
exp);cY1.ConstantFolding();if(!cY1
tD1||!isEvenInteger(cY1
x41)){goto
cannot_adopt_mul;}
}
}
}
eR1
mul,value,exponent,c51
xK);}
else
cannot_adopt_mul:
#endif
{if(mul.x22(value,exponent)==CollectionSetBase::n91)c61}
}
}
yP1
ConstantFolding_MulGrouping(eR{bool
xK
eW3
bool
should_regenerate
eW3
eS1
mul;eR1
mul
cU3,nA1
1)),c51
xK);typedef
std::pair<yZ2,eK>eT1;typedef
std::multimap<fphash_t,eT1>cK1;cK1
iJ;for(x93
eS1::xS1
j=mul.iN.y73
j!=mul.iN.end();++j
nQ&value=j
eC2.value
tX&eD2
j
eC2.xT2;if(j
eC2.e0)exponent
iW1;const
fphash_t
eU1=exponent.GetHash();x93
cK1::y83
i=iJ.xU2
eU1);for(;i!=iJ.cX1
eU1;++i)if(i
eC2.first
xI
exponent)){if(!exponent
tD1||!cZ1
x41,e62
1)))c61
i
eC2.second
c52
value);goto
skip_b;}
iJ.y43,std::make_pair(eU1,std::make_pair(exponent,eK
iZ1(1),value))));skip_b:;}
#ifdef FP_MUL_COMBINE_EXPONENTS
ConstantExponentCollection
xF
e01;for(x93
cK1::y83
j,i=iJ.y73
i!=iJ.end();i=j){j=i;++j;eT1&list=i
eC2;x32.lV1
eD2
list.first
x41;if(!(exponent==xG1)e01.MoveToSet_Unique(exponent,list
iB2
iJ.erase(i);}
}
if(e01.iA2)c61
#endif
if(should_regenerate
nQ
before=tree;before.lD1
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"Before ConstantFolding_MulGrouping: "
iL2
before
lE2"\n"
;
#endif
tree.DelParams();for(x93
cK1::y83
i=iJ.y73
i!=iJ.end();++i){eT1&list=i
eC2;
#ifndef FP_MUL_COMBINE_EXPONENTS
x32.lV1
eD2
list.first
x41;if(exponent==xG1
yA1
if(cZ1
nH2
tree.AddParamsMove(list
iB2
yA1}
}
#endif
yZ2
mul;mul
iH
cMul);mul
iF1
list
iB2
mul
iW1;if(xK&&list.first
tD1){x32
x41==e62
1)/e62
3)nQ
cbrt;cbrt
iH
cCbrt);cbrt.e8
cbrt.lT2
cbrt);yA1
cV
0.5)nQ
sqrt;sqrt
iH
cSqrt);sqrt.e8
sqrt.lT2
sqrt);yA1
cV-0.5)nQ
rsqrt;rsqrt
iH
cRSqrt);rsqrt.e8
rsqrt.lT2
rsqrt);yA1
cV-1)nQ
inv;inv
iH
cInv);inv.e8
inv.lT2
inv);yA1}
}
yZ2
pow;pow
iH
cPow);pow.e8
pow
c91
list.first);pow.lT2
pow);}
#ifdef FP_MUL_COMBINE_EXPONENTS
iJ.clear(iY1
0;a<iB
eJ3;++a)l14
eD2
iB[a
eH3;if(cZ1
nH2
tree.AddParamsMove(iB[a]iB2
yA1}
yZ2
mul;mul
iH
cMul);mul
iF1
iB[a]iB2
mul
iW1
tX
pow;pow
iH
cPow);pow.e8
pow.yA
exponent));pow.lT2
pow);}
#endif
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"After ConstantFolding_MulGrouping: "
cR
#endif
return!tree
xI
before);eE3
t23
yP1
ConstantFolding_AddGrouping(eR{bool
should_regenerate
eW3
eS1
add;x51{if(xW
a)nC==cMul)yA1
if(add.x22(xW
a))==CollectionSetBase::n91)c61}
cY2
lU2
l51);size_t
tD=0;x51{iT1&mulgroup=xW
a);if
eH2
nC==cMul){e11
cA1
eZ3{if
eH2
l8
b)tD1)yA1
x93
eS1::xS1
c=add.FindIdenticalValueTo
eH2
l8
b));if(add.Found(c)nQ
tmp
eH2
y8
CloneTag());tmp.i91
b);tmp
iW1;add.AddCollectionTo(tmp,c);c61
goto
done_a;}
}
lU2[a]=true;tD+=1;done_a:;}
}
if(tD>0){if(tD>1){eO3
std::pair<yZ2,size_t> >nY;std::multimap<fphash_t,size_t>eV1;bool
lB3
eW3
x51
if(lU2[a]){e11
xW
a).eZ3{iT1&p=xW
a)l8
b);const
fphash_t
p_hash=p.GetHash();for(std::multimap<fphash_t,size_t>::const_iterator
i=eV1.xU2
p_hash);i!=eV1.cX1
p_hash;++i){if(nY[i
eC2
eH3
xI
p)){nY[i
eC2
eI3+=1;lB3=true;goto
found_mulgroup_item_dup;}
}
nY
c52
std::make_pair(p,size_t(1)));eV1.insert(std::make_pair(p_hash,nY
eJ3-1));found_mulgroup_item_dup:;}
}
if(lB3
nQ
eE2;{size_t
max=0
eK3
p=0;p<nY
eJ3;++p)if(nY[p
eI3<=1)nY[p
eI3=0;else{nY[p
eI3*=nY[p
eH3
nP2;if(nY[p
eI3>max){eE2=nY[p
eH3;max=nY[p
eI3;}
}
}
yZ2
group_add;group_add
iH
cAdd);
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"Duplicate across some trees: "
iL2
eE2
lE2" in "
cR
#endif
x51
if(lU2[a])e11
xW
a).eZ3
if(eE2
xI
xW
a)l8
b))nQ
tmp(xW
a)y8
CloneTag());tmp.i91
b);tmp
iW1;group_add
c91
tmp);lU2[a]eW3
xY3
group_add
iW1
tX
group;group
iH
cMul);group
c91
eE2);group
c91
group_add);group
iW1;add.x22(group);c61}
}
x51
if(lU2[a]){if(add.x22(xW
a))==CollectionSetBase::n91)c61}
}
if(should_regenerate){
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"Before ConstantFolding_AddGrouping: "
cR
#endif
tree.DelParams();for(x93
eS1::xS1
j=add.iN.y73
j!=add.iN.end();++j
nQ&value=j
eC2.value
tX&coeff=j
eC2.xT2;if(j
eC2.e0)coeff
iW1;if(coeff
tD1){lG3
coeff
x41,xG1)yA1
lG3
coeff
x41
nH2
tree
c91
value);yA1}
}
yZ2
mul;mul
iH
cMul);mul
c91
value);mul
c91
coeff);mul.e33
yB
e8}
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"After ConstantFolding_AddGrouping: "
cR
#endif
lD2
eE3
t23}
l33{using
l33
FUNCTIONPARSERTYPES;using
t5
yO1
bool
ConstantFolding_IfOperations(iM3(tree.GetOpcode()==cIf||tree.GetOpcode()==cAbsIf);for(;;){if(c93==cNot){tD2
cIf);xW
0).eF2
0)c43
xW
1).swap(xW
2));}
iP1
xW
0)cQ1{tD2
t03;xW
0).eF2
0)c43
xW
1).swap(xW
2));}
else
yV3
lY1
0),t72==t03
tN1
tree.eF2
1));nZ
l03
tree.eF2
2));nZ
lZ1
if(c93==cIf||c93==cAbsIf
nQ
cond=xW
0)tX
lC3;lC3
tE2==cIf?cNotNot:cAbsNotNot);lC3
xX2
1));ConstantFolding(lC3)tX
lD3;lD3
tE2==cIf?cNotNot:cAbsNotNot);lD3
xX2
2));ConstantFolding(lD3);if(lC3
tD1||lD3
tD1
nQ
eT;eT
tE2);eT
xX2
1));eT
tX3
eT.nJ
2));eT
iW1
tX
else_tree
tS1
tE2)tS1
xX2
2))tS1.nJ
1))tS1.nJ
2))tS1
cY
cond
nC
xT1
0,cond
l8
0)yB
nB1
1,eT
yB
nB1
2,else_tree
yQ2}
if(xW
1)nC==xW
2)nC&&(xW
1)nC==cIf||xW
1)nC==t03
nQ&x53=xW
1)tX&leaf2=xW
2);if(x53
l8
0)nN2
0))&&x61
1))||x53
l8
2)nN2
2)))nQ
eT;eT
nO2;eT
tY3
eT
xY2
1));eT
xZ2
1));eT
iW1
tX
else_tree
tS1
nO2
tS1.nJ
0))tS1
xY2
2))tS1
xZ2
2))tS1
cY
x53
nC
xT1
0
cB3
0)yB
nB1
1,eT
yB
nB1
2,else_tree
yQ2
if
x61
1))&&x53
l8
2)nN2
2))nQ
eU;eU
nO2;eU
c91
xW
0));eU
xY2
0));eU
xZ2
0));eU
cY
x53
nC
yB
nB1
0,eU
xT1
2
cB3
2)xT1
1
cB3
1)yQ2
if
x61
2))&&x53
l8
2)nN2
1))nQ
eI2;eI2
iH
leaf2
nC==cIf?cNot:cC3);eI2
xZ2
0));eI2
iW1
tX
eU;eU
nO2;eU
c91
xW
0));eU
xY2
0));eU
c91
eI2);eU
cY
x53
nC
yB
nB1
0,eU
xT1
2
cB3
2)xT1
1
cB3
1)yQ2}
yZ2&xV=xW
1)tX&y4=xW
2);if(xV
xI
y4)){tree.eF2
1)yQ2
const
OPCODE
op1=xV
nC;const
OPCODE
op2=y4
nC;y93
op2){if(xV.e91
1
nQ
lO
0));y02
0))c71()n4
if(xV.e91
2&&y4.e91
2){if(xV
l8
0)xI
y4
l8
0))nQ
param0=xV
l8
0)tX
lO
1));y02
1))c71
t8
param0)n4
if(xV
l8
1)xI
y4
l8
1))nQ
param1=xV
l8
1)tX
lO
0));y02
0))c71
t8
n81
c81
param1
yQ2}
y93
yA3
cMul
lV2
cAnd
lV2
cOr
lV2
cAbsAnd
lV2
cAbsOr
lV2
cMin
lV2
cMax){eK
lE3;c2{for
iZ1
b=y4.l91
b-->0;){if(xV
lF3
y4
l8
b))){if(lE3
eU3){xV.t13
lD1}
lE3
c52
xV
n93
y4.i91
b);xV.i91
a
e02}
}
if(!lE3
eU3){xV
iW1;y4
iW1
l7
op1
yB
SetParamsMove(lE3)n4}
}
y93
yA3
cMul||(op1==cAnd&&IsLogicalValue(y4))||(op1==cOr&&IsLogicalValue(y4))){c2
if(xV
lF3
y4)){xV.lD1
xV.i91
a);xV
iW1
tX
cL1=y4;y4=tE
op1==yA3
cOr
t52
op1
c81
cL1)n4}
if((op1==cAnd
lV2
cOr)&&op2==cNotNot
nQ&lH3=y4
l8
0);c2
if(xV
lF3
lH3)){xV.lD1
xV.i91
a);xV
iW1
tX
cL1=lH3;y4=tE
op1==cOr
t52
op1
c81
cL1)n4}
if(op2==cAdd||op2==cMul||(op2==cAnd&&IsLogicalValue(xV))||(op2==cOr&&IsLogicalValue(xV))){for
iZ1
a=y4
iY
y4
lF3
xV)){y4.t13
i91
a);y4
iW1
tX
cM1=xV;xV=tE
op2==cAdd||op2==cOr
t52
op2
c81
cM1)n4}
if((op2==cAnd||op2==cOr)&&op1==cNotNot
nQ&lI3=xV
l8
0
iY1
y4
iY
y4
lF3
lI3)){y4.t13
i91
a);y4
iW1
tX
cM1=lI3;xV=tE
op2==cOr
t52
op2
c81
cM1)n4
eE3
t23}
#include <limits>
l33{using
l33
FUNCTIONPARSERTYPES;using
t5
yO1
int
maxFPExponent()yP
std::numeric_limits
xF::max_exponent;}
yP1
x91
Value_t
base,Value_t
exponent){if(base<xG1
lD2
lG3
base,xG1||lJ3
base,e62
1))cQ
return
exponent>=e62
maxFPExponent
xF())/fp_log2(base);}
yP1
ConstantFolding_PowOperations(iM3(tree.GetOpcode()==cPow);nP&&xW
1).lV1
const_value=t33
lR,xW
i8);nU
const_value)i1
t23
if(eI1
lJ3
xW
i8
nH2
tree.eF2
0)yQ2
nP&&lJ3
lR
nH2
nU
1)i1
t23
nP&&xW
1)nC==cMul){bool
y12=false
iK2
lW2=lR
tX
mulgroup=xW
1
iY1
mulgroup
iY
mulgroup
l8
a).lV1
imm=mulgroup
l8
a)x41;{if(x91
lW2,imm))break
iK2
lX2=t33
lW2,imm);lG3
lX2,xG1)break;if(!y12){y12=true;cA1
lD1}
lW2=lX2;cA1
i91
a
e02}
if(y12){cA1
e33);
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"Before pow-mul change: "
cR
#endif
xW
0).Become(eA1
lW2));xW
1).Become
eH2);
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"After pow-mul change: "
cR
#endif
}
}
if(eI1
c93==cMul)l14
lY2=xW
i8
iK2
y22=1.0;bool
y12=false
tX&mulgroup=xW
0
iY1
mulgroup
iY
mulgroup
l8
a).lV1
imm=mulgroup
l8
a)x41;{if(x91
imm,lY2))break
iK2
eW1=t33
imm,lY2);lG3
eW1,xG1)break;if(!y12){y12=true;cA1
lD1}
y22*=eW1;cA1
i91
a
e02}
if(y12){cA1
e33)tX
cZ3;cZ3
iH
cPow);cZ3
iF1
tree.lJ2));cZ3.yC2
yB
eL3
cMul
c81
cZ3
yB
AddParam(eA1
y22)yQ2}
if(xW
0)i12
eI1
xW
0)l8
1).lV1
a=xW
0)l8
i8
iK2
b=xW
i8
iK2
c=a*b;if(isEvenInteger(a)&&!isEvenInteger(c)nQ
lK3;lK3
iH
cAbs);lK3.nJ
0)c43
lK3.e33
yB
nB1
0,lK3);}
else
tree.SetParam(0,xW
0)l8
0)xT1
1,eA1
c));eE3
t23}
l33{using
l33
FUNCTIONPARSERTYPES;using
t5;e72
l5{enum
eJ2{MakeFalse=0,MakeTrue=1,t62=2,lN3=3,MakeNotNotP0=4,MakeNotNotP1=5,MakeNotP0=6,MakeNotP1=7,xJ=8}
;enum
lZ2{Never=0,Eq0=1,Eq1=2,yB3=3,yC3=4}
;eJ2
if_identical;eJ2
n02
4];e72{eJ2
what:4;lZ2
when:4;}
l02,l12,l22,l32
yO1
eJ2
Analyze
iR1
xC3
eX1)const{if(a
xI
b
y03
if_identical;yD3
tY
a);yD3
p1=CalculateResultBoundaries(b);if(p0
e61&&p1
yI){if(p0
nL3<p1
yM&&n02
0]iC
0];if(p0
nL3<=p1
yM&&n02
1]iC
1];}
if
cN3
p1
e61){if(p0
yM>p1
nL3&&n02
2]iC
2];if(p0
yM>=p1
nL3&&n02
3]iC
3];}
if(IsLogicalValue(a)){if(l02
nN3
l02.when,p1
y03
l02.what;if(l22
nN3
l22.when,p1
y03
l22.what;}
if(IsLogicalValue(b)){if(l12
nN3
l12.when,p0
y03
l12.what;if(l32
nN3
l32.when,p0
y03
l32.what;eE3
xJ;}
cB1
bool
TestCase(lZ2
when,const
yD3&p){if(!p
yI||!p
e61
cQ
e23
when){case
Eq0:nO3==e62
0.0)t43==p
yM
t83
Eq1:nO3==e62
1.0)t43==p
nL3
t83
yB3:nO3>xH1
t43<=e62
1)t83
yC3:nO3>=xH1
l61
1);yF3;eE3
t23}
;l33
RangeComparisonsData{static
const
l5
Data[6]={{l5
lL3
i3
xJ,l5::i3
xJ
lF1
Eq1
nT1
Eq1
nU1
Eq0
nV1
Eq0}
}
,{l5::x42
lM3
xJ,l5
lM3
xJ
lF1
Eq0
nT1
Eq0
nU1
Eq1
nV1
Eq1}
}
,{l5::x42
lM3
t62,l5::i3
MakeFalse
nU1
yB3
nT1
yC3
yY,{l5
lL3
xJ,l5
lM3
i3
lN3
nU1
yC3
nT1
yB3
yY,{l5::x42::i3
i3
MakeTrue,l5::t62
lF1
yC3
nV1
yB3
yY,{l5
lL3
i3
lN3,l5::xJ,l5
nC1
lF1
yB3
nV1
yC3
yY}
;}
yP1
ConstantFolding_Comparison(eR{using
l33
RangeComparisonsData;assert(tree.GetOpcode()>=cEqual&&tree.GetOpcode()<=cGreaterOrEq);e23
Data[t72-cEqual].Analyze(xW
0),xW
1))){case
l5::MakeFalse:nU
0);nZ
l5
nC1:nU
1
nY3
lN3:tD2
cEqual
nY3
t62:tD2
i51
nY3
MakeNotNotP0:tD2
cNotNot
yB
i91
1
nY3
MakeNotNotP1:tD2
cNotNot
yB
i91
0
nY3
MakeNotP0:tD2
cNot
yB
i91
1
nY3
MakeNotP1:tD2
cNot
yB
i91
0
nY3
xJ:;}
if(xW
1)tD1)e23
c93){case
cAsin:lM
fp_sin(xW
iP2
cAcos:lM
fp_cos(xW
i8))yB
eL3
t72==cLess?cGreater:t72==cLessOrEq?cGreaterOrEq:t72==cGreater?cLess:t72==cGreaterOrEq?cLessOrEq:t72);nZ
cAtan:lM
fp_tan(xW
iP2
cLog:lM
fp_exp(xW
iP2
cSinh:lM
fp_asinh(xW
iP2
cTanh:if(fp_less(fp_abs(xW
i8)nH2
lM
fp_atanh(xW
i8))yQ2
break;yF3
xY3
return
t23}
#include <list>
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
using
l33
FUNCTIONPARSERTYPES;l33{
#ifdef DEBUG_SUBSTITUTIONS
yE
double
d){union{double
d;uint_least64_t
h
c42
d=d;lR1
h
nX1
#ifdef FP_SUPPORT_FLOAT_TYPE
yE
float
f){union{float
f;uint_least32_t
h
c42
f=f;lR1
h
nX1
#endif
#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
yE
long
double
ld){union{long
double
ld;e72{uint_least64_t
a;x83
short
b;}
s
c42
ld=ld;lR1
s.b<<data.s.a
nX1
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
yE
long
ld){o<<"("
<<std::hex<<ld
nX1
#endif
#endif
}
t5{lN
nE)){}
lN
const
Value_t&i
y8
xD3
nE
i
x52
#ifdef __GXX_EXPERIMENTAL_CXX0X__
lN
Value_t&&i
y8
xD3
nE
tF2
i)x52
#endif
lN
x83
v
y8
VarTag
nE
iU2,v
x52
lN
iU1
o
y8
OpcodeTag
nE
o
x52
lN
iU1
o,x83
f
y8
FuncOpcodeTag
nE
o,f
x52
lN
const
eX1
y8
CloneTag
nE*b.data)){}
tK1
yZ2::~c02(){}
lB
ReplaceWithImmed(tX1{
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"Replacing "
iL2*this);if(IsImmed())OutFloatHex(cV2,GetImmed()lE2" with const value "
<<i;OutFloatHex(cV2,i
lE2"\n"
;
#endif
data=new
xK2
xF(i);}
tK1
e72
ParamComparer{iD2()iR1
xC3
eX1)const{if(a
nP2!=b
nP2)return
a
nP2<b
nP2
i1
a.GetHash()<b.GetHash();}
}
;xL1
xK2
xF::Sort(yY3
Opcode){case
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
i51:std::sort(iR2
iJ2
iR2
end(),ParamComparer
xF());lC
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
break;yF3
xY3}
lB
AddParam
iR1
param){y1
c52
param);}
lB
AddParamMove(yZ2&param){y1
c52
yZ2());y1.back().swap(param);}
lB
SetParam
iZ1
which,const
eX1)nW1
which
iQ2
y1[which]=b;}
lB
nB1
size_t
which,eX1)nW1
which
iQ2
y1[which
tQ3
b);}
lB
AddParams
iS1
nK){y1.insert(y1.end(),n12.iJ2
n12.end());}
lB
AddParamsMove(nK){size_t
endpos=y1
eJ3,added=n12
eJ3;y1
xF3
endpos+added,yZ2())eK3
p=0;p<added;++p)y1[endpos+p
tQ3
n12[p]);}
lB
AddParamsMove(nK,size_t
n22)nW1
n22
iQ2
i91
n22);AddParamsMove(tO1}
lB
SetParams
iS1
nK){eK
tmp(tO1
y1.swap(tmp);}
lB
SetParamsMove(nK){y1.swap(tO1
n12.clear();}
#ifdef __GXX_EXPERIMENTAL_CXX0X__
lB
SetParams(eK&&n12){SetParamsMove(tO1}
#endif
lB
i91
size_t
index){eK&Params=y1;
#ifdef __GXX_EXPERIMENTAL_CXX0X__
iR2
erase(iR2
begin()+index);
#else
t63
index].data=0
eK3
p=index;p+1<iS2;++p)t63
p].data.UnsafeSetP(&*t63
p+1
iQ2
t63
iS2-1].data.UnsafeSetP(0);iR2
resize(iS2-1);
#endif
}
lB
DelParams(){y1.clear();}
yP1
yZ2::IsIdenticalTo
iS1
eX1)const{if(&*data==&*b.data)return
true
i1
data->IsIdenticalTo(*b.data);}
yP1
xK2
xF::IsIdenticalTo
iS1
xK2
xF&b)const{if(Hash!=b.Hash
cQ
if(Opcode!=b.Opcode
cQ
e23
Opcode){case
cImmed:return
lJ3
Value,b.Value)t83
iU2:return
lL2==b.lK2
case
cFCall:case
cPCall:if(lL2!=b.lL2
cQ
break;yF3
xY3
if(iS2!=b.iS2
cQ
for
iZ1
a=0;a<iS2;++a){if(!iT2
xI
b.iT2)cQ}
lD2}
lB
Become
iS1
eX1){if(&b!=this&&&*data!=&*b.data){DataP
tmp=b.data;lD1
data.swap(tmp);}
}
lB
CopyOnWrite(){if(GetRefCount()>1)data=new
xK2
xF(*data);}
tK1
yZ2
yZ2::GetUniqueRef(){if(GetRefCount()>1)return
yZ2(*this,CloneTag())i1*this;}
i4):yS
cNop
iR3),n8
i4
const
xK2&b):yS
b.Opcode
iR3
b.Value),lL2(b.cO1,t53
b.Params),Hash(b.Hash),Depth(b.Depth
eB1
b.lM2){}
i4
tX1:yS
cImmed
iR3
i),n8
#ifdef __GXX_EXPERIMENTAL_CXX0X__
i4
xK2
xF&&b):yS
b.Opcode
iR3
tF2
b.Value)),lL2(b.cO1,t53
tF2
b.Params)),Hash(b.Hash),Depth(b.Depth
eB1
b.lM2){}
i4
Value_t&&i):yS
cImmed
iR3
tF2
i)),n8
#endif
i4
iU1
o):yS
o
iR3),n8
i4
iU1
o,x83
f):yS
o
iR3),lL2(f),t53),Hash(),Depth(1
eB1
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
l33
FUNCTIONPARSERTYPES;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
l33{xL1
tZ1
nO,std
c9&done,std::ostream&o){x51
tZ1
xW
a),done,o);std::ostringstream
buf
iL2
tree,buf);done[nJ2].insert(buf.str());}
}
#endif
t5{
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
xL1
DumpHashes(cE){std
c9
done;tZ1
tree,done,o);for(std
c9::const_iterator
i=done.y73
i!=done.end();++i){const
std::set<eM3>&flist=i
eC2;if(flist
eJ3!=1)o<<"ERROR - HASH COLLISION?\n"
;for(std::set<eM3>::const_iterator
j=flist.y73
j!=flist.end();++j){o<<'['<<std::hex<<i->first.hash1<<','<<i->first.hash2<<']'<<std::dec;o<<": "
<<*j<<"\n"
;}
}
}
xL1
DumpTree(cE){e32
iF3;e23
t72){case
cImmed:o<<tree
x41
t73
iU2:o<<"Var"
<<(tree.GetVar()-iU2)t73
cAdd:iF3"+"
;lC
cMul:iF3"*"
;lC
cAnd:iF3"&"
;lC
cOr:iF3"|"
;lC
cPow:iF3"^"
;break;yF3
iF3;o<<FP_GetOpcodeName(t72);cZ2
cFCall||t72==cPCall)o<<':'<<tree.GetFuncNo();}
o<<'(';if
l51<=1&&sep2[1])o<<(sep2+1)<<' ';x51{if(a>0)o<<' '
iL2
xW
a),o);if(a+1<tree.GetParamCount())o<<sep2;}
o<<')';}
xL1
DumpTreeWithIndent(cE,const
eM3&indent){o<<'['<<std::hex<<(void*)(&tree.lJ2))<<std::dec<<','<<tree.GetRefCount()<<']';o<<indent<<'_';e23
t72){case
cImmed:o<<"cImmed "
<<tree
x41;o<<'\n'
t73
iU2:o<<"VarBegin "
<<(tree.GetVar()-iU2);o<<'\n'
i1;yF3
o<<FP_GetOpcodeName(t72);cZ2
cFCall||t72==cPCall)o<<':'<<tree.GetFuncNo();o<<'\n';}
x51{eM3
ind=indent
eK3
p=0;p<ind
eJ3;p+=2)if(ind[p]=='\\')ind[p]=' ';ind+=(a+1<tree.GetParamCount())?" |"
:" \\"
;DumpTreeWithIndent(xW
a),o,ind);}
o<<std::flush;}
#endif
}
#endif
using
l33
lA1;using
l33
FUNCTIONPARSERTYPES;
#include <cctype>
l33
lA1{x83
ParamSpec_GetDepCode
iS1
eA2&b
yY3
b.first){case
ParamHolder:{cM*s=(cM*)b.second
i1
s->depcode;}
case
SubFunction:{cN*s=(cN*)b.second
i1
s->depcode;}
yF3
xY3
return
0;}
xL1
DumpParam
iS1
eA2&yA2
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
;x83
y32
0;e23
lX3.first){case
NumConstant:{const
ParamSpec_NumConstant
xF&param=*iS1
ParamSpec_NumConstant
xF*tF1;using
l33
FUNCTIONPARSERTYPES;o.precision(12);o<<param.constvalue;xY3
case
ParamHolder:{cM&yU3
cM*tF1;o<<ParamHolderNames[param.index];y32
param.constraints;xY3
case
SubFunction:{cN&yU3
cN*tF1;y32
param.constraints;yF
GroupFunction){if
yZ3
lG1==cNeg){o<<"-"
;n2}
iP1
param.lG1==cInv){o<<"/"
;n2}
else{eM3
opcode=FP_GetOpcodeName((iU1)param.lG1).substr(1
iY1
0;a<opcode
eJ3;++a)opcode[a]=(char)std::toupper(opcode[a]);o<<opcode<<"( "
;n2
o<<" )"
;}
}
else{o<<'('<<FP_GetOpcodeName((iU1)param.lG1)<<' ';yF
PositionalParams)o<<'[';yF
SelectedParams)o<<'{';n2
if
yZ3
data.n1!=0)o<<" <"
<<iH3.n1<<'>';yF
PositionalParams)o<<"]"
;yF
SelectedParams)o<<"}"
;o<<')';}
xY3
c73(ImmedConstraint_Value(constraints&ValueMask)){case
ValueMask:lC
Value_AnyNum:lC
x62:o<<"@E"
;lC
Value_OddInt:o<<"@O"
;lC
i11:o<<"@I"
;lC
Value_NonInteger:o<<"@F"
;lC
eY1:o<<"@L"
;yV3
ImmedConstraint_Sign(constraints&SignMask)){case
SignMask:lC
Sign_AnySign:lC
nD1:o<<"@P"
;lC
eZ1:o<<"@N"
;yV3
ImmedConstraint_Oneness(constraints&OnenessMask)){case
OnenessMask:lC
Oneness_Any:lC
Oneness_One:o<<"@1"
;lC
Oneness_NotOne:o<<"@M"
;yV3
ImmedConstraint_Constness(constraints&ConstnessMask)){case
ConstnessMask:lC
i01:if(lX3.first==ParamHolder){cM&yU3
cM*tF1;if
yZ3
index<2)xY3
o<<"@C"
;lC
Constness_NotConst:o<<"@V"
;lC
Oneness_Any:xY3}
xL1
DumpParams
lT1
paramlist,x83
count,std::ostream&o){for
lT1
a=0;a<count;++a){if(a>0)o<<' ';const
eA2&param=ParamSpec_Extract
xF(paramlist,a);DumpParam
xF(param,o);x83
depcode=ParamSpec_GetDepCode(param);if(depcode!=0)o<<"@D"
<<depcode;}
}
}
#include <algorithm>
using
l33
lA1;using
l33
FUNCTIONPARSERTYPES;l33{cM
plist_p[37]={{2,0,0x0}
nT
0,0x4}
nT
nD1,0x0}
nT
eZ1|Constness_NotConst,0x0}
nT
Sign_NoIdea,0x0}
nT
eY1,0x0}
,{3,Sign_NoIdea,0x0}
,{3,0,0x0}
,{3,eY1,0x0}
,{3,0,0x8}
,{3,Value_OddInt,0x0}
,{3,Value_NonInteger,0x0}
,{3,x62,0x0}
,{3,nD1,0x0}
,{0,eZ1|lV{0,lV{0,nD1|lV{0,x62|lV{0,i01,0x1}
,{0,i11|nD1|lV{0,i21
i01,0x1}
,{0,i21
lV{0,Oneness_One|lV{0,eY1|lV{1,lV{1,x62|lV{1,i21
lV{1,i11|lV{1,nD1|lV{1,eZ1|lV{6,0,0x0}
,{4,0,0x0}
,{4,i11,0x0}
,{4,lV{4,0,0x16}
,{5,0,0x0}
,{5,lV}
yO1
e72
plist_n_container{static
const
ParamSpec_NumConstant
xF
plist_n[20];}
yO1
const
ParamSpec_NumConstant
xF
plist_n_container
xF::plist_n[20]={{e62-2
i5-1
i5-0.5
i5-0.25
i5
0
tG2
fp_const_deg_to_rad
xF(tG2
fp_const_einv
xF(tG2
fp_const_log10inv
xF(i5
0.5
tG2
fp_const_log2
xF(i5
1
tG2
fp_const_log2inv
xF(i5
2
tG2
fp_const_log10
xF(tG2
fp_const_e
xF(tG2
fp_const_rad_to_deg
xF(tG2-fp_const_pihalf
xF(),xU1{xH1,xU1{fp_const_pihalf
xF(),xU1{fp_const_pi
xF(),xU1}
;cN
plist_s[517]={{{1,15,tH2
398,tH2
477,tH2
15,cNeg,GroupFunction,0}
,i01,0x1
lR3
15,y42
24,y42
465,y42
466,y42
498,cInv,lT
2,327995
x9
l0
2,48276
x9
l6
260151
x9
l6
470171
x9
l6
169126
x9
l6
48418
x9
l6
1328
x9
l6
283962
x9
l6
169275
x9
l6
39202
x9
l6
283964
x9
l6
283973
x9
l6
476619
x9
l6
296998
x9
l6
47
x9
SelectedParams,0}
,0,0x4
nH
161839
x9
l6
25036
x9
l6
35847
x9
l6
60440
x9
l6
30751
x9
l6
270599
x9
l6
60431
x9
l6
259119
x9
l6
183474
x9
l6
332066
x9
l6
7168
x9
l6
197632
x9
l6
291840
x9
l6
283648
x9
l6
238866
x9
l6
239902
x9
l6
31751
x9
l6
244743
x9
l6
384022
x9
SelectedParams,0}
,0,0x4
nH
385262
x9
l6
386086
x9
l6
393254
x9
SelectedParams,0}
,0,0x5
nH
393254
x9
l6
386095
x9
l6
387312
x9
l6
18662
x9
l6
61670
x9
l6
387397
x9
l6
247855
x9
SelectedParams,0}
,0,0x1
nH
342063
x9
l6
297007
x9
l6
15820
x9
l6
393263
x9
l6
393263
x9
SelectedParams,0}
,0,0x5
nH
161847
x9
l6
258103
x9
l6
249073
x9
l6
249076
x9
iO
0,0
x9
nF
0,0
i31
1,45
x9
nF
1,53
x9
nF
1,54
x9
nF
1,55
x9
nF
1,56
x9
nF
1,26
x9
nF
1,259
eK2
0x16
lR3
253
x9
nF
1,272
i31
1,323
eK2
0x16
lR3
0
x9
nF
1,21
x9
nF
1,447
eK2
0x4
lR3
449
eK2
0x4
lR3
0
eK2
0x4
lR3
0
tJ
2}
,0,0x4
lR3
15
x9
nF
1,24
tJ
2}
,0,0x0
nH
58392
i31
0,0
tJ
1}
,nD1,0x0
nH
24591
lS3
33807
lS3
48143
lS3
285720
lS3
290840
lS3
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
321551
iV2
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
296976,tQ1
324623,l1
0x14
nH
332815,l1
0x10}
,{{3,7340056,tQ1
289092,l9
92176
iV2
337935
e31
7340060
lO3
7340176,l9
338959
e31
7340061
iV2
7206,l9
x63
l9
357414,l9
368678,l9
370745,l1
0x7}
,{{3,7340177,l9
39277,tQ1
426398
lO3
40272286
iV2
490910
lO3
40336798
iV2
50600,l9
426462
iV2
490974
iV2
370726,l1
0x6
nH
371750,l1
0x6
nH
428070
e31
40336862
iV2
38378,l9
50671
e31
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
39338,tQ1
436262,l9
508966,l9
39409,tQ1
296998,tQ1
35847,l9
15,tQ1
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
e31
67579935,l9
39237,l9
58768,l9
62924,l9
121856,l9
15760
e31
64009216,l1
0x0}
,{{0,0,xL
0,0,iW
2,eC1
2,eD1
3,eC1
3,eD1
38,xL
1,38,iW
14,xL
1,57,xL
1,16,eL2
0x0
nH
471103,eL2
0x1
lR3
303,xL
1,323,yH3
0x0
nH
471363,eL2
0x16
lR3
293,eC1
294,eD1
295,xL
1,296,iW
400,xL
1,0,xL
1,460,xL
1,465,xL
1,16,eL2
0x1
lR3
57,yH3
0x1
lR3
0,iW
21,xL
1,15,eL2
0x0
nH
24591,xL
1,24,iW
517,yH3
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
50176,lF,178176,t21
tF3
283648,lF,19456,lF,27648,lF,89088,lF,86016,lF,488448,lF,14342,lF,58375,lF,46147
xX
46151,lF,284679,lF,7183,lF,46159
xX
38993
xX
50265,lF,50249,lF,283808,lF,284835,lF,24822,lF,10240,lF,11264,lF,7170,lF,x63
lF,17408,lF,164864,lF,237568,lF,242688,t21
0x14
nH
476160,lF,25607,lF,121871,lF,50252,lF,39374,lF,50183,lF,7192,lF,121887,lF,252979,lF,46155,lF,38919,lF,50267,lF,50268,lF,50253,lF,46190,lF,50295,lF,7563,t21
0x10
nH
416811,lF,416819,lF,40046,lF,46191
xX
415795,lF,40047
xX
415787,lF,39015,t21
0x5
nH
39326
xX
39326,lF,39332,t21
0x5
nH
39333
cC2
50590
xX
50590,lF,39338
xX
39338,lF,39335,t21
0x5
nH
15786
xX
146858,lF,39372,lF,39379,lF,39380,lF,39390
xX
50654
xX
50654,lF,24,t21
0x6
nH
62,lF,24,lF,62,t21
0x6
nH
43,lF,43
xX
51,lF,51
xX
50269,lF,50176
xX
50270,lF,39159,lF,39183
xX
7168
xX
31744,lF,99328,lF,31746,lF,100376,lF,39409
xX
39411
xX
39411,lF,39420,lF,39420
xX
15,lF,39025,t21
0x5
nH
39422,lF,16384,lF,62853,lF,15360,lF,15
cC2
16,lF,7183
cC2
7172
tW2
y71,nD1,0x0
nH
24591
tW2
lT
2,50200
tW2
lT
2,63521
tW2
lT
2,62500
tW2
lT
2,50453
tW2
lT
2,62488
tW2
lT
1,0,t93
7,t93
194,t93
0,cAcos
tA3
cAcosh
tA3
cAsin
tA3
cAsinh
nS
119,cAsinh
tA3
cAtan
t11
306176,cAtan2
t11
x63
cAtan2
tA3
cAtanh
nS
246,cCeil
tA3
cCeil
tB3
0,cD2
0,cCos
tB3
7,cD2
91,cD2
92,cD2
119,cD2
236,cD2
255,cD2
214,iW2
236,iW2
464,iW2
0,cCosh
tB3
0,iW2
0,cExp
nS
7,cExp
nS
91,cExp
tA3
yI3
7,yI3
91,yI3
246,cFloor
tA3
cFloor,lA
0x4
nH
309540,cHypot
t11
316708,cHypot
t11
316724,cHypot,l0
3,32513024,eM2
34627584
iL
31493120,eM2
89213952
iL
149042176
iL
246647808
iL
301234176
iL
494360576
iL
498558976
iL
62933520
iL
62933520,eM2
62933526
iL
62933526,eM2
24670208
iL
579378176
iL
573578240
iL
32513024
iL
566254592
iL
7900160
iL
588822528,cIf
nS
119,cInt
nS
246,tI2
0,tI2
7,tI2
31,tI2
194,tI2
363,tI2
15,cLog,lT
1,24,cLog,lT
1,0,cLog10
tA3
cLog2
t11
x63
cMax
t11
35847,cMax
t11
30751,cMax
tA3
cMax,AnyParams,1}
,0,0x4
nH
x63
cMin
t11
35847,cMin
t11
30751,cMin
tA3
cMin,AnyParams,1}
,0,0x4
nH
24591,cMin,lT
1,0,x72
7,x72
91,x72
92,x72
119,x72
149,x72
231,cSin,lA
0x5
lR3
246,x72
255,x72
254,x72
0,cSin
tB3
273,cSin,lA
0x1
lR3
214,y52
231,tD3
lA
0x5
lR3
246,y52
254,y52
255,y52
464,y52
0,cSinh
tB3
0,y52
15,cSqrt,lT
1,0,cTan
tA3
cTan
tB3
115,cTan
tB3
116
tC3
231
tC3
246
tC3
273
tC3
254
tC3
255,cTan
tA3
y62
0,cTanh
tB3
213,y62
231,y62
246,y62
254,y62
255,y62
0,cTrunc
t11
15384,cSub,lT
2,15384,cDiv,lT
2,476626,cDiv,lT
2,121913,tJ2
x63
lP3
tF3
x63
tJ2
31744,lP3
0x20
nH
31751,lP3
0x24
nH
31751,tJ2
121913,i51
t11
x63
cLess,lA
tF3
41984,cLess,lA
0x4
nH
41984,cLess
t11
7,cLess
t11
x63
cLessOrEq
t11
296182,cLessOrEq
t11
7168
tV2,lA
tF3
41984
tV2,lA
0x4
nH
41984
tV2
t11
7
tV2
t11
x63
y9
l0
2,296182,cGreaterOrEq
tA3
n32
245,n32
7,n32
550,n32
553,n32
554,n32
556,n32
31,n32
559,n32
15,n32
560,cNot
t11
7706,lT3
x63
lT3
35847,lT3
30751,lT3
463903,lT3
466975,cAnd,iO
0,0,cAnd,nF
2,x63
e03
7706,e03
35847,e03
463903,e03
466975,e03
30751,cOr,iO
1,0,n42
91,n42
131,n42
245,n42
215,n42
246,cDeg
nS
246,cRad
t11
x63
cAbsAnd,l6
x63
cAbsOr,iO
1,0,cC3
tA3
cAbsNotNot,l0
3,32513024,yT3
lA
0x0}
,}
;}
l33
lA1{const
Rule
grammar_rules[262]={{ProduceNewTree,17,1,0,{1,0,cAbs,eN2
409,{1,146,cAtan,eN2
403
nT
1324,cAtan2,eN2
405
nT
307201,cAtan2
eL
253174
nT
255224,cAtan2
eL
259324
nT
257274,cAtan2,eN2
152,{1,252,cCeil
i7
486,{1,68
yH2
482,{1,122
yH2
483,{1,124
yH2
151,{1,125
yH2
419,{1,123
yH2
0,{1,403,cCos,l2
2,1,246,{1,252,cCos,l2
18,1,0,{1,400
yH2
301,{1,406,cCosh,l2
2,1,246,{1,252,cCosh,l2
18,1,0,{1,400,cCosh
i7
458,{1,121,cFloor,eN2
150,{1,252,cFloor,tG3
156,{3,7382016,eF
549,{3,8430592,eF
556,{3,8436736,eF
157,{3,42998784,eF
550,{3,42999808,eF
562,{3,43039744,eF
557,{3,49291264,eF
538,{3,49325056,eF
469,{3,1058318,eF
473,{3,1058324,eF
473,{3,9438734,eF
469,{3,9438740,cIf,l2
0,3,32542225,{3,36732434,cIf,l2
0,3,32542231,{3,36732440,cIf,yR3
573,{3,32513026,cIf,yR3
515,{3,455505423,cIf,yR3
515,{3,433506837,cIf
i7
78,{1,256,cLog
i7
69,{1,258,cLog
i7
404,{1,72,cLog
i7
159,{1,147,cLog,l2
0,1,0
nT
487425,cMax
l3
16,1,445
nT
yJ3
cMax
l3
0,1,0
nT
483329,cMin
l3
16,1,446
nT
yJ3
cMin,c3
0,1,153
nT
24832
tW2
tG3
153
nT
25854
tW2
tG3
154
nT
129039
tW2
xV1
32055
tW2
xV1
32056
tW2
xV1
32057
tW2
l2
0,2,166288
nT
32137
tW2
xV1
33082
tW2
l2
0,2,7168
nT
12688
tW2
l2
0,2,7434
nT
12553,tX2
435
nT
46146,tX2
436
nT
46154,tX2
437
nT
46150,tX2
169
nT
83983,tX2
168
nT
130082,tX2
175
nT
133154
tW2
tI3
476160
nT
471055
tW2
tI3
274432
nT
273423
tW2
tI3
251904
nT
266274
tW2
tI3
251904
nT
263186,tX2
171,{1,252,cE2
421,{1,68,cE2
151,{1,122,cE2
419,{1,124,cE2
170,{1,125,cE2
482,{1,123,cE2
0,{1,405,cE2
172,{1,252,cSinh
i7
328,{1,404,cSinh
i7
173,{1,252,tJ3
0,{1,408,tJ3
176,{1,410,tJ3
177,{1,252,cTanh,l2
0,1,442
nT
449551,i6
1,441
nT
yJ3
i6
1,167
nT
268549,i6
1,181
nT
276749,i6
1,180
nT
276500,i6
2,190770
nT
189622,i6
2,194748
nT
193723,i6
2,202943
nT
196795,i6
2,59699
nT
298148,i6
2,59714
nT
325815,i6
2,59724
nT
343224
x9
c3
2,1,337,{1,333
tJ
1
tF
336,{1,338
tJ
1}
}
,{ReplaceParams,2,1,340
nT
1363
nR
342
nT
1365
nR
463
nT
472524
nR
47
nT
356711
nR
349
nT
200751
nR
360
nT
199727
nR
480
nT
207053
nR
481
nT
208077
nR
417
nT
211144
nR
209
nT
211145
nR
418
nT
215240
nR
212
nT
212329
nR
204
nT
373097
nR
211
nT
372944
nR
217
nT
201944
nR
221
nT
223448
nR
367
nT
508329
nR
219
nT
508126
nR
224
nT
225705
nR
223
nT
225776
nR
365
nT
230825
nR
426
nT
377057
nR
497
nT
377054
nR
497
nT
204201
nR
426
nT
375280
nR
224
nT
375006,cAdd
l3
2,2,407781
nT
233698,cAdd
l3
2,2,59763
nT
233842,i6
1,372
nT
1397,lW1
95
nT
24705,lW1
96
nT
24708,lW1
444
nT
449551,lW1
443
nT
yJ3
lW1
100
nT
101750,lW1
108
nT
106821,lW1
105
nT
103749
lU3
0,2,110607
nT
108869
lU3
0,2,107535
nT
109893,lJ
0
tF
112
nT
111634,cMul,SelectedParams,0
tF
567,{1,52,lJ
1
tF
568,{1,42,lJ
1}
}
,{ReplaceParams,2,1,467
nT
45516
xA
356
nT
51555
xA
468
nT
49612
xA
357
nT
47459
xA
429
nT
438699
xA
432
nT
441774
xA
486
nT
498726
xA
494
nT
504870
xA
382
nT
435579
xA
497
nT
435709
xA
426
nT
508287
xA
414
nT
500092
xA
499
nT
352744
xA
345
nT
367092
xA
381
nT
425318
xA
478
nT
425460
xA
47
nT
512501
xA
505
nT
355817
xA
47
nT
516598
xA
507
nT
518182
xA
508
nT
358896
xA
351
nT
388605
xA
511
nT
360939
xA
503
nT
354788
xA
514
nT
525350
xA
510
nT
394342
xA
386
nT
351346
lU3
2,2,363004
nT
361968
lU3
16,1,117
nT
1157
lU3
16,1,118
nT
1158
lU3
16,1,402
nT
411024
lU3
16,2,58768
nT
1472
lU3
16,2,15760
nT
1474
lU3
17,1,0,{1,400
lU3
17,1,57,{1,14,lJ
0}
}
,{ProduceNewTree,4,1,538
nT
41,lQ3
yK3
0
nT
5167,lQ3
cH
41984
e22
lQ3
cH
tU
lQ3
cH
eV
lQ3
cH
eW
cEqual
cN1
24849,cEqual
eL
tV
cEqual
eL
n52
281873,cEqual
eL
iP
cEqual
eL
lH1
lQ3
yK3
562
nT
41,i51,yK3
538
nT
5167,i51,cH
41984
e22
i51,cH
tU
i51,cH
eV
i51,cH
eW
i51
cN1
24849,i51
eL
tV
i51
eL
n52
281873,i51
eL
iP
i51
eL
lH1
i51,cH
tU
yN3
eV
yN3
eW
cLess,eN2
571
nT
46080,cLess
cN1
24832,cLess
eL
xW1
cLess
eL
tV
cLess
eL
n52
tK3
cLess
eL
nY1
cLess
eL
iP
cLess
eL
lH1
cLess,yL3
562
e22
yN3
tU
n62
cH
eV
n62
cH
eW
n62
eN2
565
nT
409615,cLessOrEq
cN1
24832,cLessOrEq
eL
xW1
cLessOrEq
eL
tV
cLessOrEq
eL
n52
tK3
cLessOrEq
eL
nY1
cLessOrEq
eL
iP
cLessOrEq
eL
lH1
n62
yL3
562
nT
409647,n62
cH
tU
cF2
eV
cF2
eW
cGreater,eN2
539
nT
409615
tV2
cN1
24832
tV2
eL
xW1
cGreater
eL
tV
cGreater
eL
n52
281856
tV2
eL
nY1
cGreater
eL
iP
cGreater
eL
lH1
cGreater,yL3
538
nT
409647,cF2
tU
y9
cH
eV
y9
cH
eW
y9
eN2
572
nT
46080,y9
nQ2
529654
nT
24832,y9
nQ2
xW1
y9
nQ2
tV
y9
nQ2
n52
tK3
y9
nQ2
nY1
y9
nQ2
iP
y9
nQ2
lH1
y9
yL3
538
e22
y9
yK3
519,{1,137,cNot,yR3
571,{1,2,cNot,l2
0,1,452
nT
yJ3
lV3
0,2,537097,{3,547892744,cAnd,c3
16,1,566,{1,5,cAnd,AnyParams,1}
}
,{ReplaceParams,16,1,569
nT
13314,lV3
16,1,544
nT
553498,lV3
16,1,546
nT
462369,lV3
16,1,548
nT
466465,lV3
0,1,457
nT
yJ3
xX1
570
nT
13314,xX1
563
nT
8197,xX1
541
nT
553498,xX1
542
nT
462369,xX1
543
nT
466465,xX1
564
nT
143365,cOr,c3
4,1,525,{1,137,yS3
yR3
572,{1,2,yS3
l4
17,1,0,{1,0,yS3
eN2
537,{1,256,cAbsNotNot,c3
18,1,531,{1,254,cAbsNotNot,c3
0,1,572,{3,43039744,yT3
tG3
571,{3,49325056,yT3
yR3
454,{3,32513586,yT3
l2
16,3,32542225,{3,36732434,yT3
y71}
,}
;e72
grammar_optimize_abslogical_type{y3
9
cT
grammar_optimize_abslogical_type
grammar_optimize_abslogical={9,{34,192,228,238,242,247,254,260,261}
}
;}
e72
grammar_optimize_ignore_if_sideeffects_type{y3
59
cT
grammar_optimize_ignore_if_sideeffects_type
grammar_optimize_ignore_if_sideeffects={59,{0,20,21,22,23,24,25,26,cO
l42
78,cP
cU
e72
grammar_optimize_nonshortcut_logical_evaluation_type{y3
56
cT
grammar_optimize_nonshortcut_logical_evaluation_type
grammar_optimize_nonshortcut_logical_evaluation={56,{0,25,cO
l42
78,cP
241,243,244,245,246,248,249,250,251,252,253,255,256,257,258,259}
}
;}
e72
grammar_optimize_recreate_type{y3
22
cT
grammar_optimize_recreate_type
grammar_optimize_recreate={22,{18,55,56,57,80,81,82,83,84,85,117,118,120,121,130,131,132,133,134,135,136,137}
}
;}
e72
grammar_optimize_round1_type{y3
125
cT
grammar_optimize_round1_type
grammar_optimize_round1={125,{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,19,25,cO
37,38,l42
45,46,47,48,49,50,51,52,53,54,58,59,60,61,62,63,64,65,66,67,68,69,70,71,78,79,80,81,82,83,84,85,86,87,88,93,94,95,96,97,98,99,100,101,117,118,119,120,121,122,123,124,125,126,127,128,129,138,160,161,162,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,cU
e72
grammar_optimize_round2_type{y3
103
cT
grammar_optimize_round2_type
grammar_optimize_round2={103,{0,15,16,17,25,cO
39,40,l42
45,46,47,48,49,50,51,52,53,54,59,60,72,73,78,79,86,87,88,89,90,91,92,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,119,122,123,124,125,126,127,128,139,159,160,161,162,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,cU
e72
grammar_optimize_round3_type{y3
79
cT
grammar_optimize_round3_type
grammar_optimize_round3={79,{74,75,76,77,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,170,171,172,173,174,175,176,177,181,182,183,184,185,186,187,188,189,190,191,193,194,195,196,197,198,199,201,202,203,205,206,207,208,209,210,211,213,214,215,217,218,219,220,221,222,223,225,226,227,229,230,231,232,233,234,235}
}
;}
e72
grammar_optimize_round4_type{y3
12
cT
grammar_optimize_round4_type
grammar_optimize_round4={12,{18,55,56,57,130,131,132,133,134,135,136,137}
}
;}
e72
grammar_optimize_shortcut_logical_evaluation_type{y3
53
cT
grammar_optimize_shortcut_logical_evaluation_type
grammar_optimize_shortcut_logical_evaluation={53,{0,25,cO
l42
78,cP
cU}
l33
lA1{tK1
eA2
ParamSpec_Extract
lT1
paramlist,lJ1){index=(paramlist>>(index*10))&1023;if(index>=57)return
eA2(SubFunction,cG2
plist_s[index-57]);if(index>=37)return
eA2(NumConstant,cG2
plist_n_container
xF::plist_n[index-37])i1
eA2(ParamHolder,cG2
plist_p[index]);}
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <stdio.h>
#include <algorithm>
#include <map>
#include <sstream>
using
l33
FUNCTIONPARSERTYPES;using
l33
lA1;using
t5;using
nG1;l33{nN1
It,x93
T,x93
Comp>t31
MyEqualRange(It
first,It
last,const
T&val,Comp
comp){size_t
len=last-first;while(len>0){size_t
xR3
len/2;It
x23(first);x23+=half;if(comp(*x23,val)){first=x23;++first;len=len-half-1;}
iP1
comp(val,*x23)){len=half;}
else{It
left(first);{It&eO2=left;It
last2(x23);size_t
len2=last2-eO2;while(len2>0){size_t
half2=len2/2;It
eN3(eO2);eN3+=half2;if(comp(*eN3,val)){eO2=eN3;++eO2;len2=len2-half2-1;}
else
len2=half2;}
}
first+=len;It
right(++x23);{It&eO2=right;It&last2=first;size_t
len2=last2-eO2;while(len2>0){size_t
half2=len2/2;It
eN3(eO2);eN3+=half2;if(comp(val,*eN3))len2=half2;else{eO2=eN3;++eO2;len2=len2-half2-1;}
}
eE3
t31(left,right);}
eE3
t31(first,first);}
tK1
e72
OpcodeRuleCompare{iD2()iS1
nP1
x83
y72)const{const
Rule&rule=grammar_rules[y72]i1
t72<rule
yR2.subfunc_opcode;}
iD2()lT1
y72,const
eR
const{const
Rule&rule=grammar_rules[y72]i1
rule
yR2.subfunc_opcode<t72;}
}
yO1
bool
TestRuleAndApplyIfMatch(t41
nP1
bool
c4{MatchInfo
xF
info;n31
found(false,cZ());if((rule.lI1
LogicalContextOnly)&&!c4{tR1
if(nB
IsIntType
xF::nS3){if(rule.lI1
NotForIntegers)tR1
cJ3
rule.lI1
OnlyForIntegers)tR1
if(nB
IsComplexType
xF::nS3){if(rule.lI1
NotForComplex)tR1
cJ3
rule.lI1
OnlyForComplex)tR1
for(;;){
#ifdef DEBUG_SUBSTITUTIONS
#endif
found=TestParams(rule
yR2
cU3,found.specs,info,true);if(found.found)break;if(!&*found.specs){fail:;
#ifdef DEBUG_SUBSTITUTIONS
DumpMatch
yG2,false);
#endif
return
t23}
#ifdef DEBUG_SUBSTITUTIONS
DumpMatch
yG2,true);
#endif
SynthesizeRule
yG2
yQ2}
nG1{yP1
ApplyGrammar
iS1
Grammar&tK2,nP1
bool
c4{if(tree.GetOptimizedUsing()==&tK2){
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"Already optimized:  "
iL2
tree
lE2"\n"
<<std::flush;
#endif
return
t23
if(true){bool
changed
eW3
e23
t72){case
cNot:case
cNotNot:case
cAnd:case
cOr:for
iZ1
a=0
tV3
true))yS1
lC
cIf:case
cAbsIf:if(ApplyGrammar(tK2,xW
0),t72==cIf))yS1
for
iZ1
a=1
tV3
c4)yS1
break;yF3
for
iZ1
a=0
tV3
false))yS1}
if(changed){tree.Mark_Incompletely_Hashed(yQ2}
typedef
const
x83
short*lW3;std::pair<lW3,lW3>range=MyEqualRange(tK2.rule_list,tK2.rule_list+tK2.rule_count
cU3,OpcodeRuleCompare
xF());eO3
x83
short>rules;rules.xH3
range.second-range.first);for
y0
if(IsLogisticallyPlausibleParamsMatch(tM1
yR2
cU3))rules
c52*r);}
range.first=!rules
eU3?&rules[0]:0;range.second=!rules
eU3?&rules[rules
eJ3-1]+1:0;if(range.first!=range.second){
#ifdef DEBUG_SUBSTITUTIONS
if(range.first!=range.second)cH2"Input ("
<<FP_GetOpcodeName(t72)<<")["
<<tree.GetParamCount()<<"]"
;if(c4
cV2<<"(Logical)"
;x83
first=l52,prev=l52;e32
sep=", rules "
;for
y0
if(first==l52)first=prev=*r;iP1*r==prev+1)prev=*r;else
cH2
sep<<first;sep=","
;if(prev!=first)cV2<<'-'<<prev;first=prev=*r;}
}
if(first!=l52)cH2
sep<<first;if(prev!=first)cV2<<'-'<<prev;}
cV2<<": "
iL2
tree
lE2"\n"
<<std::flush;}
#endif
bool
changed
eW3
for
y0
#ifndef DEBUG_SUBSTITUTIONS
if(!IsLogisticallyPlausibleParamsMatch(tM1
yR2
cU3))yA1
#endif
if(TestRuleAndApplyIfMatch(tM1
cU3,c4){yS1
xY3}
if(changed){
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"Changed."
<<std::endl;cV2<<"Output: "
iL2
tree
lE2"\n"
<<std::flush;
#endif
tree.Mark_Incompletely_Hashed(yQ2}
tree.SetOptimizedUsing(&tK2)i1
t23
yP1
ApplyGrammar
iS1
void*p,FPoptimizer_CodeTree::eR
yP
ApplyGrammar(*iS1
Grammar*)p
cU3);}
xL1
ApplyGrammars(FPoptimizer_CodeTree::eR{
#ifdef DEBUG_SUBSTITUTIONS
cV2<<iB3"grammar_optimize_round1\n"
;
#endif
n6
grammar_optimize_round1
cU3
nW
#ifdef DEBUG_SUBSTITUTIONS
cV2<<iB3"grammar_optimize_round2\n"
;
#endif
n6
grammar_optimize_round2
cU3
nW
#ifdef DEBUG_SUBSTITUTIONS
cV2<<iB3"grammar_optimize_round3\n"
;
#endif
n6
grammar_optimize_round3
cU3
nW
#ifndef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
#ifdef DEBUG_SUBSTITUTIONS
cV2<<iB3"grammar_optimize_nonshortcut_logical_evaluation\n"
;
#endif
n6
grammar_optimize_nonshortcut_logical_evaluation
cU3
nW
#endif
#ifdef DEBUG_SUBSTITUTIONS
cV2<<iB3"grammar_optimize_round4\n"
;
#endif
n6
grammar_optimize_round4
cU3
nW
#ifdef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
#ifdef DEBUG_SUBSTITUTIONS
cV2<<iB3"grammar_optimize_shortcut_logical_evaluation\n"
;
#endif
n6
grammar_optimize_shortcut_logical_evaluation
cU3
nW
#endif
#ifdef FP_ENABLE_IGNORE_IF_SIDEEFFECTS
#ifdef DEBUG_SUBSTITUTIONS
cV2<<iB3"grammar_optimize_ignore_if_sideeffects\n"
;
#endif
n6
grammar_optimize_ignore_if_sideeffects
cU3
nW
#endif
#ifdef DEBUG_SUBSTITUTIONS
cV2<<iB3"grammar_optimize_abslogical\n"
;
#endif
n6
grammar_optimize_abslogical
cU3
nW
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
l33
FUNCTIONPARSERTYPES;using
l33
lA1;using
t5;using
nG1;l33{yP1
TestImmedConstraints
lT1
bitmask,const
eR{e23
bitmask&ValueMask){case
Value_AnyNum:case
ValueMask:lC
x62:if(GetEvennessInfo
iS3
n72
Value_OddInt:if(GetEvennessInfo
iS3
x82
i11:if(GetIntegerInfo
iS3
n72
Value_NonInteger:if(GetIntegerInfo
iS3
x82
eY1:if(!IsLogicalValue(tree)cQ
nE1
SignMask){case
Sign_AnySign:lC
nD1:if(l81
n72
eZ1:if(l81
x82
Sign_NoIdea:if(l81
Unknown
cQ
nE1
OnenessMask){case
Oneness_Any:case
OnenessMask:lC
Oneness_One:if(!xY1
if(!lJ3
fp_abs(tree
x41),e62
1))cQ
lC
Oneness_NotOne:if(!xY1
lG3
fp_abs(tree
x41),e62
1))cQ
nE1
ConstnessMask){case
Constness_Any:lC
i01:if(!xY1
lC
Constness_NotConst:if(xY1
xY3
lD2}
tP1
x83
extent,x83
nbits,x93
eP2=x83
int>e72
nbitmap{private:static
const
x83
bits_in_char=8;static
const
x83
eQ2=(sizeof(eP2)*bits_in_char)/nbits;eP2
data[(extent+eQ2-1)/eQ2];e43
void
inc(lJ1,int
by=1){data[pos(index)]+=by*eP2(1<<y82);y81
void
dec(lJ1){inc(index,-1);}
int
get(lJ1
c01(data[pos(index)]>>y82)&mask()iT3
pos(lJ1)yP
index/eQ2
iT3
shift(lJ1)yP
nbits*(index%eQ2)iT3
mask()yP(1<<nbits)-1
iT3
mask(lJ1)yP
mask()<<y82;}
}
;e72
Needs{int
SubTrees:8;int
Others:8;int
l62:8;int
Immeds:8;nbitmap<iU2,2>SubTreesDetail;Needs(){std::memset(this,0,sizeof(*this));}
Needs
iS1
Needs&b){std::memcpy(this,&b,sizeof(b));}
Needs&eM1=iS1
Needs&b){std::memcpy(this,&b,sizeof(b))i1*this;}
}
yO1
Needs
CreateNeedList_uncached(eX&tL2{Needs
yQ1;for
lT1
a=0;a<params
y92;++a){const
eA2&lX3=ParamSpec_Extract
xF(params.param_list,a);e23
lX3.first){case
SubFunction:{cN&yU3
cN*tF1;yF
GroupFunction)++yQ1.Immeds;else{++iO2;assert(iH3.subfunc_opcode<VarBegin);yQ1.SubTreesDetail.inc
yZ3
lG1);}
++yQ1.l62;xY3
case
NumConstant:case
ParamHolder:++iN2;++yQ1.l62;xY3
eE3
yQ1;}
tK1
Needs&CreateNeedList(eX&tL2{typedef
std::map<eX*,Needs>eE1;static
eE1
yV1;eE1::y83
i=yV1.xU2&tL2;if(i!=yV1.cX1&tL2
return
i
eC2
i1
yV1.y43,std::make_pair(&params,CreateNeedList_uncached
xF(tL2))eC2;}
tK1
yZ2
CalculateGroupFunction
iS1
eA2&yA2
const
tC
info
yY3
lX3.first){case
NumConstant:{const
ParamSpec_NumConstant
xF&param=*iS1
ParamSpec_NumConstant
xF*tF1
i1
CodeTreeImmed
yZ3
constvalue);}
case
ParamHolder:{cM&yU3
cM*tF1
i1
info.GetParamHolderValueIfFound
yZ3
index);}
case
SubFunction:{cN&yU3
cN*tF1
tX
nS3;nS3
iH
param.lG1);nS3.lJ2).xH3
iH3
y92);for
lT1
a=0;a<iH3
y92;++a
nQ
tmp(CalculateGroupFunction
nF1
iH3.param_list,a),info));nS3
c91
tmp);}
nS3
iW1
i1
nS3;}
eE3
yZ2();}
}
nG1{yP1
IsLogisticallyPlausibleParamsMatch(eX&params,const
eR{Needs
yQ1(CreateNeedList
xF(tL2);size_t
tL3=yC
if(tL3<size_t(yQ1.l62))tM2
for
iZ1
a=0;a<tL3;++a){x83
opcode=xW
a)nC;e23
opcode){case
cImmed:if(yQ1.Immeds>0)--yQ1.Immeds;else--iN2;lC
iU2:case
cFCall:case
cPCall:--iN2;break;yF3
assert(opcode<VarBegin);if(iO2>0&&yQ1.SubTreesDetail.get(opcode)>0){--iO2;yQ1.SubTreesDetail.dec(opcode);}
else--iN2;}
}
if(yQ1.Immeds>0||iO2>0||iN2>0)tM2
if(params.match_type!=AnyParams){if(0||iO2<0||iN2<0)tM2}
lD2}
tK1
n31
TestParam
iS1
eA2&yA2
const
nP1
const
cZ&start_at,tC
info
yY3
lX3.first){case
NumConstant:{const
ParamSpec_NumConstant
xF&param=*iS1
ParamSpec_NumConstant
xF*tF1;if(!xY1
Value_t
imm=tree
x41;switch
yZ3
modulo){case
Modulo_None:lC
Modulo_Radians:imm=cK3
imm,yK
imm<xG1
imm
yN
if(imm>fp_const_pi
xF())imm-=fp_const_twopi
xF(e02
return
lJ3
imm,param.constvalue);}
case
ParamHolder:{cM&yU3
cM*tF1;if(!x3
return
info.SaveOrTestParamHolder
yZ3
index
cU3);}
case
SubFunction:{cN&yU3
cN*tF1;yF
GroupFunction){if(!x3
yZ2
xZ1=CalculateGroupFunction(yA2
info);
#ifdef DEBUG_SUBSTITUTIONS
DumpHashes(xZ1
lE2*iS1
void**)&xZ1
x41;cV2<<"\n"
;cV2<<*iS1
void**)&tree
x41;cV2<<"\n"
;DumpHashes(tree
lE2"Comparing "
iL2
xZ1
lE2" and "
iL2
tree
lE2": "
;cV2<<(xZ1
xI
tree)?"true"
:"false"
lE2"\n"
;
#endif
return
xZ1
xI
tree);}
cJ3!eK1){if(!x3
if(t72!=param.lG1
cQ
eE3
TestParams
yZ3
data
cU3,start_at,info,false);}
}
eE3
t23
tK1
e72
l31
x92
MatchInfo
xF
info;l31():start_at(),info(){}
}
yO1
class
MatchPositionSpec_PositionalParams:nQ1
l31
xF>{e43
iX2
MatchPositionSpec_PositionalParams
iZ1
n):eF1
l31
xF>(n){}
}
;e72
l72
x92
l72():start_at(){}
}
;class
yG:nQ1
l72>{e43
x83
trypos;iX2
yG
iZ1
n):eF1
l72>(n),trypos(0){}
}
yO1
n31
TestParam_AnyWhere
iS1
eA2&yA2
const
nP1
const
cZ&start_at,tC
info,cY2&used,bool
iX1{xP<yG>x8;x83
xA2(yG*)eK1;a=x8->trypos;goto
retry_anywhere_2;cO2
yG
l51);a=0;}
tM3
yC++a){if(used[a])yA1
retry_anywhere
yW3
TestParam(yA2
xW
a)xS3);yP3
used[a]=true;if(iX1
info.SaveMatchedParamIndex(a);x8->trypos=a
i1
n31(true,&*x8);}
}
retry_anywhere_2
yX3
goto
retry_anywhere;}
eE3
t23
tK1
e72
yD1
x92
MatchInfo
xF
info;cY2
used;iX2
yD1
iZ1
tL3):start_at(),info(),used(tL3){}
}
yO1
class
MatchPositionSpec_AnyParams:nQ1
yD1
xF>{e43
iX2
MatchPositionSpec_AnyParams
iZ1
n,size_t
m):eF1
yD1
xF>(n,yD1
xF(m)){}
}
yO1
n31
TestParams(eX&nN,const
nP1
const
cZ&start_at,tC
info,bool
iX1{if(nN.match_type!=AnyParams){if(xT!=tree.GetParamCount()cQ}
if(!IsLogisticallyPlausibleParamsMatch(nN
cU3))tM2
e23
nN.match_type){case
PositionalParams:{xP<cF>x8;x83
xA2(cF*)eK1;a=xT-1;goto
lK1;cO2
cF(xT);a=0;}
tM3
xT;++a){(yQ3=info;retry_positionalparams
yW3
TestParam(cW
a),xW
a)xS3);yP3
yA1}
}
lK1
yX3
info=(yQ3;goto
retry_positionalparams;}
if(a>0){--a;goto
lK1;}
info=c12
i1
t23
if(iX1
tN3
info.SaveMatchedParamIndex(a)i1
n31(true,&*x8);}
case
SelectedParams:case
AnyParams:{xP<t6>x8;cY2
used
l51);eO3
x83>iY2(xT);eO3
x83>yB2(xT);tN3{const
eA2
lX3=cW
a);iY2[a]=ParamSpec_GetDepCode(lX3);}
{x83
b=0;tN3
if(iY2[a]!=0)yB2[b++]=a;tN3
if(iY2[a]==0)yB2[b++]=a;}
x83
xA2(t6*)eK1;if(xT==0){a=0;goto
retry_anyparams_4;}
a=xT-1;goto
eG1;cO2
t6(xT
cU3.GetParamCount());a=0;if(xT!=0){c12=info;(*x8)[0].used=used;}
}
tM3
xT;++a){if(a>0){(yQ3=info;(yO3
used=used;}
retry_anyparams
yW3
TestParam_AnyWhere
xF(cW
yB2[a])cU3
xS3,used,iX1;yP3
yA1}
}
eG1
yX3
info=(yQ3;used=(yO3
used;goto
retry_anyparams;}
eH1:if(a>0){--a;goto
eG1;}
info=c12
i1
t23
retry_anyparams_4:if(nN.n1!=0){if(!TopLevel||!info.HasRestHolder(nN.n1)){eK
cI2;cI2.reserve
l51);for
lT1
b=0;b<yC++b){if(used[b])yA1
cI2
c52
xW
b));used[b]=true;if(iX1
info.SaveMatchedParamIndex(b);}
if(!info.SaveOrTestRestHolder(nN.n1,cI2)){goto
eH1;}
}
else{eA3&cI2=info.GetRestHolderValues(nN.n1
iY1
0;a<cI2
eJ3;++a){bool
found
eW3
for
lT1
b=0;b<yC++b){if(used[b])yA1
if(cI2[a]xI
xW
b))){used[b]=true;if(iX1
info.SaveMatchedParamIndex(b);found=true;xY3}
if(!found){goto
eH1;}
}
}
eE3
n31(true,xT?&*x8:0);}
case
GroupFunction:xY3
return
t23}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
#include <algorithm>
#include <assert.h>
using
t5;using
nG1;l33{tK1
yZ2
yW1
iS1
eA2&yA2
tC
info,bool
inner=true
yY3
lX3.first){case
NumConstant:{const
ParamSpec_NumConstant
xF&param=*iS1
ParamSpec_NumConstant
xF*tF1
i1
CodeTreeImmed
yZ3
constvalue);}
case
ParamHolder:{cM&yU3
cM*tF1
i1
info.GetParamHolderValue
yZ3
index);}
case
SubFunction:{cN&yU3
cN*tF1
tX
tree;tD2
param.lG1);for
lT1
a=0;a<iH3
y92;++a
nQ
nparam=yW1
nF1
iH3.param_list,a
l04
true
c81
nparam);}
if
yZ3
data.n1!=0){eK
trees(info.GetRestHolderValues
yZ3
data.n1)yB
AddParamsMove(trees);if
l51==1){assert(tree.GetOpcode()==cAdd||tree.GetOpcode()==cMul||tree.GetOpcode()==cMin||tree.GetOpcode()==cMax||tree.GetOpcode()==cAnd||tree.GetOpcode()==cOr||tree.GetOpcode()==cAbsAnd||tree.GetOpcode()==cAbsOr);tree.eF2
0));}
else
if
l51==0
yY3
t72){case
cAdd:case
cOr:tree=nA1
0));lC
cMul:case
cAnd:tree=nA1
1));yF3
xY3}
}
if(inner)tree
iW1
i1
tree;}
eE3
yZ2();}
}
nG1{xL1
SynthesizeRule(t41
nP1
tC
info
yY3
rule.ruletype){case
ProduceNewTree:{tree.Become(yW1
nF1
rule.t51
0
l04
false)e02
case
ReplaceParams:yF3{eO3
x83>list=info.GetMatchedParamIndexes();std::sort(list.iJ2
list.end()iY1
list
eJ3;a-->0;)tree.i91
list[a]);for
lT1
a=0;a<rule.repl_param_count;++a
nQ
nparam=yW1
nF1
rule.t51
a
l04
true
c81
nparam);}
xY3}
}
}
#endif
#ifdef DEBUG_SUBSTITUTIONS
#include <sstream>
#include <cstring>
using
l33
FUNCTIONPARSERTYPES;using
l33
lA1;using
t5;using
nG1;l33
lA1{xL1
DumpMatch(t41
const
nP1
const
tC
info,bool
DidMatch,std::ostream&o){DumpMatch
yG2,DidMatch?iX3"match"
:iX3"mismatch"
,o);}
xL1
DumpMatch(t41
const
nP1
const
tC
info,e32
tO3,std::ostream&o){static
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
;o<<tO3<<" (rule "
<<(&rule-grammar_rules)<<")"
<<":\n  Pattern    : "
;{eA2
tmp;tmp.first=SubFunction;ParamSpec_SubFunction
tmp2;tmp2.data=rule
yR2;tmp.second=cG2
tmp2;DumpParam
xF(tmp,o);}
o<<"\n  Replacement: "
;DumpParams
xF(rule.t51
rule.repl_param_count,o)iV3
o<<"  Tree       : "
iL2
tree,o)iV3
if(!std::strcmp(tO3,iX3"match"
))DumpHashes(tree,o
iY1
0;a<info.cA
eJ3;++a){if(!info.cA[a].nX2)yA1
o<<"           "
<<ParamHolderNames[a]<<" = "
iL2
info.cA[a],o)iV3}
e11
info.lQ
eJ3;++b){if(!iW3
eH3)yA1
for
iZ1
a=0;a<iW3
eI3
eJ3;++a){o<<"         <"
<<b<<"> = "
iL2
iW3
eI3[a],o);o<<std::endl;}
}
o<<std::flush;}
}
#endif
#include <list>
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
using
l33
FUNCTIONPARSERTYPES;l33{yP1
MarkIncompletes(FPoptimizer_CodeTree::eR{if(tree.Is_Incompletely_Hashed())lD2
bool
l82
eW3
x51
l82|=MarkIncompletes(xW
a));if(l82)tree.Mark_Incompletely_Hashed()i1
l82;}
xL1
FixIncompletes(FPoptimizer_CodeTree::eR{if(tree.Is_Incompletely_Hashed()){x51
FixIncompletes(xW
a)yB
e33);}
}
}
t5{lB
Sort()iZ3
Sort();}
lB
e33
bool
constantfolding){if(constantfolding)ConstantFolding(*this);else
Sort();data
xC
tK1
e72
cB
c03
Value_t
c13
yE1=0;
#if 0
long
double
value=Value;e9=crc32::calc(iS1
x83
char*)&value,sizeof(value));key^=(key<<24);
#elif 0
union{e72{x83
char
filler1[16]iK2
v;x83
char
filler2[16];}
buf2;e72{x83
char
filler3[sizeof
cX3)+16-sizeof(xE1)];e9;}
buf1;}
data;memset(&data,0,sizeof(data));data.buf2.v=Value;e9=data.buf1.key;
#else
int
exponent
iK2
nR2=std::frexp(Value,&tU2
e9=lT1(exponent+0x8000)&0xFFFF);if(nR2<0){nR2=-nR2;key=key^0xFFFF;}
else
key+=0x10000;nR2-=yF2;key<<=39;key|=n51(nR2+nR2)*e62
1u<<31))<<8;
#endif
lP
#ifdef FP_SUPPORT_COMPLEX_NUMBERS
nN1
T
cJ2
std::complex<T> >c03
std::complex<T>c13
cB<T>::lY3
cK2,Value.real());nB
fphash_t
temp;cB<T>::lY3
temp,Value.imag());yE1^=temp.hash2;cK2.hash2^=temp.hash1;}
}
;
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
tP1
cJ2
long>yD
long
Value){e9=Value;lP
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
tP1
cJ2
GmpInt>c03
GmpInt
c13
e9=Value.toInt();lP
#endif
xL1
xK2
xF::Recalculate_Hash_NoRecursion(){fphash_t
cK2(n51
Opcode)<<56,Opcode*iE3(0x1131462E270012B));Depth=1;e23
Opcode){case
cImmed:{cB
xF::lY3
cK2,Value
e02
case
iU2:{yE1|=n51
cO1<<48
lL1((n51
cO1)*11)^iE3(0x3A83A83A83A83A0);xY3
case
cFCall:case
cPCall:{yE1|=n51
cO1<<48
lL1((~n51
cO1)*7)^3456789;}
yF3{size_t
t61=0
e73
0;a<iS2;++a){if(iT2
nP2>t61)t61=iT2
nP2;yE1+=((iT2.tN2
hash1*(a+1))>>12)lL1
iT2.tN2
hash1
lL1(3)*iE3(0x9ABCD801357);cK2.hash2*=iE3(0xECADB912345)lL1(~iT2.tN2
hash2)^4567890;}
Depth+=t61;}
}
if(Hash!=cK2){Hash=cK2;lM2=0;}
}
lB
tE3{MarkIncompletes(*this);FixIncompletes(*this);}
}
#endif
#include <cmath>
#include <list>
#include <cassert>
#ifdef FP_SUPPORT_OPTIMIZER
using
l33
FUNCTIONPARSERTYPES;l33{using
t5
yO1
bool
x01
const
nP1
long
count,const
xK1::SequenceOpCode
xF&eQ,xK1::eP3&synth,size_t
max_bytecode_grow_length);static
const
e72
SinCosTanDataType{OPCODE
whichopcode;OPCODE
inverse_opcode;enum{eR2,nS2,nZ1,lM1}
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
,{cL2{tD3
cCosh,cM2,{tD3
cNop,{cL2
cNop,cCosh}
}
,{cCosh,cNop,{tD3
cL2
cNop}
}
,{cNop,cTanh,{cCosh,tD3
cM2,{cNop,tD3{cNop,cTanh,cCosh,cNop}
}
,{cNop,cCosh,{cTanh,tD3
cM2}
;}
t5{lB
SynthesizeByteCode(eO3
x83>&eN1,std
xM3
xF&Immed,size_t&stacktop_max){
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"Making bytecode for:\n"
;l21
#endif
while(RecreateInversionsAndNegations()){
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"One change issued, produced:\n"
;l21
#endif
tE3;using
nG1;using
l33
lA1;const
void*g=cG2
grammar_optimize_recreate;while(ApplyGrammar(*iS1
Grammar*)g,*this)){tE3;}
}
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"Actually synthesizing, after recreating inv/neg:\n"
;l21
#endif
xK1::eP3
synth;SynthesizeByteCode(synth,false
iM2.Pull(eN1,Immed,stacktop_max);}
lB
SynthesizeByteCode(xK1::eP3&synth,bool
MustPopTemps)const{y01*this))yP;}
for
iZ1
a=0;a<12;++a){const
SinCosTanDataType&data=SinCosTanData[a];if(data.whichopcode!=cNop)iC1!=data.whichopcode)yA1
c02
lZ3;lZ3.xA1);lZ3
iH
data.inverse_opcode);lZ3.yC2);y01
lZ3)){synth.tQ
else
iC1!=cInv)yA1
if(GetParam(0)nC!=data.inverse_opcode)yA1
y01
GetParam(0))){synth.tQ
size_t
found[4];e11
4;++b){c02
tree;if(data.tP3]==cNop){tD2
cInv);c02
n03;n03.xA1);n03
iH
data.tP3^2]);n03.yC2
c81
n03);}
else{tree.xA1
yB
eL3
data.tP3]);}
tree.yC2);found[b]=synth.xK3
tree);}
if(found[data.eR2
e3
nS2
tG
eR2]iM2.yT
nS2]);yB1
cDiv
tW
eR2
e3
lM1
tG
eR2]iM2.yT
lM1]);yB1
cMul
tW
nZ1
e3
lM1
tG
nZ1]iM2.yT
lM1]);yB1
cRDiv
tW
nZ1
e3
nS2
tG
nZ1]iM2.yT
nS2]);yB1
cMul,2,1
iM2.tQ
size_t
n_subexpressions_synthesized=SynthCommonSubExpressions(c92
e23
l92{case
iU2:synth.PushVar(GetVar());lC
cImmed:t71
GetImmed());lC
cAdd:case
cMul:case
cMin:case
cMax:case
cAnd:case
cOr:case
cAbsAnd:case
cAbsOr:iC1==cMul){bool
c23
eW3
yH
lX1
tD1&&isLongInteger(lX1
x41)){c11=makeLongInteger(lX1
x41);c02
tmp(*this,x93
c02::CloneTag());tmp.i91
a);tmp
iW1;if(x01
tmp,value,xK1::tU1
xF::AddSequence,synth,MAX_MULI_BYTECODE_LENGTH)){c23=true;xY3}
}
if(c23)xY3
int
yF1=0;cY2
done(GetParamCount(),false);c02
iQ;iQ
iH
l92;for(;;){bool
found
eW3
yH
done[a])yA1
if(synth.IsStackTop(lX1)){found=true;done[a]=true;lX1.n7
iQ
cS
lX1);if(++yF1>1){yU
2);iQ.yC2)yC1
iQ);yF1=yF1-2+1;}
}
}
if(!found)xY3
yH
done[a])yA1
lX1.n7
iQ
cS
lX1);if(++yF1>1){yU
2);iQ.yC2)yC1
iQ);yF1=yF1-2+1;}
}
if(yF1==0
yY3
l92{case
cAdd:case
cOr:case
cAbsOr:t71
0);lC
cMul:case
cAnd:case
cAbsAnd:t71
1);lC
cMin:case
cMax:t71
0);break;yF3
xY3++yF1;}
assert(n_stacked==1);xY3
case
cPow:{tW1
p0
i42
0);tW1
p1
i42
1);if(!p1
tD1||!isLongInteger
cR3||!x01
p0,makeLongInteger
cR3,xK1::tU1
xF::MulSequence,synth,MAX_POWI_BYTECODE_LENGTH)){p0.n7
p1.n7
yU
2);e21
cIf:case
cAbsIf:{x93
xK1::eP3::IfData
ifdata;GetParam(0)e93
SynthIfStep1(ifdata,l92;GetParam(1)e93
SynthIfStep2(ifdata);GetParam(2)e93
SynthIfStep3(ifdata
e02
case
cFCall:case
cPCall:{for
iZ1
a=0;a<l91++a)lX1.n7
yU
lT1)GetParamCount());yB1
0x80000000u|GetFuncNo(),0,0
e02
yF3{for
iZ1
a=0;a<l91++a)lX1.n7
yU
lT1)GetParamCount()e02}
synth.StackTopIs(*this);if(MustPopTemps&&n_subexpressions_synthesized>0){size_t
top=synth.GetStackTop(iM2.DoPopNMov(top-1-n_subexpressions_synthesized,top-1);}
}
}
l33{yP1
x01
const
nP1
long
count,const
xK1::SequenceOpCode
xF&eQ,xK1::eP3&synth,size_t
max_bytecode_grow_length){if
cW3!=0){xK1::eP3
backup=synth;tree.n7
size_t
bytecodesize_backup=synth.GetByteCodeSize();xK1::x01
count
e63
size_t
bytecode_grow_amount=synth.GetByteCodeSize()-bytecodesize_backup;if(bytecode_grow_amount>max_bytecode_grow_length){synth=backup
i1
t23
lD2}
else{xK1::x01
count,eQ,synth
yQ2}
}
#endif
#include <cmath>
#include <cassert>
#ifdef FP_SUPPORT_OPTIMIZER
using
l33
FUNCTIONPARSERTYPES;l33{using
t5;
#define FactorStack std xM3
const
e72
PowiMuliType{x83
opcode_square;x83
opcode_cumulate;x83
opcode_invert;x83
opcode_half;x83
opcode_invhalf;}
iseq_powi={cSqr,cMul,cInv,cSqrt,cRSqrt}
,iseq_muli={l52
x9
cNeg,l52,l52}
yO1
Value_t
cP1
const
PowiMuliType&c33,const
std
nA3,n82&stack)eQ3
1);while(IP<limit){if(n13
c33.opcode_square){if(!t42
eT3
2;yZ
opcode_invert){nS3=-nS3;yZ
opcode_half){if(nS3>xH1&&isEvenInteger(eT3
yF2;yZ
opcode_invhalf){if(nS3>xH1&&isEvenInteger(eT3
e62-0.5);++IP;yA1}
size_t
xB2=IP
iK2
lhs(1);if(n13
cFetch){lJ1=n43;if(index<y2||size_t(index-y2)>=lB2){IP=xB2;xY3
lhs=stack[index-y2];goto
yD2;}
if(n13
cDup){lhs=nS3;goto
yD2;yD2:t91
nS3);++IP
iK2
subexponent=cP1
c33
iK
if(IP>=limit||eN1[IP]!=c33.opcode_cumulate){IP=xB2;xY3++IP;stack.pop_back();nS3+=lhs*subexponent;yA1}
xY3
return
nS3;}
tK1
Value_t
ParsePowiSequence
iS1
std
nA3){n82
stack;t91
e62
1))i1
cP1
iseq_powi
iK}
tK1
Value_t
ParseMuliSequence
iS1
std
nA3){n82
stack;t91
e62
1))i1
cP1
iseq_muli
iK}
tK1
class
CodeTreeParserData{e43
iX2
CodeTreeParserData(bool
k_powi):stack(),clones(),keep_powi(k_powi){}
void
Eat
iZ1
tL3,OPCODE
opcode
nQ
xM;xM
iH
opcode);eK
params=Pop(tL3
yR1
tL2;if(!keep_powi)e23
opcode){case
cTanh:{yZ2
sinh,cosh;sinh
iH
cSinh);sinh
cS
xM
c43
sinh
iW1;cosh
iH
cCosh);cosh
c91
xM
c43
cosh
iW1
tX
pow;pow
iH
cPow);pow
c91
cosh);pow.yA
nU2));pow
iW1;xM
iH
cMul);xM.nB1
0,sinh);xM
c91
pow
e02
case
cTan:{yZ2
sin,cos;sin
iH
cSin);sin
cS
xM
c43
sin
iW1;cos
iH
cCos);cos
c91
xM
c43
cos
iW1
tX
pow;pow
iH
cPow);pow
c91
cos);pow.yA
nU2));pow
iW1;xM
iH
cMul);xM.nB1
0,sin);xM
c91
pow
e02
case
cPow:{iT1&p0=xM
l8
0);iT1&p1=xM
l8
1);if(p1
nC==cAdd){eK
mulgroup(p1.GetParamCount()iY1
0;a<p1.l91++a
nQ
pow;pow
iH
cPow);pow
cS
p0);pow
cS
p1
n93
pow
iW1;mulgroup[a
tQ3
pow);}
xM
iH
cMul
yR1
mulgroup);}
xY3
yF3
xY3
xM.e33!keep_powi);lA2,false);
#ifdef DEBUG_SUBSTITUTIONS
iA1<<tL3<<", "
<<FP_GetOpcodeName(opcode)<<"->"
<<FP_GetOpcodeName(xM
nC)<<": "
iK3
xM);iB1
xM);
#endif
t91
xM
cN2
EatFunc
iZ1
tL3,OPCODE
opcode,x83
funcno
nQ
xM=CodeTreeFuncOp
xF(opcode,funcno);eK
params=Pop(tL3
yR1
tL2;xM.yC2);
#ifdef DEBUG_SUBSTITUTIONS
iA1<<tL3<<", "
iK3
xM);iB1
xM);
#endif
lA2);t91
xM
cN2
AddConst(yL1
nQ
xM=CodeTreeImmed(value);lA2);Push(xM
cN2
AddVar
lT1
varno
nQ
xM=CodeTreeVar
xF(varno);lA2);Push(xM
cN2
SwapLastTwoInStack(){tA1
1
tQ3
tA1
2]cN2
Dup(){Fetch(lB2-1
cN2
Fetch
iZ1
which){Push(stack[which]);}
nN1
T>void
Push(T
tree){
#ifdef DEBUG_SUBSTITUTIONS
cV2<<iK3
tree);iB1
tree);
#endif
t91
tree
cN2
PopNMov
iZ1
target,size_t
source){stack[target]=stack[source];stack
xF3
target+1);}
yZ2
yE2{clones.clear()tX
nS3(stack.back());stack
xF3
lB2-1)i1
nS3;}
eK
Pop
iZ1
n_pop){eK
nS3(n_pop);for
lT1
n=0;n<n_pop;++n)nS3[n
tQ3
tA1
n_pop+n]);
#ifdef DEBUG_SUBSTITUTIONS
for
iZ1
n=n_pop;n-->0;){iA1
iL2
nS3[n]);iB1
nS3[n]);}
#endif
stack
xF3
lB2-n_pop)i1
nS3;}
size_t
GetStackTop(c01
lB2;}
private:void
FindClone(yZ2&,bool=true)yP;}
private:eK
stack;std::multimap<fphash_t,yZ2>clones;bool
keep_powi;private:CodeTreeParserData
iS1
CodeTreeParserData&);CodeTreeParserData&eM1=iS1
CodeTreeParserData&);}
yO1
e72
IfInfo{yZ2
eS2
tX
thenbranch;size_t
endif_location;IfInfo():eS2(),thenbranch(),endif_location(){}
}
;}
t5{lB
GenerateFrom
iS1
x93
FunctionParserBase
xF::Data&x33,bool
keep_powi){eK
xH2;xH2.xH3
x33.mVariablesAmount);for
lT1
n=0;n<x33.mVariablesAmount;++n){xH2
c52
CodeTreeVar
xF(n+iU2));}
GenerateFrom(x33,xH2,keep_powi);}
lB
GenerateFrom
iS1
x93
FunctionParserBase
xF::Data&x33,const
x5&xH2,bool
keep_powi){const
eO3
x83>&eN1=x33.mByteCode;const
std
xM3
xF&Immed=x33.mImmed;
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"ENTERS GenerateFrom()\n"
;
#endif
CodeTreeParserData
xF
sim(keep_powi);eO3
IfInfo
xF>eJ
eK3
IP=0,DP=0;;++IP){tO2:while(!eJ
eU3&&(eJ.eA==IP||(IP<n41&&n13
cJump&&eJ
e81.nX2))){c02
elsebranch=sim.yE2
nA2
eJ.back().eS2)nA2
eJ
e81)nA2
elsebranch)yR
3,cIf);eJ.pop_back();}
if(IP>=n41)break;x83
opcode=eN1[IP];if(eS3
cSqr||opcode==cDup||eS3
cInv&&!IsIntType
xF::nS3)||opcode==cNeg||opcode==cSqrt||opcode==cRSqrt||opcode==cFetch)){size_t
was_ip=IP
iK2
eD2
ParsePowiSequence
xF(eN1,IP,eJ
eU3?n41:eJ.eA,sim.xH
1);if(exponent!=e62
1.0)){xG
exponent
iU3);goto
tO2;}
if
eS3
cDup||opcode==cFetch||opcode==cNeg)l14
xT2=ParseMuliSequence
xF(eN1,IP,eJ
eU3?n41:eJ.eA,sim.xH
1);if(xT2!=e62
1.0)){xG
xT2)yR
2,cMul);goto
tO2;}
}
IP=was_ip;}
if(nB2>=iU2){lJ1=opcode-iU2
nA2
xH2[index]);}
else{e23
nB2){case
cIf:case
cAbsIf:{eJ
xF3
eJ
eJ3+1);c02
res(sim.yE2);eJ.back().eS2.swap(res);eJ.eA=n41;IP+=2;yA1}
case
cJump:{c02
res(sim.yE2);eJ
e81.swap(res);eJ.eA=eN1[IP+1]+1;IP+=2;yA1}
case
cImmed:xG
Immed[DP++]);lC
cDup:sim.Dup();lC
cNop:lC
cFCall:{x83
funcno=n43;assert(funcno<fpdata.mFuncPtrs.size());x83
params=x33.mFuncPtrs
n53
mParams;sim.EatFunc(params,nB2,funcno
e02
case
cPCall:{x83
funcno=n43;assert(funcno<fpdata.iI3.size());const
FunctionParserBase
xF&p=*x33.iI3
n53
mParserPtr;x83
params=x33.iI3
n53
mParams;x5
paramlist=sim.Pop(tL2;c02
tP2;tP2.GenerateFrom(*p.mData,paramlist)nA2
tP2
e02
case
cInv:xG
1
c32
cDiv);lC
cNeg
nC2
cNeg);break;xG
0
c32
cSub);lC
cSqr:xG
2
xN1
cSqrt:xG
yF2
xN1
cRSqrt:xG
e62-0.5)xN1
cCbrt:xG
e62
1)/e62
3)xN1
cDeg:xG
fp_const_rad_to_deg
xF())yG1
cRad:xG
fp_const_deg_to_rad
xF())yG1
cExp:iR)goto
yM3;xG
fp_const_e
xF()c32
cPow);lC
cExp2:iR)goto
yM3;xG
2.0
c32
cPow);lC
cCot
nC2
cTan);iR)x1
cCsc
nC2
cSin);iR)x1
cSec
nC2
cCos);iR)x1
cInt:
#ifndef __x86_64
iR)n23
1,cInt
e02
#endif
xG
yF2)tQ2
yR
1,cFloor);lC
cLog10
nC2
c53
fp_const_log10inv
xF())yG1
cLog2
nC2
c53
fp_const_log2inv
xF())yG1
cQ3:n33
c53
fp_const_log2inv
xF())yR
3,cMul);lC
cHypot:xG
2
iU3);tS3
xG
2
iU3)tQ2;xG
yF2
xN1
cSinCos:sim.Dup()yR
1,cSin);n33
cCos);lC
cSinhCosh:sim.Dup()yR
1,cSinh);n33
cCosh);lC
cRSub:tS3
case
cSub:iR)n23
2,cSub
e02
xG-1)yR
2,cMul)tQ2;lC
cRDiv:tS3
case
cDiv:iR||IsIntType
xF::nS3)n23
2,cDiv
e02
xG-1
iU3)yG1
cAdd:case
cMul:case
cMod:case
cPow:case
cEqual:case
cLess:case
cGreater:case
i51:case
cLessOrEq:case
cGreaterOrEq:case
cAnd:case
cOr:case
cAbsAnd:case
cAbsOr:sim.Eat(2,tB1
lC
cNot:case
cNotNot:case
cC3:case
cAbsNotNot
nC2
tB1
lC
cFetch:sim.Fetch(n43);lC
cPopNMov:{x83
stackOffs_target=n43;x83
stackOffs_source=n43;sim.PopNMov(stackOffs_target,stackOffs_source
e02
yF3
yM3:;x83
funcno=opcode-cAbs;assert(funcno<FUNC_AMOUNT);const
FuncDefinition&func=Functions[funcno]yR
func.params,tB1
xY3}
}
Become(sim.yE2);
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"Produced tree:\n"
;l21
#endif
}
}
#endif
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
#include <assert.h>
#define FP_MUL_COMBINE_EXPONENTS
l33{using
l33
FUNCTIONPARSERTYPES;using
t5
yO1
static
void
AdoptChildrenWithSameOpcode(eR{
#ifdef DEBUG_SUBSTITUTIONS
bool
nT2
eW3
#endif
for
l41
if(xW
a)nC==t72){
#ifdef DEBUG_SUBSTITUTIONS
if(!nT2)cH2"Before assimilation: "
cR
nT2=true;}
#endif
tree.AddParamsMove(xW
a).GetUniqueRef().lJ2),a);}
#ifdef DEBUG_SUBSTITUTIONS
if(nT2)cH2"After assimilation:   "
cR}
#endif
}
}
t5{xL1
ConstantFolding(eR{tree.Sort();
#ifdef DEBUG_SUBSTITUTIONS
void*c63=0;cV2<<"["
<<(&c63)<<"]Runs ConstantFolding for: "
cR
DumpHashes(tree
lE2
std::flush;
#endif
if(false){redo:;tree.Sort();
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"["
<<(&c63)<<"]Re-runs ConstantFolding: "
cR
DumpHashes(tree);
#endif
}
if(t72!=cImmed){yD3
p=CalculateResultBoundaries(tree);if(p
yI&&p
e61&&p
yM==p
tG1{nU
p
yM);nD}
if(false){ReplaceTreeWithOne:nU
e62
1));goto
do_return;ReplaceTreeWithZero:nU
xG1;goto
do_return;ReplaceTreeWithParam0:
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"Before replace: "
;cV2<<std::hex<<'['
eJ1
hash1<<','
eJ1
hash2<<']'<<std::dec
cR
#endif
tree.eF2
0));
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"After replace: "
;cV2<<std::hex<<'['
eJ1
hash1<<','
eJ1
hash2<<']'<<std::dec
cR
#endif
cJ
c73(t72){case
cImmed:lC
iU2:lC
cAnd:case
cAbsAnd
cK
bool
c5
eW3
for
l41{if(!c83
a)))c5=true;cA3
a),t72==cAbsAnd)){case
xX3
t7
iF2:iZ);lC
lZ1
c73
l51){case
0:iD
1:tD2
t72==cAnd?cNotNot:cAbsNotNot);cJ
yF3
cZ2
cAnd||!c5)if(ConstantFolding_AndLogic(tB2
e21
cOr:case
cAbsOr
cK
bool
c5
eW3
for
l41{if(!c83
a)))c5=true;cA3
a),t72==cAbsOr)tN1
iD
l03
iZ);lC
lZ1
c73
l51){case
0
t7
1:tD2
t72==cOr?cNotNot:cAbsNotNot);cJ
yF3
cZ2
cOr||!c5)if(ConstantFolding_OrLogic(tB2
e21
cNot:case
cC3:{x83
n61
0;e23
c93){case
cEqual:n61
i51;lC
i51:n61
cEqual;lC
cLess:n61
cGreaterOrEq;lC
cGreater:n61
cLessOrEq;lC
cLessOrEq:n61
cGreater;lC
cGreaterOrEq:n61
cLess;lC
cNotNot:n61
cNot;lC
cNot:n61
cNotNot;lC
cC3:n61
cAbsNotNot;lC
cAbsNotNot:n61
cC3;break;yF3
xY3
if(opposite){tD2
OPCODE(opposite)yB
SetParamsMove(xW
0).GetUniqueRef().lJ2));cJ
c73(lY1
0)cU3
cQ1){case
iF2
t7
l03
iD
lZ1
cZ2
cNot&&GetPositivityInfo(xW
0))==iF2)tD2
cC3);if(c93==cIf||c93==cAbsIf
nQ
iZ2=xW
0);iT1&ifp1=iZ2
l8
1);iT1&ifp2=iZ2
l8
2);if(ifp1
nC
iY3
ifp1
cQ1
cD3
ifp1
nC==cNot?cNotNot:cAbsNotNot);tR2
c43
n63)cE3
nO2;p2
tT3
iI
if(ifp2
nC
iY3
ifp2
cQ1
cD3
t72);tR2);n63)cE3
iH
ifp2
nC==cNot?cNotNot:cAbsNotNot);p2
tT3
l8
0)iI
e21
cNotNot:case
cAbsNotNot:{if(c83
0)))n73
cA3
0),t72==cAbsNotNot)){case
xX3
t7
iF2:iD
lZ1
cZ2
cNotNot&&GetPositivityInfo(xW
0))==iF2)tD2
cAbsNotNot);if(c93==cIf||c93==cAbsIf
nQ
iZ2=xW
0);iT1&ifp1=iZ2
l8
1);iT1&ifp2=iZ2
l8
2);if(ifp1
nC
iY3
ifp1
cQ1{tree.SetParam(0,iZ2
l8
0)yB
AddParam(ifp1)cE3
nO2;p2
tT3
iI
if(ifp2
nC
iY3
ifp2
cQ1
cD3
t72);tR2);n63
yB
AddParam(ifp2
yB
eL3
iZ2
nC);cJ}
e21
cIf:case
cAbsIf:{if(ConstantFolding_IfOperations(tB2
xY3
case
cMul:{NowWeAreMulGroup:;AdoptChildrenWithSameOpcode(tree)iK2
nH1=e62
1);size_t
lC2=0;bool
nI1
eW3
x51{if(!xW
tC1
yA1
Value_t
immed=xW
a)x41;if(immed==xG1
goto
ReplaceTreeWithZero;nH1*=immed;++lC2;}
if(lC2>1||(lC2==1&&lJ3
nH1,e62
1))))nI1=true;if(nI1){
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"cMul: Will add new "
iJ3
nH1<<"\n"
;
#endif
for
l41
if(xW
tC1{
#ifdef DEBUG_SUBSTITUTIONS
cV2<<" - For that, deleting "
iJ3
xW
a)x41;cV2<<"\n"
;
#endif
tU3(!lJ3
nH1,e62
1)))tree
cS
eA1
nH1));c73
l51){case
0:iD
1:n73
yF3
if(ConstantFolding_MulGrouping(tB2
if(ConstantFolding_MulLogicItems(tB2
e21
cAdd
cK
Value_t
nD2=0.0;size_t
lC2=0;bool
nI1
eW3
x51{if(!xW
tC1
yA1
Value_t
immed=xW
a)x41;nD2+=immed;++lC2;}
if(lC2>1||(lC2==1&&nD2==xG1)nI1=true;if(nI1){
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"cAdd: Will add new "
iJ3
nD2<<"\n"
;cV2<<"In: "
cR
#endif
for
l41
if(xW
tC1{
#ifdef DEBUG_SUBSTITUTIONS
cV2<<" - For that, deleting "
iJ3
xW
a)x41;cV2<<"\n"
;
#endif
tU3(!(nD2==e62
0.0)))tree
cS
eA1
nD2));c73
l51){case
0
t7
1:n73
yF3
if(ConstantFolding_AddGrouping(tB2
if(ConstantFolding_AddLogicItems(tB2
e21
cMin
cK
size_t
yI2=0;yD3
e4;x51{while(a+1<tree.GetParamCount()&&xW
a)xI
xW
a+1)))iZ+1);yE3
max.cF3(!e4
e61||(p
tG1<e4
tG1){e4
nL3=p
nL3;e4
e61=true;yI2=a;}
}
if(e4
e61)for
l41{yE3
min.cF3
a!=yI2&&p
yM>=e4
tG1
tU3
l51==1){n73
e21
cMax
cK
size_t
yI2=0;yD3
eY;x51{while(a+1<tree.GetParamCount()&&xW
a)xI
xW
a+1)))iZ+1);yE3
min.cF3(!eY
yI||p
yM>eY
yM)){eY
yM=p
yM;eY
yI=true;yI2=a;}
}
if(eY
yI){for
l41{yE3
max.cF3
a!=yI2&&(p
tG1<eY
yM){iZ);}
}
}
if
l51==1){n73
e21
cEqual:case
i51:case
cLess:case
cGreater:case
cLessOrEq:case
cGreaterOrEq:if(ConstantFolding_Comparison(tB2
lC
cAbs:{yD3
tY
xW
0));if
iE
n73
if
iF{tD2
cMul
yB
yA
e62
1)));goto
NowWeAreMulGroup;}
if(c93==cMul){iT1&p=xW
0);eK
n83;eK
cP2
e73
0;a<p.l91++a){tY
p
n93
if
iE{n83
c52
p
n93}
if
iF{cP2
c52
p
n93}
}
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"Abs: mul group has "
<<n83
eJ3<<" pos, "
<<cP2
eJ3<<"neg\n"
;
#endif
if(!n83
eU3||!cP2
eU3){
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"AbsReplace-Before: "
iL2
tree
lE2"\n"
<<std::flush;DumpHashes
nV2);
#endif
yZ2
e13;e13
iH
cMul
iY1
0;a<p.l91++a){tY
p
n93
if(iE||iF){}
else
e13
cS
p
n93}
e13
iW1
tX
nB3;nB3
iH
cAbs);nB3
c91
e13);nB3
iW1
tX
iG
cMul);mulgroup
c91
nB3);cA1
AddParamsMove(n83);if(!cP2
eU3){if(cP2
eJ3%2)cA1
yA
nU2));cA1
AddParamsMove(cP2);}
tree.Become
eH2);
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"AbsReplace-After: "
;DumpTree
nV2
lE2"\n"
<<std::flush;DumpHashes
nV2);
#endif
goto
NowWeAreMulGroup;}
}
xY3
#define HANDLE_UNARY_CONST_FUNC(funcname) nP){nU funcname(lR));nD
case
cLog:i83(fp_log);if(c93==cPow
nQ
pow=xW
0);if(GetPositivityInfo(pow
l8
0))==iF2)iV1
t2
yB
lU
if(GetEvennessInfo(pow
l8
1))==iF2)iV1()tX
abs;abs
iH
cAbs);abs
c91
pow
c43
abs.Rehash
t2);pow.nB1
0,abs
yB
lU}
iP1
c93==cAbs
nQ
pow=xW
0)l8
0);if(pow
nC==cPow)iV1()tX
abs;abs
iH
cAbs);abs
c91
pow
c43
abs.Rehash
t2);pow.nB1
0,abs
yB
lU}
lC
cAcosh:i83(fp_acosh);lC
cAsinh:i83(fp_asinh);lC
cAtanh:i83(fp_atanh);lC
cAcos:i83(fp_acos);lC
cAsin:i83(fp_asin);lC
cAtan:i83(fp_atan);lC
cCosh:i83(fp_cosh);lC
cSinh:i83(fp_sinh);lC
cTanh:i83(fp_tanh);lC
cSin:i83(fp_sin);lC
cCos:i83(fp_cos);lC
cTan:i83(fp_tan);lC
cCeil:if(n5
i83(fp_ceil);lC
cTrunc:if(n5
i83(fp_trunc);lC
cFloor:if(n5
i83(fp_floor);lC
cInt:if(n5
i83(fp_int);lC
cCbrt:i83(fp_cbrt);lC
cSqrt:i83(fp_sqrt);lC
cExp:i83(fp_exp);lC
cLog2:i83(fp_log2);lC
cLog10:i83(fp_log10);lC
cQ3:if
lI
fp_log2(lR)*xW
nC3
cArg:i83(fp_arg);lC
cConj:i83(fp_conj);lC
cImag:i83(fp_imag);lC
cReal:i83(fp_real);lC
cPolar:if
lI
fp_polar
e51
lC
cMod:if
lI
fp_mod
e51
lC
cAtan2:{yD3
tY
xW
yP2
p1=tZ
1));nP&&lJ3
lR,xG1){if(p1
e61&&(p1
tG1<xG1{nU
fp_const_pi
cH3
if(p1
yI&&p1
yM>=tE1
xG1;nD}
if(eI1
lJ3
xW
i8,xG1){if(p0
e61&&(p0
tG1<xG1{nU-fp_const_pihalf
cH3
if
cN3
p0
yM>xG1{nU
fp_const_pihalf
cH3}
if
lI
fp_atan2
e51
if((p1
yI&&p1
yM>xG1||(p1
e61&&(p1
tG1<fp_const_negativezero
xF())nQ
yJ2;yJ2
iH
cPow);yJ2
c91
xW
1));yJ2.yA
nU2));yJ2
iW1
tX
yK2;yK2
iH
cMul);yK2
c91
xW
0));yK2
c91
yJ2);yK2
cY
cAtan
yB
nB1
0,yK2
yB
i91
1);e21
cPow:{if(ConstantFolding_PowOperations(tB2
xY3
case
cDiv:nP&&eI1
xW
i8!=tE1
lR/xW
nC3
cInv:nP&&lR!=tE1
e62
1)/lR);nD
lC
cSub:if
lI
lR-xW
nC3
cNeg:nP){nU-lR);nD
lC
cRad:nP){nU
RadiansToDegrees
eG2
cDeg:nP){nU
DegreesToRadians
eG2
cSqr:nP){nU
lR*lR);nD
lC
cExp2:i83(fp_exp2);lC
cRSqrt:nP){nU
e62
1)/fp_sqrt
eG2
cCot:c22
fp_tan(lZ
cSec:c22
fp_cos(lZ
cCsc:c22
fp_sin(lZ
cHypot:if
lI
fp_hypot
e51
lC
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
cFCall:xY3
do_return:;
#ifdef DEBUG_SUBSTITUTIONS
cV2<<"["
<<(&c63)<<"]Done ConstantFolding, result: "
cR
DumpHashes(tree);
#endif
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
t5{xL1
yD3::set_abs(nL
bool
has_negative=!min.known||min.val<e62);bool
has_positive=!cG3||max.val>e62);bool
crosses_axis=has_negative&&has_positive;cI1
xF
newmax;if(min.cF3
cG3)newmax.set(fp_max(i61,i71);if(crosses_axis
cI3
e62));cJ3
min.cF3
cG3
cI3
fp_min(i61,i71);iP1
min.known
cI3
i61);else
min.set(i71;}
max=newmax;}
xL1
yD3::set_neg(){std::swap(min,max);min.val=-min.val;max.val=-max.val;}
yP1
IsLogicalTrueValue
iS1
yD3&p,bool
abs){if(nB
IsIntType
xF::nS3){if(p
yI&&p
yM>=e62
1))lD2
if(!abs&&p
e61
t43<=nU2)lD2}
cJ3
p
yI&&p
yM>=yF2)lD2
if(!abs&&p
e61
t43<=e62-0.5))lD2
eE3
t23
yP1
IsLogicalFalseValue
iS1
yD3&p,bool
abs){if(nB
IsIntType
xF::nS3){if(abs)eG3
e61
l61
1);else
eG3
yI&&p
e61&&p
yM>nU2
l61
1);}
cJ3
abs)eG3
e61
l61
0.5);else
eG3
yI&&p
e61&&p
yM>e62-0.5)l61
0.5);}
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
l33
FUNCTIONPARSERTYPES;using
t5;t5{tK1
yD3
CalculateResultBoundaries
iS1
eR
#ifdef DEBUG_SUBSTITUTIONS_extra_verbose
{using
l33
FUNCTIONPARSERTYPES;yD3
tmp=CalculateResultBoundaries_do(tree
lE2"Estimated boundaries: "
;if(tmp
yI)cV2<<tmp
yM;else
cV2<<"-inf"
;cV2<<" .. "
;if(tmp
e61)cV2<<tmp
nL3;else
cV2<<"+inf"
;cV2<<": "
iL2
tree
lE2
std::endl
i1
tmp;}
tK1
yD3
yZ2::CalculateResultBoundaries_do
iS1
eR
#endif
{iS
yH1(-fp_const_pihalf
xF(),fp_const_pihalf
xF());iS
pi_limits(-fp_const_pi
xF(),fp_const_pi
xF());iS
abs_pi_limits(xH1,fp_const_pi
xF());iS
plusminus1_limits(nU2,e62
1));using
l33
std;e23
t72){case
cImmed:nM
tree
x41
cU3
x41)t83
cAnd:case
cAbsAnd:case
cOr:case
cAbsOr:case
cNot:case
cC3:case
cNotNot:case
cAbsNotNot:case
cEqual:case
i51:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:{nM
xH1,e62
1));}
case
cAbs:tW3
set_abs();tK
cLog:tW3
tS2
fp_log);nD3
fp_log);tK
cLog2:tW3
tS2
fp_log2);nD3
fp_log2);tK
cLog10:tW3
tS2
fp_log10);nD3
fp_log10);tK
cAcosh:tW3
min.tL1
cGreaterOrEq
tH1
fp_acosh
tR3
cGreaterOrEq
tH1
fp_acosh);tK
cAsinh:tW3
min.set(fp_asinh);tO
set(fp_asinh);tK
cAtanh:tW3
min.n3-1),fp_atanh
tR3
cLess
tH1
fp_atanh);tK
cAcos:lD
nM(tO
cF3(tO
val)<e62
1))?fp_acos(tO
val):xH1,(m
yI&&(m
yM)>=nU2)?fp_acos(m
yM):fp_const_pi
xF());}
case
cAsin:tW3
min.n3-1),fp_asin,yH1
yM
tR3
cLess
tH1
fp_asin,yH1
tG1;tK
cAtan:tW3
min.set(fp_atan,yH1
yM);tO
set(fp_atan,yH1
tG1;tK
cAtan2:{nP&&lJ3
lR,xG1)yP
abs_pi_limits;}
if(eI1
lJ3
xW
i8,xG1)yP
yH1;eE3
pi_limits;}
case
cSin:lD
bool
x11=!m
yI||!tO
known||(tO
val-m
yM)>=(yK
x11)tM
Value_t
min=cK3
m
yM,yK
min<xG1
min
yN
Value_t
max=cK3
tO
val,yK
max<xG1
max
yN
if(max<min)max
yN
bool
y21=(min<=fp_const_pihalf
xF()&&max>=fp_const_pihalf
xF());bool
nJ1=(min<=cD&&max>=cD);if(y21&&nJ1)tM
if(nJ1)nM
nU2,nW2
if(y21)nM
yL2
e62
1));nM
yL2
nW2}
case
cCos:lD
if(m
yI)m
yM+=fp_const_pihalf
xF();if(xI2
tO
val+=fp_const_pihalf
xF();bool
x11=!m
yI||!tO
known||(tO
val-m
yM)>=(yK
x11)tM
Value_t
min=cK3
m
yM,yK
min<xG1
min
yN
Value_t
max=cK3
tO
val,yK
max<xG1
max
yN
if(max<min)max
yN
bool
y21=(min<=fp_const_pihalf
xF()&&max>=fp_const_pihalf
xF());bool
nJ1=(min<=cD&&max>=cD);if(y21&&nJ1)tM
if(nJ1)nM
nU2,nW2
if(y21)nM
yL2
e62
1));nM
yL2
nW2}
case
cTan:{nM);}
case
cCeil:lD
tO
i81
cFloor:lD
m
yI1
tK
cTrunc:lD
m
yI1
tO
i81
cInt:lD
m
yI1
tO
i81
cSinh:tW3
min.set(fp_sinh);tO
set(fp_sinh);tK
cTanh:tW3
min.set(fp_tanh,plusminus1_limits.min);tO
set(fp_tanh,plusminus1_limits.max);tK
cCosh:lD
if(m
yI){if(xI2{if(m
yM>=xH1&&tO
val>=xG1{m
yM
xS}
iP1(m
yM)<xH1&&tO
val>=xG1
l14
tmp
xS
if(tmp>tO
val)tO
val=tmp;m
yM=e62
1);}
else{m
yM
xS
std::swap(m
yM,tO
val);}
}
cJ3
m
yM>=xG1{m.tL
m
yM=fp_cosh(m
yM);}
else{m.tL
m
yM=e62
1);}
}
}
else{m
yI=true;m
yM=e62
1);if(xI2{m
yM=fp_cosh(tO
val);m.tL}
else
m.tL}
tK
cIf:case
cAbsIf:{yD3
res1=tZ
1));yD3
res2=tZ
2));if(!res2
yI)res1
yI
eW3
iP1
res1
yI&&(res2
yM)<res1
yM)res1
yM=res2
yM;if(!res2
e61)res1.tL
iP1
res1
e61&&(res2
tG1>res1
tG1
res1
nL3=res2
nL3
i1
res1;}
case
cMin:{bool
iT
eW3
bool
iU
eW3
e83;x7
m
tZ3
m
yI)iT
iN3
yI||(m
yM)<nF3)nF3=m
yM;if(!xI2
iU
iN3
e61||(tO
val)<nG3
nH3=tO
val;}
if(iT)nK3
iU)nQ3
return
nR3
cMax:{bool
iT
eW3
bool
iU
eW3
e83;x7
m
tZ3
m
yI)iT
iN3
yI||m
yM>nF3)nF3=m
yM;if(!xI2
iU
iN3
e61||tO
val>nG3
nH3=tO
val;}
if(iT)nK3
iU)nQ3
return
nR3
cAdd:{e83(xH1,xG1;x7
item=tZ
a));if(item
yI)nF3+=item
yM;else
nK3
item
e61)nH3+=item
nL3;else
nQ3
if(!nI3&&!nJ3)xY3
if(nI3&&nJ3&&nF3>nG3
std::swap(nF3,nG3
i1
nR3
cMul:{e72
Value{enum
nT3{tT2,lG2,PlusInf}
;nT3
eB
iK2
value;Value(nT3
t):eB(t),value(0){}
Value
cX3
v):eB(tT2),value(v){}
bool
cQ2
c01
eB==lG2||(eB==tT2&&value<xG1;}
void
eM1*=iS1
Value&rhs){if(eB==tT2&&rhs.eB==tT2)value*=rhs.value;else
eB=(cQ2)!=rhs.cQ2)?lG2:PlusInf);}
iD2<iS1
Value&rhs
c01(eB==lG2&&rhs.eB!=lG2)||(eB==tT2&&(rhs.eB==PlusInf||(rhs.eB==tT2&&value<rhs.value)));}
}
;e72
yJ1{Value
yN2,yO2;yJ1():yN2(Value::PlusInf),yO2(Value::lG2){}
void
xJ2
Value
cM3,const
Value&value2){cM3*=value2;if(cM3<yN2)yN2=cM3;if(yO2<cM3)yO2=cM3;}
}
;e83(nM3
x7
item
tZ3
item
yI&&!item
e61)nM);Value
iP3=nI3?Value(nF3):xC2
lG2);Value
nU3=nJ3?Value(nG3:xC2
PlusInf);Value
nV3=item
yI?Value(item
yM):xC2
lG2);Value
nW3=item
e61?Value(item
tG1:xC2
PlusInf);yJ1
range;range.xJ2
iP3,nV3);range.xJ2
iP3,nW3);range.xJ2
nU3,nV3);range.xJ2
nU3,nW3);if(range.yN2.eB==Value::tT2)nF3=range.yN2.value;else
nK3
range.yO2.eB==Value::tT2)nH3=range.yO2.value;else
nQ3
if(!nI3&&!nJ3)xY3
if(nI3&&nJ3&&nF3>nG3
std::swap(nF3,nG3
i1
nR3
cMod:{yD3
x=tZ
yP2
y=tZ
1));if(y
e61){if(y
nL3>=xG1{if(!x
yI||(x
yM)<xG1
nM-y
nL3,y
tG1;i03
xH1,y
tG1;}
cJ3!x
e61||(x
tG1>=xG1
nM
y
nL3,-y
tG1;i03
y
nL3,fp_const_negativezero
xF());}
}
i03);}
case
cPow:{if(eI1
xW
i8==xG1{nM
nM3}
nP&&lR==xG1{nM
xH1,xG1;}
nP&&lJ3
lR
nH2
nM
nM3}
if(eI1
xW
i8>xH1&&GetEvennessInfo(xW
1))==iF2)l14
eD2
xW
i8;yD3
tmp=tZ
yP2
nS3;nI3=true;nF3=0;if(tmp
yI&&tmp
yM>=xG1
nF3=t33
tmp
yM,tU2
iP1
tmp
e61&&tmp
nL3<=xG1
nF3=t33
tmp
nL3,tU2
nQ3
if(tmp
yI&&tmp
e61){nJ3=true;nH3=fp_max(fp_abs(tmp
yM),fp_abs(tmp
tG1);nH3=t33
nH3,tU2
eE3
nS3;}
yD3
tY
xW
yP2
p1=tZ
1));TriTruthValue
p0_positivity=cN3(p0
yM)>=xG1?iF2:(p0
e61&&(p0
tG1<xH1?l03
Unknown);TriTruthValue
cR2=GetEvennessInfo(xW
1));TriTruthValue
eZ=Unknown;e23
p0_positivity
tN1
eZ=iF2;lC
l03{eZ=cR2;xY3
yF3
e23
cR2
tN1
eZ=iF2;lC
l03
lC
Unknown:{if(eI1!t42
xW
i8)&&xW
i8>=xG1{eZ=iF2;}
xY3}
c73(eZ
tN1
l14
min=xH1;if
cN3
p1
yI){min=t33
p0
yM,p1
yM);if(p0
yM<xH1&&(!p1
e61||p1
nL3>=xG1&&min>=xG1
min=xH1;}
if
cN3
p0
yM>=xH1&&p0
e61&&p1
e61)l14
max=t33
p0
nL3,p1
tG1;if(min>max)std::swap(min,max);nM
min,max);}
nM
min,false);}
case
l03{nM
false,fp_const_negativezero
xF());}
yF3{xY3
e21
cNeg:tW3
set_neg();tK
cSub:yJ
cNeg);tmp2
tX3
tmp
iH
cAdd);tmp
tY3
tmp
c91
tmp2
i0
cInv:{c02
lW-1))i0
cDiv:yJ
cInv);tmp2
tX3
tmp
iH
xC1
AddParamMove(tmp2
i0
cRad:yV
xC1
yA
fp_const_rad_to_deg
xF())i0
cDeg:yV
xC1
yA
fp_const_deg_to_rad
xF())i0
cSqr:{c02
lW
2))i0
cExp:yV
cPow);tmp.yA
fp_const_e
xF()));tmp.nJ
0)i0
cExp2:yV
cPow);tmp.yA
nX3
tmp.nJ
0)i0
cCbrt:tW3
min.set(fp_cbrt);tO
set(fp_cbrt);tK
cSqrt:lD
if(m
yI)m
yM=(m
yM)<xH1?0:fp_sqrt(m
yM);if(xI2
tO
val=(tO
val)<xH1?0:fp_sqrt(tO
val);tK
cRSqrt:{c02
lW-0.5))i0
cHypot:{yZ2
xsqr,ysqr,add,sqrt;xsqr
tY3
xsqr.yA
nX3
ysqr
tX3
ysqr.yA
nX3
xsqr
iH
cPow);ysqr
iH
cPow);add
c91
xsqr);add
c91
ysqr);add
iH
cAdd);sqrt
c91
add);sqrt
iH
cSqrt)i1
CalculateResultBoundaries(sqrt);}
case
cQ3:yJ
cLog2);tmp2
tY3
tmp
iH
cMul);tmp
c91
tmp2);tmp.nJ
1)i0
cCot:yJ
cTan)x6
lH
cSec:yJ
cCos)x6
lH
cCsc:yJ
cSin)x6
CalculateResultBoundaries(tmp);}
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
iU2:lC
cArg:case
cConj:case
cImag:case
cReal:case
cPolar:lC
cPCall:lC
cFCall:xY3
nM);}
tK1
TriTruthValue
GetIntegerInfo
iS1
eR{e23
t72){case
cImmed:return
t42
tree
x41)?iF2:xX3
t83
cFloor:case
cCeil:case
cTrunc:case
cInt:return
iF2
t83
cAnd:case
cOr:case
cNot:case
cNotNot:case
cEqual:case
i51:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:return
iF2
t83
cIf:{TriTruthValue
a=GetIntegerInfo(xW
1));TriTruthValue
b=GetIntegerInfo(xW
2));if(a==b)return
a
xO2
case
cAdd:case
cMul:{for
l41
if(GetIntegerInfo(xW
a))!=iF2)return
Unknown
i1
iF2;}
yF3
xY3
return
Unknown;}
yP1
IsLogicalValue
iS1
eR{e23
t72){case
cImmed:return
lJ3
tree
x41,xG1||lJ3
tree
x41,e62
1))t83
cAnd:case
cOr:case
cNot:case
cNotNot:case
cAbsAnd:case
cAbsOr:case
cC3:case
cAbsNotNot:case
cEqual:case
i51:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:nZ
cMul:{for
l41
if(!c83
a))cQ
lD2}
case
cIf:case
cAbsIf:yP
c83
1))&&c83
2));}
yF3
xY3
return
t23}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
l33
FUNCTIONPARSERTYPES;
#if defined(__x86_64) || !defined(FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS)
# define CBRT_IS_SLOW
#endif
#if defined(DEBUG_POWI) || defined(DEBUG_SUBSTITUTIONS)
#include <cstdio>
#endif
l33
xK1{extern
const
x83
char
powi_table[256];}
l33{using
t5
yO1
bool
IsOptimizableUsingPowi(long
immed,long
penalty=0){xK1::eP3
synth;synth.PushVar(iU2);size_t
bytecodesize_backup=synth.GetByteCodeSize();xK1::x01
immed,xK1::tU1
xF::MulSequence,c92
size_t
bytecode_grow_amount=synth.GetByteCodeSize()-bytecodesize_backup
i1
bytecode_grow_amount<size_t(MAX_POWI_BYTECODE_LENGTH-penalty);}
xL1
ChangeIntoRootChain(nP1
bool
l13,long
tY2,long
tZ2){while(tZ2>0
t82
cCbrt);t92);tmp.e33
xD2--tZ2;}
while(tY2>0
t82
cSqrt);if(l13){tmp
iH
cRSqrt);l13
eW3}
t92);tmp.e33
xD2--tY2;}
if(l13
t82
cInv);t92
xD2}
}
tK1
e72
RootPowerTable{static
const
Value_t
RootPowers[(1+4)*(1+3)];}
yO1
const
Value_t
t9(1+4)*(1+3)]={e62
1)lS
2)lS
2*2)lS
2*2*2)lS
2*2*2*2)lS
3)lS
3*2)lS
3*2*2)lS
3*2*2*2)lS
3*2*2*2*2)lS
3*3)lS
3*3*2
yG3
yG3*2
yG3*2*2)lS
3*3*3
eY2
eY2*2
eY2*2*2
eY2*2*2*2)}
;e72
PowiResolver{static
const
x83
MaxSep=4;static
x03
i13=5;typedef
int
cO3;typedef
long
x43;typedef
long
tP;e72
yS2{yS2():n_int_sqrt(0),n_int_cbrt(0),sep_list(),n01(0){}
int
n_int_sqrt;int
n_int_cbrt;int
tT1
MaxSep];tP
n01;}
yO1
static
yS2
CreatePowiResult
cX3
exponent){yS2
nS3;cO3
tA=FindIntegerFactor(tU2
if(tA==0){
#ifdef DEBUG_POWI
i02"no factor found for %Lg\n"
,(cP3);
#endif
return
nS3;}
nS3.n01=y31
exponent,tA);x43
eT2=EvaluateFactorCost(tA,0,0,0)+cC
nS3.n01);int
s_count=0;int
i23=0;int
nZ3=0;
#ifdef DEBUG_POWI
i02"orig = %Lg\n"
,(cP3);i02"plain factor = "
iL3"%ld\n"
,(int)tA,(long)eT2);
#endif
for
lT1
n_s=0;n_s<MaxSep;++n_s){int
xD=0;x43
yK1=eT2;cO3
yX1=tA;for(int
s=1;s<i13*4;++s){
#ifdef CBRT_IS_SLOW
if(s>=i13)break;
#endif
int
n_sqrt=s%i13;int
n_cbrt=s/i13;if(n_sqrt+n_cbrt>4)yA1
Value_t
lN1=exponent;lN1-=t9
s];iE1=FindIntegerFactor(lN1);if(xT2!=0){tP
xN=y31
lN1,xT2);x43
cost=EvaluateFactorCost(xT2,s_count+n_sqrt,i23+n_cbrt,nZ3+1)+cC
xN);
#ifdef DEBUG_POWI
i02"Candidate sep %u (%d*sqrt %d*cbrt)factor = "
iL3"%ld (for %Lg to %ld)\n"
,s,n_sqrt,n_cbrt,xT2,(long)cost
cX2
lN1,(long)xN);
#endif
if(cost<yK1){xD=s;yX1=xT2;yK1=cost;}
}
}
if(!xD)break;
#ifdef DEBUG_POWI
i02"CHOSEN sep %u (%d*sqrt %d*cbrt)factor = "
iL3"%ld, exponent %Lg->%Lg\n"
,xD,xD%i13,xD/i13,yX1,yK1
cX2(exponent)cX2(exponent-t9
xD]));
#endif
nS3.tT1
n_s]=xD
eX3-=t9
xD];s_count+=xD%i13;i23+=xD/i13;eT2=yK1;tA=yX1;nZ3+=1;}
nS3.n01=y31
exponent,tA);
#ifdef DEBUG_POWI
i02"resulting exponent is %ld (from exponent=%Lg, best_factor=%Lg)\n"
,nS3.n01,(cP3
cX2
tA);
#endif
while(tA%2==0){++nS3
e52;tA/=2;}
while(tA%3==0){++nS3.n_int_cbrt;tA/=3;eE3
nS3;}
private:static
x43
cC
tP
xN){static
std::map
cS2
iM;if(xN<0){x43
cost=22
i1
cost+cC-xN);}
std::map
cS2::y83
i=iM.xU2
xN);if(i!=iM.cX1
xN)return
i
eC2;std::pair
cS2
nS3(xN,0.0);x43&cost=nS3.i63
while(xN>1){int
xT2=0;if(xN<256){xT2=xK1::powi_table[xN];if(xT2&128)xT2&=127;else
xT2=0;if(xT2&64)xT2=-(xT2&63)-1;}
if(xT2){cost+=cC
xT2);xN/=xT2;yA1}
if(!(xN&1)){xN/=2;cost+=6;}
else{cost+=7;xN-=1;}
}
iM.y43,nS3)i1
cost;}
cB1
tP
y31
yL1,iE1)yP
makeLongInteger(value*e62
xT2));}
cB1
bool
yM1
yL1,iE1)l14
v=value*e62
xT2)i1
isLongInteger(v);}
cB1
cO3
FindIntegerFactor(yL1){iE1=(2*2*2*2);
#ifdef CBRT_IS_SLOW
#else
xT2*=(3*3*3);
#endif
cO3
nS3=0;if(yM1
value,xT2)){nS3=xT2;while((xT2%2)==0&&yM1
value,xT2/2))nS3=xT2/=2;while((xT2%3)==0&&yM1
value,xT2/3))nS3=xT2/=3;}
#ifdef CBRT_IS_SLOW
if(nS3==0){if(yM1
value,3
y03
3;}
#endif
return
nS3;}
static
int
EvaluateFactorCost(int
xT2,int
s,int
c,int
nmuls){x03
x13=6;
#ifdef CBRT_IS_SLOW
x03
eU2=25;
#else
x03
eU2=8;
#endif
int
nS3=s*x13+c*eU2;while(xT2%2==0){xT2/=2;nS3+=x13;}
while(xT2%3==0){xT2/=3;nS3+=eU2;}
nS3+=nmuls
i1
nS3;}
}
;}
t5{yP1
yZ2::RecreateInversionsAndNegations(bool
prefer_base2){bool
changed=false
e73
0;a<l91++a)if(lX1.RecreateInversionsAndNegations(prefer_base2))yS1
if(changed){exit_changed:Mark_Incompletely_Hashed(yQ2
e23
l92{case
cMul:{eK
nE2
tX
nF2,cR1;if(true){bool
nK1=false
iK2
xE2=0
e73
tI1
i12
tI
0)i22
tI
1)tD1){nK1=true;xE2=tI
i8;xY3}
if(nK1)l14
immeds=1.0
e73
tI1
tD1){immeds*=powgroup
x41;yN1}
for
iZ1
a=l91
a-->0;nQ&powgroup=lX1;if(powgroup
i12
tI
0)i22
tI
1)tD1
nQ&log2=tI
0);log2.lD1
log2
iH
cQ3);log2.yA
t33
immeds,e62
1)/xE2)));log2
iW1;xY3}
}
}
for
iZ1
a=tI1
i12
tI
1)tD1){iT1&exp_param=tI
1)iK2
eD2
exp_param
x41;if(cZ1,nU2)){lD1
nE2
c52
lX1
c43
yN1
iP1
exponent<xH1&&t42
exponent)nQ
iV;iV
iH
cPow);iV
cS
tI
0));iV.yA-exponent));iV
iW1;nE2
c52
iV);lD1
yN1}
iP1
powgroup
i22!nF2.nX2){nF2=tI
0);lD1
yN1
iP1
powgroup
nC==cQ3&&!cR1.nX2){cR1=powgroup;lD1
yN1}
if(!nE2
eU3){changed=true
tX
cF1;cF1
iH
cMul);cF1
iF1
nE2);cF1
iW1
tX
iG
cMul);cA1
SetParamsMove
eM
if(cA1
IsImmed()&&fp_equal
eH2
x41
nH2
eL3
cInv)tN
cF1);}
cJ3
cA1
GetDepth()>=cF1
nP2
tC2
cDiv
iG1
tN
cF1);}
else{eL3
cRDiv)tN
cF1
iG1;}
}
}
if(nF2.nX2
nQ
iG
l92;cA1
SetParamsMove
eM
while(cA1
RecreateInversionsAndNegations(prefer_base2))cA1
tE3;eL3
cQ3)tN
nF2
iG1;yS1}
if(cR1.nX2
nQ
iG
cMul);mulgroup
c91
cR1
l8
1));cA1
AddParamsMove
eM
while(cA1
RecreateInversionsAndNegations(prefer_base2))cA1
tE3;DelParams();eL3
cQ3)tN
cR1
l8
0)iG1;yS1
e21
cAdd:{eK
i32
e73
l91
a-->0;)if(eY3
cMul){nG2
y41:tX&mulgroup=cT2
for
iZ1
b=cA1
l91
b-->0;){if
eH2
l8
b).lV1
xT2=mulgroup
l8
b)x41;lG3
xT2
xU
y41;}
cA1
lD1
cA1
i91
b
tA2
iP1
lJ3
xT2,e62-2)))eN
y41;}
cA1
lD1
cA1
i91
b);cA1
yA
e62
2))tA2}
}
if(t0){cA1
tB
mulgroup);yN1}
iP1
eY3
cDiv&&!IsIntType
xF::nS3){nG2
y51:tX&cF1=cT2
if(cF1
l8
0)tD1){lG3
cF1
l8
0)x41
xU
y51;}
cF1.lD1
cF1.i91
0);cF1
iH
cInv
tA2}
if(t0)eN
y51;}
cF1.tB
cF1);yN1}
iP1
eY3
cRDiv&&!IsIntType
xF::nS3){nG2
xD1:tX&cF1=cT2
if(cF1
l8
1)tD1){lG3
cF1
l8
i8
xU
xD1;}
cF1.lD1
cF1.i91
1);cF1
iH
cInv
tA2}
if(t0)eN
xD1;}
cF1.tB
cF1);yN1}
if(!i32
eU3){
#ifdef DEBUG_SUBSTITUTIONS
i02"Will make a Sub conversion in:\n"
);fflush(stdout);l21
#endif
yZ2
yT2;yT2
iH
cAdd);yT2
iF1
i32);yT2
iW1
tX
cS1;cS1
iH
cAdd);cS1
iF1
lJ2));cS1
iW1;if(cS1
tD1&&lJ3
cS1
x41,xG1
tC2
cNeg);eC);}
cJ3
cS1
nP2==1
tC2
cRSub);eC)eV3}
iP1
yT2
nC==cAdd
tC2
cSub)eV3
eC
l8
0)iY1
1;a<yT2.l91++a
nQ
eV2;eV2
iH
cSub);eV2
iF1
lJ2));eV2.yC2)tN
eV2);eC
n93}
}
else{eL3
cSub)eV3
eC);}
}
#ifdef DEBUG_SUBSTITUTIONS
i02"After Sub conversion:\n"
);fflush(stdout);l21
#endif
e21
cPow:{iT1&p0
i42
0);iT1&p1
i42
1);if(p1
tD1){if(p1
x41!=xH1&&!t42
p1
x41)){eG
yS2
r=eG
CreatePowiResult(fp_abs
cR3);if(r.n01!=0){bool
lH2
eW3
if(p1
x41<xH1&&r.tT1
0]==0&&r
e52>0){lH2=true;}
#ifdef DEBUG_POWI
i02"Will resolve powi %Lg as powi(chain(%d,%d),%ld)"
cX2
fp_abs
cR3,r
e52,r.n_int_cbrt,r.n01);for
lT1
n=0;n<eG
MaxSep;++n){if(r
i52==0)break;int
n_sqrt=r
i52%eG
i13;int
n_cbrt=r
i52/eG
i13;i02"*chain(%d,%d)"
,n_sqrt,n_cbrt);}
i02"\n"
);
#endif
yZ2
cU2
i42
0)tX
yU2=cU2;yU2.lD1
ChangeIntoRootChain(yU2,lH2,r
e52,r.n_int_cbrt);yU2
iW1
tX
pow;if(r.n01!=1){pow
iH
cPow);pow
c91
yU2);pow.yA
e62
r.n01)));}
else
pow.swap(yU2)tX
mul;mul
iH
cMul);mul
c91
pow);for
lT1
n=0;n<eG
MaxSep;++n){if(r
i52==0)break;int
n_sqrt=r
i52%eG
i13;int
n_cbrt=r
i52/eG
i13
tX
eW2=cU2;eW2.lD1
ChangeIntoRootChain(eW2,false,n_sqrt,n_cbrt);eW2
iW1;mul
c91
eW2);}
if(p1
x41<xH1&&!lH2){mul
iW1;eL3
cInv);nB1
0,mul);i91
1);}
else{eL3
cMul);SetParamsMove(mul.lJ2));}
#ifdef DEBUG_POWI
l21
#endif
yS1
xY3}
}
if(GetOpcode()==cPow&&(!p1
tD1||!isLongInteger
cR3||!IsOptimizableUsingPowi
xF(makeLongInteger
cR3))){if(p0
tD1&&p0
x41>e62
0.0)){if(prefer_base2)l14
yV2=fp_log2(p0
x41);lG3
yV2
nH2
i91
0);}
else{n0
eA1
yV2))eX3
cS
p1
t81
nB1
t1}
eL3
cExp2);yS1}
else
l14
yV2=fp_log(p0
x41);lG3
yV2
nH2
i91
0);}
else{n0
eA1
yV2))eX3
cS
p1
t81
nB1
t1}
eL3
cExp);yS1}
}
iP1
GetPositivityInfo(p0)==iF2){if(prefer_base2
nQ
log;log
iH
cLog2);log
cS
p0);log
iW1;n0
p1)eX3
c91
log
t81
eL3
cExp2);nB1
t1
yS1}
else{yZ2
log;log
iH
cLog);log
cS
p0);log
iW1;n0
p1)eX3
c91
log
t81
eL3
cExp);nB1
t1
yS1}
}
e21
cDiv:{if(GetParam(0)tD1&&lJ3
GetParam(0)x41
nH2
eL3
cInv);i91
0);}
xY3
yF3
xY3
if(changed)goto
exit_changed
i1
changed;}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
l33
FUNCTIONPARSERTYPES;l33{using
t5;class
eF3{size_t
nL1;size_t
eH;size_t
eI;size_t
lO1;size_t
t3;size_t
t4;size_t
n71;e43
eF3():nL1(0),eH(0),eI(0),lO1(0),t3(0),t4(0),n71(0){}
void
i43
OPCODE
op){nL1+=1
i33
cCos)++eH
i33
cSin)++eI
i33
cSec)++eH
i33
cCsc)++eI
i33
cTan)++lO1
i33
cCot)++lO1
i33
cSinh)++t4
i33
cCosh)++t3
i33
cTanh)++n71;}
size_t
GetCSEscore()const{size_t
nS3=nL1
i1
nS3;}
int
NeedsSinCos()const{bool
y61=(nL1==(eH+eI+lO1));if((lO1&&(eI||eH))||(eI&&eH)){if(y61)return
1
i1
2;eE3
0;}
int
NeedsSinhCosh()const{bool
y61=(nL1==(t3+t4+n71));if((n71&&(t4||t3))||(t4&&t3)){if(y61)return
1
i1
2;eE3
0;}
size_t
MinimumDepth()const{size_t
n_sincos=std::min(eH,eI);size_t
n_sinhcosh=std::min(t3,t4);if(n_sincos==0&&n_sinhcosh==0)return
2
i1
1;}
}
yO1
class
TreeCountType:public
std::multimap<fphash_t,std::pair<eF3,yZ2> >{}
;xL1
FindTreeCounts(tJ1&nI2,const
nP1
OPCODE
xF2,bool
skip_root=false){cX
i=nI2.xU2
nJ2);if(!skip_root){bool
found
eW3
for(;i!=nI2.cX1
nJ2;++i){if(tree
xI
i
eC2.second)){i
eC2.first.i43
xF2);found=true;xY3}
if(!found){eF3
count;count.i43
xF2);nI2.y43,std::make_pair(nJ2,std::make_pair
cW3
cU3)));}
}
x51
FindTreeCounts(nI2,xW
a),t72);}
e72
yW{bool
BalanceGood;bool
FoundChild;}
yO1
yW
lP1
iT1&root,iT1&cS3{if(root
xI
cS3){yW
nS3={true,true}
i1
nS3;}
yW
nS3={true,false}
;if(root
nC==cIf||root
nC==t03{yW
cond=lP1
root
l8
0),cS3;yW
xV=lP1
root
l8
1),cS3;yW
y4=lP1
root
l8
2),cS3;if
l23||xV
yX||y4
yX){nS3
yX=true;}
nS3
eD=((xV
yX==y4
yX)||l23
i62)&&(cond
eD||(xV
yX&&y4
yX))&&(xV
eD||l23
i62)&&(y4
eD||l23
i62);}
else{bool
iH1
eW3
bool
nM1=false
eK3
b=root.GetParamCount(),a=0;a<b;++a){yW
tmp=lP1
root
l8
a),cS3;if(tmp
yX)nS3
yX=true;if(tmp
eD==false)iH1=true;iP1
tmp
yX)nM1=true;}
if(iH1&&!nM1)nS3
eD
eW3
eE3
nS3;}
yP1
cW2
iR1
i53
const
nP1
const
xK1::eP3&synth,const
tJ1&nI2){for
iZ1
b=tree.GetParamCount(),a=0;a<b;++a){iT1&leaf=xW
a);cX
synth_it;for(x93
tJ1::const_iterator
i=nI2.y73
i!=nI2.end();++i){if(i->first!=leaf.GetHash())yA1
const
eF3&occ
nY2
first;size_t
score=occ.GetCSEscore();iT1&candidate
nY2
i63
if(lF2
candidate))yA1
if(leaf
nP2<occ.MinimumDepth())yA1
if(score<2)yA1
if(lP1
i53
leaf)eD==false)continue
l43
if(cW2(i53
leaf,synth,nI2))lD2
eE3
t23
yP1
iI1
iT1&y53,iT1&expr){yY1
y53
lF3
expr))lD2
yY1
iI1
y53
l8
a),expr
y03
true
i1
t23
yP1
GoodMomentForCSE
iR1
y53,iT1&expr){if(y53
nC==cIf)lD2
yY1
y53
lF3
expr))lD2
size_t
i72=0;yY1
iI1
y53
l8
a),expr))++i72
i1
i72!=1;}
}
t5{tK1
size_t
yZ2::SynthCommonSubExpressions(xK1::xI1
const{if(e91
0)return
0;size_t
stacktop_before=synth.GetStackTop();tJ1
nI2;FindTreeCounts(nI2,*this,GetOpcode(),true);for(;;){size_t
yW2=0;
#ifdef DEBUG_SUBSTITUTIONS_CSE
cV2<<"Finding a CSE candidate, root is:"
<<std::y11*this);
#endif
cX
cs_it(nI2.end());for(cX
j=nI2.y73
j!=nI2.end();){cX
i(j++);const
eF3&occ
nY2
first;size_t
score=occ.GetCSEscore();iT1&tree
nY2
i63
#ifdef DEBUG_SUBSTITUTIONS_CSE
cV2<<"Score "
<<score<<":\n"
<<std::flush;DumpTreeWithIndent(tree);
#endif
if(lF2
tree))xZ
if(tree
nP2<occ.MinimumDepth())xZ
if(score<2)xZ
if(lP1*this
cU3)eD==false)xZ
if(cW2(*this
cU3,synth,nI2)){yA1}
if(!GoodMomentForCSE(*this
cU3))xZ
score*=tree
nP2;if(score>yW2){yW2=score;cs_it=i;}
}
if(yW2<=0){
#ifdef DEBUG_SUBSTITUTIONS_CSE
cV2<<"No more CSE candidates.\n"
<<std::flush;
#endif
xY3
iT1&tree=cs_it
eC2.i63
#ifdef DEBUG_SUBSTITUTIONS_CSE
cV2<<iX3"Common Subexpression:"
;DumpTree
xF(tree
lE2
std::endl;
#endif
#if 0
int
n11=occ.NeedsSinCos();int
i9=occ.NeedsSinhCosh()tX
i82,i92,yX2,yY2;if(n11){i82
eX2
i82
iH
cSin);i82
iW1;i92
eX2
i92
iH
cCos);i92
iW1;if(lF2
i82)||lF2
i92))i73==2){nI2.xG2
yA1}
n11=0;}
}
if(i9){yX2
eX2
yX2
iH
cSinh);yX2
iW1;yY2
eX2
yY2
iH
cCosh);yY2
iW1;if(lF2
yX2)||lF2
yY2)){if(i9==2){nI2.xG2
yA1}
i9=0;}
}
#endif
tree.SynthesizeByteCode(synth,false);nI2.xG2
#ifdef DEBUG_SUBSTITUTIONS_CSE
synth.x73
Dump<0>(lE2"Done with Common Subexpression:"
;DumpTree
xF(tree
lE2
std::endl;
#endif
#if 0
if(n11)i73==2||i9){synth.eL1}
yB1
cSinCos,1,2)yC1
i82,1)yC1
i92,0);}
if(i9)i73)synth.eL1
if(i9==2){synth.eL1}
yB1
cSinhCosh,1,2)yC1
yX2,1)yC1
yY2,0);}
#endif
eE3
xV3
stacktop_before;}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
tK1
lU1
xF::iA2{using
t5;CopyOnWrite()tX
tree;tree.GenerateFrom(*mData);FPoptimizer_Optimize::ApplyGrammars(tree);eO3
x83>cT3;std
xM3
xF
immed;size_t
stacktop_max=0;tree.SynthesizeByteCode(cT3,immed,stacktop_max);if(mData->mStackSize!=stacktop_max){mData->mStackSize=x83(stacktop_max);
#if !defined(FP_USE_THREAD_SAFE_EVAL) && \
    !defined(FP_USE_THREAD_SAFE_EVAL_WITH_ALLOCA)
mData->mStack
xF3
stacktop_max);
#endif
}
mData->mByteCode.swap(cT3);mData->mImmed.swap(immed);}
#define FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE(type) tP1>lU1<type>::iA2{}
#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
i93(MpfrFloat)
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
i93(GmpInt)
#endif
#ifdef FP_SUPPORT_COMPLEX_DOUBLE_TYPE
i93(std::complex<double>)
#endif
#ifdef FP_SUPPORT_COMPLEX_FLOAT_TYPE
i93(std::complex<float>)
#endif
#ifdef FP_SUPPORT_COMPLEX_LONG_DOUBLE_TYPE
i93(std::complex<long
double>)
#endif
#define FUNCTIONPARSER_INSTANTIATE_OPTIMIZE(type) x73 lU1<type>::iA2;
#ifndef FP_DISABLE_DOUBLE_TYPE
iA3(double)
#endif
#ifdef FP_SUPPORT_FLOAT_TYPE
iA3(float)
#endif
#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
iA3(long
double)
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
iA3(long)
#endif
#endif // FP_SUPPORT_OPTIMIZER

#endif
