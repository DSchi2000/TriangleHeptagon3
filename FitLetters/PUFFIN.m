(* ::Package:: *)

(* ::Title:: *)
(*PUFFIN*)


(* ::Section:: *)
(*Launch*)


(*Start*)
BeginPackage["PUFFIN`"];
Print["+++++ PUFFIN` +++++\n"];
Print["Practical Utilities For Feynman INtegrals"];


(* ::Section:: *)
(*Documentation*)


SetP::usage="SetP[p$_]. Sets P to p$ (deefault: 2^31-1).";
PrintP::usage="PrintPath[]. Prints the current value for p.";
PLog::usage="PLog[x_]. Computes the p-adic logarithm of x";
RationalReconstruct::usage="RationalReconstruct[a_]. Perform the rational reconstructoin of the p-adic number a modulo P.";
RankFactorisation::usage="";
PLogRankFactorisation::usage="";
PLogRankFactorisationQ::usage="";
PLogMatrix::usage="";
IrrationalReconstruct::usage="";
AlgebraicReconstructH::usage="Reconstruct multi variable, homogenous in variables";
AlgebraicReconstruct1V::usage="Reconstruct single variable, not homogenous";


SP::usage="";
rt::usage="";


Begin["`Private`"];


(* ::Section:: *)
(*SP*)


SetAttributes[SP,Orderless]

SP[0,0]:=0;
SP[{a_},{a_}]:=d;
SP[-a_,b_]:=-SP[a,b];
SP[a_,n_Integer b_]:=n SP[a,b]
SP[a_,n_Real b_]:=n SP[a,b];
SP[a_,n_Complex b_]:=n SP[a,b];
SP[a__+b_,c__+d_]:=SP[a,c]+SP[a,d]+SP[b,c]+SP[b,d];
SP[b_,c__+d_]:=SP[b,c]+SP[b,d];
SP[a__+b_,d_]:=SP[a,d]+SP[b,d];

SP/:SP[a_,{b_}]^2:=SP[a,a];
SP/:SP[a_,{b_}] GTr[x___,{b_},y___]:=GTr[x,a,y];
SP/:SP[a_,{b_}]SP[{b_},c_]:=SP[a,c];
SP/:SP[m_,{n_}] HoldPattern[Ep[a___,{n_},z___]]:=Ep[a,m,z];                               
SP/:SP[m_,{n_}] InactiveEp[a___,{n_},z___]:=InactiveEp[a,m,z];


(* ::Section:: *)
(*P-adic functions*)


P$=2^31-1;
SetP[p$_] := Module[{}, P$ = p$];
PrintP[] := Print[P$];


PLog[x_]:=1-PowerMod[x,P$-1,P$^2];
RationalReconstruct[a_]:=With[{v=LatticeReduce[{{a,1},{P$,0}}][[1]]},v[[1]]/v[[2]]];


(*PLogMatrix[var_,sub_]:=Module[{matrix},
(*to do: check if enough sub, both if all the variables are replaced and if length[var]+1 sub*)
matrix={ConstantArray[Mod[PLog[Mod[2,P$^2]]/P$,P$],Length[var]+1]};
Do[
matrix=Append[matrix,(Table[Expand[1/P$ PLog[Mod[var[[i]]/.sub[[j]],P$^2]/.{Rational[a_,b_]:> a*PowerMod[b,-1,P$^2]}],Modulus->P$],{j,Length[var]+1}])]
,{i,Length[var]}];
matrix
]*)

PLogMatrix[var_, sub_] := Module[{matrix, safeSub},
  safeSub = Select[sub,
    Function[pt, AllTrue[var, Mod[# /. pt, P$] =!= 0 &]]
  ];
  If[Length[safeSub] < Length[var] + 1,
    Print["Not enough safe points: ", Length[safeSub], " available, need ", Length[var]+1];
    Return[$Failed]
  ];
  safeSub = safeSub[[1 ;; Length[var]+1]];
  matrix = {ConstantArray[Mod[PLog[Mod[2, P$^2]]/P$, P$], Length[var]+1]};
  Do[
    matrix = Append[matrix,
      Table[
        Expand[1/P$ PLog[Mod[var[[i]] /. safeSub[[j]], P$^2] /.
          {Rational[a_, b_] :> a*PowerMod[b, -1, P$^2]}], Modulus -> P$],
        {j, Length[var]+1}]]
  , {i, Length[var]}];
  matrix
]


PLogRankFactorisation[var_,sub_]:=Module[{matrix,MatrixF,MatrixC={},includedj={},listBasis={},out={}},
(*to do: check if enough sub, both if all the variables are replaced and if length[var]+1 sub*)
matrix=PLogMatrix[var,sub];
MatrixF=DeleteCases[RowReduce[Transpose[matrix],Modulus-> P$],{0..}];

Do[
If[DeleteCases[Table[MatrixF[[j,i]],{j,Length[MatrixF]}],0]=={1}&&!MemberQ[includedj,Position[Table[MatrixF[[j,i]],{j,Length[MatrixF]}],1][[1,1]]],
MatrixC=Append[MatrixC,matrix[[i]]];
includedj=Append[includedj,Position[Table[MatrixF[[j,i]],{j,Length[MatrixF]}],1][[1,1]]]
],{i,Length[MatrixF[[1]]]}];
MatrixC=Transpose[MatrixC];
Do[
listBasis=Append[listBasis,Global`var[Position[Normal[matrix],Normal[Transpose[MatrixC][[i]]]][[1,1]]-1]];
,{i,2,Length[Transpose[MatrixC]]}];

MatrixF=Map[RationalReconstruct[#]&,MatrixF,{2}];

Do[
out=Append[out,Global`var[i]-> Times@@(Join[{2},listBasis]^(Table[MatrixF[[j,i+1]],{j,Length[listBasis]+1}]))];
,{i,Length[var]}];
out=DeleteCases[out,Rule[a__,a__]];

out
]

RankFactorisation[var_]:=Module[{matrix,MatrixF,MatrixC={},includedj={},listBasis={},out={}},
(*to do: check if enough sub, both if all the variables are replaced and if length[var]+1 sub*)
matrix=CoefficientArrays[var,Sort[Variables[var]]][[2]];
MatrixF=DeleteCases[RowReduce[Transpose[matrix]],{0..}];

Do[
If[DeleteCases[Table[MatrixF[[j,i]],{j,Length[MatrixF]}],0]=={1}&&!MemberQ[includedj,Position[Table[MatrixF[[j,i]],{j,Length[MatrixF]}],1][[1,1]]],
MatrixC=Append[MatrixC,matrix[[i]]];
includedj=Append[includedj,Position[Table[MatrixF[[j,i]],{j,Length[MatrixF]}],1][[1,1]]]
],{i,Length[MatrixF[[1]]]}];
MatrixC=Transpose[MatrixC];

Do[
listBasis=Append[listBasis,Global`var[Position[Normal[matrix],Normal[Transpose[MatrixC][[i]]]][[1,1]]]];
,{i,Length[Transpose[MatrixC]]}];

Do[
out=Append[out,Global`var[i]-> Plus@@(Table[MatrixF[[j,i]],{j,Length[MatrixF]}]*listBasis)]
,{i,Length[var]}];
(*out=DeleteCases[out,Rule[a__,a__]];*)

out
]


PLogRankFactorisationQ[var_,sub_]:=Module[{matrix,out=False},
matrix=PLogMatrix[var,sub];
out=MatrixRank[matrix,Modulus-> P$]<Length[var]+1;
out 
]

PLogRankFactorisationQ[matrix_]:=Module[{out=False},
out =MatrixRank[matrix,Modulus-> P$]<Length[matrix];
]


IrrationalReconstruct[irrational_,n_]:=Module[{nRe,nIm,fitRe,fitIm,outRe,outIm,out},
nRe=Re[n];
nIm=Im[n];
fitRe=FindIntegerNullVector[Prepend[N[irrational],nRe]];
fitIm=FindIntegerNullVector[Prepend[N[irrational],nIm]];
If[
And@@Map[Abs[#]<100&,fitRe],
outRe=-(Plus@@(fitRe[[2;;-1]]/fitRe[[1]]*irrational));
,
outRe=nRe;
];
If[
And@@Map[Abs[#]<100&,fitIm],
outIm=-(Plus@@(fitIm[[2;;-1]]/fitIm[[1]]*irrational));
,
outIm=nIm;
];
out=outRe+I outIm;
If[N[out]==N[n],
out,
Print["Error in the Log reconstruction"];
]
]


Options[AlgebraicReconstructH]={Debug->False};
AlgebraicReconstructH[fitvecIN_,pars_,vars_,exp1_,exp2_,roots_:{},rootAnsatz_:1,opts:OptionsPattern[]]:=
Module[{variables,fitvec,flattpos,PowerMatrixUsed,AnsatzNumi,AnsatzDeno,varsA,AnsatzPowersNumi,AnsatzPowersDeno,AnsatzLengthN,AnsatzLengthD,AnsatzLength,Ansatz,Solset,AnsatzTerms,SampleSize,parsCut,dvec,nusedvars,dropl,MatNumiPart,onevc,parsCutOne,AnsatzPowersDenoOne,MatDenoPart,FinalMat,sol3,solrule,fit,reconstructedFit,usedSQRTS,SqrtDegrees,rootAnsDim,rootAnsMono,varssqrt,PhaseSpacePointsRTS,MatEQN,NullSpaces,solutions,DelNumi,DelDeno,SolutionsFuncs,SqrtApearences},
	If[OptionValue[Debug]===True,
		debugFunction=EchoLabel,
		debugFunction=Function[label,Identity]
	];
	
	If[exp1==0 && exp2==0,
	Return[0];
	];
	fitvec=fitvecIN;
	If[ Length[DeleteDuplicates[fitvec]]==1,Return[DeleteDuplicates[fitvec][[1]]]];
	(*variables=invariants[[2;;]];*)
	varsA=vars;
    varssqrt=Join[vars,Table[rt[i],{i,1,Length[roots]}]];
	""//debugFunction["Building Ansatz"];
	(*todo: works without internet? wtf is this danyolo analog vibecoding?*)
	SqrtDegrees=(ResourceFunction["PolynomialTotalDegree"][#,vars]&/@roots)/2;
	rootAnsDim=Replace[MonomialList[rootAnsatz],n_Integer:>0,{1}]/.Times->Plus/.Thread[Rule[Table[rt[i],{i,1,Length[roots]}],SqrtDegrees]];
	rootAnsMono=MonomialList[rootAnsatz];
	AnsatzPowersNumi=Flatten[Table[Join[#,CoefficientRules[rootAnsMono[[sqrtPos]],Table[rt[i],{i,1,Length[roots]}]][[1,1]]]&/@FrobeniusSolve[ConstantArray[1,Length[vars]],exp1-rootAnsDim[[sqrtPos]]],{sqrtPos,1,Length[rootAnsMono]}],1];
	AnsatzPowersDeno=FrobeniusSolve[ConstantArray[1,Length[vars]],exp2];
	
	
	AnsatzLengthN=Length[AnsatzPowersNumi];
	AnsatzLengthD=Length[AnsatzPowersDeno];
	AnsatzLength=AnsatzLengthN+AnsatzLengthD;
	
	SampleSize=Length[fitvec];
	
	If[AnsatzLength>SampleSize, Print[ToString[StringForm["Not enough Samples for Ansatz. `` pts required",AnsatzLength]]];
	  Return[-1]];


	AnsatzNumi=Sum[fac[j]Product[varssqrt[[i]]^AnsatzPowersNumi[[j,i]],{i,1,Length[varssqrt]}],{j,1,Length[AnsatzPowersNumi]}];
	AnsatzDeno=Sum[fac[j+Length[AnsatzPowersNumi]]Product[vars[[i]]^AnsatzPowersDeno[[j,i]],{i,1,Length[vars]}],{j,1,Length[AnsatzPowersDeno]}];
	Ansatz=AnsatzNumi/AnsatzDeno;
	Solset=Sort[Cases[Ansatz,fac[___],Infinity]];
	
	fitvec=Drop[fitvec,SampleSize-(AnsatzLength)];
	parsCut=Drop[pars,SampleSize-(AnsatzLength)];
	""//debugFunction["Ansatz Done"];
	""//debugFunction["Mapping Phasespacepoints"];

PhaseSpacePointsRTS=Map[Table[root[roots[[sqn]]]/.Thread[Rule[vars,#]],{sqn,1,Length[roots]}]&,parsCut]/.root[x_]:>PowerMod[x,1/2,P$];
PhaseSpacePointsRTS=MapThread[Join,{parsCut,PhaseSpacePointsRTS}];
""//debugFunction["Mapping Phasespacepoints Done"];
(*	
	nusedvars=Complement[variables,varsA];
	dropl=Table[Position[variables,nusedvars[[i]]],{i,1,Length[nusedvars]}][[All,All,1]];
		
	parsCut=Transpose[Delete[Transpose[parsCut],dropl]];*)
	""//debugFunction["Constructing Matrix"];
	MatNumiPart=Parallelize@Outer[m31ExpDot,PhaseSpacePointsRTS,AnsatzPowersNumi,1];
	fitvec=Expand[-fitvec,Modulus->P$];
	onevc=Table[{1},AnsatzLengthD];
	parsCutOne=MapThread[Join,{parsCut,Transpose[{fitvec}]}];
	AnsatzPowersDenoOne=MapThread[Join,{AnsatzPowersDeno,onevc}];
	MatDenoPart=Parallelize@Outer[m31ExpDot,parsCutOne,AnsatzPowersDenoOne,1];
	MatEQN=MapThread[Join,{MatNumiPart,MatDenoPart}];
	""//debugFunction["Constructing Matrix Done"];
	""//debugFunction["Computing Nullspace"];
	NullSpaces=NullSpace[MatEQN,Modulus->P$];
	""//debugFunction[ToString[StringForm["Computing Nullspace Done, Solutions: ``",Length[NullSpaces]]]];
	If[ NullSpaces==={},
	""//debugFunction[ToString[StringForm["No Solution Found"]]];
	  Return[-1]
	];
	(*Print[NullSpaces[[1]]];*)
	solutions=Map[Thread[Rule[Solset,#]]&,NullSpaces];
	""//debugFunction["SolRules"];
	DelNumi=Position[ParallelMap[AnsatzNumi/.#&,solutions],0];
	DelDeno=Position[ParallelMap[AnsatzNumi/.#&,solutions],0];
	solutions=Delete[solutions,Join[DelNumi,DelDeno]];
	""//debugFunction[ToString[StringForm["Non trivial solutions: ``",Length[solutions]]]];
	""//debugFunction["Applying Solutions"];
	(*SolutionsFuncs=ParallelMap[Ansatz/.#&,solutions];*)
	SolutionsFuncs=Map[Ansatz/.#&,solutions];
	If[rootAnsatz===1,
		Return[SolutionsFuncs//First]
	];
	""//debugFunction["Choose Functions with min sqrt"];
	SqrtApearences=Sum[Map[countR2[#,rt[sqn]]&,SolutionsFuncs],{sqn,1,Length[roots]}];
	Return[Together[First[SolutionsFuncs[[First[Position[SqrtApearences,Min[SqrtApearences]]]]]]]];
	]


countR2[expr_,cexpr_]:=Module[{},Total[Cases[expr,cexpr->1,{0,Infinity}]]+Total[Cases[expr,Power[cexpr,n_.]:>n,{0,Infinity}]]-Total[Cases[expr,Power[cexpr,n_.]:>1,{0,Infinity}]]]


AlgebraicReconstruct1V[fitvec_,phaseSpacePoints_,var_,powernum_,powerden_,roots_:{},rootAnsatz_:1]:=Module[{sqs,AnsatzNumerator,AnsatzDenominator,Ansatz,Facs,PhaseSpacePointsRTS,SOlSts2,MatEQN,NullSpaces,solutions,SolutionsFuncs,SqrtApearences,sqn,DelNumi,DelDeno,AnsatzPowersNumi,AnsatzPowersDeno,MatNumiPart,MatDenoPart,vars,selectionList,out},
AnsatzNumerator=Expand[Sum[var^n,{n,0,powernum}]*rootAnsatz];
AnsatzDenominator=Expand[Sum[var^n,{n,0,powerden}]];

vars=Join[{var},Table[rt[i],{i,1,Length[roots]}]];
If[ Length[DeleteDuplicates[fitvec]]==1,Return[DeleteDuplicates[fitvec][[1]]]];
AnsatzPowersNumi=CoefficientRules[AnsatzNumerator,vars][[All,1]];
AnsatzPowersDeno=CoefficientRules[AnsatzDenominator,vars][[All,1]];

AnsatzNumerator=Sum[fac[j]Product[vars[[i]]^AnsatzPowersNumi[[j,i]],{i,1,Length[vars]}],{j,1,Length[AnsatzPowersNumi]}];
AnsatzDenominator=Sum[fac[j+Length[AnsatzNumerator]]Product[vars[[i]]^AnsatzPowersDeno[[j,i]],{i,1,Length[vars]}],{j,1,Length[AnsatzPowersDeno]}];
Ansatz=AnsatzNumerator/AnsatzDenominator;
Facs=Sort[Cases[Ansatz,_fac,Infinity]];

PhaseSpacePointsRTS=Map[Table[root[roots[[sqn]]]/.#,{sqn,1,Length[roots]}]&,phaseSpacePoints]/.root[x_]:>PowerMod[x,1/2,P$];
PhaseSpacePointsRTS=MapThread[Join,{phaseSpacePoints[[All,All,2]],PhaseSpacePointsRTS}];
selectionList=Join[First[Position[phaseSpacePoints[[1,All,1]],var]],Table[i,{i,Length[phaseSpacePoints[[1,All,1]]]+1,Length[phaseSpacePoints[[1,All,1]]]+Length[roots]}]];
PhaseSpacePointsRTS=PhaseSpacePointsRTS[[All,selectionList]];

MatNumiPart=Expand[Table[Table[Times@@PowerMod[PhaseSpacePointsRTS[[i]],AnsatzPowersNumi[[q]],P$],{q,1,Length[AnsatzPowersNumi]}],{i,1,Length[PhaseSpacePointsRTS]}],Modulus->P$];
AnsatzPowersNumi;
PhaseSpacePointsRTS=MapThread[Join,{PhaseSpacePointsRTS,Transpose[{Expand[-fitvec,Modulus->P$]}]}];

AnsatzPowersDeno=MapThread[Join,{AnsatzPowersDeno,Table[{1},powerden+1]}];

MatDenoPart=Table[Table[Times@@PowerMod[PhaseSpacePointsRTS[[i]],AnsatzPowersDeno[[j]],P$],{j,1,Length[AnsatzPowersDeno]}],{i,1,Length[PhaseSpacePointsRTS]}];

MatEQN=Expand[MapThread[Join,{MatNumiPart,MatDenoPart}],Modulus->P$];
NullSpaces=NullSpace[MatEQN,Modulus->P$];
solutions=Map[Thread[Rule[Facs,#]]&,NullSpaces];
DelNumi=Position[AnsatzNumerator/.solutions,0];
DelDeno=Position[AnsatzDenominator/.solutions,0];
solutions=Delete[solutions,Join[DelNumi,DelDeno]];
SolutionsFuncs=Ansatz/.solutions;
If[solutions=={},Return[Print["Error: reconstruction not possible."]]];
If[rootAnsatz===1,
Return[Factor[SolutionsFuncs//First,Modulus->P$]]
];
SqrtApearences=Sum[Map[countR2[#,rt[sqn]]&,SolutionsFuncs],{sqn,1,Length[roots]}];
Return[DeleteDuplicates[Factor[Together[First[SolutionsFuncs[[First[Position[SqrtApearences,Min[SqrtApearences]]]]]]],Modulus->P$]]];
]


(* ::Section:: *)
(*Ben: Fast evaluation Tech*)


m31Pow = With[{p = 2^31 - 1}, Compile[{{x, _Integer}, {n, _Integer}}, Module[{xn = 1}, Do[xn = Mod[xn*x, p], {n}]; xn]]];


m31ExpDot = With[{p=2^31-1}, 
    Compile[{{val, _Integer, 1}, {exp, _Integer, 1}}, 
        Module[{res = 1}, Do[res = Mod[res*m31Pow[val[[i]], exp[[i]]], p], {i,Length@exp}]; res], 
        CompilationOptions -> {"InlineExternalDefinitions" -> True}, RuntimeOptions -> "Speed"
    ]
 ];



End[];
EndPackage[];
