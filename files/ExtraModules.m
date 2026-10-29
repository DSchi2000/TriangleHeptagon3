(* ::Package:: *)

GetRandomPT[NMin_,NMax_,invariantsIn_,sqrtSetIn_]:=Module[{tPoint,tPointReplacements,sqrtSetT,roots,NMAXIT,s16Sol},    
    NMAXIT=1000000;
    For[i = 0, i < NMAXIT, i++,
		tPoint=RandomInteger[{NMin,NMax},Length[invariantsIn]];
		tPointReplacements=Thread[Rule[invariantsIn,tPoint]];
		sqrtSetT = Expand[sqrtSetIn /. tPointReplacements, Modulus -> $P];
        roots = JacobiSymbol[#, $P] & /@ sqrtSetT;
        If[Total[roots] == Length[roots],
            Return[tPointReplacements];];
];
 Return[{}];
]
GetUnivariantSlice[a_,b_,invariantsIn_,sqrtSetIn_]:=Module[{tPoint,tPointReplacements,sqrtSetT,roots,NMAXIT,s16Sol,tVal},    
    NMAXIT=1000000;
    For[i = 0, i < NMAXIT, i++,
		tVal=RandomInteger[{0,$P}];
		tPoint=Expand[a+tVal b,Modulus->$P];
		tPointReplacements=Thread[Rule[invariantsIn,tPoint]];
		tPointReplacements=Join[tPointReplacements,{tvar->tVal}];
		sqrtSetT = Expand[sqrtSetIn /. tPointReplacements, Modulus -> $P];
                   roots = JacobiSymbol[#, $P] & /@ sqrtSetT;
        If[Total[roots] == Length[roots],
            Return[tPointReplacements];];
	];
 Return[{}];
]
GetOneVarVariance[a_,pos_]:=Module[{tPoint,tPointReplacements,sqrtSetT,roots,NMAXIT,s16Sol,tVal},    
    NMAXIT=1000000;
    For[i = 0, i < NMAXIT, i++,
		tVal=RandomInteger[{0,$P}];
	  tPoint=a;
		tPoint[[pos]]=tVal;
		tPointReplacements=Thread[Rule[invariants,tPoint]];
		tPointReplacements=tPointReplacements;
		sqrtSetT = Expand[sqrtSet /. tPointReplacements, Modulus -> $P];
                   roots = JacobiSymbol[#, $P] & /@ sqrtSetT;
        If[Total[roots] == Length[roots],
            Return[tPointReplacements];];
	];
 Return[{}];
]


countR2[expr_,cexpr_]:=Module[{},Total[Cases[expr,cexpr->1,{0,Infinity}]]+Total[Cases[expr,Power[cexpr,n_.]:>n,{0,Infinity}]]-Total[Cases[expr,Power[cexpr,n_.]:>1,{0,Infinity}]]];


(*OneVarRecoRTAnsatz2[fitvec_,phaseSpacePoints_,var_,powernum_,powerden_,roots_,rootAnsatz_]:=Module[{sqs,AnsatzNumerator,AnsatzDenominator,Ansatz,Facs,PhaseSpacePointsRTS,SOlSts2,MatEQN,NullSpaces,solutions,SolutionsFuncs,SqrtApearences,sqn,DelNumi,DelDeno,AnsatzPowersNumi,AnsatzPowersDeno,MatNumiPart,MatDenoPart,vars,selectionList},
AnsatzNumerator=Expand[Sum[var^n,{n,0,powernum}]*rootAnsatz];
AnsatzDenominator=Expand[Sum[var^n,{n,0,powernum}]];

vars=Join[{var},Table[rt[i],{i,1,Length[roots]}]];
If[ Length[DeleteDuplicates[fitvec]]==1,Return[DeleteDuplicates[fitvec][[1]]]];
AnsatzPowersNumi=CoefficientRules[AnsatzNumerator,vars][[All,1]];
AnsatzPowersDeno=CoefficientRules[AnsatzDenominator,vars][[All,1]];

AnsatzNumerator=Sum[fac[j]Product[vars[[i]]^AnsatzPowersNumi[[j,i]],{i,1,Length[vars]}],{j,1,Length[AnsatzPowersNumi]}];
AnsatzDenominator=Sum[fac[j+Length[AnsatzNumerator]]Product[vars[[i]]^AnsatzPowersDeno[[j,i]],{i,1,Length[vars]}],{j,1,Length[AnsatzPowersDeno]}];
Ansatz=AnsatzNumerator/AnsatzDenominator;
Facs=Sort[Cases[Ansatz,_fac,Infinity]];

PhaseSpacePointsRTS=Map[Table[root[roots[[sqn]]]/.#,{sqn,1,Length[roots]}]&,phaseSpacePoints]/.root[x_]:>PowerMod[x,1/2,$P];
PhaseSpacePointsRTS=MapThread[Join,{phaseSpacePoints[[All,All,2]],PhaseSpacePointsRTS}];



selectionList=Join[First[Position[phaseSpacePoints[[1,All,1]],var]],Table[i,{i,Length[phaseSpacePoints[[1,All,1]]]+1,Length[phaseSpacePoints[[1,All,1]]]+Length[roots]}]];
PhaseSpacePointsRTS=PhaseSpacePointsRTS[[All,selectionList]];

MatNumiPart=Expand[Table[Table[Times@@PowerMod[PhaseSpacePointsRTS[[i]],AnsatzPowersNumi[[q]],$P],{q,1,Length[AnsatzPowersNumi]}],{i,1,Length[PhaseSpacePointsRTS]}],Modulus->$P];
AnsatzPowersNumi;
PhaseSpacePointsRTS=MapThread[Join,{PhaseSpacePointsRTS,Transpose[{Expand[-fitvec,Modulus->$P]}]}];
AnsatzPowersDeno=MapThread[Join,{AnsatzPowersDeno,Table[{1},powerden+1]}];

MatDenoPart=Table[Table[Times@@PowerMod[PhaseSpacePointsRTS[[i]],AnsatzPowersDeno[[j]],$P],{j,1,Length[AnsatzPowersDeno]}],{i,1,Length[PhaseSpacePointsRTS]}];

MatEQN=Expand[MapThread[Join,{MatNumiPart,MatDenoPart}],Modulus->$P];

NullSpaces=NullSpace[MatEQN,Modulus->$P];
solutions=Map[Thread[Rule[Facs,#]]&,NullSpaces];
DelNumi=Position[AnsatzNumerator/.solutions,0];
DelDeno=Position[AnsatzDenominator/.solutions,0];
solutions=Delete[solutions,Join[DelNumi,DelDeno]];
SolutionsFuncs=Ansatz/.solutions;
If[rootAnsatz===1,
Return[SolutionsFuncs//First]
];
SqrtApearences=Sum[Map[countR2[#,rt[sqn]]&,SolutionsFuncs],{sqn,1,Length[roots]}];
Return[Together[First[SolutionsFuncs[[First[Position[SqrtApearences,Min[SqrtApearences]]]]]]]];
]*)

OneVarRecoRTAnsatz2[fitvec_,phaseSpacePoints_,var_,powernum_,powerden_,roots_,rootAnsatz_]:=Module[{sqs,AnsatzNumerator,AnsatzDenominator,Ansatz,Facs,PhaseSpacePointsRTS,SOlSts2,MatEQN,NullSpaces,solutions,SolutionsFuncs,SqrtApearences,sqn,DelNumi,DelDeno,AnsatzPowersNumi,AnsatzPowersDeno,MatNumiPart,MatDenoPart,vars,selectionList},
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


PhaseSpacePointsRTS=Map[Table[root[roots[[sqn]]]/.#,{sqn,1,Length[roots]}]&,phaseSpacePoints]//.root[x_]:>PowerMod[x,1/2,$P];
PhaseSpacePointsRTS=MapThread[Join,{phaseSpacePoints[[All,All,2]],PhaseSpacePointsRTS}];



selectionList=Join[First[Position[phaseSpacePoints[[1,All,1]],var]],Table[i,{i,Length[phaseSpacePoints[[1,All,1]]]+1,Length[phaseSpacePoints[[1,All,1]]]+Length[roots]}]];
PhaseSpacePointsRTS=PhaseSpacePointsRTS[[All,selectionList]];

MatNumiPart=Expand[Table[Table[Times@@PowerMod[PhaseSpacePointsRTS[[i]],AnsatzPowersNumi[[q]],$P],{q,1,Length[AnsatzPowersNumi]}],{i,1,Length[PhaseSpacePointsRTS]}],Modulus->$P];
AnsatzPowersNumi;
PhaseSpacePointsRTS=MapThread[Join,{PhaseSpacePointsRTS,Transpose[{Expand[-fitvec,Modulus->$P]}]}];
AnsatzPowersDeno=MapThread[Join,{AnsatzPowersDeno,Table[{1},powerden+1]}];

MatDenoPart=Table[Table[Times@@PowerMod[PhaseSpacePointsRTS[[i]],AnsatzPowersDeno[[j]],$P],{j,1,Length[AnsatzPowersDeno]}],{i,1,Length[PhaseSpacePointsRTS]}];

MatEQN=Expand[MapThread[Join,{MatNumiPart,MatDenoPart}],Modulus->$P];


NullSpaces=NullSpace[MatEQN,Modulus->$P];
solutions=Map[Thread[Rule[Facs,#]]&,NullSpaces];
DelNumi=Position[AnsatzNumerator/.solutions,0];
DelDeno=Position[AnsatzDenominator/.solutions,0];
solutions=Delete[solutions,Join[DelNumi,DelDeno]];
SolutionsFuncs=Ansatz/.solutions;
If[rootAnsatz===1,
Return[SolutionsFuncs//First]
];
SqrtApearences=Sum[Map[countR2[#,rt[sqn]]&,SolutionsFuncs],{sqn,1,Length[roots]}];
Return[Together[First[SolutionsFuncs[[First[Position[SqrtApearences,Min[SqrtApearences]]]]]]]];
]


OneVarRecoRTSymmetric[fitvec_,phaseSpacePoints_,var_,powernum_,powerden_,roots_,rootAnsatz_]:=Module[{PhaseSpacePointsRTS2,sqs,AnsatzNumerator,AnsatzDenominator,Ansatz,Facs,PhaseSpacePointsRTS,SOlSts2,MatEQN,NullSpaces,solutions,SolutionsFuncs,SqrtApearences,sqn,DelNumi,DelDeno,AnsatzPowersNumi,AnsatzPowersDeno,MatNumiPart,MatDenoPart,vars,selectionList},
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


PhaseSpacePointsRTS=Map[Table[root[roots[[sqn]]]/.#,{sqn,1,Length[roots]}]&,phaseSpacePoints]//.root[x_]:>PowerMod[x,1/2,$P];
PhaseSpacePointsRTS=MapThread[Join,{phaseSpacePoints[[All,All,2]],PhaseSpacePointsRTS}];

PhaseSpacePointsRTS2=PhaseSpacePointsRTS;
PhaseSpacePointsRTS2[[All,-1]]=Expand[-PhaseSpacePointsRTS[[All,-1]],Modulus->$P];
PhaseSpacePointsRTS=Flatten[Transpose[{PhaseSpacePointsRTS,PhaseSpacePointsRTS2}],1]//Echo;
fitvec//Echo;


selectionList=Join[First[Position[phaseSpacePoints[[1,All,1]],var]],Table[i,{i,Length[phaseSpacePoints[[1,All,1]]]+1,Length[phaseSpacePoints[[1,All,1]]]+Length[roots]}]];
PhaseSpacePointsRTS=PhaseSpacePointsRTS[[All,selectionList]];

MatNumiPart=Expand[Table[Table[Times@@PowerMod[PhaseSpacePointsRTS[[i]],AnsatzPowersNumi[[q]],$P],{q,1,Length[AnsatzPowersNumi]}],{i,1,Length[PhaseSpacePointsRTS]}],Modulus->$P];
AnsatzPowersNumi;
PhaseSpacePointsRTS=MapThread[Join,{PhaseSpacePointsRTS,Transpose[{Expand[-fitvec,Modulus->$P]}]}];
AnsatzPowersDeno=MapThread[Join,{AnsatzPowersDeno,Table[{1},powerden+1]}];

MatDenoPart=Table[Table[Times@@PowerMod[PhaseSpacePointsRTS[[i]],AnsatzPowersDeno[[j]],$P],{j,1,Length[AnsatzPowersDeno]}],{i,1,Length[PhaseSpacePointsRTS]}];

MatEQN=Expand[MapThread[Join,{MatNumiPart,MatDenoPart}],Modulus->$P];


NullSpaces=NullSpace[MatEQN,Modulus->$P];
solutions=Map[Thread[Rule[Facs,#]]&,NullSpaces];
DelNumi=Position[AnsatzNumerator/.solutions,0];
DelDeno=Position[AnsatzDenominator/.solutions,0];
solutions=Delete[solutions,Join[DelNumi,DelDeno]];
SolutionsFuncs=Ansatz/.solutions;
If[rootAnsatz===1,
Return[SolutionsFuncs//First]
];
SqrtApearences=Sum[Map[countR2[#,rt[sqn]]&,SolutionsFuncs],{sqn,1,Length[roots]}];
Return[Together[First[SolutionsFuncs[[First[Position[SqrtApearences,Min[SqrtApearences]]]]]]]];
]


OneVarRecoRTAnsatzND[fitvec_,phaseSpacePoints_,var_,powernum_,powerden_,roots_,rootAnsatzN_,rootAnsatzD_]:=Module[{sqs,AnsatzNumerator,AnsatzDenominator,Ansatz,Facs,PhaseSpacePointsRTS,SOlSts2,MatEQN,NullSpaces,solutions,SolutionsFuncs,SqrtApearences,sqn,DelNumi,DelDeno,AnsatzPowersNumi,AnsatzPowersDeno,MatNumiPart,MatDenoPart,vars,selectionList},
AnsatzNumerator=Expand[Sum[var^n,{n,0,powernum}]*rootAnsatzN];
AnsatzDenominator=Expand[Sum[var^n,{n,0,powernum}]*rootAnsatzD];

vars=Join[{var},Table[rt[i],{i,1,Length[roots]}]];
If[ Length[DeleteDuplicates[fitvec]]==1,Return[DeleteDuplicates[fitvec][[1]]]];
AnsatzPowersNumi=CoefficientRules[AnsatzNumerator,vars][[All,1]];
AnsatzPowersDeno=CoefficientRules[AnsatzDenominator,vars][[All,1]];

AnsatzNumerator=Sum[fac[j]Product[vars[[i]]^AnsatzPowersNumi[[j,i]],{i,1,Length[vars]}],{j,1,Length[AnsatzPowersNumi]}];
AnsatzDenominator=Sum[fac[j+Length[AnsatzNumerator]]Product[vars[[i]]^AnsatzPowersDeno[[j,i]],{i,1,Length[vars]}],{j,1,Length[AnsatzPowersDeno]}];
Ansatz=AnsatzNumerator/AnsatzDenominator;
Facs=Sort[Cases[Ansatz,_fac,Infinity]];
If[ Length[Facs]>Length[fitvec],Print[ToString[StringForm["Not enough Samples for Ansatz. `` pts required",Length[Facs]]],Return[-1]]];

PhaseSpacePointsRTS=Map[Table[root[roots[[sqn]]]/.#,{sqn,1,Length[roots]}]&,phaseSpacePoints]//.root[x_]:>PowerMod[x,1/2,$P];
PhaseSpacePointsRTS=MapThread[Join,{phaseSpacePoints[[All,All,2]],PhaseSpacePointsRTS}];



selectionList=Join[First[Position[phaseSpacePoints[[1,All,1]],var]],Table[i,{i,Length[phaseSpacePoints[[1,All,1]]]+1,Length[phaseSpacePoints[[1,All,1]]]+Length[roots]}]];
PhaseSpacePointsRTS=PhaseSpacePointsRTS[[All,selectionList]];

MatNumiPart=Expand[Table[Table[Times@@PowerMod[PhaseSpacePointsRTS[[i]],AnsatzPowersNumi[[q]],$P],{q,1,Length[AnsatzPowersNumi]}],{i,1,Length[PhaseSpacePointsRTS]}],Modulus->$P];
AnsatzPowersNumi;
PhaseSpacePointsRTS=MapThread[Join,{PhaseSpacePointsRTS,Transpose[{Expand[-fitvec,Modulus->$P]}]}];
AnsatzPowersDeno=MapThread[Join,{AnsatzPowersDeno,Table[{1},Length[AnsatzPowersDeno]]}];

MatDenoPart=Expand[Table[Table[Times@@PowerMod[PhaseSpacePointsRTS[[i]],AnsatzPowersDeno[[j]],$P],{j,1,Length[AnsatzPowersDeno]}],{i,1,Length[PhaseSpacePointsRTS]}],Modulus->$P];

MatEQN=Expand[MapThread[Join,{MatNumiPart,MatDenoPart}],Modulus->$P];

NullSpaces=NullSpace[MatEQN,Modulus->$P];
solutions=Map[Thread[Rule[Facs,#]]&,NullSpaces];
DelNumi=Position[AnsatzNumerator/.solutions,0];
DelDeno=Position[AnsatzDenominator/.solutions,0];
solutions=Delete[solutions,Join[DelNumi,DelDeno]];
SolutionsFuncs=Ansatz/.solutions;
If[rootAnsatz===1,
Return[SolutionsFuncs/.tvar->t]
];
Return[DeleteDuplicates[Together[SolutionsFuncs,Modulus->$P]]/.tvar->t]
(*SqrtApearences=Sum[Map[countR2[#,rt[sqn]]&,SolutionsFuncs],{sqn,1,Length[roots]}];
Return[Together[First[SolutionsFuncs[[First[Position[SqrtApearences,Max[SqrtApearences]]]]]]]];*)
]


RationalReconstruct[a_]:=With[{v=LatticeReduce[{{a,1},{$P,0}}][[1]]},v[[1]]/v[[2]]];
Options[ReconstructFullSqrtAnsatz]={Debug->False};
ReconstructFullSqrtAnsatz[fitvecIN_,pars_,vars_,exp1_,exp2_,roots_,rootAnsatz_,opts:OptionsPattern[]]:=
Module[{variables,fitvec,flattpos,PowerMatrixUsed,AnsatzNumi,AnsatzDeno,varsA,AnsatzPowersNumi,AnsatzPowersDeno,AnsatzLengthN,AnsatzLengthD,AnsatzLength,Ansatz,Solset,AnsatzTerms,SampleSize,parsCut,dvec,nusedvars,dropl,MatNumiPart,onevc,parsCutOne,AnsatzPowersDenoOne,MatDenoPart,FinalMat,sol3,solrule,fit,reconstructedFit,usedSQRTS,SqrtDegrees,rootAnsDim,rootAnsMono,varssqrt,PhaseSpacePointsRTS,MatEQN,NullSpaces,solutions,DelNumi,DelDeno,SolutionsFuncs,SqrtApearences},
	If[OptionValue[Debug]===True,
		debugFunction=EchoLabel,
		debugFunction=Function[label,Identity]
	];
	If[ Length[DeleteDuplicates[fitvecIN]]==1,Return[DeleteDuplicates[fitvecIN][[1]]]];
	If[exp1==0 && exp2==0,
	Return[0];
	];
	fitvec=fitvecIN;
	
	(*variables=invariants[[2;;]];*)
	varsA=vars;
    varssqrt=Join[vars,Table[rt[i],{i,1,Length[roots]}]];
	""//debugFunction["Building Ansatz"];
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

PhaseSpacePointsRTS=Map[Table[root[roots[[sqn]]]/.Thread[Rule[vars,#]],{sqn,1,Length[roots]}]&,parsCut]/.root[x_]:>PowerMod[x,1/2,$P];
PhaseSpacePointsRTS=MapThread[Join,{parsCut,PhaseSpacePointsRTS}];
""//debugFunction["Mapping Phasespacepoints Done"];
(*	
	nusedvars=Complement[variables,varsA];
	dropl=Table[Position[variables,nusedvars[[i]]],{i,1,Length[nusedvars]}][[All,All,1]];
		
	parsCut=Transpose[Delete[Transpose[parsCut],dropl]];*)
	""//debugFunction["Constructing Matrix"];
	MatNumiPart=Parallelize@Outer[m31ExpDot,PhaseSpacePointsRTS,AnsatzPowersNumi,1];
	fitvec=Expand[-fitvec,Modulus->$P];
	onevc=Table[{1},AnsatzLengthD];
	parsCutOne=MapThread[Join,{parsCut,Transpose[{fitvec}]}];
	AnsatzPowersDenoOne=MapThread[Join,{AnsatzPowersDeno,onevc}];
	MatDenoPart=Parallelize@Outer[m31ExpDot,parsCutOne,AnsatzPowersDenoOne,1];
	MatEQN=MapThread[Join,{MatNumiPart,MatDenoPart}];
	""//debugFunction["Constructing Matrix Done"];
	""//debugFunction["Computing Nullspace"];
	NullSpaces=NullSpace[MatEQN,Modulus->$P];
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
	SolutionsFuncs=ParallelMap[Ansatz/.#&,solutions];
	If[rootAnsatz===1,
		Return[SolutionsFuncs//First]
	];
	""//debugFunction["Choose Functions with min sqrt"];
	SqrtApearences=Sum[Map[countR2[#,rt[sqn]]&,SolutionsFuncs],{sqn,1,Length[roots]}];
	Return[Together[First[SolutionsFuncs[[First[Position[SqrtApearences,Min[SqrtApearences]]]]]]]];
	]


(* ::Section:: *)
(*ReconstructSqrtAnsatzND*)


Options[ReconstructAlg]={Debug->False};
ReconstructAlg[fitvecIN_,pars_,vars_,exp1_,exp2_,roots_,rootAnsatzN_,rootAnsatzD_,opts:OptionsPattern[]]:=
Module[{stream,UsedTime,NumiN,DenoN,variables,NROOTS,fitvec,flattpos,NSqrtDegrees,rtsetNNR,PowerMatrixUsed,AnsatzNumi,AnsatzDeno,varsA,AnsatzPowersNumi,AnsatzPowersDeno,AnsatzLengthN,AnsatzLengthD,AnsatzLength,Ansatz,Solset,AnsatzTerms,SampleSize,parsCut,dvec,nusedvars,dropl,MatNumiPart,onevc,parsCutOne,AnsatzPowersDenoOne,MatDenoPart,FinalMat,sol3,solrule,fit,reconstructedFit,usedSQRTS,SqrtDegrees,rootAnsDim,rootAnsMono,varssqrt,PhaseSpacePointsRTS,MatEQN,NullSpaces,solutions,DelNumi,DelDeno,SolutionsFuncs,SqrtApearences},
	If[OptionValue[Debug]===True,
		debugFunction=EchoLabel,
		debugFunction=Function[label,Identity]
	];
	If[ Length[DeleteDuplicates[fitvecIN]]==1,Return[{{{DeleteDuplicates[fitvecIN][[1]],1}},fac[1]/fac[2]}]];
	If[exp1==0 && exp2==0,
	Return[0];
	];
	fitvec=fitvecIN;
	
	(*variables=invariants[[2;;]];*)
	varsA=vars;
    varssqrt=Join[vars,Table[rt[i],{i,1,Length[roots]}]];
	""//debugFunction["Building Ansatz"];
NROOTS=DeleteDuplicates[Cases[roots,_root,Infinity]]/.root[x_]:>x;
NSqrtDegrees=(ResourceFunction["PolynomialTotalDegree"][#,vars]&/@NROOTS)/2;
rtsetNNR=roots/.Table[Rule[root[NROOTS[[i]]],Power[vars[[1]],NSqrtDegrees[[i]]]],{i,1,Length[NSqrtDegrees]}];
SqrtDegrees=(ResourceFunction["PolynomialTotalDegree"][#,vars]&/@rtsetNNR)/2//Echo;
rootAnsDim=Replace[MonomialList[rootAnsatzN],n_Integer:>0,{1}]/.Times->Plus/.Thread[Rule[Table[rt[i],{i,1,Length[roots]}],SqrtDegrees]];
rootAnsMono=MonomialList[rootAnsatzN];
AnsatzPowersNumi=Flatten[Table[Join[#,CoefficientRules[rootAnsMono[[sqrtPos]],Table[rt[i],{i,1,Length[roots]}]][[1,1]]]&/@FrobeniusSolve[ConstantArray[1,Length[vars]],exp1-rootAnsDim[[sqrtPos]]],{sqrtPos,1,Length[rootAnsMono]}],1];
	
	
rootAnsDim=Replace[MonomialList[rootAnsatzD],n_Integer:>0,{1}]/.Times->Plus/.Thread[Rule[Table[rt[i],{i,1,Length[roots]}],SqrtDegrees]];
rootAnsMono=MonomialList[rootAnsatzD];
AnsatzPowersDeno=Flatten[Table[Join[#,CoefficientRules[rootAnsMono[[sqrtPos]],Table[rt[i],{i,1,Length[roots]}]][[1,1]]]&/@FrobeniusSolve[ConstantArray[1,Length[vars]],exp2-rootAnsDim[[sqrtPos]]],{sqrtPos,1,Length[rootAnsMono]}],1];
	
	
	
	AnsatzLengthN=Length[AnsatzPowersNumi];
	AnsatzLengthD=Length[AnsatzPowersDeno];
	AnsatzLength=AnsatzLengthN+AnsatzLengthD;
	""//debugFunction[ToString[StringForm["Number of Samples Needed: ``",AnsatzLength]]];
	

	SampleSize=Length[fitvec];
	
	If[AnsatzLength>SampleSize, Print[ToString[StringForm["Not enough Samples for Ansatz. `` pts required",AnsatzLength]]];
	  Return[-1]];


	AnsatzNumi=Sum[fac[j]Product[varssqrt[[i]]^AnsatzPowersNumi[[j,i]],{i,1,Length[varssqrt]}],{j,1,Length[AnsatzPowersNumi]}];
	AnsatzDeno=Sum[fac[j+Length[AnsatzPowersNumi]]Product[varssqrt[[i]]^AnsatzPowersDeno[[j,i]],{i,1,Length[varssqrt]}],{j,1,Length[AnsatzPowersDeno]}];
(*	AnsatzDeno=Sum[fac[j+Length[AnsatzPowersNumi]]Product[vars[[i]]^AnsatzPowersDeno[[j,i]],{i,1,Length[vars]}],{j,1,Length[AnsatzPowersDeno]}]//Echo;*)

	Ansatz=AnsatzNumi/AnsatzDeno;
	
	Solset=Sort[Cases[Ansatz,fac[___],Infinity]];
	
	fitvec=Drop[fitvec,SampleSize-(AnsatzLength)];
	parsCut=Drop[pars,SampleSize-(AnsatzLength)];
	""//debugFunction["Ansatz Done"];
	""//debugFunction["Mapping Phasespacepoints"];

PhaseSpacePointsRTS=Map[Table[root[roots[[sqn]]]/.Thread[Rule[vars,#]],{sqn,1,Length[roots]}]&,parsCut]//.root[x_]:>PowerMod[x,1/2,$P];
PhaseSpacePointsRTS=MapThread[Join,{parsCut,PhaseSpacePointsRTS}];
""//debugFunction["Mapping Phasespacepoints Done"];


(*	
	nusedvars=Complement[variables,varsA];
	dropl=Table[Position[variables,nusedvars[[i]]],{i,1,Length[nusedvars]}][[All,All,1]];
		
	parsCut=Transpose[Delete[Transpose[parsCut],dropl]];*)
	""//debugFunction["Constructing Matrix"];
	MatNumiPart=Parallelize@Outer[m31ExpDot,PhaseSpacePointsRTS,AnsatzPowersNumi,1];
	fitvec=Expand[-fitvec,Modulus->$P];
	onevc=Table[{1},AnsatzLengthD];
	parsCutOne=MapThread[Join,{PhaseSpacePointsRTS,Transpose[{fitvec}]}];
	AnsatzPowersDenoOne=MapThread[Join,{AnsatzPowersDeno,onevc}];
	MatDenoPart=Parallelize@Outer[m31ExpDot,parsCutOne,AnsatzPowersDenoOne,1];
	MatEQN=MapThread[Join,{MatNumiPart,MatDenoPart}];
	""//debugFunction["Constructing Matrix Done"];
	""//debugFunction["Computing Nullspace"];
(*	Print[CurrentDir];
	Print[MatEQN//Dimensions];*)
	(*Export[CurrentDir<>"EQNMAT.txt",MatEQN]//EchoTiming;*)
(*	
	stream = OpenWrite["matrix.bin", BinaryFormat -> True];
	BinaryWrite[stream, ToCharacterCode["GFPM"]]//EchoTiming;
	BinaryWrite[stream, Dimensions[MatEQN], "UnsignedInteger32"]//EchoTiming;
	BinaryWrite[stream, Flatten[MatEQN], "UnsignedInteger32"]//EchoTiming;
	Close[stream];
	
	Return[{{{1}},Ansatz}];*)
	
	UsedTime=AbsoluteTiming[NullSpaces=NullSpace[MatEQN,Modulus->$P]];
	""//debugFunction[ToString[StringForm["Computing Nullspace Done in `` s, Solutions: ``",UsedTime[[1]],Length[NullSpaces]]]];
	If[ NullSpaces==={},
	""//debugFunction[ToString[StringForm["No Solution Found"]]];
	  Return[-1]
	];
	NumiN=Expand[{MatNumiPart[[1]]} . Transpose[NullSpaces[[All,;;AnsatzLengthN]]],Modulus->$P];
    DenoN=Expand[({MatDenoPart[[1,;;-2]]}) . (Transpose[NullSpaces[[All,AnsatzLengthN+1;;-2]]]),Modulus->$P];
    If[exp2===0,DenoN={1}];
    If[exp1===0,NumiN={1}];


    DelNumi=Position[NumiN//First,0];
	DelDeno=Position[DenoN//First,0];
	
	NullSpaces=Delete[NullSpaces,Join[DelNumi,DelDeno]];
	""//debugFunction[ToString[StringForm["Non trivial solutions: ``",Length[NullSpaces]]]];
	Return[{NullSpaces,Ansatz}];
	
	(*Print[NullSpaces[[1]]];*)
(*	solutions=Map[Thread[Rule[Solset,#]]&,NullSpaces];*)
(*	""//debugFunction["SolRules"];
	DelNumi=Position[ParallelMap[AnsatzNumi/.#&,solutions],0];
	DelDeno=Position[ParallelMap[AnsatzNumi/.#&,solutions],0];
	solutions=Delete[solutions,Join[DelNumi,DelDeno]];
	""//debugFunction[ToString[StringForm["Non trivial solutions: ``",Length[solutions]]]];
	""//debugFunction["Applying Solutions"];
	SolutionsFuncs=ParallelMap[Ansatz/.#&,solutions];
	If[rootAnsatz===1,
		Return[SolutionsFuncs//First]
	];
	""//debugFunction["Choose Functions with min sqrt"];
	SqrtApearences=Sum[Map[countR2[#,rt[sqn]]&,SolutionsFuncs],{sqn,1,Length[roots]}];
	Return[Together[First[SolutionsFuncs[[First[Position[SqrtApearences,Min[SqrtApearences]]]]]]]];*)
	]


ReconstructFullSqrtAnsatzAllSolutions[fitvecIN_,pars_,vars_,exp1_,exp2_,roots_,rootAnsatz_]:=
Module[{variables,fitvec,flattpos,PowerMatrixUsed,AnsatzNumi,AnsatzDeno,varsA,AnsatzPowersNumi,AnsatzPowersDeno,NumiN,DenoN,AnsatzDenoN,AnsatzLengthN,AnsatzLengthD,AnsatzLength,Ansatz,AnsatzNumiN,Solset,AnsatzTerms,SampleSize,parsCut,dvec,nusedvars,dropl,MatNumiPart,onevc,parsCutOne,AnsatzPowersDenoOne,MatDenoPart,FinalMat,sol3,solrule,fit,reconstructedFit,usedSQRTS,SqrtDegrees,rootAnsDim,rootAnsMono,varssqrt,PhaseSpacePointsRTS,MatEQN,NullSpaces,solutions,DelNumi,DelDeno,SolutionsFuncs,SqrtApearences},
	If[exp1==0 && exp2==0,
	Return[0];
	];
	fitvec=fitvecIN;
	If[ Length[DeleteDuplicates[fitvec]]==1,Return[DeleteDuplicates[fitvec][[1]]]];
	(*variables=invariants[[2;;]];*)
	varsA=vars;
    varssqrt=Join[vars,Table[rt[i],{i,1,Length[roots]}]];
	Print["Building Ansatz"];
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
	Print["Ansatz Done"];
	Print["Mapping Phasespacepoints"];

PhaseSpacePointsRTS=Map[Table[root[roots[[sqn]]]/.Thread[Rule[vars,#]],{sqn,1,Length[roots]}]&,parsCut]/.root[x_]:>PowerMod[x,1/2,$P];
PhaseSpacePointsRTS=MapThread[Join,{parsCut,PhaseSpacePointsRTS}];
Print["Mapping Phasespacepoints Done"];
(*	
	nusedvars=Complement[variables,varsA];
	dropl=Table[Position[variables,nusedvars[[i]]],{i,1,Length[nusedvars]}][[All,All,1]];
		
	parsCut=Transpose[Delete[Transpose[parsCut],dropl]];*)
	Print["Constructing Matrix"];
	MatNumiPart=Parallelize@Outer[m31ExpDot,PhaseSpacePointsRTS,AnsatzPowersNumi,1];
	fitvec=Expand[-fitvec,Modulus->$P];
	onevc=Table[{1},AnsatzLengthD];
	parsCutOne=MapThread[Join,{parsCut,Transpose[{fitvec}]}];
	AnsatzPowersDenoOne=MapThread[Join,{AnsatzPowersDeno,onevc}];
	MatDenoPart=Parallelize@Outer[m31ExpDot,parsCutOne,AnsatzPowersDenoOne,1];
	MatEQN=MapThread[Join,{MatNumiPart,MatDenoPart}];
	Print["Constructing Matrix Done"];
	Print["Computing Nullspace"];
	NullSpaces=NullSpace[MatEQN,Modulus->$P];
	Print[ToString[StringForm["Computing Nullspace Done, Solutions: ``",Length[NullSpaces]]]];
	If[ NullSpaces==={},
	Print[ToString[StringForm["No Solution Found"]]];
	  Return[-1]
	];
    NumiN=Expand[{MatNumiPart[[1]]} . Transpose[NullSpaces[[All,;;AnsatzLengthN]]],Modulus->$P];
    DenoN=Expand[{MatDenoPart[[1,;;-2]]} . Transpose[NullSpaces[[All,AnsatzLengthN+1;;-2]]],Modulus->$P];
(*    Dimensions[MatDenoPart[[1,;;-2]]]//Echo;
    Dimensions[Transpose[NullSpaces[[All,AnsatzLengthN+1;;-2]]]]//Echo;*)
    
    DelNumi=Position[NumiN//First,0];
	DelDeno=Position[DenoN//First,0];
	
	NullSpaces=Delete[NullSpaces,Join[DelNumi,DelDeno]];
	Print[ToString[StringForm["Non trivial solutions: ``",Length[NullSpaces]]]];

	Return[{NullSpaces,Ansatz}];
	]


m31ExpDot = With[{p=2^31-1}, 
    Compile[{{val, _Integer, 1}, {exp, _Integer, 1}}, 
        Module[{res = 1}, Do[res = Mod[res*m31Pow[val[[i]], exp[[i]]], p], {i,Length@exp}]; res], 
        CompilationOptions -> {"InlineExternalDefinitions" -> True}, RuntimeOptions -> "Speed"
    ]
 ];
m31Pow = With[{p = 2^31 - 1}, Compile[{{x, _Integer}, {n, _Integer}}, Module[{xn = 1}, Do[xn = Mod[xn*x, p], {n}]; xn]]];

