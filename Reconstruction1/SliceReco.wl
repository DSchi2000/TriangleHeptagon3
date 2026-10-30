(* ::Package:: *)

(* ::Section:: *)
(*Load Files*)


If[$FrontEnd===Null,$DIR=Directory[],$DIR=NotebookDirectory[]];
parts = FileNameSplit[$DIR];
pos = First@FirstPosition[parts, "TriangleHeptagon3"];
HomeDir = FileNameJoin[parts[[;;pos ]]];
SectorDir = FileNameJoin[parts[[;;pos+1]]];
CurrentDir=$DIR;
kirapath=SectorDir;
$P=2^31-1;
PureIntegrals=ReadList[FileNameJoin[{SectorDir,"SubBasis.txt"}]];
(*PureIntegrals[[;;52]]=ConstantArray[0,52];*)
SetDirectory[CurrentDir];
NotbookDir=NotebookDirectory[];
Get[FileNameJoin[{HomeDir,"files/ExtraFunctions.m"}]];
SQRTS=Cases[PureIntegrals,_root,Infinity]//DeleteDuplicates;
(*SQRTPOWS=Map[Exponent[#,SQRTS]&,PureIntegrals];
SQRTVEC=Map[Times@@Power[SQRTS,#]&,SQRTPOWS];*)
OutPutPath="Samples/Univariate Slice";
phasespacepoints=ReadList[FileNameJoin[{CurrentDir,ToString[StringForm["``/phasespacepoints.txt",OutPutPath]]}]];
tvarlist=Flatten[ReadList[FileNameJoin[{CurrentDir,ToString[StringForm["``/tvars.txt",OutPutPath]]}]][[;;;; 1]]];
phasespacepoints=phasespacepoints[[All,2;;]][[1;; ;; 1]];
diffInvariants={mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56};
PhaseSpacePoints=Map[Thread[Rule[diffInvariants,#]]&,phasespacepoints];
PhaseSpacePoints=Table[Join[PhaseSpacePoints[[i]],{tvarlist[[i]]}],{i,1,Length[PhaseSpacePoints]}];
(*PhaseSpacePoints=Map[Thread[Rule[diffInvariants,#]]&,phasespacepoints];*)
rootSet=(Cases[PureIntegrals[[-5;;]],_root,Infinity]//DeleteDuplicates)/.root[x_]:>x;
GetLinPos[x1_,x2_,NMI_]:=First[Flatten[Position[Flatten[Table[{i,j},{i,1,NMI},{j,1,NMI}],1],{x1,x2}]]];
Get[ToString[StringJoin[HomeDir,"/files/ExtraModules.m"]]];
Get[ToString[StringJoin[HomeDir,"/files/ExtraModules.m"]]];
usedvars={mt2,s12,s123,s23,s234,s34,s345,s45,s56};
For[i=1,i<=Length[usedvars],i++,
LinearTerm[usedvars[[i]]]=ReadList[FileNameJoin[{Directory[],ToString[StringForm["``/LinTerm``_0001.txt",OutPutPath,usedvars[[i]]]]}]];
];









GenerateNumerOfSamples[num_,den_,vars_,roots_,rootAnsatzN_,rootAnsatzD_]:=Module[{NROOTS,NSqrtDegrees,rtsetNNR,SqrtDegrees,rootAnsDim,rootAnsMono,AnsatzPowersNumi,AnsatzPowersDeno,AnsatzLengthN,AnsatzLengthD,AnsatzLength},
NROOTS=DeleteDuplicates[Cases[roots,_root,Infinity]]/.root[x_]:>x;
NSqrtDegrees=(ResourceFunction["PolynomialTotalDegree"][#,diffInvariants]&/@NROOTS)/2;
rtsetNNR=roots/.Table[Rule[root[NROOTS[[i]]],Power[diffInvariants[[1]],NSqrtDegrees[[i]]]],{i,1,Length[NSqrtDegrees]}];
SqrtDegrees=(ResourceFunction["PolynomialTotalDegree"][#,diffInvariants]&/@rtsetNNR)/2;
rootAnsDim=Replace[MonomialList[rootAnsatzN],n_Integer:>0,{1}]/.Times->Plus/.Thread[Rule[Table[rt[i],{i,1,Length[roots]}],SqrtDegrees]];
rootAnsMono=MonomialList[rootAnsatzN];
AnsatzPowersNumi=Flatten[Table[Join[#,CoefficientRules[rootAnsMono[[sqrtPos]],Table[rt[i],{i,1,Length[roots]}]][[1,1]]]&/@FrobeniusSolve[ConstantArray[1,vars],num-rootAnsDim[[sqrtPos]]],{sqrtPos,1,Length[rootAnsMono]}],1];
(*AnsatzPowersNumi//Echo;*)	
	
rootAnsDim=Replace[MonomialList[rootAnsatzD],n_Integer:>0,{1}]/.Times->Plus/.Thread[Rule[Table[rt[i],{i,1,Length[roots]}],SqrtDegrees]];
rootAnsMono=MonomialList[rootAnsatzD];
AnsatzPowersDeno=Flatten[Table[Join[#,CoefficientRules[rootAnsMono[[sqrtPos]],Table[rt[i],{i,1,Length[roots]}]][[1,1]]]&/@FrobeniusSolve[ConstantArray[1,vars],den-rootAnsDim[[sqrtPos]]],{sqrtPos,1,Length[rootAnsMono]}],1];
	
	
	AnsatzLengthN=Length[AnsatzPowersNumi];
	AnsatzLengthD=Length[AnsatzPowersDeno];
	AnsatzLength=AnsatzLengthN+AnsatzLengthD;
Return[AnsatzLength]];


(* ::Section:: *)
(*Find Even Letter Pos*)



rootVector=Map[Times@@DeleteDuplicates[Cases[#,_root,Infinity]]&,PureIntegrals]/.{}->1/.root[x_]:>x;
DeFoldMat=Normal[Total[Table[ArrayReshape[LinearTerm[var][[1,1]],{Length[PureIntegrals],Length[PureIntegrals]}],{var,invariants2}]]]/.n_Integer/;(n>0):>1;


rootVectorN=rootVector/.PhaseSpacePoints[[1]];
RootRatioMat=Expand[Table[Together[rootVectorN[[i1]]/rootVectorN[[i2]]],{i1,1,Length[PureIntegrals]},{i2,1,Length[PureIntegrals]}],Modulus->$P]*DeFoldMat;
EvenPos=Position[RootRatioMat,1];
LinearizedEvenPos=Map[GetLinPos[#[[1]],#[[2]],95]&,EvenPos]


(* ::Section:: *)
(*Variable Shift Tech*)


ChangeAuxVarInDe[varOld_,varNew_,DIsij_,tpoint_]:=Module[{newaux,oldaux,newcoefs,ChangeVarpos,DIsijSUB,newcoefsSUB,posA,posB,EqnConst},
EqnConst=Sum[(Expand[D[G6,var]/.tpoint,Modulus->$P]) ds[var],{var,invariants}];

newaux=First[Solve[EqnConst==0,ds[varNew],Modulus->$P]];
oldaux=First[Solve[EqnConst==0,ds[varOld],Modulus->$P]];
newcoefs=Coefficient[newaux[[1,2]],Table[ds[var],{var,invariants2}]];
ChangeVarpos=Sort[Map[Flatten[Position[invariants2,#]]&,{varOld,varNew}]];
DIsijSUB=Delete[DIsij,ChangeVarpos];
newcoefsSUB=Delete[newcoefs,ChangeVarpos];
posA=First[Flatten[Position[invariants2,varOld]]];
posB=First[Flatten[Position[invariants2,varNew]]];

Return[Expand[(DIsijSUB+newcoefsSUB DIsij[[posB]]),Modulus->$P] . Table[ds[var],{var,Delete[invariants2,ChangeVarpos]}]+ Expand[DIsij[[posB]]newcoefs[[posA]],Modulus->$P] ds[varOld]];

]
pos=GetLinPos[34,34,95];
DIsij=Map[LinearTerm[#][[1,1,pos]]&,usedvars]
DIsij=Insert[DIsij,0,4]
tpoint=PhaseSpacePoints[[1,;;-2]]
ChangeAuxVarInDe[s16,s12,DIsij,tpoint]


(* ::Section:: *)
(*Reconstruction*)


(* ::Subsection:: *)
(*Prepare*)


avec=ReadList[CurrentDir<>"Samples/AVECS.txt"][[1]];
bvec=ReadList[CurrentDir<>"Samples/AVECS.txt"][[2]];
LinearTerm[s16]=SparseArray[ConstantArray[0,{Length[PhaseSpacePoints],1,Length[PureIntegrals]^2}]];


(* ::Subsection:: *)
(*Reco*)


Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
(*pos=GetLinPos[34,34,95];*)
pos=LinearizedEvenPos[[290]]
FitVals=Join[Table[LinearTerm[var][[All,1,pos]],{var,invariants2}]];
Integrands=Table[Integrand=Expand[UniVarReconstructFUNCS[FitVals[[varpos]],PhaseSpacePoints,tvar,20,20,{},1,1],Modulus->$P];Factor[Integrand[[1]],Modulus->$P],{varpos,1,10}]
NewIntDt=Factor[Expand[Integrands . bvec[[All,2]],Modulus->$P],Modulus->$P]


(* ::Section:: *)
(*Integrate*)


Get["/home/dschiebe/Documents/packages/trager/Trager.wl"]
TragerLogSubtractionNumeric[NewIntDt,t,{},{},Modulus->$P]


(* ::Section:: *)
(*SATISFACTORY YEAAAAH*)


Monitor[Do[
pos=LinearizedEvenPos[[position]];
FitVals=Join[Table[LinearTerm[var][[All,1,pos]],{var,invariants2}]];
Integrands=Table[Integrand=Expand[UniVarReconstructFUNCS[FitVals[[varpos]],PhaseSpacePoints,tvar,20,20,{},1,1],Modulus->$P];Factor[Integrand[[1]],Modulus->$P],{varpos,1,10}];
NewIntDt=Factor[Expand[Integrands . bvec[[All,2]],Modulus->$P],Modulus->$P];
dtIntegrands[position]=NewIntDt;
DLogs[position]=TragerLogSubtractionNumeric[NewIntDt,t,{},{},Modulus->$P];
,{position,1,Length[LinearizedEvenPos]}],{ProgressIndicator[position,{1,Length[LinearizedEvenPos]}],position}]


?ProgressIndicator


dlogList=Table[DLogs[position],{position,1,Length[LinearizedEvenPos]}]
IntegrandsSlice=Table[dtIntegrands[position],{position,1,Length[LinearizedEvenPos]}]


dlogList//DeleteDuplicates//Length



Export[CurrentDir<>"EvenIntegrands.txt",IntegrandsSlice]


SlicedLetters=DeleteDuplicates[Cases[dlogList,_Log,Infinity]]/.Log[x_]:>x
Export[CurrentDir<>"SlicedEvenLetters.txt",SlicedLetters]


EvenLetterCandidates=DeleteCases[ReadList[HomeDir<>"/ClaudeLetters/EvenLetterCandidates.txt"],Null];
EvenLetterCandidatesFACS=DeleteElements[DeleteDuplicates[Flatten[Map[FactorList[#][[All,1]]&,EvenLetterCandidates]]],{1,-1}];
Position[Factor[EvenLetterCandidatesFACS/.SliceTpoint,Modulus->$P],0]
EvenLetterCandidatesFACS=Delete[EvenLetterCandidatesFACS,79]
ELSclice=Factor[EvenLetterCandidatesFACS/.SliceTpoint,Modulus->$P];


(*Func=SlicedLetters[[15]]
TestLetter=ELSclice[[1]];*)
pos=1;
CheckLetterFactor[FuncIn_,LetterPos_]:=Module[{Num,Den,SymFactor,FoundFactor,TestFunc,Func,TestLetter},
TestLetter=ELSclice[[LetterPos]];
Func=FuncIn;
Num=Exponent[Numerator[Func],t];
Den=Exponent[Denominator[Func],t];
SymFactor=1;
FoundFactor=True;
While[FoundFactor,
TestFunc=Factor[Func*TestLetter,Modulus->$P];
If[(Exponent[Denominator[TestFunc],t]<Den) && (Exponent[Numerator[TestFunc],t]<=Num),
	SymFactor=SymFactor/Fac[LetterPos];
	Func=TestFunc;
	Den=Exponent[Denominator[Func],t];,
	FoundFactor=False
];
];
FoundFactor=True;
While[FoundFactor,
TestFunc=Factor[Func/TestLetter,Modulus->$P];
If[(Exponent[Numerator[TestFunc],t]<Num) && (Exponent[Denominator[TestFunc],t]<=Den),
	SymFactor=SymFactor*Fac[LetterPos];
	Func=TestFunc;
	Num=Exponent[Numerator[Func],t];,
	FoundFactor=False
];
];
Return[SymFactor];
]

FoundFactors=Table[Times@@Map[CheckLetterFactor[Letter,#]&,Range[Length[ELSclice]]],{Letter,SlicedLetters}];



(Exponent[Numerator[TestFunc],t]<Num) && (Exponent[Denominator[TestFunc],t]<=Den)


FoundFactors[[-1]]


Factor[SlicedLetters[[-1]]/FoundFactors[[-1]],Modulus->$P]


Exponents=Table[
NewFunc=Factor[SlicedLetters[[i]]/FoundFactors[[i]]/.Fac[n_]:>ELSclice[[n]],Modulus->$P];
{Exponent[Numerator[NewFunc],t],Exponent[Denominator[NewFunc],t]},{i,1,Length[SlicedLetters]}]


FullreResolved=Flatten[Position[Exponents,{0,0}]]
ResolvedEvenLettersPos=Map[DeleteElements[Table[If[Log[#]===(dlogList[[i]]/.n_ Log[x_]:>Log[x]),i,0],{i,1,Length[dlogList]}],{0}]&,SlicedLetters[[FullreResolved]]]
UsedLetters=DeleteDuplicates[Cases[FoundFactors[[FullreResolved]],_Fac,Infinity]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]
SolvedMatrixPos=EvenPos[[Flatten[ResolvedEvenLettersPos]]]


Export[CurrentDir<>"DlogFactors.txt",UsedLetters]
Export[CurrentDir<>"Solved.txt",SolvedMatrixPos]




FoundFactors//Length
SlicedLetters//Length
EvenPos//Length


AllEvenLettersPos=Map[DeleteElements[Table[If[Log[#]===(dlogList[[i]]/.n_ Log[x_]:>Log[x]),i,0],{i,1,Length[dlogList]}],{0}]&,SlicedLetters]
FactorsWithPositions=Flatten[Table[Table[{pos,FoundFactors[[i]]},{pos,AllEvenLettersPos[[i]]}] ,{i,1,Length[AllEvenLettersPos]}],1];

Export[CurrentDir<>"Factors.txt",EvenLetterCandidatesFACS]
Export[CurrentDir<>"EvenLetterPos.txt",FactorsWithPositions]




