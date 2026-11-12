(* ::Package:: *)

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
OutPutPath="Samples/Univariate-Slice-1";
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
rootVector=Map[Times@@DeleteDuplicates[Cases[#,_root,Infinity]]&,PureIntegrals]/.{}->1/.root[x_]:>x;
DeFoldMat=Normal[Total[Table[ArrayReshape[LinearTerm[var][[1,1]],{Length[PureIntegrals],Length[PureIntegrals]}],{var,usedvars}]]]/.n_Integer/;(n>0):>1;
rootVectorN=rootVector/.PhaseSpacePoints[[1]];
RootRatioMat=Expand[Table[Together[rootVectorN[[i1]]/rootVectorN[[i2]]],{i1,1,Length[PureIntegrals]},{i2,1,Length[PureIntegrals]}],Modulus->$P]*DeFoldMat;
EvenPos=Position[RootRatioMat,1];
LinearizedEvenPos=Map[GetLinPos[#[[1]],#[[2]],95]&,EvenPos];
Get[HomeDir<>"/4DBasis/Parametrization.m"]


Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
pos=GetLinPos[93,90,95];

avec=ReadList[CurrentDir<>"AVECS.txt"][[1]];
bvec=ReadList[CurrentDir<>"AVECS.txt"][[2]];
LinearTerm[s16]=SparseArray[ConstantArray[0,{Length[PhaseSpacePoints],1,Length[PureIntegrals]^2}]];
FitVals=Join[Table[LinearTerm[var][[All,1,pos]],{var,diffInvariants}]];
FitVals=Normal[Transpose[FitVals]];
SliceTpoint=Thread[Rule[avec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]]

(*Integrands=Table[Integrand=Expand[UniVarReconstructFUNCS[FitVals[[varpos]],PhaseSpacePoints,tvar,20,20,{},1,1],Modulus->$P];Factor[Integrand[[1]],Modulus->$P],{varpos,1,10}]
NewIntDt=Factor[Expand[Integrands . bvec[[All,2]],Modulus->$P],Modulus->$P]*)



(*avec=ReadList[HomeDir<>"/Reconstruction1/Samples/AVECS.txt"][[1]];
bvec=ReadList[HomeDir<>"/Reconstruction1/Samples/AVECS.txt"][[2]];*)
SliceTpoint=Thread[Rule[avec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]]


phasespacepointsPara=Map[ConvertMandsToPoint[#]&,PhaseSpacePoints];
phasespacepointsPara=Table[Join[phasespacepointsPara[[i]],{tvarlist[[i]]}],{i,1,Length[phasespacepointsPara]}];

FitValsPara=Table[ConvertDF[FitVals[[i]],phasespacepointsPara[[i]]],{i,1,Length[FitVals]}];


Integrands=Table[Integrand=Expand[UniVarReconstructFUNCS[FitValsPara[[All,varpos]],phasespacepointsPara,tvar,20,20,{},1,1],Modulus->$P];Factor[Integrand[[1]],Modulus->$P],{varpos,1,9}];
NewIntDt=Factor[Expand[Integrands . bvec[[All,2]],Modulus->$P],Modulus->$P]


Get["/home/dschiebe/Documents/packages/trager/Trager.wl"]
Dlog=TragerLogSubtractionNumeric[NewIntDt,t,{},{},Modulus->$P]


EvenLetterCandidates=DeleteCases[ReadList[HomeDir<>"/ClaudeLetters/EvenLetterCandidates.txt"],Null];
EvenLetterCandidatesFACS=DeleteElements[DeleteDuplicates[Flatten[Map[FactorList[#][[All,1]]&,EvenLetterCandidates]]],{1,-1}];
Position[Factor[EvenLetterCandidatesFACS/.SliceTpoint,Modulus->$P],0];
EvenLetterCandidatesFACS=Delete[EvenLetterCandidatesFACS,79];
ELSclice=Factor[EvenLetterCandidatesFACS/.SliceTpoint,Modulus->$P];


NewDlogs1=ReadList[HomeDir<>"/Reconstruction1/NewDlogs-500.txt"];
NewDlogs2=ReadList[HomeDir<>"/Reconstruction1/NewDlogs-500-2.txt"];
NewDlogs3=ReadList[HomeDir<>"/Reconstruction1/FinishedDlog-94-94.txt"];
NewDlogs4=ReadList[HomeDir<>"/Reconstruction1/FinishedDlog-95-95.txt"];
fiveptgrams=ReadList[HomeDir<>"/Reconstruction1/grams5pt.txt"];
EvenLetterCandidatesFACSNew=DeleteElements[DeleteDuplicates[Flatten[Map[FactorList[#][[All,1]]&,DeleteDuplicates[Cases[{NewDlogs1,NewDlogs2,NewDlogs3,NewDlogs4},_Log,Infinity]]/.Log[x_]:>x]]],{1,-1}];
EvenLetterCandidatesFACSNew=DeleteDuplicates[Join[EvenLetterCandidatesFACSNew,fiveptgrams]];
EvenLetterCandidatesFACSNew//Length
FancyDeleteDuplicates[delset_]:=Module[{size,i,j,delset2,factor},delset2=delset;size=Length[delset2];
For[i=1,i<=size,i++,
	For[j=i,j<=size,j++,factor=Factor[delset2[[i]]/delset2[[j]]];
		If[NumberQ[factor],delset2[[j]]=delset2[[i]];];];];
Return[DeleteDuplicates[delset2]];]
EvenLetterCandidatesFACS=FancyDeleteDuplicates[Join[EvenLetterCandidatesFACS,EvenLetterCandidatesFACSNew]];
ELSclice=Factor[EvenLetterCandidatesFACS/.SliceTpoint,Modulus->$P];
ELSclice//Length
EvenLetterCandidatesFACS//Length


EvenLetterCandidatesFACS=ttbbFamilyDlog[["EvenLetters"]];
LettersPara=Factor[EvenLetterCandidatesFACS/.KinematicRulesBackWards]


LettersP1=DeleteElements[FactorList[Times@@Numerator[LettersPara]][[All,1]],{1,-1}]
LettersP2=DeleteElements[FactorList[Times@@Denominator[LettersPara]][[All,1]],{1,-1}]
LettersInParametrization=Join[LettersP1,LettersP2]


LettersInParametrization=LettersInParametrization//DeleteDuplicates




Test=SortBy[LettersInParametrization,Length]


Test[[1;;40]]


LetterCandidatesSliced=Factor[LettersInParametrization/.SliceTpoint,Modulus->$P];


LetterCandidatesSliced//Length


LogArg=Dlog/.Log[x_]:>Factor[x,Modulus->$P]


pos=1;
CheckLetterFactor[FuncIn_,LetterPos_]:=Module[{Num,Den,SymFactor,FoundFactor,TestFunc,Func,TestLetter},
TestLetter=LetterCandidatesSliced[[LetterPos]];
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
KnownFacs=Times@@Map[CheckLetterFactor[LogArg,#]&,Range[Length[LetterCandidatesSliced]]]


Get[HomeDir<>"/4DBasis/FactorAlphabetSolver.wl"]
KnownLogArg=KnownFacs/.Fac[n_]:>LettersInParametrization[[n]]
FactorsInMandLetters=Map[DeleteDuplicates[DeleteElements[FactorList[#][[All,1]],{1,-1}]]&,LettersPara]
PossibleLetters=DeleteDuplicates[Flatten[Map[Position[FactorsInMandLetters,FactorList[KnownLogArg[[#]]][[2,1]],{2}][[All,1]]&,Range[Length[KnownLogArg]]]]]
BuildMandFactors=SolveFactorPowers[KnownLogArg,LettersPara]
Powers=BuildMandFactors[["exponents"]][[All,2]]
LogFacInMands=Factor[Times@@Power[EvenLetterCandidatesFACS,Powers]];





Export[CurrentDir<>"KnowLogFac-93-90.txt",LogFacInMands]


Expand[LogFacInMands/.PhaseSpacePoints[[1]],Modulus->$P]
Expand[KnownLogArg/.phasespacepointsPara[[1]],Modulus->$P]


LogFacInMands


Factor[LogArg/KnownFacs/.Fac[n_]:>LetterCandidatesSliced[[n]],Modulus->$P]


Together[KnownLogArg]
LogFacInMands


KnownLogArg
LogFacInMands


KnownLogArgMyVersion=Get["/home/dschiebe/Documents/PhD-Files/TriangleHeptagon3/Reconstruction-NewPara/KnowLogFac-93-90.txt"];


MyFac=Factor[KnownLogArgMyVersion/.KinematicRulesBackWards]
KnownLogArg


RemFacToReco=Factor[MyFac/KnownLogArg]
RemFacToReco//Variables


Expand[Factor[Test=Factor[D[Log[RemFacToReco],s56]]/.SliceTpoint,Modulus->$P],Modulus->$P]
Test=Factor[D[Log[RemFacToReco],s12]]


DuoVariateSlice
point1=RandomInteger[$P,9]
point2=RandomInteger[$P,9]
BiVariateSlice=point1+t Join[point2[[;;-4]],{0,0,0}]+c Join[{0,0,0,0,0,0},point2[[-3;;]]];
BiVariateSliceTpoint=Thread[Rule[paravars,BiVariateSlice]]


Factor[Test/.BiVariateSliceTpoint,Modulus->$P]


(* ::Section:: *)
(*Playground*)


Tpoint=SliceTpoint/.t->0
Factor[Test/.Tpoint[[2;;]],Modulus->$P]
Expand[Factor[Test/.Tpoint[[Delete[Range[9],9]]],Modulus->$P],Modulus->$P]

Get[CurrentDir<>"AnsatzBuilder.wl"]


?BuildRationalAnsatz


Factor[BackTransformation=Test/.KinematicRulesForward]


avec=ReadList[HomeDir<>"/Reconstruction1/Samples/AVECS.txt"][[1]];
bvec=ReadList[HomeDir<>"/Reconstruction1/Samples/AVECS.txt"][[2]];
SliceTpointMand=Thread[Rule[avec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]]

MandsTest=Factor[D[Log[Factor[RemFacToReco/.KinematicRulesForward]],mt2]]


Expand[Factor[MandsTest/.SliceTpointMand,Modulus->$P],Modulus->$P]


(* ::Section:: *)
(*Sebastian Letters*)


DlogBasis=Get["/home/dschiebe/Downloads/ttbb_family_dlog_bare.m"]


Test=Position[ttbbFamilyDlog["EvenMatrices"],{206,206}]


ttbbFamilyDlog["EvenMatrices"][[7,363]]
ttbbFamilyDlog["EvenMatrices"][[86,30]]



