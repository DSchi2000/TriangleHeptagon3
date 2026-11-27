(* ::Package:: *)

If[$FrontEnd===Null,$DIR=Directory[],$DIR=NotebookDirectory[]];
parts = FileNameSplit[$DIR];
pos = First@FirstPosition[parts, "TriangleHeptagon3"];
HomeDir = FileNameJoin[parts[[;;pos ]]];
SectorDir = FileNameJoin[parts[[;;pos+1]]];
CurrentDir=$DIR;
kirapath=SectorDir;
$P=2^31-1;
diffInvariants={mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56};
PureIntegrals=ReadList[FileNameJoin[{SectorDir,"SubBasis.txt"}]];
SetDirectory[CurrentDir];
NotbookDir=NotebookDirectory[];
Get[FileNameJoin[{HomeDir,"files/ExtraFunctions.m"}]];
SQRTS=Cases[PureIntegrals,_root,Infinity]//DeleteDuplicates;

rootSet=(Cases[PureIntegrals[[-5;;]],_root,Infinity]//DeleteDuplicates)/.root[x_]:>x;
GetLinPos[x1_,x2_,NMI_]:=First[Flatten[Position[Flatten[Table[{i,j},{i,1,NMI},{j,1,NMI}],1],{x1,x2}]]];
Get[ToString[StringJoin[HomeDir,"/files/ExtraModules.m"]]];

usedvars={mt2,s12,s123,s23,s234,s34,s345,s45,s56};

(*OutPutPath="Samples/Full-Samples-500";
phasespacepoints=ReadList[FileNameJoin[{CurrentDir,ToString[StringForm["``/phasespacepoints.txt",OutPutPath]]}]];
phasespacepoints=phasespacepoints[[All,2;;]][[1;; ;; 1]];
PhaseSpacePoints=Map[Thread[Rule[diffInvariants,#]]&,phasespacepoints];
For[i=1,i<=Length[usedvars],i++,
LinearTerm2[usedvars[[i]]]=ReadList[FileNameJoin[{Directory[],ToString[StringForm["``/LinTerm``_0001.txt",OutPutPath,usedvars[[i]]]]}]];
];*)


(* ::Section:: *)
(*Load Samples*)


OutPutPaths={"Reconstruction-NewPara/Samples/ParaSamples-10000-1"};
SplitLengths={2};
PhaseSpacePoints={};
Do[
phasespacepointsIN=ReadList[FileNameJoin[{HomeDir,ToString[StringForm["``/phasespacepoints.txt",OutPutPaths[[q]]]]}]];
phasespacepointsIN=phasespacepointsIN[[All,2;;]];
PhaseSpacePointsIN=Map[Thread[Rule[diffInvariants,#]]&,phasespacepointsIN];
PhaseSpacePoints=Join[PhaseSpacePoints,PhaseSpacePointsIN]
,{q,1,Length[OutPutPaths]}]

loadfile=1;
Do[
OutPutPath=OutPutPaths[[filepos]];
Print["Reading From: "<>OutPutPath];
Monitor[For[q=1,q<=SplitLengths[[filepos]],q++,
	For[i=1,i<=Length[usedvars],i++,
	LinearTerm[usedvars[[i]],loadfile]=ReadList[FileNameJoin[{HomeDir,ToString[StringForm["``/LinTerm``_``.txt",OutPutPath,usedvars[[i]],IntegerString[q,10,4]]]}]][[;;-1]];
	If[q<SplitLengths[[filepos]],LinearTerm[usedvars[[i]],loadfile]=LinearTerm[usedvars[[i]],loadfile][[;;-2]]];
	];
loadfile+=1;
],"Subfile: "<>ToString[q]<>" var: "<>ToString[usedvars[[i]]]];
,{filepos,1,Length[OutPutPaths]}]

(*PhaseSpacePoints=PhaseSpacePoints[[;;3000]];*)
Map[LinearTerm[#]=(Join@@Table[LinearTerm[#,q],{q,1,Total[SplitLengths]}]);&,usedvars];
LinearTerm[s16]=SparseArray[ConstantArray[0,{Length[PhaseSpacePoints],1,Length[PureIntegrals]^2}]];
LinearTerm[s12]//Length
PhaseSpacePoints//Length


Expand[LinearTerm[s12][[;;20,1,1]]*(PhaseSpacePoints[[;;20,2,2]]),Modulus->$P]


(* ::Section::Closed:: *)
(*Load Samples UVar for Debug*)


OutPutPath="Reconstruction1/Samples/UnivariateSlice-3";
phasespacepoints=ReadList[FileNameJoin[{HomeDir,ToString[StringForm["``/phasespacepoints.txt",OutPutPath]]}]];
tvarlist=Flatten[ReadList[FileNameJoin[{HomeDir,ToString[StringForm["``/tvars.txt",OutPutPath]]}]][[;;;; 1]]];
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
LinearTerm[usedvars[[i]]]=ReadList[FileNameJoin[{HomeDir,ToString[StringForm["``/LinTerm``_0001.txt",OutPutPath,usedvars[[i]]]]}]];
];
LinearTerm[s16]=SparseArray[ConstantArray[0,{Length[PhaseSpacePoints],1,Length[PureIntegrals]^2}]];
LinearTerm[s12]//Length
PhaseSpacePoints//Length


(* ::Section:: *)
(*Variable Shift Tech*)


G6=Get[HomeDir<>"/SetUp/SixPointGram.txt"];
ChangeAuxVarInDe[varOld_,varNew_,DIsij_,tpoint_]:=Module[{newaux,oldaux,newcoefs,ChangeVarpos,DIsijSUB,newcoefsSUB,posA,posB,EqnConst},
EqnConst=Sum[(Expand[D[G6,var]/.tpoint,Modulus->$P]) ds[var],{var,diffInvariants}];

newaux=First[Solve[EqnConst==0,ds[varNew],Modulus->$P]];
oldaux=First[Solve[EqnConst==0,ds[varOld],Modulus->$P]];
newcoefs=Coefficient[newaux[[1,2]],Table[ds[var],{var,diffInvariants}]];
ChangeVarpos=Sort[Map[Flatten[Position[diffInvariants,#]]&,{varOld,varNew}]];
DIsijSUB=Delete[DIsij,ChangeVarpos];
newcoefsSUB=Delete[newcoefs,ChangeVarpos];
posA=First[Flatten[Position[diffInvariants,varOld]]];
posB=First[Flatten[Position[diffInvariants,varNew]]];

Return[Expand[(DIsijSUB+newcoefsSUB DIsij[[posB]]),Modulus->$P] . Table[ds[var],{var,Delete[diffInvariants,ChangeVarpos]}]+ Expand[DIsij[[posB]]newcoefs[[posA]],Modulus->$P] ds[varOld]];

]
pos=GetLinPos[34,34,95];
DIsij=Map[LinearTerm[#][[1,1,pos]]&,diffInvariants]
(*DIsij=Insert[DIsij,0,4]*)
tpoint=PhaseSpacePoints[[1,;;-1]]
ChangeAuxVarInDe[s16,s12,DIsij,tpoint]


(* ::Section:: *)
(*Load Extras*)


rootVector=Map[Times@@DeleteDuplicates[Cases[#,_root,Infinity]]&,PureIntegrals]/.{}->1/.root[x_]:>x;
DeFoldMat=Normal[Total[Table[ArrayReshape[LinearTerm[var][[1,1]],{Length[PureIntegrals],Length[PureIntegrals]}],{var,usedvars}]]]/.n_Integer/;(n>0):>1;
rootVectorN=rootVector/.PhaseSpacePoints[[1]];
RootRatioMat=Expand[Table[Together[rootVectorN[[i1]]/rootVectorN[[i2]]],{i1,1,Length[PureIntegrals]},{i2,1,Length[PureIntegrals]}],Modulus->$P]*DeFoldMat;
EvenPos=Position[RootRatioMat,1];
(*LinearizedEvenPos=Map[GetLinPos[#[[1]],#[[2]],95]&,EvenPos]*)
OddPos=Complement[Position[DeFoldMat,1],EvenPos];

PrevReco1=ReadList[CurrentDir<>"OddLetterReco1AB.txt"];
PrevReco2=ReadList[CurrentDir<>"OddLetterReco2AB.txt"];
PrevReco3=ReadList[CurrentDir<>"OddLetterReco2AB-PARA.txt"];

PrevReco1=Delete[PrevReco1,DeleteDuplicates[Position[PrevReco1,$Failed][[All,1]]]/.n_Integer:>{n}];
PrevReco2=Delete[PrevReco2,DeleteDuplicates[Position[PrevReco2,$Failed][[All,1]]]/.n_Integer:>{n}];
PrevReco3=Delete[PrevReco3,DeleteDuplicates[Position[PrevReco3,$Failed][[All,1]]]/.n_Integer:>{n}];
files=FileNames["DLog*.txt", CurrentDir<>"HandRecoLetters"];
PrevReco4=Map[ReadList[#]&,files];

PrevReco=SortBy[Join[PrevReco1,PrevReco2,PrevReco3,PrevReco4],First]
OddPosRem=Complement[OddPos,PrevReco[[All,1]]]
VarGrid=ReadList[HomeDir<>"/Reconstruction1/VarGrid.txt"];


(* ::Section:: *)
(*Linear Dependencies*)


Unprotect[RowReduce];
Clear[RowReduce];
Protect[RowReduce];
rtfactorsMap=Map[{#,root[rootVector[[#[[1]]]]/rootVector[[#[[2]]]]]}&,OddPos];
rtfactors=DeleteDuplicates[Cases[rtfactorsMap,_root,Infinity]]//EchoFunction[Length];
SameRootFacPos=OddPos[[Position[rtfactorsMap,rtfactors[[38]]][[All,1]]]];

Print["Unknowns:"]
unknowns=Complement[SameRootFacPos,PrevReco[[All,1]]]
Print["Knowns:"]
knowns=Intersection[SameRootFacPos,PrevReco[[All,1]]]
SameRootFacPos=Join[knowns,unknowns];

LinPos=Map[GetLinPos[#[[1]],#[[2]],95]&,SameRootFacPos];
tpoint=Last[PhaseSpacePoints][[All,2]];
MatSpace=Expand[Sum[tpoint[[i]]Normal[LinearTerm[diffInvariants[[i]]][[;;Length[LinPos],1,LinPos]]],{i,1,Length[diffInvariants]}],Modulus->$P];


Rech=RowReduce[MatSpace,Modulus->$P]/.n_Integer:>RationalReconstruct[n];
M=MatSpace;
Rech//MatrixPlot
pivotCols = DeleteCases[FirstPosition[#, x_ /; x != 0] & /@ Rech, _Missing][[All, 1]]
LinDeps=Table[Rech[[All,q]] . Table[Ent[i],{i,1,Length[MatSpace]}],{q,1,Length[MatSpace]}];
basis=SameRootFacPos[[pivotCols]]//EchoFunction[Length]
LinDeps=LinDeps/.Ent[n_]:>Ent[SameRootFacPos[[pivotCols]][[n]]];
Complement[basis,PrevReco[[All,1]]]


Map[Flatten[Position[PrevReco[[All,1]],#]]&,basis]


(* ::Subsection:: *)
(*Todo*)


{{75,76}},{{90,71}},{{90,72}},{{76,75},{93,75}}


(* ::Subsection::Closed:: *)
(*Old*)


LinDepsEQN=DeleteElements[Thread[Thread[Ent[SameRootFacPos]]==LinDeps],{True}];
solset=DeleteDuplicates[Cases[LinDepsEQN,_Ent,Infinity]];
solset=Thread[Ent[SameRootFacPos]];
unsolved=Complement[solset,Thread[Ent[PrevReco[[All,1]]]]];
solsetSorted=Join[unsolved,Complement[solset,unsolved]];
(*sol=Solve[LinDepsEQN,solsetSorted][[1]]*)
Basis=DeleteDuplicates[Cases[sol[[All,2]],_Ent,Infinity]]//EchoFunction[Length]
Intersection[Basis,unsolved]//EchoFunction[Length]


(* ::Section:: *)
(*Reco Letter*)


G6=Get[HomeDir<>"/SetUp/SixPointGram.txt"];
DLogs=ReadList[CurrentDir<>"Dlogs.txt"];
IntegrandsReco=ReadList[CurrentDir<>"Integrands-Sorted.txt"];
pos1=90;
pos2=45;
RemPos=First[Flatten[Position[OddPos,{pos1,pos2}]]]
ExponentsNumi=Exponent[Numerator[IntegrandsReco[[RemPos,2]]],t]
ExponentsDeno=Exponent[Denominator[IntegrandsReco[[RemPos,2]]],t]

Entryusedvars=VarGrid[[pos1,pos2]]

rootVector[[pos1]]//Variables//Sort
rootVector[[pos2]]//Variables//Sort


(* ::Section:: *)
(*Get Denominators:*)


(* ::Subsection::Closed:: *)
(*Functions*)


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


(* ::Subsection:: *)
(*Things*)


(* ::Subsection::Closed:: *)
(*Uvar*)


OutPutPathUV="Reconstruction1/Samples/UnivariateSlice-3";
(*OutPutPath="Reconstruction1/Samples/Univariate Slice";*)
phasespacepointsUV=ReadList[FileNameJoin[{HomeDir,ToString[StringForm["``/phasespacepoints.txt",OutPutPathUV]]}]];
tvarlist=Flatten[ReadList[FileNameJoin[{HomeDir,ToString[StringForm["``/tvars.txt",OutPutPathUV]]}]][[;;;; 1]]];
phasespacepointsUV=phasespacepointsUV[[All,2;;]][[1;; ;; 1]];
phasespacepointsUV=Map[Thread[Rule[diffInvariants,#]]&,phasespacepointsUV];
phasespacepointsUV=Table[Join[phasespacepointsUV[[i]],{tvarlist[[i]]}],{i,1,Length[phasespacepointsUV]}];

EQN=Flatten[Table[{avector[i]+phasespacepointsUV[[1,-1,2]] bvector[i] ==phasespacepointsUV[[1,i,2]],avector[i]+phasespacepointsUV[[2,-1,2]] bvector[i] ==phasespacepointsUV[[2,i,2]]},{i,1,Length[phasespacepointsUV[[1,All,1]]]-1}]];
solset=Flatten[Table[{avector[i],bvector[i]},{i,1,Length[phasespacepointsUV[[1,All,1]]]-1}]];
sol=Solve[EQN,solset,Modulus->$P][[1]];
avec=Thread[Rule[phasespacepointsUV[[1,;;-2,1]],Table[avector[i],{i,1,Length[phasespacepointsUV[[1,All,1]]]-1}] ]]/.sol;
bvec=Thread[Rule[phasespacepointsUV[[1,;;-2,1]],Table[bvector[i],{i,1,Length[phasespacepointsUV[[1,All,1]]]-1}] ]]/.sol;
(*LinearTerm[s16]=SparseArray[ConstantArray[0,{Length[phasespacepointsUV],1,Length[PureIntegrals]^2}]];*)
SliceTpoint=Thread[Rule[bvec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]];
SliceTpoint=Thread[Rule[bvec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]];
Expand[SliceTpoint/.t->phasespacepointsUV[[-1,-1,2]],Modulus->$P]
phasespacepointsUV[[-1]]


(* ::Subsection::Closed:: *)
(*Letters*)


(*EvenLetterCandidates=ReadList[HomeDir<>"/Reconstruction-NewPara/EvenLetterCandidates.txt"];*)
EvenLetterCandidates=ReadList[HomeDir<>"/FitLetters/UsedEvenLetters.txt"];
ELSclice=Factor[EvenLetterCandidates/.SliceTpoint,Modulus->$P];
ELSclice//Length
EvenLetterCandidates//Length


(* ::Subsection:: *)
(*Fit Deno*)


intvar=s34
intvarpos=First[Flatten[Position[diffInvariants,intvar]]];
numpower=ExponentsNumi[[intvarpos]];
IntegrandSlice=IntegrandsReco[[RemPos]][[2,intvarpos]]
SliceDenominator=Denominator[IntegrandSlice]/.root[x_]:>1
DenoFactors=Times@@Map[CheckLetterFactor[SliceDenominator,#]&,Range[Length[EvenLetterCandidates]]]
DenominatorFromFactors=DenoFactors/.Fac[n_]:>EvenLetterCandidates[[n]]
SliceDenominator/Factor[DenominatorFromFactors/.SliceTpoint,Modulus->$P]


(* ::Section:: *)
(*Reconstruction*)


LaunchKernels[5];


Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
sqrt1=rootVector[[pos1]];
sqrt2=rootVector[[pos2]];
pos=GetLinPos[pos1,pos2,95];
sqrt1N=PowerMod[sqrt1/.PhaseSpacePoints,1/2,$P];
sqrt2N=PowerMod[sqrt2/.PhaseSpacePoints,1/2,$P];
DenominatorFromFactorsN=Expand[DenominatorFromFactors/.PhaseSpacePoints,Modulus->$P];
FitVals=Join[Table[LinearTerm[var][[All,1,pos]],{var,diffInvariants}]];
usedvars=VarGrid[[pos1,pos2]];
If[Length[usedvars]<10,newaux=Last[Complement[diffInvariants,usedvars]];
FitVals=ParallelTable[Coefficient[ChangeAuxVarInDe[s16,newaux,Transpose[FitVals][[i]],PhaseSpacePoints[[i]]],Thread[ds[diffInvariants]]],{i,1,Length[PhaseSpacePoints]}];,
FitVals=Transpose[FitVals]];
(*Integrands=Table[Integrand=Expand[UniVarReconstructFUNCS[Expand[FitVals[[All,varpos]] DenominatorFromFactorsN/(sqrt2N/sqrt1N),Modulus->$P],PhaseSpacePoints,tvar,20,20,{},1,1],Modulus->$P];Factor[Integrand[[1]],Modulus->$P],{varpos,1,10}];
Integrands=Integrands*(root[sqrt2]/root[sqrt1])/.root[1]->1/.SliceTpoint/.root[x_]:>root[Factor[x,Modulus->$P]];*)



(*intvar=s34;*)
fitvec=Expand[FitVals[[All,First[Flatten[Position[diffInvariants,intvar]]]]]*DenominatorFromFactorsN/(sqrt2N/sqrt1N),Modulus->$P];
Reco=ReconstructFUNCS[fitvec,PhaseSpacePoints[[All,Flatten[Map[Position[diffInvariants,#]&,usedvars]],2]],usedvars,numpower,0,Debug->True]
Func=BuildFunction[Reco[[1,1]],Reco[[2]]]/.n_Integer:>RationalReconstruct[n]
FullFunc=Func/DenominatorFromFactors *(root[sqrt2]/root[sqrt1])/.root[1]->1


(* ::Section:: *)
(*Integrate*)


Get["/home/daniel/Documents/packages/trager/Trager.wl"]
roots=DeleteDuplicates[Cases[FullFunc,_root,Infinity]]
rootrepls=Thread[Rule[roots,Map[rt[#]&,Range[Length[roots]]]]]
qring=Expand[Table[rt[i]^2-(roots[[i]]/.root[x_]:>x),{i,Range[Length[roots]]}]]
rtvars=Table[rt[i],{i,Range[Length[roots]]}]
invrootrepls = Reverse /@ rootrepls;
Dlog=TragerLogSubtraction[FullFunc/.rootrepls,intvar,rtvars,qring,Debug->True,Parallelize->True]
Dlog=Dlog/.invrootrepls;


Num=Collect[Numerator[First[Cases[{Dlog},_Log,Infinity]/.Log[x_]:>x]],_root,Factor];
usedroots=Times@@Cases[Num,_root,Infinity]
bvac=Factor[Coefficient[Num,usedroots]]
avac=Factor[Num/.root[x_]:>0]
Export[CurrentDir<>"HandRecoLetters/"<>ToString[StringForm["DLog-``-``.txt",pos1,pos2]],{{pos1,pos2},avac,bvac,usedroots}]
Export[CurrentDir<>"HandRecoLetters/"<>ToString[StringForm["IntegratedFunc-``-``.txt",pos1,pos2]],{{pos1,pos2},Dlog}]



(* ::Section:: *)
(*TestZone*)


Get["/home/daniel/Documents/PHD-LocalFiles/TriangleHeptagon3/SetUp/RunSimplifyOnGram.wls"]


Test=Factor[PrevReco3[[21,2]]/.KinematicRulesForward]//EchoFunction[ByteCount]//EchoTiming
simplified = SimplifyOnVarietyFactored[Test]//EchoFunction[ByteCount]//EchoTiming
simplifiedFull = SimplifyOnVarietyFactoredBest[Test]//EchoFunction[ByteCount]//EchoTiming

Expand[Test/.PhaseSpacePoints[[1]],Modulus->$P]
Expand[simplified/.PhaseSpacePoints[[1]],Modulus->$P]
Expand[simplifiedFull/.PhaseSpacePoints[[1]],Modulus->$P]


LeafCount[simplifiedFull]
LeafCount[simplified]


Export[CurrentDir<>"TestForClaude.txt",Test]





test=s16^2 s12^2
test2=Factor[Factor[test/.KinematicRulesBackWards]/.KinematicRulesForward]//EchoFunction[ByteCount];


simplified = SimplifyOnVarietyFactored[test2]



TestFacs=DeleteElements[FactorList[Test][[All,1]],{-1,1}]//EchoFunction[Length]





Exponent[TestFacs[[3]],s16]


GF2=Coefficient[G6,s16^2]
CF2=Coefficient[TestFacs[[3]],s16^2]
Q=Factor[CF2/GF2]


Factor[Factor[TestFacs[[3]]-Q G6]]TestFacs[[1]]/TestFacs[[2]]


simplified


Factor[Test+Q G6]


Factor[mt2/.KinematicRulesBackWards/.KinematicRulesForward]
