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


OutPutPaths={"Reconstruction-NewPara/Samples/ParaSamples-5000-1"};
SplitLengths={1};
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


LinearTerm[s16]


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
PrevReco=SortBy[Join[PrevReco1,PrevReco2,PrevReco3],First]
OddPosRem=Complement[OddPos,PrevReco[[All,1]]]
VarGrid=ReadList[HomeDir<>"/Reconstruction1/VarGrid.txt"];


OddPos[[1]]


G6=Get[HomeDir<>"/SetUp/SixPointGram.txt"];
DLogs=ReadList[CurrentDir<>"Dlogs.txt"];
IntegrandsReco=ReadList[CurrentDir<>"Integrands-Sorted.txt"];
pos1=14;
pos2=1;
RemPos=First[Flatten[Position[OddPos,{pos1,pos2}]]]
ExponentsNumi=Exponent[Numerator[IntegrandsReco[[RemPos,2]]],t]
ExponentsDeno=Exponent[Denominator[IntegrandsReco[[RemPos,2]]],t]

Entryusedvars=VarGrid[[pos1,pos2]]


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


IntegrandSlice=IntegrandsReco[[RemPos]][[2,7]]
SliceDenominator=Denominator[IntegrandSlice]/.root[x_]:>1


(* ::Subsection:: *)
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


(* ::Subsection:: *)
(*Letters*)


EvenLetterCandidates=ReadList[HomeDir<>"/Reconstruction-NewPara/EvenLetterCandidates.txt"];
ELSclice=Factor[EvenLetterCandidates/.SliceTpoint,Modulus->$P];
ELSclice//Length
EvenLetterCandidates//Length


DenoFactors=Times@@Map[CheckLetterFactor[SliceDenominator,#]&,Range[Length[EvenLetterCandidates]]]
DenominatorFromFactors=DenoFactors/.Fac[n_]:>EvenLetterCandidates[[n]]


(* ::Section:: *)
(*Reconstruction*)


Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
LaunchKernels[5];
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



FitVals[[4]]


FitVals[[1]]
IntegrandsReco[[1,2]]


fitvec=Expand[FitVals[[All,7]]*DenominatorFromFactorsN/(sqrt2N/sqrt1N),Modulus->$P];
fitvec=Expand[FitVals[[All,7]]*1/(sqrt2N/sqrt1N),Modulus->$P];

(*fitvec=Expand[FitVals[[All,2]]/(sqrt2N/sqrt1N),Modulus->$P]*)


Letters=ReadList["/home/daniel/Documents/PHD-LocalFiles/TriangleHeptagon3/Reconstruction-Odd-Letters/OddLetterReco1AB.txt"];
testA=Letters[[1,2]]
testB=Letters[[1,3]]
testR=Letters[[1,4]]
Dlog=Log[(testA+testB testR)/(testA-testB testR)];
TestEnt=Factor[ExpandDenominator[Factor[D[Dlog,s12]/.Derivative[1][root][x_]:>1/(2root[x])]]/.root[x_]^3:>root[x]x/.root[x_]^2:>x]
Factor[Expand[TestEnt/.SliceTpoint,Modulus->$P],Modulus->$P]/.root[x_]:>root[Factor[x,Modulus->$P]]
Test=Factor[Expand[TestEnt/.PhaseSpacePoints[[-1]],Modulus->$P],Modulus->$P]/.root[x_]:>root[Factor[x,Modulus->$P]]


Expand[Test/.root[n_]:>PowerMod[n,1/2,$P],Modulus->$P]
FitVals[[1]]


usedvars


TestEnt1=LinearTerm[s12][[All,1,1]];
Expand[TestEnt1*s12/.PhaseSpacePoints[[1]],Modulus->$P]


ReconstructFUNCS[fitvec,PhaseSpacePoints[[All,All,2]],usedvars,1,1,Debug->True]


?ReconstructFUNCS
