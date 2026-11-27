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
(*Get[FileNameJoin[{HomeDir,"files/ExtraFunctions.m"}]];*)
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


(*OutPutPaths={"Reconstruction1/Samples/Full-Samples-10000"};*)
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
VarGrid=ReadList[HomeDir<>"/Reconstruction1/VarGrid.txt"];


OddPos[[1]]


Get[HomeDir<>"/4DBasis/Parametrization.m"]
(*TestDLog=Get[HomeDir<>"/Reconstruction1/FinishedDlog-95-95.txt"];*)
TestDLog=Get[HomeDir<>"/Reconstruction-NewPara/FinishedDlog-93-93.txt"];
(*TestDLog=Factor[TestDLog/.KinematicRulesForward];*)
TestDLog/.Log[x_]:>Variables[x]





var=mt2;
EqnConstFull=Sum[D[G6,var] ds[var],{var,diffInvariants}];
diffauxvar=First[Solve[EqnConstFull==0,ds[s16]]];
sixpointdiff=Coefficient[diffauxvar[[1,2]],ds[var]];


Diff1=Factor[D[TestDLog,var]];
Diff2=Factor[D[TestDLog,s16]];
Entry=Diff1+sixpointdiff Diff2;


A=Expand[Entry/.PhaseSpacePoints[[;;20]],Modulus->$P]


linpos=GetLinPos[93,93,95]
B=LinearTerm[mt2][[;;20,1,linpos]]


Expand[A/B,Modulus->$P]/.n_Integer:>RationalReconstruct[n]





Point1=ConvertMandsToPoint[PhaseSpacePoints[[1]]]
linpos=GetLinPos[1,1,95]
DF=Table[LinearTerm[var][[1,1,linpos]],{var,diffInvariants}]
NEWDF=ConvertDF[DF,Point1]


(*NewTestEntry=D[TestDLog,c1];*)
NewTestEntry=D[Log[1/s12],s12];


Test=Expand[NewTestEntry/.Point1,Modulus->$P]


Expand[NEWDF/Test,Modulus->$P]/.n_Integer:>RationalReconstruct[n]


(* ::Section:: *)
(*Load Letters*)


EvenLetterCands=ReadList[HomeDir<>"/Reconstruction-NewPara/EvenLetterCandidates.txt"];
MissingDlog=Get[HomeDir<>"/Reconstruction-NewPara/DlogM-93-93.txt"];
MissingDlog=MissingDlog/.Log[x_]:>x;
MissingFacs=DeleteElements[FactorList[MissingDlog][[All,1]],{-1,1}];
EvenLetterCands=DeleteDuplicates[Join[EvenLetterCands,MissingFacs]]


(* ::Section:: *)
(*Linear Independence*)


ranpt=Delete[PhaseSpacePoints[[1]],4][[All,2]]


var=s34;
PhaseSpacePointsCut=PhaseSpacePoints[[;;200]];
NullSpaceMat=Table[Expand[EvenLetterCands/.PhaseSpacePointsCut[[i]],Modulus->$P],{i,1,Length[EvenLetterCands]}];

NullSpaceMats=Table[
sixpointdiff=Coefficient[diffauxvar[[1,2]],ds[var]];
Table[Expand[(D[Log[EvenLetterCands],var]+sixpointdiff D[Log[EvenLetterCands],s16])/.PhaseSpacePointsCut[[i]],Modulus->$P],{i,1,Length[EvenLetterCands]}]
,{var,DeleteElements[diffInvariants,{s16}]}];


NullSpaceMat=Expand[Sum[ranpt[[i]]NullSpaceMats[[i]],{i,1,9}],Modulus->$P];
rech=RowReduce[NullSpaceMat,Modulus->$P];
rech//MatrixPlot
MatrixRank[NullSpaceMat,Modulus->$P]
pivotCols = FirstPosition[#, x_ /; x != 0] & /@ rech;
pivotCols = DeleteCases[pivotCols, _Missing][[All, 1]];

(* independent COLUMNS of the original matrix M *)
(*basisColumns = NullSpaceMatBasisMat[[All, pivotCols]]; *)  (* 191 x 91 *)
IndependentEvenLetters=EvenLetterCands[[pivotCols]]//EchoFunction[Length];


(* ::Section:: *)
(*Fit one*)


Close[Kernels]
LaunchKernels[5];


(*Get[FileNameJoin[{HomeDir,"files/ExtraFunctions.m"}]];*)
var=s12;
EqnConstFull=Sum[D[G6,var] ds[var],{var,diffInvariants}];
diffauxvar=First[Solve[EqnConstFull==0,ds[s16]]];
sixpointdiff=Coefficient[diffauxvar[[1,2]],ds[var]];
linpos=GetLinPos[95,90,95];
fitvec=LinearTerm[var][[;;200,1,linpos]];
PhaseSpacePointsCut=PhaseSpacePoints[[;;200]];
NullSpaceMat=ParallelTable[Join[Expand[(D[Log[IndependentEvenLetters],var]+sixpointdiff D[Log[IndependentEvenLetters],s16])/.PhaseSpacePointsCut[[i]],Modulus->$P],{-fitvec[[i]]}],{i,1,176}]//EchoTiming;
res=NullSpace[NullSpaceMat,Modulus->$P]//EchoTiming;
res[[1]]/.n_Integer:>RationalReconstruct[n]


(* ::Section:: *)
(*Fit All*)


EvenPos=RemEvenPos;


IndependentEvenLettersCut;
IndependentEvenLettersCut//Length
PhaseSpacePointsCut=PhaseSpacePoints[[;;Length[IndependentEvenLettersCut]+10]];


CloseKernels[];
LaunchKernels[10];


(*var=s234;*)
Monitor[Table[
EqnConstFull=Sum[D[G6,var] ds[var],{var,diffInvariants}];
diffauxvar=First[Solve[EqnConstFull==0,ds[s16]]];
sixpointdiff=Coefficient[diffauxvar[[1,2]],ds[var]];
FIT[var]=Table[
linpos=GetLinPos[EvenPos[[q,1]],EvenPos[[q,2]],95];
fitvec=LinearTerm[var][[;;200,1,linpos]];
PhaseSpacePointsCut=PhaseSpacePoints[[;;200]];
NullSpaceMat=ParallelTable[Join[Expand[(D[Log[IndependentEvenLettersCut],var]+sixpointdiff D[Log[IndependentEvenLettersCut],s16])/.PhaseSpacePointsCut[[i]],Modulus->$P],{-fitvec[[i]]}],{i,1,80}];
res=NullSpace[NullSpaceMat,Modulus->$P];
fit=res[[1]]/.n_Integer:>RationalReconstruct[n];
{EvenPos[[q,1]],EvenPos[[q,2]],fit},{q,1,Length[EvenPos]}
],{var,DeleteElements[diffInvariants,{s16}]}],ToString[StringForm["Running ``: ``/``",var,q,Length[EvenPos]]]]


Table[FIT[var][[-2]],{var,DeleteElements[diffInvariants,{s16}]}][[All,3,-1]]



BuildDLogFits[dvartables_]:=Module[{Rem},
If[MemberQ[dvartables[[All,3,-1]],0],Return[{{dvartables[[1,1]],dvartables[[1,2]]},$Failed}]];
Return[{{dvartables[[1,1]],dvartables[[1,2]]},Flatten[Table[
Rem=DeleteDuplicates[dvartables[[All,3,i]]];
If[Length[Rem]>1,Rem=DeleteElements[Rem,{0}];If[Length[Rem]>1,Rem=$Failed,Rem],Rem]
,{i,1,Length[dvartables[[1,3]]]}
]][[;;-2]]}]
];

BuildDLogFits[Table[FIT[var][[-1]],{var,DeleteElements[diffInvariants,{s16}]}]]


DLogFitsNew=Table[BuildDLogFits[Table[FIT[var][[evenpos]],{var,DeleteElements[diffInvariants,{s16}]}]],{evenpos,1,Length[FIT[mt2]]}];


FIT[var][[1,3]]//Length


RemEvenPos=DLogFits[[Position[DLogFits,$Failed][[All,1]],1]]


DLogFits2=Delete[DLogFits,Position[DLogFits,$Failed][[All,1]]/.n_Integer:>{n}]


Position[DLogFits,$Failed][[All,1]]


Table[Total[DLogFits[[All,2,i]]^2],{i,1,Length[IndependentEvenLetters]}]


IndependentEvenLetters


DLogFits2[[All,2,-9]]


Export[CurrentDir<>"DLogFitsEven-old.txt",DLogFits]
Export[CurrentDir<>"EvenLettersInd.txt",IndependentEvenLetters]



Export[CurrentDir<>"DLogFitsEven.txt",DLogFitsNew]
Export[CurrentDir<>"UsedEvenLetters.txt",IndependentEvenLettersCut]


(* ::Section:: *)
(*Cleaning Letters:*)


DLogFits=SortBy[Join[DLogFits2,DLogFitsNew],First]


usedLetterPos=Flatten[Position[Table[Total[DLogFitsNew[[All,2,i]]^2],{i,1,Length[IndependentEvenLettersCut]}]/.n_/;(n>0):>1,1]]


usedLetterPos//Length


IndependentEvenLettersCut=IndependentEvenLetters[[usedLetterPos]];
