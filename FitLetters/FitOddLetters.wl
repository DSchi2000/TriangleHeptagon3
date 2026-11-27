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
Get["/home/daniel/Documents/PHD-LocalFiles/TriangleHeptagon3/SetUp/RunSimplifyOnGram.wls"]



(* ::Section:: *)
(*Load Letters*)


PrevReco1=ReadList[HomeDir<>"/Reconstruction-Odd-Letters/OddLetterReco1AB.txt"];
PrevReco2=ReadList[HomeDir<>"/Reconstruction-Odd-Letters/OddLetterReco2AB.txt"];
PrevReco3=ReadList[HomeDir<>"/Reconstruction-Odd-Letters/OddLetterReco2AB-PARA.txt"];

PrevReco1=Delete[PrevReco1,DeleteDuplicates[Position[PrevReco1,$Failed][[All,1]]]/.n_Integer:>{n}];
PrevReco2=Delete[PrevReco2,DeleteDuplicates[Position[PrevReco2,$Failed][[All,1]]]/.n_Integer:>{n}];
PrevReco3=Delete[PrevReco3,DeleteDuplicates[Position[PrevReco3,$Failed][[All,1]]]/.n_Integer:>{n}];
files=FileNames["DLog*.txt", HomeDir<>"/Reconstruction-Odd-Letters/HandRecoLetters"];
PrevReco4=Map[ReadList[#]&,files];

PrevReco=SortBy[Join[PrevReco1,PrevReco2,PrevReco3,PrevReco4],First]


files=FileNames["IntegratedFunc-*.txt", HomeDir<>"/Reconstruction-Odd-Letters/HandRecoLetters"];
DLogs=Map[ReadList[#]&,files];
files=FileNames["DLog*.txt", HomeDir<>"/Reconstruction-Odd-Letters/HandRecoLetters"];
PrevReco4=Map[ReadList[#]&,files];


LaunchKernels[5]


Simplified2=ParallelMap[SimplifyOnVarietyFactored[Factor[#/.KinematicRulesForward]]&,PrevReco3[[All,2]]]
Simplified3=ParallelMap[SimplifyOnVarietyFactored[Factor[#/.KinematicRulesForward]]&,PrevReco3[[All,3]]]


Map[LeafCount[#]&,Simplified2]
Map[LeafCount[#]&,Simplified3]



BigPos=Flatten[Position[Map[LeafCount[#]&,Simplified2], x_ /; x > 1000]]
Simplified2Full=ParallelMap[SimplifyOnVarietyFactoredBest[#]&,Simplified2[[BigPos]]]


Simplified2[[BigPos]]=Simplified2Full;


PrevReco3[[All,2]]=Simplified2;
PrevReco3[[All,3]]=Simplified3;
PrevReco=SortBy[Join[PrevReco1,PrevReco2,PrevReco3,PrevReco4],First]


PrevReco3[[All,1]]


(* ::Section:: *)
(*Test one*)



pos1=58;
pos2=55;
linpos=GetLinPos[pos1,pos2,95];
posInSolve=First[Flatten[Position[PrevReco,{pos1,pos2}]]];

aFac=PrevReco[[posInSolve,2]]
bFac=PrevReco[[posInSolve,3]]
rootPart=PrevReco[[posInSolve,4]]

Test=Factor[(aFac+bFac rootPart)/(aFac-bFac rootPart)];
TestletterN=Numerator[Test];
TestletterD=Denominator[Test];

TestletterN=Collect[TestletterN,_root,Factor];
TestletterD=Collect[TestletterD,_root,Factor];

(*TestletterN=Collect[Together[aFac+bFac rootPart],_root,Factor];
TestletterD=Collect[Together[aFac-bFac rootPart],_root,Factor];*)
Testletter=TestletterN/TestletterD
(Testletter/.root[x_]:>x)//Variables


(* ::Subsection:: *)
(*dlog*)



var=mt2;
EqnConstFull=Sum[D[G6,var] ds[var],{var,diffInvariants}];
diffauxvar=First[Solve[EqnConstFull==0,ds[s16]]];
sixpointdiff=Coefficient[diffauxvar[[1,2]],ds[var]];
fitvec=LinearTerm[var][[;;200,1,linpos]];


TestVec=Expand[(D[Log[Testletter],var]+sixpointdiff D[Log[Testletter],s16])/. {Derivative[1][root][x_]:>root[x]/(2x)}/.PhaseSpacePoints[[;;10]]/.root[n_]:>PowerMod[n,1/2,$P],Modulus->$P]
Expand[fitvec[[;;10]]/TestVec,Modulus->$P]/.n_Integer:>RationalReconstruct[n]


CloseKernels[];
LaunchKernels[5];


(* ::Section:: *)
(*Linear Independence*)


(* ::Subsection:: *)
(*Build Letters*)


Monitor[OddLetters=Table[
aFac=PrevReco[[i,2]];
bFac=PrevReco[[i,3]];
rootPart=PrevReco[[i,4]];
Test=Factor[(aFac+bFac rootPart)/(aFac-bFac rootPart)];
TestletterN=Numerator[Test];
TestletterD=Denominator[Test];

TestletterN=Collect[TestletterN,_root,Factor];
TestletterD=Collect[TestletterD,_root,Factor];
TestletterN/TestletterD,{i,Range[Length[PrevReco]]}],{i}]


(* ::Subsection:: *)
(*Linear Independence*)


Length[OddLetters]


(*var=mt2;*)
Monitor[Table[
sixpointdiff=Coefficient[diffauxvar[[1,2]],ds[var]];

NullSpaceMats[var]=ParallelTable[
tpointCut=Delete[PhaseSpacePoints[[i]],Flatten[Map[Position[diffInvariants,#]&,{var,s16}]]/.n_Integer:>{n}];
tpointFull=PhaseSpacePoints[[i]];
Expand[(D[Log[OddLetters/.tpointCut],var]+sixpointdiff D[Log[OddLetters],s16])/. {Derivative[1][root][x_]:>root[x]/(2x)}/.tpointFull/.root[n_]:>PowerMod[n,1/2,$P],Modulus->$P],{i,1,370}],{var,DeleteElements[diffInvariants,{s16}]}],{var}];


(*NullSpaceMats[mt2]=3*)


Expand[NullSpaceMats[s34]-NullSpaceMats[mt2],Modulus->$P]


MatrixRank[NullSpaceMat,Modulus->$P]
rech=RowReduce[NullSpaceMats[mt2],Modulus->$P];
rech//MatrixPlot
MatrixRank[NullSpaceMat,Modulus->$P]
pivotCols = FirstPosition[#, x_ /; x != 0] & /@ rech;
pivotCols = DeleteCases[pivotCols, _Missing][[All, 1]]


ranpt=PhaseSpacePoints[[-1,All,2]]
PhaseSpacePoints[[-1]]
ranpt[[4]]=0;
ranpt


NullSpaceMat=Expand[Sum[ranpt[[i]]NullSpaceMats[diffInvariants[[i]]],{i,1,10}],Modulus->$P];
rech=RowReduce[NullSpaceMat,Modulus->$P];
rech//MatrixPlot
MatrixRank[NullSpaceMat,Modulus->$P]
pivotCols = FirstPosition[#, x_ /; x != 0] & /@ rech;
pivotCols = DeleteCases[pivotCols, _Missing][[All, 1]];

(* independent COLUMNS of the original matrix M *)
(*basisColumns = NullSpaceMatBasisMat[[All, pivotCols]]; *)  (* 191 x 91 *)
IndependentOddLetters=OddLetters[[pivotCols]]//EchoFunction[Length];


(* ::Section:: *)
(*Fit One Entry*)


Monitor[Table[
sixpointdiff=Coefficient[diffauxvar[[1,2]],ds[var]];

NullSpaceMats[var]=ParallelTable[
tpointCut=Delete[PhaseSpacePoints[[i]],Flatten[Map[Position[diffInvariants,#]&,{var,s16}]]/.n_Integer:>{n}];
tpointFull=PhaseSpacePoints[[i]];
Expand[(D[Log[IndependentOddLetters/.tpointCut],var]+sixpointdiff D[Log[IndependentOddLetters],s16])/. {Derivative[1][root][x_]:>root[x]/(2x)}/.tpointFull/.root[n_]:>PowerMod[n,1/2,$P],Modulus->$P],{i,1,76}],{var,DeleteElements[diffInvariants,{s16}]}],{var}];


IndependentOddLetters[[1]]



Complement[OddPos,PrevReco[[All,1]]]


PrevReco[[All,1]]


var=mt2;
Mats=Table[
TestMat=NullSpaceMats[var];
pos1=95;
pos2=88;
linpos=GetLinPos[pos1,pos2,95];
fitvec=LinearTerm[var][[;;76,1,linpos]];
MapThread[Append, {TestMat, fitvec}],{var,DeleteElements[diffInvariants,{s16}]}];
nspace=NullSpace[Join@@Mats,Modulus->$P][[1]]/.n_Integer:>RationalReconstruct[n]



NullSpace[Join[Mat1,Mat2],Modulus->$P][[1]]/.n_Integer:>RationalReconstruct[n]


Position[nspace,1]
IndependentOddLetters[[10]]


(* ::Section:: *)
(*Fit All*)


DLogFits=Monitor[Table[
Mats=Table[
TestMat=NullSpaceMats[var];
pos1=OddPos[[q,1]];
pos2=OddPos[[q,2]];
linpos=GetLinPos[pos1,pos2,95];
fitvec=LinearTerm[var][[;;76,1,linpos]];
MapThread[Append, {TestMat, -fitvec}],{var,DeleteElements[diffInvariants,{s16}]}];
nspaces=NullSpace[Join@@Mats,Modulus->$P];
If[Length[nspaces]>0,nspace=nspaces[[1]]/.n_Integer:>RationalReconstruct[n];,nspace=$Failed,
If[nspaces[[1,-1]]===0,nspace=$Failed];
];
{{OddPos[[q,1]],OddPos[[q,2]]},nspace},{q,1,Length[OddPos]}],ToString[StringForm["Running: ``/``",q,Length[OddPos]]]]



(*var=s12;
Monitor[
Table[
FIT[var]=Table[
TestMat=NullSpaceMats[var];
pos1=OddPos[[i,1]];
pos2=OddPos[[i,2]];
linpos=GetLinPos[pos1,pos2,95];
fitvec=LinearTerm[var][[;;76,1,linpos]];
Mat=MapThread[Append, {TestMat, fitvec}];
fit=NullSpace[Mat,Modulus->$P][[1]]/.n_Integer:>RationalReconstruct[n];
{OddPos[[i,1]],OddPos[[i,2]],fit},{i,1,Length[OddPos]}],{var,DeleteElements[diffInvariants,{s16}]}],ToString[StringForm["Running ``: ``/``",var,i,Length[OddPos]]]]*)



Position[DLogFits,$Failed]
Fails=DLogFitsNew[[Position[DLogFits,$Failed][[All,1]],1]]
Complement[Fails,fails1]


Position[OddPos,{58,55}]
TestEntFit=DLogFits[[123]]





var=s56;
(*TestEnt=(FIT[var][[123]][[3]][[;;-2]]).Log[IndependentOddLetters];*)
TestEnt=TestEntFit[[2]][[;;-2]] . Log[IndependentOddLetters];
TestEntFit
EqnConstFull=Sum[D[G6,var] ds[var],{var,diffInvariants}];
diffauxvar=First[Solve[EqnConstFull==0,ds[s16]]];
sixpointdiff=Coefficient[diffauxvar[[1,2]],ds[var]];
fitvec=LinearTerm[var][[;;200,1,GetLinPos[58,55,95]]];


TestVec=Expand[(D[TestEnt,var]+sixpointdiff D[TestEnt,s16])/. {Derivative[1][root][x_]:>root[x]/(2x)}/.PhaseSpacePoints[[;;10]]/.root[n_]:>PowerMod[n,1/2,$P],Modulus->$P]
Expand[fitvec[[;;10]]/TestVec,Modulus->$P]/.n_Integer:>RationalReconstruct[n]


(* ::Subsection:: *)
(*Export*)


Export[CurrentDir<>"DLogFits-OddLetters.txt",DLogFits]
Export[CurrentDir<>"OddLettersInd.txt",IndependentOddLetters]


(* ::Section:: *)
(*Cleaning Letters:*)


DLogFitsS=Delete[DLogFits,Position[DLogFits,$Failed][[All,1]]/.n_Integer:>{n}]


DLogFitsS[[All,2]]


usedLetterPos=Flatten[Position[Table[Total[DLogFitsS[[All,2,i]]^2],{i,1,Length[IndependentOddLetters]}]/.n_/;(n>0):>1,1]]
Complement[Range[Length[IndependentOddLetters]],usedLetterPos]


(* ::Section:: *)
(*Check remaining Problems:*)


Fails


OddPos//Length





unknowns=Fails
Print["Knowns:"]
knowns=Complement[OddPos,Fails]
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
TODOS=Complement[basis,DLogFitsS[[All,1]]]//EchoFunction[Length]


TODOS=Complement[basis,DLogFitsS[[All,1]]]//EchoFunction[Length]
Export[CurrentDir<>"TODOODD.txt",TODOS];


TODOS[[1]]


Position[DLogFitsS,TODOS[[6]]]


DLogFits[[1]]


Intersection[Fails,PrevReco[[All,1]]]


Fails
