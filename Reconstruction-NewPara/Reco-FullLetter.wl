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


(* ::Section::Closed:: *)
(*Uvar Samples*)


OutPutPath="Reconstruction1/Samples/Univariate Slice";
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


(* ::Section:: *)
(*Load Samples*)


OutPutPaths={"Reconstruction1/Samples/Full-Samples-10000"};
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
LinearTerm[s12]//Length
PhaseSpacePoints//Length


(* ::Section:: *)
(*Load Extras*)


avec=ReadList[HomeDir<>"/Reconstruction1/Samples/AVECS.txt"][[1]];
bvec=ReadList[HomeDir<>"/Reconstruction1/Samples/AVECS.txt"][[2]];
SliceTpoint=Thread[Rule[avec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]]


G6=Get[HomeDir<>"/SetUp/SixPointGram.txt"];
KnowFactor=Get[CurrentDir<>"KnowLogFac-93-90.txt"];


(* ::Section::Closed:: *)
(*Dbug Extra*)


TestDlog=Get["/home/dschiebe/Documents/PhD-Files/TriangleHeptagon3/Reconstruction1/FinishedDlog-94-94.txt"]
KnowFactor=1/(s34 (mt2 s123^2 s34+mt2^2 s12 s56-mt2^2 s123 s56+mt2 s123 s345 s56-2 mt2 s12 s45 s56+mt2 s123 s45 s56-s123 s345 s45 s56+s12 s45^2 s56) (mt2^2 s12^2 s234-2 mt2^2 s12 s123 s234+mt2^2 s123^2 s234+mt2^2 s12 s123 s34-mt2^2 s123^2 s34-mt2 s12 s123 s16 s34+mt2 s123^2 s16 s34+2 mt2 s12 s123 s23 s34-mt2^2 s12 s234 s34+mt2^2 s123 s234 s34+mt2 s12 s123 s234 s34-mt2 s123^2 s234 s34-mt2^2 s123 s34^2+mt2 s123^2 s34^2+mt2 s123 s16 s34^2+s123^2 s16 s34^2-2 mt2 s123 s23 s34^2+2 mt2 s12 s123 s234 s345-2 mt2 s123^2 s234 s345+mt2 s123^2 s34 s345-s123^2 s16 s34 s345-mt2 s123 s234 s34 s345-s123^2 s234 s34 s345+s123^2 s234 s345^2-2 mt2 s12^2 s234 s45+2 mt2 s12 s123 s234 s45-mt2 s12 s123 s34 s45+s12 s123 s16 s34 s45+2 mt2 s12 s234 s34 s45-mt2 s123 s234 s34 s45+s12 s123 s234 s34 s45+mt2 s123 s34^2 s45-s123 s16 s34^2 s45-2 s12 s123 s234 s345 s45+s123 s234 s34 s345 s45+s12^2 s234 s45^2-s12 s234 s34 s45^2+mt2^2 s12 s23 s56-mt2^2 s123 s23 s56+2 mt2 s12 s16 s34 s56-2 mt2 s123 s16 s34 s56+mt2^2 s23 s34 s56-2 mt2 s12 s23 s34 s56+mt2 s123 s23 s34 s56-mt2 s12 s23 s345 s56+2 mt2 s123 s23 s345 s56-mt2 s12 s234 s345 s56+mt2 s123 s234 s345 s56+mt2 s123 s34 s345 s56+s123 s16 s34 s345 s56-mt2 s23 s34 s345 s56+s123 s23 s34 s345 s56-s123 s23 s345^2 s56-s123 s234 s345^2 s56-mt2 s12 s23 s45 s56+mt2 s12 s234 s45 s56-mt2 s123 s234 s45 s56-2 mt2 s12 s34 s45 s56+mt2 s123 s34 s45 s56-2 s12 s16 s34 s45 s56+s123 s16 s34 s45 s56-mt2 s23 s34 s45 s56-2 s12 s23 s34 s45 s56+s12 s23 s345 s45 s56+s12 s234 s345 s45 s56+s123 s234 s345 s45 s56-2 s123 s34 s345 s45 s56+s23 s34 s345 s45 s56-s12 s234 s45^2 s56+2 s12 s34 s45^2 s56-mt2 s23 s345 s56^2+s23 s345^2 s56^2+mt2 s23 s45 s56^2-s23 s345 s45 s56^2));


TestDlog


(* ::Section:: *)
(*Prepare Reco*)


var=mt2
EqnConstFull=Sum[D[G6,var] ds[var],{var,diffInvariants}];
diffauxvar=First[Solve[EqnConstFull==0,ds[s16]]];
sixpointdiff=Coefficient[diffauxvar[[1,2]],ds[var]];
pos=GetLinPos[93,90,95];

LinearTerm[s16]=SparseArray[ConstantArray[0,{Length[PhaseSpacePoints],1,Length[PureIntegrals]^2}]];
FitVals=Join[Table[LinearTerm[var][[All,1,pos]],{var,diffInvariants}]];
FitVals=Normal[Transpose[FitVals]];

logRational=1
DlogDiff=Map[Expand[logRational  D[ Log[KnowFactor/.Delete[#,Position[diffInvariants,var]]],var]/.#,Modulus->$P]&,PhaseSpacePoints];//EchoTiming
DlogDiffs16=Map[Expand[sixpointdiff*logRational  D[Log[KnowFactor/.Delete[#,Position[diffInvariants,s16]]],s16]/.#,Modulus->$P]&,PhaseSpacePoints];//EchoTiming

fitvec=FitVals[[All,First[Flatten[Position[diffInvariants,var]]]]];
fitvec=Expand[fitvec-((DlogDiff+DlogDiffs16)),Modulus->$P];



(* ::Subsection::Closed:: *)
(*Uvar*)


Integrand=Expand[UniVarReconstructFUNCS[fitvec,PhaseSpacePoints,tvar,30,30,{},1,1],Modulus->$P];


Integrand[[1]]


(* ::Section:: *)
(*Reconstruction*)


Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
Reco=ReconstructOnVariety[fitvec,PhaseSpacePoints[[All,All,2]],diffInvariants,5,6,{G6,s16},Debug->True];
(*Reco=ReconstructAlg[fitvec,PhaseSpacePoints[[All,All,2]],diffInvariants,5,6,{},1,1,Debug->True];*)


Ansatz=Reco[[2]];
ReconstructedFunc=Together[Ansatz/.Thread[Rule[Sort[Cases[Ansatz,_fac,Infinity]],Reco[[1,1]]]]/.n_Integer:>RationalReconstruct[n]]


?ReconstructAlg
