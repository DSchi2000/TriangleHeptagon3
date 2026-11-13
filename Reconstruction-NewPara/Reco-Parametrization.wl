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
LinearTerm[s12]//Length
PhaseSpacePoints//Length


(* ::Section:: *)
(*Load Extras*)


avec=ReadList[HomeDir<>"/Reconstruction1/Samples/AVECS.txt"][[1]];
bvec=ReadList[HomeDir<>"/Reconstruction1/Samples/AVECS.txt"][[2]];
SliceTpoint=Thread[Rule[avec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]]


G6=Get[HomeDir<>"/SetUp/SixPointGram.txt"];
KnowFactor=Get[CurrentDir<>"KnowLogFac-93-93.txt"];
Get[HomeDir<>"/4DBasis/Parametrization.m"]
KnowFactor=Factor[KnowFactor/.KinematicRulesBackWards]


(* ::Section:: *)
(*Prepare Reco*)


var=s12
EqnConstFull=Sum[D[G6,var] ds[var],{var,diffInvariants}];
diffauxvar=First[Solve[EqnConstFull==0,ds[s16]]];
sixpointdiff=Coefficient[diffauxvar[[1,2]],ds[var]];
pos=GetLinPos[93,93,95];

LinearTerm[s16]=SparseArray[ConstantArray[0,{Length[PhaseSpacePoints],1,Length[PureIntegrals]^2}]];
FitVals=Join[Table[LinearTerm[var][[All,1,pos]],{var,diffInvariants}]];
FitVals=Normal[Transpose[FitVals]];
phasespacepointsPara=Map[ConvertMandsToPoint[#]&,PhaseSpacePoints];
(*phasespacepointsPara=Table[Join[phasespacepointsPara[[i]],tvarlist[[i]]],{i,1,Length[phasespacepointsPara]}];*)

FitValsPara=Table[ConvertDF[FitVals[[i]],phasespacepointsPara[[i]]],{i,1,Length[FitVals]}];

logRational=1;
DlogDiff=Map[Expand[logRational  D[ Log[KnowFactor/.Delete[#,Position[paravars,var]]],var]/.#,Modulus->$P]&,phasespacepointsPara];//EchoTiming
(*DlogDiffs16=Map[Expand[sixpointdiff*logRational  D[Log[KnowFactor/.Delete[#,Position[diffInvariants,s16]]],s16]/.#,Modulus->$P]&,PhaseSpacePoints];//EchoTiming*)

fitvec=FitValsPara[[All,First[Flatten[Position[paravars,var]]]]];
fitvec=Expand[fitvec-((DlogDiff)),Modulus->$P];
Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
(*Integrand=Expand[UniVarReconstructFUNCS[fitvec,phasespacepointsPara,tvar,20,20,{},1,1],Modulus->$P];*)
(*Integrand=Expand[BivariateReconstruction[fitvec,phasespacepointsPara,{xvar,yvar},6,6],Modulus->$P]/.{xvar->x,yvar->y};
Expand[Integrand[[1]],Modulus->$P]
Factor[Integrand[[1]],Modulus->$P]*)




(* ::Subsection:: *)
(*Uvar*)


Integrand=Expand[UniVarReconstructFUNCS[fitvec,PhaseSpacePoints,tvar,30,30,{},1,1],Modulus->$P];


Integrand[[1]]


paravars


(* ::Section:: *)
(*Reconstruction*)


Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
(*res=ReconstructInhomogineus[fitvec,phasespacepointsPara,{paravars[[7;;]],1,2},{paravars[[;;6]],4,4},Debug->True]*)
res=ReconstructInhomogineus[fitvec,phasespacepointsPara,{paravars[[7;;]],2,2},{paravars[[;;6]],3,4},Debug->True]
(*Reco=ReconstructAlg[fitvec,PhaseSpacePoints[[All,All,2]],diffInvariants,5,6,{},1,1,Debug->True];*)


ReconstructedFunc=FuncR=BuildFunction[res[[1,1]],res[[2]]];
ReconstructedFunc=Factor[ReconstructedFunc/.n_Integer:>RationalReconstruct[n]]


Integrate[ReconstructedFunc,c1]


ReconstructedFunc=FuncR=BuildFunction[res[[1,1]],res[[2]]];
ReconstructedFunc=Factor[ReconstructedFunc/.n_Integer:>RationalReconstruct[n]]
Integrate[ReconstructedFunc,s12]


TestA=%85
TestB=%105


TestA/.Log[x_]:>x
TestB/.Log[x_]:>x


Factor[(TestA/.Log[x_]:>x)/(TestB/.Log[x_]:>x)]


Export[CurrentDir<>"Dlog-93-93.txt",Integrate[ReconstructedFunc,s12]]
