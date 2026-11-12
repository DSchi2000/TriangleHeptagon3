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
(*OutPutPath="Reconstruction1/Samples/Univariate Slice";*)
OutPutPath="Reconstruction-NewPara/Samples/Univariate-Slice-1";

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



(*avec=ReadList[CurrentDir<>"AVECS.txt"][[1]];
bvec=ReadList[CurrentDir<>"AVECS.txt"][[2]];
SliceTpoint=Thread[Rule[avec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]]*)


avec=ReadList[CurrentDir<>"AVECS.txt"][[1]];
bvec=ReadList[CurrentDir<>"AVECS.txt"][[2]];
SliceTpoint=Thread[Rule[avec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]]


G6=Get[HomeDir<>"/SetUp/SixPointGram.txt"];
KnowFactor=Get[CurrentDir<>"KnowLogFac-93-93.txt"];
LinearizedEvenPos=Map[GetLinPos[#[[1]],#[[2]],95]&,EvenPos];
Get[HomeDir<>"/4DBasis/Parametrization.m"]
KnowFactor=Factor[KnowFactor/.KinematicRulesBackWards]



var=c2
EqnConstFull=Sum[D[G6,var] ds[var],{var,diffInvariants}];
diffauxvar=First[Solve[EqnConstFull==0,ds[s16]]];
sixpointdiff=Coefficient[diffauxvar[[1,2]],ds[var]];
pos=GetLinPos[93,93,95];

LinearTerm[s16]=SparseArray[ConstantArray[0,{Length[PhaseSpacePoints],1,Length[PureIntegrals]^2}]];
FitVals=Join[Table[LinearTerm[var][[All,1,pos]],{var,diffInvariants}]];
FitVals=Normal[Transpose[FitVals]];
phasespacepointsPara=Map[ConvertMandsToPoint[#]&,PhaseSpacePoints];
phasespacepointsPara=Table[Join[phasespacepointsPara[[i]],{tvarlist[[i]]}],{i,1,Length[phasespacepointsPara]}];

FitValsPara=Table[ConvertDF[FitVals[[i]],phasespacepointsPara[[i]]],{i,1,Length[FitVals]}];




DlogDiff
Expand[D[Log[KnowFactor],c2]/.phasespacepointsPara[[1]],Modulus->$P]


logRational=1;
DlogDiff=Map[Expand[logRational  D[ Log[KnowFactor/.Delete[#,Position[paravars,var]]],var]/.#,Modulus->$P]&,phasespacepointsPara];//EchoTiming
(*DlogDiffs16=Map[Expand[sixpointdiff*logRational  D[Log[KnowFactor/.Delete[#,Position[diffInvariants,s16]]],s16]/.#,Modulus->$P]&,PhaseSpacePoints];//EchoTiming*)

fitvec=FitValsPara[[All,First[Flatten[Position[paravars,var]]]]];
fitvec=Expand[fitvec-((DlogDiff)),Modulus->$P];
Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
Integrand=Expand[UniVarReconstructFUNCS[fitvec,phasespacepointsPara,tvar,20,20,{},1,1],Modulus->$P];
Expand[Integrand[[1]],Modulus->$P]
Factor[Integrand[[1]],Modulus->$P]






fitvec=FitVals[[All,First[Flatten[Position[diffInvariants,var]]]]];
(*fitvec=Expand[fitvec-((DlogDiff+DlogDiffs16)),Modulus->$P];*)
Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
Integrand=Expand[UniVarReconstructFUNCS[fitvec,PhaseSpacePoints,tvar,20,20,{},1,1],Modulus->$P];
Expand[Integrand[[1]],Modulus->$P]
int=Factor[Integrand[[1]],Modulus->$P]





DlogDiff/.PhaseSpacePoints


LogDiffN=Expand[((DlogDiff+DlogDiffs16)/.PhaseSpacePoints[[1]]),Modulus->$P]//EchoTiming


Factor[KnowFactor/.SliceTpoint,Modulus->$P]


int


Denominator[KnowFactor][[3]]
Factor[Denominator[KnowFactor][[1]]/.SliceTpoint,Modulus->$P]


dlogList=ReadList[HomeDir<>"/Reconstruction1/dlogs.txt"];


TargetDlog=dlogList[[-6]]/.Log[x_]:>Factor[x,Modulus->$P]


test=Factor[D[Log[TargetDlog],t],Modulus->$P]

testDiff=Factor[D[Log[Factor[KnowFactor/.SliceTpoint,Modulus->$P]],t],Modulus->$P]



Factor[KnowFactor/.SliceTpoint,Modulus->$P]


GetNumberOfTermsNeeded[n_,num_,den_]:=Binomial[n + num-1, n-1]+Binomial[n + den-1, n-1]


GetNumberOfTermsNeeded[10,5,6]


Test=Apart[(1546264606 (974720773+1134661556 t+t^2) (313609644+759319654 t+1666975794 t^2+t^3))/((177873137+t) (249867599+t) (814611916+t) (869940708+t) (1683896177+t) (2062316089+t)),t,Modulus->$P]


Sum[Factor[Test[[i]],Modulus->$P],{i,1,Length[Test]}]


Factor[D[Log[TargetDlog/Factor[KnowFactor/.SliceTpoint,Modulus->$P]],t],Modulus->$P]


GetNumberOfTermsNeeded[10,6,7]
GetNumberOfTermsNeeded[10,7,8]


fitvec/


Expand[PhaseSpacePoints[[All,2,2]]*fitvec,Modulus->$P]


Expand[SliceTpoint[[2,2]]/.{t->PhaseSpacePoints[[1,11,2]]},Modulus->$P]


PhaseSpacePoints[[All,2,2]]





Expand[SliceTpoint/.{t->PhaseSpacePoints[[2,11,2]]},Modulus->$P]
PhaseSpacePoints[[2]]



(* ::Section:: *)
(*Short Test For claude*)


(*Test=(mt2+s16)/.KinematicRulesBackWards*)
Test2=Factor[KnowFactor/.KinematicRulesForward]


G6


Export[CurrentDir<>"TestExpressionForClaude.txt",Test2]


Test2





sol=Solve[G6==0,s12];
Factor[Test2/.sol[[1]]]



