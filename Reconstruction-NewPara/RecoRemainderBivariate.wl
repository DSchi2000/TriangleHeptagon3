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
OutPutPath="Reconstruction-NewPara/Samples/Bivariate-Slice-1";

phasespacepoints=ReadList[FileNameJoin[{HomeDir,ToString[StringForm["``/phasespacepoints.txt",OutPutPath]]}]];
tvarlist=ReadList[FileNameJoin[{HomeDir,ToString[StringForm["``/tvars.txt",OutPutPath]]}]][[;;;; 1]];
phasespacepoints=phasespacepoints[[All,2;;]][[1;; ;; 1]];
diffInvariants={mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56};
PhaseSpacePoints=Map[Thread[Rule[diffInvariants,#]]&,phasespacepoints];
PhaseSpacePoints=Table[Join[PhaseSpacePoints[[i]],tvarlist[[i]]],{i,1,Length[PhaseSpacePoints]}];
(*PhaseSpacePoints=Map[Thread[Rule[diffInvariants,#]]&,phasespacepoints];*)
rootSet=(Cases[PureIntegrals[[-5;;]],_root,Infinity]//DeleteDuplicates)/.root[x_]:>x;
GetLinPos[x1_,x2_,NMI_]:=First[Flatten[Position[Flatten[Table[{i,j},{i,1,NMI},{j,1,NMI}],1],{x1,x2}]]];
Get[ToString[StringJoin[HomeDir,"/files/ExtraModules.m"]]];
Get[ToString[StringJoin[HomeDir,"/files/ExtraModules.m"]]];
usedvars={mt2,s12,s123,s23,s234,s34,s345,s45,s56};
For[i=1,i<=Length[usedvars],i++,
LinearTerm[usedvars[[i]]]=ReadList[FileNameJoin[{HomeDir,ToString[StringForm["``/LinTerm``_0001.txt",OutPutPath,usedvars[[i]]]]}]];
];



PhaseSpacePoints[[1]]


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



var=c1
EqnConstFull=Sum[D[G6,var] ds[var],{var,diffInvariants}];
diffauxvar=First[Solve[EqnConstFull==0,ds[s16]]];
sixpointdiff=Coefficient[diffauxvar[[1,2]],ds[var]];
pos=GetLinPos[93,93,95];

LinearTerm[s16]=SparseArray[ConstantArray[0,{Length[PhaseSpacePoints],1,Length[PureIntegrals]^2}]];
FitVals=Join[Table[LinearTerm[var][[All,1,pos]],{var,diffInvariants}]];
FitVals=Normal[Transpose[FitVals]];
phasespacepointsPara=Map[ConvertMandsToPoint[#]&,PhaseSpacePoints];
phasespacepointsPara=Table[Join[phasespacepointsPara[[i]],tvarlist[[i]]],{i,1,Length[phasespacepointsPara]}];

FitValsPara=Table[ConvertDF[FitVals[[i]],phasespacepointsPara[[i]]],{i,1,Length[FitVals]}];




logRational=1;
DlogDiff=Map[Expand[logRational  D[ Log[KnowFactor/.Delete[#,Position[paravars,var]]],var]/.#,Modulus->$P]&,phasespacepointsPara];//EchoTiming
(*DlogDiffs16=Map[Expand[sixpointdiff*logRational  D[Log[KnowFactor/.Delete[#,Position[diffInvariants,s16]]],s16]/.#,Modulus->$P]&,PhaseSpacePoints];//EchoTiming*)

fitvec=FitValsPara[[All,First[Flatten[Position[paravars,var]]]]];
fitvec=Expand[fitvec-((DlogDiff)),Modulus->$P];
Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
(*Integrand=Expand[UniVarReconstructFUNCS[fitvec,phasespacepointsPara,tvar,20,20,{},1,1],Modulus->$P];*)
Integrand=Expand[BivariateReconstruction[fitvec,phasespacepointsPara,{xvar,yvar},6,6],Modulus->$P]/.{xvar->x,yvar->y};
Expand[Integrand[[1]],Modulus->$P]
Factor[Integrand[[1]],Modulus->$P]






(* ::Section:: *)
(*Build NSamples*)


GetNumberOfTermsNeeded[n_,num_]:=Binomial[n + num-1, n-1]
GetNumberOfTermsNeededMixed[nvar1_,nvar2_,n1_,n2_,d1_,d2_]:=GetNumberOfTermsNeeded[nvar1,n1]*(1+n2)^nvar2+GetNumberOfTermsNeeded[nvar1,d1]*(1+d2)^nvar2;


GetNumberOfTermsNeededMixed[6,3,3,2,4,2]
GetNumberOfTermsNeededMixed[6,3,4,1,4,2]



Get[CurrentDir<>"AnsatzBuilder.wl"]


ansatz=BuildRationalAnsatz[paravars[[;;6]],paravars[[7;;]],NumDegree->3,DenDegree->4,CMaxPower->2,Normalize->True,NumCoeff->a,DenCoeff->b]


ansatz[["Ansatz"]]//Numerator//Length
ansatz[["Ansatz"]]//Denominator//Length


1512+3402
