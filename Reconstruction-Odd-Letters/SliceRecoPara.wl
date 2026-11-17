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
OutPutPath="Reconstruction-NewPara/Samples/Univariate-Slice-2";
(*OutPutPath="Reconstruction1/Samples/Univariate Slice";*)
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

Get[HomeDir<>"/4DBasis/Parametrization.m"]
phasespacepointsPara=Map[ConvertMandsToPoint[#]&,PhaseSpacePoints];
phasespacepointsPara=Table[Join[phasespacepointsPara[[i]],{tvarlist[[i]]}],{i,1,Length[phasespacepointsPara]}];







LinearTerm[usedvars[[1]]]//Length


(* ::Section:: *)
(*Find Even Letter Pos*)


rootVector=Map[Times@@DeleteDuplicates[Cases[#,_root,Infinity]]&,PureIntegrals]/.{}->1/.root[x_]:>x;
rootVector=Factor[rootVector/.KinematicRulesBackWards];


(* ::Section:: *)
(*Reconstruction*)


EQN=Flatten[Table[{avector[i]+phasespacepointsPara[[1,-1,2]] bvector[i] ==phasespacepointsPara[[1,i,2]],avector[i]+phasespacepointsPara[[2,-1,2]] bvector[i] ==phasespacepointsPara[[2,i,2]]},{i,1,Length[phasespacepointsPara[[1,All,1]]]-1}]];
solset=Flatten[Table[{avector[i],bvector[i]},{i,1,Length[phasespacepointsPara[[1,All,1]]]-1}]];
sol=Solve[EQN,solset,Modulus->$P][[1]];
avec=Thread[Rule[phasespacepointsPara[[1,;;-2,1]],Table[avector[i],{i,1,Length[phasespacepointsPara[[1,All,1]]]-1}] ]]/.sol
bvec=Thread[Rule[phasespacepointsPara[[1,;;-2,1]],Table[bvector[i],{i,1,Length[phasespacepointsPara[[1,All,1]]]-1}] ]]/.sol
LinearTerm[s16]=SparseArray[ConstantArray[0,{Length[phasespacepointsPara],1,Length[PureIntegrals]^2}]];
SliceTpoint=Thread[Rule[bvec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]]
SliceTpoint=Thread[Rule[bvec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]];
Expand[SliceTpoint/.t->phasespacepointsPara[[-1,-1,2]],Modulus->$P]
phasespacepointsPara[[-1]]
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


PrevReco=ReadList[CurrentDir<>"OddLetterReco1AB.txt"];
OddPos=PrevReco[[DeleteDuplicates[Position[PrevReco,$Failed][[All,1]]]]][[All,1]]


(* ::Subsection:: *)
(*Reco*)


Transpose[FitVals]


Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
pos1=20;
pos2=2;
sqrt1=rootVector[[pos1]]
sqrt2=rootVector[[pos2]]
pos=GetLinPos[pos1,pos2,95];
sqrt1N=PowerMod[sqrt1/.phasespacepointsPara,1/2,$P];
sqrt2N=PowerMod[sqrt2/.phasespacepointsPara,1/2,$P];

FitVals=Transpose[Join[Table[LinearTerm[var][[All,1,pos]],{var,diffInvariants}]]];
FitValsPara=Table[ConvertDF[FitVals[[i]],phasespacepointsPara[[i]]],{i,1,Length[FitVals]}];

Integrands=Table[Integrand=Expand[UniVarReconstructFUNCS[Expand[FitValsPara[[All,varpos]]/(sqrt2N/sqrt1N),Modulus->$P],phasespacepointsPara,tvar,20,20,{},1,1],Modulus->$P];Factor[Integrand[[1]],Modulus->$P],{varpos,1,9}]
NewIntDt=Factor[Expand[Integrands . bvec[[All,2]],Modulus->$P],Modulus->$P]
NewIntDt=NewIntDt*(root[sqrt2]/root[sqrt1])/.root[1]->1/.SliceTpoint/.root[x_]:>root[Factor[x,Modulus->$P]]


PhaseSpacePoints//Length


(* ::Section:: *)
(*Integrate*)


(*Get["/home/dschiebe/Documents/packages/trager/Trager.wl"]*)
Get["/home/daniel/Documents/packages/trager/Trager.wl"]
roots=DeleteDuplicates[Cases[NewIntDt,_root,Infinity]]
rootrepls=Thread[Rule[roots,Map[rt[#]&,Range[Length[roots]]]]]
qring=Expand[Table[rt[i]^2-(roots[[i]]/.root[x_]:>x),{i,Range[Length[roots]]}],Modulus->$P]
rtvars=Table[rt[i],{i,Range[Length[roots]]}]
invrootrepls = Reverse /@ rootrepls;
NewIntDt/.rootrepls
Dlog=TragerLogSubtractionNumeric[NewIntDt/.rootrepls,t,rtvars,qring,Modulus->$P]
Dlog=Dlog/.invrootrepls;


LogArg=Dlog/.Log[x_]:>x;
TestN=Collect[Numerator[LogArg],_root,Factor,Modulus->$P]
TestD=Collect[Denominator[LogArg],_root,Factor,Modulus->$P]
OddLetterFactored=Factor[Expand[TestN*TestD,Modulus->$P]/.root[x_]^2:>x,Modulus->$P]


EvenLetterCandidates=ReadList[HomeDir<>"/Reconstruction-NewPara/EvenLetterCandidates.txt"];
ELSclice=Factor[EvenLetterCandidates/.SliceTpoint,Modulus->$P];
ELSclice//Length
EvenLetterCandidates//Length


rootsSlice=Times@@DeleteDuplicates[Cases[NewIntDt,_root,Infinity]]
Test=Factor[Coefficient[TestN,rootsSlice],Modulus->$P]
Test2=Factor[TestN/.root[x_]:>0,Modulus->$P]

(*FactorsOdd=Times@@Map[CheckLetterFactor[Test,#]&,Range[191]]*)
FactorsOdd=Times@@Map[CheckLetterFactor[Test,#]&,Range[191]]
FactorsOdd=Times@@Map[CheckLetterFactor[Test2,#]&,Range[191]]


FactorsOdd=Times@@Map[CheckLetterFactor[OddLetterFactored,#]&,Range[191]]
nfac=Expand[OddLetterFactored/FactorsOdd/.Fac[n_]:>ELSclice[[n]],Modulus->$P]//RationalReconstruct


FactorsOdd/.Fac[n_]:>EvenLetterCandidates[[n]]


EvenLetterCandidates[[2]]
EvenLetterCandidates[[10]]
r=sqrt1*sqrt2


b=Fac[6]/.Fac[n_]:>EvenLetterCandidates[[n]]
Factor[FactorsOdd/.Fac[n_]:>EvenLetterCandidates[[n]]]

sol=Solve[Discriminant[k Factor[FactorsOdd/.Fac[n_]:>EvenLetterCandidates[[n]]]+b^2 *r,s12]==0,k];
sol=Flatten[Select[sol, FreeQ[k /. #, _Symbol] &]]

a=root[Factor[k Factor[FactorsOdd/.Fac[n_]:>EvenLetterCandidates[[n]]]+b^2 *r/.sol]]
a=a/.root[x_^2]:>x
LetterTest=(a+b root[r])/(a-b root[r])
Expand[Factor[D[3 Log[LetterTest],s12]/.Derivative[1][root][x_]:>root[x]/(2x)]/.PhaseSpacePoints[[1]]/.root[x_]:>PowerMod[x,1/2,$P],Modulus->$P]
FitVals[[2,1]]


alpha=Coefficient[sqrt1,s12^2]
beta=Coefficient[sqrt1,s12]
gamma=sqrt1/.s12->0
p=b Sqrt[alpha]
q=Factor[b beta/(2p)]
aTerm=p s12+q

OddLetterAB=aTerm+b Sqrt[sqrt1];
OddLetterABM=aTerm-b Sqrt[sqrt1];
RecoFac=Factor[Expand[OddLetterAB *OddLetterABM]]
Factor[RecoFac/.SliceTpoint,Modulus->$P]


=






(* ::Section:: *)
(*SATISFACTORY YEAAAAH*)


(* ::Subsection:: *)
(*Reconstruct*)


Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
LaunchKernels[5];
Monitor[IntegrandsMat=ParallelTable[
pos1=matpos[[1]];
pos2=matpos[[2]];
sqrt1=rootVector[[pos1]];
sqrt2=rootVector[[pos2]];
pos=GetLinPos[pos1,pos2,95];
sqrt1N=PowerMod[sqrt1/.PhaseSpacePoints,1/2,$P];
sqrt2N=PowerMod[sqrt2/.PhaseSpacePoints,1/2,$P];
FitVals=Join[Table[LinearTerm[var][[All,1,pos]],{var,diffInvariants}]];
Integrands=Table[Integrand=Expand[UniVarReconstructFUNCS[Expand[FitVals[[varpos]]/(sqrt2N/sqrt1N),Modulus->$P],PhaseSpacePoints,tvar,20,20,{},1,1],Modulus->$P];Factor[Integrand[[1]],Modulus->$P],{varpos,1,10}];
Integrands=Integrands*(root[sqrt2]/root[sqrt1])/.root[1]->1/.SliceTpoint/.root[x_]:>root[Factor[x,Modulus->$P]];
{matpos,Integrands},{matpos,OddPos}],matpos];


(*NewIntDt=Factor[Expand[Integrands . bvec[[All,2]],Modulus->$P],Modulus->$P];
NewIntDt=NewIntDt*(root[sqrt2]/root[sqrt1])/.root[1]->1/.SliceTpoint/.root[x_]:>root[Factor[x,Modulus->$P]];*)


Export[CurrentDir<>"Integrands.txt",IntegrandsMat]


(* ::Subsection:: *)
(*Integrate*)


(*Integrands=IntegrandsMat[[-5,2]];
NewIntDt=Factor[Expand[Integrands . bvec[[All,2]],Modulus->$P],Modulus->$P]*)

Get["/home/daniel/Documents/packages/trager/Trager.wl"];
Monitor[IntegrandsMat=Table[
Integrands=IntegrandsMat[[pos,2]];
NewIntDt=Factor[Expand[Integrands . bvec[[All,2]],Modulus->$P],Modulus->$P];
roots=DeleteDuplicates[Cases[NewIntDt,_root,Infinity]];
rootrepls=Thread[Rule[roots,Map[rt[#]&,Range[Length[roots]]]]];
qring=Expand[Table[rt[i]^2-(roots[[i]]/.root[x_]:>x),{i,Range[Length[roots]]}],Modulus->$P];
rtvars=Table[rt[i],{i,Range[Length[roots]]}];
invrootrepls = Reverse /@ rootrepls;
NewIntDt/.rootrepls;
Dlog=TragerLogSubtractionNumeric[NewIntDt/.rootrepls,t,rtvars,qring,Modulus->$P];
Dlog=Dlog/.invrootrepls;

{OddPos[[pos]],Dlog},{pos,Range[Length[OddPos]]}],pos];


Export[CurrentDir<>"Dlogs.txt",IntegrandsMat]


DLogs=ReadList[CurrentDir<>"Dlogs.txt"];
Integrands=ReadList[CurrentDir<>"Integrands.txt"];


(* ::Section:: *)
(*Find a factors*)


EvenLetterCandidates=ReadList[HomeDir<>"/Reconstruction-NewPara/EvenLetterCandidates.txt"];
EvenLetterCandidates=DeleteDuplicates[Join[EvenLetterCandidates,aFactorsReconstructed]];
ELSclice=Factor[EvenLetterCandidates/.SliceTpoint,Modulus->$P];
ELSclice//Length
EvenLetterCandidates//Length


Monitor[ABReco=Table[
LogArgNum=Numerator[DLogs[[pos,2]]/.Log[x_]:>x];
LogArgDen=Denominator[DLogs[[pos,2]]/.Log[x_]:>x];
FSlice=Factor[Expand[LogArgNum*LogArgDen,Modulus->$P]/.root[x_]^2:>x,Modulus->$P];
sqrt1=rootVector[[DLogs[[pos,1,1]]]];
sqrt2=rootVector[[DLogs[[pos,1,2]]]];
r=sqrt1 sqrt2;
vars=Variables[r];
quadvar=vars[[First[Flatten[Position[Map[Exponent[r,#]&,vars],2]]]]];

b=Factor[Coefficient[LogArgNum,Times@@DeleteDuplicates[Cases[LogArgNum,_root,Infinity]]],Modulus->$P];
FFactors=Times@@Map[CheckLetterFactor[FSlice,#]&,Range[191]];

If[!NumberQ[Expand[FSlice/FFactors/.Fac[n_]:>ELSclice[[n]],Modulus->$P]],
FFactors=$Failed];
If[!NumberQ[b],FactorsOdd=Times@@Map[CheckLetterFactor[b,#]&,Range[191]];If[FactorsOdd==1,FactorsOdd=$Failed],FactorsOdd=b];
If[!NumberQ[Expand[b/FactorsOdd/.Fac[n_]:>ELSclice[[n]],Modulus->$P]],
FactorsOdd=$Failed];

FFull=Factor[FFactors/.Fac[n_]:>EvenLetterCandidates[[n]]];
bFull=Factor[FactorsOdd/.Fac[n_]:>EvenLetterCandidates[[n]]];

If[(FactorsOdd=!=$Failed) && (FFactors=!=$Failed),
sol=Solve[Expand[Discriminant[Expand[k FFull+bFull^2 *r/.Delete[PhaseSpacePoints[[1]],First[Flatten[Position[PhaseSpacePoints[[1,All,1]],quadvar]]]],Modulus->$P],quadvar]/.PhaseSpacePoints[[1]],Modulus->$P]==0,k,Modulus->$P];
aFull=root[Factor[-4 FFull+bFull^2 *r/.sol]];
sol=Flatten[sol/.n_Integer:>RationalReconstruct[n]][[All,2]];
sol=DeleteElements[sol,{0}];
sol=Delete[sol,Position[Denominator[sol], x_ /; x > 500]];
aFull=root[Factor[sol[[1]] FFull+bFull^2 *r]];
aFull=aFull/.root[x_^2]:>x;
,aFull=$Failed,bFull=$Failed;];
usedroots=root[sqrt1]root[sqrt2]/.root[1]->1;
{OddPos[[pos]],aFull,bFull,usedroots},{pos,Range[Length[OddPos]]}],pos]



DeleteDuplicates[Position[ABReco,$Failed][[All,1]]]//Length
ABReco//Length


ABReco=ABReco/.root[n_Integer x_^2]:>Sqrt[n]x;


Export[CurrentDir<>"OddLetterReco1AB.txt",ABReco]


DeleteDuplicates[ABReco[[DeleteDuplicates[Position[ABReco,$Failed][[All,1]]]]][[All,1,1]]]


(* ::Section::Closed:: *)
(*Debug*)


pos=223;
(*pos=5*)
LogArgNum=Numerator[DLogs[[pos,2]]/.Log[x_]:>x];
LogArgDen=Denominator[DLogs[[pos,2]]/.Log[x_]:>x];
FSlice=Factor[Expand[LogArgNum*LogArgDen,Modulus->$P]/.root[x_]^2:>x,Modulus->$P];
sqrt1=rootVector[[DLogs[[pos,1,1]]]];
sqrt2=rootVector[[DLogs[[pos,1,2]]]];
r=sqrt1 sqrt2;
vars=Variables[r];
quadvar=vars[[First[Flatten[Position[Map[Exponent[r,#]&,vars],2]]]]];

b=Factor[Coefficient[LogArgNum,Times@@DeleteDuplicates[Cases[LogArgNum,_root,Infinity]]],Modulus->$P]
FFactors=Times@@Map[CheckLetterFactor[FSlice,#]&,Range[191]]

If[!NumberQ[Expand[FSlice/FFactors/.Fac[n_]:>ELSclice[[n]],Modulus->$P]],
FFactors=$Failed];
If[!NumberQ[b],FactorsOdd=Times@@Map[CheckLetterFactor[b,#]&,Range[191]];If[FactorsOdd==1,FactorsOdd=$Failed],FactorsOdd=b];

FactorsOdd
If[!NumberQ[Expand[b/FactorsOdd/.Fac[n_]:>ELSclice[[n]],Modulus->$P]],
FactorsOdd=$Failed];

FFull=Factor[FFactors/.Fac[n_]:>EvenLetterCandidates[[n]]];
bFull=Factor[FactorsOdd/.Fac[n_]:>EvenLetterCandidates[[n]]];

If[(FactorsOdd=!=$Failed) && (FFactors=!=$Failed),
sol=Solve[Expand[Discriminant[Expand[k FFull+bFull^2 *r/.Delete[PhaseSpacePoints[[1]],First[Flatten[Position[PhaseSpacePoints[[1,All,1]],quadvar]]]],Modulus->$P],quadvar]/.PhaseSpacePoints[[1]],Modulus->$P]==0,k,Modulus->$P];
aFull=root[Factor[-4 FFull+bFull^2 *r/.sol]];
sol=Flatten[sol/.n_Integer:>RationalReconstruct[n]][[All,2]];
sol=DeleteElements[sol,{0}];
sol=Delete[sol,Position[Denominator[sol], x_ /; x > 500]];
aFull=root[Factor[sol[[1]] FFull+bFull^2 *r]];
aFull=aFull/.root[x_^2]:>x;
,aFull=$Failed,bFull=$Failed;];
usedroots=root[sqrt1]root[sqrt2]/.root[1]->1;


sol=Solve[Expand[Discriminant[Expand[k FFull+bFull^2 *r/.Delete[PhaseSpacePoints[[1]],First[Flatten[Position[PhaseSpacePoints[[1,All,1]],quadvar]]]],Modulus->$P],quadvar]/.PhaseSpacePoints[[1]],Modulus->$P]==0,k,Modulus->$P]

aFull=root[Factor[-4 FFull+bFull^2 *r/.sol]];
sol=Flatten[sol/.n_Integer:>RationalReconstruct[n]][[All,2]]
sol=DeleteElements[sol,{0}]
sol=Delete[sol,Position[Denominator[sol], x_ /; x > 500]]
aFull=root[Factor[sol[[1]] FFull+bFull^2 *r]];





aFull=root[Factor[-4 FFull+bFull^2 *r/.sol]];


aFull





2147483643//RationalReconstruct
