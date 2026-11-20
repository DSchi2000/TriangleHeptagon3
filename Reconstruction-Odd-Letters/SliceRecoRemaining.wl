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
OutPutPath="Reconstruction1/Samples/UnivariateSlice-3";
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









(*GenerateNumerOfSamples[num_,den_,vars_,roots_,rootAnsatzN_,rootAnsatzD_]:=Module[{NROOTS,NSqrtDegrees,rtsetNNR,SqrtDegrees,rootAnsDim,rootAnsMono,AnsatzPowersNumi,AnsatzPowersDeno,AnsatzLengthN,AnsatzLengthD,AnsatzLength},
NROOTS=DeleteDuplicates[Cases[roots,_root,Infinity]]/.root[x_]:>x;
NSqrtDegrees=(ResourceFunction["PolynomialTotalDegree"][#,diffInvariants]&/@NROOTS)/2;
rtsetNNR=roots/.Table[Rule[root[NROOTS[[i]]],Power[diffInvariants[[1]],NSqrtDegrees[[i]]]],{i,1,Length[NSqrtDegrees]}];
SqrtDegrees=(ResourceFunction["PolynomialTotalDegree"][#,diffInvariants]&/@rtsetNNR)/2;
rootAnsDim=Replace[MonomialList[rootAnsatzN],n_Integer:>0,{1}]/.Times->Plus/.Thread[Rule[Table[rt[i],{i,1,Length[roots]}],SqrtDegrees]];
rootAnsMono=MonomialList[rootAnsatzN];
AnsatzPowersNumi=Flatten[Table[Join[#,CoefficientRules[rootAnsMono[[sqrtPos]],Table[rt[i],{i,1,Length[roots]}]][[1,1]]]&/@FrobeniusSolve[ConstantArray[1,vars],num-rootAnsDim[[sqrtPos]]],{sqrtPos,1,Length[rootAnsMono]}],1];
(*AnsatzPowersNumi//Echo;*)	
	
rootAnsDim=Replace[MonomialList[rootAnsatzD],n_Integer:>0,{1}]/.Times->Plus/.Thread[Rule[Table[rt[i],{i,1,Length[roots]}],SqrtDegrees]];
rootAnsMono=MonomialList[rootAnsatzD];
AnsatzPowersDeno=Flatten[Table[Join[#,CoefficientRules[rootAnsMono[[sqrtPos]],Table[rt[i],{i,1,Length[roots]}]][[1,1]]]&/@FrobeniusSolve[ConstantArray[1,vars],den-rootAnsDim[[sqrtPos]]],{sqrtPos,1,Length[rootAnsMono]}],1];
	
	
	AnsatzLengthN=Length[AnsatzPowersNumi];
	AnsatzLengthD=Length[AnsatzPowersDeno];
	AnsatzLength=AnsatzLengthN+AnsatzLengthD;
Return[AnsatzLength]];*)
GetNumberOfTermsNeeded[n_,num_]:=Binomial[n + num-1, n-1]


(* ::Section:: *)
(*Find Even Letter Pos*)





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

(*OddPos=PrevReco[[DeleteDuplicates[Position[PrevReco,$Failed][[All,1]]]]][[All,1]]*)


PrevReco=SortBy[Join[PrevReco1,PrevReco2,PrevReco3],First]


OddPosRem=Complement[OddPos,PrevReco[[All,1]]]


(* ::Section::Closed:: *)
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
pos=GetLinPos[34,34,95];
DIsij=Map[LinearTerm[#][[1,1,pos]]&,usedvars]
DIsij=Insert[DIsij,0,4]
tpoint=PhaseSpacePoints[[1,;;-2]]
ChangeAuxVarInDe[s16,s12,DIsij,tpoint]


(* ::Section:: *)
(*LoadEvenLetters*)


EQN=Flatten[Table[{avector[i]+PhaseSpacePoints[[1,-1,2]] bvector[i] ==PhaseSpacePoints[[1,i,2]],avector[i]+PhaseSpacePoints[[2,-1,2]] bvector[i] ==PhaseSpacePoints[[2,i,2]]},{i,1,Length[PhaseSpacePoints[[1,All,1]]]-1}]];
solset=Flatten[Table[{avector[i],bvector[i]},{i,1,Length[PhaseSpacePoints[[1,All,1]]]-1}]];
sol=Solve[EQN,solset,Modulus->$P][[1]];
avec=Thread[Rule[PhaseSpacePoints[[1,;;-2,1]],Table[avector[i],{i,1,Length[PhaseSpacePoints[[1,All,1]]]-1}] ]]/.sol
bvec=Thread[Rule[PhaseSpacePoints[[1,;;-2,1]],Table[bvector[i],{i,1,Length[PhaseSpacePoints[[1,All,1]]]-1}] ]]/.sol
LinearTerm[s16]=SparseArray[ConstantArray[0,{Length[PhaseSpacePoints],1,Length[PureIntegrals]^2}]];
SliceTpoint=Thread[Rule[bvec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]]
SliceTpoint=Thread[Rule[bvec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]];
Expand[SliceTpoint/.t->PhaseSpacePoints[[-1,-1,2]],Modulus->$P]
PhaseSpacePoints[[-1]]
pos=1;


EvenLetterCandidates=ReadList[HomeDir<>"/Reconstruction-NewPara/EvenLetterCandidates.txt"];
(*EvenLetterCandidates=DeleteDuplicates[Join[EvenLetterCandidates,AFACS]];*)
ELSclice=Factor[EvenLetterCandidates/.SliceTpoint,Modulus->$P];
ELSclice//Length
EvenLetterCandidates//Length














(* ::Section:: *)
(*Check Problems:*)


DLogs=ReadList[CurrentDir<>"Dlogs.txt"];
IntegrandsReco=ReadList[CurrentDir<>"Integrands-Sorted.txt"];
OddPosRem
Position[OddPos,{43,35}]
VarGrid[[43,35]]





RemPos=Flatten[Map[Position[OddPos,#]&,OddPosRem]];
VarGrid=ReadList[HomeDir<>"/Reconstruction1/VarGrid.txt"];
ExponentsNumi=Exponent[Numerator[IntegrandsReco[[RemPos,2]]],t];
ExponentsNumi


Powers=Map[Min[DeleteElements[#,{-Infinity}]]&,ExponentsNumi]
Vars=Map[Length[VarGrid[[#[[1]],#[[2]]]]]&,OddPosRem]


Table[]


Table[GetNumberOfTermsNeeded[Vars[[i]],Powers[[i]]],{i,1,Length[Powers]}]


(* ::Section:: *)
(*Reconstruction*)


EQN=Flatten[Table[{avector[i]+PhaseSpacePoints[[1,-1,2]] bvector[i] ==PhaseSpacePoints[[1,i,2]],avector[i]+PhaseSpacePoints[[2,-1,2]] bvector[i] ==PhaseSpacePoints[[2,i,2]]},{i,1,Length[PhaseSpacePoints[[1,All,1]]]-1}]];
solset=Flatten[Table[{avector[i],bvector[i]},{i,1,Length[PhaseSpacePoints[[1,All,1]]]-1}]];
sol=Solve[EQN,solset,Modulus->$P][[1]];
avec=Thread[Rule[PhaseSpacePoints[[1,;;-2,1]],Table[avector[i],{i,1,Length[PhaseSpacePoints[[1,All,1]]]-1}] ]]/.sol
bvec=Thread[Rule[PhaseSpacePoints[[1,;;-2,1]],Table[bvector[i],{i,1,Length[PhaseSpacePoints[[1,All,1]]]-1}] ]]/.sol
LinearTerm[s16]=SparseArray[ConstantArray[0,{Length[PhaseSpacePoints],1,Length[PureIntegrals]^2}]];
SliceTpoint=Thread[Rule[bvec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]]
SliceTpoint=Thread[Rule[bvec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]];
Expand[SliceTpoint/.t->PhaseSpacePoints[[-1,-1,2]],Modulus->$P]
PhaseSpacePoints[[-1]]
pos=1;



(* ::Subsection:: *)
(*Reco*)


Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
pos1=43;
pos2=35;
sqrt1=rootVector[[pos1]]
sqrt2=rootVector[[pos2]]
pos=GetLinPos[pos1,pos2,95];
sqrt1N=PowerMod[sqrt1/.PhaseSpacePoints,1/2,$P];
sqrt2N=PowerMod[sqrt2/.PhaseSpacePoints,1/2,$P];
(*pos=LinearizedEvenPos[[290]]*)
FitVals=Join[Table[LinearTerm[var][[All,1,pos]],{var,diffInvariants}]];
FitVals=Table[Coefficient[ChangeAuxVarInDe[s12,s23,Transpose[FitVals][[i]],PhaseSpacePoints[[i]]],Thread[ds[diffInvariants]]],{i,1,Length[PhaseSpacePoints]}];
Integrands=Table[Integrand=Expand[UniVarReconstructFUNCS[Expand[FitVals[[All,varpos]]/(sqrt2N/sqrt1N),Modulus->$P],PhaseSpacePoints,tvar,40,40,{},1,1],Modulus->$P];Factor[Integrand[[1]],Modulus->$P],{varpos,1,10}]
NewIntDt=Factor[Expand[Integrands . bvec[[All,2]],Modulus->$P],Modulus->$P]
NewIntDt=NewIntDt*(root[sqrt2]/root[sqrt1])/.root[1]->1/.SliceTpoint/.root[x_]:>root[Factor[x,Modulus->$P]]


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
Monitor[IntegrandsMat=Table[
pos1=matpos[[1]];
pos2=matpos[[2]];
sqrt1=rootVector[[pos1]];
sqrt2=rootVector[[pos2]];
pos=GetLinPos[pos1,pos2,95];
sqrt1N=PowerMod[sqrt1/.PhaseSpacePoints,1/2,$P];
sqrt2N=PowerMod[sqrt2/.PhaseSpacePoints,1/2,$P];
FitVals=Join[Table[LinearTerm[var][[All,1,pos]],{var,diffInvariants}]];
FitVals=Join[Table[LinearTerm[var][[All,1,pos]],{var,diffInvariants}]];
usedvars=VarGrid[[pos1,pos2]];
If[Length[usedvars]<10,newaux=Last[Complement[diffInvariants,usedvars]];
FitVals=ParallelTable[Coefficient[ChangeAuxVarInDe[s16,newaux,Transpose[FitVals][[i]],PhaseSpacePoints[[i]]],Thread[ds[diffInvariants]]],{i,1,Length[PhaseSpacePoints]}];,
FitVals=Transpose[FitVals]];

Integrands=Table[Integrand=Expand[UniVarReconstructFUNCS[Expand[FitVals[[All,varpos]]/(sqrt2N/sqrt1N),Modulus->$P],PhaseSpacePoints,tvar,20,20,{},1,1],Modulus->$P];Factor[Integrand[[1]],Modulus->$P],{varpos,1,10}];
Integrands=Integrands*(root[sqrt2]/root[sqrt1])/.root[1]->1/.SliceTpoint/.root[x_]:>root[Factor[x,Modulus->$P]];
{matpos,Integrands},{matpos,OddPos}],matpos];


(*NewIntDt=Factor[Expand[Integrands . bvec[[All,2]],Modulus->$P],Modulus->$P];
NewIntDt=NewIntDt*(root[sqrt2]/root[sqrt1])/.root[1]->1/.SliceTpoint/.root[x_]:>root[Factor[x,Modulus->$P]];*)


Export[CurrentDir<>"Integrands-Sorted.txt",IntegrandsMat]


IntegrandsMat[[-100]]


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


PrevReco=ReadList[CurrentDir<>"OddLetterReco1AB.txt"];
PrevReco2=ReadList[CurrentDir<>"OddLetterReco2AB.txt"];
AFACS1=DeleteElements[DeleteDuplicates[Flatten[Map[FactorList[#][[All,1]]&,DeleteDuplicates[DeleteElements[PrevReco[[All,2]],{$Failed}]]/.root[n_Integer x_^2 y_^2]:>Sqrt[n]x y]]],{-1,1}];
AFACS2=DeleteElements[DeleteDuplicates[Flatten[Map[FactorList[#][[All,1]]&,DeleteDuplicates[DeleteElements[PrevReco2[[All,2]],{$Failed}]]/.root[n_Integer x_^2 y_^2]:>Sqrt[n]x y]]],{-1,1}];
AFACS=DeleteDuplicates[Join[AFACS1,AFACS2]];

Complement[AFACS2,AFACS1]


EvenLetterCandidates=ReadList[HomeDir<>"/Reconstruction-NewPara/EvenLetterCandidates.txt"];
EvenLetterCandidates=DeleteDuplicates[Join[EvenLetterCandidates,AFACS]];
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
FFactors=Times@@Map[CheckLetterFactor[FSlice,#]&,Range[Length[EvenLetterCandidates]]];

If[!NumberQ[Expand[FSlice/FFactors/.Fac[n_]:>ELSclice[[n]],Modulus->$P]],
FFactors=$Failed];
If[!NumberQ[b],FactorsOdd=Times@@Map[CheckLetterFactor[b,#]&,Range[Length[EvenLetterCandidates]]];If[FactorsOdd==1,FactorsOdd=$Failed],FactorsOdd=b];
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


ABReco=ABReco/.root[n_Integer x_^2]:>Sqrt[n]x/.root[n_Integer x_^2 y_^2]:>Sqrt[n]x y;


Export[CurrentDir<>"OddLetterReco2AB.txt",ABReco]


DeleteDuplicates[ABReco[[DeleteDuplicates[Position[ABReco,$Failed][[All,1]]]]][[All,1,1]]]


ABReco[[All,2]]


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
