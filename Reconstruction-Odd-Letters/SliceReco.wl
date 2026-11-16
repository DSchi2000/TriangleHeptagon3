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









GenerateNumerOfSamples[num_,den_,vars_,roots_,rootAnsatzN_,rootAnsatzD_]:=Module[{NROOTS,NSqrtDegrees,rtsetNNR,SqrtDegrees,rootAnsDim,rootAnsMono,AnsatzPowersNumi,AnsatzPowersDeno,AnsatzLengthN,AnsatzLengthD,AnsatzLength},
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
Return[AnsatzLength]];


(* ::Section:: *)
(*Find Even Letter Pos*)


rootVector=Map[Times@@DeleteDuplicates[Cases[#,_root,Infinity]]&,PureIntegrals]/.{}->1/.root[x_]:>x;
DeFoldMat=Normal[Total[Table[ArrayReshape[LinearTerm[var][[1,1]],{Length[PureIntegrals],Length[PureIntegrals]}],{var,usedvars}]]]/.n_Integer/;(n>0):>1;


rootVectorN=rootVector/.PhaseSpacePoints[[1]];
RootRatioMat=Expand[Table[Together[rootVectorN[[i1]]/rootVectorN[[i2]]],{i1,1,Length[PureIntegrals]},{i2,1,Length[PureIntegrals]}],Modulus->$P]*DeFoldMat;
EvenPos=Position[RootRatioMat,1];
(*LinearizedEvenPos=Map[GetLinPos[#[[1]],#[[2]],95]&,EvenPos]*)
OddPos=Complement[Position[DeFoldMat,1],EvenPos]


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
pos=GetLinPos[34,34,95];
DIsij=Map[LinearTerm[#][[1,1,pos]]&,usedvars]
DIsij=Insert[DIsij,0,4]
tpoint=PhaseSpacePoints[[1,;;-2]]
ChangeAuxVarInDe[s16,s12,DIsij,tpoint]


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



(* ::Subsection:: *)
(*Reco*)


Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
pos1=93;
pos2=3;
sqrt1=rootVector[[pos1]]
sqrt2=rootVector[[pos2]]
pos=GetLinPos[pos1,pos2,95];
sqrt1N=PowerMod[sqrt1/.PhaseSpacePoints,1/2,$P];
sqrt2N=PowerMod[sqrt2/.PhaseSpacePoints,1/2,$P];
(*pos=LinearizedEvenPos[[290]]*)
FitVals=Join[Table[LinearTerm[var][[All,1,pos]],{var,diffInvariants}]];
Integrands=Table[Integrand=Expand[UniVarReconstructFUNCS[Expand[FitVals[[varpos]]/(sqrt2N/sqrt1N),Modulus->$P],PhaseSpacePoints,tvar,40,40,{},1,1],Modulus->$P];Factor[Integrand[[1]],Modulus->$P],{varpos,1,10}]
NewIntDt=Factor[Expand[Integrands . bvec[[All,2]],Modulus->$P],Modulus->$P]
NewIntDt=NewIntDt*root[sqrt2]/root[sqrt1]/.root[1]->1/.SliceTpoint/.root[x_]:>root[Factor[x,Modulus->$P]]


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


(* ::Section::Closed:: *)
(*SATISFACTORY YEAAAAH*)


Monitor[Do[
pos=LinearizedEvenPos[[position]];
FitVals=Join[Table[LinearTerm[var][[All,1,pos]],{var,invariants2}]];
Integrands=Table[Integrand=Expand[UniVarReconstructFUNCS[FitVals[[varpos]],PhaseSpacePoints,tvar,20,20,{},1,1],Modulus->$P];Factor[Integrand[[1]],Modulus->$P],{varpos,1,10}];
NewIntDt=Factor[Expand[Integrands . bvec[[All,2]],Modulus->$P],Modulus->$P];
dtIntegrands[position]=NewIntDt;
DLogs[position]=TragerLogSubtractionNumeric[NewIntDt,t,{},{},Modulus->$P];
,{position,1,Length[LinearizedEvenPos]}],{ProgressIndicator[position,{1,Length[LinearizedEvenPos]}],position}]


?ProgressIndicator


dlogList=Table[DLogs[position],{position,1,Length[LinearizedEvenPos]}]
Export[CurrentDir<>"dlogs.txt",dlogList]
IntegrandsSlice=Table[dtIntegrands[position],{position,1,Length[LinearizedEvenPos]}]


dlogList//DeleteDuplicates//Length



Export[CurrentDir<>"EvenIntegrands.txt",IntegrandsSlice]


(* ::Section:: *)
(*Find DlogFactors*)


dlogList=ReadList[CurrentDir<>"dlogs.txt"];
SlicedLetters=DeleteDuplicates[Cases[dlogList,_Log,Infinity]]/.Log[x_]:>x
(*Export[CurrentDir<>"SlicedEvenLetters.txt",SlicedLetters]*)


(* ::Subsection:: *)
(*Load Claude Guesses*)


ABVEC=ReadList[CurrentDir<>"/Samples/AVECS.txt"];
avec=ABVEC[[1]];
bvec=ABVEC[[2]];
SliceTpoint=Thread[Rule[bvec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]];
EvenLetterCandidates=DeleteCases[ReadList[HomeDir<>"/ClaudeLetters/EvenLetterCandidates.txt"],Null];
EvenLetterCandidatesFACS=DeleteElements[DeleteDuplicates[Flatten[Map[FactorList[#][[All,1]]&,EvenLetterCandidates]]],{1,-1}];
Position[Factor[EvenLetterCandidatesFACS/.SliceTpoint,Modulus->$P],0];
EvenLetterCandidatesFACS=Delete[EvenLetterCandidatesFACS,79];
ELSclice=Factor[EvenLetterCandidatesFACS/.SliceTpoint,Modulus->$P];


(* ::Subsection:: *)
(*Load First 500 Set:*)


NewDlogs1=ReadList[CurrentDir<>"NewDlogs-500.txt"];
NewDlogs2=ReadList[CurrentDir<>"NewDlogs-500-2.txt"];
NewDlogs3=ReadList[CurrentDir<>"Dlog-94-94.txt"];
fiveptgrams=ReadList[CurrentDir<>"grams5pt.txt"];
EvenLetterCandidatesFACSNew=DeleteElements[DeleteDuplicates[Flatten[Map[FactorList[#][[All,1]]&,DeleteDuplicates[Cases[{NewDlogs1,NewDlogs2,NewDlogs3},_Log,Infinity]]/.Log[x_]:>x]]],{1,-1}];
EvenLetterCandidatesFACSNew=DeleteDuplicates[Join[EvenLetterCandidatesFACSNew,fiveptgrams]];
EvenLetterCandidatesFACSNew//Length
FancyDeleteDuplicates[delset_]:=Module[{size,i,j,delset2,factor},delset2=delset;size=Length[delset2];
For[i=1,i<=size,i++,
	For[j=i,j<=size,j++,factor=Factor[delset2[[i]]/delset2[[j]]];
		If[NumberQ[factor],delset2[[j]]=delset2[[i]];];];];
Return[DeleteDuplicates[delset2]];]
EvenLetterCandidatesFACS=FancyDeleteDuplicates[Join[EvenLetterCandidatesFACS,EvenLetterCandidatesFACSNew]];
ELSclice=Factor[EvenLetterCandidatesFACS/.SliceTpoint,Modulus->$P];
ELSclice//Length
EvenLetterCandidatesFACS//Length


EvenLetterCandidatesFACSNew[[-1]]


NewDlogs3


(* ::Subsection:: *)
(*Find Factors*)


(*Func=SlicedLetters[[15]]
TestLetter=ELSclice[[1]];*)
pos=1;
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
(*
FoundFactors=Table[Times@@Map[CheckLetterFactor[Letter,#]&,Range[Length[ELSclice]]],{Letter,SlicedLetters}];
*)



ELSclice[[-1]]


FoundFactors[[-2]]



Factor[Numerator[(NewDlogs3/.Log[x_]:>x)]/.SliceTpoint,Modulus->$P]
Factor[Denominator[(NewDlogs3/.Log[x_]:>x)]/.SliceTpoint,Modulus->$P]
Factor[SlicedLetters[[-2]],Modulus->$P]


Factor[D[dlogList[[-2]],t],Modulus->$P]


(* ::Section:: *)
(*Checked fixed Entries*)


Exponents=Table[
NewFunc=Factor[SlicedLetters[[i]]/FoundFactors[[i]]/.Fac[n_]:>ELSclice[[n]],Modulus->$P];
{Exponent[Numerator[NewFunc],t],Exponent[Denominator[NewFunc],t]},{i,1,Length[SlicedLetters]}]
FullreResolved=Flatten[Position[Exponents,{0,0}]]
ResolvedEvenLettersPos=Map[DeleteElements[Table[If[Log[#]===(dlogList[[i]]/.n_ Log[x_]:>Log[x]),i,0],{i,1,Length[dlogList]}],{0}]&,SlicedLetters[[FullreResolved]]]
UsedLetters=DeleteDuplicates[Cases[FoundFactors[[FullreResolved]],_Fac,Infinity]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]];
SolvedMatrixPos=EvenPos[[Flatten[ResolvedEvenLettersPos]]]




OldSolved=ReadList[CurrentDir<>"Solved.txt"]


(* ::Subsection:: *)
(*Export*)


Export[CurrentDir<>"DlogFactors-3.txt",UsedLetters]
Export[CurrentDir<>"Solved-3.txt",SolvedMatrixPos]




AllEvenLettersPos=Map[DeleteElements[Table[If[Log[#]===(dlogList[[i]]/.n_ Log[x_]:>Log[x]),i,0],{i,1,Length[dlogList]}],{0}]&,SlicedLetters]
FactorsWithPositions=Flatten[Table[Table[{EvenPos[[pos]],FoundFactors[[i]]},{pos,AllEvenLettersPos[[i]]}] ,{i,1,Length[AllEvenLettersPos]}],1];
FactorsWithPositions2=Flatten[Table[Table[{pos,FoundFactors[[i]]},{pos,AllEvenLettersPos[[i]]}] ,{i,1,Length[AllEvenLettersPos]}],1];

Export[CurrentDir<>"Factors-3.txt",EvenLetterCandidatesFACS]
Export[CurrentDir<>"EvenLetterPos1-3.txt",FactorsWithPositions]
Export[CurrentDir<>"EvenLetterPos2-3.txt",FactorsWithPositions2]







DeleteDuplicates[Cases[FoundFactors,_Fac,Infinity]]//Sort


FactorsWithPositions[[-6]]


EvenLetterCandidatesFACS[[61]]


TestBas=ReadList["/home/daniel/Documents/PHD-LocalFiles/TriangleHeptagon3/Reconstruction1/SubBasis.txt"]


DeleteDuplicates[Cases[TestBas[[93]],_root,Infinity]]
