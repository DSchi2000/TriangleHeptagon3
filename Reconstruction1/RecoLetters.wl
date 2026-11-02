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
OutPutPath="Samples/Full-Samples-500";
phasespacepoints=ReadList[FileNameJoin[{CurrentDir,ToString[StringForm["``/phasespacepoints.txt",OutPutPath]]}]];
(*tvarlist=Flatten[ReadList[FileNameJoin[{CurrentDir,ToString[StringForm["``/tvars.txt",OutPutPath]]}]][[;;;; 1]]];*)
phasespacepoints=phasespacepoints[[All,2;;]][[1;; ;; 1]];
diffInvariants={mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56};
PhaseSpacePoints=Map[Thread[Rule[diffInvariants,#]]&,phasespacepoints];
(*PhaseSpacePoints=Table[Join[PhaseSpacePoints[[i]],{tvarlist[[i]]}],{i,1,Length[PhaseSpacePoints]}];*)
(*PhaseSpacePoints=Map[Thread[Rule[diffInvariants,#]]&,phasespacepoints];*)
rootSet=(Cases[PureIntegrals[[-5;;]],_root,Infinity]//DeleteDuplicates)/.root[x_]:>x;
GetLinPos[x1_,x2_,NMI_]:=First[Flatten[Position[Flatten[Table[{i,j},{i,1,NMI},{j,1,NMI}],1],{x1,x2}]]];
Get[ToString[StringJoin[HomeDir,"/files/ExtraModules.m"]]];
Get[ToString[StringJoin[HomeDir,"/files/ExtraModules.m"]]];
usedvars={mt2,s12,s123,s23,s234,s34,s345,s45,s56};
For[i=1,i<=Length[usedvars],i++,
LinearTerm[usedvars[[i]]]=ReadList[FileNameJoin[{Directory[],ToString[StringForm["``/LinTerm``_0001.txt",OutPutPath,usedvars[[i]]]]}]];
];




(* ::Section:: *)
(*Load More Suff For Dlogging:*)


SolvedAlready=ReadList[CurrentDir<>"Solved.txt"];
FactorsWithPositions=ReadList[CurrentDir<>"EvenLetterPos.txt"];
FactorsWithPositions2=ReadList[CurrentDir<>"EvenLetterPos2.txt"];
EvenLetterCandidatesFACS=ReadList[CurrentDir<>"Factors.txt"];
DtIntegrands=ReadList[CurrentDir<>"EvenIntegrands.txt"];
dlogList=ReadList[CurrentDir<>"dlogs.txt"];
VarGrid=ReadList[CurrentDir<>"VarGrid.txt"];
Nterms=ReadList[CurrentDir<>"TermsNeededPowers.txt"];
AnsatzPowers=ReadList[CurrentDir<>"RemainingPowers.txt"];
ABVEC=ReadList[CurrentDir<>"/Samples/AVECS.txt"];
avec=ABVEC[[1]]
bvec=ABVEC[[2]]
SliceTpoint=Thread[Rule[bvec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]];





(* ::Section:: *)
(*Variable Shifting*)


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
G6=Get[HomeDir<>"/SetUp/SixPointGram.txt"];
pos=GetLinPos[5,4,95];
DIsij=Map[LinearTerm[#][[1,1,pos]]&,usedvars]
DIsij=Insert[DIsij,0,4]
tpoint=PhaseSpacePoints[[1]]
ChangeAuxVarInDe[s16,s12,DIsij,tpoint]


Nterms
AnsatzPowers[[evenpos]]


(* ::Section:: *)
(*Recontruct 1*)


DIsij=Table[Insert[Map[LinearTerm[#][[psp,1,pos]]&,usedvars],0,4],{psp,1,Length[PhaseSpacePoints]}];
DIsij=Table[Coefficient[ChangeAuxVarInDe[s16,s12,DIsij[[psp]],PhaseSpacePoints[[psp]]],Table[ds[var],{var,diffInvariants}]],{psp,1,Length[PhaseSpacePoints]}];


var=s123;
vars=VarGrid[[50,35]]

fitvec=DIsij[[All,First[Flatten[Position[diffInvariants,var]]]]];

testpos=Flatten[Position[FactorsWithPositions[[All,1]],{50,35}]]//First
FactorsWithPositions2[[testpos]]
evenpos=FactorsWithPositions2[[testpos]][[1]];
DtIntegrands[[evenpos]]
Expand[DtIntegrands[[evenpos]],Modulus->$P]
dlogList[[evenpos]]
logRational=dlogList[[evenpos]]/.n_ Log[x_]:>n/.Log[x_]:>1
DlogDiff=Factor[Expand[logRational  D[Log[FactorsWithPositions2[[testpos]][[2]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]],var]],Modulus->$P]
fitvec=Expand[fitvec-(DlogDiff/.PhaseSpacePoints),Modulus->$P];
Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
Reco=ReconstructFUNCS[fitvec,PhaseSpacePoints[[All,Flatten[Map[Position[diffInvariants,#]&,vars]],2]],vars,4,5,Debug->True]
Ansatz=Reco[[2]];
ReconstructedFunc=Together[Ansatz/.Thread[Rule[Sort[Cases[Ansatz,_fac,Infinity]],Reco[[1,1]]]]/.n_Integer:>RationalReconstruct[n]]


(* ::Section:: *)
(*Integrate*)


Get["/home/daniel/Documents/packages/trager/Trager.wl"]
(*Get["/home/dschiebe/Documents/packages/trager/Trager.wl"]*)


TragerLogSubtraction[ReconstructedFunc,var,Debug->True]


(* ::Section:: *)
(*Residual Power Funcs*)


GetResidualPowers[matpos_]:=Module[{testpos,evenpos,logRational,DlogDiff,NewRecoFunc,NumPow,DenPow},
testpos=Flatten[Position[FactorsWithPositions[[All,1]],matpos]]//First;
evenpos=FactorsWithPositions2[[testpos]][[1]];
logRational=dlogList[[evenpos]]/.n_ Log[x_]:>n/.Log[x_]:>1;
DlogDiff=Factor[Expand[logRational  D[Log[FactorsWithPositions2[[testpos]][[2]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]/.SliceTpoint],t]]/.SliceTpoint,Modulus->$P];
NewRecoFunc=Factor[DtIntegrands[[evenpos]]-DlogDiff,Modulus->$P];
NumPow=Exponent[Numerator[NewRecoFunc],t];
DenPow=Exponent[Denominator[NewRecoFunc],t];
Return[{matpos,NumPow,DenPow}]

]


(* ::Section:: *)
(*Get Remaining Powers:*)


GetNumberOfTermsNeeded[n_,num_,den_]:=Binomial[n + num, n]+Binomial[n + den, n]
RecoPowers=Map[GetResidualPowers[#]&,Complement[FactorsWithPositions[[All,1]],SolvedAlready]]

NAnsaetze=Table[{RecoPowers[[i,1]],GetNumberOfTermsNeeded[Length[VarGrid[[RecoPowers[[i,1,1]],RecoPowers[[i,1,2]]]]],RecoPowers[[i,2]],RecoPowers[[i,3]]]},{i,1,Length[RecoPowers]}]


(* ::Section:: *)
(*BackUp*)


AnsatzPowers[[1]]


testpos=Flatten[Position[FactorsWithPositions[[All,1]],{5,4}]]//First
FactorsWithPositions2[[testpos]]
evenpos=FactorsWithPositions2[[testpos]][[1]];
DtIntegrands[[evenpos]]
Expand[DtIntegrands[[evenpos]],Modulus->$P]
dlogList[[evenpos]]
logRational=dlogList[[evenpos]]/.n_ Log[x_]:>n/.Log[x_]:>1
DlogDiff=Factor[Expand[logRational  D[Log[FactorsWithPositions2[[testpos]][[2]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]/.SliceTpoint],t]]/.SliceTpoint,Modulus->$P]
NewRecoFunc=Factor[DtIntegrands[[evenpos]]-DlogDiff,Modulus->$P];
NumPow=Exponent[Numerator[NewRecoFunc],t]
DenPow=Exponent[Denominator[NewRecoFunc],t]



