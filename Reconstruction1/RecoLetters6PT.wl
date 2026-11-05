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

OutPutPath="Samples/Full-Samples-500";
phasespacepoints=ReadList[FileNameJoin[{CurrentDir,ToString[StringForm["``/phasespacepoints.txt",OutPutPath]]}]];
phasespacepoints=phasespacepoints[[All,2;;]][[1;; ;; 1]];
PhaseSpacePoints=Map[Thread[Rule[diffInvariants,#]]&,phasespacepoints];
For[i=1,i<=Length[usedvars],i++,
LinearTerm2[usedvars[[i]]]=ReadList[FileNameJoin[{Directory[],ToString[StringForm["``/LinTerm``_0001.txt",OutPutPath,usedvars[[i]]]]}]];
];




(* ::Section:: *)
(*Load Samples*)


OutPutPaths={"Samples/Full-Samples-10000"};
SplitLengths={2};
PhaseSpacePoints={};
Do[
phasespacepointsIN=ReadList[FileNameJoin[{CurrentDir,ToString[StringForm["``/phasespacepoints.txt",OutPutPaths[[q]]]]}]];
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
	LinearTerm[usedvars[[i]],loadfile]=ReadList[FileNameJoin[{Directory[],ToString[StringForm["``/LinTerm``_``.txt",OutPutPath,usedvars[[i]],IntegerString[q,10,4]]]}]][[;;-1]];
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
(*Load More Suff For Dlogging:*)


SolvedAlready=ReadList[CurrentDir<>"Solved-3.txt"];
FactorsWithPositions=ReadList[CurrentDir<>"EvenLetterPos1-3.txt"];
FactorsWithPositions2=ReadList[CurrentDir<>"EvenLetterPos2-3.txt"];
EvenLetterCandidatesFACS=ReadList[CurrentDir<>"Factors-3.txt"];
DtIntegrands=ReadList[CurrentDir<>"EvenIntegrands.txt"];
dlogList=ReadList[CurrentDir<>"dlogs.txt"];
VarGrid=ReadList[CurrentDir<>"VarGrid.txt"];
Nterms=ReadList[CurrentDir<>"TermsNeededPowers-3.txt"];
AnsatzPowers=ReadList[CurrentDir<>"RemainingPowers-3.txt"];
ABVEC=ReadList[CurrentDir<>"/Samples/AVECS.txt"];
avec=ABVEC[[1]];
bvec=ABVEC[[2]];
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


(* ::Section:: *)
(*Recontruct 1*)


CloseKernels[];
LaunchKernels[20];


pos1=95;
pos2=95;
pos=GetLinPos[pos1,pos2,95];
DIsij=Table[Insert[Map[LinearTerm[#][[psp,1,pos]]&,usedvars],0,4],{psp,1,Length[PhaseSpacePoints]}]//EchoTiming;
(*DIsij=Table[Coefficient[ChangeAuxVarInDe[s16,s12,DIsij[[psp]],PhaseSpacePoints[[psp]]],Table[ds[var],{var,diffInvariants}]],{psp,1,Length[PhaseSpacePoints]}]//EchoTiming;*)
setup=AnsatzPowers[[Flatten[Position[Nterms[[All,1]],{pos1,pos2}]]//First]]
Nterms[[Flatten[Position[Nterms[[All,1]],{pos1,pos2}]]//First]]
vars=VarGrid[[setup[[1,1]],setup[[1,2]]]]//Sort
(*vars={s16,mt2}*)


var=s56;
(mt2^2 s12 s23-mt2 s12 s16 s23+2 mt2 s12 s16 s234-2 mt2 s123 s16 s234+mt2^2 s23 s234-mt2 s12 s23 s234-mt2 s12 s16 s34+2 mt2 s123 s16 s34+s12 s16^2 s34-2 s123 s16^2 s34-mt2^2 s23 s34-mt2 s12 s23 s34+mt2 s16 s23 s34-s12 s16 s23 s34-mt2 s16 s234 s34-s12 s16 s234 s34+2 s123 s16 s234 s34-mt2 s23 s234 s34+s12 s23 s234 s34+mt2 s16 s34^2-2 s123 s16 s34^2-s16^2 s34^2+mt2 s23 s34^2+s16 s23 s34^2-mt2 s12 s23 s345+s12 s16 s23 s345-mt2 s12 s234 s345+2 mt2 s123 s234 s345-s12 s16 s234 s345+2 s123 s16 s234 s345-2 mt2 s23 s234 s345-s12 s23 s234 s345+mt2 s234^2 s345+s12 s234^2 s345-2 s123 s234^2 s345-2 mt2 s123 s34 s345+2 s123 s16 s34 s345+mt2 s23 s34 s345-s16 s23 s34 s345-mt2 s234 s34 s345+2 s123 s234 s34 s345+2 s16 s234 s34 s345-s23 s234 s34 s345-2 s123 s234 s345^2+s23 s234 s345^2-s234^2 s345^2-mt2 s12 s234 s45-s12 s16 s234 s45-mt2 s234^2 s45+s12 s234^2 s45+mt2 s12 s34 s45-s12 s16 s34 s45+2 mt2 s234 s34 s45-s12 s234 s34 s45-s16 s234 s34 s45-mt2 s34^2 s45+s16 s34^2 s45+2 s12 s234 s345 s45+s234^2 s345 s45-s234 s34 s345 s45+mt2 s16 s23 s56+s16^2 s34 s56-s16 s23 s34 s56-mt2 s23 s345 s56-s16 s23 s345 s56-s16 s234 s345 s56+2 s23 s234 s345 s56-s16 s34 s345 s56-s23 s34 s345 s56+s23 s345^2 s56+s234 s345^2 s56+s16 s234 s45 s56-s16 s34 s45 s56-s234 s345 s45 s56+s34 s345 s45 s56)/(-mt2 s12^2 s23+mt2 s12 s123 s23+mt2 s12^2 s234-2 mt2 s12 s123 s234+mt2 s123^2 s234+mt2 s12 s123 s34-mt2 s123^2 s34-2 s12 s123 s16 s34+2 s123^2 s16 s34+2 mt2 s12 s23 s34-mt2 s123 s23 s34+s12 s123 s23 s34-mt2 s12 s234 s34+mt2 s123 s234 s34+s12 s123 s234 s34-s123^2 s234 s34-mt2 s123 s34^2+s123^2 s34^2+2 s123 s16 s34^2-mt2 s23 s34^2-s123 s23 s34^2-s12 s123 s23 s345+s12 s123 s234 s345-s123^2 s234 s345-s123^2 s34 s345+s123 s23 s34 s345-2 s123 s234 s34 s345+s12^2 s23 s45-s12^2 s234 s45+s12 s123 s234 s45+s12 s123 s34 s45-2 s12 s23 s34 s45+s12 s234 s34 s45+s123 s234 s34 s45-s123 s34^2 s45+s23 s34^2 s45+mt2 s12 s23 s56-mt2 s123 s23 s56+2 s12 s16 s34 s56-2 s123 s16 s34 s56+mt2 s23 s34 s56-2 s12 s23 s34 s56+s123 s23 s34 s56+s12 s23 s345 s56+s123 s23 s345 s56-s12 s234 s345 s56+s123 s234 s345 s56+s123 s34 s345 s56+s23 s34 s345 s56-2 s12 s23 s45 s56+s12 s234 s45 s56-s123 s234 s45 s56-2 s12 s34 s45 s56+s123 s34 s45 s56-2 s23 s34 s45 s56-s23 s345 s56^2+s23 s45 s56^2)
setup[[2]]
setup[[3]]

fitvec=DIsij[[All,First[Flatten[Position[diffInvariants,var]]]]];

testpos=Flatten[Position[FactorsWithPositions[[All,1]],setup[[1]]]]//First
FactorsWithPositions2[[testpos]]
evenpos=FactorsWithPositions2[[testpos]][[1]];
EqnConstFull=Sum[D[G6,var] ds[var],{var,diffInvariants}];
diffauxvar=First[Solve[EqnConstFull==0,ds[s16]]];
sixpointdiff=Coefficient[diffauxvar[[1,2]],ds[var]]


DtIntegrands[[evenpos]]
Expand[DtIntegrands[[evenpos]],Modulus->$P]
dlogList[[evenpos]]
logRational=dlogList[[evenpos]]/.n_ Log[x_]:>n/.Log[x_]:>1
logRational Log[FactorsWithPositions2[[testpos]][[2]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]]
DlogDiff=Factor[logRational  D[Log[LogArg/Fac[167] FactorsWithPositions2[[testpos]][[2]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]],var],Modulus->$P]//EchoTiming;
DlogDiffs16=Factor[sixpointdiff*logRational  D[ Log[LogArg/Fac[167] FactorsWithPositions2[[testpos]][[2]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]],s16],Modulus->$P]//EchoTiming;
fitvec=DIsij[[All,First[Flatten[Position[diffInvariants,var]]]]];
fitvec=Expand[fitvec[[;;10]]-ParallelMap[Expand[(DlogDiff+DlogDiffs16)/.#,Modulus->$P]&,PhaseSpacePoints[[;;10]]],Modulus->$P]//EchoTiming;
fitvec
(*Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
Reco=ReconstructOnVariety[fitvec,PhaseSpacePoints[[All,Flatten[Map[Position[diffInvariants,#]&,vars]],2]],vars,5,6,{G6,s16},Debug->True];
Ansatz=Reco[[2]];
ReconstructedFunc=Together[Ansatz/.Thread[Rule[Sort[Cases[Ansatz,_fac,Infinity]],Reco[[1,1]]]]/.n_Integer:>RationalReconstruct[n]]*)


s23


|


(* ::Section:: *)
(*Integrate*)


(*Get["/home/daniel/Documents/packages/trager/Trager.wl"]*)
Get["/home/dschiebe/Documents/packages/trager/Trager.wl"]
Dlog=TragerLogSubtraction[ReconstructedFunc,var,{s16},{G6},Parallelize->True,Debug->True]
Factor[Cases[{Dlog},_Log,Infinity]/.Log[x_]:>x]


Dlog1=Dlog;
Export[CurrentDir<>"Dlogs23-95-95.txt",Dlog1]


LogArg=Numerator[Dlog1/.Log[x_]:>x]


(*Dlog9494=Log[LogArg/Fac[167] FactorsWithPositions2[[testpos]][[2]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]]*)
Export[CurrentDir<>"FinishedDlog-94-94.txt",Dlog9494]


Dlog9595=Log[LogArg/Fac[167] FactorsWithPositions2[[testpos]][[2]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]]
Export[CurrentDir<>"FinishedDlog-95-95.txt",Dlog9595]


D[Log[Factor[LogArg/Fac[167] FactorsWithPositions2[[testpos]][[2]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]]],mt2]


Factor[Exp[Dlog9494]/.SliceTpoint,Modulus->$P]


Factor[Exp[dlogList[[-2]]]/.SliceTpoint,Modulus->$P]
(*dlogList[[-2]]*)


DlogFacs=DeleteElements[Join[FactorList[Numerator[Factor[Exp[Dlog9494]]]][[All,1]],
FactorList[Denominator[Factor[Exp[Dlog9494]]]][[All,1]]],{-1,1}]


DlogFacsSlice=Factor[DlogFacs/.SliceTpoint,Modulus->$P]


Factor[DlogFacsSlice[[1]]^2/Times@@DlogFacsSlice[[{2,3,4}]],Modulus->$P]


R=Factor[D[dlogList[[-2]],t],Modulus->$P]


var=t
L=DlogFacsSlice[[4]]/.x_^2:>x
B2 = Factor[Denominator[R]/L, Modulus -> $P]
rem = PolynomialRemainder[Numerator[R] - c*D[L, var]*B2, L, var];
csol = Solve[CoefficientList[rem, var] == 0, c,Modulus->$P]


Denominator[R]





DlogFacsSlice[[4]]/.x_^2:>x


B2


R


DlogFacs[[4]]


Expand[Factor[DlogDiff+DlogDiffs16,Modulus->$P]/.PhaseSpacePoints[[1]],Modulus->$P]
