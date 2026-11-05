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
OutPutPath="Samples/Univariate Slice";
phasespacepoints=ReadList[FileNameJoin[{CurrentDir,ToString[StringForm["``/phasespacepoints.txt",OutPutPath]]}]];
tvarlist=Flatten[ReadList[FileNameJoin[{CurrentDir,ToString[StringForm["``/tvars.txt",OutPutPath]]}]][[;;;; 1]]];
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
LinearTerm[usedvars[[i]]]=ReadList[FileNameJoin[{Directory[],ToString[StringForm["``/LinTerm``_0001.txt",OutPutPath,usedvars[[i]]]]}]];
];


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





(* ::Section::Closed:: *)
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





AnsatzPowers


(* ::Section:: *)
(*Recontruct 1*)


pos1=94;
pos2=94;
pos=GetLinPos[pos1,pos2,95];
DIsij=Table[Insert[Map[LinearTerm[#][[psp,1,pos]]&,usedvars],0,4],{psp,1,Length[PhaseSpacePoints]}]//EchoTiming;
(*DIsij=Table[Coefficient[ChangeAuxVarInDe[s16,s12,DIsij[[psp]],PhaseSpacePoints[[psp]]],Table[ds[var],{var,diffInvariants}]],{psp,1,Length[PhaseSpacePoints]}]//EchoTiming;*)
setup=AnsatzPowers[[Flatten[Position[Nterms[[All,1]],{pos1,pos2}]]//First]]
Nterms[[Flatten[Position[Nterms[[All,1]],{pos1,pos2}]]//First]]
vars=VarGrid[[setup[[1,1]],setup[[1,2]]]]//Sort

(*var=mt2;*)
setup[[2]]
setup[[3]]

(*fitvec=DIsij[[All,First[Flatten[Position[diffInvariants,var]]]]];*)

testpos=Flatten[Position[FactorsWithPositions[[All,1]],setup[[1]]]]//First
FactorsWithPositions2[[testpos]]
evenpos=FactorsWithPositions2[[testpos]][[1]];
var=s23

EqnConstFull=Sum[D[G6,var] ds[var],{var,diffInvariants}];
diffauxvar=First[Solve[EqnConstFull==0,ds[s16]]];
sixpointdiff=Coefficient[diffauxvar[[1,2]],ds[var]];

DtIntegrands[[evenpos]];
Expand[DtIntegrands[[evenpos]],Modulus->$P];
dlogList[[evenpos]];
logRational=dlogList[[evenpos]]/.n_ Log[x_]:>n/.Log[x_]:>1
logRational Log[FactorsWithPositions2[[testpos]][[2]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]]
DlogDiff=Factor[Expand[logRational  D[ Log[1/Fac[190] FactorsWithPositions2[[testpos]][[2]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]],var]],Modulus->$P];
DlogDiffs16=Factor[sixpointdiff*Expand[logRational  D[Log[ 1/Fac[190] FactorsWithPositions2[[testpos]][[2]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]],s16]],Modulus->$P];

(*DlogDiff=Factor[Expand[logRational  D[ Log[1/Fac[190]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]],var]],Modulus->$P];
DlogDiffs16=Factor[sixpointdiff*Expand[logRational  D[Log[ 1/Fac[190]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]],s16]],Modulus->$P];*)

fitvec=DIsij[[All,First[Flatten[Position[diffInvariants,var]]]]];
fitvec=Expand[fitvec-((DlogDiff+DlogDiffs16)/.PhaseSpacePoints),Modulus->$P];
Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
Integrand=Expand[UniVarReconstructFUNCS[fitvec,PhaseSpacePoints,tvar,20,20,{},1,1],Modulus->$P];
Expand[Integrand[[1]],Modulus->$P]
Factor[Integrand[[1]],Modulus->$P]

(*vars={s16,mt2}*)


D[ Log[1/Fac[190]*FactorsWithPositions2[[testpos]][[2]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]],s23]



NewDlogs3=ReadList[CurrentDir<>"Dlog-94-94.txt"];
Factor[(D[NewDlogs3[[1]],s23]+sixpointdiff D[NewDlogs3[[1]],s16])/.SliceTpoint,Modulus->$P]


(* ::Input:: *)
(**)


TESTF=Factor[(-mt2 s12 s123 s34+s12 s123 s16 s34-mt2 s123 s234 s34-s12 s123 s234 s34+mt2 s123 s34^2-s123 s16 s34^2+s123 s234 s34 s345+mt2 s12 s23 s56-mt2 s12 s234 s56+mt2 s123 s234 s56+2 mt2 s12 s34 s56-mt2 s123 s34 s56-2 s12 s16 s34 s56+s123 s16 s34 s56+mt2 s23 s34 s56+2 s12 s23 s34 s56-s12 s23 s345 s56+s12 s234 s345 s56-s123 s234 s345 s56-s23 s34 s345 s56-mt2 s23 s56^2+s23 s345 s56^2)/(-mt2 s12^2 s23 s234+mt2 s12 s123 s23 s234+mt2 s12^2 s234^2-2 mt2 s12 s123 s234^2+mt2 s123^2 s234^2+mt2 s12 s123 s23 s34-s12 s123 s16 s23 s34+2 mt2 s12 s123 s234 s34-2 mt2 s123^2 s234 s34-s12 s123 s16 s234 s34+s123^2 s16 s234 s34+mt2 s12 s23 s234 s34+mt2 s123 s23 s234 s34+2 s12 s123 s23 s234 s34+mt2 s123^2 s34^2-s123^2 s16 s34^2-mt2 s123 s23 s34^2+s123 s16 s23 s34^2-s12 s123 s23 s234 s345+s12 s123 s234^2 s345-s123^2 s234^2 s345+s123^2 s234 s34 s345-s123 s23 s234 s34 s345+s12^2 s23 s234 s45-s12^2 s234^2 s45+s12 s123 s234^2 s45-s12 s123 s234 s34 s45-s12 s23 s234 s34 s45-mt2 s12 s23^2 s56+2 mt2 s12 s23 s234 s56-2 mt2 s123 s23 s234 s56-4 mt2 s12 s23 s34 s56+2 mt2 s123 s23 s34 s56+2 s12 s16 s23 s34 s56-s123 s16 s23 s34 s56-mt2 s23^2 s34 s56-2 s12 s23^2 s34 s56+s12 s23^2 s345 s56-s12 s23 s234 s345 s56+2 s123 s23 s234 s345 s56-s123 s23 s34 s345 s56+s23^2 s34 s345 s56-s12 s23 s234 s45 s56+2 s12 s23 s34 s45 s56+mt2 s23^2 s56^2-s23^2 s345 s56^2)/.SliceTpoint,Modulus->$P]


dlogList[[evenpos]]/.Log[x_]:>Factor[x,Modulus->$P]
NewDlogs3[[1]]/.SliceTpoint/.Log[x_]:>Factor[x,Modulus->$P]





(* ::Section:: *)
(*Integrate*)


(*Get["/home/daniel/Documents/packages/trager/Trager.wl"]*)
Get["/home/dschiebe/Documents/packages/trager/Trager.wl"]
Dlog=TragerLogSubtraction[ReconstructedFunc,var,Parallelize->True,Debug->True]
Factor[Cases[{Dlog},_Log,Infinity]/.Log[x_]:>x]


(* ::Section:: *)
(*Comments:*)


(* ::Subsubsection:: *)
(*Current Highest Reco Done: 10-4/5 (94-94s23)*)


(* ::Section:: *)
(*Problems*)


Fac[65]/(Fac[3] Fac[139])
EvenLetterCandidatesFACS[[3]]
EvenLetterCandidatesFACS[[65]]
EvenLetterCandidatesFACS[[139]]


Test=Factor[dlogList[[-2]]/.Log[x_]:>x,Modulus->$P]


TA=Factor[EvenLetterCandidatesFACS[[3]]/.SliceTpoint,Modulus->$P]
TB=Factor[EvenLetterCandidatesFACS[[65]]/.SliceTpoint,Modulus->$P]
TC=Factor[EvenLetterCandidatesFACS[[139]]/.SliceTpoint,Modulus->$P]



Expand[Test,Modulus->$P]
TEST2=Expand[Test*TA TC/TB,Modulus->$P]
TEST2=Factor[Test*TA TC/TB,Modulus->$P]


Factor[D[Log[TEST2],t],Modulus->$P]


(642363787 (831395188+t) (78407223+2042854634 t+1486215360 t^2+t^3))/((753518641+t) (847954056+t) (963332568+t) (1621757903+t) (1772277964+t))
(642363787 (831395188+t) (78407223+2042854634 t+1486215360 t^2+t^3))/((753518641+t) (847954056+t) (963332568+t) (1621757903+t) (1772277964+t))
