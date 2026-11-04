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


SolvedAlready=ReadList[CurrentDir<>"Solved-2.txt"];
FactorsWithPositions=ReadList[CurrentDir<>"EvenLetterPos1-2.txt"];
FactorsWithPositions2=ReadList[CurrentDir<>"EvenLetterPos2-2.txt"];
EvenLetterCandidatesFACS=ReadList[CurrentDir<>"Factors-2.txt"];
DtIntegrands=ReadList[CurrentDir<>"EvenIntegrands.txt"];
dlogList=ReadList[CurrentDir<>"dlogs.txt"];
VarGrid=ReadList[CurrentDir<>"VarGrid.txt"];
Nterms=ReadList[CurrentDir<>"TermsNeededPowers-2.txt"];
AnsatzPowers=ReadList[CurrentDir<>"RemainingPowers-2.txt"];
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
LaunchKernels[5];


pos1=95;
pos2=95;
pos=GetLinPos[pos1,pos2,95];
DIsij=Table[Insert[Map[LinearTerm[#][[psp,1,pos]]&,usedvars],0,4],{psp,1,Length[PhaseSpacePoints]}]//EchoTiming;
(*DIsij=Table[Coefficient[ChangeAuxVarInDe[s16,s234,DIsij[[psp]],PhaseSpacePoints[[psp]]],Table[ds[var],{var,diffInvariants}]],{psp,1,Length[PhaseSpacePoints]}]//EchoTiming;*)
setup=AnsatzPowers[[Flatten[Position[Nterms[[All,1]],{pos1,pos2}]]//First]]
Nterms[[Flatten[Position[Nterms[[All,1]],{pos1,pos2}]]//First]]
vars=VarGrid[[setup[[1,1]],setup[[1,2]]]]//Sort
(*vars={s16,mt2}*)


var=mt2;
setup[[2]]
setup[[3]]

fitvec=DIsij[[All,First[Flatten[Position[diffInvariants,var]]]]];

testpos=Flatten[Position[FactorsWithPositions[[All,1]],setup[[1]]]]//First
FactorsWithPositions2[[testpos]]
evenpos=FactorsWithPositions2[[testpos]][[1]];


DtIntegrands[[evenpos]]
Expand[DtIntegrands[[evenpos]],Modulus->$P]
dlogList[[evenpos]]
logRational=dlogList[[evenpos]]/.n_ Log[x_]:>n/.Log[x_]:>1
logRational Log[FactorsWithPositions2[[testpos]][[2]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]]
DlogDiff=Factor[Expand[logRational  D[Log[FactorsWithPositions2[[testpos]][[2]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]],var]],Modulus->$P]
fitvec=Expand[fitvec-(DlogDiff/.PhaseSpacePoints),Modulus->$P];
Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
Reco=ReconstructFUNCS[fitvec,PhaseSpacePoints[[All,Flatten[Map[Position[diffInvariants,#]&,vars]],2]],vars,5,6,Debug->True];
Ansatz=Reco[[2]];
ReconstructedFunc=Together[Ansatz/.Thread[Rule[Sort[Cases[Ansatz,_fac,Infinity]],Reco[[1,1]]]]/.n_Integer:>RationalReconstruct[n]]


(* ::Section:: *)
(*Integrate*)


(*Get["/home/daniel/Documents/packages/trager/Trager.wl"]*)
Get["/home/dschiebe/Documents/packages/trager/Trager.wl"]
Dlog=TragerLogSubtraction[ReconstructedFunc,var,Parallelize->True,Debug->True]
Factor[Cases[{Dlog},_Log,Infinity]/.Log[x_]:>x]


(* ::Section:: *)
(*Automation*)


CloseKernels[];
$DistributedContexts = None;
kernels=LaunchKernels[20];


Get[ToString[StringJoin[HomeDir,"/files/ExtraRecoModules.m"]]];
Get["/home/dschiebe/Documents/packages/trager/Trager.wl"]
GetNumberOfTermsNeeded[n_,num_,den_]:=Binomial[n + num-1, n-1]+Binomial[n + den-1, n-1]


GetDLogs[pos1_,pos2_]:=Module[{pos,DIsij,setup,vars,newauxvar,fitvec,testpos,evenpos,logRational,DlogDiff,Reco,Ansatz,ReconstructedFunc,Dlog,Dlogs},
pos=GetLinPos[pos1,pos2,95];

setup=AnsatzPowers[[Flatten[Position[Nterms[[All,1]],{pos1,pos2}]]//First]];

If[setup[[2]]-setup[[3]]===-2, setup[[2]]=setup[[2]]+1];
vars=VarGrid[[setup[[1,1]],setup[[1,2]]]]//Sort;

If[ Length[vars]===10,newauxvar=s16;,
newauxvar=Last[Complement[diffInvariants,vars]];
];

DIsij=Table[Insert[Map[LinearTerm[#][[psp,1,pos]]&,usedvars],0,4],{psp,1,Length[PhaseSpacePoints]}];
DIsij=Table[Coefficient[ChangeAuxVarInDe[s16,newauxvar,DIsij[[psp]],PhaseSpacePoints[[psp]]],Table[ds[var],{var,diffInvariants}]],{psp,1,Length[PhaseSpacePoints]}];


Dlogs=Table[
fitvec=DIsij[[All,First[Flatten[Position[diffInvariants,dvar]]]]];

testpos=Flatten[Position[FactorsWithPositions[[All,1]],setup[[1]]]]//First;
evenpos=FactorsWithPositions2[[testpos]][[1]];



logRational=dlogList[[evenpos]]/.n_ Log[x_]:>n/.Log[x_]:>1;
DlogDiff=Factor[Expand[logRational  D[Log[FactorsWithPositions2[[testpos]][[2]]/.Fac[n_]:>EvenLetterCandidatesFACS[[n]]],dvar]],Modulus->$P];
fitvec=Expand[fitvec-(DlogDiff/.PhaseSpacePoints),Modulus->$P];

Reco=ReconstructFUNCS[fitvec,PhaseSpacePoints[[All,Flatten[Map[Position[diffInvariants,#]&,vars]],2]],vars,setup[[2]],setup[[3]],Debug->False];
Ansatz=Reco[[2]];
ReconstructedFunc=Together[Ansatz/.Thread[Rule[Sort[Cases[Ansatz,_fac,Infinity]],Reco[[1,1]]]]/.n_Integer:>RationalReconstruct[n]];
(*ReconstructedFunc//Echo;*)
TragerLogSubtraction[ReconstructedFunc,dvar,Parallelize->True,Debug->False],{dvar,vars}];
Return[Dlogs/.Log[x_]:>Log[Factor[x]]];
]



Monitor[NewDLogs=Table[
matpos=AnsatzPowers[[it1,1]];
setup=AnsatzPowers[[Flatten[Position[Nterms[[All,1]],{matpos[[1]],matpos[[2]]}]]//First]];
vars=VarGrid[[matpos[[1]],matpos[[2]]]]//Sort;
If[GetNumberOfTermsNeeded[Length[vars],setup[[2]],setup[[3]]]<=500,
GetDLogs[matpos[[1]],matpos[[2]]],0],
{it1,1,Length[AnsatzPowers[[All,1]]]}],{ProgressIndicator[it1,{1,Length[AnsatzPowers[[All,1]]]}],AnsatzPowers[[it1,1]]}];
Export[CurrentDir<>"NewDlogs-500-2.txt",NewDLogs]


(* ::Section:: *)
(*Problems*)


DeleteDuplicates[Position[NewDLogs,$Failed][[All,1]]]
Export[CurrentDir<>"NewDlogs-500.txt",NewDLogs]


failedEntries=DeleteDuplicates[Position[NewDLogs,$Failed][[All,1]]]


NewDLogs[[70]]


AnsatzPowers[[222]]


Flatten[Position[AnsatzPowers[[All,2]]-AnsatzPowers[[All,3]],-2]]



GetDLogs[80,1]


Monitor[NewDLogsFailed=Table[
matpos=AnsatzPowers[[it1,1]];
setup=AnsatzPowers[[Flatten[Position[Nterms[[All,1]],{matpos[[1]],matpos[[2]]}]]//First]];
vars=VarGrid[[matpos[[1]],matpos[[2]]]]//Sort;
If[GetNumberOfTermsNeeded[Length[vars],setup[[2]],setup[[3]]]<=500,
GetDLogs[matpos[[1]],matpos[[2]]],0],
{it1,failedEntries}],{ProgressIndicator[it1,{1,Length[AnsatzPowers[[All,1]]]}],AnsatzPowers[[it1,1]]}];


NewDLogsFailed


Expand[D[Log[(906811812+1417006012 t+872189447 t^2+386053936 t^3+1484787178 t^4+40849427 t^5+t^6)/(57055600+196815237 t+720038392 t^2+1147743761 t^3+13266943 t^4+1085638404 t^5+266323194 t^6+1836178903 t^7+t^8)],t],Modulus->$P]
