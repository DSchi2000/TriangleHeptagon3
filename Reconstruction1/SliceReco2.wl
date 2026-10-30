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




(* ::Section:: *)
(*Load More Suff For Dlogging:*)


SolvedAlready=ReadList[CurrentDir<>"Solved.txt"];
FactorsWithPositions=ReadList[CurrentDir<>"EvenLetterPos.txt"];
FactorsWithPositions2=ReadList[CurrentDir<>"EvenLetterPos2.txt"];
EvenLetterCandidatesFACS=ReadList[CurrentDir<>"Factors.txt"];
DtIntegrands=ReadList[CurrentDir<>"EvenIntegrands.txt"];
dlogList=ReadList[CurrentDir<>"dlogs.txt"];



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
(*GetVariablesOfIntegral*)


ClearAll[HPT3Momenta];
HPT3Momenta[HPT3[ns__]] := HPT3Momenta[{ns}]
HPT3Momenta[powers_List] := Module[
  {edges, legOf, verts, parent, find, unite, groups, legsOfGroup},
  edges = {
    1 -> {"B", "X1"}, 2 -> {"X1", "T"}, 3 -> {"B", "T"},
    4 -> {"T", "Y1"}, 5 -> {"Y1", "Y2"}, 6 -> {"Y2", "Y3"},
    7 -> {"Y3", "Y4"}, 8 -> {"Y4", "Y5"}, 9 -> {"Y5", "B"}
  };
  legOf = <|"X1" -> p[1], "Y1" -> p[2], "Y2" -> p[3],
            "Y3" -> p[4], "Y4" -> p[5], "Y5" -> p[6]|>;
  verts = {"B", "X1", "T", "Y1", "Y2", "Y3", "Y4", "Y5"};

  parent = AssociationThread[verts -> verts];
  find[v_] := If[parent[v] === v, v, parent[v] = find[parent[v]]];
  unite[a_, b_] := With[{ra = find[a], rb = find[b]},
    If[ra =!= rb, parent[ra] = rb]];
  (* power <= 0 means the propagator is absent (pinched, or a pure
     numerator dot) -- either way its two endpoints merge *)
  Do[
    If[powers[[e]] <= 0, unite @@ (e /. edges)],
    {e, 1, 9}
  ];
  groups = GatherBy[verts, find];
  legsOfGroup[grp_] := Total[Lookup[legOf, grp, Nothing]];
  Select[legsOfGroup /@ groups, # =!= 0 &]
]


Sectors=Map[Last[Sort[DeleteDuplicates[Cases[#,_HPT3,Infinity]/.{n_Integer/;n<0:>0,n_Integer/;n>1:>1}]]]&,PureIntegrals];
MomentaConfigurations=Map[HPT3Momenta[#]&,Sectors];
VariableIntegrals=Map[Variables[Det[Table[2Expand[ vec1 vec2/.p4cons]/.mandelstamScalarProducts/.kinematics,{vec1,#[[;;-2]]},{vec2,#[[;;-2]]}]]]&,MomentaConfigurations]
VarGrid=Table[DeleteDuplicates[Join[vars1,vars2]],{vars1,VariableIntegrals},{vars2,VariableIntegrals}];


|


(* ::Section:: *)
(*Get Remaining Powers:*)


GetNumberOfTermsNeeded[n_,num_,den_]:=Binomial[n + num, n]+Binomial[n + den, n]
RecoPowers=Map[GetResidualPowers[#]&,Complement[FactorsWithPositions[[All,1]],SolvedAlready]]

NAnsaetze=Table[{RecoPowers[[i,1]],GetNumberOfTermsNeeded[Length[VarGrid[[RecoPowers[[i,1,1]],RecoPowers[[i,1,2]]]]],RecoPowers[[i,2]],RecoPowers[[i,3]]]},{i,1,Length[RecoPowers]}]


(* ::Section:: *)
(*BackUp*)


testpos=Flatten[Position[FactorsWithPositions[[All,1]],{91,82}]]//First
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


