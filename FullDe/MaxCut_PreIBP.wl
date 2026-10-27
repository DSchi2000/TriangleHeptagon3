(* ::Package:: *)

SetDirectory[ParentDirectory[NotebookDirectory[]]];

If[$FrontEnd===Null,$DIR=Directory[],$DIR=NotebookDirectory[]];
parts = FileNameSplit[$DIR];
pos = First@FirstPosition[parts, "TriangleHeptagon3"];
HomeDir = FileNameJoin[parts[[;;pos ]]];
SectorDir = FileNameJoin[parts[[;;pos+1]]];
CurrentDir=$DIR;
kirapath=SectorDir;

Get[FileNameJoin[{Directory[],"files/HexaBoxModules.m"}]];
Get[FileNameJoin[{Directory[],"files/DifferentialOperator.m"}]];
(*sqrtSet=ReadList[FileNameJoin[{NotebookDirectory[],"Roots.txt"}]];*)
sqrtSet={}
$P=2^31-1;
PureIntegrals=ReadList[FileNameJoin[{NotebookDirectory[],"SubBasis.txt"}]];
SetDirectory[NotebookDirectory[]];

SectorLaporta=ConstantArray[0,95];
(*MaxCutInts={141,142,143};
*)
MaxCutInts=Range[95];
(*SectorOnline= HPT3[1,1,1,1,1,1,1,1,0,0,0,0,0];*)
AllSectors=True;

OffCutTest[integral_]:=Module[{sector},
sector=integral/.{-1->0,-2->0,-3->0,2->1,3->1};
If[sector===SectorOnline,Return[integral]];
Return[0];
]

(*MaxCutInts=Range[35];*)
PureIntegralsMaxCut=SectorLaporta;
PureIntegralsMaxCut[[MaxCutInts]]=PureIntegrals[[MaxCutInts]];

sqrtSet=(Cases[PureIntegralsMaxCut,_root,Infinity]//DeleteDuplicates)/.root[x_]:>x;

If[!AllSectors,
PureIntegralsMaxCut[[MaxCutInts]]=PureIntegralsMaxCut[[MaxCutInts]]/.HPT3[n___]:>OffCutTest[HPT3[n]];
]
changeAuxVar[S16];
tPoint=GetRandomPT[0,$P]

epsrepl={(ep->(4-d)/2)/.tPoint};
PureIntegralN=Expand[PureIntegralsMaxCut/.tPoint/.epsrepl//.root[n_]:>PowerMod[n,1/2,$P],Modulus->$P]//EchoTiming;
diffOp=getDifferentialOperator[tPoint];

loopMomentumScalarProductsT=Expand[loopMomentumScalarProducts/.tPoint,Modulus->$P];
mandelstamScalarProductsT=Expand[mandelstamScalarProducts/.tPoint,Modulus->$P];
DeleteList=Delete[{1,2,3,4,5,6,7,8,9,10,11},Position[tPoint[[All,1]],InternalAuxvar]];

Monitor[diffInts=Table[
ComputeDIFast[PureIntegralsMaxCut[[MaxCutInts[[intpos]]]],diffOp,var,tPoint],{var,1,Length[invariantsnodUpperCase]},{intpos,1,Length[MaxCutInts]}],Row[{ProgressIndicator[(var-1)*Length[MaxCutInts]+intpos,{1,Length[MaxCutInts]*10}],"  ",var,"/10 ",intpos,"/",Length[MaxCutInts]}]]//EchoTiming;

If[!AllSectors,
diffInts=diffInts/.HPT3[n___]:>OffCutTest[HPT3[n]];
];

diffIntsMC=Join[ConstantArray[0,{10,Length[PureIntegralsMaxCut]-Length[MaxCutInts]}],diffInts,2];
ToReduce=Cases[diffIntsMC,_HPT3,Infinity]//DeleteDuplicates//EchoFunction[Length];
A=If[MemberQ[(DeleteDuplicates[ToReduce/.n_Integer/;n<0:>0/.n_Integer/;n>1:>1]),HPT3[0,1,1,1,0,1,1,1,0,0,0,0,0]],
ToReduce=Join[ToReduce,{HPT3[-2,1,1,1,0,1,1,1,0,0,-2,0,0]}];
];

Export[FileNameJoin[{NotebookDirectory[],"KIRA/ReduList.txt"}],ToReduce];
GenerateSectorsFromReductionList[FileNameJoin[{NotebookDirectory[],"KIRA/ReduList.txt"}]];

ToReduce2=Cases[PureIntegralN,_HPT3,Infinity]//DeleteDuplicates//EchoFunction[Length];
Export[FileNameJoin[{NotebookDirectory[],"KIRA/UReduList.txt"}],ToReduce2];
GenerateSectorsFromReductionList[FileNameJoin[{NotebookDirectory[],"KIRA/UReduList.txt"}]];

GenerateSectorsFromReductionList[FileNameJoin[{NotebookDirectory[],"KIRA/ExtraRelationsReduList.txt"}]];

DeleteOldKiraRun[]
(*WriteNumericsFile[tPoint]*)
WriteJobFile[{"UReduList","ReduList","ExtraRelationsReduList"}]
(*WriteJobFile[{"UReduList","ReduList"}]*)

DeleteOldKiraRun[];
tpointlist={};
(*tpoint=GetRandomPT[0,$P];*)
Monitor[For[Npos=1,Npos<=20,Npos++,
tPoint[[1,2]]=RandomInteger[$P];
tpointlist=Join[tpointlist,{tPoint[[All,2]]}];
If[Npos==1,
WriteNumericsFile[tPoint];,
AddToNumericsFile[tPoint]
];
];
,ProgressIndicator[Npos,{1,20}]]








tPoint=GetRandomPT[0,$P]


rtset=DeleteDuplicates[Cases[PureIntegrals,_root,Infinity]]/.root[x_]:>x;



rtN=Expand[rtset/.tPoint,Modulus->$P]
TestA=Table[Expand[rt1/rt2,Modulus->$P],{rt1,rtN},{rt2,rtN}];


sqrtSet={}
tPoint=GetRandomPT[0,$P]
rtN=Expand[rtset/.tPoint,Modulus->$P]
TestB=Table[Expand[rt1/rt2,Modulus->$P],{rt1,rtN},{rt2,rtN}];


TestA/TestB//MatrixForm


rtset


PointGeneratorPath="/home/dschiebe/Documents/PhD-Files/FullHexaBoxDE"<>"/SQRT-Generator";

SetDirectory[CurrentDir]
Command=PointGeneratorPath<>"/gen_6pt --roots "<>SectorDir<>"/Roots.txt --gram "<>PointGeneratorPath<>"/SixPointGram.txt --npoints 1 --auxvar "<>ToString[s16]<>" --out points.txt --threads "<>ToString[1];
Run[Command];
PTS=Import[CurrentDir<>"/points.txt","Table"];
PTS = Select[PTS, !StringStartsQ[ToString[#[[1]]], "#"] &][[All,;;-1]];
invariants={d,mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56}
PTSRPL=Map[Thread[Rule[invariants,#]]&,PTS]


rtset=DeleteDuplicates[Cases[PureIntegrals,_root,Infinity]]/.root[x_]:>x;
(*sqrtSet=rtset[[{1,2,3,4,5,6,7,8,9,10,11,12,13,15}]];
sqrtSet=rtset[[{13,14}]];
tPoint=GetRandomPT[0,$P]*)


rtset[[-3]]
rtset[[-4]]

