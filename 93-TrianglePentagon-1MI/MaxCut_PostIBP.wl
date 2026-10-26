(* ::Package:: *)

SetDirectory[ParentDirectory[NotebookDirectory[]]];

If[$FrontEnd===Null,$DIR=Directory[],$DIR=NotebookDirectory[]];
parts = FileNameSplit[$DIR];
pos = First@FirstPosition[parts, "TriangleHeptagon3"];
HomeDir = FileNameJoin[parts[[;;pos ]]];
SectorDir = FileNameJoin[parts[[;;pos+1]]];
CurrentDir=$DIR;
kirapath=SectorDir;

DoExtraRelations=True;

Get[FileNameJoin[{Directory[],"files/HexaBoxModules.m"}]];
Get[FileNameJoin[{Directory[],"files/DifferentialOperator.m"}]];
(*sqrtSet=ReadList[FileNameJoin[{NotebookDirectory[],"Roots.txt"}]];*)
sqrtSet={}
$P=2^31-1;
PureIntegrals=ReadList[FileNameJoin[{NotebookDirectory[],"SubBasis.txt"}]];
SectorOnline= HPT3[1,1,1,1,0,1,1,1,1,0,0,0,0];
(*PureIntegrals=PureIntegrals/.root[mt2^2 s12^2-2 mt2 s12^2 s16+s12^2 s16^2+2 mt2^2 s12 s234-2 mt2 s12^2 s234-2 mt2 s12 s16 s234-2 s12^2 s16 s234+mt2^2 s234^2-2 mt2 s12 s234^2+s12^2 s234^2-2 mt2^2 s12 s34+4 mt2 s12 s16 s34-2 s12 s16^2 s34-2 mt2^2 s234 s34+2 mt2 s12 s234 s34+2 mt2 s16 s234 s34+2 s12 s16 s234 s34+mt2^2 s34^2-2 mt2 s16 s34^2+s16^2 s34^2-2 mt2 s12 s234 s345+2 s12 s16 s234 s345-2 mt2 s234^2 s345-2 s12 s234^2 s345+2 mt2 s234 s34 s345-2 s16 s234 s34 s345+s234^2 s345^2+2 mt2 s12 s16 s56-2 s12 s16^2 s56-2 mt2 s16 s234 s56+2 s12 s16 s234 s56+2 mt2 s16 s34 s56-4 s12 s16 s34 s56-2 s16^2 s34 s56-2 mt2 s12 s345 s56+2 s12 s16 s345 s56+2 mt2 s234 s345 s56+2 s12 s234 s345 s56+2 s16 s234 s345 s56-2 mt2 s34 s345 s56+2 s16 s34 s345 s56-2 s234 s345^2 s56+s16^2 s56^2-2 s16 s345 s56^2+s345^2 s56^2]->-root[mt2^2 s12^2-2 mt2 s12^2 s16+s12^2 s16^2+2 mt2^2 s12 s234-2 mt2 s12^2 s234-2 mt2 s12 s16 s234-2 s12^2 s16 s234+mt2^2 s234^2-2 mt2 s12 s234^2+s12^2 s234^2-2 mt2^2 s12 s34+4 mt2 s12 s16 s34-2 s12 s16^2 s34-2 mt2^2 s234 s34+2 mt2 s12 s234 s34+2 mt2 s16 s234 s34+2 s12 s16 s234 s34+mt2^2 s34^2-2 mt2 s16 s34^2+s16^2 s34^2-2 mt2 s12 s234 s345+2 s12 s16 s234 s345-2 mt2 s234^2 s345-2 s12 s234^2 s345+2 mt2 s234 s34 s345-2 s16 s234 s34 s345+s234^2 s345^2+2 mt2 s12 s16 s56-2 s12 s16^2 s56-2 mt2 s16 s234 s56+2 s12 s16 s234 s56+2 mt2 s16 s34 s56-4 s12 s16 s34 s56-2 s16^2 s34 s56-2 mt2 s12 s345 s56+2 s12 s16 s345 s56+2 mt2 s234 s345 s56+2 s12 s234 s345 s56+2 s16 s234 s345 s56-2 mt2 s34 s345 s56+2 s16 s34 s345 s56-2 s234 s345^2 s56+s16^2 s56^2-2 s16 s345 s56^2+s345^2 s56^2];*)
SectorLaporta=ReadList[FileNameJoin[{SectorDir,"SectorLaportaBasis.txt"}]];
(*MaxCutInts={141,142,143};*)
MaxCutInts=Range[53];
AllSectors=True;
(*MaxCutInts=Range[35];*)
usedauxvar=S12;
changeAuxVar[usedauxvar];

PureIntegralsMaxCut=SectorLaporta;
PureIntegralsMaxCut[[MaxCutInts]]=PureIntegrals[[MaxCutInts]];
If[!AllSectors,
PureIntegralsMaxCut[[MaxCutInts]]=PureIntegralsMaxCut[[MaxCutInts]]/.HPT3[n___]:>OffCutTest[HPT3[n]];
]

SetDirectory[NotebookDirectory[]];
Get[FileNameJoin[{HomeDir,"files/ExtraFunctions.m"}]];
NotbookDir=NotebookDirectory[];
SetDirectory[HomeDir];
CloseKernels[];
kernels=LaunchKernels[20];
Print["Loading Modules"]
(*ParallelEvaluate[Get[FileNameJoin[{Directory[],"files/HexaBoxModules.m"}]],kernels];
ParallelEvaluate[Get[FileNameJoin[{Directory[],"files/DifferentialOperator.m"}]],kernels];
ParallelEvaluate[Get[FileNameJoin[{HomeDir,"files/ExtraFunctions.m"}]],kernels];*)

ParallelEvaluate[changeAuxVar[usedauxvar],kernels];
ParallelEvaluate[SetDirectory[NotbookDir];,kernels];
Print["Loading BASIS"]
(*ParallelEvaluate[PureIntegrals=ReadList[FileNameJoin[{Directory[],"SubBasis.txt"}]];,kernels];*)
Print["MC BASIS"]
(*ParallelEvaluate[PureIntegralsMaxCut=SectorLaporta; PureIntegralsMaxCut[[MaxCutInts]]=PureIntegrals[[MaxCutInts]];,kernels];*)
Print["Loading Reduction Table"]
ParallelEvaluate[Ureduce=ReadList[FileNameJoin[{Directory[],"KIRA/UReduList.txt"}]];,kernels];

Print["Loading Roots"]
ParallelEvaluate[sqrtSet=ReadList[FileNameJoin[{Directory[],"Roots.txt"}]],kernels];
Print["Loading DE Code"]
(*ParallelEvaluate[Get[FileNameJoin[{HomeDir,"files/ExtraFunctions.m"}]],kernels];*)
SetDirectory[NotebookDirectory[]];

parallelOut=ParallelTable[ComputeDIMC[i],{i,1,20}]//EchoTiming;

dHPT3Matrices=parallelOut[[All,2]];
tpoints=parallelOut[[All,1]];
OutPutPath="ParralelD10";
SetOptions[OpenAppend,PageWidth->Infinity];
DeleteDirectory[FileNameJoin[{NotebookDirectory[],OutPutPath}],DeleteContents->True];
CreateDirectory[FileNameJoin[{NotebookDirectory[],OutPutPath}]];
For[q=1,q<=Length[dHPT3Matrices],q++,
PutAppend[tpoints[[q]],FileNameJoin[{NotebookDirectory[],ToString[StringForm["``/parameters",OutPutPath]]}]];
];
For[i=1,i<=Length[invariantsnod],i++,
For[q=1,q<=Length[dHPT3Matrices],q++,
PutAppend[dHPT3Matrices[[q,i]],FileNameJoin[{NotebookDirectory[],ToString[StringForm["``/DEmatrix_d``",OutPutPath,invariants[[i+1]]]]}]];
];
]

OutPutPath="ParralelD10";
SetDirectory[NotebookDirectory[]];
$P=2^31-1;
diffInvariants={mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56};
ClearAll[dHPT3Matrices];
For[i=1,i<=10,i++,			
dHPT3Matrices[diffInvariants[[i]]]=ReadList[FileNameJoin[{Directory[],ToString[StringForm["``/DEmatrix_d``",OutPutPath,diffInvariants[[i]]]]}]];
]
phasespacepoints=ReadList[FileNameJoin[{Directory[],ToString[StringForm["``/parameters",OutPutPath]]}]];
RationalReconstruct[a_]:=With[{v=LatticeReduce[{{a,1},{$P,0}}][[1]]},v[[1]]/v[[2]]];

BasisLength=Length[MaxCutInts];

colors={-1->White,0->Red,1->Green,2->Blue,3->Yellow,4->Orange,5->RGBColor["#e09f7d"],6->RGBColor["#EF5D60"],7->RGBColor["#EC4067"],8->RGBColor["#A01A7D"],9->RGBColor["#703399"],10->RGBColor["#0B050F"]};
CutsA=1;
CutsB=BasisLength;
For[i=1,i<=10,i++,
MatReco[diffInvariants[[i]]]=ParallelTable[Quiet[Expand[powerReconstruct[phasespacepoints[[All,1]],dHPT3Matrices[diffInvariants[[i]]][[All,pos]]]/.{t->4-2\[Epsilon]},Modulus->$P]],{pos,1,BasisLength^2}];
(*NumiPowers=Exponent[Numerator[MatReco[diffInvariants[[i]]]],\[Epsilon]]/.{-Infinity->-1};
DenoPowers=Exponent[Denominator[MatReco[diffInvariants[[i]]]],\[Epsilon]];
NumiPowersMat=NumiPowers//Normal;
DenoPowersMat=DenoPowers//Normal;
P1=MatrixPlot[ArrayReshape[NumiPowersMat,{BasisLength,BasisLength}][[CutsA;;CutsB,CutsA;;CutsB]],ColorRules->colors,PlotLabel->"Numerator",ImageSize->Medium];
P2=MatrixPlot[ArrayReshape[DenoPowersMat,{BasisLength,BasisLength}][[CutsA;;CutsB,CutsA;;CutsB]],ColorRules->colors,PlotLabel->"Denominator",ImageSize->Medium];
L1=SwatchLegend[Values@colors[[2;;]],Keys@colors[[2;;]],LegendFunction->"Frame",LegendLayout->"Column",LegendLabel->diffInvariants[[i]]];*)
];
Factor[ArrayReshape[MatReco[mt2],{BasisLength,BasisLength}],Modulus->$P]//MatrixForm
Factor[ArrayReshape[MatReco[s12],{BasisLength,BasisLength}],Modulus->$P]//MatrixForm
Factor[ArrayReshape[MatReco[s123],{BasisLength,BasisLength}],Modulus->$P]//MatrixForm
Factor[ArrayReshape[MatReco[s16],{BasisLength,BasisLength}],Modulus->$P]//MatrixForm
Factor[ArrayReshape[MatReco[s23],{BasisLength,BasisLength}],Modulus->$P]//MatrixForm
Factor[ArrayReshape[MatReco[s234],{BasisLength,BasisLength}],Modulus->$P]//MatrixForm
Factor[ArrayReshape[MatReco[s34],{BasisLength,BasisLength}],Modulus->$P]//MatrixForm
Factor[ArrayReshape[MatReco[s345],{BasisLength,BasisLength}],Modulus->$P]//MatrixForm
Factor[ArrayReshape[MatReco[s45],{BasisLength,BasisLength}],Modulus->$P]//MatrixForm
Factor[ArrayReshape[MatReco[s56],{BasisLength,BasisLength}],Modulus->$P]//MatrixForm







(Expand[Sum[Expand[phasespacepoints[[6, q + 1]] ArrayReshape[MatReco[diffInvariants[[q]]], {BasisLength, BasisLength}], Modulus -> $P], {q, 1, 10}], Modulus -> $P])/.n_Integer:>RationalReconstruct[n] // MatrixForm





(* ::Section:: *)
(*Intense Debug trash*)


Get[FileNameJoin[{HomeDir,"files/ExtraFunctions.m"}]];
SectorLaporta=ReadList[FileNameJoin[{SectorDir,"SectorLaportaBasis.txt"}]];
PureIntegralsMaxCut=SectorLaporta;
PureIntegralsMaxCut[[MaxCutInts]]=PureIntegrals[[MaxCutInts]];
Ureduce=ReadList[FileNameJoin[{Directory[],"KIRA/UReduList.txt"}]];
SUBLAPORTA=ComputeDIMC[1];
Position[SUBLAPORTA,HPT3[-1,1,1,1,0,1,1,1,0,0,-2,0,0]]
(*Cases[PureIntegrals[[-3;;]],_HPT3,Infinity]//DeleteDuplicates//Length
Ureduce//Length*)


SUBLAPORTA


PureIntegrals[[1]]
