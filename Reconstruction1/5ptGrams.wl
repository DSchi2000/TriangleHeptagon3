(* ::Package:: *)

If[$FrontEnd===Null,$DIR=Directory[],$DIR=NotebookDirectory[]];
parts = FileNameSplit[$DIR];
pos = First@FirstPosition[parts, "TriangleHeptagon3"];
HomeDir = FileNameJoin[parts[[;;pos ]]];
SectorDir = FileNameJoin[parts[[;;pos+1]]];
CurrentDir=$DIR;
SetDirectory[ParentDirectory[NotebookDirectory[]]];
Get[FileNameJoin[{Directory[],"files/HexaBoxModules.m"}]];
Get[FileNameJoin[{Directory[],"files/DifferentialOperator.m"}]];
SetDirectory[ParentDirectory[NotebookDirectory[]]];
Get[FileNameJoin[{Directory[],"files/HexaBoxModules.m"}]];
Get[FileNameJoin[{Directory[],"files/DifferentialOperator.m"}]];

CanoGramVec={p[1],p[2],p[3],p[4]};
Momenta={p[1],p[2],p[3],p[4],p[5],p[6]};
FivePMomenta=Subsets[Momenta,{4}];

TwoLineGrams=Table[Together[16 Det[Table[Expand[vec1 vec2/.p4cons]/.kinematics/.mandelstamScalarProducts,{vec1,momentas},{vec2,CanoGramVec}]]],{momentas,FivePMomenta}]
Grams=Table[Together[16 Det[Table[Expand[vec1 vec2/.p4cons]/.kinematics/.mandelstamScalarProducts,{vec1,momentas},{vec2,momentas}]]],{momentas,FivePMomenta}]


Export[CurrentDir<>"grams5pt.txt",Join[Grams,TwoLineGrams]]
