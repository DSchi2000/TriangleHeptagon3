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


D4Repls={p[5]->c1 p[1]+c2 p[2]+c3 p[3]+ p[4]};


vecs=Table[p[i],{i,1,5}]
TestMAt=Expand[Expand[Table[vec1 vec2,{vec1,vecs},{vec2,vecs}]/.D4Repls]/.mandelstamScalarProducts]


vecs=Table[p[i],{i,1,4}]
TestMAt=Expand[Expand[Table[vec1 vec2,{vec1,vecs},{vec2,vecs}]]/.mandelstamScalarProducts]


Together[Det[2*TestMAt]]


EQN1=s45-(Expand[(p[4]+p[5])^2/.D4Repls]/.mandelstamScalarProducts)
EQN2=s345-(Expand[(p[3]+p[4]+p[5])^2/.D4Repls]/.mandelstamScalarProducts)
EQN3=s16-(Expand[(p[1]+p[6])^2/.p4cons/.D4Repls]/.mandelstamScalarProducts)


Solve[{EQN1==0,EQN2==0,EQN3==0},{c1,c2,c3}]//FullSimplify
