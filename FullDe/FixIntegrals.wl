(* ::Package:: *)

(* ::Section:: *)
(*Load Packages*)


If[$FrontEnd===Null,$DIR=Directory[],$DIR=NotebookDirectory[]];
parts = FileNameSplit[$DIR];
pos = First@FirstPosition[parts, "TriangleHeptagon3"];
HomeDir = FileNameJoin[parts[[;;pos ]]];
SectorDir = FileNameJoin[parts[[;;pos+1]]];
CurrentDir=$DIR;
SetDirectory[ParentDirectory[NotebookDirectory[]]];
Get[FileNameJoin[{Directory[],"files/HexaBoxModules.m"}]];
Get[FileNameJoin[{Directory[],"files/DifferentialOperator.m"}]];


(* ::Section:: *)
(*Compute Mus*)


GetMuIJ[vecs_,i_,j_]:=Collect[Expand[Det[Table[Expand[(vec1 vec2)/.p4cons],{vec1,Join[{k[i]},vecs]},{vec2,Join[{k[j]},vecs]}]/.loopMomentumScalarProducts/.mandelstamScalarProducts/.kinematics]],PR[___],Together]/Expand[2Det[Table[Expand[(vec1 vec2)/.p4cons],{vec1,vecs},{vec2,vecs}]/.mandelstamScalarProducts/.kinematics]];
GetGram[vecs_]:=Together[Expand[Det[2*Table[Expand[(vec1 vec2)/.p4cons],{vec1,vecs},{vec2,vecs}]/.mandelstamScalarProducts/.kinematics]]];


Together[GetGram[{p[1],p[2],p[3],p[6]}]];
Together[GetMuIJ[{p[1],p[2],p[3],p[6]},1,2]];


(* ::Section:: *)
(*Compute LS*)


Get1LLeadingSing[vecs_,intMass_]:=Module[{ModCayley,Y,OneLoopLS},
SetAttributes[Y,Orderless];
Table[Y[i,i]=2intMass[[i]] ,{i,1,Length[intMass]}];
Table[Y[i,j]=intMass[[i]]+intMass[[j]]-Sum[vecs[[q]],{q,i,j-1}]^2 ,{i,1,Length[intMass]},{j,i+1,Length[intMass]}];
Y[0,0]=0;
Table[Y[0,i]=1,{i,1,Length[intMass]}];
If[EvenQ[Length[vecs]],
ModCayley=Expand[Expand[Table[Y[i,j],{i,1,Length[intMass]},{j,1,Length[intMass]}]/.p4cons]/.kinematics/.mandelstamScalarProducts],
ModCayley=Expand[Expand[Table[Y[i,j],{i,0,Length[intMass]},{j,0,Length[intMass]}]/.p4cons]/.kinematics/.mandelstamScalarProducts]
];
OneLoopLS=root[Factor[Det[ModCayley]]]/.{root[x_^2]:>x,root[-x_^2]:>x};
Return[OneLoopLS]
]


Get1LLeadingSing[{p[1]+p[2],p[3]+p[4],p[5]+p[6]},{0,0,mt2}]


FixedBasis=ReadList[CurrentDir<>"/CanonicalizedToDo.wl"];
FixedBasis=Collect[Expand[FixedBasis],{_root,_HPT3},Factor]
Export[CurrentDir<>"/BasisFixed.txt",FixedBasis]




