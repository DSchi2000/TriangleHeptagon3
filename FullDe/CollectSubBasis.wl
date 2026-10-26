(* ::Package:: *)

(* ::Section:: *)
(*Find Basis*)


If[$FrontEnd===Null,$DIR=Directory[],$DIR=NotebookDirectory[]];
parts = FileNameSplit[$DIR];
pos = First@FirstPosition[parts, "TriangleHeptagon3"];
HomeDir = FileNameJoin[parts[[;;pos ]]];
SectorDir = FileNameJoin[parts[[;;pos+1]]];
CurrentDir=$DIR;
PureIntegralsOld=ReadList[FileNameJoin[{HomeDir,"files/SortetBasisFull.txt"}]];
MissingPos=Flatten[Position[Map[Length[Variables[#]]&,PureIntegralsOld],1]]


MissingIntList=PureIntegralsOld[[MissingPos]];
Export[CurrentDir<>"/ToDO.txt",MissingIntList];


PureIntegralsOld//Length


(* ::Section:: *)
(*Extra Stuff*)


If[$FrontEnd===Null,$DIR=Directory[],$DIR=NotebookDirectory[]];
parts = FileNameSplit[$DIR];
pos = First@FirstPosition[parts, "TriangleHeptagon3"];
HomeDir = FileNameJoin[parts[[;;pos ]]];
SectorDir = FileNameJoin[parts[[;;pos+1]]];
CurrentDir=$DIR;
SetDirectory[ParentDirectory[NotebookDirectory[]]];
HXBISPs=({k[1]+p[2],k[1]+p[3],k[1]+p[4],k[2]+p[6]}/.p4cons)^2
HXBISPs=Expand[Expand[HXBISPs]/.loopMomentumScalarProducts/.kinematics/.mandelstamScalarProducts]
Get[FileNameJoin[{Directory[],"files/HexaBoxModules.m"}]];
Get[FileNameJoin[{Directory[],"files/DifferentialOperator.m"}]];






(* ::Section:: *)
(*85-86"*)


GramVec1={k[1],p[1],p[2],p[3],p[4]+p[5]};
GramVec2={k[2],p[1],p[2],p[3],p[4]+p[5]};
mu12=Collect[Expand[Det[Table[Expand[(GramVec1[[i]] GramVec2[[j]])/.p4cons],{i,1,5},{j,1,5}]/.loopMomentumScalarProducts/.mandelstamScalarProducts/.kinematics]],PR[___],Together]/Expand[2Det[Table[Expand[(GramVec1[[i]] GramVec2[[j]])/.p4cons],{i,2,5},{j,2,5}]/.mandelstamScalarProducts/.kinematics]];
mu11=Collect[Expand[Det[Table[Expand[(GramVec1[[i]] GramVec1[[j]])/.p4cons],{i,1,5},{j,1,5}]/.loopMomentumScalarProducts/.mandelstamScalarProducts/.kinematics]],PR[___],Together]/Expand[2Det[Table[Expand[(GramVec1[[i]] GramVec1[[j]])/.p4cons],{i,2,5},{j,2,5}]/.mandelstamScalarProducts/.kinematics]];
mu22=Collect[Expand[Det[Table[Expand[(GramVec2[[i]] GramVec2[[j]])/.p4cons],{i,1,5},{j,1,5}]/.loopMomentumScalarProducts/.mandelstamScalarProducts/.kinematics]],PR[___],Together]/Expand[2Det[Table[Expand[(GramVec2[[i]] GramVec2[[j]])/.p4cons],{i,2,5},{j,2,5}]/.mandelstamScalarProducts/.kinematics]];
Variables[mu12]//Sort
Variables[mu11]//Sort
Variables[mu22]//Sort
tr5=root[Expand[16Det[Table[Expand[(GramVec1[[i]] GramVec2[[j]])/.p4cons],{i,2,5},{j,2,5}]/.mandelstamScalarProducts/.kinematics]]];



FixedBasis=ReadList[CurrentDir<>"/BasisFixed.txt"];


test=Cases[FixedBasis,_root,Infinity]/.root[x_]:>x;


(* ::Section:: *)
(*Build SubBasis*)


SetDirectory[ParentDirectory[NotebookDirectory[]]];
PureIntegralsOld=ReadList[FileNameJoin[{Directory[],"files/SortetBasisFull.txt"}]];
PureIntegralsNew=PureIntegralsOld;
Get[FileNameJoin[{Directory[],"files/HexaBoxModules.m"}]];
$P=2^31-1;
SelectInts=Range[95];
SelectInts//Length
SubPureInts=PureIntegralsNew[[SelectInts]];
SubPureInts[[MissingPos]]=FixedBasis;

SubPureInts=SubPureInts/.root[x_]:>root[Factor[x]];
SubPureInts=SubPureInts/.root[x_^2]:>x;
SubPureInts=SubPureInts/.root[x_]:>root[Expand[x]];
SubPureInts=SubPureInts/.signrepls;

PureIntegralsOld/.{n_Integer/;n>1:>1,n_Integer/;n<0:>0};
Sectors=Map[Last[Sort[Cases[{#},_HPT3,Infinity]/.{n_Integer/;n>1:>1,n_Integer/;n<0:>0}]]&,SubPureInts];
SectorsOld=Map[Last[Sort[Cases[{#},_HPT3,Infinity]/.{n_Integer/;n>1:>1,n_Integer/;n<0:>0}]]&,PureIntegralsOld[[SelectInts]]];
Permutation=FindPermutation[Sectors,SectorsOld];
SubPureInts=Permute[SubPureInts,Permutation];

PureIntegralsNew[[SelectInts]]=SubPureInts;
Newsqrts=(Cases[SubPureInts,root[___],Infinity]//DeleteDuplicates)/.root[x_]:>x;
Newsqrts//Length

Oldsqrts=(Cases[PureIntegralsOld,root[___],Infinity]//DeleteDuplicates)/.root[x_]:>x;
sqrtSetNew=Join[Oldsqrts,Newsqrts]//DeleteDuplicates;
Export[FileNameJoin[{NotebookDirectory[],"Roots.txt"}],Newsqrts];
Export[FileNameJoin[{NotebookDirectory[],"SubBasis.txt"}],SubPureInts];


sqrtSetNew[[-1]]
sqrtSetNew[[-2]]
sqrtSetNew[[-3]]
sqrtSetNew[[-4]]
sqrtSetNew[[-5]]
sqrtSetNew[[-6]]
sqrtSetNew[[-7]]


signrepls=Table[root[sqrtSetNew[[i]]]->root[-sqrtSetNew[[i]]],{i,{-1,-2,-3,-4,-5,-6}}];




signrepls


Cases[SubPureInts[[-11]],_root,Infinity]/.signrepls


sqrtSetNew[[-3]]


(* ::Section:: *)
(*Export Basis*)


Export[FileNameJoin[{ParentDirectory[NotebookDirectory[]],"/files/SortetBasisFull.txt"}],PureIntegralsNew];


(* ::Section:: *)
(*PlotDiagram*)


plotDiagram[SubPureInts[[1]]]
(*BasisChangeInt=Collect[HXBBasis[[18]]/.mandchange/.hxbtohtp1,{_root,_HPT3}]//Factor*)


SubPureInts[[2]]//Factor






Table[
(Sort[Cases[{PureIntegralsOld[[SelectInts]][[i]]},_HPT3,Infinity]/.{n_Integer/;n>1:>1,n_Integer/;n<0:>0}]//Last)-
(Sort[Cases[SubPureInts[[i]],_HPT3,Infinity]/.{n_Integer/;n>1:>1,n_Integer/;n<0:>0}]//Last),{i,1,38}]


SelectInts
selectedints={1,3,4,7,10,11,12,13,16,17,18,19,20,24,25,26,27,32,33,34,35,40,42,43,44,45,48,49,50,51,52,56,57,58,59,64,65,66,70,71,72,73,74,75,80,81,82,83,84,85,90,91,92,93,94}


Length[Join[SelectInts,selectedints]//DeleteDuplicates]



Position[Map[Length[Variables[#]]&,PureIntegralsNew],1]//Length


Sectors=Map[Last[Sort[Cases[{#},_HPT3,Infinity]/.{n_Integer/;n>1:>1,n_Integer/;n<0:>0}]]&,SubPureInts];
SectorsOld=Map[Last[Sort[Cases[{#},_HPT3,Infinity]/.{n_Integer/;n>1:>1,n_Integer/;n<0:>0}]]&,PureIntegralsOld[[SelectInts]]];


Complement[SectorsOld,Sectors]


Position[Sectors,HPT3[0,0,1,0,0,0,0,0,1,0,0,0,0]]


SubPureInts[[5]]
HXBBasis[[5]]
DeleteDuplicates[Cases[HXBBasis,_hxb,Infinity]/.{n_Integer/;n>1:>1,n_Integer/;n<0:>0}]


TBASIS=ReadList["/home/dschiebe/Documents/PhD-Files/FullHexaBoxDE/files/PureBasisHexabox.txt"];



DeleteDuplicates[Cases[TBASIS[[206]],_hxb,Infinity]/.{n_Integer/;n>1:>1,n_Integer/;n<0:>0}]


HXBBasis[[5]]


SubPureInts[[1]]


Get1LLeadingSing[{p[1],p[2],p[3],p[4]+p[5],p[6]},{0,0,0,0,mt2}]


Get1LLeadingSing[{p[1]+p[2],p[3],p[4]+p[5],p[6]},{0,0,0,mt2}]
