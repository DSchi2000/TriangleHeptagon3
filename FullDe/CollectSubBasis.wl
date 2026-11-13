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
Get[FileNameJoin[{Directory[],"files/HexaBoxModules.m"}]];
HXBISPs=({k[1]+p[2],k[1]+p[3],k[1]+p[4],k[2]+p[6]}/.p4cons)^2
HXBISPs=Expand[Expand[HXBISPs]/.loopMomentumScalarProducts/.kinematics/.mandelstamScalarProducts]
Get[FileNameJoin[{Directory[],"files/DifferentialOperator.m"}]];






(* ::Section:: *)
(*Fix 5 Point Roots*)


CanoGramVec={p[1],p[2],p[3],p[4]};
Momenta={p[1],p[2],p[3],p[4],p[5],p[6]};
FivePMomenta=Subsets[Momenta,{4}];

TwoLineGrams=Table[Together[16 Det[Table[Expand[vec1 vec2/.p4cons]/.kinematics/.mandelstamScalarProducts,{vec1,momentas},{vec2,CanoGramVec}]]],{momentas,FivePMomenta}];
Grams=Table[Together[16 Det[Table[Expand[vec1 vec2/.p4cons]/.kinematics/.mandelstamScalarProducts,{vec1,momentas},{vec2,momentas}]]],{momentas,FivePMomenta}];
Gramrepls=Table[Rule[root[Grams[[i]]],TwoLineGrams[[i]]/root[Grams[[1]]]],{i,1,Length[Grams]}];


Tpoints={d->2,mt2->1867986130,s12->400105616,s123->1077685476,s16->1404758163,s23->1154837077,s234->1712394270,s34->934781953,s345->1748259880,s45->661755111,s56->1091626529}
RationalReconstruct[a_]:=With[{v=LatticeReduce[{{a,1},{$P,0}}][[1]]},v[[1]]/v[[2]]];
Gramrepls=Table[Rule[root[Grams[[i]]],TwoLineGrams[[i]]/root[Grams[[1]]]],{i,1,Length[Grams]}];
TESTVEC=Expand[Gramrepls/.Tpoints/.root[x_]:>PowerMod[x,1/2,$P],Modulus->$P]
BranchChoiceVec=Expand[TESTVEC[[All,1]]/TESTVEC[[All,2]],Modulus->$P]/.n_Integer:>RationalReconstruct[n]
Gramrepls[[All,2]]=Gramrepls[[All,2]]*BranchChoiceVec;


(* ::Section::Closed:: *)
(*85-86"*)


GramVec1={k[1],p[1],p[2]+p[3],p[4],p[5]};
GramVec2={k[2],p[1],p[2]+p[3],p[4],p[5]};
mu12=Collect[Expand[Det[Table[Expand[(GramVec1[[i]] GramVec2[[j]])/.p4cons],{i,1,5},{j,1,5}]/.loopMomentumScalarProducts/.mandelstamScalarProducts/.kinematics]],PR[___],Together]/Expand[2Det[Table[Expand[(GramVec1[[i]] GramVec2[[j]])/.p4cons],{i,2,5},{j,2,5}]/.mandelstamScalarProducts/.kinematics]];
mu11=Collect[Expand[Det[Table[Expand[(GramVec1[[i]] GramVec1[[j]])/.p4cons],{i,1,5},{j,1,5}]/.loopMomentumScalarProducts/.mandelstamScalarProducts/.kinematics]],PR[___],Together]/Expand[2Det[Table[Expand[(GramVec1[[i]] GramVec1[[j]])/.p4cons],{i,2,5},{j,2,5}]/.mandelstamScalarProducts/.kinematics]];
mu22=Collect[Expand[Det[Table[Expand[(GramVec2[[i]] GramVec2[[j]])/.p4cons],{i,1,5},{j,1,5}]/.loopMomentumScalarProducts/.mandelstamScalarProducts/.kinematics]],PR[___],Together]/Expand[2Det[Table[Expand[(GramVec2[[i]] GramVec2[[j]])/.p4cons],{i,2,5},{j,2,5}]/.mandelstamScalarProducts/.kinematics]];
Variables[mu12]//Sort
Variables[mu11]//Sort
Variables[mu22]//Sort
tr5=root[Expand[Det[2*Table[Expand[(GramVec1[[i]] GramVec2[[j]])/.p4cons],{i,2,5},{j,2,5}]/.mandelstamScalarProducts/.kinematics]]];

TestBasis= ep^3 tr5 Collect[ Expand[Together[mu22] HPT3[0, 1, 2, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0]],{_root,_HPT3},Factor];


(* ::Section::Closed:: *)
(*2,5,6,12,13*)


prefac=ep^2 (1-2 ep)(1-3 ep)
Int[2]=prefac*PureIntegralsOld[[2]]
Int[5]=prefac*PureIntegralsOld[[5]]
Int[6]=prefac*PureIntegralsOld[[6]]
Int[12]=prefac*PureIntegralsOld[[12]]
Int[13]=prefac*PureIntegralsOld[[13]]

(*TestBasis= ep^3 tr5 Collect[ Expand[Together[mu22] HPT3[0, 1, 2, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0]],{_root,_HPT3},Factor];*)


Together[mu22];





test=Cases[FixedBasis,_root,Infinity]/.root[x_]:>x;


Position[FixedBasis,HPT3[0,1,2,0,0,1,1,1,1,0,0,0,0]]


SubPureInts[[86]]


plotDiagram[HPT3[1,0,2,1,0,1,1,1,0,0,0,0,0]]


test=Cases[SubPureInts[[87]],_HPT3,Infinity]//DeleteDuplicates
test=Cases[SubPureInts[[88]],_HPT3,Infinity]//DeleteDuplicates


(* ::Section:: *)
(*Build SubBasis*)


(*FixedBasis=ReadList[CurrentDir<>"/CanonicalizedToDo.wl"];
FixedBasis=Collect[Expand[FixedBasis],{_root,_HPT3},Factor]
Export[CurrentDir<>"/BasisFixed.txt",FixedBasis]

FixedBasis=ReadList[CurrentDir<>"/BasisFixed.txt"];*)


SetDirectory[ParentDirectory[NotebookDirectory[]]];
PureIntegralsOld=ReadList[FileNameJoin[{Directory[],"files/SortetBasisFull.txt"}]];
PureIntegralsNew=PureIntegralsOld;
Get[FileNameJoin[{Directory[],"files/HexaBoxModules.m"}]];
$P=2^31-1;
SelectInts=Range[95];
SelectInts//Length
SubPureInts=PureIntegralsNew[[SelectInts]];
SubPureInts=SubPureInts/.Gramrepls;
(*SubPureInts[[MissingPos]]=FixedBasis;

SubPureInts=SubPureInts/.root[x_]:>root[Factor[x]];
SubPureInts=SubPureInts/.root[x_^2]:>x;
SubPureInts=SubPureInts/.root[x_]:>root[Expand[x]];
SubPureInts=SubPureInts/.signrepls;

SubPureInts[[2]]=Int[2];
SubPureInts[[5]]=Int[5];
SubPureInts[[6]]=Int[6];
SubPureInts[[12]]=Int[12];
SubPureInts[[13]]=Int[13];*)

PureIntegralsOld/.{n_Integer/;n>1:>1,n_Integer/;n<0:>0};
Sectors=Map[Last[Sort[Cases[{#},_HPT3,Infinity]/.{n_Integer/;n>1:>1,n_Integer/;n<0:>0}]]&,SubPureInts];
SectorsOld=Map[Last[Sort[Cases[{#},_HPT3,Infinity]/.{n_Integer/;n>1:>1,n_Integer/;n<0:>0}]]&,PureIntegralsOld[[SelectInts]]];
Permutation=FindPermutation[Sectors,SectorsOld];
SubPureInts=Permute[SubPureInts,Permutation];

(*SubPureInts[[85]]=TestBasis;*)

PureIntegralsNew[[SelectInts]]=SubPureInts;
Newsqrts=(Cases[SubPureInts,root[___],Infinity]//DeleteDuplicates)/.root[x_]:>x;
Newsqrts//Length

Oldsqrts=(Cases[PureIntegralsOld,root[___],Infinity]//DeleteDuplicates)/.root[x_]:>x;
sqrtSetNew=Join[Oldsqrts,Newsqrts]//DeleteDuplicates;
Export[FileNameJoin[{NotebookDirectory[],"Roots.txt"}],Newsqrts];
Export[FileNameJoin[{NotebookDirectory[],"SubBasis.txt"}],SubPureInts];


(* ::Subsection::Closed:: *)
(*Roots Fixing*)


FlipPos=Flatten[Position[sqrtSetNew[[All,1]]/.Thread[Rule[invariants,1]],-1]]
signrepls=Table[root[sqrtSetNew[[i]]]->root[-sqrtSetNew[[i]]],{i,FlipPos}];




Cases[SubPureInts[[-11]],_root,Infinity]/.signrepls


sqrtSetNew


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


SubPureInts[[79]]
SubPureInts[[80]]

