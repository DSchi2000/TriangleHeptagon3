(* ::Package:: *)

(* ::Section:: *)
(*Find Basis*)


If[$FrontEnd===Null,$DIR=Directory[],$DIR=NotebookDirectory[]];
parts = FileNameSplit[$DIR];
pos = First@FirstPosition[parts, "TriangleHeptagon3"];
HomeDir = FileNameJoin[parts[[;;pos ]]];
SectorDir = FileNameJoin[parts[[;;pos+1]]];
CurrentDir=$DIR;

SetDirectory[ParentDirectory[NotebookDirectory[]]];
Get[FileNameJoin[{Directory[],"files/HexaBoxModules.m"}]];
Get[FileNameJoin[{Directory[],"files/DifferentialOperator.m"}]];
ExtraRelationsSolset={HPT3[1,0,1,1,1,1,1,1,0,0,0,0,0],HPT3[0,1,1,0,1,1,1,1,1,0,0,0,0],HPT3[1,1,1,0,1,1,1,1,0,0,0,0,0]};
$P=2^31-1;
PureIntegrals=ReadList[FileNameJoin[{Directory[],"files/SortetBasisFull.txt"}]]/.root[x_]:>1;
Ureduce=ReadList[FileNameJoin[{Directory[],"FirstDeFull/KIRA/UReduList.txt"}]];
SetDirectory[NotebookDirectory[]];
samplepoint=1;
SubReduList=Get[FileNameJoin[{ParentDirectory[Directory[]],"FirstDeFull/KIRA/results/HPT3/numerics_2147483647_0.m"}]];
(*SubReduList=Get[ToString[StringForm["/home/daniel/lanayru/PhD_Files/FullHexaBoxDE/FullDE/KIRA/results/HPT3/numerics_2147483647_0.m",IntegerString[samplepoint,10,5]]]];*)
tPoint=Thread[Rule[invariants,SubReduList[[2,1]]]]
reductionRules=SubReduList[[2,2]];
epsrepl={(ep->(4-d)/2)/.tPoint};
PureIntegralN=Expand[PureIntegrals/.tPoint/.epsrepl/.root[n_]:>PowerMod[n,1/2,$P],Modulus->$P]//EchoTiming;
ExtraRelations=GetExtraRelationReplacements[reductionRules,tPoint];
PureIntegralN2=Expand[PureIntegralN/.reductionRules,Modulus->$P];
PureIntegralN3=Expand[PureIntegralN2/.ExtraRelations,Modulus->$P];

NonZeros=Cases[reductionRules,_HPT3,Infinity]//DeleteDuplicates;
ZeroInts=Complement[Ureduce,NonZeros];
ZeroIntsRepl=Thread[Rule[ZeroInts,0]];
PureIntegralN4=Expand[PureIntegralN3/.ExtraRelations/.ZeroIntsRepl,Modulus->$P];
Cases[PureIntegralN4,_HPT3,Infinity]//DeleteDuplicates//Length
SubLaporta=Cases[PureIntegralN4,_HPT3,Infinity]//DeleteDuplicates;

PureIntegralN4=Expand[PureIntegralN4/.epsrepl,Modulus->$P];

(*BasisChangMat=Table[CoefficientArrays[PureIntegralN4[[pos]],SubLaporta],{pos,1,Length[PureIntegralN]}];
BasisChangMat=Normal[BasisChangMat[[All,2]]];
BasisChangMatInv=Inverse[BasisChangMat,Modulus->$P];*)

RemReduList=Get[FileNameJoin[{NotebookDirectory[],"KIRATestFolder/results/HPT3/numerics_2147483647_0.m"}]][[2,2]];
RemLaportaBasis=Cases[(Cases[RemReduList[[All,2]],_HPT3,Infinity]//DeleteDuplicates)/.ExtraRelations,_HPT3,Infinity]//DeleteDuplicates
RemLaportaBasis//Length

Sectors=(RemLaportaBasis/.n_Integer/;n>1:>1/.n_Integer/;n<0:>0)//DeleteDuplicates;
FullSectors=Map[Last[Sort[DeleteDuplicates[Cases[{#},_HPT3,Infinity]/.n_Integer/;n>1:>1/.n_Integer/;n<0:>0]]]&,PureIntegralN4];
Sort[Map[Position[FullSectors,#]&,Sectors]//Flatten]//DeleteDuplicates//EchoFunction[Length]

Export[FileNameJoin[{NotebookDirectory[],"SectorLaportaBasis.txt"}],RemLaportaBasis];


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

GramVec1={k[1],p[1],p[2],p[3],p[4]};
GramVec2={k[2],p[1],p[2],p[3],p[4]};
mu12=Collect[Expand[Det[Table[Expand[(GramVec1[[i]] GramVec2[[j]])/.p4cons],{i,1,5},{j,1,5}]/.loopMomentumScalarProducts/.mandelstamScalarProducts/.kinematics]],PR[___],Together]/Expand[2Det[Table[Expand[(GramVec1[[i]] GramVec2[[j]])/.p4cons],{i,2,5},{j,2,5}]/.mandelstamScalarProducts/.kinematics]];
mu11=Collect[Expand[Det[Table[Expand[(GramVec1[[i]] GramVec1[[j]])/.p4cons],{i,1,5},{j,1,5}]/.loopMomentumScalarProducts/.mandelstamScalarProducts/.kinematics]],PR[___],Together]/Expand[2Det[Table[Expand[(GramVec1[[i]] GramVec1[[j]])/.p4cons],{i,2,5},{j,2,5}]/.mandelstamScalarProducts/.kinematics]];
mu22=Collect[Expand[Det[Table[Expand[(GramVec2[[i]] GramVec2[[j]])/.p4cons],{i,1,5},{j,1,5}]/.loopMomentumScalarProducts/.mandelstamScalarProducts/.kinematics]],PR[___],Together]/Expand[2Det[Table[Expand[(GramVec2[[i]] GramVec2[[j]])/.p4cons],{i,2,5},{j,2,5}]/.mandelstamScalarProducts/.kinematics]];
Variables[mu12]//Sort
Variables[mu11]//Sort
Variables[mu22]//Sort
tr5=root[Expand[16Det[Table[Expand[(GramVec1[[i]] GramVec2[[j]])/.p4cons],{i,2,5},{j,2,5}]/.mandelstamScalarProducts/.kinematics]]];





HXBBasis=ReadList[FileNameJoin[{NotebookDirectory[],"HXB-Basis.txt"}]];
DeleteDuplicates[Cases[HXBBasis,_hxb,Infinity]];
HXBProps={k[1]^2,(k[1]+p[1])^2,(k[1]+p[1]+p[2])^2,(-k[1]+k[2])^2,(k[2]+p[1]+p[2])^2,(k[2]+p[1]+p[2]+p[3])^2,(k[2]+p[1]+p[2]+p[3]+p[4])^2,-mt2+(k[2]+p[1]+p[2]+p[3]+p[4]+p[5])^2,k[2]^2,(k[1]+p[3])^2,(k[1]+p[4])^2,(k[1]+p[5])^2,(k[2]+p[1])^2};
HXBPropsKShift=HXBProps(*/.{k[1]->k[1]-p[1],k[2]->k[2]-p[1]}*)
HXBPropsPShift=HXBPropsKShift(*/.{p[n_]:>p[n-1]}/.p[0]->p[6]*)
HXBISPs=Expand[Expand[HXBPropsPShift/.p4cons]/.loopMomentumScalarProducts/.kinematics/.mandelstamScalarProducts]

HXBTOHPT3[hxbINT_]:=Expand[Times@@Power[HXBISPs,-hxbINT/.hxb->List]HPT3[0,0,0,0,0,0,0,0,0,0,0,0,0]]
(*mandchange={s12->s16,s123->s345,s16->s56,s23->s12,s234->s123,s34->s23,s345->s234,s45->s34,s56->s45}*)

HPT3BASIS=HXBBasis/.hxb[x___]:>HXBTOHPT3[hxb[x]];
(*Etra relations in HXB integral family need to remove 5*)
HPT3BASIS[[5]]=HPT3[1,0,1,0,1,0,0,0,0,0,0,0,0];


HPT3BASIS


(* ::Section:: *)
(*Build SubBasis*)


SetDirectory[ParentDirectory[NotebookDirectory[]]];
PureIntegralsOld=ReadList[FileNameJoin[{Directory[],"files/SortetBasisFull.txt"}]];
PureIntegralsNew=PureIntegralsOld;
Get[FileNameJoin[{Directory[],"files/HexaBoxModules.m"}]];
$P=2^31-1;
SelectInts={1,3,4,7,10,11,14,15,16,17,18,22,25,26,27,28,33,34,35,37,38,39,40,41,45,46,47,48,49,50,55,56,57,58,63,64,65,66,71,72,73,74,75,76,81,82,83,84,89,90,91,92,93};
SelectInts//Length
SubPureInts=PureIntegralsNew[[SelectInts]];
SubPureInts=HPT3BASIS;
SubPureInts=SubPureInts/.root[x_]:>root[Factor[x]];
SubPureInts=SubPureInts/.root[x_^2]:>x;
SubPureInts=SubPureInts/.root[x_]:>root[Expand[x]];

PureIntegralsOld/.{n_Integer/;n>1:>1,n_Integer/;n<0:>0};
Sectors=Map[Last[Sort[Cases[{#},_HPT3,Infinity]/.{n_Integer/;n>1:>1,n_Integer/;n<0:>0}]]&,SubPureInts];
SectorsOld=Map[Last[Sort[Cases[{#},_HPT3,Infinity]/.{n_Integer/;n>1:>1,n_Integer/;n<0:>0}]]&,PureIntegralsOld[[SelectInts]]];
Permutation=FindPermutation[Sectors,SectorsOld];
SubPureInts=Permute[SubPureInts,Permutation];
PureIntegralsNew[[SelectInts]]=SubPureInts;
Newsqrts=(Cases[SubPureInts,root[___],Infinity]//DeleteDuplicates)/.root[x_]:>x;
Oldsqrts=(Cases[PureIntegralsOld,root[___],Infinity]//DeleteDuplicates)/.root[x_]:>x;
sqrtSetNew=Join[Oldsqrts,Newsqrts]//DeleteDuplicates;
Export[FileNameJoin[{NotebookDirectory[],"Roots.txt"}],Newsqrts];
Export[FileNameJoin[{NotebookDirectory[],"SubBasis.txt"}],SubPureInts];



(* ::Section:: *)
(*Export Basis*)


Export[FileNameJoin[{ParentDirectory[NotebookDirectory[]],"/files/SortetBasisFull.txt"}],PureIntegralsNew];


(* ::Section:: *)
(*PlotDiagram*)


plotDiagram[HPT3[1,1,1,0,1,1,0,1,0,0,0,0,0]]
(*BasisChangeInt=Collect[HXBBasis[[18]]/.mandchange/.hxbtohtp1,{_root,_HPT3}]//Factor*)








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
