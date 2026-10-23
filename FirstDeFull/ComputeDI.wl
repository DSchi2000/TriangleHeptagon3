(* ::Package:: *)

(* ::Section:: *)
(*Load Basis*)


If[$FrontEnd===Null,$DIR=Directory[],$DIR=NotebookDirectory[]];
parts = FileNameSplit[$DIR];
pos = First@FirstPosition[parts, "TriangleHeptagon3"];
HomeDir = FileNameJoin[parts[[;;pos ]]];
SectorDir = FileNameJoin[parts[[;;pos+1]]];
CurrentDir=$DIR;
kirapath=SectorDir;

Get[FileNameJoin[{HomeDir,"files/HexaBoxModules.m"}]];
Get[FileNameJoin[{HomeDir,"files/DifferentialOperator.m"}]];
$P=2^31-1;
PureIntegrals=ReadList[FileNameJoin[{HomeDir,"files/SortetBasisFull.txt"}]];
rootSet=(Cases[PureIntegrals,_root,Infinity]//DeleteDuplicates)/.root[x_]:>x;
sqrtSet=rootSet;


(*rootSet=sqrtSet;
sqrtSet={};
PhaseSpacePoints=Table[GetRandomPT[0,$P],{q,1,20}];*)
(*RTMAT=Expand[rootSet/.PhaseSpacePoints,Modulus->$P];
MatrixRank[RTMAT,Modulus->$P]
RowReduce[RTMAT,Modulus->$P]//MatrixPlot*)


(* ::Section::Closed:: *)
(*Fix 5pt Grams*)


GramVec={p[1],p[2],p[3],p[4]};
AllMomenta={p[1],p[2],p[3],p[4],p[5],p[6]};
GramVecs=Subsets[AllMomenta,{4}];
GramVecs=DeleteCases[GramVecs,GramVec];
CanoGram=Expand[16Det[Table[Expand[( GramVec[[i]] GramVec[[j]])/.p4cons],{i,1,4},{j,1,4}]/.mandelstamScalarProducts/.kinematics]];
OtherGrams=Table[Expand[16Det[Table[Expand[( GramVecs[[pos]][[i]] GramVecs[[pos]][[j]])/.p4cons],{i,1,4},{j,1,4}]/.mandelstamScalarProducts/.kinematics]],{pos,1,Length[GramVecs]}];
TwoLineGrams=Table[Expand[16Det[Table[Expand[( GramVec[[i]] GramVecs[[pos]][[j]])/.p4cons],{i,1,4},{j,1,4}]/.mandelstamScalarProducts/.kinematics]],{pos,1,Length[GramVecs]}];
Table[Position[rootSet,OtherGrams[[i]]],{i,1,14}]
rootSet//Length
G5Repls=Table[Rule[root[OtherGrams[[i]]],TwoLineGrams[[i]]/root[CanoGram]],{i,1,14}];
PureIntegrals=PureIntegrals/.G5Repls;
sqrtSet=(Cases[PureIntegrals,_root,Infinity]//DeleteDuplicates)/.root[x_]:>x;
sqrtSet//Length
sqrtSet


(* ::Section:: *)
(*Get Phasespace Point*)


Get[FileNameJoin[{HomeDir,"files/DifferentialOperator.m"}]];
SetDirectory[NotebookDirectory[]];
tPoint=GetRandomPT[0,$P]
epsrepl={(ep->(4-d)/2)/.tPoint};
diffOp=getDifferentialOperator[tPoint];
PureIntegralN=Expand[PureIntegrals/.tPoint/.epsrepl/.root[n_]:>PowerMod[n,1/2,$P],Modulus->$P]//EchoTiming;
checkDifferentialOperator[tPoint,diffOp];


(* ::Section:: *)
(*Compute DI Pre IBP*)


loopMomentumScalarProductsT=Expand[loopMomentumScalarProducts/.tPoint,Modulus->$P];
mandelstamScalarProductsT=Expand[mandelstamScalarProducts/.tPoint,Modulus->$P];
DeleteList=Delete[{1,2,3,4,5,6,7,8,9,10,11},Position[tPoint[[All,1]],s16]];
Monitor[diffInts=Table[ComputeDIFast[PureIntegrals[[intpos]],diffOp,var,tPoint],{var,1,Length[invariantsnodUpperCase]},{intpos,1,Length[PureIntegrals]}],Row[{ProgressIndicator[(var-1)*Length[PureIntegrals]+intpos,{1,Length[PureIntegrals]*10}],"  ",var,"/10 ",intpos,"/",Length[PureIntegrals]}]]//EchoTiming;


(* ::Section:: *)
(*KIRA Setup*)


ToReduce2


ToReduce=Cases[diffInts,_HPT3,Infinity]//DeleteDuplicates//EchoFunction[Length];
Export[FileNameJoin[{NotebookDirectory[],"KIRA/ReduList.txt"}],ToReduce];
GenerateSectorsFromReductionList[FileNameJoin[{NotebookDirectory[],"KIRA/ReduList.txt"}]];
ToReduce2=Cases[PureIntegralN,_HPT3,Infinity]//DeleteDuplicates//EchoFunction[Length];
Export[FileNameJoin[{NotebookDirectory[],"KIRA/UReduList.txt"}],ToReduce2];
GenerateSectorsFromReductionList[FileNameJoin[{NotebookDirectory[],"KIRA/UReduList.txt"}]];
GenerateSectorsFromReductionList[FileNameJoin[{NotebookDirectory[],"KIRA/ExtraRelationsReduList.txt"}]];
DeleteOldKiraRun[];
WriteNumericsFile[tPoint];
WriteJobFile[{"ExtraRelationsReduList","UReduList","ReduList"}]
(*WriteJobFile[{"ExtraRelationsReduList","UReduList"}];*)


(* ::Section:: *)
(*Rotation Matrices*)


reductionRules=Get[FileNameJoin[{Directory[],"KIRA/results/HPT3/numerics_2147483647_0.m"}]][[2,2]];
(*reductionRules=ReductionRules=GetReductionTableFromList[SubReduList[[2,2]]]//EchoTiming;*)
Ureduce=ReadList[FileNameJoin[{NotebookDirectory[],"KIRA/UReduList.txt"}]];
ExtraRelations=GetExtraRelationReplacements[reductionRules,tPoint]//EchoFunction[Length];
PureIntegralN2=Expand[PureIntegralN/.reductionRules,Modulus->$P];
PureIntegralN3=Expand[PureIntegralN2/.ExtraRelations,Modulus->$P];

NonZeros=Cases[reductionRules,_HPT3,Infinity]//DeleteDuplicates;
ZeroInts=Complement[Ureduce,NonZeros];
ZeroIntsRepl=Thread[Rule[ZeroInts,0]];
PureIntegralN4=Expand[PureIntegralN3/.ExtraRelations/.ZeroIntsRepl,Modulus->$P];
Cases[PureIntegralN4,_HPT3,Infinity]//DeleteDuplicates//Length
SubLaporta=Cases[PureIntegralN4,_HPT3,Infinity]//DeleteDuplicates;

PureIntegralN4=Expand[PureIntegralN4/.epsrepl,Modulus->$P];

BasisChangMat=Table[CoefficientArrays[PureIntegralN4[[pos]],SubLaporta],{pos,1,Length[PureIntegralN]}];
BasisChangMat=Normal[BasisChangMat[[All,2]]];
BasisChangMatInv=Inverse[BasisChangMat,Modulus->$P];
BasisChangMat//MatrixPlot
BasisChangMatInv//MatrixPlot






(* ::Section:: *)
(*Apply IBPs*)


diffIntsIBREDUCED=Expand[diffInts/.reductionRules/.ExtraRelations,Modulus->$P];//EchoTiming
NonReduced=Cases[diffIntsIBREDUCED,_hxb,Infinity]//DeleteDuplicates//EchoFunction[Length];
NullInts2=Complement[NonReduced,SubLaporta];
diffIntsIBREDUCED=diffIntsIBREDUCED/.Thread[Rule[NullInts2,0]];
dhxbMatrices=Table[Normal[Coefficient[diffIntsIBREDUCED[[i,j]],SubLaporta]],{i,1,Length[invariants]-1},{j,1,Length[SubLaporta]}];
Mplot1=MatrixPlot[Expand[dhxbMatrices[[1]] . BasisChangMatInv,Modulus->$P]];
Mplot2=MatrixPlot[Expand[dhxbMatrices[[2]] . BasisChangMatInv,Modulus->$P]];
Mplot3=MatrixPlot[Expand[dhxbMatrices[[3]] . BasisChangMatInv,Modulus->$P]];
Mplot4=MatrixPlot[Expand[dhxbMatrices[[4]] . BasisChangMatInv,Modulus->$P]];
Mplot5=MatrixPlot[Expand[dhxbMatrices[[5]] . BasisChangMatInv,Modulus->$P]];
Mplot6=MatrixPlot[Expand[dhxbMatrices[[6]] . BasisChangMatInv,Modulus->$P]];
Mplot7=MatrixPlot[Expand[dhxbMatrices[[7]] . BasisChangMatInv,Modulus->$P]];
Mplot8=MatrixPlot[Expand[dhxbMatrices[[8]] . BasisChangMatInv,Modulus->$P]];
Mplot9=MatrixPlot[Expand[dhxbMatrices[[9]] . BasisChangMatInv,Modulus->$P]];
Mplot10=MatrixPlot[Expand[dhxbMatrices[[10]] . BasisChangMatInv,Modulus->$P]];
GraphicsGrid[{{Mplot1,Mplot2,Mplot3,Mplot4,Mplot5},{Mplot6,Mplot7,Mplot8,Mplot9,Mplot10}}]





