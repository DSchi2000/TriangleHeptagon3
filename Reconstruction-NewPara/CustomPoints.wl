(* ::Package:: *)

If[$FrontEnd===Null,$DIR=Directory[],$DIR=NotebookDirectory[]];
parts = FileNameSplit[$DIR];
pos = First@FirstPosition[parts, "TriangleHeptagon3"];
HomeDir = FileNameJoin[parts[[;;pos ]]];
SectorDir = FileNameJoin[parts[[;;pos+1]]];
CurrentDir=$DIR;
kirapath=SectorDir;
$P=2^31-1;
Get[HomeDir<>"/4DBasis/Parametrization.m"]
diffInvariants={mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56};
roots=ReadList[SectorDir<>"/Roots.txt"];
rootsPara=Factor[roots/.KinematicRulesBackWards]
Export[CurrentDir<>"rootsPara.txt",rootsPara]


(* ::Section:: *)
(*Gen a b*)


NParrallel=10
NSamples=1
PointGeneratorPath=HomeDir<>"/SQRT-Generator";
Command=PointGeneratorPath<>"/gen_param --roots "<>CurrentDir<>"rootsPara.txt --npoints "<>ToString[NSamples]<>" --slice slice_a.txt slice_b.txt"<>" --out points.txt --threads "<>ToString[NParrallel];
