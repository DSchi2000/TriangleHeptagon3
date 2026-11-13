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
PointGeneratorPath=HomeDir<>"/SQRT-Generator";


(* ::Section::Closed:: *)
(*UniVariate-Slice*)


NParrallel=10
NSamples=2;
PointGeneratorPath=HomeDir<>"/SQRT-Generator";
Command=PointGeneratorPath<>"/gen_param --roots "<>CurrentDir<>"rootsPara.txt --npoints "<>ToString[NSamples]<>" --out points.txt --threads "<>ToString[NParrallel];
SetDirectory[CurrentDir];
Run[Command]
PTS=Import[CurrentDir<>"/points.txt","Table"];
PTS = Select[PTS, !StringStartsQ[ToString[#[[1]]], "#"] &][[All,;;-1]];

PTSRPL=Map[Thread[Rule[paravars,#]]&,PTS];

avec=PTSRPL[[1]];
bvec=PTSRPL[[2]];
Export["AVECS.txt",{avec,bvec}]

Export["slice_a.txt", {StringRiffle[ToString /@ avec[[All,2]], " "]}, "Lines"];
Export["slice_b.txt", {StringRiffle[ToString /@ bvec[[All,2]], " "]}, "Lines"];


NParrallel=10
NSamples=50;
Command=PointGeneratorPath<>"/gen_param --roots "<>CurrentDir<>"rootsPara.txt --npoints "<>ToString[NSamples]<>" --slice slice_a.txt slice_b.txt"<>" --out points.txt --threads "<>ToString[NParrallel];
Run[Command];

PTS=Import[CurrentDir<>"/points.txt","Table"];
PTS = Select[PTS, !StringStartsQ[ToString[#[[1]]], "#"] &][[All,;;-1]];

PTSRPL=Map[Thread[Rule[paravars,#]]&,PTS[[All,;;-2]]];
tvars=Map[Thread[Rule[{tvar},#]]&,PTS[[All,-1]]];
OutPutPath="Samples";
(*tpointlist=ParallelTable[ GetTouple[EpsRecoLenth,PTSRPL[[i]]],{i,1,NSamples}];*)
Export[ToString[StringForm["``/``/tvars.txt",CurrentDir,OutPutPath]],tvars];


mandpoints=Map[ConvertPointToMands[#]&,PTSRPL];
mandpoints=Map[Join[{d->RandomInteger[$P]},#]&,mandpoints];


G6=Get[HomeDir<>"/SetUp/SixPointGram.txt"];
Expand[G6/.mandpoints,Modulus->$P]


Export[ToString[StringForm["``/``/CustomPoints.txt",CurrentDir,OutPutPath]],mandpoints];





mandpoints[[1]]


(* ::Section::Closed:: *)
(*Bivariate Slice*)


NParrallel=10
NSamples=2;
PointGeneratorPath=HomeDir<>"/SQRT-Generator";
Command=PointGeneratorPath<>"/gen_param --roots "<>CurrentDir<>"rootsPara.txt --npoints "<>ToString[NSamples]<>" --out points.txt --threads "<>ToString[NParrallel];
SetDirectory[CurrentDir];
Run[Command]
PTS=Import[CurrentDir<>"/points.txt","Table"];
PTS = Select[PTS, !StringStartsQ[ToString[#[[1]]], "#"] &][[All,;;-1]];

PTSRPL=Map[Thread[Rule[paravars,#]]&,PTS];

avec=PTSRPL[[1]];
bvec=PTSRPL[[2]];
bvec[[Range[6],2]]=ConstantArray[0,6];
cvec=PTSRPL[[2]];
cvec[[Range[3]+6,2]]=ConstantArray[0,3];
Export["BiVariateSliceParameters.txt",{avec,bvec}]

Export["slice_a.txt", {StringRiffle[ToString /@ avec[[All,2]], " "]}, "Lines"];
Export["slice_b.txt", {StringRiffle[ToString /@ bvec[[All,2]], " "]}, "Lines"];
Export["slice_c.txt", {StringRiffle[ToString /@ cvec[[All,2]], " "]}, "Lines"];


NParrallel=10
NSamples=400;
Command=PointGeneratorPath<>"/gen_param --roots "<>CurrentDir<>"rootsPara.txt --npoints "<>ToString[NSamples]<>" --slice2 slice_a.txt slice_b.txt slice_c.txt"<>" --out points.txt --threads "<>ToString[NParrallel];
Run[Command];

PTS=Import[CurrentDir<>"/points.txt","Table"];
PTS = Select[PTS, !StringStartsQ[ToString[#[[1]]], "#"] &][[All,;;-1]];

PTSRPL=Map[Thread[Rule[paravars,#]]&,PTS[[All,;;-3]]];
tvars=Map[Thread[Rule[{xvar,yvar},#]]&,PTS[[All,-2;;]]];
OutPutPath="Samples";
(*tpointlist=ParallelTable[ GetTouple[EpsRecoLenth,PTSRPL[[i]]],{i,1,NSamples}];*)
Export[ToString[StringForm["``/``/tvars.txt",CurrentDir,OutPutPath]],tvars];


mandpoints=Map[ConvertPointToMands[#]&,PTSRPL];
mandpoints=Map[Join[{d->RandomInteger[$P]},#]&,mandpoints];
G6=Get[HomeDir<>"/SetUp/SixPointGram.txt"];
Expand[G6/.mandpoints,Modulus->$P]


Export[ToString[StringForm["``/``/CustomPoints.txt",CurrentDir,OutPutPath]],mandpoints];





mandpoints[[1]]
mandpoints[[2]]



(* ::Section:: *)
(*FullSamples*)


NParrallel=10
NSamples=10000;
SetDirectory[CurrentDir];
Command=PointGeneratorPath<>"/gen_param --roots "<>CurrentDir<>"rootsPara.txt --npoints "<>ToString[NSamples]<>" --out points.txt --threads "<>ToString[NParrallel];
Run[Command]

PTS=Import[CurrentDir<>"/points.txt","Table"];
PTS = Select[PTS, !StringStartsQ[ToString[#[[1]]], "#"] &][[All,;;-1]];
PTSRPL=Map[Thread[Rule[paravars,#]]&,PTS[[All,;;-1]]];
(*tvars=Map[Thread[Rule[{xvar,yvar},#]]&,PTS[[All,-2;;]]];
OutPutPath="Samples";
(*tpointlist=ParallelTable[ GetTouple[EpsRecoLenth,PTSRPL[[i]]],{i,1,NSamples}];*)
Export[ToString[StringForm["``/``/tvars.txt",CurrentDir,OutPutPath]],tvars];*)


mandpoints=Map[ConvertPointToMands[#]&,PTSRPL];
mandpoints=Map[Join[{d->RandomInteger[$P]},#]&,mandpoints];
G6=Get[HomeDir<>"/SetUp/SixPointGram.txt"];
Expand[G6/.mandpoints,Modulus->$P]//DeleteDuplicates


OutPutPath="Samples";
Export[ToString[StringForm["``/``/CustomPoints.txt",CurrentDir,OutPutPath]],mandpoints];
