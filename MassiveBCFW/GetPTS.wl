(* ::Package:: *)

Quiet[Get[NotebookDirectory[]<>"MassiveBCFWSlice.wl"]]


res=BuildMassiveSlice[RandomInteger[$P], False];
avec=Normal[res["a"]][[All;;10]]
bvec=Normal[res["b"]][[All;;10]]



G6=Get[ParentDirectory[NotebookDirectory[]]<>"/SetUp/SixPointGram.txt"];
newpoint=Expand[avec[[All,2]]+RandomInteger[$P]*bvec[[All,2]],Modulus->$P]
tpoint=Thread[Rule[ToExpression[avec[[All,1]]],newpoint]]
Solve[Expand[G6/.tpoint[[{1,2,3,5,6,7,8,9,10}]],Modulus->$P]==0,s16,Modulus->$P]
Expand[G6/.tpoint,Modulus->$P]


(* ::Section:: *)
(*Univariate*)


PointGeneratorPath=ParentDirectory[NotebookDirectory[]]<>"/SQRT-Generator";
SectorDir=ParentDirectory[NotebookDirectory[]]<>"/FullDe";
GramDir=ParentDirectory[NotebookDirectory[]]<>"/SetUp";

Export["slice_a.txt", {StringRiffle[ToString /@ Join[{RandomInteger[$P]},avec[[All,2]]], " "]}, "Lines"];
Export["slice_b.txt", {StringRiffle[ToString /@ Join[{RandomInteger[$P]},bvec[[All,2]]], " "]}, "Lines"];
NParrallel=10;
NSamples=10;

SetDirectory[NotebookDirectory[]]
(*Command=PointGeneratorPath<>"/gen_6pt --roots "<>SectorDir<>"/Roots.txt --gram "<>GramDir<>"/SixPointGram.txt --npoints 1 --auxvar "<>ToString[s16]<>" --out points.txt --threads "<>ToString[1];*)
Command=PointGeneratorPath<>"/gen_6pt --roots "<>SectorDir<>"/Roots.txt --gram "<>GramDir<>"/SixPointGram.txt --npoints "<>ToString[NSamples]<>" --auxvar "<>ToString[s16]<>" --slice slice_a.txt slice_b.txt"<>" --out points.txt --threads "<>ToString[NParrallel]<>" --gram-resolved";
Print[Command];

Run[Command];
PTS=Import[NotebookDirectory[]<>"/points.txt","Table"];
PTS = Select[PTS, !StringStartsQ[ToString[#[[1]]], "#"] &][[All,;;-1]];
invariants={d,mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56}
PTSRPL=Map[Thread[Rule[invariants,#]]&,PTS[[All,;;-2]]]
tvars=Map[Thread[Rule[{tvar},#]]&,PTS[[All,-1]]]
Expand[G6/.PTSRPL,Modulus->$P]
Expand[(PTSRPL[[2,2;;,2]]-avec[[All,2]])/bvec[[All,2]],Modulus->$P]


(* ::Section:: *)
(*Free Samples*)


Command=PointGeneratorPath<>"/gen_6pt --roots "<>SectorDir<>"/Roots.txt --gram "<>GramDir<>"/SixPointGram.txt --npoints "<>ToString[NSamples]<>" --auxvar "<>ToString[s16]<>" --out points.txt --threads "<>ToString[NParrallel];
Print[Command];

Run[Command];
PTS=Import[NotebookDirectory[]<>"/points.txt","Table"];
PTS = Select[PTS, !StringStartsQ[ToString[#[[1]]], "#"] &][[All,;;-1]];
invariants={d,mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56}
PTSRPL=Map[Thread[Rule[invariants,#]]&,PTS[[All,;;-1]]]
Expand[G6/.PTSRPL,Modulus->$P]



