(* ::Package:: *)

(* ::Section:: *)
(*Initialization Block*)


Get["IniBlocks.m"]
Get[HomeDir<>"/MassiveBCFW/MassiveBCFWSlice.wl"]
PointGeneratorPath=HomeDir<>"/SQRT-Generator";
GramDir=HomeDir<>"/SetUp";

GetTouple[n_,tpoint_]:=Module[{tpoints},
tpoints=Table[tpoint,{n}];
tpoints[[All,1,2]]=RandomInteger[$P,n];
Return[tpoints];
 ]  


(* ::Section:: *)
(*Generate Samples*)


SetDirectory[CurrentDir];

Print["Start Kernels"];

LaunchKernels[NParrallel];

Print["Generate Points"];

FullTimeUsed=AbsoluteTiming[If[SampleMode=="Slice",Print["Sampling Mode Slice: "<>ToString[SliceVar]];
If[reset===True,
Print["Generate Initial Vectors"];
Command=PointGeneratorPath<>"/gen_6pt --roots "<>SectorDir<>"/Roots.txt --gram "<>PointGeneratorPath<>"/SixPointGram.txt --npoints 1 --auxvar "<>ToString[AuxVar]<>" --out points.txt --threads "<>ToString[NParrallel];
Run[Command];
PTS=Import[CurrentDir<>"/points.txt","Table"];
PTS = Select[PTS, !StringStartsQ[ToString[#[[1]]], "#"] &][[All,;;-1]];

PTSRPL=Map[Thread[Rule[invariants,#]]&,PTS];
aVec=PTSRPL[[1]];
Export["AVECS.txt",{aVec}];,
aVec=Get["AVECS.txt"]
];

tpoint = {d, mt2, s12, s123, s16, s23, s234, s34, s345, s45, s56}/.aVec;
freeVals = Delete[tpoint, First[Flatten[Position[invariants,AuxVar]]]]; (* drop s16 *)
Print[freeVals];

Export["slice_a.txt", {StringRiffle[ToString /@ freeVals, " "]}, "Lines"];
SlicePos=Flatten[Position[DeleteElements[invariants,{AuxVar}],SliceVar]];
bvec=ConstantArray[0,10];
bvec[[SlicePos]]=1;
Export["slice_b.txt", {StringRiffle[ToString /@ bvec, " "]}, "Lines"];

Command=PointGeneratorPath<>"/gen_6pt --roots "<>SectorDir<>"/Roots.txt --gram "<>PointGeneratorPath<>"/SixPointGram.txt --npoints "<>ToString[NSamples]<>" --auxvar "<>ToString[AuxVar]<>" --slice slice_a.txt slice_b.txt"<>" --out points.txt --threads "<>ToString[NParrallel];
Print[Command];
Run[Command];

PTS=Import[CurrentDir<>"/points.txt","Table"];
PTS = Select[PTS, !StringStartsQ[ToString[#[[1]]], "#"] &][[All,;;-2]];

PTSRPL=Map[Thread[Rule[invariants,#]]&,PTS];

tpointlist=ParallelTable[ GetTouple[EpsRecoLenth,PTSRPL[[i]]],{i,1,NSamples}];


];
If[SampleMode=="UVar",
Print["Sampling Mode UnivariateSlice: "<>ToString[t]];
If[reset===True,
Print["Generate Initial Vectors"];
res=BuildMassiveSlice[RandomInteger[$P], False];
avec=Normal[res["a"]][[All;;10]]//Echo;
bvec=Normal[res["b"]][[All;;10]]//Echo;
avec[[All,1]]=ToExpression[avec[[All,1]]];
bvec[[All,1]]=ToExpression[bvec[[All,1]]];

;
Export["AVECS.txt",{avec,bvec}];,
aVec=Get["AVECS.txt"]
];


Export["slice_a.txt", {StringRiffle[ToString /@ Join[{RandomInteger[$P]},avec[[All,2]]], " "]}, "Lines"];
Export["slice_b.txt", {StringRiffle[ToString /@ Join[{RandomInteger[$P]},bvec[[All,2]]], " "]}, "Lines"];


Command=PointGeneratorPath<>"/gen_6pt --roots "<>SectorDir<>"/Roots.txt --gram "<>GramDir<>"/SixPointGram.txt --npoints "<>ToString[NSamples]<>" --auxvar "<>ToString[s16]<>" --slice slice_a.txt slice_b.txt"<>" --out points.txt --threads "<>ToString[NParrallel]<>" --gram-resolved";
Print[Command];
Run[Command];
(*Quit[];*)
PTS=Import[CurrentDir<>"/points.txt","Table"];
PTS = Select[PTS, !StringStartsQ[ToString[#[[1]]], "#"] &][[All,;;-1]];

PTSRPL=Map[Thread[Rule[invariants,#]]&,PTS[[All,;;-2]]];
tvars=Map[Thread[Rule[{tvar},#]]&,PTS[[All,-1]]];

tpointlist=ParallelTable[ GetTouple[EpsRecoLenth,PTSRPL[[i]]],{i,1,NSamples}];
Export[ToString[StringForm["``/``/tvars.txt",CurrentDir,OutPutPath]],tvars];


];
If[SampleMode=="Full",Print["3"]];]

tpointlist=Flatten[tpointlist,1][[All]];
Export[ToString[StringForm["``/``/phasespacepoints.txt",CurrentDir,OutPutPath]],tpointlist[[All,All,2]]];



Print["Write Points"];

(*Export["univar.txt",parList[[All,2]]//DeleteDuplicates];*)
SetDirectory[SectorDir];

DeleteOldKiraRun[];
For[Npos=1,Npos<=NSamples*EpsRecoLenth,Npos++,
If[Npos==1,
WriteNumericsFile[tpointlist[[Npos]]];,
AddToNumericsFile[tpointlist[[Npos]]];
];
];

SetDirectory[CurrentDir];
Export["TIME.txt",{"Generate Points",FullTimeUsed[[1]]}];
Run[ToString[StringForm["cd `` && python3 mailsend.py",CurrentDir]]];
Print["Done"];
Quit[];
