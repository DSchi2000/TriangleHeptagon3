(* ::Package:: *)

(* ::Title:: *)
(*Initialization Block*)


Get["IniBlocks.m"]

batchSitze=5000;
(*NParrallel=96;*)
(*NParrallel=10;
NSamples=100;
*)


(* ::Title:: *)
(*Generate Differential Equations*)


NotbookDir=CurrentDir;
kirapath=SectorDir;
SetDirectory[HomeDir];
If[FileExistsQ[FileNameJoin[{kirapath,"KIRA/splitres/0000/phasespace_0000001.m"}]],

AlreadyCompPTS={};
AlreadyCompPTS=Flatten[Table[If[FileExistsQ[FileNameJoin[{NotbookDir,ToString[StringForm["``/tmp-files/DonePts_K``",OutPutPath,IntegerString[i,10,3]]]}]],
ReadList[FileNameJoin[{NotbookDir,ToString[StringForm["``/tmp-files/DonePts_K``",OutPutPath,IntegerString[i,10,3]]]}]],{0}],{i,1,1000}]];
AlreadyCompPTS=DeleteCases[AlreadyCompPTS,0];
Print[AlreadyCompPTS];

If[AlreadyCompPTS==={},StartingPoint=1,
AllExpected=Range[Max[AlreadyCompPTS]];
StartingPoint=Min[Complement[AllExpected,AlreadyCompPTS]];
If[Complement[AllExpected,AlreadyCompPTS]==={},StartingPoint=Max[AlreadyCompPTS]];
StartingPoint=StartingPoint-Mod[StartingPoint,batchSitze]+1;
];

Print[ToString[StringForm["Start at ``",StartingPoint]]];

CloseKernels[];
Print[ToString[StringForm["Launching `` Kernels",NParrallel]]];
(*kernels=LaunchKernels[NParrallel];*)
Print["Loading Files on Kernels"];
Ureduce=ReadList[FileNameJoin[{kirapath,"KIRA/UReduList.txt"}]];
changeAuxVar[ToExpression[Capitalize[ToString[auxvar]]]];
(*ParallelEvaluate[Get[FileNameJoin[{Directory[],"files/HexaBoxModules.m"}]],kernels];
ParallelEvaluate[Get[FileNameJoin[{Directory[],"files/DifferentialOperator.m"}]],kernels];
ParallelEvaluate[Ureduce=ReadList[FileNameJoin[{kirapath,"KIRA/UReduList.txt"}]];,kernels];
ParallelEvaluate[changeAuxVar[ToExpression[Capitalize[ToString[auxvar]]]],kernels];
ParallelEvaluate[SetDirectory[NotbookDir];,kernels];
ParallelEvaluate[PureIntegrals=ReadList[FileNameJoin[{SectorDir,"SubBasis.txt"}]];,kernels];
ParallelEvaluate[PureIntegralsMaxCut=SectorLaporta;
			     PureIntegralsMaxCut[[MaxCutInts]]=PureIntegrals[[MaxCutInts]];,kernels];
ParallelEvaluate[Get[FileNameJoin[{HomeDir,"files/ExtraFunctions.m"}]],kernels];
ParallelEvaluate[printOut=True;,kernels];*)
SetDirectory[NotbookDir];

Print["Compute DEs"];

SampleSize=EpsRecoLenth*NSamples;
CreateDirectory[FileNameJoin[{NotbookDir,OutPutPath,"tmp-files"}]];

(*ParallelEvaluate[SetOptions[OpenAppend,PageWidth->Infinity];,kernels];*)



(*Do[
ParallelDo[
tpoint=Get[ToString[StringForm["``/KIRA/splitres/``/phasespace_``.m",kirapath,IntegerString[Floor[(samplepointpos-1)/1000],10,4],IntegerString[samplepointpos,10,7]]]][[2,1]];
ClearAll[tpoint];
ClearSystemCache[];
,
{samplepointpos,batchpos,batchpos+batchSitze-1}];
Print[ToString[StringForm["##### ``-`` Done!  #####",batchpos,batchpos+batchSitze-1]]];
Print[ToString[StringForm["Current Size in Memory: `` GB",N[MemoryInUse[]/10^9,3]]]];
,{batchpos,1,SampleSize,batchSitze}];*)




runtime=AbsoluteTiming[
Do[
Print[ToString[StringForm["##### Start ``-``! #####",batchpos,Min[{batchpos+batchSitze-1,SampleSize}]]]];
Print[ToString[StringForm["Launching `` Kernels",NParrallel]]];
kernels=LaunchKernels[NParrallel];
ParallelEvaluate[SetOptions[OpenAppend,PageWidth->Infinity];,kernels];
ParallelDo[
If[MemberQ[AlreadyCompPTS,samplepointpos]===False,
tmpres=ComputeDISIJ[samplepointpos];
PutAppend[tmpres[[1]],FileNameJoin[{NotbookDir,ToString[StringForm["``/tmp-files/Parameters_K``",OutPutPath,IntegerString[$KernelID,10,3]]]}]];
Table[PutAppend[tmpres[[2,varpos]],FileNameJoin[{NotbookDir,ToString[StringForm["``/tmp-files/DeMats-``_K``",OutPutPath,SubdVars[[varpos]],IntegerString[$KernelID,10,3]]]}]];,{varpos,1,Length[SubdVars]}];
PutAppend[samplepointpos,FileNameJoin[{NotbookDir,ToString[StringForm["``/tmp-files/DonePts_K``",OutPutPath,IntegerString[$KernelID,10,3]]]}]];
ClearAll[tmpres];
ClearSystemCache[];
];
,
{samplepointpos,batchpos,Min[{batchpos+batchSitze-1,SampleSize}]}];
Print[ToString[StringForm["##### ``-`` Done!  #####",batchpos,Min[{batchpos+batchSitze-1,SampleSize}]]]];
Print[ToString[StringForm["Current Size in Memory: `` GB",N[MemoryInUse[]/10^9,3]]]];
CloseKernels[];
Print["Closing Kernels! for restart (Memory Control)"];
If[N[MemoryInUse[]/10^9,3]>200,
Print["MEMORY OVERFLOW"];
Quit[];
]
,{batchpos,StartingPoint,SampleSize,batchSitze}];

Print["Load tmp Parameters"];
IniPameters=ReadList[ToString[StringForm["``/``/phasespacepoints.txt",CurrentDir,OutPutPath]]];

PhaseSpacePointsLoad=DeleteCases[Flatten[Table[If[FileExistsQ[FileNameJoin[{NotbookDir,ToString[StringForm["``/tmp-files/Parameters_K``",OutPutPath,IntegerString[i,10,3]]]}]],
ReadList[FileNameJoin[{NotbookDir,ToString[StringForm["``/tmp-files/Parameters_K``",OutPutPath,IntegerString[i,10,3]]]}]],{0}],{i,1,1000}],1],0];

Print["Find Permutation to Initial Parameters"];
Perm=FindPermutation[PhaseSpacePointsLoad,IniPameters];


(*PtsPos=Flatten[Map[Position[PhaseSpacePointsLoad,#]&,IniPameters[[;;5]]]];
Table[DeLoad[invariants[[varpos]]]=DeLoad[invariants[[varpos]]][[PtsPos]],{varpos,2,Length[invariants]}];
PhaseSpacePointsLoad=PhaseSpacePointsLoad[[PtsPos]];
IniPameters=IniPameters[[;;5]];*)



Print["Load DE-Evaluations"];
Table[
runtime2=AbsoluteTiming[DeLoad[SubdVars[[varpos]]]=DeleteCases[Flatten[Table[If[FileExistsQ[FileNameJoin[{NotbookDir,ToString[StringForm["``/tmp-files/Parameters_K``",OutPutPath,IntegerString[i,10,3]]]}]],
ReadList[FileNameJoin[{NotbookDir,ToString[StringForm["``/tmp-files/DeMats-``_K``",OutPutPath,SubdVars[[varpos]],IntegerString[i,10,3]]]}]],{Null}],{i,1,1000}],1],Null];];
Print[ToString[StringForm["Done `` in ``s",SubdVars[[varpos]],runtime2[[1]]]]];
If[N[MemoryInUse[]/10^9,3]>200,
Print["MEMORY OVERFLOW"];
Quit[];];,{varpos,1,Length[SubdVars]}];
Table[DeLoad[SubdVars[[varpos]]]=Permute[DeLoad[SubdVars[[varpos]]],Perm];,{varpos,1,Length[SubdVars]}];
];

Print[ToString[StringForm["Time: ``s",runtime[[1]]]]];
Print[ToString[StringForm["Saving DEs in Chuncks of ``",EpsRecoLenth]]];

SetOptions[OpenAppend,PageWidth->Infinity];
CreateDirectory[FileNameJoin[{NotbookDir,OutPutPath}]];

(*Table[
Export[FileNameJoin[{NotbookDir,ToString[StringForm["``/DEmatrix_d``.txt",OutPutPath,invariants[[varpos]]]]}],DeLoad[invariants[[varpos]]]],{varpos,2,Length[invariants]}];
*)

Table[
Print[ToString[StringForm["Exporting ``/``",psp,SubdVars[[varpos]]]]];
numname=IntegerString[Floor[(psp-1)/1000],10,4];
dir=FileNameJoin[{NotbookDir,ToString[StringForm["``/tmp-files/DEMatsOut_``",OutPutPath,numname]]}];
If[!DirectoryQ[dir],
    CreateDirectory[dir, CreateIntermediateDirectories -> True];
];
Export[FileNameJoin[{NotbookDir,ToString[StringForm["``/tmp-files/DEMatsOut_``/DEmatrix_d``-``.txt",OutPutPath,numname,SubdVars[[varpos]],IntegerString[psp,10,7]]]}],DeLoad[SubdVars[[varpos]]][[EpsRecoLenth*psp-EpsRecoLenth+1;;EpsRecoLenth psp]]];
Export[FileNameJoin[{NotbookDir,ToString[StringForm["``/tmp-files/DEMatsOut_``/Parameters-``.txt",OutPutPath,numname,IntegerString[psp,10,7]]]}],IniPameters[[EpsRecoLenth*psp-EpsRecoLenth+1;;EpsRecoLenth psp]]];
,
{varpos,1,Length[SubdVars]},{psp,1,NSamples}];



Print["Done"];
SetDirectory[CurrentDir];
Export["TIME.txt",{"Build Differential Equation",runtime[[1]]}];
Run[ToString[StringForm["cd `` && python3 mailsend.py",CurrentDir]]];
Print["Done"];

,Print["NO IBP REDUCTION"];
SetDirectory[CurrentDir];
Export["TIME.txt",{"NO IBP REDUCTION",0}];
Run[ToString[StringForm["cd `` && python3 mailsend.py",CurrentDir]]];
]






(*AlreadyCompPTS=DeleteCases[Range[12340],8902];
AllExpected=Range[Max[AlreadyCompPTS]];
StartingPoint=Min[Complement[AllExpected,AlreadyCompPTS]];
batchSitze=1000;
StartingPoint=StartingPoint-Mod[StartingPoint,batchSitze];*)


(*EpsRecoLenth*psp-EpsRecoLenth+1+EpsRecoLenth-1//FullSimplify*)
