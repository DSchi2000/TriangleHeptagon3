(* ::Package:: *)

(* Correctness + speed check: ComputeDISIJ (ExtraFunctions.m, unmodified) vs
   ComputeDISIJHoisted (HoistedModules.m, this folder) on real sample points.
   Does not modify any existing file. *)

ScriptDir=Directory[];
ProjectRoot=ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[FileNameJoin[{ProjectRoot,"Reconstruction1","Samples"}]];
Get["IniBlocks.m"]
Get[FileNameJoin[{ProjectRoot,"DiffSpeedupExploration","HoistedModules.m"}]];

NotbookDir=CurrentDir;
kirapath=SectorDir;
SetDirectory[HomeDir];
Ureduce=ReadList[FileNameJoin[{kirapath,"KIRA/UReduList.txt"}]];
changeAuxVar[ToExpression[Capitalize[ToString[auxvar]]]];
SetDirectory[NotbookDir];

samplepoints={1,2,3,4,5};

oldTotal=0;newTotal=0;mismatches=0;

Do[
  tOld=First[AbsoluteTiming[resOld=ComputeDISIJ[sp];]];
  tNew=First[AbsoluteTiming[resNew=ComputeDISIJHoisted[sp];]];
  oldTotal+=tOld;newTotal+=tNew;
  match=(resOld===resNew);
  If[!match,mismatches++];
  Print["Sample point ",sp,"  tOld=",tOld,"s  tNew=",tNew,"s  ratio=",N[tOld/tNew],"x  match=",match];
,{sp,samplepoints}];

Print["============================================"];
Print["Sample points tested: ",Length[samplepoints],"   Mismatches: ",mismatches];
Print["Total time, ComputeDISIJ (original)        : ",oldTotal," s"];
Print["Total time, ComputeDISIJHoisted (new file) : ",newTotal," s"];
Print["Speedup                                    : ",N[oldTotal/newTotal],"x"];

SetDirectory[ScriptDir];
