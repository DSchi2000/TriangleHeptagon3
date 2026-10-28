(* ::Package:: *)

(* ::Title:: *)
(*Initialization Block*)


Get["IniBlocks.m"]
Print["Deleting tmp-files: ",FileNameJoin[{CurrentDir,ToString[StringForm["``/tmp-files",OutPutPath]]}]];
DeleteDirectory[FileNameJoin[{CurrentDir,ToString[StringForm["``/tmp-files",OutPutPath]]}],DeleteContents -> True];
