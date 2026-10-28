(* ::Package:: *)

(* ::Title:: *)
(*Initialization Block*)


Get["IniBlocks.m"]


(* ::Title:: *)
(*Generate Differential Equations*)


NotbookDir=CurrentDir;
kirapath=SectorDir;
If[StringContainsQ[kirapath,"dschiebe"],
If[StringContainsQ[kirapath,"Documents"],
kirabinpathlanayru="/usr/local",
kirabinpathlanayru="/home/dschiebe/installs/kira3/bin"
],
kirabinpathlanayru="/home/daniel/Documents/packages/kira/bin"
];


SetDirectory[FileNameJoin[{kirapath,"KIRA"}]];
runtime=AbsoluteTiming[Run[ToString[StringForm["``/bin/kira --parallel=`` jobs.yaml",kirabinpathlanayru,KIRACores]]];
			   Run["python3 SplitNumericalOutput.py -i './results/HPT3/numerics_2147483647_0.m' -o './splitres'"];];


SetDirectory[CurrentDir];
Export["TIME.txt",{"IBP Reduction",runtime[[1]]}];
Run[ToString[StringForm["cd `` && python3 mailsend.py",CurrentDir]]];
Print["Done"];
