(* ::Package:: *)

(* Thin CLI wrapper around SimplifyOnGram.wl's SimplifyFilesOnGram -- all the
   actual batch logic lives there so it behaves identically whether you call
   it from this script or from a notebook (Get["SetUp/SimplifyOnGram.wl"]
   then SimplifyFilesOnGram[{"Reconstruction-NewPara/NewArg9494.txt"}]).

   Usage:
     math -script SetUp/RunSimplifyOnGram.wls path/to/NewArg9393.txt [more files...]

   Options (environment variables, since this is a plain positional CLI):
     SIMPLIFY_VAR   - variable to reduce in (default: s16)
     SIMPLIFY_GENS  - alternate constraint file instead of G6 (default: G6)
*)

origDir = Directory[];
SetDirectory[DirectoryName[$InputFileName]];
Get["SimplifyOnGram.wl"];
SetDirectory[origDir];

varName = Environment["SIMPLIFY_VAR"];
reduceVar = If[StringQ[varName], ToExpression[varName], s16];

gensFile = Environment["SIMPLIFY_GENS"];
gens = If[StringQ[gensFile], {Get[gensFile]}, {}];

files = Drop[$CommandLine, First[FirstPosition[$CommandLine, "-script"]] + 1];
If[files === {},
  Print["Usage: math -script SetUp/RunSimplifyOnGram.wls file1.txt [file2.txt ...]"];
  Exit[1];
];

SimplifyFilesOnGram[files, gens, reduceVar];
