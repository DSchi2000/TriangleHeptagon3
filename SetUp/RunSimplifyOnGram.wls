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
     SIMPLIFY_BEST  - if set (to anything), use SimplifyOnVarietyFactoredBest
                      instead: tries every invariant as pivot for each
                      factor and keeps whichever gives the smallest
                      polynomial reduction, iterating until no further
                      improvement. Ignores SIMPLIFY_VAR. Slower, but finds
                      reductions a single fixed pivot can miss entirely.
*)

origDir = Directory[];
SetDirectory[DirectoryName[$InputFileName]];
Get["SimplifyOnGram.wl"];
SetDirectory[origDir];

(* $CommandLine only contains "-script" <thisfile> <args...> when actually
   launched as `math -script ...` from a shell. If this file is instead Get'd
   from a notebook or another kernel (e.g. by mistake, or out of habit), there
   is no "-script" token and no file list to read off $CommandLine -- in that
   case just stop here having loaded the package, rather than crashing on a
   Drop[...] with a missing position. Call SimplifyFilesOnGram[...] (or
   SimplifyOnVarietyFactored[...]) directly instead; see SimplifyOnGram.wl. *)
scriptPos = FirstPosition[$CommandLine, "-script", Missing["NotFound"]];

If[MissingQ[scriptPos],

  Print["RunSimplifyOnGram.wls was loaded without `math -script ...` (e.g. ",
    "from a notebook), so there are no command-line files to read. The ",
    "package is loaded -- call its functions directly, e.g.:\n",
    "  SimplifyFilesOnGram[{\"Reconstruction-NewPara/NewArg9494.txt\"}]\n",
    "  SimplifyOnVarietyFactored[NewArg]"]

  ,

  varName = Environment["SIMPLIFY_VAR"];
  reduceVar = If[StringQ[varName], ToExpression[varName], s16];

  gensFile = Environment["SIMPLIFY_GENS"];
  gens = If[StringQ[gensFile], {Get[gensFile]}, {}];

  files = Drop[$CommandLine, First[scriptPos] + 1];
  If[files === {},
    Print["Usage: math -script SetUp/RunSimplifyOnGram.wls file1.txt [file2.txt ...]"];
    Exit[1];
  ];

  If[StringQ[Environment["SIMPLIFY_BEST"]],
    SimplifyFilesOnGramBest[files, gens],
    SimplifyFilesOnGram[files, gens, reduceVar]
  ];
];
