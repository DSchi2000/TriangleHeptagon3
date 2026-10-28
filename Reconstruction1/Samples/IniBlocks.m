(* ::Package:: *)

(* ::Title:: *)
(*Settings*)


NSamples=200;
AuxVar=s16;
$P=2^31-1;
(*OutPutPath="Samples-1";*)
OutPutPath="Univariate Slice";
SliceVar=s16;
NParrallel=20;
KIRACores=20;
EpsRecoLenth=1;
MaxCutInts=Range[95];
reset=True;
printOut=True;
AllSectors=True;
SampleMode="UVar"; (*Slice,UVar,Full*)
SubdVars={mt2,s12,s123,s23,s234,s34,s345,s45,s56};
DoExtraRelations=True;


(*SectorOnline= hxb[1,1,1,1,1,1,1,1,0,0,0,0,0];*)


(* ::Title:: *)
(*Initialization Block *)


HomeDir=ParentDirectory[ParentDirectory[Directory[]]];
SectorDir=ParentDirectory[Directory[]];
CurrentDir=Directory[];
SetDirectory[HomeDir];
Get[FileNameJoin[{HomeDir,"files/HexaBoxModules.m"}]];
Get[FileNameJoin[{HomeDir,"files/DifferentialOperator.m"}]];
sqrtSet=ReadList[FileNameJoin[{SectorDir,"Roots.txt"}]];

PureIntegrals=ReadList[FileNameJoin[{SectorDir,"SubBasis.txt"}]];
SectorLaporta=ReadList[FileNameJoin[{SectorDir,"SectorLaportaBasis.txt"}]];
PureIntegralsMaxCut=SectorLaporta;
PureIntegralsMaxCut[[MaxCutInts]]=PureIntegrals[[MaxCutInts]];
If[!AllSectors,
PureIntegralsMaxCut[[MaxCutInts]]=PureIntegralsMaxCut[[MaxCutInts]]/.hxb[n___]:>OffCutTest[hxb[n]];
]
sqrtSet=(Cases[PureIntegralsMaxCut,_root,Infinity]//DeleteDuplicates)/.root[x_]:>x;

SetDirectory[CurrentDir];
Get[FileNameJoin[{HomeDir,"files/ExtraFunctions.m"}]];


    
auxvar=AuxVar;
