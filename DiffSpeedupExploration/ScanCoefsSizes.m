(* ::Package:: *)

(* Scan Coefs expression sizes (LeafCount) across a broad set of real sample
   points and all (var,intpos) pairs, to check whether any real cases in
   this sector cross the ~150-300 LeafCount threshold found in
   BenchDualDiff.m's scaling test (above which DualEval beats symbolic D[]).
   Does not modify any existing file. *)

ScriptDir=Directory[];
ProjectRoot=ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[FileNameJoin[{ProjectRoot,"Reconstruction1","Samples"}]];
Get["IniBlocks.m"]

NotbookDir=CurrentDir;
kirapath=SectorDir;
SetDirectory[HomeDir];
Ureduce=ReadList[FileNameJoin[{kirapath,"KIRA/UReduList.txt"}]];
changeAuxVar[ToExpression[Capitalize[ToString[auxvar]]]];
SetDirectory[NotbookDir];

ComputeDIFastExtract[toDiffInt_,difOp_,diffVarPos_,tPoint_]:=Module[{epsrepl,toDiffIntClean,diffOpPart,DeleteList,mandelstamScalarProductsT,loopMomentumScalarProductsT,DCases,DellCoeffs,delcases,Coefs},
loopMomentumScalarProductsT=Expand[loopMomentumScalarProducts/.tPoint,Modulus->$P];
mandelstamScalarProductsT=Expand[Join[mandelstamScalarProducts,kinematics]/.tPoint,Modulus->$P];
DeleteList=Delete[{1,2,3,4,5,6,7,8,9,10,11},Position[tPoint[[All,1]],InternalAuxvar]];
epsrepl={(ep->(4-d)/2)/.tPoint};

toDiffIntClean=toDiffInt/.tPoint[[DeleteCases[DeleteList,diffVarPos+1]]]/.epsrepl;

diffOpPart=Collect[Coefficient[difOp,ds[invariantsnodUpperCase[[diffVarPos]]]],del[___],Together];

DCases=Cases[{toDiffIntClean},HPT3[___],Infinity]//DeleteDuplicates;
If[Length[DCases]==0,
DCases={HPT3[0,0,0,0,0,0,0,0,0,0,0,0,0]};
toDiffIntClean=toDiffIntClean HPT3[0,0,0,0,0,0,0,0,0,0,0,0,0];
];
Coefs=Collect[Coefficient[toDiffIntClean,DCases],root[___],Together];
delcases=Cases[diffOpPart,del[___],Infinity]//Sort;
DellCoeffs=Coefficient[diffOpPart,delcases];
Return[{Coefs,delcases,DellCoeffs,tPoint}];
]

SubdVarLIST=Flatten[Map[Position[invariantsnod,#]&,SubdVars]];

(* which sample points actually have split files on disk *)
existing=Select[Join[Range[1,500],Range[1001,1500]],
  FileExistsQ[ToString[StringForm["``/KIRA/splitres/``/phasespace_``.m",kirapath,IntegerString[Floor[(#-1)/1000],10,4],IntegerString[#,10,7]]]]&];
Print["Sample points available on disk: ",Length[existing]];

samplepoints=existing[[Round[Subdivide[1,Length[existing],19]]]]//DeleteDuplicates;
Print["Testing sample points: ",samplepoints];

allLeafCounts={};
maxEntry={-1,Null,Null,Null}; (* {leafcount,sp,var,intpos} *)

Do[
  SubReduList=Get[ToString[StringForm["``/KIRA/splitres/``/phasespace_``.m",kirapath,IntegerString[Floor[(sp-1)/1000],10,4],IntegerString[sp,10,7]]]];
  tPoint=Thread[Rule[invariants,SubReduList[[2,1]]]];
  diffOp=getDifferentialOperator[tPoint];
  Do[
    tPointCUT=tPoint[[Complement[Range[11],Flatten[Map[Position[invariants,#]&,{invariantsnod[[var]],InternalAuxvar}]]]]];
    Do[
      DiffIntegral=(PureIntegralsMaxCut[[MaxCutInts[[intpos]]]]/.tPointCUT)/.kinematicFunctions;
      res=Check[ComputeDIFastExtract[DiffIntegral,diffOp,var,tPoint],$Failed];
      If[res=!=$Failed,
        lcs=LeafCount/@res[[1]];
        allLeafCounts=Join[allLeafCounts,lcs];
        thismax=Max[lcs];
        If[thismax>maxEntry[[1]],maxEntry={thismax,sp,var,intpos}];
      ];
    ,{intpos,1,Length[MaxCutInts]}]
  ,{var,SubdVarLIST}];
  Print["Sample point ",sp," done. Running max LeafCount so far: ",maxEntry[[1]]," (sp=",maxEntry[[2]],",var=",maxEntry[[3]],",intpos=",maxEntry[[4]],")"];
,{sp,samplepoints}];

Print["============================================"];
Print["Total Coefs entries scanned: ",Length[allLeafCounts]];
Print["LeafCount: min=",Min[allLeafCounts],"  max=",Max[allLeafCounts],"  mean=",N[Mean[allLeafCounts]],"  median=",N[Median[allLeafCounts]]];
Print["Percentiles: 90%=",N[Quantile[allLeafCounts,0.90]],"  99%=",N[Quantile[allLeafCounts,0.99]],"  99.9%=",N[Quantile[allLeafCounts,0.999]]];
Print["Count above 150 (crossover region): ",Length[Select[allLeafCounts,#>150&]]," / ",Length[allLeafCounts]];
Print["Count above 300: ",Length[Select[allLeafCounts,#>300&]]];
Print["Largest single entry: ",maxEntry[[1]]," at sample point ",maxEntry[[2]],", var=",invariantsnod[[maxEntry[[3]]]],", intpos=",maxEntry[[4]]];

SetDirectory[ScriptDir];
