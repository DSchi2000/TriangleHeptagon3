(* ::Package:: *)

(* Standalone benchmark: compares the existing symbolic D[] derivative step
   (HexaBoxModules.m:305) against the new DualEval mod-p forward-mode AD
   (DualDiff.m, this folder), on real Coefs expressions pulled from real sample
   points. Does not modify any existing file -- ComputeDIFastExtract below
   is a verbatim copy of ComputeDIFast's body up to (and including) the
   point where Coefs/delcases/DellCoeffs are built, under a new name. *)

ScriptDir=Directory[];
ProjectRoot=ParentDirectory[DirectoryName[$InputFileName]];
SetDirectory[FileNameJoin[{ProjectRoot,"Reconstruction1","Samples"}]];
Get["IniBlocks.m"]
Get[FileNameJoin[{ProjectRoot,"DiffSpeedupExploration","DualDiff.m"}]];

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

(* ---- collect real test cases: Coefs/delcases for several real sample points ---- *)
samplepoints={1,2,3};
intposSample=Range[1,Length[MaxCutInts]]; (* all integrals, for a representative sample *)

testcases={};
Do[
  SubReduList=Get[ToString[StringForm["``/KIRA/splitres/``/phasespace_``.m",kirapath,IntegerString[Floor[(sp-1)/1000],10,4],IntegerString[sp,10,7]]]];
  tPoint=Thread[Rule[invariants,SubReduList[[2,1]]]];
  diffOp=getDifferentialOperator[tPoint];
  Do[
    tPointCUT=tPoint[[Complement[Range[11],Flatten[Map[Position[invariants,#]&,{invariantsnod[[var]],InternalAuxvar}]]]]];
    Do[
      DiffIntegral=(PureIntegralsMaxCut[[MaxCutInts[[intpos]]]]/.tPointCUT)/.kinematicFunctions;
      res=Check[ComputeDIFastExtract[DiffIntegral,diffOp,var,tPoint],$Failed];
      If[res=!=$Failed,AppendTo[testcases,Append[res,invariantsnod[[var]]]]];
    ,{intpos,intposSample}]
  ,{var,SubdVarLIST}];
  Print["Sample point ",sp," done, ",Length[testcases]," test cases so far."];
,{sp,samplepoints}];

Print["Collected ",Length[testcases]," (Coefs,delcases,DellCoeffs,tPoint) test cases from real data."];

(* ---- complexity diagnostics: how big are real Coefs lists actually? ---- *)
leafCounts=Flatten[Map[LeafCount,#[[1]]]&/@testcases];
Print["Coefs entry LeafCount: min=",Min[leafCounts],"  max=",Max[leafCounts],"  mean=",N[Mean[leafCounts]]];

(* ---- correctness + timing: old symbolic D[] vs new DualEval ----
   Coefs is itself a LIST (one entry per distinct HPT3 basis integral --
   Coefficient[expr,{a,b,...}] returns a list), so both D[] and DualEval
   are applied across that list.
   Coefs is NOT univariate: DeleteList (HexaBoxModules.m:285) always keeps
   InternalAuxvar (s16) symbolic in addition to whichever invariant this
   ComputeDIFast call differentiates against ("ownInvariant" below) -- so
   Coefs is a function of exactly those two symbols. delcases walks every
   direction appearing in the differential-operator ansatz, which includes:
     - momentum directions del[p[j]]            -> Coefs provably has no
       dependence on raw momenta, derivative is always 0 (skipped).
     - del[ownInvariant] or del[s16]             -> genuine partial
       derivative: substitute the OTHER symbol's known numeric value first
       (valid since partial differentiation treats it as a constant either
       way), then the expression is univariate and DualEval applies.
     - del[anyOtherInvariant]                    -> Coefs has no dependence
       on it either, derivative is always 0 (skipped). *)
oldTotal=0;newTotal=0;mismatches=0;ncomparisons=0;

Do[
  {Coefs,delcases,DellCoeffs,tPointTC,ownInvariant}=tc;
  Do[
    var2=delcases[[dv]]/.del[x_]:>x;
    relevant=Head[var2]===Symbol&&MemberQ[{ownInvariant,InternalAuxvar},var2];
    If[relevant,
      otherVar=First[Complement[{ownInvariant,InternalAuxvar},{var2}]];
      otherVal=otherVar/.tPointTC;
      x0=var2/.tPointTC;
      CoefsUnivar=Coefs/.otherVar->otherVal;

      tOld=First[AbsoluteTiming[
        dPrefacOld=Expand[(D[Coefs,var2]/.tPointTC)/.root'[x_]:>1/(2 root[x])//.root[x_]:>PowerMod[x,1/2,$P],Modulus->$P];
      ]];
      tNew=First[AbsoluteTiming[
        dualRes=Map[DualEval[#,var2,x0,$P]&,CoefsUnivar];
        dPrefacNew=dualRes[[All,2]];
      ]];

      oldTotal+=tOld;newTotal+=tNew;ncomparisons++;
      If[dPrefacOld=!=dPrefacNew,
        mismatches++;
        Print["MISMATCH  var=",var2,"  old=",dPrefacOld,"  new=",dPrefacNew];
      ];
    ];
  ,{dv,1,Length[delcases]}];
,{tc,testcases}];

Print["============================================"];
Print["Comparisons: ",ncomparisons,"   Mismatches: ",mismatches];
Print["Total time, symbolic D[] approach : ",oldTotal," s"];
Print["Total time, DualEval approach     : ",newTotal," s"];
Print["Speedup                           : ",N[oldTotal/newTotal],"x"];

SetDirectory[ScriptDir];
