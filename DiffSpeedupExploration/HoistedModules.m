(* ::Package:: *)

(* Hoisted variant of ComputeDIFast/ComputeDISIJ (HexaBoxModules.m / ExtraFunctions.m).
   Does not modify either existing file.

   ComputeDISIJ calls ComputeDIFast once per (var,intpos) pair -- 9*95=855 times
   per sample point -- always with the same tPoint, and with the same difOp/var
   for all 95 intpos in a row. But inside ComputeDIFast, several quantities only
   depend on tPoint (loopMomentumScalarProductsT, mandelstamScalarProductsT,
   DeleteList, epsrepl) or only on (difOp,var) (the substitution list built from
   DeleteList, diffOpPart, delcases, DellCoeffs, and the driver's own tPointCUT) --
   none of those depend on toDiffInt/intpos. The original recomputes all of them
   855 times; this splits ComputeDIFast into:
     - PrecomputeTPointData : once per sample point
     - PrecomputeVarData    : once per var (9x instead of 855x)
     - ComputeDIFastHoisted : the remaining per-(var,intpos) work
   ComputeDISIJHoisted is a verbatim copy of ComputeDISIJ wired to use these,
   with the driver's own dead kinematicFunctionsN computation (computed but never
   used by ComputeDISIJ -- DiffIntegral is built from kinematicFunctions, not
   kinematicFunctionsN) dropped. *)

PrecomputeTPointData[tPoint_]:=Module[{loopMomentumScalarProductsT,mandelstamScalarProductsT,DeleteList,epsrepl},
loopMomentumScalarProductsT=Expand[loopMomentumScalarProducts/.tPoint,Modulus->$P];
mandelstamScalarProductsT=Expand[Join[mandelstamScalarProducts,kinematics]/.tPoint,Modulus->$P];
DeleteList=Delete[{1,2,3,4,5,6,7,8,9,10,11},Position[tPoint[[All,1]],InternalAuxvar]];
epsrepl={(ep->(4-d)/2)/.tPoint};
Return[<|"tPoint"->tPoint,"loopMomentumScalarProductsT"->loopMomentumScalarProductsT,"mandelstamScalarProductsT"->mandelstamScalarProductsT,"DeleteList"->DeleteList,"epsrepl"->epsrepl|>];
];

PrecomputeVarData[difOp_,diffVarPos_,tPoint_,tPointData_]:=Module[{diffOpPart,delcases,DellCoeffs,substList,tPointCUT},
substList=tPoint[[DeleteCases[tPointData["DeleteList"],diffVarPos+1]]];
tPointCUT=tPoint[[Complement[Range[11],Flatten[Map[Position[invariants,#]&,{invariantsnod[[diffVarPos]],InternalAuxvar}]]]]];
diffOpPart=Collect[Coefficient[difOp,ds[invariantsnodUpperCase[[diffVarPos]]]],del[___],Together];
delcases=Cases[diffOpPart,del[___],Infinity]//Sort;
DellCoeffs=Coefficient[diffOpPart,delcases];
Return[<|"substList"->substList,"tPointCUT"->tPointCUT,"delcases"->delcases,"DellCoeffs"->DellCoeffs|>];
];

ComputeDIFastHoisted[toDiffInt_,tPointData_,varData_]:=Module[{toDiffIntClean,DCases,delInt1,dInt,dPrefac,CoefsT,FinalD,Coefs,ResInt,tPoint,loopMomentumScalarProductsT,mandelstamScalarProductsT,epsrepl,delcases,DellCoeffs},
tPoint=tPointData["tPoint"];
loopMomentumScalarProductsT=tPointData["loopMomentumScalarProductsT"];
mandelstamScalarProductsT=tPointData["mandelstamScalarProductsT"];
epsrepl=tPointData["epsrepl"];
delcases=varData["delcases"];
DellCoeffs=varData["DellCoeffs"];

toDiffIntClean=toDiffInt/.varData["substList"]/.epsrepl;

DCases=Cases[{toDiffIntClean},HPT3[___],Infinity]//DeleteDuplicates;
If[Length[DCases]==0,
DCases={HPT3[0,0,0,0,0,0,0,0,0,0,0,0,0]};
toDiffIntClean=toDiffIntClean HPT3[0,0,0,0,0,0,0,0,0,0,0,0,0];
];
Coefs=Collect[Coefficient[toDiffIntClean,DCases],root[___],Together];
ResInt=Sum[delInt1=Expand[delHPT3[delcases[[delvar]]/.del[x_]:>x]DCases];
dInt=Expand[Expand[DellCoeffs[[delvar]]*delInt1]/.loopMomentumScalarProductsT/.mandelstamScalarProductsT/.tPoint,Modulus->$P];
dPrefac=Expand[((DellCoeffs[[delvar]]D[Coefs,delcases[[delvar]]/.del[x_]:>x])/.tPoint)/.root'[x_]:>1/(2 root[x])//.root[x_]:>PowerMod[x,1/2,$P],Modulus->$P];
CoefsT=Expand[Coefs/.tPoint//.root[x_]:>PowerMod[x,1/2,$P],Modulus->$P];
FinalD=Expand[Total[Expand[dInt*CoefsT+dPrefac*DCases,Modulus->$P]]];
FinalD=FinalD/.n_Integer HPT3[arg___]:>Mod[n,$P]HPT3[arg];
FinalD
,{delvar,1,Length[delcases]}];
If[DCases=={HPT3[0,0,0,0,0,0,0,0,0,0,0,0,0]},
ResInt=ResInt/.HPT3[0,0,0,0,0,0,0,0,0,0,0,0,0]->1;
];
Return[ResInt/.n_Integer HPT3[arg___]:>Mod[n,$P]HPT3[arg]];
];

(* Verbatim copy of ExtraFunctions.m's ComputeDISIJ, with the diffInts Table
   restructured so PrecomputeTPointData runs once and PrecomputeVarData runs
   once per var instead of once per (var,intpos). Everything else (basis-change
   matrix, IBP reduction, extra relations, output format) is untouched. *)
ComputeDISIJHoisted[samplepoint_]:=Module[{SubReduList,diffIntsMC,ReducedInts,DiffIntegral,ReductionRules,ExtraRelations,tpoint,diffop,kinematicsT,dHPT3All,ReducedDiffrential,NULLS,dHPT3Matrices,dHPT3MatricesSaveFormat,reductionRules,PureIntegralN,diffOp,eps,tPoint,repl,PureIntegralN2,Basis,NullInts,BasisChangMat,BasisChangMatInv,diffInts,diffIntsIBREDUCED,NonReduced,NullInts2,epsrepl,ZeroInts,ZeroIntsRepl,PureIntegralN3,PureIntegralN4,SubLaporta,NonZeros,tPointDataH,varDataH},
If[FileExistsQ[FileNameJoin[{kirapath,"KIRA/splitres/phasespace_00001.m"}]],
	SubReduList=Get[ToString[StringForm["``/KIRA/splitres/phasespace_``.m",kirapath,IntegerString[samplepoint,10,5]]]];,
	SubReduList=Get[ToString[StringForm["``/KIRA/splitres/``/phasespace_``.m",kirapath,IntegerString[Floor[(samplepoint-1)/1000],10,4],IntegerString[samplepoint,10,7]]]];
];
tPoint=Thread[Rule[invariants,SubReduList[[2,1]]]];
reductionRules=ReductionRules=SubReduList[[2,2]];

epsrepl={(ep->(4-d)/2)/.tPoint};
PureIntegralN=Expand[PureIntegralsMaxCut/.tPoint/.(Expand[kinematicFunctions/.tPoint,Modulus->$P])/.epsrepl//.root[n_]:>PowerMod[n,1/2,$P],Modulus->$P];
diffOp=getDifferentialOperator[tPoint];

ExtraRelations={};
ReducedInts=DeleteDuplicates[Cases[reductionRules,_HPT3,Infinity]];
If [DoExtraRelations,
	ExtraRelations=GetExtraRelationReplacements[reductionRules,tPoint];
];
PureIntegralN2=Expand[PureIntegralN/.reductionRules,Modulus->$P];
PureIntegralN3=Expand[PureIntegralN2/.ExtraRelations,Modulus->$P];

NonZeros=Cases[reductionRules,_HPT3,Infinity]//DeleteDuplicates;
ZeroInts=Complement[Ureduce,NonZeros];
ZeroIntsRepl=Thread[Rule[ZeroInts,0]];
PureIntegralN4=Expand[PureIntegralN3/.ExtraRelations/.ZeroIntsRepl,Modulus->$P];

SubLaporta=Cases[PureIntegralN4,_HPT3,Infinity]//DeleteDuplicates;

BasisChangMat=Table[CoefficientArrays[PureIntegralN4[[pos]],SubLaporta],{pos,1,Length[PureIntegralN]}];
BasisChangMat=Normal[BasisChangMat[[All,2]]];
BasisChangMatInv=Inverse[BasisChangMat,Modulus->$P];

SubdVarLIST=Flatten[Map[Position[invariantsnod,#]&,SubdVars]];

tPointDataH=PrecomputeTPointData[tPoint];
diffInts=Table[
varDataH=PrecomputeVarData[diffOp,var,tPoint,tPointDataH];
Table[
DiffIntegral=(PureIntegralsMaxCut[[MaxCutInts[[intpos]]]]/.varDataH["tPointCUT"])/.kinematicFunctions;
ComputeDIFastHoisted[DiffIntegral,tPointDataH,varDataH]
,{intpos,1,Length[MaxCutInts]}]
,{var,SubdVarLIST}];

If[!AllSectors,
diffInts=diffInts/.HPT3[n___]:>OffCutTest[HPT3[n]];
];

diffIntsMC=Join[ConstantArray[0,{10,Length[PureIntegralsMaxCut]-Length[MaxCutInts]}],diffInts,2];

diffIntsIBREDUCED=Expand[diffIntsMC/.reductionRules/.ExtraRelations,Modulus->$P];
NonReduced=Cases[diffIntsIBREDUCED,_HPT3,Infinity];
NullInts2=Complement[NonReduced,SubLaporta];
diffIntsIBREDUCED=diffIntsIBREDUCED/.Thread[Rule[NullInts2,0]];

dHPT3Matrices=Table[Normal[CoefficientArrays[diffIntsIBREDUCED[[i,j]],SubLaporta]/.{0}->{0,ConstantArray[0,Length[PureIntegralN]]}][[2]],{i,1,Length[SubdVarLIST]},{j,1,Length[SubLaporta]}];

For[i=1,i<=Length[SubdVarLIST],i++,dHPT3Matrices[[i]]=Expand[dHPT3Matrices[[i]] . BasisChangMatInv,Modulus->$P][[MaxCutInts,MaxCutInts]];];

dHPT3MatricesSaveFormat=Table[SparseArray[Flatten[dHPT3Matrices[[i]]]],{i,1,Length[SubdVarLIST]}];
If[printOut,Print[ToString[StringForm["Finished `` (hoisted).",samplepoint]]]];
Return[{tPoint[[All,2]],dHPT3MatricesSaveFormat}];
];
