(* ::Package:: *)

pathDO=DirectoryName[$InputFileName];
kinematics=ReadList[FileNameJoin[{HomeDir,"SetUp/mandelstam_kinematics.txt"}]];
mandelstamScalarProducts=ReadList[FileNameJoin[{HomeDir,"SetUp/mandelstam_scalarproducts.txt"}]];
mandelstamScalarProducts=Join[kinematics,mandelstamScalarProducts]
invariantsnod={mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56};
invariantsnodUpperCase={Mt2,S12,S123,S16,S23,S234,S34,S345,S45,S56};

p4conv={p[6]->-(p[1]+p[2]+p[3]+p[4]+p[5])};
$P=2^31-1;


(* ::Title:: *)
(*Differential Operator Ansatz*)


diffOpAnsatz[expr_] := Module[{v1},
	  v1 = Sum[da[j, i]*p[j]*D[expr, p[i]], {i, 1, 6}, {j, 1, 4}];
	  v1 = v1 + Sum[ds[invariantsnodUpperCase[[i]]]*D[expr, invariantsnod[[i]]], {i, 1, Length[invariantsnod]}];
	  Return[v1];
	  ];


diffOpAnsatz2 =Sum[da[j, i]*p[j]*del[p[i]], {i, 1, 6}, {j, 1, 4}] + Sum[ds[invariantsnodUpperCase[[i]]]*del[invariantsnod[[i]]], {i, 1, Length[invariantsnod]}];


(* ::Title:: *)
(*Constraints*)


operatorConstraints={
s12-(p[1]+p[2])^2,
s23-(p[2]+p[3])^2,
s34-(p[3]+p[4])^2,
s45-(p[4]+p[5])^2,
s56-(p[5]+p[6])^2,
s16-(p[1]+p[6])^2,
s123-(p[1]+p[2]+p[3])^2,
s234-(p[2]+p[3]+p[4])^2,
s345-(p[3]+p[4]+p[5])^2,
(p[1])^2,
(p[2])^2,
(p[3])^2,
(p[4])^2,
mt2-(p[5])^2,
mt2-(p[6])^2
};
diffOpConst=Expand[Thread[diffOpAnsatz[operatorConstraints]]/.p4conv];
p4diffopConst=Table[Expand[p[i]diffOpAnsatz[p[1]+p[2]+p[3]+p[4]+p[5]+p[6]]/.p4conv],{i,1,4}];
EQNSystemDiffOp=Join[diffOpConst,p4diffopConst];
solsetDiffOp=Cases[EQNSystemDiffOp,_da,Infinity]//DeleteDuplicates;
solsetDiffOp={da[1,1],da[1,5],da[1,6],da[2,1],da[2,2],da[2,4],da[2,5],da[2,6],da[3,1],da[3,4],da[3,5],da[3,6],da[4,1],da[4,2],da[4,3],da[4,4],da[4,5],da[4,6]}
(*solsetDiffOp={da[1,1],da[1,5],da[1,6],da[2,1],da[3,3],da[2,4],da[3,2],da[1,3],da[3,1],da[3,4],da[3,5],da[1,4],da[4,1],da[4,2],da[4,3],da[2,3],da[4,5],da[1,2]}*)
solsetDiffOp=Join[solsetDiffOp,{ds[S16]}]
InternalAuxvar=S16;
changeAuxVar[var_] := Module[{},
	(* Keep solsetDiffOp to just the 18 "da" entries, regardless of history --
	   getDifferentialOperator reconstructs the correct aux-specific solve list itself
	   from InternalAuxvar every time it's called, so solsetDiffOp must never be left
	   carrying a *stale* ds[...] entry (e.g. from a previous changeAuxVar call) that
	   silently desyncs it from InternalAuxvar. The old S16->var substitution only ever
	   worked once: after the first switch away from S16, the literal symbol S16 is gone
	   from solsetDiffOp, so later calls (including switching back to S16) were silent
	   no-ops, leaving solsetDiffOp and InternalAuxvar out of sync. *)
	solsetDiffOp=Cases[solsetDiffOp,_da];
	InternalAuxvar=ToExpression[Decapitalize[ToString[var]]]
];


(* ::Title:: *)
(*Solve Linear System*)


getDifferentialOperator[tpoint_] := Module[{EQNSystemDiffOpT,solDiffOp,kinematicsT,mandelstamScalarProductsT,
	NullsetDiffOp,NullRulesDiffOp,differentialOperator,SixPointGramLocal,totalDiffG6,momentaLocal,
	auxUpper,auxLower,posAux,auxEqnPos,EQNSystemDiffOpFixedT,solsetWithAux},

	(* identify the auxiliary invariant in both upper/lower case, regardless of how changeAuxVar left InternalAuxvar *)
	posAux=Position[invariantsnod,InternalAuxvar];
	If[Length[posAux]>0,
		auxLower=InternalAuxvar;
		auxUpper=invariantsnodUpperCase[[posAux[[1,1]]]],
		posAux=Position[invariantsnodUpperCase,InternalAuxvar];
		auxUpper=InternalAuxvar;
		auxLower=invariantsnod[[posAux[[1,1]]]]
	];

	kinematicsT=Expand[kinematics/.tpoint,Modulus->$P];
	mandelstamScalarProductsT=Expand[mandelstamScalarProducts/.tpoint,Modulus->$P];
	EQNSystemDiffOpT=Expand[EQNSystemDiffOp/.kinematicsT/.mandelstamScalarProductsT,Modulus->$P];

	(* The actual D=4 (SixPointGram=0) tangency condition at this phase-space point. *)
	momentaLocal={p[1],p[2],p[3],p[4],p[5]};
	SixPointGramLocal=Expand[Det[Table[momentaLocal[[i]] momentaLocal[[j]],{i,1,5},{j,1,5}]/.kinematics/.mandelstamScalarProducts/.loopMomentumScalarProducts]];
	totalDiffG6=Expand[Sum[ds[invariantsnodUpperCase[[i]]]*(D[SixPointGramLocal,invariantsnod[[i]]]/.tpoint),{i,Length[invariantsnod]}],Modulus->$P];

	(* Replace the auxiliary invariant's own (wrong) momentum-only definition inside
	   operatorConstraints with this true G6 condition, THEN solve the full system (18 da's +
	   ds[auxUpper]) together. Just excluding ds[auxUpper] from the solve and patching it in
	   afterwards is NOT enough: the 19-equation system restricted to 18 "da" unknowns is still
	   one equation short regardless of which ds[X] you leave out, so Solve surfaces that missing
	   equation as a hidden condition on some OTHER ds[X] (whichever one happens to pivot) --
	   silently breaking that invariant's own calibration instead. Swapping in the G6 condition
	   in place of the auxiliary's defining equation keeps the system a clean, fully-determined
	   19 equations, so every one of the other 9 directions comes out both tangent to G6 AND
	   correctly calibrated (d(own definition)/d(self)=1, d(others)/d(self)=0). *)
	auxEqnPos=Select[Range[Length[operatorConstraints]],Coefficient[operatorConstraints[[#]],auxLower]===1&];
	If[Length[auxEqnPos]!=1,
		Print["getDifferentialOperator: could not uniquely locate ",auxLower,"'s defining equation in operatorConstraints (found at ",auxEqnPos,"). This auxiliary-variable choice isn't supported -- pick one of the nine S-type invariants, not mt2 (which appears in two mass-shell equations)."];
		Return[$Failed]
	];
	EQNSystemDiffOpFixedT=EQNSystemDiffOpT;
	EQNSystemDiffOpFixedT[[auxEqnPos[[1]]]]=totalDiffG6;

	solsetWithAux=Join[DeleteCases[solsetDiffOp,ds[auxUpper]],{ds[auxUpper]}];
	solDiffOp=Solve[Thread[EQNSystemDiffOpFixedT==0],solsetWithAux,Modulus->$P][[1]];
	NullsetDiffOp=Cases[solDiffOp[[All,2]],_da,Infinity]//DeleteDuplicates;
	NullRulesDiffOp=Thread[Rule[NullsetDiffOp,0]];
	solDiffOp=Join[solDiffOp/.NullRulesDiffOp,NullRulesDiffOp];
	differentialOperator=diffOpAnsatz2/.solDiffOp;

	Return[differentialOperator];
];


(* ::Title:: *)
(*Test Differential Operator*)


checkDifferentialOperator[tpoint_,diffop_] := Module[{diffopS,invariantsnod2,InvariantsConstraintsCheck,InvariantsMandsCheck,CheckMatrix,kinematicsT,mandelstamScalarProductsT},
	(*invariantsnod2={s12,s23,s34,s45,s56,s16,s123,s234,s345,mt2};*)
	invariantsnod2={mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56};
	kinematicsT=Expand[kinematics/.tpoint,Modulus->$P];
	mandelstamScalarProductsT=Expand[mandelstamScalarProducts/.tpoint,Modulus->$P];
	InvariantsConstraintsCheck={(p[1]+p[2])^2,(p[2]+p[3])^2,(p[3]+p[4])^2,(p[4]+p[5])^2,(p[5]+p[6])^2,(p[1]+p[6])^2,(p[1]+p[2]+p[3])^2,(p[2]+p[3]+p[4])^2,(p[3]+p[4]+p[5])^2,(p[1])^2,
	(p[2])^2,(p[3])^2,(p[4])^2,(p[5])^2,(p[6])^2};
	InvariantsMandsCheck={s12,s23,s34,s45,s56,s16,s123,s234,s345,0,0,0,mt2,0,mt2};
	For[i=1,i<=Length[invariantsnod2],i++,
		diffopS[invariantsnod2[[i]]]=Coefficient[Expand[diffop],ds[invariantsnodUpperCase[[i]]]];
	];
	CheckMatrix=Table[diffopS[invariantsnod2[[j]]]/.del[x_]:>D[InvariantsConstraintsCheck[[i]],x],{i,1,Length[InvariantsConstraintsCheck]},{j,1,Length[invariantsnod2]}];
	CheckMatrix=Expand[Expand[CheckMatrix/.p4conv]/.kinematicsT/.mandelstamScalarProductsT,Modulus->$P];
	MatrixPlot[CheckMatrix,GridLines->Automatic,FrameTicks -> 
   {{Transpose[ArrayReshape[Append[Range[1,15],InvariantsMandsCheck]//Flatten,{2,15}]], Transpose[ArrayReshape[Append[Range[1,15],InvariantsMandsCheck]//Flatten,{2,15}]]}, {Transpose[ArrayReshape[Append[Range[1,Length[invariantsnod2]],invariantsnod2]//Flatten,{2,Length[invariantsnod2]}]], Transpose[ArrayReshape[Append[Range[1,Length[invariantsnod2]],invariantsnod2]//Flatten,{2,Length[invariantsnod2]}]]}}]//Echo;

(*	loopMomentumScalarProducts=Get[FileNameJoin[{Directory[],"files/loopmomentum_scalarproducts.txt"}]];*)
	momenta={p[1],p[2],p[3],p[4],p[5]};
	SixPointGram=Expand[Det[Table[momenta[[i]] momenta[[j]],{i,1,5},{j,1,5}]/.kinematics/.mandelstamScalarProducts/.loopMomentumScalarProducts]];
	
	SixPointGram=Expand[Det[Table[momenta[[i]] momenta[[j]],{i,1,5},{j,1,5}]/.kinematics/.mandelstamScalarProducts/.loopMomentumScalarProducts]];	
	Print[StringForm["\[CapitalDelta]6=``",Expand[SixPointGram/.tpoint,Modulus->$P]]];
	For[i=1,i<=Length[invariantsnod2],i++,
		dSixPointGram=Expand[((diffopS[invariantsnod2[[i]]]/.del[x_]:>D[SixPointGram,x])/.p4conv//Expand)/.kinematics/.mandelstamScalarProducts/.tpoint,Modulus->$P];	
		Print[StringForm["\!\(\*SubscriptBox[\(\[PartialD]\), \(``\)]\)\[CapitalDelta]6=``",invariantsnod2[[i]],dSixPointGram]];
	];
	
]
