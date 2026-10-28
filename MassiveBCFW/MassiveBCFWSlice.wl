(* ::Package:: *)

(* ::Title:: *)
(*Massive multi-line holomorphic BCFW slice (legs 5,6 massive, mt2=p5^2=p6^2) -- TriangleHeptagon3*)


(* ::Text:: *)
(**)


(*Builds a univariate slice s(t)=a+b t for the invariants
	{mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56}
(matching ../SetUp/mandelstam_scalarproducts.txt: legs 1,2,3,4 massless,
legs 5,6 massive with p5^2=p6^2=mt2 -- same external kinematics as the ttbb
FullHexaBoxDE project, just a different loop topology/propagator basis here)
that satisfies the 6-point Gram constraint G(p1,p2,p3,p4,p5)=0 identically in t.

Method (see arXiv:2412.19884 eq.(4.30) and arXiv:2110.07541 eq.(4.13)-(4.15)):
an all-line holomorphic shift lambda_k -> lambda_k + t c_k eta (lambda-tilde fixed)
keeps every line exactly null and, once momentum conservation is imposed on the c_k,
keeps the full momentum-conserving configuration exactly 4-dimensional for every t.
The 6-point Gram constraint is then automatic: any 5 vectors spanning a <=4D space
have a singular Gram matrix, with no masslessness assumption anywhere.

Each massive leg (5,6) is represented as the sum of two null constituents,
p5 = q5a+q5b, p6 = q6a+q6b, so mt2 = <5a 5b>[5a 5b] is just one more spinor
bracket -- holding it fixed, or letting it run linearly in t, is then a matter
of one extra *linear* constraint on the shift coefficients c_k (see
GenerateShiftCoefficients below). No square root ever appears: every output
invariant is an exact degree-1 polynomial in t.

All arithmetic is mod the Mersenne prime $P=2^31-1. If you need several primes
(CRT) for one slice, derive the base point / eta / c_k from FIXED RATIONAL
seeds and reduce mod each prime, rather than re-randomizing per prime --
otherwise different primes sample unrelated rational points and cannot be
CRT-combined.*)


(* ::Section:: *)
(*Configuration*)


$P = 2^31 - 1;

(* Base value of mt2 at t=0 (and, in the FixMasses->True case, for all t). *)
MT2SEED = 123456789;

(* True  -> Case A: mt2 held exactly fixed along the slice (build several such
            slices at different MT2SEED values to resolve the mt2-dependence
            separately, as in arXiv:2608.03746 eq.(5.4)).
   False -> Case B (default, recommended): mt2 itself runs linearly in t,
            shared between leg 5 and leg 6, exactly like every other entry
            of diffInvariants. This is the more generic/useful slice since it
            probes the mt2-direction together with the Mandelstam directions
            in a single line, matching the philosophy of arXiv:2412.19884. *)
FixMasses = False;

(*NSamples = 2000;*)
OutputFolder = "UnivariantSlice-1";

SeedRandom[];


(* ::Section:: *)
(*Finite-field spinor-helicity primitives*)


(* Both brackets have the same antisymmetric form; kept as two names for readability. *)
Ang[{l1_,l2_},{m1_,m2_}] := Mod[l1 m2 - l2 m1, $P];
Sq[{l1_,l2_},{m1_,m2_}]  := Mod[l1 m2 - l2 m1, $P];
ModDiv[a_,b_,p_:$P] := Mod[a*PowerMod[b,-1,p], p];
BiSpinor[{l1_,l2_},{t1_,t2_}] := {{l1 t1,l1 t2},{l2 t1,l2 t2}};
RandSpin[] := RandomInteger[{1,$P-1},2];

(* A "line" is {lambda,lambdaTilde}. A composite (massive) momentum is a list
   of lines, e.g. p5 = {q5a,q5b}. Dot2[A,B] = 2 pA.pB for composites A,B. *)
Dot2[A_List,B_List] := Mod[Sum[Ang[a[[1]],b[[1]]]*Sq[b[[2]],a[[2]]], {a,A},{b,B}], $P];
Mass2[lines_List] := ModDiv[Dot2[lines,lines],2];
Sinv[listsOfLines__List] := Mass2[Join[listsOfLines]];
CompositeMat[lines_List] := Mod[Sum[BiSpinor[ln[[1]],ln[[2]]], {ln,lines}], $P];
GramDet[lines_List] := Mod[Det[Table[Dot2[lines[[i]],lines[[j]]], {i,Length[lines]}, {j,Length[lines]}]], $P];


(* ::Section:: *)
(*Base point (t=0): legs 1,2,3,4 massless, legs 5,6 massive with p^2=mt2SeedVal*)


(* Split a target mass into two null constituents sharing a chosen lambda_a:
   solves lambdaTilde_b from <a b>[b a] == targetMass (one linear equation). *)
SplitMassPair[lamA_,ltamA_,lamB_,targetMass_] := Module[{A,K,x0,y0},
	A = Ang[lamA,lamB];
	K = ModDiv[targetMass, A];
	x0 = RandomInteger[{1,$P-1}];
	y0 = ModDiv[Mod[x0*ltamA[[2]] - K, $P], ltamA[[1]]];
	{x0,y0}
];

(* Complete overall momentum conservation p3+p4 == Qmat with p3,p4 null, given
   a chosen lambda3: solves lambdaTilde3 from det(Qmat - lambda3 (x)lambdaTilde3)=0
   (also exactly one linear equation), then reads off p4 from the resulting
   rank-1 remainder. *)
CompleteConservation[Qmat_,lam3_] := Module[{a,b,c,d,rhs,coefX,coefY,x3,y3,ltam3,M3,M4,lam4,ltam4},
	{{a,b},{c,d}} = Qmat;
	rhs = Mod[a d - b c, $P];
	coefX = Mod[d*lam3[[1]] - b*lam3[[2]], $P];
	coefY = Mod[a*lam3[[2]] - c*lam3[[1]], $P];
	x3 = RandomInteger[{1,$P-1}];
	y3 = ModDiv[Mod[rhs - coefX*x3, $P], coefY];
	ltam3 = {x3,y3};
	M3 = BiSpinor[lam3,ltam3];
	M4 = Mod[Qmat - M3, $P];
	lam4 = {M4[[1,1]], M4[[2,1]]};
	ltam4 = {1, ModDiv[M4[[1,2]], M4[[1,1]]]};
	{ltam3, lam4, ltam4}
];

GenerateBasePoint[mt2Val_] := Module[
	{lam1,ltam1,lam2,ltam2,lam3,ltam3,lam4,ltam4,
	 lam5a,ltam5a,lam5b,ltam5b,lam6a,ltam6a,lam6b,ltam6b,
	 Qmat,res},
	lam1=RandSpin[]; ltam1=RandSpin[];
	lam2=RandSpin[]; ltam2=RandSpin[];
	lam5a=RandSpin[]; ltam5a=RandSpin[]; lam5b=RandSpin[];
	ltam5b = SplitMassPair[lam5a,ltam5a,lam5b,mt2Val];
	lam6a=RandSpin[]; ltam6a=RandSpin[]; lam6b=RandSpin[];
	ltam6b = SplitMassPair[lam6a,ltam6a,lam6b,mt2Val];
	lam3 = RandSpin[];
	(* all-outgoing convention: p1+p2+p3+p4+p5+p6=0  =>  p3+p4 = Q *)
	Qmat = Mod[-(CompositeMat[{{lam1,ltam1}}] + CompositeMat[{{lam2,ltam2}}]
	            + CompositeMat[{{lam5a,ltam5a},{lam5b,ltam5b}}]
	            + CompositeMat[{{lam6a,ltam6a},{lam6b,ltam6b}}]), $P];
	res = CompleteConservation[Qmat,lam3];
	{ltam3,lam4,ltam4} = res;
	<|
		1 -> {lam1,ltam1}, 2 -> {lam2,ltam2}, 3 -> {lam3,ltam3}, 4 -> {lam4,ltam4},
		"5a" -> {lam5a,ltam5a}, "5b" -> {lam5b,ltam5b},
		"6a" -> {lam6a,ltam6a}, "6b" -> {lam6b,ltam6b}
	|>
];


(* ::Section:: *)
(*Shift coefficients: all-line holomorphic shift lambda_k -> lambda_k + t c_k eta*)


KEYS8 = {1,2,3,4,"5a","5b","6a","6b"};

(* Solve a homogeneous linear system (list of expressions in cvars, each ==0)
   by picking all but Length[eqs] of the unknowns randomly and solving the rest. *)
SolveShiftSystem[eqs_List, cvars_List] := Module[
	{nEq,freeVars,solveVars,freeVals,subsFree,eqsExpand,constVec,coeffMat,solVals},
	nEq = Length[eqs];
	freeVars = cvars[[;;-nEq-1]];
	solveVars = cvars[[-nEq;;]];
	freeVals = RandomInteger[{1,$P-1}, Length[freeVars]];
	subsFree = Thread[freeVars -> freeVals];
	eqsExpand = Expand[eqs /. subsFree];
	{constVec,coeffMat} = CoefficientArrays[eqsExpand, solveVars];
	coeffMat = Mod[Normal[coeffMat], $P];
	constVec = Mod[Normal[constVec], $P];
	solVals = LinearSolve[coeffMat, Mod[-constVec,$P], Modulus->$P];
	Join[subsFree, Thread[solveVars -> solVals]]
];

GenerateShiftCoefficients[baseLines_Association, eta_, fixMasses_:False] := Module[
	{lam,ltam,cvars,momCons1,momCons2,sq5,sq6,massEqs,eqs,subs},
	lam[k_] := baseLines[k][[1]];
	ltam[k_] := baseLines[k][[2]];
	cvars = {c1,c2,c3,c4,c5a,c5b,c6a,c6b};
	momCons1 = Sum[cvars[[i]]*ltam[KEYS8[[i]]][[1]], {i,8}];
	momCons2 = Sum[cvars[[i]]*ltam[KEYS8[[i]]][[2]], {i,8}];
	sq5 = Sq[ltam["5a"],ltam["5b"]];
	sq6 = Sq[ltam["6a"],ltam["6b"]];
	massEqs = If[fixMasses,
		{c5a*Ang[eta,lam["5b"]] + c5b*Ang[lam["5a"],eta],
		 c6a*Ang[eta,lam["6b"]] + c6b*Ang[lam["6a"],eta]},
		{sq5*(c5a*Ang[eta,lam["5b"]] + c5b*Ang[lam["5a"],eta])
		 - sq6*(c6a*Ang[eta,lam["6b"]] + c6b*Ang[lam["6a"],eta])}
	];
	eqs = Join[{momCons1,momCons2}, massEqs];
	subs = SolveShiftSystem[eqs, cvars];
	AssociationThread[KEYS8, cvars /. subs]
];


(* ::Section:: *)
(*Evaluating the slice*)


LineAtT[baseLines_,coeffs_,eta_,k_,tval_] := {
	Mod[baseLines[k][[1]] + tval*coeffs[k]*eta, $P],
	baseLines[k][[2]]
};

InvariantsAtT[baseLines_,coeffs_,eta_,tval_] := Module[{p1,p2,p3,p4,p5,p6},
	p1 = {LineAtT[baseLines,coeffs,eta,1,tval]};
	p2 = {LineAtT[baseLines,coeffs,eta,2,tval]};
	p3 = {LineAtT[baseLines,coeffs,eta,3,tval]};
	p4 = {LineAtT[baseLines,coeffs,eta,4,tval]};
	p5 = {LineAtT[baseLines,coeffs,eta,"5a",tval], LineAtT[baseLines,coeffs,eta,"5b",tval]};
	p6 = {LineAtT[baseLines,coeffs,eta,"6a",tval], LineAtT[baseLines,coeffs,eta,"6b",tval]};
	<|
		"mt2" -> Mass2[p5], "mt2b" -> Mass2[p6],
		"s12" -> Sinv[p1,p2], "s23" -> Sinv[p2,p3], "s34" -> Sinv[p3,p4],
		"s45" -> Sinv[p4,p5], "s56" -> Sinv[p5,p6], "s16" -> Sinv[p1,p6],
		"s123" -> Sinv[p1,p2,p3], "s234" -> Sinv[p2,p3,p4], "s345" -> Sinv[p3,p4,p5],
		"gram15" -> GramDet[{p1,p2,p3,p4,p5}]
	|>
];

(* diffInvariants order -- matches ../SetUp/mandelstam_scalarproducts.txt's variable set *)
diffInvariants = {mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56};
InvKeysOrdered = {"mt2","s12","s123","s16","s23","s234","s34","s345","s45","s56"};


(* ::Section:: *)
(*Callable entry point*)


(* BuildMassiveSlice[mt2Val, fixMasses] generates a fresh random slice and
   returns an Association:
     "a", "b"        -> the a_k,b_k vectors (Associations over InvKeysOrdered
                        plus "mt2b"), with invariant_k(t) = a_k + b_k t
     "basePoint", "eta", "coeffs" -> the underlying data, in case you want to
                        evaluate InvariantsAtT[basePoint,coeffs,eta,tval]
                        at more t-values later without rebuilding the slice.
   Call it directly, e.g.:
     res = BuildMassiveSlice[123456789, False];
     res["a"]   (* aVec *)
     res["b"]   (* bVec *) *)
BuildMassiveSlice[mt2Val_:MT2SEED, fixMasses_:FixMasses] := Module[
	{bp,et,cf,i0,i1,aV,bV},
	bp = GenerateBasePoint[mt2Val];
	et = RandSpin[];
	cf = GenerateShiftCoefficients[bp, et, fixMasses];
	i0 = InvariantsAtT[bp,cf,et,0];
	i1 = InvariantsAtT[bp,cf,et,1];
	aV = Association[# -> i0[#] & /@ Append[InvKeysOrdered,"mt2b"]];
	bV = Association[# -> Mod[i1[#]-i0[#],$P] & /@ Append[InvKeysOrdered,"mt2b"]];
	<|"a"->aV, "b"->bV, "basePoint"->bp, "eta"->et, "coeffs"->cf|>
];


(* ::Section:: *)
(*Build one slice and self-test it*)


(*sliceResult = BuildMassiveSlice[MT2SEED, FixMasses];
basePoint = sliceResult["basePoint"];
eta = sliceResult["eta"];
coeffs = sliceResult["coeffs"];
aVec = sliceResult["a"];
bVec = sliceResult["b"];

Print["a vector (value at t=0): ", aVec];
Print["b vector (slope):        ", bVec];

(* Self-test: linearity in t (a 3rd, random t must match the line through 0,1)
   and the Gram constraint vanishing identically -- both are consequences of
   the construction, not something separately imposed, so any failure here
   means a coding bug, not a physics issue. *)
Module[{tCheck,invCheck,ok,predicted},
	tCheck = RandomInteger[{2,$P-1}];
	invCheck = InvariantsAtT[basePoint,coeffs,eta,tCheck];
	ok = True;
	Do[
		predicted = Mod[aVec[k] + tCheck*bVec[k], $P];
		If[predicted =!= invCheck[k], ok = False;
			Print["LINEARITY FAILURE on ", k, ": predicted ", predicted, " got ", invCheck[k]]];
	, {k, InvKeysOrdered}];
	If[invCheck["gram15"] =!= 0, ok = False; Print["GRAM CONSTRAINT FAILURE: ", invCheck["gram15"]]];
	If[FixMasses,
		If[invCheck["mt2"] =!= MT2SEED, ok = False; Print["mt2 not fixed: ", invCheck["mt2"]]],
		If[invCheck["mt2"] =!= invCheck["mt2b"], ok = False; Print["mt2/mt2b mismatch: ", invCheck["mt2"]," vs ",invCheck["mt2b"]]]
	];
	If[ok, Print["Self-test passed: all invariants linear in t, Gram(1,2,3,4,5)=0 identically, mt2 handled as requested (FixMasses=",FixMasses,")."]];
];*)


(* ::Section:: *)
(*Export*)


(*(* Column layout matches diffInvariants={mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56};
   column 1 is the slice parameter t itself. One row per sample. *)

SetDirectory[DirectoryName[$InputFileName]];
If[!DirectoryQ[OutputFolder], CreateDirectory[OutputFolder]];

(* a_k, b_k (one "key -> value" per line) so the slice can be reproduced/checked
   later without rerunning the random generation. *)
aKeysSaved = Append[InvKeysOrdered,"mt2b"];
Export[FileNameJoin[{OutputFolder,"aVector.txt"}],
	StringRiffle[(ToString[#]<>" -> "<>ToString[aVec[#]])& /@ aKeysSaved, "\n"], "Text"];
Export[FileNameJoin[{OutputFolder,"bVector.txt"}],
	StringRiffle[(ToString[#]<>" -> "<>ToString[bVec[#]])& /@ aKeysSaved, "\n"], "Text"];

tSamples = RandomInteger[{1,$P-1}, NSamples]; (* collision prob. negligible for NSamples << Sqrt[$P] *)

ppLines = Table[
	inv = InvariantsAtT[basePoint,coeffs,eta,tv];
	Join[{tv}, inv /@ InvKeysOrdered]
, {tv,tSamples}];

Export[FileNameJoin[{OutputFolder,"phasespacepoints.txt"}],
	StringRiffle[ToString[#,InputForm]& /@ ppLines, "\n"], "Text"];
Export[FileNameJoin[{OutputFolder,"tvars.txt"}],
	StringRiffle[("tvar -> "<>ToString[#])& /@ tSamples, "\n"], "Text"];

Print["Wrote ", Length[tSamples], " samples (one row each) to ", FileNameJoin[{OutputFolder,"phasespacepoints.txt"}], " and tvars.txt"];
Print["Wrote a,b vectors to ", FileNameJoin[{OutputFolder,"aVector.txt"}], " and ", FileNameJoin[{OutputFolder,"bVector.txt"}]];
Print["basePoint, eta, coeffs are in-session if you want to generate more samples of the SAME slice (same a,b) without rebuilding it."];*)
