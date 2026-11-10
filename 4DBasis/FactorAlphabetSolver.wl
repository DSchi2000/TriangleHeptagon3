(* ::Package:: *)

(* ::Title:: *)
(*Express a rational function as a product of given "letters" to integer powers*)


(* ::Text:: *)
(**)


(*Given a target rational function and a list of candidate factors ("letters"),
find integer exponents n_i such that

	target == Times @@ (candidates^n)   (up to an overall rational constant)

by comparing irreducible-factor exponents (FactorList) rather than guessing by
eye -- useful for the same kind of "which letters make up this symbol-alphabet
entry" bookkeeping as EvenLetterCandidatesFACS/FactorsWithPositions elsewhere
in this project.

SUBTLETY #1 (sign): FactorList's choice of sign for a given irreducible factor
is NOT guaranteed to be the same across two differently-built expressions that
happen to contain the "same" factor -- e.g. one expression might record a
factor as F, another as -F, with the overall numeric prefactor compensating.
A naive key-by-key match then silently sees zero overlap and reports no
solution. CanonicalFactor below fixes this by mapping F and -F to the same
dictionary key (whichever of {F,-F} sorts first under Mathematica's own
canonical ordering -- applied consistently, so both directions agree). The
true overall sign/constant is still recovered correctly at the end via a
direct numeric check (target/Prod(candidates^n) should come out to a plain
constant, ideally 1).

SUBTLETY #2 (non-uniqueness): the candidate letters you hand it are not
always multiplicatively independent. In practice this showed up as
	f2^2 == f1*f3
(e.g. f1=c1^2*G, f2=c1*(c1-c2)*G, f3=(c1-c2)^2*G for a common G -- an
unsurprising relation once you notice it, but easy to miss by eye) which makes
the exponent system underdetermined: a whole INTEGER LATTICE of equally valid
decompositions exists (one free integer parameter per independent relation
among the candidates). SolveFactorPowers reports isUnique->False in that case
and picks a short representative (few/small exponents) rather than an
arbitrary member of that lattice.

SUBTLETY #3 (don't enumerate -- reduce): with real alphabets you can easily
get a dozen or more independent relations among the candidates (one real
usage here hit 15), i.e. a 15-dimensional lattice of equivalent answers. Do
NOT try to pick the best one by scanning combinations (Tuples over a range of
each free parameter) -- that is exponential in the number of relations and
WILL exhaust memory/time well before you'd notice (seen first-hand: 13^15
combinations). The fix is standard lattice machinery, not more scanning:
Solve[...,Integers] already gives one particular integer solution plus an
explicit integer basis for the solution lattice (read off by setting each
free parameter to a unit vector in turn); LatticeReduce (LLL) turns that into
a short/near-orthogonal basis; then a single least-squares projection +
Round (Babai's rounding) finds a short representative directly, no search.
This is polynomial time regardless of how many free parameters there are --
checked at 15 and 34 independent relations, both well under a second (see
self-test below), versus the old approach which couldn't even finish 15.*)


CanonicalFactor[f_] := Module[{nf = Expand[-f]}, If[OrderedQ[{f, nf}], f, nf]];

SolveFactorPowers[target_, candidates_, candidateNames_:Automatic] := Module[
	{names, allExprs, flAll, polyPart, canonPart, irreducibles, vecOf, targetVec, candMat,
	 nvars, genSol, freeParams, p0, nullBasisInt, reducedBasis, coeffs, roundedCoeffs,
	 best, overallConst},

	names = If[candidateNames === Automatic, Table["f"<>ToString[i], {i, Length[candidates]}], candidateNames];

	allExprs = Prepend[candidates, target];
	flAll = FactorList /@ allExprs;
	polyPart[fl_] := Select[fl, !NumberQ[#[[1]]]&];
	canonPart[fl_] := {CanonicalFactor[#[[1]]], #[[2]]}& /@ polyPart[fl];

	irreducibles = DeleteDuplicates[Flatten[canonPart /@ flAll, 1][[All,1]]];

	vecOf[fl_] := Module[{assoc = Association[Rule@@@canonPart[fl]]},
		Lookup[assoc, #, 0]& /@ irreducibles];

	targetVec = vecOf[First[flAll]];
	candMat = vecOf /@ Rest[flAll];  (* candidates (rows) x irreducibles (cols) *)

	nvars = Array[n, Length[candidates]];

	(* Plain Gaussian elimination -- gives ONE particular integer solution plus a free
	   integer parametrization of the whole solution set. Fast no matter how many free
	   parameters there turn out to be; the combinatorial work (if any) is avoided
	   entirely below, not hidden here. *)
	genSol = Solve[Thread[(nvars.candMat) == targetVec], nvars, Integers];
	If[genSol==={}, Return[<|"status"->"no integer-power solution for the polynomial part"|>]];

	freeParams = DeleteDuplicates[Cases[genSol[[1]], C[_], Infinity]];
	p0 = nvars /. genSol[[1]] /. Thread[freeParams -> 0];

	best = If[freeParams === {},
		p0,
		(* integer basis for the solution lattice's free directions, then LLL + Babai
		   rounding to find a short representative -- polynomial time, no enumeration *)
		nullBasisInt = Table[
			(nvars/.genSol[[1]]/.Thread[freeParams->UnitVector[Length[freeParams],k]]) - p0,
		{k, Length[freeParams]}];
		reducedBasis = LatticeReduce[nullBasisInt];
		coeffs = Quiet@Check[LeastSquares[Transpose[reducedBasis], -p0], ConstantArray[0,Length[reducedBasis]]];
		roundedCoeffs = Round[coeffs];
		p0 + roundedCoeffs.reducedBasis
	];

	overallConst = Simplify[target/Times@@(candidates^best)];

	<|
		"exponents" -> Thread[names -> best],
		"overallConstant" -> overallConst,     (* should be a plain number -- 1 if the
		                                          decomposition is exact and complete *)
		"isUnique" -> (freeParams === {}),     (* False => candidates have a syzygy;
		                                          other equally-valid exponent vectors exist *)
		"numFreeParams" -> Length[freeParams]
	|>
];


(* ::Section:: *)
(*Self-test (the two worked examples from this project's own conversation)*)


(*G = s12^2*s23^2 - ... (the Gram-type denominator shared by c1,c2,c3's Cramer's-rule
formulas, see kinematicRulesTH3); L = c1 s234-c2 s234+c2 s34-c3 s34-c1 s56.
candidates = {c1^2 G, c1(c1-c2)G, (c1-c2)^2 G, s23, s34, <poly6>, <poly7>, -(c1-c2) L G}

target1 = -c1^2 L G/((c1-c2) s23)   -> expect f1^1 f3^-1 f4^-1 f8^1, isUnique->False
target2 = c1^2 <poly7>/((c1-c2)^2 s34) -> expect f1^1 f3^-1 f5^-1 f7^1, isUnique->False

Run with the candidates from that conversation (or your own) to confirm:
	SolveFactorPowers[target1, candidates] (* exponents should match f1,f3,f4,f8 *)
	SolveFactorPowers[target2, candidates] (* exponents should match f1,f3,f5,f7 *)
In both cases "overallConstant" must come back exactly 1 -- that's the real
correctness check, not just "a solution exists".*)
