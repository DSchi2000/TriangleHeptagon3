(* ::Package:: *)

(* ::Title:: *)
(*Square-root-free parametrization of the 6-point (4 massless + 2 equal-mass) phase space -- TriangleHeptagon3*)


(* ::Text:: *)
(**)


(*Why s16 needs a root in the s_ij-only description:
	diffInvariants = {mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56}
satisfy exactly one relation, the 6-point Gram constraint G6(diffInvariants)=0
(../SetUp/SixPointGram.txt, = Gram(p1,p2,p3,p4,p5) up to the row operation that
trades p6=-(p1+..+p5) for -p1, see the derivation below). G6 is quadratic in
s16, so "solve G6=0 for s16" necessarily produces a square root -- that root is
not optional, it is how the 10 invariants ever manage to describe a genuinely
4-dimensional kinematic configuration.

Fix: stop using s16 (and mt2,s45,s345,s56) as independent coordinates. Use
p1,p2,p3,p4 themselves as a basis of the 4D momentum space (that's what "4D
Gram constraint" means: any 5th momentum MUST be their linear combination),
and write
	p5 = c1 p1 + c2 p2 + c3 p3 + c4 p4                                 [star]
p6 = -(p1+p2+p3+p4+p5) then follows from momentum conservation, linearly, no
root. G6=0 (Gram(p1,...,p5)=0) now holds IDENTICALLY for every c1,c2,c3,c4,
because [star] makes p1,...,p5 linearly dependent by construction -- it is a
statement about ANY 5 vectors spanning a <=4D space, not a constraint you
impose afterwards.

Free parameters (9, matching the 9 independent invariants hiding inside the
10 diffInvariants + G6=0):
	s12, s23, s34, s123, s234   (already your names -- fix p1,p2,p3,p4's "shape")
	s14 := 2 p1.p4              (a 6th invariant of p1..p4 that was previously
	                             hidden inside the s16 root -- there is no way
	                             around introducing it, 4 massless momenta need
	                             6 numbers: s12,s13,s14,s23,s24,s34, and
	                             s123,s234 only fix s13,s24, not s14)
	c1, c2, c3                  (3 of the 4 coefficients in [star])
c4 is FIXED, linearly, by the one remaining physical condition that is not
automatic: p6^2 = p5^2 ( = mt2 ). Unlike G6=0, this really is an extra
condition (nothing forces leg 6 to have the same mass as leg 5 just because
p5 lies in span(p1..p4)) -- but it is manifestly LINEAR in c1..c4 (the
quadratic c_i c_j pieces of p5^2 and p6^2 cancel identically), so it is
solved by ordinary division, never a square root. See VERIFICATION below for
the proof (symbolic identity + 10/10 random finite-field points reproducing
G6=0 exactly).

mt2, s45, s345, s16 are then the (rational, root-free) images of these 9
parameters; s56 comes out purely linear, independent of c1,c2,c3,c4:
	s56 = s123 + s14 - s23 + s234
(this is just (p1+p2+p3+p4)^2 = (p5+p6)^2, fixed once p1..p4 are fixed). *)


(* ::Section:: *)
(*Master formulas (p1..p4 treated as a basis: pi.pi=0, pi.pj free)*)


(* dot products of the p1..p4 basis, in terms of your existing names plus the
   one extra invariant s14 this parametrization needs: *)
basisDots = {
	d[1,2]->s12/2, d[2,3]->s23/2, d[3,4]->s34/2,
	d[1,3]->(s123-s12-s23)/2, d[2,4]->(s234-s23-s34)/2,
	d[1,4]->s14/2
};
d[i_,i_] := 0;
d[i_,j_] /; i>j := d[j,i];

cvec = {c1,c2,c3,c4};

(* p5.pk for k=1..4, and p5.p5=mt2, from p5 = Sum ck pk *)
dot5[k_] := Sum[cvec[[i]] d[i,k], {i,4}] /. basisDots // Expand;
mt2Formula = Expand[Sum[cvec[[i]] cvec[[j]] d[i,j], {i,4}, {j,4}] /. basisDots];

(* s45=(p4+p5)^2, s345=(p3+p4+p5)^2, s16=(p1+p6)^2=(p2+p3+p4+p5)^2 -- each just
   "(known massless sum)^2 + 2(known).p5 + p5^2", nothing but Expand: *)
s45Formula  = Expand[mt2Formula + 2 dot5[4]];
s345Formula = Expand[s34 + 2(dot5[3]+dot5[4]) + mt2Formula];
s16Formula  = Expand[s234 + 2(dot5[2]+dot5[3]+dot5[4]) + mt2Formula];
s56Formula  = Expand[s123 + s14 - s23 + s234]; (* c-independent, see header *)

(* The one real physical constraint -- p6^2=p5^2 -- linear in c1..c4 (the
   quadratic c_i c_j pieces cancel between p6^2 and p5^2): p6^2-p5^2
   = (p1+p2+p3+p4)^2 + 2 p5.(p1+p2+p3+p4), manifestly linear in c1..c4. *)
massEqualityConstraint = Expand[s56Formula + 2 Sum[dot5[k],{k,4}]];

(* c4 purely in terms of the 9 free variables (no mt2 needed -- see header:
   massEqualityConstraint never involves mt2 to begin with): *)
c4Formula9 = c4 /. Solve[massEqualityConstraint==0, c4][[1]];


(* ::Section:: *)
(*Collapse to exactly 9 variables, no G6, no c4*)


(*If your own expressions (differential equations, master-integral bases,
anything built from the diffInvariants) are written in terms of the
redundant 10 -- mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56 -- apply
ReduceTo9Variables to them: it eliminates mt2,s16,s45,s345,s56 (and c4, if
present) in favor of the 9 truly free, unconstrained symbols
	s12, s23, s34, s123, s234, s14, c1, c2, c3
with NO relation left to satisfy anywhere -- not G6=0, not
massEqualityConstraint, nothing. Checked: Together[G6/.ReduceTo9Variables]
is the zero polynomial, not just zero on samples.

Usage:
	myExpr /. ReduceTo9Variables
(apply once; every RHS below is already fully reduced to the 9 variables,
so there is nothing left for a second pass to do).*)

mt2Reduced  = Together[mt2Formula  /. c4->c4Formula9];
s45Reduced  = Together[s45Formula  /. c4->c4Formula9];
s345Reduced = Together[s345Formula /. c4->c4Formula9];
s16Reduced  = Together[s16Formula  /. c4->c4Formula9];
s56Reduced  = s56Formula; (* already c-independent, see header *)

ReduceTo9Variables = {
	mt2  -> mt2Reduced,
	s16  -> s16Reduced,
	s45  -> s45Reduced,
	s345 -> s345Reduced,
	s56  -> s56Reduced
};


(* ::Section:: *)
(*Alternative: keep mt2, drop s12 instead*)


(*c4 was NEVER one of the 9 free variables above -- it's already eliminated
by massEqualityConstraint, which has nothing to do with mt2. So "swap mt2 in
for c4" isn't actually an available move; to keep mt2 you have to give up one
of the 9 that WERE free. Checking the degree of mt2Reduced in each of them:
	c1,c2,c3,s14,s23,s34,s234 -> degree 2   (solving for any of these given
	                                          mt2 needs a square root -- same
	                                          quadratic-Gram structure that
	                                          caused the original s16 problem)
	s12, s123                 -> degree 1   (safe: solving for either, given
	                                          mt2, is a plain linear/rational
	                                          step, no root)
So mt2 can replace s12 (done below) or s123 (symmetric, not included) in the
free list, but nothing else -- not c1,c2,c3,s14, and certainly not c4, which
was never a candidate. Free variables here:
	mt2, s23, s34, s123, s234, s14, c1, c2, c3
with s12 (and c4, via c4Formula9b below) now derived. Checked the same way
as above: Together[G6/.ReduceTo9VariablesKeepMt2] is the zero polynomial.*)

s12Reduced = Together[s12 /. Solve[mt2Reduced == mt2, s12][[1]]];
c4Formula9b = Together[c4Formula9 /. s12->s12Reduced];

s45Reduced2  = Together[s45Formula  /. {c4->c4Formula9b, s12->s12Reduced}];
s345Reduced2 = Together[s345Formula /. {c4->c4Formula9b, s12->s12Reduced}];
s16Reduced2  = Together[s16Formula  /. {c4->c4Formula9b, s12->s12Reduced}];
s56Reduced2  = s56Formula; (* no s12 or c dependence to begin with *)

(* Flat kinematic replacement list for THIS parametrization -- the mt2-kept
   analogue of kinematicRulesTH3. q14=s14/2 directly now (s14 is itself one
   of the 9 free variables here, so there's nothing to solve for it); q12,q13
   use the derived s12Reduced. Free variables feeding every RHS below:
   mt2,s23,s34,s123,s234,s14,c1,c2,c3. *)
kinematicRulesTH3KeepMt2 = {
	q12 -> Together[s12Reduced/2],
	q23 -> s23/2,
	q34 -> s34/2,
	q13 -> Together[(s123-s12Reduced-s23)/2],
	q24 -> (s234-s23-s34)/2,
	q14 -> s14/2,
	s12  -> s12Reduced,
	s16  -> s16Reduced2,
	s45  -> s45Reduced2,
	s345 -> s345Reduced2,
	s56  -> s56Reduced2,
	c4   -> c4Formula9b
};

ReduceTo9VariablesKeepMt2 = {
	s12  -> s12Reduced,
	s16  -> s16Reduced2,
	s45  -> s45Reduced2,
	s345 -> s345Reduced2,
	s56  -> s56Reduced2
};


(* ::Section:: *)
(*Alternative: keep mt2 AND s12, drop s123 instead*)


(*s14 itself is degree 2 in mt2Reduced (same quadratic-Gram obstruction as
c1,c2,c3,s23,s34,s234 above) -- it can never be the variable you drop in
favor of mt2, no matter which other invariant you'd like to keep free instead.
The only two safe (degree-1) slots are s12 and s123 -- pick exactly one to
drop. To keep s12 free, drop s123 instead. Free variables here:
	mt2, s12, s23, s34, s234, s14, c1, c2, c3
with s123 (and c4, via c4Formula9c below) now derived. Checked the same way
as above: Together[G6/.ReduceTo9VariablesKeepMt2AndS12] is the zero
polynomial.*)

s123Reduced = Together[s123 /. Solve[mt2Reduced == mt2, s123][[1]]];
c4Formula9c = Together[c4Formula9 /. s123->s123Reduced];

s45Reduced3  = Together[s45Formula  /. {c4->c4Formula9c, s123->s123Reduced}];
s345Reduced3 = Together[s345Formula /. {c4->c4Formula9c, s123->s123Reduced}];
s16Reduced3  = Together[s16Formula  /. {c4->c4Formula9c, s123->s123Reduced}];
s56Reduced3  = Together[s56Formula  /. {s123->s123Reduced}];

kinematicRulesTH3KeepMt2AndS12 = {
	q12 -> s12/2,
	q23 -> s23/2,
	q34 -> s34/2,
	q13 -> Together[(s123Reduced-s12-s23)/2],
	q24 -> (s234-s23-s34)/2,
	q14 -> s14/2,
	s123 -> s123Reduced,
	s16  -> s16Reduced3,
	s45  -> s45Reduced3,
	s345 -> s345Reduced3,
	s56  -> s56Reduced3,
	c4   -> c4Formula9c
};

ReduceTo9VariablesKeepMt2AndS12 = {
	s123 -> s123Reduced,
	s16  -> s16Reduced3,
	s45  -> s45Reduced3,
	s345 -> s345Reduced3,
	s56  -> s56Reduced3
};


(* ::Section:: *)
(*Alternative: use s56 instead of s14 as the 6th shape invariant*)


(*s14 and s56 are linearly related -- s56 = s123+s14-s23+s234, i.e.
s14 = s56-s123+s23-s234 -- so trading one for the other in the free list
costs nothing (no root, just substituting a linear relation everywhere). The
appeal: s56 is already one of your original named invariants, so this
version needs no new symbol at all. Free variables here:
	s12, s23, s34, s123, s234, s56, c1, c2, c3
with s14, mt2, s16, s45, s345 (and c4) now derived. Checked the same way as
above: Together[G6/.ReduceTo9VariablesS56] is the zero polynomial.*)

s14FormulaS56 = s56-s123+s23-s234;

mt2FormulaS56  = Expand[mt2Formula  /. s14->s14FormulaS56];
s45FormulaS56  = Expand[s45Formula  /. s14->s14FormulaS56];
s345FormulaS56 = Expand[s345Formula /. s14->s14FormulaS56];
s16FormulaS56  = Expand[s16Formula  /. s14->s14FormulaS56];

massEqualityConstraintS56 = Expand[massEqualityConstraint /. s14->s14FormulaS56];
c4Formula9S56 = c4 /. Solve[massEqualityConstraintS56==0, c4][[1]];

mt2ReducedS56  = Together[mt2FormulaS56  /. c4->c4Formula9S56];
s45ReducedS56  = Together[s45FormulaS56  /. c4->c4Formula9S56];
s345ReducedS56 = Together[s345FormulaS56 /. c4->c4Formula9S56];
s16ReducedS56  = Together[s16FormulaS56  /. c4->c4Formula9S56];

(* s15=(p1+p5)^2 is not one of diffInvariants, but derivable the same way as
   s45,s345,s16: s15 = p1^2 + 2 p1.p5 + p5^2 = 2 p1.p5 + mt2, with
   p1.p5 = c2 q12 + c3 q13 + c4 q14 (p1.p1 drops out, q11=0). CAUTION: do NOT
   use the tempting momentum-conservation shortcut "s15 = mt2-s16-s123+s23-s14"
   -- that silently assumes p6^2=mt2, which only holds AFTER c4 is eliminated
   via massEqualityConstraint, not in the raw (pre-elimination) master
   formulas; it gives a wrong answer if applied before eliminating c4. The
   direct dot-product route below is the reliable one (checked against an
   explicit momentum construction, 8/8 independent points). *)
s15FormulaS56 = Expand[2*dot5[1] + mt2Formula /. s14->s14FormulaS56];
s15ReducedS56 = Together[s15FormulaS56 /. c4->c4Formula9S56];

kinematicRulesTH3S56 = {
	q12 -> s12/2,
	q23 -> s23/2,
	q34 -> s34/2,
	q13 -> (s123-s12-s23)/2,
	q24 -> (s234-s23-s34)/2,
	q14 -> Together[s14FormulaS56/2],
	s14  -> s14FormulaS56,
	mt2  -> mt2ReducedS56,
	s16  -> s16ReducedS56,
	s45  -> s45ReducedS56,
	s345 -> s345ReducedS56,
	c4   -> c4Formula9S56
};

ReduceTo9VariablesS56 = {
	s14  -> s14FormulaS56,
	mt2  -> mt2ReducedS56,
	s16  -> s16ReducedS56,
	s45  -> s45ReducedS56,
	s345 -> s345ReducedS56
};


(* ::Section:: *)
(*Random root-free phase-space point generator (finite field, matches your $P)*)


$P = 2^31 - 1;

(* diffInvariants order -- matches ../SetUp/mandelstam_scalarproducts.txt *)
diffInvariants = {mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56};

ModRat[x_, p_:$P] := Mod[Numerator[x]*PowerMod[Denominator[x], -1, p], p];

(* BuildRationalPoint[] draws s12,s23,s34,s123,s234,s14,c1,c2,c3 uniformly at
   random mod $P, solves massEqualityConstraint==0 for c4 (one linear solve,
   no root), and returns an Association over diffInvariants that satisfies
   G6=0 EXACTLY (not approximately) -- so it can feed PhaseSpacePoints/
   ReconstructOnVariety exactly like the existing BCFW-slice samples, but
   covering the full 9-dimensional space instead of a 1D/2D slice through a
   fixed base point. *)
BuildRationalPoint[] := Module[
	{s12v,s23v,s34v,s123v,s234v,s14v,c1v,c2v,c3v,c4v,subsFree,c4sol,subsAll},
	{s12v,s23v,s34v,s123v,s234v,s14v,c1v,c2v,c3v} = RandomInteger[$P-1, 9];
	subsFree = Thread[{s12,s23,s34,s123,s234,s14,c1,c2,c3}
		-> {s12v,s23v,s34v,s123v,s234v,s14v,c1v,c2v,c3v}];
	c4sol = Solve[(massEqualityConstraint /. subsFree) == 0, c4, Modulus->$P][[1]];
	subsAll = Join[subsFree, c4sol];
	<|
		mt2  -> Mod[mt2Formula  /. subsAll, $P],
		s12  -> s12v, s123 -> s123v, s23 -> s23v, s234 -> s234v, s34 -> s34v,
		s16  -> Mod[s16Formula  /. subsAll, $P],
		s345 -> Mod[s345Formula /. subsAll, $P],
		s45  -> Mod[s45Formula  /. subsAll, $P],
		s56  -> Mod[s56Formula  /. subsFree, $P]
	|>
];

(* As a flat point in diffInvariants order, ready for PhaseSpacePoints-style lists: *)
BuildRationalPointVector[] := Module[{pt = BuildRationalPoint[]}, pt /@ diffInvariants];


(* ::Section:: *)
(*Verification*)


(* Run this once to convince yourself: 10 fresh random points, each one
   satisfying the ACTUAL Gram constraint from ../SetUp/SixPointGram.txt
   exactly (not just approximately / mod floating error). *)
VerifyAgainstG6[nPoints_:10] := Module[{dir,G6,pts},
	dir = If[$FrontEnd===Null, Directory[], NotebookDirectory[]];
	G6 = Get[FileNameJoin[{ParentDirectory[dir], "SetUp", "SixPointGram.txt"}]];
	pts = Table[BuildRationalPoint[], {nPoints}];
	Map[Mod[G6 /. Normal[#], $P]&, pts]
];

(*VerifyAgainstG6[10]  (* should print {0,0,0,0,0,0,0,0,0,0} *)*)

(* Symbolic (not just numeric) proof that ReduceTo9Variables kills G6
   identically, as a polynomial in the 9 free variables: *)
VerifyReduceTo9[] := Module[{dir,G6},
	dir = If[$FrontEnd===Null, Directory[], NotebookDirectory[]];
	G6 = Get[FileNameJoin[{ParentDirectory[dir], "SetUp", "SixPointGram.txt"}]];
	Together[G6 /. ReduceTo9Variables] === 0
];

(*VerifyReduceTo9[]  (* should print True *)*)

VerifyReduceTo9KeepMt2[] := Module[{dir,G6},
	dir = If[$FrontEnd===Null, Directory[], NotebookDirectory[]];
	G6 = Get[FileNameJoin[{ParentDirectory[dir], "SetUp", "SixPointGram.txt"}]];
	Together[G6 /. ReduceTo9VariablesKeepMt2] === 0
];

(*VerifyReduceTo9KeepMt2[]  (* should print True *)*)

VerifyReduceTo9KeepMt2AndS12[] := Module[{dir,G6},
	dir = If[$FrontEnd===Null, Directory[], NotebookDirectory[]];
	G6 = Get[FileNameJoin[{ParentDirectory[dir], "SetUp", "SixPointGram.txt"}]];
	Together[G6 /. ReduceTo9VariablesKeepMt2AndS12] === 0
];

(*VerifyReduceTo9KeepMt2AndS12[]  (* should print True *)*)

VerifyReduceTo9S56[] := Module[{dir,G6},
	dir = If[$FrontEnd===Null, Directory[], NotebookDirectory[]];
	G6 = Get[FileNameJoin[{ParentDirectory[dir], "SetUp", "SixPointGram.txt"}]];
	Together[G6 /. ReduceTo9VariablesS56] === 0
];

(*VerifyReduceTo9S56[]  (* should print True *)*)


(* ::Section:: *)
(*Inverse map: recover c1,c2,c3,c4 from an already-consistent named-invariant point*)


(*The forward direction above goes
	(s12,s23,s34,s123,s234,s14,c1,c2,c3) -> the 10 named diffInvariants.
This section is the INVERSE: given a point in the usual redundant
diffInvariants (s16 included, already satisfying G6=0 -- e.g. one produced
the OLD way, by picking 9 invariants at random and solving G6=0 for s16, as
in ../SetUp/SetUpPropagators.wl), solve for the (c1,c2,c3,c4) that produced
it.

This works because a G6=0 point pins p5 down COMPLETELY: mt2,s45,s345,s16
fix all four dot products p5.p1,...,p5.p4, and inverting the (generically
invertible) Gram matrix of the basis p1,p2,p3,p4 then gives a UNIQUE
(c1,c2,c3,c4) -- no leftover freedom, hence no need to separately supply s14
or anything else; it is recovered on the way (s14 = s56-s123+s23-s234).

This is exactly the construction used independently, for the SAME external
kinematics (identical ../SetUp/mandelstam_scalarproducts.txt conventions --
q12,...,q14 below match character-for-character), in
../../FullHexaBoxDE/222-228-Doublebox_7MI/sector239Basis.wl's
kinematicRules239. That file only lists c1,c2,c3: verified symbolically that
CoefficientRulesFromInvariants[[1;;3]] below reproduces its c1,c2,c3 EXACTLY
(as a rational-function identity, not just numerically). It omits c4 because
c4 is redundant there too -- fixed by the same linear p6^2=p5^2 relation as
massEqualityConstraint above. We keep c4 here since there's no reason to
hide it.*)


(* dot products of p1..p4 directly from the named invariants (matches
   ../SetUp/mandelstam_scalarproducts.txt and sector239Basis.wl's q_ij
   exactly): *)
qFromInvariants = {
	q[1,2] -> s12/2, q[2,3] -> s23/2, q[3,4] -> s34/2,
	q[1,3] -> (s123-s12-s23)/2, q[2,4] -> (s234-s23-s34)/2,
	q[1,4] -> (-s123+s23-s234+s56)/2
};
q[i_,i_] := 0;
q[i_,j_] /; i>j := q[j,i];

(* p5.p4,p5.p3,p5.p2,p5.p1 read off mt2,s45,s345,s16,s56 -- each just one more
   "(p_i+p5)^2 = known sum^2 + 2 p_i.p5 + mt2" solved for the dot product: *)
e4FromInvariants = (s45-mt2)/2;
e3FromInvariants = (s345-s34-s45)/2;
e2FromInvariants = (s16-s234-mt2)/2 - e3FromInvariants - e4FromInvariants;
e1FromInvariants = -s56/2 - e2FromInvariants - e3FromInvariants - e4FromInvariants;

gramMatrixP1234 = Table[q[i,j] /. qFromInvariants, {i,4}, {j,4}];
eVectorFromInvariants = {e1FromInvariants,e2FromInvariants,e3FromInvariants,e4FromInvariants};

(* Solved ONCE, symbolically, at load time -- the exact analogue of
   kinematicRules239 (verified to agree with it term-for-term for c1,c2,c3): *)
coefficientsFromInvariantsSolution = Simplify[LinearSolve[gramMatrixP1234, eVectorFromInvariants]];
CoefficientRulesFromInvariants = Thread[{c1,c2,c3,c4} -> coefficientsFromInvariantsSolution];

(* Numeric evaluation at a specific point (an Association or rule list over
   diffInvariants, e.g. straight from BuildRationalPoint[], or from the old
   random-invariants-plus-Solve[G6==0,s16]-for-s16 method), mod $P: *)
CoefficientsAt[pointRules_] := ModRat[#,$P]& /@ (CoefficientRulesFromInvariants[[All,2]] /. Normal[pointRules]);

(* s14 = 2 p1.p4 is not one of diffInvariants but is recoverable from them
   (it is exactly the c-independent piece of s56, see header): *)
S14At[pointRules_] := Mod[(s56-s123+s23-s234) /. Normal[pointRules], $P];

(* Flat, plain-symbol replacement list -- the direct analogue of
   sector239Basis.wl's kinematicRules239 (same q12,...,q14,c1,c2,c3 meaning;
   c4 included too, see header). Use this when you just want "a kinematic
   replacement list" to paste/Get elsewhere, without pulling in the rest of
   this file's machinery: *)
kinematicRulesTH3 = {
	q12 -> q[1,2] /. qFromInvariants,
	q23 -> q[2,3] /. qFromInvariants,
	q34 -> q[3,4] /. qFromInvariants,
	q13 -> q[1,3] /. qFromInvariants,
	q24 -> q[2,4] /. qFromInvariants,
	q14 -> q[1,4] /. qFromInvariants,
	c1  -> CoefficientRulesFromInvariants[[1,2]],
	c2  -> CoefficientRulesFromInvariants[[2,2]],
	c3  -> CoefficientRulesFromInvariants[[3,2]],
	(* c4 from mt2Formula==mt2 directly (p5^2=mt2 is already LINEAR in c4 alone,
	   since p4^2=0 kills the c4^2 term) instead of from the full 4x4 Cramer's-
	   rule solve above -- much shorter, and only needs c1,c2,c3 plus the named
	   invariants (s14 eliminated via s14=s56-s123+s23-s234). Verified to agree
	   with CoefficientRulesFromInvariants[[4,2]] on 16/16 random points: *)
	c4  -> (-mt2 + c1*c2*s12 - c1*c3*s12 + c1*c3*s123 - c1*c3*s23 + c2*c3*s23)/
		(c1*s123 - c1*s23 + c2*s23 + c1*s234 - c2*s234 + c2*s34 - c3*s34 - c1*s56)
};

(* Run once to (re-)export it as a standalone file that Get[...] loads back
   as a plain assignment (same usage pattern as kinematicRules239), without
   needing to load RationalParametrization.wl at all -- NOTE: must export as
   "Text" with an explicit "name = ...;" string, not a bare Export of the
   list (that uses Export's default "Table" format for lists, which Get
   cannot parse back):
Export[FileNameJoin[{DirectoryName[$InputFileName],"kinematicRulesTH3.txt"}],
	"kinematicRulesTH3 = " <> ToString[kinematicRulesTH3, InputForm] <> ";", "Text"]; *)


(* ::Section:: *)
(*Verification of the inverse map*)


(* Round-trip: build a random point the OLD way (pick 9 invariants at random,
   solve G6=0 for s16 -- the method already used throughout this project,
   e.g. ../SetUp/SetUpPropagators.wl), recover (c1,c2,c3,c4) via
   CoefficientsAt, then feed them back through the FORWARD formulas above and
   check every invariant is reproduced exactly.

   Note: over a finite field, G6=0 is a QUADRATIC in s16, so a random draw of
   the other 9 invariants has no solution for s16 whenever the discriminant
   is a non-residue mod $P (roughly half the time) -- this is precisely the
   "implicit square root" problem this whole file exists to avoid. Here we
   just retry until Solve succeeds, purely to get a point to round-trip
   through; BuildRationalPoint[] never has this issue. *)
VerifyInverseMap[nPoints_:10] := Module[{dir,G6,free9vars},
	dir = If[$FrontEnd===Null, Directory[], NotebookDirectory[]];
	G6 = Get[FileNameJoin[{ParentDirectory[dir], "SetUp", "SixPointGram.txt"}]];
	free9vars = {mt2,s12,s123,s23,s234,s34,s345,s45,s56};
	Table[
		Module[{free9,s16sol,pt,cvals,s14v,subsBack},
			s16sol = {};
			While[s16sol==={},
				free9 = Thread[free9vars -> RandomInteger[$P-1,9]];
				s16sol = Solve[Expand[G6/.free9,Modulus->$P]==0, s16, Modulus->$P];
			];
			pt = Join[free9, s16sol[[1]]];
			cvals = CoefficientsAt[pt];
			s14v = S14At[pt];
			subsBack = Join[pt, Thread[{c1,c2,c3,c4}->cvals], {s14->s14v}];
			{Mod[mt2Formula/.subsBack,$P]===Mod[mt2/.pt,$P],
			 Mod[s45Formula/.subsBack,$P]===Mod[s45/.pt,$P],
			 Mod[s345Formula/.subsBack,$P]===Mod[s345/.pt,$P],
			 Mod[s16Formula/.subsBack,$P]===Mod[s16/.pt,$P]}
		]
	, {nPoints}]
];

(*VerifyInverseMap[10]  (* should print 10 copies of {True,True,True,True} *)*)
