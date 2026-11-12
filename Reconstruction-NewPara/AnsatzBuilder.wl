(* ::Package:: *)

(* ===================================================================== *)
(* AnsatzBuilder.wl                                                       *)
(*                                                                        *)
(* Builds a rational-function ansatz Num/Den where:                      *)
(*   - the "s"-type variables (mass dimension 1) appear homogeneously,   *)
(*     i.e. every monomial in a given polynomial has the SAME total      *)
(*     degree in the s-variables. The exponent vectors are generated     *)
(*     with FrobeniusSolve[{1,1,...,1}, degree] -- the standard way to   *)
(*     enumerate all non-negative integer solutions of                   *)
(*         e1 + e2 + ... + en == degree                                  *)
(*     which is exactly the set of degree-"degree" monomial exponents   *)
(*     in n variables.                                                   *)
(*   - the "c"-type variables (dimensionless) appear INDEPENDENTLY, each *)
(*     capped at a given max power (e.g. "at most squared" -> 0,1,2),    *)
(*     with no constraint tying their powers together (non-homogeneous). *)
(*                                                                        *)
(* ===================================================================== *)

ClearAll[SMonomialExponents, CMonomialExponents, AnsatzMonomials, BuildAnsatzPoly, BuildRationalAnsatz];

(* Exponent vectors for the homogeneous (s-type) part, degree `deg`
   in the variables `svars` -- same format FrobeniusSolve returns:
   a list of non-negative integer vectors summing to deg. *)
SMonomialExponents[svars_List, deg_Integer] :=
  FrobeniusSolve[ConstantArray[1, Length[svars]], deg]

(* Exponent vectors for the non-homogeneous (c-type) part: each
   variable's power ranges independently over 0..maxPow. *)
CMonomialExponents[cvars_List, maxPow_Integer] :=
  Tuples[Range[0, maxPow], Length[cvars]]

(* Turn a list of exponent vectors into actual monomials in `vars`. *)
MonomialsFromExponents[vars_List, exps_List] :=
  Times @@@ (vars^# & /@ exps)

(* Full monomial list for one polynomial (numerator OR denominator):
   every s-monomial of exact degree sDeg, times every c-monomial with
   each power in 0..cMaxPow. *)
AnsatzMonomials[svars_List, cvars_List, sDeg_Integer, cMaxPow_Integer] :=
  Module[{sMons, cMons},
    sMons = MonomialsFromExponents[svars, SMonomialExponents[svars, sDeg]];
    cMons = MonomialsFromExponents[cvars, CMonomialExponents[cvars, cMaxPow]];
    Flatten[Outer[Times, sMons, cMons]]
  ]

(* Build Sum[coeff[i] * monomial_i] with a fresh coefficient head. *)
BuildAnsatzPoly[svars_List, cvars_List, deg_Integer, cMaxPow_Integer, coeffHead_Symbol] :=
  Module[{mons},
    mons = AnsatzMonomials[svars, cvars, deg, cMaxPow];
    coeffHead /@ Range[Length[mons]] . mons
  ]

(* Build the full rational ansatz Num/Den.
   Options:
     "NumDegree"  -> total s-degree of the numerator
     "DenDegree"  -> total s-degree of the denominator
     "CMaxPower"  -> max power of each c-variable (independent)
     "Normalize"  -> True/False, fix the first denominator coefficient to 1
                     to remove the overall rescaling redundancy between
                     Num and Den (default: True)
     "NumCoeff"/"DenCoeff" -> symbols used for numerator/denominator
                     coefficients (default: a, b)
*)
BuildRationalAnsatz[svars_List, cvars_List, OptionsPattern[{
    "NumDegree" -> 2, "DenDegree" -> 3, "CMaxPower" -> 2,
    "Normalize" -> True, "NumCoeff" -> a, "DenCoeff" -> b}]] :=
  Module[{numMons, denMons, numPoly, denPoly, nNum, nDen, a1, b1},
    a1 = OptionValue["NumCoeff"]; b1 = OptionValue["DenCoeff"];
    numMons = AnsatzMonomials[svars, cvars, OptionValue["NumDegree"], OptionValue["CMaxPower"]];
    denMons = AnsatzMonomials[svars, cvars, OptionValue["DenDegree"], OptionValue["CMaxPower"]];
    nNum = Length[numMons]; nDen = Length[denMons];
    numPoly = (a1 /@ Range[nNum]) . numMons;
    denPoly = If[TrueQ[OptionValue["Normalize"]],
      1 + (b1 /@ Range[2, nDen]) . denMons[[2 ;;]],
      (b1 /@ Range[nDen]) . denMons
    ];
    <|
      "Numerator" -> numPoly,
      "Denominator" -> denPoly,
      "Ansatz" -> numPoly/denPoly,
      "NumeratorMonomials" -> numMons,
      "DenominatorMonomials" -> denMons,
      "NumTerms" -> nNum,
      "DenTerms" -> If[TrueQ[OptionValue["Normalize"]], nDen - 1, nDen],
      "TotalFreeCoefficients" -> nNum + If[TrueQ[OptionValue["Normalize"]], nDen - 1, nDen]
    |>
  ]


(* ===================================================================== *)
(* Example usage, matching the case we worked out:                       *)
(*   s-vars: s12,s123,s23,s234,s34,s56  (6 vars, mass dim 1)              *)
(*   c-vars: c1,c2,c3                   (3 vars, dimensionless)           *)
(*   numerator homogeneous degree 3 in s, denominator degree 4 in s       *)
(*   each c_i up to power 2, independently, in both Num and Den           *)
(* ===================================================================== *)

(*
svars = {s12, s123, s23, s234, s34, s56};
cvars = {c1, c2, c3};

result = BuildRationalAnsatz[svars, cvars,
  "NumDegree" -> 3, "DenDegree" -> 4, "CMaxPower" -> 2, "Normalize" -> True];

result["NumTerms"]                 (* 1512 *)
result["DenTerms"]                 (* 3401, i.e. 3402 - 1 after normalization *)
result["TotalFreeCoefficients"]    (* 4913 *)
result["Ansatz"]                   (* the full Num/Den expression with a[i], b[i] unknowns *)
*)
