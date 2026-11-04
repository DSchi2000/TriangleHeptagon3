(* ::Package:: *)

(* ::Title:: *)
(*Dual-number (forward-mode) differentiation mod a prime*)

(* Drop-in replacement for the pattern
     Expand[(D[expr,var]/.pt)/.root'[x_]:>1/(2 root[x])//.root[x_]:>PowerMod[x,1/2,$P],Modulus->$P]
   when `expr` is a function of a single symbol `var` (built from Integer/Rational
   constants, Plus, Times, Power[_,_Integer], and the formal square-root head root[_]).
   DualEval[expr,var,x0,p] returns {value,derivative} of expr at var->x0, both mod p,
   computed in one pass with no intermediate symbolic expression ever built -- so there
   is no Together/quotient-rule blow-up regardless of how `expr` is structured. *)


dAdd[{a1_,d1_},{a2_,d2_},p_]:={Mod[a1+a2,p],Mod[d1+d2,p]};
dMul[{a1_,d1_},{a2_,d2_},p_]:={Mod[a1*a2,p],Mod[a1*d2+a2*d1,p]};
dPow[{a_,d_},n_Integer,p_]:={PowerMod[a,n,p],Mod[n*PowerMod[a,n-1,p]*d,p]};
dSqrt[{a_,d_},p_]:=Module[{s},s=PowerMod[a,1/2,p];{s,Mod[d*PowerMod[2 s,-1,p],p]}];

DualEval::badsym="Unexpected symbol `1` encountered (expected only the differentiation variable).";
DualEval::badhead="Unexpected head `1` in subexpression `2` -- DualEval only supports Plus, Times, Power[_,_Integer], root[_], and numeric literals.";

DualEval[expr_,var_Symbol,x0_,p_Integer]:=Module[{recEval},
	recEval[e_Integer]:={Mod[e,p],0};
	recEval[e_Rational]:={Mod[Numerator[e]*PowerMod[Denominator[e],-1,p],p],0};
	recEval[e_Symbol]:=If[e===var,{Mod[x0,p],1},
		(Message[DualEval::badsym,e];Throw[$Failed])];
	recEval[Plus[a_,b__]]:=Fold[dAdd[#1,#2,p]&,recEval[a],recEval/@{b}];
	recEval[Times[a_,b__]]:=Fold[dMul[#1,#2,p]&,recEval[a],recEval/@{b}];
	recEval[Power[b_,n_Integer]]:=dPow[recEval[b],n,p];
	recEval[Power[b_,Rational[1,2]]]:=dSqrt[recEval[b],p];
	recEval[root[a_]]:=dSqrt[recEval[a],p];
	recEval[e_]:=(Message[DualEval::badhead,Head[e],e];Throw[$Failed]);
	Catch[recEval[expr]]
];
