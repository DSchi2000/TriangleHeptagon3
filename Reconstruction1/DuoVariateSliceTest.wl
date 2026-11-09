(* ::Package:: *)

Quiet[Get[ParentDirectory[NotebookDirectory[]]<>"/MassiveBCFW/MassiveBCFWSlice.wl"]]
G6=Get[ParentDirectory[NotebookDirectory[]]<>"/SetUp/SixPointGram.txt"];

res=BuildMassiveSlice2D[RandomInteger[$P], False];
avec=Normal[res["a"]][[All;;10]]
bvec=Normal[res["b"]][[All;;10]]
cvec=Normal[res["c"]][[All;;10]]


newpoint=Expand[avec[[All,2]]+RandomInteger[$P]*bvec[[All,2]]+RandomInteger[$P]*cvec[[All,2]],Modulus->$P]
tpoint=Thread[Rule[ToExpression[avec[[All,1]]],newpoint]]
Expand[G6/.tpoint,Modulus->$P]


SliceRepls=Thread[Rule[ToExpression[avec[[All,1]]],Expand[avec[[All,2]]+x*bvec[[All,2]]+y*cvec[[All,2]],Modulus->$P]]]
Expand[G6/.SliceRepls,Modulus->$P]


Facs=ReadList["/home/dschiebe/Documents/PhD-Files/TriangleHeptagon3/Reconstruction1/Factors-3.txt"];


Factor[Facs[[-2]]]
Factor[Facs[[-2]]/.SliceRepls,Modulus->$P]


possible6ptfacs={mt2^2 s23^2+2 mt2 s16 s23 s34-2 mt2 s23^2 s34+s16^2 s34^2-2 s16 s23 s34^2+s23^2 s34^2-2 mt2 s23^2 s345+2 mt2 s23 s234 s345-4 mt2 s23 s34 s345+2 s16 s23 s34 s345-2 s23^2 s34 s345-2 s16 s234 s34 s345+2 s23 s234 s34 s345+s23^2 s345^2-2 s23 s234 s345^2+s234^2 s345^2-2 mt2 s23 s234 s45+2 mt2 s23 s34 s45-4 s16 s23 s34 s45+2 s16 s234 s34 s45+2 s23 s234 s34 s45-2 s16 s34^2 s45-2 s23 s34^2 s45+2 s23 s234 s345 s45-2 s234^2 s345 s45+2 s23 s34 s345 s45+2 s234 s34 s345 s45+s234^2 s45^2-2 s234 s34 s45^2+s34^2 s45^2,-mt2^2 s12 s23+mt2^2 s123 s23-mt2 s12 s16 s34+mt2 s123 s16 s34+mt2^2 s23 s34+mt2 s12 s23 s34-2 mt2 s123 s23 s34+mt2 s16 s34^2-s123 s16 s34^2-mt2 s23 s34^2+s123 s23 s34^2+mt2 s12 s23 s345-2 mt2 s123 s23 s345-mt2 s12 s234 s345+mt2 s123 s234 s345-2 mt2 s123 s34 s345+s123 s16 s34 s345-mt2 s23 s34 s345-2 s123 s23 s34 s345-mt2 s234 s34 s345+s123 s234 s34 s345+s123 s23 s345^2-s123 s234 s345^2+mt2 s12 s23 s45+mt2 s12 s234 s45-mt2 s123 s234 s45+mt2 s12 s34 s45+mt2 s123 s34 s45+s12 s16 s34 s45-2 s123 s16 s34 s45-mt2 s23 s34 s45+s12 s23 s34 s45+mt2 s234 s34 s45-2 s12 s234 s34 s45+s123 s234 s34 s45-mt2 s34^2 s45-s123 s34^2 s45-s16 s34^2 s45-s23 s34^2 s45-s12 s23 s345 s45+s12 s234 s345 s45+s123 s234 s345 s45+s123 s34 s345 s45+s23 s34 s345 s45+s234 s34 s345 s45-s12 s234 s45^2-s12 s34 s45^2-s234 s34 s45^2+s34^2 s45^2+mt2 s23 s345 s56-s16 s34 s345 s56+s23 s34 s345 s56-s23 s345^2 s56+s234 s345^2 s56-mt2 s23 s45 s56+s16 s34 s45 s56+s23 s34 s45 s56+s23 s345 s45 s56-2 s234 s345 s45 s56+s34 s345 s45 s56+s234 s45^2 s56-s34 s45^2 s56,mt2^2 s12^2-2 mt2^2 s12 s123+mt2^2 s123^2-2 mt2^2 s12 s34+2 mt2^2 s123 s34+2 mt2 s12 s123 s34-2 mt2 s123^2 s34+mt2^2 s34^2-2 mt2 s123 s34^2+s123^2 s34^2+2 mt2 s12 s123 s345-2 mt2 s123^2 s345-2 mt2 s123 s34 s345-2 s123^2 s34 s345+s123^2 s345^2-2 mt2 s12^2 s45+2 mt2 s12 s123 s45+4 mt2 s12 s34 s45-2 mt2 s123 s34 s45+2 s12 s123 s34 s45-2 mt2 s34^2 s45-2 s123 s34^2 s45-2 s12 s123 s345 s45+2 s123 s34 s345 s45+s12^2 s45^2-2 s12 s34 s45^2+s34^2 s45^2-2 mt2 s12 s345 s56+2 mt2 s123 s345 s56-2 mt2 s34 s345 s56+2 s123 s34 s345 s56-2 s123 s345^2 s56+2 mt2 s12 s45 s56-2 mt2 s123 s45 s56+2 mt2 s34 s45 s56-4 s12 s34 s45 s56+2 s123 s34 s45 s56+2 s12 s345 s45 s56+2 s123 s345 s45 s56+2 s34 s345 s45 s56-2 s12 s45^2 s56-2 s34 s45^2 s56+s345^2 s56^2-2 s345 s45 s56^2+s45^2 s56^2,s23,s34,-mt2^2 s234+mt2^2 s34-mt2 s16 s34+mt2 s23 s34+mt2 s234 s345+mt2 s234 s45-mt2 s34 s45+s16 s34 s45-s234 s345 s45,mt2 s123^2 s34+mt2^2 s12 s56-mt2^2 s123 s56+mt2 s123 s345 s56-2 mt2 s12 s45 s56+mt2 s123 s45 s56-s123 s345 s45 s56+s12 s45^2 s56,-mt2^2 s12^2 s234+2 mt2^2 s12 s123 s234-mt2^2 s123^2 s234-mt2^2 s12 s123 s34+mt2^2 s123^2 s34+mt2 s12 s123 s16 s34-mt2 s123^2 s16 s34-2 mt2 s12 s123 s23 s34+mt2^2 s12 s234 s34-mt2^2 s123 s234 s34-mt2 s12 s123 s234 s34+mt2 s123^2 s234 s34+mt2^2 s123 s34^2-mt2 s123^2 s34^2-mt2 s123 s16 s34^2-s123^2 s16 s34^2+2 mt2 s123 s23 s34^2-2 mt2 s12 s123 s234 s345+2 mt2 s123^2 s234 s345-mt2 s123^2 s34 s345+s123^2 s16 s34 s345+mt2 s123 s234 s34 s345+s123^2 s234 s34 s345-s123^2 s234 s345^2+2 mt2 s12^2 s234 s45-2 mt2 s12 s123 s234 s45+mt2 s12 s123 s34 s45-s12 s123 s16 s34 s45-2 mt2 s12 s234 s34 s45+mt2 s123 s234 s34 s45-s12 s123 s234 s34 s45-mt2 s123 s34^2 s45+s123 s16 s34^2 s45+2 s12 s123 s234 s345 s45-s123 s234 s34 s345 s45-s12^2 s234 s45^2+s12 s234 s34 s45^2-mt2^2 s12 s23 s56+mt2^2 s123 s23 s56-2 mt2 s12 s16 s34 s56+2 mt2 s123 s16 s34 s56-mt2^2 s23 s34 s56+2 mt2 s12 s23 s34 s56-mt2 s123 s23 s34 s56+mt2 s12 s23 s345 s56-2 mt2 s123 s23 s345 s56+mt2 s12 s234 s345 s56-mt2 s123 s234 s345 s56-mt2 s123 s34 s345 s56-s123 s16 s34 s345 s56+mt2 s23 s34 s345 s56-s123 s23 s34 s345 s56+s123 s23 s345^2 s56+s123 s234 s345^2 s56+mt2 s12 s23 s45 s56-mt2 s12 s234 s45 s56+mt2 s123 s234 s45 s56+2 mt2 s12 s34 s45 s56-mt2 s123 s34 s45 s56+2 s12 s16 s34 s45 s56-s123 s16 s34 s45 s56+mt2 s23 s34 s45 s56+2 s12 s23 s34 s45 s56-s12 s23 s345 s45 s56-s12 s234 s345 s45 s56-s123 s234 s345 s45 s56+2 s123 s34 s345 s45 s56-s23 s34 s345 s45 s56+s12 s234 s45^2 s56-2 s12 s34 s45^2 s56+mt2 s23 s345 s56^2-s23 s345^2 s56^2-mt2 s23 s45 s56^2+s23 s345 s45 s56^2};


SlicedFactors=Factor[possible6ptfacs/.SliceRepls,Modulus->$P]
SortBy[Map[FactorList[#][[All,1]]&,SlicedFactors]//Flatten//Tally,Last]


P=8*mt2^2*s123*s23*s56 + 125*mt2^2*s345^2*s45 - 96*mt2*s12*s123^2*s45 - 76*mt2*s12*s123*s34*s45 - 36*mt2*s12*s23*s234*s345 + 90*mt2*s123*s16^2*s34 + 106*mt2*s16*s23*s45^2 + 120*mt2*s16*s234*s345*s56 - mt2*s23*s234*s34*s56 - 53*mt2*s234*s34^2*s345 + 108*mt2*s345*s45*s56^2 - 117*s12*s123*s234*s34^2 + 2*s12*s16*s23*s234*s345 - 96*s12*s16*s23*s34^2 + 54*s12*s23*s34*s345*s56 + 36*s12*s234^2*s34*s345 - 69*s123*s234^2*s56^2 + 46*s23*s34*s345*s45^2


P//Variables


Factor[P/.SliceTpoint,Modulus->$P]
Factor[P/.SliceRepls,Modulus->$P]



