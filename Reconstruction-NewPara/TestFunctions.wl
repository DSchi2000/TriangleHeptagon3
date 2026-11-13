(* ::Package:: *)

Get["/home/dschiebe/Documents/PhD-Files/TriangleHeptagon3/files/ExtraRecoModules.m"]
$P=2^31-1


Func=(c1^2 s12 s23+c1 s23 s56)/((c1-c3) s12^2+s23 s12)


vars=Func//Variables
tpoints=Table[Thread[Rule[vars,RandomInteger[$P,Length[vars]]]],{i,1,1000}];
fitvec=Expand[Func/.tpoints,Modulus->$P];


ClearAll[ReconstructInhomogineus];
Get["/home/dschiebe/Documents/PhD-Files/TriangleHeptagon3/files/ExtraRecoModules.m"]
res=ReconstructInhomogineus[fitvec,tpoints,{{c1,c3},2,1},{{s12,s23,s56},2,2},Debug->True]


Get["/home/dschiebe/Documents/PhD-Files/TriangleHeptagon3/files/ExtraRecoModules.m"]
FuncR=BuildFunction[res[[1,1]],res[[2]]]


Factor[FuncR,Modulus->$P]
Factor[Func,Modulus->$P]


tpoints[[-1]]
tpoints[[-1,All,2]]


TestVec1={{0,0,2},{0,1,1},{0,2,0},{1,0,1},{1,1,0},{2,0,0}}
TestVec2=Map[Exponent[#,{c1,c2}]&,MonomialList[Expand[Times@@Table[Sum[c^n,{n,0,2}],{c,{c1,c2}}]]]]
Flatten[Outer[Join, TestVec1, TestVec2, 1, 1],{1,2}]


tpoints[[All,All,2]]//Variables
tpoints[[1,All,1]]



