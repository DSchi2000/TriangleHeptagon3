(* ::Package:: *)

If[$FrontEnd===Null,$DIR=Directory[],$DIR=NotebookDirectory[]];
parts = FileNameSplit[$DIR];
pos = First@FirstPosition[parts, "TriangleHeptagon3"];
HomeDir = FileNameJoin[parts[[;;pos ]]];
SectorDir = FileNameJoin[parts[[;;pos+1]]];
CurrentDir=$DIR;
kirapath=SectorDir;
$P=2^31-1;
diffInvariants={mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56};
(*PureIntegrals=ReadList[FileNameJoin[{SectorDir,"SubBasis.txt"}]];*)
SetDirectory[CurrentDir];
NotbookDir=NotebookDirectory[];
Get[FileNameJoin[{HomeDir,"files/ExtraFunctions.m"}]];
Get[FileNameJoin[{HomeDir,"files/HexaBoxModules.m"}]];

GetLinPos[x1_,x2_,NMI_]:=First[Flatten[Position[Flatten[Table[{i,j},{i,1,NMI},{j,1,NMI}],1],{x1,x2}]]];
Get[ToString[StringJoin[HomeDir,"/files/ExtraModules.m"]]];


Get[NotebookDirectory[]<>"kinematicRulesTH3.txt"]


kinematicRulesTH3


p5repls={p[5]->c1 p[1]+c2 p[2]+c3 p[3]+c4 p[4]}


fiveVectorList={p[1]+p[6],p[3]+p[4]+p[5],p[4]+p[5],p[5]}





NewParaVrars=Variables[NewRepls[[All,2]]][[;;-2]]


Complement[diffInvariants,Test//Variables]


NewParaVrars


sol=Solve[Factor[G6/.NewRepls]==0,s56]


sol[[1]]


tpointNewPara=Thread[Rule[NewParaVrars,RandomInteger[$P,9]]];
sol=Solve[Factor[G6/.NewRepls/.tpointNewPara]==0,s56,Modulus->$P]
TpointOG=SortBy[Join[Join[Expand[NewRepls/.tpointNewPara/.sol[[1]],Modulus->$P],sol[[1]]],tpointNewPara[[4;;]]],First]
Expand[Expand[kinematicRulesTH3/.TpointOG,Modulus->$P][[7;;9,2]]-tpointNewPara[[;;3,2]],Modulus->$P]


TpointOG=SortBy[Join[Join[Expand[NewRepls/.tpointNewPara/.sol[[2]],Modulus->$P],sol[[2]]],tpointNewPara[[4;;]]],First]
Expand[Expand[kinematicRulesTH3/.TpointOG,Modulus->$P][[7;;9,2]]-tpointNewPara[[;;3,2]],Modulus->$P]


Expand[G6/.TpointOG,Modulus->$P]





Get[NotebookDirectory[]<>"kinematicRulesTH3KeepMt2.txt"]


tpointNewPara2=kinematicRulesTH3KeepMt2[[All,2]]//Variables
NewParaVrars





kinematicRulesTH3KeepMt2[[7]]
kinematicRulesTH3KeepMt2[[8]]
kinematicRulesTH3KeepMt2[[9]]
kinematicRulesTH3KeepMt2[[10]]
kinematicRulesTH3KeepMt2[[11]]


extratestrepl=Solve[(Expand[(p[1]+p[4])^2]/.kinematics/.mandelstamScalarProducts)==s14,s56]


tpointNewPara=Thread[Rule[tpointNewPara2,RandomInteger[$P,9]]]


oldpoint=SortBy[Join[tpointNewPara[[4;;]],Expand[kinematicRulesTH3KeepMt2[[7;;11]]/.tpointNewPara,Modulus->$P]],First]



Expand[G6/.oldpoint,Modulus->$P]


Expand[kinematicRulesTH3[[-4;;-2]]/.oldpoint,Modulus->$P]


Test=Collect[Expand[(p[5]/.p4cons/.p5repls)^2]/.kinematics/.mandelstamScalarProducts,diffInvariants]


Test/.extratestrepl//Expand


(* ::Section:: *)
(*{s12, s23, s34, s123, s234, s56, c1, c2, c3}*)


Get[NotebookDirectory[]<>"kinematicRulesTH3S56.txt"]
cvals=Table[M1=M;M1[[All,i]]=b;Together[Det[M1]/Det[M]],{i,1,4}]


tpointNewPara2=kinematicRulesTH3S56[[All,2]]//Variables
tpointNewPara=Thread[Rule[tpointNewPara2,RandomInteger[$P,9]]]
oldpoint=SortBy[Join[tpointNewPara[[4;;]],Expand[kinematicRulesTH3S56[[8;;11]]/.tpointNewPara,Modulus->$P]],First]
Expand[G6/.oldpoint,Modulus->$P]



Expand[kinematicRulesTH3/.oldpoint,Modulus->$P]


crepls=kinematicRulesTH3[[-4;;-2]];
p5repls
p4cons/.p5repls


EQA=Expand[(p[5]/.p5repls)^2]/.kinematics/.mandelstamScalarProducts
EQB=Expand[(p[6]/.p4cons/.p5repls)^2]/.kinematics/.mandelstamScalarProducts


Coefficient[EQA,c4]//EchoFunction[Variables]
Coefficient[EQB,c4]//EchoFunction[Variables]


Coefficient[EQA,c3]//EchoFunction[Variables]
Coefficient[EQB,c3]//EchoFunction[Variables]


solc4=First[Solve[EQA==EQB,c3]]
Factor[EQA/.solc4]


solc4=First[Solve[EQA==EQB,c3]]
Factor[EQA/.solc4]


s123
s12
s23

s234
s34
s56


momenta={p[1],p[2],p[3],p[4]}
M=Table[vec1 vec2,{vec1,momenta},{vec2,momenta}]/.kinematics/.mandelstamScalarProducts
b=Table[p[5] vec1, {vec1,momenta}]/.kinematics/.mandelstamScalarProducts


Together[kinematicRulesTH3[[-1]]]


cvals[[4]]


(* ::Section:: *)
(*Build Parametrization:*)


variables=diffInvariants
variablesNew={s12, s23, s34, s123, s234, s56, c1, c2, c3}
Intersection[variablesNew,variables]

cvals=Thread[Rule[{c1,c2,c3,c4},Table[M1=M;M1[[All,i]]=b;Together[Det[M1]/Det[M]],{i,1,4}]]];

p5repls={p[5]->c1 p[1]+c2 p[2]+c3 p[3]+c4 p[4]}
EQA=Expand[(p[5]/.p5repls)^2]/.kinematics/.mandelstamScalarProducts
EQB=Expand[(p[6]/.p4cons/.p5repls)^2]/.kinematics/.mandelstamScalarProducts

solc4=First[Solve[EQA==EQB,c4]]
sol10thmand=First[Solve[Factor[EQA/.solc4]==mt2,mt2]]


fiveVectorList={p[1]+p[6],p[3]+p[4]+p[5],p[4]+p[5]}
newParaMands=Expand[(fiveVectorList/.p4cons/.p5repls)^2]/.kinematics/.mandelstamScalarProducts
c4sol=Solve[Test[[4]]==mt2,c4][[1]]
solmands=Thread[Rule[{s16,s345,s45},Together[newParaMands/.solc4]]]


KinematicRulesForward=Join[Thread[Rule[{s12, s23, s34, s123, s234, s56},{s12, s23, s34, s123, s234, s56}]],cvals[[;;-2]]]

KinematicRulesBackWards=SortBy[Join[Thread[Rule[{s12, s23, s34, s123, s234, s56},{s12, s23, s34, s123, s234, s56}]],sol10thmand,solmands],First]


(* ::Subsection:: *)
(*Test Parametrization*)


tpointNewPara=Thread[Rule[variablesNew,RandomInteger[$P,9]]]

mt2repls=Expand[KinematicRulesBackWards[[1]]/.tpointNewPara,Modulus->$P]
Tpoint=Thread[Rule[diffInvariants,Expand[diffInvariants/.KinematicRulesBackWards/.tpointNewPara/.mt2repls,Modulus->$P]]]

Expand[G6/.Tpoint,Modulus->$P]



(* ::Section:: *)
(*Derivative*)


func=s16^2+s16 s12
EQN=Expand[-1 Sum[D[G6,var]ds[ToExpression[Capitalize[ToString[var]]]],{var,diffInvariants}]/.Tpoint,Modulus->$P]
sixpfacs=Coefficient[Solve[EQN==0,ds[S16],Modulus->$P][[1,1,2]],Table[ds[ToExpression[Capitalize[ToString[var]]]],{var,diffInvariants}]]

DiffNumbers=Expand[Table[D[func,diffInvariants[[i]]]+sixpfacs[[i]]D[func,s16],{i,1,Length[sixpfacs]}]/.Tpoint,Modulus->$P];
DiffNumbers[[4]]=0;
DiffNumbers


KinematicRulesBackWards[[All,1]]


J=Expand[Table[D[KinematicRulesBackWards[[i,2]],var],{i,{1,2,3,5,6,7,8,9,10}},{var,variablesNew}]/.tpointNewPara,Modulus->$P]
f=DiffNumbers[[Flatten[Map[Position[diffInvariants,#]&,{mt2,s12,s123,s23,s234,s34,s345,s45,s56}]]]]
NewDiff=Expand[Transpose[J] . f,Modulus->$P]


func2=Together[func/.KinematicRulesBackWards];
DiffNumbers2=Expand[Table[D[func2,variablesNew[[i]]],{i,1,Length[variablesNew]}]/.tpointNewPara,Modulus->$P]


(* ::Section:: *)
(*Uni Variate slice:*)


tpointNewPara1=Thread[Rule[variablesNew,RandomInteger[$P,9]]]
tpointNewPara2=Thread[Rule[variablesNew,RandomInteger[$P,9]]]

SliceReplsNew=Thread[Rule[variablesNew,tpointNewPara1[[All,2]]+t tpointNewPara2[[All,2]]]]

ABVEC=ReadList[HomeDir<>"/Reconstruction1/Samples/AVECS.txt"];
avec=ABVEC[[1]];
bvec=ABVEC[[2]];
SliceTpoint=Thread[Rule[bvec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]];


possible6ptfacs={mt2^2 s23^2+2 mt2 s16 s23 s34-2 mt2 s23^2 s34+s16^2 s34^2-2 s16 s23 s34^2+s23^2 s34^2-2 mt2 s23^2 s345+2 mt2 s23 s234 s345-4 mt2 s23 s34 s345+2 s16 s23 s34 s345-2 s23^2 s34 s345-2 s16 s234 s34 s345+2 s23 s234 s34 s345+s23^2 s345^2-2 s23 s234 s345^2+s234^2 s345^2-2 mt2 s23 s234 s45+2 mt2 s23 s34 s45-4 s16 s23 s34 s45+2 s16 s234 s34 s45+2 s23 s234 s34 s45-2 s16 s34^2 s45-2 s23 s34^2 s45+2 s23 s234 s345 s45-2 s234^2 s345 s45+2 s23 s34 s345 s45+2 s234 s34 s345 s45+s234^2 s45^2-2 s234 s34 s45^2+s34^2 s45^2,-mt2^2 s12 s23+mt2^2 s123 s23-mt2 s12 s16 s34+mt2 s123 s16 s34+mt2^2 s23 s34+mt2 s12 s23 s34-2 mt2 s123 s23 s34+mt2 s16 s34^2-s123 s16 s34^2-mt2 s23 s34^2+s123 s23 s34^2+mt2 s12 s23 s345-2 mt2 s123 s23 s345-mt2 s12 s234 s345+mt2 s123 s234 s345-2 mt2 s123 s34 s345+s123 s16 s34 s345-mt2 s23 s34 s345-2 s123 s23 s34 s345-mt2 s234 s34 s345+s123 s234 s34 s345+s123 s23 s345^2-s123 s234 s345^2+mt2 s12 s23 s45+mt2 s12 s234 s45-mt2 s123 s234 s45+mt2 s12 s34 s45+mt2 s123 s34 s45+s12 s16 s34 s45-2 s123 s16 s34 s45-mt2 s23 s34 s45+s12 s23 s34 s45+mt2 s234 s34 s45-2 s12 s234 s34 s45+s123 s234 s34 s45-mt2 s34^2 s45-s123 s34^2 s45-s16 s34^2 s45-s23 s34^2 s45-s12 s23 s345 s45+s12 s234 s345 s45+s123 s234 s345 s45+s123 s34 s345 s45+s23 s34 s345 s45+s234 s34 s345 s45-s12 s234 s45^2-s12 s34 s45^2-s234 s34 s45^2+s34^2 s45^2+mt2 s23 s345 s56-s16 s34 s345 s56+s23 s34 s345 s56-s23 s345^2 s56+s234 s345^2 s56-mt2 s23 s45 s56+s16 s34 s45 s56+s23 s34 s45 s56+s23 s345 s45 s56-2 s234 s345 s45 s56+s34 s345 s45 s56+s234 s45^2 s56-s34 s45^2 s56,mt2^2 s12^2-2 mt2^2 s12 s123+mt2^2 s123^2-2 mt2^2 s12 s34+2 mt2^2 s123 s34+2 mt2 s12 s123 s34-2 mt2 s123^2 s34+mt2^2 s34^2-2 mt2 s123 s34^2+s123^2 s34^2+2 mt2 s12 s123 s345-2 mt2 s123^2 s345-2 mt2 s123 s34 s345-2 s123^2 s34 s345+s123^2 s345^2-2 mt2 s12^2 s45+2 mt2 s12 s123 s45+4 mt2 s12 s34 s45-2 mt2 s123 s34 s45+2 s12 s123 s34 s45-2 mt2 s34^2 s45-2 s123 s34^2 s45-2 s12 s123 s345 s45+2 s123 s34 s345 s45+s12^2 s45^2-2 s12 s34 s45^2+s34^2 s45^2-2 mt2 s12 s345 s56+2 mt2 s123 s345 s56-2 mt2 s34 s345 s56+2 s123 s34 s345 s56-2 s123 s345^2 s56+2 mt2 s12 s45 s56-2 mt2 s123 s45 s56+2 mt2 s34 s45 s56-4 s12 s34 s45 s56+2 s123 s34 s45 s56+2 s12 s345 s45 s56+2 s123 s345 s45 s56+2 s34 s345 s45 s56-2 s12 s45^2 s56-2 s34 s45^2 s56+s345^2 s56^2-2 s345 s45 s56^2+s45^2 s56^2,s23,s34,-mt2^2 s234+mt2^2 s34-mt2 s16 s34+mt2 s23 s34+mt2 s234 s345+mt2 s234 s45-mt2 s34 s45+s16 s34 s45-s234 s345 s45,mt2 s123^2 s34+mt2^2 s12 s56-mt2^2 s123 s56+mt2 s123 s345 s56-2 mt2 s12 s45 s56+mt2 s123 s45 s56-s123 s345 s45 s56+s12 s45^2 s56,-mt2^2 s12^2 s234+2 mt2^2 s12 s123 s234-mt2^2 s123^2 s234-mt2^2 s12 s123 s34+mt2^2 s123^2 s34+mt2 s12 s123 s16 s34-mt2 s123^2 s16 s34-2 mt2 s12 s123 s23 s34+mt2^2 s12 s234 s34-mt2^2 s123 s234 s34-mt2 s12 s123 s234 s34+mt2 s123^2 s234 s34+mt2^2 s123 s34^2-mt2 s123^2 s34^2-mt2 s123 s16 s34^2-s123^2 s16 s34^2+2 mt2 s123 s23 s34^2-2 mt2 s12 s123 s234 s345+2 mt2 s123^2 s234 s345-mt2 s123^2 s34 s345+s123^2 s16 s34 s345+mt2 s123 s234 s34 s345+s123^2 s234 s34 s345-s123^2 s234 s345^2+2 mt2 s12^2 s234 s45-2 mt2 s12 s123 s234 s45+mt2 s12 s123 s34 s45-s12 s123 s16 s34 s45-2 mt2 s12 s234 s34 s45+mt2 s123 s234 s34 s45-s12 s123 s234 s34 s45-mt2 s123 s34^2 s45+s123 s16 s34^2 s45+2 s12 s123 s234 s345 s45-s123 s234 s34 s345 s45-s12^2 s234 s45^2+s12 s234 s34 s45^2-mt2^2 s12 s23 s56+mt2^2 s123 s23 s56-2 mt2 s12 s16 s34 s56+2 mt2 s123 s16 s34 s56-mt2^2 s23 s34 s56+2 mt2 s12 s23 s34 s56-mt2 s123 s23 s34 s56+mt2 s12 s23 s345 s56-2 mt2 s123 s23 s345 s56+mt2 s12 s234 s345 s56-mt2 s123 s234 s345 s56-mt2 s123 s34 s345 s56-s123 s16 s34 s345 s56+mt2 s23 s34 s345 s56-s123 s23 s34 s345 s56+s123 s23 s345^2 s56+s123 s234 s345^2 s56+mt2 s12 s23 s45 s56-mt2 s12 s234 s45 s56+mt2 s123 s234 s45 s56+2 mt2 s12 s34 s45 s56-mt2 s123 s34 s45 s56+2 s12 s16 s34 s45 s56-s123 s16 s34 s45 s56+mt2 s23 s34 s45 s56+2 s12 s23 s34 s45 s56-s12 s23 s345 s45 s56-s12 s234 s345 s45 s56-s123 s234 s345 s45 s56+2 s123 s34 s345 s45 s56-s23 s34 s345 s45 s56+s12 s234 s45^2 s56-2 s12 s34 s45^2 s56+mt2 s23 s345 s56^2-s23 s345^2 s56^2-mt2 s23 s45 s56^2+s23 s345 s45 s56^2};


DLog=Log[(-mt2^2 s12^2 s234+2 mt2^2 s12 s123 s234-mt2^2 s123^2 s234-mt2^2 s12 s123 s34+mt2^2 s123^2 s34+mt2^2 s12 s234 s34-mt2^2 s123 s234 s34-mt2 s12 s123 s234 s34+mt2 s123^2 s234 s34+mt2^2 s123 s34^2-mt2 s123^2 s34^2-2 mt2 s12 s123 s234 s345+2 mt2 s123^2 s234 s345-mt2 s123^2 s34 s345+mt2 s123 s234 s34 s345+s123^2 s234 s34 s345-s123^2 s234 s345^2+2 mt2 s12^2 s234 s45-2 mt2 s12 s123 s234 s45+mt2 s12 s123 s34 s45-2 mt2 s12 s234 s34 s45+mt2 s123 s234 s34 s45-s12 s123 s234 s34 s45-mt2 s123 s34^2 s45+2 s12 s123 s234 s345 s45-s123 s234 s34 s345 s45-s12^2 s234 s45^2+s12 s234 s34 s45^2+mt2 s12 s234 s345 s56-mt2 s123 s234 s345 s56-mt2 s123 s34 s345 s56+s123 s234 s345^2 s56-mt2 s12 s234 s45 s56+mt2 s123 s234 s45 s56+2 mt2 s12 s34 s45 s56-mt2 s123 s34 s45 s56-s12 s234 s345 s45 s56-s123 s234 s345 s45 s56+2 s123 s34 s345 s45 s56+s12 s234 s45^2 s56-2 s12 s34 s45^2 s56)/(-mt2^2 s12^2 s234+2 mt2^2 s12 s123 s234-mt2^2 s123^2 s234-mt2^2 s12 s123 s34+mt2^2 s123^2 s34+mt2 s12 s123 s16 s34-mt2 s123^2 s16 s34-2 mt2 s12 s123 s23 s34+mt2^2 s12 s234 s34-mt2^2 s123 s234 s34-mt2 s12 s123 s234 s34+mt2 s123^2 s234 s34+mt2^2 s123 s34^2-mt2 s123^2 s34^2-mt2 s123 s16 s34^2-s123^2 s16 s34^2+2 mt2 s123 s23 s34^2-2 mt2 s12 s123 s234 s345+2 mt2 s123^2 s234 s345-mt2 s123^2 s34 s345+s123^2 s16 s34 s345+mt2 s123 s234 s34 s345+s123^2 s234 s34 s345-s123^2 s234 s345^2+2 mt2 s12^2 s234 s45-2 mt2 s12 s123 s234 s45+mt2 s12 s123 s34 s45-s12 s123 s16 s34 s45-2 mt2 s12 s234 s34 s45+mt2 s123 s234 s34 s45-s12 s123 s234 s34 s45-mt2 s123 s34^2 s45+s123 s16 s34^2 s45+2 s12 s123 s234 s345 s45-s123 s234 s34 s345 s45-s12^2 s234 s45^2+s12 s234 s34 s45^2-mt2^2 s12 s23 s56+mt2^2 s123 s23 s56-2 mt2 s12 s16 s34 s56+2 mt2 s123 s16 s34 s56-mt2^2 s23 s34 s56+2 mt2 s12 s23 s34 s56-mt2 s123 s23 s34 s56+mt2 s12 s23 s345 s56-2 mt2 s123 s23 s345 s56+mt2 s12 s234 s345 s56-mt2 s123 s234 s345 s56-mt2 s123 s34 s345 s56-s123 s16 s34 s345 s56+mt2 s23 s34 s345 s56-s123 s23 s34 s345 s56+s123 s23 s345^2 s56+s123 s234 s345^2 s56+mt2 s12 s23 s45 s56-mt2 s12 s234 s45 s56+mt2 s123 s234 s45 s56+2 mt2 s12 s34 s45 s56-mt2 s123 s34 s45 s56+2 s12 s16 s34 s45 s56-s123 s16 s34 s45 s56+mt2 s23 s34 s45 s56+2 s12 s23 s34 s45 s56-s12 s23 s345 s45 s56-s12 s234 s345 s45 s56-s123 s234 s345 s45 s56+2 s123 s34 s345 s45 s56-s23 s34 s345 s45 s56+s12 s234 s45^2 s56-2 s12 s34 s45^2 s56+mt2 s23 s345 s56^2-s23 s345^2 s56^2-mt2 s23 s45 s56^2+s23 s345 s45 s56^2)]


DLogArg=DLog/.Log[x_]:>x


DLogArg/possible6ptfacs[[1]]
DLogArg/possible6ptfacs[[1]]*possible6ptfacs[[2]]


Factor[DLogArg/.KinematicRulesBackWards/.SliceReplsNew,Modulus->$P]

Factor[possible6ptfacs[[1]]/.KinematicRulesBackWards/.SliceReplsNew,Modulus->$P]
Factor[possible6ptfacs[[2]]/.KinematicRulesBackWards/.SliceReplsNew,Modulus->$P]


Factor[DLogArg/.SliceTpoint,Modulus->$P]

Factor[possible6ptfacs[[1]]/.SliceTpoint,Modulus->$P]
Factor[possible6ptfacs[[2]]/.SliceTpoint,Modulus->$P]


SlicedFactors=Factor[possible6ptfacs/.SliceTpoint,Modulus->$P]
SortBy[Map[FactorList[#][[All,1]]&,SlicedFactors]//Flatten//Tally,Last]


SlicedFactors=Factor[possible6ptfacs/.KinematicRulesBackWards/.SliceReplsNew,Modulus->$P]
SortBy[Map[FactorList[#][[All,1]]&,SlicedFactors]//Flatten//Tally,Last]


SlicedFactors[[1]]
SlicedFactors[[6]]


Factor[possible6ptfacs[[1]]/.KinematicRulesBackWards]
Factor[possible6ptfacs[[2]]/.KinematicRulesBackWards]
Factor[possible6ptfacs[[3]]/.KinematicRulesBackWards]
Factor[possible6ptfacs[[4]]/.KinematicRulesBackWards]
Factor[possible6ptfacs[[7]]/.KinematicRulesBackWards]
Factor[possible6ptfacs[[8]]/.KinematicRulesBackWards]


Factor[(s12^2 s23^2-2 s12^2 s23 s234+2 s12 s123 s23 s234+s12^2 s234^2-2 s12 s123 s234^2+s123^2 s234^2+2 s12 s123 s23 s34-2 s12 s23^2 s34+2 s12 s123 s234 s34-2 s123^2 s234 s34+2 s12 s23 s234 s34+2 s123 s23 s234 s34+s123^2 s34^2-2 s123 s23 s34^2+s23^2 s34^2-2 s12 s23^2 s56+2 s12 s23 s234 s56-2 s123 s23 s234 s56-4 s12 s23 s34 s56+2 s123 s23 s34 s56-2 s23^2 s34 s56+s23^2 s56^2)/.KinematicRulesBackWards/.SliceReplsNew,Modulus->$P]


Expand[possible6ptfacs[[6]]/.KinematicRulesBackWards/.SliceReplsNew,Modulus->$P]
Factor[possible6ptfacs[[6]]/.KinematicRulesBackWards]/.Thread[Rule[variablesNew,t]]


NewDlog=(possible6ptfacs[[1]]possible6ptfacs[[7]]*(s123-s23))/(possible6ptfacs[[3]]possible6ptfacs[[5]])


slicedLogArg=Factor[NewDlog/.KinematicRulesBackWards/.SliceReplsNew,Modulus->$P]


NewFacs=DeleteElements[FactorList[Times@@possible6ptfacs/.KinematicRulesBackWards][[All,1]],{1}]
slicedNewFacs=Factor[NewFacs/.KinematicRulesBackWards/.SliceReplsNew,Modulus->$P]
slicedNewFacs//Length


Expand[slicedLogArg,Modulus->$P]


Expand[Factor[slicedLogArg/slicedNewFacs[[1]]^2,Modulus->$P],Modulus->$P]
Expand[Factor[slicedLogArg/(slicedNewFacs[[1]]^2slicedNewFacs[[8]]),Modulus->$P],Modulus->$P]
Expand[Factor[slicedLogArg/(slicedNewFacs[[1]]^2slicedNewFacs[[8]])*(slicedNewFacs[[2]]^2slicedNewFacs[[4]]),Modulus->$P],Modulus->$P]


A=Expand[Factor[Factor[(NewFacs[[1]]^2NewFacs[[8]])/(NewFacs[[2]]^2NewFacs[[4]])]/.KinematicRulesForward]/.Tpoint,Modulus->$P]
B=Expand[NewDlog/.Tpoint,Modulus->$P]


(*Factor[(NewFacs[[1]]^2NewFacs[[8]])/(NewFacs[[2]]^2NewFacs[[4]])]/.tpointNewPara*)
Expand[Factor[NewDlog/.KinematicRulesBackWards]/.tpointNewPara,Modulus->$P]
Expand[Factor[NewDlog]/.Tpoint,Modulus->$P]


A=Factor[D[Log[slicedLogArg],t],Modulus->$P]

B=Factor[D[Log[Factor[(NewFacs[[1]]^2NewFacs[[8]])/(NewFacs[[2]]^2NewFacs[[4]])]/.SliceReplsNew],t],Modulus->$P]


Factor[A-B,Modulus->$P]


RecoFac=Factor[(NewFacs[[1]]^2NewFacs[[8]])/(NewFacs[[2]]^2NewFacs[[4]])]

Factor[NewDlog/.KinematicRulesBackWards]

FacsNewPara=Factor[possible6ptfacs/.KinematicRulesBackWards];


FacsNewPara[[8]]


Expand[(possible6ptfacs[[1]]possible6ptfacs[[8]])/(possible6ptfacs[[3]]possible6ptfacs[[4]])/.Tpoint,Modulus->$P]


Expand[(possible6ptfacs[[1]]possible6ptfacs[[7]])/(possible6ptfacs[[3]]possible6ptfacs[[5]])/.Tpoint,Modulus->$P]


Expand[NewDlog/(s123-s23)/.Tpoint,Modulus->$P]





(* ::Section:: *)
(*Solve Func:*)


Get[CurrentDir<>"FactorAlphabetSolver.wl"]


NewDlog=(possible6ptfacs[[1]]possible6ptfacs[[7]]possible6ptfacs[[2]]^2 )/(possible6ptfacs[[3]]possible6ptfacs[[5]]^2  possible6ptfacs[[6]])
RecoFac=Factor[NewDlog/.KinematicRulesBackWards]


Res=SolveFactorPowers[RecoFac,FacsNewPara]


Factor[NewDlog]

TestLog2=Times@@Power[possible6ptfacs,Res[["exponents"]][[All,2]]]


Factor[NewDlog/.SliceTpoint,Modulus->$P]
Factor[TestLog2/.SliceTpoint,Modulus->$P]


TestLog2


(* ::Section:: *)
(*Reconstruction Test:*)


tvars=Table[tvar=RandomInteger[$P],{i,1,100}];
phasespacepoints=Table[Expand[SliceReplsNew/.t->tvars[[q]],Modulus->$P],{q,1,100}];

mt2repls=Expand[KinematicRulesBackWards[[1]]/.phasespacepoints,Modulus->$P];
phasespacepointsMands=Table[Thread[Rule[diffInvariants,Expand[diffInvariants/.KinematicRulesBackWards/.phasespacepoints[[i]]/.mt2repls[[i]],Modulus->$P]]],{i,1,100}];


func=Log[TestLog2];

DiffNumbers=Table[EQN=Expand[-1 Sum[D[G6,var]ds[ToExpression[Capitalize[ToString[var]]]],{var,diffInvariants}]/.phasespacepointsMands[[q]],Modulus->$P];
sixpfacs=Coefficient[Solve[EQN==0,ds[S16],Modulus->$P][[1,1,2]],Table[ds[ToExpression[Capitalize[ToString[var]]]],{var,diffInvariants}]];
Expand[Table[D[func,diffInvariants[[i]]]+sixpfacs[[i]]D[func,s16],{i,1,Length[sixpfacs]}]/.phasespacepointsMands[[q]],Modulus->$P],{q,1,100}];


NewDiff=Table[J=Expand[Table[D[KinematicRulesBackWards[[i,2]],var],{i,{1,2,3,5,6,7,8,9,10}},{var,variablesNew}]/.phasespacepoints[[q]],Modulus->$P];
f=DiffNumbers[[q]][[Flatten[Map[Position[diffInvariants,#]&,{mt2,s12,s123,s23,s234,s34,s345,s45,s56}]]]];
Expand[Transpose[J] . f,Modulus->$P],{q,1,100}];


(* ::Subsection:: *)
(*Test*)


NewDiff[[1]]


Expand[Map[D[func/.KinematicRulesBackWards,#]&,variablesNew]/.phasespacepoints[[1]],Modulus->$P]


phasespacepointsT=Table[Join[phasespacepoints[[q]],{tvar->tvars[[q]]}],{q,1,100}];





(* ::Subsection:: *)
(*Reco*)


ClearAll[tvar]
res=UniVarReconstructFUNCS[NewDiff[[All,1]],phasespacepointsT,tvar,20,20,{},1,1]


Factor[res[[1]],Modulus->$P]
Expand[res[[1]],Modulus->$P]


(* ::Input:: *)
(*(570823070 (1662747555+t) (710764459+1145249092 t+t^2) (1856782689+1742146491 t+372830817 t^2+t^3) 1(1165154409+1963078353 t+1164922417 t^2+1356487780 t^3+1059259609 t^4+1052767728 t^5+1182337599 t^6+1960635784 t^7+t^8))/((910208014+t) (1275985884+t) (1491081406+109681491 t+t^2) (29223145+683595211 t+t^2) (1650960977+1696743172 t+869430400 t^2+t^3) (478983085+216750222 t+1186479120 t^2+t^3) (1935502469+1190678798 t+1422393615 t^2+t^3))*)


Factor[D[func/.KinematicRulesBackWards,variablesNew[[1]]]/.SliceReplsNew,Modulus->$P]


Expand[Denominator[Factor[D[func,mt2]]]][[1]]/.{Thread[Rule[invariants,t]]}


(* ::Section:: *)
(*Extra Function*)


ConvertPointToMands[tpoint_]:=Thread[Rule[diffInvariants,Expand[{mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56}/.KinematicRulesBackWards/.tpoint/.Expand[KinematicRulesBackWards[[1]]/.tpoint,Modulus->$P],Modulus->$P]]]
ConvertDF[dfunc_,tpoint_]:=Module[{J,f},
J=Expand[Table[D[KinematicRulesBackWards[[i,2]],var],{i,{1,2,3,5,6,7,8,9,10}},{var,{s12,s23,s34,s123,s234,s56,c1,c2,c3}}]/.tpoint,Modulus->$P];
f=dfunc[[Flatten[Map[Position[{mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56},#]&,{mt2,s12,s123,s23,s234,s34,s345,s45,s56}]]]];
Expand[Transpose[J] . f,Modulus->$P]
]


diffInvariants
variablesNew


ConvertPointToMands[phasespacepoints[[5]]]
ConvertDF[DiffNumbers[[5]],phasespacepoints[[5]]]
NewDiff[[5]]


DiffNumbers[[All,4]]=0;
DiffNumbers


KinematicRulesForward


variablesNew


roots=ReadList[HomeDir<>"/Reconstruction1/Roots.txt"];


Expand[roots/.phasespacepointsMands[[1]],Modulus->$P]/.n_Integer:>PowerMod[n,1/2,$P]


Factor[roots/.KinematicRulesBackWards]//Denominator


JacobiSymbol[Expand[Factor[roots[[3]]/.KinematicRulesBackWards]/.phasespacepoints,Modulus->$P],$P]


variablesNew


TestGen[]:=Module[{},
For[i=1,i<=10000,i++,
tvars=RandomInteger[$P];
phasespacepoints=Expand[SliceReplsNew/.t->tvars,Modulus->$P];
If[Count[JacobiSymbol[Expand[Factor[roots/.KinematicRulesBackWards]/.phasespacepoints,Modulus->$P],$P],1]==15,Return[tvars]]]];


Length[roots]


TestGen[]


testpoint=Expand[SliceReplsNew/.t->108833581,Modulus->$P]


PowerMod[Expand[roots/.ConvertPointToMands[testpoint],Modulus->$P],1/2,$P]


rootsPara=Factor[roots/.KinematicRulesBackWards]


Export[CurrentDir<>"rootsPara.txt",rootsPara]


Export[CurrentDir<>"slice_a.txt", {StringRiffle[ToString /@ tpointNewPara1[[All,2]], " "]}, "Lines"];
Export[CurrentDir<>"slice_b.txt", {StringRiffle[ToString /@ tpointNewPara2[[All,2]], " "]}, "Lines"];
NParrallel=10
NSamples=100
PointGeneratorPath=HomeDir<>"/SQRT-Generator";
Command=PointGeneratorPath<>"/gen_param --roots "<>CurrentDir<>"rootsPara.txt --npoints "<>ToString[NSamples]<>" --slice slice_a.txt slice_b.txt"<>" --out points.txt --threads "<>ToString[NParrallel];


Command
Run[Command]


Directory[]
