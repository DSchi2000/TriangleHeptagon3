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


KinematicRulesForward=Join[cvals[[;;-2]]]

KinematicRulesBackWards=SortBy[Join[sol10thmand,solmands],First]


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


K=Together[Table[D[KinematicRulesForward[[i,2]],var],{i,1,3},{var,{mt2,s345,s45}}]]


KN=Expand[K/.Tpoint,Modulus->$P]
KNINV=Inverse[KN,Modulus->$P]
f=DiffNumbers[[Flatten[Map[Position[diffInvariants,#]&,{mt2,s345,s45}]]]]
Expand[KNINV . f,Modulus->$P]


func2=Together[func/.KinematicRulesBackWards]
DiffNumbers2=Expand[Table[D[func,diffInvariants[[i]]]+sixpfacs[[i]]D[func,s16],{i,1,Length[sixpfacs]}]/.Tpoint,Modulus->$P];
