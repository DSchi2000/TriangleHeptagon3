(* ::Package:: *)

(* ::Section:: *)
(*MandelStams and Kinematics*)


kinematics={p[1]^2->0,p[2]^2->0,p[3]^2->0,p[4]^2->0,p[5]^2->mt2,p[6]->mt2};
p4cons={p[6]->-p[1]-p[2]-p[3]-p[4]-p[5]};
mandelstamsList={s12-(p[1]+p[2])^2,s23-(p[2]+p[3])^2,
s34-(p[3]+p[4])^2,s45-(p[4]+p[5])^2,s56-(p[5]+p[6])^2,
s16-(p[1]+p[6])^2,s123-(p[1]+p[2]+p[3])^2,s234-(p[2]+p[3]+p[4])^2,
s345-(p[3]+p[4]+p[5])^2,mt2-p[6]^2};
SetAttributes[sp,Orderless]
mandlistexpanded=Expand[mandelstamsList/.p4cons]/.kinematics/.p[n_]p[m_]:>sp[n,m];
mandrepls=First[Solve[Thread[mandlistexpanded==0],DeleteDuplicates[Cases[mandlistexpanded,_sp,Infinity]]]/.sp[n_,m_]:>p[n]p[m]]
Export[NotebookDirectory[]<>"mandelstam_scalarproducts.txt",mandrepls];
Export[NotebookDirectory[]<>"mandelstam_kinematics.txt",kinematics];
invariants=Variables[mandrepls[[All,2]]]



(* ::Subsection:: *)
(*6PointGram*)


vecs={p[6],p[2],p[3],p[4],p[5]};
G6=Together[16*Det[Table[Expand[vec1 vec2/.p4cons]/.kinematics/.mandrepls,{vec1,vecs},{vec2,vecs}]]];
Export[NotebookDirectory[]<>"SixPointGram.txt",mandrepls];
$P=2^31-1
tpoint=Thread[Rule[DeleteElements[invariants,{s16}],RandomInteger[$P,9]]];
tpoint=Sort[Join[tpoint,{Solve[Expand[G6/.tpoint,Modulus->$P]==0,s16,Modulus->$P][[1,1]]}]]


TA=Together[G6/.mt2->0];
sixpg2= Get["/home/dschiebe/Documents/PhD-Files/FullHexaBoxDE/229-231-PentaBox_3MI/Reconstruction/SixPointGram.txt"];
TB=Together[16*sixpg2/.mt2->0];
TA-TB


(* ::Section:: *)
(*Propagators*)


propagatorsList={-PR[1]+k[1]^2,
-PR[2]+(k[1]+p[1])^2,
-PR[3]+(k[2]-k[1])^2,
-PR[4]+(k[2]+p[1])^2,
-PR[5]+(k[2]+p[1]+p[2])^2,
-PR[6]+(k[2]+p[1]+p[2]+p[3])^2,
-PR[7]+(k[2]+p[1]+p[2]+p[3]+p[4])^2,
-PR[8]+(k[2]+p[1]+p[2]+p[3]+p[4]+p[5])^2-mt2,
-PR[9]+k[2]^2
};
ISP={-PR[10]+(k[1]+p[2])^2,
-PR[11]+(k[1]+p[3])^2,
-PR[12]+(k[1]+p[4])^2,
-PR[13]+(k[1]+p[5])^2
};
propagatorsList=Join[propagatorsList,ISP];
propagatorsListExpanded=Expand[propagatorsList/.p4cons]/.kinematics/.mandrepls/.{k[n_]p[m_]:>sp[k[n],p[m]],k[n_]k[m_]:>sp[k[n],k[m]],k[n_]^2:>sp[k[n],k[n]]};
proprepls=First[Solve[Thread[propagatorsListExpanded==0],DeleteDuplicates[Cases[propagatorsListExpanded,_sp,Infinity]]]/.sp[n_,m_]:>n m]
Export[NotebookDirectory[]<>"loopmomentum_scalarproducts.txt",proprepls];



(* ::Subsection:: *)
(*Plot*)


PropVectorForm=Sqrt[propagatorsList[[;;9]]/.{PR[n_]:>0,mt2->0}]/.Sqrt[x_^2]:>x/.{k[1]->k1,k[2]->k2,p[1]->p1,p[2]->p2,p[3]->p3,p[4]->p4,p[5]->p5}
SetDirectory["/home/dschiebe/Documents/PhD-Files/FullHexaBoxDE/files"];
Get[FileNameJoin[{Directory[],"KinematicGraphs.m"}]];
buildKinematicGraph[PropVectorForm,{k1,k2},{p1,p2,p3,p4,p5}]
SetDirectory[NotebookDirectory[]];


(* ::Section:: *)
(*KIRA*)


StringRiffle[Join[{d},tpoint[[All,1]]]]
StringRiffle[Join[{RandomInteger[$P]},tpoint[[All,2]]]]



vars=invariants
Kiravars="";
For[i=1,i<=Length[vars],i++,
Kiravars=StringJoin[Kiravars,ToString[StringForm["     - [``, 2]\n",vars[[i]]]]]
];

KIRAKINEMATICS="";
For[i=1,i<=Length[kinematics],i++,
KIRAKINEMATICS=StringJoin[KIRAKINEMATICS,ToString[StringForm["     - [[p``,p``], ``]\n",kinematics[[i,1]]/.p[n_]^2:>n,kinematics[[i,1]]/.p[n_]^2:>n,ToString[CForm[kinematics[[i,2]]]]]]]
];
For[i=1,i<=Length[mandrepls],i++,
KIRAKINEMATICS=StringJoin[KIRAKINEMATICS,ToString[StringForm["     - [[p``,p``], ``]\n",mandrepls[[i,1,1]]/.p[n_]:>n,mandrepls[[i,1,2]]/.p[n_]:>n,StringDrop[ToString[CForm[mandrepls[[i,2]]]],-1]]]]
];
Kiravars
KIRAKINEMATICS

PropVectorForm=Sqrt[propagatorsList/.{PR[n_]:>0,mt2->0}]/.Sqrt[x_^2]:>x/.{k[1]->k1,k[2]->k2,p[1]->p1,p[2]->p2,p[3]->p3,p[4]->p4,p[5]->p5};
masses=Map[MemberQ[#//Variables,mt2]&,propagatorsList]/.{False->0,True->mt2};

Kiraprops="";
For[i=1,i<=Length[PropVectorForm],i++,
Kiraprops=StringJoin[Kiraprops,ToString[StringForm["      - [ \"`` \", ``]\n",PropVectorForm[[i]],masses[[i]]]]]
];
Kiraprops



(* ::Subsection:: *)
(*Load Kira*)


IBPTables=Get[NotebookDirectory[]<>"KIRAINI/results/HPT3/numerics_2147483647_0.m"];
IBPTablesRepls=IBPTables[[2,2]];
tpoint=Thread[Rule[invariants,IBPTables[[2,1,2;;]]]]





MIs=DeleteDuplicates[Cases[IBPTables[[2,2,All,2]],_HPT3,Infinity]];
MIs//Length


Propdiff[n_,q_,integral_]:=Module[{newint},
	newint=(integral/.{HPT3->List});
	newint[[n]]-=q;
	Return[newint/.{List->HPT3}];
]
HPT3/:Power[PR[n_],q_]*HPT3[n1_,n2_,n3_,n4_,n5_,n6_,n7_,n8_,n9_,n10_,n11_,n12_,n13_]:=Propdiff[n,q,HPT3[n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13]];
HPT3/:PR[n_]*HPT3[n1_,n2_,n3_,n4_,n5_,n6_,n7_,n8_,n9_,n10_,n11_,n12_,n13_]:=Propdiff[n,1,HPT3[n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13]];
HPT3/: delHPT3[x_]*HPT3[n1_,n2_,n3_,n4_,n5_,n6_,n7_,n8_,n9_,n10_,n11_,n12_,n13_] := Sum[-{n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13}[[i]]*PR[i]^-1 HPT3[n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13]*D[Pr[i], x], {i, 1, 13}];

vecs1={k[1],p[1],p[2],p[3],p[4]};
vecs2={p[5],p[1],p[2],p[3],p[4]};

G2L1=Together[64*Det[Table[Expand[vec1 vec2/.p4cons]/.kinematics/.mandrepls/.proprepls,{vec1,vecs1},{vec2,vecs2}]]];

vecs3={k[2],p[1],p[2],p[3],p[4]};
vecs2={p[5],p[1],p[2],p[3],p[4]};
G2L2=Together[64*Det[Table[Expand[vec1 vec2/.p4cons]/.kinematics/.mandrepls/.proprepls,{vec1,vecs3},{vec2,vecs2}]]];


EQN1=Expand[Expand[Expand[(Join[G2L1]/.Pr->PR)*MIs]/.tpoint,Modulus->$P]/.IBPTablesRepls,Modulus->$P];
EQN2=Expand[Expand[Expand[(Join[G2L2]/.Pr->PR)*MIs]/.tpoint,Modulus->$P]/.IBPTablesRepls,Modulus->$P];
EQN=Join[EQN1,EQN2];


EQNNEW=Expand[Delete[EQN,Position[Expand[Thread[EQN==0],Modulus->$P],True]]//EchoFunction[Length],Modulus->$P];
SUBMIS=Cases[EQN,_HPT3,Infinity]//DeleteDuplicates;
(*solset={HPT3[1,0,1,1,1,1,1,1,0,0,0,0,0],HPT3[0,1,1,0,1,1,1,1,1,0,0,0,0],HPT3[1,1,1,0,1,1,1,1,0,0,0,0,0]}*)
(*solset={HPT3[1,0,1,1,1,1,1,1,0,0,0,0,0],HPT3[0,1,1,0,1,1,1,1,1,0,0,0,0],HPT3[1,1,1,0,1,1,1,1,0,0,0,0,0]}*)
nullints=Complement[SUBMIS,MIs];
nullintsrpls=Thread[Rule[nullints,0]];
EQNNEW=EQNNEW/.nullintsrpls;
EQN=EQN/.nullintsrpls;
SUBMIS=Cases[EQN,_HPT3,Infinity]//DeleteDuplicates;
sol=First[Solve[Thread[EQNNEW==0],SUBMIS,Modulus->$P]//EchoTiming]
solset=sol[[All,1]]
sol//Length


EQNNEW=Expand[Delete[EQN,Position[Expand[Thread[EQN==0],Modulus->$P],True]]//EchoFunction[Length],Modulus->$P];
EQN1A=Expand[(Join[G2L1]/.Pr->PR)*MIs];
EQN2A=Expand[(Join[G2L2]/.Pr->PR)*MIs];
EQNA=Join[EQN1A,EQN2A];
EQNA2=Expand[Delete[EQNA,Position[Expand[Thread[EQN==0],Modulus->$P],True]]//EchoFunction[Length]]
Export[NotebookDirectory[]<>"LinearSystem.txt",EQNA2];



Export[NotebookDirectory[]<>"ExtraRelationsReduList.txt",DeleteDuplicates[Cases[EQNA2,_HPT3,Infinity]]]


MISNew=Cases[MIs/.sol,_HPT3,Infinity]//DeleteDuplicates//EchoFunction[Length]


MIs/.n_Integer/;(n<0):>0//Tally
MISNew/.sol/.n_Integer/;(n<0):>0//Tally
Export[NotebookDirectory[]<>"MIs4D.txt",MISNew];
Export[NotebookDirectory[]<>"MIsFull.txt",MIs];


Sectors=DeleteDuplicates[MISNew/.sol/.n_Integer/;(n<0):>0];
PropVectorForm=Sqrt[propagatorsList[[;;9]]/.{PR[n_]:>0,mt2->0}]/.Sqrt[x_^2]:>x/.{k[1]->k1,k[2]->k2,p[1]->p1,p[2]->p2,p[3]->p3,p[4]->p4,p[5]->p5}
SetDirectory["/home/dschiebe/Documents/PhD-Files/FullHexaBoxDE/files"];
Get[FileNameJoin[{Directory[],"KinematicGraphs.m"}]];
buildKinematicGraph[PropVectorForm[[Flatten[Position[(Sectors[[-10]]/.HPT3:>List)[[;;9]],1]]]],{k1,k2},{p1,p2,p3,p4,p5}]
SetDirectory[NotebookDirectory[]];


Export[NotebookDirectory[]<>"Sectors.txt",DeleteDuplicates[MISNew/.sol/.n_Integer/;(n<0):>0]]


sectorPos=Map[Position[MISNew/.n_Integer/;(n<0):>0,#]&,Sectors]


GenerateTableEntry[subsectorpos_]:=Module[{TableSting,MIsequence,MICount,MIposition},
MIsequence=StringJoin[ToString/@(Sectors[[subsectorpos]]/.HPT3->List)];
MICount=Length[sectorPos[[subsectorpos]]];
If[MICount===1,
MIposition=ToString[sectorPos[[subsectorpos,1,1]]],
MIposition=ToString[sectorPos[[subsectorpos,1,1]]]<>"-"<>ToString[sectorPos[[subsectorpos,-1,1]]];
];


TableSting="{}&"<>MIposition<>"&"<>MIsequence<>"& "<>ToString[MICount]<>" &  {\\tabfigure{height=0.14\\textheight}{HPT3/"<>MIsequence<>".pdf}}&\makecell[l]{}\\\\ \hline";
Return[TableSting];
]


TablePos=7;
(TablePos-1)*9+1
(TablePos-1)*9+9
Export[NotebookDirectory[]<>"LateXTable.txt",Table[GenerateTableEntry[i],{i,(TablePos-1)*9+1,(TablePos-1)*9+9}]]


Export[NotebookDirectory[]<>"LateXTable.txt",Table[GenerateTableEntry[i],{i,64,64}]]


Sectors//Length
