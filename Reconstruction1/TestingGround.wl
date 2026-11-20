(* ::Package:: *)

If[$FrontEnd===Null,$DIR=Directory[],$DIR=NotebookDirectory[]];
parts = FileNameSplit[$DIR];
pos = First@FirstPosition[parts, "TriangleHeptagon3"];
HomeDir = FileNameJoin[parts[[;;pos ]]];
SectorDir = FileNameJoin[parts[[;;pos+1]]];
CurrentDir=$DIR;
kirapath=SectorDir;
$P=2^31-1;
PureIntegrals=ReadList[FileNameJoin[{SectorDir,"SubBasis.txt"}]];
(*PureIntegrals[[;;52]]=ConstantArray[0,52];*)
SetDirectory[CurrentDir];
NotbookDir=NotebookDirectory[];
Get[FileNameJoin[{HomeDir,"files/ExtraFunctions.m"}]];
SQRTS=Cases[PureIntegrals,_root,Infinity]//DeleteDuplicates;
(*SQRTPOWS=Map[Exponent[#,SQRTS]&,PureIntegrals];
SQRTVEC=Map[Times@@Power[SQRTS,#]&,SQRTPOWS];*)
OutPutPath="Reconstruction1/Samples/UnivariateSlice-3";
(*OutPutPath="Reconstruction1/Samples/Univariate Slice";*)
phasespacepoints=ReadList[FileNameJoin[{HomeDir,ToString[StringForm["``/phasespacepoints.txt",OutPutPath]]}]];
tvarlist=Flatten[ReadList[FileNameJoin[{HomeDir,ToString[StringForm["``/tvars.txt",OutPutPath]]}]][[;;;; 1]]];
phasespacepoints=phasespacepoints[[All,2;;]][[1;; ;; 1]];
diffInvariants={mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56};
PhaseSpacePoints=Map[Thread[Rule[diffInvariants,#]]&,phasespacepoints];
PhaseSpacePoints=Table[Join[PhaseSpacePoints[[i]],{tvarlist[[i]]}],{i,1,Length[PhaseSpacePoints]}];
(*PhaseSpacePoints=Map[Thread[Rule[diffInvariants,#]]&,phasespacepoints];*)
rootSet=(Cases[PureIntegrals[[-5;;]],_root,Infinity]//DeleteDuplicates)/.root[x_]:>x;
GetLinPos[x1_,x2_,NMI_]:=First[Flatten[Position[Flatten[Table[{i,j},{i,1,NMI},{j,1,NMI}],1],{x1,x2}]]];
Get[ToString[StringJoin[HomeDir,"/files/ExtraModules.m"]]];
Get[ToString[StringJoin[HomeDir,"/files/ExtraModules.m"]]];
usedvars={mt2,s12,s123,s23,s234,s34,s345,s45,s56};
For[i=1,i<=Length[usedvars],i++,
LinearTerm[usedvars[[i]]]=ReadList[FileNameJoin[{HomeDir,ToString[StringForm["``/LinTerm``_0001.txt",OutPutPath,usedvars[[i]]]]}]];
];


EQN=Flatten[Table[{avector[i]+PhaseSpacePoints[[1,-1,2]] bvector[i] ==PhaseSpacePoints[[1,i,2]],avector[i]+PhaseSpacePoints[[2,-1,2]] bvector[i] ==PhaseSpacePoints[[2,i,2]]},{i,1,Length[PhaseSpacePoints[[1,All,1]]]-1}]];
solset=Flatten[Table[{avector[i],bvector[i]},{i,1,Length[PhaseSpacePoints[[1,All,1]]]-1}]];
sol=Solve[EQN,solset,Modulus->$P][[1]];
avec=Thread[Rule[PhaseSpacePoints[[1,;;-2,1]],Table[avector[i],{i,1,Length[PhaseSpacePoints[[1,All,1]]]-1}] ]]/.sol
bvec=Thread[Rule[PhaseSpacePoints[[1,;;-2,1]],Table[bvector[i],{i,1,Length[PhaseSpacePoints[[1,All,1]]]-1}] ]]/.sol
LinearTerm[s16]=SparseArray[ConstantArray[0,{Length[PhaseSpacePoints],1,Length[PureIntegrals]^2}]];
SliceTpoint=Thread[Rule[bvec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]]
SliceTpoint=Thread[Rule[bvec[[All,1]],avec[[All,2]]+t bvec[[All,2]]]];
Expand[SliceTpoint/.t->PhaseSpacePoints[[-1,-1,2]],Modulus->$P]
PhaseSpacePoints[[-1]]
pos=1;


Get[HomeDir<>"/files/HexaBoxModules.m"]


vector1={p[1],p[2],p[3],p[4]}
vector2={p[2],p[3],p[4],p[5]}

Gram=Det[Expand[Expand[Table[2 vec1 vec2,{vec1,vector1},{vec2,vector1}]]/.mandelstamScalarProducts/.kinematics]]
Gram2=Det[Expand[Expand[Table[2 vec1 vec2,{vec1,vector2},{vec2,vector2}]]/.mandelstamScalarProducts/.kinematics]]


Factor[Gram/.SliceTpoint,Modulus->$P]


Test1=D[Log[Gram/Gram2],s12]
Test2=D[Log[Gram/Gram2],s16]
G6=Get[HomeDir<>"/SetUp/SixPointGram.txt"];
EqnConstFull=Sum[D[G6,var] ds[var],{var,diffInvariants}];
diffauxvar=First[Solve[EqnConstFull==0,ds[s16]]];
var=s12
sixpointdiff=Coefficient[diffauxvar[[1,2]],ds[var]]
Test=Factor[Test1+sixpointdiff Test2]



Factor[Test/.SliceTpoint,Modulus->$P]
Factor[Test/.RandomSlice,Modulus->$P]


RandomSlice=Thread[Rule[diffInvariants,RandomInteger[$P,10]+RandomInteger[$P,10]t]]


Letters=ReadList["/home/daniel/Documents/PHD-LocalFiles/TriangleHeptagon3/Reconstruction-NewPara/EvenLetterCandidates.txt"]


Test=Factor[Letters/.SliceTpoint,Modulus->$P];


SortBy[Flatten[Map[FactorList[#,Modulus->$P][[All,1]]&,Test]]//Tally,Last]


Position[Test,338898196+t]


Letters[[158]]
Letters[[165]]


Factor[Letters[[165]]/.SliceTpoint,Modulus->$P]



Factor[D[Log[Letters[[165]]/Letters[[158]]],s12]/.SliceTpoint,Modulus->$P]


Factor[D[Log[Letters[[165]]/Letters[[158]]],s12]]
