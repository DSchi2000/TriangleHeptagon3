(* ::Package:: *)

(* ::Title:: *)
(*Initialization Block*)


Get["IniBlocks.m"]
batchsize=5000;


(* ::Title:: *)
(*Extra Modules*)


$P=2^31-1;
GetEpsForAnsatz[pos_,var_,Dlength_,AnsatzTerms_]:=Module[{phasespacepointsperm,SolSet,Ansatz,MatPerm,Mat2,ep1,ep2,M1Mat,M0Mat,dVec,solution,fitvecs,recMat,recMatD,recMatN,inmatpos},
phasespacepointsperm=phasespacepoints;
MatPerm=dhxbMatrices[var];
dVec=phasespacepointsperm[[pos;;pos+Dlength-1]][[All,1]];
dVec=Expand[(4-dVec)/2,Modulus->$P];
fitvecs=MatPerm[[pos;;pos+Dlength-1]];
Ansatz=Sum[fac[i]AnsatzTerms[[i]],{i,1,Length[AnsatzTerms]}];
SolSet=Table[fac[i],{i,1,Length[AnsatzTerms]}];
solution=Table[eqn=Thread[Expand[fitvecs[[All,inmatpos]]==Map[Ansatz/.{x->#}&,dVec],Modulus->$P]];Sort[First[Solve[eqn,SolSet,Modulus->$P]]][[All,2]],{inmatpos,1,Length[MaxCutInts]^2}];
Return[solution]
 ]
GetEpsForFast[pos_,var_,Dlength_,AnsatzTerms_]:=Module[{Mat2,ns,res,NullRes,Mat,constpos,AnstatzFunc,Mats,epsVec,dvec,dir,numname},
numname=IntegerString[Floor[(pos-1)/1000],10,4];
dir=FileNameJoin[{Directory[],ToString[StringForm["``/tmp-files/DEMatsOut_``",OutPutPath,numname]]}];
dvec=ReadList[FileNameJoin[{dir,ToString[StringForm["Parameters-``.txt",IntegerString[pos,10,7]]]}]][[All,1]];
epsVec=Expand[(4-dvec)/2,Modulus->$P];
Mats=ReadList[FileNameJoin[{dir,ToString[StringForm["DEmatrix_d``-``.txt",var,IntegerString[pos,10,7]]]}]];
If[AnsatzTerms==={x},
Return[{Expand[Normal[Mats]/epsVec[[1]],Modulus->$P]}]
];
AnstatzFunc[eps_]:=AnsatzTerms/.x->eps;
constpos=First[Flatten[Position[AnsatzTerms,1]]];
Mat=Expand[Map[Join[AnstatzFunc[#]]&,epsVec],Modulus->$P];
NullRes=ConstantArray[0,Length[AnsatzTerms]+1];
res=Table[
If[Mats[[1,intpos]]===0,
NullRes,
Mat2=Expand[MapThread[Join,{Mat,Transpose[{-Mats[[All,intpos]]}]}],Modulus->$P];
ns=First[NullSpace[Mat2,Modulus->$P]];
ns=Expand[ns/ns[[-1]],Modulus->$P]],{intpos,1,Length[Mats[[1]]]}];
Return[res];
]          


(* ::Title:: *)
(*Reconstruct Eps Dependence*)


SetDirectory[CurrentDir];
diffInvariants={mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56};
ClearAll[dhxbMatrices];
Print["Loading Matrices"];
kernels=LaunchKernels[NParrallel];
Print["Reconstructing Eps"];

runtime=AbsoluteTiming[Table[
Mats=ParallelTable[GetEpsForFast[psp,var,EpsRecoLenth,{x}],{psp,1,NSamples}]//EchoTiming;
q=1;
Table[
(*Export[FileNameJoin[{Directory[],ToString[StringForm["``/ConstTerm``_``.txt",OutPutPath,var,IntegerString[q,10,4]]]}],SparseArray/@Mats[[i;;Min[{i+batchsize,NSamples}],All,1]]];*)
Export[FileNameJoin[{Directory[],ToString[StringForm["``/LinTerm``_``.txt",OutPutPath,var,IntegerString[q,10,4]]]}],SparseArray/@Mats[[i;;Min[{i+batchsize,NSamples}],All,1]]];
q+=1;
,{i,1,NSamples,batchsize}];
Print[ToString[StringForm["`` Done",var]]],
{var,SubdVars}];]


SetDirectory[CurrentDir];
Export["TIME.txt",{"Reconstruct eps Dependence",runtime[[1]]}];
Run[ToString[StringForm["cd `` && python3 mailsend.py",CurrentDir]]];
Print["Done"];
(*$P=2^31-1;
(*Expand[-1/2,M*)odulus->$P]
$P=2^31-1;
Expand[-1/3,Modulus->$P]*)
