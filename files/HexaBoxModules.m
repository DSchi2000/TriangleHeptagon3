(* ::Package:: *)

pathHBM=DirectoryName[$InputFileName];


kinematics=ReadList[FileNameJoin[{HomeDir,"SetUp/mandelstam_kinematics.txt"}]];
p4cons={p[6]->-p[1]-p[2]-p[3]-p[4]-p[5]};
loopMomentumScalarProducts=ReadList[FileNameJoin[{HomeDir,"SetUp/loopmomentum_scalarproducts.txt"}]];
mandelstamScalarProducts=ReadList[FileNameJoin[{HomeDir,"SetUp/mandelstam_scalarproducts.txt"}]];
kinematicFunctions=ReadList[FileNameJoin[{pathHBM,"kinematicFunction.txt"}]];
LaportaBasis=ReadList[FileNameJoin[{HomeDir,"SetUp/MIs4D.txt"}]];
LinSysAna=ReadList[FileNameJoin[{HomeDir,"SetUp/LinearSystem.txt"}]];
ExtraRelationIntList=DeleteDuplicates[Cases[LinSysAna,_HPT3,Infinity]];
allMIs=ReadList[FileNameJoin[{HomeDir,"SetUp/MIsFull.txt"}]];
invariants={d,mt2,s12,s123,s16,s23,s234,s34,s345,s45,s56};


auxvar=s16;


(* ::Title:: *)
(*Generate Random PhaseSpace Point*)


momenta={p[1],p[2],p[3],p[4],p[5]};
SixPointGram=Expand[Det[Table[momenta[[i]] momenta[[j]],{i,1,5},{j,1,5}]/.kinematics/.mandelstamScalarProducts/.loopMomentumScalarProducts]];
(*Export[NotebookDirectory[]<>"SixPointGram.txt",SixPointGram];*)
(*sqrtSet={tr5Exp,sqrtG3Exp,sqrtG3DB}/.Sqrt[x_]:>x*)
GetRandomPT[NMin_,NMax_]:=Module[{nestedroots,tPoint,tPointReplacements,sqrtSetT,roots,NMAXIT,s16Sol,nestedrootsT},    
    NMAXIT=1000000;
    nestedroots=DeleteDuplicates[Cases[sqrtSet,_root,Infinity]]/.root[x_]:>x;
    For[i = 0, i < NMAXIT, i++,
		tPoint=RandomInteger[{NMin,NMax},Length[DeleteElements[invariants,{auxvar}]]];
		tPointReplacements=Thread[Rule[DeleteElements[invariants,{auxvar}],tPoint]];
		s16Sol=Solve[Expand[SixPointGram/.tPointReplacements,Modulus->$P]==0,auxvar,Modulus->$P];
		If[Length[s16Sol]>0, 
		nestedrootsT = Expand[nestedroots /. Join[tPointReplacements,{s16Sol[[1,1]]}], Modulus -> $P];
        roots = JacobiSymbol[#, $P] & /@ nestedrootsT;
        If[Total[roots] == Length[roots],
            sqrtSetT = Expand[sqrtSet /. Join[tPointReplacements,{s16Sol[[1,1]]}]//.root[n_]:>PowerMod[n,1/2,$P], Modulus -> $P];
	        roots = JacobiSymbol[#, $P] & /@ sqrtSetT;
	        If[Total[roots] == Length[roots],
	            Return[SortBy[Join[tPointReplacements,{s16Sol[[1,1]]}],First]];];
            ];
		];
 ];
 Return[{}];
]
GetUnivariantSlice5PT[a_,b_]:=Module[{tPoint,tPointReplacements,sqrtSetT,roots,NMAXIT,s16Sol,tVal,nestedroots,nestedrootsT},    
    NMAXIT=1000000;
    nestedroots=DeleteDuplicates[Cases[sqrtSet,_root,Infinity]]/.root[x_]:>x;
    For[i = 0, i < NMAXIT, i++,
		tVal=RandomInteger[{0,$P}];
		tPoint=Expand[a+tVal b,Modulus->$P];
		tPointReplacements=Thread[Rule[DeleteElements[invariants,{auxvar}],tPoint]];
		tPointReplacements=Join[tPointReplacements,{tvar->tVal}];
		s16Sol=Solve[Expand[SixPointGram/.tPointReplacements,Modulus->$P]==0,auxvar,Modulus->$P];
		If[Length[s16Sol]>0, 
		nestedrootsT = Expand[nestedroots /. Join[tPointReplacements,{s16Sol[[1,1]]}], Modulus -> $P];
        roots = JacobiSymbol[#, $P] & /@ nestedrootsT;
        If[Total[roots] == Length[roots],
            sqrtSetT = Expand[sqrtSet /. Join[tPointReplacements,{s16Sol[[1,1]]}]//.root[n_]:>PowerMod[n,1/2,$P], Modulus -> $P];
	        roots = JacobiSymbol[#, $P] & /@ sqrtSetT;
	        If[Total[roots] == Length[roots],
	            Return[SortBy[Join[tPointReplacements,{s16Sol[[1,1]]}],First]];];
            ];
		];
 ];
 Return[{}];
]
GetOneVarVariance5PT[a_,svar_]:=Module[{tPoint,tPointReplacements,sqrtSetT,roots,NMAXIT,s16Sol,tVal,varpos,nestedrootsT,nestedroots},    
    NMAXIT=1000000;
    nestedroots=DeleteDuplicates[Cases[sqrtSet,_root,Infinity]]/.root[x_]:>x;
(*    auxvar//Echo;
    a[[All,1]]//Echo;
    Delete[a[[All,1]],Position[a[[All,1]],auxvar]]//Echo;
    Return[{}];*)
    For[i = 0, i < NMAXIT, i++,
        varpos=Position[Delete[a[[All,1]],Position[a[[All,1]],auxvar]],svar][[1,1]];
		tVal=RandomInteger[{0,$P}];
		tPoint=a[[All,2]];
		tPoint[[varpos]]=tVal;
		tPointReplacements=Thread[Rule[DeleteElements[invariants,{auxvar}],tPoint]];
		s16Sol=Solve[Expand[SixPointGram/.tPointReplacements,Modulus->$P]==0,auxvar,Modulus->$P];
		If[Length[s16Sol]>0, 
		nestedrootsT = Expand[nestedroots /. Join[tPointReplacements,{s16Sol[[1,1]]}], Modulus -> $P];
        roots = JacobiSymbol[#, $P] & /@ nestedrootsT;
        If[Total[roots] == Length[roots],
            sqrtSetT = Expand[sqrtSet /. Join[tPointReplacements,{s16Sol[[1,1]]}]//.root[n_]:>PowerMod[n,1/2,$P], Modulus -> $P];
	        roots = JacobiSymbol[#, $P] & /@ sqrtSetT;
	        If[Total[roots] == Length[roots],
	            Return[SortBy[Join[tPointReplacements,{s16Sol[[1,1]]}],First]];];
            ];
		];
 ];
 Return[{}];
]


(* ::Title:: *)
(*KIRA SET UP*)


DeleteOldKiraRun[]:=Module[{},
	Quiet[DeleteFile[FileNameJoin[{Directory[],"KIRA/kira.log"}]]];
	Quiet[DeleteDirectory[FileNameJoin[{Directory[],"KIRA/tmp"}],DeleteContents -> True]];
	Quiet[DeleteDirectory[FileNameJoin[{Directory[],"KIRA/sectormappings"}],DeleteContents -> True]];
	Quiet[DeleteDirectory[FileNameJoin[{Directory[],"KIRA/results"}],DeleteContents -> True]];
	Quiet[DeleteDirectory[FileNameJoin[{Directory[],"KIRA/splitres"}],DeleteContents -> True]];
];
WriteNumericsFile[tPointReplacements_]:=Module[{vars,vals,file,varsString,valsString},
	Quiet[DeleteFile[FileNameJoin[{Directory[],"KIRA/numerics"}]]];
	vars=tPointReplacements[[All,1]];
	vals=tPointReplacements[[All,2]];
	varsString=StringDrop[ToString[StringJoin[Table[StringJoin[ToString[vars[[i]]]," "],{i,1,Length[vars]}]]],-1];
	valsString=StringDrop[ToString[StringJoin[Table[StringJoin[ToString[vals[[i]]]," "],{i,1,Length[vals]}]]],-1];	
	file=OpenAppend[FileNameJoin[{Directory[],"KIRA/numerics"}]];
	WriteString[file,ToString[StringForm["prime ``\n",$P]]];
	WriteString[file,ToString[StringForm["``\n",varsString]]];
	WriteString[file,ToString[StringForm["``\n",valsString]]];
	Close[file];
];
AddToNumericsFile[tPointReplacements_]:=Module[{vars,vals,file,varsString,valsString},
	vars=tPointReplacements[[All,1]];
	vals=tPointReplacements[[All,2]];
	varsString=StringDrop[ToString[StringJoin[Table[StringJoin[ToString[vars[[i]]]," "],{i,1,Length[vars]}]]],-1];
	valsString=StringDrop[ToString[StringJoin[Table[StringJoin[ToString[vals[[i]]]," "],{i,1,Length[vals]}]]],-1];		
	file=OpenAppend[FileNameJoin[{Directory[],"KIRA/numerics"}]];
	WriteString[file,ToString[StringForm["``\n",valsString]]];
	Close[file];
];


GenerateSectorsFromReductionList[filePath_]:=Module[{ReduList,ReduListT,Sectors,SectorInts,r,s,d,sectorString,file,pos,sectorpath,directorypath},
ReduList=ReadList[filePath];
ReduListT=ReduList/.HPT3[n1_,n2_,n3_,n4_,n5_,n6_,n7_,n8_,n9_,n10_,n11_,n12_,n13_]:>{n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13};
Sectors=ReduListT/.( n_Integer:>0/;n<0)/.( n_Integer:>1/;n>1) //DeleteDuplicates;
directorypath=StringDrop[filePath,Last[Last[StringPosition[filePath,"/"]]];;];
sectorpath=StringJoin[{directorypath,"/",FileBaseName[filePath],"_sector"}];
Quiet[DeleteFile[sectorpath]];

For[pos=1,pos<=Length[Sectors],pos++,
	Sectors[[pos]];
	SectorInts=ReduListT[[Sequence@@#]]&/@Position[ReduListT/.( n_Integer:>0/;n<0)/.( n_Integer:>1/;n>1) ,Sectors[[pos]]];
	SectorInts;
	r=Max[Table[Total[Select[SectorInts[[i]],Positive]],{i,1,Length[SectorInts]}]];
	s=Max[-Table[Total[Select[SectorInts[[i]],Negative]],{i,1,Length[SectorInts]}]];
	d=Max[Total[Thread[((SectorInts-1)/.( n_Integer:>0/;n<0))]]];
	sectorString=StringForm["b``",ToString[StringJoin[Table[ToString[Sectors[[pos,i]]],{i,1,Length[Sectors[[pos]]]}]]]];
	If[ToString[sectorString]=="b1011001100000", s=3;];
	If[ToString[sectorString]=="b1011010100000", s=3;];

	file=OpenAppend[sectorpath];
	WriteString[file,ToString[StringForm["      - {topologies: [HPT3], sectors: [``],r: ``,s: ``,d: ``}\n",sectorString,r,s,d]]];
	Close[file];
];
];
WriteJobFile[redulists_]:=Module[{filePath,directorypath,sectorpath,lists,redufilesectors,readline,i},
	filePath=FileNameJoin[{Directory[],"KIRA/jobs.yaml"}];
	directorypath=StringDrop[filePath,Last[Last[StringPosition[filePath,"/"]]];;];
	DeleteFile[filePath];
	file=OpenAppend[filePath];
	WriteString[file,"jobs:\n"];
	WriteString[file," - reduce_sectors:\n"];
	WriteString[file,"    reduce:\n"];
	For[lists=1,lists<=Length[redulists],lists++,
		redufilesectors=OpenRead[StringJoin[{directorypath,"/",redulists[[lists]],"_sector"}]];
		For[i=1,i<=10000,i++,
			readline=ReadLine[redufilesectors];
			If[readline===EndOfFile,
				Break[];
				Close[redufilesectors];,
				WriteString[file,StringJoin[readline,"\n"]];
			];
			];
				
	];
	WriteString[file,"    select_integrals:\n"];
	WriteString[file,"     select_mandatory_list:\n"];
	For[lists=1,lists<=Length[redulists],lists++,
		WriteString[file,ToString[StringForm["      - [``.txt]\n",redulists[[lists]]]]];	
	];
	WriteString[file,"    run_initiate: true\n"];
	WriteString[file,"    numerical_points: numerics\n"];
	WriteString[file,"    run_firefly: true\n"];
];


(* ::Title:: *)
(*Solve LinearSystem from Extra Relations*)


GetExtraRelationReplacementsFromPath[numericalpath_,tPointReplacements_]:=Module[{ReduList,LinSys,Reduction,subset,Remove1L,solset,linsysBasis,sortedBasis,LinSysMatrix,NS,RemBasisN,RemBasis,SubsInts,D4ExtraRelations},
	LinSys=LinSysAna/.tPointReplacements//EchoFunction[Length];
	Reduction=GetReductionTableFromFile[numericalpath];
	Return[GetExtraRelationReplacements[Reduction,tPointReplacements]];
];
GetExtraRelationReplacements[Reduction_,tPointReplacements_]:=Module[{ReduList,LinSys,subset,Remove1L,solset,linsysBasis,sortedBasis,LinSysMatrix,NS,RemBasisN,RemBasis,SubsInts,D4ExtraRelations},
	LinSys=LinSysAna/.tPointReplacements;
	LinSys=Expand[LinSys/.Reduction,Modulus->$P];
	subset=Cases[LinSys,HPT3[___],Infinity]//DeleteDuplicates;
	Remove1L=Thread[Rule[Complement[subset,allMIs],0]];
	LinSys=LinSys/.Remove1L;
	solset={HPT3[1,0,1,1,1,1,1,1,0,0,0,0,0],HPT3[0,1,1,0,1,1,1,1,1,0,0,0,0],HPT3[1,1,1,0,1,1,1,1,0,0,0,0,0]};
	linsysBasis=Cases[LinSys,_HPT3,Infinity]//DeleteDuplicates;
	sortedBasis=Join[solset,Complement[linsysBasis,solset]];
	LinSysMatrix=Table[Coefficient[LinSys[[j]],sortedBasis[[i]]],{j,1,Length[LinSys]},{i,1,Length[linsysBasis]}];
	NS=NullSpace[LinSysMatrix,Modulus->$P]; 
	RemBasisN=Table[Last[Position[NS[[i]],1]],{i,1,Length[NS]}];
	RemBasis=sortedBasis[[Sequence@@#]]&/@RemBasisN;
	SubsInts=Delete[sortedBasis,RemBasisN];
	D4ExtraRelations=Thread[Rule[SubsInts,Table[Total[RemBasis*NS[[All,i]]],{i,1,Length[SubsInts]}]]];
	Return[D4ExtraRelations]
];
GetReductionTableFromFile[numericalpath_]:=Module[{ReduList,Reduction},
	ReduList=Get[numericalpath];
	Reduction=Table[Rule[ReduList[[2,2,pos,1]],Sum[ReduList[[2,2,pos,i,1]]*ReduList[[2,2,pos,i,2]],{i,2,Length[ReduList[[2,2,pos]]]}]],{pos,1,Length[ReduList[[2,2]]]}];
	Return[Reduction];
];
GetReductionTableFromList[ReduSublist_]:=Module[{Reduction},
	Reduction=Table[Rule[ReduSublist[[pos,1]],Sum[ReduSublist[[pos,i,1]]*ReduSublist[[pos,i,2]],{i,2,Length[ReduSublist[[pos]]]}]],{pos,1,Length[ReduSublist]}];
	Return[Reduction];
];


(* ::Title:: *)
(*Propagators*)


Pr[1]=(k[1])^2;
Pr[2]=(k[1]+p[1])^2;
Pr[3]=((k[2]-k[1])^2);
Pr[4]=(k[2]+p[1])^2;
Pr[5]=(k[2]+p[1]+p[2])^2;
Pr[6]=(k[2]+p[1]+p[2]+p[3])^2;
Pr[7]=(k[2]+p[1]+p[2]+p[3]+p[4])^2;
Pr[8]=(k[2]+p[1]+p[2]+p[3]+p[4]+p[5])^2-mt2;
Pr[9]=(k[2])^2;
Pr[10]=(k[1]+p[2])^2;
Pr[11]=(k[1]+p[3])^2;
Pr[12]=(k[1]+p[4])^2;
Pr[13]=(k[1]+p[5])^2;

InternalMometa={k[1],k[2]}
ExternalMometa={p[1],p[2],p[3],p[4],p[5]}


(* ::Title:: *)
(*HPT3 Definitions*)


Propdiff[n_,q_,integral_]:=Module[{newint},
	newint=(integral/.{HPT3->List});
	newint[[n]]-=q;
	Return[newint/.{List->HPT3}];
]
HPT3/:Power[PR[n_],q_]*HPT3[n1_,n2_,n3_,n4_,n5_,n6_,n7_,n8_,n9_,n10_,n11_,n12_,n13_]:=Propdiff[n,q,HPT3[n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13]];
HPT3/:PR[n_]*HPT3[n1_,n2_,n3_,n4_,n5_,n6_,n7_,n8_,n9_,n10_,n11_,n12_,n13_]:=Propdiff[n,1,HPT3[n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13]];
HPT3/: delHPT3[x_]*HPT3[n1_,n2_,n3_,n4_,n5_,n6_,n7_,n8_,n9_,n10_,n11_,n12_,n13_] := Sum[-{n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13}[[i]]*PR[i]^-1 HPT3[n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13]*D[Pr[i], x], {i, 1, 13}];


applyDifferentialOperatorPart[DiffInt_,difvar_]:=Module[{apearingintegrals,v1,constfactor,i,posconst,constDiff},
apearingintegrals=Cases[{DiffInt},_HPT3,Infinity]//DeleteDuplicates;
constfactor=CoefficientRules[DiffInt,Join[apearingintegrals,{HPT3[0,0,0,0,0,0,0,0,0,0,0,0,0]}]];
posconst=Position[constfactor,Table[0,{i,1,Length[apearingintegrals]+1}]];
constDiff=0;
If[Length[posconst]>0,
	constfactor=constfactor[[posconst[[1,1]],2]];
	constDiff=D[constfactor,difvar];
];
v1=Sum[D[Coefficient[Expand[DiffInt],apearingintegrals[[i]]],difvar]apearingintegrals[[i]]+Coefficient[Expand[DiffInt],apearingintegrals[[i]]]*delHPT3[difvar]*apearingintegrals[[i]] ,{i,1,Length[apearingintegrals]}];
v1+=constDiff;
Return[v1];
];


applyDifferentialOperator[diffop_,DiffInt_]:=Module[{},
	Return[diffop/.del[x_]:>applyDifferentialOperatorPart[DiffInt,x]];
];


ComputeDIFast[toDiffInt_,difOp_,diffVarPos_,tPoint_]:=Module[{epsrepl,toDiffIntClean,diffOpPart,DeleteList,mandelstamScalarProductsT,loopMomentumScalarProductsT,DCases,DellCoeffs,delcases,delInt1,dInt,dPrefac,CoefsT,FinalD,Coefs,ResInt},
loopMomentumScalarProductsT=Expand[loopMomentumScalarProducts/.tPoint,Modulus->$P];
mandelstamScalarProductsT=Expand[Join[mandelstamScalarProducts,kinematics]/.tPoint,Modulus->$P];
DeleteList=Delete[{1,2,3,4,5,6,7,8,9,10,11},Position[tPoint[[All,1]],InternalAuxvar]];
epsrepl={(ep->(4-d)/2)/.tPoint};

(*toDiffIntClean=Collect[toDiffInt/.tPoint[[DeleteCases[DeleteList,diffVarPos+1]]]/.epsrepl,{root[___],HPT3[___]},Together,Modulus->$P];*)


toDiffIntClean=toDiffInt/.tPoint[[DeleteCases[DeleteList,diffVarPos+1]]]/.epsrepl;

diffOpPart=Collect[Coefficient[difOp,ds[invariantsnodUpperCase[[diffVarPos]]]],del[___],Together];

DCases=Cases[{toDiffIntClean},HPT3[___],Infinity]//DeleteDuplicates;
If[Length[DCases]==0,
DCases={HPT3[0,0,0,0,0,0,0,0,0,0,0,0,0]};
toDiffIntClean=toDiffIntClean HPT3[0,0,0,0,0,0,0,0,0,0,0,0,0];
];
Coefs=Collect[Coefficient[toDiffIntClean,DCases],root[___],Together];
delcases=Cases[diffOpPart,del[___],Infinity]//Sort;
DellCoeffs=Coefficient[diffOpPart,delcases];
ResInt=Sum[delInt1=Expand[delHPT3[delcases[[delvar]]/.del[x_]:>x]DCases];
dInt=Expand[Expand[DellCoeffs[[delvar]]*delInt1]/.loopMomentumScalarProductsT/.mandelstamScalarProductsT/.tPoint,Modulus->$P];
dPrefac=Expand[((DellCoeffs[[delvar]]D[Coefs,delcases[[delvar]]/.del[x_]:>x])/.tPoint)/.root'[x_]:>1/(2 root[x])//.root[x_]:>PowerMod[x,1/2,$P],Modulus->$P];
CoefsT=Expand[Coefs/.tPoint//.root[x_]:>PowerMod[x,1/2,$P],Modulus->$P];
FinalD=Expand[Total[Expand[dInt*CoefsT+dPrefac*DCases,Modulus->$P]]];
FinalD=FinalD/.n_Integer HPT3[arg___]:>Mod[n,$P]HPT3[arg];
FinalD
,{delvar,1,Length[delcases]}];
If[DCases=={HPT3[0,0,0,0,0,0,0,0,0,0,0,0,0]},
ResInt=ResInt/.HPT3[0,0,0,0,0,0,0,0,0,0,0,0,0]->1;
];
Return[ResInt/.n_Integer HPT3[arg___]:>Mod[n,$P]HPT3[arg]];
]



(* ::Title:: *)
(*Plotting*)


plotDiagram[diagram_]:=Module[{oldDir,diagramPLot,PropagatorsPlot},
	oldDir=pathHBM;
	(*SetDirectory["/home/daniel/lanayru/PhD_Files/FullHexaBoxDE/files"];*)
	SetDirectory[HomeDir<>"/files"];
	Get[FileNameJoin[{pathHBM,"KinematicGraphs.m"}]];
	diagramPLot=diagram/.{HPT3->List};
	diagramPLot=diagramPLot/.( n_Integer:>0/;n<0)/.( n_Integer:>1/;n>1);
	PropagatorsPlot=Table[Pr[i],{i,1,9}]/.mt2->0/.x_^2:>x/.p[n_]:>ToExpression[StringJoin[{"p",ToString[n]}]]/.k[n_]:>ToExpression[StringJoin[{"k",ToString[n]}]];
	PropagatorsPlot=PropagatorsPlot[[Position[diagramPLot,1][[All,1]]]];
	SetDirectory[oldDir];
	Return[buildKinematicGraph[PropagatorsPlot,{k1,k2},{p1,p2,p3,p4,p5}]];
];
