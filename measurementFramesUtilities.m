BeginPackage["measurementFramesUtilities`", {"QM`"}]



Begin["`Private`"]


opDot[op1_, op2_] := Tr@Dot[ConjugateTranspose@op1, op2];
operatorBasis[dim_] := Partition[#, dim] & /@ IdentityMatrix[dim^2];

xLikeTracelessMatrices[dim_] := Table[
    Normal@SparseArray[{{i, j} -> 1, {j, i} -> 1}, {dim, dim}],
    {i, dim}, {j, i + 1, dim}
] / Sqrt @ 2 // Flatten[#, 1] &;

yLikeTracelessMatrices[dim_] := Table[
    Normal@SparseArray[{{i, j} -> -I, {j, i} -> I}, {dim, dim}],
    {i, dim}, {j, i + 1, dim}
] / Sqrt @ 2 // Flatten[#, 1] &;

zLikeTracelessMatrices[dim_] := Table[
    Normal@SparseArray[{
            {j_, j_} :> 1 /; j <= i,
            {i + 1, i + 1} :> -i
        },
        {dim, dim}
    ] / Sqrt[i + i^2],
    {i, 1, dim - 1}
];

(* this produces an orthonormal basis of Hermitian operators (for a single qubit, it's just the Paulis) *)
identityAndTracelessOpBasis[dim_] := {
    IdentityMatrix @ dim / Sqrt @ dim,
    Sequence @@ (xLikeTracelessMatrices @ dim),
    Sequence @@ (yLikeTracelessMatrices @ dim),
    Sequence @@ (zLikeTracelessMatrices @ dim)
};


vectorizeOperator[op_] := Table[
    opDot[op, basisElement],
    {basisElement, identityAndTracelessOpBasis[Length@op]}
];

devectorizeOperator[opVec_] := Sum[
    opVec[[i]] * identityAndTracelessOpBasis[Sqrt @ Length @ opVec][[i]],
    {i, Length@opVec}
];

(* computes the canonical frame operator in the standard basis of matrices |i><j| *)


gramMatrix[listOfOperators_] := Table[
   opDot[op1, op2],
   {op1, listOfOperators},
   {op2, listOfOperators}
];


(* compute the non-rescaled canonical frame operator and represents \
it in the standard basis of matrices |i><j| *)
frameOperator[frame_, "flattenedBasis"] := Table[
    Sum[
        opDot[opI, frameElement] opDot[frameElement, opJ],
        {frameElement, frame}
    ],
    {opI, operatorBasis[Length@frame[[1]]]},
    {opJ, operatorBasis[Length@frame[[1]]]}
];

(* computes the non-rescaled canonical frame operator but represents \
in the operator basis containing identity and traceless Hermitian \
matrices *)
frameOperator[frame_] := Table[
    Sum[
        opDot[opI, frameElement] opDot[frameElement, opJ],
        {frameElement, frame}
    ],
    {opI, identityAndTracelessOpBasis[Length@frame[[1]]]},
    {opJ, identityAndTracelessOpBasis[Length@frame[[1]]]}
];


dualMeasurementFrame[operatorFrame_List] := With[{
        frameOpInv = Inverse @ frameOperator @ operatorFrame
    },
    Table[
        devectorizeOperator @ Dot[frameOpInv, vectorizeOperator @ povmOp],
        {povmOp, operatorFrame}
    ]
];

(* now compute the one rescaled with unbiased prior *)
rescaledCanonicalFrameOperator[povm_] := With[{dim = Length@First@povm},
    frameOperator[Sqrt[dim] * # / Sqrt @ Tr @ # & /@ povm]
];


canonicalEstimator[povm_] := With[{
        dim = Length@First@povm,
        frameOpInv = Inverse @ rescaledCanonicalFrameOperator @ povm
    },
    Table[
        dim * Dot[frameOpInv, Flatten @ povmOp] / Tr @ povmOp,
        {povmOp, povm}
    ]
];

nonrescaledEstimator[povm_] := With[{
        frameOpInv = Inverse @ frameOperator @ povm
    },
    Table[
        Dot[frameOpInv, Flatten@povmOp],
        {povmOp, povm}
    ]
];


sicPOVM = QStateToDensityMatrix @ # / 2 & /@ Normalize /@ {
    {1, 0},
    {1, Sqrt[2]},
    {1, Sqrt[2] Exp[2 Pi I/3]},
    {1, Sqrt[2] Exp[4 Pi I/3]}
};



End[]

EndPackage[]