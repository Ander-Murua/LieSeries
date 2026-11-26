(* ::Package:: *)

(* :Title:   LieSeries                                                            *)
(* :Author:  Ander Murua                                                          *)
(* :Version: 1.0.0                                                                *)
(* :Summary:
     Manipulation of Lie series in the free Lie algebra on two generators,
     expressed in generalized Hall bases (in the sense of Viennot).
     Includes the Baker-Campbell-Hausdorff (BCH) series, its symmetric
     version, exp(ad_A)B and the Zassenhaus terms, and a Hall basis adapted to
     the Lie algebra of classical mechanics.                                    *)
(* :References:
     A. Murua, The Hopf algebra of rooted trees, free Lie algebras, and Lie series,
       Found. Comput. Math. 6 (2006) 387-426; https://doi.org/10.1007/s10208-003-0111-0
     F. Casas, A. Murua, An efficient algorithm for computing the Baker-Campbell-Hausdorff
       series and some of its applications, J. Math. Phys. 50 (2009) 033513;
       https://doi.org/10.1063/1.3078418
     R. I. McLachlan, A. Murua, The Lie algebra of classical mechanics,
       J. Comput. Dyn. 6 (2019) 345-360; https://doi.org/10.3934/jcd.2019017;
       arXiv:1905.07554 (Hall basis "PDegree", LieAlgebra -> "ClassicalMechanics") *)
(* :History:
     2007-2008  Original version (LieSeries.m), A. Murua.
     1.0.0      Package context, basis stamp for Lie series, built-in Hall-order
                registry, Hall basis "PDegree", Lie algebra of classical mechanics,
                Lie series stored by their Hall coefficients only (the bracket
                uses the integer structure constants of the Hall basis),
                bug fixes. See CHANGELOG.md.                                     *)

BeginPackage["LieSeries`"];

Unprotect @@ Names["LieSeries`*"];
ClearAll @@ Names["LieSeries`*"];

(* ------------------------------------------------------------------------- *)
(* Public symbols                                                              *)
(* ------------------------------------------------------------------------- *)

LieSeriesInitialize::usage =
"LieSeriesInitialize[{\"X\",\"Y\"}, n] sets up the Hall basis of the free Lie algebra on the generators \"X\", \"Y\" up to degree n. \
All Lie series created afterwards are truncated at degree n. \
Calling it again creates a new basis; Lie series created under the previous basis become invalid.";

LieSeries::usage =
"LieSeries[\"X\"] gives the Lie series of the generator \"X\" (any generator name given to LieSeriesInitialize). \
LieSeries[0] gives the zero Lie series. \
Lie series can be added, multiplied by scalars, differentiated (D) and integrated (Integrate) with respect to parameters, and bracketed with LieBracket.";

LieSeriesQ::usage =
"LieSeriesQ[expr] gives True if expr is a Lie series valid in the current basis.";

LieBracket::usage =
"LieBracket[A, B] gives the Lie bracket [A, B] of two Lie series, truncated at $MaxDegree. \
LieBracket[A, B, n] keeps terms of degree up to n; LieBracket[A, B, {n}] gives only the terms of degree n. \
For square matrices of equal size, LieBracket[A, B] gives A.B - B.A. For other arguments it stays unevaluated and is displayed as [A,B].";

GenericLieSeries::usage =
"GenericLieSeries[a] gives the generic Lie series with coefficients a[i] in the Hall basis, up to $MaxDegree. \
GenericLieSeries[a, n] stops at degree n. \
With LieAlgebra -> \"ClassicalMechanics\", the coefficients a[1], a[2], ... are those of the basis of L(A,B)/I.";

HallCoefficients::usage =
"HallCoefficients[A] gives the list of coefficients of the Lie series A in the Hall basis, up to $MaxDegree, in the order of HallBasisElements. \
HallCoefficients[A, n] gives the coefficients up to degree n, and HallCoefficients[A, {n}] those of degree n. With LieAlgebra -> \"ClassicalMechanics\", it refers to the basis of L(A,B)/I given by ClassicalMechanicsBasisElements.";

NonZeroHallCoefficients::usage =
"NonZeroHallCoefficients[A] gives a list of rules HallElement -> coefficient for the nonzero coefficients of A in the Hall basis. \
NonZeroHallCoefficients[A, n] and NonZeroHallCoefficients[A, {n}] restrict to degrees up to n, or to degree n. With LieAlgebra -> \"ClassicalMechanics\", it refers to the basis of L(A,B)/I given by ClassicalMechanicsBasisElements.";

ZeroLieSeriesEquations::usage =
"ZeroLieSeriesEquations[A] gives the equations (conjunction of 0 == c) expressing that all Hall coefficients c of A vanish up to $MaxDegree. \
ZeroLieSeriesEquations[A, n] and ZeroLieSeriesEquations[A, {n}] restrict to degrees up to n, or to degree n. With LieAlgebra -> \"ClassicalMechanics\", it refers to the basis of L(A,B)/I given by ClassicalMechanicsBasisElements.";

LieSeriesSubstitution::usage =
"LieSeriesSubstitution[A, {\"X\" -> a, \"Y\" -> b}] replaces the generators of the Lie series A by a and b, computing the brackets with LieBracket (or with the function given by the option LieBracketFunction). \
a and b can be Lie series, square matrices, or any objects for which the bracket function is defined. \
LieSeriesSubstitution[A, rules, n] uses only the terms of A of degree up to n. With LieAlgebra -> \"ClassicalMechanics\", the full element of L(A,B) is substituted, so the result is exact even if a and b do not satisfy [b,[b,[b,a]]] = 0.";

LieSeriesMap::usage =
"LieSeriesMap[f, A] applies f to each coefficient of the Lie series A in the Hall basis and gives the resulting Lie series: if A = Sum h_i E_i, with h_i the Hall coefficients (as given by HallCoefficients), the result is Sum f[h_i] E_i. \
In the free Lie algebra, LieSeriesMap[f, A] coincides with f /@ A, which applies f to the internal vector of Hall coefficients. With LieAlgebra -> \"ClassicalMechanics\", f is applied to the coefficients in the basis of L(A,B)/I, and the result is the element of L(A,B) with those coefficients on that basis and zero on the other Hall elements.";

LieSeriesProj::usage =
"LieSeriesProj[A, n] gives the terms of the Lie series A of degree up to n. LieSeriesProj[A, {n}] gives the terms of degree n.";

LieSeriesSet::usage =
"LieSeriesSet[A, B, {n}] gives the Lie series that coincides with A except in degree n, where it coincides with B.";

BCH::usage =
"BCH[A, B] gives the Baker-Campbell-Hausdorff series log(exp(A) exp(B)) truncated at $MaxDegree. \
A and B can be Lie series, square matrices, or other objects (see LieSeriesSubstitution). \
BCH[A, B, n] uses only the terms of degree up to n; for Lie series the result is also truncated at degree n.";

SymmetricBCH::usage =
"SymmetricBCH[A, B] gives the symmetric BCH series log(exp(A/2) exp(B) exp(A/2)) truncated at $MaxDegree. \
SymmetricBCH[A, B, n] uses only the terms of degree up to n.";

$BCHLieSeries::usage =
"$BCHLieSeries is the Lie series log(exp(X) exp(Y)) in the generators X, Y of the current basis. It is computed on first use.";

$SymmetricBCHLieSeries::usage =
"$SymmetricBCHLieSeries is the Lie series log(exp(X/2) exp(Y) exp(X/2)) in the generators X, Y of the current basis. It is computed on first use.";

ExpAd::usage =
"ExpAd[A, B] gives exp(ad_A) B = B + [A,B] + [A,[A,B]]/2! + ... for Lie series A, B, truncated at $MaxDegree. \
ExpAd[A, B, n] keeps only the terms ad_A^k B / k! with k < n.";

ZassenhausTerms::usage =
"ZassenhausTerms[n] gives the list {C2, ..., Cn} of Zassenhaus terms, with exp(X+Y) = exp(X) exp(Y) exp(C2) exp(C3) ..., where Ck is homogeneous of degree k.";

HallBasisElements::usage =
"HallBasisElements[n] gives the Hall basis elements of degree up to n, as nested LieBracket expressions in the generator names. \
HallBasisElements[{n}] gives those of degree n. The order is the one used by HallCoefficients. With LieAlgebra -> \"ClassicalMechanics\", it refers to the basis of L(A,B)/I given by ClassicalMechanicsBasisElements.";

HallDimension::usage =
"HallDimension[n] gives the number of Hall basis elements of degree n, i.e. the dimension of the homogeneous component of degree n of the free Lie algebra. With LieAlgebra -> \"ClassicalMechanics\", it refers to the basis of L(A,B)/I given by ClassicalMechanicsBasisElements.";

ClassicalMechanicsBasisElements::usage =
"ClassicalMechanicsBasisElements[n] gives the basis of the Lie algebra of classical mechanics L(A,B)/I up to degree n, where I is the ideal generated by [B,[B,[B,A]]] and A is the first generator (kinetic energy), as the elements of the Hall basis \"PDegree\" that survive in the quotient; the other ones belong to I. \
ClassicalMechanicsBasisElements[{n}] gives those of degree n. \
It can be used with LieAlgebra -> \"ClassicalMechanics\", where it coincides with HallBasisElements, or in the free Lie algebra with HallOrder -> \"PDegree\". \
Reference: R. I. McLachlan, A. Murua, The Lie algebra of classical mechanics, J. Comput. Dyn. 6 (2019) 345-360, https://doi.org/10.3934/jcd.2019017.";

ClassicalMechanicsCoefficients::usage =
"ClassicalMechanicsCoefficients[A] gives the coefficients of the Lie series A in the basis of the Lie algebra of classical mechanics given by ClassicalMechanicsBasisElements, i.e. the coefficients of the image of A in L(A,B)/I. \
ClassicalMechanicsCoefficients[A, n] and ClassicalMechanicsCoefficients[A, {n}] restrict to degrees up to n, or to degree n. \
It can be used with LieAlgebra -> \"ClassicalMechanics\", where it coincides with HallCoefficients, or in the free Lie algebra with HallOrder -> \"PDegree\".";

HallOrders::usage =
"HallOrders[] gives an association of the built-in Hall orders that can be chosen with the option HallOrder of LieSeriesInitialize, with a short description of each.";

LieSeriesBasisInfo::usage =
"LieSeriesBasisInfo[] gives an association describing the current basis: generators, maximal degree, Lie algebra (\"Free\" or \"ClassicalMechanics\"), Hall order, bracket convention and dimensions.";

$MaxDegree::usage =
"$MaxDegree is the maximal degree fixed by the last call to LieSeriesInitialize.";

$FormatMaxDegree::usage =
"$FormatMaxDegree is the maximal degree of the terms shown when a Lie series is displayed. It can be reset by the user; terms of higher degree are kept but indicated by \[Ellipsis].";

(* Options *)

HallOrder::usage =
"HallOrder is an option of LieSeriesInitialize that selects the generalized Hall basis. Possible values are the keys of HallOrders[]. \
The default, Automatic, gives \"Standard\" with LieAlgebra -> \"Free\" and \"PDegree\" with LieAlgebra -> \"ClassicalMechanics\".";

LieAlgebra::usage =
"LieAlgebra is an option of LieSeriesInitialize. LieAlgebra -> \"Free\" (default) works in the free Lie algebra L(X,Y). \
LieAlgebra -> \"ClassicalMechanics\" works in the Lie algebra of classical mechanics L(A,B)/I, where I is the ideal generated by [B,[B,[B,A]]], the first generator A being the kinetic energy and the second B the potential energy (McLachlan & Murua 2019). \
Lie series are still computed as elements of L(A,B) in the Hall basis \"PDegree\", so substitutions (LieSeriesSubstitution, BCH with matrices, ...) are always exact, but HallBasisElements, HallCoefficients, NonZeroHallCoefficients, ZeroLieSeriesEquations, HallDimension, GenericLieSeries, LieSeriesMap and the display of Lie series refer to L(A,B)/I.";

Reversed::usage =
"Reversed is an option of LieSeriesInitialize. With Reversed -> True (default) a tree with left subtree t1 and right subtree t2 represents the bracket [t2, t1]; with Reversed -> False it represents [t1, t2].";

CheckLessFunction::usage =
"CheckLessFunction is an option of LieSeriesInitialize. With CheckLessFunction -> True the Hall order is checked to be a strict total order (asymmetric, total and transitive) on the Hall trees. The condition that each Hall tree is greater than its right factor is always checked.";

FormatMaxDegree::usage =
"FormatMaxDegree is an option of LieSeriesInitialize giving the initial value of $FormatMaxDegree.";

LieBracketFunction::usage =
"LieBracketFunction is an option of LieSeriesSubstitution, BCH and SymmetricBCH (and of LieSeriesInitialize, as default value) giving the function used to bracket the substituted objects.";

SparseArrays::usage =
"SparseArrays is an option of LieSeriesInitialize. With SparseArrays -> True the generators and the zero Lie series are stored as sparse arrays.";

BCHLieSeries::usage =
"BCHLieSeries is an option of LieSeriesInitialize. With BCHLieSeries -> True, $BCHLieSeries is computed during initialization.";

SymmetricBCHLieSeries::usage =
"SymmetricBCHLieSeries is an option of LieSeriesInitialize. With SymmetricBCHLieSeries -> True, $SymmetricBCHLieSeries is computed during initialization.";

(* Messages *)

LieSeries::noinit = "No basis has been set up yet. Call LieSeriesInitialize first.";
LieSeries::stale  = "A Lie series created before the last call to LieSeriesInitialize was used. It is not valid in the current basis.";
LieSeries::cmbasis = "This function requires LieAlgebra -> \"ClassicalMechanics\" or HallOrder -> \"PDegree\" in LieSeriesInitialize.";
LieSeries::deg    = "The degree `1` must be between 1 and $MaxDegree = `2`.";
LieSeriesInitialize::args  = "LieSeriesInitialize must be called as LieSeriesInitialize[{\"X\",\"Y\"}, n] with two distinct generator names and a positive integer n.";
LieSeriesInitialize::order = "Unknown Hall order `1`. The available orders are `2`.";
LieSeriesInitialize::lieal = "Unknown value `1` of LieAlgebra. The possible values are \"Free\" and \"ClassicalMechanics\".";
LieSeriesInitialize::cmorder = "LieAlgebra -> \"ClassicalMechanics\" requires HallOrder -> \"PDegree\" (or Automatic), not `1`.";
LieSeriesInitialize::sparse = "The value of the option SparseArrays must be True or False, not `1`.";
LieSeriesInitialize::order1 = "The Hall order `1` gives the same value for both comparisons of trees {j,k} = `2`.";
LieSeriesInitialize::order2 = "The Hall order `1` does not give a Boolean value for the comparison of trees {j,k} = `2`.";
LieSeriesInitialize::order3 = "The Hall order `1` does not satisfy the Hall condition (right factor smaller than the tree) for the tree `2`.";
LieSeriesInitialize::order4 = "The Hall order `1` is not transitive on the Hall trees: after sorting them with it, the trees `2` are not in increasing order.";
LieSeriesSubstitution::gens = "The substitution rules must give values for the generators `1`.";
ZassenhausTerms::range = "The argument must be an integer between 2 and $MaxDegree = `1`.";


Begin["`Private`"];

ClearAll["LieSeries`Private`*"];

(* ------------------------------------------------------------------------- *)
(* Registry of built-in Hall orders                                            *)
(*                                                                             *)
(* To add a new built-in basis: write an order function less[j, k] on tree     *)
(* indices (it may use l, r, degree, treeWord, tallDegree below) and add an    *)
(* entry to $hallOrders. Nothing else needs to change.                         *)
(*                                                                             *)
(* less determines the decomposition t = l[t] o r[t] of each tree with more    *)
(* than one vertex: r[t] is the greatest subtree attached to the root, grafted *)
(* onto the root of l[t]. The Hall trees are 1 (X), 2 (Y) and, recursively,    *)
(* the trees t with l[t], r[t] Hall trees and less[r[t], l[t]] (for Hall l[t], *)
(* these are the pairs (j, k) with less[k, j] and not less[k, r[j]]). They     *)
(* form a generalized Hall set (Hall-Viennot basis) provided that              *)
(*   (R1) less is a strict total order on all trees, and                       *)
(*   (R2) less[r[t], t] for every Hall tree t of degree > 1, i.e. each Hall    *)
(*        tree is greater than all its proper subtrees.                        *)
(* (R1) is needed on all trees, not only on Hall trees, because less also      *)
(* gives the normal form of the auxiliary trees (see dot in buildTables). A    *)
(* simple sufficient condition for (R2) is less[r[t], t] for every tree t,     *)
(* which holds e.g. if less compares the degree first. (R2) is always checked, *)
(* right after the Hall trees are built; (R1) is checked on the Hall trees     *)
(* with CheckLessFunction -> True.                                             *)
(*                                                                             *)
(* References: Reutenauer 1993, Ch. 4; Murua 2006, Remark 17; Casas-Murua      *)
(* 2009, Sec. II.C.                                                            *)
(* ------------------------------------------------------------------------- *)

indexLess[j_Integer, k_Integer] := j < k;

(* equal words (only possible for non-Hall trees) are ordered by index, so that the
   order is total on all trees *)
lyndonLess[j_Integer, k_Integer] :=
  With[{wj = treeWord[j], wk = treeWord[k]}, If[wj === wk, j < k, OrderedQ[{wk, wj}]]];

tallLess[j_Integer, k_Integer] :=
  tallDegree[j] < tallDegree[k] || (tallDegree[j] == tallDegree[k] && j < k);

leftLess[j_Integer, k_Integer] :=
  degree[j] < degree[k] || (degree[j] == degree[k] && (l[j] < l[k] || (l[j] == l[k] && j < k)));

rightLess[j_Integer, k_Integer] :=
  degree[j] < degree[k] || (degree[j] == degree[k] && (r[j] < r[k] || (r[j] == r[k] && j < k)));

(* P-degree (McLachlan & Murua, 2019): polynomial degree in the momenta when the first
   generator is the kinetic energy (degree 2) and the second the potential energy (degree 0);
   a bracket lowers the degree by one, and brackets of two elements of degree 0 vanish. *)
pdegreeLess[j_Integer, k_Integer] :=
  pdegree[k] > pdegree[j] || (pdegree[k] == pdegree[j] && k > j);

$hallOrders = <|
  "Standard" -> <|"Less" -> indexLess,
     "Description" -> "Trees compared by their construction index: first by degree, then in order of generation."|>,
  "Lyndon"   -> <|"Less" -> lyndonLess,
     "Description" -> "Lyndon basis: trees compared through the lexicographic order of their foliage words."|>,
  "Tall"     -> <|"Less" -> tallLess,
     "Description" -> "Trees compared by tall degree t (t = 1 for generators, t = max(t(left), t(right)+1)), then by index."|>,
  "Left"     -> <|"Less" -> leftLess,
     "Description" -> "Trees compared by degree, then by the index of the left subtree, then by index."|>,
  "Right"    -> <|"Less" -> rightLess,
     "Description" -> "Trees compared by degree, then by the index of the right subtree, then by index."|>,
  "PDegree" -> <|"Less" -> pdegreeLess,
     "Description" -> "Trees compared by P-degree (2 for the first generator, 0 for the second, max(p1 + p2 - 1, 0) for a bracket), then by index (McLachlan & Murua 2019). A subset of this basis is a basis of the Lie algebra of classical mechanics; see LieAlgebra and ClassicalMechanicsBasisElements."|>
|>;

HallOrders[] := #["Description"] & /@ $hallOrders;

(* ------------------------------------------------------------------------- *)
(* Global state of the current basis (private)                                 *)
(* ------------------------------------------------------------------------- *)

$initialized = False;
$id = 0;              (* stamp of the current basis; every Lie series carries it *)
$gens = {};
$reversed = True;
$orderName = None;
$lieAlgebra = "Free";
$quotient = False;   (* True with LieAlgebra -> "ClassicalMechanics" *)
$bracketDefault = LieBracket;
$bchCache = None;
$sbchCache = None;

(* Lie series valid in the current basis carry the current stamp *)
cur[id_] := id === $id;
staleIdQ[id_] := id =!= $id;
staleQ[args_Hold] := ! FreeQ[args, LieSeries[id_Integer, _] /; id =!= $id];

degOK[n_Integer] := If[1 <= n <= $MaxDegree, True, Message[LieSeries::deg, n, $MaxDegree]; False];
degOK[_] := False;

(* ------------------------------------------------------------------------- *)
(* Initialization                                                              *)
(* ------------------------------------------------------------------------- *)

Options[LieSeriesInitialize] = {
  HallOrder -> Automatic, LieAlgebra -> "Free", Reversed -> True,
  CheckLessFunction -> False, FormatMaxDegree -> Automatic,
  LieBracketFunction -> LieBracket, SparseArrays -> False,
  BCHLieSeries -> False, SymmetricBCHLieSeries -> False};

(* unknown options are reported before any OptionValue is evaluated (OptionValue would
   issue OptionValue::nodef for them) *)
unknownOptions[opts_List] := FilterRules[opts, Except[Options[LieSeriesInitialize]]];

LieSeriesInitialize[{x_String, y_String}, maxdeg_Integer, opts : OptionsPattern[]] /;
    x =!= y && maxdeg >= 1 && unknownOptions[{opts}] =!= {} :=
  (Message[LieSeriesInitialize::optx, First[unknownOptions[{opts}]], HoldForm[LieSeriesInitialize]];
   $Failed);

LieSeriesInitialize[{x_String, y_String}, maxdeg_Integer, opts : OptionsPattern[]] /;
    x =!= y && maxdeg >= 1 && unknownOptions[{opts}] === {} :=
  Module[{orderName, sparse = OptionValue[SparseArrays], lieal = OptionValue[LieAlgebra]},
    (* choice of the Hall order *)
    orderName = Which[
      OptionValue[HallOrder] =!= Automatic, OptionValue[HallOrder],
      lieal === "ClassicalMechanics", "PDegree",
      True, "Standard"];
    Which[
      ! MemberQ[{"Free", "ClassicalMechanics"}, lieal],
        Message[LieSeriesInitialize::lieal, lieal]; $Failed,
      ! MemberQ[Keys[$hallOrders], orderName],
        Message[LieSeriesInitialize::order, orderName, Keys[$hallOrders]]; $Failed,
      lieal === "ClassicalMechanics" && orderName =!= "PDegree",
        Message[LieSeriesInitialize::cmorder, orderName]; $Failed,
      ! MemberQ[{True, False}, sparse],
        Message[LieSeriesInitialize::sparse, sparse]; $Failed,
      True,
        $lieAlgebra = lieal; $quotient = (lieal === "ClassicalMechanics");
        initialize[x, y, maxdeg, orderName, sparse,
          OptionValue[Reversed],
          OptionValue[LieBracketFunction], OptionValue[FormatMaxDegree],
          OptionValue[CheckLessFunction],
          OptionValue[BCHLieSeries], OptionValue[SymmetricBCHLieSeries]]]];

LieSeriesInitialize[___] := (Message[LieSeriesInitialize::args]; $Failed);

initialize[x_, y_, maxdeg_, orderName_, sparse_, reversed_, bracket_, fmt_,
    check_, bch_, sbch_] := (
    (* new basis *)
    $id = $id + 1;
    $gens = {x, y};
    $orderName = orderName;
    $reversed = TrueQ[reversed];
    $bracketDefault = bracket;
    $sparse = sparse;
    $bchCache = None; $sbchCache = None;
    Unprotect[$MaxDegree]; $MaxDegree = maxdeg; Protect[$MaxDegree];
    $FormatMaxDegree = If[IntegerQ[fmt] && fmt >= 1, fmt, maxdeg];

    buildTables[maxdeg, $hallOrders[orderName]["Less"], orderName, TrueQ[check]];
    $initialized = True;

    If[TrueQ[bch], bchSeries[]];
    If[TrueQ[sbch], symmetricBCHSeries[]];
  );

(* Builds all the tables of the basis. Trees are numbered by positive integers.
   Trees 1 and 2 are the generators. A tree i of degree > 1 has left subtree l[i]
   and right subtree r[i] (r = 0 for generators). Trees 1 .. $nHall are the Hall
   trees up to degree maxdeg, sorted by degree: f[n] .. f[n+1]-1 have degree n.
   During initialization, trees $nHall+1 .. nTrees are auxiliary (non-Hall) trees,
   used only to compute the structure constants of the bracket; they are discarded
   afterwards.

   Representation of Lie series. With E_i the Hall element i as a nested bracket, a Lie
   series A = Sum h_i E_i is stored as the vector (h_1, ..., h_$nHall) of its Hall
   coefficients. (The value of A at the Hall tree i in the Lie coalgebra of trees is
   sigma_i h_i, with sigma_i the symmetry factor of the tree; trees are used only during
   initialization.) See the references in the header. *)
(* Stage 1: the Hall trees up to degree maxdeg for the order less, with their P-degree and
   survival in the Lie algebra of classical mechanics. Defines l, r, degree, f, hallQ,
   treeWord, tallDegree, pdegree, cmSurvivesQ, $treeCounter, $nHall, and dot on Hall trees. *)
hallTrees[maxdeg_, less_, orderName_, check_] :=
  Block[{i, j, k, m, n, lag1, lag2},
    Clear[l, r, degree, f, hallQ, treeWord, tallDegree, pdegree, cmSurvivesQ, dot];

    l[1] = 1; l[2] = 2; r[1] = 0; r[2] = 0; degree[1] = 1; degree[2] = 1;
    pdegree[1] = 2; pdegree[2] = 0;
    pdegree[i_Integer] := pdegree[i] = Max[pdegree[l[i]] + pdegree[r[i]] - 1, 0];
    (* Hall trees that survive in the Lie algebra of classical mechanics: no subtree is a
       bracket of two elements of P-degree 0 *)
    cmSurvivesQ[i_Integer] := cmSurvivesQ[i] = r[i] == 0 ||
       (pdegree[l[i]] + pdegree[r[i]] > 0 && cmSurvivesQ[l[i]] && cmSurvivesQ[r[i]]);
    f[1] = 1;
    hallQ[1] = True; hallQ[2] = True;
    treeWord[i_Integer] := treeWord[i] =
      If[degree[i] == 1, ToString[i], treeWord[l[i]] <> treeWord[r[i]]];
    tallDegree[i_Integer] := tallDegree[i] =
      If[degree[i] == 1, 1, Max[tallDegree[l[i]], tallDegree[r[i]] + 1]];

    (* Hall trees, degree by degree *)
    $treeCounter = 3;
    Do[
      f[n] = $treeCounter;
      Do[
        m = n - degree[k];
        Do[
          If[(r[j] == 0 || ! less[k, r[j]]) && less[k, j] && hallQ[j] && hallQ[k],
            i = $treeCounter;
            l[i] = j; r[i] = k; degree[i] = n; dot[j, k] = i;
            hallQ[i] = True;
            $treeCounter++],
          {j, f[m], f[m + 1] - 1}],
        {k, 1, f[n] - 1}],
      {n, 2, maxdeg}];
    f[maxdeg + 1] = $treeCounter;
    $nHall = $treeCounter - 1;
    hallQ[i_Integer /; i > $nHall] = False;

    (* condition (R2), always checked: every Hall tree i of degree > 1 is greater than its
       right factor r[i] *)
    Do[If[! TrueQ[less[r[i], i]], Message[LieSeriesInitialize::order3, orderName, i]],
      {i, 3, $nHall}];

    (* condition (R1) on the Hall trees, with CheckLessFunction -> True: less is a strict total
       order. First, for every pair exactly one of less[j, k], less[k, j] is True; then, after
       sorting the Hall trees with less, every pair of the sorted list is in increasing order
       (so that less is the order of that list, hence transitive). *)
    If[check,
      Do[
        lag1 = less[j, k]; lag2 = less[k, j];
        If[(lag1 && ! lag2) || (! lag1 && lag2), Null,
          Message[LieSeriesInitialize::order1, orderName, {j, k}],
          Message[LieSeriesInitialize::order2, orderName, {j, k}]],
        {j, 2, $nHall}, {k, j - 1}];
      Module[{sorted = Sort[Range[$nHall], less[#1, #2] &], a, b},
        Catch[
          Do[
            If[! TrueQ[less[sorted[[a]], sorted[[b]]]],
              Message[LieSeriesInitialize::order4, orderName, {sorted[[a]], sorted[[b]]}];
              Throw[Null]],
            {a, 1, $nHall - 1}, {b, a + 1, $nHall}]]]];
  ];

(* Number of Hall elements and of elements surviving in the Lie algebra of classical
   mechanics, by degree up to maxdeg, for the basis "PDegree". Builds only the
   Hall trees (no bracket structure), so it is fast enough for high degrees, and it does not
   modify the current basis. Used by the tests. *)
cmDimensions[maxdeg_Integer] :=
  Block[{l, r, degree, f, hallQ, treeWord, tallDegree, pdegree, cmSurvivesQ, dot,
      $treeCounter, $nHall},
    hallTrees[maxdeg, pdegreeLess, "PDegree", False];
    {Table[f[n + 1] - f[n], {n, maxdeg}],
     Table[Count[Range[f[n], f[n + 1] - 1], _?cmSurvivesQ], {n, maxdeg}]}];

(* Brackets in tree coordinates (used in step 2 of buildTables). Tree coordinates are
   Associations tree -> value. cutsByTrunk[a] is the list of {t, b} for the cuts (a, b) of
   the trees t considered (trunk a, pruned branch b, with multiplicity).

   cutProducts[x, y, cutsByTrunk] is the list of pairs {t, x(a) y(b)}, one for each cut (a, b)
   of each tree t; only trunks a in the support of x are visited. *)
cutProducts[x_, y_, cutsByTrunk_] :=
  Join @@ KeyValueMap[
    Function[{a, xa}, Map[{#[[1]], xa Lookup[y, #[[2]], 0]} &, Lookup[cutsByTrunk, a, {}]]],
    x];

(* treeBracket[x, y, cutsByTrunk] gives, for the trees t whose cuts are in cutsByTrunk,
     t -> Sum_{(a,b) in S(t)} (x(a) y(b) - y(a) x(b)),
   keeping only nonzero values (-y negates all the values of y). *)
treeBracket[x_, y_, cutsByTrunk_] :=
  Select[
    GroupBy[
      Join[cutProducts[x, y, cutsByTrunk], cutProducts[-y, x, cutsByTrunk]],
      First -> Last, Total],
    # =!= 0 &];

(* Stage 2 (after stage 1): auxiliary trees, structure constants of the bracket in the Hall
   basis, and removal of the auxiliary trees.

   E_i denotes the Hall element i as a nested commutator: E_i = [E_l[i], E_r[i]], or
   E_i = [E_r[i], E_l[i]] if $reversed. Lie series are stored as vectors of coefficients in
   the basis E_i (the Hall coefficients).

   A single-edge cut of a tree t is a pair (a, b): the trunk a, which contains the root of t,
   and the pruned branch b. bp[t] lists the cuts of t as triples {t, a, b}, with multiplicity.
   For a Hall tree t every pruned branch b is a Hall tree, but for an auxiliary tree b need
   not be. A Lie series A is also determined by its values A(t) at the trees (its tree
   coordinates): its value at the Hall tree h is sigma_h times its coefficient of E_h, and
       [A, B](t) = eps Sum_{(a,b) in bp[t]} (A(a) B(b) - A(b) B(a)),
   with eps = -1 if $reversed and eps = +1 otherwise.
   e[u] (an Association tree -> integer) gives the tree coordinates of E_u at the trees of
   degree |u|. For both values of eps,
       e[u](t) = Sum_{(a,b) in bp[t]} (e[l[u]](a) e[r[u]](b) - e[l[u]](b) e[r[u]](a)).
   The matrix mE, with mE[[a, u]] = e[u](a), is built once all e[u] are known. Since
   e[v](h) = sigma_v delta(h, v) at the Hall trees h, the structure constants
       [E_u, E_v] = Sum_s C[s, u, v] E_s,
       C[s, u, v] = (eps/sigma_s) Sum_{(a,b) in bp[s]}
                      (e[u](a) sigma_v delta(b, v) - sigma_u delta(b, u) e[v](a)),
   are integers and antisymmetric in (u, v); only the entries with u < v are stored, in
   bracketRules[n] and bracketTensor[n] (n = |s|). The degree-n part of [A, B] then has s-th
   coefficient Sum_{u<v} C[s, u, v] (A_u B_v - A_v B_u).
   Only the stage-1 tables (Hall trees) are kept. *)
buildTables[maxdeg_, less_, orderName_, check_] :=
  Block[{i, n, u, eps, nTrees, nCuts, treeDegrees, cutsByTrunk, e,
      mE, sigmaVector, epsSigmaInvD, cnt, ts, q, cn, t0, t1, t2, t3, tt, tCoords = 0, tConsts = 0},
    ClearAll[degree1, kappa, sigma, degreeMask, hallSet, visibleSet, bracketRules,
      bracketTensor, hallExpr, bp];
    t0 = AbsoluteTime[];
    hallTrees[maxdeg, less, orderName, check];
    t1 = AbsoluteTime[];

    (* dot[j, k]: index of the tree obtained by grafting k onto the root of j, creating
       auxiliary trees when needed. Normal form: tree i = l[i] o r[i], with r[i] the
       greatest (for less) of the subtrees attached to the root of i, so that every tree
       has a single index. Hall trees (stage 1) are already in this normal form. *)
    dot[j_Integer, k_Integer] :=
      (dot[j, k] = dot[dot[l[j], k], r[j]]) /; r[j] =!= 0 && less[k, r[j]];
    dot[j_Integer, k_Integer] :=
      Module[{ii = $treeCounter},
        dot[j, k] = ii; l[ii] = j; r[ii] = k; degree[ii] = degree[j] + degree[k];
        $treeCounter++;
        ii] /; r[j] === 0 || ! less[k, r[j]];

    (* bp[i]: the single-edge cuts {i, a, b} of the tree i. Computing bp may create new
       auxiliary trees, so we loop until no new tree appears. *)
    bp[1] = {}; bp[2] = {};
    bp[i_Integer] := bp[i] = Join[
       {{i, l[i], r[i]}},
       Map[{i, dot[#[[2]], r[i]], #[[3]]} &, bp[l[i]]],
       Map[{i, dot[l[i], #[[2]]], #[[3]]} &, bp[r[i]]]];
    i = 3;
    While[i < $treeCounter, bp[i]; i++];
    nTrees = $treeCounter - 1;
    treeDegrees = Table[degree[i], {i, nTrees}];
    nCuts = Total[Map[Length[bp[#]] &, Range[nTrees]]];
    t2 = AbsoluteTime[];

    (* kappa and symmetry factors sigma of Hall trees (used only below) *)
    kappa[1] = 1; kappa[2] = 1; sigma[1] = 1; sigma[2] = 1;
    kappa[i_Integer] := kappa[i] = If[r[l[i]] == r[i], kappa[l[i]] + 1, 1];
    sigma[i_Integer] := sigma[i] = kappa[i] sigma[l[i]] sigma[r[i]];

    (* tree coordinates e[u] of the Hall elements (step 2), degree by degree *)
    eps = If[$reversed, -1, 1];
    e[1] = <|1 -> 1|>; e[2] = <|2 -> 1|>;
    Do[
      tt = AbsoluteTime[];
      (* cuts {t, b} of the trees t of degree n, indexed by their trunk a *)
      cutsByTrunk = GroupBy[Join @@ Map[bp, Pick[Range[nTrees], treeDegrees, n]],
        (#[[2]] &) -> (#[[{1, 3}]] &)];
      (* E_u = [E_l[u], E_r[u]] (or [E_r[u], E_l[u]] if $reversed) in tree coordinates, for
         the Hall trees u of degree n; the formula of treeBracket is the same for both values
         of eps *)
      Do[e[u] = treeBracket[e[l[u]], e[r[u]], cutsByTrunk], {u, f[n], f[n + 1] - 1}];
      tCoords += AbsoluteTime[] - tt,
      {n, 2, maxdeg}];

    (* structure constants (step 3), with sparse arrays. mE[[a, w]] = e[w](a): the columns
       of mE are the tree coordinates of the Hall elements. For the Hall trees s of degree n,
       ts[[s, b, a]] = sigma_b (multiplicity of the cut (a, b) of s); the branch b is a Hall
       tree. Then
         q[[s, v, u]] = (ts . mE)[[s, v, u]] = sigma_v Sum_a (multiplicity of (a, v)) e[u](a),
         C[[s, u, v]] = (eps/sigma_s) (q[[s, v, u]] - q[[s, u, v]]).
       Transpose[x, {1, 3, 2}] exchanges the last two indices of a three-index array. *)
    tt = AbsoluteTime[];
    mE = SparseArray[
      Join @@ Table[KeyValueMap[({#1, w} -> #2) &, e[w]], {w, $nHall}], {nTrees, $nHall}];
    sigmaVector = Table[sigma[i], {i, $nHall}];
    epsSigmaInvD = DiagonalMatrix[SparseArray[eps/sigmaVector]];
    bracketRules[1] = {};
    Do[
      (* cuts {s, a, b} of the Hall trees s of degree n, stored as {s, b, a} with value
         sigma_b (multiplicity); built from the whole lists of keys and counts at once *)
      cnt = Counts[Join @@ Map[bp, Range[f[n], f[n + 1] - 1]]];
      ts = SparseArray[
        Keys[cnt][[All, {1, 3, 2}]] -> Values[cnt] sigmaVector[[Keys[cnt][[All, 3]]]],
        {$nHall, $nHall, nTrees}];
      q = ts . mE;
      cn = epsSigmaInvD . (Transpose[q, {1, 3, 2}] - q);
      bracketRules[n] = Select[Most[ArrayRules[cn]], #[[1, 2]] < #[[1, 3]] && #[[2]] =!= 0 &],
      {n, 2, maxdeg}];
    tConsts = AbsoluteTime[] - tt;
    bracketTensor[1] = SparseArray[{}, {$nHall, $nHall, $nHall}];
    Do[
      bracketTensor[n] = If[bracketRules[n] === {},
        SparseArray[{}, {$nHall, $nHall, $nHall}],
        SparseArray[bracketRules[n], {$nHall, $nHall, $nHall}]],
      {n, 2, maxdeg}];

    t3 = AbsoluteTime[];

    (* discard the auxiliary trees *)
    ClearAll[bp, dot];
    Do[With[{j = i},
        l[j] =.; r[j] =.; degree[j] =.;
        (* values cached by the order functions for auxiliary trees *)
        Quiet[treeWord[j] =.; tallDegree[j] =.; pdegree[j] =.]],
      {i, $nHall + 1, nTrees}];
    $treeCounter = $nHall + 1;

    (* tables on the Hall trees *)
    $degreeVector = Table[degree[i], {i, $nHall}];
    degreeMask[{n_Integer}] := degreeMask[{n}] = Map[If[# == n, 1, 0] &, $degreeVector];
    degreeMask[n_Integer] := degreeMask[n] = Map[If[# <= n, 1, 0] &, $degreeVector];
    hallSet[{n_Integer}] := hallSet[{n}] = Range[f[n], f[n + 1] - 1];
    hallSet[n_Integer] := hallSet[n] = Range[1, f[n + 1] - 1];
    (* Hall elements visible to the user: all of them in L(X,Y), the surviving ones in
       L(A,B)/I with LieAlgebra -> "ClassicalMechanics" *)
    visibleSet[spec_] := visibleSet[spec] =
      If[$quotient, Select[hallSet[spec], cmSurvivesQ], hallSet[spec]];
    $visibleMask = Table[If[! $quotient || cmSurvivesQ[i], 1, 0], {i, $nHall}];

    (* degree in the first generator *)
    degree1[1] = 1; degree1[2] = 0;
    degree1[i_Integer] := degree1[i] = degree1[l[i]] + degree1[r[i]];
    $degree1Vector = Table[degree1[i], {i, $nHall}];

    (* generators and zero as coefficient vectors *)
    $genVectors = If[$sparse,
      {SparseArray[{{1} -> 1}, {$nHall}], SparseArray[{{2} -> 1}, {$nHall}]},
      {Normal[SparseArray[{{1} -> 1}, {$nHall}]], Normal[SparseArray[{{2} -> 1}, {$nHall}]]}];
    $zeroVector = If[$sparse, SparseArray[{}, {$nHall}], ConstantArray[0, $nHall]];

    (* Hall basis elements as nested LieBracket expressions *)
    hallExpr[1] = $gens[[1]]; hallExpr[2] = $gens[[2]];
    hallExpr[i_Integer] := hallExpr[i] =
      If[$reversed,
        LieBracket[hallExpr[r[i]], hallExpr[l[i]]],
        LieBracket[hallExpr[l[i]], hallExpr[r[i]]]];

    (* time (in seconds) of each phase of the last initialization, and sizes; read by
       Benchmarks/LieSeriesTimings.wl *)
    $buildTimings = <|
      "HallTrees" -> t1 - t0, "AuxiliaryTreesAndCuts" -> t2 - t1,
      "TreeCoordinates" -> tCoords, "StructureConstants" -> tConsts,
      "TensorsAndCleanup" -> AbsoluteTime[] - t3 + (t3 - t2 - tCoords - tConsts),
      "Total" -> AbsoluteTime[] - t0,
      "HallElements" -> $nHall, "Trees" -> nTrees, "Cuts" -> nCuts,
      "StructureConstantEntries" -> Total[Table[Length[bracketRules[k]], {k, 2, maxdeg}]]|>;
  ];

(* ------------------------------------------------------------------------- *)
(* Constructors and predicates                                                 *)
(* ------------------------------------------------------------------------- *)

LieSeries[name_String] /; ! $initialized := (Message[LieSeries::noinit]; $Failed);
LieSeries[name_String] /; MemberQ[$gens, name] :=
  LieSeries[$id, $genVectors[[First[Flatten[Position[$gens, name]]]]]];
LieSeries[0] /; ! $initialized := (Message[LieSeries::noinit]; $Failed);
LieSeries[0] := LieSeries[$id, $zeroVector];

LieSeriesQ[s_] := MatchQ[s, LieSeries[_?cur, c_ /; Length[c] === $nHall]];

(* ------------------------------------------------------------------------- *)
(* Vector-space structure and operations acting on the coefficients            *)
(* ------------------------------------------------------------------------- *)

notSeriesQ[a_] := FreeQ[a, LieSeries];

LieSeries /: (a_?notSeriesQ) LieSeries[id_?cur, c_] := LieSeries[id, a c];
LieSeries /: LieSeries[id_?cur, c1_] + LieSeries[id_, c2_] := LieSeries[id, c1 + c2];

LieSeries /: D[LieSeries[id_?cur, c_], args__] := LieSeries[id, D[c, args]];
LieSeries /: Integrate[LieSeries[id_?cur, c_], args__] := LieSeries[id, Integrate[c, args]];
LieSeries /: Expand[LieSeries[id_?cur, c_], args___] := LieSeries[id, Expand[c, args]];
LieSeries /: Collect[LieSeries[id_?cur, c_], args__] := LieSeries[id, Collect[c, args]];
LieSeries /: Coefficient[LieSeries[id_?cur, c_], args__] := LieSeries[id, Coefficient[c, args]];
LieSeries /: Simplify[LieSeries[id_?cur, c_], args___] := LieSeries[id, Simplify[c, args]];
LieSeries /: FullSimplify[LieSeries[id_?cur, c_], args___] := LieSeries[id, FullSimplify[c, args]];
LieSeries /: Together[LieSeries[id_?cur, c_], args___] := LieSeries[id, Together[c, args]];
LieSeries /: Factor[LieSeries[id_?cur, c_], args___] := LieSeries[id, Factor[c, args]];
LieSeries /: N[LieSeries[id_?cur, c_], args___] := LieSeries[id, N[c, args]];
LieSeries /: Normal[LieSeries[id_?cur, c_]] := LieSeries[id, Normal[c]];
LieSeries /: SparseArray[LieSeries[id_?cur, c_]] := LieSeries[id, SparseArray[c]];
(* f /@ A applies f to the coefficient vector, as in earlier versions *)
LieSeries /: Map[fn_, LieSeries[id_?cur, c_]] := LieSeries[id, fn[c]];

(* ------------------------------------------------------------------------- *)
(* Lie bracket                                                                 *)
(* ------------------------------------------------------------------------- *)

LieBracket[LieSeries[id_?cur, c1_], LieSeries[id_, c2_], {n_Integer?degOK}] :=
  LieSeries[id, (bracketTensor[n] . c2) . c1 - (bracketTensor[n] . c1) . c2];

LieBracket[a : LieSeries[_?cur, _], b : LieSeries[_?cur, _], n_Integer?degOK] :=
  If[n < 2, LieSeries[0], Plus @@ Table[LieBracket[a, b, {m}], {m, 2, n}]];

LieBracket[a : LieSeries[_?cur, _], b : LieSeries[_?cur, _]] := LieBracket[a, b, $MaxDegree];

LieBracket[a_?MatrixQ, b_?MatrixQ] /;
    Dimensions[a] === Dimensions[b] === Reverse[Dimensions[b]] := a . b - b . a;

(* infix notation A \[Diamond] B for the bracket of two Lie series *)
LieSeries /: Diamond[a : LieSeries[_?cur, _], b : LieSeries[_?cur, _]] := LieBracket[a, b];

Format[LieBracket[a_, b_]] := Row[{"[", a, ",", b, "]"}];

(* ------------------------------------------------------------------------- *)
(* Coefficients and generic series                                             *)
(* ------------------------------------------------------------------------- *)

GenericLieSeries[coef_Symbol] := GenericLieSeries[coef, $MaxDegree];
GenericLieSeries[coef_Symbol, n_Integer?degOK] /; $initialized :=
  LieSeries[$id, hallVector[Map[coef, Range[Length[visibleSet[n]]]], visibleSet[n]]];

(* coefficient vector of the Lie series Sum_k h[[k]] E_idx[[k]] (zero on the other Hall
   elements) *)
hallVector[h_List, idx_List] :=
  Module[{v = ConstantArray[0, $nHall]}, v[[idx]] = h; v];

HallCoefficients[s : LieSeries[_?cur, _]] := HallCoefficients[s, $MaxDegree];
HallCoefficients[LieSeries[_?cur, c_], n : (_Integer?degOK | {_Integer?degOK})] :=
  hallCoefficients[c, visibleSet[n]];

(* coefficients of the coefficient vector c on the Hall trees idx *)
hallCoefficients[c_, idx_] := Normal[c][[idx]];

HallBasisElements[n : (_Integer?degOK | {_Integer?degOK})] /; $initialized := Map[hallExpr, visibleSet[n]];

HallDimension[n_Integer?degOK] /; $initialized := Length[visibleSet[{n}]];

NonZeroHallCoefficients[s : LieSeries[_?cur, _]] := NonZeroHallCoefficients[s, $MaxDegree];
NonZeroHallCoefficients[s : LieSeries[_?cur, _], n : (_Integer | {_Integer})] :=
  Select[Thread[HallBasisElements[n] -> HallCoefficients[s, n]], ! PossibleZeroQ[Last[#]] &];

ZeroLieSeriesEquations[s : LieSeries[_?cur, _]] := ZeroLieSeriesEquations[s, $MaxDegree];
ZeroLieSeriesEquations[s : LieSeries[_?cur, _], n : (_Integer | {_Integer})] :=
  And @@ Map[(0 == #) &, HallCoefficients[s, n]];

cmBasisQ[] := If[$orderName === "PDegree", True, Message[LieSeries::cmbasis]; False];

(* these use the Hall trees of L(A,B) directly, so they work in both modes *)
ClassicalMechanicsBasisElements[n : (_Integer?degOK | {_Integer?degOK})] /; $initialized && cmBasisQ[] :=
  Map[hallExpr, Select[hallSet[n], cmSurvivesQ]];

ClassicalMechanicsCoefficients[s : LieSeries[_?cur, _]] := ClassicalMechanicsCoefficients[s, $MaxDegree];
ClassicalMechanicsCoefficients[LieSeries[_?cur, c_], n : (_Integer?degOK | {_Integer?degOK})] /; cmBasisQ[] :=
  hallCoefficients[c, Select[hallSet[n], cmSurvivesQ]];

LieSeriesMap[fn_, s : LieSeries[_?cur, _]] :=
  LieSeries[$id, hallVector[Map[fn, HallCoefficients[s]], visibleSet[$MaxDegree]]];

(* ------------------------------------------------------------------------- *)
(* Projections                                                                 *)
(* ------------------------------------------------------------------------- *)

LieSeriesProj[LieSeries[id_?cur, c_], n : (_Integer?degOK | {_Integer?degOK})] :=
  LieSeries[id, c degreeMask[n]];

LieSeriesSet[LieSeries[id_?cur, c1_], LieSeries[id_, c2_], {n_Integer?degOK}] :=
  With[{mask = degreeMask[{n}]}, LieSeries[id, c1 (1 - mask) + c2 mask]];

(* divides the coefficient of each tree by its degree in the first generator *)
integrateX[LieSeries[id_?cur, c_]] :=
  LieSeries[id, Normal[c] Map[If[# == 0, 1, 1/#] &, $degree1Vector]];

(* ------------------------------------------------------------------------- *)
(* Substitution                                                                *)
(* ------------------------------------------------------------------------- *)

Options[LieSeriesSubstitution] = {LieBracketFunction -> Automatic};

LieSeriesSubstitution[s : LieSeries[_?cur, _], rules : {(_String -> _) ..}, opts : OptionsPattern[]] :=
  LieSeriesSubstitution[s, rules, $MaxDegree, opts];

LieSeriesSubstitution[LieSeries[_?cur, c_], rules : {(_String -> _) ..}, n_Integer?degOK,
    OptionsPattern[]] :=
  Module[{values, bracket, hs, coef, nz, value},
    values = Map[Replace[#, Append[rules, _ -> Missing[]]] &, $gens];
    bracket = OptionValue[LieBracketFunction];
    If[bracket === Automatic, bracket = $bracketDefault];
    value[1] = values[[1]]; value[2] = values[[2]];
    value[i_Integer] := value[i] =
      If[$reversed, bracket[value[r[i]], value[l[i]]], bracket[value[l[i]], value[r[i]]]];
    hs = hallSet[n];
    coef = Normal[c][[hs]];
    nz = Flatten[Position[coef, x_ /; x =!= 0, {1}, Heads -> False]];
    If[nz === {}, nz = Range[Length[hs]]];
    If[MemberQ[values, Missing[]],
      Message[LieSeriesSubstitution::gens, $gens]; $Failed,
      Plus @@ (coef[[nz]] Map[value, hs[[nz]]])]
  ];

(* ------------------------------------------------------------------------- *)
(* BCH series                                                                  *)
(* ------------------------------------------------------------------------- *)

bchSeries[] := If[$bchCache =!= None, $bchCache,
  $bchCache = Module[{Z, S, W, Wn, k, n, XX = LieSeries[$gens[[1]]], YY = LieSeries[$gens[[2]]]},
    Z = SparseArray[XX + YY];
    S = SparseArray[1/2 (XX - YY)];
    W[0] = SparseArray[XX + YY];
    W[k_] = SparseArray[LieSeries[0]];
    Do[
      Do[Wn[k] = LieBracket[Z, W[k - 1], {n}]/k;
         W[k] = LieSeriesSet[W[k], Wn[k], {n}],
        {k, n - 1}];
      Z = LieSeriesSet[Z,
        1/n (LieBracket[S, Z, {n}] + Sum[BernoulliB[k] Wn[k], {k, 2, n - 1}]), {n}],
      {n, 2, $MaxDegree}];
    Z]];

symmetricBCHSeries[] := If[$sbchCache =!= None, $sbchCache,
  $sbchCache = Module[{Z, W, Wn, k, n},
    W[0] = SparseArray[LieSeries[$gens[[1]]]];
    Z = SparseArray[W[0] + LieSeries[$gens[[2]]]];
    W[k_] = SparseArray[LieSeries[0]];
    If[$MaxDegree >= 2, W[1] = LieBracket[Z, W[0], {2}]];
    Do[
      Do[Wn[k] = LieBracket[Z, W[k - 1], {n}]/k;
         W[k] = LieSeriesSet[W[k], Wn[k], {n}],
        {k, n - 1}];
      Z = LieSeriesSet[Z, integrateX[Sum[BernoulliB[k] Wn[k], {k, 2, n - 1}]], {n}],
      {n, 3, $MaxDegree}];
    Z]];

$BCHLieSeries := If[$initialized, bchSeries[], Message[LieSeries::noinit]; $Failed];
$SymmetricBCHLieSeries := If[$initialized, symmetricBCHSeries[], Message[LieSeries::noinit]; $Failed];

Options[BCH] = Options[SymmetricBCH] = {LieBracketFunction -> Automatic};

projIfSeries[s : LieSeries[_?cur, _], n_] := LieSeriesProj[s, n];
projIfSeries[other_, _] := other;

BCH[a_, b_, n_Integer?degOK, opts : OptionsPattern[]] /; $initialized :=
  projIfSeries[LieSeriesSubstitution[bchSeries[], {$gens[[1]] -> a, $gens[[2]] -> b}, n, opts], n];
BCH[a_, b_, opts : OptionsPattern[]] /; $initialized :=
  LieSeriesSubstitution[bchSeries[], {$gens[[1]] -> a, $gens[[2]] -> b}, opts];

SymmetricBCH[a_, b_, n_Integer?degOK, opts : OptionsPattern[]] /; $initialized :=
  projIfSeries[LieSeriesSubstitution[symmetricBCHSeries[], {$gens[[1]] -> a, $gens[[2]] -> b}, n, opts], n];
SymmetricBCH[a_, b_, opts : OptionsPattern[]] /; $initialized :=
  LieSeriesSubstitution[symmetricBCHSeries[], {$gens[[1]] -> a, $gens[[2]] -> b}, opts];

(* ------------------------------------------------------------------------- *)
(* exp(ad_A) B and Zassenhaus terms                                            *)
(* ------------------------------------------------------------------------- *)

ExpAd[a : LieSeries[_?cur, _], b : LieSeries[_?cur, _]] := ExpAd[a, b, $MaxDegree];
ExpAd[a : LieSeries[_?cur, _], b : LieSeries[_?cur, _], n_Integer /; n >= 1] :=
  Module[{W, k},
    W[0] = b;
    Do[W[k] = LieBracket[a, W[k - 1]]/k, {k, n - 1}];
    b + Sum[W[k], {k, n - 1}]];

ZassenhausTerms[n_Integer] /; $initialized && 2 <= n <= $MaxDegree :=
  Module[{xx = SparseArray[LieSeries[$gens[[1]]]], yy = SparseArray[LieSeries[$gens[[2]]]], aa, cc, k},
    aa = BCH[-yy, BCH[-xx, xx + yy]];
    cc[2] = LieSeriesProj[aa, {2}];
    Do[aa = BCH[-cc[k - 1], aa]; cc[k] = LieSeriesProj[aa, {k}], {k, 3, n}];
    Table[cc[k], {k, 2, n}]];
ZassenhausTerms[n_] /; $initialized := (Message[ZassenhausTerms::range, $MaxDegree]; $Failed);

(* ------------------------------------------------------------------------- *)
(* Information about the basis                                                 *)
(* ------------------------------------------------------------------------- *)

LieSeriesBasisInfo[] /; ! $initialized := (Message[LieSeries::noinit]; $Failed);
LieSeriesBasisInfo[] := <|
  "Generators" -> $gens,
  "MaxDegree" -> $MaxDegree,
  "LieAlgebra" -> $lieAlgebra,
  "HallOrder" -> $orderName,
  "Reversed" -> $reversed,
  "Dimensions" -> Table[Length[visibleSet[{n}]], {n, $MaxDegree}],
  "HallElements" -> $nHall,
  "BasisStamp" -> $id|>;

(* ------------------------------------------------------------------------- *)
(* Series from an earlier basis                                                *)
(* ------------------------------------------------------------------------- *)

Scan[
  Function[fun,
    fun[args___] /; staleQ[Hold[args]] := (Message[LieSeries::stale]; $Failed)],
  {LieBracket, HallCoefficients, NonZeroHallCoefficients, ZeroLieSeriesEquations,
   LieSeriesSubstitution, LieSeriesMap, ClassicalMechanicsCoefficients, LieSeriesProj, LieSeriesSet, BCH, SymmetricBCH, ExpAd}];

(* ------------------------------------------------------------------------- *)
(* Display                                                                     *)
(*                                                                             *)
(* A Lie series with coefficient vector c is displayed as                      *)
(*     Sum_i c[[i]] [Hall element i],  over Hall trees i,                      *)
(* up to degree $FormatMaxDegree. In the notebook front end the display can    *)
(* be copied and evaluated again as input.                                     *)
(* ------------------------------------------------------------------------- *)

Format[zphi[i_Integer]] :=
  If[i <= 2, $gens[[i]],
    If[$reversed,
      Row[{"[", zphi[r[i]], ",", zphi[l[i]], "]"}],
      Row[{"[", zphi[l[i]], ",", zphi[r[i]], "]"}]]];

displayForm[c_] :=
  Module[{k, nh, d, rest},
    k = If[IntegerQ[$FormatMaxDegree], Max[1, Min[$FormatMaxDegree, $MaxDegree]], $MaxDegree];
    nh = f[k + 1] - 1;
    d = (Normal[c][[1 ;; nh]] $visibleMask[[1 ;; nh]]) . Map[zphi, Range[nh]];
    rest = Normal[c][[nh + 1 ;; $nHall]] $visibleMask[[nh + 1 ;; $nHall]];
    If[Select[rest, # =!= 0 &, 1] =!= {}, Row[{d, " + \[Ellipsis]"}], d]];

LieSeries /: MakeBoxes[s : LieSeries[_?cur, c_], fmt : (StandardForm | TraditionalForm)] :=
  With[{d = displayForm[c]}, ToBoxes[Interpretation[d, s], fmt]];
Format[LieSeries[_?cur, c_], OutputForm] := displayForm[c];

LieSeries /: MakeBoxes[LieSeries[id_Integer?staleIdQ, c_], fmt : (StandardForm | TraditionalForm)] :=
  ToBoxes[Style[Row[{"LieSeries[\[LeftGuillemet]basis ", id, ", no longer valid\[RightGuillemet]]"}], Gray], fmt];

End[];

Protect @@ DeleteCases[Names["LieSeries`*"], "$FormatMaxDegree" | "LieSeries`$FormatMaxDegree"];

EndPackage[];
