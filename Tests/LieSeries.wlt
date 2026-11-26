(* Tests for the LieSeries package.
   Run in Mathematica with
       TestReport["path/to/Tests/LieSeries.wlt"]
   after the package has been loaded (Needs["LieSeries`"]).
   Most tests check mathematical identities rather than particular outputs. *)

(* ---------------------------------------------------------------- helpers *)

seriesEqual[a_, b_] := Union[Simplify[HallCoefficients[a - b]]] === {0};

(* exact exponential and logarithm of nilpotent matrices (M^d = 0 for d x d strictly triangular) *)
nilExp[m_] := Module[{d = Length[m], p = IdentityMatrix[Length[m]], s = IdentityMatrix[Length[m]]},
   Do[p = p . m/k; s = s + p, {k, 1, d - 1}]; s];
nilLog[u_] := Module[{d = Length[u], n = u - IdentityMatrix[Length[u]], p, s},
   p = IdentityMatrix[d]; s = 0 n;
   Do[p = p . n; s = s + (-1)^(k + 1) p/k, {k, 1, d - 1}]; s];

(* 6 x 6 strictly upper triangular matrices: products of 6 factors vanish,
   so series truncated at degree 5 are exact *)
mA = {{0, 1, 2, 0, 1, 3}, {0, 0, -1, 1, 0, 2}, {0, 0, 0, 2, -1, 1}, {0, 0, 0, 0, 1, -2}, {0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 0, 0}};
mB = {{0, -1, 1, 1, 0, 2}, {0, 0, 2, 0, 1, -1}, {0, 0, 0, 1, 1, 0}, {0, 0, 0, 0, -1, 1}, {0, 0, 0, 0, 0, 2}, {0, 0, 0, 0, 0, 0}};

wittDimension[n_] := Total[MoebiusMu[#] 2^(n/#) & /@ Divisors[n]]/n;

$orders = {"Standard", "Lyndon", "Tall", "Left", "Right", "PDegree"};

(* ------------------------------------------------- dimensions (Witt formula) *)

VerificationTest[
  Union @ Flatten @ Table[
    LieSeriesInitialize[{"X", "Y"}, 6, HallOrder -> ord, Reversed -> rev];
    Table[HallDimension[n] == wittDimension[n], {n, 1, 6}],
    {ord, $orders}, {rev, {True, False}}],
  {True},
  TestID -> "Witt-dimensions"]

(* --------------------------------------- Lie algebra axioms on random series *)

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 5];
  A = GenericLieSeries[a] /. a[i_] :> i - 3;
  B = GenericLieSeries[b] /. b[i_] :> (-1)^i (i + 1);
  C1 = GenericLieSeries[c] /. c[i_] :> Mod[i, 4] - 1;
  {seriesEqual[LieBracket[A, B], -LieBracket[B, A]],
   seriesEqual[LieBracket[A, LieBracket[B, C1]] + LieBracket[B, LieBracket[C1, A]] + LieBracket[C1, LieBracket[A, B]], LieSeries[0]]},
  {True, True},
  TestID -> "Antisymmetry-Jacobi"]

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 5, HallOrder -> "Lyndon", Reversed -> False];
  HallCoefficients[GenericLieSeries[a]] === Table[a[i], {i, Length[HallBasisElements[5]]}],
  True,
  TestID -> "Generic-series-coefficients"]

(* ---------------------------------------------------------------- BCH *)

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 4];
  NonZeroHallCoefficients[$BCHLieSeries],
  {"X" -> 1, "Y" -> 1, LieBracket["X", "Y"] -> 1/2,
   LieBracket["X", LieBracket["X", "Y"]] -> 1/12, LieBracket["Y", LieBracket["X", "Y"]] -> -1/12,
   LieBracket["Y", LieBracket["X", LieBracket["X", "Y"]]] -> -1/24},
  TestID -> "BCH-degree-4-Standard"]

VerificationTest[
  Union @ Table[
    LieSeriesInitialize[{"X", "Y"}, 5, HallOrder -> ord, Reversed -> (ord =!= "Lyndon")];
    BCH[mA, mB] == nilLog[nilExp[mA] . nilExp[mB]],
    {ord, $orders}],
  {True},
  TestID -> "BCH-matrices-all-orders"]

VerificationTest[
  LieSeriesInitialize[{"P", "Q"}, 5];
  {SymmetricBCH[mA, mB] == nilLog[nilExp[mA/2] . nilExp[mB] . nilExp[mA/2]],
   seriesEqual[SymmetricBCH[-LieSeries["P"], SymmetricBCH[LieSeries["P"], LieSeries["Q"]]], LieSeries["Q"]]},
  {True, True},
  TestID -> "SymmetricBCH-other-generator-names"]

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 5];
  X = LieSeries["X"]; Y = LieSeries["Y"];
  seriesEqual[BCH[X + Y, -X], LieSeriesSubstitution[$BCHLieSeries, {"X" -> X + Y, "Y" -> -X}]],
  True,
  TestID -> "BCH-series-arguments"]

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 5];
  X = LieSeries["X"]; Y = LieSeries["Y"];
  seriesEqual[BCH[BCH[X, Y], -Y], X],
  True,
  TestID -> "BCH-inverse"]

(* ---------------------------------------------------- exp(ad) and Zassenhaus *)

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 5];
  LieSeriesSubstitution[ExpAd[LieSeries["X"], LieSeries["Y"]], {"X" -> mA, "Y" -> mB}] ==
    nilExp[mA] . mB . nilExp[-mA],
  True,
  TestID -> "ExpAd-matrices"]

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 5];
  nilExp[mA + mB] == Dot @@ Join[{nilExp[mA], nilExp[mB]},
    Map[nilExp[LieSeriesSubstitution[#, {"X" -> mA, "Y" -> mB}]] &, ZassenhausTerms[5]]],
  True,
  TestID -> "Zassenhaus-matrices"]

(* ------------------------------------------------- projections and setters *)

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 5];
  Z = $BCHLieSeries;
  {seriesEqual[LieSeriesProj[Z, 3] + LieSeriesProj[Z, {4}] + LieSeriesProj[Z, {5}], Z],
   seriesEqual[LieSeriesSet[Z, LieSeries[0], {2}], Z - LieSeriesProj[Z, {2}]]},
  {True, True},
  TestID -> "Proj-Set"]

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 4];
  A = GenericLieSeries[a] /. a[i_] :> i - 4;
  {HallCoefficients[LieSeriesMap[g, A]] === Map[g, HallCoefficients[A]],
   seriesEqual[LieSeriesMap[Identity, A], A],
   HallCoefficients[LieSeriesMap[Abs, A]] === Abs[HallCoefficients[A]]},
  {True, True, True},
  TestID -> "LieSeriesMap"]

(* ------------------------------- Lie algebra of classical mechanics *)
(* With HallOrder -> "PDegree", the elements given by ClassicalMechanicsBasisElements
   should form a basis of L(A,B)/I, I the ideal generated by [B,[B,[B,A]]], and the other Hall
   elements should span I (McLachlan & Murua 2019). The two tests below prove this up to
   degree 8, and check the number of elements up to degree 20:
   - the coefficients on the surviving elements of all elements of a spanning set of I vanish,
     so the surviving elements are linearly independent in L(A,B)/I;
   - their number equals dim L(A,B)/I in each degree (Table 5 of the paper, third column). *)

$cmDimensionsTable5 = {2, 1, 2, 2, 4, 5, 10, 14, 25, 39, 69, 110, 193, 320, 555, 938, 1630, 2786, 4852, 8370};

VerificationTest[
  LieSeries`Private`cmDimensions[20],
  {Table[wittDimension[n], {n, 1, 20}], $cmDimensionsTable5},
  TestID -> "ClassicalMechanics-dimensions"]

VerificationTest[
  Module[{a, b, spanI},
    LieSeriesInitialize[{"A", "B"}, 8, HallOrder -> "PDegree", Reversed -> False];
    a = LieSeries["A"]; b = LieSeries["B"];
    (* spanI[n] spans the component of degree n of the ideal I *)
    spanI[4] = {LieBracket[b, LieBracket[b, LieBracket[b, a]]]};
    Do[spanI[n] = Flatten[Table[LieBracket[g, y, {n}], {g, {a, b}}, {y, spanI[n - 1]}]], {n, 5, 8}];
    {Table[Length[ClassicalMechanicsBasisElements[{n}]], {n, 1, 8}],
     Complement[HallBasisElements[4], ClassicalMechanicsBasisElements[4]],
     Union[Flatten[Table[ClassicalMechanicsCoefficients[y, {n}], {n, 4, 8}, {y, spanI[n]}]]]}],
  (* the only discarded element up to degree 4 is [[[A,B],B],B] = -[B,[B,[B,A]]] *)
  {Take[$cmDimensionsTable5, 8], {LieBracket[LieBracket[LieBracket["A", "B"], "B"], "B"]}, {0}},
  TestID -> "ClassicalMechanics-ideal"]

(* ------------- working in L(A,B)/I with LieAlgebra -> "ClassicalMechanics" *)

VerificationTest[
  LieSeriesInitialize[{"A", "B"}, 8, LieAlgebra -> "ClassicalMechanics"];
  {Table[HallDimension[n], {n, 1, 8}], HallBasisElements[8] === ClassicalMechanicsBasisElements[8],
   LieSeriesBasisInfo[]["LieAlgebra"], LieSeriesBasisInfo[]["HallOrder"], LieSeriesBasisInfo[]["Dimensions"]},
  {Take[$cmDimensionsTable5, 8], True, "ClassicalMechanics", "PDegree", Take[$cmDimensionsTable5, 8]},
  TestID -> "Quotient-basis"]

VerificationTest[
  Module[{cfree, A, B, Z, xI},
    (* coefficients in L(A,B)/I computed in the free Lie algebra with the basis "PDegree" *)
    LieSeriesInitialize[{"A", "B"}, 6, HallOrder -> "PDegree", Reversed -> False];
    cfree = {ClassicalMechanicsCoefficients[$BCHLieSeries], ClassicalMechanicsCoefficients[$SymmetricBCHLieSeries]};
    (* the same coefficients, computed in quotient mode *)
    LieSeriesInitialize[{"A", "B"}, 6, LieAlgebra -> "ClassicalMechanics", Reversed -> False];
    A = LieSeries["A"]; B = LieSeries["B"];
    xI = LieBracket[A, LieBracket[B, LieBracket[B, LieBracket[B, A]]]];
    {{HallCoefficients[$BCHLieSeries], HallCoefficients[$SymmetricBCHLieSeries]} === cfree,
     ClassicalMechanicsCoefficients[$BCHLieSeries] === HallCoefficients[$BCHLieSeries],
     Union[HallCoefficients[xI]], ZeroLieSeriesEquations[xI], NonZeroHallCoefficients[xI]}],
  {True, True, {0}, True, {}},
  TestID -> "Quotient-coefficients"]

VerificationTest[
  (* substitution uses the full element of L(A,B): exact for matrices that do not satisfy [b,[b,[b,a]]] = 0 *)
  LieSeriesInitialize[{"A", "B"}, 5, LieAlgebra -> "ClassicalMechanics"];
  {LieBracket[mB, LieBracket[mB, LieBracket[mB, mA]]] =!= 0 mA,
   BCH[mA, mB] == nilLog[nilExp[mA] . nilExp[mB]]},
  {True, True},
  TestID -> "Quotient-substitution"]

VerificationTest[
  LieSeriesInitialize[{"A", "B"}, 5, LieAlgebra -> "ClassicalMechanics"];
  Module[{U = GenericLieSeries[a] /. a[i_] :> i - 5},
    {HallCoefficients[GenericLieSeries[a]] === Array[a, Length[HallBasisElements[5]]],
     HallCoefficients[LieSeriesMap[g, U]] === Map[g, HallCoefficients[U]],
     HallCoefficients[LieSeriesMap[Abs, U]] === Abs[HallCoefficients[U]]}],
  {True, True, True},
  TestID -> "Quotient-GenericLieSeries-LieSeriesMap"]

VerificationTest[
  {LieSeriesInitialize[{"A", "B"}, 3, LieAlgebra -> "ClassicalMechanics", HallOrder -> "Lyndon"],
   LieSeriesInitialize[{"A", "B"}, 3, LieAlgebra -> "Nilpotent"]},
  {$Failed, $Failed},
  {LieSeriesInitialize::cmorder, LieSeriesInitialize::lieal},
  TestID -> "Quotient-wrong-options"]

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 3];
  Head[ClassicalMechanicsCoefficients[LieSeries["X"]]],
  ClassicalMechanicsCoefficients,
  {LieSeries::cmbasis},
  TestID -> "ClassicalMechanics-wrong-basis"]

(* ------------------------------------------------------- substitution *)

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 3];
  br[f_, g_] := f D[g, x] - g D[f, x];
  Simplify[LieSeriesSubstitution[LieBracket[LieSeries["X"], LieSeries["Y"]], {"X" -> x^2, "Y" -> 1/x},
      LieBracketFunction -> br] - br[x^2, 1/x]],
  0,
  TestID -> "Substitution-custom-bracket"]

VerificationTest[
  LieBracket[{{1, 2}, {3, 4}}, {{0, 1}, {1, 0}}],
  {{1, 2}, {3, 4}} . {{0, 1}, {1, 0}} - {{0, 1}, {1, 0}} . {{1, 2}, {3, 4}},
  TestID -> "Matrix-bracket"]

(* -------------------------------------------- fixes with respect to 2008 version *)

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 3];
  Head[LieSeries["X"] LieSeries["Y"]],
  Times,
  TestID -> "No-product-of-series"]

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 3];
  X = LieSeries["X"]; Y = LieSeries["Y"];
  Length[{X, Y}],
  2,
  TestID -> "Lists-of-series-are-lists"]


VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 3];
  X = LieSeries["X"]; Y = LieSeries["Y"];
  X \[Diamond] Y === LieBracket[X, Y],
  True,
  TestID -> "Diamond-notation"]

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 3];
  Xold = LieSeries["X"];
  LieSeriesInitialize[{"X", "Y"}, 4];
  {LieSeriesQ[Xold], LieSeriesQ[LieSeries["X"]], LieBracket[Xold, LieSeries["Y"]]},
  {False, True, $Failed},
  {LieSeries::stale},
  TestID -> "Stale-series"]

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 3, HallOrder -> "NoSuchOrder"],
  $Failed,
  {LieSeriesInitialize::order},
  TestID -> "Unknown-order"]

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 3, LessFunction -> LyndonLess],
  $Failed,
  {LieSeriesInitialize::optx},
  TestID -> "Unknown-option"]

(* ------------------------------------- structure constants beyond degree 8 *)
(* The BCH data used below stop at degree 8. These tests check the bracket, i.e. the
   structure constants of the Hall basis, and the BCH series up to degree 12. They use
   13 x 13 strictly upper triangular integer matrices: products of 13 factors vanish, so
   Lie series truncated at degree 12 are mapped exactly. *)

SeedRandom[20261125];
randomNil[d_] := Table[If[j > i, RandomInteger[{-3, 3}], 0], {i, d}, {j, d}];
m12A = randomNil[13]; m12B = randomNil[13];
randomSeries[s_Symbol] := GenericLieSeries[s] /. s[_] :> RandomInteger[{-5, 5}];

VerificationTest[
  Union @ Table[
    LieSeriesInitialize[{"X", "Y"}, 12, HallOrder -> cfg[[1]], Reversed -> cfg[[2]]];
    Module[{A = randomSeries[a], B = randomSeries[b], sub},
      sub[s_] := LieSeriesSubstitution[s, {"X" -> m12A, "Y" -> m12B}];
      sub[LieBracket[A, B]] == LieBracket[sub[A], sub[B]]],
    {cfg, {{"Standard", True}, {"Standard", False}, {"Lyndon", False}, {"PDegree", True}}}],
  {True},
  TestID -> "Bracket-matrices-degree-12"]

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, 12];
  {BCH[m12A, m12B] == nilLog[nilExp[m12A] . nilExp[m12B]],
   SymmetricBCH[m12A, m12B] == nilLog[nilExp[m12A/2] . nilExp[m12B] . nilExp[m12A/2]]},
  {True, True},
  TestID -> "BCH-matrices-degree-12"]

VerificationTest[
  (* stored entries {s, u, v} -> C[s, u, v] of [E_u, E_v] = Sum_s C[s, u, v] E_s:
     u < v, |s| = |u| + |v| = n, nonzero integers *)
  Union @ Flatten @ Table[
    LieSeriesInitialize[{"X", "Y"}, 10, HallOrder -> ord, Reversed -> False];
    Table[
      Map[
        Function[rule, With[{s = rule[[1, 1]], u = rule[[1, 2]], v = rule[[1, 3]], w = rule[[2]]},
          u < v && IntegerQ[w] && w != 0 &&
            LieSeries`Private`degree[s] == n ==
              LieSeries`Private`degree[u] + LieSeries`Private`degree[v]]],
        LieSeries`Private`bracketRules[n]],
      {n, 2, 10}],
    {ord, $orders}],
  {True},
  TestID -> "Structure-constants-canonical"]

VerificationTest[
  (* [E_1, E_203] = -E_341 - E_377 - E_404 with E_i = [E_l, E_r] (checked by expansion
     in words); a first draft of the Hall representation gave -2 E_404 *)
  LieSeriesInitialize[{"X", "Y"}, 11, Reversed -> False];
  Module[{e},
    e[i_] := GenericLieSeries[a] /. a[k_] :> Boole[k == i];
    Select[Thread[Range[Length[HallBasisElements[11]]] -> HallCoefficients[LieBracket[e[1], e[203]]]],
      Last[#] =!= 0 &]],
  {341 -> -1, 377 -> -1, 404 -> -1},
  TestID -> "Bracket-regression-degree-11"]

VerificationTest[
  (* the order functions are also used for the normal form of the auxiliary trees *)
  Union @ Table[
    LieSeriesInitialize[{"X", "Y"}, 8, HallOrder -> ord, CheckLessFunction -> True];
    LieSeriesBasisInfo[]["HallOrder"] === ord,
    {ord, $orders}],
  {True},
  TestID -> "Hall-orders-checked"]

VerificationTest[
  (* an invalid order (reversed index order) violates the Hall condition that each Hall tree
     is greater than its right factor: LieSeriesInitialize::order3 must be issued. The order
     is added to the private registry only for this test. *)
  Module[{res},
    LieSeries`Private`$hallOrders["ReverseIndex"] =
      <|"Less" -> (#1 > #2 &), "Description" -> "invalid order, used only in the tests"|>;
    res = Quiet[
      Check[LieSeriesInitialize[{"X", "Y"}, 4, HallOrder -> "ReverseIndex"]; "no message",
        "order3", LieSeriesInitialize::order3],
      {LieSeriesInitialize::order3, General::stop}];
    KeyDropFrom[LieSeries`Private`$hallOrders, "ReverseIndex"];
    res],
  "order3",
  TestID -> "Invalid-order-detected"]

(* ------------------------------ comparison with the published BCH tables *)
(* Tables of F. Casas and A. Murua, J. Math. Phys. 50 (2009) 033513, available at
   https://www.ehu.eus/ccwmuura/bch.html . Put the files
       bchHall20.dat  bchLyndon20.dat  sbchHall19.dat  sbchLyndon19.dat
   (complete, or truncated to their first rows) in the folder Tests/data.
   Each row is  i, i', i'', p, q  with E_1 = X, E_2 = Y, E_i = [E_i', E_i''] and
   coefficient p/q; rows are sorted by degree. The series sum p/q E_i is rebuilt
   with LieBracket and compared with the package, so the test does not depend on
   the Hall order or bracket convention used in the tables. *)

$tableDegree = 8;   (* uses the first 2+1+2+3+6+9+18+30 = 71 rows of each file *)
$dataDir = FileNameJoin[{If[StringQ[$TestFileName], DirectoryName[$TestFileName], Directory[]], "data"}];

tableRows[file_, nrows_] := Module[{st = OpenRead[file], data},
   data = ReadList[st, Number, 5 nrows]; Close[st]; Partition[data, 5]];

tableSeries[rows_] := Module[{e, z = LieSeries[0]},
   e[1] = LieSeries["X"]; e[2] = LieSeries["Y"];
   Do[If[row[[1]] >= 3, e[row[[1]]] = LieBracket[e[row[[2]]], e[row[[3]]]]];
      z = z + row[[4]]/row[[5]] e[row[[1]]],
     {row, rows}];
   z];

tableCheck[name_, series_] := Module[{file = FileNameJoin[{$dataDir, name}]},
   If[! FileExistsQ[file], "missing data file " <> file,
     seriesEqual[tableSeries[tableRows[file, Length[HallBasisElements[$tableDegree]]]], series]]];

VerificationTest[
  LieSeriesInitialize[{"X", "Y"}, $tableDegree];
  {tableCheck["bchHall20.dat", $BCHLieSeries], tableCheck["bchLyndon20.dat", $BCHLieSeries],
   tableCheck["sbchHall19.dat", $SymmetricBCHLieSeries], tableCheck["sbchLyndon19.dat", $SymmetricBCHLieSeries]},
  {True, True, True, True},
  TestID -> "Published-BCH-tables"]
