(* ::Package:: *)

(* Timings of the initialization of LieSeries and of the BCH series, to compare
   successive versions of Kernel/LieSeries.wl.

   Usage (in a fresh kernel):

     Get["path/to/Benchmarks/LieSeriesTimings.wl"];

     configs = {{"Standard", 10}, {"Standard", 12}, {"Standard", 14},
                {"Lyndon", 12}, {"PDegree", 12}};

     (* time one version of the kernel and append the results to a file *)
     res = LieSeriesTimings["path/to/Kernel/LieSeries.wl", "1.0 draft 25 Nov", configs];
     SaveTimings[res, "path/to/timings.wl"];

     (* after modifying the kernel, time it again with another label *)
     res = LieSeriesTimings["path/to/Kernel/LieSeries.wl", "1.0 sparse step 3", configs];
     SaveTimings[res, "path/to/timings.wl"];

     (* compare: one row per configuration, one column per label *)
     CompareTimings["path/to/timings.wl"]                       (* initialization time *)
     CompareTimings["path/to/timings.wl", "StructureConstants"]  (* one phase *)
     CompareTimings["path/to/timings.wl", "BCH"]
     CompareFingerprints["path/to/timings.wl"]

   Each version of the kernel is loaded again (Get) before each configuration, so
   different versions can be timed in the same session; all of them use the context
   LieSeries`. For version 0.9 pass "InitializeOptions" -> {LieSeries`GenericLieSeries -> False}
   (in 0.9 the default True precomputes all Hall elements); 0.9 records no phase timings.

   Recorded for each configuration:
     "Init"           time of LieSeriesInitialize (seconds)
     "BCH", "SymmetricBCH"   time of the first evaluation of $BCHLieSeries,
                      $SymmetricBCHLieSeries (with "BCH" -> True)
     "Phases"         LieSeries`Private`$buildTimings (1.0): times of the phases of the
                      initialization ("HallTrees", "AuxiliaryTreesAndCuts",
                      "TreeCoordinates", "StructureConstants", "TensorsAndCleanup", "Total")
                      and sizes ("HallElements", "Trees", "Cuts", "StructureConstantEntries")
     "StructureConstantsHash"   hash of the stored structure constants (meaningful
                      between versions with the Hall representation)
     "BCHHash"        hash of the Hall coefficients of $BCHLieSeries (comparable between
                      all versions)
     "TensorBytes"    ByteCount of the bracket tensors kept after initialization
     "SessionMaxMemory"   MaxMemoryUsed[] (cumulative over the session) *)

BeginPackage["LieSeriesTimings`"];

LieSeriesTimings::usage =
"LieSeriesTimings[kernelFile, label, {{order, degree}, ...}] loads kernelFile and, for each configuration, times LieSeriesInitialize[{\"X\", \"Y\"}, degree, HallOrder -> order] and the BCH series. It returns a list of associations.";
SaveTimings::usage =
"SaveTimings[results, file] appends the results of LieSeriesTimings to file.";
CompareTimings::usage =
"CompareTimings[file] shows the initialization times stored in file, one row per configuration and one column per label. CompareTimings[file, key] shows another quantity: \"BCH\", \"SymmetricBCH\", \"TensorBytes\", or a phase or size recorded in \"Phases\".";
CompareFingerprints::usage =
"CompareFingerprints[file] shows, for each configuration, whether all labels give the same structure constants and the same BCH coefficients.";

Options[LieSeriesTimings] = {"BCH" -> True, "InitializeOptions" -> {}};

Begin["`Private`"];

LieSeriesTimings[file_String, label_String, configs : {{_String, _Integer} ..}, OptionsPattern[]] :=
  With[{bch = TrueQ[OptionValue["BCH"]], initOpts = OptionValue["InitializeOptions"]},
  Map[
    Function[cfg,
      Module[{ord = cfg[[1]], deg = cfg[[2]], tinit, tbch, tsbch, bt, rec},
        Get[file];
        ClearSystemCache[];
        tinit = First[AbsoluteTiming[
          LieSeries`LieSeriesInitialize[{"X", "Y"}, deg, LieSeries`HallOrder -> ord,
            Sequence @@ initOpts]]];
        bt = LieSeries`Private`$buildTimings;
        {tbch, tsbch} = If[bch,
          {First[AbsoluteTiming[LieSeries`$BCHLieSeries;]],
           First[AbsoluteTiming[LieSeries`$SymmetricBCHLieSeries;]]},
          {Missing["NotTimed"], Missing["NotTimed"]}];
        rec = <|
          "Label" -> label, "Order" -> ord, "Degree" -> deg,
          "Init" -> tinit, "BCH" -> tbch, "SymmetricBCH" -> tsbch,
          "Phases" -> If[AssociationQ[bt], bt, Missing["NotRecorded"]],
          "StructureConstantsHash" ->
            Hash[Table[Sort[LieSeries`Private`bracketRules[k]], {k, 2, deg}]],
          "BCHHash" -> If[bch,
            Hash[LieSeries`HallCoefficients[LieSeries`$BCHLieSeries]], Missing["NotTimed"]],
          "TensorBytes" -> ByteCount[Table[LieSeries`Private`bracketTensor[k], {k, 1, deg}]],
          "SessionMaxMemory" -> MaxMemoryUsed[],
          "MathematicaVersion" -> $VersionNumber, "Date" -> DateString[]|>;
        Print[label, " | ", ord, " ", deg, " | init ", NumberForm[tinit, {8, 2}], " s",
          If[NumberQ[tbch], Row[{" | BCH ", NumberForm[tbch, {8, 2}], " s"}], ""]];
        rec]],
    configs]];

SaveTimings[results_List, file_String] :=
  Module[{old = If[FileExistsQ[file], Get[file], {}]},
    Put[Join[old, results], file];
    file];

value[rec_Association, key_String] :=
  Which[
    KeyExistsQ[rec, key], rec[key],
    AssociationQ[rec["Phases"]] && KeyExistsQ[rec["Phases"], key], rec["Phases"][key],
    True, Missing["NotAvailable"]];

show[x_?NumberQ] := If[IntegerQ[x], x, NumberForm[x, {8, 2}]];
show[x_] := x;

CompareTimings[file_String, key_String : "Init"] :=
  Module[{data = Get[file], labels, cfgs, last},
    labels = DeleteDuplicates[data[[All, "Label"]]];
    cfgs = DeleteDuplicates[Map[{#["Order"], #["Degree"]} &, data]];
    (* the last record of each label and configuration *)
    last[lab_, cfg_] := Module[{recs = Select[data,
        #["Label"] === lab && {#["Order"], #["Degree"]} === cfg &]},
      If[recs === {}, "", show[value[Last[recs], key]]]];
    Grid[
      Prepend[
        Map[Join[{#[[1]], #[[2]]}, Table[last[lab, #], {lab, labels}]] &, cfgs],
        Join[{"order", "degree"}, labels]],
      Frame -> All, Alignment -> {{Left, Right, {Right}}}, Background -> {None, {LightGray}}]];

CompareFingerprints[file_String] :=
  Module[{data = Get[file], cfgs},
    cfgs = DeleteDuplicates[Map[{#["Order"], #["Degree"]} &, data]];
    Grid[
      Prepend[
        Map[Function[cfg, Module[{recs = Select[data, {#["Order"], #["Degree"]} === cfg &], c, b},
            (* structure constants: only versions that record phases (Hall representation) *)
            c = DeleteDuplicates[Select[recs, AssociationQ[#["Phases"]] &][[All, "StructureConstantsHash"]]];
            b = DeleteDuplicates[DeleteMissing[recs[[All, "BCHHash"]]]];
            {cfg[[1]], cfg[[2]], Length[recs],
             Which[Length[c] == 0, "-", Length[c] == 1, "same", True, Style["DIFFERENT", Red]],
             Which[Length[b] == 0, "-", Length[b] == 1, "same", True, Style["DIFFERENT", Red]]}]],
          cfgs],
        {"order", "degree", "runs", "structure constants", "BCH coefficients"}],
      Frame -> All, Background -> {None, {LightGray}}]];

End[];

EndPackage[];
