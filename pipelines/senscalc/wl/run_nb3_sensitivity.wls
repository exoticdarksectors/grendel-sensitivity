#!/usr/bin/env wolframscript
(* Notebook 3 (ALP-fermion sensitivity), headless: tabulated N_events and sensitivity regions.
   wolframscript -file run_nb3_sensitivity.wls <SensCalc tree> <experiment (acceptance export name)> <extra thresholds, comma-separated> [nKernels]
   Settings: production description "Revised (2501.04525)"; all production channels available for the facility;
   the notebook's own Submit default gives N_ev = 2.3 at the default exposure; extra thresholds are computed afterwards
   with the notebook's SensitivityBlock (same exposure, same channels).
   Output: <tree>/Tabulated Nevents/ALP-fermion/<exp>/*.dat and <tree>/Sensitivity domains/ALP-fermion/<exp>/Sensitivity_*.json|mx *)
{Drv`root, Drv`exp, Drv`extra} = $ScriptCommandLine[[2 ;; 4]];
Drv`extra = ToExpression /@ StringSplit[Drv`extra, ","];
Drv`nk = If[Length[$ScriptCommandLine] >= 5, ToExpression[$ScriptCommandLine[[5]]], 4];
Drv`root = AbsoluteFileName[Drv`root];
Get[FileNameJoin[{DirectoryName[$InputFileName], "senscalc_env.wl"}]];  (* repro: FeynCalc path, gcc *)
LaunchKernels[Drv`nk];
SCEnv`subkernelEnv[];  (* repro *)
Print["subkernels: ", Length[Kernels[]], "  ", DateString[], "  experiment: ", Drv`exp, "  extra thresholds: ", Drv`extra];

Emu`$answers = {
   {"selectionDialogExplanationOneOption", "Select the regime of launching the notebook", "Tabulated Nevents+sensitivity"},
   {"selectionDialogExplanationOneOption", "Select the description of the ALP-fermion production", "Revised (2501.04525)"},
   {"selectionDialogn", "Select the experiments for which the sensitivity will be computed", {Drv`exp}},
   {"selectionDialogn", "Select the experiments:", Function[l, If[MemberQ[l, Drv`exp], {Drv`exp}, Abort[]]]},
   {"dropdownDialog", "Calculate uncertainties of the sensitivity region", False}};
Emu`$inline = {
   {"You have not provided the Production probabilities", Emu`Ret["Continue anyway"]},
   {"Enter the number of proton collisions", Emu`Press["Submit"]}};

Get[FileNameJoin[{DirectoryName[$InputFileName], "senscalc_headless.wl"}]];
CheckAbort[Emu`runMaster[FileNameJoin[{Drv`root, "3. ALP-fermion sensitivity.nb"}], {"launcher"}], Print["EMU| ABORTED"]];
Print["ProductionList: ", Global`ProductionList, "  ProductionListBaseline: ", Global`ProductionListBaseline];
Print["Npot: ", Global`NpotVal[Drv`exp], "  NevMin: ", Global`NevMinVal[Drv`exp], "  SelectedProduction: ", Global`SelectedProduction[Drv`exp]];
Do[Print["extra threshold ", t]; CheckAbort[Global`SensitivityBlock[Drv`exp, N[t], Global`NpotVal[Drv`exp], Global`SelectedProduction[Drv`exp]], Print["EMU| ABORTED extra ", t]],
  {t, Drv`extra}];
Print["files: ", Column[FileNames["*", FileNameJoin[{Drv`root, #}], 3] & /@ {"Tabulated Nevents", "Sensitivity domains"}]];
Print["done ", DateString[]];
