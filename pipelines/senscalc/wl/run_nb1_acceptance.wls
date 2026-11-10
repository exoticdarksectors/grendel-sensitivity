#!/usr/bin/env wolframscript
(* Notebook 1 (acceptance) for ALP-fermion, headless.
   wolframscript -file run_nb1_acceptance.wls <SensCalc tree> <experiment> <export name> <cuts: none|ship|etot100> [nKernels]
   cuts: none    = ">=2" detectable products, no kinematic cut
         ship    = SHiP recommended selection (docs/recommended-selections.md): E > 1 GeV for every particle category, IP < 2.5 m
         etot100 = total energy of the detected pair > 100 GeV (charged, neutral and mixed pairs)
   Output: <tree>/Acceptances/<export name>/Acceptance_<export name>_for_ALP-fermion_decayDescription_Default.m *)
{Drv`root, Drv`target, Drv`exportName, Drv`cuts} = $ScriptCommandLine[[2 ;; 5]];
Drv`nk = If[Length[$ScriptCommandLine] >= 6, ToExpression[$ScriptCommandLine[[6]]], 4];
Drv`root = AbsoluteFileName[Drv`root];
Get[FileNameJoin[{DirectoryName[$InputFileName], "senscalc_env.wl"}]];  (* repro: FeynCalc path, gcc *)
LaunchKernels[Drv`nk];
SCEnv`subkernelEnv[];  (* repro *)
Print["subkernels: ", Length[Kernels[]], "  ", DateString[]];

Drv`cutCode = Switch[Drv`cuts,
   "none", Hold[Null],
   "ship", Hold[Global`eElCutLower = Global`e\[Mu]CutLower = Global`eHadrCutLower = Global`e\[Gamma]CutLower = Global`eh0CutLower = 1.;
     Global`IPcut = 2.5],
   "etot100", Hold[Global`EtotElCutLower = Global`Etot\[Mu]CutLower = Global`EtotHadrCutLower = Global`EtotMixedCutLower =
      Global`Etot\[Gamma]\[Gamma]CutLower = Global`Etoth0h0CutLower = Global`EtotMixedNeutralCutLower = 100.],
   _, Print["unknown cuts ", Drv`cuts]; Exit[1]];

Emu`$answers = {
   {"selectionDialog", "Select the experiments:", Function[l, 
      If[MemberQ[l, Drv`target], {Drv`target}, Take[SortBy[Select[l, StringStartsQ[Drv`target, #] &], -StringLength[#] &], UpTo[1]]]]},
   {"selectionDialogExplanation", "setup(s):", Function[l,  If[MemberQ[l, Drv`target], {Drv`target}, Abort[]]]},
   {"selectionDialogExplanation", "Select LLPs for which acceptances will be computed", {"ALP-fermion"}},
   {"combinedDialog", "", Function[info, Print["combinedDialog info (LLP, exp, descriptions, default channels, acceptance types): ", info];
      {Drv`exportName, "Keep", ">=2", "Hadronized products", {"Default"}, {"Null"}}]}};
Emu`$inline = {
   {"Impose the cuts on decay products", Emu`PressThen["Submit", Drv`cutCode]}};

Get[FileNameJoin[{DirectoryName[$InputFileName], "senscalc_headless.wl"}]];
CheckAbort[Emu`runMaster[FileNameJoin[{Drv`root, "1. Acceptances.nb"}], {"launcher"}], Print["EMU| ABORTED"]];
Print["cuts actually used: ", {Global`eElCutLower, Global`e\[Mu]CutLower, Global`eHadrCutLower, Global`e\[Gamma]CutLower, Global`eh0CutLower,
   Global`IPcut, Global`EtotElCutLower, Global`EtotMixedCutLower, Global`Etot\[Gamma]\[Gamma]CutLower}];
Print["files: ", FileNames["*", FileNameJoin[{Drv`root, "Acceptances"}], 2]];
Print["done ", DateString[]];
