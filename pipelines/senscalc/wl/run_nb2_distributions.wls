#!/usr/bin/env wolframscript
(* Notebook 2 (LLP angle-energy distributions) for ALP-fermion, headless.
   wolframscript -file run_nb2_distributions.wls <SensCalc tree> <facility> <modes, comma-separated> <nsim> [nKernels]
   e.g. ... LHC B,Eta,EtaPr,KS,Omega,RhoCh,Bremsstrahlung 1000000 4
   Output: <tree>/spectra/New physics particles spectra/ALP-fermion/DoubleDistr_ALP-fermion_<mode>_<facility>.m
   (Fragmentation and Drell-Yan are pregenerated in the repository and are not recomputed.) *)
{Drv`root, Drv`facility, Drv`modes, Drv`nsim} = $ScriptCommandLine[[2 ;; 5]];
Drv`modes = StringSplit[Drv`modes, ","];
Drv`nsim = ToExpression[Drv`nsim];
Drv`nk = If[Length[$ScriptCommandLine] >= 6, ToExpression[$ScriptCommandLine[[6]]], 4];
Drv`root = AbsoluteFileName[Drv`root];
Get[FileNameJoin[{DirectoryName[$InputFileName], "senscalc_env.wl"}]];  (* repro: FeynCalc path, gcc *)
LaunchKernels[Drv`nk];
SCEnv`subkernelEnv[];  (* repro *)
Print["subkernels: ", Length[Kernels[]], "  ", DateString[], "  modes: ", Drv`modes, "  nsim: ", Drv`nsim];

Emu`$answers = {
   {"selectionDialogExplanationOneOption", "Select the facility for which tabulated distributions", Function[l, If[MemberQ[l, Drv`facility], Drv`facility, Abort[]]]},
   {"selectionDialogExplanation", "Select one or more LLPs for", {"ALP-fermion"}},
   {"selectionDialogExplanation", "Select production channels for ALP-fermion", Function[l,
      If[! SubsetQ[l, Drv`modes], Print["EMU| requested modes not offered: ", Complement[Drv`modes, l]]; Abort[]]; Drv`modes]}};
Emu`$inline = {
   {"Number of simulated events per channel", Emu`Ret[Drv`nsim]}};

Get[FileNameJoin[{DirectoryName[$InputFileName], "senscalc_headless.wl"}]];
CheckAbort[Emu`runMaster[FileNameJoin[{Drv`root, "2. LLP distribution.nb"}], {"launcher"}], Print["EMU| ABORTED"]];
Print["files: ", Column[FileNames["*", FileNameJoin[{Drv`root, "spectra", "New physics particles spectra", "ALP-fermion"}], 1]]];
Print["done ", DateString[]];
