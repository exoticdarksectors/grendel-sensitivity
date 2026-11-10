#!/usr/bin/env wolframscript
Get[FileNameJoin[{DirectoryName[$InputFileName], "senscalc_env.wl"}]];  (* repro: FeynCalc path, gcc *)
Print["Version: ", $Version, " | ", $VersionNumber, ".", $ReleaseNumber];
Print["ProcessorCount: ", $ProcessorCount, "  Kernels license: ", $MaxLicenseProcesses, " subkernels: ", $MaxLicenseSubprocesses];
Needs["CCompilerDriver`"];
Print["DefaultCCompiler: ", DefaultCCompiler[]];
cf = Compile[{{x, _Real}}, x*x, CompilationTarget -> "C"];
Print["Compile C test: ", cf[2.]];
t = AbsoluteTiming[Quiet[Needs["FeynCalc`"]]];
Print["FeynCalc load (s): ", t[[1]], "  version: ", Quiet[FeynCalc`$FeynCalcVersion]];
