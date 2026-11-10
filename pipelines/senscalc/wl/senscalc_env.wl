(* ::Package:: *)
(* senscalc_env.wl -- environment of the reproduction, loaded first by every driver in pipelines/senscalc/wl.
   This is the only addition to the campaign's drivers (campaigns/bc10/faser2_senscalc_20261025/scripts): on macOS
   FeynCalc sat in $UserBaseDirectory/Applications and clang was found automatically; here
     * FeynCalc 10.2.1 lives in $REPRO/env/wolfram/Applications (01_fetch_sources.sh), which is prepended to $Path
       of the main kernel and, through Emu`subkernelEnv[], of every parallel subkernel;
     * the C compiler for CompilationTarget -> "C" is the system gcc (/usr/bin/gcc), not the conda gcc that
       env/activate.sh puts first on PATH.
   No SensCalc code and no physics setting is touched here. *)

Begin["SCEnv`"];
$repro = Environment["REPRO"];
If[! StringQ[$repro], Print["SCEnv: REPRO is not set; source env/activate.sh first"]; Exit[2]];
$apps = FileNameJoin[{$repro, "env", "wolfram", "Applications"}];
If[! MemberQ[$Path, $apps], PrependTo[$Path, $apps]];

Needs["CCompilerDriver`"];
Needs["CCompilerDriver`GCCCompiler`"];
If[FileExistsQ["/usr/bin/gcc"],
  CCompilerDriver`$CCompiler = {"Compiler" -> CCompilerDriver`GCCCompiler`GCCCompiler, "CompilerInstallation" -> "/usr/bin",
     "CompilerName" -> "gcc"}];

(* call after LaunchKernels[n] *)
subkernelEnv[] := With[{apps = $apps}, ParallelEvaluate[If[! MemberQ[$Path, apps], PrependTo[$Path, apps]]; Null]];

Print["SCEnv| ", $Version, " | ", $SystemID, " | licence processes ", $MaxLicenseProcesses, ", subprocesses ",
  $MaxLicenseSubprocesses, " | C compiler ", Quiet[CCompilerDriver`DefaultCCompiler[]], " | FeynCalc dir ",
  FileExistsQ[FileNameJoin[{$apps, "FeynCalc", "FeynCalc.m"}]]];
End[];
