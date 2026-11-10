#!/usr/bin/env wolframscript
(* Print the metadata of a SensCalc acceptance export: wolframscript -file inspect_acceptance.wls <Acceptance_*.m> *)
f = $ScriptCommandLine[[2]];
d = Import[f, "MX"];
Print["file: ", f];
Print["length: ", Length[d], "  facility: ", d[[1]], "  ECAL: ", d[[2]], "  table dims: ", Dimensions[d[[3]]]];
Print["particle numbers: ", TableForm[Transpose[d[[4]]]]];
Print["decay description: ", d[[6]]];
masses = DeleteDuplicates[d[[3]][[All, 1]]];
Print["acceptance mass grid (", Length[masses], "): ", masses];
brvis = d[[5]] /. Global`mLLP -> #1 &;
ctau = d[[7]];
Print["Br_vis(m), ctau(m, y=1) [m]: "];
Do[Print[m, "  ", Quiet[N[d[[5]] /. Global`mLLP -> m]], "  ", Quiet[N[ctau /. {Global`mLLP -> m, Global`y -> 1.}]]], {m, {0.05, 0.3, 0.5, 1., 1.5, 2., 3., 4.}}];
Print["free symbols of ctau: ", Union[Cases[ctau, _Symbol?(Context[#] === "Global`" &), Infinity]]];
