#!/usr/bin/env wolframscript
(* Print the structure of an MX file: wolframscript -file dump_mx.wls <file.mx> [depth]
   For each level: head, length / dimensions, and short InputForm of the leaves. *)
f = $ScriptCommandLine[[2]];
depth = If[Length[$ScriptCommandLine] >= 3, ToExpression[$ScriptCommandLine[[3]]], 3];
d = Import[f, "MX"];
Print["file: ", f, "  SHA256: ", IntegerString[FileHash[f, "SHA256"], 16, 64], "  $SystemID: ", $SystemID];
show[x_, lvl_, path_] := Module[{s},
   s = ToString[Short[x, 3], InputForm, PageWidth -> 200];
   Print[StringRepeat["  ", Length[path]], path, "  head=", Head[x], If[ListQ[x], "  len=" <> ToString[Length[x]] <>
       "  dims=" <> ToString[Dimensions[x]], ""], "  ", StringTake[s, UpTo[160]]];
   If[ListQ[x] && lvl < depth, Do[show[x[[i]], lvl + 1, Append[path, i]], {i, Min[Length[x], 8]}]]];
show[d, 0, {}];
