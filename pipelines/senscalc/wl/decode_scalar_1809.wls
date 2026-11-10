#!/usr/bin/env wolframscript
(* Decode SensCalc's dark-scalar lifetime table: wolframscript -file decode_scalar_1809.wls <SensCalc tree> <out dir>
   phenomenology/Scalar/decay widths/ctauS.mx is a list of {description, {{m_GeV, ctau_m at theta^2 = 1}, ...}};
   codes/Acceptances/Scalar.nb uses ctau(m, theta2) = 10^Interpolation[Log10[table], InterpolationOrder -> 1][Log10[m]]/theta2
   with the first description as "Default". Writes ctauS_<i>.csv (mass_GeV, ctau_m) per description and
   ctauS_descriptions.txt; also dumps the structure of BrRatioScalar.mx to BrRatioScalar_structure.txt. *)
{root, out} = $ScriptCommandLine[[2 ;; 3]];
If[! DirectoryQ[out], CreateDirectory[out, CreateIntermediateDirectories -> True]];
dir = FileNameJoin[{root, "phenomenology", "Scalar", "decay widths"}];
c = Import[FileNameJoin[{dir, "ctauS.mx"}], "MX"];
Print["ctauS.mx: ", Length[c], " descriptions; SHA256 ", IntegerString[FileHash[FileNameJoin[{dir, "ctauS.mx"}], "SHA256"], 16, 64]];
Export[FileNameJoin[{out, "ctauS_descriptions.txt"}], StringRiffle[ToString[#[[1]], InputForm] & /@ c, "\n"], "Text"];
Do[
  t = N[c[[i, 2]]];
  Print["  ", i, ": ", ToString[c[[i, 1]], InputForm], "  dims ", Dimensions[t], "  m range ", MinMax[t[[All, 1]]]];
  Export[FileNameJoin[{out, "ctauS_" <> ToString[i] <> ".csv"}], Prepend[t, {"mass_GeV", "ctau_m"}], "CSV"],
  {i, Length[c]}];
b = Import[FileNameJoin[{dir, "BrRatioScalar.mx"}], "MX"];
s = {"BrRatioScalar.mx: head " <> ToString[Head[b]] <> ", length " <> ToString[Length[b]]};
Do[AppendTo[s, ToString[i] <> ": " <> StringTake[ToString[Short[b[[i]], 6], InputForm, PageWidth -> 250], UpTo[1500]]], {i, Min[Length[b], 6]}];
Export[FileNameJoin[{out, "BrRatioScalar_structure.txt"}], StringRiffle[s, "\n"], "Text"];
Print["done"];
