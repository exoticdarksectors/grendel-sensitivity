#!/usr/bin/env wolframscript
(* List code cells of a notebook with section path and CellTags: wolframscript -file nbcells.wls nb out.tsv *)
{in, out} = $ScriptCommandLine[[2 ;; 3]];
nb = Get[in];
secs = {"Title", "Chapter", "Section", "Subsection", "Subsubsection"};
path = <||>; rows = {}; i = 0;
txt[t_] := StringTake[StringReplace[StringJoin[Cases[{t}, _String, Infinity]], {"\n" -> " ", "\t" -> " "}], UpTo[70]];
walk[c : Cell[t_, sty_String, opts___]] := Which[
   MemberQ[secs, sty], path[sty] = txt[t]; Do[KeyDropFrom[path, s], {s, Drop[secs, First@FirstPosition[secs, sty]]}],
   MemberQ[{"Input", "Code"}, sty], i++;
     AppendTo[rows, {i, sty, ToString[Lookup[{opts}, CellTags, ""], InputForm], ToString[Lookup[{opts}, InitializationCell, ""]],
        StringRiffle[Values[KeyTake[path, secs]], " > "], StringTake[StringReplace[ToString[t, InputForm], {"\n" -> " ", "\t" -> " "}], UpTo[90]]}]];
walk[Cell[CellGroupData[l_List, ___]]] := Scan[walk, l];
walk[Notebook[l_List, ___]] := Scan[walk, l];
walk[_] := Null;
walk[nb];
Export[out, Prepend[rows, {"i", "style", "tags", "init", "path", "start"}], "TSV"];
Print[Length[rows], " cells"];
