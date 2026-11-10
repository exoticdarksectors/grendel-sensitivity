#!/usr/bin/env wolframscript
(* Convert notebooks to readable code: each Input/Code cell is parsed with MakeExpression and printed in InputForm.
   wolframscript -file nb2code.wls <outdir> <nb1> [<nb2> ...] *)
out = $ScriptCommandLine[[2]]; nbs = $ScriptCommandLine[[3 ;;]];
fmt[HoldComplete[CompoundExpression[a__]]] := StringRiffle[List @@ (ToString[#, InputForm, PageWidth -> 160] & /@ Map[HoldForm, HoldComplete[a]]), ";\n"];
fmt[HoldComplete[a_]] := ToString[HoldForm[a], InputForm, PageWidth -> 160];
fmt[HoldComplete[a__]] := StringRiffle[List @@ (ToString[#, InputForm, PageWidth -> 160] & /@ Map[HoldForm, HoldComplete[a]]), "\n"];
fmt[x_] := "(* PARSE FAILED *) " <> ToString[x, InputForm, PageWidth -> 160];
clean[b_] := b //. {RowBox[{x___, "\[IndentingNewLine]", y___}] :> RowBox[{x, y}], RowBox[{x___, "\n", y___}] :> RowBox[{x, y}], RowBox[{x___, "", y___}] :> RowBox[{x, y}], RowBox[{x_}] :> x};
cellTxt[Cell[BoxData[b_], "Input" | "Code", ___]] := Module[{bb = clean[If[ListQ[b], RowBox[Riffle[b, "\n"]], b]], e},
   e = Quiet[MakeExpression[bb, StandardForm]]; fmt[e]];
cellTxt[Cell[t_, sty : ("Title" | "Chapter" | "Section" | "Subsection" | "Subsubsection"), ___]] :=
  "(* ==== [" <> sty <> "] " <> StringJoin[Cases[{t}, _String, Infinity]] <> " ==== *)";
Do[Module[{rel = StringReplace[FileNameTake[f], {" " -> "_"}], cells},
   cells = Cases[Get[f], c : Cell[_, "Input" | "Code" | "Title" | "Chapter" | "Section" | "Subsection" | "Subsubsection", ___] :> Quiet[cellTxt[c]], Infinity];
   cells = Select[cells, StringQ];
   Export[FileNameJoin[{out, StringReplace[FileNameTake[FileNameDrop[f, -1]], " " -> "_"] <> "__" <> rel <> ".txt"}], StringRiffle[cells, "\n\n(* ---- *)\n"] <> "\n", "Text", CharacterEncoding -> "UTF8"];
   Print[Length[cells], " cells <- ", f]], {f, nbs}];
