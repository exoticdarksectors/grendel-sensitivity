#!/usr/bin/env wolframscript
(* Print the Input/Code cells of a notebook as linear text (no front end needed), keeping line breaks:
   wolframscript -file nb2txt.wls in.nb out.txt *)
{in, out} = $ScriptCommandLine[[2 ;; 3]];
(* strings may embed linear box syntax (\!\(...\)); flatten it so Export does not draw it in 2D *)
lin[s_String] := FixedPoint[StringReplace[#, {
     RegularExpression["(Super|Sub)scriptBox\\[([^\\[\\],]*), ([^\\[\\],]*)\\]"] :> "$2" <> If["$1" == "Super", "^", "_"] <> "$3",
     RegularExpression["RowBox\\[\\{([^\\[\\]\\{\\}]*)\\}\\]"] :> StringReplace["$1", ", " -> ""]}] &,
   StringReplace[s, {"\\!\\(\\*" -> "", "\\(" -> "", "\\)" -> "", "\\*" -> ""}]];
bx[s_String] := Module[{t = ToString[s, InputForm, PageWidth -> Infinity]},
   If[StringFreeQ[t, "\\!"], s, lin[StringReplace[StringTake[t, {2, -2}], "\\\"" -> "\""]]]];
bx[RowBox[l_List]] := StringJoin[bx /@ l];
bx[SubscriptBox[a_, b_, ___]] := bx[a] <> "_" <> bx[b];
bx[SuperscriptBox[a_, b_, ___]] := bx[a] <> "^" <> bx[b];
bx[SubsuperscriptBox[a_, b_, c_, ___]] := bx[a] <> "_" <> bx[b] <> "^" <> bx[c];
bx[FractionBox[a_, b_, ___]] := "(" <> bx[a] <> ")/(" <> bx[b] <> ")";
bx[SqrtBox[a_, ___]] := "Sqrt[" <> bx[a] <> "]";
bx[RadicalBox[a_, n_, ___]] := "(" <> bx[a] <> ")^(1/" <> bx[n] <> ")";
bx[OverscriptBox[a_, b_, ___]] := bx[a] <> "^(over " <> bx[b] <> ")";
bx[UnderscriptBox[a_, b_, ___]] := bx[a] <> "_(under " <> bx[b] <> ")";
bx[UnderoverscriptBox[a_, b_, c_, ___]] := bx[a] <> "_(" <> bx[b] <> ")^(" <> bx[c] <> ")";
bx[GridBox[rows_List, ___]] := "{" <> StringRiffle[("{" <> StringRiffle[bx /@ #, ", "] <> "}") & /@ rows, ", "] <> "}";
bx[(TagBox | InterpretationBox | StyleBox | FormBox | AdjustmentBox | ButtonBox | FrameBox | PaneBox | ItemBox)[a_, ___]] := bx[a];
bx[TemplateBox[args_List, tag_, ___]] := tag <> "[" <> StringRiffle[bx /@ args, ", "] <> "]";
bx[DynamicBox[___] | GraphicsBox[___] | Graphics3DBox[___]] := "<graphics>";
bx[x_] := ToString[x, InputForm];
cells = Cases[Get[in], Cell[BoxData[b_], "Input" | "Code", ___] :> bx[b], Infinity];
Export[out, StringRiffle[cells, "\n\n(* ---- *)\n"] <> "\n", "Text"];
Print[Length[cells], " cells -> ", out];
