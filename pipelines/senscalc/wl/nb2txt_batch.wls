#!/usr/bin/env wolframscript
(* Convert all notebooks below <srcdir> to text in <outdir> (Input/Code cells only), one kernel.
   wolframscript -file nb2txt_batch.wls <srcdir> <outdir> *)
{src, out} = $ScriptCommandLine[[2 ;; 3]];
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
(* also keep section/text titles so the structure is readable *)
cellTxt[Cell[BoxData[b_], "Input" | "Code", ___]] := bx[b];
cellTxt[Cell[t_String, sty : ("Title" | "Chapter" | "Section" | "Subsection" | "Subsubsection"), ___]] := "(* ==== [" <> sty <> "] " <> t <> " ==== *)";
cellTxt[Cell[TextData[t_], sty : ("Title" | "Chapter" | "Section" | "Subsection" | "Subsubsection"), ___]] := "(* ==== [" <> sty <> "] " <> ToString[Cases[{t}, _String, Infinity] // StringJoin] <> " ==== *)";
nbs = FileNames["*.nb", src, Infinity];
Do[
  Module[{rel = StringReplace[StringDrop[f, StringLength[src]], {"/" -> "__", " " -> "_"}], cells},
   cells = Cases[Get[f], c : Cell[_, "Input" | "Code" | "Title" | "Chapter" | "Section" | "Subsection" | "Subsubsection", ___] :> Quiet[cellTxt[c]], Infinity];
   cells = Select[cells, StringQ];
   Export[FileNameJoin[{out, StringTrim[rel, "__"] <> ".txt"}], StringRiffle[cells, "\n\n(* ---- *)\n"] <> "\n", "Text"];
   Print[Length[cells], " cells <- ", f]], {f, nbs}];
