#!/usr/bin/env wolframscript
(* Print the InputForm text of selected code cells (by running index among Input/Code cells), using the FE-like splitter.
   wolframscript -file cellcode.wls nb i1 [i2 ...] *)
in = $ScriptCommandLine[[2]]; idx = ToExpression /@ $ScriptCommandLine[[3 ;;]];
nl = "\n" | "\[IndentingNewLine]";
parsesQ[p_List] := Quiet[ToExpression[RowBox[p], StandardForm, HoldComplete]] =!= $Failed;
pieces[RowBox[l_List]] := Module[{raw = DeleteCases[SplitBy[l, MatchQ[nl]], {nl ..}], out = {}, buf = {}},
   Do[buf = If[buf === {}, r, Join[buf, {"\n"}, r]]; If[parsesQ[buf], AppendTo[out, buf]; buf = {}], {r, raw}];
   If[buf =!= {}, AppendTo[out, buf]]; out];
pieces[b_] := {{b}};
txt[b_] := StringRiffle[(Module[{e = Quiet[ToExpression[RowBox[#], StandardForm, HoldComplete]]},
       If[e === $Failed, "(*UNPARSED*) " <> ToString[RowBox[#], InputForm], ToString[e /. HoldComplete[x___] :> HoldForm[CompoundExpression[x]], InputForm, PageWidth -> 150]]] &) /@ pieces[b], "\n"];
cells = Cases[Get[in], Cell[bd_BoxData, "Input" | "Code", ___] :> bd, Infinity];
Do[Print["(* ======== cell ", i, " ======== *)"];
  Replace[cells[[i]], {BoxData[l_List] :> Scan[Print[txt[#]] &, DeleteCases[l, "" | "\n" | "\[IndentingNewLine]"]], BoxData[b_] :> Print[txt[b]]}], {i, idx}];
