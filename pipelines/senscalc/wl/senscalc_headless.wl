(* ::Package:: *)
(* senscalc_headless.wl -- run SensCalc master notebooks without a front end.

   The SensCalc notebooks (1. Acceptances, 2. LLP distribution, 3. <model> sensitivity) are driven by dialogs and by
   front-end calls (NotebookFind/SelectionEvaluate on cell tags, NotebookRead of sections).  This file emulates exactly
   those front-end pieces in a plain kernel and answers every dialog from a preset table; the SensCalc code itself is
   evaluated cell by cell as written (cells split at top-level line breaks, as the front end does).  No SensCalc
   physics code is modified.

   Usage (from a .wls driver):
     Emu`$root = "<SensCalc tree>";  Emu`$answers = {...};  Emu`$inline = {...};
     Get["senscalc_headless.wl"];  Emu`runMaster["<root>/1. Acceptances.nb", {"launcher"}]

   Dialog answers:
     Emu`$answers : list of  {kind, phraseSubstring, answer}   kind = "selectionDialog" | "selectionDialogn" |
                     "selectionDialogExplanation" | "selectionDialogExplanationOneOption" | "dropdownDialog" | "choicesDialog" | "combinedDialog"
                     answer = value, or Function[list, ...] applied to the offered list.
     Emu`$inline  : list of  {substring of the dialog's strings, action}  for DialogInput calls written inline;
                     action = Emu`Ret[value] | Emu`Press[buttonLabel] | Emu`PressThen[buttonLabel, Hold[code]]
   Any dialog without a preset answer is printed and the run is aborted. *)

Begin["Emu`"];

$logStream = "stdout";
logs[args___] := WriteString[$logStream, "EMU| ", StringJoin[If[StringQ[#], #, ToString[#, InputForm, PageWidth -> Infinity]] & /@ {args}], "\n"];
short[s_String, n_: 400] := If[StringLength[s] > n, StringTake[s, n] <> "...", s];
short[x_, n_: 400] := short[ToString[x, InputForm, PageWidth -> Infinity], n];

If[! ValueQ[$answers], $answers = {}];
If[! ValueQ[$inline], $inline = {}];
$dirStack = {};
$cellCounter = 0;
$t0 = AbsoluteTime[];

(* ---------- parsing cells like the front end ---------- *)
nl = "\n" | "\[IndentingNewLine]";
toHeld[p_List] := Quiet[ToExpression[RowBox[p], StandardForm, HoldComplete]];
parsesQ[p_List] := toHeld[p] =!= $Failed;
pieces[RowBox[l_List]] := Module[{raw = DeleteCases[SplitBy[l, MatchQ[nl]], {nl ..}], out = {}, buf = {}},
   Do[buf = If[buf === {}, r, Join[buf, {"\n"}, r]]; If[parsesQ[buf], AppendTo[out, buf]; buf = {}], {r, raw}];
   If[buf =!= {}, AppendTo[out, buf]]; out];
pieces[b_] := {{b}};
(* fallback for pieces the kernel parser rejects: drop comment boxes and nested indenting newlines *)
cleanBoxes[b_] := b //. {RowBox[{x___, RowBox[{"(*", ___, "*)"}], y___}] :> RowBox[{x, y}],
    RowBox[{x___, "(*", ___, "*)", y___}] :> RowBox[{x, y}], RowBox[{x___, "\[IndentingNewLine]" | "", y___}] :> RowBox[{x, y}]};
evalPiece[p_List] := Module[{h = toHeld[p]},
   If[h === $Failed, h = Quiet[ToExpression[cleanBoxes[RowBox[p]], StandardForm, HoldComplete]];
     If[h === $Failed, logs["PARSE FAILED (skipped): ", short[RowBox[p], 600]]; Return[Null]]; logs["parsed after cleaning comment boxes"]];
   (* ResourceFunction is locked, so the progress monitor of the notebooks is replaced by Identity here *)
   h = h /. HoldPattern[ResourceFunction["MonitorProgress"]] :> Identity;
   (* Monitor is autoloaded on first use (replacing our DownValues) and then waits for a front end: drop the monitor *)
   h = h /. HoldPattern[Monitor[e_, ___]] :> e;
   ReleaseHold[h]];
evalBoxes[b_] := Scan[evalPiece, pieces[b]];
evalBoxData[BoxData[l_List]] := Scan[evalBoxes, DeleteCases[l, "" | "\n" | "\[IndentingNewLine]"]];
evalBoxData[BoxData[b_]] := evalBoxes[b];
evalBoxData[_] := Null;

(* ---------- sub-notebooks (NotebookEvaluate) ---------- *)
evalNB[file_String] := Module[{f = AbsoluteFileName[file], cells},
   If[f === $Failed, logs["NotebookEvaluate: missing file ", file]; Return[$Failed]];
   cells = Cases[Get[f], Cell[bd_BoxData, "Input" | "Code", ___] :> bd, Infinity];
   Block[{$dirStack = Append[$dirStack, DirectoryName[f]]}, Scan[evalBoxData, cells]];
   installOverrides[]; Null];

(* ---------- master notebook: cells with tags, sections ---------- *)
$master = <||>;
loadMaster[file_String] := Module[{nb = Get[file], cells = {}, secs = <||>, cur = None, walk},
   walk[Cell[CellGroupData[l_List, ___]]] := Scan[walk, l];
   walk[Cell[t_, "Section", ___]] := (cur = StringJoin[Cases[{t}, _String, Infinity]]);
   walk[Cell[t_, "Title" | "Chapter", ___]] := (cur = None);
   walk[Cell[bd_BoxData, "Input" | "Code", opts___]] := Module[{tags = Flatten[{CellTags /. {opts} /. CellTags -> {}}]},
      AppendTo[cells, {tags, bd}]; If[cur =!= None, secs[cur] = Append[Lookup[secs, cur, {}], bd]]];
   walk[_] := Null;
   Scan[walk, First[nb]];
   $master = <|"file" -> file, "dir" -> DirectoryName[file], "cells" -> cells, "sections" -> secs|>;
   logs["loaded master ", FileNameTake[file], ": ", Length[cells], " code cells, ", Length[secs], " sections"]];

evalMasterCell[i_Integer] := Module[{t = AbsoluteTime[]},
   logs["cell ", i, " tags=", $master["cells"][[i, 1]], " start (t=", Round[AbsoluteTime[] - $t0], " s, mem=", Round[MemoryInUse[]/2.^30, 0.01], " GB)"];
   Block[{$dirStack = {$master["dir"]}}, evalBoxData[$master["cells"][[i, 2]]]];
   installOverrides[];
   logs["cell ", i, " done in ", Round[AbsoluteTime[] - t, 0.1], " s"]];
evalTagged[tag_String] := Module[{idx = Flatten[Position[$master["cells"], {tags_List /; MemberQ[tags, tag], _}, {1}, Heads -> False]]},
   logs["SelectionEvaluate tag ", tag, " -> cells ", idx];
   Scan[evalMasterCell, idx]];
evalSection[name_String] := Module[{bds = Lookup[$master["sections"], name, Missing[]]},
   If[MissingQ[bds], logs["section not found: ", name]; Return[$Failed]];
   logs["evaluating section ", name, " (", Length[bds], " cells)"];
   Block[{$dirStack = {$master["dir"]}}, Scan[evalBoxData, bds]]; installOverrides[]; Null];

(* ---------- dialogs ---------- *)
phraseString[p_] := Module[{q = p /. {Row[l_List, ___] :> StringJoin[ToString /@ l], Style[x_, ___] :> x, TextCell[x_, ___] :> x}},
   If[StringQ[q], q, ToString[q]]];
answer[kind_String, phrase_, list_] := Module[{ph = phraseString[phrase], hit, a},
   hit = SelectFirst[$answers, #[[1]] === kind && StringContainsQ[ph, #[[2]]] &, Missing[]];
   If[MissingQ[hit], logs["UNANSWERED DIALOG ", kind, " | ", ph, " | options: ", short[list, 2000]]; Abort[]];
   a = If[Head[hit[[3]]] === Function, hit[[3]][list], hit[[3]]];
   logs["DIALOG ", kind, " | ", short[ph, 200], " | options: ", short[list, 1500], " | answer: ", short[a, 600]];
   If[kind =!= "combinedDialog" && ListQ[list] && ListQ[a] && ! SubsetQ[list, a], logs["WARNING: answer not among the options"]];
   a];

SetAttributes[pressIn, HoldFirst];
pressIn[b_, label_String] := Catch[Module[{acts},
     acts = Cases[Hold[b], Button[label, act_, ___] :> Hold[act], Infinity];
     If[acts === {}, acts = Cases[{b}, Button[label, act_, ___] :> Hold[act], Infinity]];
     If[acts === {}, logs["no button ", label, " in dialog"]; Abort[]];
     ReleaseHold[First[acts]]; Null], "EmuDialog"];
dialogInline[held_Hold] := Module[{strs, key, hit, act},
   strs = StringRiffle[Select[Cases[held, _String, Infinity], StringLength[#] > 3 &], " | "];
   (* dialogs built from previously assigned (local) symbols carry no strings in held form: look at the evaluated form *)
   If[strs === "", strs = StringRiffle[Select[Cases[{ReleaseHold[held]}, _String, Infinity], StringLength[#] > 3 &], " | "]];
   hit = SelectFirst[$inline, StringContainsQ[strs, #[[1]]] &, Missing[]];
   If[MissingQ[hit], logs["UNANSWERED INLINE DIALOG: ", short[strs, 1500]]; Abort[]];
   act = hit[[2]];
   logs["INLINE DIALOG | ", short[strs, 300], " | action: ", short[act, 300]];
   Switch[Head[act],
    Ret, First[act],
    Press | PressThen,
    Module[{r},
     r = Replace[held, {
        Hold[{inits___}, body_, ___] :> Module[{inits}, pressIn[body, act[[1]]]],
        Hold[DynamicModule[{inits___}, body_, ___], ___] :> Module[{inits}, pressIn[body, act[[1]]]],
        Hold[body_, ___] :> pressIn[body, act[[1]]]}];
     If[Head[act] === PressThen, ReleaseHold[act[[2]]]];
     logs["  dialog returned: ", short[r, 300]]; r],
    _, logs["bad inline action ", short[act]]; Abort[]]];

installOverrides[] := (
   Global`selectionDialog[list_, phrase_] := answer["selectionDialog", phrase, list];
   Global`selectionDialogn[list_, phrase_, n_] := answer["selectionDialogn", phrase, list];
   Global`selectionDialogExplanation[list_, d_, phrase_, w_] := answer["selectionDialogExplanation", phrase, list];
   Global`selectionDialogExplanationOneOption[list_, d_, phrase_, w_] := answer["selectionDialogExplanationOneOption", phrase, list];
   Global`dropdownDialog[list_, phrase_] := answer["dropdownDialog", phrase, list];
   Global`choicesDialog[llp_, spec_List] := answer["choicesDialog", llp, spec];
   Global`infoDialog[phrase_] := (logs["INFO DIALOG: ", short[phraseString[phrase], 1000]]; "");
   If[DownValues[Global`combinedDialog] =!= {},
    Global`combinedDialog[llp_, exp_] := answer["combinedDialog", exp, {llp, exp, Global`listDecayDescriptions[llp],
        Global`ProcessesListRefined[llp, Global`ECALoptionGivenExperiment], Keys[DownValues[Global`explanationDecayAcceptance]][[All, 1, 1]]}]];
   (* notebook 2: evaluate a production section by its title instead of reading the open notebook *)
   If[KeyExistsQ[$master, "sections"] && StringContainsQ[$master["file"], "LLP distribution"],
    Global`inp = (# -> # &) /@ Keys[$master["sections"]];
    Global`SectionEvaluator[section_] := evalSection[section]];
   );

(* ---------- front-end functions ---------- *)
(* the system definitions of these functions have patterns more specific than ours (e.g. NotebookDirectory[]),
   so their DownValues are replaced wholesale; LaunchKernels[]/CloseKernels[] and DynamicModule get a rule prepended *)
SetAttributes[replaceDV, HoldAll];
replaceDV[f_Symbol, rules_List] := (Unprotect[f]; DownValues[f] = rules; Protect[f]);
SetAttributes[prependDV, HoldAll];
prependDV[f_Symbol, rules_List] := (Unprotect[f]; DownValues[f] = Join[rules, DownValues[f]]; Protect[f]);
replaceDV[NotebookDirectory, {HoldPattern[NotebookDirectory[___]] :> If[$dirStack === {}, $master["dir"], Last[$dirStack]]}];
replaceDV[NotebookEvaluate, {HoldPattern[NotebookEvaluate[f_String, ___]] :> evalNB[f]}];
replaceDV[EvaluationNotebook, {HoldPattern[EvaluationNotebook[___]] :> "EmuNotebook"}];
replaceDV[NotebookFind, {HoldPattern[NotebookFind[_, tag_String, All, CellTags, ___]] :> ($selTag = tag; Null)}];
replaceDV[SelectionEvaluate, {HoldPattern[SelectionEvaluate[___]] :> evalTagged[$selTag]}];
replaceDV[FrontEndTokenExecute, {HoldPattern[FrontEndTokenExecute[___]] :> Null}];
replaceDV[EmitSound, {HoldPattern[EmitSound[___]] :> Null}];
replaceDV[NotebookDelete, {HoldPattern[NotebookDelete[___]] :> Null}];
replaceDV[SelectionMove, {HoldPattern[SelectionMove[___]] :> Null}];
replaceDV[NotebookRead, {HoldPattern[NotebookRead[___]] :> Null}];
replaceDV[Cells, {HoldPattern[Cells[___]] :> {}}];
replaceDV[CellPrint, {HoldPattern[CellPrint[___]] :> Null}];
replaceDV[NotebookWrite, {HoldPattern[NotebookWrite[___]] :> Null}];
replaceDV[ParentNotebook, {HoldPattern[ParentNotebook[___]] :> Null}];
replaceDV[CreateDialog, {HoldPattern[CreateDialog[___]] :> (logs["CreateDialog called (ignored)"]; Null)}];
replaceDV[MessageDialog, {HoldPattern[MessageDialog[m_, ___]] :> (logs["MESSAGE DIALOG: ", short[phraseString[m], 1000]]; Null)}];
replaceDV[ChoiceDialog, {HoldPattern[ChoiceDialog[args___]] :> (logs["UNANSWERED ChoiceDialog: ", short[{args}, 1000]]; Abort[])}];
replaceDV[DialogInput, {HoldPattern[DialogInput[args___]] :> dialogInline[Hold[args]]}];
replaceDV[DialogReturn, {HoldPattern[DialogReturn[x_]] :> Throw[x, "EmuDialog"], HoldPattern[DialogReturn[]] :> Throw[Null, "EmuDialog"]}];
replaceDV[Monitor, {HoldPattern[Monitor[expr_, ___]] :> expr}];
prependDV[DynamicModule, {HoldPattern[DynamicModule[vars_List, body_, ___]] :> Module[vars, body]}];
prependDV[LaunchKernels, {HoldPattern[LaunchKernels[]] :> (logs["LaunchKernels[] ignored; kernels: ", Length[Kernels[]]]; Kernels[])}];
prependDV[CloseKernels, {HoldPattern[CloseKernels[]] :> (logs["CloseKernels[] ignored"]; {})}];

(* ---------- run a master notebook ---------- *)
(* "launcher" = cell 1, which evaluates the tagged cells itself through NotebookFind/SelectionEvaluate *)
runMaster[file_String, what_List] := Module[{},
   loadMaster[file];
   Do[If[w === "launcher", evalMasterCell[1], If[IntegerQ[w], evalMasterCell[w], evalTagged[w]]], {w, what}];
   logs["runMaster finished, total ", Round[AbsoluteTime[] - $t0], " s"]];

End[];
