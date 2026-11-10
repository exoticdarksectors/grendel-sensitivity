#!/usr/bin/env wolframscript
(* Decode a pinned SensCalc fermionic-ALP decay model into reviewable formats. *)

ClearAll[
  die, render, plainString, canonicalDecayName, normalizedProducts, sha256,
  checkedImport, evaluateAtMass, writeJSON
];

die[message_] := (
  WriteString[$Output, "error: " <> ToString[message, InputForm] <> "\n"];
  Exit[1]
);

render[expression_] := ToString[expression, InputForm, PageWidth -> Infinity];
plainString[expression_] := If[StringQ[expression], expression, render[expression]];
canonicalDecayName[expression_] := Switch[plainString[expression],
  "2e", "ee",
  "2\[Mu]", "mumu",
  "2\[Tau]", "tautau",
  "2\[Gamma]", "gammagamma",
  "Non-hadronic, total", "nonhadronic_total",
  "Hadronic, total", "hadronic_total",
  "Total", "total",
  _, Null
];
normalizedProducts[products_] := plainString /@
  DeleteCases[Flatten[{products}], Null | "Null"];
sha256[path_] := IntegerString[FileHash[path, "SHA256"], 16, 64];

checkedImport[path_] := Module[{value},
  value = Quiet[Check[Import[path, "MX"], $Failed]];
  If[value === $Failed, die["failed to import MX file " <> path]];
  value
];

evaluateAtMass[expression_, mass_] := Module[{value},
  value = expression /. ((symbol_Symbol /;
       MemberQ[{"ma", "mLLP", "mN"}, SymbolName[Unevaluated[symbol]]]) :> mass);
  If[MatchQ[value, _Function | _InterpolatingFunction], value = value[mass]];
  value = Quiet[Check[N[value], Indeterminate]];
  If[!NumericQ[value] || !TrueQ[Im[N[value]] == 0],
    die["non-numeric branching ratio at mass " <> ToString[mass] <>
      ": " <> render[expression]]
  ];
  Re[N[value]]
];

writeJSON[path_, value_] := Module[{result},
  result = Quiet[Check[Export[path, value, "RawJSON"], $Failed]];
  If[result === $Failed, die["failed to write " <> path]]
];

If[!MemberQ[{3, 4}, Length[$ScriptCommandLine]],
  die["usage: export_senscalc_2501.wls SENSCALC_ROOT OUTPUT_DIR [DECAY_MODEL]"]
];

sensCalcRoot = ExpandFileName[$ScriptCommandLine[[2]]];
outputDir = ExpandFileName[$ScriptCommandLine[[3]]];
decayModel = If[Length[$ScriptCommandLine] == 4,
  $ScriptCommandLine[[4]], "2501"];
If[!DirectoryQ[outputDir], CreateDirectory[outputDir, CreateIntermediateDirectories -> True]];

modelSpecs = <|
  "2501" -> <|
    "file_suffix" -> "Default",
    "phenomenology" -> "arXiv:2501.04525",
    "hashes" -> <|
      "widths" -> "8e438677df92ed3354e7004a667604bd665b182d0ba2780f0ec2262c5caf96d1",
      "branching_ratios" -> "def3041c57db36f6074deab0d373a4fb893453649ac0fca43ec45b2fc6d4e484",
      "matrix_elements" -> "f959c2257fa349e5af6966795db8cbf0da2e3ff8d097b5fccc3316d270ff17ed",
      "acceptance_notebook" -> "d44282084ffe2565c5811254e33265f2575ac80f4bc4b168ff841fda491d8a66"
    |>
  |>,
  "2310_structural" -> <|
    "file_suffix" -> "2310.03524",
    "phenomenology" -> "arXiv:2310.03524",
    "hashes" -> <|
      "widths" -> "8b355aa2496fed8da6f760124cf2b465ce467df88fce70fa710787a5f7efa46a",
      "branching_ratios" -> "66137470087b4db7f8fef6d98a54065867e94cbe3dc052be46cab9ae69a738ed",
      "matrix_elements" -> "89ef54601fc995418e9967d8269862e27cb86852a847bff1669e8948aa12dfca",
      "acceptance_notebook" -> "d44282084ffe2565c5811254e33265f2575ac80f4bc4b168ff841fda491d8a66"
    |>
  |>
|>;
If[!KeyExistsQ[modelSpecs, decayModel],
  die["unknown decay model " <> decayModel]
];
modelSpec = modelSpecs[decayModel];
fileSuffix = modelSpec["file_suffix"];

sourceFiles = <|
  "widths" -> FileNameJoin[{
    "phenomenology", "ALP-fermion", "decay widths",
    "Widths-model-ALP-fermion-scale-1000.-GeV-" <> fileSuffix <> ".m"
  }],
  "branching_ratios" -> FileNameJoin[{
    "phenomenology", "ALP-fermion", "decay widths",
    "Br-ratios-SensCalc-model-ALP-fermion-scale-1000.-GeV-" <> fileSuffix <> ".m"
  }],
  "matrix_elements" -> FileNameJoin[{
    "phenomenology", "ALP-fermion", "decay widths",
    "Matrix-elements-squared-model-ALP-fermion-scale-1000.-GeV-" <> fileSuffix <> ".m"
  }],
  "acceptance_notebook" -> FileNameJoin[{
    "codes", "Acceptances", "ALP-fermion.nb"
  }]
|>;
expectedHashes = modelSpec["hashes"];
sourcePaths = Map[FileNameJoin[{sensCalcRoot, #}] &, sourceFiles];

KeyValueMap[
  Function[{key, path},
    If[!FileExistsQ[path], die["missing source file " <> path]];
    If[sha256[path] =!= expectedHashes[key], die["SHA-256 mismatch for " <> path]]
  ],
  sourcePaths
];

widthRaw = checkedImport[sourcePaths["widths"]];
If[!ListQ[widthRaw] || Length[widthRaw] < 4 || !ListQ[widthRaw[[4]]],
  die["unexpected widths MX structure"]
];
widthPayload = widthRaw[[4]];
If[Length[widthPayload] < 4, die["widths payload has no numeric rows"]];
widthSourceHeaders = widthPayload[[1]];
widthCanonicalNames = canonicalDecayName /@ widthSourceHeaders;
widthSourceMetadata = Take[widthPayload, 3];
widthDataRaw = N[Drop[widthPayload, 3]];
If[!MatrixQ[widthDataRaw, NumericQ], die["widths payload is not a numeric matrix"]];
If[Length[widthSourceHeaders] =!= Length[First[widthDataRaw]],
  die["width header and data-column counts differ"]
];

widthColumnIds = Join[
  {"mass_GeV"},
  Table["width_coefficient_" <> IntegerString[index, 10, 3],
    {index, Length[First[widthDataRaw]] - 1}]
];
requiredWidthNames = {
  "ee", "mumu", "tautau", "gammagamma", "nonhadronic_total",
  "hadronic_total", "total"
};
If[!And @@ (MemberQ[widthCanonicalNames, #] & /@ requiredWidthNames),
  die["required canonical width columns are missing"]
];
widthDataBNT = (Prepend[Rest[#]/4., First[#]] &) /@ widthDataRaw;
Export[
  FileNameJoin[{outputDir, "widths_raw_senscalc.csv"}],
  Prepend[widthDataRaw, widthColumnIds],
  "CSV"
];
Export[
  FileNameJoin[{outputDir, "widths_bnt.csv"}],
  Prepend[widthDataBNT, widthColumnIds],
  "CSV"
];
writeJSON[
  FileNameJoin[{outputDir, "widths_metadata.json"}],
  <|
    "column_ids" -> widthColumnIds,
    "source_headers" -> (plainString /@ widthSourceHeaders),
    "source_headers_input_form" -> (render /@ widthSourceHeaders),
    "columns" -> Table[
      <|
        "source_index" -> index,
        "id" -> widthColumnIds[[index]],
        "source_header" -> plainString[widthSourceHeaders[[index]]],
        "canonical_name" -> widthCanonicalNames[[index]]
      |>,
      {index, Length[widthColumnIds]}
    ],
    "source_metadata_rows_input_form" -> (render /@ widthSourceMetadata),
    "raw_coefficient_definition" -> "Gamma = (gY/(2 vH))^2 * coefficient",
    "bnt_coefficient_definition" -> "Gamma = (1/f_BNT)^2 * coefficient",
    "conversion_raw_to_bnt" -> 0.25
  |>
];

branchingRaw = checkedImport[sourcePaths["branching_ratios"]];
If[!ListQ[branchingRaw] || Length[branchingRaw] < 3,
  die["unexpected branching-ratio MX structure"]
];
branchingProcesses = branchingRaw[[1]];
branchingProducts = branchingRaw[[2]];
branchingExpressions = branchingRaw[[3]];
If[!SameQ[
    Length[branchingProcesses], Length[branchingProducts],
    Length[branchingExpressions]
  ],
  die["branching process, product, and expression counts differ"]
];

channelIds = Table[
  "channel_" <> IntegerString[index, 10, 3],
  {index, Length[branchingProcesses]}
];
channelCanonicalNames = canonicalDecayName /@ branchingProcesses;
If[!And @@ (MemberQ[channelCanonicalNames, #] & /@
    {"ee", "mumu", "tautau", "gammagamma"}),
  die["required canonical branching channels are missing"]
];
massGrid = DeleteDuplicates[widthDataRaw[[All, 1]]];
branchingData = Table[
  Prepend[evaluateAtMass[#, mass] & /@ branchingExpressions, mass],
  {mass, massGrid}
];
Export[
  FileNameJoin[{outputDir, "branching_ratios.csv"}],
  Prepend[branchingData, Join[{"mass_GeV"}, channelIds]],
  "CSV"
];
writeJSON[
  FileNameJoin[{outputDir, "decay_channels.json"}],
  Table[
    <|
      "source_index" -> index,
      "id" -> channelIds[[index]],
      "canonical_name" -> channelCanonicalNames[[index]],
      "process_name" -> plainString[branchingProcesses[[index]]],
      "decay_products" -> normalizedProducts[branchingProducts[[index]]],
      "process_input_form" -> render[branchingProcesses[[index]]],
      "decay_products_input_form" -> render[branchingProducts[[index]]],
      "branching_expression_input_form" -> render[branchingExpressions[[index]]]
    |>,
    {index, Length[channelIds]}
  ]
];

matrixRaw = checkedImport[sourcePaths["matrix_elements"]];
If[!ListQ[matrixRaw] || Length[matrixRaw] < 3,
  die["unexpected matrix-element MX structure"]
];
matrixProcesses = matrixRaw[[1]];
matrixMetadata = matrixRaw[[2]];
matrixExpressions = matrixRaw[[3]];
If[!SameQ[
    Length[matrixProcesses], Length[matrixMetadata],
    Length[matrixExpressions]
  ],
  die["matrix-element process, metadata, and expression counts differ"]
];
matrixIds = Table[
  "matrix_element_" <> IntegerString[index, 10, 3],
  {index, Length[matrixProcesses]}
];
writeJSON[
  FileNameJoin[{outputDir, "matrix_elements.json"}],
  Table[
    <|
      "source_index" -> index,
      "id" -> matrixIds[[index]],
      "canonical_name" -> canonicalDecayName[matrixProcesses[[index]]],
      "process_name" -> plainString[matrixProcesses[[index]]],
      "decay_products" -> normalizedProducts[matrixMetadata[[index]]],
      "process_input_form" -> render[matrixProcesses[[index]]],
      "metadata_input_form" -> render[matrixMetadata[[index]]],
      "expression_input_form" -> render[matrixExpressions[[index]]],
      "expression_c_form" -> ToString[
        CForm[matrixExpressions[[index]]], PageWidth -> Infinity
      ]
    |>,
    {index, Length[matrixIds]}
  ]
];

outputFiles = {
  "widths_raw_senscalc.csv", "widths_bnt.csv", "widths_metadata.json",
  "branching_ratios.csv", "decay_channels.json", "matrix_elements.json"
};
writeJSON[
  FileNameJoin[{outputDir, "EXPORT_MANIFEST.json"}],
  <|
    "senscalc_repository" -> "https://github.com/maksymovchynnikov/SensCalc",
    "senscalc_tag" -> "v1.4",
    "senscalc_commit" -> "4f7210b6e566f5510e2c7269adec3744a6369880",
    "phenomenology" -> modelSpec["phenomenology"],
    "decay_model" -> decayModel,
    "wolfram_version" -> $Version,
    "wolfram_system_id" -> $SystemID,
    "exported_utc" -> DateString[TimeZoneConvert[Now, 0], "ISODateTime"],
    "source_sha256" -> Map[sha256, sourcePaths],
    "output_sha256" -> AssociationMap[
      sha256[FileNameJoin[{outputDir, #}]] &,
      outputFiles
    ]
  |>
];

WriteString[$Output, "Exported " <> ToString[Length[outputFiles]] <>
  " SensCalc data files to " <> outputDir <> "\n"];
