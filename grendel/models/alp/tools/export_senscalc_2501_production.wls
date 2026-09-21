#!/usr/bin/env wolframscript
(* Decode the pinned SensCalc 2501.04525 LHC production inputs. *)

ClearAll[
  die, render, plainString, sha256, checkedImport, numericAtMass, writeJSON
];

die[message_] := (
  WriteString[$Output, "error: " <> ToString[message, InputForm] <> "\n"];
  Exit[1]
);

render[expression_] := ToString[expression, InputForm, PageWidth -> Infinity];
plainString[expression_] := If[StringQ[expression], expression, render[expression]];
sha256[path_] := IntegerString[FileHash[path, "SHA256"], 16, 64];

checkedImport[path_] := Module[{value},
  value = Quiet[Check[Import[path, "MX"], $Failed]];
  If[value === $Failed, die["failed to import MX file " <> path]];
  value
];

numericAtMass[expression_, mass_] := Module[{value},
  value = expression /. ((symbol_Symbol /;
       MemberQ[{"ma", "mLLP"}, SymbolName[Unevaluated[symbol]]]) :> mass);
  If[MatchQ[value, _Function | _InterpolatingFunction], value = value[mass]];
  value = Quiet[Check[N[value], Indeterminate]];
  If[!NumericQ[value] || !TrueQ[Abs[Im[N[value]]] < 10^-10],
    die["non-numeric production input at mass " <> ToString[mass] <>
      ": " <> render[expression]]
  ];
  Re[N[value]]
];

writeJSON[path_, value_] := Module[{result},
  result = Quiet[Check[Export[path, value, "RawJSON"], $Failed]];
  If[result === $Failed, die["failed to write " <> path]]
];

If[Length[$ScriptCommandLine] != 3,
  die["usage: export_senscalc_2501_production.wls SENSCALC_ROOT OUTPUT_DIR"]
];

sensCalcRoot = ExpandFileName[$ScriptCommandLine[[2]]];
outputDir = ExpandFileName[$ScriptCommandLine[[3]]];
If[!DirectoryQ[outputDir], CreateDirectory[outputDir, CreateIntermediateDirectories -> True]];

sourceFiles = <|
  "sensitivity_notebook" -> "3. ALP-fermion sensitivity.nb",
  "generic_notebook" -> FileNameJoin[{"codes", "generic.nb"}],
  "production_notebook" -> FileNameJoin[{
    "codes", "LLP distribution", "prod-pheno-ALP-fermion.nb"
  }],
  "experiments" -> FileNameJoin[{"codes", "experiments.nb"}],
  "coefficients" -> FileNameJoin[{
    "phenomenology", "ALP-fermion", "Production probabilities",
    "Coefficients-ALP-fermion-Lambda=1.-TeV.mx"
  }],
  "fragmentation_probability" -> FileNameJoin[{
    "phenomenology", "ALP-fermion", "Production probabilities",
    "ProductionProbability-Fragmentation-ALP-fermion.mx"
  }],
  "light_meson_decays" -> FileNameJoin[{
    "phenomenology", "ALP-fermion", "Production probabilities",
    "BrRatios-Msquared-LightMesonDecays-ALP-fermion-Lambda=1.-TeV.mx"
  }],
  "bremsstrahlung_lhc" -> FileNameJoin[{
    "phenomenology", "ALP-fermion", "Production probabilities",
    "ProdProb_ALP-fermion_Bremsstrahlung-AP_LHC.m"
  }],
  "drell_yan_lhc" -> FileNameJoin[{
    "phenomenology", "ALP-fermion", "Production probabilities",
    "sigmaDrellYan_ALP-fermion_LHC.txt"
  }],
  "fragmentation_lhc_grid" -> FileNameJoin[{
    "spectra", "New physics particles spectra", "ALP-fermion", "Pregenerated",
    "DoubleDistr_ALP-fermion_Fragmentation_LHC.m"
  }],
  "eta_lhc_spectrum" -> FileNameJoin[{
    "spectra", "SM particles", "DoubleDistr_LHC_Eta.txt"
  }],
  "eta_prime_lhc_spectrum" -> FileNameJoin[{
    "spectra", "SM particles", "DoubleDistr_LHC_EtaPr.txt"
  }],
  "omega_lhc_spectrum" -> FileNameJoin[{
    "spectra", "SM particles", "DoubleDistr_LHC_Omega.txt"
  }],
  "rho_charged_lhc_spectrum" -> FileNameJoin[{
    "spectra", "SM particles", "DoubleDistr_LHC_RhoCh.txt"
  }],
  "ks_lhc_spectrum" -> FileNameJoin[{
    "spectra", "SM particles", "DoubleDistr_LHC_KS.txt"
  }]
|>;

expectedHashes = <|
  "sensitivity_notebook" -> "48883fb7ecce956bd73965709ecf44d0cae286e637067c672c1566d736072ef3",
  "generic_notebook" -> "95cc528268e41b273182377275bf9eeb00d8cc1cb92d149ddeaa9e72967f0d4a",
  "production_notebook" -> "e1a7b91a4e82421081d165561a1b7f23e90aefe451efca32a5a3459b1ac0d0f7",
  "experiments" -> "3e60809105f49d8c5274de9a352cd0063bac9cb8feba31188e5e49ca6e468bdf",
  "coefficients" -> "e6885569ee8d33fb5310723ba2f020bf8b4e3506c2b5ad258b5230f9aa981f43",
  "fragmentation_probability" -> "4804205cfb14c0d40322fcabb7ec65843023ee44e91f6918603bca0e8835e9eb",
  "light_meson_decays" -> "c7468dcf8863a8ac182c2720d34d12b2dd12693700aec032460cd5d678985f71",
  "bremsstrahlung_lhc" -> "435022c25a2bb908ffd5c4b12a5645a0d5aa7f6914c49aab5e33fe0b8b60e3d2",
  "drell_yan_lhc" -> "12e013f6df9a65c07514e75d8f1c47d53b53ebc72b21624ab13eb70afaed185d",
  "fragmentation_lhc_grid" -> "43116651d526358c446af54a16d43ca9b2233fb42ff7ab4bae6946fe21477a79",
  "eta_lhc_spectrum" -> "a7e2323c721ce2f6b4b34571b9eaa956ce2cbb5ec8d5c437a2482cafdbfbf39c",
  "eta_prime_lhc_spectrum" -> "dfe394e4901f9d3ddf1e2e3881abe2f8fe550a1cbd0bc7a759f87b79fdf8882d",
  "omega_lhc_spectrum" -> "86f8a67260340de8091936a602eed44f33c5d78200afbaca5b5943fb5849eb02",
  "rho_charged_lhc_spectrum" -> "b2265c828c8b98762e6b7d99973ed6a97f5f66b01c2cab87708236f9ba4738ff",
  "ks_lhc_spectrum" -> "53aafc2e551fed9d01698ba37c7f230e6c5e66c30cc41e73976cd06dd21784ef"
|>;

sourcePaths = Map[FileNameJoin[{sensCalcRoot, #}] &, sourceFiles];
KeyValueMap[
  Function[{key, path},
    If[!FileExistsQ[path], die["missing source file " <> path]];
    If[sha256[path] =!= expectedHashes[key], die["SHA-256 mismatch for " <> path]]
  ],
  sourcePaths
];

lowMassGrid = N[Join[
  Range[0.22, 1.00, 0.02], Range[1.05, 2.00, 0.05]
]];
islandMassGrid = N[Join[
  Range[0.22, 1.00, 0.02], Range[1.05, 3.00, 0.05],
  Range[3.10, 4.60, 0.10], {4.65, 4.70, 4.75}
]];
fPiGeV = 0.093;
bntCoefficientScaleGeV2 = fPiGeV^2/4.;

coefficientsRaw = checkedImport[sourcePaths["coefficients"]];
If[!ListQ[coefficientsRaw] ||
   !And @@ ((ListQ[#] && Length[#] >= 2 &) /@ coefficientsRaw),
  die["unexpected production-coefficient MX structure"]
];
coefficientNames = plainString /@ coefficientsRaw[[All, 1]];
coefficientExpressions = coefficientsRaw[[All, 2]];
coefficientIds = Table[
  "production_coefficient_" <> IntegerString[index, 10, 3],
  {index, Length[coefficientNames]}
];
coefficientData = Table[
  Prepend[
    bntCoefficientScaleGeV2 *
      (numericAtMass[#, mass] & /@ coefficientExpressions),
    mass
  ],
  {mass, islandMassGrid}
];
Export[
  FileNameJoin[{outputDir, "production_coupling_coefficients_bnt.csv"}],
  Prepend[coefficientData, Join[{"mass_GeV"}, coefficientIds]],
  "CSV"
];
writeJSON[
  FileNameJoin[{
    outputDir, "production_coupling_coefficients_metadata.json"
  }],
  <|
    "coefficient_definition" ->
      "P factor = (1/f_BNT)^2 * coefficient_bnt_GeV2",
    "senscalc_definition" ->
      "P factor = (Fpi*gY/(2*vh))^2 * coefficient_raw",
    "fpi_GeV" -> fPiGeV,
    "raw_to_bnt_coefficient_factor_GeV2" -> bntCoefficientScaleGeV2,
    "rows" -> Table[
      <|
        "source_index" -> index,
        "id" -> coefficientIds[[index]],
        "name" -> coefficientNames[[index]],
        "expression_input_form" -> render[coefficientExpressions[[index]]]
      |>,
      {index, Length[coefficientIds]}
    ]
  |>
];

fragmentationRaw = checkedImport[sourcePaths["fragmentation_probability"]];
If[!ListQ[fragmentationRaw] ||
   !And @@ ((ListQ[#] && Length[#] >= 2 &) /@ fragmentationRaw),
  die["unexpected fragmentation-probability MX structure"]
];
fragmentationFacilities = plainString /@ fragmentationRaw[[All, 1]];
fragmentationExpressions = fragmentationRaw[[All, 2]];
fragmentationIds = Table[
  "fragmentation_" <> IntegerString[index, 10, 3],
  {index, Length[fragmentationFacilities]}
];
fragmentationData = Table[
  Prepend[
    bntCoefficientScaleGeV2 *
      (numericAtMass[#, mass] & /@ fragmentationExpressions),
    mass
  ],
  {mass, lowMassGrid}
];
Export[
  FileNameJoin[{
    outputDir, "fragmentation_probability_coefficients_bnt.csv"
  }],
  Prepend[fragmentationData, Join[{"mass_GeV"}, fragmentationIds]],
  "CSV"
];
writeJSON[
  FileNameJoin[{outputDir, "fragmentation_probability_metadata.json"}],
  <|
    "coefficient_definition" ->
      "P_fragmentation = (1/f_BNT)^2 * coefficient_bnt_GeV2",
    "senscalc_definition" ->
      "P_fragmentation = (Fpi*gY/(2*vh))^2 * coefficient_raw",
    "fpi_GeV" -> fPiGeV,
    "raw_to_bnt_coefficient_factor_GeV2" -> bntCoefficientScaleGeV2,
    "rows" -> Table[
      <|
        "source_index" -> index,
        "id" -> fragmentationIds[[index]],
        "facility" -> fragmentationFacilities[[index]],
        "expression_input_form" -> render[fragmentationExpressions[[index]]]
      |>,
      {index, Length[fragmentationIds]}
    ]
  |>
];

lightRaw = checkedImport[sourcePaths["light_meson_decays"]];
If[!ListQ[lightRaw] ||
   !And @@ ((ListQ[#] && Length[#] >= 2 &) /@ lightRaw),
  die["unexpected light-meson MX structure"]
];
lightProcesses = plainString /@ lightRaw[[All, 1]];
lightBranchingExpressions = lightRaw[[All, 2]];
lightMatrixExpressions = (If[Length[#] >= 3, #[[3]], 0.] &) /@ lightRaw;
lightIds = Table[
  "light_meson_" <> IntegerString[index, 10, 3],
  {index, Length[lightProcesses]}
];
lightData = Table[
  Prepend[
    bntCoefficientScaleGeV2 *
      (numericAtMass[#, mass] & /@ lightBranchingExpressions),
    mass
  ],
  {mass, lowMassGrid}
];
Export[
  FileNameJoin[{outputDir, "meson_branching_coefficients_bnt.csv"}],
  Prepend[lightData, Join[{"mass_GeV"}, lightIds]],
  "CSV"
];
writeJSON[
  FileNameJoin[{outputDir, "meson_decay_channels.json"}],
  <|
    "coefficient_definition" ->
      "BR(meson -> a + X) = (1/f_BNT)^2 * coefficient_bnt_GeV2",
    "senscalc_definition" ->
      "BR(meson -> a + X) = (Fpi*gY/(2*vh))^2 * coefficient_raw",
    "fpi_GeV" -> fPiGeV,
    "raw_to_bnt_coefficient_factor_GeV2" -> bntCoefficientScaleGeV2,
    "rows" -> Table[
      <|
        "source_index" -> index,
        "id" -> lightIds[[index]],
        "process_name" -> lightProcesses[[index]],
        "branching_expression_input_form" ->
          render[lightBranchingExpressions[[index]]]
      |>,
      {index, Length[lightIds]}
    ]
  |>
];
writeJSON[
  FileNameJoin[{outputDir, "meson_decay_matrix_elements.json"}],
  Table[
    <|
      "source_index" -> index,
      "id" -> lightIds[[index]],
      "process_name" -> lightProcesses[[index]],
      "expression_input_form" -> render[lightMatrixExpressions[[index]]],
      "expression_c_form" -> ToString[
        CForm[lightMatrixExpressions[[index]]], PageWidth -> Infinity
      ]
    |>,
    {index, Length[lightIds]}
  ]
];

fragmentationGrid = checkedImport[sourcePaths["fragmentation_lhc_grid"]];
If[!MatrixQ[fragmentationGrid, NumericQ] || Length[First[fragmentationGrid]] =!= 4,
  die["unexpected pregenerated fragmentation-grid structure"]
];
Export[
  FileNameJoin[{outputDir, "fragmentation_lhc_grid.csv"}],
  Prepend[N[fragmentationGrid], {"mass_GeV", "theta_rad", "energy_GeV", "density"}],
  "CSV"
];

drellYan = Quiet[Check[Import[sourcePaths["drell_yan_lhc"], "Table"], $Failed]];
If[drellYan === $Failed || !MatrixQ[drellYan, NumericQ] ||
   Length[First[drellYan]] =!= 4,
  die["unexpected Drell-Yan table structure"]
];
Export[
  FileNameJoin[{
    outputDir, "drell_yan_lhc_cross_section_coefficients.csv"
  }],
  Prepend[N[drellYan], {
    "mass_GeV", "sigma_coefficient_central_pb",
    "sigma_coefficient_lower_pb", "sigma_coefficient_upper_pb"
  }],
  "CSV"
];
writeJSON[
  FileNameJoin[{outputDir, "drell_yan_lhc_metadata.json"}],
  <|
    "table_definition" ->
      "sigma_DY = cG_eff^2 * sigma_coefficient",
    "probability_definition" ->
      "P_DY = sigma_DY / sigma_pNucleon",
    "senscalc_sigma_pNucleon_pb" ->
      "51e9 * A_target^(-0.29); see codes/generic.nb",
    "status" ->
      "audit-only; no transverse GRENDEL acceptance is implied"
  |>
];

bremsstrahlungRaw = checkedImport[sourcePaths["bremsstrahlung_lhc"]];
writeJSON[
  FileNameJoin[{outputDir, "bremsstrahlung_lhc_metadata.json"}],
  <|
    "dimensions" -> Dimensions[bremsstrahlungRaw],
    "input_form" -> render[bremsstrahlungRaw],
    "status" -> "audit-only; no transverse GRENDEL acceptance is implied"
  |>
];

outputFiles = {
  "production_coupling_coefficients_bnt.csv",
  "production_coupling_coefficients_metadata.json",
  "fragmentation_probability_coefficients_bnt.csv",
  "fragmentation_probability_metadata.json",
  "meson_branching_coefficients_bnt.csv", "meson_decay_channels.json",
  "meson_decay_matrix_elements.json", "fragmentation_lhc_grid.csv",
  "drell_yan_lhc_cross_section_coefficients.csv",
  "drell_yan_lhc_metadata.json", "bremsstrahlung_lhc_metadata.json"
};
writeJSON[
  FileNameJoin[{outputDir, "PRODUCTION_EXPORT_MANIFEST.json"}],
  <|
    "senscalc_repository" -> "https://github.com/maksymovchynnikov/SensCalc",
    "senscalc_tag" -> "v.1.3.3",
    "senscalc_commit" -> "0bca050633aae16e148d47f21840fa07ff4b8724",
    "phenomenology" -> "arXiv:2501.04525",
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
  " SensCalc production data files to " <> outputDir <> "\n"];
