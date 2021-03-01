Result files (.res).
LEGEND:
- dcr: Digital Clinical Record, synonym of Electronic Health Record;
- lab: laboratory, used instead of Laboratory Results;
- goHome: Early Discharge;
- MV: Mechanical Ventilation, used earlier instead of Oxygen Support intensification;
- bio: biological agents used, used earlier instead of Rescue Therapy.

INSTRUCTIONS (ipython):
- navigate to corresponding folder (current one) using cd (change directory);
- type: "R = load(key)" [ENTER]; this will upload all data to variable R and print screen information;
- use: "R.sets" to check all sections (property "sets" of structures enumerates all branches);
- refer to "key" defined in code to read sets, e.g., god = goHome dcr; godl = goHome dcr lab (dcr subset using only lab), etc.

SETS EXPLAINED (order of listing):
- god:	EHR Early Discharge;
- godc:	EHR Early Discharge, CoV2 confirmed;
- godl: EHR LR subset Early Discharge (confirmatory set);
- gol:	LR Early Discharge;
- golc: LR Early Discharge, CoV2 confirmed;
- icdc: EHR ICU admission, CoV2 confirmed;
- icdlc:EHR LR subset ICU admission, CoV2 confirmed (confirmatory set);
- iclc:	EHR ICU admission, CoV2 confirmed;
- mvdc:	EHR Oxygen Support intensification, CoV2 confirmed;
- mvdlc:EHR LR subset Oxygen Support intensification, CoV2 confirmed;
- mvlc:	LR Oxygen Support intensification, CoV2 confirmed.