# SPARK data analysis

The data contained within this folder has been downloaded in January 2025 from the [co-add database](https://spark.co-add.org/). Usage of the data remains restricted to the associated Terms of Use.

## File summary

### SPARK Accumulation Data: 
Accumulation of compound inside the cells, helping elucidate the cell membrane permeability of compounds. Relevant columns:
* 'SMILES'
* 'Curated & Transformed Accumulation Data: Accumulated Compound (nmol/1E12 CFU)'
* 'Curated & Transformed Accumulation Data: Accumulated Compound (pg/1E9 CFU)'
* 'Curated & Transformed Accumulation Data: Std Err (mass)'
* 'Curated & Transformed Accumulation Data: Accumulated Compound (ng/mL)'
* 'Curated & Transformed Accumulation Data: Accumulated Compound (μM)'
* 'Curated & Transformed Accumulation Data: Species'
Statistics:
* Total rows: 1649
* Unique SMILES: 464
* Species: Acinetobacter baumannii(3)/Bacillus subtilis(28)/Escherichia coli(927)/Mycobacterium smegmatis(28)/Pseudomonas aeruginosa(34)/Staphylococcus aureus(11)

### SPARK Data Merck & Kyorin Contribution:
The contributed data comes from a single experiment: GyrAB/ParCE inhibitors. Data is reported in MIC at 18-24h (not specified if MIC50 or MIC90). All dat available in MIC (uM). Relevant columns:
* 'SMILES'
* 'Curated & Transformed MIC Data: MIC (µg/mL)'
* 'Curated & Transformed MIC Data: MIC (in µM) (µM)'
* 'Curated & Transformed MIC Data: Species'
Statistics:
* Total rows: 2075
* Unique SMILES: 249
* Species: Acinetobacter baumannii (252)/Escherichia coli(500)/Pseudomonas aeruginosa(1000)/Staphylococcus aureus(249)
Abaumanni data_
* Strains: IID876 (wild type)
* Replicates: some data is available in replicates for the same strain.
* The 40uM cut-off does not make sense as it cannot be calculated for many SMILES

### Data Quave Lab {Emory University} Publications
Contains information related to two publications, Chassagne et al 2020, Front. Pharmacol. and Porras et al 2020, Chem. Rev. Both report plants with antimicrobial activity, The Porras et al 2020 has some SMILES, the Chassagne study does not. Relevant columns:
* 'SMILES'
* 'Porras 2020 Chem. Rev. XX:XXXX: MIC value (µg/mL)'
* 'Porras 2020 Chem. Rev. XX:XXXX: MIC (in µM) (µM)'
* 'Porras 2020 Chem. Rev. XX:XXXX: Species'
Statistics:
* Total rows: 6446
* Unique SMILES: 450
* Species: over 100
Abaumannii data:
* Only three molecules contain MIC data for Abaumannii, no strain reported

### SPARK Data Novartis Contribution
Data is related to three contributed experiments:
* Novartis Efflux Panel April 2019 Contribution: screening of nearly 200 antibiotics, many of which have been approved for human use, against 22 bacterial strains to identify the common properties of drugs effective at evading efflux pumps. Antibacterial activity was assessed using a broth microdilution assay (MIC in uM)
* Novartis GyrAB April 2019 Contribution: high-throughput density-based optical assay that detects bacterial cell growth following treatment with compounds. Results are reported as the compound concentration resulting in 50% growth inhibition compared to the 100% inhibition control. The 50% Growth Inhibitory Concentration (GIC50) was calculated using a dose response curve.
* Novartis LpxADK Contribution: This is a high-throughput density-based optical assay that detects bacterial cell growth following treatment with compounds. Results are reported as the compound concentration resulting in 50% growth inhibition compared to the 100% inhibition control. The 50% Growth Inhibitory Concentration (GIC50) was calculated using a dose response curve.
For the GyrAB and LpxADK we need to obtain the values from the GIC (ug/ml) column and for the Efflux panel we will use the MIC (uM). 

Statistics:
* Total rows: 13757
* Unique SMILES: 1727
* Species (MIC):  Acinetobacter baumannii(555)/Burkholderia cepacia(3)/Enterococcus faecalis(862)/Escherichia coli(5469)/Klebsiella pneumoniae(1021)/Neisseria gonorrhoeae(44)/Pseudomonas aeruginosa(3434)Staphylococcus aureus(278)/Streptococcus pneumoniae(1230)/Escherichia coli1(140)/Pseudomonas aeruginosa(25)
* Species (IC50): Escherichia coli(140)/Pseudomonas aeruginosa(25)

Abaumannii data:
* Most data comes from the Efflux panel
* 3 strains (NB48015 (clinical isolate), NB48015-XYL003 (lab isolate), NB48015-XYL005 (lab isolate)). All molecules in the efflux panel have been screened against the three strains, which might differ in some characteristics as for Ampicillin the MIC values are, respectively, 91.6, 22.9, 11.4  uM. To avoid confounding factors only the clinical isolate will be used. All data for Abaumannii is available in MIC uM.


### SPARK Data Achaogen Contribution
Data corresponds to an study on LpxC inhibitors. Several assays are run in parallel, including MIC, IC50 and even PK for some molecules.  IC50 data is only available for K.pneumoniae and Paeurginosa. MIC Data is separated in two experiments that need to be analysed individually: Curated and Achaogen LpxC.
Relevant columns:
* 'SMILES'
* 'Curated & Transformed MIC Data: MIC (µg/mL)'
* 'Curated & Transformed MIC Data: MIC (in µM) (µM)'
* 'Curated & Transformed MIC Data: Species'
* 'Achaogen LpxC MIC Data: MIC (µg/mL)'
* 'Achaogen LpxC MIC Data: MIC (in µM) (µM)'
* 'Achaogen LpxC MIC Data: Species'
* 'Curated & Transformed IC50 Data: IC50 (in uM) (uM)'
* 'Curated & Transformed IC50 Data: Species'
* 'Achaogen Protein Binding Data: Human Fraction Free (%)'
* 'Achaogen Protein Binding Data: Mouse Fraction Free (%)'
* 'Achaogen Protein Binding Data: Rat Fraction Free (%)'
* 'Achaogen Cytotox Data: IC50 (uM)'
Statistics:
* Total rows: 58749
* Unique SMILES: 1873
* Assays (unique SMILES): MIC (1861) / IC50(488) / Human PB (332) / Mouse (55) / Rat (316) / Cytotox (65) / LpxC MIC50 (80) / LpxC MIC90 (80) /
* Species MIC: Acinetobacter baumannii(1536) / Burkholderia thailandensis(726) / Escherichia coli(9356) / Klebsiella pneumoniae(3059) / Pseudomonas aeruginosa(28963) / Staphylococcus aureus(1713) / Yersinia enterocolitica(1406) 
* Species IC50:  Klebsiella pneumoniae(1) / Pseudomonas aeruginosa(535)
Abaumannii data:
* Molecules in the same row can correspond to more than one strain, needs to be carefully analysed. 
* Relevante experiment sets than will be individually analysed (by species and strain) are: Curated & Transformed MIC Data / Achaogen LpxC MIC Data / Curated and Transformed IC50 Data
* Correlation between Achaogen LpxC MIC Data and  Curated and Transformed IC50 Data is confusing as the Species is different from one another

### SPARK Data CO-ADD Contribution
CoADD data, contains two experiments: MIC (18h incubation) and Single Point Inhibition (18h Inhibition). The single point inhibition is reported as % of inhibition at a specific concentration, the concentrations vary in range from 2 to 300 ug/ml. Relevant columns:
* 'SMILES'
* 'Curated & Transformed MIC Data: MIC (µg/mL)'
* 'Curated & Transformed MIC Data: MIC (in µM) (µM)'
* 'Curated & Transformed MIC Data: Species'
* 'CO-ADD Single Concentration Data (Inhibition): Inhibition % (%)'
* 'CO-ADD Single Concentration Data (Inhibition): Inhibition Concentration (in µM) (µM)'
* 'CO-ADD Single Concentration Data (Inhibition): Species'
Statistics:
* Total rows: 82229
* Unique SMILES: 24090
* Species (MIC): Acinetobacter baumannii(1460) / Escherichia coli(1265) / Klebsiella pneumoniae(1265) / Pseudomonas aeruginosa(2448)
* Species (Single Point): Acinetobacter baumannii(23159) / Escherichia coli(5350) / Klebsiella pneumoniae(5349) / Pseudomonas aeruginosa(40969)
Abaumannii data:
* MIC data related to a single strain (ATCC 19606)
* IC50 data related to a single strain (ATCC 19606), how do we unify the cut-offs across concentrations?

### SPARK IC50 Data
There are two broad categories of data: Curated & Transformed or Extracted and Uploaded. For some molecules we have information in both, from others only in one of them. The results are found either in IC50(uM) or an alternate measure like Ki
Relevant columns:
* 'SMILES'
* 'Curated & Transformed IC50 Data: IC50 (uM) (uM)'
* 'Curated & Transformed IC50 Data: Alternate Inhibition Reporting'
* 'Curated & Transformed IC50 Data: Alternate Inhibition Reporting Measurement (µg/mL)'
* 'Curated & Transformed IC50 Data: Species'
* 'Extracted & Uploaded IC50 Data: IC50 (in uM)'
* 'Extracted & Uploaded IC50 Data: Alternate Inhibition Reporting',
* 'Extracted & Uploaded IC50 Data: Alternate Inhibition Reporting Measurement (ug/mL)'
* 'Extracted & Uploaded IC50 Data: Species'
Statistics:
* Total rows: 5462
* Unique SMILES: 2539
* Species Curated: Acinetobacter baumannii(51) / Aquifex aeolicus(7) / Caulobacter crescentus(24) / Escherichia coli(1647) / Escherichia coli, Staphylococcus aureus(15) / Francisella tularensis(66) / Haemophilus influenzae(14) / Klebsiella pneumoniae(29) / Mycobacterium tuberculosis(210) / Pseudomonas aeruginosa(990) / Salmonella typhimurium(4) / Staphylococcus aureus(184) / Streptococcus pneumoniae(42)
* Species Extracted: Acinetobacter baumannii(120) / Aquifex aeolicus(7) / Bacillus subtilis(23) / Citrobacter freundii(1) / Enterobacter cloacae(30) / Enterococcus faecium(30) / Escherichia coli(1842) / Escherichia coli + Staphylococcus aureus(15) / Francisella tularensis(53) / Francisella tularensis subsp. Tularensis(13) / Haemophilus influenzae(14) / Klebsiella pneumoniae(63) / Mycobacterium tuberculosis(124) / Pseudomonas aeruginosa(526) / Staphylococcus aureus(235) / Staphylococcus aureus subsp. aureus(29) / Streptococcus pneumoniae(42)
Abaumannii data:
* Data refers to two studies only
* IC50 indicates inhibition of a specific target (purified). Assayed targets: GyrA, ParE, OXA-24

### SPARK MIC Data
The MIC file contains two subsets of data: Curated & Transformed and Extracted & Uploaded. Instead of the actual MIC data, we sometimes have GIC and GIC % Inhibition (Growth Inhibition Concentration and % of inhibition at that concentration) Relevant columns:
* 'SMILES'
* 'Curated & Transformed MIC Data: MIC (in µM) (µM)'
* 'Curated & Transformed MIC Data: GIC (µg/mL)'
* 'Curated & Transformed MIC Data: GIC % inhibition'
* 'Curated & Transformed MIC Data: Species'
* 'Extracted & Uploaded MIC Data: MIC (in µM) (µM)'
* 'Extracted & Uploaded MIC Data: GIC (Growth Inhibition Concentration) (µg/mL)'
* 'Extracted & Uploaded MIC Data: GIC % Inhibition (%)'
* 'Extracted & Uploaded MIC Data: Species'
Statistics:
* Total rows: 137572
* Unique SMILES: 24581
* Species: more than 60
Abaumannii data:
* The data in Curated & Transformed corresponds mostly to the data in individual files (with an additional 345 SMILES) and the data in Extracted contains 415 additional SMILES, and only a few overlaps with the individual files.


### SPARK Data Compounds & Physicochemical Properties
Information about Physicochemical properties of compounds, no bioactivity available
Statistics:
* Total rows: 162663
* Unique SMILES: 157021

# Data cleaning process
We will start from preprocessed files where rows without SMILES or SMILES that cannot be standardised are eliminated. This will save time in removing duplicates etc that are not parsable from the start.

## MIC
Only the columns containing MIC data will be taken into account. MIC data will be represented as uM (all files represent MIC both in ug/mL and uM). Datasets will be cleaned up individually (by contributor source) and then merged together.
The general steps are:
* Keep only pathogen of interest
* Binarise data based on specified cut-offs
* Merge replicates inside each dataset. The merging is done based on the binarised data.
