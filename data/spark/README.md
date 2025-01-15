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
MIC data for several species. Not specified if MIC50 or MIC90. Relevant columns:
* 'SMILES'
* 'Curated & Transformed MIC Data: MIC (in µM) (µM)'
* 'Curated & Transformed Accumulation Data: Species'
Statistics:
* Total rows: 2075
* Unique SMILES: 249
* Species: Acinetobacter baumannii (252)/Escherichia coli(500)/Pseudomonas aeruginosa(1000)/Staphylococcus aureus(249)

### Data Quave Lab {Emory University} Publications
Contains information related to two publications, Chassagne et al 2020, Front. Pharmacol. and Porras et al 2020, Chem. Rev. Both report plants with antimicrobial activity, The Porras et al 2020 has sme SMILES. Relevant columns:
* 'SMILES'
* 'Porras 2020 Chem. Rev. XX:XXXX: MIC (in µM) (µM)'
* 'Porras 2020 Chem. Rev. XX:XXXX: Species'
Statistics:
* Total rows: 6446
* Unique SMILES: 450
* Species: over 100

### SPARK Data Novartis Contribution
Contains information for MIC and IC50.Relevant columns:
* 'SMILES'
* 'Curated & Transformed MIC Data: MIC (in µM) (µM)'
* 'Curated & Transformed MIC Data: Species'
* 'Curated & Transformed IC50 Data: IC50 (uM) (uM)'
* 'Curated & Transformed IC50 Data: Species'
Statistics:
* Total rows: 13757
* Unique SMILES: 1727
* Species (MIC):  Acinetobacter baumannii(555)/Burkholderia cepacia(3)/Enterococcus faecalis(862)/Escherichia coli(5469)/Klebsiella pneumoniae(1021)/Neisseria gonorrhoeae(44)/Pseudomonas aeruginosa(3434)Staphylococcus aureus(278)/Streptococcus pneumoniae(1230)/Escherichia coli1(140)/Pseudomonas aeruginosa(25)
* Species (IC50): Escherichia coli(140)/Pseudomonas aeruginosa(25)

### SPARK Data Achaogen Contribution
Data on MIC, IC50, Protein Binding, Citotxicity (HUVEC), as well as data related to the LpxC target (MIC50, MIC90, PK). May be all the compounds are potential LpxC inhibitors. Relevant columns:
* 'SMILES'
* 'Curated & Transformed MIC Data: MIC (in µM) (µM)'
* 'Curated & Transformed MIC Data: Species'
* 'Curated & Transformed IC50 Data: IC50 (in uM) (uM)'
* 'Curated & Transformed IC50 Data: Species'
* 'Achaogen Protein Binding Data: Human Fraction Free (%)'
* 'Achaogen Protein Binding Data: Mouse Fraction Free (%)'
* 'Achaogen Protein Binding Data: Rat Fraction Free (%)'
* 'Achaogen Cytotox Data: IC50 (uM)'
* 'Achaogen LpxC MIC50/90 data: MIC50 (ug/mL)'
* 'Achaogen LpxC MIC50/90 data: MIC90 (ug/mL)'
Statistics:
* Total rows: 58749
* Unique SMILES: 1873
* Assays (unique SMILES): MIC (1861) / IC50(488) / Human PB (332) / Mouse (55) / Rat (316) / Cytotox (65) / LpxC MIC50 (80) / LpxC MIC90 (80) /
* Species MIC: Acinetobacter baumannii(1536) / Burkholderia thailandensis(726) / Escherichia coli(9356) / Klebsiella pneumoniae(3059) / Pseudomonas aeruginosa(28963) / Staphylococcus aureus(1713) / Yersinia enterocolitica(1406) 
* Species IC50:  Klebsiella pneumoniae(1) / Pseudomonas aeruginosa(535)

### SPARK Data CO-ADD Contribution
CoADD data, MIC and Single Point Inhibition. Relevant columns:
* 'SMILES'
* 'Curated & Transformed MIC Data: MIC (in µM) (µM)'
* 'Curated & Transformed MIC Data: Species'
* 'CO-ADD Single Concentration Data (Inhibition): Inhibition % (%)'
* 'CO-ADD Single Concentration Data (Inhibition): Species'
Statistics:
* Total rows: 82229
* Unique SMILES: 24090
* Species (MIC): Acinetobacter baumannii(1460) / Escherichia coli(1265) / Klebsiella pneumoniae(1265) / Pseudomonas aeruginosa(2448)
* Species (Single Point): Acinetobacter baumannii(23159) / Escherichia coli(5350) / Klebsiella pneumoniae(5349) / Pseudomonas aeruginosa(40969)

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

### SPARK MIC Data
Like the above file, it contains a Curated and an Extracted Section for MIC data. Instead of the actual MIC data, we sometimes have GIC and GIC % Inhibition (Growth Inhibition Concentration and % of inhibition at that concentration) Relevant columns:
* 'SMILES'
* 'Curated & Transformed MIC Data: MIC (in µM) (µM)'
* 'Curated & Transformed MIC Data: GIC (µg/mL)'
* 'Curated & Transformed MIC Data: GIC % inhibition'
* 'Curated & Transformed MIC Data: Species'
* 'Extracted & Uploaded MIC Data: MIC (in µM) (µM)'
* 'Extracted & Uploaded MIC Data: GIC (Growth Inhibition Concentration) (µg/mL)'
* 'Extracted & Uploaded MIC Data: GIC % Inhibition (%)'
* 'Curated & Transformed MIC Data: Species'
Statistics:
* Total rows: 137572
* Unique SMILES: 24581
* Species: more than 60

### SPARK Data Compounds & Physicochemical Properties
Information about Physicochemical properties of compounds, no bioactivity available
Statistics:
* Total rows: 162663
* Unique SMILES: 157021