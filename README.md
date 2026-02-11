# Genomic Landscape of Pharmacogenomic Variants: Insights from WGS-based Annotation using the ClinPGX Database

This study implemented an integrated pharmacogenomic (PGx) pipeline using PharmCAT and ClinPGX to interpret Whole Genome Sequencing (WGS) data from 15 individuals. By annotating diplotypes with literature-based phenotypes, we identified critical actionable variants. Significant findings included the VKORC1 rs9923231 TT genotype in all samples, necessitating reduced warfarin dosage due to high sensitivity. In pediatric contexts, NUDT15 intermediate metabolizers (*1/*2, *1/*3) were flagged for reduced mercaptopurine dosage to mitigate toxicity risks. Additionally, CYP2C19 variants (*1/*3) were associated with clopidogrel resistance, and CYP3A5 poor metabolizers (*3/*3) with decreased tacrolimus clearance. This workflow effectively translates WGS data into precise clinical recommendations, demonstrating its utility for personalized prescribing.

# Sample
# python3 ~/work/codes/run_pharmcat_stargazer.py <VCF Path> <WorkDir> <Output Prefix>
python3 ~/work/codes/run_pharmcat_stargazer.py ../Be-5.g.vcf.gz $PWD test
