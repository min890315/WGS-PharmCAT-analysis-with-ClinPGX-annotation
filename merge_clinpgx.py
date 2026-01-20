
import pandas as pd
import os
import re
import glob
import sys

def parse_rs_diplotype(diplotype):
    """
    Parses 'rs9923231 reference (C)/rs9923231 variant (T)' format.
    Returns (rs_ids, alleles_string) or None if not matching format.
    
    Rules:
    1. Extract 'rs{digits}' -> match with Variant/Haplotypes.
    2. Extract () chars -> combine (e.g. CT).
    3. If multiple groups separated by OR or comma, combine alleles with ' + '.
    """
    # Check if this is an RS-style string
    if 'rs' not in diplotype:
        return None
    
    # Split by OR or comma if present, but be careful not to split inside the main structure if not intended.
    # User said: "OR" or "," for same form multiple times.
    # Example structure implies diplotype might be "RS_PART1 OR RS_PART2"
    
    parts = re.split(r'\s+OR\s+|\s*,\s*', diplotype)
    
    combined_alleles = []
    all_rs_ids = set()
    
    for part in parts:
        # Expected format: text with rsID and (X) ... / ... text with rsID and (Y)
        # We need to capture the letters in ().
        matches = re.findall(r'\((?P<base>[A-Za-z]+)\)', part)
        if not matches:
            return None
        
        # Form alleles like 'CT'
        # Usually size 2 for diplotype?
        allele_str = "".join(matches)
        combined_alleles.append(allele_str)
        
        # Extract rsIDs
        rs_matches = re.findall(r'rs\d+', part)
        all_rs_ids.update(rs_matches)
        
    final_allele_str = " + ".join(combined_alleles)
    return list(all_rs_ids), final_allele_str

def process_star_diplotype(diplotype):
    """
    Parses '*4/*6 OR *4/*35' format.
    Returns list of possible alleles to match.
    """
    # Split by OR
    parts = diplotype.split(' OR ')
    clean_parts = [p.strip() for p in parts]
    return clean_parts

def parse_g_diplotype(diplotype):
    """
    Parses 'Reference/g.38499696C>T' format.
    Returns 'CT'.
    """
    # Regex to capture the > notation
    # The user example is "Reference/g.numberC>T" -> "CT"
    # We look for g.<digits><Base1>><Base2>
    match = re.search(r'g\.\d+(?P<ref>[A-Za-z]+)>(?P<alt>[A-Za-z]+)', diplotype)
    if match:
        return match.group('ref') + match.group('alt')
    return None

def check_alleles_match(target, candidate):
    """
    Checks if candidate matches target, including reverse order for 2-char strings (e.g. CT == TC).
    """
    if target == candidate:
        return True
    
    # Requirement 2: Check swapped order
    # Assuming standard 2-char alleles like "CT"
    if len(target) == 2 and target == candidate[::-1]:
        return True
    
    return False

def main():
    # 1. Find report files
    if len(sys.argv) > 1:
        base_dir = sys.argv[1]
    else:
        # Fallback or default
        base_dir = '/home/jslinkspark/work/pharmcat_test'

    report_files = glob.glob(os.path.join(base_dir, '**', '*report.tsv'), recursive=True)
    if not report_files:
        print("No report.tsv files found.")
        return

    # 2. Load Annotation
    ann_file = os.path.join('/home/jslinkspark/work/pharmcat_test/forannot/var_drug_ann.tsv')
    if not os.path.exists(ann_file):
        print(f"Annotation file not found: {ann_file}")
        return
        
    ann_df = pd.read_csv(ann_file, sep='\t')
    
    # 3. Remove "Variant Annotation ID" column
    if "Variant Annotation ID" in ann_df.columns:
        ann_df = ann_df.drop(columns=["Variant Annotation ID"])
        
    # Process each report file
    for r_file in report_files:
        print(f"Processing {r_file}...")
        
        # Read report.tsv, skipping first row (header is on 2nd row)
        try:
            # Check file content first to ensure we aren't misinterpreting
            # pd.read_csv with header=1 means skip row 0, use row 1 as header.
            report_df = pd.read_csv(r_file, sep='\t', header=1)
        except Exception as e:
            print(f"Error reading {r_file}: {e}")
            continue

        merged_rows = []
        
        # Iterate over report rows
        for idx, row in report_df.iterrows():
            gene = row.get('Gene')
            source_diplo = row.get('Source Diplotype')
            
            if pd.isna(gene) or pd.isna(source_diplo):
                # If essential info missing, just keep the report row? 
                # User says "merging", usually assumes keeping the left row and adding nans if no match
                # Step 7 says "left join". So we keep report row.
                # But Step 4, 5, 6 describe filtering conditions for the RIGHT side matching.
                # So we essentially do a left join where the right side is filtered specifically for each left row.
                
                # We can replicate "left join" behavior by finding matches. 
                # If matches found, duplicate left row for each match. 
                # If no matches, keep left row with NaNs on right.
                
                merged_rows.append(row.to_dict()) # No merge data
                continue
                
            # Filter annotation by Gene first
            gene_matches = ann_df[ann_df['Gene'] == gene]
            
            if gene_matches.empty:
                merged_rows.append(row.to_dict())
                continue
            
            matched_ann_rows = pd.DataFrame()
            
            # Logic for matching
            # Check for 'rs' type or '*' type
            # Note: Some diplotypes might be simple text like "Reference/Reference". 
            # We need to handle that? User Rules 5 & 6 cover specific complex cases.
            # "Reference" usually implies wild type or star alleles like *1.
            # But usually PharmCAT output "Reference/Reference" for gene with no mutations?
            # Or "Reference/c.11266C>G".
            
            is_rs = 'rs' in str(source_diplo)
            is_star = '*' in str(source_diplo)
            
            # Check for g. notation e.g. Reference/g....
            # User condition: "Reference/g.38499696C>T" style string
            is_g = 'Reference/g.' in str(source_diplo) if not pd.isna(source_diplo) else False

            if is_g:
                # Rule 1 (New)
                g_allele = parse_g_diplotype(source_diplo)
                if g_allele:
                    # Filter gene_matches
                    matches = []
                    for _, ann_row in gene_matches.iterrows():
                        ann_alleles = str(ann_row.get('Alleles', ''))
                        if check_alleles_match(g_allele, ann_alleles):
                            matches.append(ann_row)
                    
                    if matches:
                        matched_ann_rows = pd.DataFrame(matches)

            elif is_rs:
                # Rule 5
                parsed = parse_rs_diplotype(source_diplo)
                if parsed:
                    rs_ids, allele_target = parsed
                    # Filter gene_matches
                    matches = []
                    for _, ann_row in gene_matches.iterrows():
                        ann_haplo = str(ann_row.get('Variant/Haplotypes', ''))
                        ann_alleles = str(ann_row.get('Alleles', ''))
                        
                        rs_match_found = False
                        for rs in rs_ids:
                            if rs in ann_haplo:
                                rs_match_found = True
                                break
                        
                        if rs_match_found:
                            # Rule 2 (New): Check match with swapping
                            if check_alleles_match(allele_target, ann_alleles):
                                matches.append(ann_row)
                    
                    if matches:
                        matched_ann_rows = pd.DataFrame(matches)
            
            elif is_star:
                # Rule 6
                candidates = process_star_diplotype(source_diplo)
                matches = []
                for _, ann_row in gene_matches.iterrows():
                    ann_alleles = str(ann_row.get('Alleles', ''))
                    # Check exact match against any candidate?
                    # Star alleles like *1/*1 usually don't need swapping of characters "*1" -> "1*", merging order *1/*2 vs *2/*1 might matter?
                    # Logic says: "split ... if contained in Alleles column"
                    # User didn't ask for swapping here explicitly, but it makes sense for *1/*2 vs *2/*1.
                    # However, rule 6 says "included in", implying the string from report is found in Alleles.
                    # Usually Alleles column has *1/*1.
                    if ann_alleles in candidates:
                         matches.append(ann_row)
                
                if matches:
                    matched_ann_rows = pd.DataFrame(matches)

            # Deduplication for matches (Rule 6 mentions it for star alleles, but good practice generally)
            if not matched_ann_rows.empty:
                # "Source Diplotype", "Variant/Haplotypes", "Alleles" used for deduplication?
                # The right side (ann) rows might be duplicates.
                # "Duplicate removal based on Source Diplotype, Variant/Haplotypes, Alleles"
                # But Source Diplotype is from Left. 
                # So basically for this specific report row, we want unique right-side entires defined by Var/Hap + Alleles.
                matched_ann_rows = matched_ann_rows.drop_duplicates(subset=['Variant/Haplotypes', 'Alleles'])
                
                # Perform the merge (Left Join expanded)
                for _, m_row in matched_ann_rows.iterrows():
                    # Combine row + m_row
                    combined = row.to_dict().copy()
                    for col, val in m_row.to_dict().items():
                        # Rule 8: Remove Gene from filtered ann (Right side).
                        # We already merged based on gene, so we don't need to add 'Gene' again.
                        # But pandas would suffix it.
                        if col == 'Gene':
                            continue
                        combined[col] = val
                    merged_rows.append(combined)
            else:
                # No match found in annotation, but it's a left join
                merged_rows.append(row.to_dict())

        # Construct DataFrame
        final_df = pd.DataFrame(merged_rows)
        
        # Output filename
        # "name is original report.tsv from .report.tsv to .ClinPGX_annot.tsv"
        new_name = r_file.replace('.report.tsv', '.ClinPGX_annot.tsv')
        
        final_df.to_csv(new_name, sep='\t', index=False)
        print(f"Saved to {new_name}")

if __name__ == "__main__":
    main()
