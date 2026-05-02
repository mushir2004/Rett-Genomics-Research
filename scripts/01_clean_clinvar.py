import pandas as pd

def clean_clinvar_data():
    print(" Starting Data Cleaning Process...")

    # 1. Load the raw wall of text
    # Note: \t means tab-separated, which is what ClinVar gave us
    raw_file_path = "../data/raw/clinvar_mecp2_raw.txt"
    try:
        df = pd.read_csv(raw_file_path, sep='\t')
        print("\n ACTUAL COLUMNS IN FILE:\n", df.columns.tolist(), "\n")
        print(f" Loaded {len(df)} total variants from ClinVar.")
    except FileNotFoundError:
        print(" Error: Could not find the file. Make sure it's in data/raw/")
        return

    # 2. Filter for Rett Syndrome ONLY
    # We drop blank conditions and look for 'Rett' (case-insensitive)
    df = df.dropna(subset=['Condition(s)'])
    rett_mask = df['Condition(s)'].str.contains('Rett', case=False, na=False)
    df_rett = df[rett_mask]
    print(f" Filtered down to {len(df_rett)} variants linked to Rett Syndrome.")

    # 3. Separate the 'Knowns' (Gold Standard) from the 'Unknowns' (VUS)
    # Using the updated ClinVar column: 'Germline classification'
    pathogenic = df_rett[df_rett['Germline classification'].str.startswith('Pathogenic', na=False)]
    benign = df_rett[df_rett['Germline classification'].str.startswith('Benign', na=False)]
    
    # Combine them into our training set
    gold_standard = pd.concat([pathogenic, benign])
    
    print(f" Found {len(pathogenic)} Pathogenic and {len(benign)} Benign variants.")
    print(f" Gold Standard dataset size: {len(gold_standard)} variants.")

    # 4. Save the cleaned data to the processed folder
    output_path = "../data/processed/mecp2_gold_standard.csv"
    
    # Update the columns we keep to match the new headers
    columns_to_keep = ['Name', 'Gene(s)', 'Protein change', 'Germline classification', 'Germline review status']
    gold_standard[columns_to_keep].to_csv(output_path, index=False)
    
    print(f" Clean dataset saved successfully to {output_path}")

if __name__ == "__main__":
    clean_clinvar_data()