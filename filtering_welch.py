import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def load_data():
    """
    Loads all necessary files for analysis.
    """
    try:
        # Load abundance data
        abundance_data = pd.read_csv('abundance_data.csv')
        print(f"✓ abundance_data.csv loaded: {len(abundance_data)} species")
        
        # Load metadata
        metadata_files = ['metadata.tsv', 'metadata.txt', 'metadata.csv']
        metadata = None
        
        for metadata_file in metadata_files:
            try:
                if metadata_file.endswith('.tsv'):
                    metadata = pd.read_csv(metadata_file, sep='\t')
                else:
                    metadata = pd.read_csv(metadata_file)
                print(f"✓ {metadata_file} loaded: {len(metadata)} samples")
                break
            except FileNotFoundError:
                continue
        
        if metadata is None:
            raise FileNotFoundError("No metadata file found")
        
        # Load results from previous analysis
        parameter_results = pd.read_csv('parameter_weighted_analysis.csv')
        print(f"✓ parameter_weighted_analysis.csv loaded: {len(parameter_results)} species")
        
        # Load top 30 species list
        top30_species = extract_top30_species('top30.txt')
        print(f"✓ top30.txt processed: {len(top30_species)} species extracted")
        
        return abundance_data, metadata, parameter_results, top30_species
        
    except FileNotFoundError as e:
        print(f"❌ Error loading files: {e}")
        return None, None, None, None

def extract_top30_species(top30_file):
    """
    Extracts species names from top30.txt file.
    """
    try:
        with open(top30_file, 'r', encoding='utf-8') as f:
            content = f.read()
        
        species_names = set()
        lines = content.split('\n')
        
        for line in lines:
            # Look for lines that start with numbers (rankings)
            if line.strip() and (line.strip()[0].isdigit() or line.strip().startswith(' ')):
                # Try to extract species name (after number and before |)
                if '|' in line:
                    species_part = line.split('|')[0].strip()
                    # Remove ranking number
                    if '. ' in species_part:
                        species_name = species_part.split('. ', 1)[1].strip()
                        species_names.add(species_name)
        
        return list(species_names)
        
    except FileNotFoundError:
        print("⚠️ top30.txt file not found. Continuing without TOP30 filter.")
        return []

def get_sample_mapping(metadata):
    """
    Creates mapping between sample-id and groups/conditions.
    """
    sample_mapping = {}
    for _, row in metadata.iterrows():
        sample_id = row['sample-id']
        sample_mapping[sample_id] = {
            'group': row['group'],
            'condition': row['condition']
        }
    
    return sample_mapping

def apply_welch_test(data1, data2):
    """
    Applies Welch's t-test (unequal variances).
    """
    if len(data1) < 2 or len(data2) < 2:
        return 1.0
    
    try:
        t_stat, p_value = stats.ttest_ind(data1, data2, equal_var=False)
        return p_value
    except:
        return 1.0

def calculate_welch_stats(abundance_data, metadata):
    """
    Calculates statistics using Welch's t-test for all comparisons.
    """
    taxonomy_cols = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    sample_cols = [col for col in abundance_data.columns if col not in taxonomy_cols]
    
    # Sample mapping
    sample_mapping = get_sample_mapping(metadata)
    
    # Separate samples by condition
    condition_a_samples = []
    condition_b_samples = []
    
    for sample_id in sample_cols:
        if sample_id in sample_mapping:
            if sample_mapping[sample_id]['condition'] == 'A':
                condition_a_samples.append(sample_id)
            elif sample_mapping[sample_id]['condition'] == 'B':
                condition_b_samples.append(sample_id)
    
    # Get unique groups
    groups = metadata['group'].unique()
    
    results = []
    
    print(f"Calculating Welch's t-test statistics for {len(abundance_data)} species...")
    print("Using Welch's t-test (unequal variances) for all comparisons")
    
    for idx, (_, row) in enumerate(abundance_data.iterrows()):
        if idx > 0 and idx % 100 == 0:
            print(f"  Processed {idx} species...")
        
        species_name = row['Species']
        
        species_data = {
            'Phylum': row['Phylum'],
            'Class': row['Class'],
            'Order': row['Order'],
            'Family': row['Family'],
            'Genus': row['Genus'],
            'Species': species_name
        }
        
        # Add sample values
        for sample_id in sample_cols:
            if sample_id in abundance_data.columns:
                species_data[sample_id] = row[sample_id]
        
        # Calculate global statistics
        all_values = pd.to_numeric(row[sample_cols], errors='coerce').fillna(0)
        condition_a_values = pd.to_numeric(row[condition_a_samples], errors='coerce').fillna(0)
        condition_b_values = pd.to_numeric(row[condition_b_samples], errors='coerce').fillna(0)
        
        # Means and standard deviations
        species_data['MEAN_TOTAL'] = all_values.mean()
        species_data['Mean_before'] = condition_a_values.mean()
        species_data['Mean_after'] = condition_b_values.mean()
        species_data['SD_bef'] = condition_a_values.std()
        species_data['SD_after'] = condition_b_values.std()
        
        # Global Welch's t-test
        global_p = apply_welch_test(condition_a_values, condition_b_values)
        species_data['Global-t-test'] = global_p
        
        # Tests by group
        significant_tests = 0
        
        for group in groups:
            # Samples from this group
            group_a_samples = [s for s in sample_cols if s in sample_mapping and 
                             sample_mapping[s]['group'] == group and sample_mapping[s]['condition'] == 'A']
            group_b_samples = [s for s in sample_cols if s in sample_mapping and 
                             sample_mapping[s]['group'] == group and sample_mapping[s]['condition'] == 'B']
            
            if group_a_samples and group_b_samples:
                group_a_values = pd.to_numeric(row[group_a_samples], errors='coerce').fillna(0)
                group_b_values = pd.to_numeric(row[group_b_samples], errors='coerce').fillna(0)
                
                # Welch's t-test for this group
                group_p = apply_welch_test(group_a_values, group_b_values)
                species_data[f'{group}t-test'] = group_p
                
                if group_p < 0.05:
                    significant_tests += 1
                
            else:
                species_data[f'{group}t-test'] = 1.0
        
        # Determine color based on Welch's tests
        global_significant = species_data['Global-t-test'] < 0.05
        
        if global_significant:
            species_data['Cor'] = 'vermelho'
        elif significant_tests > 1:
            species_data['Cor'] = 'laranja'
        elif significant_tests == 1:
            species_data['Cor'] = 'amarelo'
        else:
            species_data['Cor'] = 'cinza'
        
        results.append(species_data)
    
    print(f"✓ Welch's t-test statistics calculated for all {len(results)} species")
    return pd.DataFrame(results)

def apply_filters(stats_df, parameter_results, top30_species):
    """
    Applies filters to generate filtrado1.csv file.
    """
    print("Applying filters with Welch's t-test results...")
    
    # Merge with parameter weighted score results
    if 'Parameter_Weighted_Score' not in stats_df.columns:
        merge_cols = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        stats_df = stats_df.merge(
            parameter_results[merge_cols + ['Parameter_Weighted_Score']],
            on=merge_cols,
            how='left'
        )
        stats_df['Parameter_Weighted_Score'] = stats_df['Parameter_Weighted_Score'].fillna(0)
    
    print(f"Total species before filters: {len(stats_df)}")
    
    # Filter 1: Total mean >= 0.005 (0.5%)
    filter1 = stats_df['MEAN_TOTAL'] >= 0.005
    print(f"Species with total mean >= 0.5%: {filter1.sum()}")
    
    # Filter 2: At least one significant Welch's test OR in TOP30
    t_test_cols = [col for col in stats_df.columns if col.endswith('t-test')]
    significant_t_tests = (stats_df[t_test_cols] < 0.05).any(axis=1)
    in_top30 = stats_df['Species'].isin(top30_species)
    
    filter2 = significant_t_tests | in_top30
    print(f"Species with significant Welch's t-test: {significant_t_tests.sum()}")
    print(f"Species in TOP30: {in_top30.sum()}")
    print(f"Species meeting statistical/TOP30 criteria: {filter2.sum()}")
    
    # Apply both filters
    final_filter = filter1 & filter2
    filtered_df = stats_df[final_filter].copy()
    
    print(f"✓ Species after applying all filters: {len(filtered_df)}")
    
    # Order columns in desired sequence
    taxonomy_cols = ['Cor', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Parameter_Weighted_Score']
    
    # Identify t-test columns
    t_test_cols = [col for col in filtered_df.columns if col.endswith('t-test')]
    
    sample_cols = [col for col in filtered_df.columns if col not in taxonomy_cols + 
                   ['MEAN_TOTAL', 'Mean_before', 'Mean_after', 'SD_bef', 'SD_after'] + t_test_cols]
    stats_cols = ['MEAN_TOTAL', 'Mean_before', 'Mean_after', 'SD_bef', 'SD_after']
    
    # Final column order
    column_order = taxonomy_cols + sample_cols + t_test_cols + stats_cols
    
    # Reorder columns (only existing ones)
    existing_columns = [col for col in column_order if col in filtered_df.columns]
    filtered_df = filtered_df[existing_columns]
    
    return filtered_df

def generate_summary_report(filtered_df, original_count):
    """
    Generates summary report of Welch's t-test filtering.
    """
    report_lines = []
    
    report_lines.append("=" * 80)
    report_lines.append("WELCH'S T-TEST FILTERING ANALYSIS REPORT")
    report_lines.append("=" * 80)
    
    report_lines.append(f"Species in original dataset: {original_count}")
    report_lines.append(f"Species after Welch's t-test filtering: {len(filtered_df)}")
    report_lines.append(f"Percentage retained: {len(filtered_df)/original_count*100:.1f}%")
    
    report_lines.append(f"\nStatistical test used: Welch's t-test (unequal variances)")
    report_lines.append(f"Test selection rationale: Robust for unequal variances, no normality assumption")
    
    # Statistics by color
    color_counts = filtered_df['Cor'].value_counts()
    report_lines.append(f"\nDistribution by significance pattern:")
    color_meanings = {
        'vermelho': 'Global significance (Red)',
        'laranja': 'Multiple groups significant (Orange)', 
        'amarelo': 'Single group significant (Yellow)',
        'cinza': 'Not significant (Gray)'
    }
    
    for color, meaning in color_meanings.items():
        count = color_counts.get(color, 0)
        report_lines.append(f"  {meaning}: {count} species")
    
    # Top species by significance
    global_sig = filtered_df[filtered_df['Global-t-test'] < 0.05].sort_values('Parameter_Weighted_Score', ascending=False)
    report_lines.append(f"\nTop 10 globally significant species (Welch's test):")
    for i, (_, row) in enumerate(global_sig.head(10).iterrows(), 1):
        report_lines.append(f"  {i:2d}. {row['Species'][:50]:<50} | PWS: {row['Parameter_Weighted_Score']:+.4f} | p: {row['Global-t-test']:.4f}")
    
    # Statistics summary
    report_lines.append(f"\nFiltering criteria applied:")
    report_lines.append(f"  1. Mean abundance >= 0.5% across all samples")
    report_lines.append(f"  2. At least one significant Welch's t-test (p < 0.05) OR in TOP30 ranking")
    
    mean_abundance = filtered_df['MEAN_TOTAL'].mean()
    report_lines.append(f"\nFiltered species statistics:")
    report_lines.append(f"  Mean abundance: {mean_abundance*100:.3f}%")
    report_lines.append(f"  Abundance range: {filtered_df['MEAN_TOTAL'].min()*100:.3f}% - {filtered_df['MEAN_TOTAL'].max()*100:.3f}%")
    
    # Save report
    with open('welch_filtering_report.txt', 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    # Display on screen
    for line in report_lines:
        print(line)

def main():
    print("Starting Welch's t-test filtering analysis...")
    print("=" * 80)
    
    # 1. Load all necessary data
    abundance_data, metadata, parameter_results, top30_species = load_data()
    
    if any(data is None for data in [abundance_data, metadata, parameter_results]):
        print("❌ Error: Could not load all necessary files!")
        return
    
    original_count = len(abundance_data)
    
    # 2. Calculate statistics using Welch's t-test
    print(f"\nCalculating statistics using Welch's t-test...")
    stats_df = calculate_welch_stats(abundance_data, metadata)
    
    # 3. Apply filters
    print(f"\nApplying filters...")
    filtered_df = apply_filters(stats_df, parameter_results, top30_species)
    
    # 4. Save filtered file
    output_file = 'filtrado1.csv'
    filtered_df.to_csv(output_file, index=False)
    print(f"✓ Filtered file saved: {output_file}")
    
    # 5. Generate summary report
    generate_summary_report(filtered_df, original_count)
    print(f"✓ Welch's t-test report saved: welch_filtering_report.txt")
    
    print(f"\n{'='*80}")
    print("WELCH'S T-TEST ANALYSIS COMPLETED!")
    print(f"Generated files:")
    print(f"  - {output_file} ({len(filtered_df)} species - Welch's t-test)")
    print(f"  - welch_filtering_report.txt (detailed analysis report)")
    print("="*80)

if __name__ == '__main__':
    main()
