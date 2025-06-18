import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def load_abundance_data(file):
    """
    Loads the abundance file with taxonomic and sample data.
    """
    return pd.read_csv(file)

def load_metadata(file):
    """
    Loads the metadata file with sample information.
    Expected format: sample-id, group, condition, parameter
    """
    try:
        # Try to load as TSV first, then as CSV
        if file.endswith('.tsv'):
            metadata = pd.read_csv(file, sep='\t')
        else:
            metadata = pd.read_csv(file)
        
        # Check if it has the required columns
        required_cols = ['sample-id', 'group', 'condition', 'parameter']
        missing_cols = [col for col in required_cols if col not in metadata.columns]
        
        if missing_cols:
            raise ValueError(f"Required columns missing in metadata: {missing_cols}")
        
        return metadata
    
    except FileNotFoundError:
        print(f"❌ Error: Metadata file {file} not found!")
        print("Creating example file 'metadata_example.tsv'...")
        create_example_metadata()
        return None

def create_example_metadata():
    """
    Creates an example metadata file.
    """
    example_data = {
        'sample-id': ['F3A1', 'F3A2', 'F3A3', 'F3B1', 'F3B2', 'F3B3',
                      'F4A1', 'F4A2', 'F4A3', 'F4B1', 'F4B2', 'F4B3',
                      'F5A1', 'F5A2', 'F5A3', 'F5B1', 'F5B2', 'F5B3',
                      'F6A1', 'F6A2', 'F6A3', 'F6B1', 'F6B2', 'F6B3'],
        'group': ['F3']*6 + ['F4']*6 + ['F5']*6 + ['F6']*6,
        'condition': ['A', 'A', 'A', 'B', 'B', 'B'] * 4,
        'parameter': [3500]*6 + [9500]*6 + [7500]*6 + [11000]*6
    }
    
    example_df = pd.DataFrame(example_data)
    example_df.to_csv('metadata_example.tsv', sep='\t', index=False)
    print("✓ File 'metadata_example.tsv' created as template.")

def calculate_parameter_factors(metadata):
    """
    Calculates weighting factors based on the parameter provided in metadata.
    """
    # Get parameter value for each group
    group_parameters = {}
    for group in metadata['group'].unique():
        group_data = metadata[metadata['group'] == group]
        # Assumes parameter is the same for all samples in the group
        parameter_value = group_data['parameter'].iloc[0]
        group_parameters[group] = parameter_value
    
    mean_parameter = np.mean(list(group_parameters.values()))
    
    parameter_factors = {}
    print(f"Calculating weighting factors based on parameter:")
    for group, param_value in group_parameters.items():
        factor = param_value / mean_parameter
        parameter_factors[group] = factor
        print(f"{group}: {param_value} | Factor: {factor:.4f} | Weight: {((factor-1)*100):+.1f}%")
    
    return parameter_factors, mean_parameter, group_parameters

def get_sample_mapping(metadata):
    """
    Creates mapping between sample-id and groups/conditions.
    """
    sample_mapping = {}
    for _, row in metadata.iterrows():
        sample_id = row['sample-id']
        sample_mapping[sample_id] = {
            'group': row['group'],
            'condition': row['condition'],
            'parameter': row['parameter']
        }
    
    return sample_mapping

def calculate_species_variation(df, metadata):
    """
    Calculates mean variation (After - Before) for each species in each group.
    """
    taxonomy_cols = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    sample_cols = [col for col in df.columns if col not in taxonomy_cols]
    
    # Create sample mapping
    sample_mapping = get_sample_mapping(metadata)
    
    # Get unique groups
    groups = metadata['group'].unique()
    
    results = []
    
    print(f"Processing {len(df)} species across {len(groups)} groups...")
    
    for idx, (_, row) in enumerate(df.iterrows()):
        if idx > 0 and idx % 100 == 0:
            print(f"  Processed {idx} species...")
            
        species_name = row['Species']
        species_data = {
            'Species': species_name,
            'Phylum': row['Phylum'],
            'Class': row['Class'],
            'Order': row['Order'],
            'Family': row['Family'],
            'Genus': row['Genus']
        }
        
        for group in groups:
            # Find samples from this group in conditions A and B
            before_cols = []
            after_cols = []
            
            for sample_id in sample_cols:
                if sample_id in sample_mapping:
                    sample_info = sample_mapping[sample_id]
                    if sample_info['group'] == group:
                        if sample_info['condition'] == 'A':
                            before_cols.append(sample_id)
                        elif sample_info['condition'] == 'B':
                            after_cols.append(sample_id)
            
            if before_cols and after_cols:
                # Mean of samples before and after
                before_values = pd.to_numeric(row[before_cols], errors='coerce')
                after_values = pd.to_numeric(row[after_cols], errors='coerce')
                
                before_mean = before_values.mean()
                after_mean = after_values.mean()
                
                # Replace NaN with 0
                before_mean = 0 if pd.isna(before_mean) else before_mean
                after_mean = 0 if pd.isna(after_mean) else after_mean
                
                # Absolute and relative variation
                variation_abs = after_mean - before_mean
                
                # Relative variation (avoid division by zero)
                if before_mean > 0:
                    variation_rel = (after_mean - before_mean) / before_mean
                else:
                    # If before was 0 and after > 0, consider as appearance (100% increase)
                    variation_rel = 1.0 if after_mean > 0 else 0.0
                
                species_data[f'{group}_before_mean'] = before_mean
                species_data[f'{group}_after_mean'] = after_mean
                species_data[f'{group}_variation_abs'] = variation_abs
                species_data[f'{group}_variation_rel'] = variation_rel
                species_data[f'{group}_present'] = 1 if (before_mean > 0 or after_mean > 0) else 0
            else:
                # Species not present in this group
                species_data[f'{group}_before_mean'] = 0
                species_data[f'{group}_after_mean'] = 0
                species_data[f'{group}_variation_abs'] = 0
                species_data[f'{group}_variation_rel'] = 0
                species_data[f'{group}_present'] = 0
        
        results.append(species_data)
    
    print(f"✓ All {len(results)} species processed")
    return pd.DataFrame(results)

def calculate_parameter_weighted_score(variations_df, parameter_factors, groups):
    """
    Calculates the parameter-weighted score for each species.
    """
    results = []
    
    print(f"Calculating scores for {len(variations_df)} species...")
    
    for idx, (_, row) in enumerate(variations_df.iterrows()):
        if idx > 0 and idx % 100 == 0:
            print(f"  Scores calculated for {idx} species...")
            
        species_scores = []
        species_weights = []
        groups_with_data = []
        
        for group in groups:
            # Check if species has data in this group
            if row[f'{group}_present'] == 1:
                variation = row[f'{group}_variation_rel']
                parameter_factor = parameter_factors[group]
                
                # Weighted score: variation × parameter factor
                weighted_score = variation * parameter_factor
                species_scores.append(weighted_score)
                species_weights.append(parameter_factor)
                groups_with_data.append(group)
        
        # Process the species
        if species_scores:
            # Final score: weighted average of scores from groups where species is present
            final_score = np.average(species_scores, weights=species_weights)
            
            # Additional statistics
            n_groups_present = len(species_scores)
            max_score = max(species_scores)
            min_score = min(species_scores)
            std_score = np.std(species_scores) if len(species_scores) > 1 else 0
            
        else:
            # Species with no data in any group
            final_score = 0
            n_groups_present = 0
            max_score = 0
            min_score = 0
            std_score = 0
        
        result = {
            'Species': row['Species'],
            'Phylum': row['Phylum'],
            'Class': row['Class'],
            'Order': row['Order'],
            'Family': row['Family'],
            'Genus': row['Genus'],
            'Parameter_Weighted_Score': final_score,
            'N_Groups_Present': n_groups_present,
            'Max_Group_Score': max_score,
            'Min_Group_Score': min_score,
            'Score_Std': std_score,
            'Groups_With_Data': ', '.join(groups_with_data) if groups_with_data else 'None'
        }
        
        # Add individual scores per group
        for group in groups:
            if row[f'{group}_present'] == 1:
                variation = row[f'{group}_variation_rel']
                weighted_score = variation * parameter_factors[group]
                result[f'{group}_Score'] = weighted_score
                result[f'{group}_Variation_Rel'] = variation
                result[f'{group}_Before_Mean'] = row[f'{group}_before_mean']
                result[f'{group}_After_Mean'] = row[f'{group}_after_mean']
            else:
                result[f'{group}_Score'] = 0
                result[f'{group}_Variation_Rel'] = 0
                result[f'{group}_Before_Mean'] = 0
                result[f'{group}_After_Mean'] = 0
        
        results.append(result)
    
    print(f"✓ Scores calculated for all {len(results)} species")
    return pd.DataFrame(results)

def interpret_results(results_df):
    """
    Interprets results and adds classifications.
    """
    results_df = results_df.copy()
    
    # Classify parameter impact
    def classify_impact(score):
        if score > 0.1:
            return "High Positive Impact"
        elif score > 0.05:
            return "Moderate Positive Impact"
        elif score > 0:
            return "Low Positive Impact"
        elif score > -0.05:
            return "Low Negative Impact"
        elif score > -0.1:
            return "Moderate Negative Impact"
        else:
            return "High Negative Impact"
    
    results_df['Impact_Classification'] = results_df['Parameter_Weighted_Score'].apply(classify_impact)
    
    # Sort by score (most impactful first)
    results_df = results_df.sort_values('Parameter_Weighted_Score', ascending=False)
    
    return results_df

def generate_summary_stats(results_df, group_parameters, output_summary=True):
    """
    Generates summary statistics and saves to file.
    """
    summary_lines = []
    
    summary_lines.append("=" * 80)
    summary_lines.append("STATISTICAL SUMMARY OF ANALYSIS")
    summary_lines.append("=" * 80)
    
    # Group parameters information
    summary_lines.append("\nGROUP PARAMETERS:")
    for group, param_value in group_parameters.items():
        summary_lines.append(f"{group}: {param_value}")
    
    total_species = len(results_df)
    positive_impact = len(results_df[results_df['Parameter_Weighted_Score'] > 0])
    negative_impact = len(results_df[results_df['Parameter_Weighted_Score'] < 0])
    
    summary_lines.append(f"\nTotal species analyzed: {total_species}")
    summary_lines.append(f"Species with positive impact: {positive_impact} ({positive_impact/total_species*100:.1f}%)")
    summary_lines.append(f"Species with negative impact: {negative_impact} ({negative_impact/total_species*100:.1f}%)")
    
    summary_lines.append(f"\nTOP 30 SPECIES WITH HIGHEST POSITIVE IMPACT:")
    top_positive = results_df[results_df['Parameter_Weighted_Score'] > 0].head(30)
    if len(top_positive) > 0:
        for i, (_, row) in enumerate(top_positive.iterrows(), 1):
            summary_lines.append(f"{i:2d}. {row['Species'][:50]:<50} | Score: {row['Parameter_Weighted_Score']:+.4f}")
    else:
        summary_lines.append("No species with positive impact found.")
    
    summary_lines.append(f"\nTOP 30 SPECIES WITH HIGHEST NEGATIVE IMPACT:")
    negative_species = results_df[results_df['Parameter_Weighted_Score'] < 0]
    if len(negative_species) > 0:
        top_negative = negative_species.tail(30).iloc[::-1]  # Reverse to show most negative first
        for i, (_, row) in enumerate(top_negative.iterrows(), 1):
            summary_lines.append(f"{i:2d}. {row['Species'][:50]:<50} | Score: {row['Parameter_Weighted_Score']:+.4f}")
    else:
        summary_lines.append("No species with negative impact found.")
    
    summary_lines.append(f"\nSPECIES WITH NO DATA (Score = 0):")
    zero_species = results_df[results_df['Parameter_Weighted_Score'] == 0]
    summary_lines.append(f"Total: {len(zero_species)} species")
    
    # Save to file
    if output_summary:
        with open('top30.txt', 'w', encoding='utf-8') as f:
            f.write('\n'.join(summary_lines))
        print("✓ TOP 30 species summary saved in: top30.txt")
    
    # Display on screen
    for line in summary_lines:
        print(line)

def main():
    print("Starting analysis with external parameters...")
    print("=" * 80)
    
    # 1. Load abundance data
    input_file = 'abundance_data.csv'
    try:
        df = load_abundance_data(input_file)
        print(f"✓ Abundance data loaded: {len(df)} species, {len(df.columns)-6} samples")
    except FileNotFoundError:
        print(f"❌ Error: File {input_file} not found!")
        return
    
    # 2. Load metadata
    metadata_files = ['metadata.tsv', 'metadata.txt', 'metadata.csv']
    metadata = None
    
    for metadata_file in metadata_files:
        try:
            metadata = load_metadata(metadata_file)
            if metadata is not None:
                print(f"✓ Metadata loaded from: {metadata_file}")
                break
        except:
            continue
    
    if metadata is None:
        print("❌ No valid metadata file found!")
        print("Use 'metadata_example.tsv' file as template.")
        return
    
    # 3. Calculate parameter factors
    print(f"\nCalculating weighting factors:")
    parameter_factors, mean_param, group_parameters = calculate_parameter_factors(metadata)
    print(f"Mean parameter value: {mean_param:.2f}")
    
    # 4. Check correspondence between samples
    sample_cols = [col for col in df.columns if col not in ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]
    metadata_samples = set(metadata['sample-id'].tolist())
    df_samples = set(sample_cols)
    
    missing_in_metadata = df_samples - metadata_samples
    missing_in_df = metadata_samples - df_samples
    
    if missing_in_metadata:
        print(f"\n⚠️  Samples in abundance_data.csv but not in metadata: {missing_in_metadata}")
    if missing_in_df:
        print(f"⚠️  Samples in metadata but not in abundance_data.csv: {missing_in_df}")
    
    # 5. Calculate species variations  
    print(f"\nCalculating abundance variations per species and group...")
    groups = metadata['group'].unique()
    variations_df = calculate_species_variation(df, metadata)
    print(f"✓ Variations calculated for {len(variations_df)} species across {len(groups)} groups")
    
    # 6. Calculate weighted scores
    print(f"\nCalculating parameter-weighted scores...")
    results_df = calculate_parameter_weighted_score(variations_df, parameter_factors, groups)
    print(f"✓ Scores calculated for {len(results_df)} species")
    
    # 7. Interpret results
    results_df = interpret_results(results_df)
    
    # 8. Save COMPLETE results
    output_file = 'parameter_weighted_analysis.csv'
    results_df.to_csv(output_file, index=False)
    print(f"✓ COMPLETE results saved in: {output_file}")
    print(f"✓ All {len(results_df)} species included in CSV file")
    
    # 9. Generate statistical summary and TOP 30 file
    generate_summary_stats(results_df, group_parameters)
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETED!")
    print(f"Output files:")
    print(f"  - {output_file} (complete results)")
    print(f"  - top30.txt (top 30 species summary)")
    print("="*80)

if __name__ == '__main__':
    main()
