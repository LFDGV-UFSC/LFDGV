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

def test_normality(data):
    """
    Tests normality using Shapiro-Wilk test.
    Returns p-value of the test.
    """
    if len(data) < 3:
        return np.nan  # Shapiro-Wilk needs at least 3 observations
    
    # Remove zeros and very small values that may cause problems
    data_clean = data[data > 1e-10]
    
    if len(data_clean) < 3:
        return np.nan
    
    try:
        stat, p_value = stats.shapiro(data_clean)
        return p_value
    except:
        return np.nan

def test_equal_variances(data1, data2):
    """
    Tests equality of variances using Levene's test.
    Returns p-value of the test.
    """
    if len(data1) < 2 or len(data2) < 2:
        return np.nan
    
    try:
        stat, p_value = stats.levene(data1, data2)
        return p_value
    except:
        return np.nan

def choose_best_test(data1, data2):
    """
    Chooses the best statistical test based on:
    1. Normality (Shapiro-Wilk)
    2. Equality of variances (Levene)
    
    Returns: test_name, p_value, details
    """
    
    if len(data1) < 2 or len(data2) < 2:
        return "Insufficient_data", 1.0, {"reason": "Less than 2 observations per group"}
    
    # Test normality
    norm_p1 = test_normality(data1)
    norm_p2 = test_normality(data2)
    
    # Consider normal if p > 0.05 (or cannot test)
    is_normal1 = pd.isna(norm_p1) or norm_p1 > 0.05
    is_normal2 = pd.isna(norm_p2) or norm_p2 > 0.05
    
    both_normal = is_normal1 and is_normal2
    
    details = {
        "shapiro_p1": norm_p1,
        "shapiro_p2": norm_p2,
        "normal1": is_normal1,
        "normal2": is_normal2,
        "both_normal": both_normal
    }
    
    if both_normal:
        # Test equality of variances
        levene_p = test_equal_variances(data1, data2)
        equal_vars = pd.isna(levene_p) or levene_p > 0.05
        
        details["levene_p"] = levene_p
        details["equal_vars"] = equal_vars
        
        if equal_vars:
            # Student's t-test (equal variances)
            try:
                t_stat, p_value = stats.ttest_ind(data1, data2, equal_var=True)
                return "Student_t_test", p_value, details
            except:
                return "Test_failed", 1.0, details
        else:
            # Welch's t-test (unequal variances)
            try:
                t_stat, p_value = stats.ttest_ind(data1, data2, equal_var=False)
                return "Welch_t_test", p_value, details
            except:
                return "Test_failed", 1.0, details
    else:
        # Non-parametric test (Mann-Whitney U)
        try:
            stat, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
            return "Mann_Whitney_U", p_value, details
        except:
            return "Test_failed", 1.0, details

def apply_conservative_test(data1, data2, test_type):
    """
    Applies the specified conservative test.
    """
    if len(data1) < 2 or len(data2) < 2:
        return 1.0
    
    try:
        if test_type == "Student_t_test":
            t_stat, p_value = stats.ttest_ind(data1, data2, equal_var=True)
        elif test_type == "Welch_t_test":
            t_stat, p_value = stats.ttest_ind(data1, data2, equal_var=False)
        elif test_type == "Mann_Whitney_U":
            stat, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
        else:
            return 1.0
        
        return p_value
    except:
        return 1.0

def determine_conservative_test(abundance_data, metadata):
    """
    Determines which conservative test to use based on analysis of all species.
    """
    print("Determining conservative test based on global analysis...")
    
    taxonomy_cols = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    sample_cols = [col for col in abundance_data.columns if col not in taxonomy_cols]
    sample_mapping = get_sample_mapping(metadata)
    
    # Separate samples by condition
    condition_a_samples = [s for s in sample_cols if s in sample_mapping and sample_mapping[s]['condition'] == 'A']
    condition_b_samples = [s for s in sample_cols if s in sample_mapping and sample_mapping[s]['condition'] == 'B']
    groups = metadata['group'].unique()
    
    # Counters for decision
    any_non_normal = False
    any_unequal_vars = False
    
    # Sample a subset of species for analysis (to save time)
    sample_size = min(100, len(abundance_data))  # Analyze up to 100 species to decide
    sample_indices = np.random.choice(len(abundance_data), sample_size, replace=False)
    
    print(f"Analyzing {sample_size} species to determine conservative test...")
    
    for idx in sample_indices:
        row = abundance_data.iloc[idx]
        
        # Global analysis
        condition_a_values = pd.to_numeric(row[condition_a_samples], errors='coerce').fillna(0)
        condition_b_values = pd.to_numeric(row[condition_b_samples], errors='coerce').fillna(0)
        
        # Test global normality
        norm_p1 = test_normality(condition_a_values)
        norm_p2 = test_normality(condition_b_values)
        
        is_normal1 = pd.isna(norm_p1) or norm_p1 > 0.05
        is_normal2 = pd.isna(norm_p2) or norm_p2 > 0.05
        
        if not (is_normal1 and is_normal2):
            any_non_normal = True
        
        # If both normal, test variances
        if is_normal1 and is_normal2:
            levene_p = test_equal_variances(condition_a_values, condition_b_values)
            if not (pd.isna(levene_p) or levene_p > 0.05):
                any_unequal_vars = True
        
        # Group analysis
        for group in groups:
            group_a_samples = [s for s in sample_cols if s in sample_mapping and 
                             sample_mapping[s]['group'] == group and sample_mapping[s]['condition'] == 'A']
            group_b_samples = [s for s in sample_cols if s in sample_mapping and 
                             sample_mapping[s]['group'] == group and sample_mapping[s]['condition'] == 'B']
            
            if group_a_samples and group_b_samples:
                group_a_values = pd.to_numeric(row[group_a_samples], errors='coerce').fillna(0)
                group_b_values = pd.to_numeric(row[group_b_samples], errors='coerce').fillna(0)
                
                # Test group normality
                norm_p1 = test_normality(group_a_values)
                norm_p2 = test_normality(group_b_values)
                
                is_normal1 = pd.isna(norm_p1) or norm_p1 > 0.05
                is_normal2 = pd.isna(norm_p2) or norm_p2 > 0.05
                
                if not (is_normal1 and is_normal2):
                    any_non_normal = True
                
                # If both normal, test variances
                if is_normal1 and is_normal2:
                    levene_p = test_equal_variances(group_a_values, group_b_values)
                    if not (pd.isna(levene_p) or levene_p > 0.05):
                        any_unequal_vars = True
    
    # Decide conservative test
    if any_non_normal:
        conservative_test = "Mann_Whitney_U"
        reason = "Non-normal data detected in some species"
    elif any_unequal_vars:
        conservative_test = "Welch_t_test"
        reason = "Unequal variances detected in some species"
    else:
        conservative_test = "Student_t_test"
        reason = "Normal data with equal variances in all analyzed species"
    
    print(f"✓ Conservative test chosen: {conservative_test}")
    print(f"  Reason: {reason}")
    
    return conservative_test

def calculate_comprehensive_stats(abundance_data, metadata):
    """
    Calculates comprehensive statistics with individual and conservative tests.
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
    
    # Determine conservative test
    conservative_test = determine_conservative_test(abundance_data, metadata)
    
    results = []
    sup_stats = []
    
    print(f"\nCalculating comprehensive statistics for {len(abundance_data)} species...")
    
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
        
        # Data for sup_statistics
        sup_data = {
            'Species': species_name,
            'Phylum': row['Phylum'],
            'Class': row['Class'],
            'Order': row['Order'],
            'Family': row['Family'],
            'Genus': row['Genus'],
            'Conservative_Test_Type': conservative_test
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
        
        # Add to sup_statistics
        sup_data['Global_Mean_A'] = condition_a_values.mean()
        sup_data['Global_SD_A'] = condition_a_values.std()
        sup_data['Global_Mean_B'] = condition_b_values.mean()
        sup_data['Global_SD_B'] = condition_b_values.std()
        sup_data['Global_Mean_Total'] = all_values.mean()
        
        # Global test - INDIVIDUAL (best test chosen)
        test_name, p_value, details = choose_best_test(condition_a_values, condition_b_values)
        species_data['Global-t-test'] = p_value
        
        # Global test - CONSERVATIVE
        conservative_p = apply_conservative_test(condition_a_values, condition_b_values, conservative_test)
        species_data['Global-t-test-conservative'] = conservative_p
        
        # Add details to sup_statistics
        sup_data['Global_Best_Test_Type'] = test_name
        sup_data['Global_Best_P_Value'] = p_value
        sup_data['Global_Conservative_P_Value'] = conservative_p
        sup_data['Global_Shapiro_A'] = details.get('shapiro_p1', np.nan)
        sup_data['Global_Shapiro_B'] = details.get('shapiro_p2', np.nan)
        sup_data['Global_Normal_A'] = details.get('normal1', False)
        sup_data['Global_Normal_B'] = details.get('normal2', False)
        sup_data['Global_Levene_P'] = details.get('levene_p', np.nan)
        sup_data['Global_Equal_Vars'] = details.get('equal_vars', False)
        
        # Tests by group
        significant_tests_individual = 0
        significant_tests_conservative = 0
        
        for group in groups:
            # Samples from this group
            group_a_samples = [s for s in sample_cols if s in sample_mapping and 
                             sample_mapping[s]['group'] == group and sample_mapping[s]['condition'] == 'A']
            group_b_samples = [s for s in sample_cols if s in sample_mapping and 
                             sample_mapping[s]['group'] == group and sample_mapping[s]['condition'] == 'B']
            
            if group_a_samples and group_b_samples:
                group_a_values = pd.to_numeric(row[group_a_samples], errors='coerce').fillna(0)
                group_b_values = pd.to_numeric(row[group_b_samples], errors='coerce').fillna(0)
                
                # Individual test (best for this species/group)
                test_name, p_value, details = choose_best_test(group_a_values, group_b_values)
                species_data[f'{group}t-test'] = p_value
                if p_value < 0.05:
                    significant_tests_individual += 1
                
                # Conservative test
                conservative_p = apply_conservative_test(group_a_values, group_b_values, conservative_test)
                species_data[f'{group}t-test-conservative'] = conservative_p
                if conservative_p < 0.05:
                    significant_tests_conservative += 1
                
                # Add to sup_statistics
                sup_data[f'{group}_Mean_A'] = group_a_values.mean()
                sup_data[f'{group}_SD_A'] = group_a_values.std()
                sup_data[f'{group}_Mean_B'] = group_b_values.mean()
                sup_data[f'{group}_SD_B'] = group_b_values.std()
                sup_data[f'{group}_Best_Test_Type'] = test_name
                sup_data[f'{group}_Best_P_Value'] = p_value
                sup_data[f'{group}_Conservative_P_Value'] = conservative_p
                sup_data[f'{group}_Shapiro_A'] = details.get('shapiro_p1', np.nan)
                sup_data[f'{group}_Shapiro_B'] = details.get('shapiro_p2', np.nan)
                sup_data[f'{group}_Normal_A'] = details.get('normal1', False)
                sup_data[f'{group}_Normal_B'] = details.get('normal2', False)
                sup_data[f'{group}_Levene_P'] = details.get('levene_p', np.nan)
                sup_data[f'{group}_Equal_Vars'] = details.get('equal_vars', False)
                
            else:
                species_data[f'{group}t-test'] = 1.0
                species_data[f'{group}t-test-conservative'] = 1.0
                
                # Fill sup_statistics with default values
                sup_data[f'{group}_Mean_A'] = 0
                sup_data[f'{group}_SD_A'] = 0
                sup_data[f'{group}_Mean_B'] = 0
                sup_data[f'{group}_SD_B'] = 0
                sup_data[f'{group}_Best_Test_Type'] = "No_data"
                sup_data[f'{group}_Best_P_Value'] = 1.0
                sup_data[f'{group}_Conservative_P_Value'] = 1.0
                sup_data[f'{group}_Shapiro_A'] = np.nan
                sup_data[f'{group}_Shapiro_B'] = np.nan
                sup_data[f'{group}_Normal_A'] = False
                sup_data[f'{group}_Normal_B'] = False
                sup_data[f'{group}_Levene_P'] = np.nan
                sup_data[f'{group}_Equal_Vars'] = False
        
        # Determine color based on INDIVIDUAL tests
        global_significant = species_data['Global-t-test'] < 0.05
        
        if global_significant:
            species_data['Cor'] = 'vermelho'
        elif significant_tests_individual > 1:
            species_data['Cor'] = 'laranja'
        elif significant_tests_individual == 1:
            species_data['Cor'] = 'amarelo'
        else:
            species_data['Cor'] = 'cinza'
        
        # Determine color based on CONSERVATIVE tests
        global_significant_cons = species_data['Global-t-test-conservative'] < 0.05
        
        if global_significant_cons:
            species_data['Color'] = 'vermelho'
        elif significant_tests_conservative > 1:
            species_data['Color'] = 'laranja'
        elif significant_tests_conservative == 1:
            species_data['Color'] = 'amarelo'
        else:
            species_data['Color'] = 'cinza'
        
        results.append(species_data)
        sup_stats.append(sup_data)
    
    print(f"✓ Statistics calculated for all {len(results)} species")
    return pd.DataFrame(results), pd.DataFrame(sup_stats)

def apply_filters(stats_df, parameter_results, top30_species, use_conservative=False):
    """
    Applies filters to generate filtered.csv file.
    """
    filter_type = "conservative" if use_conservative else "individual"
    print(f"Applying filters with {filter_type} tests...")
    
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
    
    # Filter 2: At least one significant test OR in TOP30
    if use_conservative:
        # For conservative test, look for columns that ended with 't-test-conservative'
        # but are now renamed to 't-test'
        t_test_cols_for_filter = [col for col in stats_df.columns if col.endswith('t-test-conservative')]
    else:
        # For individual test, use columns that end with 't-test' but not 't-test-conservative'
        t_test_cols_for_filter = [col for col in stats_df.columns if 
                                 col.endswith('t-test') and not col.endswith('t-test-conservative')]
    
    significant_t_tests = (stats_df[t_test_cols_for_filter] < 0.05).any(axis=1)
    in_top30 = stats_df['Species'].isin(top30_species)
    
    filter2 = significant_t_tests | in_top30
    print(f"Species with significant t-test: {significant_t_tests.sum()}")
    print(f"Species in TOP30: {in_top30.sum()}")
    print(f"Species meeting statistical/TOP30 criteria: {filter2.sum()}")
    
    # Apply both filters
    final_filter = filter1 & filter2
    filtered_df = stats_df[final_filter].copy()
    
    print(f"✓ Species after applying all filters: {len(filtered_df)}")
    
    # Select appropriate columns based on test type
    if use_conservative:
        # For conservative file, use conservative columns
        # Remove individual test columns (non-conservative)
        individual_test_cols = [col for col in filtered_df.columns if 
                               col.endswith('t-test') and not col.endswith('t-test-conservative')]
        if individual_test_cols:
            filtered_df = filtered_df.drop(columns=individual_test_cols)
        
        # Remove Cor column (individual) and keep Color (conservative)
        if 'Cor' in filtered_df.columns:
            filtered_df = filtered_df.drop(columns=['Cor'])
        
        # Rename conservative test columns to standard format
        conservative_test_cols = [col for col in filtered_df.columns if col.endswith('t-test-conservative')]
        rename_dict = {}
        for col in conservative_test_cols:
            new_name = col.replace('t-test-conservative', 't-test')
            rename_dict[col] = new_name
        if rename_dict:
            filtered_df = filtered_df.rename(columns=rename_dict)
            
        # Define color column for ordered columns
        cor_col = 'Color'
            
    else:
        # For individual file, use individual columns
        # Remove all conservative columns
        conservative_cols = [col for col in filtered_df.columns if 'conservative' in col.lower() or col == 'Color']
        if conservative_cols:
            filtered_df = filtered_df.drop(columns=conservative_cols)
            
        # Define color column for ordered columns  
        cor_col = 'Cor'
    
    # Order columns in desired sequence
    taxonomy_cols = [cor_col, 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Parameter_Weighted_Score']
    
    # Identify t-test columns (now standardized)
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

def generate_summary_report(filtered_df_individual, filtered_df_conservative, sup_stats_df, original_count):
    """
    Generates summary report of filtering including comparison between methods.
    """
    report_lines = []
    
    report_lines.append("=" * 80)
    report_lines.append("HYBRID FILTERING AND STATISTICAL ANALYSIS REPORT")
    report_lines.append("=" * 80)
    
    report_lines.append(f"Species in original dataset: {original_count}")
    report_lines.append(f"Species after filtering (individual tests): {len(filtered_df_individual)}")
    report_lines.append(f"Species after filtering (conservative test): {len(filtered_df_conservative)}")
    report_lines.append(f"Percentage retained (individual): {len(filtered_df_individual)/original_count*100:.1f}%")
    report_lines.append(f"Percentage retained (conservative): {len(filtered_df_conservative)/original_count*100:.1f}%")
    
    # Conservative test used
    conservative_test = sup_stats_df['Conservative_Test_Type'].iloc[0]
    report_lines.append(f"\nConservative test applied: {conservative_test}")
    
    # Statistics by color - Individual
    color_counts_ind = filtered_df_individual['Cor'].value_counts()
    report_lines.append(f"\nDistribution by significance (INDIVIDUAL TESTS):")
    for color in ['vermelho', 'laranja', 'amarelo', 'cinza']:
        count = color_counts_ind.get(color, 0)
        report_lines.append(f"  {color.capitalize()}: {count} species")
    
    # Statistics by color - Conservative
    color_counts_cons = filtered_df_conservative['Color'].value_counts()
    report_lines.append(f"\nDistribution by significance (CONSERVATIVE TEST):")
    for color in ['vermelho', 'laranja', 'amarelo', 'cinza']:
        count = color_counts_cons.get(color, 0)
        report_lines.append(f"  {color.capitalize()}: {count} species")
    
    # Summary of test types used
    test_types = sup_stats_df['Global_Best_Test_Type'].value_counts()
    report_lines.append(f"\nTypes of individual tests used (Global):")
    for test_type, count in test_types.items():
        report_lines.append(f"  {test_type}: {count} species")
    
    # Normality statistics
    normal_both = sup_stats_df['Global_Normal_A'] & sup_stats_df['Global_Normal_B']
    report_lines.append(f"\nNormality test (Global):")
    report_lines.append(f"  Both conditions normal: {normal_both.sum()} species")
    report_lines.append(f"  At least one non-normal: {(~normal_both).sum()} species")
    
    # Comparison between methods
    individual_species = set(filtered_df_individual['Species'])
    conservative_species = set(filtered_df_conservative['Species'])
    
    only_individual = individual_species - conservative_species
    only_conservative = conservative_species - individual_species
    both_methods = individual_species & conservative_species
    
    report_lines.append(f"\nComparison between methods:")
    report_lines.append(f"  Species detected by both: {len(both_methods)}")
    report_lines.append(f"  Only by individual tests: {len(only_individual)}")
    report_lines.append(f"  Only by conservative test: {len(only_conservative)}")
    
    # Save report
    with open('filtering_report.txt', 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    # Display on screen
    for line in report_lines:
        print(line)

def main():
    print("Starting hybrid statistical filtering analysis...")
    print("=" * 80)
    
    # 1. Load all necessary data
    abundance_data, metadata, parameter_results, top30_species = load_data()
    
    if any(data is None for data in [abundance_data, metadata, parameter_results]):
        print("❌ Error: Could not load all necessary files!")
        return
    
    original_count = len(abundance_data)
    
    # 2. Calculate comprehensive statistics with individual and conservative tests
    print(f"\nCalculating hybrid statistics (individual + conservative)...")
    stats_df, sup_stats_df = calculate_comprehensive_stats(abundance_data, metadata)
    
    # 3. Save supplementary statistics file
    sup_stats_file = 'sup_statistics.csv'
    sup_stats_df.to_csv(sup_stats_file, index=False)
    print(f"✓ Supplementary statistics saved: {sup_stats_file}")
    
    # 4. Apply filters for both methods
    print(f"\nApplying filters with individual tests...")
    filtered_df_individual = apply_filters(stats_df, parameter_results, top30_species, use_conservative=False)
    
    print(f"\nApplying filters with conservative test...")
    filtered_df_conservative = apply_filters(stats_df, parameter_results, top30_species, use_conservative=True)
    
    # 5. Save filtered files
    output_file_individual = 'filtrado.csv'
    output_file_conservative = 'filtrado_conservativetest.csv'
    
    filtered_df_individual.to_csv(output_file_individual, index=False)
    filtered_df_conservative.to_csv(output_file_conservative, index=False)
    
    print(f"✓ Filtered file (individual) saved: {output_file_individual}")
    print(f"✓ Filtered file (conservative) saved: {output_file_conservative}")
    
    # 6. Generate comparative summary report
    generate_summary_report(filtered_df_individual, filtered_df_conservative, sup_stats_df, original_count)
    print(f"✓ Comparative report saved: filtering_report.txt")
    
    print(f"\n{'='*80}")
    print("HYBRID ANALYSIS COMPLETED!")
    print(f"Generated files:")
    print(f"  - {output_file_individual} ({len(filtered_df_individual)} species - individual tests)")
    print(f"  - {output_file_conservative} ({len(filtered_df_conservative)} species - conservative test)")
    print(f"  - {sup_stats_file} (detailed statistics for all species)")
    print(f"  - filtering_report.txt (detailed comparative report)")
    print("="*80)

if __name__ == '__main__':
    main()
