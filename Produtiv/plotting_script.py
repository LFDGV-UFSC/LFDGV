import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

def get_parameter_name():
    """
    Requests parameter name from user with validation.
    """
    while True:
        parameter = input("My parameter for correlation is: ").strip()
        
        # Check if it's only one word (no spaces)
        if len(parameter.split()) == 1 and parameter != "":
            return parameter
        else:
            print("Just one word... Eg.: Productivity, or Phosphorus...")

def detect_color_column(df):
    """
    Automatically detects the color column in the DataFrame.
    """
    if 'Color' in df.columns:
        return 'Color'
    elif 'Cor' in df.columns:
        return 'Cor'
    else:
        raise KeyError("Could not find 'Color' or 'Cor' column in file")

def prepare_dataframe(filepath):
    """
    Loads and prepares DataFrame for plotting.
    """
    print(f"Loading {filepath}...")
    
    # Check if file exists
    if not os.path.exists(filepath):
        print(f"⚠️ File {filepath} not found!")
        return None
    
    # Load CSV
    df = pd.read_csv(filepath, sep=None, engine='python')
    
    # Detect color column
    color_col = detect_color_column(df)
    
    # Map color values (Portuguese/English) to matplotlib
    color_map = {
        'amarelo': 'yellow',
        'vermelho': 'red',
        'cinza': 'gray',
        'laranja': 'orange',
        # if already in English, pass through
        'yellow': 'yellow',
        'red': 'red',
        'gray': 'gray',
        'orange': 'orange'
    }
    
    df['mpl_color'] = df[color_col].map(color_map)
    
    # Check for unmapped colors
    if df['mpl_color'].isnull().any():
        missing = set(df[color_col][df['mpl_color'].isnull()])
        raise ValueError(f"Found unmapped colors in {filepath}: {missing}")
    
    # Sort by score (descending for highest positive impact on top)
    df = df.sort_values('Parameter_Weighted_Score', ascending=False)
    
    print(f"✓ {filepath} loaded: {len(df)} species")
    return df

def create_plot(df, parameter_name, output_filename, plot_title_suffix):
    """
    Creates a horizontal bar chart for the DataFrame.
    """
    # Prepare data for plotting
    species = df['Species']
    scores = df['Parameter_Weighted_Score'].astype(float)
    colors = df['mpl_color']
    
    # Create the plot
    plt.figure(figsize=(10, max(6, len(df)*0.3)))
    plt.barh(species, scores, color=colors)
    
    # Labels and title
    plt.xlabel(f'{parameter_name}-Weighted Score')
    plt.ylabel('Species')
    plt.title(f'Top Contributing Bacterial Species to {parameter_name} {plot_title_suffix}')
    
    # Custom legend
    legend_elems = [
        Line2D([0], [0], color='red', lw=6, label='Global'),
        Line2D([0], [0], color='orange', lw=6, label='Multiple Groups'),
        Line2D([0], [0], color='yellow', lw=6, label='One Group'),
        Line2D([0], [0], color='gray', lw=6, label='Not Significant'),
    ]
    plt.legend(handles=legend_elems, title='Significance Pattern', loc='upper right')
    
    # Layout adjustments
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Plot saved as {output_filename} with {len(df)} species.")

def main():
    print("=" * 60)
    print("BACTERIAL IMPACT PLOTS GENERATOR")
    print("=" * 60)
    
    # 1. Request parameter name
    parameter_name = get_parameter_name()
    print(f"✓ Parameter selected: {parameter_name}")
    
    # 2. Define file names
    file1 = 'filtrado.csv'
    file2 = 'filtrado_conservativetest.csv'
    
    output1 = f'{parameter_name}_weighted_barplot.png'
    output2 = f'{parameter_name}_weighted_barplot2.png'
    
    # 3. Process file 1 (individual tests)
    print(f"\n--- Processing {file1} (Individual Tests) ---")
    df1 = prepare_dataframe(file1)
    
    if df1 is not None:
        create_plot(df1, parameter_name, output1, "(Individual Tests)")
    else:
        print(f"❌ Error processing {file1}")
    
    # 4. Process file 2 (conservative test)
    print(f"\n--- Processing {file2} (Conservative Test) ---")
    df2 = prepare_dataframe(file2)
    
    if df2 is not None:
        create_plot(df2, parameter_name, output2, "(Conservative Test)")
    else:
        print(f"❌ Error processing {file2}")
    
    # 5. Final summary
    print(f"\n{'='*60}")
    print("SUMMARY OF GENERATED PLOTS:")
    print(f"{'='*60}")
    
    if df1 is not None:
        print(f"✅ {output1}")
        print(f"   - {len(df1)} species (optimized individual tests)")
        
        # Distribution by color
        color_dist1 = df1['mpl_color'].value_counts()
        print(f"   - Distribution: ", end="")
        colors_summary = []
        if 'red' in color_dist1.index:
            colors_summary.append(f"{color_dist1['red']} global")
        if 'orange' in color_dist1.index:
            colors_summary.append(f"{color_dist1['orange']} multiple")
        if 'yellow' in color_dist1.index:
            colors_summary.append(f"{color_dist1['yellow']} single")
        if 'gray' in color_dist1.index:
            colors_summary.append(f"{color_dist1['gray']} n.s.")
        print(", ".join(colors_summary))
    
    if df2 is not None:
        print(f"✅ {output2}")
        print(f"   - {len(df2)} species (conservative Mann-Whitney U test)")
        
        # Distribution by color
        color_dist2 = df2['mpl_color'].value_counts()
        print(f"   - Distribution: ", end="")
        colors_summary = []
        if 'red' in color_dist2.index:
            colors_summary.append(f"{color_dist2['red']} global")
        if 'orange' in color_dist2.index:
            colors_summary.append(f"{color_dist2['orange']} multiple")
        if 'yellow' in color_dist2.index:
            colors_summary.append(f"{color_dist2['yellow']} single")
        if 'gray' in color_dist2.index:
            colors_summary.append(f"{color_dist2['gray']} n.s.")
        print(", ".join(colors_summary))
    
    # 6. Comparison between methods (if both worked)
    if df1 is not None and df2 is not None:
        common_species = set(df1['Species']) & set(df2['Species'])
        only_individual = set(df1['Species']) - set(df2['Species'])
        only_conservative = set(df2['Species']) - set(df1['Species'])
        
        print(f"\n📊 COMPARISON BETWEEN METHODS:")
        print(f"   - Species detected by both: {len(common_species)}")
        print(f"   - Only by individual tests: {len(only_individual)}")
        print(f"   - Only by conservative test: {len(only_conservative)}")
        
        if len(common_species) > 0:
            print(f"   - Agreement rate: {len(common_species)/len(set(df1['Species']) | set(df2['Species']))*100:.1f}%")
    
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
