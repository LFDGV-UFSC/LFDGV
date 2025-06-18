import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# 1) Load filtered CSV
df = pd.read_csv('filtrado1.csv', sep=None, engine='python')  # try to detect separator automatically

# 2) Detect color column name
if 'Color' in df.columns:
    color_col = 'Color'
elif 'Cor' in df.columns:
    color_col = 'Cor'
else:
    raise KeyError("Could not find 'Color' or 'Cor' column in filtrado1.csv")

# 3) Map color values (Portuguese) to matplotlib names
color_map = {
    'amarelo': 'yellow',
    'vermelho': 'red',
    'cinza': 'gray',
    'laranja': 'orange',
    # if already in English, pass through
    'yellow':   'yellow',
    'red':      'red',
    'gray':     'gray',
    'orange':   'orange'
}

df['mpl_color'] = df[color_col].map(color_map)
if df['mpl_color'].isnull().any():
    missing = set(df[color_col][df['mpl_color'].isnull()])
    raise ValueError(f"Found unmapped colors: {missing}")

# 4) Sort by score (descending for highest positive impact on top)
df = df.sort_values('Parameter_Weighted_Score', ascending=False)

# 5) Prepare data for plotting
species = df['Species']
scores  = df['Parameter_Weighted_Score'].astype(float)
colors  = df['mpl_color']

# 6) Create the plot
plt.figure(figsize=(8, max(6, len(df)*0.3)))
plt.barh(species, scores, color=colors)
plt.xlabel('Parameter-Weighted Score')
plt.ylabel('Species')
plt.title('Top Contributing Bacterial Species (Welch\'s t-test)')

# 7) Custom legend
legend_elems = [
    Line2D([0], [0], color='red',    lw=6, label='Global'),
    Line2D([0], [0], color='orange', lw=6, label='Multiple Groups'),
    Line2D([0], [0], color='yellow', lw=6, label='One Group'),
    Line2D([0], [0], color='gray',   lw=6, label='Not Significant'),
]
plt.legend(handles=legend_elems, title='Significance Pattern', loc='upper right')
plt.tight_layout()
plt.savefig('parameter_weighted_barplot_welch.png', dpi=300)
plt.close()
print(f"✅ Plot saved as parameter_weighted_barplot_welch.png with {len(df)} species.")
