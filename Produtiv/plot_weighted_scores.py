import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# 1) Carrega o CSV filtrado
df = pd.read_csv('filtrado.csv', sep=None, engine='python')  # tenta detectar sep automaticamente

# 2) Detecta o nome da coluna de cor
if 'Color' in df.columns:
    color_col = 'Color'
elif 'Cor' in df.columns:
    color_col = 'Cor'
else:
    raise KeyError("Não encontrei coluna 'Color' nem 'Cor' em filtrado.csv")

# 3) Mapeia valores de cor (em português) para nomes matplotlib
color_map = {
    'amarelo': 'yellow',
    'vermelho': 'red',
    'cinza': 'gray',
    'laranja': 'orange',
    # caso haja inglês, já passamos direto
    'yellow':   'yellow',
    'red':      'red',
    'gray':     'gray',
    'orange':   'orange'
}
df['mpl_color'] = df[color_col].map(color_map)
if df['mpl_color'].isnull().any():
    missing = set(df[color_col][df['mpl_color'].isnull()])
    raise ValueError(f"Encontrei cores não mapeadas: {missing}")

# 4) Ordena pelo score (descendente para maior impacto positivo no topo)
df = df.sort_values('Parameter_Weighted_Score', ascending=False)

# 5) Prepara dados para o plot
species = df['Species']
scores  = df['Parameter_Weighted_Score'].astype(float)
colors  = df['mpl_color']

# 6) Cria o gráfico
plt.figure(figsize=(8, max(6, len(df)*0.3)))
plt.barh(species, scores, color=colors)
plt.xlabel('Productivity-Weighted Score')
plt.ylabel('Species')
plt.title ('Top Contributing Bacterial Species to Garlic Yield')

# 7) Legenda customizada
legend_elems = [
    Line2D([0], [0], color='red',    lw=6, label='Global'),
    Line2D([0], [0], color='orange', lw=6, label='Multiple Farms'),
    Line2D([0], [0], color='yellow', lw=6, label='One Farm'),
    Line2D([0], [0], color='gray',   lw=6, label='Not Significant'),
]
plt.legend(handles=legend_elems, title='Significance Pattern', loc='upper right')

plt.tight_layout()
plt.savefig('productivity_weighted_barplot.png', dpi=300)
plt.close()

print(f"✅ Plot saved as productivity_weighted_barplot.png with {len(df)} species.")

