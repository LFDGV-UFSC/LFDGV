import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

def get_parameter_name():
    """
    Solicita o nome do parâmetro ao usuário com validação.
    """
    while True:
        parameter = input("My parameter for correlation is: ").strip()
        
        # Verifica se é apenas uma palavra (sem espaços)
        if len(parameter.split()) == 1 and parameter != "":
            return parameter
        else:
            print("Just one word... Eg.: Productivity, or Phosphorus...")

def detect_color_column(df):
    """
    Detecta automaticamente a coluna de cor no DataFrame.
    """
    if 'Color' in df.columns:
        return 'Color'
    elif 'Cor' in df.columns:
        return 'Cor'
    else:
        raise KeyError("Não encontrei coluna 'Color' nem 'Cor' no arquivo")

def prepare_dataframe(filepath):
    """
    Carrega e prepara o DataFrame para plotagem.
    """
    print(f"Carregando {filepath}...")
    
    # Verifica se o arquivo existe
    if not os.path.exists(filepath):
        print(f"⚠️ Arquivo {filepath} não encontrado!")
        return None
    
    # Carrega o CSV
    df = pd.read_csv(filepath, sep=None, engine='python')
    
    # Detecta coluna de cor
    color_col = detect_color_column(df)
    
    # Mapeia valores de cor (português/inglês) para matplotlib
    color_map = {
        'amarelo': 'yellow',
        'vermelho': 'red',
        'cinza': 'gray',
        'laranja': 'orange',
        # caso haja inglês, já passa direto
        'yellow': 'yellow',
        'red': 'red',
        'gray': 'gray',
        'orange': 'orange'
    }
    
    df['mpl_color'] = df[color_col].map(color_map)
    
    # Verifica se há cores não mapeadas
    if df['mpl_color'].isnull().any():
        missing = set(df[color_col][df['mpl_color'].isnull()])
        raise ValueError(f"Encontrei cores não mapeadas em {filepath}: {missing}")
    
    # Ordena pelo score (descendente para maior impacto positivo no topo)
    df = df.sort_values('Parameter_Weighted_Score', ascending=False)
    
    print(f"✓ {filepath} carregado: {len(df)} espécies")
    return df

def create_plot(df, parameter_name, output_filename, plot_title_suffix):
    """
    Cria um gráfico de barras horizontais para o DataFrame.
    """
    # Prepara dados para o plot
    species = df['Species']
    scores = df['Parameter_Weighted_Score'].astype(float)
    colors = df['mpl_color']
    
    # Cria o gráfico
    plt.figure(figsize=(10, max(6, len(df)*0.3)))
    plt.barh(species, scores, color=colors)
    
    # Labels e título
    plt.xlabel(f'{parameter_name}-Weighted Score')
    plt.ylabel('Species')
    plt.title(f'Top Contributing Bacterial Species to {parameter_name} {plot_title_suffix}')
    
    # Legenda customizada
    legend_elems = [
        Line2D([0], [0], color='red', lw=6, label='Global'),
        Line2D([0], [0], color='orange', lw=6, label='Multiple Groups'),
        Line2D([0], [0], color='yellow', lw=6, label='One Group'),
        Line2D([0], [0], color='gray', lw=6, label='Not Significant'),
    ]
    plt.legend(handles=legend_elems, title='Significance Pattern', loc='upper right')
    
    # Ajustes de layout
    plt.tight_layout()
    
    # Salva o gráfico
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Plot saved as {output_filename} with {len(df)} species.")

def main():
    print("=" * 60)
    print("GERADOR DE GRÁFICOS DE IMPACTO BACTERIANO")
    print("=" * 60)
    
    # 1. Solicita o nome do parâmetro
    parameter_name = get_parameter_name()
    print(f"✓ Parâmetro selecionado: {parameter_name}")
    
    # 2. Define nomes dos arquivos
    file1 = 'filtrado.csv'
    file2 = 'filtrado_conservativetest.csv'
    
    output1 = f'{parameter_name}_weighted_barplot.png'
    output2 = f'{parameter_name}_weighted_barplot2.png'
    
    # 3. Processa arquivo 1 (testes individuais)
    print(f"\n--- Processando {file1} (Testes Individuais) ---")
    df1 = prepare_dataframe(file1)
    
    if df1 is not None:
        create_plot(df1, parameter_name, output1, "(Individual Tests)")
    else:
        print(f"❌ Erro ao processar {file1}")
    
    # 4. Processa arquivo 2 (teste conservador)
    print(f"\n--- Processando {file2} (Teste Conservador) ---")
    df2 = prepare_dataframe(file2)
    
    if df2 is not None:
        create_plot(df2, parameter_name, output2, "(Conservative Test)")
    else:
        print(f"❌ Erro ao processar {file2}")
    
    # 5. Resumo final
    print(f"\n{'='*60}")
    print("RESUMO DOS GRÁFICOS GERADOS:")
    print(f"{'='*60}")
    
    if df1 is not None:
        print(f"✅ {output1}")
        print(f"   - {len(df1)} espécies (testes individuais otimizados)")
        
        # Distribui por cor
        color_dist1 = df1['mpl_color'].value_counts()
        print(f"   - Distribuição: ", end="")
        colors_summary = []
        if 'red' in color_dist1.index:
            colors_summary.append(f"{color_dist1['red']} global")
        if 'orange' in color_dist1.index:
            colors_summary.append(f"{color_dist1['orange']} múltiplos")
        if 'yellow' in color_dist1.index:
            colors_summary.append(f"{color_dist1['yellow']} único")
        if 'gray' in color_dist1.index:
            colors_summary.append(f"{color_dist1['gray']} n.s.")
        print(", ".join(colors_summary))
    
    if df2 is not None:
        print(f"✅ {output2}")
        print(f"   - {len(df2)} espécies (teste conservador Mann-Whitney U)")
        
        # Distribui por cor
        color_dist2 = df2['mpl_color'].value_counts()
        print(f"   - Distribuição: ", end="")
        colors_summary = []
        if 'red' in color_dist2.index:
            colors_summary.append(f"{color_dist2['red']} global")
        if 'orange' in color_dist2.index:
            colors_summary.append(f"{color_dist2['orange']} múltiplos")
        if 'yellow' in color_dist2.index:
            colors_summary.append(f"{color_dist2['yellow']} único")
        if 'gray' in color_dist2.index:
            colors_summary.append(f"{color_dist2['gray']} n.s.")
        print(", ".join(colors_summary))
    
    # 6. Comparação entre métodos (se ambos funcionaram)
    if df1 is not None and df2 is not None:
        common_species = set(df1['Species']) & set(df2['Species'])
        only_individual = set(df1['Species']) - set(df2['Species'])
        only_conservative = set(df2['Species']) - set(df1['Species'])
        
        print(f"\n📊 COMPARAÇÃO ENTRE MÉTODOS:")
        print(f"   - Espécies detectadas por ambos: {len(common_species)}")
        print(f"   - Apenas por testes individuais: {len(only_individual)}")
        print(f"   - Apenas por teste conservador: {len(only_conservative)}")
        
        if len(common_species) > 0:
            print(f"   - Taxa de concordância: {len(common_species)/len(set(df1['Species']) | set(df2['Species']))*100:.1f}%")
    
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
