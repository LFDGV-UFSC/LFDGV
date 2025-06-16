import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def load_abundance_data(file):
    """
    Carrega o arquivo de abundância com dados taxonômicos e amostrais.
    """
    return pd.read_csv(file)

def load_metadata(file):
    """
    Carrega o arquivo de metadata com informações das amostras.
    Formato esperado: sample-id, group, condition, parameter
    """
    try:
        # Tenta carregar como TSV primeiro, depois como CSV
        if file.endswith('.tsv'):
            metadata = pd.read_csv(file, sep='\t')
        else:
            metadata = pd.read_csv(file)
        
        # Verifica se tem as colunas necessárias
        required_cols = ['sample-id', 'group', 'condition', 'parameter']
        missing_cols = [col for col in required_cols if col not in metadata.columns]
        
        if missing_cols:
            raise ValueError(f"Colunas obrigatórias ausentes no metadata: {missing_cols}")
        
        return metadata
    
    except FileNotFoundError:
        print(f"❌ Erro: Arquivo de metadata {file} não encontrado!")
        print("Criando arquivo de exemplo 'metadata_example.tsv'...")
        create_example_metadata()
        return None

def create_example_metadata():
    """
    Cria um arquivo de exemplo do metadata.
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
    print("✓ Arquivo 'metadata_example.tsv' criado como modelo.")

def calculate_parameter_factors(metadata):
    """
    Calcula os fatores de ponderação baseados no parâmetro fornecido no metadata.
    """
    # Obtém o valor do parâmetro para cada grupo
    group_parameters = {}
    for group in metadata['group'].unique():
        group_data = metadata[metadata['group'] == group]
        # Assume que o parâmetro é o mesmo para todas as amostras do grupo
        parameter_value = group_data['parameter'].iloc[0]
        group_parameters[group] = parameter_value
    
    mean_parameter = np.mean(list(group_parameters.values()))
    
    parameter_factors = {}
    print(f"Calculando fatores de ponderação baseados no parâmetro:")
    for group, param_value in group_parameters.items():
        factor = param_value / mean_parameter
        parameter_factors[group] = factor
        print(f"{group}: {param_value} | Fator: {factor:.4f} | Peso: {((factor-1)*100):+.1f}%")
    
    return parameter_factors, mean_parameter, group_parameters

def get_sample_mapping(metadata):
    """
    Cria mapeamento entre sample-id e grupos/condições.
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
    Calcula a variação média (Depois - Antes) para cada espécie em cada grupo.
    """
    taxonomy_cols = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    sample_cols = [col for col in df.columns if col not in taxonomy_cols]
    
    # Cria mapeamento das amostras
    sample_mapping = get_sample_mapping(metadata)
    
    # Obtém grupos únicos
    groups = metadata['group'].unique()
    
    results = []
    
    print(f"Processando {len(df)} espécies em {len(groups)} grupos...")
    
    for idx, (_, row) in enumerate(df.iterrows()):
        if idx > 0 and idx % 100 == 0:
            print(f"  Processadas {idx} espécies...")
            
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
            # Encontra amostras deste grupo nas condições A e B
            antes_cols = []
            depois_cols = []
            
            for sample_id in sample_cols:
                if sample_id in sample_mapping:
                    sample_info = sample_mapping[sample_id]
                    if sample_info['group'] == group:
                        if sample_info['condition'] == 'A':
                            antes_cols.append(sample_id)
                        elif sample_info['condition'] == 'B':
                            depois_cols.append(sample_id)
            
            if antes_cols and depois_cols:
                # Média das amostras antes e depois
                antes_values = pd.to_numeric(row[antes_cols], errors='coerce')
                depois_values = pd.to_numeric(row[depois_cols], errors='coerce')
                
                antes_mean = antes_values.mean()
                depois_mean = depois_values.mean()
                
                # Substitui NaN por 0
                antes_mean = 0 if pd.isna(antes_mean) else antes_mean
                depois_mean = 0 if pd.isna(depois_mean) else depois_mean
                
                # Variação absoluta e relativa
                variation_abs = depois_mean - antes_mean
                
                # Variação relativa (evita divisão por zero)
                if antes_mean > 0:
                    variation_rel = (depois_mean - antes_mean) / antes_mean
                else:
                    # Se antes era 0 e depois > 0, considera como aparecimento (100% de aumento)
                    variation_rel = 1.0 if depois_mean > 0 else 0.0
                
                species_data[f'{group}_antes_mean'] = antes_mean
                species_data[f'{group}_depois_mean'] = depois_mean
                species_data[f'{group}_variation_abs'] = variation_abs
                species_data[f'{group}_variation_rel'] = variation_rel
                species_data[f'{group}_present'] = 1 if (antes_mean > 0 or depois_mean > 0) else 0
            else:
                # Espécie não está presente neste grupo
                species_data[f'{group}_antes_mean'] = 0
                species_data[f'{group}_depois_mean'] = 0
                species_data[f'{group}_variation_abs'] = 0
                species_data[f'{group}_variation_rel'] = 0
                species_data[f'{group}_present'] = 0
        
        results.append(species_data)
    
    print(f"✓ Todas as {len(results)} espécies processadas")
    return pd.DataFrame(results)

def calculate_parameter_weighted_score(variations_df, parameter_factors, groups):
    """
    Calcula o score ponderado pelo parâmetro para cada espécie.
    """
    results = []
    
    print(f"Calculando scores para {len(variations_df)} espécies...")
    
    for idx, (_, row) in enumerate(variations_df.iterrows()):
        if idx > 0 and idx % 100 == 0:
            print(f"  Scores calculados para {idx} espécies...")
            
        species_scores = []
        species_weights = []
        groups_with_data = []
        
        for group in groups:
            # Verifica se a espécie tem dados neste grupo
            if row[f'{group}_present'] == 1:
                variation = row[f'{group}_variation_rel']
                parameter_factor = parameter_factors[group]
                
                # Score ponderado: variação × fator do parâmetro
                weighted_score = variation * parameter_factor
                species_scores.append(weighted_score)
                species_weights.append(parameter_factor)
                groups_with_data.append(group)
        
        # Processa a espécie
        if species_scores:
            # Score final: média ponderada dos scores dos grupos onde a espécie está presente
            final_score = np.average(species_scores, weights=species_weights)
            
            # Estatísticas adicionais
            n_groups_present = len(species_scores)
            max_score = max(species_scores)
            min_score = min(species_scores)
            std_score = np.std(species_scores) if len(species_scores) > 1 else 0
            
        else:
            # Espécie sem dados em nenhum grupo
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
        
        # Adiciona os scores individuais por grupo
        for group in groups:
            if row[f'{group}_present'] == 1:
                variation = row[f'{group}_variation_rel']
                weighted_score = variation * parameter_factors[group]
                result[f'{group}_Score'] = weighted_score
                result[f'{group}_Variation_Rel'] = variation
                result[f'{group}_Antes_Mean'] = row[f'{group}_antes_mean']
                result[f'{group}_Depois_Mean'] = row[f'{group}_depois_mean']
            else:
                result[f'{group}_Score'] = 0
                result[f'{group}_Variation_Rel'] = 0
                result[f'{group}_Antes_Mean'] = 0
                result[f'{group}_Depois_Mean'] = 0
        
        results.append(result)
    
    print(f"✓ Scores calculados para todas as {len(results)} espécies")
    return pd.DataFrame(results)

def interpret_results(results_df):
    """
    Interpreta os resultados e adiciona classificações.
    """
    results_df = results_df.copy()
    
    # Classifica o impacto no parâmetro
    def classify_impact(score):
        if score > 0.1:
            return "Alto Impacto Positivo"
        elif score > 0.05:
            return "Moderado Impacto Positivo"
        elif score > 0:
            return "Baixo Impacto Positivo"
        elif score > -0.05:
            return "Baixo Impacto Negativo"
        elif score > -0.1:
            return "Moderado Impacto Negativo"
        else:
            return "Alto Impacto Negativo"
    
    results_df['Impact_Classification'] = results_df['Parameter_Weighted_Score'].apply(classify_impact)
    
    # Ordena por score (mais impactantes primeiro)
    results_df = results_df.sort_values('Parameter_Weighted_Score', ascending=False)
    
    return results_df

def generate_summary_stats(results_df, group_parameters, output_summary=True):
    """
    Gera estatísticas resumo dos resultados e salva em arquivo.
    """
    summary_lines = []
    
    summary_lines.append("=" * 80)
    summary_lines.append("RESUMO ESTATÍSTICO DA ANÁLISE")
    summary_lines.append("=" * 80)
    
    # Informações dos grupos e parâmetros
    summary_lines.append("\nPARÂMETROS DOS GRUPOS:")
    for group, param_value in group_parameters.items():
        summary_lines.append(f"{group}: {param_value}")
    
    total_species = len(results_df)
    positive_impact = len(results_df[results_df['Parameter_Weighted_Score'] > 0])
    negative_impact = len(results_df[results_df['Parameter_Weighted_Score'] < 0])
    
    summary_lines.append(f"\nTotal de espécies analisadas: {total_species}")
    summary_lines.append(f"Espécies com impacto positivo: {positive_impact} ({positive_impact/total_species*100:.1f}%)")
    summary_lines.append(f"Espécies com impacto negativo: {negative_impact} ({negative_impact/total_species*100:.1f}%)")
    
    summary_lines.append(f"\nTOP 30 ESPÉCIES COM MAIOR IMPACTO POSITIVO:")
    top_positive = results_df[results_df['Parameter_Weighted_Score'] > 0].head(30)
    if len(top_positive) > 0:
        for i, (_, row) in enumerate(top_positive.iterrows(), 1):
            summary_lines.append(f"{i:2d}. {row['Species'][:50]:<50} | Score: {row['Parameter_Weighted_Score']:+.4f}")
    else:
        summary_lines.append("Nenhuma espécie com impacto positivo encontrada.")
    
    summary_lines.append(f"\nTOP 30 ESPÉCIES COM MAIOR IMPACTO NEGATIVO:")
    negative_species = results_df[results_df['Parameter_Weighted_Score'] < 0]
    if len(negative_species) > 0:
        top_negative = negative_species.tail(30).iloc[::-1]  # Inverte para mostrar os mais negativos primeiro
        for i, (_, row) in enumerate(top_negative.iterrows(), 1):
            summary_lines.append(f"{i:2d}. {row['Species'][:50]:<50} | Score: {row['Parameter_Weighted_Score']:+.4f}")
    else:
        summary_lines.append("Nenhuma espécie com impacto negativo encontrada.")
    
    summary_lines.append(f"\nESPÉCIES SEM DADOS (Score = 0):")
    zero_species = results_df[results_df['Parameter_Weighted_Score'] == 0]
    summary_lines.append(f"Total: {len(zero_species)} espécies")
    
    # Salva em arquivo
    if output_summary:
        with open('top30.txt', 'w', encoding='utf-8') as f:
            f.write('\n'.join(summary_lines))
        print("✓ Resumo das TOP 30 espécies salvo em: top30.txt")
    
    # Mostra na tela
    for line in summary_lines:
        print(line)

def main():
    print("Iniciando análise com parâmetros externos...")
    print("=" * 80)
    
    # 1. Carrega os dados de abundância
    input_file = 'abundance_data.csv'
    try:
        df = load_abundance_data(input_file)
        print(f"✓ Dados de abundância carregados: {len(df)} espécies, {len(df.columns)-6} amostras")
    except FileNotFoundError:
        print(f"❌ Erro: Arquivo {input_file} não encontrado!")
        return
    
    # 2. Carrega o metadata
    metadata_files = ['metadata.tsv', 'metadata.txt', 'metadata.csv']
    metadata = None
    
    for metadata_file in metadata_files:
        try:
            metadata = load_metadata(metadata_file)
            if metadata is not None:
                print(f"✓ Metadata carregado de: {metadata_file}")
                break
        except:
            continue
    
    if metadata is None:
        print("❌ Nenhum arquivo de metadata válido encontrado!")
        print("Use o arquivo 'metadata_example.tsv' como modelo.")
        return
    
    # 3. Calcula fatores do parâmetro
    print(f"\nCalculando fatores de ponderação:")
    parameter_factors, mean_param, group_parameters = calculate_parameter_factors(metadata)
    print(f"Valor médio do parâmetro: {mean_param:.2f}")
    
    # 4. Verifica correspondência entre amostras
    sample_cols = [col for col in df.columns if col not in ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]
    metadata_samples = set(metadata['sample-id'].tolist())
    df_samples = set(sample_cols)
    
    missing_in_metadata = df_samples - metadata_samples
    missing_in_df = metadata_samples - df_samples
    
    if missing_in_metadata:
        print(f"\n⚠️  Amostras em abundance_data.csv mas não no metadata: {missing_in_metadata}")
    if missing_in_df:
        print(f"⚠️  Amostras no metadata mas não em abundance_data.csv: {missing_in_df}")
    
    # 5. Calcula variações das espécies  
    print(f"\nCalculando variações de abundância por espécie e grupo...")
    groups = metadata['group'].unique()
    variations_df = calculate_species_variation(df, metadata)
    print(f"✓ Variações calculadas para {len(variations_df)} espécies em {len(groups)} grupos")
    
    # 6. Calcula scores ponderados
    print(f"\nCalculando scores ponderados pelo parâmetro...")
    results_df = calculate_parameter_weighted_score(variations_df, parameter_factors, groups)
    print(f"✓ Scores calculados para {len(results_df)} espécies")
    
    # 7. Interpreta resultados
    results_df = interpret_results(results_df)
    
    # 8. Salva resultados COMPLETOS
    output_file = 'parameter_weighted_analysis.csv'
    results_df.to_csv(output_file, index=False)
    print(f"✓ Resultados COMPLETOS salvos em: {output_file}")
    print(f"✓ Todas as {len(results_df)} espécies incluídas no arquivo CSV")
    
    # 9. Gera resumo estatístico e arquivo TOP 30
    generate_summary_stats(results_df, group_parameters)
    
    print(f"\n{'='*80}")
    print("ANÁLISE CONCLUÍDA!")
    print(f"Arquivos de saída:")
    print(f"  - {output_file} (resultados completos)")
    print(f"  - top30.txt (resumo das top 30 espécies)")
    print("="*80)

if __name__ == '__main__':
    main()
