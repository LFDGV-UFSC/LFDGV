import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def load_data():
    """
    Carrega todos os arquivos necessários para a análise.
    """
    try:
        # Carrega dados de abundância
        abundance_data = pd.read_csv('abundance_data.csv')
        print(f"✓ abundance_data.csv carregado: {len(abundance_data)} espécies")
        
        # Carrega metadata
        metadata_files = ['metadata.tsv', 'metadata.txt', 'metadata.csv']
        metadata = None
        
        for metadata_file in metadata_files:
            try:
                if metadata_file.endswith('.tsv'):
                    metadata = pd.read_csv(metadata_file, sep='\t')
                else:
                    metadata = pd.read_csv(metadata_file)
                print(f"✓ {metadata_file} carregado: {len(metadata)} amostras")
                break
            except FileNotFoundError:
                continue
        
        if metadata is None:
            raise FileNotFoundError("Nenhum arquivo de metadata encontrado")
        
        # Carrega resultados da análise anterior
        parameter_results = pd.read_csv('parameter_weighted_analysis.csv')
        print(f"✓ parameter_weighted_analysis.csv carregado: {len(parameter_results)} espécies")
        
        # Carrega lista das top 30 espécies
        top30_species = extract_top30_species('top30.txt')
        print(f"✓ top30.txt processado: {len(top30_species)} espécies extraídas")
        
        return abundance_data, metadata, parameter_results, top30_species
        
    except FileNotFoundError as e:
        print(f"❌ Erro ao carregar arquivos: {e}")
        return None, None, None, None

def extract_top30_species(top30_file):
    """
    Extrai os nomes das espécies do arquivo top30.txt.
    """
    try:
        with open(top30_file, 'r', encoding='utf-8') as f:
            content = f.read()
        
        species_names = set()
        lines = content.split('\n')
        
        for line in lines:
            # Procura por linhas que começam com números (rankings)
            if line.strip() and (line.strip()[0].isdigit() or line.strip().startswith(' ')):
                # Tenta extrair o nome da espécie (após o número e antes do |)
                if '|' in line:
                    species_part = line.split('|')[0].strip()
                    # Remove o número do ranking
                    if '. ' in species_part:
                        species_name = species_part.split('. ', 1)[1].strip()
                        species_names.add(species_name)
        
        return list(species_names)
        
    except FileNotFoundError:
        print("⚠️ Arquivo top30.txt não encontrado. Continuando sem filtro por TOP30.")
        return []

def get_sample_mapping(metadata):
    """
    Cria mapeamento entre sample-id e grupos/condições.
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
    Testa normalidade usando Shapiro-Wilk.
    Retorna p-value do teste.
    """
    if len(data) < 3:
        return np.nan  # Shapiro-Wilk precisa de pelo menos 3 observações
    
    # Remove zeros e valores muito pequenos que podem causar problemas
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
    Testa igualdade de variâncias usando teste de Levene.
    Retorna p-value do teste.
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
    Escolhe o melhor teste estatístico baseado em:
    1. Normalidade (Shapiro-Wilk)
    2. Igualdade de variâncias (Levene)
    
    Retorna: test_name, p_value, details
    """
    
    if len(data1) < 2 or len(data2) < 2:
        return "Insufficient_data", 1.0, {"reason": "Less than 2 observations per group"}
    
    # Testa normalidade
    norm_p1 = test_normality(data1)
    norm_p2 = test_normality(data2)
    
    # Considera normal se p > 0.05 (ou se não pode testar)
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
        # Testa igualdade de variâncias
        levene_p = test_equal_variances(data1, data2)
        equal_vars = pd.isna(levene_p) or levene_p > 0.05
        
        details["levene_p"] = levene_p
        details["equal_vars"] = equal_vars
        
        if equal_vars:
            # Teste t de Student (variâncias iguais)
            try:
                t_stat, p_value = stats.ttest_ind(data1, data2, equal_var=True)
                return "Student_t_test", p_value, details
            except:
                return "Test_failed", 1.0, details
        else:
            # Teste t de Welch (variâncias diferentes)
            try:
                t_stat, p_value = stats.ttest_ind(data1, data2, equal_var=False)
                return "Welch_t_test", p_value, details
            except:
                return "Test_failed", 1.0, details
    else:
        # Teste não paramétrico (Mann-Whitney U)
        try:
            stat, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
            return "Mann_Whitney_U", p_value, details
        except:
            return "Test_failed", 1.0, details

def apply_conservative_test(data1, data2, test_type):
    """
    Aplica o teste conservador especificado.
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
    Determina qual teste conservador usar baseado na análise de todas as espécies.
    """
    print("Determinando teste conservador baseado em análise global...")
    
    taxonomy_cols = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    sample_cols = [col for col in abundance_data.columns if col not in taxonomy_cols]
    sample_mapping = get_sample_mapping(metadata)
    
    # Separa amostras por condição
    condition_a_samples = [s for s in sample_cols if s in sample_mapping and sample_mapping[s]['condition'] == 'A']
    condition_b_samples = [s for s in sample_cols if s in sample_mapping and sample_mapping[s]['condition'] == 'B']
    groups = metadata['group'].unique()
    
    # Contadores para decisão
    any_non_normal = False
    any_unequal_vars = False
    
    # Amostra uma parte das espécies para análise (para economizar tempo)
    sample_size = min(100, len(abundance_data))  # Analisa até 100 espécies para decidir
    sample_indices = np.random.choice(len(abundance_data), sample_size, replace=False)
    
    print(f"Analisando {sample_size} espécies para determinar teste conservador...")
    
    for idx in sample_indices:
        row = abundance_data.iloc[idx]
        
        # Análise global
        condition_a_values = pd.to_numeric(row[condition_a_samples], errors='coerce').fillna(0)
        condition_b_values = pd.to_numeric(row[condition_b_samples], errors='coerce').fillna(0)
        
        # Testa normalidade global
        norm_p1 = test_normality(condition_a_values)
        norm_p2 = test_normality(condition_b_values)
        
        is_normal1 = pd.isna(norm_p1) or norm_p1 > 0.05
        is_normal2 = pd.isna(norm_p2) or norm_p2 > 0.05
        
        if not (is_normal1 and is_normal2):
            any_non_normal = True
        
        # Se ambos normais, testa variâncias
        if is_normal1 and is_normal2:
            levene_p = test_equal_variances(condition_a_values, condition_b_values)
            if not (pd.isna(levene_p) or levene_p > 0.05):
                any_unequal_vars = True
        
        # Análise por grupos
        for group in groups:
            group_a_samples = [s for s in sample_cols if s in sample_mapping and 
                             sample_mapping[s]['group'] == group and sample_mapping[s]['condition'] == 'A']
            group_b_samples = [s for s in sample_cols if s in sample_mapping and 
                             sample_mapping[s]['group'] == group and sample_mapping[s]['condition'] == 'B']
            
            if group_a_samples and group_b_samples:
                group_a_values = pd.to_numeric(row[group_a_samples], errors='coerce').fillna(0)
                group_b_values = pd.to_numeric(row[group_b_samples], errors='coerce').fillna(0)
                
                # Testa normalidade do grupo
                norm_p1 = test_normality(group_a_values)
                norm_p2 = test_normality(group_b_values)
                
                is_normal1 = pd.isna(norm_p1) or norm_p1 > 0.05
                is_normal2 = pd.isna(norm_p2) or norm_p2 > 0.05
                
                if not (is_normal1 and is_normal2):
                    any_non_normal = True
                
                # Se ambos normais, testa variâncias
                if is_normal1 and is_normal2:
                    levene_p = test_equal_variances(group_a_values, group_b_values)
                    if not (pd.isna(levene_p) or levene_p > 0.05):
                        any_unequal_vars = True
    
    # Decide o teste conservador
    if any_non_normal:
        conservative_test = "Mann_Whitney_U"
        reason = "Dados não normais detectados em algumas espécies"
    elif any_unequal_vars:
        conservative_test = "Welch_t_test"
        reason = "Variâncias desiguais detectadas em algumas espécies"
    else:
        conservative_test = "Student_t_test"
        reason = "Dados normais com variâncias iguais em todas as espécies analisadas"
    
    print(f"✓ Teste conservador escolhido: {conservative_test}")
    print(f"  Razão: {reason}")
    
    return conservative_test

def calculate_comprehensive_stats(abundance_data, metadata):
    """
    Calcula estatísticas completas com testes individuais e conservadores.
    """
    taxonomy_cols = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    sample_cols = [col for col in abundance_data.columns if col not in taxonomy_cols]
    
    # Mapeamento das amostras
    sample_mapping = get_sample_mapping(metadata)
    
    # Separa amostras por condição
    condition_a_samples = []
    condition_b_samples = []
    
    for sample_id in sample_cols:
        if sample_id in sample_mapping:
            if sample_mapping[sample_id]['condition'] == 'A':
                condition_a_samples.append(sample_id)
            elif sample_mapping[sample_id]['condition'] == 'B':
                condition_b_samples.append(sample_id)
    
    # Obtém grupos únicos
    groups = metadata['group'].unique()
    
    # Determina teste conservador
    conservative_test = determine_conservative_test(abundance_data, metadata)
    
    results = []
    sup_stats = []
    
    print(f"\nCalculando estatísticas abrangentes para {len(abundance_data)} espécies...")
    
    for idx, (_, row) in enumerate(abundance_data.iterrows()):
        if idx > 0 and idx % 100 == 0:
            print(f"  Processadas {idx} espécies...")
        
        species_name = row['Species']
        
        species_data = {
            'Phylum': row['Phylum'],
            'Class': row['Class'],
            'Order': row['Order'],
            'Family': row['Family'],
            'Genus': row['Genus'],
            'Species': species_name
        }
        
        # Dados para sup_statistics
        sup_data = {
            'Species': species_name,
            'Phylum': row['Phylum'],
            'Class': row['Class'],
            'Order': row['Order'],
            'Family': row['Family'],
            'Genus': row['Genus'],
            'Conservative_Test_Type': conservative_test
        }
        
        # Adiciona valores das amostras
        for sample_id in sample_cols:
            if sample_id in abundance_data.columns:
                species_data[sample_id] = row[sample_id]
        
        # Calcula estatísticas globais
        all_values = pd.to_numeric(row[sample_cols], errors='coerce').fillna(0)
        condition_a_values = pd.to_numeric(row[condition_a_samples], errors='coerce').fillna(0)
        condition_b_values = pd.to_numeric(row[condition_b_samples], errors='coerce').fillna(0)
        
        # Médias e desvios padrão
        species_data['MEAN_TOTAL'] = all_values.mean()
        species_data['Mean_before'] = condition_a_values.mean()
        species_data['Mean_after'] = condition_b_values.mean()
        species_data['SD_bef'] = condition_a_values.std()
        species_data['SD_after'] = condition_b_values.std()
        
        # Adiciona ao sup_statistics
        sup_data['Global_Mean_A'] = condition_a_values.mean()
        sup_data['Global_SD_A'] = condition_a_values.std()
        sup_data['Global_Mean_B'] = condition_b_values.mean()
        sup_data['Global_SD_B'] = condition_b_values.std()
        sup_data['Global_Mean_Total'] = all_values.mean()
        
        # Teste global - INDIVIDUAL (melhor teste escolhido)
        test_name, p_value, details = choose_best_test(condition_a_values, condition_b_values)
        species_data['Global-t-test'] = p_value
        
        # Teste global - CONSERVADOR
        conservative_p = apply_conservative_test(condition_a_values, condition_b_values, conservative_test)
        species_data['Global-t-test-conservative'] = conservative_p
        
        # Adiciona detalhes ao sup_statistics
        sup_data['Global_Best_Test_Type'] = test_name
        sup_data['Global_Best_P_Value'] = p_value
        sup_data['Global_Conservative_P_Value'] = conservative_p
        sup_data['Global_Shapiro_A'] = details.get('shapiro_p1', np.nan)
        sup_data['Global_Shapiro_B'] = details.get('shapiro_p2', np.nan)
        sup_data['Global_Normal_A'] = details.get('normal1', False)
        sup_data['Global_Normal_B'] = details.get('normal2', False)
        sup_data['Global_Levene_P'] = details.get('levene_p', np.nan)
        sup_data['Global_Equal_Vars'] = details.get('equal_vars', False)
        
        # Testes por grupo
        significant_tests_individual = 0
        significant_tests_conservative = 0
        
        for group in groups:
            # Amostras deste grupo
            group_a_samples = [s for s in sample_cols if s in sample_mapping and 
                             sample_mapping[s]['group'] == group and sample_mapping[s]['condition'] == 'A']
            group_b_samples = [s for s in sample_cols if s in sample_mapping and 
                             sample_mapping[s]['group'] == group and sample_mapping[s]['condition'] == 'B']
            
            if group_a_samples and group_b_samples:
                group_a_values = pd.to_numeric(row[group_a_samples], errors='coerce').fillna(0)
                group_b_values = pd.to_numeric(row[group_b_samples], errors='coerce').fillna(0)
                
                # Teste individual (melhor para esta espécie/grupo)
                test_name, p_value, details = choose_best_test(group_a_values, group_b_values)
                species_data[f'{group}t-test'] = p_value
                if p_value < 0.05:
                    significant_tests_individual += 1
                
                # Teste conservador
                conservative_p = apply_conservative_test(group_a_values, group_b_values, conservative_test)
                species_data[f'{group}t-test-conservative'] = conservative_p
                if conservative_p < 0.05:
                    significant_tests_conservative += 1
                
                # Adiciona ao sup_statistics
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
                
                # Preenche sup_statistics com valores padrão
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
        
        # Determina a cor baseada nos testes INDIVIDUAIS
        global_significant = species_data['Global-t-test'] < 0.05
        
        if global_significant:
            species_data['Cor'] = 'vermelho'
        elif significant_tests_individual > 1:
            species_data['Cor'] = 'laranja'
        elif significant_tests_individual == 1:
            species_data['Cor'] = 'amarelo'
        else:
            species_data['Cor'] = 'cinza'
        
        # Determina a cor baseada nos testes CONSERVADORES
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
    
    print(f"✓ Estatísticas calculadas para todas as {len(results)} espécies")
    return pd.DataFrame(results), pd.DataFrame(sup_stats)

def apply_filters(stats_df, parameter_results, top30_species, use_conservative=False):
    """
    Aplica os filtros para gerar o arquivo filtrado.csv.
    """
    filter_type = "conservador" if use_conservative else "individual"
    print(f"Aplicando filtros com testes {filter_type}...")
    
    # Merge com os resultados do parameter weighted score
    if 'Parameter_Weighted_Score' not in stats_df.columns:
        merge_cols = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        stats_df = stats_df.merge(
            parameter_results[merge_cols + ['Parameter_Weighted_Score']],
            on=merge_cols,
            how='left'
        )
        stats_df['Parameter_Weighted_Score'] = stats_df['Parameter_Weighted_Score'].fillna(0)
    
    print(f"Total de espécies antes dos filtros: {len(stats_df)}")
    
    # Filtro 1: Média total >= 0.005 (0.5%)
    filter1 = stats_df['MEAN_TOTAL'] >= 0.005
    print(f"Espécies com média total >= 0.5%: {filter1.sum()}")
    
    # Filtro 2: Pelo menos um teste significativo OU estar no TOP30
    if use_conservative:
        # Para teste conservador, procura por colunas que terminavam com 't-test-conservative'
        # mas que agora foram renomeadas para 't-test'
        t_test_cols_for_filter = [col for col in stats_df.columns if col.endswith('t-test-conservative')]
    else:
        # Para teste individual, usa colunas que terminam com 't-test' mas não 't-test-conservative'
        t_test_cols_for_filter = [col for col in stats_df.columns if 
                                 col.endswith('t-test') and not col.endswith('t-test-conservative')]
    
    significant_t_tests = (stats_df[t_test_cols_for_filter] < 0.05).any(axis=1)
    in_top30 = stats_df['Species'].isin(top30_species)
    
    filter2 = significant_t_tests | in_top30
    print(f"Espécies com teste t significativo: {significant_t_tests.sum()}")
    print(f"Espécies no TOP30: {in_top30.sum()}")
    print(f"Espécies que atendem ao critério estatístico/TOP30: {filter2.sum()}")
    
    # Aplica ambos os filtros
    final_filter = filter1 & filter2
    filtered_df = stats_df[final_filter].copy()
    
    print(f"✓ Espécies após aplicação de todos os filtros: {len(filtered_df)}")
    
    # Seleciona colunas apropriadas baseado no tipo de teste
    if use_conservative:
        # Para arquivo conservador, usa as colunas conservadoras
        # Remove colunas de testes individuais (que não são conservadores)
        individual_test_cols = [col for col in filtered_df.columns if 
                               col.endswith('t-test') and not col.endswith('t-test-conservative')]
        if individual_test_cols:
            filtered_df = filtered_df.drop(columns=individual_test_cols)
        
        # Remove a coluna Cor (individual) e mantém Color (conservador)
        if 'Cor' in filtered_df.columns:
            filtered_df = filtered_df.drop(columns=['Cor'])
        
        # Renomeia colunas de teste conservador para formato padrão
        conservative_test_cols = [col for col in filtered_df.columns if col.endswith('t-test-conservative')]
        rename_dict = {}
        for col in conservative_test_cols:
            new_name = col.replace('t-test-conservative', 't-test')
            rename_dict[col] = new_name
        if rename_dict:
            filtered_df = filtered_df.rename(columns=rename_dict)
            
        # Define cor_col para usar nas colunas ordenadas
        cor_col = 'Color'
            
    else:
        # Para arquivo individual, usa as colunas individuais
        # Remove todas as colunas conservadoras
        conservative_cols = [col for col in filtered_df.columns if 'conservative' in col.lower() or col == 'Color']
        if conservative_cols:
            filtered_df = filtered_df.drop(columns=conservative_cols)
            
        # Define cor_col para usar nas colunas ordenadas  
        cor_col = 'Cor'
    
    # Ordena as colunas na ordem desejada
    taxonomy_cols = [cor_col, 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Parameter_Weighted_Score']
    
    # Identifica colunas de teste t (agora já padronizadas)
    t_test_cols = [col for col in filtered_df.columns if col.endswith('t-test')]
    
    sample_cols = [col for col in filtered_df.columns if col not in taxonomy_cols + 
                   ['MEAN_TOTAL', 'Mean_before', 'Mean_after', 'SD_bef', 'SD_after'] + t_test_cols]
    stats_cols = ['MEAN_TOTAL', 'Mean_before', 'Mean_after', 'SD_bef', 'SD_after']
    
    # Ordem final das colunas
    column_order = taxonomy_cols + sample_cols + t_test_cols + stats_cols
    
    # Reordena as colunas (apenas as que existem)
    existing_columns = [col for col in column_order if col in filtered_df.columns]
    filtered_df = filtered_df[existing_columns]
    
    return filtered_df

def generate_summary_report(filtered_df_individual, filtered_df_conservative, sup_stats_df, original_count):
    """
    Gera relatório resumo da filtragem incluindo comparação entre métodos.
    """
    report_lines = []
    
    report_lines.append("=" * 80)
    report_lines.append("RELATÓRIO DE FILTRAGEM E ANÁLISE ESTATÍSTICA HÍBRIDA")
    report_lines.append("=" * 80)
    
    report_lines.append(f"Espécies no dataset original: {original_count}")
    report_lines.append(f"Espécies após filtragem (testes individuais): {len(filtered_df_individual)}")
    report_lines.append(f"Espécies após filtragem (teste conservador): {len(filtered_df_conservative)}")
    report_lines.append(f"Percentual mantido (individual): {len(filtered_df_individual)/original_count*100:.1f}%")
    report_lines.append(f"Percentual mantido (conservador): {len(filtered_df_conservative)/original_count*100:.1f}%")
    
    # Teste conservador utilizado
    conservative_test = sup_stats_df['Conservative_Test_Type'].iloc[0]
    report_lines.append(f"\nTeste conservador aplicado: {conservative_test}")
    
    # Estatísticas por cor - Individual
    color_counts_ind = filtered_df_individual['Cor'].value_counts()
    report_lines.append(f"\nDistribuição por significância (TESTES INDIVIDUAIS):")
    for color in ['vermelho', 'laranja', 'amarelo', 'cinza']:
        count = color_counts_ind.get(color, 0)
        report_lines.append(f"  {color.capitalize()}: {count} espécies")
    
    # Estatísticas por cor - Conservador
    color_counts_cons = filtered_df_conservative['Color'].value_counts()
    report_lines.append(f"\nDistribuição por significância (TESTE CONSERVADOR):")
    for color in ['vermelho', 'laranja', 'amarelo', 'cinza']:
        count = color_counts_cons.get(color, 0)
        report_lines.append(f"  {color.capitalize()}: {count} espécies")
    
    # Resumo dos tipos de testes utilizados
    test_types = sup_stats_df['Global_Best_Test_Type'].value_counts()
    report_lines.append(f"\nTipos de testes individuais utilizados (Global):")
    for test_type, count in test_types.items():
        report_lines.append(f"  {test_type}: {count} espécies")
    
    # Estatísticas de normalidade
    normal_both = sup_stats_df['Global_Normal_A'] & sup_stats_df['Global_Normal_B']
    report_lines.append(f"\nTeste de normalidade (Global):")
    report_lines.append(f"  Ambas condições normais: {normal_both.sum()} espécies")
    report_lines.append(f"  Pelo menos uma não normal: {(~normal_both).sum()} espécies")
    
    # Comparação entre métodos
    individual_species = set(filtered_df_individual['Species'])
    conservative_species = set(filtered_df_conservative['Species'])
    
    only_individual = individual_species - conservative_species
    only_conservative = conservative_species - individual_species
    both_methods = individual_species & conservative_species
    
    report_lines.append(f"\nComparação entre métodos:")
    report_lines.append(f"  Espécies detectadas por ambos: {len(both_methods)}")
    report_lines.append(f"  Apenas por testes individuais: {len(only_individual)}")
    report_lines.append(f"  Apenas por teste conservador: {len(only_conservative)}")
    
    # Salva o relatório
    with open('filtering_report.txt', 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    # Mostra na tela
    for line in report_lines:
        print(line)

def main():
    print("Iniciando análise de filtragem estatística híbrida...")
    print("=" * 80)
    
    # 1. Carrega todos os dados necessários
    abundance_data, metadata, parameter_results, top30_species = load_data()
    
    if any(data is None for data in [abundance_data, metadata, parameter_results]):
        print("❌ Erro: Não foi possível carregar todos os arquivos necessários!")
        return
    
    original_count = len(abundance_data)
    
    # 2. Calcula estatísticas completas com testes individuais e conservadores
    print(f"\nCalculando estatísticas híbridas (individuais + conservador)...")
    stats_df, sup_stats_df = calculate_comprehensive_stats(abundance_data, metadata)
    
    # 3. Salva arquivo de estatísticas suplementares
    sup_stats_file = 'sup_statistics.csv'
    sup_stats_df.to_csv(sup_stats_file, index=False)
    print(f"✓ Estatísticas suplementares salvas: {sup_stats_file}")
    
    # 4. Aplica filtros para ambos os métodos
    print(f"\nAplicando filtros com testes individuais...")
    filtered_df_individual = apply_filters(stats_df, parameter_results, top30_species, use_conservative=False)
    
    print(f"\nAplicando filtros com teste conservador...")
    filtered_df_conservative = apply_filters(stats_df, parameter_results, top30_species, use_conservative=True)
    
    # 5. Salva os arquivos filtrados
    output_file_individual = 'filtrado.csv'
    output_file_conservative = 'filtrado_conservativetest.csv'
    
    filtered_df_individual.to_csv(output_file_individual, index=False)
    filtered_df_conservative.to_csv(output_file_conservative, index=False)
    
    print(f"✓ Arquivo filtrado (individual) salvo: {output_file_individual}")
    print(f"✓ Arquivo filtrado (conservador) salvo: {output_file_conservative}")
    
    # 6. Gera relatório resumo comparativo
    generate_summary_report(filtered_df_individual, filtered_df_conservative, sup_stats_df, original_count)
    print(f"✓ Relatório comparativo salvo: filtering_report.txt")
    
    print(f"\n{'='*80}")
    print("ANÁLISE HÍBRIDA CONCLUÍDA!")
    print(f"Arquivos gerados:")
    print(f"  - {output_file_individual} ({len(filtered_df_individual)} espécies - testes individuais)")
    print(f"  - {output_file_conservative} ({len(filtered_df_conservative)} espécies - teste conservador)")
    print(f"  - {sup_stats_file} (estatísticas detalhadas de todas as espécies)")
    print(f"  - filtering_report.txt (relatório comparativo detalhado)")
    print("="*80)

if __name__ == '__main__':
    main()
