#!/usr/bin/env python3
"""
Pipeline Principal para Análise de Correlação Bacteriana
Executa todos os scripts em sequência e organiza as saídas.
"""

import subprocess
import sys
import os
import shutil
import glob
from datetime import datetime

def print_header(message):
    """Imprime cabeçalho formatado."""
    print("\n" + "="*80)
    print(f" {message}")
    print("="*80)

def print_step(step_num, message):
    """Imprime número e descrição do passo."""
    print(f"\n[PASSO {step_num}] {message}")
    print("-" * (len(message) + 15))

def run_script(script_name, description):
    """
    Executa um script Python e aguarda sua conclusão.
    Retorna True se sucesso, False se erro.
    """
    print(f"🚀 Executando: {script_name}")
    print(f"📋 Descrição: {description}")
    
    start_time = datetime.now()
    
    try:
        # Executa o script e aguarda conclusão
        result = subprocess.run([sys.executable, script_name], 
                              capture_output=False,  # Mostra output em tempo real
                              text=True, 
                              check=True)
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        print(f"✅ {script_name} concluído com sucesso!")
        print(f"⏱️  Tempo de execução: {duration:.1f} segundos")
        
        return True
        
    except subprocess.CalledProcessError as e:
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        print(f"❌ Erro ao executar {script_name}")
        print(f"⏱️  Tempo até erro: {duration:.1f} segundos")
        print(f"🔥 Código de erro: {e.returncode}")
        
        return False
        
    except FileNotFoundError:
        print(f"❌ Arquivo {script_name} não encontrado!")
        return False

def check_required_files():
    """Verifica se os arquivos necessários existem."""
    required_files = [
        'abundance_data.csv',
        'metadata.tsv'
    ]
    
    missing_files = []
    for file in required_files:
        if not os.path.exists(file):
            # Verifica alternativas para metadata
            if file == 'metadata.tsv':
                alternatives = ['metadata.txt', 'metadata.csv']
                found_alternative = False
                for alt in alternatives:
                    if os.path.exists(alt):
                        found_alternative = True
                        break
                if not found_alternative:
                    missing_files.append(file + ' (ou metadata.txt/metadata.csv)')
            else:
                missing_files.append(file)
    
    return missing_files

def check_required_scripts():
    """Verifica se os scripts necessários existem."""
    required_scripts = [
        'pws.py',
        'filtering5.py', 
        'plotting_script.py'
    ]
    
    missing_scripts = []
    for script in required_scripts:
        if not os.path.exists(script):
            missing_scripts.append(script)
    
    return missing_scripts

def create_output_directory():
    """Cria o diretório de saída."""
    output_dir = 'out_pwd'
    
    if os.path.exists(output_dir):
        print(f"📁 Diretório {output_dir} já existe")
        # Pergunta se quer sobrescrever
        response = input(f"Deseja limpar o diretório {output_dir}? (s/n): ").lower().strip()
        if response in ['s', 'sim', 'y', 'yes']:
            shutil.rmtree(output_dir)
            os.makedirs(output_dir)
            print(f"🗑️  Diretório {output_dir} limpo e recriado")
        else:
            print(f"📂 Usando diretório {output_dir} existente")
    else:
        os.makedirs(output_dir)
        print(f"📁 Diretório {output_dir} criado")

def move_output_files():
    """Move todos os arquivos de saída para o diretório out_pwd."""
    output_dir = 'out_pwd'
    
    # Lista de arquivos de saída esperados
    output_files = [
        # Saídas do pws.py
        'parameter_weighted_analysis.csv',
        'top30.txt',
        
        # Saídas do filtering5.py
        'filtrado.csv',
        'filtrado_conservativetest.csv',
        'sup_statistics.csv',
        'filtering_report.txt',
        
        # Saídas do plotting_script.py (dinâmicas)
        '*_weighted_barplot.png',
        '*_weighted_barplot2.png'
    ]
    
    moved_files = []
    missing_files = []
    
    print("📦 Movendo arquivos de saída...")
    
    for pattern in output_files:
        if '*' in pattern:
            # Usa glob para padrões com wildcard
            matching_files = glob.glob(pattern)
            for file in matching_files:
                if os.path.exists(file):
                    dest_path = os.path.join(output_dir, os.path.basename(file))
                    shutil.move(file, dest_path)
                    moved_files.append(file)
                    print(f"  ✅ {file} → {output_dir}/")
        else:
            # Arquivo específico
            if os.path.exists(pattern):
                dest_path = os.path.join(output_dir, pattern)
                shutil.move(pattern, dest_path)
                moved_files.append(pattern)
                print(f"  ✅ {pattern} → {output_dir}/")
            else:
                missing_files.append(pattern)
    
    print(f"\n📊 RESUMO DA TRANSFERÊNCIA:")
    print(f"   ✅ Arquivos movidos: {len(moved_files)}")
    if missing_files:
        print(f"   ⚠️  Arquivos não encontrados: {len(missing_files)}")
        for file in missing_files:
            print(f"      - {file}")
    
    # Cria um arquivo de log com timestamp
    log_file = os.path.join(output_dir, 'pipeline_log.txt')
    with open(log_file, 'w') as f:
        f.write(f"Pipeline executado em: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Arquivos gerados:\n")
        for file in moved_files:
            f.write(f"  - {file}\n")
        if missing_files:
            f.write(f"\nArquivos esperados mas não encontrados:\n")
            for file in missing_files:
                f.write(f"  - {file}\n")
    
    print(f"📝 Log salvo em: {log_file}")

def main():
    """Função principal do pipeline."""
    
    print_header("PIPELINE DE ANÁLISE BACTERIANA - PWS")
    print("🧬 Execução automatizada dos scripts de análise")
    print("📅 Iniciado em:", datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    
    # Passo 0: Verificações iniciais
    print_step(0, "VERIFICAÇÕES INICIAIS")
    
    # Verifica arquivos de entrada
    missing_files = check_required_files()
    if missing_files:
        print("❌ Arquivos de entrada ausentes:")
        for file in missing_files:
            print(f"   - {file}")
        print("\nPor favor, certifique-se de que os arquivos estão no diretório atual.")
        return False
    
    print("✅ Todos os arquivos de entrada encontrados")
    
    # Verifica scripts
    missing_scripts = check_required_scripts()
    if missing_scripts:
        print("❌ Scripts ausentes:")
        for script in missing_scripts:
            print(f"   - {script}")
        print("\nPor favor, certifique-se de que os scripts estão no diretório atual.")
        return False
    
    print("✅ Todos os scripts encontrados")
    
    # Cria diretório de saída
    create_output_directory()
    
    # Passo 1: Análise Principal (PWS)
    print_step(1, "ANÁLISE PRINCIPAL - PWS")
    success = run_script('pws.py', 'Análise de correlação com parâmetros externos')
    
    if not success:
        print("💥 Pipeline interrompido devido a erro no passo 1")
        return False
    
    # Passo 2: Filtragem Estatística
    print_step(2, "FILTRAGEM ESTATÍSTICA AVANÇADA")
    success = run_script('filtering5.py', 'Filtragem com testes estatísticos individuais e conservadores')
    
    if not success:
        print("💥 Pipeline interrompido devido a erro no passo 2")
        return False
    
    # Passo 3: Geração de Gráficos
    print_step(3, "GERAÇÃO DE GRÁFICOS")
    success = run_script('plotting_script.py', 'Criação de gráficos de barras comparativos')
    
    if not success:
        print("💥 Pipeline interrompido devido a erro no passo 3")
        return False
    
    # Passo 4: Organização das Saídas
    print_step(4, "ORGANIZAÇÃO DAS SAÍDAS")
    move_output_files()
    
    # Resumo final
    print_header("PIPELINE CONCLUÍDO COM SUCESSO! 🎉")
    
    print("📂 Todos os arquivos de saída foram organizados em: out_pwd/")
    print("📋 Resumo dos arquivos gerados:")
    print("   📊 parameter_weighted_analysis.csv - Resultados completos da análise")
    print("   📝 top30.txt - Top 30 espécies mais impactantes")
    print("   🔬 filtrado.csv - Espécies filtradas (testes individuais)")
    print("   🔬 filtrado_conservativetest.csv - Espécies filtradas (teste conservador)")
    print("   📈 sup_statistics.csv - Estatísticas detalhadas")
    print("   📋 filtering_report.txt - Relatório de filtragem")
    print("   🎨 *_weighted_barplot.png - Gráficos de impacto")
    print("   📝 pipeline_log.txt - Log de execução")
    
    print(f"\n🕐 Pipeline finalizado em: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("🎯 Pronto para análise dos resultados!")
    
    return True

if __name__ == '__main__':
    try:
        success = main()
        if success:
            sys.exit(0)  # Sucesso
        else:
            sys.exit(1)  # Erro
    except KeyboardInterrupt:
        print("\n\n🛑 Pipeline interrompido pelo usuário")
        sys.exit(1)
    except Exception as e:
        print(f"\n💥 Erro inesperado: {e}")
        sys.exit(1)
