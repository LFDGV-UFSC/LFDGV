#!/usr/bin/env python3
"""
Main Pipeline for Bacterial Correlation Analysis
Executes all scripts in sequence and organizes outputs.
"""

import subprocess
import sys
import os
import shutil
import glob
from datetime import datetime

def print_header(message):
    """Prints formatted header."""
    print("\n" + "="*80)
    print(f" {message}")
    print("="*80)

def print_step(step_num, message):
    """Prints step number and description."""
    print(f"\n[STEP {step_num}] {message}")
    print("-" * (len(message) + 12))

def run_script(script_name, description):
    """
    Executes a Python script and waits for completion.
    Returns True if success, False if error.
    """
    print(f"🚀 Executing: {script_name}")
    print(f"📋 Description: {description}")
    
    start_time = datetime.now()
    
    try:
        # Execute script and wait for completion
        result = subprocess.run([sys.executable, script_name], 
                              capture_output=False,  # Show output in real time
                              text=True, 
                              check=True)
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        print(f"✅ {script_name} completed successfully!")
        print(f"⏱️  Execution time: {duration:.1f} seconds")
        
        return True
        
    except subprocess.CalledProcessError as e:
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        print(f"❌ Error executing {script_name}")
        print(f"⏱️  Time until error: {duration:.1f} seconds")
        print(f"🔥 Error code: {e.returncode}")
        
        return False
        
    except FileNotFoundError:
        print(f"❌ File {script_name} not found!")
        return False

def check_required_files():
    """Checks if required files exist."""
    required_files = [
        'abundance_data.csv',
        'metadata.tsv'
    ]
    
    missing_files = []
    for file in required_files:
        if not os.path.exists(file):
            # Check alternatives for metadata
            if file == 'metadata.tsv':
                alternatives = ['metadata.txt', 'metadata.csv']
                found_alternative = False
                for alt in alternatives:
                    if os.path.exists(alt):
                        found_alternative = True
                        break
                if not found_alternative:
                    missing_files.append(file + ' (or metadata.txt/metadata.csv)')
            else:
                missing_files.append(file)
    
    return missing_files

def check_required_scripts():
    """Checks if required scripts exist."""
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
    """Creates output directory."""
    output_dir = 'out_pwd'
    
    if os.path.exists(output_dir):
        print(f"📁 Directory {output_dir} already exists")
        # Ask if want to overwrite
        response = input(f"Do you want to clean the {output_dir} directory? (y/n): ").lower().strip()
        if response in ['y', 'yes', 's', 'sim']:
            shutil.rmtree(output_dir)
            os.makedirs(output_dir)
            print(f"🗑️  Directory {output_dir} cleaned and recreated")
        else:
            print(f"📂 Using existing {output_dir} directory")
    else:
        os.makedirs(output_dir)
        print(f"📁 Directory {output_dir} created")

def move_output_files():
    """Moves all output files to out_pwd directory."""
    output_dir = 'out_pwd'
    
    # List of expected output files
    output_files = [
        # Outputs from pws.py
        'parameter_weighted_analysis.csv',
        'top30.txt',
        
        # Outputs from filtering5.py
        'filtrado.csv',
        'filtrado_conservativetest.csv',
        'sup_statistics.csv',
        'filtering_report.txt',
        
        # Outputs from plotting_script.py (dynamic)
        '*_weighted_barplot.png',
        '*_weighted_barplot2.png'
    ]
    
    moved_files = []
    missing_files = []
    
    print("📦 Moving output files...")
    
    for pattern in output_files:
        if '*' in pattern:
            # Use glob for wildcard patterns
            matching_files = glob.glob(pattern)
            for file in matching_files:
                if os.path.exists(file):
                    dest_path = os.path.join(output_dir, os.path.basename(file))
                    shutil.move(file, dest_path)
                    moved_files.append(file)
                    print(f"  ✅ {file} → {output_dir}/")
        else:
            # Specific file
            if os.path.exists(pattern):
                dest_path = os.path.join(output_dir, pattern)
                shutil.move(pattern, dest_path)
                moved_files.append(pattern)
                print(f"  ✅ {pattern} → {output_dir}/")
            else:
                missing_files.append(pattern)
    
    print(f"\n📊 TRANSFER SUMMARY:")
    print(f"   ✅ Files moved: {len(moved_files)}")
    if missing_files:
        print(f"   ⚠️  Files not found: {len(missing_files)}")
        for file in missing_files:
            print(f"      - {file}")
    
    # Create log file with timestamp
    log_file = os.path.join(output_dir, 'pipeline_log.txt')
    with open(log_file, 'w') as f:
        f.write(f"Pipeline executed on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Generated files:\n")
        for file in moved_files:
            f.write(f"  - {file}\n")
        if missing_files:
            f.write(f"\nExpected but not found files:\n")
            for file in missing_files:
                f.write(f"  - {file}\n")
    
    print(f"📝 Log saved in: {log_file}")

def main():
    """Main pipeline function."""
    
    print_header("BACTERIAL ANALYSIS PIPELINE - PWS")
    print("🧬 Automated execution of analysis scripts")
    print("📅 Started at:", datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    
    # Step 0: Initial checks
    print_step(0, "INITIAL CHECKS")
    
    # Check input files
    missing_files = check_required_files()
    if missing_files:
        print("❌ Missing input files:")
        for file in missing_files:
            print(f"   - {file}")
        print("\nPlease ensure the files are in the current directory.")
        return False
    
    print("✅ All input files found")
    
    # Check scripts
    missing_scripts = check_required_scripts()
    if missing_scripts:
        print("❌ Missing scripts:")
        for script in missing_scripts:
            print(f"   - {script}")
        print("\nPlease ensure the scripts are in the current directory.")
        return False
    
    print("✅ All scripts found")
    
    # Create output directory
    create_output_directory()
    
    # Step 1: Main Analysis (PWS)
    print_step(1, "MAIN ANALYSIS - PWS")
    success = run_script('pws.py', 'Correlation analysis with external parameters')
    
    if not success:
        print("💥 Pipeline interrupted due to error in step 1")
        return False
    
    # Step 2: Statistical Filtering
    print_step(2, "ADVANCED STATISTICAL FILTERING")
    success = run_script('filtering5.py', 'Filtering with individual and conservative statistical tests')
    
    if not success:
        print("💥 Pipeline interrupted due to error in step 2")
        return False
    
    # Step 3: Plot Generation
    print_step(3, "PLOT GENERATION")
    success = run_script('plotting_script.py', 'Creation of comparative bar charts')
    
    if not success:
        print("💥 Pipeline interrupted due to error in step 3")
        return False
    
    # Step 4: Output Organization
    print_step(4, "OUTPUT ORGANIZATION")
    move_output_files()
    
    # Final summary
    print_header("PIPELINE COMPLETED SUCCESSFULLY! 🎉")
    
    print("📂 All output files have been organized in: out_pwd/")
    print("📋 Summary of generated files:")
    print("   📊 parameter_weighted_analysis.csv - Complete analysis results")
    print("   📝 top30.txt - Top 30 most impactful species")
    print("   🔬 filtrado.csv - Filtered species (individual tests)")
    print("   🔬 filtrado_conservativetest.csv - Filtered species (conservative test)")
    print("   📈 sup_statistics.csv - Detailed statistics")
    print("   📋 filtering_report.txt - Filtering report")
    print("   🎨 *_weighted_barplot.png - Impact plots")
    print("   📝 pipeline_log.txt - Execution log")
    
    print(f"\n🕐 Pipeline finished at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("🎯 Ready for results analysis!")
    
    return True

if __name__ == '__main__':
    try:
        success = main()
        if success:
            sys.exit(0)  # Success
        else:
            sys.exit(1)  # Error
    except KeyboardInterrupt:
        print("\n\n🛑 Pipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n💥 Unexpected error: {e}")
        sys.exit(1)
