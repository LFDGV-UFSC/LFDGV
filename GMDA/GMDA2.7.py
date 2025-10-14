import tkinter as tk
from tkinter import filedialog, messagebox, ttk, scrolledtext
from Bio import SeqIO, Seq, Phylo
from Bio.SeqUtils import gc_fraction, MeltingTemp
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
import openpyxl
from openpyxl.utils import get_column_letter
from openpyxl.styles import Font
import re
import os
import math
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.cluster.hierarchy import cophenet
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pandas as pd
from sklearn.metrics import pairwise_distances
from collections import defaultdict
import threading
import time
from packaging import version
from Bio import __version__ as bio_version
from io import StringIO
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

try:
    from statistics import mode, multimode
except ImportError:
    from statistics import mode
    def multimode(data):
        from collections import Counter
        count = Counter(data)
        max_count = max(count.values())
        return [x for x, c in count.items() if c == max_count]

from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import string
from sklearn.manifold import MDS
import matplotlib.patches as mpatches
from scipy.spatial import procrustes
from scipy.stats import pearsonr
from scipy.stats import permutation_test
import warnings
from scipy.spatial.distance import cdist

# Configuração do tema moderno
MODERN_THEME = {
    "primary": "#2C3E50",
    "secondary": "#34495E",
    "accent": "#3498DB",
    "success": "#27AE60",
    "warning": "#F39C12",
    "danger": "#E74C3C",
    "light": "#ECF0F1",
    "dark": "#2C3E50",
    "text": "#2C3E50",
    "text_light": "#7F8C8D",
    "bg": "#F8F9FA",
    "card_bg": "#FFFFFF"
}

def mantel(matrix1, matrix2, method='pearson', permutations=999):
    from scipy.stats import pearsonr
    import numpy as np
    
    vec1 = matrix1.flatten()
    vec2 = matrix2.flatten()
    
    corr_original, _ = pearsonr(vec1, vec2)
    
    corr_permutations = []
    for _ in range(permutations):
        np.random.shuffle(vec2)
        corr_perm, _ = pearsonr(vec1, vec2)
        corr_permutations.append(corr_perm)
    
    p_value = np.mean(np.abs(corr_permutations) >= np.abs(corr_original))
    
    return corr_original, p_value, np.std(corr_permutations)

warnings.filterwarnings('ignore')

class ModernButton(tk.Button):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)
        self.configure(
            relief="flat",
            bd=0,
            padx=20,
            pady=10,
            font=("Segoe UI", 10),
            cursor="hand2"
        )

class ModernEntry(tk.Entry):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)
        self.configure(
            relief="flat",
            bd=1,
            highlightthickness=1,
            font=("Segoe UI", 10),
            highlightcolor=MODERN_THEME["accent"],
            highlightbackground="#CCCCCC"
        )

class ModernLabel(tk.Label):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)
        self.configure(
            font=("Segoe UI", 10),
            bg=MODERN_THEME["bg"],
            fg=MODERN_THEME["text"]
        )

class ModernFrame(tk.Frame):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)
        self.configure(bg=MODERN_THEME["bg"])

class CardFrame(tk.Frame):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)
        self.configure(
            bg=MODERN_THEME["card_bg"],
            relief="raised",
            bd=1,
            padx=15,
            pady=15
        )

class GMDA:
    def __init__(self, root):
        self.root = root
        self.root.title("GMDA - Genetic and Morphometric Data Analysis")
        self.root.geometry("1200x800")
        self.root.configure(bg=MODERN_THEME["bg"])
        
        # Configurar estilo moderno para ttk
        self.setup_styles()
        
        self.top3_results = []
        self.morpho_results = []
        self.anova_results = []
        self.tukey_results = []
        self.effect_size_results = []
        self.input_file_directory = ""
        self.replicates_var = tk.IntVar(value=3)
        
        self.genetic_linkage_matrix = None
        self.genetic_distance_matrix = None
        self.genetic_labels = None
        self.morpho_linkage_matrix = None
        self.morpho_distance_matrix = None
        self.morpho_labels = None
        self.genetic_pcoa_coords = None
        self.morpho_pca_coords = None
        
        self.setup_ui()
        self.setup_notebook()
    
    def setup_styles(self):
        style = ttk.Style()
        style.theme_use('clam')
        
        # Configurar cores do notebook
        style.configure("TNotebook", background=MODERN_THEME["bg"])
        style.configure("TNotebook.Tab", 
                       font=("Segoe UI", 10, "bold"),
                       padding=[20, 10],
                       background=MODERN_THEME["light"],
                       foreground=MODERN_THEME["text"])
        style.map("TNotebook.Tab", 
                 background=[("selected", MODERN_THEME["accent"])],
                 foreground=[("selected", "white")])
        
        # Configurar progressbar
        style.configure("TProgressbar",
                       thickness=20,
                       background=MODERN_THEME["success"])
    
    def setup_ui(self):
        # Header moderno
        header_frame = ModernFrame(self.root)
        header_frame.pack(fill="x", padx=20, pady=(20, 10))
        
        self.title_label = tk.Label(header_frame, 
                                  text="GMDA 2.7", 
                                  font=("Segoe UI", 28, "bold"), 
                                  bg=MODERN_THEME["bg"], 
                                  fg=MODERN_THEME["primary"])
        self.title_label.pack(pady=(0, 5))
        
        self.subtitle_label = tk.Label(header_frame, 
                                     text="Genetic and Morphometric Data Analysis Platform", 
                                     font=("Segoe UI", 12), 
                                     bg=MODERN_THEME["bg"], 
                                     fg=MODERN_THEME["text_light"])
        self.subtitle_label.pack(pady=(0, 10))
        
        # Status bar
        self.status_frame = ModernFrame(self.root)
        self.status_frame.pack(fill="x", side="bottom", padx=20, pady=5)
        
        self.status_label = tk.Label(self.status_frame, 
                                   text="Ready", 
                                   font=("Segoe UI", 9),
                                   bg=MODERN_THEME["bg"],
                                   fg=MODERN_THEME["text_light"])
        self.status_label.pack(side="left")
        
        self.version_label = tk.Label(self.status_frame,
                                    text="v2.0.0",
                                    font=("Segoe UI", 9),
                                    bg=MODERN_THEME["bg"],
                                    fg=MODERN_THEME["text_light"])
        self.version_label.pack(side="right")
    
    def update_status(self, message):
        self.status_label.config(text=message)
        self.root.update_idletasks()
    
    def setup_notebook(self):
        self.notebook = ttk.Notebook(self.root)
        
        self.genetic_frame = ModernFrame(self.notebook)
        self.morpho_frame = ModernFrame(self.notebook)
        self.correlation_frame = ModernFrame(self.notebook)
        
        self.notebook.add(self.genetic_frame, text="🧬 Genetic Analysis")
        self.notebook.add(self.morpho_frame, text="📊 Morphometric Analysis")
        self.notebook.add(self.correlation_frame, text="🔗 Correlation Analysis")
        
        self.notebook.pack(expand=1, fill="both", padx=20, pady=10)
        
        self.setup_genetic_tab()
        self.setup_morpho_tab()
        self.setup_correlation_tab()
    
    def setup_genetic_tab(self):
        # Main container
        main_container = ModernFrame(self.genetic_frame)
        main_container.pack(fill="both", expand=True, padx=10, pady=10)
        
        # Input section
        input_card = CardFrame(main_container)
        input_card.pack(fill="x", pady=(0, 15))
        
        input_title = ModernLabel(input_card, text="Input Files", font=("Segoe UI", 12, "bold"))
        input_title.pack(anchor="w", pady=(0, 15))
        
        # Database file input
        db_frame = ModernFrame(input_card)
        db_frame.pack(fill="x", pady=8)
        
        label_database = ModernLabel(db_frame, text="Reference Database:", font=("Segoe UI", 10, "bold"))
        label_database.grid(row=0, column=0, sticky="w", padx=(0, 10))
        
        self.entry_database = ModernEntry(db_frame, width=50)
        self.entry_database.grid(row=0, column=1, padx=5, sticky="ew")
        
        btn_browse_db = ModernButton(db_frame, text="Browse", 
                                   command=lambda: self.load_file(self.entry_database),
                                   bg=MODERN_THEME["accent"], fg="white")
        btn_browse_db.grid(row=0, column=2, padx=5)
        
        # Query file input
        query_frame = ModernFrame(input_card)
        query_frame.pack(fill="x", pady=8)
        
        label_query = ModernLabel(query_frame, text="Query Samples:", font=("Segoe UI", 10, "bold"))
        label_query.grid(row=0, column=0, sticky="w", padx=(0, 10))
        
        self.entry_query = ModernEntry(query_frame, width=50)
        self.entry_query.grid(row=0, column=1, padx=5, sticky="ew")
        
        btn_browse_query = ModernButton(query_frame, text="Browse", 
                                      command=lambda: self.load_file(self.entry_query),
                                      bg=MODERN_THEME["accent"], fg="white")
        btn_browse_query.grid(row=0, column=2, padx=5)
        
        # Options
        options_frame = ModernFrame(input_card)
        options_frame.pack(fill="x", pady=10)
        
        self.bootstrap_var_gen = tk.BooleanVar()
        check_bootstrap = tk.Checkbutton(options_frame, text="Use Bootstrap", 
                                       variable=self.bootstrap_var_gen, 
                                       bg=MODERN_THEME["card_bg"], 
                                       font=("Segoe UI", 10),
                                       selectcolor=MODERN_THEME["light"])
        check_bootstrap.pack(anchor="w")
        
        # Configure grid weights
        db_frame.columnconfigure(1, weight=1)
        query_frame.columnconfigure(1, weight=1)
        
        # Action buttons
        action_card = CardFrame(main_container)
        action_card.pack(fill="x", pady=(0, 15))
        
        action_title = ModernLabel(action_card, text="Analysis Tools", font=("Segoe UI", 12, "bold"))
        action_title.pack(anchor="w", pady=(0, 15))
        
        btn_frame = ModernFrame(action_card)
        btn_frame.pack(fill="x")
        
        btn_top3 = ModernButton(btn_frame, 
                              text="🧬 Generate Similarity Matrix & Dendrogram", 
                              command=self.calculate_genetic, 
                              bg=MODERN_THEME["success"], fg="white", width=30)
        btn_top3.grid(row=0, column=0, padx=5, pady=5, sticky="ew")
        
        btn_save = ModernButton(btn_frame, 
                              text="💾 Save Top 3 Results", 
                              command=self.save_genetic_results, 
                              bg=MODERN_THEME["accent"], fg="white", width=30)
        btn_save.grid(row=0, column=1, padx=5, pady=5, sticky="ew")
        
        btn_genetic_indices = ModernButton(btn_frame, 
                                         text="📈 Genetic Indices Analysis", 
                                         command=self.calculate_genetic_indices, 
                                         bg=MODERN_THEME["primary"], fg="white", width=30)
        btn_genetic_indices.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky="ew")
        
        btn_frame.columnconfigure(0, weight=1)
        btn_frame.columnconfigure(1, weight=1)
        
        # Results section
        results_card = CardFrame(main_container)
        results_card.pack(fill="both", expand=True)
        
        results_title = ModernLabel(results_card, text="Results", font=("Segoe UI", 12, "bold"))
        results_title.pack(anchor="w", pady=(0, 10))
        
        # Text area with scrollbar
        text_frame = ModernFrame(results_card)
        text_frame.pack(fill="both", expand=True)
        
        text_frame.columnconfigure(0, weight=1)
        text_frame.rowconfigure(0, weight=1)
        
        self.results_text_gen = tk.Text(text_frame, 
                                      bg="white", 
                                      fg=MODERN_THEME["text"],
                                      font=("Consolas", 10),
                                      relief="flat",
                                      padx=10,
                                      pady=10)
        
        scrollbar = ttk.Scrollbar(text_frame, orient="vertical", command=self.results_text_gen.yview)
        self.results_text_gen.configure(yscrollcommand=scrollbar.set)
        
        self.results_text_gen.grid(row=0, column=0, sticky="nsew")
        scrollbar.grid(row=0, column=1, sticky="ns")
    
    def setup_morpho_tab(self):
        # Main container
        main_container = ModernFrame(self.morpho_frame)
        main_container.pack(fill="both", expand=True, padx=10, pady=10)
        
        # Input section
        input_card = CardFrame(main_container)
        input_card.pack(fill="x", pady=(0, 15))
        
        input_title = ModernLabel(input_card, text="Morphometric Data Input", font=("Segoe UI", 12, "bold"))
        input_title.pack(anchor="w", pady=(0, 15))
        
        # File input
        file_frame = ModernFrame(input_card)
        file_frame.pack(fill="x", pady=8)
        
        label_file = ModernLabel(file_frame, text="Data File:", font=("Segoe UI", 10, "bold"))
        label_file.grid(row=0, column=0, sticky="w", padx=(0, 10))
        
        self.entry_file_morpho = ModernEntry(file_frame, width=50)
        self.entry_file_morpho.grid(row=0, column=1, padx=5, sticky="ew")
        
        btn_browse = ModernButton(file_frame, text="Browse", 
                                command=lambda: self.load_file(self.entry_file_morpho),
                                bg=MODERN_THEME["accent"], fg="white")
        btn_browse.grid(row=0, column=2, padx=5)
        
        # Parameters
        params_frame = ModernFrame(input_card)
        params_frame.pack(fill="x", pady=10)
        
        label_replicates = ModernLabel(params_frame, text="Replicates per sample:", font=("Segoe UI", 10))
        label_replicates.grid(row=0, column=0, sticky="w", padx=(0, 10))
        
        self.entry_replicates = ModernEntry(params_frame, textvariable=self.replicates_var, width=10)
        self.entry_replicates.grid(row=0, column=1, sticky="w", padx=5)
        
        self.bootstrap_var_morpho = tk.BooleanVar()
        check_bootstrap = tk.Checkbutton(params_frame, text="Use Bootstrap for Statistics", 
                                       variable=self.bootstrap_var_morpho, 
                                       bg=MODERN_THEME["card_bg"],
                                       font=("Segoe UI", 10),
                                       selectcolor=MODERN_THEME["light"])
        check_bootstrap.grid(row=0, column=2, sticky="w", padx=20)
        
        file_frame.columnconfigure(1, weight=1)
        
        # Action buttons
        action_card = CardFrame(main_container)
        action_card.pack(fill="x", pady=(0, 15))
        
        action_title = ModernLabel(action_card, text="Analysis Tools", font=("Segoe UI", 12, "bold"))
        action_title.pack(anchor="w", pady=(0, 15))
        
        btn_frame = ModernFrame(action_card)
        btn_frame.pack(fill="x")
        
        btn_analyze = ModernButton(btn_frame, 
                                 text="📊 Analyze Data", 
                                 command=self.calculate_morpho, 
                                 bg=MODERN_THEME["success"], fg="white", width=30)
        btn_analyze.grid(row=0, column=0, padx=5, pady=5, sticky="ew")
        
        btn_save = ModernButton(btn_frame, 
                              text="💾 Save Results", 
                              command=self.save_morpho_results, 
                              bg=MODERN_THEME["accent"], fg="white", width=30)
        btn_save.grid(row=0, column=1, padx=5, pady=5, sticky="ew")
        
        btn_frame.columnconfigure(0, weight=1)
        btn_frame.columnconfigure(1, weight=1)
        
        # Results section
        results_card = CardFrame(main_container)
        results_card.pack(fill="both", expand=True)
        
        results_title = ModernLabel(results_card, text="Analysis Results", font=("Segoe UI", 12, "bold"))
        results_title.pack(anchor="w", pady=(0, 10))
        
        # Text area with scrollbar
        text_frame = ModernFrame(results_card)
        text_frame.pack(fill="both", expand=True)
        
        text_frame.columnconfigure(0, weight=1)
        text_frame.rowconfigure(0, weight=1)
        
        self.results_text_morpho = tk.Text(text_frame, 
                                         bg="white", 
                                         fg=MODERN_THEME["text"],
                                         font=("Consolas", 10),
                                         relief="flat",
                                         padx=10,
                                         pady=10)
        
        scrollbar = ttk.Scrollbar(text_frame, orient="vertical", command=self.results_text_morpho.yview)
        self.results_text_morpho.configure(yscrollcommand=scrollbar.set)
        
        self.results_text_morpho.grid(row=0, column=0, sticky="nsew")
        scrollbar.grid(row=0, column=1, sticky="ns")
    
    def setup_correlation_tab(self):
        # Main container
        main_container = ModernFrame(self.correlation_frame)
        main_container.pack(fill="both", expand=True, padx=10, pady=10)
        
        # Info card
        info_card = CardFrame(main_container)
        info_card.pack(fill="x", pady=(0, 15))
        
        info_label = ModernLabel(info_card, 
                               text="💡 First run genetic and morphometric analyses to generate dendrograms and PCAs", 
                               font=("Segoe UI", 10),
                               fg=MODERN_THEME["warning"])
        info_label.pack(pady=10)
        
        # Analysis card
        analysis_card = CardFrame(main_container)
        analysis_card.pack(fill="x", pady=(0, 15))
        
        analysis_title = ModernLabel(analysis_card, text="Correlation Analysis", font=("Segoe UI", 12, "bold"))
        analysis_title.pack(anchor="w", pady=(0, 15))
        
        # Buttons
        btn_frame = ModernFrame(analysis_card)
        btn_frame.pack(fill="x", pady=10)
        
        btn_mantel_dendro = ModernButton(btn_frame, 
                                       text="🌳 Mantel Test - Dendrograms", 
                                       command=self.calculate_mantel_test_dendrograms, 
                                       bg=MODERN_THEME["primary"], fg="white", width=25)
        btn_mantel_dendro.grid(row=0, column=0, padx=5, pady=5, sticky="ew")
        
        btn_mantel_pca = ModernButton(btn_frame, 
                                    text="📈 Mantel Test - PCA/PCoA", 
                                    command=self.calculate_mantel_test_pca, 
                                    bg=MODERN_THEME["primary"], fg="white", width=25)
        btn_mantel_pca.grid(row=0, column=1, padx=5, pady=5, sticky="ew")
        
        btn_frame.columnconfigure(0, weight=1)
        btn_frame.columnconfigure(1, weight=1)
        
        # Progress bar
        self.progress_bar = ttk.Progressbar(analysis_card, mode='indeterminate')
        self.progress_bar.pack(fill="x", pady=10)
        
        # Results section
        results_card = CardFrame(main_container)
        results_card.pack(fill="both", expand=True)
        
        results_title = ModernLabel(results_card, text="Correlation Results", font=("Segoe UI", 12, "bold"))
        results_title.pack(anchor="w", pady=(0, 10))
        
        # Text area with scrollbar
        text_frame = ModernFrame(results_card)
        text_frame.pack(fill="both", expand=True)
        
        text_frame.columnconfigure(0, weight=1)
        text_frame.rowconfigure(0, weight=1)
        
        self.results_text_corr = tk.Text(text_frame, 
                                       bg="white", 
                                       fg=MODERN_THEME["text"],
                                       font=("Consolas", 10),
                                       relief="flat",
                                       padx=10,
                                       pady=10)
        
        scrollbar = ttk.Scrollbar(text_frame, orient="vertical", command=self.results_text_corr.yview)
        self.results_text_corr.configure(yscrollcommand=scrollbar.set)
        
        self.results_text_corr.grid(row=0, column=0, sticky="nsew")
        scrollbar.grid(row=0, column=1, sticky="ns")

    def load_file(self, entry_widget):
        filename = filedialog.askopenfilename(
            title="Select File",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        if filename:
            entry_widget.delete(0, tk.END)
            entry_widget.insert(0, filename)
            self.input_file_directory = os.path.dirname(filename)
            self.update_status(f"Loaded: {os.path.basename(filename)}")

    @staticmethod
    def isanumber(a):
        try:
            float(a)
            return not np.isnan(float(a))
        except (ValueError, TypeError):
            return False

    def calcular_similaridade(self, query, reference):
        total_similarity = 0
        valid_loci_count = 0

        for q, r in zip(query, reference):
            if self.isanumber(q) and self.isanumber(r) and not (np.isnan(float(q)) or np.isnan(float(r))):
                q = float(q)
                r = float(r)
                loco_similarity = max(0, 1 - abs(q - r) / r)
                total_similarity += loco_similarity
                valid_loci_count += 1

        return total_similarity / valid_loci_count if valid_loci_count > 0 else 0

    def calculate_genetic(self):
        self.update_status("Calculating genetic similarity...")
        self.top3_results = []
        db_filename = self.entry_database.get().strip()
        query_filename = self.entry_query.get().strip()

        if not db_filename or not query_filename:
            messagebox.showwarning("Warning", "Please select both files.")
            self.update_status("Ready")
            return

        try:
            db_df = pd.read_excel(db_filename, engine='openpyxl', header=1)
            query_df = pd.read_excel(query_filename, engine='openpyxl', header=1)
        except Exception as e:
            messagebox.showerror("Error", f"Error reading files:\n{e}")
            self.update_status("Error reading files")
            return

        if db_df.empty or query_df.empty:
            messagebox.showerror("Error", "One or both files are empty.")
            self.update_status("Empty files")
            return

        db_name_col = db_df.columns[0]
        query_name_col = query_df.columns[0]

        self.results_text_gen.delete(1.0, tk.END)
        self.results_text_gen.insert(tk.END, "🧬 GENETIC SIMILARITY ANALYSIS\n")
        self.results_text_gen.insert(tk.END, "=" * 50 + "\n\n")

        for idx_query, query_row in query_df.iterrows():
            query_name = query_row[query_name_col]
            if pd.isna(query_name) or str(query_name).strip() == '':
                continue

            query_data = query_row.drop(query_name_col)
            similarities = {}

            for idx_ref, ref_row in db_df.iterrows():
                ref_name = ref_row[db_name_col]
                if pd.isna(ref_name) or str(ref_name).strip() == '':
                    continue
                ref_data = ref_row.drop(db_name_col)
                similarity = self.calcular_similaridade(query_data, ref_data)
                similarities[ref_name] = similarity

            if not similarities:
                continue

            top3 = sorted(similarities.items(), key=lambda x: x[1], reverse=True)[:3]
            self.results_text_gen.insert(tk.END, f"🔍 Query sample: {query_name}\n")
            for ref_name, score in top3:
                self.results_text_gen.insert(tk.END, f"   📌 {ref_name}: {score * 100:.2f}%\n")
                self.top3_results.append((query_name, ref_name, score))

            self.top3_results.append(("", "", ""))

        try:
            self.generate_genetic_matrix_and_dendrogram(query_df, db_df, query_name_col, db_name_col, self.bootstrap_var_gen.get())
        except Exception as e:
            self.results_text_gen.insert(tk.END, f"\n❌ Error generating matrix/dendrogram: {e}\n")

        self.update_status("Genetic analysis completed")

    def generate_genetic_matrix_and_dendrogram(self, query_df, db_df, query_name_col, db_name_col, bootstrap_enabled):
        combined_df = pd.concat([query_df, db_df], ignore_index=True)
        combined_df = combined_df.dropna(subset=[query_name_col if query_name_col in combined_df.columns else db_name_col])
        
        if combined_df.empty:
            self.results_text_gen.insert(tk.END, "❌ Error: No valid samples found.\n")
            return

        names = combined_df.iloc[:, 0].tolist()
        data = combined_df.iloc[:, 1:]

        valid_indices = [i for i, name in enumerate(names) if pd.notna(name) and str(name).strip() != '']
        names = [names[i] for i in valid_indices]
        data = data.iloc[valid_indices]

        if len(names) < 2:
            self.results_text_gen.insert(tk.END, "❌ Error: At least 2 samples required.\n")
            return

        matrix = pd.DataFrame(index=names, columns=names, dtype=float)
        np.fill_diagonal(matrix.values, 1.0)

        for i in range(len(data)):
            for j in range(i + 1, len(data)):
                sim = self.calcular_similaridade(data.iloc[i], data.iloc[j])
                matrix.iat[i, j] = sim
                matrix.iat[j, i] = sim

        matrix_filename = os.path.join(self.input_file_directory, "Genetic_similarity_matrix.xlsx")
        try:
            with pd.ExcelWriter(matrix_filename, engine='openpyxl') as writer:
                title_df = pd.DataFrame({"": ["GMDA - Genetic Similarity Matrix"]})
                title_df.to_excel(writer, sheet_name='Similarity Matrix', index=False, header=False)
                matrix.to_excel(writer, sheet_name='Similarity Matrix', startrow=2)
            self.results_text_gen.insert(tk.END, f"\n💾 Similarity matrix saved as: {matrix_filename}\n")
        except Exception as e:
            self.results_text_gen.insert(tk.END, f"\n❌ Error saving matrix: {e}\n")

        distance_matrix = 1 - matrix.astype(float).values
        np.fill_diagonal(distance_matrix, 0)
        
        try:
            condensed_dist = squareform(distance_matrix, checks=False)
            linkage_matrix = linkage(condensed_dist, method='average')
            
            self.genetic_linkage_matrix = linkage_matrix
            self.genetic_distance_matrix = condensed_dist
            self.genetic_labels = names
            
            if bootstrap_enabled:
                self.results_text_gen.insert(tk.END, "🔄 Calculating bootstrap supports...\n")
                self.root.update_idletasks()
                supports = self.calculate_bootstrap_support(data, linkage_matrix)
                self.plot_genetic_dendrogram(linkage_matrix, names, supports)
            else:
                self.plot_genetic_dendrogram(linkage_matrix, names)
                
        except Exception as e:
            self.results_text_gen.insert(tk.END, f"❌ Error computing dendrogram: {e}\n")

    def calculate_bootstrap_support(self, data, linkage_matrix, replicates=100):
        n = data.shape[0]
        clusters_suporte = defaultdict(int)

        original_clusters = self.get_cluster_members(linkage_matrix, n)

        for _ in range(replicates):
            try:
                sampled_cols = np.random.choice(data.shape[1], size=data.shape[1], replace=True)
                data_boot = data.iloc[:, sampled_cols]

                matrix = pd.DataFrame(index=data_boot.index, columns=data_boot.index, dtype=float)
                np.fill_diagonal(matrix.values, 1.0)
                
                for i1 in range(len(data_boot)):
                    for i2 in range(i1 + 1, len(data_boot)):
                        sim = self.calcular_similaridade(data_boot.iloc[i1], data_boot.iloc[i2])
                        matrix.iat[i1, i2] = sim
                        matrix.iat[i2, i1] = sim

                dist = 1 - matrix.values
                condensed_dist = squareform(dist, checks=False)
                linkage_boot = linkage(condensed_dist, method='average')
                bootstrap_clusters = self.get_cluster_members(linkage_boot, n)

                for cluster_id, members in original_clusters.items():
                    if any(members == b_members for b_members in bootstrap_clusters.values()):
                        clusters_suporte[cluster_id] += 1

            except Exception:
                continue

        for c in clusters_suporte:
            clusters_suporte[c] = clusters_suporte[c] / replicates * 100

        return clusters_suporte

    def get_cluster_members(self, linkage_matrix, n_samples):
        clusters = {}
        for i, (c1, c2, _, _) in enumerate(linkage_matrix):
            c1, c2 = int(c1), int(c2)
            members = set()
            if c1 < n_samples:
                members.add(c1)
            else:
                members.update(clusters.get(c1, set()))
            if c2 < n_samples:
                members.add(c2)
            else:
                members.update(clusters.get(c2, set()))
            clusters[i + n_samples] = members
        return clusters

    def plot_genetic_dendrogram(self, linkage_matrix, labels, supports=None):
        plt.figure(figsize=(12, 8))
        dendro = dendrogram(linkage_matrix, labels=labels, leaf_rotation=90, leaf_font_size=10)

        if supports:
            icoord = np.array(dendro['icoord'])
            dcoord = np.array(dendro['dcoord'])
            
            for i, d in zip(icoord, dcoord):
                x = 0.5 * (i[1] + i[2])
                y = d[1]
                
                for cluster_idx, row in enumerate(linkage_matrix):
                    if abs(row[2] - y) < 1e-5:
                        support = supports.get(cluster_idx + len(labels), 0)
                        if support >= 50:
                            plt.text(x, y, f"{int(support)}%", 
                                   ha='center', va='bottom',
                                   fontsize=8, color='red')
                        break

        title = "GMDA - Genetic UPGMA Dendrogram"
        if supports:
            title += " (with bootstrap)"
            
        plt.title(title, fontsize=14)
        plt.ylabel("Dissimilarity (1-S)", fontsize=14)
        plt.tight_layout()
        
        filename = "Genetic_Dendrogram"
        if supports:
            filename += "_with_bootstrap"
            
        plot_path = os.path.join(self.input_file_directory, f"{filename}.png")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        self.results_text_gen.insert(tk.END, f"📊 Dendrogram saved as: {plot_path}\n")

    def save_genetic_results(self):
        if not self.top3_results:
            messagebox.showinfo("No data", "No results available to save.")
            return

        df_results = pd.DataFrame(self.top3_results, columns=["Query", "Reference", "Similarity"])
        df_results["Similarity (%)"] = pd.to_numeric(df_results["Similarity"], errors='coerce') * 100
        df_results.drop(columns=["Similarity"], inplace=True)

        save_path = filedialog.asksaveasfilename(
            initialdir=self.input_file_directory,
            defaultextension=".tsv",
            filetypes=[("TSV file", "*.tsv"), ("All files", "*.*")],
            title="Save genetic similarity results as..."
        )
        if not save_path:
            return

        try:
            with open(save_path, 'w', encoding='utf-8') as f:
                f.write("GMDA - Genetic Similarity Results\n\n")
                df_results.to_csv(f, sep="\t", index=False)
            self.results_text_gen.insert(tk.END, f"\n💾 Results saved as TSV: {save_path}\n")
        except Exception as e:
            self.results_text_gen.insert(tk.END, f"❌ Error saving TSV: {e}\n")
            return

        pdf_path = save_path.rsplit('.', 1)[0] + ".pdf"
        try:
            with PdfPages(pdf_path) as pdf:
                lines_per_page = 45
                lines = ["GMDA - Genetic Similarity Results", ""]
                
                for idx, row in df_results.iterrows():
                    if pd.isna(row["Query"]) or str(row["Query"]).strip() == "":
                        lines.append("")
                        continue
                    line = f"{row['Query']} ↔ {row['Reference']}: {row['Similarity (%)']:.2f}%"
                    lines.append(line)

                for i in range(0, len(lines), lines_per_page):
                    fig, ax = plt.subplots(figsize=(8.5, 11))
                    ax.axis('off')
                    ax.text(0, 1, "\n".join(lines[i:i + lines_per_page]), 
                           va='top', fontsize=10, family='monospace')
                    pdf.savefig(fig, bbox_inches='tight')
                    plt.close(fig)
                    
            self.results_text_gen.insert(tk.END, f"📄 Results saved as PDF: {pdf_path}\n")
        except Exception as e:
            self.results_text_gen.insert(tk.END, f"❌ Error saving PDF: {e}\n")

    def calculate_genetic_indices(self):
        self.update_status("Calculating genetic indices...")
        query_filename = self.entry_query.get().strip()
        if not query_filename:
            messagebox.showwarning("Warning", "Please select the query samples file.")
            self.update_status("Ready")
            return

        try:
            query_df = pd.read_excel(query_filename, engine='openpyxl', header=1)
        except Exception as e:
            messagebox.showerror("Error", f"Error reading file:\n{e}")
            self.update_status("Error reading file")
            return

        if query_df.empty:
            messagebox.showerror("Error", "The query file is empty.")
            self.update_status("Empty file")
            return

        query_name_col = query_df.columns[0]
        query_data = query_df.drop(columns=[query_name_col])

        sample_results = []
        pop_A, pop_Ae, pop_He, pop_Ho, pop_F = [], [], [], [], []
        
        for idx, row in query_df.iterrows():
            sample = row[query_name_col]
            if pd.isna(sample) or str(sample).strip() == '':
                continue
                
            alleles = []
            heterozygotes = 0
            total_loci = 0
            
            for i in range(0, len(row)-1, 2):
                if i+1 < len(row) and i+2 < len(row):
                    locus1 = row[i+1]
                    locus2 = row[i+2]
                    
                    if self.isanumber(locus1) and self.isanumber(locus2):
                        locus1 = float(locus1)
                        locus2 = float(locus2)
                        alleles.append((i//2, locus1))
                        alleles.append((i//2, locus2))
                        if locus1 != locus2:
                            heterozygotes += 1
                        total_loci += 1
            
            if not alleles:
                sample_results.append({
                    'Sample': sample,
                    'A': 0,
                    'Ho': 0.0
                })
                continue
                
            unique_alleles = set(alleles)
            A = len(unique_alleles)
            Ho = heterozygotes / total_loci if total_loci > 0 else 0
            
            sample_results.append({
                'Sample': sample,
                'A': A,
                'Ho': round(Ho, 3)
            })
            
            allele_counts = defaultdict(int)
            for allele in alleles:
                allele_counts[allele] += 1
                
            freqs = np.array([count / len(alleles) for count in allele_counts.values()])
            p_squared_sum = np.sum(freqs ** 2)
            
            Ae = 1 / p_squared_sum if p_squared_sum > 0 else 0
            He = 1 - p_squared_sum
            F = (He - Ho) / He if He > 0 else 0

            pop_A.append(A)
            pop_Ae.append(Ae)
            pop_He.append(He)
            pop_Ho.append(Ho)
            pop_F.append(F)

        if not sample_results:
            messagebox.showwarning("Warning", "No valid samples found.")
            self.update_status("No valid samples")
            return

        mean_A = round(np.mean(pop_A), 3) if pop_A else 0
        mean_Ae = round(np.mean(pop_Ae), 3) if pop_Ae else 0
        mean_He = round(np.mean(pop_He), 3) if pop_He else 0
        mean_Ho = round(np.mean(pop_Ho), 3) if pop_Ho else 0
        mean_F = round(np.mean(pop_F), 3) if pop_F else 0

        df_samples = pd.DataFrame(sample_results)
        df_population = pd.DataFrame({
            'A': [mean_A],
            'Ae': [mean_Ae],
            'Ho': [mean_Ho],
            'He': [mean_He],
            'F': [mean_F]
        })

        tsv_path = os.path.join(self.input_file_directory, "Genetic_indices.tsv")
        pdf_path = os.path.join(self.input_file_directory, "Genetic_indices.pdf")

        try:
            with open(tsv_path, 'w', encoding='utf-8') as f:
                f.write("GMDA - Genetic Indices\n\n")
                f.write("Samples\n")
                df_samples.to_csv(f, sep="\t", index=False)
                f.write("\nPopulation\n")
                df_population.to_csv(f, sep=" ", index=False)
        except Exception as e:
            messagebox.showerror("Error", f"Error saving TSV:\n{e}")
            self.update_status("Error saving results")
            return

        try:
            with PdfPages(pdf_path) as pdf:
                fig, ax = plt.subplots(figsize=(8.5, 11))
                ax.axis('tight')
                ax.axis('off')
                               
                table = ax.table(cellText=df_samples.values,
                               colLabels=df_samples.columns,
                               loc='center',
                               cellLoc='center')
                table.auto_set_font_size(False)
                table.set_fontsize(10)
                table.scale(1, 1.5)
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
                
                fig, ax = plt.subplots(figsize=(8.5, 3))
                ax.axis('tight')
                ax.axis('off')
                ax.set_title("Population", fontsize=12, pad=20)
                
                table = ax.table(cellText=df_population.values,
                               colLabels=df_population.columns,
                               loc='center',
                               cellLoc='center')
                table.auto_set_font_size(False)
                table.set_fontsize(10)
                table.scale(1, 1.5)
                
                for j in range(len(df_population.columns)):
                    table[(1, j)].set_facecolor('#E6E6FA')
                    table[(1, j)].set_text_props(weight='bold')
                
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
                
            self.results_text_gen.insert(tk.END, f"\n💾 Genetic indices saved as:\n{tsv_path}\n{pdf_path}\n")
        except Exception as e:
            self.results_text_gen.insert(tk.END, f"❌ Error saving PDF:\n{e}\n")
            return

        try:
            query_names = query_df.iloc[:, 0].tolist()
            query_data = query_df.iloc[:, 1:]
            
            valid_indices = [i for i, name in enumerate(query_names) if pd.notna(name) and str(name).strip() != '']
            query_names = [query_names[i] for i in valid_indices]
            query_data = query_data.iloc[valid_indices]
            
            if len(query_names) < 2:
                self.results_text_gen.insert(tk.END, "\n⚠️ At least 2 Query samples needed for dendrogram and PCoA\n")
                return
            
            matrix = pd.DataFrame(index=query_names, columns=query_names, dtype=float)
            np.fill_diagonal(matrix.values, 1.0)

            for i in range(len(query_data)):
                for j in range(i + 1, len(query_data)):
                    sim = self.calcular_similaridade(query_data.iloc[i], query_data.iloc[j])
                    matrix.iat[i, j] = sim
                    matrix.iat[j, i] = sim

            distance_matrix = 1 - matrix.astype(float).values
            np.fill_diagonal(distance_matrix, 0)
            
            condensed_dist = squareform(distance_matrix, checks=False)
            linkage_matrix = linkage(condensed_dist, method='average')
            
            self.genetic_linkage_matrix = linkage_matrix
            self.genetic_distance_matrix = condensed_dist
            self.genetic_labels = query_names
            
            if self.bootstrap_var_gen.get():
                self.results_text_gen.insert(tk.END, "\n🔄 Calculating bootstrap supports for Query samples...\n")
                self.root.update_idletasks()
                supports = self.calculate_bootstrap_support(query_data, linkage_matrix)
                self.plot_query_dendrogram(linkage_matrix, query_names, supports)
            else:
                self.plot_query_dendrogram(linkage_matrix, query_names)
            
            self.genetic_pcoa_coords = self.generate_query_pcoa(distance_matrix, query_names)
            
        except Exception as e:
            self.results_text_gen.insert(tk.END, f"\n❌ Error in Query-only analysis: {e}\n")

        self.update_status("Genetic indices analysis completed")

    def plot_query_dendrogram(self, linkage_matrix, labels, supports=None):
        plt.figure(figsize=(12, 8))
        dendro = dendrogram(linkage_matrix, labels=labels, leaf_rotation=90, leaf_font_size=14)

        if supports:
            icoord = np.array(dendro['icoord'])
            dcoord = np.array(dendro['dcoord'])
            
            for i, d in zip(icoord, dcoord):
                x = 0.5 * (i[1] + i[2])
                y = d[1]
                
                for cluster_idx, row in enumerate(linkage_matrix):
                    if abs(row[2] - y) < 1e-5:
                        support = supports.get(cluster_idx + len(labels), 0)
                        if support >= 50:
                            plt.text(x, y, f"{int(support)}%", 
                                   ha='center', va='bottom',
                                   fontsize=8, color='red')
                        break

        title = "GMDA - Query Samples UPGMA Dendrogram"
        if supports:
            title += " (with bootstrap)"
            
        plt.title(title, fontsize=14)
        plt.ylabel("Dissimilarity (1-S)", fontsize=14)
        plt.tight_layout()
        
        filename = "Query_Dendrogram"
        if supports:
            filename += "_with_bootstrap"
            
        plot_path = os.path.join(self.input_file_directory, f"{filename}.png")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        self.results_text_gen.insert(tk.END, f"📊 Query dendrogram saved as: {plot_path}\n")

    def generate_query_pcoa(self, distance_matrix, labels):
        try:
            mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
            coords = mds.fit_transform(distance_matrix)
            
            total_inertia = np.sum(mds.dissimilarity_matrix_ ** 2)
            inertia_per_component = []
            for i in range(2):
                inertia = np.sum((coords[:, i] - coords[:, i].mean()) ** 2)
                inertia_per_component.append(inertia)
            
            var_explained = [inertia / total_inertia * 100 for inertia in inertia_per_component]
            
            plt.figure(figsize=(10, 8))
            
            colors = plt.cm.get_cmap('tab10', len(labels))
            
            for i, (x, y) in enumerate(coords):
                plt.scatter(x, y, color=colors(i), s=100, alpha=0.7)
                plt.text(x, y, labels[i], fontsize=10, ha='center', va='bottom')
                        
            plt.xlabel(f"PCo1 ({var_explained[0]:.1f}%)", fontsize=14)
            plt.ylabel(f"PCo2 ({var_explained[1]:.1f}%)", fontsize=14)
            plt.title("GMDA - Query Samples PCoA", fontsize=14)
            plt.grid(True, alpha=0.3, linestyle='--')
            plt.axhline(0, color='gray', linewidth=0.5)
            plt.axvline(0, color='gray', linewidth=0.5)
            plt.tight_layout()
            
            pcoa_path = os.path.join(self.input_file_directory, "Query_PCoA.png")
            plt.savefig(pcoa_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.results_text_gen.insert(tk.END, f"📈 Query PCoA plot saved as: {pcoa_path}\n")
            self.results_text_gen.insert(tk.END, f"📊 Variance explained: PCo1 = {var_explained[0]:.1f}%, PCo2 = {var_explained[1]:.1f}%\n")
            
            return coords
            
        except Exception as e:
            self.results_text_gen.insert(tk.END, f"❌ Error generating PCoA: {e}\n")
            return None

    def calculate_morpho(self):
        self.update_status("Performing morphometric analysis...")
        filename = self.entry_file_morpho.get().strip()
        if not filename:
            messagebox.showwarning("Warning", "Please select a data file.")
            self.update_status("Ready")
            return

        try:
            df = pd.read_excel(filename, engine='openpyxl', header=0)
        except Exception as e:
            messagebox.showerror("Error", f"Error reading file:\n{e}")
            self.update_status("Error reading file")
            return

        if df.empty:
            messagebox.showerror("Error", "The file is empty.")
            self.update_status("Empty file")
            return

        if df.columns[0] != 'Amostra':
            messagebox.showerror("Error", "A primeira coluna deve ser 'Amostra'")
            self.update_status("Invalid file format")
            return

        try:
            n_replicates = int(self.replicates_var.get())
            if n_replicates < 1:
                raise ValueError("Number of replicates must be at least 1")
        except ValueError as e:
            messagebox.showerror("Error", f"Invalid number of replicates: {e}")
            self.update_status("Invalid replicates")
            return

        variable_names = df.columns[1:]
        sample_names = []
        measurements = defaultdict(list)

        for i in range(0, len(df), n_replicates):
            sample_name = df.iloc[i, 0]
            if pd.isna(sample_name):
                continue
                
            sample_names.append(str(sample_name))
            
            for var in variable_names:
                replicates = []
                for j in range(n_replicates):
                    if i+j < len(df):
                        value = df.iloc[i+j][var]
                        if self.isanumber(value):
                            replicates.append(float(value))
            
                measurements[var].append(replicates)

        self.results_text_morpho.delete(1.0, tk.END)
        self.results_text_morpho.insert(tk.END, "📊 MORPHOMETRIC ANALYSIS RESULTS\n")
        self.results_text_morpho.insert(tk.END, "="*50 + "\n\n")

        self.morpho_results = []
        self.anova_results = []
        self.tukey_results = []
        self.effect_size_results = []
        
        mean_values = {}
        for var in variable_names:
            mean_values[var] = {}
            for i, sample in enumerate(sample_names):
                if measurements[var][i]:
                    mean_values[var][sample] = round(np.mean(measurements[var][i]), 2)
                else:
                    mean_values[var][sample] = np.nan

        for var in variable_names:
            data_for_anova = []
            group_labels = []
            
            for i, sample_measures in enumerate(measurements[var]):
                if sample_measures:
                    data_for_anova.extend(sample_measures)
                    group_labels.extend([sample_names[i]] * len(sample_measures))
            
            if not data_for_anova or len(set(group_labels)) < 2:
                continue
                
            stats_result = self.calculate_statistics(data_for_anova)
            
            effect_sizes = {}
            ci_lower = {}
            ci_upper = {}
            
            if len(measurements[var]) >= 2:
                for i in range(len(measurements[var])):
                    for j in range(i+1, len(measurements[var])):
                        if measurements[var][i] and measurements[var][j]:
                            group1 = measurements[var][i]
                            group2 = measurements[var][j]
                            
                            n1, n2 = len(group1), len(group2)
                            mean1, mean2 = np.mean(group1), np.mean(group2)
                            pooled_std = np.sqrt(((n1-1)*np.std(group1)**2 + (n2-1)*np.std(group2)**2) / (n1+n2-2))
                            d = (mean1 - mean2) / pooled_std
                            
                            se = np.sqrt((n1 + n2)/(n1*n2) + d**2/(2*(n1+n2)))
                            ci_low = d - 1.96 * se
                            ci_high = d + 1.96 * se
                            
                            pair_name = f"{sample_names[i]}-{sample_names[j]}"
                            effect_sizes[pair_name] = round(d, 2)
                            ci_lower[pair_name] = round(ci_low, 2)
                            ci_upper[pair_name] = round(ci_high, 2)
            
            if effect_sizes:
                for pair, d in effect_sizes.items():
                    self.effect_size_results.append({
                        'Trait': var,
                        'Comparison': pair,
                        "Cohen's d": d,
                        'CI 95% Lower': ci_lower[pair],
                        'CI 95% Upper': ci_upper[pair]
                    })
            
            self.results_text_morpho.insert(tk.END, f"📏 Trait: {var}\n")
            self.results_text_morpho.insert(tk.END, f"  • Mean: {stats_result['mean']:.2f} ± {stats_result['mean_std']:.2f}\n")
            self.results_text_morpho.insert(tk.END, f"  • Median: {stats_result['median']:.2f} ± {stats_result['median_std']:.2f}\n")
            
            if isinstance(stats_result['mode'], str):
                self.results_text_morpho.insert(tk.END, f"  • Mode: {stats_result['mode']}\n")
            else:
                self.results_text_morpho.insert(tk.END, f"  • Mode: {stats_result['mode']:.2f} ± {stats_result['mode_std']:.2f}\n")
            
            if effect_sizes:
                self.results_text_morpho.insert(tk.END, "  • Effect sizes (Cohen's d) between groups:\n")
                for pair, d in effect_sizes.items():
                    ci_low = ci_lower[pair]
                    ci_high = ci_upper[pair]
                    self.results_text_morpho.insert(tk.END, f"    ↳ {pair}: d = {d:.2f} (95% CI: {ci_low:.2f} to {ci_high:.2f})\n")
            
            p_value = None
            try:
                if len(measurements[var]) > 1:
                    groups = []
                    for sample_measures in measurements[var]:
                        if sample_measures:
                            groups.append(sample_measures)
                    
                    if len(groups) >= 2:
                        f_val, p_val = stats.f_oneway(*groups)
                        anova_result = {'F': f_val, 'p': p_val}
                    else:
                        anova_result = None
                else:
                    anova_result = None
                    
                if anova_result:
                    p_value = anova_result['p']
                    self.results_text_morpho.insert(tk.END, f"  • ANOVA p-value: {p_value:.4f}\n")
                    
                    if p_value < 0.05:
                        self.results_text_morpho.insert(tk.END, "    ↳ ✅ Significant differences found (p < 0.05)\n")
                        
                        tukey = pairwise_tukeyhsd(
                            endog=data_for_anova,
                            groups=group_labels,
                            alpha=0.05
                        )
                        
                        tukey_df = pd.DataFrame(data=tukey._results_table.data[1:], 
                                               columns=tukey._results_table.data[0])
                        
                        groups = self._assign_tukey_groups(tukey_df, sample_names)
                        
                        tukey_formatted = {
                            'Trait': var,
                            'Results': []
                        }
                        
                        for sample in sample_names:
                            tukey_formatted['Results'].append({
                                'Sample': sample,
                                'Mean': mean_values[var][sample],
                                'Group': groups.get(sample, '?')
                            })
                        
                        self.tukey_results.append(tukey_formatted)
                        
                        self.results_text_morpho.insert(tk.END, "    ↳ Tukey HSD grouping:\n")
                        for result in tukey_formatted['Results']:
                            self.results_text_morpho.insert(tk.END, 
                                f"      {result['Sample']}: {result['Mean']:.2f}{result['Group']}\n")
                        
                        self.anova_results.append({
                            'Trait': var,
                            'ANOVA_p_value': p_value,
                            'Tukey_results': tukey_df.to_string()
                        })
                    else:
                        self.results_text_morpho.insert(tk.END, "    ↳ ❌ No significant differences found (p ≥ 0.05)\n")
                        self.anova_results.append({
                            'Trait': var,
                            'ANOVA_p_value': p_value,
                            'Tukey_results': "Not applicable (p ≥ 0.05)"
                        })
                else:
                    self.results_text_morpho.insert(tk.END, "  • ANOVA: Not enough groups for analysis\n")
                    self.anova_results.append({
                        'Trait': var,
                        'ANOVA_p_value': "N/A",
                        'Tukey_results': "Not enough groups for analysis"
                    })
            except Exception as e:
                self.results_text_morpho.insert(tk.END, f"  • ❌ Error in ANOVA: {str(e)}\n")
                self.anova_results.append({
                    'Trait': var,
                    'ANOVA_p_value': "Error",
                    'Tukey_results': f"Error: {str(e)}"
                })
            
            self.results_text_morpho.insert(tk.END, f"  • Count: {stats_result['count']}\n\n")
            
            self.morpho_results.append({
                'Trait': var,
                'Count': stats_result['count'],
                'Mean': f"{stats_result['mean']:.2f} ± {stats_result['mean_std']:.2f}",
                'Median': f"{stats_result['median']:.2f} ± {stats_result['median_std']:.2f}",
                'Mode': stats_result['mode'],
                'ANOVA_p_value': f"{p_value:.4f}" if p_value is not None else "N/A"
            })

        try:
            data_for_clustering = []
            for var in variable_names:
                var_means = []
                for sample_measures in measurements[var]:
                    if sample_measures:
                        var_means.append(np.mean(sample_measures))
                    else:
                        var_means.append(np.nan)
                data_for_clustering.append(var_means)
            
            df_clustering = pd.DataFrame(data_for_clustering, index=variable_names, columns=sample_names).T
            
            df_clustering = df_clustering.dropna()
            
            if len(df_clustering) > 1:
                morpho_result = self.generate_morpho_dendrogram(
                    df_clustering.values,
                    df_clustering.index.tolist(),
                    bootstrap=self.bootstrap_var_morpho.get()
                )
                
                if morpho_result:
                    self.morpho_linkage_matrix = morpho_result['linkage_matrix']
                    self.morpho_distance_matrix = morpho_result['distance_matrix']
                    self.morpho_labels = df_clustering.index.tolist()
                
                pca_result = self.generate_pca(df_clustering, df_clustering.index.tolist())
                if pca_result:
                    self.morpho_pca_coords = pca_result['pca_result']
            else:
                self.results_text_morpho.insert(tk.END, "\n⚠️ Not enough valid samples for dendrogram and PCA\n")
            
        except Exception as e:
            self.results_text_morpho.insert(tk.END, f"\n❌ Error in analysis: {e}\n")

        self.update_status("Morphometric analysis completed")

    def _assign_tukey_groups(self, tukey_df, sample_names):
        groups = defaultdict(set)
        
        for _, row in tukey_df.iterrows():
            group1 = row['group1']
            group2 = row['group2']
            reject = row['reject']
            
            if not reject:
                found = False
                for g in groups.values():
                    if group1 in g or group2 in g:
                        g.add(group1)
                        g.add(group2)
                        found = True
                        break
                if not found:
                    groups[len(groups)] = {group1, group2}
        
        all_samples_in_groups = set()
        for g in groups.values():
            all_samples_in_groups.update(g)
        
        for sample in sample_names:
            if sample not in all_samples_in_groups:
                groups[len(groups)] = {sample}
        
        letters = string.ascii_lowercase
        group_letters = {}
        for i, (_, samples) in enumerate(groups.items()):
            if i < len(letters):
                letter = letters[i]
            else:
                letter = f"z{i-len(letters)+1}"
            
            for sample in samples:
                group_letters[sample] = letter
        
        return group_letters

    def save_morpho_results(self):
        if not self.morpho_results:
            messagebox.showinfo("No data", "No results available to save.")
            return

        df_stats = pd.DataFrame(self.morpho_results)
        df_stats = df_stats[['Trait', 'Count', 'Mean', 'Median', 'Mode', 'ANOVA_p_value']]
        
        if self.tukey_results:
            sample_data = {}
            for result in self.tukey_results:
                trait = result['Trait']
                for res in result['Results']:
                    sample = res['Sample']
                    if sample not in sample_data:
                        sample_data[sample] = {'Sample': sample}
                    sample_data[sample][f"{trait}_Mean"] = f"{res['Mean']:.2f}{res['Group']}"
            
            df_tukey = pd.DataFrame.from_dict(sample_data, orient='index')
            df_tukey = df_tukey.sort_index()
            df_tukey = df_tukey.reset_index(drop=True)
            
            cols = ['Sample'] + [col for col in df_tukey.columns if col != 'Sample']
            df_tukey = df_tukey[cols]
        else:
            df_tukey = pd.DataFrame()

        if self.effect_size_results:
            df_effect_sizes = pd.DataFrame(self.effect_size_results)
            df_effect_sizes = df_effect_sizes[['Trait', 'Comparison', "Cohen's d", 'CI 95% Lower', 'CI 95% Upper']]
        else:
            df_effect_sizes = pd.DataFrame(columns=['Trait', 'Comparison', "Cohen's d", 'CI 95% Lower', 'CI 95% Upper'])

        save_path = filedialog.asksaveasfilename(
            initialdir=self.input_file_directory,
            defaultextension=".xlsx",
            filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")],
            title="Save morphometric results as..."
        )
        if not save_path:
            return

        try:
            with pd.ExcelWriter(save_path, engine='openpyxl') as writer:
                title_df = pd.DataFrame({"": ["GMDA - Morphometric Statistics"]})
                title_df.to_excel(writer, sheet_name='Statistics', index=False, header=False)
                df_stats.to_excel(writer, sheet_name='Statistics', startrow=2, index=False)
                
                if not df_tukey.empty:
                    title_df = pd.DataFrame({"": ["GMDA - Tukey HSD Grouping (Mean±Group)"]})
                    title_df.to_excel(writer, sheet_name='Tukey_HSD', index=False, header=False)
                    df_tukey.to_excel(writer, sheet_name='Tukey_HSD', startrow=2, index=False)
                
                title_df = pd.DataFrame({"": ["GMDA - Effect Sizes and Confidence Intervals"]})
                title_df.to_excel(writer, sheet_name='Effect_Sizes', index=False, header=False)
                df_effect_sizes.to_excel(writer, sheet_name='Effect_Sizes', startrow=2, index=False)
            
            self.results_text_morpho.insert(tk.END, f"\n💾 Results saved as Excel: {save_path}\n")
        except Exception as e:
            self.results_text_morpho.insert(tk.END, f"❌ Error saving Excel: {e}\n")

    def calculate_statistics(self, data):
        if not data:
            return {
                'count': 0,
                'mean': "N/A",
                'mean_std': "N/A",
                'median': "N/A",
                'median_std': "N/A",
                'mode': "N/A",
                'mode_std': "N/A"
            }
        
        data = np.array(data)
        stats_result = {
            'count': len(data),
            'mean': round(np.mean(data), 2),
            'mean_std': round(np.std(data), 2),
            'median': round(np.median(data), 2),
            'median_std': round(np.std(data - np.median(data)), 2)
        }
        
        try:
            modes = multimode(data)
            if len(modes) == 1:
                stats_result['mode'] = round(modes[0], 2)
                stats_result['mode_std'] = round(np.std(data - modes[0]), 2)
            else:
                stats_result['mode'] = f"multiple"
                stats_result['mode_std'] = "N/A"
        except Exception:
            stats_result['mode'] = "N/A"
            stats_result['mode_std'] = "N/A"
            
        return stats_result

    def generate_morpho_dendrogram(self, data, labels, bootstrap=False, replicates=100):
        try:
            scaler = StandardScaler()
            data_scaled = scaler.fit_transform(data)
            
            dist_matrix = pdist(data_scaled, metric='euclidean')
            linkage_matrix = linkage(dist_matrix, method='average')
            
            plt.figure(figsize=(12, 8))
            dendro = dendrogram(
                linkage_matrix,
                labels=labels,
                orientation='top',
                leaf_rotation=90,
                leaf_font_size=16,
                color_threshold=0.7*np.max(linkage_matrix[:,2]))
            
            if bootstrap:
                supports = self.calculate_bootstrap_support(
                    pd.DataFrame(data), 
                    linkage_matrix,
                    replicates=replicates
                )
                
                icoord = np.array(dendro['icoord'])
                dcoord = np.array(dendro['dcoord'])
                
                for i, d in zip(icoord, dcoord):
                    x = 0.5 * (i[1] + i[2])
                    y = d[1]
                    
                    for cluster_idx, row in enumerate(linkage_matrix):
                        if abs(row[2] - y) < 1e-5:
                            support = supports.get(cluster_idx + len(labels), 0)
                            if support >= 50:
                                plt.text(x, y, f"{int(support)}%", 
                                       ha='center', va='bottom',
                                       fontsize=8, color='red')
                            break
            
            title = "GMDA - Morphometric UPGMA Dendrogram"
            if bootstrap:
                title += f" (Bootstrap n={replicates})"
                
            plt.title(title, fontsize=14, pad=20)
            plt.xlabel("Samples", fontsize=12)
            plt.ylabel("Euclidean Distance", fontsize=12)
            plt.grid(True, alpha=0.3, linestyle='--')
            plt.tight_layout()
            
            filename = "Morphometric_Dendrogram"
            if bootstrap:
                filename += "_with_bootstrap"
                
            plot_path = os.path.join(self.input_file_directory, f"{filename}.png")
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.results_text_morpho.insert(tk.END, f"\n📊 Dendrogram saved as: {plot_path}\n")
            
            return {
                'linkage_matrix': linkage_matrix,
                'distance_matrix': dist_matrix,
                'bootstrap_supports': supports if bootstrap else None,
                'plot_path': plot_path
            }
            
        except Exception as e:
            plt.close()
            self.results_text_morpho.insert(tk.END, f"\n❌ Error generating dendrogram: {e}\n")
            return None

    def generate_pca(self, data, sample_names):
        try:
            scaler = StandardScaler()
            data_scaled = scaler.fit_transform(data)
            
            pca = PCA()
            pca_result = pca.fit_transform(data_scaled)
            explained_var = pca.explained_variance_ratio_
            
            plt.figure(figsize=(10, 8))
            
            plt.scatter(pca_result[:, 0], pca_result[:, 1], alpha=0.7, s=100, c='steelblue')
            
            for i, name in enumerate(sample_names):
                plt.annotate(str(name), (pca_result[i, 0], pca_result[i, 1]),
                             xytext=(5, 5), textcoords='offset points', fontsize=16)
            
            loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
            for i, var in enumerate(data.columns):
                plt.arrow(0, 0, loadings[i, 0]*3, loadings[i, 1]*3,
                          color='red', alpha=0.7, head_width=0.1)
                plt.text(loadings[i, 0]*3.2, loadings[i, 1]*3.2, var, 
                         color='red', ha='center', va='center', fontsize=10)
            
            pc1_var = explained_var[0] * 100
            pc2_var = explained_var[1] * 100
            
            plt.xlabel(f'PC1 ({pc1_var:.1f}%)', fontsize=16)
            plt.ylabel(f'PC2 ({pc2_var:.1f}%)', fontsize=16)
            plt.title('GMDA - Principal Component Analysis', fontsize=14)
            plt.grid(True, alpha=0.3)
            plt.axhline(0, color='gray', linewidth=0.5)
            plt.axvline(0, color='gray', linewidth=0.5)
            plt.tight_layout()
            
            pca_path = os.path.join(self.input_file_directory, "Morphometric_PCA.png")
            plt.savefig(pca_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.results_text_morpho.insert(tk.END, f"\n📈 PCA plot saved as: {pca_path}\n")
            
            self.results_text_morpho.insert(tk.END, "\n📊 PCA Variance Explained:\n")
            for i, var in enumerate(explained_var):
                self.results_text_morpho.insert(tk.END, f"  PC{i+1}: {var*100:.1f}%\n")
            
            return {
                'pca_result': pca_result,
                'explained_variance': explained_var,
                'loadings': loadings,
                'plot_path': pca_path
            }
            
        except Exception as e:
            plt.close()
            self.results_text_morpho.insert(tk.END, f"\n❌ Error generating PCA: {e}\n")
            return None

    def calculate_mantel_test_dendrograms(self):
        self.update_status("Running Mantel test on dendrograms...")
        if (self.genetic_distance_matrix is None or self.morpho_distance_matrix is None or
            self.genetic_labels is None or self.morpho_labels is None):
            messagebox.showwarning("Warning", "Please run both genetic and morphometric analyses first to generate distance matrices.")
            self.update_status("Ready")
            return
        
        self.progress_bar.start(10)
        
        genetic_samples = set(self.genetic_labels)
        morpho_samples = set(self.morpho_labels)
        common_samples = list(genetic_samples & morpho_samples)
        
        if len(common_samples) < 3:
            messagebox.showerror("Error", f"Not enough common samples between datasets. Found: {len(common_samples)}")
            self.progress_bar.stop()
            self.update_status("Not enough common samples")
            return
        
        genetic_indices = [i for i, label in enumerate(self.genetic_labels) if label in common_samples]
        morpho_indices = [i for i, label in enumerate(self.morpho_labels) if label in common_samples]
        
        genetic_square = squareform(self.genetic_distance_matrix)
        morpho_square = squareform(self.morpho_distance_matrix)
        
        genetic_common = genetic_square[np.ix_(genetic_indices, genetic_indices)]
        morpho_common = morpho_square[np.ix_(morpho_indices, morpho_indices)]
        
        genetic_condensed = squareform(genetic_common, checks=False)
        morpho_condensed = squareform(morpho_common, checks=False)
        
        corr, p_value, z_score = self.perform_mantel_test(genetic_condensed, morpho_condensed)
        
        self.save_mantel_results(genetic_condensed, morpho_condensed, corr, p_value, z_score, common_samples, "Dendrograms")
        
        self.results_text_corr.delete(1.0, tk.END)
        self.results_text_corr.insert(tk.END, "🔗 MANTEL TEST RESULTS (DENDROGRAMS)\n")
        self.results_text_corr.insert(tk.END, "="*50 + "\n\n")
        self.results_text_corr.insert(tk.END, f"📊 Mantel correlation coefficient: {corr:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📈 P-value: {p_value:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📋 Z-score: {z_score:.4f}\n")
        self.results_text_corr.insert(tk.END, f"🔢 Number of samples: {len(common_samples)}\n")
        
        if p_value < 0.05:
            self.results_text_corr.insert(tk.END, "✅ The correlation is statistically significant (p < 0.05)\n")
        else:
            self.results_text_corr.insert(tk.END, "❌ The correlation is not statistically significant (p ≥ 0.05)\n")
        
        self.progress_bar.stop()
        self.update_status("Mantel test completed")

    def calculate_mantel_test_pca(self):
        self.update_status("Running Mantel test on PCA/PCoA...")
        if (self.genetic_pcoa_coords is None or self.morpho_pca_coords is None or
            self.genetic_labels is None or self.morpho_labels is None):
            messagebox.showwarning("Warning", "Please run both genetic and morphometric analyses first to generate PCoA/PCA coordinates.")
            self.update_status("Ready")
            return
        
        self.progress_bar.start(10)
        
        genetic_samples = set(self.genetic_labels)
        morpho_samples = set(self.morpho_labels)
        common_samples = list(genetic_samples & morpho_samples)
        
        if len(common_samples) < 3:
            messagebox.showerror("Error", f"Not enough common samples between datasets. Found: {len(common_samples)}")
            self.progress_bar.stop()
            self.update_status("Not enough common samples")
            return
        
        genetic_indices = [i for i, label in enumerate(self.genetic_labels) if label in common_samples]
        morpho_indices = [i for i, label in enumerate(self.morpho_labels) if label in common_samples]
        
        genetic_coords_common = self.genetic_pcoa_coords[genetic_indices]
        morpho_coords_common = self.morpho_pca_coords[morpho_indices]
        
        genetic_dist = pdist(genetic_coords_common, metric='euclidean')
        morpho_dist = pdist(morpho_coords_common, metric='euclidean')
        
        corr, p_value, z_score = self.perform_mantel_test(genetic_dist, morpho_dist)
        
        self.save_mantel_pca_results(genetic_dist, morpho_dist, corr, p_value, z_score, common_samples)
        
        self.results_text_corr.delete(1.0, tk.END)
        self.results_text_corr.insert(tk.END, "🔗 MANTEL TEST RESULTS (PCA/PCoA)\n")
        self.results_text_corr.insert(tk.END, "="*50 + "\n\n")
        self.results_text_corr.insert(tk.END, f"📊 Mantel correlation coefficient: {corr:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📈 P-value: {p_value:.4f}\n")
        self.results_text_corr.insert(tk.END, f"📋 Z-score: {z_score:.4f}\n")
        self.results_text_corr.insert(tk.END, f"🔢 Number of samples: {len(common_samples)}\n")
        
        if p_value < 0.05:
            self.results_text_corr.insert(tk.END, "✅ The correlation is statistically significant (p < 0.05)\n")
        else:
            self.results_text_corr.insert(tk.END, "❌ The correlation is not statistically significant (p ≥ 0.05)\n")
        
        self.progress_bar.stop()
        self.update_status("Mantel test completed")
    
    def save_mantel_pca_results(self, dist1, dist2, corr, p_value, z_score, sample_names):
        save_path = os.path.join(self.input_file_directory, "Mantel Test PCA_PCoA.pdf")
        
        try:
            with PdfPages(save_path) as pdf:
                fig, ax = plt.subplots(figsize=(8.5, 11))
                ax.axis('off')
                
                text_content = [
                    "GMDA - Mantel Test Analysis (PCA/PCoA)",
                    "=" * 50,
                    "",
                    f"📊 Mantel correlation coefficient: {corr:.4f}",
                    f"📈 P-value: {p_value:.4f}",
                    f"📋 Z-score: {z_score:.4f}",
                    f"🔢 Number of samples: {len(sample_names)}",
                    "",
                    "Sample names:",
                    ", ".join(sample_names)
                ]
                
                ax.text(0.1, 0.9, "\n".join(text_content), transform=ax.transAxes, 
                       verticalalignment='top', fontsize=12)
                
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
                
                fig, ax = plt.subplots(figsize=(8, 8))
                ax.scatter(dist1, dist2, alpha=0.7)
                
                z = np.polyfit(dist1, dist2, 1)
                p = np.poly1d(z)
                ax.plot(dist1, p(dist1), "r--", alpha=0.8)
                
                ax.set_xlabel("PCoA distance", fontsize=12)
                ax.set_ylabel("PCA distance", fontsize=12)
                ax.set_title(f"Mantel Test PCA/PCoA\n(r = {corr:.3f}, p = {p_value:.3f})", fontsize=14)
                ax.grid(True, alpha=0.3)
                
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
                
            self.results_text_corr.insert(tk.END, f"\n💾 Results saved as PDF: {save_path}\n")
            
        except Exception as e:
            self.results_text_corr.insert(tk.END, f"❌ Error saving results: {str(e)}\n")
    
    def perform_mantel_test(self, dist1, dist2, permutations=1000):
        matrix1 = squareform(dist1)
        matrix2 = squareform(dist2)
        
        result = mantel(matrix1, matrix2, method='pearson', permutations=permutations)
        
        return result[0], result[1], result[2]

    def save_mantel_results(self, dist1, dist2, corr, p_value, z_score, sample_names, analysis_type):
        default_filename = f"Mantel Test {analysis_type}.pdf"
        save_path = os.path.join(self.input_file_directory, default_filename)
        
        try:
            with PdfPages(save_path) as pdf:
                fig, ax = plt.subplots(figsize=(8.5, 11))
                ax.axis('off')
                
                text_content = [
                    f"GMDA - Mantel Test Analysis ({analysis_type})",
                    "=" * 50,
                    "",
                    f"📊 Mantel correlation coefficient: {corr:.4f}",
                    f"📈 P-value: {p_value:.4f}",
                    f"📋 Z-score: {z_score:.4f}",
                    f"🔢 Number of samples: {len(sample_names)}",
                    "",
                    "Sample names:",
                    ", ".join(sample_names)
                ]
                
                ax.text(0.1, 0.9, "\n".join(text_content), transform=ax.transAxes, 
                       verticalalignment='top', fontsize=12)
                
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
                
                fig, ax = plt.subplots(figsize=(8, 8))
                ax.scatter(dist1, dist2, alpha=0.7)
                
                z = np.polyfit(dist1, dist2, 1)
                p = np.poly1d(z)
                ax.plot(dist1, p(dist1), "r--", alpha=0.8)
                
                ax.set_xlabel(f"{analysis_type.split('/')[0]} distance", fontsize=12)
                ax.set_ylabel(f"{analysis_type.split('/')[-1]} distance", fontsize=12)
                ax.set_title(f"Mantel Test ({analysis_type})\n(r = {corr:.3f}, p = {p_value:.3f})", fontsize=14)
                ax.grid(True, alpha=0.3)
                
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
                
            self.results_text_corr.insert(tk.END, f"\n💾 Results saved as PDF: {save_path}\n")
            
        except Exception as e:
            self.results_text_corr.insert(tk.END, f"❌ Error saving results: {str(e)}\n")
    
if __name__ == "__main__":
    root = tk.Tk()
    app = GMDA(root)
    root.mainloop()