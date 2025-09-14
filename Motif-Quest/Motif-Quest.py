import streamlit as st
from Bio import SeqIO
import re
import pandas as pd
from io import StringIO
import base64
import altair as alt

# Add these functions here
def get_base64_of_bin_file(bin_file):
    with open(bin_file, 'rb') as f:
        data = f.read()
    return base64.b64encode(data).decode()

def display_local_image(image_path, width=120):
    # Get base64 representation of the image
    img_base64 = get_base64_of_bin_file(image_path)
    
    # HTML to display the image
    img_html = f'''
    <div style="text-align: center;">
        <img src="data:image/png;base64,{img_base64}" width="{width}">
    </div>
    '''
    st.markdown(img_html, unsafe_allow_html=True)

# Configure custom favicon
def set_dna_favicon():
    # SVG code for a DNA icon
    dna_icon = """
    <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 512 512">
        <path fill="#0062cc" d="M448 96c0-35.3-28.7-64-64-64s-64 28.7-64 64c0 23.6 12.9 44.3 32 55.4v209.2c-19.1 11.1-32 31.8-32 55.4c0 35.3 28.7 64 64 64s64-28.7 64-64c0-23.6-12.9-44.3-32-55.4V151.4c19.1-11.1 32-31.8 32-55.4zM128 96c0-35.3-28.7-64-64-64S0 60.7 0 96c0 23.6 12.9 44.3 32 55.4v209.2c-19.1 11.1-32 31.8-32 55.4c0 35.3 28.7 64 64 64s64-28.7 64-64c0-23.6-12.9-44.3-32-55.4V151.4c19.1-11.1 32-31.8 32-55.4zM319.9 351.6c-19.6-10.8-41.9-16.9-65.9-16.9c-29.7 0-57.3 9.5-79.9 25.5c-6.3 4.5-12.3 9.5-17.8 15.1c-5.4-5.5-11.4-10.5-17.8-15.1c-22.6-16-50.2-25.5-79.9-25.5c-2.5 0-4.9 .1-7.3 .2C59.5 326 66.3 310.6 66.3 296c0-40.8-35.5-86.1-36.9-87.9c-4.6-5.7-13-6.6-18.7-1.9c-5.7 4.6-6.6 13-1.9 18.7c.3 .4 32.5 40.5 32.5 71.1c0 41.2-40.6 70.1-43.2 72.2c-5.3 4.1-6.2 11.8-2.1 17.1c4.1 5.3 11.8 6.2 17.1 2.1c.7-.5 1.6-1.3 2.7-2.1l0 0 0 0c41.2-31.2 124.2-7.9 124.2 56.7c0 .1-.1 .2-.1 .2c0 .1 0 .1 0 .2c0 6.5 5.2 11.8 11.7 12l.3 0c6.5 0 11.8-5.3 11.8-11.8c0-19.7 11.7-38 30.1-55.5c14.7-14.1 34.4-26.2 56.2-33.5c-2.2 4.9-3.4 10.3-3.4 16c0 23.4 19 42.4 42.4 42.4s42.4-19 42.4-42.4s-19-42.4-42.4-42.4c-12.1 0-22.9 5.1-30.7 13.2zm42.4 121.1L320 480l-42.3-7.4 42.3-7.4 42.3-7.4-42.3-7.4-42.3-7.4 42.3-7.4 42.3-7.4-42.3-7.4-42.3-7.4 42.3-7.4 42.3-7.4-42.3-7.4-42.3-7.4 42.3-7.4 42.3-7.4-42.3-7.4zm-212.8 0L0 480l42.3-7.4L0 465.1l0-7.4 42.3-7.4L0 442.8l0-7.4 42.3-7.4L0 420.6l0-7.4 42.3-7.4L0 398.3l0-7.4 42.3-7.4L0 376.1l0-7.4 42.3-7.4L0 353.8l0-7.4 42.3-7.4 42.3-7.4 42.3 7.4 42.3 7.4-42.3 7.4-42.3 7.4z"/>
    </svg>
    """
    # Convert SVG to base64
    b64_icon = base64.b64encode(dna_icon.encode()).decode()
    
    # Create the favicon
    favicon_html = f'''
    <link rel="icon" href="data:image/svg+xml;base64,{b64_icon}" type="image/svg+xml">
    '''
    
    # Apply the favicon
    st.markdown(favicon_html, unsafe_allow_html=True)

def generate_pattern(motif_length, min_repeats):
    """
    Generates a regex pattern to find repetitions of any nucleotide sequence.
    Ex: for motif_length=1, it searches for (A|T|G|C){min_repeats,}
    for motif_length=2, it searches for (AT|GC|...){min_repeats,}
    """
    # For a 1 nucleotide motif, we can look for any repeating nucleotide
    if motif_length == 1:
        # Search for A{X,}, T{X,}, G{X,}, C{X,}
        return r'[ATGCatgc]{1}' + '{%d,}' % min_repeats
    else:
        # For motifs of 2 or more nucleotides, we use a capture group and back-reference
        # This means: capture N nucleotides (the motif) and see if this same motif repeats M times
        # r"([ATGCatgc]{" + str(motif_length) + r"})" + "{" + str(min_repeats -1) + r",}"
        return r'([ATGCatgc]{' + str(motif_length) + r'})' + r'{\1,' + str(min_repeats) + r'}'

def find_microsatellites(sequence, min_repeats, max_distance):
    """
    Finds microsatellites (SSRs) in the DNA sequence.
    """
    results = []
    
    for repeat_type, repeat_info in min_repeats.items():
        motif_length = repeat_info['motif_length']
        min_repeats_value = repeat_info['min_repeats']
        
            # Regex pattern to find repetitions of motif
            # Ex: for motif of 1, ([ATGCatgc]) followed by {min_repeats_value,}
            # For motif > 1, ([ATGCatgc]{motif_length}) followed by \1{min_repeats_value-1,}
            # This ensures that the same motif is repeated.
        if motif_length == 1:
            # For motifs of 1, we look for repetitions of A, T, C or G
            # Ex: AAAAAAAAAAAAA or CCCCCCCCCCCC
            pattern = r'([ATGCatgc])\1{' + str(min_repeats_value - 1) + r',}'
        else:
            # For motifs of 2 or more nucleotides, such as (AT)(AT)(AT)
            pattern = r'([ATGCatgc]{' + str(motif_length) + r'})' + r'\1{' + str(min_repeats_value - 1) + r',}'
        
        try:
            matches = re.finditer(pattern, sequence, re.IGNORECASE)
            for match in matches:
                start, end = match.span()
                repeat_unit = match.group(1) # The repeating motif (e.g. 'AG', 'TGT')
                full_repeat_sequence = match.group(0) # The full SSR sequence (e.g. 'AGAGAG')
                
                results.append({
                    "type": repeat_type,
                    "start": start + 1,
                    "end": end,
                    "repeat_unit": repeat_unit,
                    "length": len(full_repeat_sequence),
                    "sequence": full_repeat_sequence,
                    "num_repeats": len(full_repeat_sequence) // len(repeat_unit) if len(repeat_unit) > 0 else 0,
                    "category": "Simple"
                })
        except Exception as e:
            st.warning(f"Error processing pattern for {repeat_type}: {e}")

    results.sort(key=lambda x: x['start'])
    composite_results = find_composite_ssrs(results, max_distance)
    
    results.extend(composite_results)

    return results

def find_composite_ssrs(ssrs, max_distance):
    """
    Identifies composite SSRs based on maximum distance.
    """
    composite_ssrs = []
    # Creates a copy so as not to modify the original list during iteration
    # ssrs_copy = sorted(ssrs, key=lambda x: x['start']) # Ensures it is sorted
    
    # It is already sorted before calling this function
    for i in range(len(ssrs) - 1):
        current_ssr = ssrs[i]
        next_ssr = ssrs[i + 1]
        
        # Distance is calculated between the end of the current SSR and the start of the next SSR
        distance = next_ssr["start"] - current_ssr["end"] -1 # Distance between the last nucleotide of one and the first of the other
        
        # Consider composite SSRs if the distance between them is <= max_distance and is positive (non-overlapping)
        if distance >= 0 and distance <= max_distance:
            composite_ssrs.append({
                "type": "composite",
                "start": current_ssr["start"],
                "end": next_ssr["end"],
                "components": [current_ssr, next_ssr],
                "distance": distance, # Distance between SSRs
                "repeat_unit": f"{current_ssr['repeat_unit']}-{next_ssr['repeat_unit']}", # Combination of reasons
                "length": next_ssr['end'] - current_ssr['start'] + 1 # Total length of compound
            })
    return composite_ssrs


def export_results(results, file_format="custom", project_name=None):
    """
    Exports the results in TXT format.
    """
    if file_format == "custom":
        output_lines = []
        
        # Add header if there is a project name
        if project_name:
            output_lines.append(f"Project: {project_name}\n")

        for seq in results:
            output_lines.append(f"ID: {seq['id']}\n")
            
            # Sort simple and compound results for more logical output
            sorted_individual_results = sorted([r for r in seq["results"] if r["type"] != "composite"], key=lambda x: x['start'])
            sorted_composite_results = sorted([r for r in seq["results"] if r["type"] == "composite"], key=lambda x: x['start'])

            if sorted_individual_results:
                output_lines.append("--- Simple SSRs ---\n")
                for result in sorted_individual_results:
                    output_lines.append(f"  Type: {result['type'].capitalize()}\n")
                    output_lines.append(f"  Start: {result['start']}\n")
                    output_lines.append(f"  End: {result['end']}\n")
                    output_lines.append(f"  Repeat Unit: {result['repeat_unit']}\n")
                    output_lines.append(f"  Length: {result['length']}\n")
                    output_lines.append(f"  Category: Simple\n")
                    output_lines.append("\n") # Blank line to separate entries

            if sorted_composite_results:
                output_lines.append("--- Composite SSRs ---\n")
                for result in sorted_composite_results:
                    output_lines.append(f"  Type: Composite\n")
                    output_lines.append(f"  Start: {result['start']}\n")
                    output_lines.append(f"  End: {result['end']}\n")
                    output_lines.append(f"  Components:\n")
                    for comp in result["components"]:
                        output_lines.append(f"    - Type: {comp['type'].capitalize()}, Unit: {comp['repeat_unit']}, Pos: {comp['start']}-{comp['end']}\n")
                    output_lines.append(f"  Distance Between Components: {result['distance']}\n")
                    output_lines.append(f"  Total Length: {result['length']}\n")
                    output_lines.append("\n") 
            
            output_lines.append("\n" + "="*50 + "\n\n") # Separator between sequences

        return "".join(output_lines)
   
    elif file_format == "txt":
        # Convert the result list in a single DataFrame to TXT (string format)
        all_results_data = []
        for seq_result in results:
            seq_id = seq_result['id']
            for item in seq_result['results']:
                row = {
                    "Project": project_name if project_name else "",
                    "ID": seq_id,
                    "Type": item['type'].capitalize(),
                    "Start": item['start'],
                    "End": item['end'],
                    "Repeat Unit": item['repeat_unit'],
                    "Length": item['length'],
                    "Category": item['type'].capitalize() if item['type'] != 'composite' else 'Composite',
                    "Details": ""
                }
                if item['type'] == 'composite':
                    components_str = "; ".join([f"{comp['type'].capitalize()}:{comp['repeat_unit']}[{comp['start']}-{comp['end']}]" for comp in item['components']])
                    row["Details"] = f"Distance: {item['distance']}, Components: {components_str}"
                else:
                    row["Details"] = f"Repeats: {item['num_repeats']}"
                all_results_data.append(row)
        
        if not all_results_data:
            return "No SSRs found with the current parameters."
            
        df_export = pd.DataFrame(all_results_data)
        return df_export.to_string(index=False) # Returns formatted string to TXT

    return ""


# Main Streamlit layout
st.set_page_config(
    layout="wide",
    page_title="Motif-Quest",
    page_icon="🧬"  # DNA emoji as default icon (can be replaced with custom icon)
)

# Apply custom favicon
set_dna_favicon()

# Add custom CSS to style navigation and sidebars
st.markdown("""
<style>
    /* Main color scheme */
    :root {
        --main-green: #006400;
        --main-blue: #1e88e5;
        --main-red: #006400;
        --light-green: #CCFFCC;
        --light-blue: #e3f2fd;              
        --light-red: #ffebee;
        --dark-green: #0d5d24;
        --dark-blue: #0d47a1;
        --dark-red: #b71c1c;
    }
    
     /* Scrollbar Customization */
    ::-webkit-scrollbar {
        width: 8px;
        height: 8px;
    }

    ::-webkit-scrollbar-track {
        background: var(--light-green);
        border-radius: 4px;
    }

    ::-webkit-scrollbar-thumb {
        background: var(--main-green);
        border-radius: 4px;
    }

    ::-webkit-scrollbar-thumb:hover {
        background: var(--dark-green);
    }
                
    /* Style for headers */
    h1, h2, h3 {
        color: var(--main-green) !important;
    }
    
    h2 {
        text-align: center !important;
        color: var(--main-green) !important; /* Keep color if already set */
    }
    
    /* Style for navigation */
    .nav-option {
        background-color: var(--light-blue);
        border: 1px solid #ddd;
        border-radius: 5px;
        padding: 8px 15px;
        margin-bottom: 8px;
        text-align: center;
        cursor: pointer;
        transition: background-color 0.3s;
    }
    .nav-option:hover {
        background-color: var(--main-blue);
        color: white;
    }
    .nav-option.active {
        background-color: var(--dark-blue);
        color: white;
        font-weight: bold;
    }
    .stRadio > div {
        display: none;
    }
    .navigation-title {
        font-weight: bold;
        font-size: 1.1em;
        margin-bottom: 10px;
        color: var(--main-green);
        text-align: center;
    }
    
    /* CSS to raise the settings panel */
    .parameter-config {
        margin-top: -40px;
        background-color: var(--light-green);
        padding: 15px;
        border-radius: 8px;
        border-left: 4px solid var(--main-green);
    }
    
    /* Center the title */
    .main-title {
        text-align: center;
        margin-bottom: 20px;
        color: var(--main-green);
    }
    
    /* Reduce the size of the configuration title */
    .config-header {
        font-size: 1.4em;
        margin-bottom: 10px;
        font-weight: bold;
        color: var(--main-green);
    }
    
    /* Style for file format selector */
    .file-format-select {
        margin-top: 10px;
        margin-bottom: 15px;
        background-color: var(--light-blue);
        padding: 10px;
        border-radius: 5px;
    }
    
    /* Style for sidebars */
    [data-testid="stSidebar"] strong {
        color: var(--main-green) !important;
    }
    
    /* Style for the right sidebar */
    .right-sidebar {
        position: fixed;
        top: 0;
        right: 0;
        width: 17%;
        height: 100vh;
        background-color: var(--light-green);
        padding: 1rem;
        border-left: 1px solid var(--main-green);
        overflow-y: auto;
        z-index: 999;
    }
    
    /* Style for main content */
    .main-content {
        width: 83%;
        padding-right: 17%;
    }
    
    /* Style for format selection buttons */
    .format-option {
        display: inline-block;
        background-color: var(--light-blue);
        border: 1px solid #ddd;
        border-radius: 5px;
        padding: 8px 15px;
        margin-right: 10px;
        margin-bottom: 10px;
        text-align: center;
        cursor: pointer;
        transition: background-color 0.3s;
    }
    .format-option:hover {
        background-color: var(--main-blue);
        color: white;
    }
    .format-option.selected {
        background-color: var(--main-blue);
        color: white;
        font-weight: bold;
    }
    
    /* Style for format buttons together */
    .format-buttons {
        display: flex;
        gap: 10px;
        margin-bottom: 15px;
    }
    
    /* Styling for buttons */
    .stButton > button {
        background-color: var(--main-blue);
        color: white;
        border: none;
        border-radius: 5px;
        padding: 8px 16px;
        font-weight: bold;
        transition: background-color 0.3s;
    }
    .stButton > button:hover {
        background-color: var(--dark-blue);
    }
    
    /* Analysis button */
    button[data-testid="baseButton-secondary"] {
        background-color: var(--main-green);
        color: white;
    }
    button[data-testid="baseButton-secondary"]:hover {
        background-color: var(--dark-green);
    }
    
    /* Download button */
    button[data-testid="baseButton-primary"] {
        background-color: var(--main-red);
        color: white;
    }
    button[data-testid="baseButton-primary"]:hover {
        background-color: var(--dark-red);
    }
    
    /* Style for expanding sections */
    .streamlit-expanderHeader {
        background-color: var(--light-blue);
        border-radius: 5px;
        color: var(--dark-blue);
    }
    
    /* Style for results tables */
    table {
        width: 100%;
        border-collapse: collapse;
    }
    th {
        background-color: var(--main-green);
        color: white;
        padding: 8px;
    }
    td {
        padding: 8px;
        border-bottom: 1px solid #ddd;
    }
    tr:nth-child(even) {
        background-color: var(--light-green);
    }
    tr:hover {
        background-color: #ddd;
    }
    
    /* Style for inputs */
    input, textarea, select {
    border: 1px solid var(--main-blue) !important;
    border-radius: 5px !important;
    background-color: var(--light-green) !important;
    color: #000000 !important; /* Add this line to set the text color */
    }
    input:focus, textarea:focus, select:focus {
    border: 2px solid var(--main-green) !important;
    box-shadow: 0 0 5px rgba(23, 156, 59, 0.5) !important;
    background-color: var(--light-green) !important;
    color: #000000 !important; /* Add this line too to keep the color when focusing */
    }
    
    /* Style for notices and messages */
    .stAlert {
        border-left: 5px solid var(--main-red) !important;
    }  

    /* ESSENTIAL: Ensures column width is fixed */
    .stDataFrame table {
        table-layout: fixed !important;
        width: 100% !important;        
    }
            
    /* "Index" column in results table */
    .stDataFrame td:nth-child(1), .stDataFrame th:nth-child(1) {
        width: 6px !important; /* Adjust this value as needed. (ex: 50px, 70px) */
        min-width: 6px !important;
        max-width: 10px !important;
        white-space: nowrap; /* For numbers, you generally don't want to break */
        overflow: hidden;
        text-overflow: ellipsis;
    }    
    
    /* "Type" Column (second column) */
    .stDataFrame td:nth-child(2), .stDataFrame th:nth-child(2) {
        width: 30px !important; /* Increased to be visible */
        min-width: 20px !important;
        max-width: 40px !important;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
    }

    /* "Start" Column in results table */
    .stDataFrame td:nth-child(3), .stDataFrame th:nth-child(3) {
        width: 10px !important; /* Adjust this value as needed (e.g., 50px, 70px) */
        min-width: 10px !important;
        max-width: 10px !important;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
    }

    /* "End" Column in results table */
    .stDataFrame td:nth-child(4), .stDataFrame th:nth-child(4) {
        width: 10px !important; /* Adjust this value as needed (e.g., 50px, 70px) */
        min-width: 10px !important;
        max-width: 10px !important;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
    }

    /* "Repeat Unit" Column in results table */
    .stDataFrame td:nth-child(5), .stDataFrame th:nth-child(5) {
        width: 80px !important; /* Adjust this value as needed (e.g., 80px, 120px) */
        min-width: 80px !important; /* Minimum width to ensure it doesn't get too small */
        max-width: 120px !important; /* Maximum width to prevent it from growing too large */
        white-space: nowrap; /* Prevents line breaks if you prefer not to wrap short units */
        overflow: hidden; /* Hides extra content if it doesn't fit */
        text-overflow: ellipsis;
    }

    /* "Length" Column in results table */
    .stDataFrame td:nth-child(6), .stDataFrame th:nth-child(6) {
        width: 15px !important; /* Adjust this value as needed (e.g., 50px, 70px) */
        min-width: 15px !important;
        max-width: 40px !important;
        white-space: nowrap; /* For numbers, wrapping is generally not desired */
        overflow: hidden;
        text-overflow: ellipsis;
    }

    /* "Category" Column in results table */
    .stDataFrame td:nth-child(7), .stDataFrame th:nth-child(7) {
        width: 200px !important; /* Adjust this value as needed (e.g., 50px, 70px) */
        min-width: 200px !important;
        max-width: 200px !important;
        white-space: nowrap; /* For numbers, wrapping is generally not desired */
        overflow: hidden;
        text-overflow: ellipsis;
    }

    /* "Details" Column in results table */
    .stDataFrame td:nth-child(8), .stDataFrame th:nth-child(8) {
        width: 400px !important;
        min-width: 400px !important; /* Minimum width for details */
        max-width: 400px !important; /* Maximum width for details */
        white-space: normal !important; /* ESSENTIAL: Allows line breaks for long details */
        word-wrap: break-word;
        overflow-wrap: break-word;
        overflow: visible; /* Ensures content is visible even if it exceeds the initial max-width */
        text-overflow: clip; /* Remove ellipsis if you want to see all the text with breaks */
    }

    /* Style for the "ID:" prefix */
    .id-prefix {
        color: var(--main-blue) !important; /* Or the color you want, for example, #FF0000 for red */
    }

    /* Style for the frame around the sequence ID */
    .id-frame {
        border: 3px solid var(--main-green); /* Green border */
        border-radius: 8px; /* Rounded corners */
        padding: 1px; /* Internal spacing between the border and the text */
        margin-bottom: 6px; /* Spacing below each ID frame */
        background-color: var(--light-green); /* Soft background for the frame */
        text-align: center; /* Centers the text inside the frame */
    }
</style>
""", unsafe_allow_html=True)

    # With the left sidebar
with st.sidebar:
    # Display local image (replace with the actual path to your image)
    # Make sure to have an image file in the same directory as your script
    display_local_image('lfdgv.png', width=180)
    
    st.markdown("")

    # Define navigation options
    options = ["Analysis", "Usage & Parameters", "Output", "About"]
    
    # Check if there is already an active page in the session
    if 'current_page' not in st.session_state:
        st.session_state.current_page = "Analysis"
    
    # Create a custom button for each option
    for option in options:
        if st.sidebar.button(option, key=f"nav_{option}", 
                            use_container_width=True,
                            type="secondary" if st.session_state.current_page != option else "primary"):
            st.session_state.current_page = option
            st.rerun()
    
    page = st.session_state.current_page

    # Parameter configuration in the sidebar (only when on the Analysis page)
    if page == "Analysis":
        st.markdown("---")

        display_local_image('dna.png', width=80)

        
        st.markdown('<h3 style="text-align: center;">SSR Search Configuration</h3>', unsafe_allow_html=True)
        
        
        # Initialize session state for dynamic fields
        if 'dynamic_repeats' not in st.session_state:
            st.session_state.dynamic_repeats = []
        
        # Default repeats dictionary
        default_repeats = {
            'mono': {'motif_length': 1, 'min_repeats': 12},
            'di': {'motif_length': 2, 'min_repeats': 6},
            'tri': {'motif_length': 3, 'min_repeats': 5},
            'tetra': {'motif_length': 4, 'min_repeats': 4},
            'penta': {'motif_length': 5, 'min_repeats': 4},
            'hexa': {'motif_length': 6, 'min_repeats': 4}
        }
        
        # Default input fields
        min_repeats = {}
        for repeat_type, repeat_info in default_repeats.items():
            with st.expander(f"**SSR motif length ({repeat_info['motif_length']})**"):
                min_repeats[repeat_type] = {
                    'motif_length': repeat_info['motif_length'],
                    'min_repeats': st.number_input(
                        f"Minimum Repeats", 
                        min_value=0,
                        value=repeat_info['min_repeats'],
                        key=f"default_{repeat_type}"
                    )
                }

        st.markdown("")

        # Dynamic fields after hexanucleotides
        st.markdown('<h3 style="text-align: center;">+ Add New Fields</h3>', unsafe_allow_html=True)

        
        # Function to add a new repeat field
        def add_dynamic_repeat():

            # Determine the next motif size
            if st.session_state.dynamic_repeats:
                next_motif_length = max(r['motif_length'] for r in st.session_state.dynamic_repeats) + 1
            else:
                next_motif_length = 7
            
            # Check if a field with this motif size already exists
            existing_motif = [r for r in st.session_state.dynamic_repeats if r['motif_length'] == next_motif_length]
            
            if not existing_motif:
                st.session_state.dynamic_repeats.append({
                    'motif_length': next_motif_length,
                    'min_repeats': 4  # Default value
                })

        # Button to add a new field
        st.button("+ Add New Nucleotide Field", on_click=add_dynamic_repeat)

        # Show dynamic fields
        for idx, repeat_info in enumerate(st.session_state.dynamic_repeats):
            with st.expander(f"Repeats of {repeat_info['motif_length']} nucleotide(s)"):
                
                # Use a unique key for each input field
                input_key = f"dynamic_repeats_{repeat_info['motif_length']}"
                
                # Add default value using the current session state value
                new_min_repeats = st.number_input(
                    f"Minimum Repeats (Motif of {repeat_info['motif_length']} nucleotide(s))", 
                    min_value=0,
                    value=repeat_info['min_repeats'],
                    key=input_key
                )
                
                # Update the value in the session state
                st.session_state.dynamic_repeats[idx]['min_repeats'] = new_min_repeats
                
                # Function to remove the field
                def remove_dynamic_repeat(index):
                    st.session_state.dynamic_repeats.pop(index)
                
                remove_key = f"remove_{repeat_info['motif_length']}"
                if st.button(f"Remove Field of {repeat_info['motif_length']} nucleotide(s)", 
                            key=remove_key, 
                            on_click=remove_dynamic_repeat, 
                            args=(idx,)):
                    st.rerun()

        # Add dynamic fields to the repeats dictionary
        for repeat_info in st.session_state.dynamic_repeats:
            min_repeats[f"custom_{repeat_info['motif_length']}"] = repeat_info

        st.markdown("---")
        
        # Maximum distance between composite SSRs
        max_distance = st.number_input(
            "**Maximum Distance Between Composite SSRs**", 
            min_value=0, 
            value=100
        )
        
        # Save settings in session state for access during analysis
        st.session_state.min_repeats = min_repeats
        st.session_state.max_distance = max_distance

        # Add version information to left sidebar
    st.markdown("---")

    display_local_image('ufsc.png')

    st.markdown("")
    st.markdown("")
    st.markdown("")
    st.markdown("")
    st.markdown("")
    st.markdown("")
    st.markdown("")
    st.markdown("")
    st.markdown("")
    st.markdown("")
    st.markdown("**Motif-Quest**")
    st.markdown("*Versão: v1.0*")
    st.markdown("*Data: 2025-03-28*")

display_local_image('abc.png', width=1150)

st.markdown("---")

# Content based on selected page
if page == "Analysis":
    # Main area (content)

    st.markdown("")

    # Project Name above Upload File
    st.markdown('<p style="font-size: 1em; font-weight: bold;">Project Name (optional)</p>', unsafe_allow_html=True)
    project_name = st.text_input("project_name_input", placeholder="Enter your project name", label_visibility="collapsed")  
    
    st.markdown("")
    
    # FASTA File Input
    st.markdown('<p style="font-size: 1em; font-weight: bold;">Upload a FASTA file (optional)</p>', unsafe_allow_html=True)
    uploaded_file = st.file_uploader("fasta_file_uploader", type=["fasta", "fa"], label_visibility="collapsed")

    st.markdown("")

    # Campo de texto para colar a sequência FASTA
    st.markdown('<p style="font-size: 1em; font-weight: bold;">Or paste FASTA sequence(s) here (optional)</p>', unsafe_allow_html=True)
    fasta_text = st.text_area("fasta_text_area", height=200, help="Paste your DNA sequences in FASTA format.", label_visibility="collapsed")
   
    # Initialize session state to file format if none exists
    if 'output_format' not in st.session_state:
        st.session_state.output_format = "Custom" 

    st.markdown("")
    
    # Analysis button
    run_button = st.button("**Run Analysis**")
    
    # Process data if button is clicked
    if run_button:
        sequences = []

        # Process FASTA file if provided
        if uploaded_file is not None:
            try:
                # Decode the file contents to a string
                fasta_content = uploaded_file.getvalue().decode('utf-8')
                fasta_io = StringIO(fasta_content)
                sequences.extend(list(SeqIO.parse(fasta_io, "fasta")))
            except Exception as e:
                st.error(f"Error processing FASTA file: {e}")

        # Process pasted text if provided
        if fasta_text.strip():
            try:
                fasta_io = StringIO(fasta_text)
                sequences.extend(list(SeqIO.parse(fasta_io, "fasta")))
            except Exception as e:
                st.error(f"Error processing FASTA text: {e}")

        # Check if there are sequences to process
        if not sequences:
            st.warning("No FASTA sequence provided. Please upload a file or paste a sequence.")
            # Optional: Clear previous results if no new sequence is provided
            st.session_state.analysis_results = None
        else:
            # Process each sequence
            processed_sequences = []
            for record in sequences:
                sequence = str(record.seq)
                results = find_microsatellites(sequence, 
                          st.session_state.min_repeats, 
                          st.session_state.max_distance)
                processed_sequences.append({
                    "id": record.id,
                    "results": results,
                    "sequence_length": len(sequence)
                })
            # STORE RESULTS IN SESSION STATE
            st.session_state.analysis_results = processed_sequences
    
    #  SHOW RESULTS AND DOWNLOAD BUTTONS IF THERE ARE RESULTS IN SESSION STATE
    if 'analysis_results' in st.session_state and st.session_state.analysis_results:
        st.header("Results")
        
        # Iterates over the results stored in the session state
        for seq in st.session_state.analysis_results: 
            st.markdown(
                f'<div class="id-frame">'
                f'<h2 style="margin: 0; padding: 0;"><span class="id-prefix">ID:</span> {seq["id"]}</h2>'
                f'</div>',
                unsafe_allow_html=True
            )
            
            # Sort the results for display
            sorted_display_results = sorted(seq["results"], key=lambda x: x['start'])

            if not sorted_display_results:
                st.info(f"No SSRs found for sequence {seq['id']} with the current parameters.")
                continue

            # Display results in a DataFrame for better visualization in the interface
            display_data = []
            for result in sorted_display_results:
                if result["type"] == "composite":
                    components_str = "; ".join([f"{comp['type'].capitalize()}:{comp['repeat_unit']}[{comp['start']}-{comp['end']}]" for comp in result['components']])
                    display_data.append({
                        "Type": "Composite",
                        "Start": result['start'],
                        "End": result['end'],
                        "Repeat Unit": result['repeat_unit'],
                        "Length": result['length'],
                        "Category": "Composite",
                        "Details": f"Distance: {result['distance']}, Components: {components_str}"
                    })
                else:
                    display_data.append({
                        "Type": result['type'].capitalize(),
                        "Start": result['start'],
                        "End": result['end'],
                        "Repeat Unit": result['repeat_unit'],
                        "Length": result['length'],
                        "Category": "Simple",
                        "Details": f"Repeats: {result['num_repeats']}"
                    })
            
            df_display = pd.DataFrame(display_data)
            st.dataframe(df_display, width=1200, height=300) # Display as interactive table         
        
        # START OF GRAPH GENERATION SECTION
        st.markdown("---")
        st.markdown("""
            <div class="id-frame">
                <h3 style="margin: 0; padding: 0;">SSR Position Visualization</h3>
            </div>
            """, unsafe_allow_html=True)
        st.markdown("")

        # Iterate over each sequence to create an individual position graph
        for seq_data in st.session_state.analysis_results:
            seq_id = seq_data['id']
            sequence_length = seq_data.get('sequence_length', 0)
            
            # 1. Only get results that are NOT of type 'composite'
            ssrs_simples_para_plotar = [
                result for result in seq_data['results'] if result['type'] != 'composite'
            ]

            # 2. If there are no simple SSRs, inform the user and skip to the next sequence
            if not ssrs_simples_para_plotar:
                st.info(f"No simple SSRs found to display in the sequence {seq_id}.")
                st.markdown("---")
                continue

            # 3. 
            dados_grafico = []
            for ssr in ssrs_simples_para_plotar:
                dados_grafico.append({
                    'Start': ssr['start'],
                    'End': ssr['end'],
                    'Type': ssr['type'].capitalize(),
                    'Repeat Unit': ssr['repeat_unit']
                })
                         
            df_grafico = pd.DataFrame(dados_grafico)

            # 4. Create the chart with only the simple SSRs
            chart_simples = alt.Chart(df_grafico).mark_rect(height=20, color='#1f77b4').encode(
                x=alt.X('Start',
                        scale=alt.Scale(domain=[1, sequence_length]),
                        axis=alt.Axis(title="Position in the Sequence (bp)")),
                x2='End',
                y=alt.Y('Type', title="Type of SSR"),
                tooltip=['Start', 'End', 'Type', 'Repeat Unit']
            ).properties(
                title=f"SSR positions in {seq_id}"
            ).interactive()

            st.altair_chart(chart_simples, use_container_width=True)
            st.markdown("---")

        # END OF GRAPH GENERATION SECTION

        # Download button for TXT format
        txt_data = export_results(st.session_state.analysis_results, "txt", project_name)
        st.download_button(
            label="Download Results (TXT)",
            data=txt_data,
            file_name=f"motif_quest_results_{project_name or 'default'}.txt",
            mime="text/plain",
            key="download_txt_results"
        )
                
elif page == "Usage & Parameters":
    st.markdown("""
    ### How to Use Motif-Quest
    
    1. **Upload or Paste Sequence**: 
       - Upload a FASTA file or paste DNA sequences directly into the text area.
    
    2. **Configure Parameters**:
       - Set minimum repeat counts for different motif lengths
       - Add custom motif length fields if needed
       - Adjust maximum distance for composite SSRs
    
    3. **Run Analysis**:
       - Click "Run Analysis" to process your sequences
    
    ---
                    
    ### Parameters Explanation
    - **Motif Length**: Number of nucleotides in the repeated sequence
    - **Minimum Repeats**: Minimum number of times a motif must be repeated
    - **Composite SSR Distance**: Maximum distance between two simple SSRs to be considered a composite SSR
    """, unsafe_allow_html=True)

elif page == "Output":
    st.markdown("""
    ### Results Columns
    - **Project**: Name of the project (if provided)
    - **ID**: Sequence identifier
    - **Type**: Type of microsatellite (mono, di, tri, etc.)
    - **Start**: Starting position of the microsatellite
    - **End**: Ending position of the microsatellite
    - **Repeat Unit**: The actual repeated sequence
    - **Length**: Total length of the microsatellite
    - **Category**: Simple or Composite SSR

    ---
                
    ### Export Options
    - CSV: Comma-Separated Values format
    - TXT: Plaintext table format
    - Custom: Formato personalizado com detalhes de cada SSR, como solicitado.
    """, unsafe_allow_html=True)

elif page == "About":
    
    # Layout with two columns
    col1, col2 = st.columns([3, 1])
    
    with col1:
        st.markdown("""
        ### What is Motif-Quest?
        Motif-Quest is a bioinformatics tool designed to identify and analyze microsatellites (Simple Sequence Repeats or SSRs) in DNA sequences. Microsatellites are tandem repeats of short DNA motifs (typically 1-6 bp) that are abundant in eukaryotic genomes and have important applications in genetics, breeding, and evolutionary studies.

        ---            
                    
        ### Key Features
        - Identification of perfect microsatellites with motifs from 1-6 nucleotides (mono- to hexa-nucleotides)
        - Support for custom motif lengths beyond hexanucleotides
        - Detection of compound/composite microsatellites
        - Customizable repeat thresholds for different motif lengths
        - User-friendly interface with flexible input methods
        - Comprehensive results with detailed position and pattern information
        - Export options in CSV and TXT formats

        ---            

        ### Applications
        - Marker development for genetic mapping
        - Population genetics studies
        - Genome-wide microsatellite distribution analysis
        - Evolutionary studies
        - Variety identification and fingerprinting

        ---            

        ### Developed By
        **Cassiano Eric de Carvalho Marques** Undergraduate Student in Agronomy  
        Federal University of Santa Catarina (UFSC) - Florianópolis  
        Research Intern at LFDGV (Laboratory of Plant Development Physiology and Genetics)

        """, unsafe_allow_html=True)

    
    with col2:
        
        st.write("""
        ### Contact
        For questions, suggestions, or bug reports, please contact:  
        cassiomirm1@gmail.com
        """)
    
    # Citation information
    st.subheader("How to Cite")
    st.write("""
    If you use Motif-Quest in your research, please cite:
    
    Marques, C.E.C. (2025). Motif-Quest: A comprehensive tool for microsatellite identification and analysis. LFDGV, Federal University of Santa Catarina, Florianópolis, Brazil.
    """)