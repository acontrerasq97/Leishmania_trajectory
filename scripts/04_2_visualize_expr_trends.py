from pathlib import Path
import anndata
import matplotlib.pyplot as plt

# ============================================================
# CONFIGURATION
# ============================================================
BASE_DIR = Path("./leishmania_trajectory")
TRENDS_DIR = BASE_DIR / "results" / "wot" / "expression_trends"


# Scan all h5ad files available in the trends folder
trend_files = sorted(list(TRENDS_DIR.glob("*.h5ad")))

# Dictionary to hold file listings categorized by their source day
day_mapping = defaultdict(list)

print("-> Parsing WOT file metadata naming conventions...")

for path in trend_files:
    # Regular expression to extract the source day from 'expression_trends_day_{day}_{index}.h5ad'
    match = re.search(r"expression_trends_day_(\d+)_", path.name)
    if match:
        source_day = int(match.group(1))
        day_mapping[source_day].append(path)

# ============================================================
# Print Mapping Summary Report
# ============================================================
print("\n========================================================")
print("WOT EXPERIMENTAL TIMELINE MAPPING REPORT")
print("========================================================\n")

for day in sorted(day_mapping.keys()):
    files_in_day = day_mapping[day]
    print(f"EXPERIMENTAL DAY {day}: Total of {len(files_in_day)} sub-trajectories detected")
    print("--------------------------------------------------------")
    
    for idx, path in enumerate(files_in_day):
        # Quick-load headers to show data shape without slowing down RAM
        trend_meta = anndata.read_h5ad(path, backed="r")
        # Extract the specific population trajectory index identifier from filename
        trajectory_idx = path.stem.split("_")[-1]
        
        print(f"  • Sub-trajectory index [{trajectory_idx}]: Shape = {trend_meta.shape} | File: {path.name}")
    print("\n")


# ============================================================
# CONFIGURATION
# ============================================================
BASE_DIR = Path("./leishmania_trajectory")
TRENDS_DIR = BASE_DIR / "results" / "wot" / "expression_trends"

# Define the single specific population index you want to track across time
# Change this number to isolate other specific lineages (e.g., 0, 1, 2, etc.)
TARGET_INDEX = 0 

# Target gene of interest for molecular kinetic profiling
gene = "LmjF.01.0410"

# Scan the target directory for available h5ad files
trend_files = sorted(list(TRENDS_DIR.glob("*.h5ad")))

# Lists to compile the filtered chronological data vectors
filtered_trends = []
extracted_days = []

print(f"-> Filtering files for unique population index [{TARGET_INDEX}] across timeline...")

for path in trend_files:
    # Match the exact trailing index pattern to isolate the specific population
    # Looks for matches like '_day_{X}_{TARGET_INDEX}.h5ad'
    if path.stem.endswith(f"_{TARGET_INDEX}"):
        # Regex to safely extract the experimental day from the filename matrix
        match = re.search(r"expression_trends_day_(\d+)_", path.name)
        if match:
            source_day = int(match.group(1))
            
            # Load the anndata trend file into memory
            trend = anndata.read_h5ad(path)
            
            # Verify if the chosen gene resides in the variables slice of this matrix
            if gene in trend.var_names:
                filtered_trends.append(trend)
                extracted_days.append(source_day)
                print(f"  • Matched: Day {source_day} -> {path.name} | Shape: {trend.shape}")

# Sort vectors synchronously to ensure chronological day alignment on the X-axis
sort_indices = np.argsort(extracted_days)
extracted_days = [extracted_days[idx] for idx in sort_indices]
filtered_trends = [filtered_trends[idx] for idx in sort_indices]

# ============================================================
# Plot Isolated Population Index Kinetic Trend
# ============================================================
if len(filtered_trends) == 0:
    print(f"❌ Error: No trend files found matching target index [{TARGET_INDEX}] containing gene {gene}")
else:
    plt.figure(figsize=(8, 5), facecolor="white")
    
    # Initialize lists to gather point values for a continuous single line trajectory
    x_timeline = []
    y_expression = []
    
    for trend, day in zip(filtered_trends, extracted_days):
        # Extract the continuous expression intensity values matrix
        X_data = trend[:, gene].X
        flat_values = np.asarray(X_data).flatten()
        
        # Calculate the mean value for this specific timeline step coordinate
        mean_val = np.mean(flat_values)
        
        x_timeline.append(day)
        y_expression.append(mean_val)
        
    # Plot the consolidated single-population lineage curve
    plt.plot(
        x_timeline,
        y_expression,
        color="#2c7fb8",    # Clean high-contrast cyan/blue tone
        marker="o",
        linewidth=3,
        markersize=8,
        label=f"Isolated Lineage Index: {TARGET_INDEX}"
    )
    
    plt.xlabel("Experimental Timeline (Days)", fontsize=11)
    plt.ylabel("WOT-weighted log-expression", fontsize=11)
    plt.title(f"Single Lineage Kinetic Profile: {gene} (Index {TARGET_INDEX})", fontsize=12, fontweight="bold", pad=12)
    
    plt.xticks(x_timeline)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend(frameon=True, facecolor="white", edgecolor="#e0e0e0")
    plt.tight_layout()
    
    plt.show()
    print(f"🎉 Plot successfully rendered for single lineage index [{TARGET_INDEX}] across experimental days!")