import pysam
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns # Optional, for nicer styling, install with 'pip install seaborn'

# Set consistent, professional styling (optional, but makes papers look better)
try:
    sns.set_theme(style="whitegrid")
except ImportError:
    plt.style.use('seaborn-v0_8-whitegrid') # Fallback style

# --- CONFIGURATION ---
BAM_FILE = "results/alignment_sorted.bam"
OUTPUT_DIR = "results/figures/"
MECP2_CHROM = "chrX"
# Focus just on the MECP2 region for metrics
MECP2_START = 153287263
MECP2_END = 153363188
DELETION_START = MECP2_START + 20000 # Matches simulation
DELETION_END = MECP2_START + 25000

# ------------------------------------------------------------------------------
# Figure 1: Mapping Quality (MapQ) Distribution
# Shows how confident BWA is in placing the reads.
# ------------------------------------------------------------------------------
def plot_mapq_distribution(samfile):
    print(" Plotting Figure 1: Advanced Mapping Quality (MapQ)...")
    
    mapq_scores = []
    
    for read in samfile.fetch(MECP2_CHROM, MECP2_START, MECP2_END):
        if not read.is_unmapped and not read.is_secondary:
            mapq_scores.append(read.mapping_quality)
            
    if not mapq_scores:
        print(" No valid reads found in MECP2 region.")
        return

    plt.figure(figsize=(10, 6))
    
    # 1. Draw Professional "Quality Zones" (Background Shading)
    plt.axvspan(-2, 29.5, color='#ffcccc', alpha=0.4, label='Low Quality Zone (<30)')
    plt.axvspan(29.5, 62, color='#ccffcc', alpha=0.4, label='High Quality Zone (>=30)')
    
    # 2. Plot the data using a Seaborn histogram with a smooth KDE density curve
    sns.histplot(mapq_scores, bins=range(0, 62), kde=True, color='#2c3e50', edgecolor='black', stat='density', alpha=0.7)
    
    # 3. Calculate dynamic statistics
    total_reads = len(mapq_scores)
    high_qual = sum(1 for q in mapq_scores if q >= 30)
    percent_high = (high_qual / total_reads) * 100 if total_reads > 0 else 0
    
    # 4. Add a floating statistics box
    textstr = f'Total Mapped Reads: {total_reads:,}\nHigh Quality (>=30): {percent_high:.1f}%'
    props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray')
    plt.gca().text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=12,
            verticalalignment='top', bbox=props)

    plt.title("Genomic Mapping Quality (MapQ) Density in MECP2 Region", fontsize=15, fontweight='bold')
    plt.xlabel("Phred Quality Score (MapQ)", fontsize=13)
    plt.ylabel("Read Density", fontsize=13)
    plt.xlim(-2, 62)
    
    # Move legend outside the main data area
    plt.legend(fontsize=11, loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "fig1_mapq_distribution.png"), dpi=300)
    print(f" Figure 1 saved.")

# ------------------------------------------------------------------------------

def plot_read_length(samfile):
    print(" Plotting Figure 2: Advanced Read Length Distribution...")
    
    lengths = []
    
    for read in samfile.fetch(MECP2_CHROM, MECP2_START, MECP2_END):
        if not read.is_unmapped:
            lengths.append(read.query_length)
            
    if not lengths:
        return

    plt.figure(figsize=(10, 6))
    
    # 1. Use a discrete plot for exact base-pair lengths
    sns.histplot(lengths, discrete=True, color='#3498db', edgecolor='black', alpha=0.8)
    
    # 2. Calculate exact metrics
    mean_len = np.mean(lengths)
    std_len = np.std(lengths)
    min_len = min(lengths)
    max_len = max(lengths)
    
    # 3. Add a floating summary statistics box
    textstr = f'Insert Size Metrics\n{"-"*20}\nMean Length: {mean_len:.1f} bp\nStd Deviation: {std_len:.2f} bp\nMin Length: {min_len} bp\nMax Length: {max_len} bp'
    props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray')
    plt.gca().text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=12,
            verticalalignment='top', bbox=props, family='monospace')
    
    plt.title("Post-Alignment Read Length Distribution (MECP2)", fontsize=15, fontweight='bold')
    plt.xlabel("Simulated Read Length (Base Pairs)", fontsize=13)
    plt.ylabel("Frequency (Read Count)", fontsize=13)
    
    # Add a target line
    plt.axvline(x=150, color='#e74c3c', linestyle='--', linewidth=2, label="Target Length (150bp)")
    plt.legend(fontsize=11, loc='upper right')
    
    plt.xlim(130, 170) # Zoom in slightly to make the graph look cleaner
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "fig2_read_length_check.png"), dpi=300)
    print(f" Figure 2 saved.")
# ------------------------------------------------------------------------------
# Figure 3: Genomic Coverage Heatmap (A professional summary)
# Shows the exact location of the deletion visually.
# ------------------------------------------------------------------------------
def plot_coverage_heatmap(samfile):
    print(" Plotting Figure 3: Coverage Heatmap (Summary View)...")
    
    # We create low-resolution bins to create a heatmap feel
    # Divide the 75kbp gene into 100 bins
    num_bins = 100
    gene_length = MECP2_END - MECP2_START
    bin_size = gene_length // num_bins
    
    # Calculate coverage per bin using pileup
    depths_per_pos = np.zeros(gene_length)
    for pileupcolumn in samfile.pileup(MECP2_CHROM, MECP2_START, MECP2_END, truncate=True):
        idx = pileupcolumn.pos - MECP2_START
        if 0 <= idx < gene_length:
            depths_per_pos[idx] = pileupcolumn.n
            
    # Reshape the data into 100 bins and average
    # (Note: we trim the end slightly to make it divisible by num_bins)
    depths_trimmed = depths_per_pos[:(num_bins * bin_size)]
    binned_depths = depths_trimmed.reshape(num_bins, bin_size).mean(axis=1)
    
    # Format the data for a heatmap (requires a 2D array, we'll make a 1x100 matrix)
    heatmap_data = binned_depths.reshape(1, -1)
    
    plt.figure(figsize=(16, 3))
    
    # Create the heatmap. We use a coolwarm color scheme to show contrast.
    # Reds = high coverage, Blues = low coverage
    img = plt.imshow(heatmap_data, cmap='coolwarm', aspect='auto', interpolation='nearest',
                     extent=[MECP2_START, MECP2_END, 0, 1])
    
    # Styling
    ax = plt.gca()
    ax.get_yaxis().set_visible(False) # Hide the Y-axis as it's not meaningful
    plt.colorbar(img, orientation='horizontal', label='Average Normalized Read Depth', pad=0.3)
    
    plt.title(f"Genomic Coverage Heatmap: Validating MECP2 Structural Deletion", fontsize=10, fontweight='bold')
    plt.xlabel("Genomic Coordinates (chrX)", fontsize=13)

    # Highlight the expected gap
    plt.axvline(x=DELETION_START, color='black', linestyle='--', alpha=0.7)
    plt.axvline(x=DELETION_END, color='black', linestyle='--', alpha=0.7)
    plt.annotate('DELETION ZONE', xy=((DELETION_START + DELETION_END) / 2, 0.5), xytext=((DELETION_START + DELETION_END) / 2, 1.2),
                 color='black', weight='bold', ha='center')

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "fig3_coverage_heatmap.png"), dpi=300)
    print(f" Figure 3 saved.")

# ------------------------------------------------------------------------------
# Main Execution
# ------------------------------------------------------------------------------
if __name__ == "__main__":
    # 1. Create the directory for the figures
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # 2. Open the BAM file and run the plotting functions
    try:
        # Pysam handles BAM reading and region fetching very efficiently
        sam = pysam.AlignmentFile(BAM_FILE, "rb")
        print(f" Processing {BAM_FILE} for advanced metrics...")
        
        plot_mapq_distribution(sam)
        plot_read_length(sam)
        plot_coverage_heatmap(sam)
        
        sam.close()
        print(f"\n SUCCESS! Check your 'results/figures/' folder for the 3 new research paper figures.")
        
    except (IOError, ValueError) as e:
        print(f" Error: Could not open {BAM_FILE}. Ensure you ran 'run_alignment.sh' first!")
        print(f"Detailed Error: {e}")