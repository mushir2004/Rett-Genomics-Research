import pysam

def find_structural_variant(bam_file, chrom, start, end):
    print(f" Analyzing Read Alignments in {chrom}:{start}-{end}...")
    
    samfile = pysam.AlignmentFile(bam_file, "rb")
    
    # We will look for 'CIGAR' strings that contain 'N' or 'D' (Deletions)
    # or reads that are 'Split' across a large distance.
    found_mutation = False

    for read in samfile.fetch(chrom, start, end):
        if not read.is_unmapped:
            # CIGAR 'D' means a deletion relative to the reference
            cigar_stats = read.get_cigar_stats()[0]
            deletion_len = cigar_stats[2] # Index 2 is the 'D' operator (Deletion)
            
            if deletion_len > 1000: # We are looking for our 5000bp hit
                print(f" ALERT: Structural Variant (SV) Found!")
                print(f" Read ID: {read.query_name}")
                print(f" Deletion Size: ~{deletion_len} bp")
                print(f" CIGAR String: {read.cigarstring}")
                found_mutation = True
                break # We found it, no need to spam the terminal

    if not found_mutation:
        # Fallback: Check for a sudden drop in mean coverage
        print("⚠️ No explicit CIGAR deletion found. Checking for coverage gaps...")
        coverage = samfile.count_coverage(chrom, start, end)
        total_cov = [sum(base) for base in zip(*coverage)]
        
        # If the average coverage is near zero, it's a deletion
        avg_cov = sum(total_cov) / len(total_cov)
        if avg_cov < 0.1:
            print(f" ALERT: Possible Large Deletion detected (Average Coverage: {avg_cov:.4f})")
        else:
            print(f" Region appears normal (Average Coverage: {avg_cov:.2f})")

    samfile.close()

if __name__ == "__main__":
    BAM_PATH = "results/alignment_sorted.bam"
    find_structural_variant(BAM_PATH, "chrX", 153287263, 153363188)