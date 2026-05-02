import os
import random

# MECP2 Gene Coordinates (hg38) on chrX
CHR = "chrX"
START = 153287263
END = 153363188

def simulate_patient_data():
    print(" Phase 2 (Realistic): Simulating 'Messy' Illumina Data...")
    
    # Create necessary folders
    os.makedirs("data/patient", exist_ok=True)
    
    # 1. Extract healthy MECP2 region (unchanged)
    # This is our baseline "Control"
    os.system(f"samtools faidx data/reference/chrX.fa {CHR}:{START}-{END} > data/patient/healthy_mecp2.fa")
    
    print(" Healthy MECP2 sequence extracted.")
    
    # 2. Extract and create the Mutation (unchanged)
    with open("data/patient/healthy_mecp2.fa", "r") as f:
        lines = f.readlines()
        sequence = "".join(line.strip() for line in lines[1:])
    
    # Simulate the 5,000bp deletion
    deletion_start = 20000
    deletion_end = 25000
    mutated_full_sequence = sequence[:deletion_start] + sequence[deletion_end:]
    
    # 3. Simulate "Illumina" Short-Reads (THE NEW PART!)
    print(f"  Shredding patient DNA into 50,000 short, random reads...")
    
    num_reads = 50000     # Generate a lot of data for realistic coverage
    read_length = 150    # Standard Illumina read length
    
    with open("data/patient/patient_mecp2_shredded.fa", "w") as f:
        seq_len = len(mutated_full_sequence)
        
        for i in range(num_reads):
            # Pick a random starting point for the read
            start_pos = random.randint(0, seq_len - read_length - 1)
            end_pos = start_pos + read_length
            
            # Extract the short read and write it to the file
            read_seq = mutated_full_sequence[start_pos:end_pos]
            f.write(f">read_{i}\n")
            f.write(read_seq + "\n")
            
    print(f"  MUTATION SIMULATED: Deleted 5000bp, shredded into 50k reads.")
    print(" Saved as data/patient/patient_mecp2_shredded.fa")

if __name__ == "__main__":
    # Ensure raw reference index exists (important for extracting healthy MECP2)
    if not os.path.exists("data/reference/chrX.fa.fai"):
        print(" Missing data/reference/chrX.fa.fai index. Please run samtools faidx data/reference/chrX.fa first!")
    else:
        simulate_patient_data()