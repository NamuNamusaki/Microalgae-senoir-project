import os
import csv
import glob
import sys

# ==============================================================================
# 1. การตั้งค่าโปรเจกต์ (PROJECT SETUP)
# ==============================================================================

BASE_DIR = os.getcwd()
# โฟลเดอร์ "ปลายทาง" ที่จะสร้าง Links
REF_DIR = os.path.join(BASE_DIR, "reference_data") 
# โฟลเดอร์ "ต้นทาง" ที่เก็บผล Genomics
GENOMICS_DIR = os.path.join(BASE_DIR, "../Genomics") 

# --- กำหนดไฟล์แผนที่ ---
SAMPLE_SHEET_FILE = os.path.join(BASE_DIR, "samples.csv")
GENOME_MAP_FILE = os.path.join(BASE_DIR, "genome_map.csv")

# ==============================================================================
# 2. ฟังก์ชันหลัก (Main "Bridge" Function)
# ==============================================================================

def main():
    """
    รัน "สะพาน" เพื่ออ่าน 'genome_map.csv' และ 'samples.csv' 
    และสร้าง symbolic links สำหรับไฟล์ FASTA และ GFF ที่จำเป็น
    """
    print("="*70)
    print("STEP 1: Preparing Reference Files (Linking)...")
    print("="*70)
    
    os.makedirs(REF_DIR, exist_ok=True)
    
    # 1. อ่าน genome_map.csv เก็บเป็น Dictionary
    genome_map = {}
    try:
        with open(GENOME_MAP_FILE, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                genome_map[row['species_name']] = row['assembly_accession']
    except FileNotFoundError:
        print(f"❌ ERROR: Genome map file not found at {GENOME_MAP_FILE}")
        sys.exit(1)
        
    # 2. หาสปีชีส์ทั้งหมดที่ต้องใช้ (แบบไม่ซ้ำกัน) จาก samples.csv
    unique_species = set()
    try:
        with open(SAMPLE_SHEET_FILE, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                unique_species.add(row['species_name'])
    except FileNotFoundError:
        print(f"❌ ERROR: Sample sheet not found at {SAMPLE_SHEET_FILE}")
        sys.exit(1)

    print(f"Found {len(unique_species)} unique species to prepare: {unique_species}")

    # 3. วนลูปสร้าง Links
    for species_name in unique_species:
        if species_name not in genome_map:
            print(f"  [✗] WARNING: No 'assembly_accession' found for '{species_name}' in {GENOME_MAP_FILE}. Skipping.")
            continue
            
        accession = genome_map[species_name]
        print(f"  Processing {species_name} (Accession: {accession})...")

        # 4. กำหนด Path (ตามโครงสร้างโฟลเดอร์ Genomics ล่าสุด)
        source_fasta_dir = os.path.join(GENOMICS_DIR, "Data", species_name)
        source_gff_dir = os.path.join(GENOMICS_DIR, "Result", "AUGUSTUS", species_name)

        # 5. ค้นหาไฟล์ (ใช้ glob)
        fasta_files = glob.glob(os.path.join(source_fasta_dir, "*.fna")) + \
                      glob.glob(os.path.join(source_fasta_dir, "*.fa")) + \
                      glob.glob(os.path.join(source_fasta_dir, "*.fasta"))
        source_fasta = fasta_files[0] if fasta_files else None

        gff_files = glob.glob(os.path.join(source_gff_dir, "*.gff")) + \
                    glob.glob(os.path.join(source_gff_dir, "*.gff3"))
        source_gff = gff_files[0] if gff_files else None

        # 6. กำหนดไฟล์ปลายทาง
        dest_fasta = os.path.join(REF_DIR, f"{species_name}.fa")
        dest_gff = os.path.join(REF_DIR, f"{species_name}.gff3")

        # 7. สร้าง Link (FASTA)
        if source_fasta:
            if os.path.islink(dest_fasta): os.remove(dest_fasta)
            os.symlink(os.path.abspath(source_fasta), dest_fasta)
            print(f"    [✓] Linked FASTA: {os.path.basename(source_fasta)}")
        else:
            print(f"    [✗] WARNING: FASTA file not found in {source_fasta_dir}")
            
        # 8. สร้าง Link (GFF)
        if source_gff:
            if os.path.islink(dest_gff): os.remove(dest_gff)
            os.symlink(os.path.abspath(source_gff), dest_gff)
            print(f"    [✓] Linked GFF: {os.path.basename(source_gff)}")
        else:
            print(f"    [✗] WARNING: GFF file not found in {source_gff_dir}")

    print("\nReference preparation complete.")

# --- รันสคริปต์ ---
if __name__ == "__main__":
    main()