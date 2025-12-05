import os
import sys
import subprocess
import glob
import re
import multiprocessing
from math import floor 

# ==============================================================================
# --- CONFIGURATION ---
# --- แก้ไขค่าต่างๆ ในส่วนนี้ให้ตรงกับโปรเจกต์ของคุณ ---
# ==============================================================================

# 1. ระบุ Path ไปยังไดเรกทอรีหลักที่เก็บโฟลเดอร์ของทุกสปีชีส์
BASE_DIR = "/home_sbi_cold/salilthip.pray/Senior/Genomics/data" 
RESULT_BASE_DIR = "/home_sbi_cold/salilthip.pray/Senior/Genomics/result"


# 2. ระบุ Path ของไดเรกทอรีที่จะใช้เก็บผลลัพธ์ (สคริปต์จะสร้างให้ถ้ายังไม่มี)
QUAST_OUTPUT_DIR = os.path.join(RESULT_BASE_DIR, "QUAST_results")
BUSCO_OUTPUT_DIR = os.path.join(RESULT_BASE_DIR, "BUSCO_results")
AUGUSTUS_OUTPUT_DIR = os.path.join(RESULT_BASE_DIR, "AUGUSTUS_results")
PROTEIN_OUTPUT_DIR = os.path.join(RESULT_BASE_DIR, "Proteins_faa") # เพิ่มสำหรับเก็บไฟล์โปรตีน
CDS_OUTPUT_DIR = os.path.join(RESULT_BASE_DIR, "CDS_fasta") #
DIAMOND_OUTPUT_DIR = os.path.join(RESULT_BASE_DIR, "DIAMOND_results")
EGGNOG_OUTPUT_DIR = os.path.join(RESULT_BASE_DIR, "EGGNOG_results")

# 3. การตั้งค่า CPU& Parallel
PARALLEL_JOBS = 4    #ทำงานละประมาณ 4-5 core
TOTAL_CPU_CORE = 18  # จำนวน CPU threads ที่จะใช้

# 4. การตั้งค่า BUSCO
# (สำคัญมาก) ระบุ lineage ที่จะใช้ เช่น embryophyta_odb10, eukaryota_odb10, etc.
# ดูลิสต์ทั้งหมดได้โดยการรัน `busco --list-datasets`
BUSCO_LINEAGE_MAP = {
    "aurantiochytrium_limacinum" : "eukaryota_odb10",
    "chlorella_sorokiniana" : "chlorophyta_odb10",
}

# 5. การตั้งค่า AUGUSTUS
# (สำคัญมาก) สร้าง mapping ระหว่าง "ชื่อโฟลเดอร์" กับ "ชื่อโมเดลสปีชีส์ของ AUGUSTUS"
# คุณต้องหาชื่อโมเดลที่เหมาะสมกับสปีชีส์ของคุณ (เช่น human, arabidopsis, fly)
# Key คือชื่อโฟลเดอร์, Value คือชื่อโมเดลของ AUGUSTUS
AUGUSTUS_SPECIES_MAP = { #Choose the closest with our species
    "aurantiochytrium_limacinum": "generic",
    "chlorella_sorokiniana" : "chlamydomonas"

}
# 6. ตั้งค่า BLAST & DIAMOND
# โหลดและจัดเก็บฐานข้อมูล ??? -> จะเตรียมก่อนหน้าหรือจะโหลดมาทีเดียว
# (สำคัญ) ระบุ Path ไปยังไฟล์ฐานข้อมูลของ DIAMOND ที่สร้างด้วย 'diamond makedb'
DIAMOND_DB_PATH = "/path/to/your/database.dmnd"

# 7. ตั้งค่า EGGNoG
EGGNOG_DATA_DIR = "/path/to/eggnog-mapper/data"

# ==============================================================================
# --- SCRIPT LOGIC ---
# --- ไม่จำเป็นต้องแก้ไขโค้ดด้านล่างนี้ ---
# ==============================================================================

def run_command(command, log_file):
    """ฟังก์ชันสำหรับรัน command line และจัดการ error/logging"""
    try:
        # เขียน log ทันทีว่าเริ่มทำ
        with open(log_file, 'a') as log:
            log.write(f"COMMAND: {' '.join(command)}\n{'='*30}\n")
        
        # ใช้ Popen เพื่อให้สามารถเขียน stdout/stderr ลงไฟล์ log ได้แบบ real-time
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding='utf-8')
        
        with open(log_file, 'a') as log:
            for line in process.stdout:
                # พิมพ์ output ออกหน้าจอ (สำหรับดูความคืบหน้าแบบ real-time)
                sys.stdout.write(line)
                # เขียนลง log file
                log.write(line)
        
        process.wait()

        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, command)
        
        print(f"  > Command completed successfully for log: {log_file}\n")

    except FileNotFoundError:
        error_msg = f"  [ERROR] Command not found: {command[0]}. Is it installed and in your PATH?"
        print(error_msg)
        with open(log_file, 'a') as log:
            log.write(error_msg + "\n")
        raise
    except subprocess.CalledProcessError as e:
        error_msg = f"  [ERROR] Command failed with exit code {e.returncode}. Check log: {log_file}"
        print(error_msg)
        with open(log_file, 'a') as log:
            log.write(error_msg + "\n")
        raise
    except Exception as e:
        error_msg = f"  [ERROR] An unexpected error occurred: {e}. Check log: {log_file}"
        print(error_msg)
        with open(log_file, 'a') as log:
            log.write(error_msg + "\n")
        raise

def find_genome_file(species_dir):
    """ค้นหาไฟล์ genome ในโฟลเดอร์ (รองรับ .fa, .fna, .fasta)"""
    for ext in ("*.fasta", "*.fa", "*.fna"):
        files = glob.glob(os.path.join(species_dir, ext))
        if files:
            return files[0] # คืนค่าไฟล์แรกที่เจอ
    return None

def extract_seq(gff_file, protein_out,cds_out):
    try:
        print(f"  > Reading from: {gff_file}")
        with open(gff_file, 'r') as f:
            content = f.read()

        gene_blocks = content.split('# end gene')
        all_proteins_fasta = []
        all_cds_fasta = []

        for block in gene_blocks:
            gene_id_match = re.search(r'# start gene (\S+)', block)
            if not gene_id_match:
                continue
            gene_id = gene_id_match.group(1)

            # สกัด Protein และทำความสะอาด
            protein_match = re.search(r'# protein sequence = \[(.*?)\]', block, re.DOTALL)
            if protein_match:
                raw_sequence = protein_match.group(1)
                cleaned_sequence = re.sub(r'[\s#$]', '', raw_sequence)
                all_proteins_fasta.append(f">{gene_id}\n{cleaned_sequence}\n")
            
            # สกัด CDS และทำความสะอาด
            cds_match = re.search(r'# coding sequence = \[(.*?)\]', block, re.DOTALL)
            if cds_match:
                raw_sequence = cds_match.group(1)
                cleaned_sequence = re.sub(r'[\s#$]', '', raw_sequence)
                all_cds_fasta.append(f">{gene_id}\n{cleaned_sequence}\n")

        # เขียนไฟล์ Protein
        with open(protein_out, 'w') as f_prot:
            f_prot.write("".join(all_proteins_fasta))
        print(f"  > Created protein file: {protein_out}")

        # เขียนไฟล์ CDS
        with open(cds_out, 'w') as f_cds:
            f_cds.write("".join(all_cds_fasta))
        print(f"  > Created CDS file: {cds_out}\n")
        return True

    except Exception as e:
        print(f"  [ERROR] Could not process file {gff_file}: {e}\n")
        return False
    
def run_species_pipeline(species_name_config_tuple):
    """
    ฟังก์ชันนี้คือ Pipeline ทั้งหมดสำหรับ 1 สปีชีส์
    มันจะถูกเรียกโดย multiprocessing.Pool
    """
    species_name, config = species_name_config_tuple
    
    # ดึงค่า settings จาก config dictionary
    base_dir = config["BASE_DIR"]
    cpus_per_job = config["CPUS_PER_JOB"]
    
    # สร้าง log file หลักสำหรับสปีชีส์นี้
    # เราจะรวม log ของทุกขั้นตอนไว้ในไฟล์เดียวเพื่อให้ง่ายต่อการ debug
    species_log_dir = os.path.join(config["RESULT_BASE_DIR"], "00_Logs")
    os.makedirs(species_log_dir, exist_ok=True)
    main_log_file = os.path.join(species_log_dir, f"{species_name}_pipeline.log")

    print(f"--- 🚀 STARTING: {species_name} (Using {cpus_per_job} threads) ---")

    try:
        # --- เตรียมไฟล์เบื้องต้น ---
        species_dir = os.path.join(base_dir, species_name)
        genome_file = find_genome_file(species_dir)
        if not genome_file:
            print(f"  [WARNING] No genome file found for {species_name}. Skipping.")
            return (species_name, "Skipped - No Genome File")
        base_name = os.path.splitext(os.path.basename(genome_file))[0]

        # --- Step 1: Run QUAST ---
        print(f"  [{species_name}] Step 1: Running QUAST...")
        quast_output_path = os.path.join(config["QUAST_OUTPUT_DIR"], species_name)
        if not os.path.exists(os.path.join(quast_output_path, "report.txt")):
            os.makedirs(quast_output_path, exist_ok=True)
            command = [
                "quast.py",
                "--output-dir", quast_output_path,
                "--threads", str(cpus_per_job),
                genome_file
            ]
            run_command(command, main_log_file)
        else:
            print(f"  [{species_name}] QUAST output already exists. Skipping.")

        # --- Step 3: Run BUSCO ---
        print(f"  [{species_name}] Step 3: Running BUSCO...")
        if species_name not in config["BUSCO_LINEAGE_MAP"]:
            print(f"  [WARNING] No BUSCO lineage defined for '{species_name}'. Skipping BUSCO.")
        else:
            busco_lineage = config["BUSCO_LINEAGE_MAP"][species_name]
            augustus_model = config["AUGUSTUS_SPECIES_MAP"].get(species_name, "generic")
            busco_output_name = species_name
            
            command = [
                "busco",
                "-i", genome_file,
                "-o", busco_output_name,
                "-l", busco_lineage,
                "-m", "genome",
                "-c", str(cpus_per_job),
                "--out_path", config["BUSCO_OUTPUT_DIR"],
                "--augustus_species", augustus_model,
                "--force"
            ]
            run_command(command, main_log_file)

        # --- Step 4: Run AUGUSTUS ---
        print(f"  [{species_name}] Step 4: Running AUGUSTUS...")
        if species_name not in config["AUGUSTUS_SPECIES_MAP"]:
            print(f"  [WARNING] No AUGUSTUS species model defined for '{species_name}'. Skipping AUGUSTUS.")
            return (species_name, "Skipped - No AUGUSTUS map")
        
        augustus_model = config["AUGUSTUS_SPECIES_MAP"][species_name]
        # แก้ไข: สร้างโฟลเดอร์ย่อยสำหรับผลลัพธ์ AUGUSTUS เพื่อความเป็นระเบียบ
        augustus_species_dir = os.path.join(config["AUGUSTUS_OUTPUT_DIR"], species_name)
        os.makedirs(augustus_species_dir, exist_ok=True)
        augustus_output_file = os.path.join(augustus_species_dir, f"{base_name}.gff")

        command = [
            "augustus",
            "--species", augustus_model,
            "--outfile", augustus_output_file,
            "--gff3", "on", "--UTR", "off", "--uniqueGeneId", "true",
            "--noInFrameStop", "true", "--codingseq", "on", "--protein", "on",
            genome_file
        ]
        run_command(command, main_log_file)

        # --- Step 4.5: Extracting Protein Seq. ---
        print(f"  [{species_name}] Step 4.5: Extracting Sequences...")
        protein_output_file = os.path.join(config["PROTEIN_OUTPUT_DIR"], f"{species_name}_proteins.faa")
        cds_output_file = os.path.join(config["CDS_OUTPUT_DIR"], f"{species_name}_cds.fna")
        
        if not os.path.exists(augustus_output_file):
            print(f"  [WARNING] AUGUSTUS GFF file not found at '{augustus_output_file}'. Skipping extraction.")
            return (species_name, "Failed - AUGUSTUS GFF missing")
            
        extract_seq(augustus_output_file, protein_output_file, cds_output_file)

        # --- Step 5: Run DIAMOND ---
        print(f"  [{species_name}] Step 5: Running DIAMOND...")
        if not os.path.exists(protein_output_file) or os.path.getsize(protein_output_file) == 0:
            print(f"  [WARNING] Protein file not found or empty for '{species_name}'. Skipping DIAMOND.")
        else:
            diamond_species_dir = os.path.join(config["DIAMOND_OUTPUT_DIR"], species_name)
            os.makedirs(diamond_species_dir, exist_ok=True)
            output_diamond_file = os.path.join(diamond_species_dir, f"{species_name}_diamond.tsv")
            
            command = [
                "diamond", "blastp",
                "-d", config["DIAMOND_DB_PATH"],
                "-q", protein_output_file,
                "-o", output_diamond_file,
                "-p", str(cpus_per_job),
                "-k", "1",
                "-e", "1e-5",
                "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle"
            ]
            run_command(command, main_log_file)

        # --- Step 6: Run EggNOG-mapper ---
        print(f"  [{species_name}] Step 6: Running EggNOG-mapper...")
        if not os.path.exists(protein_output_file) or os.path.getsize(protein_output_file) == 0:
            print(f"  [WARNING] Protein file not found or empty for '{species_name}'. Skipping EggNOG.")
        else:
            eggnog_species_dir = os.path.join(config["EGGNOG_OUTPUT_DIR"], species_name)
            os.makedirs(eggnog_species_dir, exist_ok=True)
            output_prefix = os.path.join(eggnog_species_dir, species_name)
            
            command = [
                "emapper.py",
                "-i", protein_output_file,
                "-o", output_prefix,
                "--output_dir", eggnog_species_dir,
                "--data_dir", config["EGGNOG_DATA_DIR"],
                "--cpu", str(cpus_per_job),
                "-m", "diamond",
                "--force"
            ]
            run_command(command, main_log_file)
            
        print(f"--- ✅ FINISHED: {species_name} ---")
        return (species_name, "Success")

    except Exception as e:
        # ดักจับ error ที่อาจเกิดขึ้นและไม่ได้ถูกจัดการโดย run_command
        error_msg = f"--- ❌ FAILED: {species_name} with critical error: {e} ---"
        print(error_msg)
        with open(main_log_file, 'a') as log:
            log.write(f"\n{error_msg}\n")
        return (species_name, f"Failed - {e}")

def main():
    """ฟังก์ชันหลักในการรัน Pipeline"""
    print("Starting Parallel Genomics Analysis Pipeline...")
    print(f"Running up to {PARALLEL_JOBS} species in parallel.")

    # --- 1. คำนวณการจัดสรร CPU ---
    if TOTAL_CPU_CORE < PARALLEL_JOBS:
        print(f"[WARNING] TOTAL_CPU_THREADS ({TOTAL_CPU_CORE}) is less than PARALLEL_JOBS ({PARALLEL_JOBS}).")
        print("          Setting CPUs per job to 1.")
        cpus_per_job = 1
    else:
        cpus_per_job = floor(TOTAL_CPU_CORE / PARALLEL_JOBS)
    
    print(f"System will use max {TOTAL_CPU_CORE} threads.")
    print(f"Each of the {PARALLEL_JOBS} parallel jobs will be allocated {cpus_per_job} threads.")
    
    # --- 2. ค้นหาสปีชีส์ทั้งหมด (เหมือนเดิม) ---
    try:
        all_dirs = [d for d in os.listdir(BASE_DIR) if os.path.isdir(os.path.join(BASE_DIR, d))]
        # กรองเอาโฟลเดอร์ผลลัพธ์ออกไป
        output_folder_names = {
            os.path.basename(d) for d in 
            [QUAST_OUTPUT_DIR, BUSCO_OUTPUT_DIR, AUGUSTUS_OUTPUT_DIR, 
             PROTEIN_OUTPUT_DIR, CDS_OUTPUT_DIR, DIAMOND_OUTPUT_DIR, 
             EGGNOG_OUTPUT_DIR, os.path.join(RESULT_BASE_DIR, "00_Logs")]
        }
        # กรองชื่อโฟลเดอร์ที่เป็นชื่อเดียวกับโฟลเดอร์ Input (กรณี BASE_DIR = RESULTS_BASE_DIR)
        input_folder_names = {os.path.basename(RESULT_BASE_DIR)}
        
        exclude_folders = output_folder_names.union(input_folder_names)
        
        species_list = [s for s in all_dirs if s not in exclude_folders]

        if not species_list:
            print(f"[ERROR] No species directories found in {BASE_DIR}")
            print(f"  (Note: Excluding folders named: {exclude_folders})")
            sys.exit(1)
        
        print(f"Found {len(species_list)} species to process: {', '.join(species_list)}\n")
    except FileNotFoundError:
        print(f"[ERROR] The base directory '{BASE_DIR}' does not exist.")
        sys.exit(1)

    # --- 3. สร้าง Config Dictionary ---
    # เราจะส่ง dictionary นี้ไปยังทุกๆ worker process
    config = {
        "BASE_DIR": BASE_DIR,
        "RESULT_BASE_DIR": RESULT_BASE_DIR,
        "QUAST_OUTPUT_DIR": QUAST_OUTPUT_DIR,
        "BUSCO_OUTPUT_DIR": BUSCO_OUTPUT_DIR,
        "AUGUSTUS_OUTPUT_DIR": AUGUSTUS_OUTPUT_DIR,
        "PROTEIN_OUTPUT_DIR": PROTEIN_OUTPUT_DIR,
        "CDS_OUTPUT_DIR": CDS_OUTPUT_DIR,
        "DIAMOND_OUTPUT_DIR": DIAMOND_OUTPUT_DIR,
        "EGGNOG_OUTPUT_DIR": EGGNOG_OUTPUT_DIR,
        "BUSCO_LINEAGE_MAP": BUSCO_LINEAGE_MAP,
        "AUGUSTUS_SPECIES_MAP": AUGUSTUS_SPECIES_MAP,
        "DIAMOND_DB_PATH": DIAMOND_DB_PATH,
        "EGGNOG_DATA_DIR": EGGNOG_DATA_DIR,
        "CPUS_PER_JOB": cpus_per_job
    }

    # --- 4. เตรียม Tasks สำหรับ Pool ---
    # สร้าง list ของ tuples ที่จะส่งให้ worker
    # แต่ละ tuple คือ (species_name, config)
    tasks_to_run = [(species_name, config) for species_name in species_list]

    # --- 5. รัน Pool ---
    print("="*50)
    print(f"Starting Process Pool... (Processing {len(tasks_to_run)} tasks)")
    print("="*50)

    # ใช้ with-statement เพื่อให้แน่ใจว่า pool จะถูกปิดอย่างถูกต้อง
    with multiprocessing.Pool(processes=PARALLEL_JOBS) as pool:
        # .map() จะส่ง task (tuple) ไปยัง 'run_species_pipeline' ทีละตัว
        # และรอจนกว่าทุกอย่างจะเสร็จสิ้น
        results = pool.map(run_species_pipeline, tasks_to_run)

    # --- 6. สรุปผลลัพธ์ ---
    print("="*50)
    print("All tasks completed.")
    print("="*50)
    
    success_count = 0
    failed_count = 0
    for species, status in results:
        print(f"  - {species}: {status}")
        if status == "Success":
            success_count += 1
        else:
            failed_count += 1
            
    print("\n--- Pipeline Summary ---")
    print(f"Total Species:   {len(results)}")
    print(f"Succeeded:       {success_count}")
    print(f"Failed/Skipped:  {failed_count}")
    print("Pipeline finished successfully!")


# --- 7. Entry Point (สำคัญมากสำหรับ multiprocessing) ---
if __name__ == "__main__":
    # สร้างไดเรกทอรีสำหรับเก็บผลลัพธ์ทั้งหมด *ก่อน* ที่จะเริ่ม
    # เพื่อป้องกันไม่ให้ processes หลายตัวพยายามสร้างพร้อมกัน
    print("Creating output directories...")
    os.makedirs(QUAST_OUTPUT_DIR, exist_ok=True)
    os.makedirs(BUSCO_OUTPUT_DIR, exist_ok=True)
    os.makedirs(AUGUSTUS_OUTPUT_DIR, exist_ok=True)
    os.makedirs(PROTEIN_OUTPUT_DIR, exist_ok=True)
    os.makedirs(CDS_OUTPUT_DIR, exist_ok=True)
    os.makedirs(DIAMOND_OUTPUT_DIR, exist_ok=True)
    os.makedirs(EGGNOG_OUTPUT_DIR, exist_ok=True)
    os.makedirs(os.path.join(RESULT_BASE_DIR, "00_Logs"), exist_ok=True)
    
    # เริ่มการทำงานหลัก
    main()