import os
import subprocess
import csv
import multiprocessing
import sys
import shutil

# ==============================================================================
# 1. การตั้งค่าโปรเจกต์ (PROJECT SETUP)
# ==============================================================================

# --- กำหนดจำนวน "คนงาน" ที่จะรันพร้อมกัน ---
# นี่คือจำนวน SRA ID ที่จะประมวลผลพร้อมกันในแต่ละขั้นตอน
NUM_PARALLEL_JOBS = 4

# --- กำหนด Path หลัก ---
BASE_DIR = os.getcwd()
OUTPUT_DIR = os.path.join(BASE_DIR, "analysis_output")
REF_DIR = os.path.join(BASE_DIR, "reference_data")
SAMPLE_SHEET_FILE = os.path.join(BASE_DIR, "samples.csv")
ADAPTER_FILE_PATH = os.path.join(REF_DIR, "TruSeq3-SE.fa") 

# --- สร้าง Directories หลัก (สำหรับเก็บผลลัพธ์) ---
os.makedirs(OUTPUT_DIR, exist_ok=True)
for subdir in ["sra", "fastq_raw", "fastqc_raw", "fastq_trimmed"]:
    os.makedirs(os.path.join(OUTPUT_DIR, subdir), exist_ok=True)

# ==============================================================================
# 2. ฟังก์ชันช่วยรันคำสั่ง (HELPER FUNCTION)
# ==============================================================================

def execute_command(command_list, description, sra_id):
    """ฟังก์ชันรันคำสั่งพร้อม Logging"""
    log_prefix = f"[{sra_id}]"
    print(f"\n{log_prefix} 🚀 Starting: {description}...")
    print(f"{log_prefix}    Command: {' '.join(command_list)}")

    try:
        # ใช้ capture_output=True เพื่อไม่ให้ log ของทุก job ปนกันมั่วในหน้าจอหลัก
        result = subprocess.run(command_list, check=True, text=True, 
                                executable=None, capture_output=True, timeout=3600) # 1 hour timeout
        print(f"{log_prefix} ✅ Finished: {description} successfully.")
        return
    except subprocess.CalledProcessError as e:
        print(f"❌ ERROR in '{description}' for {sra_id}: {e}")
        # พิมพ์ 5 บรรทัดสุดท้ายของ Stderr เพื่อ Debug
        print(f"{log_prefix} STDERR: ...\n" + "\n".join(e.stderr.splitlines()[-5:]))
        raise e
    except subprocess.TimeoutExpired:
        print(f"❌ TIMEOUT: '{description}' for {sra_id} took too long.")
        raise Exception(f"Timeout on {sra_id}")

# ==============================================================================
# 3. ฟังก์ชัน "คนงาน" (WORKER FUNCTIONS) - รันแบบขนาน
# ==============================================================================

def run_qc_step(job_tuple):
    """
    คนงานสำหรับขั้นตอนที่ 1: Download, FastQC, Trimmomatic
    """
    sra_id, species_name = job_tuple
    try:
        # --- 1. กำหนด Path ---
        sra_path = os.path.join(OUTPUT_DIR, "sra")
        raw_fastq_path = os.path.join(OUTPUT_DIR, "fastq_raw")
        raw_fastqc_path = os.path.join(OUTPUT_DIR, "fastqc_raw")
        trimmed_fastq_path = os.path.join(OUTPUT_DIR, "fastq_trimmed")
        
        raw_fastq = os.path.join(raw_fastq_path, f"{sra_id}.fastq")
        trimmed_fastq = os.path.join(trimmed_fastq_path, f"{sra_id}_trimmed.fastq")

        # --- 2. Acquisition ---
        if not os.path.exists(raw_fastq):
            cmd_prefetch = ["prefetch", sra_id, "-O", sra_path]
            execute_command(cmd_prefetch, "Downloading", sra_id)
            sra_file = os.path.join(sra_path, sra_id, f"{sra_id}.sra")
            
            if not os.path.exists(sra_file):
                 sra_file = os.path.join(sra_path, f"{sra_id}.sra")
            
            cmd_dump = ["fastq-dump", "--outdir", raw_fastq_path, "--split-files", sra_file]
            execute_command(cmd_dump, "Converting to FASTQ", sra_id)
        
        # --- 3. QC Check ---
        cmd_fastqc = ["fastqc", raw_fastq, "-o", raw_fastqc_path]
        execute_command(cmd_fastqc, "Running FastQC", sra_id)

        # --- 4. QC Trim ---
        cmd_trim = [
            "trimmomatic", "SE",
            "-threads", "2",
            raw_fastq,
            trimmed_fastq,
            f"ILLUMINACLIP:{ADAPTER_FILE_PATH}:2:30:10",
            "LEADING:3",
            "TRAILING:3",
            "SLIDINGWINDOW:4:15",
            "MINLEN:36"
        ]
        execute_command(cmd_trim, "Trimming adapters", sra_id)
        
        return (sra_id, "QC_Success")
    except Exception as e:
        return (sra_id, f"QC_Failed: {e}")

def run_align_step(job_tuple):
    """
    คนงานสำหรับขั้นตอนที่ 3: STAR Alignment
    """
    sra_id, species_name = job_tuple
    try:
        # --- 1. กำหนด Path ---
        trimmed_fastq = os.path.join(OUTPUT_DIR, "fastq_trimmed", f"{sra_id}_trimmed.fastq")
        species_output_dir = os.path.join(OUTPUT_DIR, species_name)
        star_index_dir = os.path.join(REF_DIR, f"{species_name}_star_index")
        bam_path = os.path.join(species_output_dir, "bam_files")
        os.makedirs(bam_path, exist_ok=True)
        
        # --- 2. Alignment ---
        output_prefix = os.path.join(bam_path, f"{sra_id}_")
        cmd_star_align = [
            "STAR",
            "--runThreadN", "4",
            "--genomeDir", star_index_dir,
            "--readFilesIn", trimmed_fastq,
            "--outFileNamePrefix", output_prefix,
            "--outSAMtype", "BAM", "SortedByCoordinate"
        ]
        execute_command(cmd_star_align, "Aligning reads (STAR)", sra_id)
        
        return (sra_id, "Align_Success")
    except Exception as e:
        return (sra_id, f"Align_Failed: {e}")

def run_quantify_step(job_tuple):
    """
    คนงานสำหรับขั้นตอนที่ 4: htseq-count
    """
    sra_id, species_name = job_tuple
    try:
        # --- 1. กำหนด Path ---
        species_output_dir = os.path.join(OUTPUT_DIR, species_name)
        gff_file_path = os.path.join(REF_DIR, f"{species_name}.gff3")
        bam_path = os.path.join(species_output_dir, "bam_files")
        counts_path = os.path.join(species_output_dir, "counts_htseq")
        os.makedirs(counts_path, exist_ok=True)

        bam_file = os.path.join(bam_path, f"{sra_id}_Aligned.sortedByCoord.out.bam")
        count_file = os.path.join(counts_path, f"{sra_id}_counts.txt")

        # --- 2. Quantification ---
        cmd_htseq = [
            "htseq-count",
            "-f", "bam",
            "-r", "pos",
            "-s", "no",
            "--idattr=ID",
            bam_file,
            gff_file_path
        ]

        result = execute_command(cmd_htseq, "Counting reads (htseq-count)", sra_id)
        with open(count_file, 'w') as f:
            f.write(result.stdout)        
        return (sra_id, "Quant_Success")
    
    except Exception as e:
        return (sra_id, f"Quant_Failed: {e}")

# ==============================================================================
# 4. ฟังก์ชัน "ผู้จัดการ" (MANAGER FUNCTIONS)
# ==============================================================================

def build_star_indices(unique_species):
    """
    ขั้นตอนที่ 2: สร้าง STAR Index (ทีละตัว) เพื่อป้องกันการชนกัน
    """
    print("\n" + "="*70)
    print("STEP 2: Building STAR Indices (Sequential)...")
    print("="*70)
    
    for species_name in unique_species:
        print(f"--- Checking Index for: {species_name} ---")
        genome_fasta_path = os.path.join(REF_DIR, f"{species_name}.fa")
        gff_file_path = os.path.join(REF_DIR, f"{species_name}.gff3")
        star_index_dir = os.path.join(REF_DIR, f"{species_name}_star_index")
        
        if not os.path.exists(genome_fasta_path) or not os.path.exists(gff_file_path):
            print(f"  [✗] WARNING: Missing {species_name}.fa or .gff3 in {REF_DIR}. Skipping index build.")
            continue
            
        if not os.path.exists(os.path.join(star_index_dir, "SA")):
            print(f"  [i] Index not found. Building...")
            os.makedirs(star_index_dir, exist_ok=True)
            cmd_star_index = [
                "STAR",
                "--runThreadN", str(NUM_PARALLEL_JOBS * 2),
                "--runMode", "genomeGenerate",
                "--genomeDir", star_index_dir,
                "--genomeFastaFiles", genome_fasta_path,
                "--sjdbGTFfile", gff_file_path,
                "--sjdbOverhang", "99"
            ]

            print(f"  🚀 Starting: Building STAR index for {species_name}...")
            # รันแบบ List (shell=False โดยอัตโนมัติ)
            subprocess.run(cmd_star_index, check=True, text=True, 
                           capture_output=True)
            print(f"  [✓] Finished: Index for {species_name} built successfully.")

            try:
                # รัน Index build โดยตรง (ไม่ผ่าน worker)
                print(f"  🚀 Starting: Building STAR index for {species_name}...")
                subprocess.run(cmd_star_index, shell=True, check=True, text=True, 
                               executable='/bin/bash', capture_output=True)
                print(f"  [✓] Finished: Index for {species_name} built successfully.")
            except subprocess.CalledProcessError as e:
                print(f"  ❌ ERROR building index for {species_name}: {e.stderr}")
        else:
            print("  [✓] Index already exists.")

def main():
    
    # 1. เตรียม "รายชื่องาน" (Job List)
    print("="*70)
    print("Reading Job List from samples.csv...")
    print("="*70)
    
    jobs = []
    unique_species = set()
    try:
        with open(SAMPLE_SHEET_FILE, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                jobs.append((row['sra_id'], row['species_name']))
                unique_species.add(row['species_name'])
    except FileNotFoundError:
        print(f"❌ ERROR: Sample sheet not found at {SAMPLE_SHEET_FILE}")
        sys.exit(1)
        
    if not jobs:
        print("❌ ERROR: No jobs found in 'samples.csv'.")
        sys.exit(1)
        
    print(f"Found {len(jobs)} total SRA samples to process across {len(unique_species)} species.")
    
    # --- เริ่มต้น Pool ---
    pool = multiprocessing.Pool(processes=NUM_PARALLEL_JOBS)
    
    # 2. ขั้นตอนที่ 1: QC (ขนาน)
    print("\n" + "="*70)
    print(f"STEP 1: Running QC (Parallel Jobs: {NUM_PARALLEL_JOBS})...")
    print("="*70)
    qc_results = pool.map(run_qc_step, jobs)
    
    # 3. ขั้นตอนที่ 2: สร้าง Index (ลำดับ)
    build_star_indices(unique_species)
    
    # 4. ขั้นตอนที่ 3: Alignment (ขนาน)
    print("\n" + "="*70)
    print(f"STEP 3: Running Alignment (Parallel Jobs: {NUM_PARALLEL_JOBS})...")
    print("="*70)
    align_results = pool.map(run_align_step, jobs)

    # 5. ขั้นตอนที่ 4: Quantification (ขนาน)
    print("\n" + "="*70)
    print(f"STEP 4: Running Quantification (Parallel Jobs: {NUM_PARALLEL_JOBS})...")
    print("="*70)
    quant_results = pool.map(run_quantify_step, jobs)
    
    # --- ปิด Pool ---
    pool.close()
    pool.join()
    
    # 6. สรุปผลลัพธ์
    print("\n" + "="*70)
    print("🎉🎉🎉 All Pipeline Stages Finished 🎉🎉🎉")
    print("="*70)
    
    # --- สร้าง Dictionary จากผลลัพธ์เพื่อง่ายต่อการค้นหา ---
    qc_status = dict(qc_results)
    align_status = dict(align_results)
    quant_status = dict(quant_results)
    
    all_sra_ids = [job[0] for job in jobs]
    failures = []
    success_count = 0
    
    print("--- Final Job Status Summary ---")
    print(f"{'SRA ID':<12} | {'QC':<12} | {'Alignment':<12} | {'Quantify':<12}")
    print("-" * 54)

    for sra_id in all_sra_ids:
        # ดึงสถานะ
        qc_stat = qc_status.get(sra_id, 'N/A')
        align_stat = align_status.get(sra_id, 'N/A')
        quant_stat = quant_status.get(sra_id, 'N/A')

        # ตรวจสอบว่าสำเร็จทุกขั้นตอนหรือไม่
        is_qc_success = "Success" in qc_stat
        is_align_success = "Success" in align_stat
        is_quant_success = "Success" in quant_stat

        job_failed = False
        
        # ตรวจสอบความล้มเหลวทีละขั้นตอน
        # เราใช้ 'elif' เพราะถ้า QC ล้มเหลว, Alignment และ Quantify ก็ไม่ควรรัน (หรือจะล้มเหลวตาม)
        if not is_qc_success:
            failures.append((sra_id, "QC", qc_stat))
            job_failed = True
        elif not is_align_success:
            failures.append((sra_id, "Alignment", align_stat))
            job_failed = True
        elif not is_quant_success:
            failures.append((sra_id, "Quantify", quant_stat))
            job_failed = True

        # พิมพ์สรุปสถานะในตาราง
        if not job_failed:
            success_count += 1
            print(f"{sra_id:<12} | {'Success':<12} | {'Success':<12} | {'Success':<12}")
        else:
            qc_print = "Success" if is_qc_success else "FAILED"
            # ถ้า QC ล้มเหลว, Alignment จะยังไม่ถูกรัน
            align_print = "Success" if is_align_success else ("FAILED" if is_qc_success else "Not Run")
            quant_print = "Success" if is_quant_success else ("FAILED" if is_align_success else "Not Run")
            print(f"{sra_id:<12} | {qc_print:<12} | {align_print:<12} | {quant_print:<12}")


    print("-" * 54)
    print(f"\nOverall Summary: {success_count} / {len(jobs)} samples processed successfully.")
    
    # พิมพ์รายละเอียดของ SRA ID ที่ล้มเหลว
    if failures:
        print("\n--- 🔥 Failed Samples Details 🔥 ---")
        for sra_id, stage, status in failures:
            # ตัดข้อความ error ให้สั้นลง
            error_message = str(status).split('\n')[0] # เอาแค่บรรทัดแรกของ error
            print(f"  SRA ID: {sra_id}")
            print(f"  Stage : {stage}")
            print(f"  Error : {error_message}...")
            print("-" * 30)

# --- รันสคริปต์ ---
if __name__ == "__main__":
    main()