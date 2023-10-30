import os
import subprocess
import sys

GAS_dir = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\GAS_motifs\\ALL_motifs\\FIMO_gainOfFunction"

for dirname in sorted(os.listdir(GAS_dir)):
    dir1 = os.path.join(GAS_dir, dirname)
    if os.path.isdir(dir1):
        print(dir1)

        #check if dir exists
        results_FIRST = os.path.join(dir1,"results_FIRST")
        results_SECOND = os.path.join(dir1,"results_SECOND")

        print(results_FIRST)
        cmd = [
            "docker", "run",
            "-v",
            r"C:\Users\hoffmannmd\OneDrive - National Institutes of Health\00_PROJECTS\GAS_motifs\ALL_motifs\FIMO_gainOfFunction\{}:/data".format(
                dirname),
            "memesuite/memesuite",
            "fimo",
            "--thresh", "0.05",
            "--no-qvalue",
            "--max-stored-scores", "900000000",
            "--norc",
            "--oc", "/data/results_FIRST",
            "/data/FIMO_input_GAS_motifs_matrix_FIRST.txt",
            "/data/{}.fa".format(dirname)
        ]
        with open("C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\ALL_motifs\\FIMO_gainOfFunction\\output.txt", "w") as f:
            subprocess.run(cmd, stdout=f, stderr=f)

        print(results_SECOND)
        cmd = [
            "docker", "run",
            "-v",
            r"C:\Users\hoffmannmd\OneDrive - National Institutes of Health\00_PROJECTS\GAS_motifs\ALL_motifs\FIMO_gainOfFunction\{}:/data".format(
                dirname),
            "memesuite/memesuite",
            "fimo",
            "--thresh", "0.05",
            "--no-qvalue",
            "--max-stored-scores", "900000000",
            "--norc",
            "--oc", "/data/results_SECOND",
            "/data/FIMO_input_GAS_motifs_matrix_SECOND.txt",
            "/data/{}.fa".format(dirname)
        ]
        with open("C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\ALL_motifs\\FIMO_gainOfFunction\\output.txt", "w") as f:
            subprocess.run(cmd, stdout=f, stderr=f)

    print("X")