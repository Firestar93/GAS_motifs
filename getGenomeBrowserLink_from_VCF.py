import os

import openpyxl
import pandas as pd

#!/bin/python3

import requests
import json
import dataclasses
from openpyxl.workbook import Workbook
@dataclasses.dataclass
class SNP:
    rsid: str
    chromosome: str
    start: int
    end: int
    ref: str
    alt: str

    def __init__(self, line: str):
        splitted = line.split("\t")
        self.rsid = splitted[0]
        self.chromosome = splitted[1].replace("chr", "")
        self.start = int(splitted[3])
        self.end = int(splitted[4])

        self.ref = splitted[2][0:splitted[2].find("/")]
        self.alt = splitted[2][splitted[2].find("/") + 1:]


def build_config(snps):
    #snps = [SNP(line) for line in open(snp_path, "r").readlines() if not line.startswith("RSID")]

    config = {
        "genome": "hg38",
        "investigatedGenes": [],
        "investigatedDiseases": [],
        "manualTissues": [],
        "removedTissues": [],
        "investigatedSNPs": [
            snp for snp in snps
        ]
    }

    return {"config": config}


def create_link(config: {}):
    response = requests.post(
        "https://api.epistasis-disease-atlas.com/save_request_config", json = config
    )

    return "https://exbio.wzw.tum.de/genome-browser?uid=" + response.text.strip('"')


def create_link_snps(snps):
    config = build_config(snps)

    with open("C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\found_SNPs\\PBMC_TCell\\config.json", "w") as f:
        json.dump(config, f)

    return create_link(config)

if __name__ == "__main__":

    input_file = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\found_SNPs\\PBMC_TCell\\SNPs_in_highAcetylation_AND_immuneGenes_AND_STATsignals.vcf"
    df_snps = pd.read_csv(input_file, sep='\t', comment="#", header=None)
    snps = df_snps[2].tolist()

    snps_A = snps[:len(snps) // 2]
    snps_B = snps[len(snps) // 2:]

    print(create_link_snps(snps_B))

    print("X")



    """
    input_dir="C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\found_SNPs\\old\\SurvivingSNPs_AFTERH3K27ac.bed"
    df_snps = pd.read_csv(input_dir, sep='\t')
    df_snps = df_snps.drop_duplicates(subset=['POS'])
    snps = df_snps["POS"]

    snps_excel = snps
    snps_excel=snps_excel.to_frame()
    snps_excel.to_excel("/home/markus/Dropbox/UNI/NIH_Projects/GAS_motifs_SNPs/found_SNPs/snps.xlsx")

    snps = set(snps)

    snps = list(snps)
    snps=sorted(snps)

    snps_A = snps[:len(snps) // 2]
    snps_B = snps[len(snps) // 2:]

    print(snps)

    print("X")
    #print(create_link_snps(snps_B))
    """



print("X")