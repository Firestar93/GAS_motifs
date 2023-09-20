#!/bin/python3

import requests
import json
import dataclasses

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


def build_config(snp_path: str):
    snps = [SNP(line) for line in open(snp_path, "r").readlines() if not line.startswith("RSID")]

    config = {
        "genome": "hg38",
        "investigatedGenes": [],
        "investigatedDiseases": [],
        "manualTissues": [],
        "removedTissues": [],
        "investigatedSNPs": [
            snp.rsid for snp in snps
        ]
    }

    return {"config": config}


def create_link(config: {}):
    response = requests.post(
        "https://api.epistasis-disease-atlas.com/save_request_config", json = config
    )

    return "https://exbio.wzw.tum.de/genome-browser?uid=" + response.text.strip('"')


def create_link_snps(snp_path: str):
    config = build_config(snp_path)

    with open("config.json", "w") as f:
        json.dump(config, f)

    return create_link(config)

if __name__ == "__main__":
    print(create_link_snps("/Users/chhatralasd/Downloads/found_snp_motif/chr6_snp_destroy_motif.tsv"))