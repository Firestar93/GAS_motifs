import requests
import csv
import xml.etree.ElementTree as ET
import time

# NCBI requires an email parameter in queries
EMAIL = "youremail@example.com"


def fetch_clinvar_significance(rs_id):
    """
    Query ClinVar for a given rsID and return the first reported clinical significance.
    If none found, return 'N/A'.
    """

    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    # 1) esearch: find relevant ClinVar IDs for the variant
    esearch_params = {
        "db": "clinvar",
        "term": rs_id,
        "retmax": 10,  # adjust as needed (some rsIDs might have multiple hits)
        "email": "markus.hoffmann@nih.gov"
    }

    search_response = requests.get(base_url + "esearch.fcgi", params=esearch_params)
    search_root = ET.fromstring(search_response.text)

    # Extract the first ClinVar ID (if any)
    ids = search_root.find("IdList")
    if not ids or len(list(ids)) == 0:
        return "N/A"  # No ClinVar records for this rsID

    # We'll just take the first ID for brevity
    record_id = ids[0].text

    # 2) esummary: get details about that ClinVar record
    esummary_params = {
        "db": "clinvar",
        "id": record_id,
        "retmode": "xml",
        "email": EMAIL
    }

    summary_response = requests.get(base_url + "esummary.fcgi", params=esummary_params)
    summary_root = ET.fromstring(summary_response.text)

    # Typically the clinical significance is found in a <Item Name="clinical_significance"> node
    significance_elem = summary_root.find(".//Item[@Name='clinical_significance']")
    if significance_elem is None or not significance_elem.text:
        return "N/A"

    return significance_elem.text.strip()


###############################################################################
# FULL VARIANT LIST: rsID + Gene
# If you have more or fewer variants, edit this list accordingly.
###############################################################################
variants = [
    {"rsID": "rs762389787", "Gene": "CSF3R"},
    {"rsID": "rs1209949613", "Gene": "CD101"},
    {"rsID": "rs566449782", "Gene": "RORC"},
    {"rsID": "rs762957571", "Gene": "ARHGEF2"},
    {"rsID": "rs1320629886", "Gene": "LBR"},
    {"rsID": "rs1257658099", "Gene": "IL7R"},
    {"rsID": "rs1412706856", "Gene": "ERBIN"},
    {"rsID": "rs556295198", "Gene": "ERBIN"},
    {"rsID": "rs1166900155", "Gene": "GCNT4"},
    {"rsID": "rs538911235", "Gene": "MSH3"},
    {"rsID": "rs1211652644", "Gene": "SLC12A2"},
    {"rsID": "rs1300469777", "Gene": "SLC12A2"},
    {"rsID": "rs1244482467", "Gene": "SLC12A2"},
    {"rsID": "rs908390275", "Gene": "TCF7"},
    {"rsID": "rs1313579135", "Gene": "TCF7"},
    {"rsID": "rs1407404362", "Gene": "TCF7"},
    {"rsID": "rs1485439074", "Gene": "ANKHD1"},
    {"rsID": "rs898478536", "Gene": "ITK"},
    {"rsID": "rs770307514", "Gene": "ITK"},
    {"rsID": "rs1340591740", "Gene": "IL12B"},
    {"rsID": "rs1054504837", "Gene": "DOCK2"},
    {"rsID": "rs539020058", "Gene": "PRR7"},
    {"rsID": "rs1328443888", "Gene": "DUSP22"},
    {"rsID": "rs964783222", "Gene": "WRNIP1"},
    {"rsID": "rs1290995737", "Gene": "TMEM14C"},
    {"rsID": "rs1453228287", "Gene": "TMEM14C"},
    {"rsID": "rs148371811", "Gene": "JARID2"},
    {"rsID": "rs1356982733", "Gene": "GNL1"},
    {"rsID": "rs1480616251", "Gene": "DHX16"},
    {"rsID": "rs1460296291", "Gene": "TAPBP"},
    {"rsID": "rs940616305", "Gene": "CGAS"},
    {"rsID": "rs1403816118", "Gene": "FYN"},
    {"rsID": "rs1484558508", "Gene": "VNN1"},
    {"rsID": "rs56273545", "Gene": "IFNGR1"},
    {"rsID": "rs191964705", "Gene": "IFNGR1"},
    {"rsID": "rs999557026", "Gene": "IFNGR1"},
    {"rsID": "rs1243423407", "Gene": "CITED2"},
    {"rsID": "rs889130510", "Gene": "CITED2"},
    {"rsID": "rs1208491606", "Gene": "EZR"},
    {"rsID": "rs1218894838", "Gene": "PMS2"},
    {"rsID": "rs926571814", "Gene": "RAC1"},
    {"rsID": "rs1286087613", "Gene": "AHR"},
    {"rsID": "rs1328762768", "Gene": "SKAP2"},
    {"rsID": "rs1273816947", "Gene": "ANLN"},
    {"rsID": "rs17496067", "Gene": "CDK13"},
    {"rsID": "rs891335123", "Gene": "DBNL"},
    {"rsID": "rs891975351", "Gene": "DBNL"},
    {"rsID": "rs562693929", "Gene": "LAT2"},
    {"rsID": "rs1242612897", "Gene": "SLC25A40"},
    {"rsID": "rs1339419515", "Gene": "BPGM"},
    {"rsID": "rs1341587730", "Gene": "BPGM"},
    {"rsID": "rs528194510", "Gene": "CNOT4"},
    {"rsID": "rs1262912023", "Gene": "KIF13B"},
    {"rsID": "rs1328504999", "Gene": "ASH2L"},
    {"rsID": "rs1023605113", "Gene": "RIPK2"},
    {"rsID": "rs1185793126", "Gene": "RPL30"},
    {"rsID": "rs1049184903", "Gene": "RPL30"},
    {"rsID": "rs146201851", "Gene": "KLF10"},
    {"rsID": "rs370669851", "Gene": "JAK2"},
    {"rsID": "rs1047339704", "Gene": "CD274"},
    {"rsID": "rs141486767", "Gene": "CD274"},
    {"rsID": "rs559384562", "Gene": "CD274"},
    {"rsID": "rs1389497816", "Gene": "RRAGA"},
    {"rsID": "rs570659761", "Gene": "C9orf72"},
    {"rsID": "rs1340280311", "Gene": "RIGI"},
    {"rsID": "rs183974004", "Gene": "SIT1"},
    {"rsID": "rs1164988424", "Gene": "TRIM14"},
    {"rsID": "rs552879920", "Gene": "SLC46A2"},
    {"rsID": "rs527767092", "Gene": "SLC46A2"},
    {"rsID": "rs760432279", "Gene": "SLC46A2"},
    {"rsID": "rs1166368265", "Gene": "TRIM32"},
    {"rsID": "rs1343459382", "Gene": "ABL1"},
    {"rsID": "rs1408230109", "Gene": "VIM"},
    {"rsID": "rs1318433111", "Gene": "VIM"},
    {"rsID": "rs1354909126", "Gene": "ITGB1"},
    {"rsID": "rs1288081410", "Gene": "SHLD2"},
    {"rsID": "rs1343658469", "Gene": "C12orf4"},
    {"rsID": "rs1156445361", "Gene": "C12orf4"},
    {"rsID": "rs537252651", "Gene": "CD27"},
    {"rsID": "rs746108614", "Gene": "PHB2"},
    {"rsID": "rs782801340", "Gene": "C1RL"},
    {"rsID": "rs1235159557", "Gene": "OLR1"},
    {"rsID": "rs1353819412", "Gene": "LRRK2"},
    {"rsID": "rs369581626", "Gene": "LMBR1L"},
    {"rsID": "rs372908103", "Gene": "ZNF385A"},
    {"rsID": "rs1052974951", "Gene": "TESPA1"},
    {"rsID": "rs548115646", "Gene": "RAB5B"},
    {"rsID": "rs11176083", "Gene": "IRAK3"},
    {"rsID": "rs889545139", "Gene": "IRAK3"},
    {"rsID": "rs1274171486", "Gene": "LYZ"},
    {"rsID": "rs969252028", "Gene": "EIF2B1"},
    {"rsID": "rs1162946510", "Gene": "RB1"},
    {"rsID": "rs1489272927", "Gene": "TRIM13"},
    {"rsID": "rs934803990", "Gene": "GPR183"},
    {"rsID": "rs1234876633", "Gene": "GPR183"},
    {"rsID": "rs1333515728", "Gene": "ANG"},
    {"rsID": "rs1320083106", "Gene": "RNASE3"},
    {"rsID": "rs1323186370", "Gene": "RNASE3"},
    {"rsID": "rs138171691", "Gene": "RNASE3"},
    {"rsID": "rs1016224188", "Gene": "CEBPE"},
    {"rsID": "rs1023572517", "Gene": "CNIH1"},
    {"rsID": "rs1236319671", "Gene": "GCH1"},
    {"rsID": "rs1298660313", "Gene": "ZFP36L1"},
    {"rsID": "rs917109425", "Gene": "BATF"},
    {"rsID": "rs1210326020", "Gene": "YY1"},
    {"rsID": "rs1225476024", "Gene": "RCOR1"},
    {"rsID": "rs570008008", "Gene": "KLF13"},
    {"rsID": "rs1298690347", "Gene": "TYRO3"},
    {"rsID": "rs1483085124", "Gene": "PPIB"},
    {"rsID": "rs1263545453", "Gene": "CSK"},
    {"rsID": "rs74025847", "Gene": "CSK"},
    {"rsID": "rs957501412", "Gene": "ISG20"},
    {"rsID": "rs541463494", "Gene": "LRRK1"},
    {"rsID": "rs565093000", "Gene": "MEFV"},
    {"rsID": "rs910130021", "Gene": "SOCS1"},
    {"rsID": "rs1285344328", "Gene": "IL4R"},
    {"rsID": "rs139142497", "Gene": "IL4R"},
    {"rsID": "rs1007157254", "Gene": "LAT"},
    {"rsID": "rs1232784476", "Gene": "MYL11"},
    {"rsID": "rs965766271", "Gene": "NOD2"},
    {"rsID": "rs1173850428", "Gene": "TRADD"},
    {"rsID": "rs938025290", "Gene": "PSMB10"},
    {"rsID": "rs1465699524", "Gene": "PSMB10"},
    {"rsID": "rs970365552", "Gene": "EXOSC6"},
    {"rsID": "rs1426850182", "Gene": "RNF166"},
    {"rsID": "rs1211118328", "Gene": "RNF166"},
    {"rsID": "rs566682313", "Gene": "CXCL16"},
    {"rsID": "rs1227693178", "Gene": "C1QBP"},
    {"rsID": "rs1377617251", "Gene": "EIF5A"},
    {"rsID": "rs917179073", "Gene": "GPS2"},
    {"rsID": "rs1015657660", "Gene": "GPS2"},
    {"rsID": "rs771380306", "Gene": "VAMP2"},
    {"rsID": "rs1427979634", "Gene": "PIK3R5"},
    {"rsID": "rs747210934", "Gene": "SUPT6H"},
    {"rsID": "rs73285000", "Gene": "RFFL"},
    {"rsID": "rs981669204", "Gene": "RFFL"},
    {"rsID": "rs565692092", "Gene": "NR1D1"},
    {"rsID": "rs1428746001", "Gene": "RARA"},
    {"rsID": "rs1291868776", "Gene": "TMEM106A"},
    {"rsID": "rs574137262", "Gene": "HEXIM1"},
    {"rsID": "rs898372560", "Gene": "SP2"},
    {"rsID": "rs34167845", "Gene": "CD300A"},
    {"rsID": "rs1262896680", "Gene": "CD300E"},
    {"rsID": "rs1239109578", "Gene": "UNC13D"},
    {"rsID": "rs567920860", "Gene": "SPHK1"},
    {"rsID": "rs181657872", "Gene": "JMJD6"},
    {"rsID": "rs912159055", "Gene": "ACTG1"},
    {"rsID": "rs949644659", "Gene": "ACTG1"},
    {"rsID": "rs1315780565", "Gene": "COLEC12"},
    {"rsID": "rs138606888", "Gene": "PTPN2"},
    {"rsID": "rs932943623", "Gene": "MKNK2"},
    {"rsID": "rs113402573", "Gene": "SPPL2B"},
    {"rsID": "rs1389583407", "Gene": "VAV1"},
    {"rsID": "rs144739782", "Gene": "TMED1"},
    {"rsID": "rs540787753", "Gene": "ACP5"},
    {"rsID": "rs990431105", "Gene": "KLF2"},
    {"rsID": "rs952750220", "Gene": "KLF2"},
    {"rsID": "rs571421696", "Gene": "JAK3"},
    {"rsID": "rs143632289", "Gene": "FFAR2"},
    {"rsID": "rs1268829036", "Gene": "FFAR2"},
    {"rsID": "rs1327438586", "Gene": "BCL3"},
    {"rsID": "rs948542112", "Gene": "CLPTM1"},
    {"rsID": "rs1029399914", "Gene": "BRD1"},
    {"rsID": "rs1209008068", "Gene": "BID"},
    {"rsID": "rs1327557545", "Gene": "KLRB1"},
    {"rsID": "rs1164186254", "Gene": "KLRB1"},
    {"rsID": "rs1207563094", "Gene": "TNFRSF13C"},
    {"rsID": "rs1415066403", "Gene": "CCDC134"},
    {"rsID": "rs1416479337", "Gene": "RHBDD3"},
    {"rsID": "rs1313866100", "Gene": "LGALS1"},
    {"rsID": "rs953877230", "Gene": "GTPBP1"},
    {"rsID": "rs545618650", "Gene": "GTPBP1"},
    {"rsID": "rs39519520", "Gene": "ATF4"},
    {"rsID": "rs1010340846", "Gene": "CCDC134"},
    {"rsID": "rs1415066403", "Gene": "CCDC134"},
    {"rsID": "rs1207563094", "Gene": "TNFRSF13C"},
]


def main():
    output_filename = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\GAS_motifs\\00_PUBLICATION_FOR_SUBMISSION\\CURRENTVERSION\\BMC_Genomics\\variants_with_clinvar.csv"
    fieldnames = ["rsID", "Gene", "ClinVar_Sig"]

    with open(output_filename, mode="w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        count = 0
        for variant in variants:
            count=count+1
            if count%10 == 0:
                print(count)
            rs_id = variant["rsID"]
            gene = variant["Gene"]
            time.sleep(1)
            clinvar_sig = fetch_clinvar_significance(rs_id)

            writer.writerow({
                "rsID": rs_id,
                "Gene": gene,
                "ClinVar_Sig": clinvar_sig
            })

    print(f"Done. CSV saved to {output_filename}")


if __name__ == "__main__":
    main()
