import pandas as pd

# get ENA IDs and FTP links
mapping = pd.read_csv("data/CDnr_to_ENA_ms_Volkova_2020.csv")
file_report_p1 = pd.read_csv("data/filereport_read_run_PRJEB2730_tsv.txt", sep="\t")
file_report_p2 = pd.read_csv("data/filereport_read_run_PRJEB4769_tsv.txt", sep="\t")
file_report = pd.concat([file_report_p1, file_report_p2])

merged = mapping.merge(
    file_report,
    left_on="Accession number",
    right_on="secondary_sample_accession",
)
merged["aligned_ftp"] = merged["submitted_ftp"].apply(lambda f: f.split(";")[0].replace("#", "%23"))
merged["is_bam"] = merged["submitted_ftp"].apply(lambda f: f.split(";")[0].endswith(".bam"))
merged.to_csv("data/combined_file_manifest.tsv", sep="\t", index=False)
