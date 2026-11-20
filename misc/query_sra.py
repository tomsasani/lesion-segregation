import xml.etree.ElementTree as ET
from Bio import Entrez
import csv
import pandas as pd

manifest = pd.read_csv("data/combined_file_manifest.tsv", sep="\t")
print (manifest)

# Always provide your email to NCBI when using E-utilities
Entrez.email = "tomsasani1493@gmail.com"

# Replace with your actual list of hundreds of accessions
accessions = manifest["experiment_accession"].unique()#["ERX529749", "SRX123456"] 

output_file = "sra_designs.csv"

with open(output_file, mode="w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["Accession", "Design_Description"])

    # Fetching in batches of 200 is recommended to avoid timing out NCBI's servers
    batch_size = 200
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i+batch_size]

        # Join the accessions with commas
        id_list = ",".join(batch)

        try:
            # Fetch the data from the SRA database in XML format
            handle = Entrez.efetch(db="sra", id=id_list, retmode="xml")
            xml_data = handle.read()
            handle.close()

            # Parse the XML string
            root = ET.fromstring(xml_data)

            # Look for every EXPERIMENT block in the returned XML
            for experiment in root.findall(".//EXPERIMENT"):
                acc = experiment.get("accession")

                # Find the DESIGN_DESCRIPTION element
                design_desc_node = experiment.find(".//DESIGN_DESCRIPTION")
                design_desc = design_desc_node.text.strip() if design_desc_node is not None and design_desc_node.text else "No description provided"

                writer.writerow([acc, design_desc])

        except Exception as e:
            print(f"Error fetching batch starting at index {i}: {e}")

print(f"Finished! Data written to {output_file}")