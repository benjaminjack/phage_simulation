from Bio import Entrez, SeqIO

Entrez.email = "benjamin.r.jack@gmail.com"

# Download T7 wild-type and E. coli K12 genbank records
handle = Entrez.efetch(db="nuccore", id=["NC_001604"], rettype="gb", retmode="text")

records = SeqIO.parse(handle, "genbank")

outfile = open("T7_ref.fasta", "w")

for record in records:

    # Grab whole sequence for extracting records later
    seq = str(record.seq)

    for feature in record.features:
        # Grab coding, rRNA, and tRNA sequences
        if feature.type == "CDS":
            # Record some informationa about the sequence for the FASTA header
            if "gene" in feature.qualifiers:
                gene_name = feature.qualifiers["gene"][0]
            elif "product" in feature.qualifiers:
                gene_name = feature.qualifiers["product"][0]

            # Everything should have a locus tag
            if "locus_tag" in feature.qualifiers:
                id = feature.qualifiers["locus_tag"][0]

            # Grab protein ID if applicable
            prot_id = "NA"
            if "protein_id" in feature.qualifiers:
                prot_id = feature.qualifiers["protein_id"][0]

            # Construct a string in FASTA format
            out = ">"+str(record.id)+"|"+str(id)+"|"+str(prot_id)+"|"+str(feature.type)+"|"+str(gene_name)+"\n"
            out += seq[feature.location.start.position:feature.location.end.position]

            out += "\n"

            # Write FASTA sequence to file
            outfile.write(out)

outfile.close()

handle.close()
