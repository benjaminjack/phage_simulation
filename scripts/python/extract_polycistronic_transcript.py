from Bio import Entrez, SeqIO

Entrez.email = "benjamin.r.jack@gmail.com"

# Download T7 wild-type and E. coli K12 genbank records
handle = Entrez.efetch(db="nuccore", id=["NC_001604"], rettype="gb", retmode="text")

records = SeqIO.parse(handle, "genbank")

outfile = open("T7polycistron.fasta", "w")

# Set starting and ending genes for polycistronic transcript
start_gene = "T7p44"
end_gene = "T7p46"

for record in records:

    for feature in record.features:
        # Grab coding sequences
        if feature.type == "CDS":
            # Everything should have a locus tag
            if "locus_tag" in feature.qualifiers:
                id = feature.qualifiers["locus_tag"][0]
            # Check to see if this is the gene that we are starting or ending
            # with
            if id == start_gene:
                start = feature.location.start.position
            elif id == end_gene:
                end = feature.location.end.position
            else:
                continue


# Grab sequence based on ranges colected above
dna = record.seq[start:end]
# Write out the file
out = ">" + start_gene + ":" + end_gene + "\n" + str(dna)
outfile.write(out)
outfile.close()
handle.close()
