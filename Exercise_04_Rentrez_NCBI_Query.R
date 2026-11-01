#===============================================================
# EXERCISE 04 — QUERYING NCBI DATABASES WITH RENTREZ
# From Genome to Drug | 2025
# Università degli Studi di Perugia
#===============================================================
# Package: rentrez
# Goal:    Learn to programmatically search and retrieve
#          biological sequences from NCBI databases using R —
#          the same way bioinformaticians do it in pipelines,
#          instead of downloading files manually from the web.
#===============================================================
# Databases used:
#   protein  → amino acid sequences (RefSeq)
#   nuccore  → nucleotide sequences (mRNA, genomic DNA)
#   gene     → gene information and metadata
#===============================================================

#---------------------------------------------------------------
# INSTALL AND LOAD REQUIRED PACKAGES
#---------------------------------------------------------------
install.packages("rentrez")
library(rentrez)

#---------------------------------------------------------------
# PART 1 — SEARCH THE NCBI PROTEIN DATABASE
#---------------------------------------------------------------
# We will search for the TP53 protein sequence in Homo sapiens,
# filtered to RefSeq entries only (curated, high-quality records).

# Q1. Run the search below. How many results does it return?
#     What does the [Gene Name], [Organism] and [Filter] syntax
#     mean in an NCBI query?

tp53_protein_search <- entrez_search(
  db   = "protein",
  term = "TP53[Gene Name] AND Homo sapiens[Organism] AND RefSeq[Filter]"
)

# Inspect the search result
tp53_protein_search         # Summary of the search object
tp53_protein_search$count   # Total number of results found
tp53_protein_search$ids     # List of accession IDs returned

# Q2. How many IDs were returned in $ids?
#     Is this the same as $count? Why might they differ?
#     Hint: entrez_search() returns a limited number of IDs
#     by default — check the retmax parameter.

#---------------------------------------------------------------
# PART 2 — FETCH THE PROTEIN SEQUENCE
#---------------------------------------------------------------
# Now retrieve the actual amino acid sequence for the first hit.
# rettype = "fasta" returns the sequence in FASTA format.

tp53_protein_sequence <- entrez_fetch(
  db      = "protein",
  id      = tp53_protein_search$ids[1],
  rettype = "fasta"
)

# Print the sequence to the console
cat("TP53 Protein Sequence:\n")
cat(tp53_protein_sequence)

# Q3. What is the accession number in the FASTA header?
#     Does it start with NP_ or XP_?
#     What does each prefix mean?

# Q4. How many amino acids does this sequence contain?
#     Hint: count the characters excluding the header line.
#     You can use nchar() on the sequence string.

# Q5. Is this the canonical isoform of TP53 (isoform a)?
#     How can you verify this from the FASTA header alone?

#---------------------------------------------------------------
# PART 3 — SAVE THE SEQUENCE TO A FILE
#---------------------------------------------------------------
# Save the retrieved sequence as a FASTA file.
# This is the standard format used by BLAST, ClustalOmega,
# and most bioinformatics tools.

write(tp53_protein_sequence, "TP53_protein_sequence.fasta")
cat("\nProtein sequence saved to: TP53_protein_sequence.fasta\n")

# Q6. Open the saved file in a text editor.
#     What are the two parts of a FASTA file?
#     (Header line and sequence block — describe each.)

#---------------------------------------------------------------
# PART 4 — FETCH THE mRNA SEQUENCE
#---------------------------------------------------------------
# Now retrieve the mRNA (nucleotide) sequence for TP53.
# We search the nuccore database instead of protein.

tp53_mrna_search <- entrez_search(
  db   = "nuccore",
  term = "TP53[Gene Name] AND Homo sapiens[Organism] AND mRNA[Filter] AND RefSeq[Filter]"
)

tp53_mrna_search$count
tp53_mrna_search$ids

# Fetch the first mRNA sequence
tp53_mrna_sequence <- entrez_fetch(
  db      = "nuccore",
  id      = tp53_mrna_search$ids[1],
  rettype = "fasta"
)

cat("TP53 mRNA Sequence:\n")
cat(tp53_mrna_sequence)

# Save to file
write(tp53_mrna_sequence, "TP53_mrna_sequence.fasta")
cat("\nmRNA sequence saved to: TP53_mrna_sequence.fasta\n")

# Q7. What prefix does the mRNA accession number start with?
#     (NM_ or XM_?) What does each prefix mean?

# Q8. The mRNA sequence is longer than the protein sequence.
#     Why? What additional regions does an mRNA contain
#     beyond the coding sequence (CDS)?
#     Hint: think about 5'UTR, CDS, and 3'UTR.

#---------------------------------------------------------------
# PART 5 — FETCH GENE METADATA
#---------------------------------------------------------------
# Retrieve structured information about the TP53 gene itself
# from the NCBI Gene database.

tp53_gene_search <- entrez_search(
  db   = "gene",
  term = "TP53[Gene Name] AND Homo sapiens[Organism]"
)

tp53_gene_search$ids

# Fetch gene summary in XML format
tp53_gene_summary <- entrez_summary(
  db = "gene",
  id = tp53_gene_search$ids[1]
)

# Inspect available fields
names(tp53_gene_summary)

# Print key metadata
cat("Gene name:       ", tp53_gene_summary$name, "\n")
cat("Full name:       ", tp53_gene_summary$description, "\n")
cat("Chromosome:      ", tp53_gene_summary$chromosome, "\n")
cat("Gene ID:         ", tp53_gene_summary$uid, "\n")
cat("Genomic location:", tp53_gene_summary$genomicsummary, "\n")

# Q9. On which chromosome is TP53 located?
#     What is its official NCBI Gene ID?

# Q10. Run the same gene search but change the organism
#      to "Mus musculus". What is the mouse TP53 Gene ID?
#      Is it the same as the human one?

#---------------------------------------------------------------
# PART 6 — FETCH SEQUENCES FOR MULTIPLE SPECIES
#---------------------------------------------------------------
# Instead of fetching one sequence at a time, we can loop
# over multiple species and retrieve all sequences at once.
# This is the basis of automated comparative genomics pipelines.

species_list <- c(
  "Homo sapiens",
  "Mus musculus",
  "Danio rerio",
  "Gallus gallus"
)

# Loop over each species and retrieve the TP53 protein sequence
for (sp in species_list) {
  
  cat("\n--- Fetching TP53 protein for:", sp, "---\n")
  
  search_result <- entrez_search(
    db   = "protein",
    term = paste0("TP53[Gene Name] AND ", sp, "[Organism] AND RefSeq[Filter]")
  )
  
  if (search_result$count > 0) {
    seq <- entrez_fetch(
      db      = "protein",
      id      = search_result$ids[1],
      rettype = "fasta"
    )
    cat(seq)
    
    # Save each sequence to its own file
    filename <- paste0("TP53_", gsub(" ", "_", sp), ".fasta")
    write(seq, filename)
    cat("Saved to:", filename, "\n")
    
  } else {
    cat("No RefSeq protein found for", sp, "\n")
  }
}

# Q11. Did all four species return a result?
#      Were any sequences missing? What does this tell you
#      about RefSeq annotation coverage across species?

# Q12. Open the four saved FASTA files.
#      Are all sequences the same length?
#      Which species has the most similar sequence to human?
#      How could you confirm this quantitatively?
#      Hint: think about what tool you used in Exercise 03.

#---------------------------------------------------------------
# BONUS — LINK DATABASES WITH entrez_link()
#---------------------------------------------------------------
# NCBI databases are interconnected. You can navigate from
# one database to another — for example, from a gene to
# its associated PubMed publications.

tp53_links <- entrez_link(
  dbfrom = "gene",
  id     = tp53_gene_search$ids[1],
  db     = "pubmed"
)

# How many PubMed publications are linked to TP53?
length(tp53_links$links$gene_pubmed)

# Fetch the titles of the first 5 linked papers
top5_papers <- entrez_summary(
  db = "pubmed",
  id = tp53_links$links$gene_pubmed[1:5]
)

# Print titles
for (paper in top5_papers) {
  cat("Title:", paper$title, "\n\n")
}

# BONUS Q1. How many PubMed papers are linked to TP53?
#           What does this number tell you about the
#           importance of this gene in biomedical research?

# BONUS Q2. Read the titles of the top 5 papers.
#           What are the main research topics covered?
#           Are they focused on basic biology, clinical
#           applications, or both?

#===============================================================
# END OF EXERCISE 04
#===============================================================
