# ============================================================
#  CTLA-4 Homology Exercise — complete script
#  Species: Human, Mouse, Rabbit, Pig
# ============================================================

# --- 0. Load libraries --------------------------------------
library(Biostrings)
library(msa)
library(ape)
library(pheatmap)
library(bio3d)


# --- 1. Fetch sequences from UniProt -----------------------
fetch_seq <- function(acc) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", acc, ".fasta")
  readAAStringSet(url)
}

human  <- fetch_seq("P16410")
mouse  <- fetch_seq("P09793")
rabbit <- fetch_seq("P42072")
pig    <- fetch_seq("Q9MYX7")

ctla4 <- c(human, mouse, rabbit, pig)
names(ctla4) <- c("Human", "Mouse", "Rabbit", "Pig")

# Check lengths — should all be ~223
width(ctla4)


# --- 2. Multiple sequence alignment -------------------------
aln <- msa(ctla4, method = "ClustalOmega")
print(aln, show = "complete")


# --- 3. Pairwise % identity matrix --------------------------
aln_mat <- as.matrix(aln)
n       <- nrow(aln_mat)
nms     <- rownames(aln_mat)
pid_mat <- matrix(0, n, n, dimnames = list(nms, nms))

for (i in 1:n) {
  for (j in 1:n) {
    same  <- sum(aln_mat[i,] == aln_mat[j,] & aln_mat[i,] != "-")
    total <- sum(aln_mat[i,] != "-" | aln_mat[j,] != "-")
    pid_mat[i,j] <- round(same / total * 100, 1)
  }
}

print(pid_mat)

# Heatmap
pheatmap(pid_mat,
         display_numbers = TRUE,
         number_format   = "%.1f%%",
         color           = colorRampPalette(c("#e6f1fb", "#185fa5"))(40),
         main            = "CTLA-4 pairwise % identity",
         border_color    = NA)


# --- 4. Phylogenetic tree -----------------------------------
aln_seqinr <- msaConvert(aln, type = "seqinr::alignment")
dist_mat   <- seqinr::dist.alignment(aln_seqinr, matrix = "identity")

tree <- nj(dist_mat)

plot(tree,
     type        = "phylogram",
     tip.color   = "#185fa5",
     edge.color  = "#444441",
     edge.width  = 2,
     cex         = 1.2,
     main        = "CTLA-4 — NJ tree (4 mammals)")
add.scale.bar(cex = 0.8)

# Save as Newick (open in iTOL: itol.embl.de)
write.tree(tree, file = "ctla4_4species.nwk")

# ============================================================
#  Done! Check online tools for 3D visualization:
#  Human (crystal):  rcsb.org/3d-view/1I85
#  Mouse (AlphaFold): alphafold.ebi.ac.uk/entry/P09793
#  Rabbit (AlphaFold): alphafold.ebi.ac.uk/entry/P42072
#  Pig (AlphaFold):  alphafold.ebi.ac.uk/entry/Q9MYX7
# ============================================================
