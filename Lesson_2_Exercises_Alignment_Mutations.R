# ============================================================
#  R EXERCISES — Alignment, Sequence Comparison & Mutations
#  From Genome to Drug | Lesson II
#  Università degli Studi di Perugia
# ============================================================
#
#  PACKAGES NEEDED (run once):
#  install.packages("BiocManager")
#  BiocManager::install("Biostrings")
#  install.packages("ggplot2")
#
#  HOW TO USE:
#  Run each exercise block by block.
#  Read the comments before running each line.
#  Answer the questions in the console output.
# ============================================================

library(Biostrings)
library(ggplot2)

cat("Packages loaded. Let's start!\n\n")


# ============================================================
# EXERCISE 1 — Your first sequence object
# ============================================================
cat("=== EXERCISE 1: Creating a DNA sequence ===\n\n")

# In R, DNA sequences are stored as special objects
# DNAString() creates a single sequence
# DNAStringSet() creates a collection of sequences

seq1 <- DNAString("ATGCAGTACCTGAATCGG")

cat("Our sequence:\n")
print(seq1)

cat("\nLength:", length(seq1), "bases\n")
cat("Complement:\n")
print(complement(seq1))
cat("\nReverse complement (as in the other strand):\n")
print(reverseComplement(seq1))

cat("\n--- QUESTION ---\n")
cat("If this is the coding strand, what is the template strand?\n")
cat("Hint: the template strand is the reverse complement.\n\n")


# ============================================================
# EXERCISE 2 — Comparing two sequences manually
# ============================================================
cat("=== EXERCISE 2: Match, mismatch and gap by eye ===\n\n")

# Two short fragments from CTLA4 — human vs mouse signal peptide
human <- "MACLGFQRHKAQLNLATR"
mouse <- "MACLGFQRHKAQLNLATR"

# Let's introduce one difference to simulate the +49A/G SNP
# Position 17: T (Threonine) in human → A (Alanine) in mutant
human_wt  <- "MACLGFQRHKAQLNLATR"
human_mut <- "MACLGFQRHKAQLNLAAR"
#                              ^
#                        position 17: T → A

cat("Wild-type:  ", human_wt,  "\n")
cat("Mutant:     ", human_mut, "\n\n")

# Compare position by position
wt_v  <- strsplit(human_wt,  "")[[1]]
mut_v <- strsplit(human_mut, "")[[1]]
match <- ifelse(wt_v == mut_v, "|", "X")

cat("Position-by-position comparison:\n")
cat("WT:   ", paste(wt_v,  collapse=" "), "\n")
cat("      ", paste(match, collapse=" "), "\n")
cat("MUT:  ", paste(mut_v, collapse=" "), "\n\n")

n_match    <- sum(wt_v == mut_v)
n_mismatch <- sum(wt_v != mut_v)
pct_id     <- round(n_match / length(wt_v) * 100, 1)

cat("Matches:    ", n_match,    "\n")
cat("Mismatches: ", n_mismatch, "\n")
cat("% Identity: ", pct_id,     "%\n\n")

cat("--- QUESTION ---\n")
cat("This is the +49A/G SNP. The sequences are", pct_id, "% identical.\n")
cat("Yet this single mismatch reduces CTLA-4 surface\n")
cat("expression by ~30% and increases autoimmune risk.\n")
cat("Does % identity alone tell the full story? Why not?\n\n")


# ============================================================
# EXERCISE 3 — Pairwise alignment with Biostrings
# ============================================================
cat("=== EXERCISE 3: Pairwise alignment ===\n\n")

# Now we use a real alignment algorithm
# pairwiseAlignment() implements Smith-Waterman (local)
# or Needleman-Wunsch (global) with BLOSUM62

seq_human <- AAString("MACLGFQRHKAQLNLATRTWPCTLLFFLLFIPVFCK")
seq_mouse <- AAString("MACLGFQRHKAQLNLATRTWPCTLLFFLLFIPVFCK")

# Use a slightly different mouse sequence to see mismatches
seq_mouse2 <- AAString("MACLGFQRHKAQLNLATRTWPCTLLFFLLFTPAFCK")
#                                                        ^^
#                                        two differences here

cat("Human: ", as.character(seq_human),  "\n")
cat("Mouse: ", as.character(seq_mouse2), "\n\n")

# Global alignment (Needleman-Wunsch)
aln <- pairwiseAlignment(
  seq_human, seq_mouse2,
  substitutionMatrix = "BLOSUM62",
  gapOpening         = 10,
  gapExtension       = 0.5,
  type               = "global"
)

cat("=== Alignment result ===\n")
cat("Score:      ", score(aln), "\n")
cat("% Identity: ", round(pid(aln), 1), "%\n\n")

cat("Aligned sequences:\n")
cat("Human: ", as.character(pattern(aln)), "\n")
cat("Mouse: ", as.character(subject(aln)), "\n\n")

cat("--- QUESTION ---\n")
cat("Where are the differences between human and mouse?\n")
cat("Are they in the signal peptide or another region?\n\n")


# ============================================================
# EXERCISE 4 — Global vs local alignment
# ============================================================
cat("=== EXERCISE 4: Global vs local alignment ===\n\n")

cat("Global = aligns from start to end (Needleman-Wunsch)\n")
cat("Local  = finds best matching region (Smith-Waterman)\n\n")

# Sequence A: full CTLA-4 signal peptide
seq_A <- AAString("MACLGFQRHKAQLNLATRTWPCTLLFFLLFIPVFCK")

# Sequence B: only a fragment (middle portion)
seq_B <- AAString("LATRTWPCTLLFFLLFIPVF")

cat("Sequence A (full):    ", as.character(seq_A), "\n")
cat("Sequence B (fragment):", as.character(seq_B), "\n\n")

# Global alignment
aln_global <- pairwiseAlignment(seq_A, seq_B,
  substitutionMatrix = "BLOSUM62",
  gapOpening = 10, gapExtension = 0.5,
  type = "global")

# Local alignment
aln_local <- pairwiseAlignment(seq_A, seq_B,
  substitutionMatrix = "BLOSUM62",
  gapOpening = 10, gapExtension = 0.5,
  type = "local")

cat("--- GLOBAL alignment ---\n")
cat("Score:      ", score(aln_global), "\n")
cat("% Identity: ", round(pid(aln_global), 1), "%\n")
cat("Pattern:    ", as.character(pattern(aln_global)), "\n")
cat("Subject:    ", as.character(subject(aln_global)), "\n\n")

cat("--- LOCAL alignment ---\n")
cat("Score:      ", score(aln_local), "\n")
cat("% Identity: ", round(pid(aln_local), 1), "%\n")
cat("Pattern:    ", as.character(pattern(aln_local)), "\n")
cat("Subject:    ", as.character(subject(aln_local)), "\n\n")

cat("--- QUESTION ---\n")
cat("Why does local alignment give a higher % identity here?\n")
cat("When would you use global vs local alignment in practice?\n\n")


# ============================================================
# EXERCISE 5 — BLOSUM62 scoring matrix
# ============================================================
cat("=== EXERCISE 5: Understanding BLOSUM62 ===\n\n")

data(BLOSUM62)

cat("BLOSUM62 tells us how similar two amino acids are.\n")
cat("Positive score = conservative substitution (tolerated)\n")
cat("Negative score = radical substitution (potentially harmful)\n\n")

# Check scores for several substitutions
substitutions <- list(
  "Thr -> Ala (CTLA4 +49G SNP)" = c("T", "A"),
  "Leu -> Ile (conservative)"   = c("L", "I"),
  "Gly -> Trp (radical)"        = c("G", "W"),
  "Val -> Val (identity)"       = c("V", "V"),
  "Tyr -> Cys (binding domain)" = c("Y", "C")
)

cat(sprintf("%-35s  %s\n", "Substitution", "BLOSUM62 score"))
cat(strrep("-", 50), "\n")

for(name in names(substitutions)){
  pair  <- substitutions[[name]]
  score <- BLOSUM62[pair[1], pair[2]]
  interpretation <- ifelse(score > 0, "conservative",
                    ifelse(score == 0, "neutral", "radical"))
  cat(sprintf("%-35s  %3d  (%s)\n", name, score, interpretation))
}

cat("\n--- QUESTION ---\n")
cat("The CTLA4 +49G SNP (Thr->Ala) has score 0 (neutral).\n")
cat("Yet it reduces surface expression by 30%.\n")
cat("What does this tell us about computational predictions?\n\n")


# ============================================================
# EXERCISE 6 — Visualise the BLOSUM62 matrix as a heatmap
# ============================================================
cat("=== EXERCISE 6: BLOSUM62 heatmap ===\n\n")

library(reshape2)

# Convert matrix to long format for ggplot
b_df <- melt(BLOSUM62)
colnames(b_df) <- c("AA1", "AA2", "Score")

p <- ggplot(b_df, aes(x = AA2, y = AA1, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low      = "#F96167",
    mid      = "#F4F7FB",
    high     = "#00A896",
    midpoint = 0,
    name     = "Score"
  ) +
  labs(
    title    = "BLOSUM62 substitution matrix",
    subtitle = "Green = conservative | Red = radical | White = neutral",
    x        = "Amino acid",
    y        = "Amino acid"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", color = "#1A2B4A"),
    plot.subtitle = element_text(color = "#6B7280"),
    axis.text     = element_text(size = 8),
    panel.grid    = element_blank()
  )

print(p)

cat("Heatmap generated.\n\n")
cat("--- QUESTION ---\n")
cat("Look at the heatmap. Which amino acids are most\n")
cat("similar to each other (green clusters)?\n")
cat("Hint: look at V, L, I, M — what do they have in common?\n\n")


# ============================================================
# EXERCISE 7 — Human vs Mouse CTLA-4 comparison
# ============================================================
cat("=== EXERCISE 7: Human vs Mouse CTLA-4 ===\n\n")

# Signal peptide region of both species
human_ctla4 <- AAString("MACLGFQRHKAQLNLATRTWPCTLLFFLLFIPVFCK")
mouse_ctla4 <- AAString("MACLGFQRHKAQLNLATRTWPCTLLFFLLFTPAFCK")

cat("Human CTLA-4 signal peptide:\n", as.character(human_ctla4), "\n")
cat("Mouse CTLA-4 signal peptide:\n", as.character(mouse_ctla4), "\n\n")

aln_hm <- pairwiseAlignment(
  human_ctla4, mouse_ctla4,
  substitutionMatrix = "BLOSUM62",
  gapOpening = 10, gapExtension = 0.5,
  type = "global"
)

cat("Alignment score:   ", score(aln_hm), "\n")
cat("% Identity:        ", round(pid(aln_hm), 1), "%\n\n")

# Show where differences are
h_chars <- strsplit(as.character(pattern(aln_hm)), "")[[1]]
m_chars <- strsplit(as.character(subject(aln_hm)), "")[[1]]
symbols <- ifelse(h_chars == m_chars, "|", "*")

cat("Human: ", paste(h_chars,  collapse=" "), "\n")
cat("       ", paste(symbols,  collapse=" "), "\n")
cat("Mouse: ", paste(m_chars,  collapse=" "), "\n\n")
cat("(* = difference between species)\n\n")

n_diff <- sum(symbols == "*")
cat("Differences found:", n_diff, "positions\n\n")

cat("--- QUESTION ---\n")
cat("The full protein is ~76% identical between human and mouse.\n")
cat("Based on what you see here in the signal peptide region:\n")
cat("Is the signal peptide more or less conserved than the full protein?\n")
cat("What evolutionary pressure might explain this?\n\n")


# ============================================================
# EXERCISE 8 — Clinical scenario: score a mystery mutation
# ============================================================
cat("=== EXERCISE 8: Clinical scenario ===\n\n")

cat("A melanoma patient has a somatic mutation in CTLA-4:\n")
cat("position 104, Tyrosine (Y) → Cysteine (C)\n")
cat("This position is inside the Ig V-set domain —\n")
cat("exactly where ipilimumab binds.\n\n")
cat("Question for the oncologist:\n")
cat("'Will ipilimumab still work for this patient?'\n\n")

data(BLOSUM62)

# Step 1: score the mutation
wt_aa  <- "Y"   # Tyrosine — wild type
mut_aa <- "C"   # Cysteine — patient mutation
score  <- BLOSUM62[wt_aa, mut_aa]

cat("--- Step 1: BLOSUM62 score ---\n")
cat("Wild-type:  Tyrosine  (Y) — large, aromatic, forms H-bonds\n")
cat("Mutant:     Cysteine  (C) — small, forms disulfide bonds\n")
cat("BLOSUM62 score:", score, "\n\n")

if(score < 0){
  cat("Interpretation: RADICAL substitution\n")
  cat("Very different size and chemistry — likely disruptive\n\n")
} else {
  cat("Interpretation: conservative or neutral substitution\n\n")
}

# Step 2: visualise the two amino acids
aa_props <- data.frame(
  Property  = c("Molecular weight (Da)",
                "Side chain volume (A3)",
                "Hydrophobicity"),
  Tyrosine  = c(181.2, 193, -1.3),
  Cysteine  = c(121.2,  86,  2.5)
)

cat("--- Step 2: Amino acid properties ---\n")
print(aa_props)

aa_long <- reshape2::melt(aa_props, id.vars = "Property",
                          variable.name = "Residue",
                          value.name    = "Value")

aa_long <- aa_long |>
  dplyr::group_by(Property) |>
  dplyr::mutate(Value_norm = (Value - min(Value)) /
                             (max(Value) - min(Value) + 0.001)) |>
  dplyr::ungroup()

p2 <- ggplot(aa_long,
             aes(x = Property, y = Value_norm, fill = Residue)) +
  geom_col(position = "dodge", width = 0.55) +
  scale_fill_manual(
    values = c("Tyrosine" = "#0D6E8A", "Cysteine" = "#F96167"),
    name   = "Amino acid"
  ) +
  labs(
    title    = "Clinical scenario: Y104C in CTLA-4",
    subtitle = "Tyrosine (WT) vs Cysteine (patient mutation) — property comparison",
    x        = "",
    y        = "Normalised value"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", color = "#1A2B4A"),
    plot.subtitle = element_text(color = "#6B7280"),
    axis.text.x   = element_text(angle = 20, hjust = 1)
  )

print(p2)

cat("\n--- Step 3: Clinical interpretation ---\n")
cat("BLOSUM62 score:", score, "(radical substitution)\n")
cat("Location: Ig V-set domain — ipilimumab binding site\n")
cat("Chemistry: large aromatic → small sulfhydryl\n\n")
cat("Conclusion: HIGH RISK that ipilimumab will NOT bind\n")
cat("Recommendation: test binding experimentally before treating\n")
cat("Alternative: anti-PD-1 therapy (nivolumab, pembrolizumab)\n\n")

cat("--- FINAL QUESTION ---\n")
cat("This is exactly the kind of analysis done in\n")
cat("personalised oncology today — before deciding\n")
cat("which checkpoint inhibitor to use for a patient.\n")
cat("What additional experiments would you run to confirm\n")
cat("your computational prediction?\n\n")


# ============================================================
# SUMMARY
# ============================================================
cat("============================================================\n")
cat("  EXERCISES COMPLETE — what you practised today\n")
cat("============================================================\n\n")

cat("Ex 1: Creating DNA sequence objects in R\n")
cat("Ex 2: Match, mismatch and gap — position by position\n")
cat("Ex 3: Pairwise alignment with BLOSUM62\n")
cat("Ex 4: Global vs local alignment — when to use each\n")
cat("Ex 5: BLOSUM62 scoring — conservative vs radical\n")
cat("Ex 6: BLOSUM62 heatmap — visualising amino acid similarity\n")
cat("Ex 7: Human vs mouse CTLA-4 — cross-species comparison\n")
cat("Ex 8: Clinical scenario — mutation impact on drug binding\n\n")

cat("All exercises used CTLA-4 — the gene behind ipilimumab,\n")
cat("the first immune checkpoint inhibitor (FDA 2011),\n")
cat("and the Nobel Prize in Medicine 2018.\n\n")

cat("============================================================\n")
