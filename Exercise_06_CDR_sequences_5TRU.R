# ============================================================
#  Extract CDR sequences from 5TRU (ipilimumab + CTLA-4)
#  Then open in ChimeraX for proper 3D visualization
# ============================================================

# --- 0. Load libraries --------------------------------------
library(bio3d)


# --- 1. Download structure ----------------------------------
pdb <- read.pdb("5TRU")

# Check chains
unique(pdb$atom$chain)
# H = heavy chain, L = light chain, A = CTLA-4


# --- 2. Extract CDR sequences (Kabat numbering) -------------

# Function to get sequence of a selected region
get_seq <- function(pdb, chain, resno_range) {
  sel <- atom.select(pdb,
                     chain = chain,
                     resno = resno_range,
                     elety = "CA")
  res <- pdb$atom$resid[sel$atom]
  aa  <- aa321(res)        # convert 3-letter to 1-letter code
  paste(aa, collapse = "")
}

# Heavy chain CDRs
CDR_H1 <- get_seq(pdb, "H", 31:35)
CDR_H2 <- get_seq(pdb, "H", 50:65)
CDR_H3 <- get_seq(pdb, "H", 95:102)

# Light chain CDRs
CDR_L1 <- get_seq(pdb, "L", 24:34)
CDR_L2 <- get_seq(pdb, "L", 50:56)
CDR_L3 <- get_seq(pdb, "L", 89:97)

# --- Print all CDRs -----------------------------------------
cat("=== Ipilimumab CDR sequences (from PDB 5TRU) ===\n\n")
cat("Heavy chain:\n")
cat("  CDR-H1 (26-35) :", CDR_H1, "\n")
cat("  CDR-H2 (50-65) :", CDR_H2, "\n")
cat("  CDR-H3 (95-102):", CDR_H3, "\n\n")
cat("Light chain:\n")
cat("  CDR-L1 (24-34) :", CDR_L1, "\n")
cat("  CDR-L2 (50-56) :", CDR_L2, "\n")
cat("  CDR-L3 (89-97) :", CDR_L3, "\n")


# --- 3. Get full heavy and light chain sequences ------------
get_full_seq <- function(pdb, chain) {
  sel <- atom.select(pdb, chain = chain, elety = "CA")
  res <- pdb$atom$resid[sel$atom]
  # keep unique residues in order
  resno <- pdb$atom$resno[sel$atom]
  idx   <- !duplicated(resno)
  aa    <- aa321(res[idx])
  paste(aa, collapse = "")
}

seq_H <- get_full_seq(pdb, "H")
seq_L <- get_full_seq(pdb, "L")
seq_A <- get_full_seq(pdb, "A")

cat("\n=== Full sequences ===\n")
cat("Heavy chain :", nchar(seq_H), "aa\n", seq_H, "\n\n")
cat("Light chain :", nchar(seq_L), "aa\n", seq_L, "\n\n")
cat("CTLA-4      :", nchar(seq_A), "aa\n", seq_A, "\n\n")


# --- 4. Save PDB for ChimeraX ------------------------------
# Save the full structure
write.pdb(pdb, file = "5TRU_ipilimumab.pdb")
cat("PDB saved as: 5TRU_ipilimumab.pdb\n")
cat("Working directory:", getwd(), "\n")


# --- 5. Generate ChimeraX commands -------------------------
# Save these commands as a .cxc script — open in ChimeraX
chimerax_script <- '
# ============================================================
#  ChimeraX script — ipilimumab + CTLA-4 (5TRU)
# ============================================================

# Open structure
open 5TRU

# 1. Start clean — cartoon style for everything
hide all
show cartoon

# 2. Color the 3 chains clearly
color /H #185fa5     
color /L #1d9e75     
color /A #993c1d     

# 3. Highlight the 6 CDR loops
color /H:26-35  #f5c842    
color /H:50-65  #ef9f27    
color /H:95-102 #e24b4a    
color /L:24-34  #85b7eb    
color /L:50-56  #378add    
color /L:89-97  #0c447c    

# 4. Make CDR loops thicker
cartoon style /H:26-35   xsection round width 1.8 thick 1.8
cartoon style /H:50-65   xsection round width 1.8 thick 1.8
cartoon style /H:95-102  xsection round width 1.8 thick 1.8
cartoon style /L:24-34   xsection round width 1.8 thick 1.8
cartoon style /L:50-56   xsection round width 1.8 thick 1.8
cartoon style /L:89-97   xsection round width 1.8 thick 1.8

# 5. CTLA-4 as transparent surface so you see the interface
surface /A
transparency /A 60 target s

# 6. Zoom to the binding interface
view /H:95-102 /L:89-97 /A
'

writeLines(chimerax_script, "ipilimumab_CDR.cxc")
cat("ChimeraX script saved as: ipilimumab_CDR.cxc\n")
cat("Open ChimeraX -> File -> Open -> ipilimumab_CDR.cxc\n")
