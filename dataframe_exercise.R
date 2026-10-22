# ============================================================
#  EXERCISE — Understanding Data Frames in R
#  Applied to RNA-seq sample metadata
# ============================================================
# Tip: press Ctrl+Enter (Windows) or Cmd+Enter (Mac)
#      to run one line at a time in RStudio.


# ------------------------------------------------------------
# PART 1: What is a data frame?
# ------------------------------------------------------------

# A data frame is a TABLE — just like an Excel spreadsheet.
#   - Each ROW    = one observation  (here: one sample)
#   - Each COLUMN = one variable     (name, condition, counts...)
#   - Columns can be DIFFERENT types (text, numbers, categories)

# Let's build a simple RNA-seq sample table from scratch:

sample_table <- data.frame(
  sample    = c("Ctrl_1",   "Ctrl_2",   "Ctrl_3",
                "Treat_1",  "Treat_2",  "Treat_3"),
  condition = c("Control",  "Control",  "Control",
                "Treated",  "Treated",  "Treated"),
  batch     = c("B1", "B2", "B1", "B2", "B1", "B2"),
  total_reads = c(12000000, 15000000, 11000000,
                  18000000, 16000000, 17000000),
  stringsAsFactors = FALSE
)

# Print the whole table
sample_table


# ------------------------------------------------------------
# PART 2: Exploring the data frame
# ------------------------------------------------------------

nrow(sample_table)      # how many samples (rows)?
ncol(sample_table)      # how many variables (columns)?
dim(sample_table)       # both at once: rows x columns
colnames(sample_table)  # column names
str(sample_table)       # structure: type of each column
summary(sample_table)   # quick summary of every column


# ------------------------------------------------------------
# PART 3: Accessing data
# ------------------------------------------------------------

# Access a whole column with $
sample_table$condition
sample_table$total_reads

# Access one specific row (one sample)
sample_table[1, ]       # first row
sample_table[4, ]       # fourth row

# Access one specific cell [row, column]
sample_table[2, "total_reads"]   # total reads of sample 2

# Filter: keep only Treated samples
sample_table[sample_table$condition == "Treated", ]


# ------------------------------------------------------------
# PART 4: Adding a new column
# ------------------------------------------------------------

# Calculate reads in millions (easier to read)
sample_table$reads_millions <- sample_table$total_reads / 1e6
sample_table


# ------------------------------------------------------------
# PART 5: Connecting to DESeq2
# ------------------------------------------------------------

# In DESeq2, this table is called "colData".
# DESeq2 needs the condition column as a FACTOR
# so it knows which group is the reference (Control).

sample_table$condition <- factor(sample_table$condition,
                                 levels = c("Control", "Treated"))

# Control is now the reference level (= denominator of fold change)
levels(sample_table$condition)

# Set row names = sample names (required by DESeq2)
rownames(sample_table) <- sample_table$sample
sample_table

# This table is ready to be passed to DESeq2 as colData:
#
# dds <- DESeqDataSetFromMatrix(
#   countData = count_matrix,   # genes x samples
#   colData   = sample_table,   # <-- our data frame!
#   design    = ~ condition
# )
