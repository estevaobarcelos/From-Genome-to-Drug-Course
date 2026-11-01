#===============================================================
# EXERCISE 01 — R DATA STRUCTURES
# From Genome to Drug | 2025
# Università degli Studi di Perugia
#===============================================================
# Goal: Understand how R stores and handles different types
#       of data — the foundation of all bioinformatics
#       analysis in R.
#===============================================================

#---------------------------------------------------------------
# PART 1 — DATA TYPES: typeof(), mode(), class()
#---------------------------------------------------------------
# R has three functions to inspect how data is stored.
# They look similar but answer different questions:
#
#   typeof() → How does R store the object internally?
#   mode()   → What is the general storage type? (more generic)
#   class()  → How does R treat the object for programming?
#              (relevant for object-oriented programming / OOP)

# Case 1: An integer
x <- 10L  # The "L" suffix tells R to store this as an integer
typeof(x)  # "integer" — stored as integer internally
mode(x)    # "numeric" — integers and doubles both fall under numeric
class(x)   # "integer" — treated as integer at the OOP level

# Case 2: A decimal number (double)
y <- 10.5
typeof(y)  # "double"  — stored as double-precision floating point
mode(y)    # "numeric"
class(y)   # "numeric"

# Case 3: A string (character)
s <- "Hello, it's me"
typeof(s)  # "character"
mode(s)    # "character"
class(s)   # "character"
# All three return "character" — no ambiguity here.

# Case 4: A factor
# A factor is used to represent categorical variables in R.
# It stores a limited set of distinct values called levels.
# Examples: colors, gender, treatment group, disease status.

f <- factor(c("red", "green", "blue", "yellow", "purple", "black", "orange"))
h <- factor(c("Male", "Female"))

typeof(f)  # "integer"  — factors are stored internally as integers
mode(f)    # "numeric"  — same logic
class(f)   # "factor"   — treated as a categorical variable

# Case 5: A data frame
df <- data.frame(a = c(1, 2, 3), b = c("x", "y", "z"))
typeof(df)  # "list"       — a data frame is a list of columns internally
mode(df)    # "list"
class(df)   # "data.frame" — treated as a table

# Quick reference:
# typeof() → internal storage (low-level)
# mode()   → generic storage description
# class()  → OOP type (high-level, affects how functions behave)


#---------------------------------------------------------------
# PART 2 — DATA STRUCTURES
#---------------------------------------------------------------

#---------------------------------------------------------------
# 1. Vector
#---------------------------------------------------------------
# A vector is the most basic data structure in R.
# All elements must be of the same type.

# Integer vector from 1 to 10
v <- c(1:10)
class(v)
typeof(v)
mode(v)

# Sum all elements
sum_v <- sum(v)
print(sum_v)

# Extract the third element
print(v[3])

# Double vector (contains at least one decimal → becomes numeric)
v <- c(1, 2, 3, 4.5, 5)
class(v)
typeof(v)
mode(v)

# Logical vector
v <- c(TRUE, FALSE, TRUE)
class(v)
typeof(v)
mode(v)
# Tip: TRUE = 1, FALSE = 0 — you can sum logical vectors!
sum(c(TRUE, FALSE, TRUE, TRUE))

# Complex vector
v <- c(2+3i, 5-4i, 1+0i)
class(v)
typeof(v)
mode(v)
Re(v)  # Real part
Im(v)  # Imaginary part

# Character vector (strings)
v <- c("red", "green", "blue")
class(v)
typeof(v)
mode(v)


#---------------------------------------------------------------
# 2. Matrix
#---------------------------------------------------------------
# A matrix is a 2D structure where all elements are the same type.
# Rows and columns are indexed with [row, column].

# Create a 3x3 matrix with numbers 1 to 9
m <- matrix(1:9, nrow = 3, ncol = 3)
class(m)
typeof(m)
mode(m)

# Convert matrix to data frame
df <- as.data.frame(m)
df

# Extract the third row
print(m[3, ])

# Sum of each column
col_sums <- colSums(m)
print(col_sums)


#---------------------------------------------------------------
# 3. Array
#---------------------------------------------------------------
# An array extends a matrix to more than 2 dimensions.
# Think of it as stacked matrices (e.g., 3D data).

# Create a 2x2x2 array with numbers 1 to 8
arr <- array(1:8, dim = c(2, 2, 2))
class(arr)
typeof(arr)
mode(arr)

# Extract element at position [1, 2, 2]
print(arr[1, 2, 2])

# Sum of all elements
sum_arr <- sum(arr)
print(sum_arr)


#---------------------------------------------------------------
# 4. Data Frame
#---------------------------------------------------------------
# A data frame is the most common structure for biological data.
# Each column can have a different type (numeric, character, etc.)
# Think of it as a spreadsheet or a clinical data table.

# Create a data frame with Name and Age columns
df <- data.frame(Name = c("Anna", "Luca", "Marco"), Age = c(22, 30, 25))
df
class(df)
typeof(df)
mode(df)

# Add a new row
df <- rbind(df, data.frame(Name = "Giulia", Age = 28))

# Select patients older than 25
df_over_25 <- df[df$Age > 25, ]
print(df_over_25)


#---------------------------------------------------------------
# 5. List
#---------------------------------------------------------------
# A list is the most flexible structure in R.
# It can hold elements of different types and dimensions.
# Useful for storing complex results (e.g., DESeq2 output).

# Create a list with a numeric vector, a string and a matrix
my_list <- list(
  numbers = c(1, 2, 3, 4),
  text    = "Hello R",
  matrix  = matrix(1:4, nrow = 2)
)
class(my_list)

# Extract the second element by name
print(my_list$text)

# Modify the first element
my_list$numbers <- c(10, 20, 30, 40)
print(my_list)


#---------------------------------------------------------------
# 6. Factor
#---------------------------------------------------------------
# Factors are essential for categorical variables in statistics
# and plotting (e.g., treatment group, vital status, gender).

# Create a character vector of colors
colors <- c("red", "blue", "green", "red")
class(colors)

# Convert to factor
factor_colors <- factor(colors)

# Show the levels (distinct categories)
print(levels(factor_colors))


#---------------------------------------------------------------
# 7. Tibble
#---------------------------------------------------------------
# A tibble is a modern version of a data frame from the tidyverse.
# It is cleaner, prints better, and handles data types more safely.

library(tidyverse)

# Create a tibble with Name and Age
tb <- tibble(Name = c("Anna", "Luca", "Marco"), Age = c(22, 30, 25))
class(tb)

# Filter people older than 25
tb_over_25 <- tb[tb$Age > 25, ]
print(tb_over_25)


#---------------------------------------------------------------
# 8. Function
#---------------------------------------------------------------
# Functions are first-class objects in R — you can store them,
# pass them as arguments, and inspect their class.

# Define a simple function
square <- function(x) {
  return(x^2)
}
class(square)

# Use the function
print(square(5))

#===============================================================
# END OF EXERCISE 01
#===============================================================
