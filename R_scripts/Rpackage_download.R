if (!require(BiocManager)) install.packages('BiocManager', repos='https://www.stats.bris.ac.uk/R/')

if (!require(optparse)) install.packages('optparse', repos='https://www.stats.bris.ac.uk/R/')

if (!require(tidyverse)) install.packages('tidyverse', repos='https://www.stats.bris.ac.uk/R/')

if (!require(magrittr)) install.packages('magrittr', repos='https://www.stats.bris.ac.uk/R/')

if (!require(readr)) install.packages('readr', repos='https://www.stats.bris.ac.uk/R/')

if (!require(purrr)) install.packages('purrr', repos='https://www.stats.bris.ac.uk/R/')

if (!require(readxl)) install.packages('readxl', repos='https://www.stats.bris.ac.uk/R/')

if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")

if (!requireNamespace("rtracklayer", quietly = TRUE)) BiocManager::install("rtracklayer")

if (!requireNamespace("GenomicRanges", quietly = TRUE)) BiocManager::install("GenomicRanges")