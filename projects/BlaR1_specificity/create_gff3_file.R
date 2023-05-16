# Load libraries
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))


# === Functions ================================================================

create_granges <- function(file_path) {
    # Generates a GRanges object from the csv data frame.
    
    # Reads the csv file
    # --> Converts to data frame
    annotation_data <- read.csv(file = file_path, sep = ',')
    
    # Creates a GRanges object from the csv data frame
    gr <- GRanges(seqnames = annotation_data$seqname,
                  source = annotation_data$source,
                  ranges = IRanges(start = annotation_data$start,
                                   end = annotation_data$end),
                  score = annotation_data$score,
                  strand = annotation_data$strand,
                  domain = annotation_data$domain)
    
    return(gr)
}


# === Script ===================================================================

# Reads cmd line arguments
args <- commandArgs(trailingOnly = TRUE)

# Only one arguments can be provided
# --> File path
if (length(args) == 0) {
    
    print("Error: No arguments provided...")
    
} else if (length(args) == 1) {
    
    file_path <- args
    # Normalizes the path for every operating system
    file_path <- normalizePath(path = file_path)
    
    # Creates the GRanges
    gr <- create_granges(file_path = file_path)
    
    # Create export file
    file_name <- paste0(strsplit(basename(args), '\\.')[[1]][1], ".gff3")
    export_path <- file.path(dirname(args), file_name)
    # Export with gff3 file format
    export(gr, con = export_path, format = 'GFF3')
    
} else {
    
    print("Error: Passed too many arguments...")
    
}
