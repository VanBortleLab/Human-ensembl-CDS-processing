# title: "MCB290: Codon Frequency Analysis"

## Section 1. Retrieve Coding Sequence Information (Annotation)



codon_table <- c(
  TTT="Phe", TTC="Phe", TTA="Leu", TTG="Leu",
  CTT="Leu", CTC="Leu", CTA="Leu", CTG="Leu",
  ATT="Ile", ATC="Ile", ATA="Ile", ATG="Met",
  GTT="Val", GTC="Val", GTA="Val", GTG="Val",
  
  TCT="Ser", TCC="Ser", TCA="Ser", TCG="Ser",
  CCT="Pro", CCC="Pro", CCA="Pro", CCG="Pro",
  ACT="Thr", ACC="Thr", ACA="Thr", ACG="Thr",
  GCT="Ala", GCC="Ala", GCA="Ala", GCG="Ala",
  
  TAT="Tyr", TAC="Tyr", TAA="Stop", TAG="Stop",
  CAT="His", CAC="His", CAA="Gln", CAG="Gln",
  AAT="Asn", AAC="Asn", AAA="Lys", AAG="Lys",
  GAT="Asp", GAC="Asp", GAA="Glu", GAG="Glu",
  
  TGT="Cys", TGC="Cys", TGA="Stop", TGG="Trp",
  CGT="Arg", CGC="Arg", CGA="Arg", CGG="Arg",
  AGT="Ser", AGC="Ser", AGA="Arg", AGG="Arg",
  GGT="Gly", GGC="Gly", GGA="Gly", GGG="Gly"
)



# so we got our cds downloaded, now let's go and 
# 1: break into triplets
# 2: divisible by 3

# THIS CODE READS A FASTA FILE AND TURNS IT INTO THE RDS FILE
# uses the library (biostrings)

library(Biostrings)
library(stringr)

# Read the FASTA file
fasta_path <- "Homo_sapiens.GRCh38.cds.all.fa"
fasta <- readDNAStringSet(fasta_path)

headers  <- names(fasta)          # full header line for each record, minus ">"
cds_seqs <- as.character(fasta)   # full CDS sequence for each record, as a string

# Parse the header fields

# TranscriptID.TranscriptIDv is always the first whitespace delimited token
transcript_full <- str_extract(headers, "^\\S+")
TranscriptID <- str_split_fixed(transcript_full, "\\.", 2)[, 1]  # base ID only
TranscriptIDv <- transcript_full # full ID.version

# gene:GeneID.GeneIDv
gene_full <- str_remove(str_extract(headers, "gene:\\S+"), "^gene:")
GeneID <- str_split_fixed(gene_full, "\\.", 2)[, 1]  # base ID only
GeneIDv <- gene_full    # full ID.version

# gene_symbol:GeneName
gene_name <- str_remove(str_extract(headers, "gene_symbol:\\S+"), "^gene_symbol:")

# Make the data frame
cds <- data.frame(
  GeneName      = gene_name,
  GeneID        = GeneID,
  GeneIDv       = GeneIDv,
  TranscriptID  = TranscriptID,
  TranscriptIDv = TranscriptIDv,
  CDS           = cds_seqs,
  stringsAsFactors = FALSE
)

cat("Rows parsed:", nrow(cds), "\n")
cat("Records missing GeneName:", sum(is.na(cds$GeneName)), "\n")

# save as the RDS file
saveRDS(cds, "Homo_sapiens_GRCh38_cds.rds")


## Section 2. Break up CDS into Individual Codons 


# CDS is similar to a human dictionary
# It's the overall sequences of codons that we need to filter

CodSeqs = cds[1,6] # row 1, column 6
NCodons = nchar(CodSeqs)/3 # it counts # of characters in a string

# so we got our cds downloaded, now let's go and 
# 1: break into triplets
# 2: divisible by 3

# now we split it into codons if it is divisible by 3

start = seq(1, nchar(CodSeqs), 3)
# starts from 1 letter, then up to how many letters there are
# and it / it by 3 (start, end, by 3). 
# this is a vector with start positions only! so 1,4,7,10...

end = seq(3, nchar(CodSeqs), 3)

# so now we have the start of the codons (1,3,7..)
# and now we have the end... (3,6,9..)

codons= substring(CodSeqs, start, end)
codons

# substring actually exports your code!


# now for ATG 
# substr(x, start, stop) extracts part of a string, so...

ATGCodon = substr(CodSeqs, 1, 3) == "ATG"

# also we can't use an if statement since R takes true/false statements one at a time 
# since cds is a spreadsheet vs line one row it won't work

# now for stop codons: TAA, TAG, TGA

endCodon = substr(CodSeqs, nchar(CodSeqs) - 2, nchar(CodSeqs))

# so it outputs codon sequences, counts them, and then it does total - 2, and then up to the end
# so like 9-2 = 7, then 7 -> 9 7,8,9 and then let's check if it's = to a stop codon

#now let's set a value to the end codons

endCodon = c("TGA", "CCC", "TAG", "AAA")

endCodon %in% c("TAA", "TAG", "TGA")

# is each codon on the left IN the list on the right?? 

# now we need to apply this to ALL of the data set

# now we split it into codons if it is divisible by 3

start = seq(1, nchar(CodSeqs), 3)

# starts from 1 letter, then up to how many letters there are
# and it / it by 3 (start, end, by 3). 
# this is a vector with start positions only! so 1,4,7,10...

end = seq(3, nchar(CodSeqs), 3)

# so now we have the start of the codons (1,3,7..)
# and now we have the end... (3,6,9..)

codons= substring(CodSeqs, start, end)
codons

# substring actually exports your code!


# now for ATG 
# substr(x, start, stop) extracts part of a string, so...

ATGCodon = substr(CodSeqs, 1, 3) == "ATG"

# also we can't use an if statement since R takes true/false statements one at a time 
# since cds is a spreadsheet vs line one row it won't work

# now for stop codons: TAA, TAG, TGA

endCodon = substr(CodSeqs, nchar(CodSeqs) - 2, nchar(CodSeqs))
# so it outputs codon sequences, counts them, and then it does total - 2, and then up to the end
# so like 9-2 = 7, then 7 -> 9 7,8,9 and then let's check if it's = to a stop codon

# now let's set a value to the end codons

endCodon = c("TGA", "CCC", "TAG", "AAA")

endCodon %in% c("TAA", "TAG", "TGA")

# is each codon on the left IN the list on the right?? that's what % mean
# % on left and right 

# now we need to apply this to ALL of the data set
# so we have multiple columns, but only the 6th one is the data that we need

cds_seqs = cds[,6]
# so all data in the 6th column

Correcting_CDS = function(seq) {
  nchar(seq) %% 3 == 0 &&
    substr(seq, 1, 3) == "ATG" &&
    substr(seq, nchar(seq) - 2, nchar(seq)) %in% c("TAA", "TAG", "TGA")
}

# so we name the function
# seq = character string which represents CDS
# function checks if the whole cds follows the rules
# nchar(seq) counts the number of nucleotides in the sequence
# nchar(seq) %% 3 == 0 this is a really cool code
# this checks for the REMAINDER to see if it's equal to 0 when we divide total 
# number of codons by 3. if it is divisible then it has a remainder of 0! 
# we have %% on the left and right side. this means if the left is true, continue to the right, next function
# so the the /3 was true, so we go to the next statement ATG
# since ATG was then right, we go to the ending codon and we finish with curly brackets.



## Section 3. Determine Codon Frequency (%)


# next, connect our command with the CDS

validCDS = sapply(cds_seqs, Correcting_CDS)

# sapply rund correct cds on each cds sequence
# we check if it's true or false

# we made it into a vector of true/false values kinda
# we take vectpr from the right side and store it in the validCDS

CleanCDS = cds_seqs[validCDS]

# [] select elements from a vector
# im telling R give me cds seqs where the value of validCDS is true 
# so it won't print false statements
# cds seq has all 6th column, while correct cds reads one line at a time
# so now we have to loop it 

ALLcodons = unlist(lapply(CleanCDS, function(s) {
  starts = seq(1, nchar(s), 3)
  ends = seq(3, nchar(s), 3)
  substring(s, starts, ends)
}))

# unlist means it combines all codons into one chunky vector
# lapply, as name implies applies function to each element
# function(seq) is an anonymous function 
# in a nutshell lapply needs to apply a function, function(seq) is it
# it means for this sequence, do these steps
# for each DNA seq, split it into codons



# we need t get the codon frequencies to see which amino acids show up with those codons 


## Section 4. Plot Codon Frequency (%) Patterns


codon_counts = table(ALLcodons)
codon_counts = sort(codon_counts, decreasing = TRUE)
print(codon_counts)

# table tallies all codons and creates a frequency table of codons
# orders the table so the most frequent codons come first

barplot(codon_counts,
  main = "Codon Sequence Counts",
  ylab = "Count",
  col = "skyblue",
  las = 2,       # Rotate x-axis labels vertically so it can be readable
  cex.names = .4)  # Make labels smaller to fit


## Section 5. Plot Amino Acid Fequency (%) Patterns


# NOW THE AMINO ACIDS

# we will be using the codon and amino acid chart from the very top of the coding 
# I tried to use bio strings but it kept telling me how it was not 
# loaded properly

ALLcodons = toupper(ALLcodons)
AAs = codon_table[ALLcodons]

AA_counts = sort(table(AAs), decreasing = TRUE)

barplot(AA_counts,
  main = "Amino Acid Counts",
  ylab = "Count",
  col = "orange",
  las = 2,
  cex.names = .7)

# it kept telling me that my codon chart didn't exist 
# so I had to put it AT THE VERY TOP
# after that it finally worked

# Leucine is the most common amino acid!



## Section 6. Summary / Conclusions and Calculate Relative Synonymous Codon Usage (RSCU)

# so we know how each amino acid uses different codons
# I am going to find out which amino acid uses which codons the most!
# I also learned that RSCU is a measurement to describe codon bias apparently
# so we are going to be talking about RSCU and which codons are preferred 
# more for a certain amino acid
# RSCU = usage of a codon / usage of all codons 
# if it's = 1 it is used equally, R>1 is over represented for a certain codon
# <1 under represented like a rare codon

# let's count codons and how much each codon appears
# we get a table with codon counts

codon_counts = table(ALLcodons)

# now we are going to build a dataframe by putting each codon to 
# its amino acid
# build codon_df (Codon, Count, AA)
# table of whose names are the codons
# so the names return to names like GTC GCC etc
# this is my codon column in the data frame
# we then convert the values from table to numeric vector
# the AA looks up which amino acid each codon = to within my table
# so like GCT -> Ala, GCC -> Ala
# strings makes sure that text as regular characters and not factors 
# ^ apparently that is an issue in older versions of R and it when I was having
# some issues chat told me to add this along other code lines and it worked so I kept it

codon_df = data.frame(
  Codon = names(codon_counts),
  Count = as.numeric(codon_counts),
  AA = codon_table[names(codon_counts)],
  stringsAsFactors = FALSE
)

# calculate RSCU for each codon
# we are counting how many codons belong to each amino acid (AA)

syn_counts = table(codon_df$AA)

# total codon counts per amino acid
# this sums all codon usages for each amino acid
# we will use it for the observed/expected
total_AA_counts = tapply(codon_df$Count, codon_df$AA, sum)


# RSCU formula
# how often this specific codon appears that we are looking at
# total number of codons in that amino acid
# how many codons are synonymous (code for the same amino acod)


codon_df$RSCU = codon_df$Count/
  (total_AA_counts[codon_df$AA] / syn_counts[codon_df$AA])

# split codons into groups by amino acid

split_codons = split(codon_df, codon_df$AA)

# set up plot layout (4 rows × 6 columns) small plot

rows = 4
cols = 6
# 24 plots

par(mfrow = c(rows, cols),
mar = c(3, 3, 2, 1), #inner margin
oma = c(0, 0, 3, 0)) #outer margin

# plot RSCU barplots for each amino acid
                                     
for (aa in names(split_codons)) {
codon_data = split_codons[[aa]]

# aa vecomes each amino acid name
# codon data is that amino acid's codon, counts, and their RSCU

  
barplot(
  height = codon_data$RSCU,
  names.arg = codon_data$Codon,
  main = paste("RSCU:", aa),
  ylab = "RSCU",
  col = "pink",
  las = 2,
  ylim = c(0, max(codon_data$RSCU, 1.5)),
  cex.names = 0.6
)
  

abline(h = 1, col = "purple", lty = 2)
}

# overall title
mtext("Relative Synonymous Codon Usage (RSCU) per Amino Acid",
      outer = TRUE, cex = 1.5)


# in our results, we can see which codons get more preference in the amino acids!
# very cool
