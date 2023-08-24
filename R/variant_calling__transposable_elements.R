library(data.table)
library(dplyr)

#sample <- "test"

if (length(commandArgs(trailingOnly = TRUE)) == 0) {
  stop("Please provide a sample name as a command-line argument.")
}

sample <- commandArgs(trailingOnly = TRUE)[1]

df <-
  paste0(
    "data/variantCalling/transposableElements/",
    sample, ".vcf"
  ) %>%
  fread() %>%
  as.data.frame() %>%
  rename(
    ID = UUID,
    '#CHROM' = Chrom,
    FILTER = Filter,
    POS = Start
  ) %>%
  mutate(
    REF = "N",
    QUAL = ".",
    FORMAT = ".",
    ALT = ifelse(Inversion, "<INV>", "<INS>"),
    INFO = paste0(
      "SVTYPE=",
      ifelse(Inversion, "INV", "INS"),
      ";SVLEN=",
      LengthIns,
      ";END=",
      End,
      ";SVCALLERS=TLDR"
    )
  ) %>%
  filter(FILTER == "PASS") %>%
  relocate('#CHROM')


fwrite(
  df,
  paste0(
    "data/variantCalling/transposableElements/",
    sample, ".vcf"
  ), sep = "\t"
)