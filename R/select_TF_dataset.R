require(tidyverse)
require(stringr)

rawdf <- read_tsv("https://www.encodeproject.org/metadata/type=Experiment&assay_slims=DNA+binding&assembly=hg19&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&target.investigated_as=transcription+factor&biosample_type=immortalized+cell+line&limit=all/metadata.tsv")

expDF <- read_tsv("https://www.encodeproject.org/report.tsv?type=Experiment&assay_title=ChIP-seq&assembly=hg19",
                  skip=1)


df <- rawdf %>% 
  rename(
    cell = `Biosample term name`,
    TF = `Experiment target`,
    format = `File format`
    ) %>% 
  mutate(TF = str_replace(TF, "-human", "")) %>% 
  left_join(expDF, by = c("Experiment accession" = "Accession")) %>% 
  select(TF, cell, format, Assembly, Description, matches("Lab"), `File download URL`, matches("treatment"), everything())


df %>%
  filter(Assembly == "hg19") %>% 
  filter(TF == "NR3C1") %>%
  filter(format == "bam") %>%
  filter(cell == "A549") %>% 
  filter(`Output type` == "alignments") %>% 
  filter(Lab.x == "ENCODE Processing Pipeline") %>% 
  View()


# testing 
id = "ENCSR000BJT"

rawdf %>% filter(`Experiment accession` == "ENCSR000BJT") %>% 
  select(`Experiment target`)

expDF %>% filter(`Accession` == "ENCSR000BJT") %>% 
  select(`Target label`, Description, everything())

