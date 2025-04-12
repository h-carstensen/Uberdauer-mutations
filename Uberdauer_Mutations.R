##### Analyzing csu60 enhancer WGS data #####
### Edited 2025-03-25 ### 

##### Load Libraries #####
library(tidyverse)
library(here)

##### Read in data #####
WGS_SNP <- read_csv(here("WGS_uberdauer_all_SNP.csv")) #Need to change location on other computers

WGS_indel_1 <- read_csv(here("WGS_uberdauer_all_indels_1.csv"))
WGS_indel_2 <- read_csv(here("WGS_uberdauer_all_indels_2.csv")) #Indels file was too large for GitHub - split into two

WGS_indel <- full_join(WGS_indel_1, WGS_indel_2) #Joining the two indel files

##### How to see unique values ######

unique(WGS_version_1$ExonicFunc.refGene)

#### How to count number of distinct values

Data_frame %>%
  summarise(count = n_distinct(Column))

##### Filtering SNP sheet 

WGS_SNP_filtered <- WGS_SNP %>% 
  filter(Func.refGene != "intergenic", Func.refGene != "downstream", Func.refGene != "upstream",
         Func.refGene != "upstream;downstream", Func.refGene != "UTR3", Func.refGene != "UTR5",
         Func.refGene != "UTR5;UTR3") %>% #Filtering out all these values. 
  mutate(Gene = Gene.refGene) %>%  #Adding a column called Gene so the sheet retains a column with all gene names
  pivot_wider(names_from = Genome, values_from = Gene.refGene) %>% #Pivoting wider to make columns for genotypes, lining up identical mutations
  filter(!complete.cases(csu60)) %>%  
  filter(!complete.cases(csu60_2023)) %>% #Filtering out mutations that occur in both versions of csu60
  unite(col = "Ref_Alt", c(Ref,Alt), sep = "_", remove = FALSE) %>%  #Combine the reference and alternate columns
  filter(Ref_Alt == "G_A" | Ref_Alt == "C_T" | Ref_Alt == "T_C" | Ref_Alt == "A_G") %>% #Filter for just EMS mutations
  arrange(Chr,Start) %>% #Arrange by chromosome and start site, to line up multiple rows of the same gene
  filter(Func.refGene == "exonic" | Func.refGene == "splicing" | Func.refGene == "exonic;splicing" ) %>% #Filtering to just these values, removing intronic
  filter(ExonicFunc.refGene == "nonsynonymous SNV" | ExonicFunc.refGene == "." | ExonicFunc.refGene == "stopgain" | 
           ExonicFunc.refGene == "unknown" | ExonicFunc.refGene == "stoploss")  #Filtering to just these values, removing synonymous SNV


##### Filtering indel sheet

WGS_indel_filtered <- WGS_indel %>%
  filter(ExonicFunc.refGene != ".") %>%  #Filtering out exonic function "." leaving ones that affect coding regions
  unite(col = "Start_end", c(Genome,Start,End), sep = "_", remove = FALSE) %>% #Creating a column with mutant_start_end, which can be used to remove duplicate rows in the data
  distinct(Start_end, .keep_all = TRUE) %>%  #Removing duplicate rows in the data
  select(!(Start_end)) %>% #Now removing the Start_end column, because it messes up pivoting wider in the next step
  mutate(Gene = Gene.refGene) %>%  #Adding a column called Gene so the sheet retains a column with all gene names
  pivot_wider(names_from = Genome, values_from = Gene.refGene) %>% #pivoting wider to create a column with each mutant, lining up the same exact mutation in a row
  unite(col = "Ref_Alt", c(Ref,Alt), sep = "_", remove = FALSE) %>%  #Making this column that exists in the SNP sheet so the sheets can be easily joined later
  filter(!complete.cases(csu60)) %>%  
  filter(!complete.cases(csu60_2023))  #Filtering out mutations that occur in both versions of csu60

#### Joining the SNP and indel dataframes, arranging by chromosomal position

WGS_indel_SNP <- full_join(WGS_indel_filtered, WGS_SNP_filtered) %>%
  arrange(Chr,Start)

WGS_indel_filtered %>% 
  write_csv(here("List_1.csv"))


###### Filtering to genes with two or more rows (two or more unique mutations)

WGS_Indel_SNP_2ormore <- WGS_indel_SNP %>%
  group_by(Gene) %>%
  filter(n()>1)  #### works to filter to two or more

WGS_Indel_SNP_2ormore %>% 
  write_csv(here("List_2.csv"))


##### Filtering to genes with three or more rows (three or more unique mutations)

WGS_Indel_SNP_3ormore <- WGS_indel_SNP %>%
  group_by(Gene) %>%
  filter(n()>2)

WGS_Indel_SNP_3ormore %>% 
  write_csv(here("List_3.csv"))
