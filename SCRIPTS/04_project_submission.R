#library(tidyverse)
#library(dplyr)
#library(ggpubr)
#library(ggrepel)
#library(data.table)
#library(VariantAnnotation)
#library(GenomicAlignments)
#library(GenomicFeatures)
#library(karyoploteR)
#library(BSgenome)
#library(BSgenome.Dmelanogaster.UCSC.dm6)
#library(sigfit)
#library(emmeans)

#making list of files
cohort_vcfs <- list.files(
  path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/ALL_VCF/ALL_VCF_UNZIP",
  pattern = ".vcf$",
  full.names = TRUE
)

#one pipeline to read in, obtain fix and info and combine
vcf_process <- function(x){
  
  #extract sample id
  sample_id <- tools::file_path_sans_ext(basename(x))
  
  #read VCF
  cohort_read_in <- readVcf(x, genome = "dm6")
  
  #FIX
  cohort_fix <- data.frame(
    CHROM = as.character(seqnames(cohort_read_in)),
    POS = start(cohort_read_in),
    ID = names(cohort_read_in),
    REF = as.character(ref(cohort_read_in)),
    ALT = sapply(alt(cohort_read_in), function(y) paste(as.character(y), collapse = ",")),
    QUAL = qual(cohort_read_in),
    FILTER = as.character(fixed(cohort_read_in)$FILTER),
    stringsAsFactors = FALSE
  )
  
  #INFO
  cohort_info <- as.data.frame(
    info(cohort_read_in)
  )
  
  #combine and add sample id
  cbind(cohort_fix, cohort_info) |>
    mutate(SAMPLE_ID = sample_id)
}

#apply function to all VCFs
cohort_fullset <- lapply(cohort_vcfs, vcf_process)

#bind resulting data frame by common rows
cohort_fullset <- bind_rows(cohort_fullset)


#renaming sample id
cohort_fullset$SAMPLE_ID <- substr(cohort_fullset$SAMPLE_ID, 9, 18)

#input columns for sex, temperature, and developmental stage
#generating a dataframe for sex, temperature, and developmental stage
#reading in the cohort metadata
cohort_metadata <- read.csv("C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/Drosophila_cohort_metadata - Sequenced.csv", na.strings = ".")
cohort_metadata <- cohort_metadata |>
  as.data.frame() |>
  #selecting desired columns from the metadata
  dplyr::select(5,10,13,14,15,21,22,23)
#renaming for uniformity and ease of joining the dataframes
colnames(cohort_metadata)[1] <- "SAMPLE_ID"
colnames(cohort_metadata)[2] <- "DEV_STAGE"
colnames(cohort_metadata)[3] <- "SEX"
colnames(cohort_metadata)[4] <- "AGE"
colnames(cohort_metadata)[5] <- "TEMP"
colnames(cohort_metadata)[6] <- "SANGER_SAMPLE_ID"


#removing unwanted samples
cohort_metadata <- cohort_metadata[!grepl("Grandarents", cohort_metadata$DEV_STAGE),]

#cleaning the DEV_STAGE column values for downstream analysis
cohort_metadata$DEV_STAGE <- ifelse(
  grepl("embryo", cohort_metadata$DEV_STAGE),
  "EMBRYO",
  cohort_metadata$DEV_STAGE
)
cohort_metadata$DEV_STAGE <- ifelse(
  grepl("1st", cohort_metadata$DEV_STAGE),
  "1ST_INSTAR",
  cohort_metadata$DEV_STAGE
)
cohort_metadata$DEV_STAGE <- ifelse(
  grepl("2st", cohort_metadata$DEV_STAGE),
  "2ST_INSTAR",
  cohort_metadata$DEV_STAGE
)
cohort_metadata$DEV_STAGE <- ifelse(
  grepl("3st", cohort_metadata$DEV_STAGE),
  "3ST_INSTAR",
  cohort_metadata$DEV_STAGE
)
cohort_metadata$DEV_STAGE <- ifelse(
  grepl("Early_Pupa", cohort_metadata$DEV_STAGE),
  "EARLY_PUPA",
  cohort_metadata$DEV_STAGE
)
cohort_metadata$DEV_STAGE <- ifelse(
  grepl("Pupa", cohort_metadata$DEV_STAGE),
  "PUPA",
  cohort_metadata$DEV_STAGE
)
cohort_metadata$DEV_STAGE <- ifelse(
  grepl("Adult", cohort_metadata$DEV_STAGE),
  "ADULT",
  cohort_metadata$DEV_STAGE
)

#adding the metadata to the fullset
setDT(cohort_fullset)
setDT(cohort_metadata)
setkey(cohort_fullset, SAMPLE_ID)
setkey(cohort_metadata, SAMPLE_ID)
cohort_fullset <- cohort_metadata[cohort_fullset]

#cleaning the age column for analysis

#remove white space
cohort_fullset <- cohort_fullset |>
  mutate(
    age_clean = tolower(trimws(AGE)))

#note whether unit is days or hours
cohort_fullset <- cohort_fullset |>
  mutate(
    unit = ifelse(str_detect(age_clean, "d"), "days", "hours"),
    #extract raw ages
    numbers = str_extract_all(age_clean, "\\d+"))

#note start and end of age range
cohort_fullset <- cohort_fullset |>
  mutate(
    start = as.numeric(sapply(numbers, '[', 1)),
    end = as.numeric(sapply(numbers, '[', 2)))

#if no end value, start = end
cohort_fullset$end[is.na(cohort_fullset$end)] <- cohort_fullset$start[is.na(cohort_fullset$end)]

#convert all day ages to hours
cohort_fullset <- cohort_fullset |>
  mutate(
    start_hours = ifelse(unit == "days",
                         start*24 + 216,
                         start),
    end_hours = ifelse(unit == "days",
                       end * 24 + 216,
                       end))

#paste new age range in hours
cohort_fullset <- cohort_fullset |>
  mutate(
    age_hours_range = paste0(start_hours, "-", end_hours)
  )

#calculate midpoint of age range
cohort_fullset$midpoint <- (cohort_fullset$start_hours + cohort_fullset$end_hours)/2  

#drop unecessary columns
cohort_fullset <- dplyr::select(cohort_fullset, -age_clean)
cohort_fullset <- cohort_fullset |>
  dplyr::select(-unit, -numbers, -start, -end, -start_hours, -end_hours)

#amending the CHROM portion of the data from NCBI chromosome names to standard Chromosome numbers
cohort_fullset <- cohort_fullset |>
  
  #adapting the CHROM column
  mutate(
    
    #any adaptation stored in the same column
    CHROM = case_when(
      
      CHROM == "NC_004354.4" ~ "chrX",
      CHROM == "NT_033779.5" ~ "chr2L",
      CHROM == "NT_033778.4" ~ "chr2R",
      CHROM == "NT_037436.4" ~ "chr3L",
      CHROM == "NT_033777.3" ~ "chr3R",
      CHROM == "NC_004353.4" ~ "chr4",
      CHROM == "NC_024512.1" ~ "chrY",
      CHROM == "NC_024511.2" ~ "chrM"
    )
  )

#analysing mutation burden data
#creating a list of files
mut_burd_tsvs <- list.files(
  path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/MUTATION_BURDEN_TSVS",
  pattern = "mut_burden.tsv$",
  full.names = TRUE
)

#function for reading in all tsvs
mut_tsv_reader <- function(x) {
  read.delim(file = x, nrows = 13, header = TRUE) |>
    as.data.frame() |>
    mutate(SAMPLE_ID = basename(x))
}

#reading in all tsvs
mutation_data_tsvs <- bind_rows(lapply(mut_burd_tsvs, mut_tsv_reader))

mutation_data_tsvs$type <- sub("\\..*$", "", rownames(mutation_data_tsvs))

#first subset to retain correct burdens, uci and lci for all samples
mutation_burden_corrected <- subset.data.frame(mutation_data_tsvs, mutation_data_tsvs$type == "corrected")

#clean sample id
mutation_burden_corrected$SAMPLE_ID <- substr(mutation_burden_corrected$SAMPLE_ID, 9, 18)

#removing NA
mutation_burden_corrected <- drop_na(mutation_burden_corrected)
#lost sample 103

#remove "corrected" column
mutation_burden_corrected <- mutation_burden_corrected[,-7]


#plotting mutation burden with confidence intervals for all samples
ggplot(data = mutation_burden_corrected,
       mapping = aes(x = SAMPLE_ID, y = burden)) +
  geom_boxplot()+
  geom_errorbar(aes(ymin = burden_lci, ymax = burden_uci), color = "red", width = 0.2)+
  theme_minimal()+
  labs(
    x = "Sample ID",
    y = "Mutation Burden",
    title = "Cohort Mutation Burden (PASSED SNVS)"
  )+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

#remove sample 103b from cohort fullset
cohort_fullset <- subset(cohort_fullset, cohort_fullset$SAMPLE_ID != "FFLYD0103b")

#incorporating mutation burden data into fullset
setDT(cohort_fullset)
setDT(mutation_burden_corrected)
setkey(cohort_fullset, SAMPLE_ID)
setkey(mutation_burden_corrected, SAMPLE_ID)
cohort_fullset <- mutation_burden_corrected[cohort_fullset]


#incorporating QC metrics
#creating a file list to read in all efficiency tsvs
effi_total_unclean <- list.files(
  #using list files function which allows specify path
  path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/QC_ANALYSIS/QC_DATA/QC_total/",
  #pattern allows importing of all files ending tsv
  pattern = "effi.tsv$",
  #fullnames enables us to retain the directory to create a source file column to keep track of data
  full.names = TRUE
)


#defining a function to clean all tsvs at once -> these tsvs are in horizontal format and directly reading them will read the column names as data instead
QC_total_cleaning <- function(files) {
  #reads in the files as a table, splitting the data into separate columns
  read.table(files, header = FALSE)|>
    #transposes the data so that the first column becomes rownames (horizontal to vertical)
    column_to_rownames(var = "V1")|>
    t()|>
    #makes the data read into a data frame
    as.data.frame()|>
    #retaining source file for sample ID
    mutate(SANGER_SAMPLE_ID = basename(files)) |>
    #retains all original data types and retains data frame format as t() changes all data types to be the same in a matrix
    type.convert(as.is = TRUE)
}

#running all files through the function
effi_total_clean <- lapply(effi_total_unclean, QC_total_cleaning)

#binding rows to make one compact df
effi_total_clean <- do.call(rbind, effi_total_clean)

#uniform sanger ids
effi_total_clean$SANGER_SAMPLE_ID <- substr(effi_total_clean$SANGER_SAMPLE_ID, 1, 16)



#repeating for meancov data
#creating a file list to read in all efficiency tsvs at once
mean_cov_total <- list.files(
  #using list files function which allows specify path, pattern
  path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/QC_ANALYSIS/QC_DATA/QC_total/",
  #pattern allows importing of all files ending tsv
  pattern = "meancov.tsv$",
  #fullnames enables us to retain the directory to create a source file column to keep track of data
  full.names = TRUE
)

#different df function as no need to clean here
mean_cov_readin <- function(apart) {
  read.delim(file = apart, nrows = 1, header = FALSE) |>
    as.data.frame()|>
    type.convert(as.is = TRUE)
}

#running all files through the function
mean_cov_total_df <- lapply(mean_cov_total, mean_cov_readin)

#binding rows to make one compact df
mean_cov_total_df <- do.call(rbind, mean_cov_total_df)

#standardise column names
colnames(mean_cov_total_df)[1] <- "SANGER_SAMPLE_ID"

colnames(mean_cov_total_df)[2] <- "MEAN_COVERAGE"

#combine all QC metrics
setDT(effi_total_clean)
setDT(mean_cov_total_df)
setkey(effi_total_clean, SANGER_SAMPLE_ID)
setkey(mean_cov_total_df, SANGER_SAMPLE_ID)
total_qc <- mean_cov_total_df[effi_total_clean]

#remove 103b from qc
total_qc <- subset(total_qc, total_qc$SANGER_SAMPLE_ID != "6416STDY14810179")

#incorporate into cohort fullset
setDT(total_qc)
setDT(cohort_fullset)
setkey(total_qc, SANGER_SAMPLE_ID)
setkey(cohort_fullset, SANGER_SAMPLE_ID)
cohort_fullset <- cohort_fullset[total_qc]

#plotting QC metrics

#function for including an upper or lower recommended threshold line
upper_threshold <- function(u, upper_label, offset, offset2) {
  list(geom_hline(
    yintercept= u,
    linetype = "dashed",
    color = "black"),
    annotate("text",
             x = offset,
             y = offset2,
             label = upper_label,
             hjust = 1.05,
             vjust = -1)
  )
}

lower_threshold <- function(l, lower_label, offset) {
  list(geom_hline(
    yintercept= l,
    linetype = "dashed",
    color = "black"),
    annotate("text",
             x = offset,
             y = l,
             label = lower_label,
             hjust = 1.05,
             vjust = -1)
  )
}


#DUPLICATE RATE
#ADD THRESHOLD VALUE -> 81%, lower 65%, upper 90%, optimal 75-76
cohort_fullset$DEV_STAGE <- factor(cohort_fullset$DEV_STAGE, levels = c("EMBRYO", "1ST_INSTAR", "2ST_INSTAR", "3ST_INSTAR", "EARLY_PUPA", "PUPA", "ADULT"))
ggplot(cohort_fullset,
       mapping = aes(x = reorder(SAMPLE_ID, DUPLICATE_RATE), y = DUPLICATE_RATE, colour = DEV_STAGE)) +
  geom_point()+
  scale_color_brewer(palette = "Set2")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "DUPLICATE RATE DATA FROM TOTAL QC", x = "SAMPLE_ID", y = "DUPLICATE RATE")+
  upper_threshold(0.81, "EFFICIENCY OPTIMUM", 20,0.81)+
  lower_threshold(0.755, "EMPIRICAL OPTIMUM", 57)




#EFFICIENCY
#ADD THRESHOLD VALUE -> 0.07
#adjusting upper threshold function so that label is seen

ggplot(cohort_fullset,
       mapping = aes(x = reorder(SAMPLE_ID, EFFICIENCY), y = EFFICIENCY, colour = DEV_STAGE)) +
  geom_point()+
  scale_color_brewer(palette = "Set2")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "EFFICIENCY DATA FROM TOTAL QC", x = "SAMPLE_ID", y = "EFFICIENCY")+
  upper_threshold(0.07, "MAXIMISED EFFICIENCY", 55, 0.069)


#GC BOTH
ggplot(cohort_fullset,
       mapping = aes(x = reorder(SAMPLE_ID, GC_BOTH), y = GC_BOTH, colour = DEV_STAGE)) +
  geom_point()+
  scale_color_brewer(palette = "Set2")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "GC_BOTH DATA FROM TOTAL QC", x = "SAMPLE_ID", y = "GC_BOTH")+
  upper_threshold(0.42, "GENOME AVERAGE", 15, 0.42)


#GC SINGLE
ggplot(cohort_fullset,
       mapping = aes(x = reorder(SAMPLE_ID, GC_SINGLE), y = GC_SINGLE, colour = DEV_STAGE)) +
  geom_point()+
  scale_color_brewer(palette = "Set2")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "GC_SINGLE DATA FROM TOTAL QC", x = "SAMPLE_ID", y = "GC_SINGLE")+
  upper_threshold(0.42, "GENOME AVERAGE", 15, 0.42)



#F-EFF
ggplot(cohort_fullset,
       mapping = aes(x = reorder(SAMPLE_ID, `F-EFF`), y = `F-EFF`, colour = DEV_STAGE)) +
  geom_point()+
  scale_color_brewer(palette = "Set2")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "DROP-OUT DATA FROM TOTAL QC", x = "SAMPLE_ID", y = "DROP-OUT")+
  upper_threshold(0.3, "UPPER THRESHOLD", 55, 0.3) +
  lower_threshold(0.1, "LOWER THRESHOLD", 55)

#MEANCOV
ggplot(cohort_fullset,
       mapping = aes(x = reorder(SAMPLE_ID, MEAN_COVERAGE), y = MEAN_COVERAGE, colour = DEV_STAGE)) +
  geom_point()+
  scale_color_brewer(palette = "Set2")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "MEAN COVERAGE DATA FROM TOTAL QC", x = "SAMPLE_ID", y = "MEAN COVERAGE")

#PLOTTED VAF AGAINST MUTATION BURDEN TO ASCERTAIN IF HIGH VAF CONTAINING SAMPLES CORRELATED WITH HIGH BURDEN, FOUND NO SIGNIFICANT SO LEFT IN
ggplot(cohort_snv_pass_check,
       mapping = aes(x = DUPLEX_VAF, y = burden)) +
  geom_point()+
  labs(title = "MUTATION BURDEN AGAINST DUPLEX VAF (PASSED SNVs)", x = "DUPLEX VAF", y = "MUTATION BURDEN")


#so total QC
#obtain fullset
#remove 96b, 103b, 105b and 147b
#filter for snvs which passed
#plot
#note overarchingly that coverage varies largely
#check if high vaf samples correlate with high mut burd?

#removing sample 96b
cohort_fullset <- subset(cohort_fullset, cohort_fullset$SAMPLE_ID != "FFLYD0096b")

#removing sample 105b
cohort_fullset <- subset(cohort_fullset, cohort_fullset$SAMPLE_ID != "FFLYD0105b")

#removing sample 147b
cohort_fullset <- subset(cohort_fullset, cohort_fullset$SAMPLE_ID != "FFLYD0147b")


#obtaining data for rate of somatic mutation
#calculate burden/age in hours
cohort_fullset <- cohort_fullset |> mutate(rate = burden/midpoint)
cohort_fullset <- cohort_fullset |> mutate(rate_year = rate/8760)
cohort_fullset <- cohort_fullset |> mutate(lwr_hour = burden_lci/midpoint)
cohort_fullset <- cohort_fullset |> mutate(lwr_year = lwr_hour/8760)
cohort_fullset <- cohort_fullset |> mutate(upr_hour = burden_uci/midpoint)
cohort_fullset <- cohort_fullset |> mutate(upr_year = upr_hour/8760)


#SUBSETTING
#generating PASSED SNVs subset
cohort_snv_pass <- subset(cohort_fullset, cohort_fullset$TYPE == "snv")
cohort_snv_pass <- subset(cohort_snv_pass, cohort_snv_pass$FILTER == "PASS")

#949 observations total: ALL PASSED SNVs, 96b, 103b, 105b, 147b out

#generating subsets for temperature
cohort_snv_pass_normal <- subset(cohort_snv_pass, cohort_snv_pass$TEMP == "normal")

#445 observations total: PASSED SNVs, NORMAL TEMPERATURE, 96b, 103b, 105b, 147b out


cohort_snv_pass_cold <- subset(cohort_snv_pass, cohort_snv_pass$TEMP == "cold")

#294 observations total: PASSED SNVs, COLD TEMPERATURE, 96b, 103b, 105b, 147b out


cohort_snv_pass_hot <- subset(cohort_snv_pass, cohort_snv_pass$TEMP == "hot")

#210 observations total: PASSED SNVs, HOT TEMPERATURE, 96b, 103b, 105b, 147b out

#generating subsets for sex - only assigned in adults
cohort_snv_pass_adults <- subset(cohort_snv_pass, cohort_snv_pass$DEV_STAGE == "ADULT")

#446 observations total: PASSED SNVs, ADULTS ONLY, 96b, 103b, 105b, 147b out

cohort_snv_pass_adult_f <- subset(cohort_snv_pass_adults, cohort_snv_pass_adults$SEX == "Female")

#298 observations total: PASSED SNVs, FEMALE ADULTS ONLY, 96b, 103b, 105b, 147b out


cohort_snv_pass_adult_m <- subset(cohort_snv_pass_adults, cohort_snv_pass_adults$SEX == "Male")

#148 observations total: PASSED SNVs, FEMALE ADULTS ONLY, 96b, 103b, 105b, 147b out

#sanity check, male + female subsets = total adult obs

#proceeding with mutation burden graphs

#mutation burden by age, normal temperature
#445 observations total: PASSED SNVs, NORMAL TEMPERATURE, 96b, 103b, 105b, 147b out

ggplot(data = cohort_snv_pass_normal,
       mapping = aes(x = midpoint, y = burden, colour = DEV_STAGE)) +
  geom_errorbar(aes(ymin = burden_lci, ymax = burden_uci), color = "black", width = 0.2)+
  geom_point()+
  labs(title = "SNV MUTATION BURDEN BY AGE (HRS), NORMAL TEMPERATURE", x = "AGE (HRS)", y = "MUTATION BURDEN")

#mutation burden by age, cold temperature
#294 observations total: PASSED SNVs, COLD TEMPERATURE, 96b, 103b, 105b, 147b out

ggplot(data = cohort_snv_pass_cold,
       mapping = aes(x = midpoint, y = burden, colour = DEV_STAGE)) +
  geom_errorbar(aes(ymin = burden_lci, ymax = burden_uci), color = "black", width = 0.2)+
  geom_point()+
  labs(title = "SNV MUTATION BURDEN BY AGE (HRS), COLD TEMPERATURE", x = "AGE (HRS)", y = "MUTATION BURDEN")

#mutation burden by age, hot temperature
#210 observations total: PASSED SNVs, HOT TEMPERATURE, 96b, 103b, 105b, 147b out

ggplot(data = cohort_snv_pass_hot,
       mapping = aes(x = midpoint, y = burden, colour = DEV_STAGE)) +
  geom_errorbar(aes(ymin = burden_lci, ymax = burden_uci), color = "black", width = 0.2)+
  geom_point()+
  labs(title = "SNV MUTATION BURDEN BY AGE (HRS), HOT TEMPERATURE", x = "AGE (HRS)", y = "MUTATION BURDEN")

#mutation burden by age, all temperatures
#949 observations total: ALL PASSED SNVs, 96b, 103b, 105b, 147b out
ggplot(data = cohort_snv_pass,
       mapping = aes(x = midpoint, y = burden)) +
  geom_errorbar(aes(ymin = burden_lci, ymax = burden_uci), color = "black", width = 0.2)+
  geom_point(data = cohort_snv_pass |> filter(TEMP == "cold"), color = "blue")+
  geom_point(data = cohort_snv_pass |> filter(TEMP == "normal"), color = "green")+
  geom_point(data = cohort_snv_pass |> filter(TEMP == "hot"), color = "red")+
  scale_x_continuous(
    n.breaks = 15
  )+
  labs(
    x = "Age (hrs)",
    y = "Mutation Burden",
    title = "Cohort Mutation Burden (PASSED SNVS, ALL TEMPERATURES)"
  )

#generating a facet wrap

ggplot(data = cohort_snv_pass,
       mapping = aes(x = midpoint, y = burden, colour = DEV_STAGE)) +
  geom_errorbar(aes(ymin = burden_lci, ymax = burden_uci), color = "black", width = 0.2)+
  geom_point()+
  facet_grid(~factor(TEMP, levels=c('cold', 'normal', 'hot')))+
  scale_x_continuous(
    n.breaks = 15
  ) +
  labs(
    x = "Age (hrs)",
    y = "Mutation Burden",
    title = "SNV MUTATION BURDEN BY AGE (hrs), ALL TEMPERATURES"
  )


#mutation rate by age, normal temperature
#445 observations total: PASSED SNVs, NORMAL TEMPERATURE, 96b, 103b, 105b, 147b out
cohort_snv_pass_normal$DEV_STAGE <- factor(cohort_snv_pass_normal$DEV_STAGE, levels = c("EMBRYO", "1ST_INSTAR", "2ST_INSTAR", "3ST_INSTAR", "EARLY_PUPA", "PUPA","ADULT"))


ggplot(data = cohort_snv_pass_normal,
       mapping = aes(x = midpoint, y = rate_year, colour = DEV_STAGE)) +
  geom_errorbar(aes(ymin = lwr_year, ymax = upr_year), color = "black", width = 0.2)+
  geom_point()+
  labs(title = "SNVs PER YEAR BY AGE, NORMAL TEMPERATURE", x = "AGE (HRS)", y = "MUTATION RATE (SNVs per year)")

#mutation rate by age, cold temperature
#294 observations total: PASSED SNVs, COLD TEMPERATURE, 96b, 103b, 105b, 147b out
cohort_snv_pass_cold$DEV_STAGE <- factor(cohort_snv_pass_cold$DEV_STAGE, levels = c("EMBRYO", "1ST_INSTAR", "2ST_INSTAR", "3ST_INSTAR", "EARLY_PUPA", "PUPA","ADULT"))

ggplot(data = cohort_snv_pass_cold,
       mapping = aes(x = midpoint, y = rate_year, colour = DEV_STAGE)) +
  geom_errorbar(aes(ymin = lwr_year, ymax = upr_year), color = "black", width = 0.2)+
  geom_point()+
  labs(title = "SNVs PER YEAR BY AGE, COLD TEMPERATURE", x = "AGE (HRS)", y = "MUTATION RATE (SNVs PER YEAR)")

#mutation rate by age, hot temperature
#210 observations total: PASSED SNVs, HOT TEMPERATURE, 96b, 103b, 105b, 147b out

ggplot(data = cohort_snv_pass_hot,
       mapping = aes(x = midpoint, y = rate_year, colour = DEV_STAGE)) +
  geom_point()+
  labs(title = "SNV MUTATION BURDEN BY AGE (HRS), HOT TEMPERATURE", x = "AGE (HRS)", y = "MUTATION BURDEN")

#mutation rate by age, all temperatures
#949 observations total: ALL PASSED SNVs, 96b, 103b, 105b, 147b out
ggplot(data = cohort_snv_pass,
       mapping = aes(x = midpoint, y = rate_year)) +
  geom_point(data = cohort_snv_pass |> filter(TEMP == "cold"), color = "blue")+
  geom_point(data = cohort_snv_pass |> filter(TEMP == "normal"), color = "green")+
  geom_point(data = cohort_snv_pass |> filter(TEMP == "hot"), color = "red")+
  scale_x_continuous(
    n.breaks = 15
  )+
  labs(
    x = "Age (hrs)",
    y = "Mutation Burden",
    title = "Cohort Mutation Burden (PASSED SNVS, ALL TEMPERATURES)"
  )

#generate a facet wrap

ggplot(data = cohort_snv_pass,
       mapping = aes(x = midpoint, y = rate_year, colour = DEV_STAGE)) +
  geom_point()+
  facet_grid(~factor(TEMP, levels=c('cold', 'normal', 'hot')))+
  scale_x_continuous(
    n.breaks = 15
  ) +
  labs(
    x = "Age (hrs)",
    y = "Mutation Burden",
    title = "SNV MUTATION BURDEN BY AGE (hrs), ALL TEMPERATURES"
  )



#446 observations total: PASSED SNVs, ADULTS ONLY, 96b, 103b, 105b, 147b out

cohort_snv_pass_adults$TEMP <- factor(cohort_snv_pass_adults$TEMP, levels = c("hot", "normal", "cold"))

ggplot(data = cohort_snv_pass_adults,
       mapping = aes(x = midpoint, y = burden, colour = TEMP)) +
  geom_errorbar(aes(ymin = burden_lci, ymax = burden_uci), color = "black", width = 0.2)+
  geom_point()+
  facet_wrap(facets = vars(SEX))+
  labs(
    x = "Age (hrs)",
    y = "Mutation Burden",
    title = "ADULT SNV MUTATION BURDEN BY AGE (hrs), COLOUR BY TEMP, FACET BY SEX"
  )

#linear regression accounting for temperature covariate
lin_reg_sex_temp <- lm(burden ~ SEX + midpoint + TEMP, data = cohort_snv_pass_adults)

#assess difference between male and female
anova(lin_reg_sex_temp)
#significant: 152.8084 < 2.2e-16
#pairwise

emmeans(lin_reg_sex_temp, pairwise ~ SEX, adjust = "bonferroni")
#female to male: 
#estimated difference of adjusted means:5.19e-09
#adjusted p: <0.0001

#significance testing pre vs post hatch (mostly equal groups)
#splitting between larval, pupal, adult

#normal only initially, let's see
#if interesting, we'll put all together and colour for temp

#dev stage long first
cohort_snv_pass_normal <- cohort_snv_pass_normal |>
  
  #adapting the DEV_STAGE column
  mutate(
    
    #any adaptation stored in the same column
    DEV_STAGE_LONG = case_when(
      
      DEV_STAGE == "EMBRYO" ~ "PRE-PUPAL",
      DEV_STAGE == "1ST_INSTAR" ~ "PRE-PUPAL",
      DEV_STAGE == "2ST_INSTAR" ~ "PRE-PUPAL",
      DEV_STAGE == "3ST_INSTAR" ~ "PRE-PUPAL",
      DEV_STAGE == "EARLY_PUPA" ~ "PUPAL",
      DEV_STAGE == "PUPA" ~ "PUPAL",
      DEV_STAGE == "ADULT" ~ "ADULT"
    )
  )

cohort_snv_pass_normal$DEV_STAGE_LONG <- factor(cohort_snv_pass_normal$DEV_STAGE_LONG, levels = c("PRE-PUPAL", "PUPAL", "ADULT"))


#adjust for multiple testing#
cohort_snv_pass_normal_prepup <- cohort_snv_pass_normal |> filter(DEV_STAGE_LONG == "PRE-PUPAL") 
cor_test_normal_prepup <- cor.test(cohort_snv_pass_normal_prepup$midpoint, cohort_snv_pass_normal_prepup$burden, method = "pearson", conf.level = 0.95)
p.adjust(cor_test_normal_prepup$p.value, method = "bonferroni")
print(cor_test_normal_prepup$estimate)

#adjusted p value: 1.914258e-85
#r value: 0.9516227 

cohort_snv_pass_normal_pup <- cohort_snv_pass_normal |> filter(DEV_STAGE_LONG == "PUPAL") 
cor_test_normal_pup <- cor.test(cohort_snv_pass_normal_pup$midpoint, cohort_snv_pass_normal_pup$burden, method = "pearson", conf.level = 0.95)
p.adjust(cor_test_normal_pup$p.value, method = "bonferroni")
print(cor_test_normal_pup$estimate)

#adjusted p value: 9.462844e-14
#r value: -0.7484617 

cohort_snv_pass_normal_adult <- cohort_snv_pass_normal |> filter(DEV_STAGE_LONG == "ADULT") 
cor_test_normal_adult <- cor.test(cohort_snv_pass_normal_adult$midpoint, cohort_snv_pass_normal_adult$burden, method = "pearson", conf.level = 0.95)
p.adjust(cor_test_normal_adult$p.value, method = "bonferroni")
print(cor_test_normal_adult$estimate)

#adjusted p value: 3.198561e-25
#r value: 0.6361954 

cohort_snv_pass_normal$DEV_STAGE <- factor(cohort_snv_pass_normal$DEV_STAGE, levels = c("EMBRYO", "2ST_INSTAR", "3ST_INSTAR", "EARLY_PUPA", "PUPA", "ADULT"))

ggplot(data = cohort_snv_pass_normal,
       mapping = aes(x = midpoint, y = burden, colour = DEV_STAGE)) +
  geom_errorbar(aes(ymin = burden_lci, ymax = burden_uci), color = "black", width = 0.2)+
  geom_point()+
  facet_wrap(facets = vars(DEV_STAGE_LONG), scales = "free_x")+
  labs(title = "SNV MUTATION BURDEN BY AGE (HRS), NORMAL TEMPERATURE", x = "AGE (HRS)", y = "MUTATION BURDEN")


#kruskal wallis test
kruskal.test(burden ~ DEV_STAGE_LONG, data = cohort_snv_pass_normal)
#p value 0.0004721, proceed with pairwise wilcox test
pairwise.wilcox.test(cohort_snv_pass_normal$burden, cohort_snv_pass_normal$DEV_STAGE_LONG,
                     p.adjust.method = "BH")
#significant difference between pre-pupal and pupal: adjusted p is 0.1981
#significant difference between pupal and adult: 0.0016
#diff prepupal adult: 0.0026    

#add r and p values by textbox

#again with all temps

#dev stage long first
cohort_snv_pass <- cohort_snv_pass |>
  
  #adapting the DEV_STAGE column
  mutate(
    
    #any adaptation stored in the same column
    DEV_STAGE_LONG = case_when(
      
      DEV_STAGE == "EMBRYO" ~ "PRE-PUPAL",
      DEV_STAGE == "1ST_INSTAR" ~ "PRE-PUPAL",
      DEV_STAGE == "2ST_INSTAR" ~ "PRE-PUPAL",
      DEV_STAGE == "3ST_INSTAR" ~ "PRE-PUPAL",
      DEV_STAGE == "EARLY_PUPA" ~ "PUPAL",
      DEV_STAGE == "PUPA" ~ "PUPAL",
      DEV_STAGE == "ADULT" ~ "ADULT"
    )
  )

cohort_snv_pass$DEV_STAGE_LONG <- factor(cohort_snv_pass$DEV_STAGE_LONG, levels = c("PRE-PUPAL", "PUPAL", "ADULT"))


#adjust for multiple testing
#account for temperature as a covariate, linear regression
cohort_snv_pass_prepup$midpoint <- as.numeric(cohort_snv_pass_prepup$midpoint)
cohort_snv_pass_prepup <- cohort_snv_pass |> filter(DEV_STAGE_LONG == "PRE-PUPAL") 
cor_test_pass_prepup <- lm(burden ~ midpoint + TEMP, data = cohort_snv_pass_prepup)


#p value: <2e-16 #accounted for temp
#slope: 2.508e-10 

cohort_snv_pass_pup <- cohort_snv_pass |> filter(DEV_STAGE_LONG == "PUPAL") 
cor_test_pass_pup <- lm(burden ~ midpoint + TEMP, data = cohort_snv_pass_pup)
summary(cor_test_pass_pup)

#p value: <2e-16
#slope: -9.525e-11 

cohort_snv_pass_adult <- cohort_snv_pass |> filter(DEV_STAGE_LONG == "ADULT") 
cor_test_pass_adult <- lm(burden ~ midpoint + TEMP, data = cohort_snv_pass_adult)
summary(cor_test_pass_adult)

#p value: < 2e-16
#slope: 4.628e-12

#applying multiple testing
p_vals_pass <- c(2e-16, 2e-16, 2e-16)
p_vals_adj <- p.adjust(p_vals_pass, method = "bonferroni")
p_vals_adj

#adjusted: 6e-16

#assessing difference across stages accounting to temperature and age
stage_model <- lm(burden ~ DEV_STAGE_LONG + TEMP + midpoint, data = cohort_snv_pass)
anova(stage_model)
#significant: 68.2098 <2e-16
#pairwise

emmeans(stage_model, pairwise ~ DEV_STAGE_LONG, adjust = "bonferroni")
#significant difference pre-pup to pup
#estimate difference of adjusted means:4.82e-09
#adjusted p value:<0.0001

#significant difference pre-pup to adult
#estimate difference of adjusted means:1.11e-08
#adjusted p value:<0.0001

#significant difference pup to adult
#estimate difference of adjusted means:6.29e-09
#adjusted p value: <0.0001



cohort_snv_pass$DEV_STAGE <- factor(cohort_snv_pass$DEV_STAGE, levels = c("EMBRYO", "1ST_INSTAR", "2ST_INSTAR", "3ST_INSTAR", "EARLY_PUPA", "PUPA", "ADULT"))
cohort_snv_pass$TEMP <- factor(cohort_snv_pass$TEMP, levels = c("hot", "normal", "cold"))

ggplot(data = cohort_snv_pass,
       mapping = aes(x = midpoint, y = burden, colour = TEMP)) +
  geom_errorbar(aes(ymin = burden_lci, ymax = burden_uci), color = "black", width = 0.2)+
  geom_point()+
  facet_wrap(facets = vars(DEV_STAGE_LONG), scales = "free_x")+
  labs(title = "SNV MUTATION BURDEN BY AGE (HRS), ALL TEMPERATURES", x = "AGE (HRS)", y = "MUTATION BURDEN")

#add r and p values by textbox

#all temps genomic distribution
#now we first generate GRanges for plotting
cohort_snv_pass_gr <- GRanges(
  seqnames = cohort_snv_pass$CHROM,
  ranges = IRanges(start = cohort_snv_pass$POS, end = cohort_snv_pass$POS),
  sample = cohort_snv_pass$SAMPLE_ID,
  stage = cohort_snv_pass$DEV_STAGE,
  variant_type = cohort_snv_pass$TYPE,
  coverage = cohort_snv_pass$DUPLEX_COV,
  vaf = cohort_snv_pass$DUPLEX_VAF,
  temp = cohort_snv_pass$TEMP
)

#blank karyploter graph
kp_all <- plotKaryotype(genome = "dm6", plot.type = 4)

stage_levels <- unique(cohort_snv_pass_gr$stage)
cols <- setNames(rainbow(length(stage_levels)), stage_levels)
vaf_cols <- ifelse(cohort_snv_pass_gr$vaf > 0.01, "red", "black")
temp_levels <- unique(cohort_snv_pass_gr$temp)
shapes <- c(16,17,15)
names(shapes) <- temp_levels
vaf_scaled <- cohort_snv_pass_gr$vaf / max(cohort_snv_pass_gr$vaf, na.rm = TRUE)


legend(
  "topright",
  legend = stage_levels,
  col = cols,
  pch = 16,
  title = "Dev Stage",
  cex = 1
)
legend(
  "topright",
  legend = c("Low VAF", "Medium VAF", "High VAF"),
  pt.cex = c(0.05,0.5,1),
  pch = 16,
  title = "VAF (*10)",
  cex = 1,
  inset = c(0, 0.8))

legend(
  "topright",
  legend = temp_levels,
  col = "black",
  pch = shapes,
  title = "Temperature",
  cex = 1
)
kpAxis(kp_all,
       ymin = 0,
       ymax = max(cohort_snv_pass_gr$coverage)
       )
kpAddLabels(
  kp_all,
  labels = "Coverage (Scaled)",
  side = "left",
  cex = 3,
  srt = 90)
kpAddMainTitle(kp_all,
               "Genomic Distribution of All Passed SNVs",
               cex = 1)
kpPoints(
  kp_all,
  data = cohort_snv_pass_gr,
  y = cohort_snv_pass_gr$coverage,
  ymin = 0,
  ymax = max(cohort_snv_pass_gr$coverage),
  cex = 0.6,
  col = vaf_cols,
  pch = shapes
)

#another example with cov < 100 removed
#generating subset
cohort_snv_pass_cov <- subset(cohort_snv_pass, cohort_snv_pass$DUPLEX_COV > 100)

#now we first generate GRanges for plotting
cohort_snv_pass_cov_gr <- GRanges(
  seqnames = cohort_snv_pass_cov$CHROM,
  ranges = IRanges(start = cohort_snv_pass_cov$POS, end = cohort_snv_pass_cov$POS),
  sample = cohort_snv_pass_cov$SAMPLE_ID,
  stage = cohort_snv_pass_cov$DEV_STAGE,
  variant_type = cohort_snv_pass_cov$TYPE,
  coverage = cohort_snv_pass_cov$DUPLEX_COV,
  vaf = cohort_snv_pass_cov$DUPLEX_VAF,
  temp = cohort_snv_pass_cov$TEMP
)
min(cohort_snv_pass_cov_gr$coverage)

kp_all_cov <- plotKaryotype(genome = "dm6", plot.type = 4)

stage_levels <- unique(cohort_snv_pass_cov_gr$stage)
cols <- setNames(rainbow(length(stage_levels)), stage_levels)
temp_levels <- unique(cohort_snv_pass_cov_gr$temp)
shapes <- c(16,17,15)
names(shapes) <- temp_levels


legend(
  "topright",
  legend = stage_levels,
  col = cols,
  pch = 16,
  title = "Dev Stage",
  cex = 0.6
)
legend(
  "topright",
  legend = c("Low VAF", "Medium VAF", "High VAF"),
  pt.cex = c(0.05,0.5,1),
  pch = 16,
  title = "% VAF",
  cex = 0.6,
  inset = c(0, 0.3))

legend(
  "topright",
  legend = temp_levels,
  col = "black",
  pch = shapes,
  title = "Temperature",
  cex = 0.6,
  inset = c(0, 0.5)
)
kpAxis(kp_all,
       ymin = 100,
       ymax = 920
       )
kpAddLabels(
  kp_all,
  labels = "Coverage",
  side = "left",
  cex = 3,
  srt = 90)
kpAddMainTitle(kp_all,
               "Genomic Distribution of All Passed SNVs (Coverage > 100)",
               cex = 1)
kpPoints(
  kp_all,
  data = cohort_snv_pass_cov_gr,
  y = cohort_snv_pass_cov_gr$coverage,
  ymin = 100,
  ymax = 920,
  cex = cohort_snv_pass_cov_gr$vaf*100,
  col = cols,
  pch = shapes
)


#the following code is adapted from source courtesy of Martin Santamarina Garcia

#all temperatures
# ---- Read reference genome ----

dnaREF <- BSgenome.Dmelanogaster.UCSC.dm6

cohort_tri <- cohort_snv_pass
cohort_tri <- dplyr::select(cohort_tri,6,7,14,15,17,18,22)
# ---- Extract upstream & downstream sequence context from reference
cohort_tri$contextREF <- apply(cohort_tri, 1, function(x) get_context(as.character(x["CHROM"]), as.integer(x["POS"]), dnaREF))
cohort_tri$contextREF_5   = substr(cohort_tri$contextREF, 1, 1)
cohort_tri$contextREF_pos = substr(cohort_tri$contextREF, 2, 2)
cohort_tri$contextREF_3   = substr(cohort_tri$contextREF, 3, 3)

cohort_tri$contextALT<-paste0(cohort_tri$contextREF_5,cohort_tri$ALT,cohort_tri$contextREF_3)


### --- Normalize to pyrimidine context ----------------------------------
cohort_tri$contexREF_normalised<-ifelse(cohort_tri$REF %in% c("A","G"), reverseComplement(DNAStringSet(cohort_tri$contextREF)), cohort_tri$contextREF)
cohort_tri$contextALT_normalised<-ifelse(cohort_tri$REF %in% c("A","G"), reverseComplement(DNAStringSet(cohort_tri$contextALT)), cohort_tri$contextALT)

cohort_tri$mutation_normalised<-paste0(cohort_tri$contexREF_normalised,">",cohort_tri$contextALT_normalised)
table(cohort_tri$mutation_normalised)


### Set mutations object for plotting
data("cosmic_signatures_v2")
mutations <- rep(0, length(colnames(cosmic_signatures_v2))) # 96
names(mutations) <- colnames(cosmic_signatures_v2)
mutations

### Add mutation counts
tbl <- as.integer(table(cohort_tri$mutation_normalised))
names(tbl) <- names(table(cohort_tri$mutation_normalised))

mutations[names(tbl)] <- tbl


### Plot raw spectra
plot_spectrum(mutations, name = "SNV Spectra - Full Cohort", pdf_path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/SNV_SPECTRA_FINAL.pdf")

#seeing high oxidative signatures, doing by dev stage
cohort_tri_adult <- subset(cohort_tri, cohort_tri$DEV_STAGE == "ADULT")

### Set mutations object for plotting
data("cosmic_signatures_v2")
mutations_adult <- rep(0, length(colnames(cosmic_signatures_v2))) # 96
names(mutations_adult) <- colnames(cosmic_signatures_v2)
mutations_adult

### Add mutation counts
tbl_adult <- as.integer(table(cohort_tri_adult$mutation_normalised))
names(tbl_adult) <- names(table(cohort_tri_adult$mutation_normalised))

mutations_adult[names(tbl_adult)] <- tbl_adult


### Plot raw spectra
plot_spectrum(mutations_adult, name = "SNV Spectra - Adult", pdf_path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/SNV_SPECTRA_ADULT_FINAL.pdf")

#embryo+larvae
cohort_tri_prepup <- subset(cohort_tri, cohort_tri$DEV_STAGE != "ADULT")
cohort_tri_prepup <- subset(cohort_tri_prepup, cohort_tri_prepup$DEV_STAGE != "PUPA")
cohort_tri_prepup <- subset(cohort_tri_prepup, cohort_tri_prepup$DEV_STAGE != "EARLY_PUPA")

### Set mutations object for plotting
data("cosmic_signatures_v2")
mutations_prepup <- rep(0, length(colnames(cosmic_signatures_v2))) # 96
names(mutations_prepup) <- colnames(cosmic_signatures_v2)
mutations_prepup

### Add mutation counts
tbl_prepup <- as.integer(table(cohort_tri_prepup$mutation_normalised))
names(tbl_prepup) <- names(table(cohort_tri_prepup$mutation_normalised))

mutations_prepup[names(tbl_prepup)] <- tbl_prepup


### Plot raw spectra
plot_spectrum(mutations_prepup, name = "SNV Spectra - Embryo + Larval", pdf_path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/SNV_SPECTRA_PRE_PUPAL_FINAL.pdf")


#pupa
cohort_tri_pupa <- subset(cohort_tri, cohort_tri$DEV_STAGE == "PUPA" | cohort_tri$DEV_STAGE == "EARLY_PUPA")

### Set mutations object for plotting
data("cosmic_signatures_v2")
mutations_pupa <- rep(0, length(colnames(cosmic_signatures_v2))) # 96
names(mutations_pupa) <- colnames(cosmic_signatures_v2)
mutations_pupa

### Add mutation counts
tbl_pupa <- as.integer(table(cohort_tri_pupa$mutation_normalised))
names(tbl_pupa) <- names(table(cohort_tri_pupa$mutation_normalised))

mutations_pupa[names(tbl_pupa)] <- tbl_pupa


### Plot raw spectra
plot_spectrum(mutations_pupa, name = "SNV Spectra - Pupa", pdf_path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/SNV_SPECTRA_PUPA_FINAL.pdf")

#comparing mutation spectra between temperature cohorts
#using normal

# ---- Read reference genome ----

dnaREF <- BSgenome.Dmelanogaster.UCSC.dm6

cohort_tri_normal <- cohort_snv_pass_normal
cohort_tri_normal <- dplyr::select(cohort_tri_normal,1,9,14,15,17,18,21)
# ---- Extract upstream & downstream sequence context from reference
cohort_tri_normal$contextREF <- apply(cohort_tri_normal, 1, function(x) get_context(as.character(x["CHROM"]), as.integer(x["POS"]), dnaREF))
cohort_tri_normal$contextREF_5   = substr(cohort_tri_normal$contextREF, 1, 1)
cohort_tri_normal$contextREF_pos = substr(cohort_tri_normal$contextREF, 2, 2)
cohort_tri_normal$contextREF_3   = substr(cohort_tri_normal$contextREF, 3, 3)

cohort_tri_normal$contextALT<-paste0(cohort_tri_normal$contextREF_5,cohort_tri_normal$ALT,cohort_tri_normal$contextREF_3)


### --- Normalize to pyrimidine context ----------------------------------
cohort_tri_normal$contexREF_normalised<-ifelse(cohort_tri_normal$REF %in% c("A","G"), reverseComplement(DNAStringSet(cohort_tri_normal$contextREF)), cohort_tri_normal$contextREF)
cohort_tri_normal$contextALT_normalised<-ifelse(cohort_tri_normal$REF %in% c("A","G"), reverseComplement(DNAStringSet(cohort_tri_normal$contextALT)), cohort_tri_normal$contextALT)

cohort_tri_normal$mutation_normalised<-paste0(cohort_tri_normal$contexREF_normalised,">",cohort_tri_normal$contextALT_normalised)
table(cohort_tri_normal$mutation_normalised)


### Set mutations object for plotting
data("cosmic_signatures_v2")
mutations_normal <- rep(0, length(colnames(cosmic_signatures_v2))) # 96
names(mutations_normal) <- colnames(cosmic_signatures_v2)
mutations_normal

### Add mutation counts
tbl_normal <- as.integer(table(cohort_tri_normal$mutation_normalised))
names(tbl_normal) <- names(table(cohort_tri_normal$mutation_normalised))

mutations_normal[names(tbl_normal)] <- tbl_normal


### Plot raw spectra
plot_spectrum(mutations_normal, name = "SNV Spectra - Normal Temperature", pdf_path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/SNV_SPECTRA_NORMAL.pdf")

#Now for cold temperature

cohort_tri_cold <- cohort_snv_pass_cold
cohort_tri_cold <- dplyr::select(cohort_tri_cold,1,9,14,15,17,18,21)
# ---- Extract upstream & downstream sequence context from reference
cohort_tri_cold$contextREF <- apply(cohort_tri_cold, 1, function(x) get_context(as.character(x["CHROM"]), as.integer(x["POS"]), dnaREF))
cohort_tri_cold$contextREF_5   = substr(cohort_tri_cold$contextREF, 1, 1)
cohort_tri_cold$contextREF_pos = substr(cohort_tri_cold$contextREF, 2, 2)
cohort_tri_cold$contextREF_3   = substr(cohort_tri_cold$contextREF, 3, 3)

cohort_tri_cold$contextALT<-paste0(cohort_tri_cold$contextREF_5,cohort_tri_cold$ALT,cohort_tri_cold$contextREF_3)


### --- Normalize to pyrimidine context ----------------------------------
cohort_tri_cold$contexREF_normalised<-ifelse(cohort_tri_cold$REF %in% c("A","G"), reverseComplement(DNAStringSet(cohort_tri_cold$contextREF)), cohort_tri_cold$contextREF)
cohort_tri_cold$contextALT_normalised<-ifelse(cohort_tri_cold$REF %in% c("A","G"), reverseComplement(DNAStringSet(cohort_tri_cold$contextALT)), cohort_tri_cold$contextALT)

cohort_tri_cold$mutation_normalised<-paste0(cohort_tri_cold$contexREF_normalised,">",cohort_tri_cold$contextALT_normalised)
table(cohort_tri_cold$mutation_normalised)


### Set mutations object for plotting
data("cosmic_signatures_v2")
mutations_cold <- rep(0, length(colnames(cosmic_signatures_v2))) # 96
names(mutations_cold) <- colnames(cosmic_signatures_v2)
mutations_cold

### Add mutation counts
tbl_cold <- as.integer(table(cohort_tri_cold$mutation_normalised))
names(tbl_cold) <- names(table(cohort_tri_cold$mutation_normalised))

mutations_cold[names(tbl_cold)] <- tbl_cold


### Plot raw spectra
plot_spectrum(mutations_cold, name = "SNV Spectra - Cold Temperature", pdf_path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/SNV_SPECTRA_COLD.pdf")


#hot temperature
cohort_tri_hot <- cohort_snv_pass_hot
cohort_tri_hot <- dplyr::select(cohort_tri_hot,1,9,14,15,17,18,21)
# ---- Extract upstream & downstream sequence context from reference
cohort_tri_hot$contextREF <- apply(cohort_tri_hot, 1, function(x) get_context(as.character(x["CHROM"]), as.integer(x["POS"]), dnaREF))
cohort_tri_hot$contextREF_5   = substr(cohort_tri_hot$contextREF, 1, 1)
cohort_tri_hot$contextREF_pos = substr(cohort_tri_hot$contextREF, 2, 2)
cohort_tri_hot$contextREF_3   = substr(cohort_tri_hot$contextREF, 3, 3)

cohort_tri_hot$contextALT<-paste0(cohort_tri_hot$contextREF_5,cohort_tri_hot$ALT,cohort_tri_hot$contextREF_3)


### --- Normalize to pyrimidine context ----------------------------------
cohort_tri_hot$contexREF_normalised<-ifelse(cohort_tri_hot$REF %in% c("A","G"), reverseComplement(DNAStringSet(cohort_tri_hot$contextREF)), cohort_tri_hot$contextREF)
cohort_tri_hot$contextALT_normalised<-ifelse(cohort_tri_hot$REF %in% c("A","G"), reverseComplement(DNAStringSet(cohort_tri_hot$contextALT)), cohort_tri_hot$contextALT)

cohort_tri_hot$mutation_normalised<-paste0(cohort_tri_hot$contexREF_normalised,">",cohort_tri_hot$contextALT_normalised)
table(cohort_tri_hot$mutation_normalised)


### Set mutations object for plotting
data("cosmic_signatures_v2")
mutations_hot <- rep(0, length(colnames(cosmic_signatures_v2))) # 96
names(mutations_hot) <- colnames(cosmic_signatures_v2)
mutations_hot

### Add mutation counts
tbl_hot <- as.integer(table(cohort_tri_hot$mutation_normalised))
names(tbl_hot) <- names(table(cohort_tri_hot$mutation_normalised))

mutations_hot[names(tbl_hot)] <- tbl_hot


### Plot raw spectra
plot_spectrum(mutations_hot, name = "SNV Spectra - Hot Temperature", pdf_path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/SNV_SPECTRA_HOT.pdf")

######## Martin's code ends here #######

#comparing dupcaller and nanoseq outputs
#call in dupcaller

#making list of files
dupcaller_vcfs <- list.files(
  path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/dupcaller",
  pattern = ".vcf$",
  full.names = TRUE
)

#one pipeline to read in, obtain fix and info and combine
dup_vcf_process <- function(x){
  
  #extract sample id
  dup_sample_id <- tools::file_path_sans_ext(basename(x))
  
  #read VCF
  dup_read_in <- readVcf(x, genome = "dm6")
  
  #FIX
  dup_fix <- data.frame(
    CHROM = as.character(seqnames(dup_read_in)),
    POS = start(dup_read_in),
    ID = names(dup_read_in),
    REF = as.character(ref(dup_read_in)),
    ALT = sapply(alt(dup_read_in), function(y) paste(as.character(y), collapse = ",")),
    QUAL = qual(dup_read_in),
    FILTER = as.character(fixed(dup_read_in)$FILTER),
    stringsAsFactors = FALSE
  )
  
  #INFO
  dup_info <- as.data.frame(
    info(dup_read_in)
  )
  
  #GT
  dup_geno <- lapply(geno(dup_read_in), as.data.frame)
  full_geno <- cbind(dup_geno$AC, dup_geno$RC, dup_geno$DP)
  colnames(full_geno) <- c("AC_T", "AC_N", "RC_T", "RC_N", "DP_T", "DP_N")
  
  
  #combine and add sample id
  cbind(dup_fix, dup_info, full_geno) |>
    mutate(SAMPLE_ID = dup_sample_id)
}

dup_fullset <- lapply(dupcaller_vcfs, dup_vcf_process)
dup_fullset <- bind_rows(dup_fullset)

dup_fullset <- dplyr::select(dup_fullset, "CHROM", "POS", "ID", "REF", "ALT", "FILTER", "TN", "AC_T", "DP_T", "SAMPLE_ID")
#EXTRACT AP AND DP FOR DUPLEX COV AND VAF CALCULATIONS!!
#ap/dp is vaf
#dp is cov
dup_fullset <- dup_fullset |> mutate(DUPLEX_VAF = AC_T/DP_T)
colnames(dup_fullset)[9] <- "DUPLEX_COV"
colnames(dup_fullset)[8] <- "NV"
colnames(dup_fullset)[7] <- "TRI"


#PUT ONTO GENOMIC DIST GRAPH THINGY

dup_fullset_pass <- subset(dup_fullset, dup_fullset$FILTER == "PASS")

#renaming sample id so it's easier to manage
dup_fullset_pass$SAMPLE_ID <- substr(dup_fullset_pass$SAMPLE_ID,11, 20)

#add in metadata
#input columns for sex, temperature, and developmental stage
#generating a dataframe for sex, temperature, and developmental stage
#reading in the cohort metadata

dup_metadata <- read.csv("C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/Drosophila_cohort_metadata - Sequenced.csv", na.strings = ".")

dup_metadata <- dup_metadata |>
  as.data.frame() |>
  #selecting desired columns from the metadata
  dplyr::select(Donor, Developmental.Stage..Original., Sex..Original., Temperature, Age, SangerSampleID)

#renaming for uniformity and ease of merging the dataframes
colnames(dup_metadata)[1] <- "SAMPLE_ID"
colnames(dup_metadata)[2] <- "DEV_STAGE"
colnames(dup_metadata)[3] <- "SEX"
colnames(dup_metadata)[4] <- "TEMP"
colnames(dup_metadata)[5] <- "AGE"
colnames(dup_metadata)[6] <- "SANGER_SAMPLE_ID"

#only retaining wanted samples
dup_metadata <- subset(dup_metadata, dup_metadata$SAMPLE_ID == "FFLYD0136b"|dup_metadata$SAMPLE_ID == "FFLYD0166b"| dup_metadata$SAMPLE_ID == "FFLYD0219b")

#cleaning the DEV_STAGE column values for downstream analysis
dup_metadata$DEV_STAGE <- ifelse(
  grepl("embryo", dup_metadata$DEV_STAGE),
  "EMBRYO",
  dup_metadata$DEV_STAGE
)
dup_metadata$DEV_STAGE <- ifelse(
  grepl("1st", dup_metadata$DEV_STAGE),
  "1ST_INSTAR",
  dup_metadata$DEV_STAGE
)
dup_metadata$DEV_STAGE <- ifelse(
  grepl("2st", dup_metadata$DEV_STAGE),
  "2ST_INSTAR",
  dup_metadata$DEV_STAGE
)
dup_metadata$DEV_STAGE <- ifelse(
  grepl("3st", dup_metadata$DEV_STAGE),
  "3ST_INSTAR",
  dup_metadata$DEV_STAGE
)
dup_metadata$DEV_STAGE <- ifelse(
  grepl("Early_Pupa", dup_metadata$DEV_STAGE),
  "EARLY_PUPA",
  dup_metadata$DEV_STAGE
)
dup_metadata$DEV_STAGE <- ifelse(
  grepl("Pupa", dup_metadata$DEV_STAGE),
  "PUPA",
  dup_metadata$DEV_STAGE
)
dup_metadata$DEV_STAGE <- ifelse(
  grepl("Adult", dup_metadata$DEV_STAGE),
  "ADULT",
  dup_metadata$DEV_STAGE
)

#adding the metadata to the fullset
setDT(dup_fullset_pass)
setDT(dup_metadata)
setkey(dup_fullset_pass, SAMPLE_ID)
setkey(dup_metadata, SAMPLE_ID)
dup_fullset_pass <- dup_metadata[dup_fullset_pass]


#cleaning the age column for analysis
dup_fullset_pass <- dup_fullset_pass |>
  mutate(
    age_clean = tolower(trimws(AGE)))

dup_fullset_pass <- dup_fullset_pass |>
  mutate(
    unit = ifelse(str_detect(age_clean, "d"), "days", "hours"),
    numbers = str_extract_all(age_clean, "\\d+"))

dup_fullset_pass <- dup_fullset_pass |>
  mutate(
    start = as.numeric(sapply(numbers, '[', 1)),
    end = as.numeric(sapply(numbers, '[', 2)))

dup_fullset_pass$end[is.na(dup_fullset_pass$end)] <- dup_fullset_pass$start[is.na(dup_fullset_pass$end)]


dup_fullset_pass <- dup_fullset_pass |>
  mutate(
    start_hours = ifelse(unit == "days",
                         start*24 + 216,
                         start),
    end_hours = ifelse(unit == "days",
                       end * 24 + 216,
                       end))

dup_fullset_pass <- dup_fullset_pass |>
  mutate(
    age_hours_range = paste0(start_hours, "-", end_hours)
  )

dup_fullset_pass$midpoint <- (dup_fullset_pass$start_hours + dup_fullset_pass$end_hours)/2  
dup_fullset_pass <- dplyr::select(dup_fullset_pass, -age_clean)
dup_fullset_pass <- dup_fullset_pass |>
  dplyr::select(-unit, -numbers, -start, -end, -start_hours, -end_hours)


#amending the CHROM portion of the data from NCBI chromosome names to standard Chromosome numbers
dup_fullset_pass <- dup_fullset_pass |>
  
  #adapting the CHROM column
  mutate(
    
    #any adaptation stored in the same column
    CHROM = case_when(
      
      CHROM == "NC_004354.4" ~ "chrX",
      CHROM == "NT_033779.5" ~ "chr2L",
      CHROM == "NT_033778.4" ~ "chr2R",
      CHROM == "NT_037436.4" ~ "chr3L",
      CHROM == "NT_033777.3" ~ "chr3R",
      CHROM == "NC_004353.4" ~ "chr4",
      CHROM == "NC_024512.1" ~ "chrY",
      CHROM == "NC_024511.2" ~ "chrM"
    )
  )

#comparing against our nanoseq variants for these samples in all analyses
#let's do a genomic distribution graph
#and then a trinucleotide spectra?
#start with that
#for genomic distribution let's structure such that common is one colour, dup only is one colour, nano only is one colour
#generate an equivalent nano df
nano_fullset_pass <- subset(cohort_snv_pass, cohort_snv_pass$SAMPLE_ID == "FFLYD0136b"|cohort_snv_pass$SAMPLE_ID == "FFLYD0166b"| cohort_snv_pass$SAMPLE_ID == "FFLYD0219b")
nano_fullset_pass <- dplyr::select(nano_fullset_pass, 6,7,8,10,9,11,14,15,16,17,18,20,22,23,25,24)
dup_fullset_pass <- dplyr::select(dup_fullset_pass, -17, -18)

#make separate GRanges objects for dup only, shared, and nano only
#look at dup only vs nano only genomic dist with cov and vaf
#look at common variants genomic dist with dup cov and then nano cov

dup_only <- anti_join(dup_fullset_pass, nano_fullset_pass, by = c("SAMPLE_ID", "CHROM", "POS", "ID", "REF", "ALT")) 
nano_only <- anti_join(nano_fullset_pass, dup_fullset_pass, by = c("SAMPLE_ID", "CHROM", "POS", "ID", "REF", "ALT"))
dup_nano_common <- inner_join(dup_fullset_pass, nano_fullset_pass, by = c("SAMPLE_ID", "CHROM", "POS", "ID", "REF", "ALT"))
#works!
nano_only$DEV_STAGE <- as.character(nano_only$DEV_STAGE)

dup_only_gr <- GRanges(
  seqnames = dup_only$CHROM,
  ranges = IRanges(start = dup_only$POS, end = dup_only$POS),
  sample = dup_only$SAMPLE_ID,
  vaf = dup_only$DUPLEX_VAF,
  coverage = dup_only$DUPLEX_COV,
  stage = dup_only$DEV_STAGE
)

nano_only_gr <- GRanges(
  seqnames = nano_only$CHROM,
  ranges = IRanges(start = nano_only$POS, end = nano_only$POS),
  sample = nano_only$SAMPLE_ID,
  vaf = nano_only$DUPLEX_VAF,
  coverage = nano_only$DUPLEX_COV,
  stage = nano_only$DEV_STAGE
)

dup_nano_common_gr <- GRanges(
  seqnames = dup_nano_common$CHROM,
  ranges = IRanges(start = dup_nano_common$POS, end = dup_nano_common$POS),
  sample = dup_nano_common$SAMPLE_ID,
  vaf_dup = dup_nano_common$DUPLEX_VAF.x,
  vaf_nano = dup_nano_common$DUPLEX_VAF.y,
  cov_dup = dup_nano_common$DUPLEX_COV.x,
  cov_nano = dup_nano_common$DUPLEX_COV.y,
  stage = dup_nano_common$DEV_STAGE.x
)



dup_nano_compare <- plotKaryotype(genome = "dm6", plot.type = 4)
sample_levels <- unique(dup_only_gr$sample)
shapes <- c(16,17,15)  
names(shapes) <- sample_levels
pch_vals_dup <- shapes[dup_only_gr$sample]
pch_vals_nano <- shapes[nano_only_gr$sample]
pch_vals_both <- shapes[dup_nano_common_gr$sample]

dist_levels <- c("DUPCALLER ONLY", "COMMON", "NANOSEQ ONLY")
dist_cols <- c("red", "green", "blue")  
names(dist_cols) <- dist_levels


legend(
  "topright",
  legend = sample_levels,
  col = "black",
  pch = shapes,
  title = "SAMPLES",
  cex = 0.6
)

legend(
  "topright",
  legend = dist_levels,
  col = dist_cols,
  title = "VARIANT DISTRIBUTION",
  pch = 16,
  cex = 0.6,
  inset = c(0,0.3)
)

kpAddMainTitle(kp_dup_nano_compare,
               "Genomic Distribution of Dupcaller vs Nanoseq Variant Calls",
               cex = 1)
kpPoints(
  karyoplot = kp_dup_nano_compare,
  data = dup_only_gr,
  y = 0.3,
  pch = pch_vals_dup,
  cex = 0.8,
  col = "red"
)

kpPoints(
  karyoplot = kp_dup_nano_compare,
  data = nano_only_gr,
  y = 0.1,
  pch = pch_vals_nano,
  cex = 0.8,
  col = "blue"
)
kpPoints(
  karyoplot = kp_dup_nano_compare,
  data = dup_nano_common_gr,
  y = 0.2,
  pch = pch_vals_both,
  cex = 0.8,
  col = "green"
)




#look at common variants genomic dist with dup cov and then nano cov


#COMPARE DUP COV VS NANO COV FOR COMMON VARIANTS
#plotting
#y = cov
#cex = vaf
#cols = dev stage
#pch = samples

common_variants <- plotKaryotype(genome = "dm6", plot.type = 4, cex = 3)

sample_levels <- unique(dup_nano_common_gr$sample)
shapes <- c(16,17,15)  
names(shapes) <- sample_levels
pch_vals_common <- shapes[dup_nano_common_gr$sample]

legend(
  "topright",
  legend = sample_levels,
  col = "black",
  pch = shapes,
  title = "SAMPLES",
  cex = 0.8
)


legend(
  "topright",
  legend = c("Low VAF", "Medium VAF", "High VAF"),
  pt.cex = c(0.05,0.5,1),
  pch = 16,
  title = "%VAF (Dup)",
  cex = 0.8,
  inset = c(0, 0.4))

kpAxis(
  karyoplot = common_variants,
  ymin = min(dup_nano_common_gr$cov_nano),
  ymax = max(dup_nano_common_gr$cov_dup),
  cex = 3
)

kpAddMainTitle(common_variants,
               "Genomic Distribution of Dupcaller and Nanoseq Common Variant Calls (NanoSeq Coverage)",
               cex = 3)

#first with nano coverage
kpPoints(
  karyoplot = common_variants,
  data = dup_nano_common_gr,
  y = dup_nano_common_gr$cov_nano,
  ymax = max(dup_nano_common_gr$cov_dup),
  ymin = min(dup_nano_common_gr$cov_nano),
  pch = pch_vals_common,
  cex = dup_nano_common_gr$vaf_nano*100,
  col = "blue"
)

#then with dup coverage
kpPoints(
  karyoplot = common_variants,
  data = dup_nano_common_gr,
  y = dup_nano_common_gr$cov_dup,
  ymax = max(dup_nano_common_gr$cov_dup),
  ymin = min(dup_nano_common_gr$cov_nano),
  pch = pch_vals_common,
  cex = dup_nano_common_gr$vaf_dup*100,
  col = "red"
)
min(dup_nano_common_gr$cov_dup)
max(dup_nano_common_gr$cov_dup)
min(dup_nano_common_gr$cov_nano)
max(dup_nano_common_gr$cov_nano)


#COMPARING TRI SPECTRA

#using Martin's code:
# ---- Define Functions ---- 
get_context <- function(chr, pos, ref_genome) {
  if (!(chr %in% names(ref_genome))) return(NA)
  # Extract 3bp region
  seq <- subseq(ref_genome[[chr]], pos-1, pos+1)
  as.character(seq)
}

# ---- Read reference genome ----

dnaREF <- BSgenome.Dmelanogaster.UCSC.dm6

duplex_tri <- dup_fullset_pass
duplex_tri <- dplyr::select(duplex_tri,1,7,8,10,11,13)
# ---- Extract upstream & downstream sequence context from reference
duplex_tri$contextREF <- apply(duplex_tri, 1, function(x) get_context(as.character(x["CHROM"]), as.integer(x["POS"]), dnaREF))
duplex_tri$contextREF_5   = substr(duplex_tri$contextREF, 1, 1)
duplex_tri$contextREF_pos = substr(duplex_tri$contextREF, 2, 2)
duplex_tri$contextREF_3   = substr(duplex_tri$contextREF, 3, 3)

duplex_tri$contextALT<-paste0(duplex_tri$contextREF_5,duplex_tri$ALT,duplex_tri$contextREF_3)


### --- Normalize to pyrimidine context ----------------------------------
duplex_tri$contexREF_normalised<-ifelse(duplex_tri$REF %in% c("A","G"), reverseComplement(DNAStringSet(duplex_tri$contextREF)), duplex_tri$contextREF)
duplex_tri$contextALT_normalised<-ifelse(duplex_tri$REF %in% c("A","G"), reverseComplement(DNAStringSet(duplex_tri$contextALT)), duplex_tri$contextALT)

duplex_tri$mutation_normalised<-paste0(duplex_tri$contexREF_normalised,">",duplex_tri$contextALT_normalised)
table(duplex_tri$mutation_normalised)


### Set mutations object for plotting
data("cosmic_signatures_v2")
mutations_dup <- rep(0, length(colnames(cosmic_signatures_v2))) # 96
names(mutations_dup) <- colnames(cosmic_signatures_v2)
mutations_dup

### Add mutation counts
tbl_dup <- as.integer(table(duplex_tri$mutation_normalised))
names(tbl_dup) <- names(table(duplex_tri$mutation_normalised))

mutations_dup[names(tbl_dup)] <- tbl_dup


### Plot raw spectra
plot_spectrum(mutations_dup, name = "SNV Spectra - Dupcaller", pdf_path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/SNV_SPECTRA_DUPCALLER.pdf")

#plotting for nano same

#using Martin's code:
# ---- Define Functions ---- 
get_context <- function(chr, pos, ref_genome) {
  if (!(chr %in% names(ref_genome))) return(NA)
  # Extract 3bp region
  seq <- subseq(ref_genome[[chr]], pos-1, pos+1)
  as.character(seq)
}

# ---- Read reference genome ----

dnaREF <- BSgenome.Dmelanogaster.UCSC.dm6

nano_tri <- nano_fullset_pass
nano_tri <- dplyr::select(nano_tri,1,7,8,10,11,13)
# ---- Extract upstream & downstream sequence context from reference
nano_tri$contextREF <- apply(nano_tri, 1, function(x) get_context(as.character(x["CHROM"]), as.integer(x["POS"]), dnaREF))
nano_tri$contextREF_5   = substr(nano_tri$contextREF, 1, 1)
nano_tri$contextREF_pos = substr(nano_tri$contextREF, 2, 2)
nano_tri$contextREF_3   = substr(nano_tri$contextREF, 3, 3)

nano_tri$contextALT<-paste0(nano_tri$contextREF_5,nano_tri$ALT,nano_tri$contextREF_3)


### --- Normalize to pyrimidine context ----------------------------------
nano_tri$contexREF_normalised<-ifelse(nano_tri$REF %in% c("A","G"), reverseComplement(DNAStringSet(nano_tri$contextREF)), nano_tri$contextREF)
nano_tri$contextALT_normalised<-ifelse(nano_tri$REF %in% c("A","G"), reverseComplement(DNAStringSet(nano_tri$contextALT)), nano_tri$contextALT)

nano_tri$mutation_normalised<-paste0(nano_tri$contexREF_normalised,">",nano_tri$contextALT_normalised)
table(nano_tri$mutation_normalised)


### Set mutations object for plotting
data("cosmic_signatures_v2")
mutations_nano <- rep(0, length(colnames(cosmic_signatures_v2))) # 96
names(mutations_nano) <- colnames(cosmic_signatures_v2)
mutations_nano

### Add mutation counts
tbl_nano <- as.integer(table(nano_tri$mutation_normalised))
names(tbl_nano) <- names(table(nano_tri$mutation_normalised))

mutations_nano[names(tbl_nano)] <- tbl_nano


### Plot raw spectra
plot_spectrum(mutations_nano, name = "SNV Spectra - NanoSeq", pdf_path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/SNV_SPECTRA_NANOSEQ.pdf")

#plotting dup only

dup_only_tri <- dup_only_kp
dup_only_tri <- dplyr::select(dup_only_tri,1,7,8,10,11)
# ---- Extract upstream & downstream sequence context from reference
dup_only_tri$contextREF <- apply(dup_only_tri, 1, function(x) get_context(as.character(x["CHROM"]), as.integer(x["POS"]), dnaREF))
dup_only_tri$contextREF_5   = substr(dup_only_tri$contextREF, 1, 1)
dup_only_tri$contextREF_pos = substr(dup_only_tri$contextREF, 2, 2)
dup_only_tri$contextREF_3   = substr(dup_only_tri$contextREF, 3, 3)

dup_only_tri$contextALT<-paste0(dup_only_tri$contextREF_5,dup_only_tri$ALT,dup_only_tri$contextREF_3)


### --- Normalize to pyrimidine context ----------------------------------
dup_only_tri$contexREF_normalised<-ifelse(dup_only_tri$REF %in% c("A","G"), reverseComplement(DNAStringSet(dup_only_tri$contextREF)), dup_only_tri$contextREF)
dup_only_tri$contextALT_normalised<-ifelse(dup_only_tri$REF %in% c("A","G"), reverseComplement(DNAStringSet(dup_only_tri$contextALT)), dup_only_tri$contextALT)

dup_only_tri$mutation_normalised<-paste0(dup_only_tri$contexREF_normalised,">",dup_only_tri$contextALT_normalised)
table(dup_only_tri$mutation_normalised)


### Set mutations object for plotting
data("cosmic_signatures_v2")
mutations_dup_only <- rep(0, length(colnames(cosmic_signatures_v2))) # 96
names(mutations_dup_only) <- colnames(cosmic_signatures_v2)
mutations_dup_only

### Add mutation counts
tbl_dup_only <- as.integer(table(dup_only_tri$mutation_normalised))
names(tbl_dup_only) <- names(table(dup_only_tri$mutation_normalised))

mutations_dup_only[names(tbl_dup_only)] <- tbl_dup_only


### Plot raw spectra
plot_spectrum(mutations_dup_only, name = "SNV Spectra - DupCaller Only", pdf_path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/SNV_SPECTRA_DUP_ONLY.pdf")


#plotting nano only
nano_only_tri <- nano_only_kp
nano_only_tri <- dplyr::select(nano_only_tri,1,7,8,10,11)
# ---- Extract upstream & downstream sequence context from reference
nano_only_tri$contextREF <- apply(nano_only_tri, 1, function(x) get_context(as.character(x["CHROM"]), as.integer(x["POS"]), dnaREF))
nano_only_tri$contextREF_5   = substr(nano_only_tri$contextREF, 1, 1)
nano_only_tri$contextREF_pos = substr(nano_only_tri$contextREF, 2, 2)
nano_only_tri$contextREF_3   = substr(nano_only_tri$contextREF, 3, 3)

nano_only_tri$contextALT<-paste0(nano_only_tri$contextREF_5,nano_only_tri$ALT,nano_only_tri$contextREF_3)


### --- Normalize to pyrimidine context ----------------------------------
nano_only_tri$contexREF_normalised<-ifelse(nano_only_tri$REF %in% c("A","G"), reverseComplement(DNAStringSet(nano_only_tri$contextREF)), nano_only_tri$contextREF)
nano_only_tri$contextALT_normalised<-ifelse(nano_only_tri$REF %in% c("A","G"), reverseComplement(DNAStringSet(nano_only_tri$contextALT)), nano_only_tri$contextALT)

nano_only_tri$mutation_normalised<-paste0(nano_only_tri$contexREF_normalised,">",nano_only_tri$contextALT_normalised)
table(nano_only_tri$mutation_normalised)


### Set mutations object for plotting
data("cosmic_signatures_v2")
mutations_nano_only <- rep(0, length(colnames(cosmic_signatures_v2))) # 96
names(mutations_nano_only) <- colnames(cosmic_signatures_v2)
mutations_nano_only

### Add mutation counts
tbl_nano_only <- as.integer(table(nano_only_tri$mutation_normalised))
names(tbl_nano_only) <- names(table(nano_only_tri$mutation_normalised))

mutations_nano_only[names(tbl_nano_only)] <- tbl_nano_only


### Plot raw spectra
plot_spectrum(mutations_nano_only, name = "SNV Spectra - NanoSeq Only", pdf_path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/SNV_SPECTRA_NANO_ONLY.pdf")



#PLOTTING COMMON ONLY

dup_nano_common_tri <- dup_nano_common_kp
dup_nano_common_tri <- dplyr::select(dup_nano_common_tri,1,7,8,10,11)
# ---- Extract upstream & downstream sequence context from reference
dup_nano_common_tri$contextREF <- apply(dup_nano_common_tri, 1, function(x) get_context(as.character(x["CHROM"]), as.integer(x["POS"]), dnaREF))
dup_nano_common_tri$contextREF_5   = substr(dup_nano_common_tri$contextREF, 1, 1)
dup_nano_common_tri$contextREF_pos = substr(dup_nano_common_tri$contextREF, 2, 2)
dup_nano_common_tri$contextREF_3   = substr(dup_nano_common_tri$contextREF, 3, 3)

dup_nano_common_tri$contextALT<-paste0(dup_nano_common_tri$contextREF_5,dup_nano_common_tri$ALT,dup_nano_common_tri$contextREF_3)


### --- Normalize to pyrimidine context ----------------------------------
dup_nano_common_tri$contexREF_normalised<-ifelse(dup_nano_common_tri$REF %in% c("A","G"), reverseComplement(DNAStringSet(dup_nano_common_tri$contextREF)), dup_nano_common_tri$contextREF)
dup_nano_common_tri$contextALT_normalised<-ifelse(dup_nano_common_tri$REF %in% c("A","G"), reverseComplement(DNAStringSet(dup_nano_common_tri$contextALT)), dup_nano_common_tri$contextALT)

dup_nano_common_tri$mutation_normalised<-paste0(dup_nano_common_tri$contexREF_normalised,">",dup_nano_common_tri$contextALT_normalised)
table(dup_nano_common_tri$mutation_normalised)


### Set mutations object for plotting
data("cosmic_signatures_v2")
mutations_common_only <- rep(0, length(colnames(cosmic_signatures_v2))) # 96
names(mutations_common_only) <- colnames(cosmic_signatures_v2)
mutations_common_only

### Add mutation counts
tbl_common_only <- as.integer(table(dup_nano_common_tri$mutation_normalised))
names(tbl_common_only) <- names(table(dup_nano_common_tri$mutation_normalised))

mutations_common_only[names(tbl_common_only)] <- tbl_common_only


### Plot raw spectra
plot_spectrum(mutations_common_only, name = "SNV Spectra - DupCaller + NanoSeq Common", pdf_path = "C:/Users/Arya Kalavath/Documents/3rd Year/PROJECT/PROJECT_CODING/VCF_ANALYSIS_FINAL/VCF_DATA/SNV_SPECTRA_NANO_DUP_COMMON.pdf")

