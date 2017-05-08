#
# Generate a PX file for submission to proteomeXchange
#

library(tidyverse)
library(readxl)
library(stringr)

# Submitter name
out <- tibble(col1='MTD', col2='submitter_name', col3='Benjamin R. Jack', col4='', col5='', col6='', col7='')
# Submitter email
out <- bind_rows(out, c(col1='MTD', col2='submitter_email', col3='benjamin.r.jack@gmail.com', col4='', col5='', col6='', col7=''))
# Submitter affiliation
out <- bind_rows(out, c(col1='MTD', col2='submitter_affiliation', col3='The University of Texas at Austin', col4='', col5='', col6='', col7=''))
# Pride login name
out <- bind_rows(out, c(col1='MTD', col2='submitter_pride_login', col3='b.jack@utexas.edu', col4='', col5='', col6='', col7=''))
# Lab head name, email, and affiliation
out <- bind_rows(out, c(col1='MTD', col2='lab_head_name', col3='Claus O. Wilke', col4='', col5='', col6='', col7=''))
out <- bind_rows(out, c(col1='MTD', col2='lab_head_email', col3='wilke@austin.utexas.edu', col4='', col5='', col6='', col7=''))
out <- bind_rows(out, c(col1='MTD', col2='lab_head_affiliation', col3='The University of Texas at Austin', col4='', col5='', col6='', col7=''))
# Project title and description
out <- bind_rows(out, c(col1='MTD', col2='project_title', col3='Reduced protein expression in a virus attenuated by codon deoptimization', col4='', col5='', col6='', col7=''))
out <- bind_rows(out, c(col1='MTD', col2='project_description', col3='A general means of viral attenuation involves the extensive recoding of synonymous codons in the viral genome. The mechanistic underpinnings of this approach remain unclear, however.  Using quantitative proteomics and RNA sequencing, we explore the molecular basis of attenuation in a strain of bacteriophage T7 whose major capsid gene was engineered to carry 182 suboptimal codons. We do not detect transcriptional effects from recoding. Proteomic observations reveal that translation is halved for the recoded major capsid gene, and a more modest reduction applies to several co-expressed downstream genes. We observe no changes in protein abundances of other co-expressed genes that are encoded upstream. Viral burst size, like capsid protein abundance, is also decreased by half. Together, these observations suggest that, in this virus, reduced translation of an essential polycistronic transcript and diminished virion assembly form the molecular basis of attenuation.', col4='', col5='', col6='', col7=''))
out <- bind_rows(out, c(col1='MTD', col2='sample_processing_protocol', col3='E. coli was grown in LB broth to a concentration of 108 cells/mL at 37C with agitation, then infected with phage at an MOI of 2.5. At 1, 5, and 9 minutes post-infection, 2 mL of bacterial suspension were
                        removed from the phage-infected cultures and pelleted in a microcentrifuge. Pellets were either flash frozen
                        in liquid nitrogen or immediately used for downstream processes. T7-infected E. coli cell pellets were resuspended in 50 mM Tris-HCl
                        pH 8.0, 10 mM DTT. 2,2,2-trifluoroethanol (Sigma) was added to 50% (v/v) final concentration and samples
                        were incubated at 56C for 45 minutes. Following incubation, iodoacetamide was added to a concentration
                        of 25 mM and samples were incubated at room temperature in the dark for 30 minutes. Samples were diluted
                        10-fold with 2 mM CaCl2, 50 mM Tris-HCl, pH 8.0. Samples were digested with trypsin (Pierce) at 37C for
                        5 hours. Digestion was quenched by adding formic acid to 1% (v/v). Tryptic peptides were bound, washed,
                        and eluted from HyperSep C18 SpinTips (Thermo Scientific). Eluted peptides were dried by speed-vac and
                        resuspended in Buffer C (5% acetonitrile, 0.1% formic acid) for analysis by LC-MS/MS.
                        For LC-MS/MS analysis, peptides were subjected to separation by C18 reverse phase chromatography on a
                        Dionex Ultimate 3000 RSLCnano UHPLC system (Thermo Scientific). Peptides were loaded onto an Acclaim
                        C18 PepMap RSLC column (Dionex; Thermo Scientific) and eluted using a 5–40% acetonitrile gradient
                        over 250 minutes at 300 nL/min flow rate. Eluted peptides were directly injected into an Orbitrap Elite
                        mass spectrometer (Thermo Scientific) by nano-electrospray and subject to data-dependent tandem mass
                        spectrometry, with full precursor ion scans (MS1) collected at 60,0000 resolution. Monoisotopic precursor
                        selection and charge-state screening were enabled, with ions of charge > +1 selected for collision-induced
                        dissociation (CID). Up to 20 fragmentation scans (MS2) were collected per MS1. Dynamic exclusion was
                        active with 45 s exclusion for ions selected twice within a 30 s window.', col4='', col5='', col6='', col7=''))
out <- bind_rows(out, c(col1='MTD', col2='data_processing_protocol', col3='Peptides were subjected to separation by C18 reverse phase chromatography on a
Dionex Ultimate 3000 RSLCnano UHPLC system (Thermo Scientific). Peptides were loaded onto an Acclaim
                        C18 PepMap RSLC column (Dionex; Thermo Scientific) and eluted using a 5–40% acetonitrile gradient
                        over 250 minutes at 300 nL/min flow rate. Eluted peptides were directly injected into an Orbitrap Elite
                        mass spectrometer (Thermo Scientific) by nano-electrospray and subject to data-dependent tandem mass
                        spectrometry, with full precursor ion scans (MS1) collected at 60,0000 resolution. Monoisotopic precursor
                        selection and charge-state screening were enabled, with ions of charge > +1 selected for collision-induced
                        dissociation (CID). Up to 20 fragmentation scans (MS2) were collected per MS1. Dynamic exclusion was
                        active with 45 s exclusion for ions selected twice within a 30 s window.
                        We assigned each peptide to a protein or protein group (in the case of ambiguous peptides which map to
                        multiple proteins) using Proteome Discoverer (Thermo Scientific) and REL606 and T7 reference proteomes
                        (NCBI: NC_012967, NC_001604.1) concatenated with a database of contaminant proteins (http://www.
                        biochem.mpg.de/5111795/maxquant).', col4='', col5='', col6='', col7=''))

# Keywords
out <- bind_rows(out, c(col1='MTD', col2='keywords', col3='Codon deoptimization, viral attenuation, codon usage', col4='', col5='', col6='', col7=''))
# Submission type
out <- bind_rows(out, c(col1='MTD', col2='submission_type', col3='PARTIAL', col4='', col5='', col6='', col7=''))
# Experiment type
out <- bind_rows(out, c(col1='MTD', col2='experiment_type', col3='[PRIDE, PRIDE:0000429, Shotgun proteomics, ]', col4='', col5='', col6='', col7=''))
# Species
out <- bind_rows(out, c(col1='MTD', col2='species', col3='[NEWT, 83333, Escherichia coli K-12, ]', col4='', col5='', col6='', col7=''))
out <- bind_rows(out, c(col1='MTD', col2='species', col3='[NEWT, 10760, Enterobacteria phage T7, ]', col4='', col5='', col6='', col7=''))
# Tissue
out <- bind_rows(out, c(col1='MTD', col2='tissue', col3='[PRIDE, PRIDE: 0000442, Tissue not applicable to dataset,]', col4='', col5='', col6='', col7=''))
# Instrument
out <- bind_rows(out, c(col1='MTD', col2='instrument', col3='[MS, MS:1001910, LTQ Orbitrap Elite,]', col4='', col5='', col6='', col7=''))
out <- bind_rows(out, c(col1='MTD', col2='modification', col3='[PRIDE, PRIDE:0000398, No PTMs are included in the dataset, ]', col4='', col5='', col6='', col7=''))

#
# File mappings and experimental factors
# 

# Insert a blank row
out <- bind_rows(out, c(col1='', col2='', col3='', col4='', col5='', col6='', col7=''))
samples <- read_excel("sample_list.xlsx") %>% 
  filter(`Data Type` == "mass-spec", `Biological replicate` != 1) %>% 
  mutate(file_type = "raw")
samples_msf <- mutate(samples, File = str_replace(File, ".raw", ".msf"), Directory = str_replace(Directory, "raw", "msf_1_4")) %>%
  mutate(file_type="search")
samples <- bind_rows(samples, samples_msf) %>%
  mutate(technical_replicate = str_match(Sample, "PA[0-9]+([a-z])")[,2]) %>%
  mutate(file_id = row_number(), file_path = paste0(Directory, "/", File)) %>%
  group_by(Sample) %>%
  mutate(file_mapping = lst(file_id)) %>%
  rowwise() %>%
  mutate(file_mapping = paste0(file_mapping[file_mapping != file_id], collapse = ",")) %>%
  mutate(fme = "FME", sme = "SME") %>%
  mutate(Strain = str_replace(Strain, "11-46", "wildtype") %>% str_replace("11-44", "recoded") %>% str_replace("11-42", "evolved")) %>%
  mutate(experimental_factor = paste0(Strain, " strain, ", `Time point`, " min after infection, biological replicate ", `Biological replicate`, ", technical replicate ", technical_replicate))

mappings_head <- tibble(col1='FMH', col2='file_id', col3='file_type', col4='file_path', col5='file_mapping', col6='', col7='')
mappings <- tibble(col1=samples$fme, col2=as.character(samples$file_id), col3=samples$file_type, col4=samples$'file_path', col5=samples$'file_mapping', col6='', col7='')

out <- bind_rows(out, mappings_head, mappings)

# Insert a blank row
out <- bind_rows(out, c(col1='', col2='', col3='', col4='', col5='', col6='', col7=''))
annot_head <- tibble(col1='SMH', col2='file_id', col3='species', col4='tissue', col5='instrument', col6='modification', col7='experimental_factor')
annots <- tibble(col1=samples$sme, col2=as.character(samples$file_id), col3='[NEWT, 10760, Enterobacteria phage T7, ]', col4='[PRIDE, PRIDE: 0000442, Tissue not applicable to dataset,]', col5='[MS, MS:1001910, LTQ Orbitrap Elite,]', col6='[PRIDE, PRIDE:0000398, No PTMs are included in the dataset, ]', col7=samples$experimental_factor)
out <- bind_rows(out, annot_head, annots)

write_tsv(out, "summary.px")