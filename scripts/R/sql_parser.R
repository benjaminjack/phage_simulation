###########################################################################
## This script processes Thermo formatted SQLite database files into 
## normalized protein abundances.
##
## Author: Benjamin R. Jack
## Email: benjamin.r.jack@gmail.com
## December 2015
###########################################################################

# Remove any variables in memory
rm(list = ls())

# Load appropriate libraries
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# Read in Thermo SQLite database file
my_db <- src_sqlite("11-46_1min_1a_20150515-2.msf")

PeptidesProteins <- tbl(my_db, "PeptidesProteins") %>%
  select(PeptideID, ProteinID)

Peptides <- tbl(my_db, "Peptides") %>%
  select(PeptideID, SpectrumID, ConfidenceLevel, Sequence) %>%
  filter(ConfidenceLevel == 3)

CustomFields <- tbl(my_db, "CustomDataFields")

CustomPeptides <- tbl(my_db, sql("SELECT FieldID, PeptideID, CAST(FieldValue as REAL) AS FieldValue FROM CustomDataPeptides"))

PEP <- left_join(CustomPeptides, CustomFields) %>%
  select(PeptideID, DisplayName, FieldValue) %>%
  collect() %>%
  spread(DisplayName, FieldValue)

ProteinAnnotations <- tbl(my_db, "ProteinAnnotations") %>%
  select(ProteinID, Description)

Events <- tbl(my_db, "Events") %>%
  select(-RT, -LeftRT, -RightRT, -SN, -FileID, -Intensity)

EventAreaAnnotations <- tbl(my_db, "EventAreaAnnotations") %>%
  select(EventID, QuanResultID)

PrecursorIonAreaSearchSpectra <- tbl(my_db, "PrecursorIonAreaSearchSpectra")

# Here is where all the scan information, including mass and charge
SpectrumHeaders <- tbl(my_db, "SpectrumHeaders") %>%
  select(SpectrumID, FirstScan, Mass, Charge, RetentionTime) %>%
  collect()

# Grab intensities
MassPeaks <- tbl(my_db, "MassPeaks") %>%
  select(MassPeakID, Intensity) %>%
  collect()

# Here are all the areas
events_joined <- inner_join(Events, EventAreaAnnotations) %>% 
  inner_join(PrecursorIonAreaSearchSpectra) %>%
  collect()  %>%
  rename(m_z = Mass)

# NOTE: mass and m/z are probably not correct right now! I have check them in more detail

# Join areas to the spectrum IDs. Not all spectrum IDs have areas!
spectra <- left_join(SpectrumHeaders, events_joined, by = c("SpectrumID" = "SearchSpectrumID")) %>%
  inner_join(MassPeaks, by = c("SpectrumID" = "MassPeakID"))

# Build protein/peptide list
mapped <- inner_join(PeptidesProteins, Peptides, by = c("PeptideID" = "PeptideID")) %>% 
  inner_join(ProteinAnnotations) %>%
  collect() %>%
  mutate(Description = str_match(Description, "^>([a-zA-Z0-9._]+)\\b")[,2]) %>%
  group_by(PeptideID) %>%
  summarize(Sequence = unique(Sequence), SpectrumID = unique(SpectrumID), Description = paste(Description, collapse = "; "))

# Join peptide info to mass/area/charge/etc.
all_data <- right_join(spectra, mapped, by=c("SpectrumID" = "SpectrumID")) %>% 
  left_join(PEP, by = c("PeptideID" = "PeptideID")) %>%
  select(PeptideID, Sequence, Description, Area, Mass, m_z, Charge, Intensity, FirstScan, PEP, `q-Value`)



# Intensities come from MassPeak table, MassPeakID == SpectrumID == SearchSpectrumID?

