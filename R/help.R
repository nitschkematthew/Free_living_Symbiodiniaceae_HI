
# Convert a DNAStringSet into a data.frame to use with tidyverse functions
DNAStringSet_to_df <- function(DNAStringSet){
  seq_df <- data.frame(names = names(DNAStringSet),
                       seqs = paste(DNAStringSet))
  return(seq_df)
}

# Convert a data.frame (seq_df) into a DNAStringSet to use with Biostrings functions
df_to_DNAStringset <- function(df, seqs = "seqs", names = "names"){
  DNAstr <- DNAStringSet(df[[seqs]])
  names(DNAstr) <- df[[names]]
  return(DNAstr)
}

# FIX THIS
# Convert data.frame (seq_df) into a DNAbin object to use with ape/insect functions
# df_to_DNAbin <- function(seq_df, seqs = "seqs", names = "names"){
#   DNA <- as.character(seq_df[[seqs]])
#   names(DNA) <- seq_df[[names]]
#   return(DNA)
# }

# Convert DNAbin to DNAStringSet to use with Biostrings functions
DNAbin_to_DNAStringSet <- function(DNAbin){
  DNAstrings <- DNAbin %>% as.character() %>% lapply(.,paste0, collapse="") %>% unlist() %>% DNAStringSet()
  return(DNAstrings)
}

# Convert DNAstringset to DNAbin
DNAStringSet_to_DNAbin <- function(DNAStringSet){
  DNAbin <- as.DNAbin(DNAStringSet)
  return(DNAbin)
}

# Convert DNAbin to data.frame (seq_df) to use with tidyverse functions
DNAbin_to_df <- function(DNAbin){
  seq_df <- data.frame(names = labels(DNAbin),
                       seqs = paste(DNAbin))
  return(seq_df)
}

# Read fasta file into df
fasta_to_df <- function(filepath){
  fasta_stringset <- readDNAStringSet(filepath)
  seq_df <- DNAStringSet_to_df(fasta_stringset) 
  return(seq_df)
}

# Write seq_df to .fasta file

# Extract accession numbers (assuming each accession finishes with a .1) from DNAstringset to use in database merging
extract_accn <- function(DNAStringSet){
  names <- names(DNAStringSet)
  accn <- stringi::stri_extract(names, regex='[^.]*')
  return(accn)
}

# Dereplicate a DNAStringSet and concatenate the names together with | as the delimiter
derep_cat_names <- function(dnaSet){
  if (!is(dnaSet, "DNAStringSet")) {
    dnaSet <- DNAStringSet(dnaSet)
  }
  if (is.null(names(dnaSet))) {
    message("No names attribute found in dnaSet object...", 
            "using artifically generated names")
    names(dnaSet) <- paste("read", 1:length(dnaSet), 
                           sep = "-")
  }
  
  seq_df <- DNAStringSet_to_df(dnaSet) %>%
    distinct(names, .keep_all = TRUE) %>%
    ungroup() %>%
    group_by(seqs) %>%
    summarise(names = paste(names, collapse = '|')) %>%
    ungroup()
  
  seq_ss <- df_to_DNAStringset(seq_df, "seqs", "names")
  return(seq_ss)
}

# Function to retrieve the query from genbank nucleotide database in fasta format using accession numbers. 
# If the accession number no longer exists in the database, include the accession in the output as 404NotFound
fetch_seqs <- function(query){
  recs <- try(entrez_fetch(db = "nuccore", id = query, rettype = "fasta"))
  if(str_detect(recs[1], "Error|Fail")){
    recs <- paste0(">", query, "_404NotFound\nAAAAAAAAAA\n\n")
    return(recs)
  }
  else{
    return(recs)
  }
}

# Run virtualPCR from insect in a loop using multiple primer pairs
virtualPCR_multiplex <- function(DNAbin, fwd_list, rev_list){
  if(length(fwd_list) != length(rev_list)){message("Error: Must have equal number of primer pair combinations")}
  else{
    PCR <- list()
    for(i in 1:length(fwd_list)){
      trim <- virtualPCR(DNAbin, up = fwd_list[i], down = rev_list[i], 
                         minfsc = 50, minrsc = 50, rcdown = FALSE, trimprimers = FALSE)
      trim <- DNAbin_to_DNAStringSet(trim)
      PCR[[i]] <- trim
      message(paste("Finished round ", i))
    }
    return(PCR)
  }
}

# Subset phyloseq object inside a function
get_sample_subsets_specific <- function(ps, factor){
  sample_subset <- sample_data(ps)[ which(sample_data(ps)$Specific == factor),]
  phy_subset <- merge_phyloseq(tax_table(ps), 
                               otu_table(ps, taxa_are_rows = FALSE),
                               refseq(ps),
                               sample_subset)
  return(phy_subset)
}

get_sample_subsets_culture <- function(ps, factor){
  sample_subset <- sample_data(ps)[ which(sample_data(ps)$Culture_genus == factor),]
  phy_subset <- merge_phyloseq(tax_table(ps), 
                               otu_table(ps, taxa_are_rows = FALSE),
                               refseq(ps),
                               sample_subset)
  return(phy_subset)
}

# Create network data from phyloseq object
ps_network_data <- function(ps, referenceDB = "", Specific_factor = "", genus = "Symbiodinium"){
  
  if(missing(Specific_factor)) {
    
    # n <- length(sample_names(ps))
    # flist <- filterfun(kOverA(n, 0))
    # a <- filter_taxa(ps, flist)
    # core <- prune_taxa(a, ps)
    
    taxlist <- as.data.frame(tax_table(ps)) %>%
      select(Species) %>%
      distinct(Species)
    
    ASVdata <- enframe(taxa_sums(ps)) %>%
      select(names = name, abundance = value) %>%
      mutate(origin = "ASV")
    
    SDB <- readDNAStringSet(referenceDB) %>%
      DNAStringSet_to_df() %>%
      filter(str_detect(names, paste(taxlist$Species, collapse = "|"))) %>%
      mutate(names = str_remove(names, paste0("Eukaryota;Dinoflagellata;Dinophyceae;Suessiales;Symbiodiniaceae;",genus,";")))
    
    SDBdata <- SDB %>%
      select(names) %>%
      mutate(abundance = mean(ASVdata$abundance),
             origin = "Reference")
    
    network_data <- rbind(SDBdata, ASVdata)
    
    return(network_data)
      
  } else {
    ps_sf <- get_sample_subsets_specific(ps, factor = Specific_factor)
    ps_culture <- get_sample_subsets_culture(ps, factor = genus)
    ps <- merge_phyloseq(ps_sf, ps_culture)
    
    # n <- length(sample_names(ps))
    # flist <- filterfun(kOverA(n, 0))
    # a <- filter_taxa(ps, flist)
    # core <- prune_taxa(a, ps)
    
    taxlist <- as.data.frame(tax_table(ps)) %>%
      select(Species) %>%
      distinct(Species)
    
    ASVdata <- enframe(taxa_sums(ps)) %>%
      select(names = name, abundance = value) %>%
      mutate(origin = "ASV")
    
    SDB <- readDNAStringSet(referenceDB) %>%
      DNAStringSet_to_df() %>%
      filter(str_detect(names, paste(taxlist$Species, collapse = "|"))) %>%
      mutate(names = str_remove(names, paste0("Eukaryota;Dinoflagellata;Dinophyceae;Suessiales;Symbiodiniaceae;",genus,";")))
    
    SDBdata <- SDB %>%
      select(names) %>%
      mutate(abundance = mean(ASVdata$abundance),
             origin = "Reference")
    
    node_data <- rbind(SDBdata, ASVdata)
    
    both_list <- psmelt(ps) %>%
      filter(Abundance > 0,
             Type != "Host") %>%
      group_by(OTU, Type) %>%
      top_n(1, Abundance) %>%
      ungroup() %>%
      group_by(OTU) %>%
      filter(n() > 1) %>%
      select(names = OTU, Type) %>%
      ungroup() %>%
      group_by(names, Type) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      group_by(names) %>%
      filter(n() > 1) %>%
      dplyr::slice(1) %>%
      select(names) %>%
      mutate(label = "Both")

      network_data <- left_join(node_data, both_list) %>%
      mutate(label = case_when(origin == "Reference" & is.na(label) ~ "Reference",
                              origin == "ASV" & is.na(label) ~ "Environment",
                              origin == "ASV" & label == "Both" ~ "Culture and Environment"))
    
    return(network_data)
  }
}

# Create sequence similarity network k-distances
ps_SSN <- function(ps, referenceDB = "", Specific_factor = "", genus = "Symbiodinium", k = 6){
  
  if(missing(Specific_factor)) {
    taxlist <- as.data.frame(tax_table(ps)) %>%
      select(Species) %>%
      distinct(Species)
    
    ASV <- refseq(ps) %>%
      DNAStringSet_to_df()
    
    SDB <- readDNAStringSet(referenceDB) %>%
      DNAStringSet_to_df() %>%
      filter(str_detect(names, paste(taxlist$Species, collapse = "|"))) %>%
      mutate(names = str_remove(names, paste0("Eukaryota;Dinoflagellata;Dinophyceae;Suessiales;Symbiodiniaceae;",genus,";")))
    
    SDB_ASV_seqs <- rbind(SDB, ASV)
    
    kdist <- SDB_ASV_seqs %>%
      df_to_DNAStringset() %>%
      DNAStringSet_to_DNAbin() %>%
      kdistance(k = k, residues = "DNA", method = "edgar") %>%
      as.matrix()
    
    long_dist <- melt(kdist)[melt(lower.tri(kdist))$value,] %>%
      mutate(value = (1-value)*100) %>%
      select(from = Var1, to = Var2, kdist = value)
    
    return(long_dist)
  }
  else{
    ps <- get_sample_subsets_specific(ps, factor = Specific_factor)
    
    taxlist <- as.data.frame(tax_table(ps)) %>%
      select(Species) %>%
      distinct(Species)
    
    ASV <- refseq(ps) %>%
      DNAStringSet_to_df()
    
    SDB <- readDNAStringSet(referenceDB) %>%
      DNAStringSet_to_df() %>%
      filter(str_detect(names, paste(taxlist$Species, collapse = "|"))) %>%
      mutate(names = str_remove(names, paste0("Eukaryota;Dinoflagellata;Dinophyceae;Suessiales;Symbiodiniaceae;",genus,";")))
    
    SDB_ASV_seqs <- rbind(SDB, ASV)
    
    kdist <- SDB_ASV_seqs %>%
      df_to_DNAStringset() %>%
      DNAStringSet_to_DNAbin() %>%
      kdistance(k = k, residues = "DNA", method = "edgar") %>%
      as.matrix()
    
    long_dist <- melt(kdist)[melt(lower.tri(kdist))$value,] %>%
      mutate(value = (1-value)*100) %>%
      select(from = Var1, to = Var2, kdist = value)
    
    return(long_dist)
  }
}

ps_SSN_all <- function(ps, referenceDB = "", k = 6){
  
  taxlist <- as.data.frame(tax_table(ps)) %>%
    select(Species) %>%
    distinct(Species)
  
  ASV <- refseq(ps) %>%
    DNAStringSet_to_df()
  
  SDB <- readDNAStringSet(referenceDB) %>%
    DNAStringSet_to_df() %>%
    filter(str_detect(names, paste(taxlist$Species, collapse = "|"))) %>%
    mutate(names = str_remove(names, "Eukaryota;Dinoflagellata;Dinophyceae;Suessiales;Symbiodiniaceae;")) %>%
    separate(names, into = c("Genus","names"), sep = ";") %>%
    select(names, seqs)
  
  SDB_ASV_seqs <- rbind(SDB, ASV)
  
  kdist <- SDB_ASV_seqs %>%
    df_to_DNAStringset() %>%
    DNAStringSet_to_DNAbin() %>%
    kdistance(k = k, residues = "DNA", method = "edgar") %>%
    as.matrix()
  
  long_dist <- melt(kdist)[melt(lower.tri(kdist))$value,] %>%
    mutate(value = (1-value)*100) %>%
    select(from = Var1, to = Var2, kdist = value)
  
  return(long_dist)
}

ps_network_data_all <- function(ps, referenceDB = ""){
  
  taxlist <- as.data.frame(tax_table(ps)) %>%
    select(Genus, Species) %>%
    distinct(Genus, Species)
  
  ASVdata <- enframe(taxa_sums(ps)) %>%
    select(names = name, abundance = value) %>%
    mutate(origin = "ASV")
  
  genuslist <- as.data.frame(tax_table(ps)) %>%
    tibble::rownames_to_column(var = "names") %>%
    select(names, Genus)
  
  ASVdata <- left_join(ASVdata, genuslist)
  
  SDB <- readDNAStringSet(referenceDB) %>%
    DNAStringSet_to_df() %>%
    filter(str_detect(names, paste(taxlist$Species, collapse = "|"))) %>%
    mutate(names = str_remove(names, "Eukaryota;Dinoflagellata;Dinophyceae;Suessiales;Symbiodiniaceae;")) %>%
    separate(names, into = c("Genus","names"), sep = ";")
  
  SDBdata <- SDB %>%
    select(names, Genus) %>%
    mutate(abundance = mean(ASVdata$abundance),
           origin = "Reference") %>%
    select(names, abundance, origin, Genus)
  
  network_data <- rbind(SDBdata, ASVdata)
  
  both_list <- psmelt(ps) %>%
    filter(Abundance > 0,
           Type != "Host") %>%
    group_by(OTU, Type) %>%
    top_n(1, Abundance) %>%
    ungroup() %>%
    group_by(OTU) %>%
    filter(n() > 1) %>%
    select(names = OTU, Type) %>%
    ungroup() %>%
    group_by(names, Type) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    group_by(names) %>%
    filter(n() > 1) %>%
    dplyr::slice(1) %>%
    select(names) %>%
    mutate(label = "Both")
  
  network_data <- left_join(network_data, both_list) %>%
    mutate(label = case_when(origin == "Reference" & is.na(label) ~ "Reference",
                             origin == "ASV" & is.na(label) ~ "Environment",
                             origin == "ASV" & label == "Both" ~ "Culture and Environment"))
  
  return(network_data)
}
