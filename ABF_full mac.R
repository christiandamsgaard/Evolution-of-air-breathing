###########################################################################
#### SCRIPT FOR THE EVOLUTION OF AIR-BREATHING AND AMPHIBOUS BEHAVIOUR ####
####             by Christian Damsgaard and Vikram Baliga              ####
###########################################################################



###########################################
#### READ ME BEFORE RUNNING THE SCRIPT ####
###########################################

## Before running the scripts, you should do the following

# 1. Make an new folder that will be the working directory and specify its path here
setwd("/Users/christiandamsgaard/Dropbox/Projects/Air-breathing review_full/")                     # Set working directory

# 2. Download, unzip and paste "Complete “all-taxon assembled” chronograms" and "PFC_short_classification.csv" from "https://fishtreeoflife.org/downloads/" into the working directory. 

# 3. Download and paste "Data Set 1.xlsx" into the working directory. To save time, you can further 
#    download the maximum clade credibility tree "mcc.nexus", and outputs from sections 1 and 3 called
#    "times_bpp_Organ.txt","counts_bpp_Organ.txt","times_bpp_LS.txt","counts_bpp_LS.txt"

# 4. Within the working directory folder make a folder called "Figures".

# 5. Install the packages specified under ## LOAD PACKAGES ## using the "Tools" --> "Install Packages..." option. 

#    To prepare running the model, first run the lines from the sections #### LOAD PACKAGES, COLOR SCHEMES, FUNCTIONS, IMPORT TREES, IMPORT DATA #### 
#    From there, you can run modelling within section 1-2 on air-breathing evolution or section 3-4 on amphibious evolution.
#    Sections 1-2 and sections 3-4 are followed by lines that code the figures. 
#    Sections 1-2 and sections 3-4 each take up ~ 10GB ram, so if your computer has < 16GB memory,
#        consider running Sections 1-2, plot the associated data, and shut down R to clear the memory. 
#        Then re-open R, and run sections 3-4, plot the associated data. 



############################################################################
#### LOAD PACKAGES, COLOR SCHEMES, FUNCTIONS, IMPORT TREES, IMPORT DATA ####
############################################################################

## LOAD PACKAGES ##
lapply(c("diversitree","vcd","ape","phyloseq","readxl","geiger","phytools","ggplot2", "phangorn","splitstackshape","plotrix","cowplot","gridExtra","RColorBrewer"), library, character.only = TRUE) 



## COLOR SCHEMES FOR ACTA PHYSIOLOGICA##
cols<-c(rgb(red=79, green=129, blue=189, maxColorValue = 255),
        rgb(red=180, green=42, blue=81, maxColorValue = 255),
        rgb(red=254, green=254, blue=218, maxColorValue = 255),
        rgb(red=127, green=127, blue=127, maxColorValue = 255),
        rgb(red=30 , green=200, blue=110, maxColorValue = 255),
        rgb(red=250, green=180, blue=140, maxColorValue = 255),
        rgb(red=000, green=000, blue=000, maxColorValue = 255),
        rgb(red=255, green=255, blue=255, maxColorValue = 255))



## CUSTOM FUNCTION ##
# Prune tree function
keep.tip                 <- function(tree,tip) {
  drop.tip(tree,setdiff(tree$tip.label,tip))
}


## IMPORT PHYLOGENIES ##
# Load basal sarcopterygean phylogeny (source Betancur-R 2017)
sarc                     <-read.tree(text="(Latimeria_chalumnae:409.4,(Neoceratodus_forsteri:241.804369,(Protopterus_aethiopicus_annectens:103.2,Lepidosiren_paradoxa:103.2)100:138.604369)100:167.595631)58;")
sarc$root.edge           <-424.8-409.4 
elas                     <-read.tree(text = "(Hydrolagus:399.3989352,(((Scyliorhinus:139.3745,Mustelus:139.3745)29:59.6255,Squalus:199)39:65.67780412,Dasyatis:264.6778041)140:134.7211311);")
elas$root.edge           <-63.00106
agna                     <-read.tree(text = "(Myxine:470.51250000,(Lampetra:16.00000000,Petromyzon:16.00000000)'14':454.51250000);")
agna$root.edge           <-144.4875

# Import 100 versions of the actinopterygean phylogeny
trees                    <-read.tree(file = "actinopt_full.trees")                                 # Import 100 BPP trees from the working directory (can be downloaded from https://fishtreeoflife.org)


# Bind sarcopterygean phylogeny to each BPP tree
for (i in 1:100) {
  
  hor <- as.data.frame(tree_layout(trees[[i]])$edgeDT)                                             # Data frame with the divergence times of all branches in phylogeny
  crown<-max(max(hor$xright)-hor$xleft)                                                            # The root node time
  trees[[i]]$root.edge<-100                                                                        # Add root note edge length
  trees[[i]]<-bind.tree(trees[[i]],sarc,position = 424.8-crown)                                    # Add the sarcopterygean (sarc) tree each BPP tree
  
  hor <- as.data.frame(tree_layout(trees[[i]])$edgeDT)                                             # Data frame with the divergence times of all branches in phylogeny
  crown<-max(max(hor$xright)-hor$xleft)                                                            # The root node time
  trees[[i]]$root.edge<-100                                                                        # Add root note edge length
  trees[[i]]<-bind.tree(trees[[i]],elas,position = 462.4-crown)                                    # Add the sarcopterygean (sarc) tree each BPP tree
  
  hor <- as.data.frame(tree_layout(trees[[i]])$edgeDT)                                             # Data frame with the divergence times of all branches in phylogeny
  crown<-max(max(hor$xright)-hor$xleft)                                                            # The root node time
  trees[[i]]$root.edge<-200                                                                        # Add root note edge length
  trees[[i]]<-bind.tree(trees[[i]],agna,position = 615-crown)                                    # Add the sarcopterygean (sarc) tree each BPP tree

  
}


# Generate a maximum clade credibility tree based on the 100 BPP trees using first 
# and second line OR by loading it using the read.tree function in line 3.
#mcc                      <-maxCladeCred(trees)                                                     # Generate a maximum clade credibility tree (mcc) tree using the 100 BPP trees
#write.tree(mcc,"mcc.nexus")                                                                        # Write the mcc tree to the working directory
mcc                      <-read.tree("mcc.nexus")                                                   # Load maximum clade credibility tree from working directory


# Prune each BPP trees to genus level #
for (i in 1:100){
  species                <- as.data.frame(trees[[i]]$tip.label)                                    # Load all species names from the species-level BPP tree
  colnames(species)      <- c("Names")            
  a                      <- cSplit(indt = species,splitCols = "Names",sep="_",type.convert=FALSE)  # Separate genus and species name
  genera                 <- unique(a$Names_1)                                                      # Vector with all genus names
  trees[[i]]             <- keep.tip(trees[[i]],trees[[i]]$tip.label[match(genera,a$Names_1)])     # Prune tree to only represent 1 species per genus
  tip.labels             <- as.data.frame(trees[[i]]$tip.label); colnames(tip.labels)<-c("Names")  # Load all species names from the genus level BPP tree
  tip.labels             <- cSplit(indt = tip.labels,splitCols="Names",sep="_",type.convert=FALSE) # Separate genus and species name
  trees[[i]]$tip.label   <- tip.labels$Names_1                                                     # Use only the genus name
  trees[[i]]             <- ladderize(trees[[i]])                                                  # Ladderize the tree
}



# Prune MCC tree to genus level and save it as ge_tree <-- i.e. genus level MCC tree
species                  <- as.data.frame(mcc$tip.label)                                           # Load all species names from the species level MCC tree
colnames(species)        <- c("Names")                                                  
a                        <- cSplit(indt = species,splitCols = "Names", sep="_", type.convert=FALSE)# Separate genus and species name
genera                   <- unique(a$Names_1)                                                      # Vector with all genus names
ge_tree                  <- keep.tip(mcc,mcc$tip.label[match(genera,a$Names_1)])                   # Prune tree to only represent 1 species per genera
tip.labels               <- as.data.frame(ge_tree$tip.label); colnames(tip.labels)<-c("Names")     # Load all species names from the genus level BPP tree
tip.labels               <- cSplit(indt = tip.labels,splitCols ="Names",sep="_",type.convert=FALSE)# Separate genus and species name
ge_tree$tip.label        <- tip.labels$Names_1                                                     # Use only the genus name as tip label
ge_tree                  <- ladderize(ge_tree)                                                     # Ladderize the tree


# add shark phylogeny from timetree.org
elas                     <-read.tree(text = "(Hydrolagus:399.3989352,(((Scyliorhinus:139.3745,Mustelus:139.3745)29:59.6255,Squalus:199)39:65.67780412,Dasyatis:264.6778041)140:134.7211311);")
hor                      <-as.data.frame(tree_layout(elas)$edgeDT);elas_crown<-max(max(hor$xright)-hor$xleft) # Get coordinates of the shark tree
elas_oste_sarc_dt        <-462.4                                                                   # Divergence time between bony fishes and sharks.
elas$root.edge           <-elas_oste_sarc_dt-elas_crown                                            # Root node length
hor                      <-as.data.frame(tree_layout(ge_tree)$edgeDT)                              # Get coordinates on the bony fish tree
oste_crown               <-max(max(hor$xright)-hor$xleft)                                          # Root node age of the bony fish tree
ge_tree$root.edge        <-462.4-424.8                                                             # Root node length of the bony fish tree
ge_tree                     <-bind.tree(ge_tree,                                                   # Merge bony fish tree and shark tree
                                     elas,
                                     position = elas_oste_sarc_dt-oste_crown)

# add agnathan phylogeny (from timetree)
agna                     <-read.tree(text = "(Myxine:470.51250000,(Lampetra:16.00000000,Petromyzon:16.00000000)'14':454.51250000);")
hor                      <-as.data.frame(tree_layout(agna)$edgeDT);agna_crown<-max(max(hor$xright)-hor$xleft) # Get coordinates of the agnatha tree
agna_elas_oste_sarc_dt   <-615                                                                     # Divergence time between agnatha and sharks.
agna$root.edge           <-agna_elas_oste_sarc_dt-agna_crown                                       # Root node length of agnatha tree
hor                      <-as.data.frame(tree_layout(ge_tree)$edgeDT)                              # Get coordinates on the shark+bony fish tree
gnat_crown               <-max(max(hor$xright)-hor$xleft)                                          # Root node age of the shark+bony fish tree
ge_tree$root.edge           <-agna_elas_oste_sarc_dt-462.4                                         # Root node length of the shark+bony fish tree
ge_tree                     <-bind.tree(ge_tree,                                                   # Merge shark+bony fish tree and agnatha tree                                                   
                                     agna,
                                     position = agna_elas_oste_sarc_dt-462.4)



## IMPORT DATA ##
# Reduce species level data (data_sp) on air-breathing to genus level and save it as ge_data (i.e. genus level data)
data_sp                  <- read_xlsx("Data Set 1.xlsx",sheet = "species")                         # Load Data Set 1
data_sp                  <- cbind(data_sp$Species[!is.na(data_sp$Species)],
                                  data_sp$Organ[!is.na(data_sp$Organ)])                            # Reduce Data Set 1 to species and air-breathing organ (Organ)
species                  <- as.data.frame(data_sp[,1]); colnames(species)<-c("Names")              # List of all air-breathing species
a                        <- cSplit(indt = species,splitCols="Names",sep="_",type.convert=FALSE)    # Separate genus and species name
data_ge                  <- as.data.frame(cbind(a$Names_1,data_sp[,2]))                            # Data frame with genus names and air-breating organs                               
colnames(data_ge)        <- c("Genus","Organ")
data_ge                  <- aggregate(Organ~Genus,data = data_ge,unique)                           # Remove identical entries 


# List of amphibious genera (source Wright Turko 2016 JEB)
amphibious<-c("Anabas","Ctenopoma","Ctenopoma","Channa","Channa","Channa","Channa","Channa","Betta","Anguilla","Anguilla","Leuresthes","Leuresthes","Porichthys","Alticus","Alticus","Alticus","Alticus","Andamia","Andamia","Andamia","Andamia","Blennius","Blennius","Blennius","Blennius","Blennius","Blennius","Blennius","Blennius","Coryphoblennius","Entomacrodus","Entomacrodus","Entomacrodus","Entomacrodus","Hypsoblennius","Istiblennius","Praealticus","Salarias","Arcos","Gobiesox","Gobiesox","Lepadogaster","Pherallodiscus","Sicyases","Tomicodon","Tomicodon","Tomicodon","Tomicodon","Dialommus","Dialommus","Bellapiscis","Forsterygion","Helcogramma","Acanthoclinus","Neoceratodus","Brycon","Erythrinus","Hoplerythrinus","Copella","Pyrrhulina","Lepidocephalichthys","Misgurnus","Misgurnus","Anableps","Anableps","Cyprinodon","Jordanella","Fundulus","Fundulus","Fundulus","Fundulus","Fundulus","Fundulus","Fundulus","Lucania","Aphyosemion","Aphyosemion","Aphyosemion","Aphyosemion","Aphyosemion","Aphyosemion","Fundulopanchax","Fundulopanchax","Gambusia","Anablepsoides","Anablepsoides","Anablepsoides","Anablepsoides","Anablepsoides","Anablepsoides","Anablepsoides","Anablepsoides","Anablepsoides","Cynodonichthys","Kryptolebias","Kryptolebias","Laimosemion","Laimosemion","Laimosemion","Galaxias","Galaxias","Galaxias","Galaxias","Galaxias","Galaxias","Galaxias","Galaxiella","Galaxiella","Neochanna","Neochanna","Neochanna","Neochanna","Neochanna","Neochanna","Bostrychus","Oxyeleotris","Chlamydogobius","Ctenogobius","Gillichthys","Gillichthys","Kelloggella","Quietula","Sicyopterus","Apocryptes","Boleophthalmus","Boleophthalmus","Boleophthalmus","Boleophthalmus","Boleophthalmus","Boleophthalmus","Periophthalmodon","Periophthalmodon","Periophthalmodon","Periophthalmus","Periophthalmus","Periophthalmus","Periophthalmus","Periophthalmus","Periophthalmus","Periophthalmus","Periophthalmus","Periophthalmus","Periophthalmus","Periophthalmus","Periophthalmus","Periophthalmus","Periophthalmus","Periophthalmus","Pseudapocryptes","Pseudapocryptes","Scartelaos","Scartelaos","Scartelaos","Scartelaos","Zappa","Parioglossus","Lepidogalaxias","Lepidosiren","Protopterus","Protopterus","Protopterus","Protopterus","Artedius","Ascelichthys","Clinocottus","Clinocottus","Clinocottus","Clinocottus","Leptocottus","Oligocottus","Oligocottus","Oligocottus","Taurulus","Girella","Apodichthys","Pholis","Pholis","Pholis","Apodichthys","Zoarces","Anoplarchus","Cebidichthys","Xiphister","Xiphister","Erpetoichthys","Polypterus","Callichthys","Channallabes","Clarias","Clarias","Clarias","Clarias","Clarias","Clarias","Clarias","Clarias","Centrochir","Platydoras","Heteropneustes","Heteropneustes","Heteropneustes","Mastacembelus","Mastacembelus","Mastacembelus","Monopterus","Monopterus","Synbranchus")
amphibious<-unique(amphibious)


# Generate named vectors for 
     # respiratory medium (AW: values A(ir) or W(ater))
     # respiratory medium as binary character (AW_bi: values 1 and 0)
     # Air-breathing organ (Organ)
     # Lifestyle (LS: values (Am)phibious and (Aq)uatic)

genus                    <- c(ge_tree$tip.label)                                                   # Generate a vector with all the genera in the phylogeny
df                       <- as.data.frame(genus)                                                   # Transform into a data frame
for (i in 1:length(genus)) {                                                                  
  df$AW[i]               <-ifelse(df$genus[i] %in% data_ge$Genus==T,                               # If the genera contains air-breathing species,              
                                  "A",                                                             # then assign air-breathing (A) as breathing mode
                                  "W")                                                             # else assign water-breathing (W) as breathing mode
  
  df$AW_bi[i]               <-ifelse(df$genus[i] %in% data_ge$Genus==T,                            # If the genera contains air-breathing species,              
                                     1,                                                            # then assign air-breathing (1) as breathing mode
                                     0)                                                            # else assign water-breathing (0) as breathing mode
  
  df$Organ[i]            <-ifelse(df$genus[i] %in% data_ge$Genus==T,                               # If the genera contains air-breathing species,              
                                  as.character(data_ge$Organ[match(df$genus[i],data_ge$Genus)]),   # then assign its air-breathing organ  
                                  "Gill")                                                          # else assign gill
  
  df$LS[i]               <-ifelse(df$genus[i] %in% amphibious==T,                                  # If the genera contains amphibious species,              
                                  "Am",                                                            # then assign amphibious (Am) as lifestyle (LS)
                                  "Aq")                                                            # else assign aquatic (Aq) as lifestyle (LS)
}     

# Generate named vectors for the four characters
AW                       <- setNames(df$AW,df$genus)
AW_bi                    <- setNames(df$AW_bi,df$genus)
Organ                    <- setNames(df$Organ,df$genus)
LS                       <- setNames(df$LS,df$genus)




# Set up a data matrix (tax) containing order, familiy, genus, and species for each species in the phylogeny
species                  <- read_xlsx("Data Set 1.xlsx",sheet = "species")$Species                # Load Data Set 1
species                  <- sub("_", "\\ ",species)                                               # Separate genus and species by space
species                  <- sub("_", "\\ ",species)                                               # Separate genus and species by space (to remove double underscores)
tax                      <- read.csv("PFC_short_classification.csv")                              # Load the fish taxonomy

# Add taxonomy data to basal sarcopterygeans, sharks and agnathans
tax2<-as.data.frame(rbind(c("Coelacanthiformes", "Latimeriidae","Latimeria","Latimeria chalumnae"),
                          c("Ceratodontiformes","Neoceratodontidae","Neoceratodus","Neoceratodus forsteri"),
                          c("Lepidosireniformes","Protopteridae","Protopterus","Protopterus aethiopicus annectens"),
                          c("Lepidosireniformes","Lepidosirenidae","Lepidosiren","Lepidosiren paradoxa"),
                          c("Chimaeriformes","Chimaeridae", "Hydrolagus","Hydrolagus affinis"),
                          c("Carcharhiniformes","Scyliorhinidae","Scyliorhinus","Scyliorhinus canicula"),
                          c("Carcharhiniformes","Triakidae","Mustelus","Mustelus asterias"),
                          c("Squaliformes","Squalidae","Squalus","Squalus acanthias"),
                          c("Myliobatiformes","Dasyatidae","Dasyatis","Dasyatis marmorata"),
                          c("Myxiniformes","Myxinidae","Myxine","Myxine glutinosa"),
                          c("Petromyzontiformes","Petromyzontidae","Lampetra","Lampetra fluviatilis"),
                          c("Petromyzontiformes","Petromyzontidae","Petromyzon","Petromyzon marinus")     
                          )
                    )
colnames(tax2)           <-colnames(tax)
tax                      <-rbind(tax2,tax) 

##################################################################
#### QUESITONS ON AIR-BREATHING AND AMPHIBIOUS FISH DIVERSITY ####
##################################################################


# How many species of air-breathing fishes?
unique(species)
length(unique(species))
# Result: 656 species


# How many genera of air-breating fishes?
genera<-tax[match(species,tax$genus.species),3]
unique(genera)
length(unique(genera))
# Result: 129 genera


# How many families of air-breating fishes?
families<-tax[match(species,tax$genus.species),2]
unique(families)
length(unique(families))
# Result: 41 families


# How many orders of air-breating fishes?
orders<-tax[match(species,tax$genus.species),1]
unique(orders)
length(unique(orders))
# Result: 22 orders 

## AMPHIBIOUS FISHES
# How many genera of amphibious fishes?
unique(amphibious)
length(unique(amphibious))
# 93 genera

# How many families of amphibious fishes?
families<-tax[match(amphibious,tax$genus),2]
unique(families)
length(unique(families))
# Result: 41 families

# How many orders of amphibious fishes?
orders<-tax[match(amphibious,tax$genus),1]
unique(orders)
length(unique(orders))
# Result: 19 orders




#######################
#### MODEL FITTING ####
#######################

## MODEL AIR-BREATHING EVOLUTION
# Note that his may take around 24 hours
ace_ER                   <- fitDiscrete(ge_tree, Organ, type="discrete", model = "ER")             # Fits an equal rates model to the data
ace_SYM                  <- fitDiscrete(ge_tree, Organ, type="discrete", model = "SYM")            # Fits a symmetric rates model to the data  
ace_ARD                  <- fitDiscrete(ge_tree, Organ, type="discrete", model = "ARD")            # Fits an all-rates-different model to the data  

# Extract likelihoods, df, AICc values and perform likelihood ratio tests
S                        <-c(ace_ER$opt$lnL,ace_ER$opt$k,ace_ER$opt$aicc,
                             ace_SYM$opt$lnL,ace_SYM$opt$k,ace_SYM$opt$aicc,
                             ace_ARD$opt$lnL,ace_ARD$opt$k,ace_ARD$opt$aicc,
                             pchisq(abs(2*(ace_ER$opt$lnL-ace_ARD$opt$lnL)) , 
                                    ace_ARD$opt$k-ace_ER$opt$k, lower.tail=FALSE), 
                             pchisq(abs(2*(ace_ER$opt$lnL-ace_SYM$opt$lnL)) , 
                                    ace_SYM$opt$k-ace_ER$opt$k, lower.tail=FALSE) ,
                             pchisq(abs(2*(ace_SYM$opt$lnL-ace_ARD$opt$lnL)) , 
                                    ace_ARD$opt$k-ace_SYM$opt$k, lower.tail=FALSE))
names(S)<-c("lnL_ER","df_ER","AICc_ER","lnL_SYM","df_SYM","AICc_SYM","lnL_ARD","df_ARD","AICc_ARD",
            "p_ERvsARD","p_ERvsSYM","p_SYMvsARD")
S
# Conclusion: ARD > SYM > ER, but ARD and SYM disregards transitions to/from lung. Use ER model for air-breathing evolution



## MODEL AMPHIBIOUS BEHAVIOUR EVOLUTION 
# Note that his may take around 24 hours
ace_ER                   <- fitDiscrete(ge_tree, LS, type="discrete", model = "ER")             # Fits an equal rates model to the data
ace_SYM                  <- fitDiscrete(ge_tree, LS, type="discrete", model = "SYM")            # Fits a symmetric rates model to the data  
ace_ARD                  <- fitDiscrete(ge_tree, LS, type="discrete", model = "ARD")            # Fits an all-rates-different model to the data  

# Extract likelihoods, df, AICc values and perform likelihood ratio tests
S                        <-c(ace_ER$opt$lnL,ace_ER$opt$k,ace_ER$opt$aicc,
                             ace_SYM$opt$lnL,ace_SYM$opt$k,ace_SYM$opt$aicc,
                             ace_ARD$opt$lnL,ace_ARD$opt$k,ace_ARD$opt$aicc,
                             pchisq(abs(2*(ace_ER$opt$lnL-ace_ARD$opt$lnL)) , 
                                    ace_ARD$opt$k-ace_ER$opt$k, lower.tail=FALSE), 
                             pchisq(abs(2*(ace_ER$opt$lnL-ace_SYM$opt$lnL)) , 
                                    ace_SYM$opt$k-ace_ER$opt$k, lower.tail=FALSE) ,
                             pchisq(abs(2*(ace_SYM$opt$lnL-ace_ARD$opt$lnL)) , 
                                    ace_ARD$opt$k-ace_SYM$opt$k, lower.tail=FALSE))
names(S)<-c("lnL_ER","df_ER","AICc_ER","lnL_SYM","df_SYM","AICc_SYM","lnL_ARD","df_ARD","AICc_ARD",
            "p_ERvsARD","p_ERvsSYM","p_SYMvsARD")
S
# Comment: ARD > SYM = ER. However, ARD predicts exceptionally high rates of character evolution. Use ER model for amphibious behaviour evolution.



###########################################################################
####                            SECTION 1                              ####
#### NUMBER AND TIMING OF TRANSITIONS BETWEEN WATER- AND AIR-BREATHING #### 
####         USING 100 BAYESIAN POSTERIOR PROBABILITY TREES            ####
####                                                                   ####
####     THERE ARE TWO OPTIONS (A OR B): A: RUNNING THE MODEL FROM     #### 
####              SCRATCH OR B: LOADING PREVIOUSLY SAVED RUN           ####
###########################################################################



### OPTION A: RUN THE MODEL FROM SCRATCH
# Loop that runs stochatic character mapping on each of the 100 BPP trees. 
# Requires around 1.5 GB of available memory, and takes around 48 hours. 
# The output of the loop is the files "counts_bpp_Organ.txt" and "times_bpp_Organ.txt" (available on GitHub)
# that contains the transition counts and transition timing, respectively, from each loop.

nsim                     <- 100                                                                    # Number of simulations per BPP tree
model                    <- "ER"                                                                   # Set model for character evolution to equal rates model
counts_bpp_Organ         <- matrix(ncol = 39)                                                      # Generate matrix to store counts of transitions
times_bpp_Organ          <- matrix(ncol = 2);colnames(times_bpp_Organ)<-c("mode","time")           # Generate matrix to store timing of transitions

for (k in 1:100) {
  ## STOCHASTIC CHARACTER MAPPING  
  SCM_bpp_Organ          <- make.simmap(trees[[k]],x = Organ,model = model,nsim = nsim)            # Generate stochatic character maps
  dSCM_bpp_Organ                   <- describe.simmap(SCM_bpp_Organ)                               # Summarize information from the stochastic character maps
  counts                 <- as.data.frame(countSimmap(SCM_bpp_Organ,states=rownames(Organ))$Tr)    # Calculate the number of transitions between states
  counts[,38]            <-counts[,3]+counts[,4]+counts[,5]+counts[,6]+counts[,7]                  # Summarize the number of transitions from gill to any ABO
  counts[,39]            <-counts[,8]+counts[,14]+counts[,20]+counts[,26]+counts[,32]              # Summarize the number of transitions from any ABO to gill
  names(counts)[38:39]   <-c("gains","losses")                                                    
  colnames(counts_bpp_Organ)   <-names(counts)
  counts_bpp_Organ             <-na.omit(rbind(counts_bpp_Organ,counts))                           # Add the count data to the counts_bpp_Organ matrix
  write.table(counts_bpp_Organ,file = paste("./counts_bpp_Organ.txt"))                             # Save the number of counts
  
  ## FIND TIME OF TRANSITIONS
  hor <- as.data.frame(tree_layout(trees[[k]])$edgeDT)                                             # Extract divergence times from all internal nodes in the tree in the vector called "hor" (as in horizontal lines in the phylogeny)
  ACE <- cbind(dSCM_bpp_Organ$ace[,1],1-dSCM_bpp_Organ$ace[,1])                                    # Extract Bayesian posterior probabilities from all internal nodes in the tree
  
  # Generate a matrix H, that finds the timing of gains and losses of air-breathing
  H<-matrix(nrow = length(hor$V1), ncol = 7)                                                       # Generate a matrix for data storage, where each row contains an edge in the tree
  for (i in 1:length(H[,1])) { H[i,1] <- max(hor$xright)-hor$xleft[i] }                            # Edge goes from this age
  for (i in 1:length(H[,1])) { H[i,2] <- max(hor$xright)-hor$xright[i] }                           # to this age
  for (i in 1:length(H[,1])) { H[i,3] <- mean(H[i,1],H[i,2]) }                                     # where age is the midpoint of that edge
  for (i in 1:length(H[,1])) { H[i,4] <- ACE[match(hor$V1[i],row.names(ACE))]}                     # Bayesian posterior probability at first point of the edge
  for (i in 1:length(H[,1])) { H[i,5] <- ACE[match(hor$V2[i],row.names(ACE))]}                     # Bayesian posterior probability at second point of the edge
  for (i in 1:length(H[,1])) { H[i,5] <- ifelse(is.na(ACE[match(hor$V2[i],row.names(ACE))])==T,    # If the branch ends in an extant species, set Bayesian posterior probability to 1 for air-breathers and to 0 for water-breathers
                                                unname(ifelse(Organ[match(trees[[k]]$tip[hor$V2[i]],names(Organ))]=="Gill",1,0)),
                                                H[i,5])} 
  for (i in 1:length(H[,1])) { H[i,6] <- ifelse(H[i,4]<0.5&H[i,5]>0.5,                             # If Bayesian posterior probability changes between <0.5 and >0.5, find the time where Bayesian posterior probability is 0.5 using linear regression (i.e. time of the origin of air-breathing). 
                                                (0.5-lm(c(H[i,4],H[i,5])~c(H[i,1],H[i,2]))$coefficients[1])/lm(c(H[i,4],H[i,5])~c(H[i,1],H[i,2]))$coefficients[2],
                                                NA)}                                              
  for (i in 1:length(H[,1])) { H[i,7] <- ifelse(H[i,4]>0.5&H[i,5]<0.5,                             # If Bayesian posterior probability changes between >0.5 and <0.5, find the time where Bayesian posterior probability is 0.5 using linear regression (i.e. time of the origin of water-breathing). 
                                                (0.5-lm(c(H[i,4],H[i,5])~c(H[i,1],H[i,2]))$coefficients[1])/lm(c(H[i,4],H[i,5])~c(H[i,1],H[i,2]))$coefficients[2],
                                                NA)}
  
  wa                     <- H[,6];   wa<-wa[!is.na(wa)]                                            # "wa" is a vector with all the times of transition from water to air. wa as in water to air (i.e. gain of air-breathing)
  aw                     <- H[,7];   aw<-aw[!is.na(aw)]                                            # "aw" is a vector with all the times of transition from air to water. aw as in air to water (i.e. loss of air-breathing)
  
  # Generate a dataframe (dataframe) with type of transition and their divergence times in first and second column
  M                      <- matrix(ncol = 2,nrow = (length(aw)+length(wa)))                        # Matrix to store data      
  for (i in 1:length(wa)) { M[i,1]<-"WA"} 
  for (i in 1:length(aw)) { j=i+length(wa);  M[j,1]<-"AW"} 
  for (i in 1:length(wa)) { M[i,2]<-wa[i]} 
  for (i in 1:length(aw)) { j=i+length(wa);  M[j,2]<-aw[i]} 
  dataframe              <- as.data.frame(M)
  colnames(dataframe)    <- c("mode","time")
  dataframe$time         <- as.numeric(as.character(dataframe$time))
  times_bpp_Organ        <-na.omit(rbind(times_bpp_Organ,dataframe))
  write.table(times_bpp_Organ,file = paste("./times_bpp_Organ.txt"))                               # Save the timing of origins
}

### OPTION B: LOADING PREVIOUSLY SAVED RUN
# Load SCM counts and times. This can be done if option A loop has already been run. 
counts_bpp_Organ         <-read.table("./counts_bpp_Organ.txt")
times_bpp_Organ          <-read.table("./times_bpp_Organ.txt")



###################################################################
####                        SECTION 2                          ####
#### INDENTIFY CLADES WHERE AIR-BREATHING EVOLVED AND WAS LOST #### 
####         USING THE MAXIMUM CLADE CREDIBILITY TREE          ####
###################################################################
# Requires around 10 GB of available memory, and takes around 6 hours for 1000 simulations

nsim                     <- 1000                                                                   # Number of simulations 
model                    <- "ER"                                                                   # Set model for character evolution to equal rates model

SCM_mcc_Organ            <- make.simmap(tree=ge_tree,x=Organ,nsim=nsim,model=model)                # Run stochastic character mapping on MCC tree
dSCM_mcc_Organ           <- describe.simmap(SCM_mcc_Organ)                                         # Summarize stochastic character mapping
dSCM_mcc_Organ


##########################################################
####          PLOTS RELATED TO SECTION 1 & 2          ####
#### REQUIRES THAT SECTION 1-2 HAS BEEN RUN OR LOADED ####
##########################################################


### FIGURE 5+S1: CIRCULAR PHYLOGENY COLOR MAPPED WITH PROBABILITY FOR AIR-BREATHING
# Prepare for plotting: Find edges where air-breathing was gained and lost
hor                      <- as.data.frame(tree_layout(ge_tree)$edgeDT)                           # Extract divergence times from all internal nodes in the tree in the vector called "hor" (as in horizontal lines in the phylogeny)
ACE                      <- cbind(dSCM_mcc_Organ$ace[,1],1-dSCM_mcc_Organ$ace[,1])               # Extract Bayesian posterior probabilities from all internal nodes in the tree

H<-matrix(nrow = length(hor$V1), ncol = 7)                                                       # Generate a matrix for data storage, where each row contains an edge in the tree
for (i in 1:length(H[,1])) { H[i,1] <- max(hor$xright)-hor$xleft[i] }                            # Edge goes from this age
for (i in 1:length(H[,1])) { H[i,2] <- max(hor$xright)-hor$xright[i] }                           # to this age
for (i in 1:length(H[,1])) { H[i,3] <- mean(H[i,1],H[i,2]) }                                     # where age is the midpoint of that edge
for (i in 1:length(H[,1])) { H[i,4] <- ACE[match(hor$V1[i],row.names(ACE))]}                     # Bayesian posterior probability at first point of the edge
for (i in 1:length(H[,1])) { H[i,5] <- ACE[match(hor$V2[i],row.names(ACE))]}                     # Bayesian posterior probability at second point of the edge
for (i in 1:length(H[,1])) { H[i,5] <- ifelse(is.na(ACE[match(hor$V2[i],row.names(ACE))])==T,    # If the edge ends in an extant species, set Bayesian posterior probability to 1 for air-breathers and to 0 for water-breathers
                                              unname(ifelse(Organ[match(ge_tree$tip[hor$V2[i]],names(Organ))]=="Gill",1,0)),
                                              H[i,5])} 
for (i in 1:length(H[,1])) { H[i,6] <- ifelse(H[i,4]<0.5&H[i,5]>0.5,                             # If Bayesian posterior probability changes between <0.5 and >0.5, find the time where Bayesian posterior probability is 0.5 (i.e. time of the origin of air-breathing). 
                                              (0.5-lm(c(H[i,4],H[i,5])~c(H[i,1],H[i,2]))$coefficients[1])/lm(c(H[i,4],H[i,5])~c(H[i,1],H[i,2]))$coefficients[2],
                                              NA)}                                              
for (i in 1:length(H[,1])) { H[i,7] <- ifelse(H[i,4]>0.5&H[i,5]<0.5,                             # If Bayesian posterior probability changes between >0.5 and <0.5, find the time where Bayesian posterior probability is 0.5 (i.e. time of the origin of water-breathing). 
                                              (0.5-lm(c(H[i,4],H[i,5])~c(H[i,1],H[i,2]))$coefficients[1])/lm(c(H[i,4],H[i,5])~c(H[i,1],H[i,2]))$coefficients[2],
                                              NA)}

# Find the branches in the tree, where air-breathing originates (wa.edges) and was lost (aw.edges)
m                      <- cbind(hor$V1[as.numeric(row.names((subset(as.data.frame(H),is.na(H[,6])==F))))],
                                hor$V2[as.numeric(row.names((subset(as.data.frame(H),is.na(H[,6])==F))))])
l                      <- cbind(hor$V1[as.numeric(row.names((subset(as.data.frame(H),is.na(H[,7])==F))))],
                                hor$V2[as.numeric(row.names((subset(as.data.frame(H),is.na(H[,7])==F))))])
g                      <- as.data.frame(ge_tree$edge)
signif(subset(as.data.frame(H),is.na(H[,7])==F),3)
aw.edges<-matrix(ncol=1,nrow=length(m[,2]))
for (i in 1:length(m[,1])) {
  aw.edges[i]          <-rownames(subset(g,g$V1==m[i,1]&g$V2==m[i,2]))
}

wa.edges<-matrix(ncol=1,nrow=length(l[,2]))
for (i in 1:length(l[,1])) {
  wa.edges[i]          <-rownames(subset(g,g$V1==l[i,1]&g$V2==l[i,2]))
}

#pdf(file = paste("./Figures/Fig. S1 for AR_2.pdf",sep = ""),width = 50,height = 50,useDingbats = F)
#plot(ge_tree,lwd = 0.5,fsize = 0.03,ftype = "i",show.tip.label = T,type = "fan",legend=F)
#edgelabels(text = length(wa.edges):1,frame="circle",bg=cols[2],edge=as.numeric(as.vector(wa.edges)))
#edgelabels(text = length(aw.edges):1,frame="circle",bg=cols[1],edge=as.numeric(as.vector(aw.edges)))
#dev.off()


# Save a density plot
obj                      <-densityMap(mergeMappedStates(tree = SCM_mcc_Organ,
                                                        old.states = c("Lung","Swimbladder","GIT","Mouth","Skin"),"ABO"),
                                    plot=FALSE,states=c("ABO","Gill"),
                                    res=1000)
obj$cols[1:1001]         <-colorRampPalette(cols[2:1], space="Lab")(1001)                          # Change color theme to Acta colors

# Fake plot used to (somehow) activate the axis function
plot(obj,lwd = 0.5,fsize = 0.1,ftype = "i",show.tip.label = F,type = "fan",legend=F)
axis_obj<-axis(1,pos=0,at=seq(max(nodeHeights(obj$tree)),100,by=-100),cex.axis=0.5,labels=FALSE)

# Plot Fig. 5
pdf(file = "./Figures/Figure 5 - air-breathing evolution.pdf",width = 19/2.54,height = 19/2.54,useDingbats = F)                # PDF saving options
plot(obj,lwd = 1,fsize = .000001,ftype = "i",show.tip.label = F,type = "fan",legend=F,no.margin=T) # Plot density map
for(i in 1:length(axis_obj)){                                                                      # Plot circular time axis
  a1<-0
  a2<-2*3.141528
  draw.arc(0,0,radius=axis_obj[i]-100,angle1 = a1,angle2 = a2,lwd=0.5,lty="dashed",col=make.transparent("black",0.5))}
edgelabels(pch=21,bg=cols[2],cex=2,edge=as.numeric(as.vector(wa.edges)))                       # Add red dot on edges where air-breathing evolved
edgelabels(pch=21,bg=cols[1],cex=2,edge=as.numeric(as.vector(aw.edges)))                       # Add red blue on edges where air-breathing was lost
dev.off()



### Fig. S1
# Large sized density map with genus, family and order names (Fig. S1)
obj_2                    <- obj                                                                    # Copy the first density map
obj_2$tree$tip.label     <- paste(obj_2$tree$tip.label," (",
                                  as.character(tax[match(obj_2$tree$tip.label,tax$genus),1]),", ",
                                  as.character(tax[match(obj_2$tree$tip.label,tax$genus),2]),")",
                                  sep = "")                                                        # Change tip labels to include genus, family and order

# Plot and save Fig. S1
pdf(file = paste("./Figures/Figure S1 - air-breathing evolution.pdf",sep = ""),width = 50,height = 50,useDingbats = F)
plot(obj_2,lwd = 1.5,fsize = .2,show.tip.label = F,type = "fan",legend=F,no.margin=T)
for(i in 1:length(axis_obj)){
  a1<-0
  a2<-2*3.141528
  draw.arc(0,0,radius=axis_obj[i]-100,angle1 = a1,angle2 = a2,lwd=0.1,lty="dashed",
           col=make.transparent("black",0.5))  
}
edgelabels(pch=21,bg=cols[2],cex=1.5,edge=as.numeric(as.vector(wa.edges)))
edgelabels(pch=21,bg=cols[1],cex=1.5,edge=as.numeric(as.vector(aw.edges)))
dev.off()





### FIGURE 6 
# Individual panels are set up as a function, so they can be arranged in a matrix below. 
# Each function has 1 option: Save F (plots the figure) or T (saves the figure)

### FIG. 6A: COMBINED HISTROGRAM OF GAINS AND LOSSES
Fig.6A<-function(save){
  
  # These lines make a data frame (df) with the type of transition in first column and the number of those transitions in the second column
  gains                  <- counts_bpp_Organ$gains                                                 # Vector with the number of gains of air-breathing
  losses                 <- counts_bpp_Organ$losses                                                # Vector with the number of losses of air-breathing
  M                      <- matrix(ncol = 2,nrow = (length(gains)+length(losses)))                 # Matrix (M) for data storage
  
  for (i in 1:length(gains)) { M[i,1]<-"WA"} 
  for (i in 1:length(losses)) { j=i+length(losses);  M[j,1]<-"AW"} 
  for (i in 1:length(gains)) { M[i,2]<-gains[i]} 
  for (i in 1:length(losses)) { j=i+length(losses);  M[j,2]<-losses[i]} 
  df                     <-as.data.frame(M)
  colnames(df)           <-c("transition","counts")
  df$counts              <-as.numeric(as.character(df$counts))
  
  # Use poisson test to calculate the most likely number of gains and losses
  aw.events              <-unname(poisson.test(round(mean(subset(df,transition=="AW")[,2])),T=1,r=0)$statistic) # The most likely number of losses
  wa.events              <-unname(poisson.test(round(mean(subset(df,transition=="WA")[,2])),T=1,r=0)$statistic) # The most likely number of gains
  
  # Plot
  ggplot(df,aes(x = counts,fill=transition))+
    geom_bar(aes(y = 100*(..count..)/sum(..count..)))+
    scale_x_continuous(limits = c(-1,1+max(counts_bpp_Organ[,1])))+
    scale_y_continuous(expand = c(0,0))+
    theme_classic()+
    theme(text = element_text(size=6),
          axis.text = element_text(size=6),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line = element_line(colour = 'black', size = 0.5),
          panel.grid.major = element_line(color = "grey", size = 0.1, linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          legend.title = element_text(face="bold"),
          legend.position = "none",
          legend.justification = c("center", "top"),
          legend.box.just = "center",
          legend.margin = margin(1,1,1,1),
          panel.background = element_rect(fill = cols[3]))+
    geom_vline(xintercept = aw.events, linetype="dotted", 
               color = "black", size=0.5)+
    geom_vline(xintercept = wa.events, linetype="dotted", 
               color = "black", size=0.5)+
    scale_fill_manual(name="Transition",
                      breaks=c("AW", "WA"),
                      labels=c("Air- to water-breathing", "Water- to air-breating"),
                      values = c("WA" = unname(cols[2]),
                                 "AW" = unname(cols[1])))+
    labs(y = expression("%"),x = expression(bold("Number of transitions")))->plot
  if(save==T) {
    ggsave(filename = "Fig. 6A.pdf",width = 6,height = 6,units = "cm",path = "./Figures/")
  }
  if(save==F){
    return(plot)
  }
}
Fig.6A(F)
Fig.6A(T)


## FIG. 6B: Timing of gains and losses
Fig.6B<-function(save){
  
  plot<-ggplot(data=times_bpp_Organ, aes(x=mode, y=time,fill=mode)) +
    
    geom_violin(colour = "black",size = 0.1,trim = F)+
    coord_cartesian(ylim=c(600,0))+
    theme_classic()+
    theme(text = element_text(size=6),
          axis.text.x = element_text(size=6),
          axis.line.x = element_line(size = 0.5),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y=element_blank(),
          legend.position = "none",
          panel.grid.major = element_line(color = "grey", size = 0.1, linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.background = element_rect(fill = cols[3]))+
    scale_fill_manual(name="Transition",
                      breaks=c("AW", "WA"),
                      labels=c("Air-breathing to Water-breathing","Water-breathing to Air-breathing"),
                      values = c("WA" = unname(cols[1]),
                                 "AW" = unname(cols[2])))+
    labs(y = expression(bold("Time before present (MYA)")))+
    scale_y_reverse()+
    coord_flip()
  plot  
  if(save==T) {
    ggsave(filename = "Fig. 6B.pdf",width = 6,height = 6,units = "cm",path = "./Figures/")
  }
  if(save==F){
    return(plot)
  }
}
Fig.6B(F)
Fig.6B(T)










## FIG. 6C: EVOLUTIONARY TRAJECTORY OF BREATHING MODE CONNECTING THE MRCA OF BONY FISHES AND AN EXTANT GENUS
# The number specify a genus in the phylogeny (see ge_tree$tip.label), where number 200 is Brienomyrus 

Fig.6C<-function(number,save){
  s<-number
  
  # Extract ancestral states for air-breathing
  ACE                    <- dSCM_mcc_Organ$ace
  ACE                    <- cbind(ACE[,1],1-ACE[,1])
  
  
  hor                    <-as.data.frame(tree_layout(ge_tree)$edgeDT)                              # Get time for each branch
  
  species_traj           <-ge_tree$tip.label[s]                                                    # The name of the genus 
  trajectory             <-Ancestors(ge_tree, match(species_traj,ge_tree$tip.label))               # Vector with the ancestors of that genus (i.e. internal node numbers)
  trajectory             <-append(x = trajectory,
                                  values = species_traj,
                                  after = length(trajectory))                                      # Add the extant genus to that vector
  l                      <-length(trajectory)                                                      # The length of the trajectory
  
  # Set up a matrix A, where each row specifies  
  # 1) node number
  # 2) the probability for air-breathing
  # 3) the geological time of that node
  # 4) the position in the trajectory
  
  A                      <-matrix(ncol = 4,nrow=l)
  colnames(A)            <-c("Node","Value","Time","Seq")
  
  A[,1]                  <-trajectory                                                              # node numbers in first column
  
  for (i in 1:l) {
    A[i,2]               <-1-unname(ACE[match(A[i,1],row.names(ACE))])                             # Find probability for air-breathing of that internal node
    A[i,3]               <-max(hor$xright)-hor$xright[match(trajectory[i],hor$V2)]                 # Extract the geological time of that internal node
    A[i,3]               <-ifelse(is.na(A[i,3])==TRUE, yes = 0, no = A[i,3])                       # For the extant genus, set time to zero
  }
  A[l,2]                 <-ifelse(unname(Organ[match(species_traj,names(Organ))])=="Gill",0,1)     # For the extant genus, set probability to 1 for air-breathers and to 0 for water breathers
  
  A                      <-as.data.frame(A)
  A$Value                <-as.numeric(as.character(A$Value));A$Time<-as.numeric(as.character(A$Time)) # Store as numeric values
  A                      <-A[order(A[,3]),]                                                        # Sort with respect to geological time
  A[,4]                  <- seq(from = 1, to = l, by = 1) 
  
  plot<-ggplot(A, aes(x = Time, y = Value,color=Value)) + 
    geom_line(color = "black",size = 0.1,linetype = "dashed")+
    scale_color_gradient(low = cols[1],high = cols[2])+
    theme_classic(10) + 
    theme(text=element_text(size=6), 
          axis.title = element_text(face = "bold"), 
          axis.line = element_blank(),
          axis.ticks = element_line(color = "grey",size = 0.1),
          legend.position="none",
          panel.grid.major = element_line(color = "grey", size = 0.1, linetype = "solid"),
          panel.background = element_rect(fill = cols[3]))+
    scale_x_reverse()+
    scale_y_continuous(limits = c(0,1))+
    geom_point(pch = 19, size = 1) + 
    labs(y = "Probability for air-breathing",
         x = "Time before present (MYA)")
  
  if(save==T) {
    ggsave(filename = "Fig. 6C.pdf",width = 6/2.54,height = 6/2.54,units = "in",path = "./Figures/")
  }
  if(save==F){
    return(plot)
  }
}

Fig.6C(200,T)
Fig.6C(200,F)







## FIGURE 7: Evolution of air-breathing in early fishes

# Build tree
ef_tree                  <- keep.tip(tree = ge_tree,tip = c("Latimeria","Protopterus","Polypterus","Acipenser","Amia","Anguilla"))

# Name the internal branches
ef_tree$node.label       <- c("Osteichthyes","Actinopterygii","Actinopteri","Neopterygii","Sarcopterygii")

# Set the Root node length
ef_tree$root.edge        <- 50

# Extract BPP for air-breathing from stochastic character mapping
ACE                      <- dSCM_mcc_Organ$ace
ACE                      <- cbind(ACE[,1],1-ACE[,1])

# Find the BPP for air-breathing in the corresponding internal branches in the pruned phylogeny
internal_pies          <- rbind(ACE[match(getMRCA(ge_tree,c("Latimeria","Protopterus")),row.names(ACE)),], # MRCA of sarcopterygeii (lobe finned fishes)
                                ACE[match(getMRCA(ge_tree,c("Latimeria","Anguilla")),row.names(ACE)),],    # MRCA of osteichthyes (bony fishes)
                                ACE[match(getMRCA(ge_tree,c("Polypterus","Anguilla")),row.names(ACE)),],   # MRCA of actinopterygii (ray finned fishes)
                                ACE[match(getMRCA(ge_tree,c("Acipenser","Anguilla")),row.names(ACE)),],    # MRCA of acipenseri
                                ACE[match(getMRCA(ge_tree,c("Amia","Anguilla")),row.names(ACE)),])         # MRCA of neopterygii
row.names(internal_pies) <- node.label


external_pies<-
  rbind(
    ACE[match(getMRCA(ge_tree,c("Erpetoichthys","Polypterus")),row.names(ACE)),],
    ACE[match(getMRCA(ge_tree,c("Psephurus", "Polyodon", "Scaphirhynchus", "Pseudoscaphirhynchus", "Huso", "Acipenser")),row.names(ACE)),],
    ACE[match(getMRCA(ge_tree,c("Amia", "Lepisosteus")),row.names(ACE)),],
    ACE[match(getMRCA(ge_tree,c("Elops","Perca")),row.names(ACE)),],
    ACE[match(getMRCA(ge_tree,c("Neoceratodus","Protopterus")),row.names(ACE)),],
    c(1,0)
  )

row.names(external_pies)<- 
  c("MRCA of bichirs and reedfishes",
    "MRCA of sturgeons and paddlefishes",
    "MRCA of bowfin and gars",
    "MRCA of teleosts",
    "MRCA of lungfishes",
    "Extant coelacanths")
external_pies


# Rename the tip labels in the pruned phylogeny
ef_tree$tip.label<-c("MRCA of bichirs and reedfishes",
                     "MRCA of sturgeons and paddlefishes",
                     "MRCA of bowfin and gars",
                     "MRCA of teleosts",
                     "MRCA of lungfishes",
                     "Extant coelacanths")

# Save figure 7 as PDF

pdf("./Figures/Figure 7 - early fishes.pdf",width = 14/2.54,height = 14/2.54,useDingbats = F)
par(mar=c(2,2,2,2))
plot(ef_tree,root.edge = T,type = "cladogram",use.edge.length = F,edge.width = 1,
     show.node.label = F,
     show.tip.label=F,
     #no.margin = T,
     label.offset = .5)
nodelabels(pie = internal_pies, piecol = cols[c(1,2)], cex = 1.7)
tiplabels(pie=external_pies,
          piecol=cols[c(1,2)],cex=1.7)
dev.off()





## FIG. S2: FREQUENCY DISTRIBUTION OF ALL TYPES OF TRANSITIONS BETWEEN AIR-BREATHING ORGANS
# Generate a function that plots each panel
Fig.S2.panel<-function(i){
  
  # Set color scheme for frequency distributions
  cols_fig2              <- c(rep(cols[2],7),rep(cols[1],30),cols[2],cols[1])
  
  j<-i
  
  df_plot<-as.data.frame(counts_bpp_Organ[,j])
  names(df_plot)<-"counts"
  events<-poisson.test(round(mean(counts_bpp_Organ[,i])),T=1,r=0)$statistic
  
  ggplot(df_plot,aes(x = counts))+
    geom_bar(col=cols_fig2[i],aes(y = 100*(..count..)/sum(..count..)))+
    scale_x_continuous(limits = c(-1,1+max(counts_bpp_Organ)))+
    scale_y_continuous(expand = c(0,0))+
    theme_classic()+
    geom_vline(xintercept = events, linetype="dotted", 
               color = "black", size=0.5)+
    theme(legend.position = "none",
          text = element_text(size=6),
          axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line = element_line(colour = 'black', size = 0.5),
          panel.grid.major = element_line(color = "grey", size = 0.1, linetype = "dashed"),
          panel.grid.major.y = element_blank(),
          panel.background = element_rect(fill = cols[3]))+
    annotate("text",size=2.1167, x = Inf, y = Inf, vjust = 1, hjust = 1,
             label = events)+
    labs(y = expression("%"),
         x = expression("Number of transitions"))
}



# Specify the plot
plot                   <- plot_grid(Fig.S2.panel(3),Fig.S2.panel(8),
                                    Fig.S2.panel(4),Fig.S2.panel(14),
                                    Fig.S2.panel(5),Fig.S2.panel(20),
                                    Fig.S2.panel(6),Fig.S2.panel(26),
                                    Fig.S2.panel(7),Fig.S2.panel(32),
                                    Fig.S2.panel(38),Fig.S2.panel(39),
                                    nrow = 6,ncol =2,align = 'hv',
                                    labels = c("GIT","GIT","Lung","Lung","Skin","Skin","Swimbladder","Swimbladder","Mouth","Mouth","Total gains","Total losses"),
                                    label_size = 6,hjust = 1,vjust = 1,label_x = 0.7,label_y = 0.95, axis = 'tb');plot

# Add x and y axis titles
y.grob                 <- textGrob("Frequency", gp=gpar(fontface="bold", fontsize=12), rot=90)
x.grob                 <- textGrob("Number of transitions", gp=gpar(fontface="bold", fontsize=12))

# Save af PDF
pdf("./Figures/Figure S2 - frequency distributions.pdf",width = 7.7/2.54,height = 7.7/2*6/2.5,useDingbats = F)
grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob))
dev.off()








###########################################################################
####                            SECTION 3                              ####
#### NUMBER AND TIMING OF TRANSITIONS BETWEEN AMPHIBIOUSNESS AND AQUA- #### 
#### TIC LIFELSTYLE USING 100 BAYESIAN POSTERIOR PROBABILITY TREES     ####
####                                                                   ####
####     THERE ARE TWO OPTIONS (A OR B): A: RUNNING THE MODEL FROM     #### 
####              SCRATCH OR B: LOADING PREVIOUSLY SAVED RUN           ####
###########################################################################



### OPTION A: RUN THE MODEL FROM SCRATCH
# Loop that runs stochatic character mapping on each of the 100 BPP trees. 
# Requires around 1.5 GB of available memory, and takes around 48 hours. 
# The output of the loop is the files "counts_bpp_LS.txt" and "times_bpp_LS.txt"
# that contains the transition counts and transition timing, respectively, from each loop.
nsim                     <- 100                                                                    # Number of simulations per BPP tree
model                    <- "ER"                                                                   # Set model for character evolution to equal rates model
counts_bpp_LS            <- matrix(ncol = 5)                                                       # Generate matrix to store counts of transitions
times_bpp_LS             <- matrix(ncol = 2);colnames(times_bpp_LS)<-c("mode","time")              # Generate matrix to store timing of transitions

for (k in 1:100) {
  ## STOCHASTIC CHARACTER MAPPING  
  SCM_bpp_LS                    <- make.simmap(trees[[k]],x = LS,model = model,nsim = nsim)        # Generate stochatic character maps
  dSCM_bpp_LS                   <- describe.simmap(SCM_bpp_LS)                                     # Summarize information from the stochastic character maps
  counts                 <- as.data.frame(countSimmap(SCM_bpp_LS,states=rownames(LS))$Tr)          # Calculate the number of transitions between states
  colnames(counts_bpp_LS)<-names(counts)
  counts_bpp_LS          <-na.omit(rbind(counts_bpp_LS,counts))                                    # Add the count data to the counts_bpp_LS matrix
  write.table(counts_bpp_LS,file = paste("./counts_bpp_LS.txt"))                                   # Save the number of counts
  
  ## FIND TIME OF TRANSITIONS
  hor <- as.data.frame(tree_layout(trees[[k]])$edgeDT)                                             # Extract divergence times from all internal nodes in the tree in the vector called "hor" (as in horizontal lines in the phylogeny)
  ACE <- cbind(dSCM_bpp_LS$ace[,1],1-dSCM_bpp_LS$ace[,1])                                          # Extract Bayesian posterior probabilities from all internal nodes in the tree
  
  # Generate a matrix H, that finds the timing of gains and losses of air-breathing
  H<-matrix(nrow = length(hor$V1), ncol = 7)                                                       # Generate a matrix for data storage, where each row contains an edge in the tree
  for (i in 1:length(H[,1])) { H[i,1] <- max(hor$xright)-hor$xleft[i] }                            # Edge goes from this age
  for (i in 1:length(H[,1])) { H[i,2] <- max(hor$xright)-hor$xright[i] }                           # to this age
  for (i in 1:length(H[,1])) { H[i,3] <- mean(H[i,1],H[i,2]) }                                     # where age is the midpoint of that edge
  for (i in 1:length(H[,1])) { H[i,4] <- ACE[match(hor$V1[i],row.names(ACE))]}                     # Bayesian posterior probability at first point of the edge
  for (i in 1:length(H[,1])) { H[i,5] <- ACE[match(hor$V2[i],row.names(ACE))]}                     # Bayesian posterior probability at second point of the edge
  for (i in 1:length(H[,1])) { H[i,5] <- ifelse(is.na(ACE[match(hor$V2[i],row.names(ACE))])==T,    # If the branch ends in an extant species, set Bayesian posterior probability to 1 for air-breathers and to 0 for water-breathers
                                                unname(ifelse(LS[match(trees[[k]]$tip[hor$V2[i]],names(LS))]=="Am",1,0)),
                                                H[i,5])} 
  for (i in 1:length(H[,1])) { H[i,6] <- ifelse(H[i,4]<0.5&H[i,5]>0.5,                             # If Bayesian posterior probability changes between <0.5 and >0.5, find the time where Bayesian posterior probability is 0.5 using linear regression (i.e. time of the origin of air-breathing). 
                                                (0.5-lm(c(H[i,4],H[i,5])~c(H[i,1],H[i,2]))$coefficients[1])/lm(c(H[i,4],H[i,5])~c(H[i,1],H[i,2]))$coefficients[2],
                                                NA)}                                              
  for (i in 1:length(H[,1])) { H[i,7] <- ifelse(H[i,4]>0.5&H[i,5]<0.5,                             # If Bayesian posterior probability changes between >0.5 and <0.5, find the time where Bayesian posterior probability is 0.5 using linear regression (i.e. time of the origin of water-breathing). 
                                                (0.5-lm(c(H[i,4],H[i,5])~c(H[i,1],H[i,2]))$coefficients[1])/lm(c(H[i,4],H[i,5])~c(H[i,1],H[i,2]))$coefficients[2],
                                                NA)}
  
  
  wa                     <- H[,6];   wa<-wa[!is.na(wa)]                                            # "wa" is a vector with all the times of transition from water to air
  aw                     <- H[,7];   aw<-aw[!is.na(aw)]                                            # "aw" is a vector with all the times of transition from air to water
  
  
  
  # Generate a dataframe (dataframe) with type of transition and their divergence times in first and second column
  M                      <- matrix(ncol = 2,nrow = (length(aw)+length(wa)))                        
  for (i in 1:length(wa)) { M[i,1]<-"WA"} 
  if (length(aw)>0){
    for (i in 1:length(aw)) { j=i+length(wa);  M[j,1]<-"AW"} 
  }
  
  for (i in 1:length(wa)) { M[i,2]<-wa[i]} 
  if (length(aw)>0){
    for (i in 1:length(aw)) { j=i+length(wa);  M[j,2]<-aw[i]}
  }
  dataframe              <- as.data.frame(M)
  colnames(dataframe)    <- c("mode","time")
  dataframe$time         <- as.numeric(as.character(dataframe$time))
  times_bpp_LS           <-na.omit(rbind(times_bpp_LS,dataframe))
  write.table(times_bpp_LS,file = paste("./times_bpp_LS.txt"))                                     # Save the timing of origins
}

### OPTION B: LOADING PREVIOUSLY SAVED RUN
# Load SCM counts and times. This can be done if option A loop has already been run. 
counts_bpp_LS            <-read.table("./counts_bpp_LS.txt")
times_bpp_LS             <-read.table("./times_bpp_LS.txt")




####################################################################
####                        SECTION 4                           ####
#### INDENTIFY CLADES WHERE AMPHIBIOUSNESS EVOLVED AND WAS LOST #### 
####         USING THE MAXIMUM CLADE CREDIBILITY TREE           ####
####################################################################



# Requires around 10 GB of available memory, and takes around 6 hours. 
nsim                     <- 1000                                                                    # Number of simulations 
model                    <- "ER"                                                                   # Set model for character evolution to equal rates model

SCM_mcc_LS               <- make.simmap(tree=ge_tree,x=LS,nsim=nsim,model=model)                   # Run stochastic character mapping on MCC tree
dSCM_mcc_LS              <- describe.simmap(SCM_mcc_LS)                                            # Summarize stochastic character mapping







##########################################################
####          PLOTS RELATED TO SECTION 3 & 4          ####
#### REQUIRES THAT SECTION 3-4 HAS BEEN RUN OR LOADED ####
##########################################################

head(dSCM_mcc_LS$ace)
### FIGURE 7+S3: CIRCULAR PHYLOGENY COLOR MAPPED WITH PROBABILITY FOR AMPHIBIOUSNESS
# Prepare for plotting: Find edges where air-breathing was gained and lost
hor                      <- as.data.frame(tree_layout(ge_tree)$edgeDT)                           # Extract divergence times from all internal nodes in the tree in the vector called "hor" (as in horizontal lines in the phylogeny)
ACE                      <- cbind(dSCM_mcc_LS$ace[,2],1-dSCM_mcc_LS$ace[,2])                      # Extract Bayesian posterior probabilities from all internal nodes in the tree

H<-matrix(nrow = length(hor$V1), ncol = 7)                                                       # Generate a matrix for data storage, where each row contains an edge in the tree
for (i in 1:length(H[,1])) { H[i,1] <- max(hor$xright)-hor$xleft[i] }                            # Edge goes from this age
for (i in 1:length(H[,1])) { H[i,2] <- max(hor$xright)-hor$xright[i] }                           # to this age
for (i in 1:length(H[,1])) { H[i,3] <- mean(H[i,1],H[i,2]) }                                     # where age is the midpoint of that edge
for (i in 1:length(H[,1])) { H[i,4] <- ACE[match(hor$V1[i],row.names(ACE))]}                     # Bayesian posterior probability at first point of the edge
for (i in 1:length(H[,1])) { H[i,5] <- ACE[match(hor$V2[i],row.names(ACE))]}                     # Bayesian posterior probability at second point of the edge
for (i in 1:length(H[,1])) { H[i,5] <- ifelse(is.na(ACE[match(hor$V2[i],row.names(ACE))])==T,    # If the edge ends in an extant species, set Bayesian posterior probability to 1 for amphibious and to 0 for aquatic
                                              unname(ifelse(LS[match(ge_tree$tip[hor$V2[i]],names(LS))]=="Aq",1,0)),
                                              H[i,5])} 
for (i in 1:length(H[,1])) { H[i,6] <- ifelse(H[i,4]<0.5&H[i,5]>0.5,                             # If Bayesian posterior probability changes between <0.5 and >0.5, find the time where Bayesian posterior probability is 0.5 (i.e. time of the origin of amphibiousness). 
                                              (0.5-lm(c(H[i,4],H[i,5])~c(H[i,1],H[i,2]))$coefficients[1])/lm(c(H[i,4],H[i,5])~c(H[i,1],H[i,2]))$coefficients[2],
                                              NA)}                                              
for (i in 1:length(H[,1])) { H[i,7] <- ifelse(H[i,4]>0.5&H[i,5]<0.5,                             # If Bayesian posterior probability changes between >0.5 and <0.5, find the time where Bayesian posterior probability is 0.5 (i.e. time of the origin of aquatic). 
                                              (0.5-lm(c(H[i,4],H[i,5])~c(H[i,1],H[i,2]))$coefficients[1])/lm(c(H[i,4],H[i,5])~c(H[i,1],H[i,2]))$coefficients[2],
                                              NA)}



# Find the branches in the tree, where amphibiousness originated (wa.edges) and was lost (aw.edges)
m                      <- cbind(hor$V1[as.numeric(row.names((subset(as.data.frame(H),is.na(H[,6])==F))))],
                                hor$V2[as.numeric(row.names((subset(as.data.frame(H),is.na(H[,6])==F))))])
l                      <- cbind(hor$V1[as.numeric(row.names((subset(as.data.frame(H),is.na(H[,7])==F))))],
                                hor$V2[as.numeric(row.names((subset(as.data.frame(H),is.na(H[,7])==F))))])
g                      <- as.data.frame(ge_tree$edge)

aw.edges<-matrix(ncol=1,nrow=length(m[,2]))
for (i in 1:length(m[,1])) {
  aw.edges[i]<-rownames(subset(g,g$V1==m[i,1]&g$V2==m[i,2]))
}

if (length(l>0)){
  wa.edges<-matrix(ncol=1,nrow=length(l[,2]))
  for (i in 1:length(l[,1])) {
    wa.edges[i]<-rownames(subset(g,g$V1==l[i,1]&g$V2==l[i,2]))
  }
}


# Save a density plot
obj_3<-densityMap(SCM_mcc_LS,res=1000)
obj_3$cols[1:1001]<-colorRampPalette(cols[2:1], space="Lab")(1001)

# Fake plot used to (somehow) active the axis function
plot(obj_3,lwd = 0.5,fsize = 0.1,ftype = "i",show.tip.label = F,type = "fan",legend=F)
axis_obj<-axis(1,pos=0,at=seq(max(nodeHeights(obj$tree)),100,by=-100),cex.axis=0.5,labels=FALSE)

# Plot Fig. 9
pdf(file = "./Figures/Figure 9 - amphibious evolution.pdf",width = 19/2.54,height = 19/2.54,useDingbats = F)                              # PDF saving options
plot(obj_3,lwd = 1,fsize = .000001,ftype = "i",show.tip.label = F,type = "fan",legend=F,no.margin=T)# Plot density map
for(i in 1:length(axis_obj)){                                                                    # Plot circular time axis
  a1<-0
  a2<-2*3.141528
  draw.arc(0,0,radius=axis_obj[i]-100,angle1 = a1,angle2 = a2,lwd=0.5,lty="dashed",col=make.transparent("black",0.5))}
if (length(l>0)){edgelabels(pch=21,bg=cols[2],cex=2,edge=as.numeric(as.vector(wa.edges)))}                       # Add red dot on edges where amphibiousness evolved
edgelabels(pch=21,bg=cols[1],cex=2,edge=as.numeric(as.vector(aw.edges)))                       # Add red blue on edges where amphibiousness was lost
dev.off()



### Fig. S3
# Large sized density map with genus, family and order names (Fig. S1)
obj_4<-obj_3                                                                                      # Copy the first density map
obj_4$tree$tip.label     <- paste(obj_4$tree$tip.label," (",
                                  as.character(tax[match(obj_4$tree$tip.label,tax$genus),1]),", ",
                                  as.character(tax[match(obj_4$tree$tip.label,tax$genus),2]),")",
                                  sep = "")                                                      # Change tip labels to include genus, family and order

# Plot and save Fig. S3
pdf(file = paste("./Figures/Figure S3 - amphbious evolution",sep = ""),width = 50,height = 50,useDingbats = F)
plot(obj_4,lwd = 1.5,fsize = .2,show.tip.label = F,type = "fan",legend=F,no.margin=T)
for(i in 1:length(axis_obj)){
  a1<-0
  a2<-2*3.141528
  draw.arc(0,0,radius=axis_obj[i]-100,angle1 = a1,angle2 = a2,lwd=0.1,lty="dashed",
           col=make.transparent("black",0.5))  
}
edgelabels(pch=21,bg=cols[2],cex=2,edge=as.numeric(as.vector(wa.edges)))
dev.off()






  
  
######################################################
#### ANSWERS TO QUESTIONS ASKED IN THE MANUSCRIPT ####
######################################################  

# How many percent of air-breathing origins occured within the last 65 MYA?
100*sum(ifelse(subset(times_bpp_Organ,mode =="AW")$time<65,1,0))/length(subset(times_bpp_Organ,mode =="AW")$time)

# How many percent of amphibious origins occured within the last 65 MYA?
100*sum(ifelse(subset(times_bpp_LS,mode =="WA")$time<65,1,0))/length(subset(times_bpp_LS,mode =="WA")$time)

## Which transitions among air-breathing organs significantly different from zero. TRUE for not different from zero. 
W                      <- cbind(colnames(counts_bpp_Organ),rep(NA,39),rep(NA,39))                        # Empty data frame for data storage
for (i in 1:39) {W[i,2]<- poisson.test(round(mean(counts_bpp_Organ[,i])),T=1,r=0)$p.value
                 W[i,3]<- poisson.test(round(mean(counts_bpp_Organ[,i])),T=1,r=0)$statistic}             # Test the hypothesis that counts are zero
as.data.frame(W) 



## Are gains and losses of amphibiousness significantly different from zero. TRUE for not different from zero. 
X                      <- cbind(colnames(counts_bpp_LS),rep(NA,5),rep(NA,5))                        # Empty data frame for data storage
for (i in 1:5) {X[i,2]<- poisson.test(round(mean(counts_bpp_LS[,i])),T=1,r=0)$p.value
X[i,3]<- poisson.test(round(mean(counts_bpp_LS[,i])),T=1,r=0)$statistic}             # Test the hypothesis that counts are zero
as.data.frame(X) 




# Is the number of origins of air-breathing higher than 67?
poisson.test(round(mean(counts_bpp_Organ$gains)),T=1,r=67)
# Comment: Yes!

# Are the number of air-breathing orgins lower than amphibious origins
t.test(counts_bpp_Organ$gains,counts_bpp_LS$Aq.Am)
# Air-breathing evolved fewer times that amphibious behaviour. 


# What is the probability for water- and air-breathing in the MRCA of bony fishes?
node<-which(rownames(dSCM_mcc_Organ$ace)==getMRCA(ge_tree,c("Neoceratodus","Amia")))
ACE                      <- cbind(dSCM_mcc_Organ$ace[node,1],1-dSCM_mcc_Organ$ace[node,1])                      # Extract Bayesian posterior probabilities from all internal nodes in the tree
ACE

# What is the probability for Amphibiousness in the MRCA of bony fishes?
ACE                      <- cbind(dSCM_mcc_LS$ace[,1],1-SCM_mcc_LS$ace[,1])                      # Extract Bayesian posterior probabilities from all internal nodes in the tree
ACE

node<-which(rownames(dSCM_mcc_LS$ace)==getMRCA(ge_tree,c("Neoceratodus","Amia")))
ACE                      <- cbind(dSCM_mcc_LS$ace[node,1],1-dSCM_mcc_LS$ace[node,1])                      # Extract Bayesian posterior probabilities from all internal nodes in the tree
ACE

# What is the probability for the gastrointestinal tract being the air-breathing organ in the most recent common ancestor of Loricariidae?
node<-which(rownames(dSCM_mcc_Organ$ace)==getMRCA(ge_tree,c("Corydoras","Microlepidogaster")))
ACE                      <- cbind(dSCM_mcc_Organ$ace[node,2],1-dSCM_mcc_Organ$ace[node,2])                      # Extract Bayesian posterior probabilities from all internal nodes in the tree
ACE
# 83.1% probability that the GIT was the ABO of the MRCA of loricariidae


###################################################
#### LINKING SPECIATION RATE TO BREATHING MODE ####
####               BiSSE model                 ####
###################################################
ge_tree<-nnls.tree(cophenetic(ge_tree),ge_tree,rooted=TRUE) # Make tree ultrametric

lik <- make.bisse(nnls.tree(cophenetic(ge_tree),ge_tree,rooted=TRUE), 
                  AW_bi)
p <- starting.point.bisse(ge_tree)
fit <- find.mle(lik, p)
round(coef(fit), 3)

lik.l <- constrain(lik, lambda1 ~ lambda0)                                                         # Test if lambda is equal for water- and air-breathers
fit.l <- find.mle(lik.l, p[argnames(lik.l)])
rbind(full=coef(fit), equal.l=coef(fit.l, TRUE))
anova(fit, equal.l=fit.l)
# Conclusion: The are different

lik.l <- constrain(lik, mu1 ~ mu0)                                                                 # Test if mu is equal for water- and air-breathers
fit.l <- find.mle(lik.l, p[argnames(lik.l)])
rbind(full=coef(fit), equal.l=coef(fit.l, TRUE))
anova(fit, equal.l=fit.l)
# Conclusion: They are not different

lik.l <- constrain(lik, q01 ~ q10)                                                                 # Test if forward transitions rate is equal to back transitions rate
fit.l <- find.mle(lik.l, p[argnames(lik.l)])
rbind(full=coef(fit), equal.l=coef(fit.l, TRUE))
anova(fit, equal.l=fit.l)
# Conclusion: They are different



##### 
# Other plots to the paper #
######
cols<-c(rgb(red=79, green=129, blue=189, maxColorValue = 255),
        rgb(red=180, green=42, blue=81, maxColorValue = 255),
        rgb(red=254, green=254, blue=218, maxColorValue = 255),
        rgb(red=127, green=127, blue=127, maxColorValue = 255),
        rgb(red=30 , green=200, blue=110, maxColorValue = 255),
        rgb(red=250, green=180, blue=140, maxColorValue = 255),
        rgb(red=000, green=000, blue=000, maxColorValue = 255),
        rgb(red=255, green=255, blue=255, maxColorValue = 255))

df1<-data.frame(PO2 = c(20,17,14,20,17,14),
               Onset = c(34,39,40.1,38.1,37.8,33.5),
               Species = c(rep("Trichopodus trichopterus",3),rep("Betta splendens",3)))

ggplot(df1,aes(x = PO2,y=Onset,col=Species))+
  scale_color_manual(values = cols[1:2])+
  geom_path()+
  theme_classic(10) + 
  theme(text=element_text(size=12), 
        axis.title = element_text(face = "bold"), 
        axis.line = element_blank(),
        axis.ticks = element_line(color = "grey",size = 0.1),
        legend.position="none",
        legend.text = element_text(face="italic"),
        legend.title = element_blank(),
        panel.grid.major = element_line(color = "grey", size = 0.1, linetype = "solid"),
        panel.background = element_rect(fill = cols[3]))+
  geom_point(pch = 19, size = 3) + 
  scale_y_continuous(limits = c(33,41))+
  scale_x_continuous(breaks = c(14,17,20))+
  annotate(geom = "text",x = 15,y = 34,label = "Betta splendens",fontface = 'italic',hjust = 0,size=12/3)+
  annotate(geom = "text",x = 15,y = 40.2,label = "Trichopodus trichopterus",fontface = 'italic',hjust = 0,size=12/3)+
  labs(y = expression("Onset of air-breathing (days)"),
       x = expression("Rearing PO"[2]*" (kPa)"))->p1;p1
ggsave(paste("./Figures/Figure 4.pdf",sep = ""),height = 7.7/2.54,width = 7.7/2.54)




df2<-data.frame(lsa=c(24.70,24.89,13.80,5.30,5.41,5.34),
                lsaerr=c(3.3,0.52,0.17,0.87,0.89,1.25),
                abo=c(0.25,0.39,0.28,0.06,0.06,0.1),
                aboerr=c(0.06,0.07,0.04,0.01,0.01,0.03),
                PO2 = c(20,17,14,20,17,14),
                Species = c(rep("Betta splendens",3),rep("Trichopodus trichopterus",3)))
labels<-c("A","A","B",rep("a",3))

ggplot(df2,aes(x = PO2,y=lsa,col=Species))+
  scale_color_manual(values = cols[1:2])+
  geom_path()+
  theme_classic(10) + 
  geom_text(aes(label=labels), vjust=-0.5, color="black",
            y=as.numeric(df2$lsa+df2$lsaerr),
            size=4)+
  theme(text=element_text(size=12), 
        axis.title = element_text(face = "bold"), 
        axis.line = element_blank(),
        axis.ticks = element_line(color = "grey",size = 0.1),
        legend.position="none",
        legend.text = element_text(face="italic"),
        legend.title = element_blank(),
        panel.grid.major = element_line(color = "grey", size = 0.1, linetype = "solid"),
        panel.background = element_rect(fill = cols[3]))+
  geom_point(pch = 19, size = 3)+
  geom_errorbar(df2,mapping=aes(x=PO2,ymin=lsa-lsaerr,ymax=lsa+lsaerr),width = 0.1)+
  scale_y_continuous(limits = c(0,30))+
  scale_x_continuous(breaks = c(14,17,20))+
  labs(y = expression("Total lamellar surface area (mm"^"2"*")"),
       x = expression("Rearing PO"[2]*" (kPa)"))->p2;p2
ggsave(paste("./Figures/Fig. 3A.pdf",sep = ""),height = 7.7/2.54,width = 7.7/2.54)


labels<-c("A","B","A",rep("a",3))
ggplot(df2,aes(x = PO2,y=abo,col=Species))+
  scale_color_manual(values = cols[1:2])+
  geom_path()+
  geom_text(aes(label=labels), vjust=-0.5, color="black",
            y=as.numeric(df2$abo+df2$aboerr),
            size=4)+
  theme_classic(10) + 
  theme(text=element_text(size=12), 
        axis.title = element_text(face = "bold"), 
        axis.line = element_blank(),
        axis.ticks = element_line(color = "grey",size = 0.1),
        legend.position="none",
        legend.text = element_text(face="italic"),
        legend.title = element_blank(),
        panel.grid.major = element_line(color = "grey", size = 0.1, linetype = "solid"),
        panel.background = element_rect(fill = cols[3]))+
  geom_point(pch = 19, size = 3)+
  geom_errorbar(df2,mapping=aes(x=PO2,ymin=abo-aboerr,ymax=abo+aboerr),width = 0.1)+
  #scale_y_continuous(limits = c(0,30))+
  scale_x_continuous(breaks = c(14,17,20))+
  labs(y = expression("Total labyrinth surface area (mm"^"2"*")"),
       x = expression("Rearing PO"[2]*" (kPa)"))->p3;p3
ggsave(paste("./Figures/Fig. 3B.pdf",sep = ""),height = 7.7/2.54,width = 7.7/2.54)



# Figure panel for Fig. 10 (makes 50 different plots - chose the one that looks best)
for (i in 1:50){
  random.tree<-
    pbtree(n=10)
  char<-setNames(object = c("B","B","A","A","A","A","A","A","A","A"),nm = random.tree$tip.label)
  
  simmap<-make.simmap(random.tree,char,nsim = 100,model = "ER")
  
  obj_5<-densityMap(simmap,legend = F)
  obj_5$cols[1:1001]         <-colorRampPalette(cols[1:2], space="Lab")(1001)                          # Change color theme to Acta colors
  
  pdf(paste("./Figures/Fig10/Figure 10_",i,".pdf"),width = 7/2.54,height = 7/2.54,useDingbats = F)
  plot(obj_5,lwd = 3,fsize = .000001,ftype = "i",show.tip.label = F,legend=F,no.margin=T,use.edge.lengths = F)# Plot density map
  dev.off()
}

