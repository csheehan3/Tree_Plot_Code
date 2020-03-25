######Get these packages or load them 
library("kdensity")
library(devtools)
library(tidyverse)
library("ggpubr")
library("readr")
library("biomaRt")
library("rlist")
######need to get rid of the genes with average value 0
BRCA_data <- read_tsv(file = "TCGA-BRCA_fpkm.tsv")
row_end <- dim(BRCA_data)[1]
col_end <- dim(BRCA_data)[2]
new_BRCA_data <- apply(BRCA_data[, 2:col_end],1,mean) %>% cbind(average_values, BRCA_data) #had to get the mean of each row for next step
cleaned_data <- filter(new_BRCA_data, new_BRCA_data$average_values!=0) #eliminate mean of zeros from the dataset
row_end <- dim(cleaned_data)[1]
col_end <- dim(cleaned_data)[2]
######
set.seed(10)
##### run 500 simulations of mean of 10 random genes correlated with another random gene
mean_spearman_values <- c()
simulation_count=0
while (simulation_count<501) {
    A_very_interesting_gene_index <- sample(1:row_end, 1, replace=TRUE)
    A_very_interesting_gene <- as.numeric(cleaned_data[A_very_interesting_gene_index, 3:col_end])
    row_indices <- sample(1:row_end, 10, replace=TRUE)
    random_rows <- cleaned_data[row_indices, 3:col_end]
    spearman_list <- c()
    for (A_gene in 1:10){ ##do the correlations for the 10 randomly chosen genes 
      A_correlation <- cor(as.numeric(random_rows[A_gene,]), A_very_interesting_gene, method = c("spearman"))
      spearman_list <- c(spearman_list, A_correlation)
      if(is.na(A_correlation)==TRUE){
        print("warning!")
        print(head(as.numeric(random_rows[A_gene,])))
        print(head(A_very_interesting_gene))
        break()
      }
    }
    mean_spearman_values <- c(mean_spearman_values, mean(spearman_list))
    simulation_count <- list.filter(mean_spearman_values, is.na(.)==FALSE) %>% length()
    print(paste("simulation count is equal to:", simulation_count))
}
Ten_gene_set_simulations <- kdensity(mean_spearman_values, kernel=c("gaussian"))
shapiro.test(mean_spearman_values)
####### This is now repeated to look at mean correlations with a set of 50 random genes
set.seed(10)
##
mean_spearman_values <- c()
simulation_count=0
while (simulation_count<501) {
  A_very_interesting_gene_index <- sample(1:row_end, 1, replace=TRUE)
  A_very_interesting_gene <- as.numeric(cleaned_data[A_very_interesting_gene_index, 3:col_end])
  row_indices <- sample(1:row_end, 50, replace=TRUE)
  random_rows <- cleaned_data[row_indices, 3:col_end]
  spearman_list <- c()
  for (A_gene in 1:50){
    A_correlation <- cor(as.numeric(random_rows[A_gene,]), A_very_interesting_gene, method = c("spearman"))
    spearman_list <- c(spearman_list, A_correlation)
    if(is.na(A_correlation)==TRUE){
      print("warning!")
      print(head(as.numeric(random_rows[A_gene,])))
      print(head(A_very_interesting_gene))
      spearman_list <- c()
      break()
    }
  }
  mean_spearman_values <- c(mean_spearman_values, mean(spearman_list))
  simulation_count <- list.filter(mean_spearman_values, is.na(.)==FALSE) %>% length()
  print(paste("simulation count is equal to:", simulation_count))
}
Fifty_gene_set_simulations <- kdensity(mean_spearman_values, kernel=c("gaussian"))
Fifty_gene_set_simulations(0.05)
shapiro.test(mean_spearman_values)
###### now for 100 genes
set.seed(10)
##
mean_spearman_values <- c()
simulation_count=0
while (simulation_count<501) {
  A_very_interesting_gene_index <- sample(1:row_end, 1, replace=TRUE)
  A_very_interesting_gene <- as.numeric(cleaned_data[A_very_interesting_gene_index, 3:col_end])
  row_indices <- sample(1:row_end, 100, replace=TRUE)
  random_rows <- cleaned_data[row_indices, 3:col_end]
  spearman_list <- c()
  for (A_gene in 1:100){
    A_correlation <- cor(as.numeric(random_rows[A_gene,]), A_very_interesting_gene, method = c("spearman"))
    spearman_list <- c(spearman_list, A_correlation)
    if(is.na(A_correlation)==TRUE){
      print("warning!")
      print(head(as.numeric(random_rows[A_gene,])))
      print(head(A_very_interesting_gene))
      spearman_list <- c()
      break()
    }
  }
  mean_spearman_values <- c(mean_spearman_values, mean(spearman_list))
  simulation_count <- list.filter(mean_spearman_values, is.na(.)==FALSE) %>% length()
  print(paste("simulation count is equal to:", simulation_count))
}
Hundred_gene_set_simulations <- kdensity(mean_spearman_values, kernel=c("gaussian"))
Hundred_gene_set_simulations(0.05)
shapiro.test(mean_spearman_values)
#####Again for 200 hundred genes
set.seed(10)
##
mean_spearman_values <- c()
simulation_count=0
while (simulation_count<501) {
  A_very_interesting_gene_index <- sample(1:row_end, 1, replace=TRUE)
  A_very_interesting_gene <- as.numeric(cleaned_data[A_very_interesting_gene_index, 3:col_end])
  row_indices <- sample(1:row_end, 200, replace=TRUE)
  random_rows <- cleaned_data[row_indices, 3:col_end]
  spearman_list <- c()
  for (A_gene in 1:200){
    A_correlation <- cor(as.numeric(random_rows[A_gene,]), A_very_interesting_gene, method = c("spearman"))
    spearman_list <- c(spearman_list, A_correlation)
    if(is.na(A_correlation)==TRUE){
      print("warning!")
      print(head(as.numeric(random_rows[A_gene,])))
      print(head(A_very_interesting_gene))
      spearman_list <- c()
      break()
    }
  }
  mean_spearman_values <- c(mean_spearman_values, mean(spearman_list))
  simulation_count <- list.filter(mean_spearman_values, is.na(.)==FALSE) %>% length()
  print(paste("simulation count is equal to:", simulation_count))
}
simulation_density <- kdensity(mean_spearman_values, kernel=c("gaussian"))
########Now we can plot the density functions created from each of the 4 types of simulations
plot(Ten_gene_set_simulations, col="gold3", ylim=c(0,20), xlim=c(-0.15,0.15), xlab="Mean Spearman Coefficient", main=NULL) 
lines(Fifty_gene_set_simulations, col="orange")
lines(Hundred_gene_set_simulations, col="red")
lines(simulation_density, col="purple")
legend(-0.1, 15, legend=c("10 genes", "50 genes", "100 genes", "200 genes"),
       col=c("gold3", "orange", "red", "purple"), lty=1:2, cex=1.1)
#######Here I add the Entrez gene name conversions to the BRCA dataset
gene_list <- c(cleaned_data[,2])
library("biomaRt")
cleaved_gene_names <- str_replace(gene_list, #Ensemble_ID gene names are modified to a more general ID name
                                  pattern = ".[0-9]+$",
                                  replacement = "")
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gene_conversion <- getBM(attributes = c('ensembl_gene_id', #query sent to biomaRt and desired Entrez names return
                                          'entrezgene_id'),
                           filters = 'ensembl_gene_id', 
                           values = cleaved_gene_names,
                           mart = mart)
n=0
entrez_dictionary <- c()
for (translat in cleaned_data$Ensembl_ID[1:row_end]){
  cleaved_translat <- str_replace(translat, pattern = ".[0-9]+$", replacement = "")
  newline <- filter(gene_conversion, gene_conversion$ensembl_gene_id==cleaved_translat)
  n <- n + 1 #####this is just a counter to let you know how its running
  if (n %% 100 == 0) {
    print(paste("Working on line:",n," of ",row_end)) #reporter
  }
  if (nrow(newline) > 1){
    newline <- newline[1, 1:2] ###picks the first for cases with multiple gene transcripts
  }
  if (length(newline$entrezgene_id)==0 || is.na(newline$entrezgene_id)==TRUE) { ###puts in NA for cases that missed a conversion
    entrez_dictionary <- c(entrez_dictionary, NA)
  } else {
    entrez_dictionary <- c(entrez_dictionary, newline$entrezgene_id) ####grabs the Entrez conversion for all other cases
  }
}
BRCA_genes_altered <- cbind(as.data.frame(entrez_dictionary), cleaned_data)
########Here I make a Tibble with the GO terms I'm interested in & grab all the Entrez gene names included in each set
GO_cellular_metabolic_process <- read_tsv(file="basket.tsv", col_names = FALSE)
GO_metabolism_additional_terms <- read_tsv(file="basket (1).tsv", col_names=FALSE)
GO_cellular_metabolic_process <- bind_rows(GO_cellular_metabolic_process, GO_metabolism_additional_terms)
GO_query_list <- c()
for (term in 1:142){
  GO_query_list <- c(GO_query_list, as.character(GO_cellular_metabolic_process[term,1]))
}
GO_metabolic_list <- c()
n=1
for (term in 1:142){
  GO_conversion <- getBM(attributes = c('hgnc_symbol',
                                        'entrezgene_id'),
                         filters = "go_parent_term", 
                         values = GO_query_list[term],
                         mart = mart)
  print(paste("completed:", GO_query_list[term]))
  if (length(GO_conversion$entrezgene_id) > 50 && length(GO_conversion$entrezgene_id) < 300){ ##I only care about sets between 50 and 300 genes in length
    GO_metabolic_list[[n]] <-  as.list(GO_conversion$entrezgene_id)
    names(GO_metabolic_list)[[n]] <- GO_query_list[term]
    n=n+1
  }
}
#####calculate the mean of the multiple correlations tests
BACH1_row <- filter(BRCA_genes_altered, BRCA_genes_altered$entrez_dictionary==571)[,4:1220] %>% as.numeric()
GO_spearman_values <- c()
mean_GO_spearman_values <- c()
for (term in 1:length(GO_metabolic_list)) {
  filtered_frame <- filter(BRCA_genes_altered, BRCA_genes_altered$entrez_dictionary %in% GO_metabolic_list[[term]])
  length_of_filter <- dim(filtered_frame)[1]
  for (n in 1:length_of_filter){
    test_row <- filtered_frame[n,4:1220] %>% as.numeric()
    correlation_value <- cor(BACH1_row, test_row, method=c("spearman"))
    GO_spearman_values <- c(GO_spearman_values, correlation_value)
  }
  mean_GO_spearman_values <- c(mean_GO_spearman_values, mean(GO_spearman_values))
  GO_spearman_values <- c() ##reset the GO_spearman_values list for the next GO term to fill up
  print(paste("completed:", names(GO_metabolic_list)[term])) #reporter
}
GO_test_tibble <- tibble("Metabolic_GO_Term" = names(GO_metabolic_list), "Mean_Spearman_Value" = mean_GO_spearman_values)
arrange(GO_test_tibble, desc(GO_test_tibble$Mean_Spearman_Value))
######now for a probability value,
######first I need to find spearman value where integration to +Inf is = 0.5
starting_value=0.03
while (integrate(function(x) simulation_density(x), lower=starting_value, upper=Inf)$value > 0.5){ #increases value until integration is equal to 0.5
  starting_value = starting_value + 0.00001
}
Probability_function = function(X1){ #Function that passes value to integration function depending on whether it is above or below the midway point
  if (X1 > starting_value){
    integrate(function(x) simulation_density(x), lower=X1, upper=Inf)$value
  } else {
    integrate(function(x) simulation_density(x), lower=-Inf, upper=X1)$value
  }
}
GO_test_tibble <- sapply(GO_test_tibble$Mean_Spearman_Value, Probability_function) %>% ##Gets the probability from each mean and adds to tibble
  cbind(GO_test_tibble, "Probability_Value"=.) 
Bunch_of_zeros <- GO_test_tibble$Mean_Spearman_Value * 0
####### This was for plotting the results with the GO term analysis and probability distribution
kdensity(mean_spearman_values, kernel=c("gaussian")) %>% 
  plot(main="GO Terms overlayed with Probability Distribution", xlab="Mean Spearman Value", xlim=c(-0.25,0.25), col="blue")
plot(x=GO_test_tibble$Mean_Spearman_Value, y=Bunch_of_zeros, col="red", xlim=c(-0.25,0.25))
####### I wanted to save the actual values to look at which terms were hits 
arrange(GO_test_tibble, desc(GO_test_tibble$Mean_Spearman_Value)) %>% 
  write.table(., file = "Mean_Spearman_GO_Table.csv")
######Performs the correlations again for the interesting terms to generate a tree plot
######Collects the spearman coefficient for each gene against bach1, the p-value, and translate Entrez number to gene name
BACH1_row <- filter(BRCA_genes_altered, BRCA_genes_altered$entrez_dictionary==571)[,4:1220] %>% as.numeric()
filtered_frame <- filter(BRCA_genes_altered, BRCA_genes_altered$entrez_dictionary %in% GO_metabolic_list[[28]])
length_of_filter <- dim(filtered_frame)[1]
GO_single_spearman_values <- c()
GO_single_p_values <- c()
for (n in 1:length_of_filter){
    single_gene_row <- filtered_frame[n,4:1220] %>% as.numeric()
    correlation_value <- cor.test(BACH1_row, single_gene_row, method=c("spearman"))
    correlation_value$estimate
    GO_single_spearman_values <- c(GO_single_spearman_values, correlation_value$estimate)
    GO_single_p_values <- c(GO_single_p_values, correlation_value$p.value)
}
hgnc_conversion <- getBM(attributes = c('entrezgene_id', 
                                        'hgnc_symbol'),
                         filters = 'entrezgene_id', 
                         values = filtered_frame$entrez_dictionary,
                         mart = mart)
Another_Random <- tibble("Gene_Name"=hgnc_conversion$hgnc_symbol, "Spearman_Coefficient"=GO_single_spearman_values, "P-Value"=GO_single_p_values)
Andother_Random <- arrange(Another_Random, Another_Random$`P-Value`)
Random_Tibble <- tibble("Gene_Name"=hgnc_conversion$hgnc_symbol, "Spearman_Coefficient"=GO_single_spearman_values, "P-Value"=GO_single_p_values)
Random_Tibble <- arrange(Random_Tibble, Random_Tibble$`P-Value`)
ROS_Tibble <- tibble("Gene_Name"=hgnc_conversion$hgnc_symbol, "Spearman_Coefficient"=GO_single_spearman_values, "P-Value"=GO_single_p_values)
ROS_Tibble <- arrange(ROS_Tibble, ROS_Tibble$`P-Value`)
head(ROS_Tibble)
Carbs_Tibble <- tibble("Gene_Name"=hgnc_conversion$hgnc_symbol, "Spearman_Coefficient"=GO_single_spearman_values, "P-Value"=GO_single_p_values)
Carbs_Tibble <- arrange(Carbs_Tibble, Carbs_Tibble$`P-Value`)
head(Carbs_Tibble)
OXPHOS_tibble <- tibble("Gene_Name"=hgnc_conversion$hgnc_symbol, "Spearman_Coefficient"=GO_single_spearman_values, "P-Value"=GO_single_p_values)
OXPHOS_tibble <- arrange(OXPHOS_tibble, OXPHOS_tibble$`P-Value`)
head(OXPHOS_tibble)
####Generate the Tree Plots for each of the interesting terms
ggplot(OXPHOS_tibble, aes(x=reorder(OXPHOS_tibble$Gene_Name, OXPHOS_tibble$`P-Value`), y=OXPHOS_tibble$Spearman_Coefficient)) + 
  geom_bar(stat="identity" , width=0.3, col="red") +
  theme_light() +
  theme(axis.text=element_text(size=3)) +
  labs(title="BACH1 correlations with OXPHOS Gene Set", y="Spearman Coefficient", x="P-Value")

ggplot(Carbs_Tibble, aes(x=reorder(Carbs_Tibble$Gene_Name, Carbs_Tibble$`P-Value`), y=Carbs_Tibble$Spearman_Coefficient)) + 
  geom_bar(stat="identity" , width=0.3, col="red") +
  theme_light() +
  theme(axis.text=element_text(size=3)) +
  labs(title="BACH1 correlations with Positive Regulation of Carbohydrate Metabolic Process Gene Set", y="Spearman Coefficient", x="P-Value")

ggplot(ROS_Tibble, aes(x=reorder(ROS_Tibble$Gene_Name, ROS_Tibble$`P-Value`), y=ROS_Tibble$Spearman_Coefficient)) +
  geom_bar(stat="identity" , width=0.3, col="red") +
  theme_light() +
  theme(axis.text=element_text(size=3)) +
  labs(title="BACH1 correlations with Positive Regulation of ROS Metabolic Process Gene Set", y="Spearman Coefficient", x="P-Value")

ggplot(Random_Tibble, aes(x=reorder(Random_Tibble$Gene_Name, Random_Tibble$`P-Value`), y=Random_Tibble$Spearman_Coefficient)) +
  geom_bar(stat="identity" , width=0.3, col="red") +
  theme_light() +
  theme(axis.text=element_text(size=3)) +
  labs(title="BACH1 correlations with Membrane Biosynthesis Gene Set", y="Spearman Coefficient", x="P-Value")

ggplot(Another_Random, aes(x=reorder(Another_Random$Gene_Name, Another_Random$`P-Value`), y=Another_Random$Spearman_Coefficient)) +
  geom_bar(stat="identity" , width=0.3, col="red") +
  theme_light() +
  theme(axis.text=element_text(size=3)) +
  labs(title="BACH1 correlations with Lipid Catabolic Process Gene Set", y="Spearman Coefficient", x="P-Value")
         
######produce a Tree plot for a simulation
set.seed(14)
##
spearman_values <- c()
P_values <- c()
row_end <- dim(BRCA_genes_altered)[1]
col_end <- dim(BRCA_genes_altered)[2]
A_very_interesting_gene_index <- sample(1:row_end, 1, replace=TRUE)
A_very_interesting_gene <- as.numeric(BRCA_genes_altered[A_very_interesting_gene_index, 4:col_end])
row_indices <- sample(1:row_end, 400, replace=TRUE)
random_rows <- BRCA_genes_altered[row_indices,]
for (A_gene in 1:400){
  A_correlation <- cor.test(as.numeric(random_rows[A_gene, 4:col_end]), A_very_interesting_gene, method = c("spearman"))
  spearman_values <- c(spearman_values, A_correlation$estimate)
  P_values <- c(P_values, A_correlation$p.value)
}
Simulation_tibble <- tibble("Gene_Name"=random_rows$entrez_dictionary, "Spearman_Coefficient"=spearman_values, "P_Value"=P_values)
Simulation_tibble <- arrange(Simulation_tibble, Simulation_tibble$`P_Value`)
ggplot(data=subset(Simulation_tibble, !is.na(Simulation_tibble$Gene_Name)), aes(x=reorder(Gene_Name, P_Value), y=Spearman_Coefficient)) +
  geom_bar(stat="identity" , width=0.3, col="red") +
  theme_light() +
  theme(axis.text=element_text(size=3)) +
  labs(title="Simulation Tree Plot", y="Spearman Coefficient", x="P-Value")



