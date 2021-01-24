
##############
### README ###
##############

# BLAST of all against all contigs was done with BLAST_allVSall.sh script
# The results are in blastout directory
# All combinations of contigs are in BLAST_comb.txt
# Positions of GENE within contigs are in NDM_positions.txt
# The BLAST output s converted into a network with blast2net.py
# The netwoek is in blast_network.net

# This script analyses the network in order to identify structural variants
# It produces 


# Set working directory
setwd(getSrcDirectory()[1])

library(data.table)
library(igraph)
library(seqinr)
library(data.tree)
library(ape)
library(ggtree)
library(ggplot2)
library(scales)
library(cowplot)
library(svglite)


#####################
# SETTING VARIABLES #
#####################

# read the network and create a graph
network <- fread("results/blast_network.net",stringsAsFactors = FALSE)
colnames(network) <- c("V1", "V2","weight_right")
network <- network[,c(1,2,3)]
g <- graph_from_data_frame(network,directed = FALSE)
#vertices_nr <- length(V(g))
#big_comp <- length(V(g))

# Get contig lengths
seq_directory="test_data/"
files <- list.files(seq_directory,full.names = T)
files <- files[grep("\\.fna$|\\.fasta$|\\.fsa$",files)]
contig_lens <- c()
for(i in files){
  fasta <- read.fasta(i)
  contig_lens <- c(contig_lens,unlist(lapply(fasta,length)))
}
contig_lens <- data.table(ID = names(contig_lens), length = contig_lens)
setkey(contig_lens,ID)

# start position for the splitting threshold
thrsh <- 0
# set splitting threshold step
step <- 100
# counter measures the position in the component_membs table
counter <- 1
# set end position
end <- 4500

# creating data frame for storing the results
#stats <- data.frame(threshold=NA, big_cluster_size=NA, components_nr=NA, singletons_nr=NA)
component_membs <- data.frame(row.names = V(g)$name)



#########################
# SPLITTING THE NETWORK #
#########################

while(thrsh<=end){
  network_thrsh <- network[network$weight>=thrsh]
  g <- graph_from_data_frame(network_thrsh,directed = FALSE)
  comp <- components(g)
  component_membs[names(comp$membership),counter] <- comp$membership
  
  # gettin contigs of sufficient length that are singletons.
  x <- contig_lens[!names(comp$membership)] 
  x <- x[x$length > thrsh+step,]$ID
  y <- max(comp$membership)+1
  if(length(x) > 0)
    component_membs[x,counter] <- y:(y+length(x)-1)
  
  thrsh <- thrsh+step
  counter <- counter+1
}




##########################
# TRACKING THE SPLITTING #
##########################

# this bit of code track how the network was split.
# The data of how the splitting happened is contained in component_membs

# vector of thresholds
thrsh <- seq(0,step*(ncol(component_membs)-1),step)

# Represent the most common structural variants 
# (i.e. struct variants which numbers are >= limit)
limit <- 5
groups <- unique(component_membs[,1])
other_counter <- 1
for (i in 2:ncol(component_membs)){
  for(group in groups){
    counter <- 1
    group_membs <- rownames(component_membs[component_membs[,i-1]==group & 
                                              !is.na(component_membs[,i]),])
    next_col <- component_membs[group_membs,i]
    next_groups <- sort(table(next_col),decreasing = TRUE)
    if(!grepl("otherVariant", group))
    # if the group already is proclaimed other variant (i.e members are < limit)
    # then this ensures the children of the group do not get another 
    # 'otherVariant' label
      otherVar_groups <- names(which(next_groups < limit))
    next_groups <- names(next_groups)
    next_groups <- next_groups[!next_groups %in% otherVar_groups]

    # naming first the next_groups
    if(length(next_groups) == 1)
      next_col[next_col == next_groups] <- group
    else if(length(next_groups) != 0){
      for(next_group in next_groups){
        new_group_name <- paste(group,counter,sep="_")
        counter <- counter + 1
        next_col[next_col==next_group] <- new_group_name
      }
    }
    
    # naming the other variants
    if(length(otherVar_groups) != 0){
      for(next_group in otherVar_groups){
        new_group_name <- paste("otherVariant",other_counter,sep="_")
        other_counter <- other_counter + 1
        next_col[next_col==next_group] <- new_group_name
      }
    }
    component_membs[group_membs,i] <- next_col
  }
  groups <- table(component_membs[!is.na(component_membs[,i]),i])
  groups <- names(sort(groups,decreasing = TRUE))
}


##################
# RENDERING TREE #
##################

# singleton is an early terminated sequence (i.e. sequence cutting short)
component_membs[is.na(component_membs)] <- "singleton" 

# component_members -> columns are positions in sequences
#                      rows are the sequences
#                      within each cell there is a name of the structural var.
#                      This name codes for the tree of structural variations:
#                      e.g. 1_2  means that this is the second structural var
#                      that emerged from the parent (root brranch 1)
#                      e.g. 1_3_2_2  means that this is the second structural var
#                      that emerged from the structural var 2 which is child of
#                      the struct var 3 which parent is the root branch.
#                      IMPORTANT: otherVariants do follow this tree structure,
#                      but they are detached from the main tree naming, as if
#                      new tree started growing on the main tree
#                      This was done for plotting purposes

write.csv(component_membs, "results/Table_of_structural_variations_raw.csv",quote = FALSE)

# before rendering the tree we need to count other variants at each position
# Counting other variants and singletons for plotting
singletons <- lapply(component_membs, function(x) {sum(x=="singleton")})
singletons <- data.frame(singletons)
singletons[2,] <- thrsh
singletons <- data.frame(t(singletons))
colnames(singletons) <- c("Count","Threshold")
otherVar <- apply(component_membs,2,function(x){
  return(sum(grepl("otherVariant",x)))})
otherVar <- data.frame(otherVar)
otherVar[,2] <- thrsh
colnames(otherVar) <- c("Count","Threshold")
# removing other variants to render a tree
component_membs <- data.frame(lapply(component_membs, function(x) {
  gsub("otherVariant.*", "otherVariant", x)}))

# samples for each nodes
nodes <- unique(unlist(component_membs))
samples <- c("root"=nrow(component_membs))

# rendering the tree goes like this
# cmpnt_membs -->LIST OF LISTS (tree as nested lists) --> PHYLO --> tree ATA FRAME
lol <- list(thrsh=0,branchHeight=0)

# percent change cutoff is to ensure preservation of the most samples on each branch.
# it aims to keep the most interesting structural variations
# if the tree starts splitting to fast (ie the percentage of branching is > 0.1)
# it will stop the process of the tree branching out on a particular branch
percent_change_cutoff <- 0.1
for(name in nodes){
  t <- table(which(component_membs==name,arr.ind=TRUE)[,2])
  if(max(t)<20){
    t2 <- max(which(component_membs==name,arr.ind=TRUE)[,2])
    samples[name] <- t[as.character(t2)]
    t <- thrsh[t2]
  }
  else{
    t2 <- 1-t/max(t)
    t2 <- max(as.numeric(names(t2[t2<percent_change_cutoff])))
    samples[name] <- t[as.character(t2)]
    t <- thrsh[t2]
  }
  t_sum <- 0
  splt <- unlist(strsplit(name,"_"))
  exprs <- ""
  for (j in 1:length(splt)){
    a <- paste('[["',paste(splt[1:j],collapse = "_"),'"]]',sep="")
    exprs <- paste(exprs,a,sep="")
    if (j!=length(splt))
      t_sum <- t_sum + eval(parse(text=paste("lol",exprs,"$branchHeight",sep="")))
  }
  exprs1 <- paste("lol",exprs,"$thrsh <-",t,sep="")
  eval(parse(text=exprs1))
  t <- t-t_sum
  exprs2 <- paste("lol",exprs,"$branchHeight <-",t,sep="")
  eval(parse(text=exprs2))
}

# making a tree
tree <- as.Node(lol)
tree <- as.phylo(tree,heightAttribute = "thrsh")
tree$edge.length <- tree$edge.length*-1
# making a tree data frame using fortify from ggtree
tree_table <- fortify(tree)
tree_table <- tree_table[!(tree_table$label %in% c("singleton","otherVariant")),]
tree_table$y <- tree_table$y-1
for(i in 1:nrow(tree_table)){
  tree_table[i,"samples"] <- samples[tree_table$label[i]]
}
tree_table[tree_table$label == "Root","y"] <- tree_table[tree_table$label == "1","y"]

write.table(tree_table,"results/Tree_of_structural_variations_table.tsv",quote = F,sep = "\t",
            row.names = F)



############
# PLOTTING #
############

p1 <- ggtree(tree_table,aes(size=samples, color=samples),branch.length = 0) + 
  theme_tree2() + scale_size(guide = FALSE) + 
  scale_x_continuous(name="Position (bp)", position = "top",
                     expand = expansion(mult = c(0.02, 0.02), add = c(0,50)),
                     limits = c(0,max(thrsh)),
                     breaks = pretty_breaks(n=10)) +
  scale_color_gradient(name = "# Samples", trans = "log",
                      expand = expansion(mult = c(0.002, 0.002)),
                      breaks = trans_breaks("log10", function(x) round(10 ^ x)), 
                      low = "#F0E68C", high = "darkred") + 
  coord_cartesian(clip = "off") +
  ggtitle("Splitting of Structural Variants") +
  theme(legend.position=c(0.1,0.8),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_line(linetype = 2,size = 0.5),
        panel.grid.minor.x = element_line(linetype = 2,size = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 0.3, vjust = 1),
        plot.margin=unit(c(1,0,0.1,1), "cm"))

# adding node labels that indicate number of samples at each node
p1 <- p1 + geom_label(tree_table, 
                      mapping = aes(x=x,y=y,label=samples),inherit.aes = FALSE)

p_single <- ggplot(singletons,aes(x=Threshold,y=Count)) + 
  geom_area(size=3,fill="#998ec3") +
  scale_x_continuous(name=NULL, position = "top",
                     expand = expansion(mult = c(0.02, 0.02), add = c(0,50)),
                     limits = c(0,max(thrsh)),
                     breaks = pretty_breaks(n=10)) + 
  scale_y_reverse(name="# Contigs Cutting Short",
                  expand = expansion(mult = c(0, 0)),
                  breaks = pretty_breaks()) +
  scale_color_manual(values = "#225ea8") + theme_classic() +
  theme(legend.position = "none", legend.title = element_blank(),
        panel.grid.major.y = element_line(linetype = 2,size = 0.5),
        panel.grid.minor.y = element_line(linetype = 2,size = 0.5),
        panel.grid.major.x = element_line(linetype = 2,size = 0.5),
        panel.grid.minor.x = element_line(linetype = 2,size = 0.5),
        axis.text = element_text(color="black"),
        axis.text.x = element_blank())

p_otherVar <- ggplot(otherVar,aes(x=Threshold,y=Count)) + 
  geom_area(size=3, fill="#bdbdbd") +
  scale_x_continuous(name=NULL, position = "top",
                     expand = expansion(mult = c(0.02, 0.02), add = c(0,50)),
                     limits = c(0,max(thrsh)),
                     breaks = pretty_breaks(n=10)) + 
  scale_y_reverse(name="# Other Structural Vatiants",
                  expand = expansion(mult = c(0, 0)),
                  breaks = pretty_breaks()) +
  scale_color_manual(values = "#225ea8") + theme_classic()+
  theme(legend.position = "none", legend.title = element_blank(),
        panel.grid.major.y = element_line(linetype = 2,size = 0.5),
        panel.grid.minor.y = element_line(linetype = 2,size = 0.5),
        panel.grid.major.x = element_line(linetype = 2,size = 0.5),
        panel.grid.minor.x = element_line(linetype = 2,size = 0.5),
        axis.text = element_text(color="black"),
        axis.text.x = element_blank())

p_cow <- plot_grid(p1,p_single,p_otherVar,align = "v",
                   axis = "t",ncol = 1,rel_heights = c(3,1,1))
  
pdf("results/Tree_of_structural_variations.pdf",width = 12,height = 10)
p_cow
dev.off()

svglite("results/Tree_of_structural_variations.svg",width = 12,height = 10)
p_cow
dev.off()


