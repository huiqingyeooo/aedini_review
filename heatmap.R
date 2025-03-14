# This script was used to create a tree and a heatmap summarizing the review of aedini phylogenies across published papers

rm(list=ls())
setwd(" ")
library(readxl); library(ggtree); library(dplyr)

#clade_assignment<-read.csv("clade_assignment.csv")
tree<-read.tree("mol_hypothesis.nwk")
dat<-read_excel("data.xlsx");nrow(dat)

#create column with combination of genus, subgenus, author, and year to obtain unique rows for matching
dat$uniq<-paste(dat$Genus, dat$Subgenus, dat$Author, dat$Year, sep="_")
uniq2<-as.data.frame(unique(dat$uniq)); colnames(uniq2)<-"uniq"; nrow(uniq2)

#left join to merge unique list with dat
#specify multiple=first to return on the first match to the unique list
dat2 <- left_join(uniq2, dat, by="uniq",multiple="first")
nrow(dat2)

# create columns for genus+subgenus and author+year
dat2$gen_subgen<-paste(dat2$Genus, dat2$Subgenus, sep="_")
dat2$ref<-paste(dat2$Author, dat2$Year)

# create pivot table
dat3<-dat2 %>% group_by(ref) %>% select(gen_subgen, Clade)%>% #first group data by unique ref column, then select taxa and clade columns
  pivot_wider(names_from = "ref", values_from ="Clade")%>% #apply pivot
  rownames_to_column()%>% #to obtain column with gen_subgen column
  left_join((as.data.frame(tree$tip.label)%>%`colnames<-`(c("gen_subgen"))),.,by="gen_subgen")%>% #extracts all gen_subgen from tip labels and matches taxa from dat3 to the full gen_subgen list
  select(!(rowname))%>% #remove column with rowname numbers
  column_to_rownames('gen_subgen')%>% #convert first column to rownames to suit format needed for heatmap
  rename_with(make.names)%>% #makes syntactically valid names for column headers, e.g. replaces spaces with .
  replace(is.na(.), 0) %>% #replace empty cells (NA) with 0
  select(order(names(.)))%>% #arrange columns by alphabetical order
  relocate(RHK.2009, .after = last_col())

# lists taxa in clades A, B and Psorophora
clade_a<-extract.clade(tree, node = findMRCA(tree, tips = c('Opifex_Opifex','Aedes_Abraedes')))$tip.label
clade_b<-extract.clade(tree, node = findMRCA(tree, tips = c('Zeugnomyia_0','Aedes_Aedes')))$tip.label
psor<-extract.clade(tree, node = findMRCA(tree, tips = tree$tip.label[grepl('Psorophora', tree$tip.label)]))$tip.label

clades<-list(clade_abp = c(clade_a, clade_b, psor),clade_a = clade_a, clade_b = clade_b, psorophora = psor)
tree.clades<-groupOTU(.data=tree, .node=clades)

# plot tree
(tree.dat<-ggtree(tree.clades, aes(color = group), lwd = 0.7) %<+% dat3)+geom_text(aes(label=node)) #impt to attach dataset to tree with %<+%
tree.dat<-ggtree(tree.clades, aes(color = group), lwd = 0.7) %<+% dat3
tree.dat2<-tree.dat %>% ggtree::rotate(99) %>% ggtree::rotate(100) %>% ggtree::rotate(101)

# add in heat map
pdf("clade_allTrees_suppl.pdf", width=8, height=10)
gheatmap(tree.dat2,dat3, width=6, offset=8, colnames=F)+
  #add colours for heatmap
  scale_fill_manual(breaks=c("0", "A", "B","Neither","Outgroup","Psorophora"), 
                  values=c("#F5F5F5","#21908CFF","#440154FF","#FDE725FF","#5DC863FF","#3B528BFF"),
                  name="Clade")+
  #add colours for tree branches
  scale_color_manual(values = c("#21908CFF","black","#440154FF","#3B528BFF"), 'Clade',
                     labels = c('clade_abp' = 'Root', 'clade_a' = 'Clade A', 'clade_b' = 'Clade B',
                                'psorophora' = 'Psorophora'))+
  geom_tiplab(size=2.5, offset=0.5, color="black")+
  scale_x_ggtree()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 5))+
  theme(aspect.ratio = 1.6)
dev.off()

##### To generate plot for molecular phylogenies only #####
# create list of gen/subgen from molecular phylogenies
mol.list<-dat2 %>% filter(Author!="RHK") %>%
  group_by(ref) %>% select(gen_subgen, Clade)%>% #first group data by unique ref column, then select taxa and clade columns
  pivot_wider(names_from = "ref", values_from ="Clade")%>% #apply pivot
  select(gen_subgen)%>% column_to_rownames('gen_subgen') #to obtain column with gen_subgen column and convert to rownames

# remove tips not present in molecular phylogenies
mol.trees<-geiger::treedata(tree, data=mol.list, warnings = T, sort = T)$phy
length(tree$tip.label); length(mol.trees$tip.label) #check number of tips

# generate matrix for molecular data only
mol.dat <-dat2 %>% filter(Author!="RHK") %>%
  group_by(ref) %>% select(gen_subgen, Clade)%>% #first group data by unique ref column, then select taxa and clade columns
  pivot_wider(names_from = "ref", values_from ="Clade")%>% #apply pivot
  column_to_rownames('gen_subgen')%>% #convert first column to rownames to suit format needed for heatmap
  rename_with(make.names)%>% #makes syntactically valid names for column headers, e.g. replaces spaces with .
  replace(is.na(.), 0) %>% #replace empty cells (NA) with 0
  select(order(names(.))) #arrange columns by alphabetical order

# define clades in mol.trees
clades<-list(clade_abp = c(clade_a, clade_b, psor),clade_a = clade_a, clade_b = clade_b, psorophora = psor)
mol.tree.clades<-groupOTU(.data=mol.trees, .node=clades)

# Plot tree
(mol.tree.dat<-ggtree(mol.tree.clades, aes(color = group), lwd = 0.7) %<+% mol.dat)+geom_text(aes(label=node))
mol.tree.dat<-ggtree(mol.tree.clades, aes(color = group), lwd = 0.7) %<+% mol.dat 
mol.tree.dat2<-mol.tree.dat %>% ggtree::rotate(69) %>% ggtree::rotate(70)%>% ggtree::rotate(71)

# add in heat map
pdf("clade_molTrees.pdf", width=8, height=10)
gheatmap(mol.tree.dat2, mol.dat, width=6, offset=9, colnames=F)+
  #add colours for heatmap
  scale_fill_manual(breaks=c("0", "A", "B","Neither","Outgroup","Psorophora"), 
                    values=c("#F5F5F5","#21908CFF","#440154FF","#FDE725FF","#5DC863FF","#3B528BFF"),
                    name="Clade")+
  #add colours for tree branches
  scale_color_manual(values = c("#21908CFF","black","#440154FF","#3B528BFF"), 'Clade',
                     labels = c('clade_abp' = 'Root', 'clade_a' = 'Clade A', 'clade_b' = 'Clade B',
                                'psorophora' = 'Psorophora'))+
  geom_tiplab(size=2.9, offset=0.5, color="black")+
  scale_x_ggtree()+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 6))+
  theme(aspect.ratio = 1.6)
dev.off()
