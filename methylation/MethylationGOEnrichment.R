## MethylationGOEnrichment.R
## Jessie Pelosi 
## Last modified 1 November, 2024 

# Read in 

library(dplyr)
library(ggplot2)
library(topGO)
library(clusterProfiler)
library(viridis)


# Read in differentially methylated sites and regions (5kb windows)

sites_hypo <- read.delim("Lygodium/v1.1.9/sites_hypoTSS-Lygodium_1Nov2024.txt", sep='')
sites_hyper <- read.delim("Lygodium/v1.1.9/sites_hyperTSS-Lygodium_1Nov2024.txt", sep='')

sites_hypo_filtered <- sites_hypo %>% 
  filter(dist.to.feature > -5000  & dist.to.feature < 5000) %>% 
  mutate(Gene = feature.name)

sites_hyper_filtered <- sites_hyper %>% 
  filter(dist.to.feature > -5000  & dist.to.feature < 5000) %>% 
  mutate(Gene = feature.name)


tiles_hypo <- read.delim("Lygodium/v1.1.9/5kbtiles_hypoTSS-Lygodium_1Nov2024.txt", sep = '')
tiles_hyper <- read.delim("Lygodium/v1.1.9/5kbtiles-hyperTSS-Lygodium_1Nov2024.txt", sep = '')

tiles_hypo_filtered <- tiles_hypo %>% 
  filter(dist.to.feature > -5000  & dist.to.feature < 5000) %>% 
  mutate(Gene = feature.name)

tiles_hyper_filtered <- tiles_hyper %>% 
  filter(dist.to.feature > -5000  & dist.to.feature < 5000) %>% 
  mutate(Gene = feature.name)


#Let's look at some basic stats: 

features <- read.delim("Lygodium/v1.1.9/diffMethFeatures-Lygo.txt")
features_tiled <- features %>% 
  filter(Dataset == "tiled")

ggplot(data = features_tiled, mapping=aes(x="", y = PercDiff, fill = feature)) + geom_col(colour = "white") +
  facet_wrap(~MethylationStatus) + coord_polar("y", start = 0) + theme_void() + 
  scale_fill_viridis_d(option = "viridis")


features_site <- features %>% 
  filter(Dataset == "site")

ggplot(data = features_site, mapping=aes(x="", y = PercDiff, fill = feature)) + geom_col(colour = "white") +
  facet_wrap(~MethylationStatus) + coord_polar("y", start = 0) + theme_void() + 
  scale_fill_viridis_d(option = "viridis")


# Are there functions that are enriched in these different gene sets? 

# Let's do a KEGG enrichment first: 
LygoBlast <- read.delim("../Lygo_Arabidopsis.blast") 

lygo_genes_AT <- LygoBlast %>% 
  dplyr::select(Gene, ArabidopsisGene, evalue) %>%
  group_by(Gene) %>% 
  arrange(evalue, .by_group = TRUE) %>% 
  filter(row_number()==1)

sites_hypo_filtered_AT <- inner_join(sites_hypo_filtered, lygo_genes_AT)

kegg <- enrichKEGG(gene = sites_hypo_filtered_AT$ArabidopsisGene, organism = "ath", pvalueCutoff = 0.05)
kegg


sites_hyper_filtered_AT <- inner_join(sites_hyper_filtered, lygo_genes_AT)
kegg <- enrichKEGG(gene = sites_hyper_filtered_AT$ArabidopsisGene, organism = "ath", pvalueCutoff = 0.05)
kegg



tiles_hypo_filtered_AT <- inner_join(tiles_hypo_filtered, lygo_genes_AT)
kegg <- enrichKEGG(gene = tiles_hypo_filtered_AT$ArabidopsisGene, organism = "ath", pvalueCutoff = 0.05)
kegg

tiles_hyper_filtered_AT <- inner_join(tiles_hyper_filtered, lygo_genes_AT)
kegg <- enrichKEGG(gene = tiles_hyper_filtered_AT$ArabidopsisGene, organism = "ath", pvalueCutoff = 0.05)
kegg


