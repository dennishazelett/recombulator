library(tidyverse)
library(GenomicRanges)
library(foreach)
library(shiny)

map_position <- function(linkage_group, position) {
  ## "takes nt position in nucleotides as argument"
  x <- position/1e6
  if (linkage_group == 'X') {
    genetic_map <- -0.0097*x^3 + 0.2996*x^2 + 1.1626*x - 1.8904
  } else if (linkage_group == '2L') {
    genetic_map <- -0.0099*x^3 + 0.2087*x^2 + 2.5183*x - 0.8057
  } else if (linkage_group == '2R') {
    genetic_map <- -0.0083*x^3 + 0.3627*x^2 - 1.5224*x + 58.252
  } else if (linkage_group == '3L') {
    genetic_map <- -0.006*x^3 + 0.1091*x^2 + 2.6663*x -1.6899
  } else if (linkage_group == '3R') {
    genetic_map <- -0.0038*x^3 + 0.233*x^2 - 1.5567*x + 50.127
  }
  return(genetic_map)
}

# genelib <- read.delim2(file = 'geneheader.txt', header = F, sep = ';', comment.char = '#')
# libsize <- dim(genelib)[1]
# fly_genes <- GRanges()
# 
# create_gene <- function(i) {
#   V1 <- str_split(genelib$V1, ' ')
#   V2 <- str_split(genelib$V2, '=')[[i]][2]
#   fbgn <- str_split(V1[[i]][1], '>')[[1]][2]
#   type <- str_split(V1[[i]][2], '=')[[1]][2]
#   linkage_group <- str_split(V2, ':')[[1]][1]
#   loc <- str_extract(str_split(V2, ':')[[1]][2], regex("[0-9]*\\.\\.[0-9]*"))
#   start <- as.numeric(str_split(loc, pattern = '\\.\\.')[[1]][1])
#   end   <- as.numeric(str_split(loc, pattern = '\\.\\.')[[1]][2])
#   if (str_detect(V2, 'complement')) {
#     strand <- '-'
#   } else { 
#     strand <- '+'  
#   }
#   gene_symbol <- str_split(genelib$V4, '=')[[i]][2]
#   
#   gene <- GRanges(
#     seqnames <- linkage_group,
#     ranges = IRanges(start = start,
#                      end = end),
#     strand = strand,
#     mcols = data.frame(fbgn = fbgn,
#                        gene_symbol = gene_symbol))
#   
#   gene
# }  
# 
# fly_genes <- foreach(i=1:libsize, .combine = c, .inorder = TRUE) %do% create_gene(i)
# 
# save(fly_genes, file = 'fly_genes.Rda')

load('fly_genes.Rda')
seqnames(fly_genes[1,])

estimateRF <- function(distance) {
  ## "estimate the recombinant fraction given two genetic map loci"
  if (distance < 20) {
    rf <- distance/100 } else {
    ## adapted from Ashburner "Drosophila A Laboratory Handbook", pp 458
    ## see Kosambi's formula
    y <- estimateRF(distance/2)
    rf <- 2*y/(1+4*y^2)
  }
  return(rf)
}

predictRF <- function(locus1, locus2) {
  # inputs are GRanges (single line)
  locus1_lg <- as.character(seqnames(locus1)) 
  locus2_lg <- as.character(seqnames(locus2)) 
  if (locus1_lg==locus2_lg) {
    locus1_prox <- map_position(locus1_lg, start(locus1))
    locus1_dist <- map_position(locus1_lg, end(locus1))
    locus2_prox <- map_position(locus2_lg, start(locus2))
    locus2_dist <- map_position(locus2_lg, end(locus2))
    a <- abs(locus1_prox-locus2_dist)
    b <- abs(locus1_dist-locus2_prox)
    if (a < b) {
      return(estimateRF(a))
    } else {
      return(estimateRF(b))
    }
  } else {return(0.5)}
}

recombine_two_genes <- function(gene1, gene2) {
  locus1 <- fly_genes[fly_genes$mcols.gene_symbol==gene1,]
  locus2 <- fly_genes[fly_genes$mcols.gene_symbol==gene2,]
  predictRF(locus1, locus2)
}

number_of_progeny <- function(crossover_prob, num_events, confidence, n=0) {
  prob_num_events_n <- pbinom(
    num_events, 
    size = num_events+n, 
    prob = crossover_prob
  )
  if (prob_num_events_n < 1-confidence) {
    return(num_events + n)
  } else {
    number_of_progeny(crossover_prob, num_events, confidence, n = num_events+n)
  }
}

graph_result <- function(gene1, gene2, confidence=0.95, nsuccesses = 3) {
  ptest <- recombine_two_genes(gene1, gene2)
  probs <- pretty(c(0.4, 0.99), n=60)
  progeny <- numeric(length(probs)) 
  for (i in 1:length(probs)) {progeny[i] <- number_of_progeny(ptest, nsuccesses, probs[i])}
  my_data <- data.frame(probs, progeny)
  minprogeny <- my_data$progeny[min(which(my_data$probs >= confidence))]
  my_graph <- ggplot(my_data, aes(x=progeny, y=probs)) + 
    geom_line() + 
    geom_hline(yintercept = confidence, lty=2, col='darkgrey') +
    geom_vline(xintercept = minprogeny, lty=3, col='red') + 
    annotate(geom="text", x=minprogeny, y= 0.7, label=paste("n =", minprogeny), color = 'red') +
    xlab("number of progeny screened") +
    ylab(paste("probability")) +
    ggtitle(paste(100 * confidence, "% chance of obtaining at least", nsuccesses, "recombinants"))
  my_graph
}


