# install binr package from CRAN
# install.packages('binr')

library('binr')
library('stringr')
library('magrittr')

digest_files = list.files(pattern = 'digests_');
enzymes <- lapply(digest_files, function(x){
  return( (
    (x %>% str_split('digests_') )[[1]][2] %>%
      str_split('.tsv')
    )[[1]][[1]])
}) %>% unlist

parseEnzTables <- function(enzyme.name) {
  assign(paste0(enzyme.name, '.table'),
         readMzTable(enzyme.name),
         envir = parent.env(environment()))
}

readMzTable <- function (enz) {
  mz.table <- read.table(paste0('digests_', enz, '.tsv'),
           sep = '\t',
           header = FALSE)
}

lapply(enzymes, parseEnzTables)
# new variables in environment: arg_c.table, glu_c.table, ...

processEnzTables <- function(enzyme.name) {
  min <- 1000
  enz.table <- get(paste0(enzyme.name, '.table'), envir = parent.env(environment()))
  enz.mass.vector <- enz.table[3]
  print(enzyme.name)
  bincount_total <- 0
  for (j in 1:5) {
    binmin <- (j-1)*100 + min
    binmax <- j*100 + min
    bincount <- sum(enz.mass.vector >= binmin & enz.mass.vector < binmax)
    binstatement <- paste0(binmin, " to ", binmax, " : ", bincount)
    print(binstatement)
    bincount_total <- bincount_total + bincount
  }
  print(paste0("Total for ", enzyme.name, " : ", bincount_total))
}


for (i in 1:length(enzymes)) {
  processEnzTables(enzymes[i])
  #print(enzymes[i])
}

for (i in 1:length(enzymes)) {
  enzyme.name <- enzymes[i]
  enz.table <- get(paste0(enzyme.name, '.table'))
  enz.seq.vector <- enz.table[6][,]
  mean_seq_length <- mean(nchar(as.character(enz.seq.vector)))
  print(paste0("Average sequence length for ", enzyme.name, " : ", mean_seq_length))
}

library(ggplot2)

max.seqlength <- 50
seqlengths <- nchar(as.character(arg_c.table[6][,]))
arg_c_lengths <- data.frame(length = seqlengths[seqlengths < max.seqlength])
seqlengths <- nchar(as.character(glu_c.table[6][,]))
glu_c_lengths <- data.frame(length = seqlengths[seqlengths < max.seqlength])
seqlengths <- nchar(as.character(lys_c.table[6][,]))
lys_c_lengths <- data.frame(length = seqlengths[seqlengths < max.seqlength])
seqlengths <- nchar(as.character(trypsin.table[6][,]))
trypsin_lengths <- data.frame(length = seqlengths[seqlengths < max.seqlength])
arg_c_lengths$enzyme <- 'arg_c'
glu_c_lengths$enzyme <- 'glu_c'
lys_c_lengths$enzyme <- 'lys_c'
trypsin_lengths$enzyme <- 'trypsin'

enz.lengths <- rbind(arg_c_lengths, glu_c_lengths, lys_c_lengths, trypsin_lengths)
ggplot(enz.lengths, aes(length, fill = enzyme)) + geom_density(alpha = 0.2) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())