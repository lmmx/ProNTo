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