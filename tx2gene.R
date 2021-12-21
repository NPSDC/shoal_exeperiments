library(tximport)
library(readr)
args = commandArgs(trailingOnly=TRUE)

library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

k <- keys(edb, keytype = "TXNAME")
tx2gene <- select(edb, k, "GENEID", "TXNAME")
tx2gene <- tx2gene[c('TXNAME', 'GENEID')]
txi <- tximport(c(args[1]), type = "salmon", tx2gene = tx2gene, ignoreTxVersion = T)

# # head(tx2gene)
# txi <- tximport(c(args[1]), type = "salmon", tx2gene = tx2gene)
# names(txi)

write.table(txi, args[2], quote=FALSE, sep="\t")