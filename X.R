setwd("~/Dropbox/GitHub/X")

# Mouse symbol, Human symbol, Gene name, Gene summary

library(mygene)
mus2hg <- read.delim("hg2mus.map", header = F, stringsAsFactors = F)
geneId <- unique(ens.map$Symbol)
table <- queryMany(geneId, scopes="symbol", species="mouse", fields = c("name", "summary"))
table <- table[!duplicated(table$query), ]
table$hg <- mus2hg$V1[match(geneId, mus2hg$V4)]
query.hg <- queryMany(table$hg, scopes="symbol", species="human", fields = c("name", "summary"))
query.hg <- query.hg[!duplicated(query.hg$query), ]
summary_hg <- query.hg$summary[match(table$hg, query.hg$query)]
table$summary <- paste("Mus:", table$summary, "Hg:", summary_hg)
table <- table[c("query", "name", "hg", "summary")]
rownames(table) <- NULL
table <- as.data.frame(table)
save(table, file = "summary.rdt")
