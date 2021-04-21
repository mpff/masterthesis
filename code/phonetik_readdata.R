# From '..phonetik/code/Rfiles/readdata.R'
pouplier <- read.csv(file = "./data/phonetik/data/pouplier.csv", header = FALSE)
pouplier <- as.matrix(pouplier)
pouplier[which(pouplier == -99)] <- NA
n <- nrow(pouplier) / 5

cols = rep(1:n, each = 5)

pdf("./results/data.pdf")
plot(0:1, range(pouplier[!is.na(pouplier)]), type = "n",
     xlab = 'time', ylab = "index")
for (i in 1:nrow(pouplier)){
  m = sum(!is.na(pouplier[i,]))
  lines((0:(ncol(pouplier)-1)) / (m-1), pouplier[i,],
        col = cols[i])
}
dev.off()

reflabels <- read.csv(file = "./data/phonetik/data/reflabels.csv", header = FALSE)
reflabels <- paste(reflabels[,1], reflabels[,2], reflabels[,3],
                   reflabels[,4], reflabels[,5], reflabels[,6],
                   reflabels[,7], reflabels[,8], sep = "")

words <- rep(reflabels, each = 5) #contains words for each row
repeats <- rep(1:5, n)  #contains repeat number for each row

library(manipulate)
manipulate(plot(pouplier[i, ], t = "l"), i = slider(1, 140))
