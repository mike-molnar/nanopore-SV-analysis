library(karyoploteR)
args = commandArgs(trailingOnly=TRUE)

cov_data_normal <- read.table(args[1], sep="\t", strip.white=TRUE)
chromosome_normal <- cov_data_normal[, "V1"]
x_normal <- cov_data_normal[, "V3"]
y_normal <- cov_data_normal[, "V5"]

cov_data_tumor <- read.table(args[2], sep="\t", strip.white=TRUE)
chromosome_tumor <- cov_data_tumor[, "V1"]
x_tumor <- cov_data_tumor[, "V3"]
y_tumor <- cov_data_tumor[, "V5"]

pdf(args[3], height=30, width=25)
pp <- getDefaultPlotParams(plot.type=1)
pp$data1outmargin <- 80
kp <- plotKaryotype(genome=args[4], chromosomes=c("chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrY"), labels.plotter = NULL, plot.params = pp, plot.type=1)
kpAddChromosomeNames(kp, cex=1.3)

kpAddLabels(kp, labels = "Copy\nCount", label.margin = 0.03, cex=1.3)
kpAxis(kp, ymin=0, ymax=4, side=1, numticks=5)
kpAxis(kp, ymin=0, ymax=4, side=2, numticks=5)
kpAbline(kp, h=c(0.25, 0.5, 0.75), lty="39", data.panel=1)
kpPlotLoess(kp, chr=chromosome_normal, x=x_normal, y=y_normal, conf.interval=NULL, span=0.01, col="blue")
kpPlotLoess(kp, chr=chromosome_tumor, x=x_tumor, y=y_tumor, conf.interval=NULL, span=0.01, col="red")
legend("right", fill = c("blue", "red"), legend = c(args[5], args[6]), bty="n", ncol=1, cex=1.7)
dev.off()