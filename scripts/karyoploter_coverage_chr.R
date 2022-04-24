library(karyoploteR)
args = commandArgs(trailingOnly=TRUE)

# usage: karyotyper_chr.R [normal coverage] [tumor coverage] [normal BAF] [tumor BAF] [normal name] [tumor name] [reference genome] [chromosome] [output filename]

print("Loading datasets...")
# load datasets
cov_data_normal <- read.table(args[1], sep="\t", strip.white=TRUE)
cov_data_tumor <- read.table(args[2], sep="\t", strip.white=TRUE)
baf_data_normal <- read.table(args[3], sep="\t", strip.white=TRUE)
baf_data_tumor <- read.table(args[4], sep="\t", strip.white=TRUE)

# extract data from files 
chromosome_cov_normal <- cov_data_normal[, "V1"]
x_cov_normal <- cov_data_normal[, "V2"]
y_cov_normal <- cov_data_normal[, "V5"]

chromosome_cov_tumor <- cov_data_tumor[, "V1"]
x_cov_tumor <- cov_data_tumor[, "V2"]
y_cov_tumor <- cov_data_tumor[, "V5"]

chromosome_baf_normal <- baf_data_normal[, "V1"]
x_baf_normal <- baf_data_normal[, "V2"]
y_baf_normal <- baf_data_normal[, "V4"]

chromosome_baf_tumor <- baf_data_tumor[, "V1"]
x_baf_tumor <- baf_data_tumor[, "V2"]
y_baf_tumor <- baf_data_tumor[, "V4"]

print("Plotting the data...")
# build the karyoplot
pdf(args[9], width=25, height=10)

print("Plotting the data...")
# build the karyoplot
kp <- plotKaryotype(genome = args[7], chromosomes=args[8], plot.type=2, labels.plotter=NULL)
kpAddCytobandLabels(kp, srt=65, cex=1.05, col="red4")
kpAddChromosomeNames(kp, cex=1.7)

#add a legend to the plot
legend("top", fill = c("blue", "red"), legend = c(args[5], args[6]), bty="n", ncol=2, cex=2.5)

# plot the copy counts
kpAddLabels(kp, labels = "Copy\nCount", r0=0.0, r1=1.0, label.margin = 0.03, cex=1.7, data.panel=1)
kpAxis(kp, r0=0.0, r1=1.0, ymin=0, ymax=4, side=1, numticks=5, data.panel=1)
kpAxis(kp, r0=0.0, r1=1.0, ymin=0, ymax=4, side=2, numticks=5, data.panel=1)
kpAbline(kp, r0=0.0, r1=1.0, h=c(0.25, 0.5, 0.75), lty="39", data.panel=1)
kpPlotLoess(kp, chr=chromosome_cov_normal, x=x_cov_normal, y=y_cov_normal, r0=0.0, r1=1.0, conf.interval=NULL, span=0.01, col="blue")
kpPlotLoess(kp, chr=chromosome_cov_tumor, x=x_cov_tumor, y=y_cov_tumor, r0=0.0, r1=1.0, conf.interval=NULL, span=0.01, col="red")

# plot the b-allele frequencies
kpAddLabels(kp, labels = "B-Allele\nFrequency\n(%)", r0=1.0, r1=0.05, label.margin = 0.03, cex=1.7, data.panel=2)    
kpAxis(kp, r0=0.45, r1=0, ymin=0, ymax=100, side=1, data.panel=2)
kpAxis(kp, r0=0.45, r1=0, ymin=0, ymax=100, side=2, data.panel=2)
kpAxis(kp, r0=1.0, r1=0.55, ymin=0, ymax=100, side=1, cex=0.95, data.panel=2)
kpAxis(kp, r0=1.0, r1=0.55, ymin=0, ymax=100, side=2, cex=0.95, data.panel=2)
kpPoints(kp, chr=chromosome_baf_normal, x=x_baf_normal, y=y_baf_normal, r0=0.45, r1=0, cex=0.07, col="blue", data.panel=2)
kpPoints(kp, chr=chromosome_baf_tumor, x=x_baf_tumor, y=y_baf_tumor, r0=1.0, r1=0.55, cex=0.07, col="red", data.panel=2)

dev.off()