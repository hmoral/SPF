#!/usr/bin/Rscript
###### Adapted from Pečnerová et al.
# https://github.com/KHanghoj/leopardpaper/
# original script:
# https://github.com/KHanghoj/leopardpaper/blob/master/mapping_and_qc/reference_filter/depth/04_quantile_thresholds.R

#Get arguments
args <- commandArgs(TRUE)
depth.dist.file <- args[1]
outname <- args[2]

depth.global <- scan(depth.dist.file)

depth.global.cumsum <- cumsum(as.numeric(depth.global))

# Find thresholds.
q99 <- depth.global.cumsum[length(depth.global.cumsum)]*0.99
q01 <- depth.global.cumsum[length(depth.global.cumsum)]*0.01
q99.threshold <- min(which(depth.global.cumsum > q99))
q01.threshold <- min(which(depth.global.cumsum > q01))

# Save on file
write.table(t(c(q01.threshold-1,q99.threshold+1)),paste0(outname,".treshold.txt"),quote=F,row.names=F,col.names=F,sep="\t")

# Plot
png(paste(outname, '-fullrange.png', sep=''), width = 600, height = 300)
#pdf(paste(plotdir, group, '-fullrange.pdf', sep=''), width = 6, height = 3)
par(mfrow=c(1,2))
plot(depth.global, type='l',
    xlab = 'Global depth (count)', ylab = 'Number of sites')
abline(v = q99.threshold, col='red', lty=2); abline(v = q01.threshold, col='red', lty=2)
legend('topright', legend = c(paste('1% thres:', q01.threshold, sep=''),
                            paste('99% thres: ', q99.threshold, sep='')),  bty = "n")
plot(depth.global.cumsum, type='l',
    xlab = 'Global depth (count)', ylab = 'Cumulative number of sites')
abline(v = q99.threshold, col='red', lty=2); abline(v = q01.threshold, col='red', lty=2)
mtext(paste('Full range - ', outname, sep=''), side = 3, line = -2, outer = TRUE)
dev.off()

# Zoomed in.
x_max <- 600
png(paste(outname, '-zoomin.png', sep=''), width = 600, height = 300)
#pdf(paste(plotdir, group, '-zoomin.pdf', sep=''), width = 6, height = 3)
par(mfrow=c(1,2))
plot(depth.global, type='l', xlim=c(0, x_max),
        xlab = 'Global depth (count)', ylab = 'Number of sites')
abline(v = q99.threshold, col='red', lty=2); abline(v = q01.threshold, col='red', lty=2)
legend('topright', legend = c(paste('1% thres:', q01.threshold, sep=''),
                                paste('99% thres: ', q99.threshold, sep='')),  bty = "n")
plot(depth.global.cumsum, type='l', xlim=c(0, x_max),
        xlab = 'Global depth (count)', ylab = 'Cumulative number of sites')
abline(v = q99.threshold, col='red', lty=2); abline(v = q01.threshold, col='red', lty=2)
abline(h = q99, col='grey', lty=2); abline(h = q01, col='grey', lty=2)
mtext(paste('Zoomed in - ', outname, sep=''), side = 3, line = -2, outer = TRUE)
dev.off()