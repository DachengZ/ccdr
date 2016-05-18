time <- c(test.50.50$time / 50, test.50.100$time / 50,
          test.100.100$time / 50, test.100.150$time / 50,
          test.100.200$time / 50,
          test.200.100$time / 20, test.200.200$time / 20,
          test.200.300$time / 20, test.200.400$time / 5)
pdf(file = "time.pdf", width = 7, height = 7)
    par(oma = c(0, 0, 0, 0))
    plot(time, xlab = "# of nodes, # of edges", ylab = "time (seconds)", main = "average running time", xaxt = "n")
    text(1:length(time), time + 0.3, labels = c("50, 50", "50, 100", "100, 100", "100, 150", "100, 200", "200, 100", "200, 200", "200, 300", "200, 400"), cex = 0.7)
dev.off()


tpr <- c(mean(test.50.50$metric[, 5]), mean(test.50.100$metric[, 5]),
         mean(test.100.100$metric[, 5]), mean(test.100.150$metric[, 5]),
         mean(test.100.200$metric[, 5]),
         mean(test.200.100$metric[, 5]), mean(test.200.200$metric[, 5]),
         mean(test.200.300$metric[, 5]), mean(test.200.400$metric[, 5])
         )
pdf(file = "tpr.pdf", width = 7, height = 7)
    par(oma = c(0, 0, 0, 0))
    plot(tpr, xlab = "# of nodes, # of edges", ylab = "TPR", main = "average TPR", xaxt = "n")
    text(1:length(tpr), tpr + 0.005, labels = c("50, 50", "50, 100", "100, 100", "100, 150", "100, 200", "200, 100", "200, 200", "200, 300", "200, 400"), cex = 0.7)
dev.off()

fdr <- c(mean(test.50.50$metric[, 6]), mean(test.50.100$metric[, 6]),
         mean(test.100.100$metric[, 6]), mean(test.100.150$metric[, 6]),
         mean(test.100.200$metric[, 6]),
         mean(test.200.100$metric[, 6]), mean(test.200.200$metric[, 6]),
         mean(test.200.300$metric[, 6]), mean(test.200.400$metric[, 6])
)
pdf(file = "fdr.pdf", width = 7, height = 7)
    par(oma = c(0, 0, 0, 0))
    plot(fdr, xlab = "# of nodes, # of edges", ylab = "fpr", main = "average FDR", xaxt = "n")
    text(1:length(fdr), fdr + 0.005, labels = c("50, 50", "50, 100", "100, 100", "100, 150",  "100, 200", "200, 100", "200, 200", "200, 300", "200, 400"), cex = 0.7)
dev.off()

