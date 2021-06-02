library(reshape2)
library(lattice)
library(gridExtra)

### Make sure to run "CME 250 Project Part 2 - psrin.R" first to populate all data fields ###

colors = c('deepskyblue','firebrick1')

#Graph balanced accuracy results

plot11 <- barchart(BalAcc ~ Resampling, data = results[c(1:8),c(2:5)], groups=Metric, ylim=c(0.5,1), main=" Logistic regression \n (glm)", par.settings=list(superpose.polygon=list(col=colors)), scales=list(x=list(rot=0,cex=0.8)))
plot12 <- barchart(BalAcc ~ Resampling, data = results[c(9:16),c(2:5)], groups=Metric, ylim=c(0.5,1), main=" Random forest \n (ranger)", par.settings=list(superpose.polygon=list(col=colors)), scales=list(x=list(rot=0,cex=0.8)))
plot13 <- barchart(BalAcc ~ Resampling, data = results[c(17:24),c(2:5)], groups=Metric, ylim=c(0.5,1), main=" Neural network \n (nnet)", par.settings=list(superpose.polygon=list(col=colors)), scales=list(x=list(rot=0,cex=0.8)))

#Graph computation time results

plot21 <- barchart(CompTime ~ Resampling, data = results[c(1:8),c(2:5)], groups=Metric, ylim=c(0,40), auto.key=list(space='bottom'), par.settings=list(superpose.polygon=list(col=colors)), scales=list(x=list(rot=0,cex=0.8)))
plot22 <- barchart(CompTime ~ Resampling, data = results[c(9:16),c(2:5)], groups=Metric, ylim=c(0,40), auto.key=list(space='bottom'), par.settings=list(superpose.polygon=list(col=colors)), scales=list(x=list(rot=0,cex=0.8)))
plot23 <- barchart(CompTime ~ Resampling, data = results[c(17:24),c(2:5)], groups=Metric, ylim=c(0,40), auto.key=list(space='bottom'), par.settings=list(superpose.polygon=list(col=colors)), scales=list(x=list(rot=0,cex=0.8)))

grid.arrange(plot11,plot12,plot13,plot21,plot22,plot23, nrow = 2, ncol=3)

#Visualize confusion matrices for top model of each model type (glm, ranger, nnet)
top_glm <- which.max(results$BalAcc[9:16])              #Index of top hyperparameters for glm
top_ranger <- which.max(results$BalAcc[9:16]) + 8       #Index of top hyperparameters for ranger
top_nnet <- which.max(results$BalAcc[17:24]) + 16       #Index of top hyperparameters for nnet

par(mfrow=c(1,3))

fourfoldplot(ConfMatrices[[top_glm]], color = c("#CC6666", "#99CC99"), conf.level = 0, margin = 1)
fourfoldplot(ConfMatrices[[top_ranger]], color = c("#CC6666", "#99CC99"), conf.level = 0, margin = 1)
fourfoldplot(ConfMatrices[[top_nnet]], color = c("#CC6666", "#99CC99"), conf.level = 0, margin = 1)

mtext(" Logistic regression \n (glm) \n ROC + smote", at=-5.2, line=-6, font=2)
mtext(" Random forest \n (ranger) \n ROC + smote", at=-2.6, line=-6, font=2)
mtext(" Neural network \n (nnet) \n ROC + up", at=0, line=-6, font=2)