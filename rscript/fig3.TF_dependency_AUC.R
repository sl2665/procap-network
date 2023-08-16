source('rscript/scatterdist.R')

# Plot correlation scatterplot along the distance
window = make.window.list(pro.near)
pro.cor  = cor.pro(window, pro, pro)

# Correlation decay curve

# Test the effect of CTCF intersection on the CDC
ChIPnames

intersect.tfcor = function(tf, filenamebase) {
  cdc.list = cdc.intersect(window, ct[[tf]], pro.cor)
  paste0("pdf/", filenamebase, "-1.pdf")
  auc1 = plot.cdc(cdc.list, 1, 2, file= paste0("pdf/", filenamebase, "-1.95percentile.pdf"), factor=tf)
  auc2 = plot.cdc(cdc.list, 1, 2, file= paste0("pdf/", filenamebase, "-2.median.pdf"), factor=tf, ylim = c(-0.1, 0.5), do.med = T)
  return(c(auc1[2]-auc1[1], auc2[2]-auc2[1]))
}

overlap.tfcor = function(tf, filenamebase) {
  cdc.list = cdc.overlap(window, ct[[tf]], pro.cor)
  paste0("pdf/", filenamebase, "-1.pdf")
  auc1 = plot.cdc.over(cdc.list, file= paste0("pdf/", filenamebase, "-1.95percentile.pdf"), factor=tf)
  auc2 = plot.cdc.over(cdc.list, file= paste0("pdf/", filenamebase, "-2.median.pdf"), factor=tf, ylim = c(-0.1, 0.5), do.med = T)
  return(c(auc1[2]-auc1[1], auc2[2]-auc2[1]))
}

# Generate intersect and overlap plots for all ChIP factors
AUC.in = data.frame(TF = c(), AUC.95 = c(), AUC.med = c())
AUC.ol = data.frame(TF = c(), AUC.95 = c(), AUC.med = c())

for(tf in ChIPnames) {
  res.in = intersect.tfcor(tf, paste0("Fig3/ChIPall/intersect-", tf))
  res.ol = overlap.tfcor(tf, paste0("Fig3/ChIPall/overlap-", tf))
  AUC.in = rbind(AUC.in, data.frame(TF = tf, AUC.95 = res.in[1], AUC.med = res.in[2]))
  AUC.ol = rbind(AUC.ol, data.frame(TF = tf, AUC.95 = res.ol[1], AUC.med = res.ol[2]))
}


