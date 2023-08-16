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
  auc1 = plot.cdc(cdc.list, 2, 2, file= paste0("pdf/", filenamebase, "-1.95percentile.pdf"), factor=tf)
  auc2 = plot.cdc(cdc.list, 2, 2, file= paste0("pdf/", filenamebase, "-2.median.pdf"), factor=tf, ylim = c(-0.1, 0.5), do.med = T)
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

for(tf in ChIPnames2) {
  res.in = intersect.tfcor(tf, paste0("Fig3/ChIPall/intersect-", tf))
  res.ol = overlap.tfcor(tf, paste0("Fig3/ChIPall/overlap-", tf))
  AUC.in = rbind(AUC.in, data.frame(TF = tf, AUC.95 = res.in[1], AUC.med = res.in[2]))
  AUC.ol = rbind(AUC.ol, data.frame(TF = tf, AUC.95 = res.ol[1], AUC.med = res.ol[2]))
}


#######################################
# Significant correlation percentage using fdr values

# Function to return all PRO-cap correlations and their fdr values
fdr.intersect = function(window, pos, ref.cor, max.cross = 5) {
  window.in = list()
  window.ex = list()
  for(i in 1:max.cross) {
    window.in[[i]] = exclude.window(window, pro[,1:2], pro[,1:2], pos, n.cross=i)
    window.ex[[i]] = diff.window(window, window.in[[i]])
  }
  pro.in.cor = list()
  pro.ex.cor = list()
  for(i in 1:max.cross) {
    pro.in.cor[[i]] = recor.pro(window.in[[i]], window, ref.cor)
    pro.ex.cor[[i]] = recor.pro(window.ex[[i]], window, ref.cor)
  } 
  # retrieve just the distance, correlation coefficient, pvalue
  fdr=p.adjust(cor2pvalue(ref.cor$y, 76), method = "fdr")
  pfunc = approxfun(ref.cor$y, fdr)
  
  fdr.in.res = list()
  for(i in 1:max.cross) {
    fdr.in.res[[i]] = data.frame(
      pos = pro.in.cor[[i]]$x,
      cor = pro.in.cor[[i]]$y,
      fdr = pfunc(pro.in.cor[[i]]$y))
  }
  return(fdr.in.res)
}


fdr.overlap = function(window, pos, ref.cor) {
  window.o1 = overlap.window(window, pro[,1:2], pro[,1:2], pos, n.over=1)
  window.o2 = overlap.window(window, pro[,1:2], pro[,1:2], pos, n.over=2)
  window.o0 = diff.window(window, window.o1)

  fdr=p.adjust(cor2pvalue(ref.cor$y, 76), method = "fdr")
  pfunc = approxfun(ref.cor$y, fdr)
  
  pro.ol.cor = list()
  pro.ol.cor[[1]] = recor.pro(window.o0, window, ref.cor)
  pro.ol.cor[[2]] = recor.pro(window.o1, window, ref.cor)
  pro.ol.cor[[3]] = recor.pro(window.o2, window, ref.cor)
  
  fdr.ol.res = list()
  for(i in 1:3) {
    fdr.ol.res[[i]] = data.frame(
      pos = pro.ol.cor[[i]]$x,
      cor = pro.ol.cor[[i]]$y,
      fdr = pfunc(pro.ol.cor[[i]]$y))
  }
  return(fdr.ol.res)
}

fdr.in.merge = function(fdr.in, cdc.in=1, cdc.ex=2) {
  fdr.in.list = 
    Reduce(rbind, fdr.in[1:cdc.in])
  fdr.ex.list = 
    Reduce(rbind, fdr.in[cdc.ex:length(fdr.in)])
  return(list(fdr.in.list, fdr.ex.list))
}

plot.fdr.dist2 = function(fdr.list, file="pdf/test.pdf",
                          xlim=c(0,1000*2^(0:8)), ylim = c(0, 1),
                          cdc.in = 1, cdc.ex = 2,
                          fdr = 0.1, main = "All tTRE", width = 1) {
  ndata1 = split(fdr.list[[1]]$fdr,cut(fdr.list[[1]]$pos/1000,xlim/1000))
  f1 = unlist(lapply(ndata1, function(x) sum(x<fdr)/length(x))) * 100
  n1 = unlist(lapply(ndata1, function(x) length(x)))
  ndata2 = split(fdr.list[[2]]$fdr,cut(fdr.list[[2]]$pos/1000,xlim/1000))
  f2 = unlist(lapply(ndata2, function(x) sum(x<fdr)/length(x))) * 100
  n2 = unlist(lapply(ndata2, function(x) length(x)))
  nxlim = length(xlim)-1
  labels = paste(xlim[1:nxlim]/1000,xlim[1:nxlim+1]/1000,sep="-")
  labels = factor(labels, levels = labels)
  plotdata = data.frame(Percentage = f1,
                        Distance = labels,
                        n = n1,
                        class = 1) %>%
    bind_rows(data.frame(Percentage = f2,
                         Distance = labels,
                         n = n2,
                         class = 2))
  
  g = ggplot(plotdata, aes(y = Percentage, x = Distance, fill = factor(class))) +
    geom_bar(stat='identity', position='dodge2', color = "black") +
    theme_bw() +
    scale_fill_manual(values = c("#ffffff", "#e0e0e0")) +
    theme(legend.position = "none") +
    ggtitle(main)

  return(list(g, plotdata))
}

pdf("pdf/Fig3/ChIPfdr/all.pdf")
fdr.data.all = data.frame(Percentage = c(),
                          Distance = c(),
                          n = c(),
                          class = c(),
                          type = c(),
                          factor = c())
for(tf in ChIPnames) {
  fdr.in.tf = fdr.in.merge(fdr.intersect(window, ct[[tf]], pro.cor))
  fdr.ol.tf = fdr.overlap(window, ct[[tf]], pro.cor)
  g1 = plot.fdr.dist2(fdr.in.tf, main = paste(tf, "intersect"))
  g2 = plot.fdr.dist2(fdr.ol.tf[c(1,2)], main = paste(tf, "overlap"))
  print(g1[[1]])
  print(g2[[1]])
  fdr.data.all = fdr.data.all %>% 
    bind_rows(data.frame(g1[[2]], type = "intersect", factor = tf)) %>%
    bind_rows(data.frame(g2[[2]], type = "overlap", factor = tf))
}
dev.off()

pdf("pdf/Fig3/FigS3.select.pdf", width = 3, height = 3)
ChIPnames2 = c("CTCF", "P300", "RAD21")
fdr.data.all = data.frame(Percentage = c(),
                          Distance = c(),
                          n = c(),
                          class = c(),
                          type = c(),
                          factor = c())
for(tf in ChIPnames2) {
  fdr.in.tf = fdr.in.merge(fdr.intersect(window, ct[[tf]][ ,-1], pro.cor))
  fdr.ol.tf = fdr.overlap(window, ct[[tf]][,-1], pro.cor)
  g1 = plot.fdr.dist2(fdr.in.tf, main = paste(tf, "intersect"))
  g2 = plot.fdr.dist2(fdr.ol.tf[c(1,2)], main = paste(tf, "overlap"))
  print(g1[[1]])
  print(g2[[1]])
  fdr.data.all = fdr.data.all %>% 
    bind_rows(data.frame(g1[[2]], type = "intersect", factor = tf)) %>%
    bind_rows(data.frame(g2[[2]], type = "overlap", factor = tf))
}
dev.off()

pdf("pdf/Fig3/FigS3.select.pdf", width = 4, height = 4)
for(tf in ChIPnames2) {
  plotdata = fdr.data.all %>%
    filter(factor == tf)
  paste(tf, "intersect")
  g1 = ggplot(plotdata %>% filter(type == "intersect"),
             aes(y = Percentage, x = Distance, color = factor(class))) +
    geom_bar(stat='identity', position='dodge2', fill = "white", size = 0.5) +
    scale_color_manual(values = color$adjust(c("#084FA0", "#E09810"), 1, 1)) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.75),
          legend.position = "none",
          axis.line = element_line(linewidth = 0),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, color = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.text.y=element_text(color="black")) +
    ggtitle(paste(tf, "intersection"))
  g2 = ggplot(plotdata %>% filter(type == "overlap"),
              aes(y = Percentage, x = Distance, color = factor(class))) +
    geom_bar(stat='identity', position='dodge2', fill = "white", size = 0.5) +
    scale_color_manual(values = color$adjust(c("#084FA0", "#E09810"), 1, 1)) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.75),
          legend.position = "none",
          axis.line = element_line(linewidth = 0),
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, color = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.text.y=element_text(color="black")) +
    ggtitle(paste(tf, "occupancy")) 
  print(g1)
  print(g2)
}
dev.off()

