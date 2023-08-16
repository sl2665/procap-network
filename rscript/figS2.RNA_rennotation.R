ranot = list()
ranot$rnaseq.geneloc = read.table("window/rnaseq.geneloc.txt", header = F,
                                  col.names = c("id", "name",
                                                "chr", "pos1", "pos2", "strand"))
ranot$rnaseq.pos = read.table("window/rnaseq.pos.old.bed", col.names = c("chr", "pos1", "pos2", "b1", "b2", "strand"))

ranot$rnaseq.pos = ranot$rnaseq.pos %>%
  left_join(ranot$rnaseq.geneloc, by = c("chr", "pos1", "pos2", "strand"))

# Retrive procap positions (all pro-cap, not just variable ones)
ranot$procap.pos = read.table("readcount/procap.ambr.norm.txt", header = F)[, 1:2]
colnames(ranot$procap.pos) = c("chr", "pos")
ranot$procap.npos = sort(conv.pos(ranot$procap.pos[,1:2]))

###
# Define active TSSs
ranot$ensG = read.table("window/ensGene-hg19.txt")[,c(3,4,5,6,13)]
colnames(ranot$ensG) = c("chr", "strand", "start", "end", "id")
ranot$ensG = ranot$ensG %>%
  mutate(ntss = ifelse(strand == "+", conv.pos(ranot$ensG[,c(1, 3)]),
                       conv.pos(ranot$ensG[,c(1, 4)])))
# find nearest upstream procap and downstream procap positions to find nearest procap positions
ranot$idx = findInterval(ranot$ensG$ntss, ranot$procap.npos)
ranot$uidx = ranot$idx
ranot$uidx[ranot$idx == 0] = 1
ranot$didx = ranot$idx + 1
ranot$didx[ranot$idx == length(ranot$procap.npos)] = length(ranot$procap.npos)
ranot$distup = ranot$ensG$ntss - ranot$procap.npos[ranot$uidx]
ranot$distdn = ranot$procap.npos[ranot$didx] - ranot$ensG$ntss
ranot$aTSS = abs(ranot$distup) < 1000 | abs(ranot$distdn) < 1000
ranot$ensG_a = ranot$ensG[ranot$aTSS, ] %>%
  mutate(distup = ranot$distup[ranot$aTSS],
         distdn = ranot$distdn[ranot$aTSS]) %>%
  mutate(dist = ifelse(distup > distdn, distdn, distup)) %>%
  group_by(id) %>%
  filter(dist == min(dist)) %>%
  ungroup()

View(ranot$ensG_a)
View(ranot$rnaseq.pos)

# Next steps 
ranot$rnaseq.ensG = ranot$rnaseq.pos %>%
  inner_join(ranot$ensG_a, by = "id") %>%
  select(-start, -end) %>%
  distinct()
View(ranot$rnaseq.ensG)
View(ranot$rnaseq.ensG)

ranot$rnaseq.pos = read.table("window/rnaseq.pos.old.bed", col.names = c("chr", "pos1", "pos2", "b1", "b2", "strand"))
View(ranot$rnaseq.pos)

ranot$rnaseq.pos = ranot$rnaseq.pos %>%
  left_join(ranot$rnaseq.ensG %>%
              mutate(chr = chr.x) %>%
              select(chr, pos1, ntss) %>%
              distinct(chr, pos1, .keep_all = T),
            by = c("chr", "pos1")) %>%
  mutate(pos1 = ntss %% 1000000000, pos2 = ntss %% 1000000000 + 1) %>%
  select(-ntss)
ranot$rnaseq.pos[is.na(ranot$rnaseq.pos)] = -1000000
View(ranot$rnaseq.pos)

write.table(ranot$rnaseq.pos, "window/rnaseq.pos.bed", quote = F, sep = "\t", col.names = F, row.names = F)
