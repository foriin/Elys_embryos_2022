library(reshape2)
library(stringr)
library(dm3)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(seqinr)
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(Biostrings)
library(ggplot2)


load("RData/emb.lam.elys.nogap.RData", verbose = T)
load("RData/nup98.damid.hmm3.nogap.RData", verbose = T)
load("RData/emb.hmm3.RData", verb = T)


elys.x.nup.npc <- GenomicRanges::intersect(elys.emb.hmm3.ng,
                                           nup98.npc.hmm3.uq.ng,
                                           ignore.strand = T)
elys.x.nup.npc <- elys.x.nup.npc[width(elys.x.nup.npc) >= 100]

elys.x.nup.nuc <- GenomicRanges::intersect(elys.emb.hmm3.ng,
                                           nup98.nuc.hmm3.uq.ng,
                                           ignore.strand = T)
elys.x.nup.nuc <- elys.x.nup.nuc[width(elys.x.nup.nuc) >= 100]


# 20 bp bins
count_at.2 <- function(dom, area = 3000, bin = 20){
  a <- resize(dom, area*2 + bin, fix = 'center')
  if (!all(grepl("^chr", seqlevels(a)))) seqlevels(a) <-
      paste0("chr", seqlevels(a))
  seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3, a)
  b <- sapply(as.character(seqs), function(x){
    cutt <- sapply(seq(from=1, to=nchar(x), by=bin),
                   function(i) substr(x, i, i+(bin - 1)))
    1 - sapply(cutt, function(x){GC(s2c(x))}, USE.NAMES = F)
  }, USE.NAMES = F)
  
  return(b)
}

rowMeans.2 <- function(df, na.rm = T){
  apply(df, 1, function(x){
    x <- x[x != 0]
    mean(x, na.rm = T)
  })
}

mean.conf <- function(vec){
  vec <- vec[vec != 0 & !is.na(vec)]
  error <- qt(0.975,df=length(vec)-1)*sd(vec)/sqrt(length(vec))
  low <- mean(vec, na.rm = T) - error
  high <- mean(vec, na.rm = T) + error
  c(low,high)
}


at.1 <- count_at.2(elys.x.nup.npc)
at.1.df <- cbind(seq(-3, 3, 0.02),
                 rowMeans.2(at.1),
                 apply(at.1, 1, mean.conf) %>% t()) %>% as_tibble() %>% 
  setNames(c("kb", "y", "ymin", "ymax"))
at.2 <- count_at.2(elys.x.nup.nuc)
at.2.df <- cbind(seq(-3, 3, 0.02),
                 rowMeans.2(at.2),
                 apply(at.2, 1, mean.conf) %>% t()) %>% as_tibble() %>% 
  setNames(c("kb", "y", "ymin", "ymax"))
rranges <- elys.emb.pr[start(elys.emb.pr) > 6100] %>% delete.het.gr()
at.4 <- count_at.2(rranges[sample(1:length(rranges), 3000)])
at.4.df <- cbind(seq(-3, 3, 0.02),
                 rowMeans.2(at.4),
                 apply(at.4, 1, mean.conf) %>% t()) %>% as_tibble() %>% 
  setNames(c("kb", "y", "ymin", "ymax"))



p1 <- ggplot(at.1.df,
             aes(x = kb, y = y))+
  # geom_line()+
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "red3", alpha = 0.4)+
  stat_smooth(geom = "line",
              method = "loess",alpha = 0.9, se = T, span = 1/5, size = 1.2,
              col = "red")+
  ylab("AT-percent per 20 bp bin")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("Elys embryo x Kc167 Nup98 NPC*")

p2 <- ggplot(at.2.df,
             aes(x = kb, y = y))+
  # geom_line()+
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.4)+
  stat_smooth(geom = "line",
              method = "loess",alpha = 0.9, se = T, span = 1/5, size = 1.2,
              col = "blue4")+
  ylab("AT-percent per 20 bp bin")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("Elys embryo x Kc167 Nup98 NUC*")

p4 <- ggplot(at.4.df,
             aes(x = kb, y = y))+
  # geom_line()+
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "burlywood", alpha = 0.4)+
  stat_smooth(geom = "line",
              method = "loess",alpha = 0.9, se = T, span = 1/5, size = 1.2,
              col = "bisque4")+
  ylab("AT-percent per 20 bp bin")+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, 1), expand = c(0.01,0.01))+
  ggtitle("Random locations")

pdf("plots/at_content.elys.x.nup.npc.nuc.lads.20bp.pdf", width = 8, height = 8)
grid.arrange(
  p1 + ylim(c(0.55, 0.61)),
  p2 + ylim(c(0.55, 0.61)),
  p4 + ylim(c(0.55, 0.61)),
  ncol = 2
)
dev.off()

# Sliding window, very slow
count_at.3 <- function(dom){
  a <- resize(dom, 6100, fix = 'center')
  if (!all(grepl("^chr", seqlevels(a)))) seqlevels(a) <-
      paste0("chr", seqlevels(a))
  seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3, a)
  b <- sapply(as.character(seqs), function(x){
    cutt <- sapply(seq(from=1, to=nchar(x), by=1),
                   function(i) substr(x, i, i+99))
    1 - sapply(cutt, function(x){GC(s2c(x))}, USE.NAMES = F)
  }, USE.NAMES = F)
  
  c <- apply(b, 1, mean, na.rm = T)
  return(c)
}

domains.1 <- subsetByOverlaps(elys)

# Heat maps
library(genomation)
whitered <- colorRampPalette(c("white", "red3"))

at.1 <- count_at.2(elys.x.nup.npc, area = 3000, bin = 300)
dim(at.1)
at.1.df <- cbind(seq(-3, 3, 0.3),
                 rowMeans.2(at.1),
                 apply(at.1, 1, mean.conf) %>% t()) %>% as_tibble() %>% 
  setNames(c("kb", "y", "ymin", "ymax"))
at.2 <- count_at.2(elys.x.nup.nuc, area = 3000, bin = 300)
at.2.df <- cbind(seq(-3, 3, 0.3),
                 rowMeans.2(at.2),
                 apply(at.2, 1, mean.conf) %>% t()) %>% as_tibble() %>% 
  setNames(c("kb", "y", "ymin", "ymax"))
rranges <- elys.emb.pr[start(elys.emb.pr) > 6100] %>% delete.het.gr()
at.4 <- count_at.2(rranges[sample(1:length(rranges), 3000)],
                   area = 3000, bin = 300)
at.4.df <- cbind(seq(-3, 3, 0.3),
                 rowMeans.2(at.4),
                 apply(at.4, 1, mean.conf) %>% t()) %>% as_tibble() %>% 
  setNames(c("kb", "y", "ymin", "ymax"))

sml.1 <- ScoreMatrixBin(target = lam.spg.spc[[1]], windows = aly.dep.1,
                        strand.aware = F,weight.col = "log2damid", bin.num = nrow(at.1))
at.1a <- t(at.1)[(apply(at.1, 2, function(x){all(x > 0.4)}) &
                    apply(at.1, 2, function(x){all(x < 0.8)})),]
S3Part(sml.1) <- at.1a[order(at.1a[,11], decreasing = T),]

sml.2 <- ScoreMatrixBin(target = lam.spg.spc[[1]], windows = aly.dep.1,
                        strand.aware = F,weight.col = "log2damid", bin.num = nrow(at.2))
at.2a <- t(at.2)[(apply(at.2, 2, function(x){all(x > 0.4)}) &
                    apply(at.2, 2, function(x){all(x < 0.8)})),]
S3Part(sml.2) <- at.2a[order(at.2a[,11], decreasing = T),]

sml.3 <- ScoreMatrixBin(target = lam.spg.spc[[1]], windows = aly.dep.1,
                        strand.aware = F,weight.col = "log2damid", bin.num = nrow(at.4))
at.4a <- t(at.4)[apply(at.4, 2, function(x){all(x > 0.4)}) &
                   apply(at.4, 2, function(x){all(x < 0.8)}),]
S3Part(sml.3) <- at.4a[order(at.4a[,11], decreasing = T),]


pdf("plots/test.pdf")
par(mfrow = c(1,3))
heatMatrix(sml.1, order = F, xcoords = c(-3000, 3000), col = whitered(10), main = "Elys x NPC")
heatMatrix(sml.2, order = F, xcoords = c(-3000, 3000), col = whitered(10), main = "Elys x NUP")
heatMatrix(sml.3, order = F, xcoords = c(-3000, 3000), col = whitered(10), main = "Random")
dev.off()

heatMatrix(test, order = F, xcoords = c(-3000, 3000), col = whitered(30) )


save(elys.x.nup.npc,
     elys.x.nup.nuc,
     file = "RData/emb.elys.hmm3.x.nup98.npc.nuc.uq.RData")
