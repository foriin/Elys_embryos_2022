library(GenomicRanges)
library(dplyr)
library(dm3)
library(rtracklayer)
library(ggplot2)
library(gridExtra)

# Function to get damid signals around TSSs
tss.surr <- function(tss.gr, signal.gr, sig.col, vic = 5){
  ovs <- findOverlaps(tss.gr, signal.gr, ignore.strand = T)
  ovs <- ovs[ovs@to > vic]
  # print(mcols(signal.gr)[[sig.col]][1:10])
  # sapply(1:length(ovs), function(x){
  #   if(strand(tss.gr[ovs@from[x]]) %>% as.vector() == "-"){
  #     rev(mcols(signal.gr)[[sig.col]][(ovs@to[x] - vic):(ovs@to[x] + vic)])
  #   } else{
  #     mcols(signal.gr)[[sig.col]][(ovs@to[x] - vic):(ovs@to[x] + vic)]
  #   }
  # })
  
  # Don't take the strandness into account
  
  sapply(1:length(ovs), function(x){
    mcols(signal.gr)[[sig.col]][(ovs@to[x] - vic):(ovs@to[x] + vic)]
    
  })
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

load("RData/elys.npc.x.chrom.colours.RData", verbose = T)
load("RData/kc.and.nrn.lam.hmm.RData", verbose = T)
load("RData/emb.hmm3.RData", verbose = T)
load("RData/kc.lam.profile.ma.RData", verbose = T)
load("RData/s2.elys.chip.bin.RData", verbose = T)
load("RData/h3k27ac.RData", verbose = T)
load("RData/emb.elys.hmm3.x.nup98.npc.nuc.uq.RData", verb = T)
seqlevels(e.elys.x.nup.npc.gr) <- sub("chr", "", seqlevels(e.elys.x.nup.npc.gr))
e.elys.x.nup.npc.gr <- reduce(e.elys.x.nup.npc.gr, min.gapwidth = 901)


# Let's make profiles of Kc167 Lamin DamID around centers of embryonic Elys
# domains that overlap with Kc167 NPC domains, LADs and inactive chromatin
# colours

# 1. Embryonic Elys domains that overlap w/ Kc NPC and LADs

domains.1 <- subsetByOverlaps(e.elys.x.nup.npc.gr[width(e.elys.x.nup.npc.gr) > 100],
                              kc.lam.hmm,
                              maxgap = 2001,
                              ignore.strand = T)
summary(width(domains.1))
# 
domains.1.c <- resize(domains.1, 1,fix = "center")

profs.1 <- list(
  kc.bin.lam.pr,
  lam.emb.pr,
  elys.emb.pr,
  s2.elys.chip.bins.2,
  h3k27ac.chip.bin
)
profs.1[[1]]@metadata <- list(
  protein = "Lam",
  title = "Emb Elys x Nup98 NPC* x Kc167 LADs\nLam DamID in Kc167",
  y_lab = "log2(Dam/Dam-Lam)",
  col = "#660000"
)
profs.1[[2]]@metadata <- list(
  protein = "Lam",
  title = "Emb Elys x Nup98 NPC* x Kc167 LADs\nLam DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Lam)",
  col = "#660000"
)
profs.1[[3]]@metadata <- list(
  protein = "Elys",
  title = "Emb Elys x Nup98 NPC* x Kc167 LADs\nElys DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Elys)",
  col = "red"
)
profs.1[[4]]@metadata <- list(
  protein = "Elys",
  title = "Emb Elys x Nup98 NPC* x Kc167 LADs\nElys ChIP-seq in S2",
  y_lab = "ChIP-seq Elys",
  col = "yellow"
)
profs.1[[5]]@metadata <- list(
  protein = "H3K27Ac",
  title = "Emb Elys x Nup98 NPC* x Kc167 LADs\nH3K27Ac ChIP-seq in S2",
  y_lab = "ChIP-seq H3K27Ac/Input",
  col = "orange"
)

ps.1 <- lapply(profs.1, function(pr){
  prof  <- tss.surr(domains.1.c, pr, vic = 20, sig.col = "score")
  dom.prof <- prof %>% rowMeans.2()  %>% 
    as.data.frame()  %>% 
    mutate(kb = seq(-6, by = 0.3, length.out = 41))
  dom.prof <- cbind(dom.prof, 
                      apply(prof, 1, mean.conf) %>% t %>% 
                        as.data.frame() )
  
  
  ggplot(dom.prof %>% 
                 setNames(c("score", "kb", "ymin", "ymax")),
               aes(x = kb, y = score))+
    # geom_line()+
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70", alpha = 0.4)+
    stat_smooth(geom = "line",
                method = "loess",alpha = 0.9,se = T, span = 1/5, size = 1.5,
                col = pr@metadata$col)+
    ylab(pr@metadata$y_lab)+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_x_continuous(breaks = seq(-6, 6, 1), expand = c(0.01,0.01))+
    ggtitle(pr@metadata$title)
})


# 2. Embryonic Elys domains that overlap w/ Kc NPC and embryonic LADs
# + S2 Elys Chip profiles
load("RData/s2.elys.chip.bin.RData", verbose = T)
domains.2 <- subsetByOverlaps(e.elys.x.nup.npc.gr[width(e.elys.x.nup.npc.gr) > 100],
lam.emb.hmm3,
maxgap = 2001,
ignore.strand = T)
summary(width(domains.2))

domains.2.c <- resize(domains.2, 1,fix = "center")

profs.2 <- list(
  kc.bin.lam.pr,
  lam.emb.pr,
  elys.emb.pr,
  s2.elys.chip.bins.2
)
profs.2[[1]]@metadata <- list(
  protein = "Lam",
  title = "Emb Elys x Nup98 NPC x embryo LADs\nLam DamID in Kc167",
  y_lab = "log2(Dam/Dam-Lam)",
  col = "#660000"
)
profs.2[[2]]@metadata <- list(
  protein = "Lam",
  title = "Emb Elys x Nup98 NPC x embryo LADs\nLam DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Lam)",
  col = "#660000"
)
profs.2[[3]]@metadata <- list(
  protein = "Elys",
  title = "Emb Elys x Nup98 NPC x embryo LADs\nElys DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Elys)",
  col = "#EE0000"
)
profs.2[[4]]@metadata <- list(
  protein = "Elys",
  title = "Emb Elys x Nup98 NPC x embryo LADs\nElys ChIP-seq in S2",
  y_lab = "ChIP-seq Elys/Input",
  col = "#EE3223"
)


ps.2 <- lapply(profs.2, function(pr){
  prof  <- tss.surr(domains.2.c, pr, vic = 20, sig.col = "score")
  dom.prof <- prof %>% rowMeans.2()  %>% 
    as.data.frame()  %>% 
    mutate(kb = seq(-6, by = 0.3, length.out = 41))
  dom.prof <- cbind(dom.prof, 
                    apply(prof, 1, mean.conf) %>% t %>% 
                      as.data.frame() )
  
  
  ggplot(dom.prof %>% 
           setNames(c("score", "kb", "ymin", "ymax")),
         aes(x = kb, y = score))+
    # geom_line()+
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70", alpha = 0.4)+
    stat_smooth(geom = "line",
                method = "loess",alpha = 0.9,se = T, span = 1/5, size = 1.5,
                col = pr@metadata$col)+
    ylab(pr@metadata$y_lab)+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_x_continuous(breaks = seq(-6, 6, 1), expand = c(0.01,0.01))+
    ggtitle(pr@metadata$title)
})


pdf("plots/Fig3b.pdf", width = 19)
grid.arrange(ps.1[[1]],
             ps.2[[2]],
             ps.1[[3]],
             ps.1[[4]],
             ps.1[[5]],
             nrow = 1)
dev.off()

# (Embryonic Elys x Nup98 NUC x Kc167 LADs)
domains.3 <- subsetByOverlaps(elys.x.nup.nuc[width(elys.x.nup.nuc) > 100],
                              kc.lam.hmm,
                              maxgap = 2001,
                              ignore.strand = T)
summary(width(domains.3))
# 
domains.3.c <- resize(domains.3, 1,fix = "center")

profs.1[[1]]@metadata <- list(
  protein = "Lam",
  title = "Emb Elys x Nup98 NUC* x Kc167 LADs\nLam DamID in Kc167",
  y_lab = "log2(Dam/Dam-Lam)",
  col = "#660000"
)
profs.1[[2]]@metadata <- list(
  protein = "Lam",
  title = "Emb Elys x Nup98 NUC* x Kc167 LADs\nLam DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Lam)",
  col = "#660000"
)
profs.1[[3]]@metadata <- list(
  protein = "Elys",
  title = "Emb Elys x Nup98 NUC* x Kc167 LADs\nElys DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Elys)",
  col = "red"
)
profs.1[[4]]@metadata <- list(
  protein = "Elys",
  title = "Emb Elys x Nup98 NUC* x Kc167 LADs\nElys ChIP-seq in S2",
  y_lab = "ChIP-seq Elys",
  col = "yellow"
)
profs.1[[5]]@metadata <- list(
  protein = "H3K27Ac",
  title = "Emb Elys x Nup98 NUC* x Kc167 LADs\nH3K27Ac ChIP-seq in S2",
  y_lab = "ChIP-seq H3K27Ac/Input",
  col = "orange"
)

ps.3 <- lapply(profs.1, function(pr){
  prof  <- tss.surr(domains.3.c, pr, vic = 20, sig.col = "score")
  dom.prof <- prof %>% rowMeans.2()  %>% 
    as.data.frame()  %>% 
    mutate(kb = seq(-6, by = 0.3, length.out = 41))
  dom.prof <- cbind(dom.prof, 
                    apply(prof, 1, mean.conf) %>% t %>% 
                      as.data.frame() )
  
  
  ggplot(dom.prof %>% 
           setNames(c("score", "kb", "ymin", "ymax")),
         aes(x = kb, y = score))+
    # geom_line()+
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70", alpha = 0.4)+
    stat_smooth(geom = "line",
                method = "loess",alpha = 0.9,se = T, span = 1/5, size = 1.5,
                col = pr@metadata$col)+
    ylab(pr@metadata$y_lab)+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_x_continuous(breaks = seq(-6, 6, 1), expand = c(0.01,0.01))+
    ggtitle(pr@metadata$title)
})

# (Embryonic Elys x Nup98 NUC x embryo LADs)
domains.4 <- subsetByOverlaps(elys.x.nup.nuc[width(elys.x.nup.nuc) > 100],
                              lam.emb.hmm3,
                              maxgap = 2001,
                              ignore.strand = T)
summary(width(domains.4))
# 
domains.4.c <- resize(domains.4, 1,fix = "center")

profs.2[[1]]@metadata <- list(
  protein = "Lam",
  title = "Emb Elys x Nup98 NUC x embryo LADs\nLam DamID in Kc167",
  y_lab = "log2(Dam/Dam-Lam)",
  col = "#660000"
)
profs.2[[2]]@metadata <- list(
  protein = "Lam",
  title = "Emb Elys x Nup98 NUC x embryo LADs\nLam DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Lam)",
  col = "#660000"
)
profs.2[[3]]@metadata <- list(
  protein = "Elys",
  title = "Emb Elys x Nup98 NUC x embryo LADs\nElys DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Elys)",
  col = "#EE0000"
)
profs.2[[4]]@metadata <- list(
  protein = "Elys",
  title = "Emb Elys x Nup98 NUC x embryo LADs\nElys ChIP-seq in S2",
  y_lab = "ChIP-seq Elys/Input",
  col = "#EE3223"
)


ps.4 <- lapply(profs.2, function(pr){
  prof  <- tss.surr(domains.4.c, pr, vic = 20, sig.col = "score")
  dom.prof <- prof %>% rowMeans.2()  %>% 
    as.data.frame()  %>% 
    mutate(kb = seq(-6, by = 0.3, length.out = 41))
  dom.prof <- cbind(dom.prof, 
                    apply(prof, 1, mean.conf) %>% t %>% 
                      as.data.frame() )
  
  
  ggplot(dom.prof %>% 
           setNames(c("score", "kb", "ymin", "ymax")),
         aes(x = kb, y = score))+
    # geom_line()+
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70", alpha = 0.4)+
    stat_smooth(geom = "line",
                method = "loess",alpha = 0.9,se = T, span = 1/5, size = 1.5,
                col = pr@metadata$col)+
    ylab(pr@metadata$y_lab)+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_x_continuous(breaks = seq(-6, 6, 1), expand = c(0.01,0.01))+
    ggtitle(pr@metadata$title)
})

pdf("plots/Fig3c.pdf", width = 19)
grid.arrange(ps.3[[1]],
             ps.4[[2]],
             ps.3[[3]],
             ps.3[[4]],
             ps.3[[5]],
             nrow = 1)
dev.off()

