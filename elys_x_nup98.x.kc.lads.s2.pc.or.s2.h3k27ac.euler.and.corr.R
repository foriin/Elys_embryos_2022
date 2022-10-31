library(dm3)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(tibble)
library(eulerr)

load("RData/h3k27ac.RData", verbose = T)
load("RData/s2.pc.chip.modencode.RData", verbose = T)
load("RData/nups.x.elys.emb.RData", verbose = T)
load("RData/kc.and.nrn.lam.hmm.RData", verbose = T)
load("RData/emb.hmm3.RData", verbose = T)

# Let's find how much do emb Elys x Kc Nup98 (NPC or nucleoplasmic) domains
# intersect with some other domains: Kc167 LADs, S2 Polycomb or S2 H3K27Ac

doms <- list(
  kc.lam.hmm,
  s2.pc.domains,
  h3k27ac.hmm3
)
names(doms) <- c("Kc LADs", "S2 Pc", "S2 H3K27Ac")
# 1. Elys x Nup98 NPC
ix <- sapply(1:3, function(x){
  x.no.y.z <- GenomicRanges::setdiff(
    doms[[x]],
    c(doms[-x][[1]], doms[-x][[2]]),
    ignore.strand = T
  )
  x.with.y.no.z <- GenomicRanges::setdiff(
    GenomicRanges::intersect(
      doms[-x][[1]],
      doms[[x]],
      ignore.strand = T
    ),
    doms[-x][[2]], ignore.strand = T
  )
  x.with.z.no.y <- GenomicRanges::setdiff(
    GenomicRanges::intersect(
      doms[-x][[2]],
      doms[[x]],
      ignore.strand = T
    ),
    doms[-x][[1]], ignore.strand = T
  )
  
  c(
    sum(width(GenomicRanges::intersect(nup.npc.x.elys.emb.1,
                                       x.no.y.z, ignore.strand = T))),
    sum(width(GenomicRanges::intersect(nup.npc.x.elys.emb.1,
                                       x.with.y.no.z, ignore.strand = T))),
    sum(width(GenomicRanges::intersect(nup.npc.x.elys.emb.1,
                                       x.with.z.no.y, ignore.strand = T)))
  )
})

colnames(ix) <- c("Kc LADs", "S2 Pc", "S2 H3K27Ac")

ix.2 <- rbind(ix, colSums(ix))/sum(width(nup.npc.x.elys.emb.1))*100
ix.all <- sum(width(GenomicRanges::intersect(
  nup.npc.x.elys.emb.1,
  Reduce(GenomicRanges::intersect, doms), ignore.strand = T
)))/sum(width(nup.npc.x.elys.emb.1))*100

fit <- euler(c("Kc_LADs" = unname(ix.2[1,1]),
               "S2_H3K27Ac" = unname(ix.2[1,3]),
               "Kc_LADs&S2_H3K27Ac" = unname(ix.2[3,1]),
               "S2_Pc" = unname(ix.2[1,2]),
               "S2_H3K27Ac&S2_Pc" = unname(ix.2[3,3]),
               "S2_Pc&Kc_LADs" = unname(ix.2[2.2]),
               "Kc_LADs&S2_H3K27Ac&S2_Pc" = ix.all))

fit$original.values <- round(fit$original.values, digits = 2)
# Save euler plot
pdf("plots/elys.npc.x.diff.doms.venn.pdf")
plot(fit, fill = c("blue", "red", "khaki"), quantities = T, main = "Percentage from emb Elys x Nup98 NPC")
dev.off()

# Elys x Nup98 Nucleoplasmic

ix.nuc <- sapply(1:3, function(x){
  x.no.y.z <- GenomicRanges::setdiff(
    doms[[x]],
    c(doms[-x][[1]], doms[-x][[2]]),
    ignore.strand = T
  )
  x.with.y.no.z <- GenomicRanges::setdiff(
    GenomicRanges::intersect(
      doms[-x][[1]],
      doms[[x]],
      ignore.strand = T
    ),
    doms[-x][[2]], ignore.strand = T
  )
  x.with.z.no.y <- GenomicRanges::setdiff(
    GenomicRanges::intersect(
      doms[-x][[2]],
      doms[[x]],
      ignore.strand = T
    ),
    doms[-x][[1]], ignore.strand = T
  )
  
  c(
    sum(width(GenomicRanges::intersect(nup.nuc.x.elys.emb.1,
                                       x.no.y.z, ignore.strand = T))),
    sum(width(GenomicRanges::intersect(nup.nuc.x.elys.emb.1,
                                       x.with.y.no.z, ignore.strand = T))),
    sum(width(GenomicRanges::intersect(nup.nuc.x.elys.emb.1,
                                       x.with.z.no.y, ignore.strand = T)))
  )
})

colnames(ix.nuc) <- c("Kc LADs", "S2 Pc", "S2 H3K27Ac")
ix.nuc.all <- sum(width(GenomicRanges::intersect(
  nup.nuc.x.elys.emb.1,
  Reduce(GenomicRanges::intersect, doms), ignore.strand = T
)))/sum(width(nup.nuc.x.elys.emb.1))*100
ix.nuc <- rbind(ix.nuc/sum(width(nup.nuc.x.elys.emb.1))*100,
                colSums(ix.nuc)/sum(width(nup.nuc.x.elys.emb.1))*100 + ix.nuc.all)


fit <- euler(c("Kc_LADs" = unname(ix.nuc[1,1]),
               "S2_H3K27Ac" = unname(ix.nuc[1,3]),
               "Kc_LADs&S2_H3K27Ac" = unname(ix.nuc[3,1]),
               "S2_Pc" = unname(ix.nuc[1,2]),
               "S2_H3K27Ac&S2_Pc" = unname(ix.nuc[3,3]),
               "S2_Pc&Kc_LADs" = unname(ix.nuc[2.2]),
               "Kc_LADs&S2_H3K27Ac&S2_Pc" = ix.nuc.all))

fit$original.values <- round(fit$original.values, digits = 2)
# Save euler plot
pdf("plots/elys.nuc.x.diff.doms.venn.pdf")
plot(fit, fill = c("blue", "red", "khaki"), quantities = T, main = "Percentage from emb Elys x Nup98 NUC")
dev.off()


# Correlation of Elys profiles with H3K27Ac in Elys x Nup NUC

e.e.nup.no.npc.pr <- elys.emb.pr %>% subsetByOverlaps(nup.npc.x.elys.emb.1,
                                                      invert = T)
e.e.nup.no.npc.pr <- e.e.nup.no.npc.pr[!is.na(e.e.nup.no.npc.pr$score)]
s2.h3k27ac.nup.nuc.pr <- h3k27ac.chip.bin %>% subsetByOverlaps(e.e.nup.no.npc.pr)

cor.test(e.e.nup.no.npc.pr$score, s2.h3k27ac.nup.nuc.pr$score)
cor(e.e.nup.no.npc.pr$score, s2.h3k27ac.nup.nuc.pr$score, use = "pairwise")

e.e.nup.npc.pr <- elys.emb.pr %>% subsetByOverlaps(nup.npc.x.elys.emb.1)
s2.h3k27ac.nup.npc.pr <- h3k27ac.chip.bin %>% subsetByOverlaps(nup.npc.x.elys.emb.1)

cor.test(e.e.nup.npc.pr$score, s2.h3k27ac.nup.npc.pr$score)
cor(e.e.nup.npc.pr$score, s2.h3k27ac.nup.npc.pr$score, use = "pairwise")

# It's not very high!