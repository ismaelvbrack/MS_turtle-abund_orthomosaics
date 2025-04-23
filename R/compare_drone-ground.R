
#### Compare the counts of drone-based orthomosaics and visual ground


library(ggplot2)
library(gridExtra)


#* Population counts
pop.counts <- read.table(file.path("data","PopulationCounts.txt"),h=T)

pop.counts$data <- factor(pop.counts$data, levels=pop.counts$data) # to dates

n.visits <- as.numeric(nrow(pop.counts)) # number of visits / days

# Organize count data
counts <- data.frame(
  date=rep(pop.counts$data, 3),
  drone=rep(pop.counts$walking + pop.counts$nesting, 3),
  ground=c(pop.counts$ground1,pop.counts$ground2,pop.counts$ground3)
)


# Fig S1a - Drone vs. Ground counts ---------------------------------------
rango <- range(counts$drone,counts$ground, na.rm=T)

fs1a <- ggplot(counts, aes(x=drone,y=ground)) +
  geom_point(size=2) +
  geom_smooth(method="lm",se=F) +
  scale_y_continuous(limits=rango) +
  scale_x_continuous(limits=rango) +
  geom_abline(intercept=0, slope=1, linewidth=1,col="gray50",linetype="dashed") +
  labs(x="Counts in the orthomosaic",y="Visual ground counts",title="a)") +
  theme_classic(base_size=16) 

# Proportion of the ground counts in comparison to drone-based counts
data.frame(
  drone=rowSums(pop.counts[,c("walking","nesting")]),
  ground=round(rowMeans(pop.counts[,c("ground1","ground2","ground3")],na.rm=T)),
  prop=paste0(round(rowMeans(pop.counts[,c("ground1","ground2","ground3")],na.rm=T) / rowSums(pop.counts[,c("walking","nesting")]) * 100,2), "%")
)

# Abundance estimates -----------------------------------------------------
samps1 <- read.table(file.path("outputs","MR5_postMCMC.txt"),h=T)

resu1 <- t(apply(samps1, 2, function(x)
  c(mean(x),quantile(x, probs=c(0.025,0.975)))
))

# Root folder with all files
fl.root <- file.path("outputs","Counts5_phi1_postMCMC")

# List files
fls <- list.files(fl.root,pattern="s")

# Import each file with the posterior samples
samps2 <- list()
for(i in 1:length(fls)){
  samps2[[i]] <- read.table(file.path(fl.root,fls[i]),h=T)
}

samps2 <- do.call(rbind, samps2) # combine all samples

length(fls)*200*3 == nrow(samps2) # check

# Summarize results
resu <- t(apply(samps2, 2, function(x) 
  c(Mean=mean(x),SE=sd(x),quantile(x, probs=c(0.025,0.975)))
))

resu <- as.data.frame(resu)

round(resu[13:24,])

counts$c.ground <- counts$ground / resu1["phi.1.",1]

counts$est.N <- round(resu[13:24,"Mean"])


# Fig S1b - Difference between drone and ground counts vs. estimated abundance -------------------------

#rango <- range(counts$est.N,counts$c.ground, na.rm=T)
fs1b <- ggplot(counts, aes(x=est.N,y=drone - ground)) +
  geom_point(size=2) +
  #scale_y_continuous(limits=rango) +
  #scale_x_continuous(limits=rango) +
  #geom_abline(intercept=0, slope=1, linewidth=1,col="gray50",linetype="dashed") +
  geom_hline(yintercept=0, linewidth=1,col="gray50",linetype="dashed") +
  labs(x="Estimated abundance",y="Orthomosaic - Ground turtle counts", title="b)") +
  theme_classic(base_size=16) 


round(sum(rowSums(pop.counts[,c("walking","nesting")], na.rm=T)))
round(sum(rowMeans(pop.counts[,c("ground1","ground2","ground3")], na.rm=T)))
#round(sum(rowMeans(matrix(counts$c.ground,ncol=3), na.rm=T)))

ggsave(file.path("ms","FigS1 Comparison Ground Counts.png"),grid.arrange(fs1a,fs1b,ncol=2),
       dpi=400,units="cm", width=28, height=14)
