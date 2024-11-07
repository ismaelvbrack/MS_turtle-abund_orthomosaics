
library(ggplot2)
library(gridExtra)
library(bayesplot)


# Step1: Mark-resight estimates -------------------------------------------
# Samples from the 1st step
samps1 <- read.table(file.path("outputs","MR5_postMCMC.txt"),h=T)

apply(samps1, 2, function(x)
  c(mean(x),quantile(x, probs=c(0.025,0.975)))
)

color_scheme_set("brightblue")

fig5 <- mcmc_areas(samps1,pars=c("omega","delta","phi.2.","phi.1.","theta"),
                   point_est="mean",prob=0.95) +
  scale_x_continuous(breaks=seq(0,1,0.2),limits=c(0,1)) +
  scale_y_discrete(labels=rev(c(expression("Nesting "(theta)),
                                expression("Availability[1] "(phi[1])),
                                expression("Availability[2] "(phi[2])),
                                expression("Mark identification "(delta)),
                                expression("Double count "(omega))
                                ))) +
  theme_bw(base_size=16) +
  theme(axis.text.y=element_text(face="bold",color="black")) +
  labs(x="Probability") 

# Raw counts --------------------------------------------------------------
pop.counts <- read.table(file.path("data","PopulationCounts.txt"),h=T)
pop.counts$data <- factor(pop.counts$data, levels=pop.counts$data)
J <- nrow(pop.counts)

figa <-
  ggplot(data=data.frame(counts=c(pop.counts$walking,pop.counts$nesting),
                         type=rep(c("walking","nesting"),each=nrow(pop.counts)),
                         date=factor(rep(pop.counts$data, 2), levels=pop.counts$data)
  ), aes(x=date,y=counts, fill=type)) +
  geom_bar(stat="identity",position="dodge",width=0.7) +
  scale_y_continuous(breaks=seq(0,2000,by=500)) +
  scale_fill_manual(values=c("darkorange2","royalblue3")) +
  labs(x="",y="Count data",fill="",title="a)") +
  theme_classic(base_size=16) +
  theme(legend.position="bottom",
        legend.text=element_text(face="bold",size=14)) +
  #theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank())


# Import and Arrange ------------------------------------------------------

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


# Total population size ---------------------------------------------------

figd <- mcmc_areas(samps2, pars=c("Ntot"),point_est="mean",prob=0.95,bw=500) +
  theme_bw(base_size=16) + 
  scale_x_continuous(breaks=seq(30000,60000,4000)) +
  coord_cartesian(ylim=c(1,1.001)) +
  labs(x="",title="Total Population Size") +
  theme(title=element_text(size=14),axis.title.y=element_blank())


# Entrants ---------------------------------------------------------------

par <- "B"

tab1 <- resu[which(rownames(resu) %in% paste0(par,".",1:J,".")),]

tab1$date <- pop.counts$data

#fige <- 
  ggplot(data=tab1,aes(x=date,y=Mean,ymin=`2.5%`,ymax=`97.5%`)) +
  geom_point(size=3) +
  geom_errorbar(linewidth=1,width=.4) +
  labs(x="",y="Number of entrants (B)",title="c)") +
  scale_y_continuous(breaks=seq(0,20000,by=4000),limits=c(0,20000)) +
  theme_classic(base_size=16) +
  #theme(axis.text.x=element_blank())
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))

# Pop. size per day ---------------------------------------------------------------

par <- "N"

tab2 <- resu[which(rownames(resu) %in% paste0(par,".",1:J,".")),]

tab2$date <- pop.counts$data

figc <- 
  ggplot(data=tab2,aes(x=as.factor(date),y=Mean,ymin=`2.5%`,ymax=`97.5%`)) +
  geom_point(size=3) +
  geom_errorbar(linewidth=1,width=.4) +
  labs(x="",y="Estimated Daily Population (Nt)",title="c)") +
  scale_y_continuous(breaks=seq(0,20000,by=4000)) +
  theme_classic(base_size=16) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text())


# Nesting & Walking -------------------------------------------------------

tab <- resu[c(grep("N.nest", rownames(resu)),
              grep("N.walk", rownames(resu))),]

tab$date <- rep(pop.counts$data,2)

tab$state <- rep(c("nesting","walking"),each=J)

tab3 <- subset(tab, state=="walking")

  ggplot(data=tab,aes(x=date,y=Mean,
                      ymin=`2.5%`,ymax=`97.5%`,
                      col=state)) +
  geom_point(size=3, position=position_dodge(0.3)) +
  geom_errorbar(linewidth=1,position=position_dodge(0.3),width=.4) +
  scale_y_continuous(breaks=seq(2000,12000,by=2000)) +
  labs(x="",y="Number of individuals") +
  theme_classic(base_size=16) +
  theme(legend.position="top") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))

# Entering and leaving
tab4 <- rbind(tab1,tab3[-J,-6])
tab4$type <- factor(rep(c("Entries","Departures"),c(J,J-1)),levels=c("Entries","Departures"))
tab4[which(tab4$type=="Departures"),"date"] <- pop.counts$data[-1]

figb <- 
  ggplot(data=tab4,aes(x=date,y=Mean,
                       ymin=`2.5%`,ymax=`97.5%`,
                       col=type)) +
  geom_point(size=2.6, position=position_dodge(0.5)) +
  geom_errorbar(linewidth=0.8,position=position_dodge(0.5),width=0.6) +
  scale_y_continuous(breaks=seq(2000,12000,by=2000)) +
  scale_color_manual(values=c("cyan4","red3")) +
  labs(x="",y="Estimated individuals",title="b)") +
  geom_vline(xintercept=seq(1.5,nrow(pop.counts)),col="gray70",linetype="dashed") +
  theme_classic(base_size=16) +
  theme(legend.position="bottom",legend.title=element_blank(),
        legend.text=element_text(face="bold",size=14)) +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        #axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text())


# -------------------------------------------------------------------------
titl1 <- ggplot() + 
  labs(title = "Mark-resight model") +
  theme(panel.background = element_blank(),
        title=element_text(size=16,face="bold"))
titl2 <- ggplot() + 
  labs(title = "Population counts model") +
  theme(panel.background = element_blank(),
        title=element_text(size=16,face="bold"))


fig1 <- grid.arrange(titl1,
                     fig5,
                     ncol=1,
                     heights=c(0.3,2)
)
ggsave("Model5_Step1-estimates.png",fig1,
       dpi=400,height=12,width=20,units="cm")



fig2 <- grid.arrange(titl2,
             figa,
             figb,
             figc,
             ncol=1,
             heights=c(0.3,2,2,2.1)
             )
ggsave("Model5_Step2-estimates.png",fig2,
               dpi=400,height=32,width=15,units="cm")

