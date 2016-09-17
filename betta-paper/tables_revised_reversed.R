mywd <- "/Users/adw96/Documents/Research w Bunge/4th paper - beta diversity/JRSSC/revisions1/"
setwd(mywd)

## Install nec software, load function
## Need to have CatchAll (and mono) available for the medium-diversity analyses
install.packages("breakaway")
require(breakaway)

## Dethlefsen data; available at
## http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0060280#s5
otu_deth <- read.csv("dethlefsen-data.csv")
rownames(otu_deth) <- otu_deth[,1]
otu_deth <- otu_deth[,2:19] ## delete classification info
freq_deth <- as.data.frame(table(otu_deth[,1]))
freq_deth <- freq_deth[freq_deth[,1]!=0,]

## Whitman data; available at the NCBI SRA under accession number SRP072444
## http://www.ncbi.nlm.nih.gov/sra
otu_whitman <- read.csv("whitman-data.csv")
rownames(otu_whitman) <- otu_whitman[,1]
otu_whitman <- otu_whitman[,-1]
otu_whitman <- as.matrix(otu_whitman)
freq_whitman <- as.data.frame(table(otu_whitman[,1]))
freq_whitman <- freq_whitman[freq_whitman[,1]!=0,]
metadata_whitman <- read.csv("whitman-metadata.csv")
attach(metadata_whitman)

#### Just do 5000 for this one
t1 <- proc.time()
my_iterations  <- 5000

#################################################################################
#################################### Table 1 ####################################
#################################################################################
real_subs_whitman <- function(freq, its, nlower, nupper) {
  outputs <- matrix(NA,nrow=its,ncol=3)
  organisms <- as.numeric(as.character(rep(freq[,1],freq[,2])))
  names(organisms) <- paste("Org",1:length(organisms),sep="")
  orgs_pr <- organisms/sum(organisms)
  ns <- round(runif(its, nlower, nupper))
  for (k in 1:its) {
    my_sub <- sample(names(organisms), ns[k], prob=orgs_pr, replace=T)
    tmp <- cbind(as.numeric(names(table(table(my_sub)))),as.numeric(table(table(my_sub))))
    run <- breakaway(tmp,answers=1,plot=0,print=0)
    baway_results <- c(run$est,run$seest)
    n_results <- sum(tmp[,2])
    try(outputs[k, 1:2] <- baway_results, silent=T)
    try(outputs[k, 3] <- n_results, silent=T)
    if (k %% 5000 == 0) cat("k is", k, "\n")
  }
  colnames(outputs) <- c("breakaway est", "breakaway se", "obs")
  return(outputs)
}

display_function_whitman <- function(simulated, outly=5) {
  tab3 <- matrix(NA,nrow=2, ncol=6); rownames(tab3) <- c("betta", "n")
  colnames(tab3) <- c("dps10, 0.01", "dps10, 0.05", "dps10, 0.10", "dps20, 0.01", "dps20, 0.05", "dps20, 0.10")
  for (dps in c(10,20)) {
    tests <- dim(simulated)[1]/dps
    p_baway <- rep(NA,tests); p_none <- rep(NA, tests)
    X <- cbind(rep(1,dps), seq(ifelse(dps==10, 10, 5), 100, length=dps))
    for (i in 1:tests) {
      range <- (i*dps-(dps-1)):(i*dps)
      mysubset_baway <- simulated[range,1:2]
      mysubset_baway[mysubset_baway[,2]/mysubset_baway[,1]*100 < outly, ] <- NA
      mysubset_obs <- simulated[range,3]
      p_baway[i] <- try(betta(mysubset_baway[,1],mysubset_baway[,2],X)$table[2,3],silent=T)
      p_none[i] <- summary(lm(mysubset_obs ~ X[,2]))$coef[2,4]
    }
    p_baway <- as.numeric(p_baway); p_none <- as.numeric(p_none)
    for (alpha_ind in 1:3) {
      alpha <- c(0.01, 0.05, 0.10)[alpha_ind]
      tab3[,3*dps/10+alpha_ind-3] <- round(c(mean(p_baway<alpha,na.rm=T), mean(p_none<alpha,na.rm=T)),3)
    }
  }
  apply(tab3, 1, function(x) cat(c(" & (", x[1], ",", x[2], ",", x[3], ") & (", x[4], ",", x[5], ",", x[6], ") \\", "\\", "\n"), sep=""))
}

ns <- apply(otu_whitman, 2, sum); range(ns) ## range is ~10,000 to 200,000

## iterations about 0.76 s/each
## seed 1 for table 1
set.seed(1)
whitman_table1a <- real_subs_whitman(freq_whitman, my_iterations, 1e4, 2e5) ## 4.5 hours
display_function_whitman(whitman_table1a, outly=5)


#################################################################################
#################################### Table 2 ####################################
#################################################################################
real_subs_deth <- function(freq, its, nlower, nupper) {
  outputs <- matrix(NA,nrow=its,ncol=3)
  organisms <- as.numeric(as.character(rep(freq[,1],freq[,2])))
  names(organisms) <- paste("Org",1:length(organisms),sep="")
  orgs_pr <- organisms/sum(organisms)
  ns <- round(runif(its, nlower, nupper))
  for (k in 1:its) {
    my_sub <- sample(names(organisms), ns[k], prob=orgs_pr,replace=T)
    tmp <- cbind(as.numeric(names(table(table(my_sub)))),as.numeric(table(table(my_sub))))
    write.table(tmp,paste(mywd,"/temp_catchall.csv",sep=""),sep=",",row.names = FALSE,col.names = FALSE)
    system("mono ./CatchAllCmdL.exe temp_catchall.csv temp_catchall 0 >/dev/null")
    catchall_results <- c(as.numeric(data.matrix(read.csv(paste(mywd,"temp_catchall/temp_catchall_BestModelsAnalysis.csv",sep=""))[1,4])),
                          as.numeric(data.matrix(read.csv(paste(mywd,"temp_catchall/temp_catchall_BestModelsAnalysis.csv",sep=""))[1,5])))
    n_results <- sum(tmp[,2])
    try(outputs[k,1:2] <- catchall_results, silent=T)
    try(outputs[k,3] <- n_results, silent=T)
    if (k %% 5000 == 0) cat("k is", k, "\n")
  }
  colnames(outputs) <- c("catchall est", "catchall se", "obs")
  return(outputs)
}

display_function_deth <- function(simulated, outly=1) {
  tab3 <- matrix(NA,nrow=2,ncol=6); rownames(tab3) <- c("betta", "n")
  colnames(tab3) <- c("dps10, 0.01", "dps10, 0.05", "dps10, 0.10", "dps20, 0.01", "dps20, 0.05", "dps20, 0.10")
  x10p <- as.factor(c(rep("A",5), rep("B",5))); X10 <- model.matrix(lm(rnorm(10)~x10p))[1:10,]
  x20p <- as.factor(c(rep("A",10), rep("B",10))); X20 <- model.matrix(lm(rnorm(20)~x20p))[1:20,]
  for (dps in c(10,20)) {
    tests <- dim(simulated)[1]/dps
    p_baway <- rep(NA,tests); p_catchall <- rep(NA,tests); p_chao <- rep(NA, tests); p_chao_reg <- rep(NA, tests); p_none <- rep(NA, tests)
    if (dps==10) {
      X <- X10; xp <- x10p
    } else {
      X <- X20; xp <- x20p
    }
    for (i in 1:tests) {
      range <- (i*dps-(dps-1)):(i*dps)
      mysubset_catchall <- simulated[range,1:2]
      mysubset_catchall[mysubset_catchall[,2]/mysubset_catchall[,1]*100 < outly, ] <- NA
      mysubset_obs <- simulated[range,3]
      p_catchall[i] <- try(betta(mysubset_catchall[,1],mysubset_catchall[,2],X)$table[2,3],silent=T)
      p_none[i] <- summary(lm(mysubset_obs ~ xp))$coef[2,4]
    }
    p_catchall <- as.numeric(p_catchall); p_none <- as.numeric(p_none)
    for (alpha_ind in 1:3) {
      alpha <- c(0.01, 0.05, 0.10)[alpha_ind]
      tab3[,3*dps/10+alpha_ind-3] <- round(c(mean(p_catchall<alpha, na.rm=T), mean(p_none<alpha,na.rm=T)),3)
    }
  }
  apply(tab3, 1, function(x) cat(c(" & (", x[1], ",", x[2], ",", x[3], ") & (", x[4], ",", x[5], ",", x[6], ") \\", "\\", "\n"), sep=""))
  cat("\n \n")
}

## iterations about 6.3s each = 33 hours
system("mkdir temp_catchall")
set.seed(2)
deth_table1a <- real_subs_deth(freq_deth, my_iterations, 1e4, 3e4)
display_function_deth(deth_table1a, outly=1)
system("rm -r temp_catchall")
system("rm -r temp_catchall.csv")



#############
display_function_whitman(whitman_table1a, outly=5)
display_function_deth(deth_table1a, outly=1)

display_function_whitman(deth_table1a, outly=1)
display_function_deth(whitman_table1a, outly=5)

#################################################################################
#################################### Figure 1 ###################################
#################################################################################
freq_list_whitman <- apply(otu_whitman, 2, function(x) as.data.frame(table(x)))
freq_list_whitman <- lapply(freq_list_whitman, function(x) x[x[,1]!=0,])

Day <- as.factor(Day)
Amdmt <- as.factor(Amdmt)

ns_whitman <- unlist(lapply(freq_list_whitman, function(x) sum(x[,2])))

breakaway_whitman <- matrix(NA,nrow=length(freq_list_whitman),ncol=4)
rownames(breakaway_whitman) <- names(freq_list_whitman)
colnames(breakaway_whitman) <- c("baway_est","baway_seest","baway_lcb","baway_ucb")
for (i in 1:length(freq_list_whitman)) {
  baway <- try(breakaway(freq_list_whitman[[i]],plot=FALSE,print=FALSE,answers=TRUE),silent=T)
  if(class(try(baway$est==0,silent=1))!="try-error") { # if it works, store it
    breakaway_whitman[i,1] <- baway$est
    breakaway_whitman[i,2] <- baway$seest
    breakaway_whitman[i,3:4] <- baway$ci
  }
}

sample_d0_a0 <- Day==0 & Amdmt == 0 & breakaway_whitman[,2]/breakaway_whitman[,1] > 0.01
ests_d0_a0 <- breakaway_whitman[sample_d0_a0,]
ns_d0_a0 <- ns_whitman[sample_d0_a0]

range(ests_d0_a0[,1]/ns_d0_a0)

betta(ests_d0_a0[,1], ests_d0_a0[,2])$table
betta(ests_d0_a0[,1], ests_d0_a0[,2])$homogeneity ## homogeneous!

#pdf("isme_fig3.pdf", height=5, width=6)
par(mfrow=c(1,1))
z_t <- qnorm(0.975, 0, 1); ylim_t1 <- 1000; ylim_t2 <- 7000
plot(0,0,type="n",bty="n",main="",
     xlab="Samples",ylab="Estimates of species richness",cex.lab=0.8,
     xlim=c(0,length(ests_d0_a0[,1])+1),ylim=c(ylim_t1,ylim_t2))
for (i in 1:length((ns_d0_a0))) {
  lines(c(rank(ns_d0_a0)[i],rank(ns_d0_a0)[i]),
        c(max(0,ests_d0_a0[i,1]-z_t*ests_d0_a0[i,2], ns_d0_a0[i]),
          min(ylim_t2,ests_d0_a0[i,1]+z_t*ests_d0_a0[i,2])))
  points(rank(ns_d0_a0)[i],ests_d0_a0[i,1],pch=16,cex=0.8)
  points(rank(ns_d0_a0)[i],ns_d0_a0[i],pch=1,cex=0.8)
}
legend(6,2000, c("Observed richness","Estimated total richness","95% interval estimate of total richness"),
       pch=c(1,16, NA), lty=c(NA, NA, 1),cex=0.7,bty="o",box.col = "grey")
#dev.off()

#################################################################################
#################################### Figure 2 ###################################
#################################################################################
freq_list_deth <- apply(otu_deth, 2, function(x) as.data.frame(table(x)))
freq_list_deth <- lapply(freq_list_deth, function(x) x[x[,1]!=0,])

setwd(mywd)
system("mkdir frequencytables")
system("cp runcatchall0.sh frequencytables")
system("cp CatchAllCmdL.exe frequencytables")
for (i in 1:length(freq_list_deth)) {
  tmp <- freq_list_deth[[i]]
  col1 <- as.numeric(tmp[,1])
  col2 <- as.numeric(tmp[,2])
  freq_list_deth[[i]] <- cbind(col1,col2)
  freqtable <- freq_list_deth[[i]]
  write.table(freqtable,paste(mywd,"/frequencytables/",names(freq_list_deth)[i],".csv",sep=""),sep=",",row.names = FALSE,col.names = FALSE)
}
setwd(paste(mywd,"frequencytables/",sep = ""))
system("./runcatchall0.sh")
setwd(mywd)
txt.sources = paste(mywd,"/",list.files(pattern="*BestModelsAnalysis.csv",recursive=T),sep="")
myread <- function(x)  {
  y <- read.csv(x,stringsAsFactors=FALSE)
  return(as.numeric(y[1,4:7]))
}
catchall_deth <- matrix(sapply(txt.sources,myread),byrow=TRUE,ncol=4)
rownames(catchall_deth) <- sapply(strsplit(txt.sources,"/"),"[",11)
colnames(catchall_deth) <- c("catchall_est","catchall_seest","catchall_lcb","catchall_ucb")

participant <- c(rep("A",8),rep("B",5),rep("C",5))
treatment <- c(rep("prior",4),rep("treatment",2),rep("after",2),
               rep(c(rep("prior",2),rep("treatment",1),rep("after",2)),2))

tmt <- as.factor(treatment); tmt <- relevel(tmt, ref="prior")
xm_r <- model.matrix(lm(rnorm(18)~tmt))
betta_random(catchall_deth[,1], catchall_deth[,2], xm_r, participant)$table

ns_d <- unlist(lapply(freq_list_deth, function(x) sum(x[,2])))

mixed_effects_d <- summary(lme4::lmer(ns_d ~ tmt + (1|participant)))$coef
mixed_effects_d
2*(1-pnorm(abs(mixed_effects_d[,3])))

#pdf("isme_fig2.pdf", height=5, width=6)
col_by_patient <- c(rep("blue",8),rep("green",5),rep("black",5))
character_by_patient <- c(rep(0,8),rep(2,5),rep(5,5))
character_by_patient2 <- c(rep(15,8),rep(17,5),rep(18,5))
par(xpd=T)
plot(0,0,type="n",main="",xlim=c(0,19),ylim=c(0,2500),bty="n",xaxt="n",xlab="",
     ylab="Estimates of species richness",cex.lab=0.8)
z <- qnorm(0.975, 0, 1)
for (i in 1:18) {
  points(i,catchall_deth[i,1],col=col_by_patient[i],pch=character_by_patient2[i],cex=0.6)
  lines(c(i,i),c(catchall_deth[i,1]-z*catchall_deth[i,2],catchall_deth[i,1]+z*catchall_deth[i,2]),col=col_by_patient[i])
  points(i,ns_d[i],col=col_by_patient[i],pch=character_by_patient[i], cex=0.5)
}
xaxs1 <- -250; xaxs2 <- xaxs1; xaxs3 <- -200; txth <- -100; txts <- 0.7; pth <- 2500
text(5,pth,"Patient A",col="blue"); text(11,pth,"Patient B",col="green");text(16,pth,"Patient C",col="black")
lines(c(0.75,4.25),c(xaxs1,xaxs1)); lines(c(0.75,0.75),c(xaxs2,xaxs3)); lines(c(4.25,4.25),c(xaxs2,xaxs3))
lines(c(4.75,6.25),c(xaxs1,xaxs1));lines(c(4.75,4.75),c(xaxs2,xaxs3));lines(c(6.25,6.25),c(xaxs2,xaxs3))
lines(c(6.75,8.25),c(xaxs1,xaxs1));lines(c(6.75,6.75),c(xaxs2,xaxs3));lines(c(8.25,8.25),c(xaxs2,xaxs3))
text(2.5,txth,"PRE",cex=txts);text(5.5,txth,"TR",cex=txts);text(7.5,txth,"POST",cex=txts)
lines(c(8.75,10.25),c(xaxs1,xaxs1));lines(c(8.75,8.75),c(xaxs2,xaxs3));lines(c(10.25,10.25),c(xaxs2,xaxs3))
lines(c(10.75,11.25),c(xaxs1,xaxs1));lines(c(10.75,10.75),c(xaxs2,xaxs3));lines(c(11.25,11.25),c(xaxs2,xaxs3))
lines(c(11.75,13.25),c(xaxs1,xaxs1));lines(c(11.75,11.75),c(xaxs2,xaxs3));lines(c(13.25,13.25),c(xaxs2,xaxs3))
text(9.5,txth,"PRE",cex=txts);text(11,txth,"TR",cex=txts);text(12.5,txth,"POST",cex=txts)
lines(c(13.75,15.25),c(xaxs1,xaxs1));lines(c(13.75,13.75),c(xaxs2,xaxs3));lines(c(15.25,15.25),c(xaxs2,xaxs3))
lines(c(15.75,16.25),c(xaxs1,xaxs1));lines(c(15.75,15.75),c(xaxs2,xaxs3));lines(c(16.25,16.25),c(xaxs2,xaxs3))
lines(c(16.75,18.25),c(xaxs1,xaxs1));lines(c(16.75,16.75),c(xaxs2,xaxs3));lines(c(18.25,18.25),c(xaxs2,xaxs3))
text(14.5,txth,"PRE",cex=txts);text(16,txth,"TR",cex=txts);text(17.5,txth,"POST",cex=txts)
legend(8, 500, c("Observed richness","Estimated total richness","95% interval estimate of total richness"),
       pch=c(1,16, NA), lty=c(NA, NA, 1),cex=0.7,col="dark grey", bty="o",box.col = "grey")
#dev.off()