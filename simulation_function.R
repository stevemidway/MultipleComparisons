

#install.packages("agricolae")
library(agricolae)
#install.packages("DescTools")
library(DescTools)
#install.packages("MHTdiscrete")
library(MHTdiscrete)
library(ggplot2)
####Simulate data

####SIMULATE DATA FUNCTION
#n is the number of samples in a group;
#mean is the mean of the distribution;
#sd is the sd of the distribution;
#sims defines the number of simulations to run
#ntests defines the number of MC tests being run
#sig_dif is the location of the means that are significantly different than the others
#so if mean=c(0,1,0) and you want to see how often the tests identify sig diffs from the mean of 1, then you would
#use sig_dif=c(2) because the mean of 1 is the second number in the vector

sim_data<-function(n=c(50,50,50,50,50), mean=c(0,0,0,0,0),
                   sd=c(1,1,1,1,1), sims=50, ntests=9, sig_dif=NA){
  #based on the number of groups, identify the number of combinations
  num_combs<-factorial(length(n))/(2*factorial(length(n)-2))
  
  #set up the array to hold the simulation results
  sig_diffs<-array(dim=c(sims, num_combs, ntests))
  type_two<-array(dim=c(sims, num_combs, ntests))
  #run the simulation loop
  for(j in 1:sims){
    #create a list for the group simulated data
    group_sims<-list()
    for(i in 1:length(n)){
      group_sims[[i]]<-rnorm(n[i], mean[i], sd[i])
    }
    
    
    #Set up data frame to hold data
    group_df<-data.frame("n_x"=unlist(group_sims),
                         "group"=rep(seq(from=1, to=length(n),
                                         by=1), times=n))
    
    ###ANOVA
    a1 <- aov(group_df$n_x ~ as.factor(group_df$group))
    
    ##TUKEY
    tukey <- TukeyHSD(x=a1, "as.factor(group_df$group)", conf.level=0.95)
    sig_diffs[j,1:num_combs,1]<-tukey$`as.factor(group_df$group)`[,4]
    
    new_list<-list()
    for(k in 1:length(sig_dif)){
      #extract the p-values that compare to significantly different groups
      new_list[[k]]<-tukey$`as.factor(group_df$group)`[grep(sig_dif[k],
                                                            rownames(tukey$`as.factor(group_df$group)`)),4]
      
      #subset the data to ensure we don't resample p-values if more than one group that is significantly different
      tukey$`as.factor(group_df$group)`<-tukey$`as.factor(group_df$group)`[-grep(sig_dif[k],
                                                                                 rownames(tukey$`as.factor(group_df$group)`)),]
    }
    all_comps<-unlist(new_list, use.names=FALSE)
    type_two[j,1:length(all_comps),1]<- all_comps
    
    ##Student-Neuman-Keuls Test
    #need to look more into what the group=FALSE is actually doing
    snk <- SNK.test(a1, "as.factor(group_df$group)", alpha = 0.05, group=FALSE)
    sig_diffs[j,1:num_combs,2]<-snk$comparison[,2]
    
    
    new_list<-list()
    for(k in 1:length(sig_dif)){
      #extract the p-values that compare to significantly different groups
      new_list[[k]]<-snk$comparison[grep(sig_dif[k],rownames(snk$comparison)),2]
      #subset the data to ensure we don't resample p-values if more than one group that is significantly different
      snk$comparison<-snk$comparison[-grep(sig_dif[k],rownames(snk$comparison)),]
    }
    all_comps<-unlist(new_list, use.names=FALSE)
    type_two[j,1:length(all_comps),2]<- all_comps
    
    ###Duncan's MRT
    dmrt<-duncan.test(a1, "as.factor(group_df$group)", group=FALSE)
    sig_diffs[j,1:num_combs,3]<-dmrt$comparison[,2]
    
    new_list<-list()
    for(k in 1:length(sig_dif)){
      #extract the p-values that compare to significantly different groups
      new_list[[k]]<-dmrt$comparison[grep(sig_dif[k],rownames(dmrt$comparison)),2]
      #subset the data to ensure we don't resample p-values if more than one group that is significantly different
      dmrt$comparison<-dmrt$comparison[-grep(sig_dif[k],rownames(dmrt$comparison)),]
    }
    all_comps<-unlist(new_list, use.names=FALSE)
    type_two[j,1:length(all_comps),3]<- all_comps
    
    #FISHER LSD test, this can have p.adjust
    lsd<-LSD.test(a1, "as.factor(group_df$group)", group=FALSE)
    sig_diffs[j,1:num_combs,4]<-lsd$comparison[,2]
    
    new_list<-list()
    for(k in 1:length(sig_dif)){
      new_list[[k]]<-lsd$comparison[grep(sig_dif[k],rownames(lsd$comparison)),2]
      
      lsd$comparison<-lsd$comparison[-grep(sig_dif[k],rownames(lsd$comparison)),]
    }
    all_comps<-unlist(new_list, use.names=FALSE)
    type_two[j,1:length(all_comps),4]<- all_comps
    
    
    #for now i will run with bonf correction as well since bonf is so common
    lsd.bonf<-LSD.test(a1, "as.factor(group_df$group)", group=FALSE, p.adj="bonferroni")
    sig_diffs[j,1:num_combs,5]<-lsd.bonf$comparison[,2]
    
    new_list<-list()
    for(k in 1:length(sig_dif)){
      new_list[[k]]<-lsd.bonf$comparison[grep(sig_dif[k],rownames(lsd.bonf$comparison)),2]
    #  
      lsd.bonf$comparison<-lsd.bonf$comparison[-grep(sig_dif[k],rownames(lsd.bonf$comparison)),]
    }
    all_comps<-unlist(new_list, use.names=FALSE)
    type_two[j,1:length(all_comps),5]<- all_comps
    
    #LSD_Sidak
    lsd.sidak<-LSD.test(a1, "as.factor(group_df$group)", group=FALSE, p.adj="none")
    combos<-rownames(lsd.sidak$comparison)
    lsd_sidak_df<-data.frame(combos, pvals=as.vector(lsd.sidak$comparison[,2]))
    no.nas<-na.omit(lsd_sidak_df)
    sidak<-MHTdiscrete::Sidak.p.adjust(no.nas$pvals, 0.05)
    sig_diffs[j,1:num_combs,6]<-sidak
    
    sidak_df<-data.frame(combos=no.nas$combos, pvals=sidak)
    new_list<-list()
    for(k in 1:length(sig_dif)){
      new_list[[k]]<-sidak_df[grep(sig_dif[k],sidak_df$combos),2]
      
      sidak_df<-sidak_df[-grep(sig_dif[k],sidak_df$combos),]
    }
    all_comps<-unlist(new_list, use.names=FALSE)
    type_two[j,1:length(all_comps),6]<- all_comps

    #T.test with bonferonni correction
    bonf_ttest<-pairwise.t.test(group_df$n_x, group_df$group, p.adj = "bonf")
    combos<-paste(rep(colnames(bonf_ttest$p.value),each=length(bonf_ttest$p.value[,1])),
                  rep(rownames(bonf_ttest$p.value), length(bonf_ttest$p.value[,1])), sep="-")
    new_df<-data.frame(combos, pvals=as.vector(bonf_ttest$p.value))
    no.nas.bonf<-na.omit(new_df)
    sig_diffs[j,1:num_combs,7]<-no.nas.bonf$pvals
    
    new_list<-list()
    for(k in 1:length(sig_dif)){
      new_list[[k]]<-no.nas.bonf[grep(sig_dif[k],no.nas.bonf$combos),2]
      
      no.nas.bonf<-no.nas.bonf[-grep(sig_dif[k],no.nas.bonf$combos),]
    }
    all_comps<-unlist(new_list, use.names=FALSE)
    type_two[j,1:length(all_comps),7]<- all_comps
    
    #T.test with Sidak correction
    ttest<-pairwise.t.test(group_df$n_x, group_df$group, p.adj = "none")
    combos<-paste(rep(colnames(ttest$p.value),each=length(ttest$p.value[,1])),
                  rep(rownames(ttest$p.value), length(ttest$p.value[,1])), sep="-")
    ttest_df<-data.frame(combos, pvals=as.vector(ttest$p.value))
    no.nas<-na.omit(ttest_df)
    sidak<-MHTdiscrete::Sidak.p.adjust(no.nas$pvals, 0.05)
    sig_diffs[j,1:num_combs,8]<-sidak
    
    sidak_df<-data.frame(combos=no.nas$combos, pvals=sidak)
    new_list<-list()
    for(k in 1:length(sig_dif)){
      new_list[[k]]<-sidak_df[grep(sig_dif[k],sidak_df$combos),2]
      
      sidak_df<-sidak_df[-grep(sig_dif[k],sidak_df$combos),]
    }
    all_comps<-unlist(new_list, use.names=FALSE)
    type_two[j,1:length(all_comps),8]<- all_comps
    
    ##Scheffe Test
    scheffe <- scheffe.test(a1, "as.factor(group_df$group)", alpha = 0.05, group=FALSE)
    sig_diffs[j,1:num_combs,9]<-scheffe$comparison[,2]
    
    new_list<-list()
    for(k in 1:length(sig_dif)){
      #extract the p-values that compare to significantly different groups
      new_list[[k]]<-scheffe$comparison[grep(sig_dif[k],rownames(scheffe$comparison)),2]
      #subset the data to ensure we don't resample p-values if more than one group that is significantly different
      scheffe$comparison<-scheffe$comparison[-grep(sig_dif[k],rownames(scheffe$comparison)),]
    }
    all_comps<-unlist(new_list, use.names=FALSE)
    type_two[j,1:length(all_comps),9]<- all_comps
    
  }
  

  test.names<-c("tukey", "snk", "dmrt", "lsd","lsd.bonf","lsd.sidak", "t.test_bonf", "t.test_sidak", #"REGW",
                "scheffe")
  groups=length(n)
  prop_type1<-matrix(ncol=2, nrow=9)
  prop_type1[,2]<-test.names
  for(i in 1:9){
    prop_type1[i,1]<-round(length(which(sig_diffs[,,i]<0.05))/(sims*groups),3)
  }
  
  true_diff_grps<-length(n)-length(sig_dif) 
  prop_type2<-matrix(ncol=2, nrow=9)
  prop_type2[,2]<-test.names
  for(i in 1:9){
    bad_vec<-as.vector(type_two[,,i][,1:true_diff_grps])
    prop_type2[i,1]<-round(length(which(bad_vec>0.05))/(sims*true_diff_grps),3)
  }
  
  output_list<-list("Tukey"=as.vector(sig_diffs[,,1]), "SNK"=as.vector(sig_diffs[,,2]),
                    "DMRT"=as.vector(sig_diffs[,,3]),#"nk"=as.vector(sig_diffs[,,4]),
                    "LSD"=as.vector(sig_diffs[,,4]),"LSD.bonf"=as.vector(sig_diffs[,,5]),
                    "LSD.sidak"=as.vector(sig_diffs[,,6]),
                    "t.test_bonf"=as.vector(sig_diffs[,,7]), "t.test_sidak"=as.vector(sig_diffs[,,8]),
                    #"REGW"=as.vector(sig_diffs[,,7]),
                    "Scheffe"=as.vector(sig_diffs[,,9]),
                    "sig_diffs"=sig_diffs,
                    "type_two"=type_two, "prop_type1"=prop_type1,
                    "prop_type2"=prop_type2)
  return(output_list)}


###BALANCED DATA
###No differences

##LOW SAMPLE SIZE

##FEW GROUPS

blf_nd<-sim_data(sims=100, n=c(10,10,10), mean=c(0,0,0), sd=c(3,3,3))

###MANY GROUPS

blm_nd<-sim_data(sims=100, n=rep(10,7), mean=rep(0,7), sd=rep(3,7))


###HIGH SAMPLE SIZE

##FEW GROUPS

bhf_nd<-sim_data(sims=100, n=c(100,100,100), mean=c(0,0,0), sd=c(3,3,3))

##MANY GROUPS

bhm_nd<-sim_data(sims=100, n=rep(100,7), mean=rep(0,7), sd=rep(3,7))

bal_type1<-as.data.frame(rbind(blf_nd$prop_type1, blm_nd$prop_type1, bhf_nd$prop_type1,
                               bhm_nd$prop_type1))
colnames(bal_type1)<-c("prop", "test")
bal_type1$type<-rep(c("lowss_few", "lowss_many", "highss_few", "highss_many"), each=9)

####Differences in mean

blf_d<-sim_data(sims=100, n=c(10,10,10), mean=c(1,0,0), sd=c(3,3,3), sig_dif=1)

blm_d<-sim_data(sims=100, n=c(10,10,10,10,10,10,10),
                mean=c(1,0,0,0,0,0,0), sd=c(3,3,3,3,3,3,3), sig_dif=1)

bhf_d<-sim_data(sims=100, n=c(100,100,100), mean=c(1,0,0), sd=c(3,3,3), sig_dif=1)

bhm_d<-sim_data(sims=100, n=c(100,100,100,100,100,100,100), 
                mean=c(1,0,0,0,0,0,0), sd=c(3,3,3,3,3,3,3), sig_dif=1)

bal_type2<-as.data.frame(rbind(blf_d$prop_type2, blm_d$prop_type2, bhf_d$prop_type2,
                               bhm_d$prop_type2))
colnames(bal_type2)<-c("prop", "test")
bal_type2$type<-rep(c("lowss_few", "lowss_many", "highss_few", "highss_many"), each=9)


###UNBALANCED DATA
###No differences

##LOW SAMPLE SIZE

##FEW GROUPS

ublf_nd<-sim_data(sims=100, n=c(5,10,15), mean=c(0,0,0), sd=c(3,3,3))

###MANY GROUPS

ublm_nd<-sim_data(sims=100, n=c(5,10,10,15,15,5,5), mean=rep(0,7), sd=rep(3,7))


###HIGH SAMPLE SIZE

##FEW GROUPS

ubhf_nd<-sim_data(sims=100, n=c(85,100,115), mean=c(0,0,0), sd=c(3,3,3))

##MANY GROUPS

ubhm_nd<-sim_data(sims=100, n=c(85,85,100,100,115,115,85), mean=rep(0,7), sd=rep(3,7))

unbal_type1<-as.data.frame(rbind(ublf_nd$prop_type1, ublm_nd$prop_type1, ubhf_nd$prop_type1,
                                 ubhm_nd$prop_type1))
colnames(unbal_type1)<-c("prop", "test")
unbal_type1$type<-rep(c("lowss_few", "lowss_many", "highss_few", "highss_many"), each=8)

###Differences

##LOW SAMPLE SIZE

##FEW GROUPS

ublf_d<-sim_data(sims=100, n=c(5,10,15), mean=c(1,0,0), sd=c(3,3,3), sig_dif=1)

###MANY GROUPS

ublm_d<-sim_data(sims=100, n=c(5,10,10,15,15,5,5), 
                 mean=c(1,0,0,0,0,0,0), sd=rep(3,7), sig_dif=1)


###HIGH SAMPLE SIZE

##FEW GROUPS

ubhf_d<-sim_data(sims=100, n=c(85,100,115), mean=c(1,0,0), sd=c(3,3,3), sig_dif=1)

##MANY GROUPS

ubhm_d<-sim_data(sims=100, n=c(85,85,100,100,115,115,85), 
                 mean=c(1,0,0,0,0,0,0), sd=rep(3,7), sig_dif=1)

unbal_type2<-as.data.frame(rbind(ublf_d$prop_type2, ublm_d$prop_type2, ubhf_d$prop_type2,
                                 ubhm_d$prop_type2))
colnames(unbal_type2)<-c("prop", "test")
unbal_type2$type<-rep(c("lowss_few", "lowss_many", "highss_few", "highss_many"), each=8)


mylist<-list(
  bal_type1 = bal_type1,
  bal_type2 = bal_type2, 
  unbal_type1 = unbal_type1,
  unbal_type2 = unbal_type2,
  blf_nd = blf_nd,
  blm_nd = blm_nd,
  bhf_nd = bhf_nd, 
  bhm_nd = bhm_nd,
  blf_d = blf_d,
  blm_d = blm_d,
  bhf_d = bhf_d, 
  bhm_d = bhm_d,
  ublf_nd = ublf_nd,
  ublm_nd = ublm_nd,
  ubhf_nd = ubhf_nd, 
  ubhm_nd = ubhm_nd,
  ublf_d = ublf_d,
  ublm_d = ublm_d,
  ubhf_d = ubhf_d, 
  ubhm_d = ubhm_d
)

save(mylist, file="mc_dat_final1.RData")
