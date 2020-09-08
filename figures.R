
load("mc_dat_final1.RData") 

mylist$bal_type1$type<-revalue(mylist$bal_type1$type, c("lowss_few"="LSFG", "lowss_many"="LSMG", "highss_few"="HSFG", "highss_many"="HSMG"))

mylist$bal_type1$test<-revalue(mylist$bal_type1$test, c("dmrt"="DMRT", "lsd"="LSD", "lsd.bonf"="LSD Bonf", "lsd.sidak"="LSD Sidák","scheffe"="Scheffe", "snk"="SNK", "t.test_bonf"="t-test Bonf", "t.test_sidak"="t-test Sidak", "tukey"="Tukey"))

mylist$unbal_type1$type<-revalue(mylist$unbal_type1$type, c("lowss_few"="LSFG", "lowss_many"="LSMG", "highss_few"="HSFG", "highss_many"="HSMG"))

mylist$unbal_type1$test<-revalue(mylist$unbal_type1$test, c("dmrt"="DMRT", "lsd"="LSD", "lsd.bonf"="LSD Bonf", "lsd.sidak"="LSD Sidák","scheffe"="Scheffe", "snk"="SNK", "t.test_bonf"="t-test Bonf", "t.test_sidak"="t-test Sidak", "tukey"="Tukey"))

g1<-ggplot(mylist$bal_type1, aes(x=type, y=as.numeric(as.character(mylist$bal_type1$prop)))) +
  geom_bar(aes(fill = test), position = "dodge", stat="identity", color="black") +  
  labs(fill = "MC Test") + 
  xlab("Treatment") + 
  ylab("Proportion type I error") +
  theme_classic(base_size = 12) +
  scale_fill_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6"), labels=c("DMRT", "LSD", "LSD Bonferroni", "LSD Sidák","Scheffé's S ", "SNK", substitute(paste(italic('t'), "-test Bonferroni")), substitute(paste(italic('t'), "-test Sidák")), "HSD")) +
  ggtitle("a) Balanced study design")

g2<-ggplot(mylist$unbal_type1, aes(x=type, y=as.numeric(as.character(mylist$unbal_type1$prop)))) +
  geom_bar(aes(fill = test), position = "dodge", stat="identity", color="black") +  
  labs(fill = "MC Test") + 
  xlab("Treatment") + 
  ylab("Proportion type I error") +
  theme_classic(base_size = 12) +
  scale_fill_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6"),labels=c("DMRT", "LSD", "LSD Bonferroni", "LSD Sidák","Scheffé's S", "SNK", substitute(paste(italic('t'), "-test Bonferroni")), substitute(paste(italic('t'), "-test Sidák")), "HSD")) +
  ggtitle("b) Unbalanced study design")

library(gridExtra)

grid.arrange(g1, g2,nrow = 2)

library(ggridges)

ridge_plot_t1<-function(plot.data){
  #can also examine ridge plots of p-values to identify if those proportions really tell the whole story
  ridge_df<-data.frame("pvals"=c(plot.data$Tukey, plot.data$SNK, plot.data$DMRT, 
                                 plot.data$LSD, plot.data$LSD.bonf, plot.data$LSD.sidak,
                                 plot.data$t.test_bonf, plot.data$t.test_sidak,
                                 plot.data$Scheffe), 
                       test=rep(c("Tukey", "SNK", "DMRT", "Fisher's LSD", "Fisher's LSD Bonferroni", "Fisher's LSD Sidak",
                                  "t-test Bonferroni", "t-test Sidák", "Scheffé's S"), each=length(plot.data$Tukey)))
  
  ggplot(ridge_df, aes(y = test, x = as.numeric(as.character(pvals))), fill=test) +
    geom_density_ridges2(aes(alpha=0.5, fill=test)) +
    xlim(0,1) +
    xlab(substitute(paste(italic('p'), "-values"))) + 
    ylab("Test") +
    geom_vline(xintercept=0.05, lwd=1, col="red", lty=2) + 
    theme(legend.position = "none") +
    scale_fill_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    scale_y_discrete(labels=c("DMRT", "LSD", "LSD Bonferroni", "LSD Sidák", "SNK", "Scheffé's S", substitute(paste(italic('t'), "-test Bonferroni")), substitute(paste(italic('t'), "-test Sidák")), "HSD"))
}


r1<-ridge_plot_t1(mylist$blf_nd)

r1<-r1+theme(axis.title.x=element_blank())+ggtitle("a) LSFG")# Low sample size, few groups")

r2<-ridge_plot_t1(mylist$blm_nd)

r2<-r2+theme(axis.title.x=element_blank(),
             axis.title.y=element_blank(),
             axis.text.y=element_blank())+ggtitle("b) LSMG")# Low sample size, many groups")

r3<-ridge_plot_t1(mylist$bhf_nd)

r3<-r3+ggtitle("c) HSFG")# High sample size, few groups")

r4<-ridge_plot_t1(mylist$bhm_nd)

r4<-r4+theme(axis.title.y=element_blank(),
             axis.text.y=element_blank())+ggtitle("d) HSMG")# High sample size, many groups")

pdf(file="ridges_type_1_error.pdf", width=8,height=6)
grid.arrange(
  grobs = list(r1,r2,r3,r4),
  widths = c(2, 1, 2, 1),
  layout_matrix = rbind(c(1, 1, 2),
                        c(3, 3, 4))
)
dev.off()

mylist$bal_type2$type<-revalue(mylist$bal_type2$type, c("lowss_few"="LSFG", "lowss_many"="LSMG", "highss_few"="HSFG", "highss_many"="HSMG"))
mylist$bal_type2$test<-revalue(mylist$bal_type2$test, c("dmrt"="DMRT", "lsd"="LSD", "lsd.bonf"="LSD.Bonf", "lsd.sidak"="LSD.Sidak",
                                                        "scheffe"="Scheffe", "snk"="S-N-K", "t.test_bonf"="T-test.Bonf", "t.test_sidak"="T-test.Sidak", "tukey"="Tukey"))
mylist$unbal_type2$type<-revalue(mylist$unbal_type2$type, c("lowss_few"="LSFG", "lowss_many"="LSMG", "highss_few"="HSFG", "highss_many"="HSMG"))
mylist$unbal_type2$test<-revalue(mylist$unbal_type2$test, c("dmrt"="DMRT", "lsd"="LSD", "lsd.bonf"="LSD.Bonf", "lsd.sidak"="LSD.Sidak",
                                                            "scheffe"="Scheffe", "snk"="S-N-K", "t.test_bonf"="T-test.Bonf", "t.test_sidak"="T-test.Sidak", "tukey"="Tukey"))


g3<-ggplot(mylist$bal_type2, aes(x=type, y=as.numeric(as.character(mylist$bal_type2$prop)))) +
  geom_bar(aes(fill = test), position = "dodge", stat="identity", color="black") +  
  labs(fill = "MC Test") + 
  xlab("Treatment") + 
  ylab("Proportion type II error") +
  theme_classic(base_size = 12) +
  scale_fill_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6"), labels=c("DMRT", "LSD", "LSD Bonferroni", "LSD Sidák","Scheffé's S", "SNK", substitute(paste(italic('t'), "-test Bonferroni")), substitute(paste(italic('t'), "-test Sidák")), "HSD")) +
  ggtitle("a) Balanced study design")

g4<-ggplot(mylist$unbal_type2, aes(x=type, y=as.numeric(as.character(mylist$unbal_type2$prop)))) +
  geom_bar(aes(fill = test), position = "dodge", stat="identity", color="black") +  
  labs(fill = "MC Test") + 
  xlab("Treatment") + 
  ylab("Proportion type II error") +
  theme_classic(base_size = 12) +
  scale_fill_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6"), labels=c("DMRT", "LSD", "LSD Bonferroni", "LSD Sidák","Scheffé's S", "SNK", substitute(paste(italic('t'), "-test Bonferroni")), substitute(paste(italic('t'), "-test Sidák")), "HSD")) +
  ggtitle("b) Unbalanced study design")

grid.arrange(g3, g4,nrow = 2)

ridge_plot_t2<-function(plot_data){
  #this type of ridge plot only examines the amount of type II error when groups with differences were accounted for
  ridge_df<-data.frame("pvals"=c(plot_data$type_two[,,1], plot_data$type_two[,,2],
                                 plot_data$type_two[,,3], plot_data$type_two[,,4],
                                 plot_data$type_two[,,5], plot_data$type_two[,,6],
                                 plot_data$type_two[,,7], plot_data$type_two[,,8],
                                 plot_data$type_two[,,9]), 
                       test=rep(c("Tukey", "S-N-K", "DMRT", "LSD", "LSD.Bonf", "LSD.Sidak",
                                  "T-test.Bonf", "T-test.Sidak", "Scheffe"), each=300))
  
  ggplot(ridge_df, aes(y = test, x = as.numeric(as.character(pvals))), fill=test) +
    geom_density_ridges2(aes(fill=test, alpha=0.5)) +
    xlim(0,1) +
    xlab(substitute(paste(italic('p'), "-values"))) + 
    ylab("Test") +
    geom_vline(xintercept=0.05, lwd=1, col="red", lty=2) + 
    theme(legend.position = "none") +
    scale_fill_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    scale_y_discrete(labels=c("DMRT", "LSD", "LSD Bonferroni", "LSD Sidák", "SNK", "Scheffé's S", substitute(paste(italic('t'), "-test Bonferroni")), substitute(paste(italic('t'), "-test Sidák")), "HSD"))
}


r1<-ridge_plot_t2(mylist$blf_d)

r1<-r1+theme(axis.title.x=element_blank())+ggtitle("a) LSFG")# Low sample size, few groups")

r2<-ridge_plot_t2(mylist$blm_d)

r2<-r2+theme(axis.title.x=element_blank(),
             axis.title.y=element_blank(),
             axis.text.y=element_blank())+ggtitle("b) LSMG")# Low sample size, many groups")

r3<-ridge_plot_t2(mylist$bhf_d)

r3<-r3+ggtitle("c) HSFG")# High sample size, few groups")

r4<-ridge_plot_t2(mylist$bhm_d)

r4<-r4+theme(axis.title.y=element_blank(),
             axis.text.y=element_blank())+ggtitle("d) HSMG")# High sample size, many groups")

pdf(file="ridges_type_2_error.pdf", width=8,height=6)
grid.arrange(
  grobs = list(r1,r2,r3,r4),
  widths = c(2, 1, 2, 1),
  layout_matrix = rbind(c(1, 1, 2),
                        c(3, 3, 4))
)
dev.off()
