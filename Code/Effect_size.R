library(lme4)
library(lmerTest)

matr<-read.csv("data/TestI_raw.csv",header=T)
matr<-matr[which(matr$rm==FALSE),]
matr$P.layer<-factor(matr$P.layer)
matr$layer<-factor(matr$layer)

##-------ef.size-----------------
aa<-matrix(ncol=8)
for (treat in c("Drought","Salt1","Salt2","Cd1","Cd2","Nutrient1","Nutrient2",
                "Leaf")){
  matr.sub<-matr[matr$treat.1=="actr"|matr$treat.1==treat,]
  for (gen in c("F1","F2","F3","F4")){
    matr.sub.gen<-matr.sub[matr.sub$Generation==gen,]
    print(gen)
    for (trait in c("flt","height","diam","biomass","fruit")){
      matr.sub.gen.tr<-matr.sub.gen
      formu<-as.formula(paste(trait, "~layer+ecotype+ecotype:treat.1+(1|P.layer)+(1|P.row)"))
      mod<-lmer(formu
                ,na.action=na.omit
                ,data=matr.sub.gen.tr
      )
      a<-cbind(treat,gen,trait,coef(summary(mod)))
      print(dim(a))
      aa<-rbind(aa,a)
    }
  }
}

for (treat in c("JA")){
  matr.sub<-matr[matr$treat.1=="JActr"|matr$treat.1==treat,]
  matr.sub$treat.1 <- relevel(factor(matr.sub$treat.1), ref="JActr")
  for (gen in c("F1","F2","F3","F4")){
    matr.sub.gen<-matr.sub[matr.sub$Generation==gen,]
    print(gen)
    for (trait in c("flt","height","diam","biomass","fruit")){
      matr.sub.gen.tr<-matr.sub.gen
      formu<-as.formula(paste(trait, "~layer+ecotype+ecotype:treat.1+(1|P.layer)+(1|P.row)"))
      mod<-lmer(formu
                ,na.action=na.omit
                ,data=matr.sub.gen.tr
      )
      a<-cbind(treat,gen,trait,coef(summary(mod)))
      print(dim(a))
      aa<-rbind(aa,a)
    }
  }
}

for (treat in c("actr1")){
  matr.sub<-matr[matr$treat=="actr2"|matr$treat==treat,]
  for (gen in c("F1","F2","F3","F4")){
    matr.sub.gen<-matr.sub[matr.sub$Generation==gen,]
    print(gen)
    for (trait in c("flt","height","diam","biomass","fruit")){
      matr.sub.gen.tr<-matr.sub.gen
      formu<-as.formula(paste(trait, "~layer+ecotype+ecotype:treat+(1|P.layer)+(1|P.row)"))
      mod<-lmer(formu
                ,na.action=na.omit
                ,data=matr.sub.gen.tr
      )
      a<-cbind(treat,gen,trait,coef(summary(mod)))
      print(dim(a))
      aa<-rbind(aa,a)
    }
  }
}

colnames(aa)<-c("treat","gen","trait","Estimate","Std.Error","df","tvalue","Pvalue")
write.table(aa[-1,],"clipboard.txt",sep="\t",col.names=T)
