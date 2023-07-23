#-----R packages-----
library(tidyverse)
library(lme4)
library(car)
library(lmerTest)

#----data prepare-----
matr1<-read.csv("data/TestI_raw.csv",header=T)
colnames(matr1)
table(matr1$Generation)
matr1<-matr1[which(matr1$rm==FALSE),]
matr1$P.layer<-factor(matr1$P.layer)
matr1$layer<-factor(matr1$layer)


#----Anc effsize-----

aa<-c()
matr.sub <- subset(matr1,treat.1!="JA")
for (gen in c("Anc")){
  matr.sub.gen<-matr.sub[matr.sub$Generation==gen,]
  print(gen)
  for (trait in c("flt","height","diam","biomass","fruit")){
    matr.sub.gen.tr<-matr.sub.gen
    formu<-as.formula(paste(trait, "~ecotype+ecotype:treat.1+(1|P.layer)+(1|P.row)"))
    mod<-lmer(formu
              ,na.action=na.omit
              ,data=matr.sub.gen.tr
    )
    aaa <- as.data.frame(coef(summary(mod)))
    aaa <- subset(aaa,grepl("treat", rownames(aaa)))
    a<-cbind(gen,trait,aaa)
    print(dim(a))
    aa<-rbind(aa,a)
  }
}

for (treat in c("JA")){
  matr.sub<-matr1[matr1$treat.1=="JActr"|matr1$treat.1==treat,]
  matr.sub$treat.1 <- relevel(factor(matr.sub$treat.1), ref="JActr")
  for (gen in c("Anc")){
    matr.sub.gen<-matr.sub[matr.sub$Generation==gen,]
    print(gen)
    for (trait in c("flt","height","diam","biomass","fruit")){
      matr.sub.gen.tr<-matr.sub.gen
      formu<-as.formula(paste(trait, "~ecotype+ecotype:treat.1+(1|P.layer)+(1|P.row)"))
      mod<-lmer(formu
                ,na.action=na.omit
                ,data=matr.sub.gen.tr
      )
      aaa <- as.data.frame(coef(summary(mod)))
      aaa <- subset(aaa,grepl("treat", rownames(aaa)))
      a<-cbind(gen,trait,aaa)
      print(dim(a))
      aa<-rbind(aa,a)
    }
  }
}

for (treat in c("actr1")){
  matr.sub<-matr1[matr1$treat=="actr2"|matr1$treat==treat,]
  for (gen in c("Anc")){
    matr.sub.gen<-matr.sub[matr.sub$Generation==gen,]
    print(gen)
    for (trait in c("flt","height","diam","biomass","fruit")){
      matr.sub.gen.tr<-matr.sub.gen
      formu<-as.formula(paste(trait, "~ecotype+ecotype:treat+(1|P.layer)+(1|P.row)"))
      mod<-lmer(formu
                ,na.action=na.omit
                ,data=matr.sub.gen.tr
      )
      aaa <- as.data.frame(coef(summary(mod)))
      aaa <- subset(aaa,grepl("treat", rownames(aaa)))
      a<-cbind(gen,trait,aaa)
      print(dim(a))
      aa<-rbind(aa,a)
    }
  }
}

write.csv(aa,"Anc_eff.csv")

#----Test I effsize-----

aa<-c()
matr.sub <- subset(matr1,treat.1!="JA")
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
    aaa <- as.data.frame(coef(summary(mod)))
    aaa <- subset(aaa,grepl("treat", rownames(aaa)))
    a<-cbind(gen,trait,aaa)
    print(dim(a))
    aa<-rbind(aa,a)
  }
}

for (treat in c("JA")){
  matr.sub<-matr1[matr1$treat.1=="JActr"|matr1$treat.1==treat,]
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
      aaa <- as.data.frame(coef(summary(mod)))
      aaa <- subset(aaa,grepl("treat", rownames(aaa)))
      a<-cbind(gen,trait,aaa)
      print(dim(a))
      aa<-rbind(aa,a)
    }
  }
}

for (treat in c("actr1")){
  matr.sub<-matr1[matr1$treat=="actr2"|matr1$treat==treat,]
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
      aaa <- as.data.frame(coef(summary(mod)))
      aaa <- subset(aaa,grepl("treat", rownames(aaa)))
      a<-cbind(gen,trait,aaa)
      print(dim(a))
      aa<-rbind(aa,a)
    }
  }
}

write.csv(aa,"TestI_effect.size.csv")
