blt<-function(triplet,sptree,gene_trees){
  result<-data.frame(treelength=character(),branchlength1=character(),branchlength2=character(),
                     interlength=character(),outgrouplength=character(),
                     type=character(),triplet=character(),outgroup=character())

  sptree_trip<-keep.tip(sptree,triplet)
  sp.outtip<-sptree_trip$tip.label[min(sptree_trip$edge[sptree_trip$edge[,1]==length(sptree_trip$tip.label)+1,2])]
  sp.sistertips<-sptree_trip$tip.label[sptree_trip$tip.label!=sp.outtip]

  for(i in 1:length(gene_trees)){
    tre<-gene_trees[[i]]

    tl<-sum(tre$edge.length)

    triplet_tre<-keep.tip(tre,triplet)

    temp<-triplet_tre$edge[,1]==max(triplet_tre$edge)
    bls<-triplet_tre$edge.length[temp]
    bl1<- bls[1]
    bl2<- bls[2]

    ol<-triplet_tre$edge.length[triplet_tre$edge[,2]==min(triplet_tre$edge[triplet_tre$edge[,1]==length(triplet_tre$tip.label)+1,2])]

    il=triplet_tre$edge.length[triplet_tre$edge[,2]==max(triplet_tre$edge[triplet_tre$edge[,1]==length(triplet_tre$tip.label)+1,2])]

    out<-triplet_tre$tip.label[min(triplet_tre$edge[triplet_tre$edge[,1]==length(triplet_tre$tip.label)+1,2])]

    type<-ifelse(out==sp.outtip,"concordant",ifelse(out==sp.sistertips[1],"discordant1","discordant2"))

    r<-data.frame(treelength=tl,branchlength1=bl1,branchlength2=bl2,interlength=il,outgrouplength=ol,
                  type=type
                  ,triplet=paste(triplet,collapse = "_")
                  ,outgroup=out)

    result<-rbind(result,r)

  }
  return(result)

}


wilcox_test<-function(result,triplet,sptree){
  counts<-as.data.frame(table(result$outgroup))
  sptree_trip<-keep.tip(sptree,triplet)
  sptre_outgroup<-sptree_trip$tip.label[min(sptree_trip$edge[sptree_trip$edge[,1]==length(sptree_trip$tip.label)+1,2])]

  concor<-sptre_outgroup
  concor_count<-counts[counts$Var1%in%concor,]

  discor<-counts$Var1[!counts$Var1%in%concor]
  discor_count<-counts[!counts$Var1%in%concor,]

  result$sbl<-(result$branchlength1+result$branchlength2)/result$treelength
  concor_sbl<-result$sbl[result$type=="concordant"]
  discor1_sbl<-result$sbl[result$type=="discordant1"]
  discor2_sbl<-result$sbl[result$type=="discordant2"]

  wilcox_test1<-wilcox.test(concor_sbl,discor1_sbl)$p.value
  wilcox_test2<-wilcox.test(discor1_sbl,discor2_sbl)$p.value
  wilcox_test3<-wilcox.test(concor_sbl,discor2_sbl)$p.value
  # c(wilcox_test1,wilcox_test2,wilcox_test3)
  #chi=chisq.test(discor_count$Freq)$p.value

  re<-data.frame(triplet=paste(triplet,collapse = "_"),
                 Concor_sisBL=mean(concor_sbl),Discor1_sisBL=mean(discor1_sbl),Discor2_sisBL=mean(discor2_sbl),
                 wilcox_ConcorDiscor1=wilcox_test1,wilcox_ConcorDiscor2=wilcox_test2,wilcox_Discor1Discor2=wilcox_test3)

  return(re)
}

