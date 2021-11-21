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


wilcox.text<-function(result,triplet){
  counts<-as.data.frame(table(result$outgroup))
  concor<-max(triplet)
  concor_count<-counts[counts$Var1%in%concor,]
  
  discor<-counts$Var1[!counts$Var1%in%concor]
  discor_count<-counts[!counts$Var1%in%concor,]
  
  result$proxy_t<-(result$branchlength1+result$branchlength2)/result$treelength
  
  concor_proxy_t<-result$proxy_t[result$type=="concordant"]
  discor1_proxy_t<-result$proxy_t[result$type=="discordant1"]
  discor2_proxy_t<-result$proxy_t[result$type=="discordant2"]
  
  wilcox_text1<-wilcox.test(concor_proxy_t,discor1_proxy_t)$p.value
  wilcox_text2<-wilcox.test(discor1_proxy_t,discor2_proxy_t)$p.value
  wilcox_text3<-wilcox.test(concor_proxy_t,discor2_proxy_t)$p.value
  chi=chisq.test(discor_count$Freq)$p.value
  
  re<-data.frame(triplet=paste(triplet,collapse = "_"),outgroup=counts$Var1,frequence=counts$Freq,chisq=chi,
                 concor_proxy_t=mean(concor_proxy_t),discor1_proxy_t=mean(discor1_proxy_t),discor2_proxy_t=mean(discor2_proxy_t),
                 wilcox_text_cd1=wilcox_text1,wilcox_text_cd2=wilcox_text2,wilcox_text_d1d2=wilcox_text3)
  
  
  return(re)
}
