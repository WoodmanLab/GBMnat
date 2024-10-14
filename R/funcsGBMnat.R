# This script contains functions used in GBMnat.Rmd.

# ---- Kaplan-Meier survival analysis ----
#' get KM medium survival time
#'
#' @param survdf data frame with survival data.
#' @param duration name of the column with the duration of the event.
#' @param variable name of the column with the variable to stratify the analysis.
#' @param censoring name of the column with the censoring information.
#' @param conversion "d2m" for days to months conversion.
#' @param timeUnit time unit for the duration column.
#'
#' @return a data frame with the medium survival time for each group.
#' @export getKM_medium
#'
#' @import dplyr survival

getKM_medium=function(survdf,duration,variable=NULL,censoring, conversion="d2m",timeUnit="month"){
  require(dplyr)
  require(survival)
  survdf=as.data.frame(survdf)

  if(conversion=="d2m"){
    survdf[[duration]]=survdf[[duration]]*12/365.25
    timeUnit="month"
  }

  if(!is.null(variable)){
    survdf[[variable]]=as.factor(survdf[[variable]])
    colnames(survdf)[[which(colnames(survdf)==variable)]]="theVariable"
  }
  colnames(survdf)[[which(colnames(survdf)==duration)]]="theDuration"
  colnames(survdf)[[which(colnames(survdf)==censoring)]]="theStatus"

  if(!is.null(variable)){
    dff=summary(survfit(Surv(theDuration, theStatus) ~ theVariable, data = survdf))$table
    names(dff)=gsub("theVariable=","",names(dff))
    # dff=data.frame(variable=names(dff),medium=dff)
    # colnames(dff)=c(variable,"medium")
  }else{
    dff=summary(survfit(Surv(theDuration, theStatus) ~ 1, data = survdf))$table
  }

  return(dff)
}

#' generate Kaplan-Meier survival curve
#'
#' @param survdf \code{data.frame()}.
#' Should at least include columns of duration (eg. survival), censoring (eg. vital status), grouping varible.
#' @param duration  \code{character()}. The name of the column of survival.
#' @param censoring \code{character()}. The name of the column of vital status.
#' @param variable \code{character()}. The name of the column of grouping variable.
#' @param conversion if "d2m". Convert duration unit from days to months.
#' @param palette \code{character()}. color palette. character vector.
#' @param OrderedVarFactor \code{character()}. order of the factors in grouping variable.
#' @param plotTitle \code{character()}. plot title.
#' @param table_base_size annotation font size.
#'
#' @return a Kaplan-Meier plot.
#' @importFrom survminer surv_fit
#'
#' @export
#'
plotKM=function(survdf,duration,censoring,variable=NULL,conversion="d2m",statOnly=F,palette=NULL,OrderedVarFactor=NULL,plotTitle=NULL,table_base_size=6,timeUnit="month",timeYLab="OS"){
  require(survminer)
  require("survival")
  require(gridExtra)
  survdf=as.data.frame(survdf)
  if(!is.null(variable)){
    survdf[[variable]]=as.factor(survdf[[variable]])
    if(is.null(palette)){
      palette=cols4all::c4a("dark24",n = nlevels(survdf[[variable]]))
    }

    if(!is.null(OrderedVarFactor)){
      survdf[[variable]]=factor(survdf[[variable]],levels = OrderedVarFactor)
      # names(palette)=paste(variable,"=",OrderedVarFactor,sep="")
      names(palette)=levels(survdf[[variable]])
    }
  }

  if(conversion=="d2m"){
    survdf[[duration]]=survdf[[duration]]*12/365.25
    timeUnit="month"
  }


  if(is.null(variable)){
    call=as.formula(sprintf("Surv(%s, %s) ~ 1",duration,censoring))
  }else{
    call=as.formula(sprintf("Surv(%s, %s) ~ %s",duration,censoring,variable))
  }
  fit=survminer::surv_fit(call,data = survdf)

  if(!is.null(variable)){
    logRankTest.oa=survdiff(call,data = survdf)$pvalue

    if(nlevels(survdf[[variable]])==2){
      tmp=summary(coxph(call, data = survdf))
      stats=setNames(
        c(logRankTest.oa,tmp$conf.int[,"exp(coef)"],tmp$conf.int[,"lower .95"],tmp$conf.int[,"upper .95"]),
        c("p.value","HR","lowerCI95","upperCI95"))
      stats.str=stats.str_=sprintf("p = %.3f \nHR = %.3f (%.3f - %.3f)",
                                   logRankTest.oa,tmp$conf.int[,"exp(coef)"],tmp$conf.int[,"lower .95"],tmp$conf.int[,"upper .95"])
    }
    if(nlevels(survdf[[variable]])>2){
      logRankTest.pw=pairwise_survdiff(call,data = survdf,p.adjust.method="none")$p.value
      HR=list()
      refU=unique(combn(levels(survdf[[variable]]),2)[1,])
      for (ref in refU){
        dataSurv=survdf[,c(duration,censoring,variable)]
        dataSurv[[variable]]=factor(dataSurv[[variable]],levels = c(ref,levels(dataSurv[[variable]])[levels(dataSurv[[variable]])!=ref]))
        tmp=summary(coxph(call, data =dataSurv))
        tmp=cbind(
          HR=tmp$coefficients[,"exp(coef)"] %>% format(digits=2),
          CI=paste(tmp$conf.int[,"lower .95"]%>% format(digits=2),tmp$conf.int[,"upper .95"]%>% format(digits=2),sep = "-"),
          `p.coxph`=tmp$coefficients[,"Pr(>|z|)"] %>% format(digits=2))
        rownames(tmp)=apply(combn(levels(dataSurv[[variable]]),2)[,1:(nlevels(dataSurv[[variable]])-1),drop=F], 2, function(x) paste(x[2],"vs", x[1],sep = ""))
        HR[[ref]]=tmp
      }
      if(nlevels(dataSurv[[variable]])>2){
        HR=Reduce(rbind,HR)
        HR=HR[!duplicated(lapply(strsplit(rownames(HR),"vs"), sort)),]
      }else{HR=HR[[1]]}
      stats=cbind(HR,`p.logrank`=na.omit(as.vector(logRankTest.pw)) %>% format(digits=2))
      stats=rbind(stats,overall=c(NA,NA,NA,logRankTest.oa%>% format(digits=2)))
      # stats.str=readr::format_delim(stats %>% as.data.frame%>% tibble::rownames_to_column("pair"),delim = "\t")
      stats.str=stats.str_=capture.output(stats %>% as.data.frame%>%knitr::kable("simple",align = 'c'))
      stats.str=paste(stats.str[stats.str!=""&!grepl("-----",stats.str)],collapse = "\n")
    }
    if(statOnly){
      # if(nlevels(survdf[[variable]])==2){
      #   return(logRankTest.oa)
      # }else if(nlevels(dataSurv[[variable]])>2){
      #   return(stats)
      return(stats)

    }else{
      ggsurv=ggsurvplot(fit, data = survdf,
                        palette = palette,
                        # pval = TRUE,
                        xlab = paste("Time",sprintf("(%s)",timeUnit)),
                        ylab = sprintf("%s (probability)",timeYLab),
                        conf.int = F,
                        risk.table = TRUE,
                        risk.table.col = "strata",
                        legend.title = variable,
                        risk.table.y.text = F,
                        legend.labs=levels(survdf[[variable]]),
                        risk.table.height = 0.25,
                        break.x.by = 10,
                        ggtheme = theme_bw())
      if(nlevels(survdf[[variable]]) %in% c(1,2)){
        ggsurv$plot=ggsurv$plot +
          annotate("text",
                   size=(table_base_size+0.4)/.pt,
                   x=range(survdf[[duration]],na.rm = T)[2]*1/4,y=0.2,
                   label=stats.str,alpha = .9,)

      }else if(nlevels(survdf[[variable]])>2){
        ggsurv$plot=ggsurv$plot +
          annotation_custom(gridExtra::tableGrob(stats,theme=ttheme_minimal(base_size = table_base_size)),
                            xmin=range(na.omit(survdf[[duration]]))[2]*2.5/5,
                            xmax=range(na.omit(survdf[[duration]]))[2]*2.5/5,
                            ymin=0.75,ymax=0.90)
      }
      if(!is.null(plotTitle)){
        ggsurv$plot=ggsurv$plot +
          labs(title = sprintf("Kaplan-Meier Curve on %s",duration))
      }

      ggsurv$plot=ggsurv$plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(legend.position=c(.8,.8))
      ggsurv$table=ggsurv$table + theme(panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank())
    }
  }else{
    ggsurv=ggsurvplot(fit, data = survdf,
                      # palette = palette,
                      # pval = TRUE,
                      xlab = paste("Time",sprintf("(%s)",timeUnit)),
                      ylab = sprintf("%s (probability)",timeYLab),
                      conf.int = F,
                      risk.table = TRUE,
                      legend=c(.8,.8),
                      risk.table.col = "strata",
                      # legend.title = variable,
                      risk.table.y.text = F,
                      # legend.labs=levels(survdf[[variable]]),
                      risk.table.height = 0.25,
                      break.x.by = 10,
                      ggtheme = theme_bw())
  }

  return(ggsurv)

}

# ---- get GRanges object for copy number analysis ----

#' @title Convert data.frame to GRanges
#'
#' @param df  data.frame. Must contain columns of chromosome, start, end, strand and other metadata columns.
#' @param genome genome version. Default is "hg19".
#' @param seqlevelsStyle style/format of seqlevels. Default is "NCBI".
#' @param simplified logical. Weather to keep only major chromosomes or include other scaffolds. Default is TRUE.
#' @param xy logical. Weather to include X and Y chromosomes. Default is FALSE.
#' @param seqnames_col column name of chromosome.
#' @param start_col column name of start position.
#' @param end_col column name of end position.
#' @param strand_col column name of strand.
#' @param meta_cols column names of metadata.
#'
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import GenomicRanges
#' @return GRanges object.
#' @export
#'
df2granges<-function(df,genome=c("hg19","hg38"),seqlevelsStyle=c("NCBI","UCSC"),simplified=TRUE,xy=FALSE,seqnames_col="chromosome",start_col="start",end_col="end",strand_col=NULL,meta_cols=NULL){
  # library(GenomicRanges)
  genome=match.arg(genome)
  if(genome=='hg19'){
    library(BSgenome.Hsapiens.UCSC.hg19)
    txdb<-BSgenome.Hsapiens.UCSC.hg19
  } else {
    library(BSgenome.Hsapiens.UCSC.hg38)
    txdb<-BSgenome.Hsapiens.UCSC.hg38
  }
  seqinfo_<-seqinfo(txdb)
  seqs=seqlevels(txdb)
  if(simplified){
    seqs=seqlevels(txdb)[1:24]
    if(!xy){
      seqs=seqlevels(txdb)[1:22]
    }
    seqinfo_<-seqinfo_[seqs]
  }
  df[[seqnames_col]]<-paste("chr",df[[seqnames_col]],sep="")
  df<-df[df[[seqnames_col]] %in% seqs,]
  if(is.null(strand_col)){
    strand="*"
  } else {
    strand=df[[strand_col]]
  }
  if(is.null(meta_cols)) {
    mcols=NULL
  } else {
    mcols=df[,meta_cols,drop=F]
  }
  granges<-GRanges(seqnames=df[[seqnames_col]],ranges=IRanges(start=df[[start_col]],end=df[[end_col]]),strand = strand,mcols=mcols,seqinfo=seqinfo_)
  colnames(mcols(granges))=meta_cols
  return(granges)
}

#' @title Bin GRanges object
#'
#' @param gr GRanges object.
#' @param window window size for binning. Default is 500000.
#'
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import GenomicRanges
#' @import IRanges
#' @return GRanges object.
#' @export

bingranges<-function(gr,window=500000){
  seqinfo_<-seqinfo(gr)
  bingenome<-lapply(1:length(seqinfo_),function(i){chromosome=seqinfo_@seqnames[i];
  start=seq(1,seqinfo_@seqlengths[i],by=window);
  end=start+window-1;
  end[length(end)]=seqinfo_@seqlengths[i];
  return(GRanges(seqnames=chromosome,ranges = IRanges(start=start,end=end),strand="*",seqinfo = seqinfo_))})
  bingenome<-do.call(c,bingenome)
  default_dict<-list("character"="","logical"=NA,"integer"=0L,"numeric"=0)
  for (col in colnames(mcols(gr))){
    mcols(bingenome)[[col]]<-default_dict[[class(mcols(gr)[[col]])]]
  }
  ind<-findOverlaps(bingenome,gr,select="first")
  mcols(bingenome)[!is.na(ind),]<-mcols(gr)[ind[!is.na(ind)],]
  return(bingenome)
}

#---- prepare unsupervised data ----

#' prepare_unsupervised_data
#'
#' @description extract qualified subset of data with most variance for unsupervised analysis
#' @param expressions \code{data.frame()}. Gene expressions table, optimally clean and filtered by prepare_clean_RNA_sample_info_and_protein_expressions.
#' @param method \code{character()}. Methods for calculation of variance, and etc, and later used by filtering. CV for coefficient of variance, DQ for quantile of both mean and variance (dual quantile), GUMBLE for gumbel distribution of mad. Default set as "MAD".
#' @param mad_top_n \code{numeric()}. Number of genes to be kept with top variance, if it is -1 keep all genes. Default set as -1.
#' @param mad_top_quantile \code{numeric()}. Percentage of genes to be kept with top mad values if mad_top_n is -1. Default set as 0.75.
#' @param cv_top_n \code{numeric()}. Number of genes to be kept with top cv values, if it is -1 keep all genes. Default set as -1,
#' @param cv_top_mean_quantile \code{numeric()}. Percentage of genes to be kept with top mad values if cv_top_n is -1. Default set as 0.5.
#' @param dq_top_mean_quantile \code{numeric()}. Genes to be kept with mean value quantile above defined dq_top_mean_quantile. Default set as 0.5.
#' @param dq_top_var_quantile \code{numeric()}. Genes to be kept with variance value quantile above defined dq_top_var_quantile. Default set as 0.5.
#' @param gumbel_p_cutoff \code{numeric()}. Genes with p value of MAD gumbel distribution less than gumbel_p_cutoff will be kept. Default set as 0.1.
#' @param remove_outlier \code{logical(1)}. Whether to remove outlier value before calculation of variance, and etc. Default set as FALSE.
#' @import goeveg ordinal
#' @return a \code{data.frame()} of a subset of expressions.
#' @export
#' @examples
#' \dontrun{
#' results<-prepare_unsupervised_data(log2_expression,method="MAD",mad_top_n=1000,remove_outlier=F)
#' }

prepare_unsupervised_data<-function(expressions,method=c("MAD","CV","DQ","GUMBEL"),mad_top_n=-1,mad_top_quantile=0.75,cv_top_n=-1,cv_top_mean_quantile=0.5,dq_top_mean_quantile=0.5,dq_top_var_quantile=0.5,gumbel_p_cutoff=0.1,remove_outlier=F){
  library(goeveg)
  library(ordinal)
  removeoutlier<-function(data){
    quartiles<-quantile(data,probs=c(0.25,0.75),na.rm=T)
    IQR<-IQR(data)
    Lower <- quartiles[1] - 1.5*IQR
    Upper <- quartiles[2] + 1.5*IQR
    data_wo_outlier <- subset(data, data > Lower & data < Upper)
    return(data_wo_outlier)
  }

  method=method[1] #XL1
  #MAD
  if(method=='MAD'){
    mads<-apply(expressions,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);mad=mad(dat,na.rm=T);return(mad)})
    if(mad_top_n==-1){
      results=expressions[mads>quantile(mads,mad_top_quantile),]
    } else{
      results=expressions[order(mads,decreasing = T)[1:mad_top_n],]
    }
  }
  #GUMBEL
  if(method=='GUMBEL'){
    mads<-apply(expressions,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);mad=mad(dat,na.rm=T);return(mad)})
    gumbel_mean=mean(mads,na.rm=T)
    gumbel_var<-var(mads,na.rm=T)
    threshold=qgumbel(1-gumbel_p_cutoff,gumbel_mean,gumbel_var)
    results=expressions[mads>=threshold,]
  }
  #CV
  if(method=="CV"){
    means=apply(expressions,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);mean=mean(dat,na.rm=T);return(mean)})
    #      rowMeans(expressions,na.rm=T)
    top_mean_filter=quantile(means,cv_top_mean_quantile,na.rm=T)
    expressions_=expressions[means>=top_mean_filter,]
    expressions_<-expressions_[rowSums(is.na(as.matrix(expressions_)))<(0.25*ncol(expressions_)),]
    cvs=apply(expressions_,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);cv=cv(dat,na.rm=T);return(cv)})
    if(cv_top_n==-1){
      results=expressions_
    } else {
      results=expressions_[order(cvs,decreasing = T)[1:min(cv_top_n,nrow(expressions_))],]
    }
  }
  #DQ
  if(method=="DQ"){
    means=apply(expressions,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);mean=mean(dat,na.rm=T);return(mean)})
    top_mean_filter=quantile(means,dq_top_mean_quantile,na.rm=T)
    expressions_=expressions[(means>=top_mean_filter & (!is.na(means))),]
    expressions_<-expressions_[rowSums(is.na(as.matrix(expressions_)))<(0.25*ncol(expressions_)),]
    vars=apply(expressions_,1,function(dat) {if(remove_outlier) dat=removeoutlier(dat);var=stats::var(dat,na.rm=T);return(var)})
    top_var_filter=quantile(vars,dq_top_var_quantile)
    results<-expressions_[(!is.na(vars)) & vars>=top_var_filter,]
  }
  results<-results[rowSums(is.na(as.matrix(results)))<(0.25*ncol(results)),]
  return(results)
}

# ---- clustering functions ----
#' estimate best number of clusters
#'
#' @description estimate the best number of cluster for clustering, using 30 indices from NbClust::NbClust
#'
#' @param data \code{data.frame()} or \code{matrix()}. Gene expression data (or matrix).
#' @param scale \code{character()}. Direction for scale data. Default ="row".
#' @param data_preTreat \code{logical()}. If TRUE, data will be PCA transformed; Components whose cumulative sum of variance reach variance_explained_cutoff will be kept. Default FALSE.
#' @param removeVar \code{numeric()}. Remove this percentage of variables based on low variance. Default=0.2.
#' @param variance_explained_cutoff \code{numeric()}. If data_preTreat is TRUE, components in PCA transformed data matrix whose cumulative sum of variance should reach this value. default=0.8.
#' @param method \code{character()}.The clustering method to be used. This should be one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans". Default = "complete".
#' @param min.nc \code{numeric()}. Minimal number of clusters, between 1 and (number of objects - 1), default=2.
#' @param plotOP \code{logical()}. Whether to output plots.
#' @param max.nc \code{numeric()}. Maximal number of clusters, between 2 and (number of objects - 1), greater or equal to min.nc. By default, max.nc=15.
#'
#' @return \code{list()} containing \code{Best.nc()} (statistics for various indices) and \code{Best.NumberofCluster()} (best number of clusters).
#' @export
#'
#' @examples
#' \dontrun{
#' result<-estimate_bestNumberofClusters(data=data,scale='none',data_preTreat = F,max.nc = 8)
#' }
estimate_bestNumberofClusters<-function(
    data,
    scale=c("row","column","none"),
    data_preTreat=T,
    removeVar=0.2,
    variance_explained_cutoff=0.8,
    method="complete",
    min.nc=2,
    max.nc=15,
    plotOP=FALSE
){
  if(ncol(data)<nrow(data)){cat("sample space is less than feature space, reduction of feature space using PCA would be recommended!\n\n")}
  scale=match.arg(scale)
  if(scale=="none"){
    data<-t(data)
  }
  if(scale=="row"){
    data<-scale(t(data))
  }
  if(scale=="column"){
    data=scale(data)
  }

  if(data_preTreat){
    pca_res<-PCAtools::pca(t(data),removeVar = removeVar)
    n_components<-min(sum(cumsum(pca_res$variance)<variance_explained_cutoff)+1,nrow(data))
    data<-scale(pca_res$rotated[,pca_res$components[1:n_components]])
    # res<-NbClust(data,  min.nc=min.nc, max.nc=max.nc, method = method, index = "all",plotOP = plotOP)
    # # bestNumberofCluster<-as.integer(names(which.max(table(as.numeric(res$Best.nc["Number_clusters",])))))
    # bestNumberofCluster=max(res$Best.partition)
    index_for_NbClust<-c("kl","ch","hartigan",'ccc',"scott","marriot","trcovw","tracew","friedman","rubin",
                         "cindex","db","silhouette","duda","pseudot2","beale","ratkowsky","ball","ptbiserial","gap",
                         "frey","mcclain","gamma","gplus","tau","dunn","hubert","sdindex",'dindex',"sdbw")
    res<-lapply(index_for_NbClust,function(ind){
      res<-tryCatch(
        {NbClust(data,method=method,index=ind,min.nc = min.nc,max.nc = max.nc,plotOP = plotOP)$Best.nc},
        error=function(e){return(c("Number_clusters"=NA,"Value_Index"=NA))})
      # if(is.null(res)){return(c("Number_clusters"=NA,"Value_Index"=NA))};
      return(res)})
    names(res)<-index_for_NbClust
    res=res[!sapply(res,is.null)]
    res<-as.data.frame(res)
    bestNumberofCluster<-as.integer(names(which.max(table(as.numeric(res["Number_clusters",])))))
  }else{
    index_for_NbClust<-c("kl","ch","hartigan",'ccc',"scott","marriot","trcovw","tracew","friedman","rubin",
                         "cindex","db","silhouette","duda","pseudot2","beale","ratkowsky","ball","ptbiserial","gap",
                         "frey","mcclain","gamma","gplus","tau","dunn","hubert","sdindex",'dindex',"sdbw")
    res<-lapply(index_for_NbClust,function(ind){
      res<-tryCatch(
        {NbClust(data,method=method,index=ind,min.nc = min.nc,max.nc = max.nc,plotOP = plotOP)$Best.nc},
        error=function(e){return(c("Number_clusters"=NA,"Value_Index"=NA))})
      # if(is.null(res)){return(c("Number_clusters"=NA,"Value_Index"=NA))};
      return(res)})
    names(res)<-index_for_NbClust
    res=res[!sapply(res,is.null)]
    res<-as.data.frame(res)
    bestNumberofCluster<-as.integer(names(which.max(table(as.numeric(res["Number_clusters",])))))}
  return(list(Best.nc=res,Best.NumberofCluster=bestNumberofCluster))
}

#' map_cluster
#' @description map cluster_id from query clusters to cluster_id from target clusters
#' @param query_clusters integer vector.
#' @param target_clusters integer vector.
#' @param bestnumberofclusters integer vector, containing mapped cluster_id from target clusters
#'
#' @return integer vector, containing mapped cluster_id from target clusters
#'
map_clusters<-function(query_clusters,target_clusters,bestnumberofclusters){
  stopifnot("query clusters should be as same length as target clusters!"=length(query_clusters)==length(target_clusters))
  cluster_table<-table(query_clusters,target_clusters)
  cluster_table<-cluster_table[order(apply(cluster_table,1,max),decreasing=T),]
  assigned_target_cluster_ids<-NULL
  extra_k=1
  for(i in rownames(cluster_table)){
    if(length(assigned_target_cluster_ids)>=length(unique(target_clusters))){
      assigned_target_cluster_ids<-c(assigned_target_cluster_ids,as.character(length(unique(target_clusters))+extra_k))
      names(assigned_target_cluster_ids)[length(assigned_target_cluster_ids)]<-i
      extra_k=extra_k+1
    }
    cluster_table<-cluster_table[,setdiff(colnames(cluster_table),assigned_target_cluster_ids),drop=F]
    if(max(cluster_table[i,])!=0){
      assigned_target_cluster_ids<-c(assigned_target_cluster_ids,colnames(cluster_table)[which.max(cluster_table[i,])])
      names(assigned_target_cluster_ids)[length(assigned_target_cluster_ids)]<-i
    } else {
      j=bestnumberofclusters+1
      assigned_target_cluster_ids<-c(assigned_target_cluster_ids,as.character(j))
      names(assigned_target_cluster_ids)[length(assigned_target_cluster_ids)]<-i
      j=j+1
    }
  }
  return(as.integer(assigned_target_cluster_ids[order(names(assigned_target_cluster_ids))]))
}

# unified cluster call (kmeans,pam,hclust,fuzzy,mclust,apclust,dbscan,mclclust,specc,kkmean,skmeans,nmf,som), not export for user
kmeansCluster<-function(data,k,...){return(stats::kmeans(t(data),centers=k,...)$cluster)}
pamCluster<-function(data,k,...){return(cluster::pam(t(data),k = k,...)$clustering)}
hclustCluster<-function(data,k,...){
  hclust.out <- stats::hclust(dist(t(data)))
  hclust_clusters<-stats::cutree(hclust.out,k=k,...)
  return(hclust_clusters)
}
fuzzyCluster<-function(data,k,...){return(cluster::fanny(t(data),k=k,...)$clustering)}
mclustCluster<-function(data,k,...){return(mclust::Mclust(t(data),G=k)$classification)}
apclustCluster<-function(data,k,...){
  apclust_clusters<-apcluster::apclusterK(s=negDistMat(r=2),x=t(data),K=k)
  apclust_clusters<-rep(1:length(apclust_clusters@clusters),times=sapply(apclust_clusters@clusters,length))[order(unlist(apclust_clusters@clusters))]
  names(apclust_clusters)<-colnames(data)
  return(apclust_clusters)
}
hdbscanCluster<-function(data,k,...){return(dbscan::hdbscan(t(data),minPts = k)$cluster+1)}
mclCluster<-function(data,k,...){
  cor_mat<-cor(data)
  mcl_clusters<-MCL::mcl(cor_mat,addLoops = T,ESM=T,allow1=T)$Cluster
  names(mcl_clusters)<-colnames(data)
  return(mcl_clusters)
}

speccCluster<-function(data,k,...){
  require(kernlab)
  specc_clusters<-kernlab::specc(t(data),centers=k,...)@.Data
  names(specc_clusters)<-colnames(data)
  return(specc_clusters)
}

kkmeansCluster<-function(data,k,...){
  require(kernlab)
  kkmeans_cluster<-kernlab::kkmeans(t(data),centers=k,...)@.Data
  names(kkmeans_cluster)<-colnames(data)
  return(kkmeans_cluster)
}

skmeansCluster<-function(data,k,...){
  skmeans_cluster<-skmeans::skmeans(x=t(data),k=k,...)$cluster
  names(skmeans_cluster)<-colnames(data)
  return(skmeans_cluster)
}

nmfCluster<-function(data,k,...){
  require(NMF)
  data<-(data-min(data,na.rm=T))/(max(data,na.rm=T)-min(data,na.rm=T))
  fit = NMF::nmf(data, rank = k, ...)
  nmf_cluster<-apply(fit@fit@H, 2, which.max)
  names(nmf_cluster)<-colnames(data)
  return(nmf_cluster)
}

somCluster<-function(data,k,...){
  kr = floor(sqrt(ncol(data)))
  somfit = kohonen::som(t(data), grid = somgrid(kr, kr, "hexagonal"), ...)
  m = somfit$codes[[1]]
  m = m[seq_len(nrow(m)) %in% somfit$unit.classif, ]
  cl = cutree(hclust(dist(m)), k)
  group = numeric(ncol(data))
  for(cl_unique in unique(cl)) {
    ind = as.numeric(gsub("V", "", names(cl)[which(cl == cl_unique)]))
    l = somfit$unit.classif %in% ind
    group[l] = cl_unique
  }
  names(group)<-colnames(data)
  return(group)
}

#' multiCluster
#'
#' @description estimate clustering using multiple clustering algorithms, such as "kmeans","pam","hclust","fuzzy","mclust","apclust","hdbscan","MCL","specc","kkmeans","SKmeans","NMF","SOM".
#' @param data \code{data.frame()} or \code{matrix()}. Gene expression data (or matrix).
#' @param scale \code{character()}. Direction for scale data. Default ="row".
#' @param bestnumberofclusters \code{integer(1)}. Optimal number of clusters, if missing, calculated by estimate_bestNumberofClusters.
#' @param data_preTreat \code{logical()}. If TRUE, data will be PCA transformed; Components whose cumulative sum of variance reach variance_explained_cutoff will be kept. Default FALSE.
#' @param removeVar \code{numeric()}. Remove this percentage of variables based on low variance. Default=0.2.
#' @param variance_explained_cutoff \code{numeric()}. If data_preTreat is TRUE, components in PCA transformed data matrix whose cumulative sum of variance should reach this value. default=0.8.
#' @param nbclust_method 	the cluster analysis method to be used. This should be one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans".
#' @param min.nc minimal number of clusters, between 1 and (number of objects - 1)
#' @param max.nc maximal number of clusters, between 2 and (number of objects - 1), greater or equal to min.nc. By default, max.nc=15.
#' @param cluster_methods \code{vector(mode="character")}. Methods used for clustering and evaluation. default=c("kmeans","pam","hclust","fuzzy","mclust","apcluster","hdbscan","MCL","specc","kkmeans","SKmeans","NMF","SOM").
#' @param base_cluster_method \code{character()}. if missing, the most correlated method will be used.
#' @importFrom stats kmeans hclust cutree
#' @importFrom cluster pam fanny
#' @importFrom mclust Mclust mclustBIC adjustedRandIndex
#' @importFrom dbscan hdbscan
#' @importFrom MCL mcl
#' @importFrom skmeans skmeans
#' @import NMF
#' @import kernlab
#' @import apcluster
#' @importFrom kohonen som somgrid
#' @importFrom ggcorrplot ggcorrplot
#'
#' @return list, containing best partition of sample, and co-clusters derived from multiple clustering algorithms.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' multiCluster(data=data,scale='none',data_preTreat = T)
#' }
#'
multiCluster<-function(
    data,
    scale=c("row","column","none"),
    bestnumberofclusters,
    data_preTreat=F,
    removeVar=0.2,
    variance_explained_cutoff=0.8,
    nbclust_method="complete",
    min.nc=2,
    max.nc=15,
    cluster_methods=c("kmeans","pam","hclust","fuzzy","mclust","apclust","hdbscan","MCL","specc","kkmeans","skmeans","nmf","som"), #
    base_cluster_method){

  scale=match.arg(scale)
  Best.partition=NULL
  result=list()
  if(missing(bestnumberofclusters)){
    bnc<-estimate_bestNumberofClusters(data,scale=scale,data_preTreat=data_preTreat,removeVar = removeVar, variance_explained_cutoff = variance_explained_cutoff,method = nbclust_method,min.nc=min.nc,max.nc=max.nc)
    bestnumberofclusters<-bnc$Best.NumberofCluster
    Best.partition=bnc$Best.nc$Best.partition
    cat("The best number of clusters is ",bestnumberofclusters,"\n\n")
    result[["Best.partition"]]=Best.partition
  }
  data_o<-data
  if(scale=="row"){
    data<-t(scale(t(data)))
  }
  if(scale=="column"){
    data=t(scale(data))
  }

  kmeans_clusters=pam_clusters=hclust_clusters=
    fuzzy_clusters=mclust_clusters=apcluster_clusters=
    hdbscan_clusters=mcl_clusters=specc_clusters=
    kkmeans_clusters= skmeans_clusters=nmf_clusters=
    som_clusters=NULL


  if("kmeans" %in% cluster_methods){
    kmeans_clusters<-kmeansCluster(data,k = bestnumberofclusters)
  }
  if("pam" %in% cluster_methods){
    pam_clusters<-pamCluster(data,k = bestnumberofclusters)
  }
  if("hclust" %in% cluster_methods){
    hclust_clusters<-hclustCluster(data,k=bestnumberofclusters)
  }
  if("fuzzy" %in% cluster_methods){
    fuzzy_clusters<-fuzzyCluster(data,k=bestnumberofclusters)
  }
  if("mclust" %in% cluster_methods){
    mclust_clusters<-mclustCluster(data,k=bestnumberofclusters)
  }

  if("apclust" %in% cluster_methods){
    apclust_clusters<-apclustCluster(data,k=bestnumberofclusters)
  }
  if("hdbscan" %in% cluster_methods){
    hdbscan_clusters<-hdbscanCluster(data,k = bestnumberofclusters)
  }
  if("MCL" %in% cluster_methods){
    mcl_clusters<-mclCluster(data,k=bestnumberofclusters)
  }
  if("specc" %in% cluster_methods){
    specc_clusters<-speccCluster(data,k=bestnumberofclusters)
  }
  if("kkmeans" %in% cluster_methods){
    kkmeans_clusters<-kkmeansCluster(data,k=bestnumberofclusters)
  }
  if("skmeans" %in% cluster_methods){
    skmeans_clusters<-skmeansCluster(data,k=bestnumberofclusters)
  }
  if("nmf" %in% cluster_methods){
    nmf_clusters<-nmfCluster(data_o,k=bestnumberofclusters)
  }
  if("som" %in% cluster_methods){

    som_clusters<-somCluster(data,k=bestnumberofclusters)
  }
  res<-rbind(kmeans_clusters,pam_clusters,hclust_clusters,fuzzy_clusters,mclust_clusters,apclust_clusters,hdbscan_clusters,mcl_clusters,specc_clusters,kkmeans_clusters,skmeans_clusters,nmf_clusters,som_clusters)
  rownames(res)<-cluster_methods
  res=res[apply(res,1,function(x) length(unique(x)))>1,] # XL fix error when methods output same cluster number for all samples (cannot set k for hdbscan) https://hdbscan.readthedocs.io/en/latest/faq.html#q-hdbscan-is-failing-to-separate-the-clusters-i-think-it-should

  if(missing(base_cluster_method)){ # XL base_cluster_method=names(which.max(colSums(cor(t(res)))))
    pairs=data.frame(utils::combn(rownames(res),2))
    rand.sim=data.frame(matrix(nrow = nrow(res),ncol = nrow(res), dimnames = list(rownames(res),rownames(res))))
    for (x in pairs){rand.sim[x[1],x[2]]=rand.sim[x[2],x[1]]=mclust::adjustedRandIndex(res[x[1],],res[x[2],])}
    base_cluster_method=names(which.max(rowSums(rand.sim,na.rm = T)))
    p=ggcorrplot::ggcorrplot(rand.sim,hc.order = TRUE,type = "lower",lab = TRUE,method = c("circle"),title = "Similarity matrix (ARI)")
    # The adjusted Rand Index (ARI) should be interpreted as follows: ARI >= 0.90 excellent recovery; 0.80 =< ARI < 0.90 good recovery; 0.65 =< ARI < 0.80 moderate recovery; ARI < 0.65 poor recovery
    result[["rand.plot"]]=p
    result[[" base_method"]]=base_cluster_method
    result[["rand.sim"]]=rand.sim
  }
  res<-as.data.frame(t(as.data.frame(lapply(as.data.frame(t(res)),function(col){map_clusters(col,res[base_cluster_method,],bestnumberofclusters = bestnumberofclusters)[col]}))))
  res_<-as.data.frame(matrix(paste("Cluster_",as.matrix(res),sep=""),nrow=nrow(res)))
  rownames(res_)<-rownames(res)
  colnames(res_)<-colnames(data)
  result[["CoClusters"]]=res_
  return(result)
}


future_consensusCluster<-function(data,scale=c("row","column","none"),method=c("bootstrap","perturb","combine"),subFeatureSize=0.8,subSampleSize=1,noise=1,cutFUN,nTimes=100,clusters=2,verbose=F,num_cores,...){
  #	A function that, given a data matrix, returns a vector of cluster assignments. Examples of functions with this behavior are cutHclust, cutKmeans, cutPam, and cutRepeatedKmeans, or kmeansCluster,pamCluster,hclustCluster,fuzzyCluster,mclustCluster,apclustCluster,hdbscanCluster,mclCluster,speccCluster,kkmeansCluster,skmeansCluster,nmfCluster,somCluster
  # library(future.apply)
  tictoc::tic("future_consensusCluster:")
  tictoc::tic("prep:")
  scale=match.arg(scale)
  if(scale=="row"){
    data<-t(scale(t(data)))
  }
  if(scale=="column"){
    data=t(scale(data))
  }
  method=match.arg(method)
  if(missing(cutFUN)){
    cutFUN=c(kmeansCluster,pamCluster,hclustCluster,fuzzyCluster,mclustCluster,apclustCluster,hdbscanCluster,mclCluster,speccCluster,kkmeansCluster,skmeansCluster,nmfCluster,somCluster)
  }else{
    if(any(cutFUN%in%c("cutHclust", "cutKmeans","cutPam","cutRepeatedKmeans"))){
      require(ClassDiscovery)
      tmp=setNames(
        c(kmeansCluster,pamCluster,hclustCluster,fuzzyCluster,mclustCluster,apclustCluster,hdbscanCluster,mclCluster,speccCluster,kkmeansCluster,skmeansCluster,nmfCluster,somCluster,cutHclust,cutKmeans,cutPam,cutRepeatedKmeans),
        c("kmeans","pam","hclust","fuzzy","mclust","apclust","hdbscan","MCL","specc","kkmeans","skmeans","nmf","som","cutHclust", "cutKmeans","cutPam","cutRepeatedKmeans"))
    }else{
      tmp=setNames(
        c(kmeansCluster,pamCluster,hclustCluster,fuzzyCluster,mclustCluster,apclustCluster,hdbscanCluster,mclCluster,speccCluster,kkmeansCluster,skmeansCluster,nmfCluster,somCluster),
        c("kmeans","pam","hclust","fuzzy","mclust","apclust","hdbscan","MCL","specc","kkmeans","skmeans","nmf","som"))
    }

    cutFUN=tmp[cutFUN]
  }
  N <- ncol(data)
  subFeatureSize <- as.integer(nrow(data)*subFeatureSize)
  subSampleSize <- as.integer(ncol(data)*subSampleSize)
  tictoc::toc()
  tictoc::tic("future_lapply:")
  stableMatch <- matrix(0, nrow = N, ncol = N,)
  rownames(stableMatch)<-colnames(data)
  colnames(stableMatch)<-colnames(data)

  if(missing(num_cores)){
    num_cores <- max(1,future::availableCores()-2)
  }
  future::plan(future::multisession, workers = num_cores)
  loopFun=function(){
    for (i1 in 1:nTimes) {
      if(method=="bootstrap"){
        tempData <- data[sample(1:nrow(data), subFeatureSize, replace = F),sample(1:ncol(data), subSampleSize, replace = F)]
      }
      if(method=="perturb"){
        tempData <- data + matrix(rnorm(N * nrow(data), 0, noise),ncol = N)
      }
      if(method=="combine"){
        tempData <- data[sample(1:nrow(data), subFeatureSize, replace = F),sample(1:ncol(data), subSampleSize, replace = F) ]
        tempData <- tempData + matrix(rnorm(ncol(tempData) * nrow(tempData), 0, noise),ncol = ncol(tempData))
      }

      if (verbose) {
        cat(paste("[", i1, "] ", nrow(tempData), " ", sep = ""))
        if (i1%%10 == 0)
          cat("\n")
      }
      for(k in clusters){
        for(cutfun in cutFUN){
          tempCut <- cutfun(tempData, k=k,...)
          tempMatch <- matrix(0, nrow = N, ncol = N)
          rownames(tempMatch)<-colnames(data)
          colnames(tempMatch)<-colnames(data)

          for (i2 in 1:k) {
            for(samp1 in names(tempCut)){
              for(samp2 in names(tempCut)){
                tempMatch[samp1, samp2] <- ifelse(tempCut[samp1]==tempCut[samp2],1,0)
              }
            }
          }
          stableMatch <- stableMatch + tempMatch
        }
      }
    }
    return(stableMatch)
  }
  stableMatch =future::value(future::future({loopFun()}))
  tictoc::toc()
  # tictoc::tic("fold:")
  # stableMatch <- future.apply:::fold(tempMatch,function(a,b){a+b})
  # tictoc::toc()

  if (verbose)
    cat("\n")
  result=stableMatch/stableMatch[row(stableMatch)==col(stableMatch)]
  tictoc::toc()
  return(result)
}


#' @title consensusCluster
#'
#' @description consensus clustering using bootstrap or(and) noise addition of expression.
#'
#' @param data \code{data.frame()} or \code{matrix()}. Gene expression data.
#' @param scale \code{character()}. Default ="row", Direction for scale data.
#' @param method \code{character()}. default="bootstrap",  method used to produce subset expression for clustering
#' @param subFeatureSize \code{numeric()}. Default=0.8. The subset fraction of features for bootstrap.
#' @param subSampleSize \code{numeric()}. Default=1. The subset fraction of samples for bootstrap.
#' @param noise \code{numeric()}. Default=1, noise with variance of noise defined unit will be added onto expression
#' @param cutFUN \code{function()}. A vector of clustering functions. Options include:
#' \itemize{
#'  \item "Methods in ClassDiscovery package" - "cutHclust", "cutKmeans", "cutPam", and "cutRepeatedKmeans"
#'  \item "Methods provided by expr" - "kmeansCluster", "pamCluster", "hclustCluster", "fuzzyCluster", "mclustCluster", "apclustCluster", "hdbscanCluster", "mclCluster", "speccCluster", "kkmeansCluster", "skmeansCluster", "nmfCluster", "somCluster"
#' }
#' @param nTimes \code{integer()}. Default=100. Times of iteration of Clustering
#' @param clusters \code{vector(mode = "integer")}. Default=2, predefined number of clusters.
#' @param verbose \code{logical()}. default \code{FALSE}, if \code{TRUE}, print detailed process information.
#' @param ... params passed on to cutFUN.
#'
#' @importFrom mclust Mclust mclustBIC adjustedRandIndex
#' @importFrom stats kmeans hclust cutree
#' @importFrom cluster pam fanny
#' @importFrom mclust Mclust mclustBIC adjustedRandIndex
#' @importFrom dbscan hdbscan
#' @importFrom MCL mcl
#' @importFrom skmeans skmeans
#' @import NMF
#' @import kernlab
#' @import apcluster
#' @return \code{data.frame()}, containing possibility of co-clustering.
#' @export
#'
#' @examples
#' \dontrun{
#' results<-future_consensusCluster(data,method='combine',clusters=4)
#' }
consensusCluster<-function(data,scale=c("row","column","none"),method=c("bootstrap","perturb","combine"),subFeatureSize=0.8,subSampleSize=1,noise=1,cutFUN,nTimes=100,clusters=2,verbose=F,...){
  #	A function that, given a data matrix, returns a vector of cluster assignments. Examples of functions with this behavior are cutHclust, cutKmeans, cutPam, and cutRepeatedKmeans, or kmeansCluster,pamCluster,hclustCluster,fuzzyCluster,mclustCluster,apclustCluster,hdbscanCluster,mclCluster,speccCluster,kkmeansCluster,skmeansCluster,nmfCluster,somCluster
  scale=match.arg(scale)
  if(scale=="row"){
    data<-t(scale(t(data)))
  }
  if(scale=="column"){
    data=t(scale(data))
  }
  method=match.arg(method)
  if(missing(cutFUN)){
    cutFUN=c(kmeansCluster,pamCluster,hclustCluster,fuzzyCluster,mclustCluster,apclustCluster,hdbscanCluster,mclCluster,speccCluster,kkmeansCluster,skmeansCluster,nmfCluster,somCluster)
  }else{
    if(any(cutFUN%in%c("cutHclust", "cutKmeans","cutPam","cutRepeatedKmeans"))){
      require(ClassDiscovery)
      tmp=setNames(
        c(kmeansCluster,pamCluster,hclustCluster,fuzzyCluster,mclustCluster,apclustCluster,hdbscanCluster,mclCluster,speccCluster,kkmeansCluster,skmeansCluster,nmfCluster,somCluster,cutHclust,cutKmeans,cutPam,cutRepeatedKmeans),
        c("kmeans","pam","hclust","fuzzy","mclust","apclust","hdbscan","MCL","specc","kkmeans","skmeans","nmf","som","cutHclust", "cutKmeans","cutPam","cutRepeatedKmeans"))
    }else{
      tmp=setNames(
        c(kmeansCluster,pamCluster,hclustCluster,fuzzyCluster,mclustCluster,apclustCluster,hdbscanCluster,mclCluster,speccCluster,kkmeansCluster,skmeansCluster,nmfCluster,somCluster),
        c("kmeans","pam","hclust","fuzzy","mclust","apclust","hdbscan","MCL","specc","kkmeans","skmeans","nmf","som"))
    }
    cutFUN=tmp[cutFUN]
  }
  N <- ncol(data)
  subFeatureSize <- as.integer(nrow(data)*subFeatureSize)
  subSampleSize <- as.integer(ncol(data)*subSampleSize)
  #  stableSampling<-
  stableMatch <- matrix(0, nrow = N, ncol = N,)
  rownames(stableMatch)<-colnames(data)
  colnames(stableMatch)<-colnames(data)
  for (i1 in 1:nTimes) {
    if(method=="bootstrap"){
      tempData <- data[sample(1:nrow(data), subFeatureSize, replace = F),sample(1:ncol(data), subSampleSize, replace = F)]
    }
    if(method=="perturb"){
      tempData <- data + matrix(rnorm(N * nrow(data), 0, noise),ncol = N)
    }
    if(method=="combine"){
      tempData <- data[sample(1:nrow(data), subFeatureSize, replace = F),sample(1:ncol(data), subSampleSize, replace = F) ]
      tempData <- tempData + matrix(rnorm(ncol(tempData) * nrow(tempData), 0, noise),ncol = ncol(tempData))
    }

    if (verbose) {
      cat(paste("[", i1, "] ", nrow(tempData), " ", sep = ""))
      if (i1%%10 == 0)
        cat("\n")
    }
    for(k in clusters){
      for(cutfun in cutFUN){
        tempCut <- tryCatch(
          {cutfun(tempData, k=k,...)},
          error=function(e){return(setNames(rep(NA,ncol(tempData)),colnames(tempData)))})

        tempMatch <- matrix(0, nrow = N, ncol = N)
        rownames(tempMatch)<-colnames(data)
        colnames(tempMatch)<-colnames(data)

        for (i2 in 1:k) {
          for(samp1 in names(tempCut)){
            for(samp2 in names(tempCut)){
              tempMatch[samp1, samp2] <- ifelse(tempCut[samp1]==tempCut[samp2] & !is.na(tempCut[samp1]==tempCut[samp2]),1,0)
            }
          }
        }
        stableMatch <- stableMatch + tempMatch
      }
    }
  }
  if (verbose)
    cat("\n")
  result=stableMatch/stableMatch[row(stableMatch)==col(stableMatch)]
  return(result)
}

# ---- dge ----
#' dge_limma
#'
#' @description differential gene expression analysis using limma algorithm.
#' @param expressions \code{data.frame()}. Gene expressions table.
#' @param is_rawcount \code{logical()}. Whether expressions is raw count or not. Default set to FALSE.
#' @param is_logged \code{logical()}. Whether provided expressions was log2 transformed or not. Default set to TRUE.
#' @param normalize \code{logical()}. Normalize provided expressions by library size. Default set to FALSE.
#' @param sample_frequency_threshold \code{numeric()}. Genes with low expression occurred in at least sample_frequency_threshold fraction of all samples will be removed. Default 0.5.
#' @param clinic_info \code{data.frame()}. Clinical information table.
#' @param ID_col \code{character()}. Column name for Sample ID (or Patient ID), should be consistent with column names of expressions.
#' @param group_col \code{character()}. Column name for group factor that differential gene expression will be conducted on.
#' @param covariate_col \code{character()}. column name for covariate factor that effects should be removed from model.
#' @param block_col \code{character()}. Column name for block factor if test is conducted on block model, such as paired test.
#' @param contrasts \code{vector(mode="character")} Specific contrasts if preferred, elements should be exactly same as group factor.
#' @param method \code{character(1)}. Method to conduct test. One of:
#' \itemize{
#'  \item "limma_trend" for non raw count expressions
#'  \item "limma_voom" for raw count expressions.
#' }
#' Default set as "limma_trend".
#'
#' @import DESeq2 statmod
#' @return \code{list()}, contains expressions, method, design, contrasts, test, and statistics of limma test.
#' @export
#'
#' @examples
#' \dontrun{
#' data(woodman)
#' clean_log2_protein_expressions=log2_expressions
#' results<-dge_limma(
#'   clean_log2_protein_expressions,
#'   clinical_info = clean_RNA_sample_info,
#'   ID_col = "Sample_ID",
#'   group_col = "APOLLO_TIMPOINT",
#'   contrasts = c("TP2-Baseline"),
#'   method ="limma_trend")
#' }
#'
#'
dge_limma<-function(expressions,is_rawcount=FALSE,is_logged=T,normalize=FALSE,sample_frequency_threshold=0.5,
                    clinic_info,
                    ID_col,
                    group_col,
                    covariate_col,
                    block_col,
                    contrasts,
                    method=c("limma_trend","limma_voom")){
  require(limma)
  require(edgeR)
  require(DESeq2)
  stopifnot("ID column was not found in clinic_info"=ID_col %in% colnames(clinic_info))
  stopifnot("group column was not found in clinic_info"=group_col %in% colnames(clinic_info))
  clinic_info[[group_col]]<-factor(clinic_info[[group_col]])
  if(!missing(covariate_col)) {
    stopifnot("covariate column was not found in clinic_info" = covariate_col %in% colnames(clinic_info))
    clinic_info[[covariate_col]]<-factor(clinic_info[[covariate_col]])
    if(any(levels(clinic_info[[covariate_col]]) %in% levels(clinic_info[[group_col]]))){clinic_info[[covariate_col]]<-factor(paste(clinic_info[[covariate_col]],"_",sep=""))}
    formula<-paste(paste("~",covariate_col,sep="+"),group_col,sep="+")
    design<-model.matrix(as.formula(formula),data=clinic_info)
    colnames(design)<-c("Base",paste(levels(clinic_info[[covariate_col]])[-1],levels(clinic_info[[covariate_col]])[1],sep="-"),paste(levels(clinic_info[[group_col]])[-1],levels(clinic_info[[group_col]])[1],sep="-"))
    if(!missing(contrasts)){
      contrast.matrix<-matrix(0,nrow=ncol(design),ncol=length(contrasts))
      colnames(contrast.matrix)<-contrasts
      rownames(contrast.matrix)<-colnames(design)
      for(contrast in contrasts){
        if(contrast %in% rownames(contrast.matrix)){
          contrast.matrix[contrast,contrast]<-1
        } else {
          contrast_<-paste(unlist(strsplit(contrast,"-")),levels(clinic_info[[group_col]])[1],sep="-")
          contrast.matrix[contrast_,contrast]<-c(1,-1)
        }
      }

    }
  } else {
    formula<-paste("~0+",group_col,sep="")
    design<-model.matrix(as.formula(formula),data=clinic_info)
    colnames(design)<-levels(clinic_info[[group_col]])
    contrast.matrix<-makeContrasts(contrasts=contrasts,levels = design)
  }
  if(!missing(block_col)) {
    stopifnot("block column was not found in clinic_info"=block_col %in% colnames(clinic_info))
    clinic_info[[block_col]]<-factor(clinic_info[[block_col]])
  }
  stopifnot("expression colnames do not match ID column of clinic_info"=all(colnames(expressions)==clinic_info[[ID_col]]))
  method=match.arg(method)
  if(is_rawcount){
    dge <- DGEList(counts=expressions)
    keep <- filterByExpr(dge, design)
    dge <- dge[keep,,keep.lib.sizes=FALSE]
    dge <- calcNormFactors(dge)
    if(method=="limma_trend"){
      logCPM <- cpm(dge, log=TRUE, prior.count=3)
      is_rawcount<-FALSE
      is_logged<-TRUE
      expressions<-logCPM
    } else if(method=="limma_voom"){
      v <- voom(dge, design, plot=TRUE, normalize="quantile")
      expressions<-v$E
      if(!missing(block_col)){
        cor <- duplicateCorrelation(v, design, block = clinic_info[[block_col]])
        v <- voom(v, design, plot = TRUE, block = clinic_info[[block_col]], correlation = cor$consensus)
        cor <- duplicateCorrelation(v, design, block = clinic_info[[block_col]])
        fit <- lmFit(v, design, block = clinic_info[[block_col]], correlation = cor$consensus)
      } else {
        fit <- lmFit(v, design)
      }
    }
  }
  if(!is_rawcount) {
    if(!is_logged){
      expressions<-log2(expressions+1)
    }
    if(normalize){
      library_size<-apply(expressions,2,sum)
      expressions<-t(t(expressions)/library_size)*mean(library_size)
    }
    expressions<-expressions[(rowSums(expressions==0))<(sample_frequency_threshold*ncol(expressions)),]
    cat("non raw_count data will be analyzed with limma\n")
    if(!missing(block_col)){
      corfit <- duplicateCorrelation(expressions,design,block=clinic_info[[block_col]])
      fit <- lmFit(expressions,design,block=clinic_info[[block_col]],correlation=corfit$consensus)
    } else {
      fit <- lmFit(expressions,design)
    }
  }
  fit<-contrasts.fit(fit,contrast.matrix)
  fit<-eBayes(fit,robust = T)
  group_mean<-list()
  for (group in levels(clinic_info[[group_col]])){
    group_mean[[paste(group,"_mean",sep="")]]<-rowMeans(expressions[,clinic_info[[ID_col]][clinic_info[[group_col]]==group]],na.rm=T)
  }
  group_mean<-as.data.frame(group_mean)
  contrast_statistics=group_mean
  for(contrast in contrasts){
    statistics=topTable(fit,coef=contrast,number=nrow(expressions))
    colnames(statistics)<-paste(contrast,colnames(statistics),sep=":")
    contrast_statistics<-cbind(contrast_statistics,statistics[rownames(contrast_statistics),])
  }
  if(length(contrasts)>=2){
    F_statistics<-topTable(fit,number = nrow(expressions))
    contrast_statistics<-cbind(group_mean,F_statistics[rownames(group_mean),-c(1:length(contrasts))],contrast_statistics[,-c(1:ncol(group_mean))])
  }
  return(list(expressions=expressions,method=method,design=design,contrast.matrix=contrast.matrix,fit=fit,statistics=contrast_statistics))
}

#---- gsea ----

#'  Gene Set Enrichment Analysis
#'
#' @param statistics statistics from limma or other differential expression analysis
#' @param output_dir output directory
#' @param logFC_col colname of logFC
#' @param pval_col colname of p-value
#' @param pval_cutoff p-value cutoff
#' @param pathways pathway gene sets
#' @param enrichment_padj_cutoff padj cutoff for enrichment
#' @param pathway_gene_table_file file name of pathway gene table
#' @param fgseatableplot_file file name of fgsea table plot
#' @param enrich_pathways up2down2 or up2down
#' @param enrichment_map whether to plot enrichment map
#' @param FC_cutoff fold change cutoff
#' @param kappa_cutoff kappa cutoff
#' @param ... other parameters
#' @import pathfindR
#' @return
#' @export
#'
gsea<-function(statistics,output_dir=".",logFC_col="logFC",pval_col="P.Value",pval_cutoff=0.05,pathways,enrichment_padj_cutoff=0.05,pathway_gene_table_file,fgseatableplot_file="fgseatableplot.tiff",enrich_pathways="up2down2",enrichment_map=T,FC_cutoff=1.5,kappa_cutoff=0.3,...){
  if(!dir.exists(output_dir)){dir.create(output_dir)}
  if(class(statistics)=="data.frame"){
    stats<-statistics[[logFC_col]]
    names(stats)<-rownames(statistics)
  }
  fgseaRes <- fgsea(pathways = pathways, stats = stats,minSize=15,maxSize=500,nperm=10000)
  fgseaRes<-fgseaRes[order(fgseaRes$NES,decreasing=T)]
  sig_pathways<-pathways[fgseaRes$pathway[fgseaRes$padj<=enrichment_padj_cutoff & (!is.na(fgseaRes$padj))]]
  if(length(sig_pathways)==0) {
    cat("no significant signaling pathway was found.")
    return(NULL)
  }
  if(sum(fgseaRes$padj<=enrichment_padj_cutoff  & (!is.na(fgseaRes$padj)))>30){
    sig_pathways<-pathways[fgseaRes$pathway[order(fgseaRes$padj)][1:30]]
  }
  pathway_gene_table=pathway_gene_table(sig_pathways,stats)
  if(missing(pathway_gene_table_file)){
    pathway_gene_table_file="pathway_gene_table.csv"
  }
  write.csv(pathway_gene_table,file.path(output_dir,pathway_gene_table_file))

  tiff(filename = file.path(output_dir,fgseatableplot_file),units="in", width=15, height=5, res=300, compression = 'lzw')
  fgseatableplot(pathways = sig_pathways,stats = stats,fgseaRes = fgseaRes)
  dev.off()

  if(!is.na(enrich_pathways)){
    up_enrich_pathways<-NULL
    down_enrich_pathways<-NULL
    enrich_pathways_<-NULL
    if(grepl("up",enrich_pathways)){
      n_up_pathways<-as.numeric(unlist(regmatches(enrich_pathways,gregexpr("(?<=up)[0-9]+",enrich_pathways,perl=T))))
      up_enrich_pathways<-head(fgseaRes$pathway,n_up_pathways)
      enrich_pathways_<-c(up_enrich_pathways,down_enrich_pathways)
    }
    if(grepl("down",enrich_pathways)){
      n_down_pathways<-as.numeric(unlist(regmatches(enrich_pathways,gregexpr("(?<=down)[0-9]+",enrich_pathways,perl=T))))
      down_enrich_pathways<-tail(fgseaRes$pathway,n_down_pathways)
      enrich_pathways_<-c(up_enrich_pathways,down_enrich_pathways)
    }
    if(is.null(enrich_pathways_)){enrich_pathways_<-enrich_pathways}
  }
  enrichpathwaysplot<-list()
  for(enrich_pathway in enrich_pathways_){
    enrich_filename= paste(enrich_pathway,".tiff",sep="")
    tiff(filename = file.path(output_dir,enrich_filename),units="in",width=6,height=4,res=300,compression='lzw')
    plot_enrichment<-plotEnrichment(pathways[[enrich_pathway]],stats)+labs(title=enrich_pathway)
    print(plot_enrichment)
    dev.off()
    enrichpathwaysplot[[enrich_pathway]]<-paste(enrich_pathway,".tiff",sep="")
  }
  enrichment_map_data<-NULL
  if(enrichment_map){
    sig_genes<-rownames(statistics)[statistics[[pval_col]]<=pval_cutoff]
    enrichment_map_data<-fgseaRes[fgseaRes$padj<=enrichment_padj_cutoff & (!is.na(fgseaRes$padj)),]
    enrichment_map_data<-data.frame(ID=paste("pathway",1:nrow(enrichment_map_data),sep="_"),
                                    Term_Description=enrichment_map_data[["pathway"]],
                                    lowest_p=enrichment_map_data[["padj"]])
    changed_genes<-sapply(enrichment_map_data[["Term_Description"]],function(path){up_regulated=na.omit(pathway_gene_table[["Gene"]][pathway_gene_table[[path]]>=FC_cutoff]);
    up_regulated<-up_regulated[up_regulated %in% sig_genes];
    down_regulated=na.omit(pathway_gene_table[["Gene"]][pathway_gene_table[[path]]<=(-FC_cutoff)]);
    down_regulated<-down_regulated[down_regulated %in% sig_genes];
    return(c("Up_regulated"=paste(up_regulated,collapse = ", "),"Down_regulated"=paste(down_regulated,collapse = ", ")))
    })
    changed_genes<-as.data.frame(t(changed_genes))
    assertthat::are_equal(enrichment_map_data$Term_Description,rownames(changed_genes))
    enrichment_map_data<-cbind(enrichment_map_data,changed_genes)
    enrichment_map_data<-enrichment_map_data[(enrichment_map_data[["Up_regulated"]]!="") | (enrichment_map_data[["Down_regulated"]]!=""),]
    if(nrow(enrichment_map_data)>=2){
      tiff(filename = file.path(output_dir,"enrichment_map.tiff"),units="in", width=8, height=8, res=300, compression = 'lzw')
      enrichment_map_data<-cluster_enriched_terms(enrichment_map_data,use_description = T,plot_clusters_graph = T,plot_dend=T,kappa_threshold=kappa_cutoff)
      dev.off()
    } else{
      cat("rows of enrichment_map_data is less than 2. therefore, no enriched terms plot will be produced.")
    }
  }

  return(list(fgseatable=fgseaRes,sig_pathways=sig_pathways,pathway_gene_table=pathway_gene_table,enrichment_map_data=enrichment_map_data,pathway_gene_table_file=pathway_gene_table_file,fgseatableplot_file=fgseatableplot_file,enrichpathwaysplot=enrichpathwaysplot))
}

pathway_gene_table<-function(pathways,stats){
  stats<-data.frame(Gene_ID=names(stats),logFC=stats)
  genes_ordered_by_freq<-names(genes_freq<-sort(table(Reduce(c,pathways))[intersect(stats[["Gene_ID"]],unique(Reduce(c,pathways)))],decreasing = T))
  results<-data.frame(Gene=genes_ordered_by_freq)
  for(pathway in names(pathways)){
    genes_of_pathway<-intersect(genes_ordered_by_freq,pathways[[pathway]])
    results[[pathway]]=NA
    results[[pathway]][match(genes_of_pathway,genes_ordered_by_freq)]<-stats[["logFC"]][match(genes_of_pathway,stats[["Gene_ID"]])]
  }
  return(results)
}

fgseatableplot<-function (pathways, stats, fgseaRes, gseaParam = 1, colwidths = c(5, 3, 3, 1.2, 1.2), render = TRUE)
{
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathways <- lapply(pathways, function(p) {
    unname(as.vector(na.omit(match(p, names(statsAdj)))))
  })
  pathways <- pathways[sapply(pathways, length) > 0]
  ps <- lapply(names(pathways), function(pn) {
    p <- pathways[[pn]]
    annotation <- fgseaRes[match(pn, fgseaRes$pathway), ]
    list(textGrob(pn, just = "right", x = unit(0.95, "npc")),
         ggplot() +
           geom_segment(aes(x = p, xend = p, y = 0, yend = statsAdj[p]), size = 0.2) +
           scale_x_continuous(limits = c(0, length(statsAdj)), expand = c(0, 0)) +
           scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
           xlab(NULL) + ylab(NULL) +
           theme(panel.background = element_blank(), axis.line = element_blank(),
                 axis.text = element_blank(), axis.ticks = element_blank(),
                 panel.grid = element_blank(), axis.title = element_blank(),
                 plot.margin = rep(unit(0, "null"), 4), panel.spacing = rep(unit(0, "null"), 4)),
         ggplot()+geom_rect(aes(xmin=0,xmax=fgseaRes$NES[fgseaRes$pathway==pn],ymin=-0.6,ymax=0.6),fill=ifelse(sign(fgseaRes$NES[fgseaRes$pathway==pn])==1,"green","red")) +
           scale_x_continuous(limits = c(min(fgseaRes$NES,na.rm=T)*1.1, max(max(fgseaRes$NES,na.rm=T)*1.1,0)), expand = c(0, 0)) +
           scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
           xlab(NULL) + ylab(NULL) +
           theme(panel.background = element_blank(), axis.line = element_blank(),
                 axis.text = element_blank(), axis.ticks = element_blank(),
                 panel.grid = element_blank(), axis.title = element_blank(),
                 plot.margin = rep(unit(0, "null"), 4), panel.spacing = rep(unit(0, "null"), 4)),
         #         textGrob(sprintf("%.2f", annotation$NES)),
         textGrob(sprintf("%.1e", annotation$pval)),
         textGrob(sprintf("%.1e", annotation$padj)))
  })
  rankPlot_segment <- ggplot() + geom_blank() + scale_x_continuous(limits = c(0, length(statsAdj)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1,1), expand = c(0, 0)) + xlab(NULL) + ylab(NULL) +
    theme(panel.background = element_blank(), axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(),axis.title = element_blank(), plot.margin = unit(c(0, 0, 0.5, 0), "npc"), panel.spacing = unit(c(0, 0, 0, 0), "npc"))

  rankPlot_rect <- ggplot() + geom_blank() + scale_x_continuous(limits = c(min(fgseaRes$NES)*1.1, max(fgseaRes$NES)*1.1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1,1), expand = c(0, 0)) + xlab(NULL) + ylab(NULL) +
    theme(panel.background = element_blank(), axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(),axis.title = element_blank(), plot.margin = unit(c(0, 0, 0.5, 0), "npc"), panel.spacing = unit(c(0, 0, 0, 0), "npc"))

  grobs <- c(lapply(c("Pathway", "Gene ranks", "NES", "pval", "padj"), textGrob), unlist(ps, recursive = FALSE), list(nullGrob(), rankPlot_segment, rankPlot_rect, nullGrob(), nullGrob()))

  grobsToDraw <- rep(colwidths != 0, length(grobs)/length(colwidths))

  p <- arrangeGrob(grobs = grobs[grobsToDraw], ncol = sum(colwidths != 0), widths = colwidths[colwidths != 0])
  if (render) {
    grid.draw(p)
  }
  else {
    p
  }
}


#---- sankey_heatmap ----
sankey_heatmap<-function(sig_expressions,scale=TRUE,pathways,keep_other=TRUE,other_color="gray",pathway_colors,sankey_width=unit(10, "cm"),line_size=1,text_size=6,...){
  if(scale){
    sig_expressions<-t(scale(t(sig_expressions)))
  }
  if(missing(pathway_colors)){
    pathway_colors<-sample(c("#006400","#00008b","#b03060","#ff4500","#ffd700","#7fff00","#00ffff","#ff00ff","#6495ed","#ffdab9"),length(pathways))
    names(pathway_colors)<-names(pathways)
  }
  pathways<-data.frame(gene=unname(unlist(pathways)),pathway=rep(names(pathways),times=sapply(pathways,length)))
  sankey_df<-pathways[pathways[["gene"]] %in% rownames(sig_expressions),]
  sankey_df<-sankey_df[!duplicated(sankey_df),]
  if(keep_other){
    sankey_df_other<-data.frame(gene=setdiff(rownames(sig_expressions),sankey_df[["gene"]]),pathway="Other")
    sankey_df<-rbind(sankey_df,sankey_df_other)
    pathway_colors<-c(pathway_colors,"Other"=other_color)
  }
  sankey_anno = HeatmapAnnotation(sankey= anno_empty(border = F, width = sankey_width),which = "row")
  ht = Heatmap(sig_expressions, right_annotation = sankey_anno,...)
  ht = draw(ht)
  h_genes = data.frame(gene=rownames(sig_expressions)[row_order(ht)])
  h_genes[["h_order"]]=nrow(h_genes):1
  sankey_df[["h_order"]]<-h_genes[["h_order"]][match(sankey_df[["gene"]],h_genes[["gene"]])]
  sankey_df<-sankey_df[order(sankey_df[["pathway"]],sankey_df[["h_order"]],decreasing=T),]
  sankey_df[["point_x"]]=0.02
  sankey_df[["point_y"]]=sankey_df[["h_order"]]
  sankey_df[["point_color"]]<-pathway_colors[sankey_df[["pathway"]]]
  sankey_df[['point_color']][which(sankey_df[["gene"]] %in% sankey_df[["gene"]][duplicated(sankey_df[["gene"]])])]<-"black"
  sankey_df[["rect_x"]]=0.95
  sankey_df[["rect_y"]]=nrow(sankey_df):1/nrow(sankey_df)*nrow(sig_expressions)
  sankey_df[["rect_width"]]=0.05
  sankey_df[["rect_height"]]=1*nrow(sig_expressions)/nrow(sankey_df)
  sankey_df[["curve_x_start"]]<-0.05
  sankey_df[["curve_x_end"]]<-0.95
  sankey_df[["curve_y_start"]]<-sankey_df[["h_order"]]
  sankey_df[["curve_y_end"]]<-sankey_df[["rect_y"]]
  sankey_df[["text_x"]]=0.94
  sankey_df<-as.data.frame(sankey_df %>% group_by(pathway) %>% mutate(text_y=mean(rect_y)))
  sankey_df[["color"]]<-pathway_colors[sankey_df[["pathway"]]]
  decorate_annotation("sankey", {
    pushViewport(viewport(xscale = c(0, 1),yscale = c(0.5, nrow(sig_expressions)+0.5)))
    sigmoid=1 / (1 + exp(-seq(-5,5,length=100)))
    sigmoid<-(sigmoid-sigmoid[1])/(sigmoid[100]-sigmoid[1])
    for(i in 1:nrow(sankey_df)){
      grid.circle(x=sankey_df[i,"point_x"], y=sankey_df[i,"point_y"], r=unit(1,"mm"),gp = gpar(col="black",fill = sankey_df[i,"point_color"],alpha=0.5), default.units = "native")
      grid.lines(x=seq(sankey_df[i,"curve_x_start"],sankey_df[i,"curve_x_end"],length=100),y=(sankey_df[i,"curve_y_start"]+sigmoid)+sigmoid*(sankey_df[i,"curve_y_end"]-sankey_df[i,"curve_y_start"]-1),gp=gpar(col= sankey_df[i,"color"],alpha=0.2,lwd=line_size),default.units = "native")
      grid.rect(x=sankey_df[i,"rect_x"],y=sankey_df[i,"rect_y"],width=sankey_df[i,"rect_width"],height=sankey_df[i,"rect_height"],gp=gpar(fill=sankey_df[i,"color"],col=NA,alpha=1),just="left",default.units = "native")
      grid.text(label=sankey_df[i,"pathway"],x=sankey_df[i,"text_x"],y=sankey_df[i,"text_y"],just = "right",gp=gpar(fontsize=text_size),default.units = "native")
    }
    popViewport()
  })
}


# ---- help functions ----

#' reassign NA
#'
#' @param x vector
#' @param y a character or a vector to replace NA in x.
#'
#' @return new vector.
#' @export
#'
reassignNA=function(x,y){
  x[is.na(x)]=y
  return(x)}
