check_low_expression<-function(expressions,low_expression_cutoff=1,low_expression_pct_outlier_factor=1.5){
  tatal_low_expression<-sum(expressions<low_expression_cutoff)
  average_pct_low_expression<-tatal_low_expression/(nrow(expressions)*ncol(expressions))*100
  cat("On average, every sample have ",average_pct_low_expression,"% low expression genes.\n\n")
  sample_pct_low_expression<-colSums(expressions<low_expression_cutoff)/nrow(expressions)*100
  outlier_sample_pct_low_expression<-sample_pct_low_expression[sample_pct_low_expression>average_pct_low_expression*low_expression_pct_outlier_factor]
  if(length(outlier_sample_pct_low_expression)==0){
    cat("distribution of low expression gene for every sample is similar.")
  } else{
    cat("sample",paste(names(outlier_sample_pct_low_expression),collapse = ","),"looks having too many low expression genes: ",paste(outlier_sample_pct_low_expression,collapse=","),"% respectively.\n\n")
  }
  return(
    list(sample_pct_low_expression=sample_pct_low_expression,
         outlier_samples=outlier_sample_pct_low_expression))
}

immu_deconvolution<-function(project_dir,gene_expressions,sample_info,ID_col,cf_algorithm='xcell',log2_transform_cf=TRUE,background=0.00001,mode="supervised",test_method=c("ttest","limma"),paired=F,group_col,contrasts,pvalue_cutoff=0.05,top_anno=NULL,...){
  require(immunedeconv)
  test_method<-match.arg(test_method)
  if(!dir.exists(project_dir)){
    dir.create(project_dir)
  }
  stopifnot(all(colnames(gene_expressions)==sample_info[[ID_col]]))
  cellular_fraction<-deconvolute(gene_expressions,cf_algorithm,tumor = T)
  cellular_fraction<-as.data.frame(cellular_fraction)
  rownames(cellular_fraction)<-cellular_fraction$cell_type
  cellular_fraction<-cellular_fraction[,-1]
  write.csv(cellular_fraction,file.path(project_dir,"cellular_fraction.csv"))
  cf<-cellular_fraction
  if(log2_transform_cf){
    log2_cellular_fraction<-log2(cellular_fraction+background)
    cf<-log2_cellular_fraction
    write.csv(log2_cellular_fraction,file.path(project_dir,"log2_cellular_fratcion.csv"))
  }

  if(mode=="supervised"){
    groups<-unlist(strsplit(contrasts,"-"))
    if(test_method=="ttest"){
      statistics=apply(cf,1,function(scf){df<-data.frame(value=as.numeric(scf),group=sample_info[[group_col]]);df<-df[df$group %in% groups,];test=t.test(value~group,data=df,paired=paired);return(test$p.value)})
    }
    if(test_method=='limma'){
      if(!log2_transform_cf){
        cf=log2(cf+background)
      }
      formula<-as.formula(paste("~0+",group_col,sep=""))
      design<-model.matrix(formula,data=sample_info)
      colnames(design)<-levels(factor(sample_info[[group_col]]))
      fit<-lmFit(cf,design)
      contrast.matrix<-makeContrasts(contrasts=contrasts,levels = design)
      fit <- contrasts.fit(fit,contrast.matrix)
      fit<-eBayes(fit)
      statistics<-topTable(fit,number=nrow(cf))
    }

    mean_<-list()
    for (group in groups){
      mean_[[paste(group,"_mean",sep="")]]<-rowMeans(cf[,sample_info[[group_col]]==group],na.rm=T)
    }
    mean_<-as.data.frame(mean_)
    if(test_method=="ttest"){
      results<-data.frame(mean_[names(statistics),],P.Value=statistics)
      results[["adj.P.Val"]]<-p.adjust(results[["P.Value"]],method='BH')
    }
    if(test_method=="limma"){
      results<-cbind(mean_[rownames(statistics),],statistics)
    }
    write.csv(results,file.path(project_dir,paste(group_col,test_method,"statistics.csv",sep="_")))
    stopifnot("rownames of cellular fraction should be same as rownames of statistics results"=all(rownames(results)==rownames(cf)))
  }
  filtered_cf<-cf[rowSums(cellular_fraction>background)>=(ncol(cf)/2),]
  if(mode=="supervised"){
    neg_log_pvalue=-log10(results[rownames(filtered_cf),][["P.Value"]])
    row_anno<-rowAnnotation(neg_log_pvalue=row_anno_barplot(neg_log_pvalue,gp=gpar(fill=ifelse(neg_log_pvalue>(-log10(pvalue_cutoff)),"red","gray")),ylim=c(0,max(c(neg_log_pvalue,-log10(pvalue_cutoff)))+0.5)),width=unit(3,"cm"),annotation_name_side="top",annotation_name_rot=0)
    h<-Heatmap(t(scale(t(filtered_cf))),name="cellular_fraction",cluster_columns = cluster_within_group(t(scale(t(filtered_cf))),sample_info[[group_col]]),show_row_names = T,show_column_dend = F,show_column_names = T,right_annotation = row_anno,top_annotation = top_anno,...)
    tiff(filename = file.path(project_dir,"immun_deconvolution_heatmap.tiff"),width = 10,height=8,units = "in",res=300,compression = "lzw")
    draw(h)
    xpos<-0.95/(max(c(neg_log_pvalue,-log10(pvalue_cutoff)))+0.5)*(-log10(pvalue_cutoff))
    decorate_annotation("neg_log_pvalue",{
      grid.lines(c(xpos,xpos),c(0,1),gp = gpar(lty = 2,col="red"))
    })
    dev.off()
    return(list(cf=cf,statistics=results))
  } else {
    h<-Heatmap(t(scale(t(filtered_cf))),name="cellular_fraction",show_row_names = T,show_column_names = T,top_annotation = top_anno,...)
    tiff(filename = file.path(project_dir,"immun_deconvolution_heatmap.tiff"),width = 10,height=8,units = "in",res=300,compression = "lzw")
    draw(h)
    dev.off()
    return(cf)
  }

}
