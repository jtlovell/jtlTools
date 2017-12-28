#' @title Parse genStat QTL model output
#' @name genstatEffect
#' @aliases genstatLod
#' @aliases manualDrop
#'
#' @description
#' \code{genstatEffect} calculates the effects from genstat models
#'
#' @param path The path to the files to be parsed
#' @param pattern The pattern of filenames to look for
#' @details parses QTL effects from genstat
#'
#' @rdname genstatEffect
#' @export
genstatEffect<-function(path, pattern = "QTL-"){
  qtlEff<-lapply(list.files(path, pattern = pattern), function(x){
    tmp<-xlsx::read.xlsx(paste0(path,x), sheetIndex = 1)
    x<-gsub("-","_",x,fixed=T)
    colnames(tmp)[1]<-"POS"
    id=tmp$POS
    ef<-tmp[,c(1,grep("EFF",colnames(tmp), fixed=T))]
    ef<-melt(ef, id.var = "POS")
    sd<-tmp[,c(1,grep("SD",colnames(tmp), fixed=T))]
    sd<-melt(sd, id.var = "POS")

    ef$variable<-qdap::mgsub(c("....",".","EFF"),c("_","",""),ef$variable, fixed=T)
    ef$SITE = splitText(ef$variable,num=2)
    ef$ParentCross = splitText(ef$variable)
    ef<-ef[,-which(colnames(ef)=="variable")]
    colnames(ef)[which(colnames(ef)=="value")]<-"effect"

    sd$variable<-qdap::mgsub(c("....",".","SD"),c("_","",""),sd$variable, fixed=T)
    sd$SITE = splitText(sd$variable,num=2)
    sd$ParentCross = splitText(sd$variable)
    sd<-sd[,-which(colnames(sd)=="variable")]
    colnames(sd)[which(colnames(sd)=="value")]<-"sd"
    out<-merge(ef,sd, by = c("SITE","POS","ParentCross"), all=T)
    out$pheno<-qdap::mgsub(c("QTL_",".xlsx"),c("",""),x,fixed=T)
    return(out)
  })
  return(do.call(rbind, qtlEff))
}

#'
#' @rdname genstatLod
#' @title genstatLod
#' @description
#'  \code{genstatLod} Extracts the -log10 P-value profile from genstat output
#'
#' @param path The path to the files to be parsed
#' @param pattern The pattern of filenames to look for
#' @param cross The R/qtl cross object to be populated.
#' @details parses LOD scores from genstat
#'
#' @export
genstatLod<-function(path, pattern = "LOD", cross){
  lods<-lapply(list.files(path,pattern = pattern), function(x){
    tmp<-xlsx::read.xlsx(paste0(path,x), sheetIndex = 1)
    if(ncol(tmp)==3) colnames(tmp)[3]<-gsub(".xlsx","",gsub("LOD-","",x, fixed=T),fixed=T)
    return(tmp)
  })
  lods<-lapply(lods,function(x) x[,c(1,3)])
  lods<-reshape::merge_recurse(lods, id = c("Chromosome"))
  colnames(lods)[1]<-c("marker.name")
  colnames(lods)<-gsub("-","_",colnames(lods), fixed=T)

  map<-pullMap(cross)
  map<-merge(map, lods, by = "marker.name")
  rownames(map)<-map$marker.name
  map$pos<-map$pos.female
  s1<-map
  s1<-s1[,c("chr","pos",colnames(lods)[-1])]
  s1<-s1[with(s1,order(chr,pos)),]
  class(s1)<-c("scanone", "data.frame")

  return(s1)
}

#'
#' @rdname manualDrop
#' @title manualDrop
#' @description
#'  \code{manualDrop} calculates CI from the -log10 P-value curves
#'
#' @param qtlMarkerList A list of marker names that represent the QTL of interest.
#' @param s1 The scanone-classed -log10 P-value curves (usually generated from genstatLod)
#'
#' @details manual calculation of confidence intervals
#'
#' @export
manualDrop<-function(s1, qtlMarkerList, drop = 1){
  out<-lapply(1:length(qtlMarkerList), function(x){
    pheno = names(qtlMarkerList)[x]
    mars = qtlMarkerList[[x]]
    om<-lapply(1:length(mars), function(y){
      qtlMarker<-mars[y]
      chr<-as.character(s1[qtlMarker,"chr"])
      tmp<-s1[s1$chr == chr,c("pos",pheno)]
      ind = which(rownames(tmp)==qtlMarker)
      up<-ind:1
      drop.up<-cumsum(-(diff(tmp[,pheno][up])))
      bound.up<-rownames(tmp)[up[min(which(drop.up>=drop))+1]]
      down<-ind:nrow(tmp)
      drop.down<-cumsum(-(diff(tmp[,pheno][down])))
      bound.down<-rownames(tmp)[down[min(which(drop.down>=drop))+1]]
      if(is.na(bound.up)){
        bound.up<-rownames(tmp)[1]
      }
      if(is.na(bound.down)){
        bound.down<-rownames(tmp)[nrow(tmp)]
      }
      return(data.frame(
        pheno = pheno, chr = chr, pos = tmp[qtlMarker,"pos"], peakLod = tmp[qtlMarker,pheno],
        drop = drop,
        marker.low = bound.up, maker.hi = bound.down,
        cM.low = tmp[bound.up,"pos"], cM.hi = tmp[bound.down,"pos"],
        stringsAsFactors=F))
    })
    return(do.call(rbind, om))
  })
  return(do.call(rbind,out))
}
