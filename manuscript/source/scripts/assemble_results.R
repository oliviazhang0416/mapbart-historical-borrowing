options(device=function(...) pdf(file=NULL, ...))
library(jsonlite)
library(ggplot2)
write_json <- function(x, path, ...) jsonlite::write_json(x, path, ..., digits=12)
script_arg <- grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)
if(length(script_arg)!=1L) stop("Run this manuscript builder with Rscript")
script_file <- normalizePath(sub("^--file=","",script_arg))
root <- normalizePath(file.path(dirname(script_file),"../../.."))
out <- file.path(root, "manuscript/source")
folders <- c("lrcbart-sim-gaussian","lrcbart-sim-survival",
             "lrcbart-sim-gaussian-single-arm","lrcbart-sim-survival-single-arm")
extract <- function(pattern, x, fallback="") {
  z <- regmatches(x,regexpr(pattern,x,perl=TRUE))
  if(length(z)) z else fallback
}
rows <- list(); region_rows <- list(); provenance <- list()
for(study in folders) {
  result_dirs <- if(grepl("survival-single-arm", study)) file.path(root,study,"res",c("n30","n200")) else file.path(root,study,"res")
  files <- unlist(lapply(result_dirs, list.files, pattern="\\.RData$", full.names=TRUE))
  for(f in files) {
    d <- readRDS(f); stopifnot(is.data.frame(d), nrow(d)==100, setequal(d$iteration,1:100))
    b <- basename(f); method <- sub("_p10.*","",b)
    scenario <- extract("sc[1-5](?:_cor-?[0-9.]+)?(?:_X[57])?(?:_d[0-9.]+)?",b)
    target <- sub("^_N","",extract("_N[^_]+",b))
    prior <- sub("^_prior","",extract("_prior[0-9.]+",b))
    w <- if(grepl("_w0.9_",b,fixed=TRUE)) "0.9" else "1"
    config <- if(method=="LRC-BART") {
      if(grepl("single-arm",study)) paste0(if(w=="0.9") "w=0.9; " else "",target) else target
    } else if(nzchar(prior)) prior else if(method=="MAP") target else ""
    nT <- if(grepl("survival-single-arm",study)) sub("^_n","",extract("_n[0-9]+",b)) else if(grepl("single-arm",study)) "200" else ""
    primary <- if(method %in% c("LRC-BART","MAP")) target=="100" && w=="1" else TRUE
    for(estimand in if("rmse_rmst"%in%names(d)) c("primary","RMST") else "primary") {
      suf <- if(estimand=="RMST") "_rmst" else ""
      val <- function(nm) d[[paste0(nm,suf)]]
      bias <- val("bias"); mse <- val("rmse"); cv<-val("coverage"); pw<-val("tp")
      stopifnot(all(is.finite(bias)),all(is.finite(mse)),all(mse>=0),all(is.finite(cv)),all(is.finite(pw)))
      r <- list(study=study,file=b,scenario=scenario,method=method,config=config,
        target=target,nT=nT,primary=primary,estimand=estimand,n=nrow(d),
        bias=mean(bias),rmse=sqrt(mean(mse)),sd=mean(val("sd")),width=mean(val("ci")),
        coverage=mean(cv),power=mean(pw),se_bias=sd(bias)/10,
        se_rmse=if(mean(mse)>0) sd(mse)/(20*sqrt(mean(mse))) else 0,
        se_coverage=sqrt(mean(cv)*(1-mean(cv))/100),
        se_power=sqrt(mean(pw)*(1-mean(pw))/100))
      if(method=="LRC-BART") for(nm in c("s0_sq","ess_prior","ess_ceiling","ess_target","ess_target_requested","ess_capped","w"))
        r[[nm]] <- mean(d[[nm]],na.rm=TRUE)
      rows[[length(rows)+1L]]<-r
    }
    if(method=="LRC-BART" && primary && grepl("^sc5_X5",scenario) && !grepl("single",study)) {
      for(region in c("Compatible","Incompatible")) for(metric in c("Probability of small discrepancy","Mean discrepancy")) {
        col <- paste0(if(metric=="Mean discrepancy") "g_" else "map_",if(region=="Compatible") "region" else "outside")
        v<-d[[col]]
        region_rows[[length(region_rows)+1L]] <- data.frame(
          outcome=if(grepl("gaussian",study)) "Gaussian" else "Survival",
          shift=as.numeric(sub(".*_d","",scenario)),region=region,metric=metric,
          estimate=mean(v,na.rm=TRUE),lower=quantile(v,.1,na.rm=TRUE),upper=quantile(v,.9,na.rm=TRUE),n=sum(is.finite(v)))
      }
    }
    provenance[[length(provenance)+1L]] <- list(path=f,md5=unname(tools::md5sum(f)),rows=nrow(d))
  }
}
write_json(rows,file.path(out,"generated/simulation_results.json"),auto_unbox=TRUE,pretty=TRUE,na="null")
write_json(provenance,file.path(out,"generated/result_manifest.json"),auto_unbox=TRUE,pretty=TRUE)
reg<-do.call(rbind,region_rows)
status <- system2(file.path(R.home("bin"),"Rscript"), shQuote(file.path(out,"scripts/build_figure1.R")))
stopifnot(status==0L)
write_json(reg,file.path(out,"generated/regional_results.json"),auto_unbox=TRUE,pretty=TRUE)

appdir<-file.path(root,"lrcbart-case-study-mm")
app<-readRDS(file.path(appdir,"res/results_table.RData"))
app$rmst_trt<-app$rmst_ctrl<-app$rmst_ucmm<-NA_real_
for(i in seq_len(nrow(app))){
 z<-readRDS(file.path(appdir,"res",app$file[i]))
 for(arm in c("trt","ctrl","ucmm")){
   nm<-paste0("rmst_",arm,"_est");v<-z[[nm]]
   if(!is.null(v)) app[i,paste0("rmst_",arm)]<-as.numeric(v[1])
 }
}
write_json(app,file.path(out,"generated/application_results.json"),dataframe="rows",auto_unbox=TRUE,pretty=TRUE,na="null")
e<-new.env();load(file.path(appdir,"data_cleaned/merged_elokrd_ucmm_n283.RData"),e);dat<-e$merged
cohort<-list()
for(k in c(1,0)){
 d<-dat[dat$trt==k,]; x<-list(cohort=if(k==1)"EloKRd" else "UCMM",n=nrow(d),
  age_median=median(d$age),age_q1=unname(quantile(d$age,.25)),age_q3=unname(quantile(d$age,.75)))
 for(v in c("male","race_Black","race_Other","hispanic","high_risk_cyto","asct","pfs_status","os_status"))
  x[[v]]<-sum(d[[v]]==1,na.rm=TRUE)
 x$race_White<-sum(d$race_Black==0 & d$race_Other==0)
 x$missing_model<-sum(!complete.cases(d[,c("age","male","race_Black","race_Other","hispanic","high_risk_cyto","asct")]))
 cohort[[length(cohort)+1]]<-x
}
write_json(cohort,file.path(out,"generated/cohort.json"),auto_unbox=TRUE,pretty=TRUE)
cat("Aggregated",length(provenance),"simulation files and",nrow(app),"application results.\n")

# Recover the historical prediction draws using the exact saved calibration settings.
# Canonical study checkpoints and result files are never overwritten.
library(Rcpp)
sourceCpp(file.path(root,"../lrcBART/clrcbart.cpp"),cacheDir=file.path(tempdir(),"manuscript_cpp"))
sourceCpp(file.path(root,"../lrcBART/cess.cpp"),cacheDir=file.path(tempdir(),"manuscript_cpp"))
oldcal<-readRDS(file.path(appdir,"res/ess/ess_pfs_n283_Hf10.RData"))
ess_data_file<-file.path(appdir,"data_cleaned/merged_elokrd_ucmm_n283.RData")
ess_outcome<-"PFS";ess_data_tag<-"n283"
ess_X_all<-as.matrix(dat[,c("age","male","race_Black","race_Other","hispanic","high_risk_cyto","asct")])
ess_time<-dat$pfs_months/12;ess_status<-dat$pfs_status;ess_trt<-dat$trt
ess_H_g<-oldcal$H_g;ess_H_f<-oldcal$H_f
ess_n_burn<-oldcal$settings$n_burn;ess_n_draw<-oldcal$settings$n_draw;ess_seed<-oldcal$settings$seed
ess_checkpoint_dir<-tempfile("manuscript_pfs_recovery_")
source(file.path(appdir,"ess_local/ess_cal.R"),local=TRUE)
recovery_error<-max(abs(ess_calibration$stage1_mu_draws-oldcal$stage1_mu_draws))
tr<-dat[dat$trt==1,]
groups<-data.frame(Age=cut(tr$age,c(-Inf,50,60,70,Inf),right=FALSE,labels=c("<50","50-59","60-69","70+")),
 Sex=ifelse(tr$male==1,"Male","Female"),
 Race=ifelse(tr$race_Black==1,"Black",ifelse(tr$race_Other==1,"Other","White")),
 Hispanic=ifelse(tr$hispanic==1,"Yes","No"),Risk=ifelse(tr$high_risk_cyto==1,"High","Standard"),
 ASCT=ifelse(tr$asct==1,"Yes","No"),stringsAsFactors=FALSE)
key<-apply(groups,1,paste,collapse="|");keys<-unique(key);gr<-list()
s0<-oldcal$s0_sq[["100"]];set.seed(460116);tau<-3*s0/rchisq(100000,3)
for(i in seq_along(keys)){
 idx<-which(key==keys[i]);Vf<-var(rowMeans(ess_stage1_f_test[,idx,drop=FALSE]))
 cg<-cess_cg(t(ess_stage1_xtest[idx,,drop=FALSE]),ess_cg_cutpoints,.5,3,10000L,12L,as.integer(1700+i))
 ess<-mean(oldcal$sigma1_sq/(Vf+5*cg*tau))
 row<-groups[idx[1],];row$n<-length(idx);row$Vf<-Vf;row$weight_sq_sum<-5*cg;row$ess<-ess
 gr[[i]]<-row
}
g<-do.call(rbind,gr);rownames(g)<-NULL
write_json(g,file.path(out,"generated/subgroup_ess.json"),dataframe="rows",auto_unbox=TRUE,pretty=TRUE)
write_json(list(age_bands=c("<50","50-59","60-69","70+"),historical_prediction_recovery_max_error=recovery_error,
 sigma1_sq=oldcal$sigma1_sq,s0_sq=s0,groups=nrow(g),tree_draws_per_group=10000,scale_draws=100000,
 historical_draws=nrow(ess_stage1_f_test),ESS_range=range(g$ess),overall_saved_calibration=oldcal$targets),
 file.path(out,"generated/subgroup_provenance.json"),auto_unbox=TRUE,pretty=TRUE)

status <- system2(file.path(R.home("bin"),"Rscript"), shQuote(file.path(out,"scripts/build_subgroup_table.R")))
stopifnot(status==0L)
cat("Subgroup table ESS:",nrow(g),"joint groups; recovered-draw maximum error",recovery_error,
    "; ESS range",range(g$ess),"\n")

# Calibration illustrations use fixed replicate 1 of Sc1, not a selected favourable dataset.
curves<-list();curve_prov<-list()
for(study in folders){
 fs<-list.files(file.path(root,study,"res","ess"),recursive=TRUE,pattern="\\.RData$",full.names=TRUE)
 fs<-fs[grepl("^ess_",basename(fs)) & grepl("sc1_",basename(fs)) & grepl("alternative_1_",basename(fs))]
 if(!length(fs)) next
 if(grepl("survival-single-arm",study)) fs<-fs[!duplicated(sub(".*(n30|n200).*","\\1",fs))] else fs<-fs[1]
 for(f in fs){
   cal<-readRDS(f)
   if(is.null(cal$ess_tau0_grid)) cal$ess_tau0_grid<-cal$ess0_grid
   lab<-switch(study, "lrcbart-sim-gaussian"="Gaussian, two-arm",
    "lrcbart-sim-survival"="Survival, two-arm","lrcbart-sim-gaussian-single-arm"="Gaussian, single-arm",
    paste("Survival, single-arm,",extract("n30|n200",f)))
   curves[[length(curves)+1]]<-data.frame(panel=lab,s=cal$grid,ess=cal$ess_tau0_grid,ceiling=cal$ceiling,method="lrcBART")
   curve_prov[[length(curve_prov)+1]]<-list(panel=lab,path=f)
   map_path<-file.path(dirname(f),paste0("map_ess_",tools::file_path_sans_ext(basename(cal$data_file)),"_checkpoint.RData"))
   if(file.exists(map_path)) {
     mp<-readRDS(map_path)
     curves[[length(curves)+1]]<-data.frame(panel=lab,s=mp$S_grid,ess=mp$ESS,ceiling=NA_real_,method="MAP")
     curve_prov[[length(curve_prov)+1]]<-list(panel=lab,path=map_path)
   }
 }
}
for(h in c(10,50)){
 cal<-readRDS(file.path(appdir,paste0("res/ess/ess_pfs_n283_Hf",h,".RData")))
 curves[[length(curves)+1]]<-data.frame(panel="Application PFS",s=cal$grid,ess=cal$ess_tau0_grid,ceiling=cal$ceiling,method=paste0("lrcBART, Hf=",h))
}
cc<-do.call(rbind,curves)
pp<-ggplot(cc,aes(s,ess,color=method))+geom_line(linewidth=.65)+
 geom_hline(data=unique(cc[is.finite(cc$ceiling),c("panel","ceiling")]),aes(yintercept=ceiling),linetype="dotted",color="grey45")+
 scale_x_log10()+facet_wrap(~panel,ncol=2,scales="free_y")+theme_bw(base_size=9)+
 scale_color_manual(values=c("lrcBART"="#176B87","MAP"="#CA6A36","lrcBART, Hf=10"="#176B87","lrcBART, Hf=50"="#735C91"),
 labels=c("lrcBART"=expression(lrcBART),"MAP"=expression(MAP),
 "lrcBART, Hf=10"=expression(lrcBART*", "*H[f]==10),
 "lrcBART, Hf=50"=expression(lrcBART*", "*H[f]==50)))+
 labs(x=expression(s[0]^2),y="Saved calibration ESS",color=NULL)+theme(legend.position="bottom")
ggsave(file.path(out,"figures/figureS1.pdf"),pp,width=7,height=8)
write_json(curve_prov,file.path(out,"generated/calibration_figure_provenance.json"),auto_unbox=TRUE,pretty=TRUE)

app_cal <- list()
for(ep in c("pfs","os")) for(h in c(10,50)) {
  cp<-file.path(appdir,paste0("res/ess/ess_",ep,"_n283_Hf",h,".RData"))
  z<-readRDS(cp)
  for(i in seq_len(nrow(z$targets))) {
    t<-as.list(z$targets[i,]);t$outcome<-toupper(ep);t$H_f<-h;t$ceiling<-z$ceiling
    t$V_f<-z$V_mu_f;t$sigma1_sq<-z$sigma1_sq;t$weight_sq_sum<-z$H_g*z$c_g
    t$whole_variance_ess<-mean(z$sigma1_sq/(z$V_mu_f+z$H_g*z$c_g*(3*t$s0_sq/qchisq((seq_len(100000)-.5)/100000,3))))
    app_cal[[length(app_cal)+1]]<-t
  }
}
write_json(app_cal,file.path(out,"generated/application_calibration.json"),auto_unbox=TRUE,pretty=TRUE)
