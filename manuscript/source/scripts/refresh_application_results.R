# Refresh only the application from current saved results; preserve simulation and subgroup summaries.
options(device=function(...) pdf(file=NULL, ...))
library(jsonlite)
library(ggplot2)
write_json <- function(x, path, ...) jsonlite::write_json(x, path, ..., digits=12)
script_arg <- grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)
if(length(script_arg)!=1L) stop("Run this manuscript builder with Rscript")
script_file <- normalizePath(sub("^--file=","",script_arg))
root <- normalizePath(file.path(dirname(script_file),"../../.."))
out <- file.path(root, "manuscript/source")
appdir<-file.path(root,"lrcbart-case-study-mm")
app<-readRDS(file.path(appdir,"res/results_table.RData"))
app$rmst_trt<-app$rmst_ctrl<-app$rmst_ucmm<-NA_real_
if("rmst_ucmm_population" %in% names(app))
 app$rmst_ucmm_population<-ifelse(app$method=="KM","UCMM",NA_character_)
for(i in seq_len(nrow(app))){
 z<-readRDS(file.path(appdir,"res",app$file[i]))
 for(arm in c("trt","ctrl")){
   nm<-paste0("rmst_",arm,"_est");v<-z[[nm]]
   if(!is.null(v)) app[i,paste0("rmst_",arm)]<-as.numeric(v[1])
 }
 if(app$method[i]=="KM") app$rmst_ucmm[i]<-as.numeric(z$rmst_ctrl_est[1])
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

stopifnot(nrow(app)==58L, all(is.finite(app$rmst_ucmm[app$method=="KM"])),
          all(is.na(app$rmst_ucmm[app$method!="KM"])))
manifest <- lapply(app$file, function(f) {
 p <- file.path(appdir,"res",f)
 list(file=f,md5=unname(tools::md5sum(p)))
})
write_json(manifest,file.path(out,"generated/application_result_manifest.json"),auto_unbox=TRUE,pretty=TRUE)
cat("Refreshed",nrow(app),"application results and saved calibration summaries; no fits rerun.\n")
