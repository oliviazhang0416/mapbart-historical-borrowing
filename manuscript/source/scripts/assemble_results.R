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
# Construct the application summary directly from saved model fits.
app_env <- new.env()
invisible(capture.output(sys.source(file.path(appdir,"summarize.R"), envir=app_env)))
app <- app_env$results_table
rm(app_env)
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
e<-new.env();load(file.path(appdir,"data_cleaned/merged_elokrd_ucmm_n230.RData"),e);dat<-e$merged
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

# Build the patient-profile PFS information table from the saved calibration.
status <- system2(file.path(R.home("bin"),"Rscript"), shQuote(file.path(out,"scripts/build_profile_table.R")))
stopifnot(status==0L)

# The calibration illustration uses fixed replicate 1 of Gaussian Sc1.
cal_path<-file.path(root,"lrcbart-sim-gaussian/res/ess/ess_data_p10_sc1_alternative_1_Hg5.RData")
cal<-readRDS(cal_path)
if(is.null(cal$ess_tau0_grid)) cal$ess_tau0_grid<-cal$ess0_grid
map_path<-file.path(dirname(cal_path),paste0("map_ess_",tools::file_path_sans_ext(basename(cal$data_file)),"_checkpoint.RData"))
mp<-readRDS(map_path)
cc<-rbind(
 data.frame(s=cal$grid,ess=cal$ess_tau0_grid,ceiling=cal$ceiling,method="lrcBART"),
 data.frame(s=mp$S_grid,ess=mp$ESS,ceiling=NA_real_,method="MAP")
)
pp<-ggplot(cc,aes(s,ess,color=method))+geom_line(linewidth=.8)+
 geom_hline(yintercept=cal$ceiling,linetype="dotted",color="grey45")+
 scale_x_log10()+theme_bw(base_size=11)+
 scale_color_manual(values=c("lrcBART"="#176B87","MAP"="#CA6A36"),
 labels=c("lrcBART"=expression(lrcBART),"MAP"=expression(MAP)))+
 labs(x=expression(s[0]^2),y="Prior ESS",color=NULL)+theme(legend.position="bottom")
ggsave(file.path(out,"figures/figureS1.pdf"),pp,width=6.5,height=4.2)
write_json(list(list(panel="Gaussian, two-arm",path=cal_path),
                list(panel="Gaussian, two-arm",path=map_path)),
           file.path(out,"generated/calibration_figure_provenance.json"),
           auto_unbox=TRUE,pretty=TRUE)

app_cal <- list()
for(ep in c("pfs","os")) for(h in c(10,50)) {
  cp<-file.path(appdir,paste0("res/ess/ess_",ep,"_n230_Hf",h,".RData"))
  z<-readRDS(cp)
  for(i in seq_len(nrow(z$targets))) {
    t<-as.list(z$targets[i,]);t$outcome<-toupper(ep);t$H_f<-h;t$ceiling<-z$ceiling
    t$V_f<-z$V_mu_f;t$sigma1_sq<-z$sigma1_sq;t$weight_sq_sum<-z$H_g*z$c_g
    t$whole_variance_ess<-mean(z$sigma1_sq/(z$V_mu_f+z$H_g*z$c_g*(3*t$s0_sq/qchisq((seq_len(100000)-.5)/100000,3))))
    app_cal[[length(app_cal)+1]]<-t
  }
}
write_json(app_cal,file.path(out,"generated/application_calibration.json"),auto_unbox=TRUE,pretty=TRUE)
