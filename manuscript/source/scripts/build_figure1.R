# Recover profile summaries from the current primary Gaussian fits, then draw
# the revision manuscript's two-panel covariate display. Study files are read-only.
library(Rcpp)
library(jsonlite)
arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
script <- normalizePath(sub("^--file=", "", arg))
root <- normalizePath(file.path(dirname(script), "../../.."))
out <- file.path(root, "manuscript/source")
study <- file.path(root, "lrcbart-sim-gaussian")
n_rep <- 40L
scenarios <- c("sc5_X5_d0.5", "sc5_X5_d1", "sc5_X5_d2")
delta <- c(.5, 1, 2)
if (!("--plot-only" %in% commandArgs(TRUE))) {
sourceCpp(file.path(root, "../lrcBART/clrcbart.cpp"), cacheDir=file.path(tempdir(), "figure1_cpp"))
# Evaluate the unchanged control-fitting block, using the saved scale and seed.
code <- readLines(file.path(study, "lrcBART.R"))
start <- tail(grep('    x_names <- paste0', code, fixed=TRUE),1L)
stopifnot(length(start)==1L)
end <- grep('    ctrl$f_test <- ctrl$f_test + y_center', code, fixed=TRUE)
fit_code <- parse(text=code[start:(end-1L)])
cache <- file.path(out, "generated/figure1_profile_cache")
dir.create(cache, showWarnings=FALSE)
profiles <- list(); errors <- numeric()
for (sc in scenarios) {
  saved <- readRDS(file.path(study, "res", paste0("LRC-BART_p10_",sc,"_N100_alternative.RData")))
  for (r in seq_len(n_rep)) {
    data_file <- file.path(study,"data",paste0("data_p10_",sc,"_alternative_",r,".RData"))
    dat <- readRDS(data_file)
    settings <- readRDS(file.path(study,"res/cache",paste0("lrcbart_trt_p10_",sc,"_alternative_",r,".RData")))$settings
    this_seed <- settings$seed
    stopifnot(settings$ntree==50L,settings$ndpost==1000L,settings$nskip==1000L)
    row <- saved[saved$iteration==r,,drop=FALSE]
    signature <- list(data_md5=unname(tools::md5sum(data_file)),
      source_md5=unname(tools::md5sum(file.path(study,"lrcBART.R"))),
      sampler_md5=unname(tools::md5sum(file.path(root,"../lrcBART/clrcbart.cpp"))),
      s0_sq=row$s0_sq, seed=this_seed+101L)
    cf <- file.path(cache,paste0(sc,"_",r,".rds"))
    recovered <- if(file.exists(cf)) readRDS(cf) else NULL
    if(is.null(recovered) || !identical(recovered$signature,signature)) {
      env <- new.env(parent=globalenv())
      list2env(list(dat=dat,p_obs=10L,ntree=50L,H_g=5L,w_prior=c(1,1),w_fixed=NA_real_,
        nskip=1000L,ndpost=1000L,this_seed=this_seed,target_index=1L,
        target_row=data.frame(s0_sq=row$s0_sq)),env)
      eval(fit_code,env)
      keep <- dat$X[dat$X[,"D"]==1,"Z"]==0
      gg <- env$ctrl$g_test[,keep,drop=FALSE]
      p <- data.frame(scenario=sc,replicate=r,
        x5=dat$X[dat$X[,"D"]==1 & dat$X[,"Z"]==0,"X5"],
        map=colMeans(abs(gg)<.5),g=colMeans(gg))
      err <- abs(mean(p$map)-row$map_mean)
      if(grepl("sc5",sc)) {
        inside <- p$x5>2
        err <- max(err,abs(mean(p$map[inside])-row$map_region),
          abs(mean(p$map[!inside])-row$map_outside),
          abs(mean(p$g[inside])-row$g_region),abs(mean(p$g[!inside])-row$g_outside))
      }
      if(!is.finite(err) || err>=1e-10) stop("Profile recovery mismatch: ",sc," replicate ",r,
        "; max error ",err,"; map ",mean(p$map)," saved ",row$map_mean,
        "; seed ",this_seed+101L)
      recovered <- list(signature=signature,profiles=p,recovery_error=err)
      saveRDS(recovered,cf)
    }
    profiles[[length(profiles)+1L]] <- recovered$profiles
    errors <- c(errors,recovered$recovery_error)
    if(r%%5==0) cat(sc, r, "of", n_rep, "recovered\n")
  }
}
res <- do.call(rbind,profiles)
brk <- seq(-.5,4.5,by=.5)
res$bin <- cut(res$x5,brk,include.lowest=TRUE)
agg <- aggregate(cbind(map,g)~scenario+bin,res,mean)
agg$n <- aggregate(map~scenario+bin,res,length)$map
agg$mid <- ((brk[-length(brk)]+brk[-1])/2)[as.integer(agg$bin)]
write_json(agg,file.path(out,"generated/figure1_binned_profiles.json"),dataframe="rows",pretty=TRUE,digits=12)
write_json(list(replicates=n_rep,scenarios=scenarios,bin_breaks=brk,
  max_recovery_error=max(errors),trial_control_profiles=nrow(res),
  profiles_in_display_range=sum(!is.na(res$bin)),source="Current Gaussian primary fits; exact saved seeds and scales"),
  file.path(out,"generated/figure1_provenance.json"),auto_unbox=TRUE,pretty=TRUE,digits=12)
} else {
  agg <- read_json(file.path(out,"generated/figure1_binned_profiles.json"),simplifyVector=TRUE)
}
agg <- agg[agg$scenario %in% scenarios,]
brk <- seq(-.5,4.5,by=.5)
cols <- setNames(c("#777777","#0072B2","#D55E00"),scenarios)
shapes <- setNames(c(16,17,15),scenarios)
pdf(file.path(out,"figures/figure1.pdf"),width=7.4,height=4.5)
par(mfrow=c(1,2),mar=c(4,4.2,4,.8),oma=c(0,0,2.8,0),mgp=c(2.6,.7,0),cex=.9)
regions <- function() {
  abline(v=2,lty=3,col="grey55")
  mtext("Historical shift",side=3,line=.45,at=.75,cex=.8)
  mtext("No historical shift",side=3,line=.45,at=3.25,cex=.8)
}
plot(NA,xlim=range(brk),ylim=c(0,1),xlab=expression(X[5]),
  ylab=expression(P(group("{",abs(g(x))<0.5,"}")~"| data")),
  main="(a) Small-discrepancy probability",cex.main=.95)
regions()
for(sc in scenarios) {
 a<-agg[agg$scenario==sc,]
 lines(a$mid,a$map,col=cols[sc],lwd=2)
 points(a$mid,a$map,col=cols[sc],pch=shapes[sc],cex=.75)
}
plot(NA,xlim=range(brk),ylim=c(-2.3,.3),xlab=expression(X[5]),
  ylab=expression(paste("Discrepancy ",g(x))),
  main="(b) Mean discrepancy and truth",cex.main=.95)
regions()
for(i in seq_along(scenarios)) {
 sc<-scenarios[i]
 segments(brk[1],-delta[i],2,-delta[i],col=cols[sc],lty=2,lwd=1.3)
}
segments(2,0,tail(brk,1),0,col="grey25",lty=2,lwd=1.3)
for(sc in scenarios) {
 a<-agg[agg$scenario==sc,]
 lines(a$mid,a$g,col=cols[sc],lwd=2)
 points(a$mid,a$g,col=cols[sc],pch=shapes[sc],cex=.75)
}
legend("bottomright",legend=c("Posterior mean","True discrepancy"),lty=c(1,2),lwd=c(2,1.3),bty="n",cex=.8)
par(fig=c(0,1,0,1),new=TRUE,mar=c(0,0,0,0),oma=c(0,0,0,0))
plot.new()
legend("top",legend=c(expression(delta==0.5),expression(delta==1),expression(delta==2)),
 col=cols,lwd=2,pch=shapes,horiz=TRUE,bty="n",title="Historical shift (Sc5)",cex=.85)
dev.off()
cat("Figure 1 written from existing profile summaries.\n")
