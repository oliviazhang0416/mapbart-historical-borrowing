# Source dgp.R first. Extend region logic in this pilot only.
in_region_original<-in_region
in_region<-function(X,region) switch(region,
 RECT75=!(X[,5]>2 & X[,7]>2),
 ONE75=X[,5]<=2+qnorm(.75),
 XOR50=(X[,5]>2)==(X[,7]>2),
 ONE50=X[,5]>2,
 in_region_original(X,region))
gen_interaction<-function(region,delta,seed) {
 d<-gen_data(4,'gaussian',region=region,delta_rwd=delta,seed=seed,n_mc=0)
 d$constants$scenario_id<-paste0(region,'_d',delta)
 d
}
