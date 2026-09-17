# Observed covariate combinations: cohort counts and unchanged PFS ESS.
library(jsonlite)
arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
script <- normalizePath(sub("^--file=","",arg))
out <- normalizePath(file.path(dirname(script),".."))
root <- normalizePath(file.path(out,"../.."))
data_file <- file.path(root,"lrcbart-case-study-mm/data_cleaned/merged_elokrd_ucmm_n283.RData")
e <- new.env();load(data_file,e);dat <- e$merged
levels_by_column <- list(Age=c("<50","50-59","60-69","70+"),
  Sex=c("Female","Male"),Race=c("White","Black","Other"),
  Hispanic=c("No","Yes"),Risk=c("Standard","High"),ASCT=c("No","Yes"))
groups <- data.frame(Age=cut(dat$age,c(-Inf,50,60,70,Inf),right=FALSE,labels=levels_by_column$Age),
  Sex=ifelse(dat$male==1,"Male","Female"),
  Race=ifelse(dat$race_Black==1,"Black",ifelse(dat$race_Other==1,"Other","White")),
  Hispanic=ifelse(dat$hispanic==1,"Yes","No"),Risk=ifelse(dat$high_risk_cyto==1,"High","Standard"),
  ASCT=ifelse(dat$asct==1,"Yes","No"),stringsAsFactors=FALSE)
stopifnot(!anyNA(groups),all(dat$trt %in% c(0,1)))
cols <- names(levels_by_column)
key <- function(d) apply(d[,cols,drop=FALSE],1,paste,collapse="|")
grid <- expand.grid(levels_by_column,stringsAsFactors=FALSE)
ranks <- lapply(cols,function(nm) match(grid[[nm]],levels_by_column[[nm]]))
grid <- grid[do.call(order,ranks),];rownames(grid)<-NULL
gkey <- key(grid);pkey <- key(groups)
grid$EloKRd_n <- as.integer(table(factor(pkey[dat$trt==1],levels=gkey)))
grid$UCMM_n <- as.integer(table(factor(pkey[dat$trt==0],levels=gkey)))
saved <- fromJSON(file.path(out,"generated/subgroup_ess.json"))
idx <- match(gkey,key(saved));grid$ess <- saved$ess[idx]
stopifnot(nrow(grid)==192L,!anyDuplicated(gkey),sum(grid$EloKRd_n)==30L,
  sum(grid$UCMM_n)==253L,sum(grid$EloKRd_n>0)==22L,
  identical(is.na(grid$ess),grid$EloKRd_n==0L),
  all(grid$EloKRd_n[!is.na(idx)]==saved$n[idx[!is.na(idx)]]),
  all(grid$ess[match(key(saved),gkey)]==saved$ess))
grid <- grid[grid$EloKRd_n>0,,drop=FALSE]
rownames(grid)<-NULL
stopifnot(nrow(grid)==22L,sum(grid$EloKRd_n)==30L,all(is.finite(grid$ess)))
write_json(grid,file.path(out,"generated/subgroup_table.json"),dataframe="rows",pretty=TRUE,na="null",digits=12)
write_json(list(rows=nrow(grid),EloKRd_n=sum(grid$EloKRd_n),UCMM_n=sum(grid$UCMM_n),
  UCMM_total_cohort_n=sum(dat$trt==0),ESS_rows=sum(!is.na(grid$ess)),UCMM_observed_combinations=sum(grid$UCMM_n>0),
  either_cohort_combinations=sum(grid$EloKRd_n+grid$UCMM_n>0),
  neither_cohort_combinations=sum(grid$EloKRd_n+grid$UCMM_n==0),
  ESS_values_unchanged=TRUE,data_md5=unname(tools::md5sum(data_file)),display_order=levels_by_column),
  file.path(out,"generated/subgroup_table_validation.json"),auto_unbox=TRUE,pretty=TRUE)
header <- paste(c("Age band","Sex","Race", "\\shortstack{Hispanic/\\\\Latino}",
  "\\shortstack{Cytogenetic\\\\risk}","ASCT","\\shortstack{EloKRd\\\\$n$}",
  "\\shortstack{UCMM\\\\$n$}","\\shortstack{Historical prior\\\\ESS}"),collapse=" & ")
note <- paste0("Only the 22 combinations observed in EloKRd are shown, ordered by the covariate columns. ",
  "Counts are observed patients within each combination. The displayed counts sum to 30 EloKRd and ",sum(grid$UCMM_n),
  " UCMM patients; other UCMM combinations are omitted from this table. ",
  "Each ESS uses the full 253-patient UCMM fit, exact trial ages, and the same all-spike calibration law and globally selected PFS scale. ",
  "Units are equivalent uncensored normal-reference controls for mean control log-PFS. Subgroup ESS values are not additive and do not measure outcome compatibility. ASCT is postinduction.")
t <- c("\\begingroup","\\setlength{\\tabcolsep}{3pt}","\\renewcommand{\\arraystretch}{1.0}",
  "\\fontsize{10}{11.5}\\selectfont","\\begin{longtable}{llllllrrr}",
  "\\caption{Observed cohort counts and historical prior ESS for joint covariate combinations observed in EloKRd (PFS).}\\label{tab:table7}\\\\",
  "\\toprule",paste0(header," \\\\ \\midrule"),"\\endfirsthead",
  "\\multicolumn{9}{l}{Table \\thetable\\ (continued)}\\\\","\\toprule",paste0(header," \\\\ \\midrule"),"\\endhead",
  "\\midrule\\multicolumn{9}{r}{Continued on next page}\\\\\\endfoot","\\bottomrule\\endlastfoot")
for(i in seq_len(nrow(grid))) {
  if(i>1 && (grid$Age[i]!=grid$Age[i-1] || grid$Sex[i]!=grid$Sex[i-1])) t<-c(t,"\\addlinespace[2pt]")
  row <- as.character(grid[i,cols])
  row[1] <- switch(row[1],"<50"="$<50$","50-59"="50--59","60-69"="60--69","70+"="$\\geq70$")
  t<-c(t,paste0(paste(c(row,grid$EloKRd_n[i],grid$UCMM_n[i],
    if(is.na(grid$ess[i])) "---" else sprintf("%.1f",grid$ess[i])),collapse=" & ")," \\\\"))
}
t<-c(t,"\\end{longtable}",paste0("\\noindent\\begin{minipage}{\\linewidth}\\footnotesize ",note,"\\end{minipage}"),"\\endgroup")
writeLines(t,file.path(out,"generated/table7.tex"))
cat("Table 7: 22 EloKRd combinations; 30 EloKRd and",sum(grid$UCMM_n),"UCMM patients shown; unchanged ESS values.\n")
