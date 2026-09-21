"""Produce manuscript tables exclusively from the audited aggregate JSON files."""
import json, math
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
GEN = BASE / 'generated'
sim = json.loads((GEN/'simulation_results.json').read_text())
app = json.loads((GEN/'application_results.json').read_text())
cohort = json.loads((GEN/'cohort.json').read_text())
cal = json.loads((GEN/'application_calibration.json').read_text())

def fmt(x, n=3):
    return '--' if x is None or not math.isfinite(x) else f'{x:.{n}f}'

def tex(s):
    return str(s).replace('_', r'\_').replace('%', r'\%').replace('&', r'\&')

scenarios = ['sc1','sc2','sc3_cor-0.5','sc3_cor0','sc3_cor0.5','sc4_d1','sc4_d2',
             'sc5_X5_d0.5','sc5_X5_d1','sc5_X5_d2','sc5_X7_d1','sc5_X7_d2']
def scen(s):
    if s in ['sc1','sc2']: return 'Sc'+s[-1]
    if s.startswith('sc3'): return r'Sc3, $\rho='+s.split('cor')[1]+'$'
    if s.startswith('sc4'): return r'Sc4, $\delta='+s.split('_d')[1]+'$'
    a=s.split('_'); return 'Sc5, $X_'+a[1][-1]+r'$, $\delta='+a[2][1:]+'$'

order = ['LRC-BART','BARTv1','BARTv2','BARTv3','LMv1','LMv2','LMv3',
         'AFTv1','AFTv2','AFTv3','hierLM','hierAFT','PSCL','MAP']
target_order = ['100','75','50','f90','f50','f25','s0min']
def sortkey(x):
    return (int(x['nT'] or 0),scenarios.index(x['scenario']),order.index(x['method']),
            '0.9' in x['config'],target_order.index(x['target']) if x['target'] else 0,x['config'])

def method(x, sensitivity=False):
    m=x['method']; c=x['config']
    if 'single-arm' in x['study']:
        if m=='BARTv2': return 'BART-CP'
        if m in ['LMv2','AFTv2']: return m[:-2]+'-CP'
    if m=='LRC-BART':
        if not sensitivity: return 'lrcBART'
        c=c.replace('w=0.9; ',r'$w=0.9$, ')
        return 'lrcBART: '+c
    if m.startswith('hier'): return ('hier. LM' if m=='hierLM' else 'hier. AFT')
    if m=='MAP': return 'MAP'+(' ('+c+')' if sensitivity else '')
    if m=='PSCL': return 'PSCL'
    for a in ['BART','LM','AFT']:
        if m.startswith(a+'v'): return a+'-'+{'1':'NP','2':'CP','3':'PP'}[m[-1]]
    return tex(m)

def retained_single_arm_result(x):
    """Exclude the discontinued s_tau^2=0.5 single-arm hierarchical runs."""
    return not (x['study'] in ['lrcbart-sim-gaussian-single-arm',
                               'lrcbart-sim-survival-single-arm']
                and x['method'] in ['hierLM', 'hierAFT']
                and x['config'] == '0.5')

def table(name, caption, headers, rows, spec, note='', landscape=False, small=True):
    if name in ['table2','table3','tableS4','tableS5']:
        note += r' The hierarchical comparator uses discrepancy-variance prior scale $s_\tau^2=0.25$.'
    elif name in ['table4','tableS6','tableS7']:
        note += r' The single-arm hierarchical comparator uses discrepancy-variance prior scale $s_\tau^2=0.05$.'
    elif name == 'table6':
        note += r' The hier. AFT comparator uses hierarchical discrepancy-variance prior scale $s_\tau^2=0.05$, where $\tau_\alpha^2,\tau_\beta^2\sim\operatorname{IG}(3/2,3s_\tau^2/2)$.'
    elif name == 'tableS10':
        note += r' Within each endpoint, the first and second hier. AFT rows use hierarchical discrepancy-variance prior scales $s_\tau^2=0.05$ and $s_\tau^2=0.5$, respectively.'
    t=[]
    if landscape: t.append(r'\begin{landscape}')
    intros={
      'tableS4':('Gaussian two-arm results','All 12 settings and 17 configurations are reported. Main Table 2 selects Sc1, Sc2, Sc3 with $\\rho=0$, Sc4 with $\\delta=1$, and Sc5 outside $X_5>2$ with $\\delta=1$.'),
      'tableS5':('Survival two-arm results',"The primary median-survival ratio is followed by the supplementary three-year RMST ratio. Each arm uses its own fitted residual standard deviation, and the control RMST denominator is floored at 0.05 (Appendix D.5)."),
      'tableS6':('Gaussian single-arm results','Sc1 and Sc2 at $n_1=200$ include all five lrcBART priors and all three comparator settings.'),
      'tableS7':('Survival single-arm results','Both trial sizes and both scenarios are reported, first for median-survival ratios and then for five-year RMST ratios with floored control denominators (Appendix D.5).'),
      'tableS10':('Complete application results','All 58 PFS and OS results are included. Similarity among lrcBART estimates does not resolve prior dependence for unobserved trial-control outcomes.')}
    if name in intros:
        title,desc=intros[name];t.extend([r'\subsection{'+title+'}',r'\noindent '+desc])
    style = r'\fontsize{8.5}{10.3}\selectfont' if small else r'\small'
    stretch='1.13'
    if name in ['table2','table4','table6']: style=r'\fontsize{10}{11.5}\selectfont';stretch='1.04'
    if name=='table3': style=r'\fontsize{10}{11.5}\selectfont';stretch='1.04'
    if name=='tableS11b': stretch='1.0';style=r'\fontsize{8.5}{9.2}\selectfont'
    if name in ['tableS4','tableS10']: stretch='1.0'
    t.extend([r'\begingroup',r'\setlength{\tabcolsep}{4pt}',r'\renewcommand{\arraystretch}{1.13}',
              r'\renewcommand{\arraystretch}{'+stretch+'}',style,
              r'\begin{longtable}{'+spec+'}',r'\caption{'+caption+r'}\label{tab:'+name+r'}\\',
              r'\toprule',' & '.join(headers)+r' \\ \midrule',r'\endfirsthead',
              r'\multicolumn{'+str(len(headers))+r'}{l}{\tablename\ \thetable\ (continued)}\\',
              r'\toprule',' & '.join(headers)+r' \\ \midrule',r'\endhead',
              r'\midrule\multicolumn{'+str(len(headers))+r'}{r}{Continued on next page}\\\endfoot',
              r'\bottomrule\endlastfoot'])
    for row in rows:
        if row=='__PAGE__': t.append(r'\newpage')
        elif isinstance(row,str): t.append(r'\multicolumn{'+str(len(headers))+r'}{l}{\strut\textit{'+row+r'}}\\*')
        else: t.append(' & '.join(row)+r' \\')
    t.append(r'\end{longtable}')
    if note: t.append(r'\noindent\begin{minipage}{\linewidth}\footnotesize '+note+r'\end{minipage}')
    t.append(r'\endgroup')
    if landscape: t.append(r'\end{landscape}')
    (GEN/(name+'.tex')).write_text('\n'.join(t)+'\n')

selected=['sc1','sc2','sc3_cor0','sc4_d1','sc5_X5_d1']
metricnote=(r'Each row summarizes 100 alternative-hypothesis replicates. Coverage is for a 95\% interval. '
            r'Empirical power is the frequency of posterior probability of benefit exceeding 0.95 (PSCL: lower 95\% confidence limit above the threshold). Gaussian decisions use a threshold of 0.5 '
            r'and survival decisions use a ratio threshold of 1. RMSE is the square root of the mean squared error. '
            r'NP denotes no historical pooling (trial controls only); CP denotes complete pooling; PP denotes source-adjusted pooling using a source indicator. PP is distinct from the hierarchical comparator. '
            r'lrcBART requests ESS 100; the MAP comparator in Table 2 also requests 100. Full configurations and Monte Carlo standard errors appear in the appendix.')
for study,num in [('lrcbart-sim-gaussian',2),('lrcbart-sim-survival',3)]:
    rr=sorted([x for x in sim if x['study']==study and x['primary'] and x['estimand']=='primary' and x['scenario'] in selected],key=sortkey)
    rows=[];prev=None
    for x in rr:
        if x['scenario']!=prev:
            if x['scenario']=='sc4_d1': rows.append('__PAGE__')
            rows.append(scen(x['scenario']));prev=x['scenario']
        rows.append([method(x)]+[fmt(x[k]) for k in ['bias','rmse','coverage','power']])
    table('table'+str(num),('Gaussian mean-difference' if num==2 else 'Median-survival-ratio')+' performance in representative two-arm scenarios.',
          ['Method','Bias','RMSE','Coverage','Power'],rows,'lrrrr',metricnote)

rows=[]
for study,ns in [('lrcbart-sim-gaussian-single-arm',['200']),('lrcbart-sim-survival-single-arm',['200'])]:
    for n in ns:
        for sc in ['sc1','sc2']:
            rows.append(('Gaussian' if 'gaussian' in study else 'Survival')+', '+scen(sc))
            rr=sorted([x for x in sim if x['study']==study and x['nT']==n and x['scenario']==sc and x['primary'] and x['estimand']=='primary' and retained_single_arm_result(x)],key=sortkey)
            rows.extend([[method(x)]+[fmt(x[k]) for k in ['bias','rmse','coverage','power']] for x in rr])
table('table4','Single-arm performance for Gaussian mean differences and median-survival ratios.',
      ['Method','Bias','RMSE','Coverage','Power'],rows,'lrrrr',
      r'Each row summarizes 100 alternative-hypothesis replicates. Coverage is for a 95\% interval; empirical power is the frequency of posterior probability of benefit exceeding 0.95, using a threshold of 0.5 for Gaussian mean differences and 1 for survival ratios. '
      r'All single-arm lrcBART primary rows use fixed $w=1$, $H_f=50$, $H_g=5$ and requested ESS 100. The hierarchical discrepancy variance has prior $\operatorname{IG}(3/2,3s_\tau^2/2)$. In single-arm studies, CP uses all historical controls with full likelihood weight and no source discrepancy; there are no observed trial controls to pool. Full configurations and Monte Carlo standard errors appear in the appendix.')

rows=[['Participants']+[str(x['n']) for x in cohort],['Age (years), median (IQR)']+[f"{x['age_median']:.1f} ({x['age_q1']:.1f}, {x['age_q3']:.1f})" for x in cohort]]
for label,key in [('Male','male'),('White','race_White'),('Black','race_Black'),('Other race','race_Other'),('Hispanic or Latino','hispanic'),('High-risk cytogenetics','high_risk_cyto'),('ASCT','asct'),('Progression or death','pfs_status'),('All-cause death','os_status')]:
    if label=='White': rows.append(['Race','',''])
    if key.startswith('race_'): label=r'\quad '+label.replace('Other race','Other')
    rows.append([label]+[f"{x[key]} ({100*x[key]/x['n']:.1f})" for x in cohort])
table('table5','Characteristics of the harmonized multiple-myeloma analysis cohorts.',
      ['Characteristic','EloKRd','UCMM'],rows,'lrr',
      r'Entries are $n$ (\%) except where indicated. IQR, interquartile range; ASCT, autologous stem-cell transplantation. '
      r'ASCT is a postinduction characteristic and is not presented as a pretreatment baseline variable. All analysis covariates were complete in both cohorts.',small=False)

def appmethod(x):
    return 'lrcBART' if x['method']=='LRC-BART' else ('hier. AFT' if x['method']=='HierAFT' else ('BART-CP' if x['method'] in ['BART','Standard BART'] else x['method']))
def ci(x,v,lo,hi): return f"{x[v]:.3f} ({x[lo]:.3f}, {x[hi]:.3f})"
rows=[]
for ep in ['PFS','OS']:
    rows.append(ep)
    aa=[x for x in app if x['outcome']==ep
        and (x['method']!='LRC-BART' or x['config']=='default' and x['target']=='100')
        and (x['method']!='HierAFT' or abs(float(x['prior'])-0.05)<1e-12)]
    aa.sort(key=lambda x: {'KM': 0, 'lrcBART': 1, 'BART-CP': 2, 'AFT-CP': 3, 'hier. AFT': 4}.get(appmethod(x), 5))
    for x in aa:
        is_km = x['method'] == 'KM'
        rows.append([appmethod(x),ci(x,'estimate','lower','upper'),ci(x,'difference','difference_lower','difference_upper'),
                     fmt(x['rmst_trt']),r'---' if is_km else fmt(x['rmst_ctrl']),
                     fmt(x['rmst_ucmm']) if is_km else r'---'])
table('table6','Five-year restricted mean survival time (RMST) in the application.',
      ['Method',r'\shortstack{RMST ratio\\(95\% interval)}',r'\shortstack{RMST difference\\(95\% interval)}',r'\shortstack{EloKRd treatment\\RMST}',r'\shortstack{EloKRd hypothetical\\control RMST}',r'\shortstack{UCMM control\\RMST}'],rows,'lrrrrr',
      r'RMST and differences are in years. For KM, EloKRd treatment RMST and UCMM control RMST are observed cohort estimates. For adjusted methods, EloKRd treatment and hypothetical control RMSTs are model-based estimates standardized to EloKRd covariates; UCMM control RMST is not reported. A dash denotes an inapplicable estimate. '
      r'KM is an unadjusted comparison with confidence intervals; model-based intervals are equal-tailed credible intervals. '
      r'Arm summaries are posterior medians, except KM estimates. The median of a ratio or difference need not equal the ratio or difference of marginal medians. '
      r'Primary lrcBART uses $H_f=10$, $H_g=5$, $w=1$, requested ESS 30.',landscape=True)

fullnote=(r'Entries with parentheses are estimates (Monte Carlo standard errors). Post. SD is the mean posterior standard deviation, '
          r'or the estimated standard error for PSCL; width is the mean 95\% interval width. '
          r'Coverage and empirical power are proportions; Bayesian decisions require posterior probability of benefit above 0.95, whereas PSCL requires its lower 95\% confidence limit to exceed the benefit threshold. Each row has 100 replicates. Gaussian power uses a benefit threshold of 0.5; survival power uses a ratio threshold of 1. '
          r'Targets 100, 75 and 50 are absolute counts; f90, f50 and f25 are fractions of the calibration ceiling. '
          r's0min fixes $s_0^2=10^{-6}$ without target selection. NP, CP and PP denote trial-only, completely pooled and source-adjusted control fits. In single-arm studies, CP uses all historical controls without a discrepancy term. The hierarchical discrepancy variance has prior $\operatorname{IG}(3/2,3s_\tau^2/2)$.')
for study,num in [('lrcbart-sim-gaussian',4),('lrcbart-sim-survival',5),('lrcbart-sim-gaussian-single-arm',6),('lrcbart-sim-survival-single-arm',7)]:
    rows=[]
    for estimand in (['primary','RMST'] if 'survival' in study else ['primary']):
        rr=sorted([x for x in sim if x['study']==study and x['estimand']==estimand and retained_single_arm_result(x)],key=sortkey)
        prev=None
        for x in rr:
            group=(x['nT'],x['scenario'])
            if group!=prev:
                if num==4 and x['scenario']=='sc5_X7_d2': rows.append('__PAGE__')
                rows.append(('Denominator-floored RMST ratio: ' if estimand=='RMST' else ('Median-survival ratio: ' if 'survival' in study else 'Mean difference: '))+
                            ('$n_1='+x['nT']+'$, ' if x['nT'] else '')+scen(x['scenario']));prev=group
            rows.append([method(x,True)]+[fmt(x[k])+' ('+fmt(x['se_'+k])+')' for k in ['bias','rmse']]+
                        [fmt(x['sd']),fmt(x['width'])]+[fmt(x[k])+' ('+fmt(x['se_'+k])+')' for k in ['coverage','power']])
    table('tableS'+str(num),'Complete '+study.replace('lrcbart-sim-','').replace('-',' ').capitalize()+' simulation results.',
          ['Method / setting','Bias (MCSE)','RMSE (MCSE)','Post. SD','Width','Coverage (MCSE)','Power (MCSE)'],rows,'lrrrrrr',fullnote,landscape=True)

rows=[]
for study in ['lrcbart-sim-gaussian','lrcbart-sim-survival','lrcbart-sim-gaussian-single-arm','lrcbart-sim-survival-single-arm']:
    prev=None
    for x in sorted([x for x in sim if x['study']==study and x['method']=='LRC-BART' and x['estimand']=='primary'],key=sortkey):
        g=(study,x['nT'],x['scenario'])
        if g!=prev:
            rows.append(study.replace('lrcbart-sim-','').replace('-',' ')+(', $n_1='+x['nT']+'$' if x['nT'] else '')+', '+scen(x['scenario']));prev=g
        rows.append([method(x,True),fmt(x.get('ess_target_requested'),1),fmt(x.get('ess_target'),1),fmt(x.get('ess_prior'),1),
                     fmt(x.get('ess_ceiling'),1),fmt(x.get('ess_capped'),2),fmt(x.get('s0_sq'),6)])
table('tableS8','Replicate-averaged lrcBART prior ESS calibration summaries.',
      ['Setting','Requested','Capped target','Curve ESS','Ceiling','Cap fraction',r'$s_0^2$'],rows,'lrrrrrr',
      r'All quantities are means over 100 datasets, except cap fraction, which is the proportion flagged as capped. Curve ESS is the mean across block-specific ESS calculations at the selected scale. '
      r'The ceiling uses the variance across all historical posterior draws. These are numerical summaries of distinct rules (Appendix B.6); neither is a posterior borrowing count. '
      r's0min bypasses target selection. A dash denotes an unavailable or inapplicable quantity.',landscape=False)

rows=[]
for ep in ['PFS','OS']:
    rows.append(ep)
    for x in sorted([x for x in app if x['outcome']==ep], key=lambda x: {'KM':0,'lrcBART':1,'BART-CP':2,'AFT-CP':3,'hier. AFT':4}.get(appmethod(x),5)):
        setting=(x['config'].replace('_',', ')+', '+x['target']) if x['method']=='LRC-BART' else '--'
        is_km = x['method']=='KM'
        rows.append([appmethod(x),tex(setting),ci(x,'estimate','lower','upper'),ci(x,'difference','difference_lower','difference_upper'),
                     fmt(x['rmst_trt']),r'---' if is_km else fmt(x['rmst_ctrl']),
                     fmt(x['rmst_ucmm']) if is_km else r'---'])
table('tableS10','Complete application results, including all lrcBART sensitivity settings.',
      ['Method','Setting',r'Ratio (95\% interval)',r'Difference (95\% interval)',r'\shortstack{EloKRd\\treatment RMST}',r'\shortstack{EloKRd hypothetical\\control RMST}',r'\shortstack{UCMM control\\RMST}'],rows,'llrrrrr',
      r'All 58 method/endpoint results are included. Default: $H_f=10$, $H_g=5$, $w=1$. Hf50 changes $H_f$ to 50; w0.9 fixes $w=0.9$. '
      r'Application labels 100, 75 and 50 request ESS 30, 22.5 and 15, respectively; f90, f50 and f25 request the corresponding fraction of the calibration ceiling. '
      r'Intervals, units and arm-column definitions follow Table 6. Contrasts use the EloKRd hypothetical control for adjusted methods and the observed UCMM control for KM.',landscape=True)

rows=[]
for ep in ['PFS','OS']:
    for h in [10,50]:
        rows.append(ep+', $H_f='+str(h)+'$')
        for x in [x for x in cal if x['outcome']==ep and x['H_f']==h]:
            rows.append([x['target_name'],fmt(x['requested'],2),fmt(x['target'],2),fmt(x['s0_sq'],6),fmt(x['ess_tau0'],2),fmt(x['whole_variance_ess'],2),fmt(x['ceiling'],2)])
table('tableS11a','Application prior ESS calibration.',
      ['Target label','Requested','Capped target',r'$s_0^2$','Block-based ESS','Full-variance ESS','Ceiling'],rows,'lrrrrrr',
      r'Full-variance ESS evaluates the expectation-of-ratios definition with the variance across all historical posterior draws and tree-weight estimate, using 100,000 chi-square quantiles. '
      r'Block-based ESS is the block average used in the existing scale selector, with 4,000 random scale draws. Differences also reflect numerical integration error; the reciprocal-variance ordering applies when the same integration draws are used. The same calibration is used at $w=1$ and $w=0.9$; it remains an all-spike reference.',landscape=True)
rows=[]
for ep in ['PFS','OS']:
    rows.append(ep)
    for x in [x for x in app if x['outcome']==ep and x['method']=='LRC-BART']:
        rows.append([tex(x['config'])+', '+x['target'],fmt(x['rhat'],4),fmt(x['ess_bulk'],0),fmt(x['ess_tail'],0)])
table('tableS11b','MCMC diagnostics for the application log RMST ratio.',
      ['Setting',r'$\widehat R$','Bulk MCMC ESS','Tail MCMC ESS'],rows,'lrrr',
      r'Four chains, each with 2,000 warm-up and 2,000 retained draws. These effective sample sizes measure MCMC efficiency and are unrelated to prior ESS. '
      r'These scalar diagnostics alone do not establish adequate exploration of all tree structures or parameters.')
reported_sim = [x for x in sim if retained_single_arm_result(x)]
reported_files = len({(x['study'], x['file']) for x in reported_sim})
print(f'Generated 13 numerical table fragments from {reported_files} reported simulation files and {len(app)} application rows.')
