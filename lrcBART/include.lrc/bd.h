#ifndef LRCBART_BD_H
#define LRCBART_BD_H

// ---------------------------------------------------------------------------
// The proposal is built on mBART/include.m/bd.h: it uses the same pointer tree,
// good-bottom-node search, ancestor cutpoint bounds, and in-place birth/death.
//
// LRC-BART MODIFICATION START
// Add ensemble-specific leaf likelihoods/minimum counts and the final
// grow-prune-change move set from the revision.
// ---------------------------------------------------------------------------
enum move_type { GROW=0, PRUNE=1, CHANGE=2 };

double move_probability(bool stump, int move, const pinfo& pi)
{
  if(stump) return move==GROW ? 1.0 : 0.0;
  if(move==GROW) return pi.p_grow;
  if(move==PRUNE) return pi.p_prune;
  return pi.p_change;
}

int pick_move(tree& t, const pinfo& pi, rn& gen)
{
  if(t.treesize()==1) return GROW;
  const double u=gen.uniform();
  if(u<pi.p_grow) return GROW;
  if(u<pi.p_grow+pi.p_prune) return PRUNE;
  return CHANGE;
}

double log_split_prior_ratio(const pinfo& pi, size_t depth,
                             double ps_left, double ps_right)
{
  const double ps=pi.alpha/std::pow(1.0+(double)depth,pi.mybeta);
  return std::log(ps)+std::log1p(-ps_left)+std::log1p(-ps_right)-std::log1p(-ps);
}

double pair_loglik(size_t nl1, size_t nl2, double sl1, double sl2,
                   size_t nr1, size_t nr2, double sr1, double sr2,
                   double sigma1_sq, double sigma2_sq, const pinfo& pi)
{
  return lh(nl1,nl2,sl1,sl2,sigma1_sq,sigma2_sq,pi)+
         lh(nr1,nr2,sr1,sr2,sigma1_sq,sigma2_sq,pi);
}

// Returns whether the proposal was accepted and reports which move was
// proposed so bart.h can preserve live acceptance diagnostics.
bool bd(tree& t, xinfo& xi, dinfo& di, pinfo& pi,
        double sigma1_sq, double sigma2_sq,
        rn& gen, int& proposed_move)
{
  proposed_move=pick_move(t,pi,gen);

  if(proposed_move==GROW) {
    tree::npv allbots;
    tree::npv goodbots;
    t.getbots(allbots);
    for(size_t i=0;i<allbots.size();i++)
      if(cansplit(allbots[i],xi)) goodbots.push_back(allbots[i]);
    if(goodbots.empty()) return false;

    tree::tree_p cand=goodbots[(size_t)std::floor(gen.uniform()*goodbots.size())];
    std::vector<size_t> goodvars;
    getgoodvars(cand,xi,goodvars);
    if(goodvars.empty()) return false;
    const size_t v=goodvars[(size_t)std::floor(gen.uniform()*goodvars.size())];
    int L=0;
    int U=(int)xi[v].size()-1;
    cand->rg(v,&L,&U);
    const size_t c=(size_t)(L+(int)std::floor(gen.uniform()*(U-L+1)));

    size_t nl1,nl2,nr1,nr2;
    double sl1,sl2,sr1,sr2;
    getsuff(t,cand,v,c,xi,di,
            nl1,nl2,sl1,sl2,nr1,nr2,sr1,sr2);
    if(!min_leaf_ok(nl1,nl2,pi) || !min_leaf_ok(nr1,nr2,pi)) return false;

    const double ll_new=pair_loglik(
      nl1,nl2,sl1,sl2,nr1,nr2,sr1,sr2,sigma1_sq,sigma2_sq,pi);
    const double ll_old=lh(
      nl1+nr1,nl2+nr2,sl1+sr1,sl2+sr2,sigma1_sq,sigma2_sq,pi);

    // Preserve mBART's exact good-bottom-node proposal count.  The child
    // split probabilities are zero when the proposed region has exhausted
    // its final admissible cutpoint.
    const double ps_child=pi.alpha/std::pow(2.0+(double)cand->depth(),pi.mybeta);
    const double ps_left=(goodvars.size()>1 || (int)c-1>=L) ? ps_child : 0.0;
    const double ps_right=(goodvars.size()>1 || U>=(int)c+1) ? ps_child : 0.0;
    const size_t nnog_old=t.nnogs();
    size_t nnog_new;
    if(cand->getp()==0) nnog_new=1;
    else if(cand->getp()->isnog()) nnog_new=nnog_old;
    else nnog_new=nnog_old+1;

    const bool old_stump=(t.treesize()==1);
    const double q_forward=move_probability(old_stump,GROW,pi);
    const double q_reverse=move_probability(false,PRUNE,pi);
    if(q_forward<=0.0 || q_reverse<=0.0 || nnog_new==0) return false;

    double log_alpha=ll_new-ll_old+
      log_split_prior_ratio(pi,cand->depth(),ps_left,ps_right);
    log_alpha+=std::log(q_reverse)-std::log(q_forward);
    log_alpha+=std::log((double)goodbots.size())-std::log((double)nnog_new);
    log_alpha=std::min(0.0,log_alpha);

    if(std::log(gen.uniform())<log_alpha) {
      cand->birthp(cand,v,c,0.0,0.0,1,1);
      return true;
    }
    return false;
  }

  tree::npv allnogs;
  t.getnogs(allnogs);
  if(allnogs.empty()) return false;
  tree::tree_p cand=allnogs[(size_t)std::floor(gen.uniform()*allnogs.size())];

  if(proposed_move==PRUNE) {
    size_t nl1,nl2,nr1,nr2;
    double sl1,sl2,sr1,sr2;
    getsuff(t,cand->getl(),cand->getr(),xi,di,
            nl1,nl2,sl1,sl2,nr1,nr2,sr1,sr2);

    const double ll_old=pair_loglik(
      nl1,nl2,sl1,sl2,nr1,nr2,sr1,sr2,sigma1_sq,sigma2_sq,pi);
    const double ll_new=lh(
      nl1+nr1,nl2+nr2,sl1+sr1,sl2+sr2,sigma1_sq,sigma2_sq,pi);

    const size_t nnog_old=t.nnogs();
    tree::npv current_leaves;
    t.getbots(current_leaves);
    size_t good_current=0;
    for(size_t l=0;l<current_leaves.size();l++)
      if(cansplit(current_leaves[l],xi)) ++good_current;
    size_t good_new=good_current;
    if(cansplit(cand->getl(),xi)) --good_new;
    if(cansplit(cand->getr(),xi)) --good_new;
    if(cansplit(cand,xi)) ++good_new;
    const bool new_stump=(cand->getp()==0);
    const double q_forward=move_probability(false,PRUNE,pi);
    const double q_reverse=move_probability(new_stump,GROW,pi);
    if(q_forward<=0.0 || q_reverse<=0.0 || good_new==0) return false;

    const double ps_child=pi.alpha/std::pow(2.0+(double)cand->depth(),pi.mybeta);
    const double ps_left=cansplit(cand->getl(),xi) ? ps_child : 0.0;
    const double ps_right=cansplit(cand->getr(),xi) ? ps_child : 0.0;

    double log_alpha=ll_new-ll_old-
      log_split_prior_ratio(pi,cand->depth(),ps_left,ps_right);
    log_alpha+=std::log(q_reverse)-std::log(q_forward);
    log_alpha+=std::log((double)nnog_old)-std::log((double)good_new);
    log_alpha=std::min(0.0,log_alpha);

    if(std::log(gen.uniform())<log_alpha) {
      cand->deathp(cand,0.0,1);
      return true;
    }
    return false;
  }

  // CHANGE: the revision changes a nog node's split rule.  The topology is
  // unchanged, and the uniform rule prior/proposal terms cancel.
  size_t onl1,onl2,onr1,onr2;
  double osl1,osl2,osr1,osr2;
  getsuff(t,cand->getl(),cand->getr(),xi,di,
          onl1,onl2,osl1,osl2,onr1,onr2,osr1,osr2);
  const double ll_old=pair_loglik(
    onl1,onl2,osl1,osl2,onr1,onr2,osr1,osr2,sigma1_sq,sigma2_sq,pi);

  std::vector<size_t> goodvars;
  getgoodvars(cand,xi,goodvars);
  if(goodvars.empty()) return false;
  const size_t old_v=cand->getv();
  const size_t old_c=cand->getc();
  const size_t v=goodvars[(size_t)std::floor(gen.uniform()*goodvars.size())];
  int L=0;
  int U=(int)xi[v].size()-1;
  cand->rg(v,&L,&U);
  const size_t c=(size_t)(L+(int)std::floor(gen.uniform()*(U-L+1)));
  cand->setrule(v,c);

  size_t nnl1,nnl2,nnr1,nnr2;
  double nsl1,nsl2,nsr1,nsr2;
  getsuff(t,cand->getl(),cand->getr(),xi,di,
          nnl1,nnl2,nsl1,nsl2,nnr1,nnr2,nsr1,nsr2);
  if(!min_leaf_ok(nnl1,nnl2,pi) || !min_leaf_ok(nnr1,nnr2,pi)) {
    cand->setrule(old_v,old_c);
    return false;
  }

  const double ll_new=pair_loglik(
    nnl1,nnl2,nsl1,nsl2,nnr1,nnr2,nsr1,nsr2,sigma1_sq,sigma2_sq,pi);
  const double log_alpha=std::min(0.0,ll_new-ll_old);
  if(std::log(gen.uniform())<log_alpha) return true;
  cand->setrule(old_v,old_c);
  return false;
}

// LRC-BART MODIFICATION END

#endif
