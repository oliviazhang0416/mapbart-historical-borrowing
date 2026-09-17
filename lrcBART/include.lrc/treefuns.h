#ifndef LRCBART_TREEFUNS_H
#define LRCBART_TREEFUNS_H

//--------------------------------------------------
// The traversal helpers are inherited from mBART/include.m/treefuns.h.
void fit(tree& t, xinfo& xi, size_t p, size_t n,
         double *x, double *s, double* fv)
{
  tree::tree_p bn;
  for(size_t i=0;i<n;i++) {
    bn=t.bn(x+i*p,xi);
    // LRC-BART MODIFICATION: both ensembles now use one scalar leaf value.
    // `s` remains in the signature so the inherited call pattern is stable.
    fv[i]=bn->gettheta();
  }
}

// ---------------------------------------------------------------------------
// LRC-BART ADDITION: expose reached g-leaf indicators for realized ESS and the
// optional all-spike diagnostic.  This uses the same traversal as fit().
// ---------------------------------------------------------------------------
void fitz(tree& t, xinfo& xi, size_t p, size_t n,
          double *x, int* zv)
{
  for(size_t i=0;i<n;i++) zv[i]=t.bn(x+i*p,xi)->getz();
}

//--------------------------------------------------
bool cansplit(tree::tree_p n, xinfo& xi)
{
   int L,U;
   bool v_found=false;
   size_t v=0;
   while(!v_found && v<xi.size()) {
      L=0;
      U=(int)xi[v].size()-1;
      n->rg(v,&L,&U);
      if(U>=L) v_found=true;
      v++;
   }
   return v_found;
}

//--------------------------------------------------
void getgoodvars(tree::tree_p n, xinfo& xi,
                 std::vector<size_t>& goodvars)
{
   goodvars.clear();
   for(size_t v=0;v<xi.size();v++) {
      int L=0;
      int U=(int)xi[v].size()-1;
      n->rg(v,&L,&U);
      if(U>=L) goodvars.push_back(v);
   }
}

#endif
