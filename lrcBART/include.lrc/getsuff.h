#ifndef LRCBART_GETSUFF_H
#define LRCBART_GETSUFF_H

// ---------------------------------------------------------------------------
// LRC-BART MODIFICATION START
// These are direct replacements for the two hetergetsuff() overloads in
// mBART/include.m/getsuff.h.  The original tree walk is retained.  The old
// two-mean/per-leaf-tau model needed weighted first and second moments;
// LRC-BART needs source counts and raw partial-residual sums for the exact f
// and g conjugate formulas.
// ---------------------------------------------------------------------------

// Sufficient statistics for a proposed grow: split observations currently in
// nx into proposed left and right children.
void getsuff(tree& x,
                tree::tree_p nx,
                size_t v, size_t c,
                xinfo& xi, dinfo& di,
                size_t& nl_1, size_t& nl_2,
                double& sl_1, double& sl_2,
                size_t& nr_1, size_t& nr_2,
                double& sr_1, double& sr_2)
{
  nl_1=nl_2=nr_1=nr_2=0;
  sl_1=sl_2=sr_1=sr_2=0.0;

  for(size_t i=0;i<di.n;i++) {
    double *xx=di.x+i*di.p;
    if(nx!=x.bn(xx,xi)) continue;

    const bool left=(xx[v] < xi[v][c]);
    if(di.source[i]==1) {
      if(left) { ++nl_1; sl_1+=di.y[i]; }
      else     { ++nr_1; sr_1+=di.y[i]; }
    } else if(di.source[i]==2) {
      if(left) { ++nl_2; sl_2+=di.y[i]; }
      else     { ++nr_2; sr_2+=di.y[i]; }
    }
  }
}

// Sufficient statistics for the two existing children of a nog node.  This
// overload is used by prune and by the before/after evaluations of change.
void getsuff(tree& x,
                tree::tree_p l, tree::tree_p r,
                xinfo& xi, dinfo& di,
                size_t& nl_1, size_t& nl_2,
                double& sl_1, double& sl_2,
                size_t& nr_1, size_t& nr_2,
                double& sr_1, double& sr_2)
{
  nl_1=nl_2=nr_1=nr_2=0;
  sl_1=sl_2=sr_1=sr_2=0.0;

  for(size_t i=0;i<di.n;i++) {
    double *xx=di.x+i*di.p;
    tree::tree_p bn=x.bn(xx,xi);

    if(bn==l) {
      if(di.source[i]==1) { ++nl_1; sl_1+=di.y[i]; }
      else if(di.source[i]==2) { ++nl_2; sl_2+=di.y[i]; }
    } else if(bn==r) {
      if(di.source[i]==1) { ++nr_1; sr_1+=di.y[i]; }
      else if(di.source[i]==2) { ++nr_2; sr_2+=di.y[i]; }
    }
  }
}

// The applicable minimum count differs by ensemble: pooled controls for f,
// RCT controls only for g.  A prior-only single-arm g forest sets n_min=0.
bool min_leaf_ok(size_t n1, size_t n2, const pinfo& pi)
{
  if(pi.kind==G_ENSEMBLE) return n1>=pi.n_min;
  return n1+n2>=pi.n_min;
}

// LRC-BART MODIFICATION END

#endif
