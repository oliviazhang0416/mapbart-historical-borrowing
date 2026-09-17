#ifndef LRCBART_INFO_H
#define LRCBART_INFO_H

// The dinfo/pinfo split is inherited from mBART/include.m/info.h.

// ---------------------------------------------------------------------------
// LRC-BART MODIFICATION START
// The original dinfo stored source membership in `s` and the original pinfo
// stored the shared-leaf inverse-gamma parameters.  LRC-BART keeps the same
// non-owning data/prior roles, but distinguishes the f and g ensembles and
// carries their final production priors.
// ---------------------------------------------------------------------------
enum ensemble_kind {
  F_ENSEMBLE = 0,
  G_ENSEMBLE = 1
};

class dinfo {
public:
  dinfo(): p(0), n(0), x(0), y(0), source(0) {}
  size_t p;
  size_t n;
  double *x;
  double *y;
  int *source; // 1 = RCT control, 2 = RWD control
};

class pinfo {
public:
  pinfo():
    pb(.5), alpha(.95), mybeta(2.0),
    p_grow(.3), p_prune(.3), p_change(.4),
    kind(F_ENSEMBLE), n_min(5),
    lambda_f_sq(1.0), tau0_sq(.01), tau1_sq(1.0), w(.5) {}

  // Kept for compatibility with the original grow/prune organization.
  double pb;
  double alpha;
  double mybeta;

  // Final revision move probabilities.
  double p_grow;
  double p_prune;
  double p_change;

  ensemble_kind kind;
  size_t n_min;

  // f leaf variance and g spike/slab parameters.
  double lambda_f_sq;
  double tau0_sq;
  double tau1_sq;
  double w;
};
// LRC-BART MODIFICATION END

#endif
