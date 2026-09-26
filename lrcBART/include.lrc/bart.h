#ifndef LRCBART_BART_H
#define LRCBART_BART_H

#include <ctime>

// ---------------------------------------------------------------------------
// This class is a direct structural descendant of mBART/include.m/bart.h: it
// still owns a vector of pointer-based trees, the common cutpoint grid, raw
// non-owning data pointers, allfit/r/ftemp buffers, and the per-tree Gibbs
// backfitting loop.
//
// LRC-BART MODIFICATION START
// One instance is designated f and another g.  The role controls source
// inclusion, leaf formulas, minimum counts, and the offset supplied by the
// other ensemble.
// ---------------------------------------------------------------------------
class bart {
public:
   bart();
   explicit bart(size_t m);
   ~bart();

   void setdata(size_t p, size_t n, double *x, double *y,
                int *source, size_t nc=100);
   void setdata(size_t p, size_t n, double *x, double *y,
                int *source, int* nc);

   void setprior(ensemble_kind kind,
                 double alpha, double beta,
                 double lambda_f_sq,
                 double tau0_sq, double tau1_sq, double w,
                 size_t n_min,
                 double p_grow=.3, double p_prune=.3)
   {
      pi.kind=kind;
      pi.alpha=alpha;
      pi.mybeta=beta;
      pi.lambda_f_sq=lambda_f_sq;
      pi.tau0_sq=tau0_sq;
      pi.tau1_sq=tau1_sq;
      pi.w=w;
      pi.n_min=n_min;
      pi.p_grow=p_grow;
      pi.p_prune=p_prune;
      pi.p_change=1.0-p_grow-p_prune;
   }

   void set_g_hyper(double tau0_sq, double tau1_sq, double w)
   {
      pi.tau0_sq=tau0_sq;
      pi.tau1_sq=tau1_sq;
      pi.w=w;
   }

   tree& gettree(size_t i) {return t[i];}
   const tree& gettree(size_t i) const {return t[i];}
   size_t getm() const {return m;}
   xinfo& getxinfo() {return xi;}
   void setxinfo(xinfo& _xi);

   // JOINT LRC-BART MODIFICATION START
   // offset is S*g + A*psi for f and f + A*psi for g in joint mode;
   // the original non-joint caller omits A*psi. Source-1 includes both arms.
   // K0/Lg/S0 accumulate for g and must be reset before every forest sweep.
   // JOINT LRC-BART MODIFICATION END
   void draw(double sigma1_sq, double sigma2_sq, const double *offset,
             int& K0, int& Lg, double& S0,
             rn& gen);

   void predict(size_t p, size_t n, double *x, double *fp);
   void predict_spike(size_t p, size_t n, double *x, int *nspike);
   double f(size_t i) const {return allfit[i];}
   const double* fitted() const {return allfit;}
   double getaccept(size_t move) const
   {
      return proposed[move]>0 ? (double)accepted[move]/proposed[move] : NA_REAL;
   }

private:
   size_t m;
   std::vector<tree> t;
   pinfo pi;
   size_t p,n;
   double *x;
   double *y;
   int *source;
   xinfo xi;
   double *allfit;
   double *r;
   double *ftemp;
   dinfo di;
   long accepted[3];
   long proposed[3];
};

//--------------------------------------------------
bart::bart():m(200),t(m),pi(),p(0),n(0),x(0),y(0),source(0),xi(),
             allfit(0),r(0),ftemp(0),di()
{
  for(int k=0;k<3;k++) accepted[k]=proposed[k]=0;
}

bart::bart(size_t im):m(im),t(m),pi(),p(0),n(0),x(0),y(0),source(0),xi(),
                       allfit(0),r(0),ftemp(0),di()
{
  for(int k=0;k<3;k++) accepted[k]=proposed[k]=0;
}

bart::~bart()
{
  delete[] allfit;
  delete[] r;
  delete[] ftemp;
}

//--------------------------------------------------
void bart::setxinfo(xinfo& _xi)
{
   const size_t pp=_xi.size();
   xi.resize(pp);
   for(size_t i=0;i<pp;i++) {
      xi[i].resize(_xi[i].size());
      for(size_t j=0;j<_xi[i].size();j++) xi[i][j]=_xi[i][j];
   }
}

//--------------------------------------------------
void bart::setdata(size_t p, size_t n, double *x, double *y,
                   int *source, size_t numcut)
{
   int* nc=new int[p];
   for(size_t i=0;i<p;i++) nc[i]=(int)numcut;
   setdata(p,n,x,y,source,nc);
   delete[] nc;
}

void bart::setdata(size_t p, size_t n, double *x, double *y,
                   int *source, int *nc)
{
   this->p=p;
   this->n=n;
   this->x=x;
   this->y=y;
   this->source=source;

   // Retain mBART's equally spaced fallback grid.  Production R workflows
   // supply the revision's explicit cutpoint grid through setxinfo().
   if(xi.empty()) {
      std::vector<double> minx(p,INFINITY), maxx(p,-INFINITY);
      for(size_t v=0;v<p;v++) {
         for(size_t i=0;i<n;i++) {
            const double value=*(x+p*i+v);
            if(value<minx[v]) minx[v]=value;
            if(value>maxx[v]) maxx[v]=value;
         }
      }
      xi.resize(p);
      for(size_t v=0;v<p;v++) {
         const double inc=(maxx[v]-minx[v])/(nc[v]+1.0);
         xi[v].resize(nc[v]);
         for(int j=0;j<nc[v];j++) xi[v][j]=minx[v]+(j+1)*inc;
      }
   }

   delete[] allfit;
   allfit=new double[n];
   predict(p,n,x,allfit);

   delete[] r;
   delete[] ftemp;
   r=new double[n];
   ftemp=new double[n];

   di.n=n;
   di.p=p;
   di.x=x;
   di.y=r;
   di.source=source;
}

//--------------------------------------------------
void bart::predict(size_t p, size_t n, double *x, double *fp)
{
   double *one=new double[n];
   for(size_t i=0;i<n;i++) fp[i]=0.0;
   for(size_t h=0;h<m;h++) {
      fit(t[h],xi,p,n,x,0,one);
      for(size_t i=0;i<n;i++) fp[i]+=one[i];
   }
   delete[] one;
}

void bart::predict_spike(size_t p, size_t n, double *x, int *nspike)
{
   int *one=new int[n];
   for(size_t i=0;i<n;i++) nspike[i]=0;
   for(size_t h=0;h<m;h++) {
      fitz(t[h],xi,p,n,x,one);
      for(size_t i=0;i<n;i++) nspike[i]+=one[i];
   }
   delete[] one;
}

//--------------------------------------------------
void bart::draw(double sigma1_sq, double sigma2_sq, const double *offset,
                int& K0, int& Lg, double& S0,
                rn& gen)
{
   for(size_t h=0;h<m;h++) {
      fit(t[h],xi,p,n,x,0,ftemp);
      for(size_t i=0;i<n;i++) {
         allfit[i]-=ftemp[i];
         // RWD rows are ignored by getsuff/drmu for g; setting their residual
         // to zero also makes that exclusion explicit at the forest level.
         if(pi.kind==G_ENSEMBLE && source[i]!=1) r[i]=0.0;
         else r[i]=y[i]-allfit[i]-(offset ? offset[i] : 0.0);
      }

      int move=-1;
      const bool didstep=bd(t[h],xi,di,pi,sigma1_sq,sigma2_sq,gen,move);
      if(move>=0 && move<3) {
         proposed[move]++;
         if(didstep) accepted[move]++;
      }

      drmu(t[h],xi,di,pi,sigma1_sq,sigma2_sq,K0,Lg,S0,gen);
      fit(t[h],xi,p,n,x,0,ftemp);
      for(size_t i=0;i<n;i++) allfit[i]+=ftemp[i];
   }
}

// LRC-BART MODIFICATION END

#endif
