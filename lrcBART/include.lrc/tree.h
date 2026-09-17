#ifndef LRCBART_TREE_H
#define LRCBART_TREE_H

#include <map>
#include <cmath>
#include <cstddef>

//--------------------------------------------------
// These aliases and the pointer-based tree organization are inherited from
// mBART/include.m/tree.h.
typedef std::vector<double> vec_d;
typedef std::vector<vec_d> xinfo;

//--------------------------------------------------
class tree {
public:
   typedef tree* tree_p;
   typedef const tree* tree_cp;
   typedef std::vector<tree_p> npv;
   typedef std::vector<tree_cp> cnpv;

   // ------------------------------------------------------------------------
   // LRC-BART MODIFICATION START
   // mBART stored theta_RCT, theta_RWD, and a leaf-specific tau.  The final
   // LRC-BART model uses one scalar contribution in either ensemble and adds
   // z only for a g leaf (z=1 spike, z=0 slab).
   // ------------------------------------------------------------------------
   tree(): theta(0.0), z(1), v(0), c(0), p(0), l(0), r(0) {}
   ~tree() { delete l; delete r; }

   void settheta(double value) { theta=value; }
   void setz(int value) { z=value; }
   double gettheta() const { return theta; }
   int getz() const { return z; }

   // Needed by the change proposal while retaining the original tree class.
   void setrule(size_t variable, size_t cut) { v=variable; c=cut; }
   bool isleaf() const { return l==0; }
   // LRC-BART MODIFICATION END

   size_t getv() const {return v;}
   size_t getc() const {return c;}
   tree_p getp() {return p;}
   tree_p getl() {return l;}
   tree_p getr() {return r;}

   size_t treesize();
   size_t nnogs();
   size_t nbots();

   // ------------------------------------------------------------------------
   // LRC-BART MODIFICATION: birth/death keep their original in-place pointer
   // mechanics; only the terminal payload changes to (theta,z).
   // ------------------------------------------------------------------------
   void birthp(tree_p np, size_t v, size_t c,
               double thetal, double thetar, int zl=1, int zr=1);
   void deathp(tree_p nb, double theta, int z=1);

   void getbots(npv& bv);
   void getnogs(npv& nv);
   void getnodes(npv& v);
   void getnodes(cnpv& v) const;

   tree_p bn(double *x,xinfo& xi);
   void rg(size_t v, int* L, int* U);

   size_t nid() const;
   size_t depth();
   char ntype();
   bool isnog();
   size_t getbadcut(size_t v);

private:
   // LRC-BART MODIFICATION: scalar leaf payload plus spike/slab indicator.
   double theta;
   int z;
   size_t v;
   size_t c;
   tree_p p;
   tree_p l;
   tree_p r;
};

std::ostream& operator<<(std::ostream&, const tree&);
std::ostream& operator<<(std::ostream& os, const tree& t)
{
  tree::cnpv nds;
  t.getnodes(nds);
  os << nds.size() << std::endl;
  for(size_t i=0;i<nds.size();i++) {
    os << nds[i]->nid() << " ";
    os << nds[i]->getv() << " ";
    os << nds[i]->getc() << " ";
    // LRC-BART MODIFICATION: serialize the active scalar and z, not the old
    // pair of source means and per-leaf tau.
    os << nds[i]->gettheta() << " ";
    os << nds[i]->getz() << std::endl;
  }
  return os;
}

// The node numbering, traversal, and cutpoint-bound code below is retained
// from mBART/include.m/tree.h.
//--------------------
size_t tree::nid() const
{
   if(!p) return 1;
   if(this==p->l) return 2*(p->nid());
   return 2*(p->nid())+1;
}
//--------------------
size_t tree::depth()
{
   if(!p) return 0;
   return (1+p->depth());
}
//--------------------
size_t tree::treesize()
{
   if(l==0) return 1;
   return 1+l->treesize()+r->treesize();
}
//--------------------
char tree::ntype()
{
   if(!p) return 't';
   if(!l) return 'b';
   if(!(l->l) && !(r->l)) return 'n';
   return 'i';
}
//--------------------
bool tree::isnog()
{
   if(!l) return false;
   return !(l->l) && !(r->l);
}
//--------------------
size_t tree::nnogs()
{
   if(!l) return 0;
   if(l->l || r->l) return l->nnogs()+r->nnogs();
   return 1;
}
//--------------------
size_t tree::nbots()
{
   if(!l) return 1;
   return l->nbots()+r->nbots();
}
//--------------------
void tree::getbots(npv& bv)
{
   if(l) {
      l->getbots(bv);
      r->getbots(bv);
   } else {
      bv.push_back(this);
   }
}
//--------------------
void tree::getnogs(npv& nv)
{
   if(l) {
      if(l->l || r->l) {
         if(l->l) l->getnogs(nv);
         if(r->l) r->getnogs(nv);
      } else {
         nv.push_back(this);
      }
   }
}
//--------------------
void tree::getnodes(npv& v)
{
   v.push_back(this);
   if(l) {
      l->getnodes(v);
      r->getnodes(v);
   }
}
void tree::getnodes(cnpv& v) const
{
   v.push_back(this);
   if(l) {
      l->getnodes(v);
      r->getnodes(v);
   }
}
//--------------------
tree::tree_p tree::bn(double *x,xinfo& xi)
{
   if(!l) return this;
   if(x[v] < xi[v][c]) return l->bn(x,xi);
   return r->bn(x,xi);
}
//--------------------
void tree::rg(size_t variable, int* L, int* U)
{
   if(!p) return;
   if(p->v == variable) {
      if(this == p->l) {
         if((int)(p->c) <= *U) *U=(int)p->c-1;
      } else {
         if((int)(p->c) >= *L) *L=(int)p->c+1;
      }
   }
   p->rg(variable,L,U);
}
//--------------------
void tree::birthp(tree_p np, size_t variable, size_t cut,
                  double thetal, double thetar, int zl, int zr)
{
   tree_p left = new tree;
   left->theta=thetal;
   left->z=zl;
   tree_p right = new tree;
   right->theta=thetar;
   right->z=zr;
   np->l=left;
   np->r=right;
   np->v=variable;
   np->c=cut;
   left->p=np;
   right->p=np;
}
//--------------------
void tree::deathp(tree_p nb, double value, int indicator)
{
   delete nb->l;
   delete nb->r;
   nb->l=0;
   nb->r=0;
   nb->v=0;
   nb->c=0;
   nb->theta=value;
   nb->z=indicator;
}
//--------------------
size_t tree::getbadcut(size_t variable)
{
   tree_p par=getp();
   if(par->getv()==variable) return par->getc();
   return par->getbadcut(variable);
}

#endif
