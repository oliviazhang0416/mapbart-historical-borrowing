

class heterbart : public bart
{
  public:
   heterbart():bart() { }
   heterbart(size_t m):bart(m) { }
   void pr();
   void draw(double *sigma, 
             std::vector<double>& allLsum,
             rn& gen, 
             int shards=1);
};

//--------------------------------------------------
void heterbart::pr()
{
   cout << "+++++heterbart object:\n";
   bart::pr();
}
//--------------------------------------------------
void heterbart::draw(double *sigma, 
                     std::vector<double>& allLsum,
                     rn& gen, 
                     int shards)
{
   size_t i=0;
   for(size_t j=0;j<m;j++) {
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) {
         allfit[k] = allfit[k]-ftemp[k];
         r[k] = y[k]-allfit[k];
      }
      if(heterbd(t[j],xi,di,pi,sigma,nv,pv,aug,gen,shards)) i++;
      heterdrmu(t[j],xi,di,pi,sigma,gen);
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
      
      hetergetdiff(t[j],allLsum[j]);
   }
   accept=i;
   
   if(dartOn) {
     if(grp) draw_s_grp(nv,lpv,theta,gen,grp,rho);
     else draw_s(nv,lpv,theta,gen);
     draw_theta0(const_theta,theta,lpv,a,b,rho,gen);
     for(size_t j=0;j<p;j++) pv[j]=::exp(lpv[j]);
   }
}

