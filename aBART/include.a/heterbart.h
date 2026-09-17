
class heterbart : public bart
{
  public:
   heterbart():bart() { }
   heterbart(size_t m):bart(m) { }
   void draw(double *sigma, std::vector<double>& allLsum,
             rn& gen, int shards=1);
};

//--------------------------------------------------
void heterbart::draw(double *sigma, 
                     std::vector<double>& allLsum,
                     rn& gen, int shards)
{
   size_t i=0;
   for(size_t j=0;j<m;j++) {
     // cout << "      Tree " << j << "      "<< endl;
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
}
