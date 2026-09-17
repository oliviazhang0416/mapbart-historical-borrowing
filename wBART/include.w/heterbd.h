
bool heterbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, 
	     std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen, int shards=1);

bool heterbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, 
	     std::vector<size_t>& nv, std::vector<double>& pv, bool aug, rn& gen, int shards)
{
   tree::npv goodbots;  //nodes we could birth at (split on)
   double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x

   if(gen.uniform() < PBx) { //do birth or death

      //--------------------------------------------------
      //draw proposal
      tree::tree_p nx; //bottom node
      size_t v,c; //variable and cutpoint
      double pr; //part of metropolis ratio from proposal and prior
      bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,nv,pv,aug,gen);

      //--------------------------------------------------
      //compute sufficient statistics
      size_t nr,nl; //counts in proposed bots
      double bl,br; //sums of weights
      double Ml, Mr; //weighted sum of y in proposed bots
      hetergetsuff(x,nx,v,c,xi,di,nl,bl,Ml,nr,br,Mr,sigma,shards);
      // cout << "--- l ---" << endl;
      // cout << "bl = " << bl << endl;
      // cout << "Ml = " << Ml << endl;
      // cout << "--- r ---" << endl;
      // cout << "br = " << br << endl;
      // cout << "Mr = " << Mr << endl;
      //--------------------------------------------------
      //compute alpha
      double alpha=0.0, lalpha=0.0;
      double lhl, lhr, lht;
      if((nl>=5) && (nr>=5)) { //cludge?
        // cout << "--- l ---" << endl;
        lhl = heterlh(bl,Ml,pi.tau);
        // cout << "lhl = " << lhl << endl;
        
        // cout << "--- r ---" << endl;
         lhr = heterlh(br,Mr,pi.tau);
         // cout << "lhr = " << lhr << endl;
         
         // cout << "--- l+r ---" << endl;
         lht = heterlh(bl+br,Ml+Mr,pi.tau);
         // cout << "lht = " << lht << endl;
         
         alpha=1.0;
         lalpha = log(pr) + (lhl+lhr-lht); 
         lalpha = std::min(0.0,lalpha);
         
         // cout << "log(lik) = " << (lhl+lhr-lht) << endl;
         // cout << "log(pr) = " << log(pr) << endl;
      }

      //--------------------------------------------------
      //try metrop
      double mul,mur; //means for new bottom nodes, left and right
      double uu = gen.uniform();
      // cout << "log(uu) = " << log(uu) << endl;
      bool dostep = (alpha > 0) && (log(uu) < lalpha);
      if(dostep) {
         mul = heterdrawnodemu(bl,Ml,pi.tau,gen);
         mur = heterdrawnodemu(br,Mr,pi.tau,gen);
         x.birthp(nx,v,c,mul,mur);
	 nv[v]++;
	 
	 // cout << "     BIRTH YES " << "    "<< endl;
         return true;
      } else {
        
        // cout << "* BIRTH NO *" << endl;
         return false;
      }
   } else {
      //--------------------------------------------------
      //draw proposal
      double pr;  //part of metropolis ratio from proposal and prior
      tree::tree_p nx; //nog node to death at
      dprop(x,xi,pi,goodbots,PBx,nx,pr,gen);

      //--------------------------------------------------
      //compute sufficient statistics
      double br,bl; //sums of weights
      double Ml, Mr; //weighted sums of y
      hetergetsuff(x, nx->getl(), nx->getr(), xi, di, bl, Ml, br, Mr, sigma, shards);

      //--------------------------------------------------
      //compute alpha
      double lhl, lhr, lht;
      lhl = heterlh(bl,Ml,pi.tau);
      lhr = heterlh(br,Mr,pi.tau);
      lht = heterlh(bl+br,Ml+Mr,pi.tau);

      double lalpha = log(pr) + (lht - lhl - lhr);
      lalpha = std::min(0.0,lalpha);

      //--------------------------------------------------
      //try metrop
      //double a,b,s2,yb;
      double mu;
      if(log(gen.uniform()) < lalpha) {
         mu = heterdrawnodemu(bl+br,Ml+Mr,pi.tau,gen);
	 nv[nx->getv()]--;
         x.deathp(nx,mu);
         return true;
      } else {
         return false;
      }
   }
}

