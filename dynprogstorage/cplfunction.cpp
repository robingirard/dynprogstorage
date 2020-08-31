/*
 * cplfunction.hpp
 *
 *  Created on: 27 avr. 2013
 *      Author: robin
 */

#ifndef cplfunction_CPP_
#define cplfunction_CPP_

#include "cplfunction.hpp"
/*
 * cplfunction.hpp
 *
 *  Created on: 16 avr. 2013
 *      Author: robin
 */

std::vector<double> diffinv_(std::vector<double> x)
{
  int l=x.size();
  std::vector<double> res(l);
  res[0]=x[0];
  for(int i = 1; i < l; ++i) {
    res[i] = x[i]+res[i-1];
  }
  return(res);
}


bool isincreasing(std::vector<double> arg){
	int length=arg.size();
	  bool res=true;
	   for (int n=1; n<(length); n++)
		  if (arg[n]<=arg[n-1]){
			  res=false;
			  break;
		  }
	  return res;
}



// here polynom is 1/2 ax^2+bx+c
double getSlope(std::pair<double,double> Coefficients,double val){
// returns the slope at val given Coefficients a and b.f
	// Coefficients are (a,b)
   if (val==-std::numeric_limits<double>::infinity()&&Coefficients.first!=0){
  	 if (Coefficients.first<0){
			 return(std::numeric_limits<double>::infinity());
		 }else{
			 return(-std::numeric_limits<double>::infinity());
		 }
	 }else{
		 if (val==std::numeric_limits<double>::infinity()&&Coefficients.first!=0){
			 if (Coefficients.first<0){
				 return(-std::numeric_limits<double>::infinity());
			 }else{
				 return(std::numeric_limits<double>::infinity());
			 }
		 }else
		 {
			if (Coefficients.first==0){
				return(Coefficients.second);
			}
			else if (Coefficients.first==-std::numeric_limits<double>::infinity()){
				if (val<0)
				{
					return(std::numeric_limits<double>::infinity());
				}else
				{
					return(-std::numeric_limits<double>::infinity());
				}

			}else if (Coefficients.first==std::numeric_limits<double>::infinity())
			{
				if (val<0)
				{
					return(-std::numeric_limits<double>::infinity());
				}else
				{
					return(std::numeric_limits<double>::infinity());
				}
			}else
			{
				return(Coefficients.first*val+Coefficients.second);
			}
		 }
	 }
}

double getVal(std::pair<double,double> Coefficients,double val){
// returns the val at val given Coefficients a and b
	 if (val==-std::numeric_limits<double>::infinity()&&Coefficients.first!=0){
		 if (Coefficients.first<0){
			 return(std::numeric_limits<double>::infinity());
		 }else{
			 return(-std::numeric_limits<double>::infinity());
		 }
	 }else{
		 if (val==std::numeric_limits<double>::infinity()&&Coefficients.first!=0){
			 if (Coefficients.first<0){
				 return(-std::numeric_limits<double>::infinity());
			 }else{
				 return(std::numeric_limits<double>::infinity());
			 }
		 }else{
			 return(Coefficients.first/2*val*val+Coefficients.second*val);
		 }
	 }
}

double getXetoile(std::pair<double,double> Coefficients){
	 if (Coefficients.first==0){
		 if (Coefficients.second==0)
		 {
			 return(0);
		 }
		 else{
			 if (Coefficients.second<0){
				 return(std::numeric_limits<double>::infinity());
			 }else{
				 return(-std::numeric_limits<double>::infinity());
			 }
		 }
	 }else{
		 return(-Coefficients.second/Coefficients.first);
	 }
}

std::pair<double,double> Slopes2Coeffs(double Slope0,double Slope1){
  // returns the a and b coefficient of 1/2 ax^2+bx+c given the slopes in zero and the slopes in 1
  // a= S1-S0
  std::pair<double,double> res;
  res.first=Slope1-Slope0;
  res.second=Slope0;
  return(res);
}

std::vector<double> GetNextBreakVal(std::map<double, double>::iterator it,
                                    std::map<double, double>::iterator itplus1,
                                    std::map<double, double>::iterator itend,
                                    std::vector<double> curval, //size 3
                                    double FirstBreakVal
)
  //curval[0] BreakPoint, curval[1] is Val, curval[2] is Slope idem pour res
{
  std::vector<double> res(3);
  // determination du BreakPoint suivant
  if (itplus1==itend){
    //BreakPoint
    res[0]=std::numeric_limits<double>::infinity();
    //slope
    res[2]=std::numeric_limits<double>::infinity();
    //Value
    if (curval[2]<0){
      res[1]=-std::numeric_limits<double>::infinity();
    }else{
      if (curval[2]>0){
        res[1]=std::numeric_limits<double>::infinity();
      }else{// curval.second==0
        if (curval[1]==std::numeric_limits<double>::infinity()&&curval[0]==-std::numeric_limits<double>::infinity()){
          res[1]=FirstBreakVal;
        }else{
          res[1]=curval[1];
        }

      }
    }
  }else{
    //BreakPoint
    res[0]=itplus1->first;
    //slope
    res[2]=curval[2]+itplus1->second;
    //Value
    if (curval[1]==std::numeric_limits<double>::infinity()&&curval[0]==-std::numeric_limits<double>::infinity()){
      res[1]=FirstBreakVal+curval[2]*itplus1->first;
    }else{
      res[1]=curval[1]+curval[2]*itplus1->first;
      res[1]=curval[1];
    }
  }
  return(res);
}





    cplfunction::~cplfunction(){
      Breakpoints_.clear();
    }

    cplfunction::cplfunction()
    	: Breakpoints_(),FirstBreakVal_(0),
    	  FirstSlopeVal_(-std::numeric_limits<double>::infinity()){}


    class nonsuitableinput : public std::exception {
      public:
      const char * what() { return "non suitable input, if it is a dict, should have Slopes  and BreakPoints fields"; }
    };

    class emptyfunc : public std::exception {
      public:
      const char * what() { return "empty function"; }
    };

    class nonincreasingslopes : public std::exception {
      public:
      const char * what() { return "non increasing slopes"; }
    };

    class nonincreasingbreakpoints : public std::exception {
      public:
      const char * what() { return "non increasing breakpoints"; }
    };


        cplfunction::cplfunction(std::vector<double> Slopes, std::vector<double> BreakPoints,double FirstBreakVal){

  		int NbSlopes=  Slopes.size();
            FirstSlopeVal_=Slopes[0];
            Breakpoints_[BreakPoints[0]]=0;
  		if (NbSlopes==BreakPoints.size()){
  			if (isincreasing(Slopes)){
  				if (isincreasing(BreakPoints)){
  					FirstSlopeVal_=Slopes[0];
  					Breakpoints_[BreakPoints[0]]=0;
  					for (int i=1; i<NbSlopes; i++)
  					{
  						Breakpoints_[BreakPoints[i]]=Slopes[i]-Slopes[i-1];
  					}
  					FirstBreakVal_= FirstBreakVal;
  				}else{
  					//Rprintf( "Error: non increasing breakpoints" ) ;
  					throw nonincreasingbreakpoints() ;
  				}
  			}else{
  				//Rprintf( "Error: non increasing Slopes" ) ;
  				throw nonincreasingslopes() ;
  			}
  		}else{
  			//Rprintf( "Error: number of Slopes must be number of breaks+1  " ) ;
  			throw nonincreasingslopes() ;
  		}
	  }



	  cplfunction::cplfunction(cplfunction const & x) :
		  Breakpoints_(x.Breakpoints_),FirstBreakVal_(x.FirstBreakVal_),
		  FirstSlopeVal_(x.FirstSlopeVal_){}

    cplfunction* cplfunction::clone() const {
        return new cplfunction(*this) ;
    }

    cplfunction::cplfunction(double uniquebreak,double val)//a single point has 1 break and infinite FirstSlopeVal_
    	:Breakpoints_(),FirstBreakVal_(val),FirstSlopeVal_(std::numeric_limits<double>::infinity())
    {
    	Breakpoints_[uniquebreak]=0;
    }

    //a half line is just 1 point and a slope
    cplfunction::cplfunction(double uniquebreak,double val,double Slope1)
    	:Breakpoints_(),
    	 FirstBreakVal_(val),FirstSlopeVal_(Slope1)
    {
    	Breakpoints_[uniquebreak]=0;
    }


    //a "V" function is 2 points and 2 slopes the first point being infinity
    cplfunction::cplfunction(double uniquebreak,double val,double Slope1, double Slope2)
		:Breakpoints_(),
		 FirstBreakVal_(val),FirstSlopeVal_(Slope1)
    {
    	Breakpoints_[-std::numeric_limits<double>::infinity()]=0;
    	Breakpoints_[uniquebreak]=Slope2-Slope1;
    }

    cplfunction::cplfunction(double breakleft,double breakright,double val,double Slope1,double Slope2)
    	:Breakpoints_(),
    	 FirstBreakVal_(val),FirstSlopeVal_(Slope1)
    {
    	Breakpoints_[breakleft]=0;
    	Breakpoints_[breakright]=Slope2-Slope1;
    }


     std::map<double,double> cplfunction::get_BreakPoints()
     {
     	int nbBreaks=Breakpoints_.size();
     	//std::cout<<"nbBreaks"<<nbBreaks<<std::endl;
     	//std::vector<double> Breakpoints(nbBreaks);
  	 	//std::vector<double> Slopes(nbBreaks);
  	 	std::map<double,double> res;
  	 	double precslope=0;
  	 	int compteur=0;
  	 	//std::cout<<"Breakpoints_: "<< Breakpoints_.size()<<std::endl;
  	 	std::map<double,double>::iterator Breakpoints_it=Breakpoints_.begin();
  	 	while(Breakpoints_it != Breakpoints_.end())
  	 	{

  	 		if (compteur==0)
  	 		{
  	 		    res[Breakpoints_it->first]=FirstSlopeVal_;
  	 			precslope=FirstSlopeVal_;
  	 		}else
  	 		{
		 		res[Breakpoints_it->first]=precslope+Breakpoints_it->second;
		 		precslope=res[Breakpoints_it->first];
  	 		}
  	 		Breakpoints_it++; compteur++;
  	 	}
        return(res);
  	 }

    void cplfunction::print_()
      {
     	int nbBreaks=Breakpoints_.size();
     	std::vector<double> Breakpoints(nbBreaks);
   	 	std::vector<double> Slopes(nbBreaks);
   	 	int compteur=0;
   	 	std::map<double,double>::iterator Breakpoints_it=Breakpoints_.begin();
   	 	if (Breakpoints_it->second!=0)
   	 	{
   			std::cout<<"Warning first Slope diff non null =  "<< Breakpoints_it->second <<", ";
   	 	}

   	 	while(Breakpoints_it != Breakpoints_.end())
   	 	{
   			std::cout<<"|"<<Breakpoints_it->first<<"|";
   			if (compteur==0)
   			{
   	 			Slopes[compteur]=FirstSlopeVal_;
   	 			std::cout<<"__"<<FirstSlopeVal_<<"__";
   			}else
   			{
   	 			Slopes[compteur]=Slopes[compteur-1]+Breakpoints_it->second;
   	 			std::cout<<"__"<<Slopes[compteur]<<"__";
   			}
   			Breakpoints_it++; compteur++;
   	 	}
   	 	std::cout<<std::endl;
   	}



/*
    cplfunction(double * twobreaks,double slope, double val){
 	   int NbSlopes=1;
 	   double Slopes [2]={slope, numeric_limits<double>::infinity()};
 	   create_cplfunction(NbSlopes,Slopes,twobreaks,val);
    }
*/
    /* cplfunction(simplefunction sfunc){
 	   int NbSlopes=2;
 	   double Slopes [2]={sfunc.leftslope_, sfunc.rightslope_};
 	   double BreakPoints [3]={-numeric_limits<double>::infinity(),sfunc.breakpoint_,numeric_limits<double>::infinity()};
 	   create_cplfunction(NbSlopes,Slopes,BreakPoints,sfunc.val_);
    };*/

    void cplfunction::operator = (cplfunction s) {
     /* Cleanup current data */
     if(this != &s) {
      Breakpoints_.clear();
      /* copy needed data, call copy constructor
       * not efficient but will call copy constructor
       * */
      Breakpoints_=s.Breakpoints_;
      FirstBreakVal_=s.FirstBreakVal_;
      FirstSlopeVal_=s.FirstSlopeVal_;
     }
    }

    void cplfunction::AddSimple(double leftslope, double rightslope, double val, double breakpoint)
    {
    	//std::cout << __FUNCTION__ << "("<<leftslope<< ","<<rightslope<<","<<breakpoint<<")"<<"in "<<std::endl;
    	//this->print();
    	std::map<double, double>::iterator i = Breakpoints_.begin();
    	FirstBreakVal_=FirstBreakVal_+val;

		if (breakpoint<=Breakpoints_.begin()->first)
		{//BreakPoint is out of the domain, on the left
			  FirstSlopeVal_=FirstSlopeVal_+rightslope;
		}else
		{
		  if (breakpoint>=Breakpoints_.rbegin()->first && Breakpoints_.rbegin()->second== std::numeric_limits<double>::infinity()){
			  FirstSlopeVal_=FirstSlopeVal_+leftslope;
		  }else
		  {
			/*here the new breakpoint is inside the domain of this*/

			FirstSlopeVal_=FirstSlopeVal_+leftslope;
			double diff=rightslope-leftslope;
			std::pair<std::map<double, double>::iterator,bool> tmp_insert=Breakpoints_.insert(pair<double, double> (breakpoint, diff));
			if (!tmp_insert.second)
			{//insert the new breakpoint if it does not exist and if it exists increment :
				double tmpval=tmp_insert.first->second;
				(*tmp_insert.first).second=tmpval+rightslope-leftslope;
			}
		  }
		}
		//std::cout << __FUNCTION__ << "out "<<std::endl;
		//this->print();
    }

    bool cplfunction::eq(cplfunction  const & cplfunction1){
 	   if (FirstBreakVal_!=cplfunction1.FirstBreakVal_){
 		   return(false);
 	   }
 	   if (FirstSlopeVal_!=cplfunction1.FirstSlopeVal_){
 		  return(false);
 	   }
 	   if (Breakpoints_.size()!=cplfunction1.Breakpoints_.size()){
 		   return(false);
 	   }else{
   		   std::map<double, double>::iterator i = Breakpoints_.begin();
				std::map<double, double>::const_iterator	   i2=cplfunction1.Breakpoints_.begin();
   		   while(i != Breakpoints_.end()) {
   			   if (i->first==i2->first&&i->second==i2->second){
   				 ++i;++i2;
   			   }else{
   				   return(false);
   			   }
   		   }
   		   return(true);
 	   }
    }


double cplfunction::evalf(double x)
{
  //std::cout << __FUNCTION__ << "out "<<std::endl;
  //this->print();

  // cas très particulier où la fonction est réduite à un point
  // Dans ce cas FirstBreakVal_ est la valeur en ce point
  if (this->is_a_point()){
    if (x==Breakpoints_.begin()->first){
      return(FirstBreakVal_);
    }else{
      return(std::numeric_limits<double>::infinity());
    }
  }


  // si x est en dehors du domaine de definition --> infini
  if (x<Breakpoints_.begin()->first){
    return(std::numeric_limits<double>::infinity());
  }
  if ((x>Breakpoints_.rbegin()->first)&&(Breakpoints_.rbegin()->second==numeric_limits<double>::infinity())){
    return(std::numeric_limits<double>::infinity());
  }

  std::map<double, double>::iterator itend=Breakpoints_.end();
  std::map<double, double>::iterator it=Breakpoints_.begin();
  // cas simples d'une fonction linéaire
  std::map<double, double>::iterator itplus1 = Breakpoints_.begin(); ++itplus1;
  double curVal=FirstBreakVal_+FirstSlopeVal_*x;
  if (itplus1==itend){//un seul breakpoint
    return(FirstBreakVal_+FirstSlopeVal_*x);
  }

  //maintenant çà n'est ni une ligne ni un point donc il y a deux breakpoints

  double curBreakPoint=it->first;
  double curSlope=FirstSlopeVal_+it->second;
  double NextBreakPoint=itplus1->first;
  double NextBreakVal=FirstBreakVal_+curSlope*itplus1->first;
  bool done=false;

  if ((curBreakPoint<x)&&(x<=NextBreakVal)) done= true; // on garde curVal=FirstBreakVal_+FirstSlopeVal_*x;
  while (!done){
    // x_k<x,  on doit avoir itplus1!=itend normelement
    ++it;
    curBreakPoint=it->first;
    curSlope=curSlope+it->second;
    ++itplus1;
    if (itplus1==itend){
      curVal=NextBreakVal+curSlope*(x-curBreakPoint);
      done= true;
    }else{
      NextBreakPoint=itplus1->first;
      if (curBreakPoint<x&x<=NextBreakVal){
        done= true;
        curVal=NextBreakVal+curSlope*(x-curBreakPoint);
      }else{
        NextBreakVal=NextBreakVal+curSlope*(NextBreakPoint-curBreakPoint);
      }
    }
  }

  //curBreakPoint<x<
  return(curVal);

  }

//  boost::python::list evalDeltaf2(double x){
//      std::vector<double> res=evalDeltaf2_(x);
//       boost::python::list res_List=stl2py(res);
//      return(res_List);
//  }

     std::vector<double> cplfunction::evalDeltaf2_(double x)
    {
      // returns (f'(x-),f'(x+))
      std::vector<double> res(2);
      res[0]=FirstSlopeVal_;
      res[1]=FirstSlopeVal_;
      std::map<double,double>::iterator it = Breakpoints_.upper_bound(x);// it est le premier élément strictement plus grand que x
      std::map<double,double>::iterator it2 = Breakpoints_.begin();
      std::map<double,double>::iterator it2plus1 = Breakpoints_.begin();
      if(it != it2) {
        // le premier élément strictement plus grand que x n'est pas la première clés (x0)
        // donc xub > x >= x0  (xub =it->first)
        if (it2->first == x){// x==x0 res[0]=-infini res[1] = FirstSlopeVal_
          res[0] = -std::numeric_limits<double>::infinity();
        }else{
          while (it2!=it){
            // xit2 <= x< xub (xub =it->first)
            if (it2->first != x) res[0]=res[0]+it2->second; // si x == xit2 on garde le précédent
            res[1]=res[1]+it2->second;
            ++it2;
          }
        }
      }else{//  x est en dessous du premier breakpoint
        res[0]=-std::numeric_limits<double>::infinity();
        res[1]=-std::numeric_limits<double>::infinity();
      }

      return(res);
    }

    double cplfunction::evalDeltaf(double x)
    {
      double res=FirstSlopeVal_;
      std::map<double,double>::iterator it = Breakpoints_.upper_bound(x);
      std::map<double,double>::iterator it2 = Breakpoints_.begin();
      if(it != it2) {
        while (it2!=it){
          res=res+it2->second;
          ++it2;
        }
      }else{
        res=-std::numeric_limits<double>::infinity();
      }
      return(res);
    }

    bool cplfunction::is_last_infinity()
    {
    	if (Breakpoints_.size()==1)
    	{
    		return(FirstSlopeVal_!=std::numeric_limits<double>::infinity());
    	}else
    	{
    		return(((Breakpoints_.rbegin()->second!=std::numeric_limits<double>::infinity()) &&
    				FirstSlopeVal_!=std::numeric_limits<double>::infinity()));
    	}
    }

    bool cplfunction::is_a_point()
    {
    	return(FirstSlopeVal_==std::numeric_limits<double>::infinity());
    }

    bool cplfunction::is_an_infinite_line()
    {
    	return(FirstSlopeVal_!=std::numeric_limits<double>::infinity() &&
    			Breakpoints_.size()==1 &&
    			(Breakpoints_.begin())->first==-std::numeric_limits<double>::infinity());
    }

    // Etoile is replaced by "Legendre"
    void cplfunction::Etoile()
    {
    //	std::cout << __FUNCTION__<< " in" <<std::endl;
    //	this->print();
		cplfunction tmp(*this);
		Breakpoints_.clear();
		bool done=false;

		if (tmp.is_a_point())
		{// a point get transformed into a line
			FirstBreakVal_=-tmp.FirstBreakVal_;
			FirstSlopeVal_=(tmp.Breakpoints_.begin())->first;
			Breakpoints_[-std::numeric_limits<double>::infinity()]=0;
			done=true;
		}

		if (tmp.is_an_infinite_line())
		{// an infinite line get transformed into a point
			FirstBreakVal_=-tmp.FirstBreakVal_;
			Breakpoints_[(tmp.Breakpoints_.begin())->second+FirstSlopeVal_]=0;
			FirstSlopeVal_=std::numeric_limits<double>::infinity();
			done=true;
		}
		if (tmp.Breakpoints_.size()==1 && !done)
		{// only one breakpoint not a line not a point : this is a half line
		 // gives a half line with 2 breaks
			FirstBreakVal_=-tmp.FirstBreakVal_;
			FirstSlopeVal_=(tmp.Breakpoints_.begin())->first;
			Breakpoints_[-std::numeric_limits<double>::infinity()]=0;
			Breakpoints_[tmp.FirstSlopeVal_]=std::numeric_limits<double>::infinity();
			done=true;
		}
		if (done)
		{
			// do nothing
		}
		else
		{
			std::map<double,double>::iterator it=tmp.Breakpoints_.begin();
			std::map<double,double>::iterator itplus1=tmp.Breakpoints_.begin(); ++itplus1;
			FirstBreakVal_=-tmp.FirstBreakVal_;
			double NewBreak,NewSlope,PastBreakPrec,NewBreakPrec;
			if ((it->first!=-std::numeric_limits<double>::infinity()))
			{// first break is not infinity : this will give a slope and a first break with inf
				FirstSlopeVal_=it->first;
				Breakpoints_[-std::numeric_limits<double>::infinity()]=0;
			}else
			{
				FirstSlopeVal_=itplus1->first;
			}
			PastBreakPrec=FirstSlopeVal_;
			NewBreakPrec=tmp.FirstSlopeVal_;

			while (itplus1!=tmp.Breakpoints_.end())
			{
				NewSlope=itplus1->first-PastBreakPrec;
				PastBreakPrec=itplus1->first;
				NewBreak=it->second+NewBreakPrec;
				NewBreakPrec=NewBreak;
				Breakpoints_[NewBreak]=NewSlope;
				++it; ++itplus1;
			}
			if (tmp.is_last_infinity())
			{
				NewBreak=it->second+NewBreakPrec;
				Breakpoints_[NewBreak]=std::numeric_limits<double>::infinity();
			}
		}
	//	std::cout << __FUNCTION__<< " out" <<std::endl;
	//	this->print();
    }

    void cplfunction::EpiSum_Withline(double lowerbound,double upperbound,double Slope)
    {
    	//std::cout << __FUNCTION__<< " in" <<std::endl;
    	//this->print();
    	bool done=false;
		if (is_a_point())
		{// a point get transformed into a line
			double uniquebreak=Breakpoints_.begin()->first;
			Breakpoints_.clear();
			FirstSlopeVal_=Slope;
			Breakpoints_[lowerbound+uniquebreak]=0;
			Breakpoints_[upperbound+uniquebreak]=std::numeric_limits<double>::infinity();
			done=true;
		}

		if (is_an_infinite_line())
		{// an infinite line get transformed into a point a point remains the same and back it is the line
			// ... no change
			done=true;
		}
		if (Breakpoints_.size()==1 && !done)
		{// only one breakpoint not a line not a point : this is a half line
		 // gives a half line with 2 breaks
			double uniquebreak=Breakpoints_.begin()->first;
			double uniqueslope=Breakpoints_.begin()->second+FirstSlopeVal_;
			Breakpoints_.clear();
			if (Slope<uniqueslope)
			{// f*+g* has two breaks and two slopes
				FirstSlopeVal_=Slope;
				Breakpoints_[lowerbound+uniquebreak]=0;
				Breakpoints_[upperbound+uniquebreak]=uniqueslope-Slope;
			}else
			{
				FirstSlopeVal_=uniqueslope;
				Breakpoints_[lowerbound+uniquebreak]=0;
			}

			done=true;
		}
		if (done)
		{
			// do nothing
		}
		else
		{// there are 2 breaks or more.

	    	std::map<double,double>::iterator it=Breakpoints_.begin();
	    	std::map<double,double>::iterator itplus1=Breakpoints_.begin();++itplus1;
	    	std::map<double,double>::iterator itend=Breakpoints_.end();
	    	double CurrentSlope=FirstSlopeVal_;
	    	double newbreakval,newbreakslope,x;
	    	bool newbreak;
	    	if (CurrentSlope<Slope)
	    	{
				while (itplus1!=itend && (CurrentSlope+itplus1->second)<Slope )
				{// for the breakpoint that have slope < Slope just shift their key "on the left" by "lowerbound"
					CurrentSlope=CurrentSlope+itplus1->second;
					x=it->first;
					const_cast<double&>(it->first) = x+lowerbound;
					++it;++itplus1;
				}

		    	if (itplus1==itend)
		    	{
		    		newbreak=false;
		    	}else
		    	{//CurrentSlope<Slope <=CurrentSlope+itplus1->second
					x=it->first;
					const_cast<double&>(it->first) = x+lowerbound;
					if ((CurrentSlope+itplus1->second)==Slope)
					{
						newbreak=false;
					}else
					{
						newbreak=true;
						newbreakval=itplus1->first+lowerbound;
						newbreakslope=Slope-CurrentSlope;
						if (itplus1->second!=std::numeric_limits<double>::infinity())
						{
							itplus1->second=itplus1->second-newbreakslope;
						}

					}

					while (itplus1!=itend)
					{
						x=itplus1->first;
						const_cast<double&>(itplus1->first) = x+upperbound;
						++itplus1;
					}
		    	}
	    	}else
	    	{// Slope were lower or equal to the first Slope
	    		//Slope <=CurrentSlope
	    		if (it->first==-std::numeric_limits<double>::infinity())
	    		{
	    			newbreak=false;
					while (it!=itend)
					{
						x=it->first;
						const_cast<double&>(it->first) = x+upperbound;
						++it;
					}
	    		}else
	    		{//it->first>-numeric_limits<double>::infinity()
	    			if (CurrentSlope==Slope)
	    			{
	    				newbreak=false;
						x=it->first;
						const_cast<double&>(it->first) = x+lowerbound;
	    				while (itplus1!=itend)
	    				{
	    					x=itplus1->first;
	    					const_cast<double&>(itplus1->first) = x+upperbound;
	    					++itplus1;
	    				}
	    			}else
	    			{//CurrentSlope>Slope and it->first>-numeric_limits<double>::infinity()
						newbreak=true;
						newbreakval=it->first+lowerbound;
						it->second=it->second+FirstSlopeVal_-Slope; // from zero to true val
						FirstSlopeVal_=Slope;
						newbreakslope=0;

						while (it!=itend)
						{
							x=it->first;
							const_cast<double&>(it->first) = x+upperbound;
							++it;
						}
	    			}
	    		}
	    	}


	    	if (newbreak)
	    	{
	    		Breakpoints_[newbreakval]=newbreakslope;
	    	}

		}
		//std::cout << __FUNCTION__<< " out" <<std::endl;
		//this->print();
    }

    double cplfunction::Argmin(){
 	 // std::cout << __FUNCTION__ << std::endl;
 	 // this->print();
 	   double res;
 	   double precslope=FirstSlopeVal_;
 	   int NbSlopes=Breakpoints_.size();

 	   if (is_a_point()){res=Breakpoints_.begin()->first;}
 	   if (is_an_infinite_line()){
 		   if (FirstSlopeVal_==0){res=0;}
 		   else if (FirstSlopeVal_<0){res=std::numeric_limits<double>::infinity();}
 		   else if (FirstSlopeVal_>0){res=-std::numeric_limits<double>::infinity();}
 	   }

	   if (precslope>0){
		   res =Breakpoints_.begin()->first;
	   }else{
		 std::map<double, double>::iterator i = Breakpoints_.begin();
		 std::map<double, double>::iterator iend = Breakpoints_.end();
		++i;
		  while(i != iend) {
				res=i->first;
				precslope=precslope+i->second;
				if (precslope>0){ break;}
				++i;
				if(precslope<0 && is_last_infinity())
				{
					res=std::numeric_limits<double>::infinity();
				}
		  }
	   }
	//  std::cout<<"res="<<res<<endl;
	   return(res);

    }

    void cplfunction::Squeeze(double leftBreak,double rightBreak)
    {
     	  // std::cout << __FUNCTION__ << "("<<leftBreak<< ","<<rightBreak<<")"<<" in "<<std::endl;
     	//   std::cout<< "left : "<<leftBreak<<", right : "<<rightBreak<<std::endl;
     	 //  std::cout<< "this left "<< Breakpoints_.begin()->first << "this right : "<<Breakpoints_.rbegin()->first<<std::endl;
     	  //std::cout<<  "this right slope : "<<Breakpoints_.rbegin()->second<<std::endl;

     	  // this->print();
     	  // Rcout<<"FirstSlopeVal_ : "<<FirstSlopeVal_<<std::endl;
    	// test for empty interval or empty intersection of function and interval
    	bool done=false;
		if (  (leftBreak>rightBreak) ||
				(Breakpoints_.rbegin()->first<leftBreak && !is_last_infinity() ) ||
				(Breakpoints_.begin()->first>rightBreak && !is_last_infinity()) )
		{
			double arg_ = 1e-7;
			if ((leftBreak-Breakpoints_.rbegin()->first<arg_))
			{
				double tmpbreak=Breakpoints_.rbegin()->first;
				Breakpoints_.clear();
				Breakpoints_[tmpbreak]=std::numeric_limits<double>::infinity();
				done=true;
			}else if (Breakpoints_.begin()->first-rightBreak<arg_){
				double tmpbreak=Breakpoints_.begin()->first;
				Breakpoints_.clear();
				Breakpoints_[tmpbreak]=std::numeric_limits<double>::infinity();
				done=true;
			}else
			{
				if (leftBreak>=rightBreak){
				//std::cout<<"leftBreak>=rightBreak"<<std::endl;
        }
				//std::cout<<"Empty function thrown in Squeeze"<<std::endl;
				//std::cout << __FUNCTION__ << "("<<leftBreak<< ","<<rightBreak<<")"<<" in "<<std::endl;
				//this->print();
				//std::cout<< "Breakpoints_.begin()->first-rightBreak"<<Breakpoints_.begin()->first-rightBreak<<std::endl;
				throw emptyfunc();
			}

		}
		else if (is_a_point()||done)
		{

		}else
		{
			/// taking care of left break
			std::map<double, double>::iterator it=Breakpoints_.begin();
			std::map<double, double>::iterator itend=Breakpoints_.end();
			if (it->first<leftBreak)
			{// something will have to be cut on the left
				while (it!=itend && it->first<leftBreak )
				{
					FirstSlopeVal_=FirstSlopeVal_+it->second;
					++it;
				}
				if (it==itend)
				{
					Breakpoints_.erase(Breakpoints_.begin(),it);
					Breakpoints_.insert(std::pair<double, double> (leftBreak, 0.0));
				}else
				{
					if (it!=Breakpoints_.begin()){
						Breakpoints_.erase(Breakpoints_.begin(),it);
					}

					if (Breakpoints_.begin()->first==leftBreak){
						FirstSlopeVal_=FirstSlopeVal_+Breakpoints_.begin()->second;
						Breakpoints_.begin()->second=0.0;
					}else{
						Breakpoints_.insert(std::pair<double, double> (leftBreak, 0.0));
					}
				}
			}

			// taking care of right break
			if (is_last_infinity())
			{
				if (rightBreak!=std::numeric_limits<double>::infinity())
				{// something has to be cut on the right
					std::pair<std::map<double, double>::iterator,bool> tmp_insert=Breakpoints_.insert(std::pair<double, double> (rightBreak, numeric_limits<double>::infinity()));
					std::map<double, double>::iterator tmp_insert_it=tmp_insert.first;
					if (!tmp_insert.second)
					{
						(tmp_insert_it)->second=std::numeric_limits<double>::infinity();
					}

					++tmp_insert_it;
					if (tmp_insert_it!=Breakpoints_.end())
					{
						Breakpoints_.erase(tmp_insert_it,Breakpoints_.end());
					}
					if (Breakpoints_.size()==1)
					{
						FirstSlopeVal_=std::numeric_limits<double>::infinity();
					}
				}
			}else
			{
				if (rightBreak!=std::numeric_limits<double>::infinity() && Breakpoints_.rbegin()->first>rightBreak )
				{// something has to be cut on the right
					std::pair<std::map<double, double>::iterator,bool> tmp_insert=Breakpoints_.insert(std::pair<double, double> (rightBreak, numeric_limits<double>::infinity()));

					std::map<double, double>::iterator tmp_insert_it=tmp_insert.first;
					if (!tmp_insert.second)
					{
						(tmp_insert_it)->second=std::numeric_limits<double>::infinity();
					}

					++tmp_insert_it;
					if (tmp_insert_it!=Breakpoints_.end())
					{
						Breakpoints_.erase(tmp_insert_it,Breakpoints_.end());
					}
					if (Breakpoints_.size()==1)
					{
						Breakpoints_.begin()->second=0;
						FirstSlopeVal_=std::numeric_limits<double>::infinity();
					}
				}
			}
		}

		//std::cout << __FUNCTION__ << "("<<leftBreak<< ","<<rightBreak<<")"<<" out "<<std::endl;
		//this->print();
		//std::cout << "FirstSlopeVal_: "<<FirstSlopeVal_<<std::endl;
     }


    std::vector<double> cplfunction::Squeeze2(double leftBreak,double rightBreak)
        {
         	   //std::cout << __FUNCTION__ << "("<<leftBreak<< ","<<rightBreak<<")"<<" in "<<std::endl;
         	  // std::cout<< "left : "<<leftBreak<<", right : "<<rightBreak<<std::endl;
         	  // std::cout<< "this left "<< Breakpoints_.begin()->first << "this right : "<<Breakpoints_.rbegin()->first<<std::endl;
         	  //std::cout<<  "this right slope : "<<Breakpoints_.rbegin()->second<<std::endl;

         	   //this->print();
         	   //std::cout<<"FirstSlopeVal_ : "<<FirstSlopeVal_<<std::endl;
        	// test for empty interval or empty intersection of function and interval
    	std::vector<double> res(2);
    	res[0]=0; res[1]=0;
    		if (  (leftBreak>rightBreak) ||
    				(Breakpoints_.rbegin()->first<leftBreak && !is_last_infinity() ) ||
    				(Breakpoints_.begin()->first>rightBreak && !is_last_infinity()) )
    		{
    			if (leftBreak>=rightBreak){
          //  std::cout<<"leftBreak>=rightBreak"<<std::endl;
          }
    			//std::cout<<"Empty function thrown in Squeeze"<<std::endl;
    			//std::cout << __FUNCTION__ << "("<<leftBreak<< ","<<rightBreak<<")"<<" in "<<std::endl;
    			//this->print();
    			throw emptyfunc();
    		}
    		else if (is_a_point())
    		{

    		}else
    		{
    			/// taking care of left break
    			std::map<double, double>::iterator it=Breakpoints_.begin();
    			std::map<double, double>::iterator itend=Breakpoints_.end();
    			if (it->first<leftBreak)
    			{// something will have to be cut on the left
    				while (it!=itend && it->first<leftBreak )
    				{
    					FirstSlopeVal_=FirstSlopeVal_+it->second;
    					++it;
    				}
    				if (it==itend)
    				{
    					Breakpoints_.erase(Breakpoints_.begin(),it);
    					Breakpoints_.insert(std::pair<double, double> (leftBreak, 0.0));
    				}else
    				{
    					if (it!=Breakpoints_.begin()){
    						Breakpoints_.erase(Breakpoints_.begin(),it);
    					}

    					if (Breakpoints_.begin()->first==leftBreak){
    						FirstSlopeVal_=FirstSlopeVal_+Breakpoints_.begin()->second;
    						Breakpoints_.begin()->second=0.0;
    					}else{
    						Breakpoints_.insert(std::pair<double, double> (leftBreak, 0.0));
    					}
    				}
    			}

    			// taking care of right break
    			if (is_last_infinity())
    			{
    				if (rightBreak!=numeric_limits<double>::infinity())
    				{// something has to be cut on the right
    					std::pair<std::map<double, double>::iterator,bool> tmp_insert=Breakpoints_.insert(std::pair<double, double> (rightBreak, numeric_limits<double>::infinity()));
    					std::map<double, double>::iterator tmp_insert_it=tmp_insert.first;
    					if (!tmp_insert.second)
    					{
    						(tmp_insert_it)->second=numeric_limits<double>::infinity();
    					}

    					++tmp_insert_it;
    					if (tmp_insert_it!=Breakpoints_.end())
    					{
    						Breakpoints_.erase(tmp_insert_it,Breakpoints_.end());
    					}
    					if (Breakpoints_.size()==1)
    					{
    						FirstSlopeVal_=numeric_limits<double>::infinity();
    					}
    				}
    			}else
    			{
    				if (rightBreak!=numeric_limits<double>::infinity() && Breakpoints_.rbegin()->first>rightBreak )
    				{// something has to be cut on the right
    					std::pair<std::map<double, double>::iterator,bool> tmp_insert=Breakpoints_.insert(std::pair<double, double> (rightBreak, numeric_limits<double>::infinity()));

    					std::map<double, double>::iterator tmp_insert_it=tmp_insert.first;
    					if (!tmp_insert.second)
    					{
    						(tmp_insert_it)->second=numeric_limits<double>::infinity();
    					}

    					++tmp_insert_it;
    					if (tmp_insert_it!=Breakpoints_.end())
    					{
    						Breakpoints_.erase(tmp_insert_it,Breakpoints_.end());
    					}
    					if (Breakpoints_.size()==1)
    					{
    						FirstSlopeVal_=numeric_limits<double>::infinity();
    					}
    				}
    			}
    		}

    		double tmp=Breakpoints_.begin()->first-leftBreak;
    		if (tmp>0)
    		{
    			res[0]=tmp;
    		}
    		tmp=rightBreak-Breakpoints_.rbegin()->first;
    		if (tmp>0)
    		{
    			res[1]=tmp;
    		}
    		return(res);
			//std::cout << __FUNCTION__ << "("<<leftBreak<< ","<<rightBreak<<")"<<" out "<<std::endl;
    		//this->print();
    		//std::cout << "FirstSlopeVal_: "<<FirstSlopeVal_<<endl;
         }


    void cplfunction::Sumf(cplfunction & cplfunction1){
   // std::cout << __FUNCTION__ <<" in "<<std::endl;
   //	this->print();
   //	cplfunction1.print();
   	//  std::cout<<"cplfunction1.FirstSlopeVal_ : "<<cplfunction1.FirstSlopeVal_<<std::endl;

   	//domaine de définition de la fonction Sumf est l'intersection des deux domaines
   	  if (cplfunction1.is_last_infinity())
   	  {
   		(*this).Squeeze(cplfunction1.Breakpoints_.begin()->first,numeric_limits<double>::infinity());
   	  }else
   	  {
  	   	(*this).Squeeze(cplfunction1.Breakpoints_.begin()->first,cplfunction1.Breakpoints_.rbegin()->first);
   	  }

   	  if (	(cplfunction1.Breakpoints_.size()==1) ||
   			  (cplfunction1.Breakpoints_.size()==2 && !cplfunction1.is_last_infinity()))
   	  {// linear function ... only one slope
   		  FirstSlopeVal_=FirstSlopeVal_+cplfunction1.FirstSlopeVal_;
   		  FirstBreakVal_=FirstBreakVal_+cplfunction1.FirstBreakVal_;
   	  }else if (is_a_point())
   	  {

   	  }else
   	  {// at least two slopes with 2 breaks or more
		std::map<double,double>::const_iterator it=cplfunction1.Breakpoints_.begin();
		++it;
		std::map<double, double>::const_iterator itplus=it;
		it=cplfunction1.Breakpoints_.begin();

		(*this).AddSimple(cplfunction1.FirstSlopeVal_,
				itplus->second+cplfunction1.FirstSlopeVal_,
				cplfunction1.FirstBreakVal_,itplus->first);

		++itplus;++it;
		while (itplus!=cplfunction1.Breakpoints_.end())
		{// enter this loop if there are more than 2 slopes (3 breaks or more)...
			(*this).AddSimple(0.0,itplus->second,0.0,itplus->first);
			++itplus;++it;
		}
      }

   	//std::cout << __FUNCTION__ <<" out "<<endl;
   	//this->print();
    }

    bool cplfunction::StartsLargerthan(cplfunction & Cplfunction1){
      cplfunction tmp1(*this);
      cplfunction tmp2(Cplfunction1);

      //domaine de définition de la fonction max est l'intersection des deux domaines
      if (tmp2.is_last_infinity())
      {//si Cplfunction1 est définie jusqu'à l'infini
        tmp1.Squeeze(tmp2.Breakpoints_.begin()->first, std::numeric_limits<double>::infinity());
      }else
      {//si Cplfunction1 est infinie au dela de son dernier breakpoint
        tmp1.Squeeze(tmp2.Breakpoints_.begin()->first,tmp2.Breakpoints_.rbegin()->first);
      }

      if (tmp1.is_last_infinity())
      {//si Cplfunction1 est définie jusqu'à l'infini
        tmp2.Squeeze(tmp1.Breakpoints_.begin()->first, std::numeric_limits<double>::infinity());
      }else
      {//si Cplfunction1 est infinie au dela de son dernier breakpoint
        tmp2.Squeeze(tmp1.Breakpoints_.begin()->first,tmp1.Breakpoints_.rbegin()->first);
      }

      if (tmp1.Breakpoints_.begin()->first == -std::numeric_limits<double>::infinity()){
        return(tmp1.FirstSlopeVal_<=tmp2.FirstSlopeVal_);
      }else{
        return(tmp1.evalf(tmp1.Breakpoints_.begin()->first)>=tmp2.evalf(tmp2.Breakpoints_.begin()->first));
      }
    }


    cplfunction cplfunction::MaxfStartsLarger_(cplfunction & Cplfunction1){
      //std::cout << __FUNCTION__ <<" in "<<std::endl;
      //this->print();
      cplfunction f1(*this);
      cplfunction f2(Cplfunction1);


      //domaine de définition de la fonction max est l'intersection des deux domaines
      if (f2.is_last_infinity())
      {//si Cplfunction1 est définie jusqu'à l'infini
        f1.Squeeze(f2.Breakpoints_.begin()->first, std::numeric_limits<double>::infinity());
      }else
      {//si Cplfunction1 est infinie au dela de son dernier breakpoint
        f1.Squeeze(f2.Breakpoints_.begin()->first,f2.Breakpoints_.rbegin()->first);
      }

      if (f1.is_last_infinity())
      {//si Cplfunction1 est définie jusqu'à l'infini
        f2.Squeeze(f1.Breakpoints_.begin()->first, std::numeric_limits<double>::infinity());
      }else
      {//si Cplfunction1 est infinie au dela de son dernier breakpoint
        f2.Squeeze(f1.Breakpoints_.begin()->first,f1.Breakpoints_.rbegin()->first);
      }


      //res.Breakpoints_.insert(f2.Breakpoints_.begin(), f2.Breakpoints_.end());
      // cas très particulier où la fonction est réduite à un point

      std::map<double, double>::iterator it_f1 = f1.Breakpoints_.begin();
      std::map<double, double>::iterator it_f1_end = f1.Breakpoints_.end();
      std::map<double, double>::iterator it_f2 = f2.Breakpoints_.begin();
      std::map<double, double>::iterator it_f2_end = f2.Breakpoints_.end();
      std::map<double, double>::iterator itplus1_f1 = f1.Breakpoints_.begin(); ++itplus1_f1;
      std::map<double, double>::iterator itplus1_f2 = f2.Breakpoints_.begin(); ++itplus1_f2;

      // initialisation
      cplfunction res;
      res.Breakpoints_[it_f1->first]=it_f1->second;
      res.FirstBreakVal_=f1.FirstBreakVal_;
      res.FirstSlopeVal_=f1.FirstSlopeVal_;
      double curx,tmp,Deltaf1, Deltaf2;
      double f1DEcurx,f2DEcurx, f1DEprecx, f2DEprecx;
      double a,b,d,e;
      bool f1islargerthanf2=true;
      double precx=-std::numeric_limits<double>::infinity();


      std::vector<double> cur_f1(3),cur_f2(3),Next_f1(3),Next_f2(3);
      //cur_f1[0] BreakPoint, cur_f1[1] is Val, cur_f1[2] is Slope
      curx=it_f1->first;
      if (curx==-std::numeric_limits<double>::infinity()){
        // on a besoin de gérer le cas où it_f1->first -std::numeric_limits<double>::infinity()
        // dans tous les cas it_f2->first=it_f1->first
        if (f1.is_a_point())  return(f1);
        if (f1.is_an_infinite_line()){
          if (f2.is_an_infinite_line()) {
            // f1 and f2 are an infinite line
            if (f1.FirstSlopeVal_==f2.FirstSlopeVal_){
              return(res);
            }else{
              res.Breakpoints_.insert(std::pair<double,double>((f2.FirstBreakVal_-f1.FirstBreakVal_)/(f1.FirstSlopeVal_-f2.FirstSlopeVal_),f2.FirstSlopeVal_-f1.FirstSlopeVal_));
              return(res);
            }
          }else{
            // f1 is an infinite line but not f2
            tmp=(f2.FirstBreakVal_-f1.FirstBreakVal_)/(f1.FirstSlopeVal_-f2.FirstSlopeVal_);
            if (tmp<itplus1_f2->first){ // f2 est devenu plus grand que f1 avant curx
              res.Breakpoints_.insert(std::pair<double,double>(tmp,f2.FirstSlopeVal_-f1.FirstSlopeVal_));
              f1islargerthanf2=false;
              curx=tmp;
            }
          }
        }else{
          //f1 is not an infinite line
          if (f2.is_an_infinite_line()){
            //f2 is an infinite line

            tmp=(f1.FirstBreakVal_-f2.FirstBreakVal_)/(f2.FirstSlopeVal_-f1.FirstSlopeVal_);
            //std::cout << "curx : " << curx << " tmp : " << tmp << std::endl;
            if (tmp<itplus1_f1->first){ // f2 est devenu plus grand que f1 avant curx
              res.Breakpoints_.insert(std::pair<double,double>(tmp,f2.FirstSlopeVal_-f1.FirstSlopeVal_));
              f1islargerthanf2=false;
              curx=tmp;
            }
            // ++it_f1;++itplus1_f1;
          }else{  // none of f1 and f2 are an infinite line
            tmp=(f1.FirstBreakVal_-f2.FirstBreakVal_)/(f2.FirstSlopeVal_-f1.FirstSlopeVal_);
            if (tmp<curx){ // f2 est devenu plus grand que f1 avant curx
              res.Breakpoints_.insert(std::pair<double,double>(tmp,f2.FirstSlopeVal_-f1.FirstSlopeVal_));
              f1islargerthanf2=false;
              curx=tmp;
            }
          }
        }
      }
      //std::cout << "curx : " << curx << std::endl;

      // maintenant curx != -infini
      cur_f1[0]=curx;   cur_f2[0]=curx;
      cur_f1[1]=f1.evalf(curx);
      cur_f2[1]=f2.evalf(curx);
      cur_f1[2]=f1.evalDeltaf(curx);  //f1.FirstSlopeVal_+it_f1->second;
      cur_f2[2]=f2.evalDeltaf(curx); // f2.FirstSlopeVal_+it_f2->second;
      // Recherche d'un point de départ

      Next_f1=GetNextBreakVal(it_f1,itplus1_f1,it_f1_end,cur_f1,f1.FirstBreakVal_);
      Next_f2=GetNextBreakVal(it_f2,itplus1_f2,it_f2_end,cur_f2,f2.FirstBreakVal_);
      //std::cout << "Next_f1[0] : " << Next_f1[0] << " Next_f1[1] : " << Next_f1[1]<<" Next_f1[2] : " <<Next_f1[2]<< std::endl;
      //std::cout << "Next_f2[0] : " << Next_f2[0] << " Next_f2[1] : " << Next_f2[1]<<" Next_f2[2] : " <<Next_f2[2]<< std::endl;



      while (Next_f1[0]!=Next_f2[0] || (Next_f1[0]!=std::numeric_limits<double>::infinity() && itplus1_f1!=it_f1_end))
        {
        //Rcout<<"entering the loop"<<endl;
        // condition de sortie :
        // Soit Next_f1[0]==Next_f2[0] ET Next_f1[0]==std::numeric_limits<double>::infinity() droite infinie
        // Soit Next_f1[0]==Next_f2[0] ET itplus1_f1==it_f1_end point final
        //Il existe un prochain breakpoint, alors on avance en le définissant comme courant
        // avancer d'un cran jusqu'au plus proche breakpoint
        f1DEprecx=f1.evalf(curx); f2DEprecx=f2.evalf(curx);
        precx=curx;
        if (Next_f1[0]<=Next_f2[0]){
          ++itplus1_f1;++it_f1;
          cur_f1=Next_f1;
          curx=cur_f1[0];
        }else{
          ++itplus1_f2;++it_f2;
          cur_f2=Next_f2;
          curx=cur_f2[0];
        }
        f1DEcurx=f1.evalf(curx); f2DEcurx=f2.evalf(curx);
        //std::cout << "f1DEcurx " <<f1DEcurx<< " f2DEcurx "<<f2DEcurx<< std::endl;
        //std::cout << "f2DEprecx " <<f2DEprecx<< " f1DEprecx "<<f1DEprecx<< std::endl;
        if (f1DEcurx>=f2DEcurx){
          if (!f1islargerthanf2){
            //std::cout << "If (1-a) " << std::endl;
            // on vient de changer de situation f1DEcurx>=f2DEcurx et f1DEprecx<=f2DEprecx
            // il faut chercher le point d'intersection et l'ajouter
            // on utilise le théorème de Thalès pour trouver le x0 de croissement
            a= f2DEprecx - f1DEprecx ;
            b= f1DEcurx - f2DEcurx ;
            d= curx-precx ;
            res.Breakpoints_.insert(std::pair<double,double>(precx+a*d/(a+b),f1.evalDeltaf(curx)-f2.evalDeltaf(precx)));
            f1islargerthanf2=true;
          }
          //Rcout << "If (1-b) " << std::endl;
          tmp=f1.evalDeltaf(curx)-f1.evalDeltaf(precx);
          //Rcout<< "(2) curx : "<<curx<<" tmp : "<< tmp<< "precx : "<<precx<<std::endl;
          if (tmp>0) res.Breakpoints_.insert(std::pair<double,double>(curx,tmp));
        }else{// f2(curx)>f1(curx)
          if (f1islargerthanf2){
            //std::cout << "If (2-a) " << std::endl;
            // on vient de changer de situation f1DEcurx<f2DEcurx et f1DEprecx>f2DEprecx
            // il faut chercher le point d'intersection et l'ajouter
            // on utilise le théorème de Thalès pour trouver le x0 de croissement
            a= f1DEprecx - f2DEprecx ;
            b= f2DEcurx - f1DEcurx ;
            d= curx-precx ;
            res.Breakpoints_.insert(std::pair<double,double>(precx+a*d/(a+b),f2.evalDeltaf(curx)-f1.evalDeltaf(precx)));
            f1islargerthanf2=false;
          }
          //std::cout << "If (2-b) " << std::endl;
          tmp=f2.evalDeltaf(curx)-f2.evalDeltaf(precx);
          if (tmp>0)  res.Breakpoints_.insert(std::pair<double,double>(curx,tmp));
        }


        Next_f1=GetNextBreakVal(it_f1,itplus1_f1,it_f1_end,cur_f1,f1.FirstBreakVal_);
        Next_f2=GetNextBreakVal(it_f2,itplus1_f2,it_f2_end,cur_f2,f2.FirstBreakVal_);


        //return(0);
      }

      if (itplus1_f1==it_f1_end){
        // il reste alors à couvrir la droite infinie restante
        Deltaf1=f1.evalDeltaf(curx); Deltaf2=f2.evalDeltaf(curx);
        //std::cout<< "Deltaf1 : "<<Deltaf1<<"Deltaf2 : "<<Deltaf2<< std::endl;
        f1DEcurx=f1.evalf(curx); f2DEcurx=f2.evalf(curx);
        if (f1islargerthanf2 && Deltaf1<Deltaf2){
          precx=curx;
          curx=curx+(f2DEcurx-f1DEcurx)/(Deltaf1-Deltaf2);
          tmp=f2.evalDeltaf(curx)-f1.evalDeltaf(precx);
          //std::cout<<"tmp : "<<tmp<<" curx :"<<curx<<std::endl;
          if (tmp>0) res.Breakpoints_.insert(std::pair<double,double>(curx,tmp));
        }

        if (!f1islargerthanf2 && Deltaf1>Deltaf2){
          precx=curx;
          curx=curx+(f2DEcurx-f1DEcurx)/(Deltaf1-Deltaf2);
            tmp=f1.evalDeltaf(curx)-f2.evalDeltaf(precx);
          if (tmp>0) res.Breakpoints_.insert(std::pair<double,double>(curx,tmp));
        }
      }else{
        // il ne reste qu'un point final avec une pente infinie
        curx=Next_f2[0]; // == Next_f1[0]
        res.Breakpoints_.insert(std::pair<double,double>(curx,std::numeric_limits<double>::infinity()));
      }



      return(res);
     }

    cplfunction cplfunction::Maxf(cplfunction & Cplfunction1){
      //f1 : tmp and f2 : Cplfunction1
      //cplfunction tmp(*this);
      //Breakpoints_.clear();

      if (this->StartsLargerthan(Cplfunction1)){
        return(this->MaxfStartsLarger_(Cplfunction1));
      }else{
        return(Cplfunction1.MaxfStartsLarger_(*this));
      }

      //std::cout << __FUNCTION__ <<" out "<<std::endl;
      //this->print();
    }



    void cplfunction::Swap(double y)
    {
		//Rcout << __FUNCTION__ << " " << y << " in "<< endl;
		//Rcout << "y-0.6" << " " << y-0.6 << " in "<< endl;
		//this->print();
 	   cplfunction tmp(*this);
 	   Breakpoints_.clear();
 	   std::map<double,double>::reverse_iterator rit = tmp.Breakpoints_.rbegin();

 	   if (tmp.is_a_point())
 	   {
 		  Breakpoints_[y-rit->first]=0;
 		  FirstSlopeVal_=numeric_limits<double>::infinity();
 	   }else
 	   {
 	 	   std::map<double,double>::reverse_iterator ritplus1 = tmp.Breakpoints_.rbegin();
 	 	   ++ritplus1;
 			double LastSlopeVal=0.;
 			if (tmp.is_last_infinity())
 			{
 				Breakpoints_[-numeric_limits<double>::infinity()]=0;
 			}else
 			{
 				Breakpoints_[y-rit->first]=0; ++rit;++ritplus1;
 			}

 			while(ritplus1 != tmp.Breakpoints_.rend()){
				Breakpoints_[y-rit->first] = rit->second;
				LastSlopeVal=LastSlopeVal+(rit->second);
				++rit;++ritplus1;
 			}
 			if (rit->first!=-numeric_limits<double>::infinity())
 			{
 				Breakpoints_[y-rit->first] = numeric_limits<double>::infinity();
 			}

 			FirstSlopeVal_=-(FirstSlopeVal_+LastSlopeVal);
 	   }
 	 // Rcout << __FUNCTION__ << " " << y << " out "<< endl;
 	  // this->print();
 	   //return(*this);
    }

    void cplfunction::Legendre()
    {
        	//Rcout << __FUNCTION__<< " in" <<endl;
        	//this->print();
    	FirstBreakVal_=-FirstBreakVal_;
		if (is_a_point())
		{// a point get transformed into a line
			FirstSlopeVal_=(Breakpoints_.begin())->first;
			Breakpoints_.clear();
			Breakpoints_[-numeric_limits<double>::infinity()]=0;
		} else if (is_an_infinite_line())
		{
			double breakval=(Breakpoints_.begin())->second+FirstSlopeVal_;
			Breakpoints_.clear();
			Breakpoints_[breakval]=0;
			FirstSlopeVal_=numeric_limits<double>::infinity();
		} else if (Breakpoints_.size()==1)
		{// only one breakpoint not a line not a point : this is a half line
		 // gives a half line with 2 breaks
			double breakval=FirstSlopeVal_;
			FirstSlopeVal_=(Breakpoints_.begin())->first;
			Breakpoints_.clear();
			Breakpoints_[-numeric_limits<double>::infinity()]=0;
			Breakpoints_[breakval]=numeric_limits<double>::infinity();
		}else{// 2 breaks or more
			double breakval;
	    	std::map<double,double>::iterator it=Breakpoints_.begin();
	    	if (it->first==-numeric_limits<double>::infinity())
	    	{
	    		breakval=FirstSlopeVal_+it->second;
	    		++it;
	    		FirstSlopeVal_=FirstSlopeVal_+it->second;
	    		it->second=0;
	    		Breakpoints_.erase(Breakpoints_.begin());
	    	}else
	    	{
	    		breakval=-numeric_limits<double>::infinity();
	    	}

		    if (is_last_infinity())
		    {
	    		double lastbreak=flip_push_left(breakval);
	    		Breakpoints_[lastbreak]=numeric_limits<double>::infinity();
	    	}
	    	else
	    	{
	    		flip_push_left(breakval);
	    	}
		}

			//Rcout << __FUNCTION__<< " out" <<endl;
			//this->print();
    }


    inline double cplfunction::flip_push_left( double left_val)
    {
	// input Breakpoints
	// x0,		x1,				x2,				...,		xn
	// f(x_0)	f(x1)-f(x0),	f(x2)-f(x1),	...,		f(xn)-f(xn-1)
	//output Breakpoints
	// left_val, 	f(x_0),		f(x1),		...,		f(xn-1),
	// x0,			x1-x0,		x2-x1,		...,		xn-xn-1,
    //	returns f(xn)
    	std::map<double,double>::iterator itbegin = Breakpoints_.begin();
    	std::map<double,double>::iterator itend = Breakpoints_.end();
    	std::map<double, double>::iterator i=itbegin;
    	//assert(left_val< itbegin->first);


		double x = FirstSlopeVal_;
		double y = i->first;
    	FirstSlopeVal_=i->first;
		i->second = 0;
		const_cast<double&>(i->first) = left_val;
		double old_second=0;
		++i;

		while ( i != itend)
		{
			double old_second = i->second+x;
			double old_first=i->first;
			i->second = i->first-y;
			const_cast<double&>(i->first) = x;
			x = old_second;
			y = old_first;
			++i;
		}

		return x;

    }


  inline void cplfunction::shift_Breakpoint( std::map<double,double>::iterator it,double shift)
  {
	  // about changing the key of an element in a map
	  // http://stackoverflow.com/questions/5743545/what-is-the-fastest-way-to-change-a-key-of-an-element-inside-stdmap
    std::swap(Breakpoints_[it->first+shift], it->second);
    Breakpoints_.erase(it);
  }

 inline void cplfunction::shift_Breakpoint( std::map<double,double>::reverse_iterator it,double shift)
  {
    double first=it->first;
    double second=it->second;
    std::map<double,double>::reverse_iterator j = it ; ++j;
    Breakpoints_.erase(j.base()); // http://stackoverflow.com/questions/1830158/how-to-call-erase-with-a-reverse-iterator
    Breakpoints_[first+shift]=second;
  }


 // end of class cplfunction definition

//}// fin namespace Linear




////////////:
////////////


cplfunction Suml(cplfunction const & cplfunction_1,cplfunction const & cplfunction_2)
{
	 cplfunction tmp1(cplfunction_1);cplfunction  tmp2(cplfunction_2);
   		if (cplfunction_2.Breakpoints_.size()>cplfunction_1.Breakpoints_.size()){
   			tmp2.Sumf(tmp1);
   			return(tmp2);
   		}else{
   			tmp1.Sumf(tmp2);
   			return(tmp1);
   		}
}




cplfunction InfConv(cplfunction const & cplFunction_1,cplfunction const & cplFunction_2){
       cplfunction tmp1=cplFunction_1,tmp2=cplFunction_2;
    	 tmp1.Etoile();
    	 tmp2.Etoile();
    	 cplfunction res=Suml(tmp1,tmp2);
    	 res.Etoile();
    	 return(res);
     }
// static void finalizer_of_cplfunction( cplfunction* ptr ){
//    ptr->cplfunction::~cplfunction();
    //printf("finalizer has been called\n");
// }

cplfunction InfConfFunct(cplfunction const & cplFunction_1,cplfunction const & cplFunction_2,double y ){
       cplfunction tmp1(cplFunction_1),tmp2(cplFunction_2);
    	 tmp2.Swap(y);
    	 cplfunction B=Suml(tmp1,tmp2);
    	 return(B);
     }






  // Destructor
  cplfunctionvec::~cplfunctionvec(){
    MycplfunctionList_.clear();
  }

  //Constructors
  cplfunctionvec::cplfunctionvec() : MycplfunctionList_(){}
  cplfunctionvec::cplfunctionvec(int i) : MycplfunctionList_(i){}
  cplfunctionvec::cplfunctionvec(cplfunctionvec const & x):MycplfunctionList_(x.MycplfunctionList_){}

   cplfunctionvec::cplfunctionvec(std::vector<double> S1,std::vector<double> S2,std::vector<double> B1,std::vector<double> B2, std::vector<double> f0)
     : MycplfunctionList_()
     {
     	  int length=S1.size();
	  std::vector<double> Slopes(2);
	  std::vector<double> BreakPoints(2);
	  //std::vector<double> ZeroVal(1);
	  for (int compteur=0; compteur<length; compteur++){
		Slopes[0]=S1[compteur];Slopes[1]=S2[compteur];
		BreakPoints[0]=B1[compteur];BreakPoints[1]=B2[compteur];
		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
		MycplfunctionList_.push_back(cplfunction(Slopes,BreakPoints,f0[compteur]));
	  }
   }

  cplfunctionvec::cplfunctionvec(std::vector<double> S1, std::vector<double> B1, std::vector<double> f0)
  : MycplfunctionList_()
  {
	  int length=S1.size();
	  std::vector<double> Slopes(1);
	  std::vector<double> BreakPoints(1);
	 // std::vector<double> ZeroVal(1);
	  for (int compteur=0; compteur<length; compteur++){
		Slopes[0]=S1[compteur];
		BreakPoints[0]=B1[compteur];
		//ZeroVal[0]=f0[compteur];
		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
		MycplfunctionList_.push_back(cplfunction(Slopes,BreakPoints,f0[compteur]));
	  }
  }

  class nonsuitableinputVec : public std::exception {
      public:
      const char * what() { return "non suitable input, if it is a dict, should have Slopes  and BreakPoints fields"; }
    };
  //Wrapper to base functions
  std::vector<cplfunction>::iterator cplfunctionvec::begin(){return(MycplfunctionList_.begin());}
  std::vector<cplfunction>::iterator cplfunctionvec::end(){return(MycplfunctionList_.end());}
  std::vector<cplfunction>::reverse_iterator cplfunctionvec::rbegin(){return(MycplfunctionList_.rbegin());}
  void cplfunctionvec::vec_set( int i,cplfunction value ) { MycplfunctionList_.at(i) = value; }
  cplfunction cplfunctionvec::vec_get( int i) { return(MycplfunctionList_.at(i)); }
  int cplfunctionvec::size(){ return(MycplfunctionList_.size()); }

  void cplfunctionvec::push_back(cplfunction func){MycplfunctionList_.push_back(func);}

  void cplfunctionvec::Max_(cplfunctionvec func1,cplfunctionvec func2){
    //assert(MycplfunctionList_.size())==this->size());
    std::vector<cplfunction>::iterator it_1=func2.begin();
    std::vector<cplfunction>::iterator it_end=func2.end();
    std::vector<cplfunction>::iterator it_2=func1.begin();
    while(it_1!=it_end){
      MycplfunctionList_.push_back(it_1->Maxf(*it_2));
      ++it_1;++it_2;
    }
  }

  void cplfunctionvec::Maxf_1Breaks_withO(std::vector<double> S1,std::vector<double> B1, std::vector<double> f0){
    //assert(MycplfunctionList_.size())==this->size());
    cplfunctionvec func(S1,B1,f0);
    std::vector<cplfunction>::iterator it_1=MycplfunctionList_.begin();
    std::vector<cplfunction>::iterator it_end=MycplfunctionList_.end();
    std::vector<cplfunction>::iterator it_2=func.begin();
    while(it_1!=it_end){
      it_1->Maxf(*it_2);
      ++it_1;++it_2;
    }
  }

  void cplfunctionvec::Maxf_2Breaks_withO(std::vector<double> S1,std::vector<double> S2,std::vector<double> B1,std::vector<double> B2, std::vector<double> f0){
    //assert(MycplfunctionList_.size())==this->size());
    cplfunctionvec func(S1,S2,B1,B2,f0);
    std::vector<cplfunction>::iterator it_1=MycplfunctionList_.begin();
    std::vector<cplfunction>::iterator it_end=MycplfunctionList_.end();
    std::vector<cplfunction>::iterator it_2=func.begin();
    while(it_1!=it_end){
      it_1->Maxf(*it_2);
      ++it_1;++it_2;
    }
  }

  // serialized push



  void cplfunctionvec::SerialPush_1Breaks_Functions(std::vector<double> S1, std::vector<double> B1)
  {
	  int length=S1.size();
	  std::vector<double> Slopes(1);
	  std::vector<double> BreakPoints(1);
	  for (int compteur=0; compteur<length; compteur++){
		Slopes[0]=S1[compteur];
		BreakPoints[0]=B1[compteur];
		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
		MycplfunctionList_.push_back(cplfunction(Slopes,BreakPoints,0.));
	  }
  }

  std::vector<double> cplfunctionvec::Evalf(std::vector<double> x)
  {
    double compteur=0; // should test the size of x
    if (x.size()!=MycplfunctionList_.size()) std::cout<<"wrong size"<<std::endl;
    std::vector<double> fx(x.size());

 	std::vector<cplfunction>::iterator itend,it = MycplfunctionList_.begin();
 	itend=MycplfunctionList_.end();
 	while ( it!=itend)
 	{
		fx[compteur]=it->evalf(x[compteur]);
		++it; ++compteur;
	}
    return(fx);
  };

  std::vector<double> cplfunctionvec::EvalDeltaf2(std::vector<double> x)
  {
    int compteur=0; // should test the size of x
    if (x.size()!=MycplfunctionList_.size()) std::cout<<"wrong size"<<std::endl;
    std::vector<double> Deltafmoins(x.size());
    std::vector<double> Deltafplus(x.size());
    std::vector<double> tmp(2);
 	std::vector<cplfunction>::iterator itend,it = MycplfunctionList_.begin();
 	itend=MycplfunctionList_.end();
 	while ( it!=itend)
 	{
     	tmp=it->evalDeltaf2_(x[compteur]);
		Deltafmoins[compteur]=tmp[0];
		Deltafplus[compteur]=tmp[1];
		++it; ++compteur;
	}

        return(Deltafmoins);
  };

  std::vector<double> cplfunctionvec::EvalDeltafMoins(std::vector<double> x)
  {
    int compteur=0; // should test the size of x
    if (x.size()!=MycplfunctionList_.size()) std::cout<<"wrong size"<<std::endl;
    std::vector<double> Deltafmoins(x.size());
    std::vector<double> Deltafplus(x.size());
    std::vector<double> tmp(2);
 	std::vector<cplfunction>::iterator itend,it = MycplfunctionList_.begin();
 	itend=MycplfunctionList_.end();
 	while ( it!=itend)
 	{
     	tmp=it->evalDeltaf2_(x[compteur]);
		Deltafmoins[compteur]=tmp[0];
		Deltafplus[compteur]=tmp[1];
		++it; ++compteur;
	}

        return(Deltafmoins);
  };

  std::vector<double> cplfunctionvec::EvalDeltafPlus(std::vector<double> x)
  {
    int compteur=0; // should test the size of x
    if (x.size()!=MycplfunctionList_.size()) std::cout<<"wrong size"<<std::endl;
    std::vector<double> Deltafmoins(x.size());
    std::vector<double> Deltafplus(x.size());
    std::vector<double> tmp(2);
 	std::vector<cplfunction>::iterator itend,it = MycplfunctionList_.begin();
 	itend=MycplfunctionList_.end();
 	while ( it!=itend)
 	{
     	tmp=it->evalDeltaf2_(x[compteur]);
		Deltafmoins[compteur]=tmp[0];
		Deltafplus[compteur]=tmp[1];
		++it; ++compteur;
	}

        return(Deltafplus);
  };


  void cplfunctionvec::Etoile()
  {
 	std::vector<cplfunction>::iterator itend,it = MycplfunctionList_.begin();
 	itend=MycplfunctionList_.end();
 	while ( it!=itend)
 	{
     	it->Etoile();
		++it;
	}
  }

//  void Swap(boost::python::list x_list)
//  {
//    int compteur=0; // should test the size of x
//    if (boost::python::len(x_list)!=MycplfunctionList_.size()) std::cout<<"wrong size"<<std::endl;
// 	std::vector<cplfunction>::iterator itend,it = MycplfunctionList_.begin();
// 	itend=MycplfunctionList_.end();
// 	while ( it!=itend)
// 	{
//     	it->Swap(boost::python::extract<double>(x_list[compteur]));
//		++it; ++compteur;
//	}
//  };


  void cplfunctionvec::SerialPush_2Breaks_Functions_withO(std::vector<double> S1,std::vector<double> S2,
                                     std::vector<double> B1,std::vector<double> B2, std::vector<double> f0)
  {
	  int length=S1.size();
	  std::vector<double> Slopes(2);
	  std::vector<double> BreakPoints(2);
	  //std::vector<double> ZeroVal(1);
	  for (int compteur=0; compteur<length; compteur++){
		Slopes[0]=S1[compteur];Slopes[1]=S2[compteur];
		BreakPoints[0]=B1[compteur];BreakPoints[1]=B2[compteur];
		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
		MycplfunctionList_.push_back(cplfunction(Slopes,BreakPoints,f0[compteur]));
	  }
  }
  void cplfunctionvec::SerialPush_2Breaks_Functions(std::vector<double> S1,std::vector<double> S2, std::vector<double> B1,std::vector<double> B2)
  {
	  int length=S1.size();
	  std::vector<double> Slopes(2);
	  std::vector<double> BreakPoints(2);
	  for (int compteur=0; compteur<length; compteur++){
		Slopes[0]=S1[compteur];Slopes[1]=S2[compteur];
		BreakPoints[0]=B1[compteur];BreakPoints[1]=B2[compteur];
		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
		MycplfunctionList_.push_back(cplfunction(Slopes,BreakPoints,0.));
	  }
  }

  // void SerialPush_nBreaks_Functions(Rcpp::NumericMatrix S, Rcpp::NumericMatrix B)
  // {
	//   int length=S.nrow(),nbbreak=S.ncol();
	//   std::vector<double> Slopes(nbbreak);
	//   std::vector<double> BreakPoints(nbbreak);
	//   for (int compteur=0; compteur<length; compteur++){
	// 	Slopes=S(compteur,_);BreakPoints=B(compteur,_);
	// 	//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
	// 	MycplfunctionList_.push_back(cplfunction(Slopes,BreakPoints,0));
	//   }
  // };
  //
  //
  //
  // void PushSlopesandBreakfromOptim(Rcpp::NumericMatrix Pmax,Rcpp::NumericMatrix pi)
  // {
  //   int length=Pmax.nrow();
  //   Rcpp::List tmpres;
  //   std::vector<double> Slopes,Breaks;
  //   for (int compteur=0; compteur<length; compteur++){
  //     tmpres=GetSlopesandBreakfromOptim_(Pmax(compteur,_),pi(compteur,_));
  //     Slopes= as<NumericVector>(tmpres["Slopes"]);
  //       Breaks= as<NumericVector>(tmpres["Breaks"]);
  //     MycplfunctionList_.push_back(cplfunction(Slopes,Breaks,0));
  //   }
  // }





//   void SerialPush_nBreaks_Market_Functions(Rcpp::NumericMatrix Pmax, Rcpp::NumericMatrix pi,
// 		  std::vector<double> NetConso,cplfunctionvec cplfunctionvecToadd)
//   {
// 	  int length=Pmax.nrow(),nbbreak=S.ncol();
// 	  std::vector<double> Slopes(nbbreak);
// 	  std::vector<double> BreakPoints(nbbreak);
// 	  Rcpp::List tmpList;
// 	  for (int compteur=0; compteur<length; compteur++){
// 		tmpList=GetSlopesandBreakfromOptim(Pmax(compteur,_),pi(compteur,_));
// 		Slopes=tmpList["Slope"];BreakPoints=tmpList["Break"];
// 		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
// 		MycplfunctionList_.push_back(cplfunction(Slopes,BreakPoints,0));
// 		MycplfunctionList_.rbegin()->Swap(NetConso[compteur]);
// 		MycplfunctionList_.rbegin()->suml(cplfunctionvecToadd[compteur]);
// 	  }
//   };


  void cplfunctionvec::SerialPush_Store_Functions(double gamma1, double gamma2, std::vector<double> Pflex, std::vector<double> Prod)
  {// gamma1 rendement stockage, gamma2 rendement d?stockage
	  int length=Pflex.size();

	  std::vector<double> Slopes2(3);
	  std::vector<double> BreakPoints2(3);
	  std::vector<double> Slopes1(2);
	  std::vector<double> BreakPoints1(2);
	  std::vector<cplfunction> res2;
	  cplfunction tmp1,tmp2,res;
	  double Prodtmp;
	  for (int compteur=0; compteur<length; compteur++){
		  //calcul de f1=(x*1(x>0)+gamma2 * x*1(x>0)-P/gamma1)_+
		  // + f2=(x*1(x>0)+gamma2 * x*1(x>0)-P/gamma1)_+
		  Prodtmp=Prod[compteur]+Pflex[compteur];
		if (Prod[compteur]<0)
		{// deux breakpoints
			Slopes2[0]=0;Slopes2[1]=gamma2;Slopes2[2]=1;
			BreakPoints2[0]=-std::numeric_limits<double>::infinity();
			BreakPoints2[1]=Prod[compteur]/gamma1; BreakPoints2[2]=0;
			res2.push_back(cplfunction(Slopes2,BreakPoints2,0));

		}else
		{
			Slopes1[0]=0;Slopes1[1]=1;
				BreakPoints1[0]=-std::numeric_limits<double>::infinity();
				BreakPoints1[1]=Prod[compteur]/gamma1;
				res2.push_back(cplfunction(Slopes1,BreakPoints1,0));
		}


		if (Prodtmp<0)
		{// deux breakpoints
			Slopes2[0]=0;Slopes2[1]=gamma2;Slopes2[2]=1;
			BreakPoints2[0]=-std::numeric_limits<double>::infinity();
			BreakPoints2[1]=Prodtmp/gamma1; BreakPoints2[2]=0;
			//(MycplfunctionList_.rbegin())->Sumf(cplfunction(Slopes2,BreakPoints2,0))
			MycplfunctionList_.push_back(cplfunction(Slopes2,BreakPoints2,0));

		}else
		{
			Slopes1[0]=0;Slopes1[1]=1;
				BreakPoints1[0]=-std::numeric_limits<double>::infinity();
				BreakPoints1[1]=Prodtmp/gamma1;
				MycplfunctionList_.push_back(cplfunction(Slopes1,BreakPoints1,0));
		}

		// res=Suml(tmp1,tmp2);
    //cplfunction tmp3=Suml(tmp2,tmp1);
    //MycplfunctionList_.push_back(tmp1);
		(MycplfunctionList_.rbegin())->Sumf(*(res2.rbegin()));
		//itend=MycplfunctionList_.rbegin();
		//itend->Sumf(tmp2);
		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
		//MycplfunctionList_.push_back(Suml(tmp1,tmp2));
	  }
  }

  void cplfunctionvec::SerialPenalize(std::vector<double> alpha,
		  std::vector<double> inf,std::vector<double> sup)
  {
	  int length=MycplfunctionList_.size();
	  //assert(alpha.size()==length);
	  std::vector<double> Slopes(2);
	  std::vector<double> BreakPoints(2);
	  std::vector<cplfunction> f;
    double zero=0;
	  cplfunction tmp1,tmp2,tmp3;
	  for (int compteur=0; compteur<length; compteur++){
		Slopes[0]=alpha[compteur];Slopes[1]=std::numeric_limits<double>::infinity();
		BreakPoints[0]=inf[compteur];BreakPoints[1]=sup[compteur];
		tmp1=MycplfunctionList_[compteur];
		//f.push_back(cplfunction(Slopes,BreakPoints,zero));
		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
		vec_set(compteur,Suml(tmp1,cplfunction(Slopes,BreakPoints,zero)));
	  }
  }


  //Optim problem solving

std::vector<double> cplfunctionvec::OptimMargInt(std::vector<double>  Pmoins,
                std::vector<double>  Pplus,
                std::vector<double>  Cmoins,
                std::vector<double>  Cplus)
{

 	//cplfunctionvec Couts =*Coutsptr;
 	int length=Pmoins.size();
 	int compteur=0;
 	std::vector<double> xEtoile(length);
 	std::vector<cplfunction> f;
 	cplfunction tmpfunc,tmpfunc2,tmpfunc3;
 	std::vector<cplfunction>::iterator itend,it = MycplfunctionList_.begin();
 	tmpfunc2=*it;
 	tmpfunc2.Squeeze(Pmoins[compteur],Pplus[compteur]);
 	f.push_back(tmpfunc2);
 	compteur++; ++it;
 	itend=MycplfunctionList_.end();
 	while ( it!=itend)
 	{
		tmpfunc=*it;
		cplfunction tmpfunc2= *(f.rbegin());
		if (tmpfunc.is_an_infinite_line())
		{
 			tmpfunc2.Squeeze(Cmoins[compteur-1],Cplus[compteur-1]);
 			tmpfunc2.EpiSum_Withline(Pmoins[compteur],Pplus[compteur],tmpfunc.FirstSlopeVal_);
 			f.push_back(tmpfunc2);
		}else
		{
 			tmpfunc.Squeeze(Pmoins[compteur],Pplus[compteur]);
 			tmpfunc.Legendre();
 			tmpfunc2.Squeeze(Cmoins[compteur-1],Cplus[compteur-1]);
 			tmpfunc2.Legendre();
 			if (tmpfunc.Breakpoints_.size()>tmpfunc2.Breakpoints_.size())
 			{
				tmpfunc.Sumf(tmpfunc2);
				tmpfunc.Legendre();
				f.push_back(tmpfunc);
 			}else
 			{
				tmpfunc2.Sumf(tmpfunc);
				tmpfunc2.Legendre();
				f.push_back(tmpfunc2);
 			}
		}
		compteur++; ++it;
 	}

 	std::vector<cplfunction>::reverse_iterator itr,itf,itfrend;
 	itr = MycplfunctionList_.rbegin();
 	itf= f.rbegin();
 	compteur=length-1;
 	tmpfunc= *(itf);  ++itf;
 	tmpfunc.Squeeze(Cmoins[compteur],Cplus[compteur]);
 	xEtoile[compteur]=tmpfunc.Argmin();
 	double z=xEtoile[compteur];

 	itfrend=f.rend();
 	while(itf!= itfrend)
 	{
		--compteur;
		tmpfunc=*itr; ++itr;
		tmpfunc2=*itf; ++itf;
		tmpfunc.Squeeze(Pmoins[compteur+1],Pplus[compteur+1]);
		tmpfunc2.Squeeze(Cmoins[compteur],Cplus[compteur]);
		tmpfunc2.Swap(z);
		if (tmpfunc.Breakpoints_.size()>tmpfunc2.Breakpoints_.size())
		{
 			tmpfunc.Sumf(tmpfunc2);
 			xEtoile[compteur]=tmpfunc.Argmin();
		}else
		{
 			tmpfunc2.Sumf(tmpfunc);
 			xEtoile[compteur]=tmpfunc2.Argmin();
		}
		z=z-xEtoile[compteur];
		xEtoile[compteur]=z;
 	}
 	double tmpval,tmpval1=0;
 	for (int i=0;i<length;i++)
 	{
		tmpval=xEtoile[i];
		xEtoile[i]=xEtoile[i]-tmpval1;
		tmpval1=tmpval;
 	}

 	return(xEtoile);

 	};

	std::vector<double> cplfunctionvec::OptimPriceMarket_(std::vector<double> Pplus, double Conso)
	{
	int length=MycplfunctionList_.size();
	int compteur=0;
	std::vector<double> xEtoile(length);


	std::vector<cplfunction> f;
	cplfunction tmpfunc,tmpfunc2,tmpfunc3;
	std::vector<cplfunction>::iterator itend,it = MycplfunctionList_.begin();
	tmpfunc2=*it;
	//Rcout<<"Compteur1 "<<compteur<<endl;
	tmpfunc2.Squeeze((double)0.,Pplus[compteur]);
	f.push_back(tmpfunc2);
	compteur++; ++it;

	itend=MycplfunctionList_.end();
	while ( it!=itend)
	{
			tmpfunc=*it;
			cplfunction tmpfunc2= *(f.rbegin());
			if (tmpfunc.is_an_infinite_line())
			{
				tmpfunc2.EpiSum_Withline((double)0.,Pplus[compteur],tmpfunc.FirstSlopeVal_);
				f.push_back(tmpfunc2);
			}else
			{
				tmpfunc.Squeeze((double)0.,Pplus[compteur]);
				tmpfunc.Legendre();
				tmpfunc2.Legendre();
				if (tmpfunc.Breakpoints_.size()>tmpfunc2.Breakpoints_.size())
				{
					tmpfunc.Sumf(tmpfunc2);
					tmpfunc.Legendre();
					f.push_back(tmpfunc);
				}else
				{
					tmpfunc2.Sumf(tmpfunc);
					tmpfunc2.Legendre();
					f.push_back(tmpfunc2);
				}
			}

			compteur++; ++it;
  	}
  	std::vector<cplfunction>::reverse_iterator itr,itf,itfrend;
     itr = MycplfunctionList_.rbegin();
     itf= f.rbegin();
     compteur=length-1;
     tmpfunc= *(itf);  ++itf;
   	xEtoile[compteur]=Conso;
   	double z=xEtoile[compteur];
	itfrend=f.rend();
	while(itf!= itfrend)
	{
			--compteur;
			tmpfunc=*itr; ++itr;
			tmpfunc2=*itf; ++itf;
			tmpfunc.Squeeze((double)0.,Pplus[compteur+1]);
			tmpfunc2.Swap(z);
			if (tmpfunc.Breakpoints_.size()>tmpfunc2.Breakpoints_.size())
			{
				tmpfunc.Sumf(tmpfunc2);
				xEtoile[compteur]=tmpfunc.Argmin();
			}else
			{
				tmpfunc2.Sumf(tmpfunc);
				xEtoile[compteur]=tmpfunc2.Argmin();
			}
			z=z-xEtoile[compteur];
			xEtoile[compteur]=z;
		}

  	 double tmpval,tmpval1=0;
  	 for (int i=0;i<length;i++){
  		 tmpval=xEtoile[i];
  		 xEtoile[i]=xEtoile[i]-tmpval1;
  		 tmpval1=tmpval;
  	 }
     return(xEtoile);
	}

cplfunctionvec cplfunctionvec::Maxf(cplfunctionvec const & func1)
{
  //assert(MycplfunctionList_.size())==this->size());
  cplfunctionvec tmp1(func1);
  cplfunctionvec  tmp2(*this);
  cplfunctionvec res;
  std::vector<cplfunction>::iterator it_1=tmp2.begin();
  std::vector<cplfunction>::iterator it_end=tmp2.end();
  std::vector<cplfunction>::iterator it_2=tmp1.begin();
  while(it_1!=it_end){
    res.push_back(it_1->Maxf(*it_2));
    ++it_1;++it_2;
  }
  return(res);
}


//boost::python::list OptimPriceStorage_(boost::python::list Prices_List,boost::python::list  Pmoins_List,
//                boost::python::list  Pplus_List,
//                boost::python::list  Cmoins_List,
//                boost::python::list  Cplus_List){
//
//      	 // convert from python
//      	 std::vector<double> Prices;
//      	 py2stl(Prices_List,Prices);
//      	 std::vector<double> Pmoins;
//      	 py2stl(Pmoins_List,Pmoins);
//      	 std::vector<double> Pplus;
//      	 py2stl(Pplus_List,Pplus);
//      	 std::vector<double> Cmoins;
//      	 py2stl(Cmoins_List,Cmoins);
//      	 std::vector<double> Cplus;
//      	 py2stl(Cplus_List,Cplus);
//	int length=Pmoins.size();
//	std::vector<double>::iterator Pmoins_it=Pmoins.begin(), Cmoins_it=Cmoins.begin();
//	std::vector<double>::iterator Pplus_it=Pplus.begin(),Cplus_it=Cplus.begin();
//	std::vector<double>::iterator Prices_it=Prices.begin(),Prices_itend = Prices.end();
//	std::vector<cplfunction> f;
//	cplfunction tmpfunc,tmpfunc3;
//
//	cplfunction tmpfunc2=cplfunction(*Pmoins_it,*Pplus_it,0,*Prices_it,numeric_limits<double>::infinity());
//	f.push_back(tmpfunc2);
//	++Pplus_it;++Pmoins_it; ++Prices_it;
//
//	while ( Prices_it!=Prices_itend)
//	{
//		//cplfunction tmpfunc2= *(f.rbegin());
//		tmpfunc2.Squeeze(*Cmoins_it,*Cplus_it);
//		tmpfunc2.EpiSum_Withline(*Pmoins_it,*Pplus_it,*Prices_it);
//		f.push_back(tmpfunc2);
//		++Prices_it;++Cmoins_it;++Cplus_it;++Pmoins_it; ++Pplus_it;
//	}
//
//	std::vector<cplfunction>::reverse_iterator itr,itf,itfend;
//	itfend=f.rend();
//	itf= f.rbegin();
//	int compteur=length-1;
//	itf->Squeeze(Cmoins[compteur],Cplus[compteur]);
//	//tmpfunc.Squeeze(Pmoins[length-1],Pplus[length-1]);
//	std::vector<double> xEtoile(length);
//	xEtoile[compteur]=itf->Argmin();  ++itf;
//	double z=xEtoile[compteur];
//	while(itf!= itfend)
//	{
//		--compteur;
//		cplfunction tmpfunc=cplfunction(Pmoins[compteur+1],Pplus[compteur+1],0,Prices[compteur+1],numeric_limits<double>::infinity());
//		itf->Squeeze(Cmoins[compteur],Cplus[compteur]);
//		itf->Swap(z);
//		itf->Sumf(tmpfunc);
//		xEtoile[compteur]=itf->Argmin();
//		++itf;
//		z=z-xEtoile[compteur];
//		xEtoile[compteur]=z;
//	}
//	double tmpval,tmpval1=0;
//	for (int i=0;i<length;i++)
//	{
//		tmpval=xEtoile[i];
//		xEtoile[i]=xEtoile[i]-tmpval1;
//		tmpval1=tmpval;
//	}
//	boost::python::list res=stl2py(xEtoile);
// 	return(res);
//	}
//




double evalf_(std::vector<double> BreakPoints, std::vector<double> Prices,double x)
{
	if (abs(x)<=50)
	{
		return(x*Prices[1]);
	}else if (x>0)
	{
		return 50*Prices[1]+(x-50)*Prices[2];
	}else
	{
		return -50*Prices[1]+(x+50)*Prices[0];
	}
}


// std::vector<double> SerialOptimPriceStorage(NumericMatrix Prices,NumericMatrix BreakPoints,NumericVector Pmoins,NumericVector Pplus,NumericVector Cmoins,NumericVector Cplus)
// {
// 	  cplfunctionvec f(Prices.nrow());
// 	  for (int compteur=0; compteur<Prices.nrow(); compteur++){
// 		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
// 		  f.vec_set(compteur,cplfunction(Prices(compteur,_),BreakPoints(compteur,_),0));
// 	  }
//
// 	  int ncases=Pmoins.size();
// 	  NumericVector benefit(ncases);
//
// 	  for (int compteur=0; compteur<ncases; compteur++){
// 		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
// 		NumericVector Pmoinstmp(Prices.nrow(),Pmoins[compteur]);
// 		NumericVector Cmoinstmp(Prices.nrow(),Cmoins[compteur]);
// 		NumericVector PPlustmp(Prices.nrow(),Pplus[compteur]);
// 		NumericVector Cplustmp(Prices.nrow(),Cplus[compteur]);
// 		NumericVector res=(f.OptimMargInt(Pmoinstmp,PPlustmp,Cmoinstmp,Cplustmp))["xEtoile"];
// 		benefit[compteur]=0;
// 		 for (int compteur2=0; compteur2<Prices.nrow(); compteur2++){
// 				//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
// 			 benefit[compteur]=benefit[compteur]+evalf_(BreakPoints(compteur2,_),Prices(compteur2,_),res[compteur2]);
// 		}
// 	  }
// 	return benefit;
// }
//


//
//
//     boost::python::list OptimPriceMarket_l_(std::vector<double> Prices,std::vector<double> Pplus, double Conso)
//     {
//         std::vector<double> Prices;
//         std::vector<double> Pplus;
//
//         int length=Prices.size();
//     std::vector<double>::iterator Pplus_it=Pplus.begin(),Prices_itend = Prices.end(),Prices_it=Prices.begin();
//     std::vector<cplfunction> f;
//     cplfunction tmpfunc,tmpfunc3;
//
//     cplfunction tmpfunc2=cplfunction((double)0.,*Pplus_it,0,*Prices_it,numeric_limits<double>::infinity());
//     f.push_back(tmpfunc2);
//     ++Pplus_it; ++Prices_it;
//     while ( Prices_it!=Prices_itend)
//     {
//     	tmpfunc2.EpiSum_Withline((double)0.,*Pplus_it,*Prices_it);
//     	f.push_back(tmpfunc2);
//     	++Prices_it; ++Pplus_it;
//     }
//     std::vector<cplfunction>::reverse_iterator itr,itf,itfend;
//     itfend=f.rend();
//     itf= f.rbegin();
//     int compteur=length-1;
//     std::vector<double> xEtoile(length);
//     xEtoile[compteur]=Conso;  ++itf;
//     double z=xEtoile[compteur];
//     while(itf!= itfend)
//     {
//     	--compteur;
//     	cplfunction tmpfunc=cplfunction((double)0.,Pplus[compteur+1],0,Prices[compteur+1],numeric_limits<double>::infinity());
//     	itf->Swap(z);
//     	itf->Sumf(tmpfunc);
//     	xEtoile[compteur]=itf->Argmin();
//     	++itf;
//     	z=z-xEtoile[compteur];
//     	xEtoile[compteur]=z;
//     }
//     double tmpval,tmpval1=0;
//     for (int i=0;i<length;i++)
//     {
//     	tmpval=xEtoile[i];
//     	xEtoile[i]=xEtoile[i]-tmpval1;
//     	tmpval1=tmpval;
//     }
//
//      boost::python::list xEtoile_list=stl2py(xEtoile);
//      return(xEtoile);
//     };


// Rcpp::NumericMatrix OptimPriceMarket_l(NumericMatrix OffresPrix,NumericMatrix Availability,NumericVector Conso)
// {
// 	int nbpasTemps=OffresPrix.nrow(),nbProd=OffresPrix.ncol();
// 	Rcpp::NumericMatrix Power(nbpasTemps,nbProd);
//
// 	for (int compteurt=0; compteurt<nbpasTemps; compteurt++){
// 		  Power(compteurt,_)=OptimPriceMarket_l_(OffresPrix(compteurt,_),Availability(compteurt,_),Conso[compteurt]);
// 	}
//
// 	return Power;
// }

cplfunctionvec Sum(cplfunctionvec const & cplfunctionvec_1,cplfunctionvec const & cplfunctionvec_2)
{
	 cplfunctionvec res(cplfunctionvec_1);cplfunctionvec  tmp2(cplfunctionvec_2);
	 cplfunctionvec  tmp1(cplfunctionvec_1);
   	    if (tmp2.size()!=tmp1.size()) std::cout<<"wrong size"<<std::endl;
 	std::vector<cplfunction>::iterator itend,it_res = res.begin();
 	std::vector<cplfunction>::iterator it_tmp2=tmp2.begin();
 	std::vector<cplfunction>::iterator it_tmp1=tmp1.begin();
 	itend=res.end();
 	while ( it_res!=itend)
 	{
     	it_res->Sumf(*it_tmp2);
		++it_res; it_tmp1++; it_tmp2++;
	}
	return(res);
}

cplfunctionvec Max(cplfunctionvec const & func1,cplfunctionvec const & func2)
{
  //assert(MycplfunctionList_.size())==this->size());
  cplfunctionvec tmp1(func1);cplfunctionvec  tmp2(func2);
  cplfunctionvec res;
  std::vector<cplfunction>::iterator it_1=tmp2.begin();
  std::vector<cplfunction>::iterator it_end=tmp2.end();
  std::vector<cplfunction>::iterator it_2=tmp1.begin();
  while(it_1!=it_end){
    res.push_back(it_1->Maxf(*it_2));
    ++it_1;++it_2;
  }
  return(res);
}
//}





#endif /* cplfunction_CPP_ */
