

#ifndef CPLFUNCTION_HPP
#define CPLFUNCTION_HPP

#include <iostream>
#include <math.h>
#include <limits>
#include <vector>
#include <map>


using namespace std;




//RCPP_EXPOSED_CLASS(cplfunction)
//RCPP_EXPOSED_CLASS(cplfunctionvec)
//RCPP_EXPOSED_CLASS(cpqfunction)
//RCPP_EXPOSED_CLASS(cpqfunctionvec)

std::vector<double> diffinv_(std::vector<double> x);
bool isincreasing(std::vector<double> arg);
double getSlope(std::pair<double,double> Coefficients,double val);
double getVal(std::pair<double,double> Coefficients,double val);
double getXetoile(std::pair<double,double> Coefficients);
std::pair<double,double> Slopes2Coeffs(double Slopes0,double Slopes1);
std::vector<double> GetNextBreakVal(std::map<double, double>::iterator it,
                                    std::map<double, double>::iterator itplus1,
                                    std::map<double, double>::iterator itend,
                                    std::vector<double> curval, //size 3
                                    double FirstBreakVal);

class cplfunction {
	// this class implements the convex continuous piecewise functions with a map
	// this allows nlog(n) sum of two such function
	// FirstSlopeVal_ is the absolute slope associated to the first breakpoint
	// Breakpoints_ is
	//
	//

    public:
    std::map<double,double> Breakpoints_; // Breakpoints_[x_0] vaut toujours 0
                                          // entre x_0 et x_1 la pente est FirstSlopeVal_
                                          // entre x_1 et x_2 la pente est FirstSlopeVal_+Breakpoints_[x_1]
                                          // ...
                                          // entre x_{n-1} et x_n la pente est FirstSlopeVal_+sum_{i=1}^{n-1}Breakpoints_[x_i]
                                          // après x_n la pente est FirstSlopeVal_+sum_{i=1}^{n}Breakpoints_[x_i]
                                          // si Breakpoints_[x_n] est infinie cela signifie que la fonction n'est pas définie après x_n
    double FirstBreakVal_; // firstbreakval : valeur en zero de l'extrapolation
                           // du premier morceau sauf si la fonction est restreinte à un point
                           // alors valeur en ce point
    double FirstSlopeVal_ ;
    // methods
    ~cplfunction();
    cplfunction();
    cplfunction(std::vector<double> Slopes, std::vector<double> BreakPoints,double FirstBreakVal);
    cplfunction(const cplfunction&);
    cplfunction* clone() const;
    cplfunction(double uniquebreak,double val);
    cplfunction(double uniquebreak,double val,double Slope1);
    cplfunction(double uniquebreak,double val,double Slope1, double Slope2);
    cplfunction(double breakleft,double breakright,double val,double Slope1,double Slope2);
    void print_();
    std::map<double,double> get_BreakPoints();
    void operator = (cplfunction s) ;
    void AddSimple(double leftslope, double rightslope, double val, double breakpoint);
    bool eq(cplfunction  const & cplfunction1);
    double evalf(double x);
    std::vector<double> evalDeltaf2_(double x);
    double evalDeltaf(double x);
    bool is_last_infinity();
    bool is_a_point();
    bool is_an_infinite_line();
    void Etoile();
    void EpiSum_Withline(double lowerbound,double upperbound,double Slope);
    double Argmin();
    void Squeeze(double leftBreak,double rightBreak);
    std::vector<double> Squeeze2(double leftBreak,double rightBreak);
    void Sumf(cplfunction & cplfunction1);
    bool StartsLargerthan(cplfunction & Cplfunction1);
    cplfunction MaxfStartsLarger_(cplfunction & Cplfunction1);
    cplfunction Maxf(cplfunction & Cplfunction1);
    void Swap(double y);
    void Legendre();
    inline double flip_push_left( double left_val);
    inline void shift_Breakpoint( std::map<double,double>::iterator it,double shift);
    inline void shift_Breakpoint( std::map<double,double>::reverse_iterator it,double shift);
};


cplfunction Suml(cplfunction const & cplfunction_1,cplfunction const & cplfunction_2);
cplfunction InfConv(cplfunction const & cplFunction_1,cplfunction const & cplFunction_2);
cplfunction InfConfFunct(cplfunction const & cplFunction_1,cplfunction const & cplFunction_2,double y );


class cplfunctionvec {

  public:
  std::vector<cplfunction> MycplfunctionList_;

   // methods
   ~cplfunctionvec();
   cplfunctionvec();
   cplfunctionvec(int i);
   cplfunctionvec(cplfunctionvec const & x);
   cplfunctionvec(std::vector<double> S1, std::vector<double> B1, std::vector<double> f0);
   cplfunctionvec(std::vector<double> S1,std::vector<double> S2,std::vector<double> B1,std::vector<double> B2, std::vector<double> f0);
   std::vector<cplfunction>::iterator begin();
   std::vector<cplfunction>::iterator end();
   std::vector<cplfunction>::reverse_iterator rbegin();
   void vec_set( int i,cplfunction value );
   cplfunction vec_get( int i);
   int size();
   void push_back(cplfunction func);
   void Max_(cplfunctionvec func1,cplfunctionvec func2);
   void Maxf_1Breaks_withO(std::vector<double> S1,std::vector<double> B1, std::vector<double> f0);
   void Maxf_2Breaks_withO(std::vector<double> S1,std::vector<double> S2,std::vector<double> B1,std::vector<double> B2, std::vector<double> f0);
   void SerialPush_1Breaks_Functions(std::vector<double> S1, std::vector<double> B1);
   void Etoile();
   void SerialPush_2Breaks_Functions_withO(std::vector<double> S1,std::vector<double> S2,std::vector<double> B1,std::vector<double> B2, std::vector<double> f0);
   void SerialPush_2Breaks_Functions(std::vector<double> S1,std::vector<double> S2, std::vector<double> B1,std::vector<double> B2);
   void SerialPush_Store_Functions(double gamma1, double gamma2, std::vector<double> Pflex, std::vector<double> Prod);
   void SerialPenalize(std::vector<double> alpha,std::vector<double> inf,std::vector<double> sup);
     std::vector<double> EvalDeltaf2(std::vector<double> x);
     std::vector<double> EvalDeltafMoins(std::vector<double> x);
     std::vector<double> EvalDeltafPlus(std::vector<double> x);
       std::vector<double> Evalf(std::vector<double> x);
   std::vector<double> OptimMargInt(std::vector<double>  Pmoins,
                std::vector<double>  Pplus,
                std::vector<double>  Cmoins,
                std::vector<double>  Cplus);
   std::vector<double> OptimPriceMarket_(std::vector<double> Pplus, double Conso);
   cplfunctionvec Maxf(cplfunctionvec const & func1);


};
#endif