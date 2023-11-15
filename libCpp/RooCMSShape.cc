/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: PhysicsTools/TagAndProbe/RooCMSShape
 *
 *
 * Authors:
 *   Nadia Adam, Princeton - neadam@princeton.edu
 *   Adam Hunt, Princeton  - ahunt@princeton.edu
 *   Kalanand Mishra, Fermilab - kalanand@fnal.gov
 *
 * Description:
 *   Defines a probability density function which has exponential decay
 *   distribution at high mass beyond the pole position (say, Z peak)
 *   but turns over (i.e., error function) at low mass due to threshold
 *   effect. We use this to model the background shape in Z->ll invariant
 *   mass.
 * History:
 *
 *
 *****************************************************************************/

#include "RooCMSShape.h"

ClassImp(RooCMSShape)

 RooCMSShape::RooCMSShape(const char *name, const char *title,
                        RooAbsReal& _x,
                        RooAbsReal& _alpha,
                        RooAbsReal& _beta,
                        RooAbsReal& _gamma,
                        RooAbsReal& _peak) :
   RooAbsPdf(name,title),
   x("x","x",this,_x),
   alpha("alpha","alpha",this,_alpha),
   beta("beta","beta",this,_beta),
   gamma("gamma","gamma",this,_gamma),
   peak("peak","peak",this,_peak)
 { }


 RooCMSShape::RooCMSShape(const RooCMSShape& other, const char* name):
   RooAbsPdf(other,name),
   x("x",this,other.x),
   alpha("alpha",this,other.alpha),
   beta("beta",this,other.beta),
   gamma("gamma",this,other.gamma),
   peak("peak",this,other.peak)
 { }



 Double_t RooCMSShape::evaluate() const
 {
  // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE

  //Double_t erf = TMath::Erfc((alpha - x) * beta);
  // Double_t erf = RooMath::erfc((alpha - x) * beta);
  Double_t erf = (1.0 + RooMath::erf((x - alpha) / beta));
  Double_t u = (x - peak)*gamma;

  if(u < -70) u = 1e20;
  else if( u>70 ) u = 0;
  else u = exp(-u);   //exponential decay
  return erf*u;
 }


 Int_t RooCMSShape::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
 {
   if (matchArgs(allVars,analVars,x)) return 1 ;
   return 0 ;
 }

 double RooCMSShape::analyticalIntegral(Int_t code, const char* rangeName) const
 {
   switch(code) {
     case 1:
     {

       Double_t max = x.max(rangeName);
       Double_t min = x.min(rangeName);

       Double_t inner = (-alpha/beta) + 0.5*beta*gamma;

       Double_t term0 = 1.0/gamma;

       Double_t term1 = exp( 0.25*gamma*( (-4.0*alpha) + (beta*beta*gamma) + (4.0*peak) ) );

       Double_t term2 = exp(gamma*(peak-min)) * RooMath::erfc((alpha-min)/beta);

       Double_t term3 = exp(gamma*(peak-max)) * RooMath::erfc((alpha-max)/beta);

       Double_t add1 = RooMath::erf(inner + (max/beta));

       Double_t add2 = RooMath::erf(inner + (min/beta));

       Double_t integral = term0 * ( (term1*(add1-add2)) + term2 - term3 );

       //std::cout << integral << std::endl;

       return integral ;
     }
   }
   assert(0) ;
   return 0 ;
 }

