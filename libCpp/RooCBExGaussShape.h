#ifndef ROO_CB_EX_GAUSS_SHAPE
#define ROO_CB_EX_GAUSS_SHAPE

#include "TROOT.h"
#include "TFitResult.h"
#include "TMath.h"

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"
// #include "RooAbsArg.h"
// #include "RooCategoryProxy.h"
// #include "RooAbsCategory.h"
// #include "Riostream.h"

class RooCBExGaussShape : public RooAbsPdf {
public:
  RooCBExGaussShape() {} ;
  RooCBExGaussShape(const char *name, const char *title, RooAbsReal& _m,
		    RooAbsReal& _m0,  RooAbsReal& _sigma,     RooAbsReal& _alpha,
		    RooAbsReal& _n,   RooAbsReal& _sigma_2,   RooAbsReal& _tailLeft
		    );
  RooCBExGaussShape(const RooCBExGaussShape& other, const char* name);
  inline virtual TObject* clone(const char* newname) const { return new RooCBExGaussShape(*this,newname);}
  inline ~RooCBExGaussShape(){}
  Double_t evaluate() const ;

  ClassDef(RooCBExGaussShape, 2)

protected:

  RooRealProxy m ;
  RooRealProxy m0 ;
  RooRealProxy sigma ;
  RooRealProxy alpha ;
  RooRealProxy n ;
  RooRealProxy sigma_2 ;
  RooRealProxy tailLeft ;

};

#endif
