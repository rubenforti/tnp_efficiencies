#include <memory>

#include "TFitResult.h"

#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"



// Double CB with RooFit
// copied from https://github.com/gdujany/chibAnalysis/blob/master/My_double_CB/My_double_CB.h
//
class My_double_CB : public RooAbsPdf {
 public:
  My_double_CB() {} ;
  My_double_CB(const char *name, const char *title,
	       RooAbsReal& _x,
	       RooAbsReal& _mu,
	       RooAbsReal& _sig,
	       RooAbsReal& _a1,
	       RooAbsReal& _n1,
	       RooAbsReal& _a2,
	       RooAbsReal& _n2);
  My_double_CB(const My_double_CB& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new My_double_CB(*this,newname); }
  inline virtual ~My_double_CB() { }

 protected:

  RooRealProxy x ;
  RooRealProxy mu ;
  RooRealProxy sig ;
  RooRealProxy a1 ;
  RooRealProxy n1 ;
  RooRealProxy a2 ;
  RooRealProxy n2 ;

  Double_t evaluate() const ;

};
