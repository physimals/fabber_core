
#include "newimage/newimage.h"
#include "warpfns/warpfns.h"

#include <ostream>
#include <string>
#include <map>
#include <newmat.h>

void UpdateDeformation(const NEWIMAGE::volume4D<float>& wholeimage, const NEWIMAGE::volume4D<float>& modelpred, int no_iter,
		       const NEWIMAGE::volume4D<float>& prevdefx,
		       const NEWIMAGE::volume4D<float>& prevdefy, const NEWIMAGE::volume4D<float>& prevdefz,
			   NEWIMAGE::volume4D<float>& finalimage,
			   NEWIMAGE::volume4D<float>& defx, NEWIMAGE::volume4D<float> & defy,NEWIMAGE::volume4D<float> & defz);

