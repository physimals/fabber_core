#include <ostream>
#include <string>
#include <map>
#include <newmat.h>
#include "newimage/newimageall.h"
#include "warpfns/warpfns.h"

void UpdateDeformation(const volume4D<float>& wholeimage, const volume4D<float>& modelpred, int no_iter,
		       const volume4D<float>& prevdefx, 
		       const volume4D<float>& prevdefy, const volume4D<float>& prevdefz,
		       volume4D<float>& finalimage, 
		       volume4D<float>& defx, volume4D<float> & defy,volume4D<float> & defz);

