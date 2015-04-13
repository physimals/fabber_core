/*  dist_mvn.h - MultiVariate Normal distribution class/structure

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include <stdexcept>
#include "assert.h"
#include "easylog.h"
#include "newimage/newimageall.h"

using namespace NEWMAT;
using namespace MISCMATHS;

class MVNDist {
public:

  // Constructors:
  MVNDist(); // unknown size -- will be fixed by first SetPrecisions/SetCovariance
  MVNDist(int dim) { len = -1; SetSize(dim); } // known size
  MVNDist(const MVNDist& from) { len = -1; *this = from; } // copy constructor
  MVNDist(const MVNDist& from1, const MVNDist& from2); // concat constructor
  MVNDist(const string filename)
    { len = -1; Load(filename); }

  void CopyFromSubmatrix(const MVNDist& from, int first, int last, 
    bool checkIndependence = true);
  MVNDist GetSubmatrix(int first, int last, bool checkIndependence = true)
    { MVNDist ret; 
      ret.CopyFromSubmatrix(*this, first, last, checkIndependence); 
      return ret; }  

  const MVNDist& operator=(const MVNDist& from);
  void SetSize(int dim);
  int GetSize() const
    { assert(len == means.Nrows() || len<0); return len; }
    

  ColumnVector means; // You shouldn't ReSize this manually -- use SetSize instead

  const SymmetricMatrix& GetPrecisions() const;
  const SymmetricMatrix& GetCovariance() const;
  void SetPrecisions(const SymmetricMatrix& from);
  void SetCovariance(const SymmetricMatrix& from);

  void Dump(const string indent = "") const { DumpTo(LOG, indent); }
  void DumpTo(ostream& out, const string indent = "") const;

  void Load(const string& filename);
  static void Load(vector<MVNDist*>& mvns, const string& filename, const NEWIMAGE::volume<float>& mask);
  static void Save(const vector<MVNDist*>& mvns, const string& filename, const NEWIMAGE::volume<float>& mask);
  
 protected:
  int len; // should only be changed explicitly

 private:
  // Mutable, because they're a cache calculated on first use --
  // to the outside world, changes here don't affect const-ness.
  mutable SymmetricMatrix precisions;
  mutable SymmetricMatrix covariance;
  mutable bool precisionsValid;
  mutable bool covarianceValid;
  // Note that you shouldn't store the references from GetPrecisions/GetCovariance
  // to use later, because they may be out of date if a Set function has been 
  // called since.  That kinda violates const-ness.. sorry. 
};

inline ostream& operator<<(ostream& out, const MVNDist& dist)
{ dist.DumpTo(out); return out; }

