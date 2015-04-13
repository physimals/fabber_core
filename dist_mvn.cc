/*  dist_mvn.cc - MultiVariate Normal distribution class/structure

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "dist_mvn.h"
#include "easyoptions.h"
#include "miscmaths/miscmaths.h"

using namespace NEWIMAGE;
using namespace Utilities;
using namespace MISCMATHS;

// Constructors

MVNDist::MVNDist()
{
  Tracer_Plus tr("MVNDist::MVNDist()");
  len = -1;
  precisionsValid = covarianceValid = false;
}

MVNDist::MVNDist(const MVNDist& from1, const MVNDist& from2)
{
  Tracer_Plus tr("MVNDist::MVNDist(from1,from2)");
  len = from1.len + from2.len;
  means = from1.means & from2.means;
  precisionsValid = false;
  
  // Always duplicate the covariances (even if this means some recalculation)
  // Otherwise if we use precisions.i(), zeros won't stay exactly zero
  covariance.ReSize(len);
  covariance = 0;
  covariance.SymSubMatrix(1, from1.len) = from1.GetCovariance();  
  covariance.SymSubMatrix(from1.len+1, from1.len+from2.len) = from2.GetCovariance();
  covarianceValid = true;

  assert(means.Nrows() == len);
}

const MVNDist& MVNDist::operator=(const MVNDist& from)
{
  // Not useful and dominates --debug-running-stack:
  // Tracer_Plus tr("MVNDist::operator=");

  assert(&from != NULL); // yes, this can happen.  References are but pointers in disguise...

  // Special case: assignment to self (is a no-op)
  if (&from == this)
    return *this;

  // Special case: assigned from an uninitialized MVNDist
  if (from.len == -1)
    {
      len = -1;
      precisionsValid = covarianceValid = false;
      // Note, might still be consuming large amounts of memory, even though
      // precisions & covariance are now inaccessible from the outside
      return *this;
    }

  assert(from.len == from.means.Nrows());
  SetSize(from.len);
  //  len = from.len;
  means = from.means;
  precisionsValid = from.precisionsValid;
  covarianceValid = from.covarianceValid;
  if (precisionsValid)
    precisions = from.precisions;
  //  else if (precisions.Nrows() != len)
  //    precisions.ReSize(len);

  if (covarianceValid)
    covariance = from.covariance;
  //  else if (covariance.Nrows() != len)
  //    covariance.ReSize(len);

  assert(means.Nrows() == len);

  return *this;
}

void MVNDist::CopyFromSubmatrix(const MVNDist& from, int first, int last, 
    bool checkIndependence)
{
    Tracer_Plus tr("MVNDist::CopyFromSubmatrix");
    len = last-first+1;
    means = from.means.Rows(first, last);
    precisionsValid = from.precisionsValid;
    covarianceValid = from.covarianceValid;
    if (precisionsValid)
        precisions = from.precisions.SymSubMatrix(first, last);
    else if (precisions.Nrows() != len)
        precisions.ReSize(len);
        
    if (covarianceValid)
        covariance = from.covariance.SymSubMatrix(first, last);
    else if (covariance.Nrows() != len)
        covariance.ReSize(len);
            
    assert(means.Nrows() == len);
    
    if (checkIndependence)
    {
        Matrix deps1 = from.GetCovariance().Rows(first, last).Columns(1, first-1);
        Matrix deps2 = from.GetCovariance().Rows(first, last).Columns(last+1, from.covariance.Ncols());
        if (!deps1.IsZero() || !deps2.IsZero() )
            throw Invalid_option("Covariance found in part of MVN that should be independent from the rest!");
    } 
    return;
}        

void MVNDist::SetSize(int dim)
{
  // Not useful and dominates --debug-running-stack:
  // Tracer_Plus tr("MVNDist::SetSize");
  if (dim<=0)
    throw RBD_COMMON::Logic_error("Can't have dim<=0\n");
    
  assert(means.Nrows() == len || len<0);    
    
  if (len != dim)
  {
    //Tracer_Plus tr("MVNDist::SetSize (actually resizing)");
    len = dim;
    means.ReSize(dim);
    precisions.ReSize(dim);
    covariance.ReSize(dim);
  }
  precisionsValid = false;
  covarianceValid = false;
  // means is also now undefined (or at least out-of-date)  

  assert(means.Nrows() == len);
  assert(precisions.Nrows() == len);
  assert(covariance.Nrows() == len);
}

// Accessors
const SymmetricMatrix& MVNDist::GetPrecisions() const
{
  Tracer_Plus tr("MVNDist::GetPrecisions");
  if (len == -1) throw Logic_error("MVN is uninitialized!\n");
  assert(means.Nrows() == len);
  if (!precisionsValid)
    {
      Tracer_Plus tr("MVNDist::GetPrecisions calculation");        
      assert(covarianceValid);
      // precisions and precisionsValid are mutable, 
      // so we can change them even in a const function
      precisions = covariance.i();
      precisionsValid = true;
    }
  assert(means.Nrows() == len);
  assert(precisions.Nrows() == len);
  return precisions;
}

const SymmetricMatrix& MVNDist::GetCovariance() const
{
  Tracer_Plus tr("MVNDist::GetCovariance");    
  if (len == -1) throw Logic_error("MVN is uninitialized!\n");
  assert(means.Nrows() == len);
  if (!covarianceValid)
    {
      Tracer_Plus tr("MVNDist::GetCovariance calculation");
      assert(precisionsValid);
      // covariance and covarianceValid are mutable, 
      // so we can change them even in a const function
      try 
        {
          covariance = precisions.i();
        } 
      catch (Exception)
        {
	  //          LOG_ERR("Oh dear, it didn't like that\nPrecisions = \n");
	  //          LOG_ERR(precisions);
	  //          LOG_ERR("Rethrowing...\n");
	  //          throw;
	  // Better behaviour: adds a tiny amount to the diagonal and tries again
	  Warning::IssueOnce("MVN precision (len==" + stringify(len) + ") was singular, adding 1e-10 to diagonal");
	  covariance = (precisions + IdentityMatrix(len)*1e-10).i();
        }
      covarianceValid = true;
    }
  assert(means.Nrows() == len);
  assert(covariance.Nrows() == len);
  return covariance;
}

void MVNDist::SetPrecisions(const SymmetricMatrix& from)
{
  Tracer_Plus tr("MVNDist::SetPrecisions");
  assert(from.Nrows() == len);
  assert(means.Nrows() == len);
  precisions = from;
  precisionsValid = true;
  covarianceValid = false;
  assert(means.Nrows() == len);
}

void MVNDist::SetCovariance(const SymmetricMatrix& from)
{
  Tracer_Plus tr("MVNDist::SetCovariance");
  //cout << from.Nrows() << " ---- " << len << endl;  
  assert(from.Nrows() == len);
  assert(means.Nrows() == len);
  covariance = from;
  covarianceValid = true;
  precisionsValid = false;
  assert(means.Nrows() == len);
}

void MVNDist::DumpTo(ostream& out, const string indent) const
{ 
  Tracer_Plus tr("MVNDist::Dump");
  out << indent << "MVNDist, with len == " << len 
       << ", precisionsValid == " << precisionsValid
       << ", covarianceValid == " << covarianceValid << endl;
  out << indent << "  Means: " << means.t();
  if (precisionsValid || covarianceValid)
    {
      out << indent << "  Covariance matrix:" << endl;
      for (int i=1; i<=len; i++)
	out << indent << "  " << GetCovariance().Row(i);
    }
  else
    out << indent << "  Covariance undefined." << endl;

  assert(means.Nrows() == len);
}

void MVNDist::Load(const string& filename)
{
    cout << "Reading MVN from file '" << filename << "'...\n";
    Matrix mat = read_vest(filename);
    
    cout << "Converting to an MVN...\n";
    // Format: [covariance means(:); means(:)' 1.0]
    const int N = mat.Nrows() - 1;
    
    if (N < 1 || mat != mat.t() || mat(N+1,N+1) != 1.0)
    {
        cout << "N == " << N << ", matrix:\n" << mat;
        throw Invalid_option("Inputted MVNs must be symmetric matrices!\nFormat = [covariance means(:); means(:)' 1.0]\n");
    }
    SetSize(N);
    means = mat.Column(len+1).Rows(1, N);
    SymmetricMatrix sym; 
    sym << mat.SubMatrix(1,N,1,N);
    SetCovariance(sym); 
    
    assert(means.Nrows() == len);
}

void MVNDist::Load(vector<MVNDist*>& mvns, const string& filename, const volume<float>& mask)
{
    Tracer_Plus tr("MVNDist::Load (static)");
    
    Matrix vols; 
    LOG_ERR("Reading MVNs from " << filename << endl);
 
    volume4D<float> input;
    read_volume4D(input,filename);
    vols=input.matrix(mask);
    
    const int nVoxels = vols.Ncols();
    for (unsigned i = 0; i < mvns.size(); i++) 
        assert(mvns[i] == NULL); // should've deleted everything first. 
        
    mvns.resize(nVoxels, NULL);
    const int nParams = ((int)sqrt(8*vols.Nrows()+1)-3)/2;
    assert( vols.Nrows() == nParams*(nParams+1)/2 + nParams+1 );
    
    SymmetricMatrix tmp(nParams);
    
//        cout << "--------\n";
//        cout << nParams << endl;
//        cout << nVoxels << endl;
//        cout << vols.Nrows() << endl;
//        cout << tmp.Nrows() << endl;
  
    assert(nVoxels > 0);
  
    for (int vox = 1; vox <= nVoxels; vox++)
    {
        assert(mvns.at(vox-1) == NULL);
        mvns[vox-1] = new MVNDist(nParams);
        //tmp << vols.Column(vox).Rows(1,nParams*(nParams+1)/2); // Doesn't work
	int index = 0;
        for (int r = 1; r <= nParams; r++)
	  for (int c = 1; c <= r; c++)
	    tmp(r,c) = vols(++index,vox);
	assert(index == nParams*(nParams+1)/2);
        mvns[vox-1]->SetCovariance(tmp);
        mvns[vox-1]->means = vols.Column(vox).Rows(nParams*(nParams+1)/2 + 1, nParams*(nParams+1)/2 + nParams);
        assert(vols(vols.Nrows(), vox) == 1);
        assert(mvns[vox-1]->means.Nrows() == mvns[vox-1]->len);
    } 
    
    assert(mvns[0] != NULL);

}

void MVNDist::Save(const vector<MVNDist*>& mvns, const string& filename, const volume<float>& mask)
{    
    Tracer_Plus tr("MVNDist::Save");
     
    // Save the MVNs in a NIFTI file as a single NIFTI_INTENT_SYMMATRIX 
    // last row/col is the means (1 in the corner).
    // Note that I'm using the 4th dim and should really be using the 5th,
    // according to the specification -- but I don't think it really matters.

    Matrix vols;

    const int nVoxels = mvns.size();
    assert(nVoxels > 0 && mvns.at(0) != NULL);   
    const int nParams = mvns.at(0)->means.Nrows();
    
//    means.ReSize(nParams, nVoxels);
//    precisions.ReSize(nParams*(nParams+1)/2, nVoxels);
    vols.ReSize(nParams*(nParams+1)/2 + nParams+1, nVoxels);
    
    ColumnVector aOne(1); aOne = 1.0;
    
    for (int vox = 1; vox <= nVoxels; vox++)
      {
        vols.Column(vox) = 
            mvns.at(vox-1)->GetCovariance().AsColumn()
            & mvns.at(vox-1)->means & aOne;
        // AsColumn uses row ordering on the lower triangular part, 
        // as NIFTI_INTENT_SYMMATRIX specifies: (1,1) (2,1) (2,2) (3,1)...
      }
    // Write the file
    volume4D<float> output(mask.xsize(),mask.ysize(),mask.zsize(),vols.Nrows());
    output.setmatrix(vols,mask);
    output.set_intent(NIFTI_INTENT_SYMMATRIX,0,0,0);
    output.setDisplayMaximumMinimum(output.max(),output.min());
    save_volume4D(output,filename);
}
