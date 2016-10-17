/*  dist_mvn.h - MultiVariate Normal distribution class/structure

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include "dataset.h"
#include "easylog.h"

/**
 * Multivariate Normal Distribution
 *
 * For a given number of parameters, the distribution
 * is defined by:
 *    - The means, one for each parameter
 *    - The covariances, or equivalently their precisions,
 *      one for each distinct pair of parameters and
 *      irrespective of order
 */
class MVNDist
{
public:

	/**
	 * Default constructor
	 *
	 * Size will be undefined -- will be fixed by first SetPrecisions/SetCovariance
	 */
	MVNDist();

	/**
	 * Create distribution of known size
	 */
	MVNDist(int dim)
	{
		m_size = -1;
		SetSize(dim);
	}

	/**
	 * Copy constructor
	 */
	MVNDist(const MVNDist& from)
	{
		m_size = -1;
		*this = from;
	}
	/**
	 * Concatentate two MVNs
	 */
	MVNDist(const MVNDist& from1, const MVNDist& from2);

	/**
	 * Create from NIFTI file
	 */
	MVNDist(const string filename)
	{
		m_size = -1;
		Load(filename);
	}

	/**
	 * Copy using a subset of another MVN distribution's parameters
	 */
	void CopyFromSubmatrix(const MVNDist& from, int first, int last, bool checkIndependence = true);

	/**
	 * Get a subset of this MVN distribution as another MVN distribution
	 */
	MVNDist GetSubmatrix(int first, int last, bool checkIndependence = true)
	{
		MVNDist ret;
		ret.CopyFromSubmatrix(*this, first, last, checkIndependence);
		return ret;
	}

	/**
	 * Assignment operator
	 */
	const MVNDist& operator=(const MVNDist& from);

	/**
	 * Set the size
	 *
	 * Will resize means and covariances/precisions to match
	 */
	void SetSize(int dim);

	/**
	 * Get the size (number of parameters)
	 *
	 * @return size or -1 if not initialized
	 */
	int GetSize() const
	{
		assert(m_size == means.Nrows() || m_size<0);
		return m_size;
	}

	/**
	 * Mean values of each parameter
	 *
	 * FIXME this is a public member which shouldn't be
	 * resized - should encapsulate it to protect
	 * against misuse
	 */
	NEWMAT::ColumnVector means;

	/**
	 * Get the precisions
	 *
	 * This is a symmetric matrix, i.e. the precision of two
	 * parameters is unaffected by what order you specify
	 * them in.
	 *
	 * Precision is the inverse of covariance, just
	 * as in 1-dimension the precision of a random
	 * variable is the reciprocol of the variance
	 */
	const NEWMAT::SymmetricMatrix& GetPrecisions() const;

	/**
	 * Get the covariances
	 *
	 * @see GetPrecisions
	 */
	const NEWMAT::SymmetricMatrix& GetCovariance() const;

	/**
	 * Set the precisions
	 *
	 * Also checks that size of matrix
	 * is correct for current number of parameters
	 *
	 * Covariances will be updated lazily on next
	 * call to GetCovariances
	 */
	void SetPrecisions(const NEWMAT::SymmetricMatrix& from);

	/**
	 * Set the covariances
	 *
	 * Also checks that size of matrix
	 * is correct for current number of parameters
	 *
	 * Precisions will be updated lazily on next
	 * call to GetPrecisions
	 */
	void SetCovariance(const NEWMAT::SymmetricMatrix& from);

	/**
	 * Dump info to the default log
	 */
	void Dump(const string indent = "") const
	{
		DumpTo(LOG, indent);
	}

	/**
	 * Dump info to the specified output stream
	 */
	void DumpTo(ostream& out, const string indent = "") const;

	/**
	 * Load from VEST file
	 *
	 * The file must contain an N+1 x N+1 symmetric matrix.
	 * This contains the covariances and the means in the
	 * following example form (C=covariance, M=mean)
	 *
	 * CCCM
	 * CCCM
	 * CCCM
	 * MMM1
	 *
	 * Precisions will be updated lazily on next
	 * call to GetPrecisions
	 */
	void Load(const string& filename);

	/**
	 * Load a per-voxel vector of MVN distributions from a NIFTI file
	 *
	 * @param mvns One MVN for each voxel
	 */
	static void Load(vector<MVNDist*>& mvns, const string& filename, FabberRunData &data);

	/**
	 * Save a per-voxel vector of MVN distributions to a NIFTI file
	 *
	 * The file will contain an N+1xN+1 symmetric matrix
	 * as described in /ref Load, however the use
	 * of NIFTI_INTENT_SYMMATRIX means that only the
	 * lower triangle data needs to be specified
	 *
	 * @param mvns One MVN for each voxel
	 */
	static void Save(const vector<MVNDist*>& mvns, const string& filename, FabberRunData &data);

protected:
	int m_size; // should only be changed explicitly

private:
	// Mutable, because they're a cache calculated on first use --
	// to the outside world, changes here don't affect const-ness.
	mutable NEWMAT::SymmetricMatrix precisions;
	mutable NEWMAT::SymmetricMatrix covariance;
	mutable bool precisionsValid;
	mutable bool covarianceValid;
	// Note that you shouldn't store the references from GetPrecisions/GetCovariance
	// to use later, because they may be out of date if a Set function has been
	// called since.  That kinda violates const-ness.. sorry.
};

inline ostream& operator<<(ostream& out, const MVNDist& dist)
{
	dist.DumpTo(out);
	return out;
}

