/*  fabber_mc.h - Motion correction class

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "fabber_mc.h"

#ifdef __FABBER_MOTION

MCobj::MCobj(FabberRunData& allData, int dof)
{
	Matrix coords = allData.GetVoxelCoords();
	std::vector<int> size(3);
	size[0] = coords.Row(1).Maximum() + 1;
	size[1] = coords.Row(2).Maximum() + 1;
	size[2] = coords.Row(3).Maximum() + 1;
	mask = volume<float>(size[0], size[1], size[2]);
	mask = 0.0;

	//cout << size[0] << ", " << size[1] << ", " << size[2] << endl;
	//cout << coords.Ncols() << endl;
	for (int v=1; v<=coords.Ncols(); v++)
	{
		mask(coords(1, v), coords(2, v), coords(3, v)) = 1.0;
	//	cout << coords(1, v) << ", " << coords(2, v) << ", " <<  coords(3, v) << endl;
	}
	Matrix datamat = allData.GetMainVoxelData();
	wholeimage.setmatrix(datamat,mask);

	// the following sets up an initial zero deformation field
	userdof=dof;
	num_iter=10;
	modelpred=wholeimage;
	modelpred=0.0f;
	if (userdof > 12)
	{
		defx=modelpred;
		defx.setROIlimits(0,2);
		defx.activateROI();
		defx=defx.ROI();
		defy=defx;
		defz=defx;
		// Unnecessary initialisations?!?
		tmpx=defx;
		tmpy=defx;
		tmpz=defx;
	}
	else
	{
		mcf.setparams("verbose",false);
	}
	affmat=IdentityMatrix(4);
	finalimage=modelpred;
}

void MCobj::run_mc(const Matrix& modelpred_mat, Matrix& finalimage_mat)
{
	modelpred.setmatrix(modelpred_mat,mask);
	if (userdof>12)
	{
		UpdateDeformation(wholeimage,modelpred,num_iter,defx,defy,defz,finalimage,tmpx,tmpy,tmpz);
		defx=tmpx;
		defy=tmpy;
		defz=tmpz;
	}
	else
	{
		// mcf.register_volumes(4D reference,4D input image,refweight,inweight,4D output image);
		affmat = mcf.register_volumes(modelpred,wholeimage,mask,mask,finalimage);

		// apply transforms to wholeimage to get finalimage (the above is a dummy)
		for (int n=0; n<wholeimage.maxt(); n++)
		{
			affine_transform(wholeimage[n],finalimage[n],affmat.Rows(n*4+1, n*4+4));
		}
	}
	finalimage_mat = finalimage.matrix(mask);
}

#endif //__FABBER_MOTION
