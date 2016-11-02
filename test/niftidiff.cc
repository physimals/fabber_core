/**
 * A simple NIFTII diff program
 */

#define EPSILON 0.01

#include <newimage/newimageall.h>
#include <newmat.h>

#include <string>

using namespace std;

int main(int argc, char **argv)
{
	if (argc != 3)
	{
		cout << "Usage: niftidiff <file1> <file2>" << endl;
		return 1;
	}
	string f1 = argv[1];
	string f2 = argv[2];

	NEWIMAGE::volume4D<float> d1;
	NEWIMAGE::volume4D<float> d2;

	try
	{
		read_volume4D(d1, f1);

	} catch (...)
	{
		cerr << "Error reading file " << f1 << endl;
		return 1;
	}
	try
	{
		read_volume4D(d2, f2);
	} catch (...)
	{
		cerr << "Error reading file " << f2 << endl;
		return 1;
	}

	if (d1.xsize() != d2.xsize())
	{
		cerr << "Dimension mismatch: (x): " << d1.xsize() << "!=" << d2.xsize() << endl;
		return 1;
	}
	if (d1.ysize() != d2.ysize())
	{
		cerr << "Dimension mismatch: (y): " << d1.ysize() << "!=" << d2.ysize() << endl;
		return 1;
	}
	if (d1.zsize() != d2.zsize())
	{
		cerr << "Dimension mismatch: (z): " << d1.zsize() << "!=" << d2.zsize() << endl;
		return 1;
	}
	if (d1.tsize() != d2.tsize())
	{
		cerr << "Dimension mismatch: (t): " << d1.tsize() << "!=" << d2.tsize() << endl;
		return 1;
	}

	for (int x = 0; x < d1.xsize(); x++)
	{
		for (int y = 0; y < d1.ysize(); y++)
		{
			for (int z = 0; z < d1.zsize(); z++)
			{
				for (int t = 0; t < d1.tsize(); t++)
				{
					if (abs(d1(x, y, z, t) - d2(x, y, z, t)) > EPSILON)
					{
						cerr << "Value mismatch (" << x << ", " << y << ", " << z << ", " << t << "): "
								<< d1(x, y, z, t) << "!=" << d2(x, y, z, t) << endl;
						return 1;
					}
				}
			}
		}
	}

	return 0;
}

