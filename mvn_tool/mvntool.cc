/* mvntool.cc - Tool for adding parameters to an MVN

 Michael Chappell, FMRIB Analysis Group

 Copyright (C) 2007 University of Oxford  */
/*  CCOPYRIGHT  */

#include "dist_mvn.h"
#include "rundata_newimage.h"

#include "armawrap/newmat.h"

#include <exception>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

using namespace std;
using NEWMAT::ColumnVector;
using NEWMAT::Matrix;
using NEWMAT::SymmetricMatrix;


/* Function declarations */
void Usage(const string &errorString = "");

int main(int argc, char **argv)
{
    try
    {
        EasyLog log;
        log.StartLog(cout);
        cout << "FABBER: MVNtool" << endl;

        FabberRunDataNewimage args;
        args.SetLogger(&log);
        args.Parse(argc, argv);

        // MVN tool uses a different argument for the main
        // data, set the 'data' option to match before we might use it
        args.Set("data", args.GetString("input"));
        args.SetExtentFromData();

        if ((argc == 1) || args.ReadBool("help"))
        {
            Usage();
            return 0;
        }

        /* parse command line arguments*/
        bool verbose = args.ReadBool("v");

        string infile;
        string outfile;
        infile = args.Read("input");
        int param = 0;
        int cparam = 0;
        outfile = args.ReadWithDefault("output", infile);

        // Unfortunately, libfabber uses the output option to decide where
        // to save files. MVNtool uses it for the output file. So no we've
        // got it, we clear the value so it won't confuse libfabber when
        // we get to save
        args.Unset("output");

        bool ins;
        bool write;

        double val = 0;
        double var = 0;
        string valimfile;
        string varimfile;
        bool bval = false;
        bool bvar = false;
        bool cvar = false;
        /* Choose what we want to do to/with the parameter - default is to read */
        ins = args.ReadBool("new");     //insert a new parameter
        write = args.ReadBool("write"); //overwrite an existing parameter

        // determine which parameter we are reading/writing etc
        string plistfile = args.ReadWithDefault("param-list", "");
        if (plistfile == "")
        { //use must have specified a parameter number
            param = convertTo<int>(args.Read("param"));

            //deal with --cvar option
            cparam = convertTo<int>(args.ReadWithDefault("cvar", "0"));
        }
        else
        {
            string paramname = args.Read("param");
            // read in parameter names
            ifstream paramFile((plistfile).c_str());
            vector<string> paramNames;

            //NOTE covariance parameter not setup for parameter lists

            string currparam;
            while (paramFile.good())
            {
                getline(paramFile, currparam);
                paramNames.push_back(currparam);
            }
            paramNames.pop_back(); //remove final empty line assocaited with eof

            //user will have specified a parameter name in the list
            string nplistfile = args.ReadWithDefault("new-param-list", "");
            if (nplistfile == "")
            { //simply deal with the parameter list
                for (unsigned int i = 0; i < paramNames.size(); i++)
                {
                    if (paramname.compare(paramNames[i]) == 0)
                        param = i + 1;
                }
                if (param == 0)
                    throw std::runtime_error("Cannot find specfied parameter name in list");
            }
            else
            { //deal with current and new parameter lists
                // if we are here we must be inserting a new parameter
                ins = true;

                //load new parameter list
                ifstream newparamFile((nplistfile).c_str());
                vector<string> newparamNames;
                while (newparamFile.good())
                {
                    getline(newparamFile, currparam);
                    newparamNames.push_back(currparam);
                }
                newparamNames.pop_back(); //remove final empty line assocaited with eof

                // check parameter name is not in old list
                for (unsigned int i = 0; i < paramNames.size(); i++)
                {
                    if (paramname.compare(paramNames[i]) == 0)
                        throw std::runtime_error(
                            "Parameter name found in parameter list for this MVN, cannot insert an identical parameter");
                }
                //find parameter in new list
                int newparam = 0;
                for (unsigned int i = 0; i < newparamNames.size(); i++)
                {
                    if (paramname.compare(newparamNames[i]) == 0)
                        newparam = i + 1;
                }
                if (newparam == 0)
                    throw std::runtime_error("Cannot find specfied parameter name in new parameter name list");
                //cout << newparam<<endl;

                if (newparam == 1)
                {
                    param = 1;
                }
                else
                {
                    //get name of previous parameter in list
                    string preparamname = newparamNames[newparam - 1 - 1];
                    //cout << preparamname << endl;
                    // find this in the old list
                    for (unsigned int i = 0; i < paramNames.size(); i++)
                    {
                        if (preparamname.compare(paramNames[i]) == 0)
                            param = i + 1;
                    }
                    if (param == 0)
                        throw std::runtime_error(
                            "Cannot complete this operation since the new list contains other parameters not present in the old list.\n You may be able to get around this by adding new parameters in a different order:\n new parameters that are sequential in the new list should be added starting with the earliest appearing.");
                    // param currently hold location of previous parameter to the one we want to insert, so increment
                    param++;
                    //cout << param <<endl;
                }

                //need to write out the updated parameter list
                vector<string> outParams;
                for (int i = 0; i < param - 1; i++)
                {
                    outParams.push_back(paramNames[i]);
                }
                outParams.push_back(paramname);
                for (unsigned int i = param - 1; i < paramNames.size(); i++)
                {
                    outParams.push_back(paramNames[i]);
                }

                string paramlistout = args.ReadWithDefault("out-param-file", "");
                ofstream outparamFile((paramlistout).c_str());
                if (paramlistout == "")
                {
                    for (unsigned int i = 0; i < outParams.size(); i++)
                    {
                        outparamFile << outParams[i] << endl;
                    }
                }
            }
        }

        // Deal with --val and --var
        if (ins | write)
        {
            if (ins & write)
            {
                throw FabberRunDataError("Cannot insert and write at same time - choose either --new or --write");
            }

            valimfile = args.ReadWithDefault("valim", "");
            varimfile = args.ReadWithDefault("varim", "");

            val = convertTo<double>(args.ReadWithDefault("val", "-1e-6"));
            var = convertTo<double>(args.ReadWithDefault("var", "-1e-6"));
        }
        else
        {
            if (outfile == infile)
                throw MandatoryOptionMissing("output");

            bval = args.ReadBool("val");
            bvar = args.ReadBool("var");
            if (cparam > 0)
                cvar = true;

            if (bval && bvar)
            {
                throw FabberRunDataError("Cannot output value and variance at same time - choose either --val or --val");
            }
            if (bval && cvar)
            {
                throw FabberRunDataError("Cannot output value and covariance at same time");
            }
            if (bvar && cvar)
            {
                throw FabberRunDataError("Cannot output variance and covariance at same time");
            }
            if (!bval && !bvar && !cvar)
            {
                throw FabberRunDataError(
                    "Please select whether you want to extract the value (--val) or the variance (--var) or a covaraince (--cvar)");
            }
        }

        if (verbose)
            cout << "Read file" << endl;
        vector<MVNDist *> vmvnin;
        MVNDist::Load(vmvnin, "input", args, &log);

        if (ins | write)
        {
            /* section to deal with writing to or inserting into an MVN*/

            // deal with the values we are going to insert
            ColumnVector inmean(vmvnin.size());
            ColumnVector invar(vmvnin.size());
            if (valimfile == "")
            {
                inmean = val;
            }
            else
            {
                inmean = args.GetVoxelData("valim").AsColumn();
            }
            if (varimfile == "")
            {
                invar = var;
            }
            else
            {
                invar = args.GetVoxelData("varim").AsColumn();
            }

            vector<MVNDist *> vmvnout(vmvnin);

            int oldsize;
            MVNDist mvnin;
            mvnin = *vmvnin[1];
            oldsize = mvnin.GetSize();

            if (ins)
            {
                if (verbose)
                    cout << "Inserting new parameter" << endl;
                if (param > oldsize + 1)
                    throw FabberRunDataError("Cannot insert parameter here, not enough parameters in existing MVN");
            }
            else
            {
                if (param > oldsize)
                    throw FabberRunDataError("Cannot edit this parameter, not enough parameters in existing MVN");
            }
            SymmetricMatrix mvncov;
            MVNDist mvnnew(1);
            mvnnew.means = 0;
            SymmetricMatrix covnew(1);
            covnew(1, 1) = 0;
            mvnnew.SetCovariance(covnew);
            MVNDist mvnout;

            /* Loop over each enrty in mvnin - each voxel! */
            for (unsigned v = 0; v < vmvnin.size(); v++)
            {
                mvnin = *vmvnin[v];

                if (ins)
                { /* insert new parameter */
                    //cout << "Add new parameter" << endl;

                    if (param - 1 >= 1)
                    {
                        MVNDist mvn1(param - 1);
                        mvn1.CopyFromSubmatrix(mvnin, 1, param - 1, 0);
                        mvnout = MVNDist(mvn1, mvnnew);
                    }
                    else
                    {
                        mvnout = mvnnew;
                    }

                    if (oldsize >= param)
                    {
                        MVNDist mvn2(oldsize + 1 - param);
                        mvn2.CopyFromSubmatrix(mvnin, param, oldsize, 0);
                        mvnout = MVNDist(mvnout, mvn2);
                    }

                    mvncov = mvnout.GetCovariance();
                }
                else
                {
                    mvnout = mvnin;
                    mvncov = mvnin.GetCovariance();
                }

                /* Set parameters mean value and variance*/
                //cout << "Set parameter mean and variance" << endl;
                mvnout.means(param) = inmean(v + 1);
                mvncov(param, param) = invar(v + 1);
                mvnout.SetCovariance(mvncov);

                vmvnout[v] = new MVNDist(mvnout);
            }

            /* Save MVN to output */
            if (verbose)
                cout << "Save file" << endl;
            MVNDist::Save(vmvnout, outfile, args);
        }

        else
        {
            /* Section to deal with reading a parameter out to an image */
            int nVoxels = vmvnin.size();
            Matrix image;
            image.ReSize(1, nVoxels);

            if (bval)
            {
                if (verbose)
                    cout << "Extracting value for parameter:" << param << endl;
                for (int vox = 1; vox <= nVoxels; vox++)
                {
                    image(1, vox) = vmvnin[vox - 1]->means(param);
                }
            }
            else if (bvar)
            {
                if (verbose)
                    cout << "Extracting variance for parameter:" << param << endl;
                for (int vox = 1; vox <= nVoxels; vox++)
                {
                    image(1, vox) = vmvnin[vox - 1]->GetCovariance()(param, param);
                }
            }
            else if (cvar)
            {
                if (verbose)
                    cout << "Extracting co-variance for parameter " << param << "with parameter" << cparam << endl;
                for (int vox = 1; vox <= nVoxels; vox++)
                {
                    image(1, vox) = vmvnin[vox - 1]->GetCovariance()(param, cparam);
                }
            }

            if (verbose)
                cout << "Writing output file" << endl;

            args.SaveVoxelData(outfile, image);
        }

        if (verbose)
            cout << "Done." << endl;

        return 0;
    }
    catch (std::runtime_error &e)
    {
        cout << e.what() << endl;
    }
    catch (const std::exception &e)
    {
        cout << e.what() << endl;
        Usage();
    }
    catch (...)
    {
        cout << "There was an error!" << endl;
    }

    return 1;
}

void Usage(const string &errorString)
{
    cout << "\nUsage: mvntool <arguments>\n"
         << "Arguments are mandatory unless they appear in [brackets].\n\n";

    cout << " --help : Prints this information." << endl
         << " --input=<MVNfile> : Name of input MVN file." << endl
         << " --param=<n> : Number of parameter to read/replace/insert." << endl
         << endl
         << " Extract behaviour (default):" << endl
         << " --output=<NIFIfile> : Name of file for output" << endl
         << "   [--val] : Write parameter value to file." << endl
         << "   [--var] : Write parameter variance to file." << endl
         << endl
         << " Writing behaviour:" << endl
         << "   [--write] : Overwrite an existing parameter" << endl
         << "   [--new] : Insert a new parameter."
         << endl
         << "[--output]=<MVNfile> : Name of file for output, overwrites input if not set" << endl
         << " --valim=<NIFITfile> : Image to write for mean of parameter." << endl
         << " --varim=<NIFITfile> : Image to write for variance of parameter." << endl
         << " --val=<mean_value>  : Mean value for parameter to be written." << endl
         << " --var=<variance>    : Variance of parameter to be written." << endl
         << endl;
}
