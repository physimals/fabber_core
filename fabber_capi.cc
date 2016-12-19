#include "fabber_capi.h"

#include "dataset.h"
#include "fabber_io.h"
#include "easylog.h"
#include "fwdmodel.h"
#include "inference.h"
#include "setup.h"

#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

#include <string.h>
#include <stdio.h>

using NEWMAT::Matrix;
using namespace std;

static int fabber_err(int code, const char *msg, char *err_buf)
{
	// Error buffer is optional
	if (!err_buf)
		return code;

	strncpy(err_buf, msg, FABBER_ERR_MAXC - 1);
	err_buf[FABBER_ERR_MAXC - 1] = '\0';
	return code;
}

void *fabber_new(int nx, int ny, int nz, const int *mask, char *err_buf)
{
	if (!mask)
	{
		fabber_err(FABBER_ERR_FATAL, "Mask is NULL", err_buf);
		return NULL;
	}
	if ((nx <= 0) || (ny <= 0) || (nz <= 0))
	{
		fabber_err(FABBER_ERR_FATAL, "Dimensions must be >0", err_buf);
		return NULL;
	}

	try
	{
		FabberSetup::SetupDefaults();
		FabberIoCarray *io = new FabberIoCarray(nx, ny, nz, mask);
		FabberRunData* rundata = new FabberRunData(io);
		return rundata;
	} catch (...)
	{
		// FIXME theoretical memory leak of io
		fabber_err(FABBER_ERR_FATAL, "Failed to allocate memory for run data", err_buf);
		return NULL;
	}
}

int fabber_set_opt(void *fab, const char *key, const char *value, char *err_buf)
{
	if (!fab)
		return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);

	if (!key || !value)
		return fabber_err(FABBER_ERR_FATAL, "Option key or value is NULL", err_buf);

	try
	{
		FabberRunData* rundata = (FabberRunData*) fab;
		if (strcmp(key, "loadmodels") == 0)
		{
			FwdModel::LoadFromDynamicLibrary(value);
		}
		else
		{
			rundata->Set(key, value);
		}
		return 0;
	} catch (exception &e)
	{
		return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
	}
}

int fabber_set_data(void *fab, const char *name, int data_size, const float *data, char *err_buf)
{
	if (!fab)
		return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
	if (!data)
		return fabber_err(FABBER_ERR_FATAL, "Data buffer is NULL", err_buf);
	if (!name)
		return fabber_err(FABBER_ERR_FATAL, "Data name is NULL", err_buf);
	if (data_size <= 0)
		return fabber_err(FABBER_ERR_FATAL, "Data size must be >0", err_buf);

	FabberRunData* rundata = (FabberRunData*) fab;
	try
	{
		((FabberIoCarray*) rundata->GetIo())->SetVoxelData(name, data_size, data);
	} catch (exception &e)
	{
		return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
	} catch (...)
	{
		return fabber_err(FABBER_ERR_FATAL, "Error setting data", err_buf);
	}
	return 0;
}

int fabber_get_data_size(void *fab, const char *name, char *err_buf)
{
	if (!fab)
		return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
	if (!name)
		return fabber_err(FABBER_ERR_FATAL, "Data name is NULL", err_buf);

	FabberRunData* rundata = (FabberRunData*) fab;
	try
	{
		Matrix data = rundata->GetVoxelData(name);
		return data.Nrows();
	} catch (DataNotFound &e)
	{
		return fabber_err(-1, "Data not found", err_buf);
	} catch (...)
	{
		return fabber_err(FABBER_ERR_FATAL, "Error getting data", err_buf);
	}
}

int fabber_get_data(void *fab, const char *name, float *data_buf, char *err_buf)
{
	if (!fab)
		return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
	if (!name)
		return fabber_err(FABBER_ERR_FATAL, "Data name is NULL", err_buf);
	if (!data_buf)
		return fabber_err(FABBER_ERR_FATAL, "Data name is NULL", err_buf);

	FabberRunData* rundata = (FabberRunData*) fab;
	try
	{
		((FabberIoCarray*) rundata->GetIo())->GetVoxelData(name, data_buf);
		return 0;
	} catch (DataNotFound &e)
	{
		return fabber_err(-1, "Data not found", err_buf);
	} catch (...)
	{
		return fabber_err(FABBER_ERR_FATAL, "Error getting data", err_buf);
	}
}

int fabber_dorun(void *fab, int log_bufsize, char *log_buf, char *err_buf)
{
	if (!fab)
		return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
	if (!log_buf)
		return fabber_err(FABBER_ERR_FATAL, "Log buffer is NULL", err_buf);
	if (!err_buf)
		return fabber_err(FABBER_ERR_FATAL, "Error buffer is NULL", err_buf);

	FabberRunData* rundata = (FabberRunData*) fab;
	stringstream log;
	try
	{
		EasyLog::CurrentLog().StartLog(log);
		rundata->Run();
		EasyLog::CurrentLog().ReissueWarnings();
	} catch (const DataNotFound& e)
	{
		EasyLog::CurrentLog().ReissueWarnings();
		LOG_ERR("Data not found:\n  " << e.what() << endl);
		fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
	} catch (const Invalid_option& e)
	{
		EasyLog::CurrentLog().ReissueWarnings();
		LOG_ERR("Invalid_option exception caught in fabber:\n  " << e.what() << endl);
		fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
	} catch (const exception& e)
	{
		EasyLog::CurrentLog().ReissueWarnings();
		LOG_ERR("STL exception caught in fabber:\n  " << e.what() << endl);
		fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
	} catch (NEWMAT::Exception& e)
	{
		EasyLog::CurrentLog().ReissueWarnings();
		LOG_ERR("NEWMAT exception caught in fabber:\n  " << e.what() << endl);
		fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
	} catch (...)
	{
		EasyLog::CurrentLog().ReissueWarnings();
		LOG_ERR("Some other exception caught in fabber!" << endl);
		fabber_err(FABBER_ERR_FATAL, "Unrecognized exception", err_buf);
	}

	EasyLog::CurrentLog().StopLog();

	if (log_bufsize < 0)
		return fabber_err(FABBER_ERR_FATAL, "Log buffer size is < 0", err_buf);
	if (log_bufsize > 0 && !log_buf)
		return fabber_err(FABBER_ERR_FATAL, "Log buffer is NULL", err_buf);

	strncpy(log_buf, log.str().c_str(), log_bufsize - 1);
	log_buf[log_bufsize - 1] = '\0';

	return 0;
}

void fabber_destroy(void *fab)
{
	if (fab)
	{
		// Get rid of registered models etc
		FabberSetup::Destroy();
		FabberRunData* rundata = (FabberRunData*) fab;
		delete rundata->GetIo();
		delete rundata;
	}
}

int fabber_get_options(void *fab, const char *key, const char *value, int out_bufsize, char *out_buf, char *err_buf)
{
	if (!fab)
		return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
	if (out_bufsize < 0)
		return fabber_err(FABBER_ERR_FATAL, "Output buffer size is < 0", err_buf);
	if (!out_buf)
		return fabber_err(FABBER_ERR_FATAL, "Output buffer is NULL", err_buf);
	if (key && !value)
		return fabber_err(FABBER_ERR_FATAL, "Key specified but no value", err_buf);

	vector<OptionSpec> options;
	string desc;
	try
	{
		if (!key || (strlen(key) == 0))
		{
			FabberRunData::GetOptions(options);
		}
		else if (strcmp(key, "model") == 0)
		{
			std::auto_ptr<FwdModel> model(FwdModel::NewFromName(value));
			desc = model->GetDescription();
			model->GetOptions(options);
		}
		else if (strcmp(key, "method") == 0)
		{
			std::auto_ptr<InferenceTechnique> method(InferenceTechnique::NewFromName(value));
			desc = method->GetDescription();
			method->GetOptions(options);
		}

		stringstream out;
		vector<OptionSpec>::iterator iter;
		for (iter = options.begin(); iter != options.end(); iter++)
		{
			out << iter->name << "\t" << iter->description << "\t" << iter->type << "\t" << iter->optional << "\t"
					<< iter->def << endl;
		}
		string outstr = out.str();

		if (outstr.size() >= out_bufsize)
		{
			return fabber_err(-1, "Buffer too small", err_buf);
		}
		strncpy(out_buf, outstr.c_str(), outstr.size());
		return 0;
	} catch (exception &e)
	{
		return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
	} catch (...)
	{
		return fabber_err(FABBER_ERR_FATAL, "Error in get_opts", err_buf);
	}
}

int fabber_get_models(void *fab, int out_bufsize, char *out_buf, char *err_buf)
{
	if (!fab)
		return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
	if (out_bufsize < 0)
		return fabber_err(FABBER_ERR_FATAL, "Output buffer size is < 0", err_buf);
	if (!out_buf)
		return fabber_err(FABBER_ERR_FATAL, "Output buffer is NULL", err_buf);

	try
	{
		vector<string> known = FwdModel::GetKnown();
		stringstream out;
		vector<string>::iterator iter;
		for (iter = known.begin(); iter != known.end(); iter++)
		{
			out << *iter << endl;
		}
		string outstr = out.str();
		if (outstr.size() >= out_bufsize)
		{
			return fabber_err(-1, "Buffer too small", err_buf);
		}
		strncpy(out_buf, outstr.c_str(), outstr.size());
		return 0;
	} catch (exception &e)
	{
		return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
	} catch (...)
	{
		return fabber_err(FABBER_ERR_FATAL, "Error in get_opts", err_buf);
	}
}

int fabber_get_methods(void *fab, int out_bufsize, char *out_buf, char *err_buf)
{
	if (!fab)
		return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
	if (out_bufsize < 0)
		return fabber_err(FABBER_ERR_FATAL, "Output buffer size is < 0", err_buf);
	if (!out_buf)
		return fabber_err(FABBER_ERR_FATAL, "Output buffer is NULL", err_buf);

	try
	{
		vector<string> known = InferenceTechnique::GetKnown();
		stringstream out;
		vector<string>::iterator iter;
		for (iter = known.begin(); iter != known.end(); iter++)
		{
			out << *iter << endl;
		}
		string outstr = out.str();
		if (outstr.size() >= out_bufsize)
		{
			return fabber_err(-1, "Buffer too small", err_buf);
		}
		strncpy(out_buf, outstr.c_str(), outstr.size());
		return 0;
	} catch (exception &e)
	{
		return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
	} catch (...)
	{
		return fabber_err(FABBER_ERR_FATAL, "Error in get_opts", err_buf);
	}
}

int fabber_get_model_params(void *fab, const char *model_name, int out_bufsize, char *out_buf, char *err_buf)
{
	if (!fab)
		return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
	if (!model_name || (strlen(model_name) == 0))
		return fabber_err(FABBER_ERR_FATAL, "Model name is NULL or empty", err_buf);
	if (out_bufsize < 0)
		return fabber_err(FABBER_ERR_FATAL, "Output buffer size is < 0", err_buf);
	if (!out_buf)
		return fabber_err(FABBER_ERR_FATAL, "Output buffer is NULL", err_buf);

	try
	{
		FabberRunData* rundata = (FabberRunData*) fab;
		std::auto_ptr<FwdModel> model(FwdModel::NewFromName(model_name));
		model->Initialize(*rundata);
		vector<string> params;
		model->NameParams(params);
		stringstream out;
		vector<string>::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			out << *iter << endl;
		}
		string outstr = out.str();
		if (outstr.size() >= out_bufsize)
		{
			return fabber_err(-1, "Buffer too small", err_buf);
		}
		strncpy(out_buf, outstr.c_str(), outstr.size());
		return 0;
	} catch (exception &e)
	{
		return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
	} catch (...)
	{
		return fabber_err(FABBER_ERR_FATAL, "Error in get_opts", err_buf);
	}
}
