/*  fabber_capi.h - Pure C API for fabber

 Used by language bindings, e.g. Python interface

 Martin Craig

 Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include "rundata.h"

// Export required symbols to Windows DLL
#ifdef _WIN32
#ifdef fabbercore_shared_EXPORTS
#define FABBER_DLL_API __declspec(dllexport)
#else
#define FABBER_DLL_API __declspec(dllimport)
#endif
#else
#define FABBER_DLL_API
#endif

extern "C" {

#define FABBER_ERR_MAXC 255
#define FABBER_ERR_FATAL -255
#define FABBER_ERR_NEWMAT -254

/**
 * Create a new context for running fabber.
 *
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return A handle to the context. This should not be used for any purpose apart from
 *         to pass to other API functions
 */
FABBER_DLL_API void* fabber_new(char* err_buf);

/**
 * Load models from a dynamic library
 *
 * @param fab Fabber context, returned by fabber_new
 * @param libpath Path to library, non NULL, not empty
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return 0 on success, <0 on failure
 */
FABBER_DLL_API int fabber_load_models(void* fab, const char* libpath, char* err_buf);

/**
 * Set the extent of the volume to be procesed.
 *
 * This must be called before any data is provided.
 *
 * @param fab Fabber context, returned by fabber_new
 * @param nx Extent in x direction
 * @param ny Extent in y direction
 * @param nz Extent in z direction
 * @param mask Array of length nx*ny*nz of 0 and 1 defining a mask region.
 *             (where mask=0, voxel is ignored, where mask!=0, voxel is included)
 *             Row-major order is used (see fabber_set_data). If NULL, data is
 *             assumed to be unmasked.
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return 0 on success, <0 on failure
 */
FABBER_DLL_API int fabber_set_extent(void* fab, int nx, int ny, int nz, const int* mask, char* err_buf);

/**
 * Destroy fabber context previously created in fabber_new
 *
 * Will not return any errors.
 *
 * @param fab Fabber context, returned by fabber_new. NULL will be ignored. Anything
 *            else will probably cause a crash.
 */
FABBER_DLL_API void fabber_destroy(void* fab);

/**
 * Set an option
 *
 * @param fab Fabber context, returned by fabber_new
 * @param key Key value, not null or empty
 * @param value Value, not null, may be empty (e.g. for boolean options)
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return 0 on success, <0 on failure
 */
FABBER_DLL_API int fabber_set_opt(void* fab, const char* key, const char* value, char* err_buf);

/**
 * Set voxel data
 *
 * @param fab Fabber context, returned by fabber_new
 * @param name Name of data. Main timeseries has name 'data'. Not null, not empty
 * @param data_size Data size in 4th dimension, i.e. number of t-slices for the main
 *                  timeseries, 1 for simple image data.
 * @param data Array of size nx*ny*nz*data_size. Will be masked using the mask supplied
 *             in fabber_new. Row-major order is used, i.e. 4th dimension varies fastest
 *             (default in C and Python Numpy)
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return 0 on success, <0 on failure
 */
FABBER_DLL_API int fabber_set_data(void* fab, const char* name, int data_size, const float* data, char* err_buf);

/**
 * Get the size of output voxel data
 *
 * @param fab Fabber context, returned by fabber_new
 * @param name Name of data. Not null, not empty
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return Size of data in 4th dimension, or <0 on error
 */
FABBER_DLL_API int fabber_get_data_size(void* fab, const char* name, char* err_buf);

/**
 * Get output voxel data
 *
 * @param fab Fabber context, returned by fabber_new
 * @param name Name of data. Not null, not empty
 * @param data_buf Data buffer of size nx*ny*nz*data_size. Data size may be obtained
 *                 from fabber_get_data_size. Row-major order is used (see fabber_set_data)
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return 0 on success, <0 on failure
 */
FABBER_DLL_API int fabber_get_data(void* fab, const char* name, float* data_buf, char* err_buf);

/**
 * Run Fabber model fitting on already configured options and data
 *
 * @param fab Fabber context, returned by fabber_new
 * @param log_bufsize Size of the log buffer. If too small, only this number of characters
 *                    will be returned.
 * @param log_buf Char buffer of size log_bufsize to receive output log
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 * @param progress_cb Function pointer which takes two integers (current voxel, total voxels). Pass NULL if not required
 *
 * @return 0 on success, <0 on failure
 */
FABBER_DLL_API int fabber_dorun(void* fab, int log_bufsize, char* log_buf, char* err_buf, void (*progress_cb)(int, int));

/**
 * Get fabber options, optionally for a specific method or model
 *
 * @param fab Fabber context, returned by fabber_new
 * @param key "method" for method option, "model" for model options, NULL or empty
 *            for general Fabber options
 * @param value Name of method or model. Ignored for general options
 * @param out_bufsize Size of the output buffer. If too small, no output is returned
 * @param out_buf Char buffer of size log_bufsize to receive output. First line will contain
 *                description of the model, method, or Fabber generally. This will be blank
 *                if no description exists, but the line will be included. The remaining lines
 *                will contain tab-separated data with each option on a single line. Order of
 *                data is option name, description, type, optional (1 or 0), default value
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return 0 on success, <0 on failure
 */
FABBER_DLL_API int fabber_get_options(void* fab, const char* key, const char* value, int out_bufsize, char* out_buf, char* err_buf);

/**
 * Get list of known models
 *
 * @param fab Fabber context, returned by fabber_new
 * @param out_bufsize Size of the output buffer. If too small, no output is returned
 * @param out_buf Char buffer of size out_bufsize to receive output. Will contain
 *                known model names separated by newlines
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return 0 on success, <0 on failure
 */
FABBER_DLL_API int fabber_get_models(void* fab, int out_bufsize, char* out_buf, char* err_buf);

/**
 * Get list of known inference methods
 *
 * @param fab Fabber context, returned by fabber_new
 * @param out_bufsize Size of the output buffer. If too small, no output is returned
 * @param out_buf Char buffer of size out_bufsize to receive output. Will contain
 *                known method names separated by newlines
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return 0 on success, <0 on failure
 */
FABBER_DLL_API int fabber_get_methods(void* fab, int out_bufsize, char* out_buf, char* err_buf);

/**
 * Get list model parameters that will be output. Note that this will depend
 * on the options specified, so must be called after all options are set
 *
 * @param fab Fabber context, returned by fabber_new
 * @param out_bufsize Size of the output buffer. If too small, no output is returned
 * @param out_buf Char buffer of size out_bufsize to receive output. Will contain
 *                model parameter names separated by newlines
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return 0 on success, <0 on failure
 */
FABBER_DLL_API int fabber_get_model_params(void* fab, int out_bufsize, char* out_buf, char* err_buf);
}
