#pragma once

extern "C"
{

#define FABBER_ERR_MAXC 255
#define FABBER_ERR_FATAL -255

/**
 * Create a new context for running fabber.
 *
 * Currently a mask must be specified
 *
 * @param nx Extent in x direction
 * @param ny Extent in y direction
 * @param nz Extent in z direction
 * @param mask Array of length nx*ny*nz of 0 and 1 defining a mask region (1=include voxel)
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return A handle to the context. This should not be used for any purpose apart from
 *         to pass to other API functions
 */
void *fabber_new(int nx, int ny, int nz, const int *mask, char *err_buf);

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
int fabber_set_opt(void *fab, const char *key, const char *value, char *err_buf);

/**
 * Set voxel data
 *
 * @param fab Fabber context, returned by fabber_new
 * @param name Name of data. Main timeseries has name 'data'. Not null, not empty
 * @param data_size Data size in 4th dimension, i.e. number of t-slices for the main
 *                  timeseries, 1 for simple image data.
 * @param data Array of size nx*ny*nz*data_size. Will be masked using the mask supplied
 *             in fabber_new
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return 0 on success, <0 on failure
 */
int fabber_set_data(void *fab, const char *name, int data_size, const float *data, char *err_buf);

/**
 * Get the size of output voxel data
 *
 * @param fab Fabber context, returned by fabber_new
 * @param name Name of data. Not null, not empty
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return Size of data in 4th dimension, or <0 on error
 */
int fabber_get_data_size(void *fab, const char *name, char *err_buf);

/**
 * Get output voxel data
 *
 * @param fab Fabber context, returned by fabber_new
 * @param name Name of data. Not null, not empty
 * @param data_buf Data buffer of size nx*ny*nz*data_size. Data size may be obtained
 *                 from fabber_get_data_size.
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return 0 on success, <0 on failure
 */
int fabber_get_data(void *fab, const char *name, float *data_buf, char *err_buf);

/**
 * Run Fabber model fitting on already configured options and data
 *
 * @param fab Fabber context, returned by fabber_new
 * @param log_bufsize Size of the log buffer. If too small, only this number of characters
 *                    will be returned.
 * @param log_buf Char buffer of size log_bufsize to receive output log
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return 0 on success, <0 on failure
 */
int fabber_dorun(void *fab, int log_bufsize, char *log_buf, char *err_buf);

/**
 * Destroy fabber context previously created in fabber_new
 */
int fabber_destroy(void *fab, char *err_buf);

/**
 * Get fabber options, optionally for a specific method or model
 *
 * @param fab Fabber context, returned by fabber_new
 * @param key "method" for method option, "model" for model options, NULL or empty
 *            for general Fabber options
 * @param value Name of method or model. Ignored for general options
 * @param out_bufsize Size of the output buffer. If too small, no output is returned
 * @param out_buf Char buffer of size log_bufsize to receive output. Will contain
 *                tab-separated data with each option on a single line. Order of data
 *                is option name, description, type, optional (1 or 0), default value
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return 0 on success, <0 on failure
 */
int fabber_get_options(const char *key, const char *value, int out_bufsize, char *out_buf, char *err_buf);

/**
 * Get list of known models
 *
 * @param fab Fabber context, returned by fabber_new
 * @param out_bufsize Size of the output buffer. If too small, no output is returned
 * @param out_buf Char buffer of size log_bufsize to receive output. Will contain
 *                known model names separated by newlines
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return 0 on success, <0 on failure
 */
int fabber_get_models(int out_bufsize, char *out_buf, char *err_buf);

/**
 * Get list of known inference methods
 *
 * @param fab Fabber context, returned by fabber_new
 * @param out_bufsize Size of the output buffer. If too small, no output is returned
 * @param out_buf Char buffer of size log_bufsize to receive output. Will contain
 *                known method names separated by newlines
 * @param err_buf Optional buffer for error message. Max message length=FABBER_ERR_MAXC
 *
 * @return 0 on success, <0 on failure
 */
int fabber_get_methods(int out_bufsize, char *out_buf, char *err_buf);

}
