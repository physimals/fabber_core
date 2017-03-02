#pragma once

#include <string>

#define V_MAJ 3
#define V_MIN 9
#define V_PAT 10
#define V_FL ""

/**
* Get current release version string
*
* @return A version string in form major.minor.patch, e.g. 1.2.3
*/
std::string fabber_release_version();

/**
 * Get string identifying the current source code revision.
 *
 * Currently this is a SHA1 hash of the current Git revision, if available
 *
 * @return identifying string, or 'unknown' if not available at build time
 */
std::string fabber_source_version();

/**
  * Get date of last commit, if available
  *
  * @return date string, or 'unknown' if not available at build time
  */
std::string fabber_source_date();
