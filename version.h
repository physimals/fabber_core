#pragma once

#include <string>

/**
 * Get string identifying the current source code revision.
 *
 * This will be of the form <tag>[-<commits>-<hash>] as returned by 
 * git describe. For a tagged release the commits and has are not included
 * and the version will normally be in the form 'vX.Y.Z'
 *
 * @return identifying string, or 'unknown' if not available at build time
 */
std::string fabber_version();

/**
  * Get date of last commit, if available
  *
  * @return date string, or 'unknown' if not available at build time
  */
std::string fabber_source_date();
