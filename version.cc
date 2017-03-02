#include "version.h"

#include <string>
#include <sstream>

std::string fabber_release_version()
{
    std::stringstream v;
    v << V_MAJ << "." << V_MIN << "." << V_PAT << V_FL;
    return v.str();
}

std::string fabber_source_version()
{
#ifdef GIT_SHA1
    return GIT_SHA1;
#else
    return "unknown";
#endif
}

std::string fabber_source_date()
{
#ifdef GIT_DATE
    return GIT_DATE;
#else
    return "unknown";
#endif
}
