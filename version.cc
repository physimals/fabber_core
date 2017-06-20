#include "version.h"

#include <string>

std::string fabber_version()
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
