/********************************************************************************
 * Copyright (c) 2021, The University Corporation for Atmospheric Research (UCAR).
 *
 * Unless noted otherwise source code is licensed under the BSD license.
 * Additional copyright and license information can be found in the LICENSE file
 * distributed with this code, or at http://mpas-dev.github.com/license.html
 ********************************************************************************/

#include <string.h>

static const char build_opts[] = 
#include "build_opts.string"

const char *get_build_opts(int *len);


/********************************************************************************
 *
 * get_build_opts
 *
 * One-line summary.
 *
 * This function does
 * something that can be described here.
 *
 ********************************************************************************/
const char *get_build_opts(int *len)
{
	if (len != NULL) {
		*len = (int)strlen(build_opts);
	}

	return build_opts;
}
