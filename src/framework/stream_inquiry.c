/*
 * Copyright (c) 2023, The University Corporation for Atmospheric Research (UCAR).
 *
 * Unless noted otherwise source code is licensed under the BSD license.
 * Additional copyright and license information can be found in the LICENSE file
 * distributed with this code, or at http://mpas-dev.github.com/license.html
 */

#include <stddef.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include "ezxml.h"

#ifdef _MPI
#include "mpi.h"
#endif

#define MSGSIZE 256


/*
 *  Interface routines for writing log messages; defined in mpas_log.F
 *  messageType_c may be any of "MPAS_LOG_OUT", "MPAS_LOG_WARN", "MPAS_LOG_ERR", or "MPAS_LOG_CRIT"
 */
void mpas_log_write_c(const char *message_c, const char *messageType_c);


/*********************************************************************************
 *
 *  Function: read_and_broadcast
 *
 *  Reads the contents of a file into a buffer in distributed-memory parallel code.
 *
 *  The buffer buf is allocated with size bufsize, which will be exactly the
 *  number of bytes in the file fname. Only the master task will actually read the
 *  file, and the contents are broadcast to all other tasks. The mpi_comm argument
 *  is a Fortran MPI communicator used to determine which task is the master task.
 *
 *  A return code of 0 indicates the file was successfully read and broadcast to
 *  all MPI tasks that belong to the communicator.
 *
 *********************************************************************************/
int read_and_broadcast(const char *fname, int mpi_comm, char **buf, size_t *bufsize)
{
	int iofd;
	int rank;
	struct stat s;
	int err;
	char msgbuf[MSGSIZE];

#ifdef _MPI
	MPI_Comm comm;

	comm = MPI_Comm_f2c((MPI_Fint)mpi_comm);
	err = MPI_Comm_rank(comm, &rank);
#else
	rank = 0;
#endif

	if (rank == 0) {
		iofd = open(fname, O_RDONLY);
		if (iofd <= 0) {
			snprintf(msgbuf, MSGSIZE, "Could not open run-time I/O config file %s", fname);
			mpas_log_write_c(msgbuf, "MPAS_LOG_ERR");
			return 1;
		}

		fstat(iofd, &s);
		*bufsize = (size_t)s.st_size;
#ifdef _MPI
		err = MPI_Bcast((void *)bufsize, (int)sizeof(size_t), MPI_BYTE, 0, comm);
#endif
	
		*buf = (char *)malloc(*bufsize);
		err = read(iofd, (void *)(*buf), *bufsize);

#ifdef _MPI
		err = MPI_Bcast((void *)(*buf), (int)(*bufsize), MPI_CHAR, 0, comm);
#endif
	}
	else {
#ifdef _MPI
		err = MPI_Bcast((void *)bufsize, (int)sizeof(size_t), MPI_BYTE, 0, comm);
#endif
		*buf = (char *)malloc(*bufsize);

#ifdef _MPI
		err = MPI_Bcast((void *)(*buf), (int)(*bufsize), MPI_CHAR, 0, comm);
#endif
	}

	return 0;
}

/********************************************************************************
 *
 * parse_streams_file
 *
 * One-line summary
 *
 * Detailed paragraph(s) about the function.
 * Detailed paragraph(s) about the function.
 * Detailed paragraph(s) about the function.
 * Detailed paragraph(s) about the function.
 *
 ********************************************************************************/
ezxml_t parse_streams_file(int mpi_comm, const char *filename)
{
	char *xml_buf;
	size_t bufsize;

	if (read_and_broadcast(filename, mpi_comm, &xml_buf, &bufsize) != 0) {
		return NULL;
	}

	return ezxml_parse_str(xml_buf, bufsize);
}

/********************************************************************************
 *
 * free_streams_file
 *
 * One-line summary
 *
 * Detailed paragraph(s) about the function.
 * Detailed paragraph(s) about the function.
 * Detailed paragraph(s) about the function.
 * Detailed paragraph(s) about the function.
 *
 ********************************************************************************/
void free_streams_file(ezxml_t xmltree)
{
	ezxml_free(xmltree);
}


/********************************************************************************
 *
 * query_streams_file
 *
 * One-line summary
 *
 * Detailed paragraph(s) about the function.
 * Detailed paragraph(s) about the function.
 * Detailed paragraph(s) about the function.
 * Detailed paragraph(s) about the function.
 *
 ********************************************************************************/
int query_streams_file(ezxml_t xmltree, const char *streamname, const char *attname, const char **attvalue)
{
	ezxml_t stream_xml;
	const char *streamID;
	const char *attval_local;

	for (stream_xml = ezxml_child(xmltree, "immutable_stream"); stream_xml; stream_xml = ezxml_next(stream_xml)) {
		streamID = ezxml_attr(stream_xml, "name");

		if (strcmp(streamID, streamname) == 0) {
	                if (attname != NULL) {
				*attvalue = ezxml_attr(stream_xml, attname);
			}
			return 1;
		}
	}

	return 0;
}

#if 0
/*********************************************************************************
 *
 *  Function: query_stream_attribute
 *
 *  Parses an XML file and searches for the stream whose name matches the 'streamname'
 *  argument; then, returns the associated attributes for that stream in
 *  the 'filename', 'ref_time', 'filename_interval', and 'io_type' arguments.
 *
 *  The fname argument provides the name of the XML file that contains the stream
 *  definitions, and mpi_comm is the Fortran MPI communicator used by MPAS.
 *
 *********************************************************************************/
void query_stream_attribute(int *mpi_comm, char *fname, char *streamname, char *attname, char *attval, int *status)
{
	char *xml_buf;
	size_t bufsize;
	ezxml_t streams;
	ezxml_t stream_xml;
	const char *streamID;
	const char *attval_local;
	int found;
	char msgbuf[MSGSIZE];

	*status = 0;

	if (read_and_broadcast(fname, mpi_comm, &xml_buf, &bufsize) != 0) {
		*status = 1;
		return;
	}

	streams = ezxml_parse_str(xml_buf, bufsize);
	if (!streams) {
		snprintf(msgbuf, MSGSIZE, "Problems encountered while parsing run-time I/O config file %s", fname);
		mpas_log_write_c(msgbuf, "MPAS_LOG_ERR");
		*status = 1;
		return;
	}	

	found = 0;
	for (stream_xml = ezxml_child(streams, "immutable_stream"); stream_xml; stream_xml = ezxml_next(stream_xml)) {
		streamID = ezxml_attr(stream_xml, "name");

		if (strcmp(streamID, streamname) == 0) {
			found = 1;
			attval_local = ezxml_attr(stream_xml, attname);
			if (attval_local == NULL) {
				*status = 1;
				return;
			}

			strcpy(attval, attval_local);
		}
	}
	if (found == 0) {
		*status = 1;
		return;
	}
}
#endif
