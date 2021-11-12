/********************************************************************************
 * Copyright (c) 2021, The University Corporation for Atmospheric Research (UCAR).
 *
 * Unless noted otherwise source code is licensed under the BSD license.
 * Additional copyright and license information can be found in the LICENSE file
 * distributed with this code, or at http://mpas-dev.github.com/license.html
 ********************************************************************************/

#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>


/********************************************************************************
 *
 * gen_string_consts
 *
 * What does this program do?
 *
 ********************************************************************************/
int main(int argc, char **argv)
{
	int len;
	int i;
	int fd_in, fd_out;
	char buf[1024];

	if (argc != 3) {
		fprintf(stderr, "\nUsage: %s <input filename> <output filename>\n\n", argv[0]);
		fprintf(stderr, "\n");
		fprintf(stderr, "       If <input_filename> is -, input will be read from stdin.\n");
		fprintf(stderr, "       If <output_filename> is -, output will be written to stdout.\n");
		return 1;
	}

	if (strcmp(argv[1], "-") != 0) {
		fd_in = open(argv[1], O_RDONLY);
	} else {
		fd_in = 0;
	}

	if (strcmp(argv[2], "-") != 0) {
		fd_out = open(argv[2], O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
	} else {
		fd_out = 1;
	}

	write(fd_out, "\"", 1);

	while ((len = read(fd_in, buf, 1024)) > 0) {
		for (i = 0; i < len; i++) {
			if (buf[i] == '"') {                  /* double quote */
				write(fd_out, "\\\"", 2);
			} else if (buf[i] == '\'') {          /* single quote */
				write(fd_out, "\\\'", 2);
			} else if (buf[i] == '\n') {          /* newline character */
				write(fd_out, "\\n", 2);
			} else if (buf[i] == '\t') {          /* tab character */
				write(fd_out, "\\t", 2);
			} else {
				write(fd_out, &buf[i], 1);
			}
		}
	}

	write(fd_out, "\";\n", 3);

	if (fd_in != 0) {
		close(fd_in);
	}

	if (fd_out != 1) {
		close(fd_out);
	}

	return 0;
}
