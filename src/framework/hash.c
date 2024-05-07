#include <stdint.h>
#include <stddef.h>
#include "hash.h"

uint64_t string_hash_64(const uint8_t *s)
{
	uint64_t hash = 5381;
	uint8_t c;

	if (s == NULL) {
		return 0ul;
	}

	while ((c = *s++)) {
		hash = ((hash << 5) + hash) ^ c; /* hash * 33 ^ c */
	}

	return hash;
}
