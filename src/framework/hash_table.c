#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "hash_table.h"
#include "hash.h"


/********************************************************************************
 *
 * hash_table_init
 *
 * One-line summary
 *
 * Full description
 *
 ********************************************************************************/
int hash_table_init(hash_table *t, size_t n)
{
	if (t == NULL || n <= 0) {
		return HT_ERR_NULL_ARGUMENT;
	}

	t->size = n;

	t->table = malloc(sizeof(hash_node *) * n);
	memset(t->table, 0, sizeof(hash_node *) * n);

	return HT_SUCCESS;
}


/********************************************************************************
 *
 * hash_table_insert
 *
 * One-line summary
 *
 * Full description
 *
 ********************************************************************************/
int hash_table_insert(hash_table *t, hash_node *node)
{
	uint64_t hash;
	hash_node *cursor;

	if (t == NULL || node == NULL) {
		return HT_ERR_NULL_ARGUMENT;
	}

	hash = string_hash_64(node->key) % t->size;

	node->next = NULL;

	if (t->table[hash] == NULL) {
		t->table[hash] = node;
	} else {
		cursor = t->table[hash];
		while (cursor->next) {
			if (!strncmp(cursor->key, node->key, KEYLEN_MAX)) {
				return HT_ERR_DUPLICATE_KEY;
			}
			cursor = cursor->next;
		}
		cursor->next = node;
	}

	return HT_SUCCESS;
}


/********************************************************************************
 *
 * hash_table_remove
 *
 * One-line summary
 *
 * Full description
 *
 ********************************************************************************/
hash_node * hash_table_remove(hash_table *t, const char *key)
{
	uint64_t hash;
	hash_node *cursor, *prev_cursor;

	if (t == NULL || key == NULL) {
		return NULL;
	}

	hash = string_hash_64(key) % t->size;

	if (t->table[hash] == NULL) {
		return NULL;
	} else {
		cursor = t->table[hash];
		if (!strncmp(cursor->key, key, KEYLEN_MAX)) {
			t->table[hash] = cursor->next;
			return cursor;
		}
		prev_cursor = cursor;
		cursor = cursor->next;
		while (cursor) {
			if (!strncmp(cursor->key, key, KEYLEN_MAX)) {
				prev_cursor->next = cursor->next;
				return cursor;
			}
			prev_cursor = cursor;
			cursor = cursor->next;
		}
		return NULL;
	}
}


/********************************************************************************
 *
 * hash_table_query
 *
 * One-line summary
 *
 * Full description
 *
 ********************************************************************************/
hash_node * hash_table_query(hash_table *t, const char *key)
{
	uint64_t hash;
	hash_node *cursor;

	if (t == NULL || key == NULL) {
		return NULL;
	}

	hash = string_hash_64(key) % t->size;

	if (t->table[hash] == NULL) {
		return NULL;
	} else {
		cursor = t->table[hash];
		while (cursor) {
			if (!strncmp(cursor->key, key, KEYLEN_MAX)) {
				return cursor;
			}
			cursor = cursor->next;
		}
		return NULL;
	}
}


/********************************************************************************
 *
 * hash_table_print
 *
 * One-line summary
 *
 * Full description
 *
 ********************************************************************************/
int hash_table_print(FILE *f, hash_table *t)
{
	int i, j, k, kmax;
	hash_node *cursor;

	if (t == NULL || f == NULL) {
		return HT_ERR_NULL_ARGUMENT;
	}

	fprintf(f, "\n****************************************\n");
        j = 0;
        kmax = 0;
	for (i = 0; i < t->size; i++) {
		cursor = t->table[i];
		fprintf(f, "table[%i]", i);
        	k = 0;
		while (cursor) {
			k++;
			j++;
			fprintf(f, " -> \"%s\"", cursor->key);
			cursor = cursor->next;
		}
		if (k > kmax) kmax = k;
		fprintf(f, "\n");
	}
	fprintf(f, "\n");
	fprintf(f, "Table size: %lu\n", t->size);
	fprintf(f, "# entries : %lu\n", j);
	fprintf(f, "Max list length : %lu\n", kmax);
	fprintf(f, "****************************************\n\n");

	return HT_SUCCESS;
}


/********************************************************************************
 *
 * hash_table_print_stdout
 *
 * One-line summary
 *
 * Full description
 *
 ********************************************************************************/
int hash_table_print_stdout(hash_table *t)
{
	return hash_table_print(stdout, t);
}


/********************************************************************************
 *
 * hash_table_destroy
 *
 * One-line summary
 *
 * Full description
 *
 ********************************************************************************/
int hash_table_destroy(hash_table *t)
{
	if (t == NULL) {
		return HT_ERR_NULL_ARGUMENT;
	}

	free(t->table);
	t->table = NULL;
	t->size = 0;

	return HT_SUCCESS;

}


/********************************************************************************
 *
 * hash_node_set_key
 *
 * One-line summary
 *
 * Full description
 *
 ********************************************************************************/
int hash_node_set_key(hash_node *n, const char *key)
{
	size_t keylen;

	if (n == NULL) {
		return HT_ERR_NULL_ARGUMENT;
	}

	keylen = strlen(key);

	if (keylen > KEYLEN_MAX) {
		return HT_ERR_KEY_TOO_LONG;
	}

	n->keylen = strlen(key);
	strncpy(n->key, key, KEYLEN_MAX);

	return HT_SUCCESS;
}


/********************************************************************************
 *
 * hash_node_set_data
 *
 * One-line summary
 *
 * Full description
 *
 ********************************************************************************/
int hash_node_set_data(hash_node *n, void *data)
{
	if (n == NULL) {
		return HT_ERR_NULL_ARGUMENT;
	}

	n->data = data;

	return HT_SUCCESS;
}


/********************************************************************************
 *
 * hash_node_get_key
 *
 * One-line summary
 *
 * Full description
 *
 ********************************************************************************/
const char * hash_node_get_key(hash_node *n)
{
	if (n == NULL) {
		return NULL;
	}

	return (const char *)n->key;
}


/********************************************************************************
 *
 * hash_node_get_data
 *
 * One-line summary
 *
 * Full description
 *
 ********************************************************************************/
void * hash_node_get_data(hash_node *n)
{
	if (n == NULL) {
		return NULL;
	}

	return n->data;
}
