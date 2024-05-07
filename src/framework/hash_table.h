#ifndef HASH_TABLE_H
#define HASH_TABLE_H

#include <stdio.h>

#define HT_SUCCESS           ( 0)
#define HT_ERR_NULL_ARGUMENT (-1)
#define HT_ERR_KEY_TOO_LONG  (-2)
#define HT_ERR_DUPLICATE_KEY (-3)

#define KEYLEN_MAX (size_t)32

typedef struct _hash_node {
	char key[KEYLEN_MAX];
	size_t keylen;
	struct _hash_node *next;
	void *data;
} hash_node;

typedef struct {
	size_t size;
	hash_node **table;
} hash_table;

int hash_table_init(hash_table *t, size_t n);
int hash_table_insert(hash_table *t, hash_node *node);
hash_node * hash_table_remove(hash_table *t, const char *key);
hash_node * hash_table_query(hash_table *t, const char *key);
int hash_table_print(FILE *f, hash_table *t);
int hash_table_print_stdout(hash_table *t);
int hash_table_destroy(hash_table *t);

int hash_node_set_key(hash_node *n, const char *key);
int hash_node_set_data(hash_node *n, void *data);
const char * hash_node_get_key(hash_node *n);
void * hash_node_get_data(hash_node *n);

#endif
