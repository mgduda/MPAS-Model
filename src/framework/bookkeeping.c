#include <stdio.h>
#include <stdlib.h>
#include "hash_table.h"


struct bookkeeper {
	hash_table ht;
};


int bookkeeping_init(struct bookkeeper **b)
{
	*b = malloc(sizeof(struct bookkeeper));
	hash_table_init(&(*b)->ht, (size_t)1000);

	return 0;
}


int bookkeeping_finalize(struct bookkeeper **b)
{
	hash_table_destroy(&(*b)->ht);
	free(*b);
	*b = NULL;

	return 0;
}


int bookkeeping_declare_requisite(struct bookkeeper *b, const char *fieldname)
{
	fprintf(stderr, "fieldname = %s\n", fieldname);
}
