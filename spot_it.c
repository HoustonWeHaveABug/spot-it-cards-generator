#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#define FREE(p) free(p), p = NULL

#define PRIME_MIN 2
#define ORDER_MAX 4096
#define POWER_MIN 1

struct polynomial_s {
	unsigned degree;
	unsigned *values;
};
typedef struct polynomial_s polynomial_t;

unsigned prime_g, order_g, *values_g = NULL, order_square_g, *mul_cache_g = NULL, *add_cache_g = NULL;
polynomial_t *pn_irreducible_g = NULL, **pn_ff_g = NULL;

void log_message(const char *, ...);
int is_prime(void);
unsigned pff_opposite(unsigned);
unsigned pff_inverse(unsigned);
polynomial_t *pn_creation(unsigned);
void pn_destruction(polynomial_t *);
int pn_reduction(polynomial_t *);
polynomial_t *pn_copy(const polynomial_t *);
void pn_min_max(const polynomial_t *, const polynomial_t *, const polynomial_t **, const polynomial_t **);
polynomial_t *pn_addition(const polynomial_t *, const polynomial_t *);
polynomial_t *pn_subtraction(const polynomial_t *, const polynomial_t *);
polynomial_t *pn_multiplication(const polynomial_t *, const polynomial_t *);
polynomial_t *pn_modulus(const polynomial_t *, const polynomial_t *);
polynomial_t *pn_gcd(const polynomial_t *, const polynomial_t *);
int pn_is_irreducible(void);
int search_irreducible_pn(unsigned);
unsigned ff_value(const unsigned *, unsigned);
int ff_values(unsigned, unsigned);
unsigned *cache_creation(void);
unsigned operation(unsigned, unsigned, unsigned *, polynomial_t *(*)(const polynomial_t *, const polynomial_t *));
void cleanup(void);

void log_message(const char *format, ...) {
va_list arguments;
	va_start(arguments, format);
	vfprintf(stderr, format, arguments);
	va_end(arguments);
	fflush(stderr);
}

int is_prime(void) {
unsigned i;
	if (prime_g > 2 && prime_g%2 == 0) return 0;
	for (i = 3; i*i <= prime_g; i += 2) if (prime_g%i == 0) return 0;
	return 1;
}

unsigned pff_opposite(unsigned value) {
	return prime_g-value;
}

unsigned pff_inverse(unsigned value) {
unsigned i;
	for (i = 1; i < prime_g && (value*i)%prime_g != 1; i++);
	return i;
}

polynomial_t *pn_creation(unsigned degree) {
polynomial_t *pn_result = NULL;
	if (!(pn_result = malloc(sizeof(polynomial_t)))) {
		log_message("pn_creation error: pn_result = malloc(%u)\n", sizeof(polynomial_t));
		return NULL;
	}
	pn_result->degree = degree;
	if (!(pn_result->values = malloc(sizeof(unsigned)*(degree+1)))) {
		log_message("pn_creation error: pn_result->values = malloc(%u)\n", sizeof(unsigned)*(degree+1));
		FREE(pn_result);
		return NULL;
	}
	return pn_result;
}

void pn_destruction(polynomial_t *pn) {
	if (pn) {
		if (pn->values) FREE(pn->values);
		FREE(pn);
	}
}

int pn_reduction(polynomial_t *pn) {
unsigned i, *tmp = NULL;
	for (i = 0; i <= pn->degree; i++) pn->values[i] %= prime_g;
	for (i = pn->degree; i > 0 && !pn->values[i]; i--);
	if (i < pn->degree) {
		pn->degree = i;
		if (!(tmp = realloc(pn->values, sizeof(unsigned)*(i+1)))) {
			log_message("pn_reduction error: tmp = realloc(pn->values, %u)\n", sizeof(unsigned)*(i+1));
			pn_destruction(pn);
			return EXIT_FAILURE;
		}
		pn->values = tmp;
	}
	return EXIT_SUCCESS;
}

polynomial_t *pn_copy(const polynomial_t *pn) {
unsigned i;
polynomial_t *pn_result = NULL;
	if (!(pn_result = pn_creation(pn->degree))) {
		log_message("pn_copy error: pn_result = pn_creation(%u)\n", pn->degree);
		return NULL;
	}
	for (i = 0; i <= pn_result->degree; i++) pn_result->values[i] = pn->values[i];
	return pn_result;
}

void pn_min_max(const polynomial_t *pn_a, const polynomial_t *pn_b, const polynomial_t **pn_min, const polynomial_t **pn_max) {
	if (pn_b->degree > pn_a->degree) {
		*pn_min = pn_a;
		*pn_max = pn_b;
	}
	else {
		*pn_min = pn_b;
		*pn_max = pn_a;
	}
}

polynomial_t *pn_addition(const polynomial_t *pn_a, const polynomial_t *pn_b) {
unsigned i;
const polynomial_t *pn_min = NULL, *pn_max = NULL;
polynomial_t *pn_result = NULL;
	pn_min_max(pn_a, pn_b, &pn_min, &pn_max);
	if (!(pn_result = pn_creation(pn_max->degree))) {
		log_message("pn_addition error: pn_result = pn_creation(%u)\n", pn_max->degree);
		return NULL;
	}
	for (i = 0; i <= pn_min->degree; i++) pn_result->values[i] = pn_a->values[i]+pn_b->values[i];
	for (i = pn_min->degree+1; i <= pn_max->degree; i++) pn_result->values[i] = pn_max->values[i];
	if (pn_reduction(pn_result) == EXIT_SUCCESS) return pn_result; else {
		log_message("pn_addition error: pn_reduction(pn_result)\n");
		return NULL;
	}
}

polynomial_t *pn_subtraction(const polynomial_t *pn_a, const polynomial_t *pn_b) {
unsigned i;
const polynomial_t *pn_min = NULL, *pn_max = NULL;
polynomial_t *pn_result = NULL;
	pn_min_max(pn_a, pn_b, &pn_min, &pn_max);
	if (!(pn_result = pn_creation(pn_max->degree))) {
		log_message("pn_subtraction error: pn_result = pn_creation(%u)\n", pn_max->degree);
		return NULL;
	}
	for (i = 0; i <= pn_min->degree; i++) pn_result->values[i] = pn_a->values[i]+pff_opposite(pn_b->values[i]);
	for (i = pn_min->degree+1; i <= pn_a->degree; i++) pn_result->values[i] = pn_a->values[i];
	for (i = pn_min->degree+1; i <= pn_b->degree; i++) pn_result->values[i] = pff_opposite(pn_b->values[i]);
	if (pn_reduction(pn_result) == EXIT_SUCCESS) return pn_result; else {
		log_message("pn_subtraction error: pn_reduction(pn_result)\n");
		return NULL;
	}
}

polynomial_t *pn_multiplication(const polynomial_t *pn_a, const polynomial_t *pn_b) {
unsigned i, j;
polynomial_t *pn_result = NULL;
	if (!(pn_result = pn_creation(pn_a->degree+pn_b->degree))) {
		log_message("pn_multiplication error: pn_result = pn_creation(%u)\n", pn_a->degree+pn_b->degree);
		return NULL;
	}
	for (i = 0; i <= pn_result->degree; i++) pn_result->values[i] = 0;
	for (i = 0; i <= pn_a->degree; i++) if (pn_a->values[i]) for (j = 0; j <= pn_b->degree; j++) if (pn_b->values[j]) pn_result->values[i+j] += pn_a->values[i]*pn_b->values[j];
	if (pn_reduction(pn_result) == EXIT_SUCCESS) return pn_result; else {
		log_message("pn_multiplication error: pn_reduction(pn_result)\n");
		return NULL;
	}
}

polynomial_t *pn_modulus(const polynomial_t *pn_a, const polynomial_t *pn_b) {
unsigned i;
polynomial_t *pn_result = NULL, *pn_high = NULL, *pn_tmp = NULL;
	if (!pn_b->degree && !pn_b->values[0]) {
		log_message("pn_modulus error: division by 0\n");
		return NULL;
	}
	if (!(pn_result = pn_copy(pn_a))) {
		log_message("pn_modulus error: pn_result = pn_copy(pn_a)\n");
		return NULL;
	}
	while ((pn_result->degree || pn_result->values[0]) && pn_result->degree >= pn_b->degree) {
		if (!(pn_high = pn_creation(pn_result->degree-pn_b->degree))) {
			log_message("pn_modulus error: pn_high = pn_creation(%u)\n", pn_result->degree-pn_b->degree);
			pn_destruction(pn_result);
			return NULL;
		}
		for (i = 0; i < pn_high->degree; i++) pn_high->values[i] = 0;
		pn_high->values[pn_high->degree] = pn_result->values[pn_result->degree]*pff_inverse(pn_b->values[pn_b->degree]);
		if (pn_reduction(pn_high) == EXIT_FAILURE) {
			log_message("pn_modulus error: pn_reduction(pn_high)\n");
			pn_destruction(pn_result);
			return NULL;
		}
		pn_tmp = pn_high;
		pn_high = pn_multiplication(pn_b, pn_tmp);
		pn_destruction(pn_tmp);
		if (!pn_high) {
			log_message("pn_modulus error: pn_high = pn_multiplication(pn_b, pn_tmp)\n");
			pn_destruction(pn_result);
			return NULL;
		}
		pn_tmp = pn_result;
		pn_result = pn_subtraction(pn_tmp, pn_high);
		pn_destruction(pn_tmp);
		pn_destruction(pn_high);
		if (!pn_result) {
			log_message("pn_modulus error: pn_result = pn_subtraction(pn_tmp, pn_high)\n");
			return NULL;
		}
	}
	return pn_result;
}

polynomial_t *pn_gcd(const polynomial_t *pn_a, const polynomial_t *pn_b) {
polynomial_t *pn_result = NULL, *pn_tmp = NULL;
	if (pn_b->degree || pn_b->values[0]) {
		if (!(pn_tmp = pn_modulus(pn_a, pn_b))) {
			log_message("pn_gcd error: pn_tmp = pn_modulus(pn_a, pn_b)\n");
			return NULL;
		}
		pn_result = pn_gcd(pn_b, pn_tmp);
		pn_destruction(pn_tmp);
	}
	else if (!(pn_result = pn_copy(pn_a))) log_message("pn_gcd error: pn_result = pn_copy(pn_a)\n");
	return pn_result;
}

int pn_is_irreducible(void) {
unsigned power, degree, i, j;
polynomial_t *pn_result = NULL, *pn_tmp = NULL;
	power = prime_g;
	degree = pn_irreducible_g->degree >> 1;
	for (i = 0; i < degree; i++) {
		if (!(pn_result = pn_creation(power))) {
			log_message("pn_is_irreducible error: pn_result = pn_creation(%u)\n", power);
			return -1;
		}
		pn_result->values[0] = 0;
		pn_result->values[1] = pff_opposite(1);
		for (j = 2; j < pn_result->degree; j++) pn_result->values[j] = 0;
		pn_result->values[pn_result->degree] = 1;
		pn_tmp = pn_result;
		pn_result = pn_modulus(pn_tmp, pn_irreducible_g);
		pn_destruction(pn_tmp);
		if (!pn_result) {
			log_message("pn_is_irreducible error: pn_result = pn_modulus(pn_tmp, pn_irreducible_g)\n");
			return -1;
		}
		pn_tmp = pn_result;
		pn_result = pn_gcd(pn_irreducible_g, pn_tmp);
		pn_destruction(pn_tmp);
		if (!pn_result) {
			log_message("pn_is_irreducible error: pn_result = pn_gcd(pn_irreducible_g, pn_tmp, %u)\n");
			return -1;
		}
		if (pn_result->degree) {
			pn_destruction(pn_result);
			return 0;
		}
		pn_destruction(pn_result);
		power *= prime_g;
	}
	return 1;
}

int search_irreducible_pn(unsigned degree) {
int result = 0;
unsigned i;
	if (degree < pn_irreducible_g->degree) {
		if (degree) i = 0; else i = 1;
		while (i < prime_g && !result) {
			pn_irreducible_g->values[degree] = i;
			if ((result = search_irreducible_pn(degree+1)) < 0) log_message("search_irreducible_pn error: result = search_irreducible_pn(%u)\n", degree+1);
			i++;
		}
	}
	else {
		pn_irreducible_g->values[degree] = 1;
		if ((result = pn_is_irreducible()) < 0) log_message("search_irreducible_pn error: result = pn_is_irreducible()\n");
	}
	return result;
}

unsigned ff_value(const unsigned *values, unsigned degree) {
unsigned result, factor, i;
	result = 0;
	factor = 1;
	for (i = 0; i < degree; i++) {
		result += values[i]*factor;
		factor *= prime_g;
	}
	result += values[i]*factor;
	return result;
}

int ff_values(unsigned degree_max, unsigned degree) {
int result = 1;
unsigned i, value, j;
	if (degree < degree_max)
		for (i = 0; i < prime_g && result; i++) {
			values_g[degree] = i;
			if (!(result = ff_values(degree_max, degree+1))) log_message("ff_values error: result = ff_values(%u, %u)\n", degree_max, degree+1);
		}
	else {
		if (degree) i = 1; else i = 0;
		while (i < prime_g && result) {
			values_g[degree] = i;
			value = ff_value(values_g, degree);
			if ((pn_ff_g[value] = pn_creation(degree_max))) for (j = 0; j <= degree_max; j++) pn_ff_g[value]->values[j] = values_g[j]; else {
				log_message("ff_values error: pn_ff_g[%u] = pn_creation(%u)\n", value, degree_max);
				result = 0;
			}
			i++;
		}
	}
	return result;
}

unsigned *cache_creation(void) {
unsigned *cache, i;
	if (!(cache = malloc(sizeof(unsigned)*order_square_g))) {
		log_message("cache_creation error: cache = malloc(%u)\n", sizeof(unsigned)*order_square_g);
		return NULL;
	}
	for (i = 0; i < order_square_g; i++) cache[i] = order_g;
	return cache;
}

unsigned operation(unsigned i, unsigned j, unsigned *cache, polynomial_t *(*pn_operation)(const polynomial_t *, const polynomial_t *)) {
unsigned icache, result;
polynomial_t *pn_result = NULL, *pn_tmp = NULL;
	icache = i*order_g+j;
	if (cache[icache] < order_g) return cache[icache];
	if (!(pn_result = pn_operation(pn_ff_g[i], pn_ff_g[j]))) {
		log_message("operation error: pn_result = pn_operation(%u, %u)\n", pn_ff_g[i], pn_ff_g[j]);
		return order_g;
	}
	pn_tmp = pn_result;
	pn_result = pn_modulus(pn_tmp, pn_irreducible_g);
	pn_destruction(pn_tmp);
	if (!pn_result) {
		log_message("operation error: pn_result = pn_modulus(pn_tmp, pn_irreducible_g)\n");
		return order_g;
	}
	result = ff_value(pn_result->values, pn_result->degree);
	pn_destruction(pn_result);
	cache[icache] = result;
	cache[j*order_g+i] = result;
	return result;
}

void cleanup(void) {
unsigned i;
	if (add_cache_g) FREE(add_cache_g);
	if (mul_cache_g) FREE(mul_cache_g);
	if (pn_ff_g) {
		for (i = 0; i < order_g; i++) pn_destruction(pn_ff_g[i]);
		FREE(pn_ff_g);
	}
	pn_destruction(pn_irreducible_g);
}

int main(void) {
int result;
unsigned power, i, j, k, l;
	fprintf(stdout, "Prime ? ");
	fflush(stdout);
	scanf("%u", &prime_g);
	if (prime_g < PRIME_MIN || prime_g > ORDER_MAX) {
		log_message("Prime is not between %u and %u\n", PRIME_MIN, ORDER_MAX);
		return EXIT_FAILURE;
	}
	if (!is_prime()) {
		log_message("%u is not a prime number\n", prime_g);
		return EXIT_FAILURE;
	}
	fprintf(stdout, "Power ? ");
	fflush(stdout);
	scanf("%u", &power);
	if (power < POWER_MIN) {
		log_message("Power is not greater than or equal to %u\n", POWER_MIN);
		return EXIT_FAILURE;
	}
	order_g = prime_g;
	for (i = 1; i < power && order_g <= ORDER_MAX; i++) order_g *= prime_g;
	if (order_g > ORDER_MAX) {
		log_message("Order is not lesser than or equal to %u\n", ORDER_MAX);
		return EXIT_FAILURE;
	}
	if (!(pn_irreducible_g = pn_creation(power))) {
		log_message("main error: pn_irreductible_g = pn_creation(%u)\n", power);
		return EXIT_FAILURE;
	}
	if (power > 1) {
		if ((result = search_irreducible_pn(0)) < 1) {
			log_message(!result ? "Irreducible polynomial not found\n":"main error: result = search_irreducible_pn(0)\n");
			cleanup();
			return EXIT_FAILURE;
		}
	}
	else {
		for (i = 0; i < power; i++) pn_irreducible_g->values[i] = 0;
		pn_irreducible_g->values[power] = 1;
	}
	if (!(pn_ff_g = malloc(sizeof(polynomial_t *)*order_g))) {
		log_message("main error: pn_ff_g = malloc(%u)\n", sizeof(polynomial_t *)*order_g);
		cleanup();
		return EXIT_FAILURE;
	}
	for (i = 0; i < order_g; i++) pn_ff_g[i] = NULL;
	for (i = 0; i < power; i++) {
		if (!(values_g = malloc(sizeof(unsigned)*(i+1)))) {
			log_message("main error: values_g = malloc(%u)\n", sizeof(unsigned *)*(i+1));
			cleanup();
			return EXIT_FAILURE;
		}
		if (!ff_values(i, 0)) {
			log_message("main error: ff_values(%u, 0)\n", i);
			FREE(values_g);
			cleanup();
			return EXIT_FAILURE;
		}
		FREE(values_g);
	}
	order_square_g = order_g*order_g;
	if (!(mul_cache_g = cache_creation())) {
		log_message("main error: mul_cache_g = cache_creation()\n");
		cleanup();
		return EXIT_FAILURE;
	}
	if (!(add_cache_g = cache_creation())) {
		log_message("main error: add_cache_g = cache_creation()\n");
		cleanup();
		return EXIT_FAILURE;
	}
	for (i = 0; i < order_g; i++)
		for (j = 0; j < order_g; j++) {
			for (k = 0; k < order_g; k++) {
				if ((l = operation(i, k, mul_cache_g, pn_multiplication)) == order_g) {
					log_message("main error: l = operation(%u, %u, mul_cache_g, pn_multiplication)\n", i, k);
					cleanup();
					return EXIT_FAILURE;
				}
				if ((l = operation(j, l, add_cache_g, pn_addition)) == order_g) {
					log_message("main error: l = operation(%u, %u, add_cache_g, pn_addition)\n", j, l);
					cleanup();
					return EXIT_FAILURE;
				}
				if (k) fprintf(stdout, " ");
				fprintf(stdout, "%u", l*order_g+k);
			}
			fprintf(stdout, " %u\n", order_square_g+i);
		}
	for (i = 0; i < order_g; i++) {
		fprintf(stdout, "%u", i);
		for (j = 1; j < order_g; j++) fprintf(stdout, " %u", j*order_g+i);
		fprintf(stdout, " %u\n", order_square_g+order_g);
	}
	fprintf(stdout, "%u", order_square_g);
	for (i = 1; i <= order_g; i++) fprintf(stdout, " %u", order_square_g+i);
	fprintf(stdout, "\n");
	fflush(stdout);
	cleanup();
	return EXIT_SUCCESS;
}
