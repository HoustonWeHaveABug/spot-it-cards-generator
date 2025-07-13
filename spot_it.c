#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#define PRIME_MIN 2
#define POWER_MIN 1

typedef struct {
	unsigned degree;
	unsigned *values;
}
polynomial_t;

unsigned prime_g, order_g, *inverses_g = NULL, *values_g = NULL, order_square_g, *mul_cache_g = NULL, *add_cache_g = NULL;
polynomial_t *pn_irreducible_g = NULL, **pn_ff_g = NULL;

static int is_prime(void);
static int search_irreducible_pn(unsigned);
static int pn_is_irreducible(void);
static polynomial_t *pn_gcd(const polynomial_t *, const polynomial_t *);
static int ff_values(unsigned, unsigned);
static unsigned *cache_creation(void);
static unsigned operation(unsigned, unsigned, unsigned *, polynomial_t *(*)(const polynomial_t *, const polynomial_t *));
static polynomial_t *pn_modulus(const polynomial_t *, const polynomial_t *);
static polynomial_t *pn_subtraction(const polynomial_t *, const polynomial_t *);
static polynomial_t *pn_multiplication(const polynomial_t *, const polynomial_t *);
static polynomial_t *pn_addition(const polynomial_t *, const polynomial_t *);
static void pn_min_max(const polynomial_t *, const polynomial_t *, const polynomial_t **, const polynomial_t **);
static unsigned pff_opposite(unsigned);
static unsigned ff_value(const unsigned *, unsigned);
static polynomial_t *pn_copy(const polynomial_t *);
static polynomial_t *pn_creation(unsigned);
static void pn_reduction(polynomial_t *);
static void log_message(const char *, ...);
static void cleanup(void);
static void pn_destruction(polynomial_t *);

int main(void) {
	int result;
	unsigned power, i;
	if (scanf("%u", &prime_g) != 1 || prime_g < PRIME_MIN) {
		log_message("Prime is not greater than or equal to %u\n", PRIME_MIN);
		return EXIT_FAILURE;
	}
	if (!is_prime()) {
		log_message("%u is not a prime number\n", prime_g);
		return EXIT_FAILURE;
	}
	if (scanf("%u", &power) != 1 || power < POWER_MIN) {
		log_message("Power is not greater than or equal to %u\n", POWER_MIN);
		return EXIT_FAILURE;
	}
	order_g = prime_g;
	for (i = POWER_MIN; i < power && order_g <= UINT_MAX/prime_g; ++i) {
		order_g *= prime_g;
	}
	if (order_g > UINT_MAX/prime_g) {
		log_message("Order is not lesser than or equal to %u\n", UINT_MAX);
		return EXIT_FAILURE;
	}
	inverses_g = malloc(sizeof(unsigned)*prime_g);
	if (!inverses_g) {
		log_message("main error: inverses_g = malloc(%u)\n", sizeof(unsigned)*prime_g);
		return EXIT_FAILURE;
	}
	for (i = 1; i < prime_g; ++i) {
		unsigned j;
		for (j = 1; j < prime_g && (i*j)%prime_g != 1; ++j);
		inverses_g[i] = j;
	}
	pn_irreducible_g = pn_creation(power);
	if (!pn_irreducible_g) {
		log_message("main error: pn_irreductible_g = pn_creation(%u)\n", power);
		cleanup();
		return EXIT_FAILURE;
	}
	result = search_irreducible_pn(0);
	if (result < 1) {
		if (!result) {
			log_message("Irreducible polynomial not found\n");
		}
		cleanup();
		return EXIT_FAILURE;
	}
	if (sizeof(polynomial_t *) > UINT_MAX/order_g) {
		log_message("Will not be able to allocate memory for pn_ff_g\n");
		cleanup();
		return EXIT_FAILURE;
	}
	pn_ff_g = malloc(sizeof(polynomial_t *)*order_g);
	if (!pn_ff_g) {
		log_message("main error: pn_ff_g = malloc(%u)\n", sizeof(polynomial_t *)*order_g);
		cleanup();
		return EXIT_FAILURE;
	}
	for (i = 0; i < order_g; ++i) {
		pn_ff_g[i] = NULL;
	}
	values_g = malloc(sizeof(unsigned)*power);
	if (!values_g) {
		log_message("main error: values_g = malloc(%u)\n", sizeof(unsigned *)*power);
		cleanup();
		return EXIT_FAILURE;
	}
	for (i = 0; i < power; ++i) {
		if (!ff_values(i, 0)) {
			free(values_g);
			cleanup();
			return EXIT_FAILURE;
		}
	}
	free(values_g);
	if (order_g*sizeof(unsigned) > UINT_MAX/order_g) {
		log_message("Will not be able to allocate memory for cache\n");
		cleanup();
		return EXIT_FAILURE;
	}
	order_square_g = order_g*order_g;
	mul_cache_g = cache_creation();
	if (!mul_cache_g) {
		log_message("main error: mul_cache_g = cache_creation()\n");
		cleanup();
		return EXIT_FAILURE;
	}
	add_cache_g = cache_creation();
	if (!add_cache_g) {
		log_message("main error: add_cache_g = cache_creation()\n");
		cleanup();
		return EXIT_FAILURE;
	}
	for (i = 0; i < order_g; ++i) {
		unsigned j;
		for (j = 0; j < order_g; ++j) {
			unsigned k;
			for (k = 0; k < order_g; ++k) {
				unsigned l = operation(i, k, mul_cache_g, pn_multiplication);
				if (l == order_g) {
					log_message("main error: l = operation(%u, %u, mul_cache_g, pn_multiplication)\n", i, k);
					cleanup();
					return EXIT_FAILURE;
				}
				l = operation(j, l, add_cache_g, pn_addition);
				if (l == order_g) {
					log_message("main error: l = operation(%u, %u, add_cache_g, pn_addition)\n", j, l);
					cleanup();
					return EXIT_FAILURE;
				}
				printf("%u ", l*order_g+k);
			}
			printf("%u\n", order_square_g+i);
		}
	}
	for (i = 0; i < order_g; ++i) {
		unsigned j;
		for (j = 0; j < order_g; ++j) {
			printf("%u ", j*order_g+i);
		}
		printf("%u\n", order_square_g+order_g);
	}
	for (i = 0; i < order_g; ++i) {
		printf("%u ", order_square_g+i);
	}
	printf("%u\n", order_square_g+order_g);
	fflush(stdout);
	cleanup();
	return EXIT_SUCCESS;
}

static int is_prime(void) {
	unsigned i;
	if (prime_g < 4) {
		return 1;
	}
	if (prime_g%2 == 0 || prime_g%3 == 0) {
		return 0;
	}
	for (i = 5; i*i <= prime_g; i += 6) {
		if (prime_g%i == 0 || prime_g%(i+2) == 0) {
			return 0;
		}
	}
	return 1;
}

static int search_irreducible_pn(unsigned degree) {
	int result;
	if (degree < pn_irreducible_g->degree) {
		unsigned i;
		result = 0;
		for (i = 0; i < prime_g && !result; ++i) {
			pn_irreducible_g->values[degree] = i;
			result = search_irreducible_pn(degree+1);
			if (result < 0) {
				log_message("search_irreducible_pn error: result = search_irreducible_pn(%u)\n", degree+1);
			}
		}
	}
	else {
		pn_irreducible_g->values[degree] = 1;
		result = pn_is_irreducible();
		if (result < 0) {
			log_message("search_irreducible_pn error: result = pn_is_irreducible()\n");
		}
	}
	return result;
}

static int pn_is_irreducible(void) {
	unsigned power = prime_g, degree = pn_irreducible_g->degree >> 1, i;
	for (i = 0; i < degree; ++i) {
		unsigned j;
		polynomial_t *pn_result = pn_creation(power), *pn_tmp;
		if (!pn_result) {
			log_message("pn_is_irreducible error: pn_result = pn_creation(%u)\n", power);
			return -1;
		}
		pn_result->values[0] = 0;
		pn_result->values[1] = pff_opposite(1);
		for (j = 2; j < pn_result->degree; ++j) {
			pn_result->values[j] = 0;
		}
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

static polynomial_t *pn_gcd(const polynomial_t *pn_a, const polynomial_t *pn_b) {
	polynomial_t *pn_result;
	if (pn_b->degree || pn_b->values[0]) {
		polynomial_t *pn_tmp = pn_modulus(pn_a, pn_b);
		if (!pn_tmp) {
			log_message("pn_gcd error: pn_tmp = pn_modulus(pn_a, pn_b)\n");
			return NULL;
		}
		pn_result = pn_gcd(pn_b, pn_tmp);
		pn_destruction(pn_tmp);
	}
	else {
		pn_result = pn_copy(pn_a);
		if (!pn_result) {
			log_message("pn_gcd error: pn_result = pn_copy(pn_a)\n");
		}
	}
	return pn_result;
}

static int ff_values(unsigned degree_max, unsigned degree) {
	unsigned i;
	if (degree < degree_max) {
		for (i = 0; i < prime_g; ++i) {
			values_g[degree] = i;
			if (!ff_values(degree_max, degree+1)) {
				log_message("ff_values error: result = ff_values(%u, %u)\n", degree_max, degree+1);
				return 0;
			}
		}
		return 1;
	}
	for (i = degree > 0; i < prime_g; ++i) {
		unsigned value, j;
		values_g[degree] = i;
		value = ff_value(values_g, degree);
		pn_ff_g[value] = pn_creation(degree_max);
		if (!pn_ff_g[value]) {
			log_message("ff_values error: pn_ff_g[%u] = pn_creation(%u)\n", value, degree_max);
			return 0;
		}
		for (j = 0; j <= degree_max; ++j) {
			pn_ff_g[value]->values[j] = values_g[j];
		}
	}
	return 1;
}

static unsigned *cache_creation(void) {
	unsigned *cache = malloc(sizeof(unsigned)*order_square_g), i;
	if (!cache) {
		log_message("cache_creation error: cache = malloc(%u)\n", sizeof(unsigned)*order_square_g);
		return NULL;
	}
	for (i = 0; i < order_square_g; ++i) {
		cache[i] = order_g;
	}
	return cache;
}

static unsigned operation(unsigned i, unsigned j, unsigned *cache, polynomial_t *(*pn_operation)(const polynomial_t *, const polynomial_t *)) {
	unsigned icache = i*order_g+j, result;
	polynomial_t *pn_result, *pn_tmp;
	if (cache[icache] < order_g) {
		return cache[icache];
	}
	pn_result = pn_operation(pn_ff_g[i], pn_ff_g[j]);
	if (!pn_result) {
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

static polynomial_t *pn_modulus(const polynomial_t *pn_a, const polynomial_t *pn_b) {
	polynomial_t *pn_result;
	if (!pn_b->degree && !pn_b->values[0]) {
		log_message("pn_modulus error: division by 0\n");
		return NULL;
	}
	pn_result = pn_copy(pn_a);
	if (!pn_result) {
		log_message("pn_modulus error: pn_result = pn_copy(pn_a)\n");
		return NULL;
	}
	while ((pn_result->degree || pn_result->values[0]) && pn_result->degree >= pn_b->degree) {
		unsigned i;
		polynomial_t *pn_high = pn_creation(pn_result->degree-pn_b->degree), *pn_tmp;
		if (!pn_high) {
			log_message("pn_modulus error: pn_high = pn_creation(%u)\n", pn_result->degree-pn_b->degree);
			pn_destruction(pn_result);
			return NULL;
		}
		for (i = 0; i < pn_high->degree; ++i) {
			pn_high->values[i] = 0;
		}
		pn_high->values[pn_high->degree] = pn_result->values[pn_result->degree]*inverses_g[pn_b->values[pn_b->degree]];
		pn_reduction(pn_high);
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

static polynomial_t *pn_subtraction(const polynomial_t *pn_a, const polynomial_t *pn_b) {
	unsigned i;
	const polynomial_t *pn_min, *pn_max;
	polynomial_t *pn_result;
	pn_min_max(pn_a, pn_b, &pn_min, &pn_max);
	pn_result = pn_creation(pn_max->degree);
	if (!pn_result) {
		log_message("pn_subtraction error: pn_result = pn_creation(%u)\n", pn_max->degree);
		return NULL;
	}
	for (i = 0; i <= pn_min->degree; ++i) {
		pn_result->values[i] = pn_a->values[i]+pff_opposite(pn_b->values[i]);
	}
	for (i = pn_min->degree+1; i <= pn_a->degree; ++i) {
		pn_result->values[i] = pn_a->values[i];
	}
	for (i = pn_min->degree+1; i <= pn_b->degree; ++i) {
		pn_result->values[i] = pff_opposite(pn_b->values[i]);
	}
	pn_reduction(pn_result);
	return pn_result;
}

static polynomial_t *pn_multiplication(const polynomial_t *pn_a, const polynomial_t *pn_b) {
	unsigned i;
	polynomial_t *pn_result = pn_creation(pn_a->degree+pn_b->degree);
	if (!pn_result) {
		log_message("pn_multiplication error: pn_result = pn_creation(%u)\n", pn_a->degree+pn_b->degree);
		return NULL;
	}
	for (i = 0; i <= pn_result->degree; ++i) {
		pn_result->values[i] = 0;
	}
	for (i = 0; i <= pn_a->degree; ++i) {
		if (pn_a->values[i]) {
			unsigned j;
			for (j = 0; j <= pn_b->degree; ++j) {
				if (pn_b->values[j]) {
					pn_result->values[i+j] += pn_a->values[i]*pn_b->values[j];
				}
			}
		}
	}
	pn_reduction(pn_result);
	return pn_result;
}

static polynomial_t *pn_addition(const polynomial_t *pn_a, const polynomial_t *pn_b) {
	unsigned i;
	const polynomial_t *pn_min, *pn_max;
	polynomial_t *pn_result;
	pn_min_max(pn_a, pn_b, &pn_min, &pn_max);
	pn_result = pn_creation(pn_max->degree);
	if (!pn_result) {
		log_message("pn_addition error: pn_result = pn_creation(%u)\n", pn_max->degree);
		return NULL;
	}
	for (i = 0; i <= pn_min->degree; ++i) {
		pn_result->values[i] = pn_a->values[i]+pn_b->values[i];
	}
	for (i = pn_min->degree+1; i <= pn_max->degree; ++i) {
		pn_result->values[i] = pn_max->values[i];
	}
	pn_reduction(pn_result);
	return pn_result;
}

static void pn_min_max(const polynomial_t *pn_a, const polynomial_t *pn_b, const polynomial_t **pn_min, const polynomial_t **pn_max) {
	if (pn_b->degree > pn_a->degree) {
		*pn_min = pn_a;
		*pn_max = pn_b;
	}
	else {
		*pn_min = pn_b;
		*pn_max = pn_a;
	}
}

static unsigned pff_opposite(unsigned value) {
	return prime_g-value;
}

static unsigned ff_value(const unsigned *values, unsigned degree) {
	unsigned factor = 1, result = values[0], i;
	for (i = 1; i <= degree; ++i) {
		factor *= prime_g;
		result += values[i]*factor;
	}
	return result;
}

static polynomial_t *pn_copy(const polynomial_t *pn) {
	unsigned i;
	polynomial_t *pn_result = pn_creation(pn->degree);
	if (!pn_result) {
		log_message("pn_copy error: pn_result = pn_creation(%u)\n", pn->degree);
		return NULL;
	}
	for (i = 0; i <= pn_result->degree; ++i) {
		pn_result->values[i] = pn->values[i];
	}
	return pn_result;
}

static polynomial_t *pn_creation(unsigned degree) {
	polynomial_t *pn_result = malloc(sizeof(polynomial_t));
	if (!pn_result) {
		log_message("pn_creation error: pn_result = malloc(%u)\n", sizeof(polynomial_t));
		return NULL;
	}
	pn_result->degree = degree;
	pn_result->values = malloc(sizeof(unsigned)*(degree+1));
	if (!pn_result->values) {
		log_message("pn_creation error: pn_result->values = malloc(%u)\n", sizeof(unsigned)*(degree+1));
		free(pn_result);
		return NULL;
	}
	return pn_result;
}

static void pn_reduction(polynomial_t *pn) {
	unsigned i;
	for (i = 0; i <= pn->degree; ++i) {
		pn->values[i] %= prime_g;
	}
	for (; pn->degree && !pn->values[pn->degree]; --pn->degree);
}

static void log_message(const char *format, ...) {
	va_list arguments;
	va_start(arguments, format);
	vfprintf(stderr, format, arguments);
	va_end(arguments);
	fflush(stderr);
}

static void cleanup(void) {
	if (add_cache_g) {
		free(add_cache_g);
	}
	if (mul_cache_g) {
		free(mul_cache_g);
	}
	if (pn_ff_g) {
		unsigned i;
		for (i = 0; i < order_g; ++i) {
			pn_destruction(pn_ff_g[i]);
		}
		free(pn_ff_g);
	}
	pn_destruction(pn_irreducible_g);
	free(inverses_g);
}

static void pn_destruction(polynomial_t *pn) {
	if (pn) {
		if (pn->values) {
			free(pn->values);
		}
		free(pn);
	}
}
