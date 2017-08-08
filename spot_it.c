#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#define FREE(p) free(p), p = NULL

#define PRIME_MIN 2
#define ORDER_MAX 4096
#define POWER_MIN 1

typedef struct {
	unsigned long degree;
	unsigned long *values;
}
polynomial_t;

unsigned long prime_g, order_g, *values_g = NULL, order_square_g, *mul_cache_g = NULL, *add_cache_g = NULL;
polynomial_t *pn_irreducible_g = NULL, **pn_ff_g = NULL;

int is_prime(void);
int search_irreducible_pn(unsigned long);
int pn_is_irreducible(void);
polynomial_t *pn_gcd(const polynomial_t *, const polynomial_t *);
int ff_values(unsigned long, unsigned long);
unsigned long *cache_creation(void);
unsigned long operation(unsigned long, unsigned long, unsigned long *, polynomial_t *(*)(const polynomial_t *, const polynomial_t *));
polynomial_t *pn_modulus(const polynomial_t *, const polynomial_t *);
polynomial_t *pn_subtraction(const polynomial_t *, const polynomial_t *);
polynomial_t *pn_multiplication(const polynomial_t *, const polynomial_t *);
polynomial_t *pn_addition(const polynomial_t *, const polynomial_t *);
void pn_min_max(const polynomial_t *, const polynomial_t *, const polynomial_t **, const polynomial_t **);
unsigned long pff_opposite(unsigned long);
unsigned long pff_inverse(unsigned long);
unsigned long ff_value(const unsigned long *, unsigned long);
void cleanup(void);
polynomial_t *pn_creation(unsigned long);
polynomial_t *pn_copy(const polynomial_t *);
int pn_reduction(polynomial_t *);
void pn_destruction(polynomial_t *);
void log_message(const char *, ...);

int main(void) {
int result;
unsigned long power, i, j, k, l;
	if (scanf("%lu", &prime_g) != 1 || prime_g < PRIME_MIN || prime_g > ORDER_MAX) {
		log_message("Prime is not between %lu and %lu\n", PRIME_MIN, ORDER_MAX);
		return EXIT_FAILURE;
	}
	if (!is_prime()) {
		log_message("%lu is not a prime number\n", prime_g);
		return EXIT_FAILURE;
	}
	if (scanf("%lu", &power) != 1 || power < POWER_MIN) {
		log_message("Power is not greater than or equal to %lu\n", POWER_MIN);
		return EXIT_FAILURE;
	}
	order_g = prime_g;
	for (i = 1; i < power && order_g <= ORDER_MAX; i++) {
		order_g *= prime_g;
	}
	if (order_g > ORDER_MAX) {
		log_message("Order is not lesser than or equal to %lu\n", ORDER_MAX);
		return EXIT_FAILURE;
	}
	pn_irreducible_g = pn_creation(power);
	if (!pn_irreducible_g) {
		log_message("main error: pn_irreductible_g = pn_creation(%lu)\n", power);
		return EXIT_FAILURE;
	}
	if (power > 1) {
		result = search_irreducible_pn(0);
		if (result < 1) {
			log_message(!result ? "Irreducible polynomial not found\n":"main error: result = search_irreducible_pn(0)\n");
			cleanup();
			return EXIT_FAILURE;
		}
	}
	else {
		for (i = 0; i < power; i++) {
			pn_irreducible_g->values[i] = 0;
		}
		pn_irreducible_g->values[power] = 1;
	}
	pn_ff_g = malloc(sizeof(polynomial_t *)*order_g);
	if (!pn_ff_g) {
		log_message("main error: pn_ff_g = malloc(%lu)\n", sizeof(polynomial_t *)*order_g);
		cleanup();
		return EXIT_FAILURE;
	}
	for (i = 0; i < order_g; i++) {
		pn_ff_g[i] = NULL;
	}
	for (i = 0; i < power; i++) {
		values_g = malloc(sizeof(unsigned long)*(i+1));
		if (!values_g) {
			log_message("main error: values_g = malloc(%lu)\n", sizeof(unsigned long *)*(i+1));
			cleanup();
			return EXIT_FAILURE;
		}
		if (!ff_values(i, 0)) {
			log_message("main error: ff_values(%lu, 0)\n", i);
			FREE(values_g);
			cleanup();
			return EXIT_FAILURE;
		}
		FREE(values_g);
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
	for (i = 0; i < order_g; i++) {
		for (j = 0; j < order_g; j++) {
			for (k = 0; k < order_g; k++) {
				l = operation(i, k, mul_cache_g, pn_multiplication);
				if (l == order_g) {
					log_message("main error: l = operation(%lu, %lu, mul_cache_g, pn_multiplication)\n", i, k);
					cleanup();
					return EXIT_FAILURE;
				}
				l = operation(j, l, add_cache_g, pn_addition);
				if (l == order_g) {
					log_message("main error: l = operation(%lu, %lu, add_cache_g, pn_addition)\n", j, l);
					cleanup();
					return EXIT_FAILURE;
				}
				if (k) {
					printf(" ");
				}
				printf("%lu", l*order_g+k);
			}
			printf(" %lu\n", order_square_g+i);
		}
	}
	for (i = 0; i < order_g; i++) {
		printf("%lu", i);
		for (j = 1; j < order_g; j++) {
			printf(" %lu", j*order_g+i);
		}
		printf(" %lu\n", order_square_g+order_g);
	}
	printf("%lu", order_square_g);
	for (i = 1; i <= order_g; i++) {
		printf(" %lu", order_square_g+i);
	}
	puts("");
	cleanup();
	return EXIT_SUCCESS;
}

int is_prime(void) {
unsigned long i;
	if (prime_g > 2 && prime_g%2 == 0) {
		return 0;
	}
	for (i = 3; i*i <= prime_g; i += 2) {
		if (prime_g%i == 0) {
			return 0;
		}
	}
	return 1;
}

int search_irreducible_pn(unsigned long degree) {
int result = 0;
unsigned long i;
	if (degree < pn_irreducible_g->degree) {
		i = degree == 0;
		while (i < prime_g && !result) {
			pn_irreducible_g->values[degree] = i;
			result = search_irreducible_pn(degree+1);
			if (result < 0) {
				log_message("search_irreducible_pn error: result = search_irreducible_pn(%lu)\n", degree+1);
			}
			i++;
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

int pn_is_irreducible(void) {
unsigned long power, degree, i, j;
polynomial_t *pn_result = NULL, *pn_tmp = NULL;
	power = prime_g;
	degree = pn_irreducible_g->degree >> 1;
	for (i = 0; i < degree; i++) {
		pn_result = pn_creation(power);
		if (!pn_result) {
			log_message("pn_is_irreducible error: pn_result = pn_creation(%lu)\n", power);
			return -1;
		}
		pn_result->values[0] = 0;
		pn_result->values[1] = pff_opposite(1);
		for (j = 2; j < pn_result->degree; j++) {
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
			log_message("pn_is_irreducible error: pn_result = pn_gcd(pn_irreducible_g, pn_tmp, %lu)\n");
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

polynomial_t *pn_gcd(const polynomial_t *pn_a, const polynomial_t *pn_b) {
polynomial_t *pn_result = NULL, *pn_tmp = NULL;
	if (pn_b->degree || pn_b->values[0]) {
		pn_tmp = pn_modulus(pn_a, pn_b);
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

int ff_values(unsigned long degree_max, unsigned long degree) {
int result = 1;
unsigned long i, value, j;
	if (degree < degree_max) {
		for (i = 0; i < prime_g && result; i++) {
			values_g[degree] = i;
			result = ff_values(degree_max, degree+1);
			if (!result) {
				log_message("ff_values error: result = ff_values(%lu, %lu)\n", degree_max, degree+1);
			}
		}
	}
	else {
		i = degree > 0;
		while (i < prime_g && result) {
			values_g[degree] = i;
			value = ff_value(values_g, degree);
			pn_ff_g[value] = pn_creation(degree_max);
			if (pn_ff_g[value]) {
				for (j = 0; j <= degree_max; j++) {
					pn_ff_g[value]->values[j] = values_g[j];
				}
			}
			else {
				log_message("ff_values error: pn_ff_g[%lu] = pn_creation(%lu)\n", value, degree_max);
				result = 0;
			}
			i++;
		}
	}
	return result;
}

unsigned long *cache_creation(void) {
unsigned long *cache = malloc(sizeof(unsigned long)*order_square_g), i;
	if (!cache) {
		log_message("cache_creation error: cache = malloc(%lu)\n", sizeof(unsigned long)*order_square_g);
		return NULL;
	}
	for (i = 0; i < order_square_g; i++) {
		cache[i] = order_g;
	}
	return cache;
}

unsigned long operation(unsigned long i, unsigned long j, unsigned long *cache, polynomial_t *(*pn_operation)(const polynomial_t *, const polynomial_t *)) {
unsigned long icache = i*order_g+j, result;
polynomial_t *pn_result = NULL, *pn_tmp = NULL;
	if (cache[icache] < order_g) {
		return cache[icache];
	}
	pn_result = pn_operation(pn_ff_g[i], pn_ff_g[j]);
	if (!pn_result) {
		log_message("operation error: pn_result = pn_operation(%lu, %lu)\n", pn_ff_g[i], pn_ff_g[j]);
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

polynomial_t *pn_modulus(const polynomial_t *pn_a, const polynomial_t *pn_b) {
unsigned long i;
polynomial_t *pn_result = NULL, *pn_high = NULL, *pn_tmp = NULL;
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
		pn_high = pn_creation(pn_result->degree-pn_b->degree);
		if (!pn_high) {
			log_message("pn_modulus error: pn_high = pn_creation(%lu)\n", pn_result->degree-pn_b->degree);
			pn_destruction(pn_result);
			return NULL;
		}
		for (i = 0; i < pn_high->degree; i++) {
			pn_high->values[i] = 0;
		}
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

polynomial_t *pn_subtraction(const polynomial_t *pn_a, const polynomial_t *pn_b) {
unsigned long i;
const polynomial_t *pn_min = NULL, *pn_max = NULL;
polynomial_t *pn_result = NULL;
	pn_min_max(pn_a, pn_b, &pn_min, &pn_max);
	pn_result = pn_creation(pn_max->degree);
	if (!pn_result) {
		log_message("pn_subtraction error: pn_result = pn_creation(%lu)\n", pn_max->degree);
		return NULL;
	}
	for (i = 0; i <= pn_min->degree; i++) {
		pn_result->values[i] = pn_a->values[i]+pff_opposite(pn_b->values[i]);
	}
	for (i = pn_min->degree+1; i <= pn_a->degree; i++) {
		pn_result->values[i] = pn_a->values[i];
	}
	for (i = pn_min->degree+1; i <= pn_b->degree; i++) {
		pn_result->values[i] = pff_opposite(pn_b->values[i]);
	}
	if (pn_reduction(pn_result) == EXIT_SUCCESS) {
		return pn_result;
	}
	else {
		log_message("pn_subtraction error: pn_reduction(pn_result)\n");
		return NULL;
	}
}

polynomial_t *pn_multiplication(const polynomial_t *pn_a, const polynomial_t *pn_b) {
unsigned long i, j;
polynomial_t *pn_result = pn_creation(pn_a->degree+pn_b->degree);
	if (!pn_result) {
		log_message("pn_multiplication error: pn_result = pn_creation(%lu)\n", pn_a->degree+pn_b->degree);
		return NULL;
	}
	for (i = 0; i <= pn_result->degree; i++) {
		pn_result->values[i] = 0;
	}
	for (i = 0; i <= pn_a->degree; i++) {
		if (pn_a->values[i]) {
			for (j = 0; j <= pn_b->degree; j++) {
				if (pn_b->values[j]) {
					pn_result->values[i+j] += pn_a->values[i]*pn_b->values[j];
				}
			}
		}
	}
	if (pn_reduction(pn_result) == EXIT_SUCCESS) {
		return pn_result;
	}
	else {
		log_message("pn_multiplication error: pn_reduction(pn_result)\n");
		return NULL;
	}
}

polynomial_t *pn_addition(const polynomial_t *pn_a, const polynomial_t *pn_b) {
unsigned long i;
const polynomial_t *pn_min = NULL, *pn_max = NULL;
polynomial_t *pn_result = NULL;
	pn_min_max(pn_a, pn_b, &pn_min, &pn_max);
	pn_result = pn_creation(pn_max->degree);
	if (!pn_result) {
		log_message("pn_addition error: pn_result = pn_creation(%lu)\n", pn_max->degree);
		return NULL;
	}
	for (i = 0; i <= pn_min->degree; i++) {
		pn_result->values[i] = pn_a->values[i]+pn_b->values[i];
	}
	for (i = pn_min->degree+1; i <= pn_max->degree; i++) {
		pn_result->values[i] = pn_max->values[i];
	}
	if (pn_reduction(pn_result) == EXIT_SUCCESS) {
		return pn_result;
	}
	else {
		log_message("pn_addition error: pn_reduction(pn_result)\n");
		return NULL;
	}
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

unsigned long pff_opposite(unsigned long value) {
	return prime_g-value;
}

unsigned long pff_inverse(unsigned long value) {
unsigned long i;
	for (i = 1; i < prime_g && (value*i)%prime_g != 1; i++);
	return i;
}

unsigned long ff_value(const unsigned long *values, unsigned long degree) {
unsigned long result, factor, i;
	result = 0;
	factor = 1;
	for (i = 0; i < degree; i++) {
		result += values[i]*factor;
		factor *= prime_g;
	}
	result += values[i]*factor;
	return result;
}

void cleanup(void) {
unsigned long i;
	if (add_cache_g) {
		FREE(add_cache_g);
	}
	if (mul_cache_g) {
		FREE(mul_cache_g);
	}
	if (pn_ff_g) {
		for (i = 0; i < order_g; i++) {
			pn_destruction(pn_ff_g[i]);
		}
		FREE(pn_ff_g);
	}
	pn_destruction(pn_irreducible_g);
}

polynomial_t *pn_creation(unsigned long degree) {
polynomial_t *pn_result = NULL;
	pn_result = malloc(sizeof(polynomial_t));
	if (!pn_result) {
		log_message("pn_creation error: pn_result = malloc(%lu)\n", sizeof(polynomial_t));
		return NULL;
	}
	pn_result->degree = degree;
	pn_result->values = malloc(sizeof(unsigned long)*(degree+1));
	if (!pn_result->values) {
		log_message("pn_creation error: pn_result->values = malloc(%lu)\n", sizeof(unsigned long)*(degree+1));
		FREE(pn_result);
		return NULL;
	}
	return pn_result;
}

polynomial_t *pn_copy(const polynomial_t *pn) {
unsigned long i;
polynomial_t *pn_result = pn_creation(pn->degree);
	if (!pn_result) {
		log_message("pn_copy error: pn_result = pn_creation(%lu)\n", pn->degree);
		return NULL;
	}
	for (i = 0; i <= pn_result->degree; i++) pn_result->values[i] = pn->values[i];
	return pn_result;
}

int pn_reduction(polynomial_t *pn) {
unsigned long i, *tmp = NULL;
	for (i = 0; i <= pn->degree; i++) {
		pn->values[i] %= prime_g;
	}
	for (i = pn->degree; i > 0 && !pn->values[i]; i--);
	if (i < pn->degree) {
		pn->degree = i;
		tmp = realloc(pn->values, sizeof(unsigned long)*(i+1));
		if (!tmp) {
			log_message("pn_reduction error: tmp = realloc(pn->values, %lu)\n", sizeof(unsigned long)*(i+1));
			pn_destruction(pn);
			return EXIT_FAILURE;
		}
		pn->values = tmp;
	}
	return EXIT_SUCCESS;
}

void pn_destruction(polynomial_t *pn) {
	if (pn) {
		if (pn->values) {
			FREE(pn->values);
		}
		FREE(pn);
	}
}

void log_message(const char *format, ...) {
va_list arguments;
	va_start(arguments, format);
	vfprintf(stderr, format, arguments);
	va_end(arguments);
	fflush(stderr);
}
