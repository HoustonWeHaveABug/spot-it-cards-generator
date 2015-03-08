This program generates the set of n^2+n+1 cards to play the Spot It (aka Dobble) game for all n power of prime.

This set of cards has below properties that are those of a projective plane of order n:
- There are n^2+n+1 cards and n^2+n+1 symbols
- Each card contains n+1 symbols
- Each symbol appears on n+1 cards
- Every two cards have exactly one symbol in common
- For every two symbols there is exactly one card that has both of them

I would have never been able to write that program without the help of below article (thank you Maxime):

http://images.math.cnrs.fr/Dobble-et-la-geometrie-finie.html

The program takes as input two numbers: a prime (p) and a power (q), n = p^q.

The computations are performed using multiplication and addition in the finite field of order n (modulus an irreducible polynomial of degree q). Computation results are cached for better performance.

If q > 1 the irreducible polynomial that will be used in the computations is determined using the Ben-Or algorithm, otherwise it will be x.

Each element in the finite field is represented by a polynomial.

The result is presented as per below:
- One line per card
- For each card the n+1 symbols are represented by a number ranging from 0 to n^2+n
