
# COMMENTS

ViaPermutationAction.gi contains code that practically solves the following problem.

Input: G = <g_1, g_2, …, g_m>  <=  S_n

Output: filter/sieve F = [f_1, …, f_k] <= G

Problem: Suppose we are given F’ = [f_1’, …, f_k’] such that there exists an x in G with f_i^x = f_i’ for all 1 <= i <= k. Then we can compute x in G, i.e. x = RepresentativeAction(G, F, F’, OnTuples);

Practical Applications: Given any phi in Inn(G), we can efficiently compute x such that phi = Inn(x) by evaluating phi k times.


Idea:

Construct the sieve F in a heuristic approach by choosing random elements in G.

Check then if F satisfies certain properties:

- For each orbit O_i of F on [1,...,n] choose a point p_i and a transversal to reach the points O_i.
- For each orbit O_i:
  - For each point p in O_i:
    - Compute cycle signature: array containing the cycle length of the point p in the permutation f_j.
  - Does there exist a unique cycle signature among the points in O_i?
    - If no, then F is not a sieve.
- F is a sieve.

