# walk.c

Estimate the entropy from a random walker.

## Format of the output log file

Col    | Meaning
-------|-----------
1      | Time step 
2      | Uncorrected estimated entropy, mean
3      | Uncorrected estimated entropy, std.
4      | True value, which is log(q) for uniform distribution
5      | Linearly corrected entropy, 1st order, mean
6      | Linearly corrected entropy, 1st order, std.
7      | Linearly corrected entropy, 2nd order, mean
8      | Linearly corrected entropy, 2nd order, std.
9      | Exponentially corrected entropy, 1st order, mean
10     | Exponentially corrected entropy, 1st order, std.
11     | Exponentially corrected entropy, 2nd order, mean
12     | Exponentially corrected entropy, 2nd order, std.


# potts.c


The pair correction has to be nonpositive
S(i, j) <= S(i) + S(j)

The block average has to be less than the entire entropy
