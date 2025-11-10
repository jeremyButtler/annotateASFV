# Goal:

Describe how to use the genMath library. You probably
  should use the system libraries.

Files: genMath.h and genMath.c

# genMath:

These are my general math functions. Most are branchless
  functions from the stanfford bithacking guide. Has
  min, max, abloxlute values (SL for long, regular int).
  Also includes and ifMax/ifMin function to find the
  maximum/minimum and copy its paired input
  (is branchless).

The min/max ifMin/ifMax exist are .h file macros. The
  other functions are in the .c file.

## highest bit

To get the high bit in an unsigned integer use
  the `hibit2_genMath(unsigned int)` function. The return
  value is the highest bit.

`signed int log2SI = hibit2_genMath(100);`

## absolute values

For absolute values for signed ints use
  the `ab_genMath(signed int)` function. The returned
  value is the absolute value.

`signed int valueSI = ab_genMath(-10);`

For absolute values for signed longs use
  the `abSL_genMath(signed long)` function. The returned
  value is the absolute value.

`signed long valueSL = ab_genMath(-1000000000000);`

# maximum and minimums

For to get maximum use the `max_genMath(value1, value2)`
  function (macro). The returned value is the maximum. You
  can also use the `ifmax_genMath()` function to simulate
  an ifstatement for a number.

- ifmax Input:
  1. first value to compare
  2. second value to compare
  3. value to return if first value is greater
  4. value to return if second value is greater
- Output:
  - Returns input 3 if input 1 is greater then input 2
  - Returns input 4 if input 1 is less then input 2

**Example**

```
signed int firstGreater = ifmax_genMath(10, 1, 'A', 'x');
   /*will return A since 10 >= 1*/
```

The `min_genMath()` and `ifmin_genMath()` fucntions are
  like the `max_genMath()` and `ifmax_genMath()`
  functions, except they find the minimum instead of the
  maximum.

