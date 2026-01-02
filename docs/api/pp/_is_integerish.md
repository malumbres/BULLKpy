# Is integerish

```{eval-rst}
.. autofunction:: bullkpy.pp._is_integerish

```

Heuristically determine whether a matrix **looks like raw count data**, i.e. contains **non-negative integers** (up to numerical tolerance).

This is a **private utility function** used internally to distinguish raw counts from normalized or transformed expression matrices.

## Purpose

In many workflows (normalization, DE, GLMs), it is important to know whether an input matrix:  
- represents raw counts (integers ≥ 0), or
- has already been normalized / transformed (floats, negatives, decimals). 

_is_integerish provides a fast, heuristic check that works for both dense and sparse matrices without scanning the entire array.

## Heuristic definition

A matrix X is considered integer-like if:   
1. All checked values are ≥ 0.   
2. All checked values are close to integers, i.e.:  
```python
|x − round(x)| < tol
```
Only a subset of values is inspected for performance.

## Parameters

**X**   
Expression matrix to test.  
Accepted types:   
- NumPy array
- SciPy sparse matrix
- Array-like convertible to NumPy.  

Both dense and sparse inputs are supported.

**max_check**  
Maximum number of entries to inspect.   
- Default: 50_000   
- Prevents scanning very large matrices.  
- If X has more elements, only the first max_check values are used.   

**tol**   
Numerical tolerance for integer checking.  
- Default: 1e-8
- Allows for tiny floating-point deviations (e.g. 1.0000000001)

## Behavior by matrix type

**Sparse matrices**   
- Only non-zero entries (X.data) are inspected
- Efficient for large, sparse count matrices
- An empty sparse matrix (all zeros) returns True

**Dense matrices**   
- Matrix is flattened
- Only the first max_check values are tested
- An empty matrix returns True

#### Return value

**True**  
Matrix appears to contain only non-negative integers (within tolerance). 

**False**    
Any of the following detected:   
- Negative values
- Non-integer values (e.g. decimals)
- Normalized / log-transformed data

## Examples

Raw counts (integer)
```python
_is_integerish(np.array([[0, 1, 2], [3, 4, 5]]))
# True
```

Floating-point normalized data
```python
_is_integerish(np.array([[0.0, 1.2], [3.7, 4.1]]))
# False
```

Log-transformed data
```python
_is_integerish(np.log1p([[0, 10, 100]]))
# False
```

Sparse count matrix
```python
from scipy.sparse import csr_matrix
X = csr_matrix([[0, 1, 0], [2, 0, 3]])
_is_integerish(X)
# True
```

## Notes
- This is a heuristic, not a strict validator
- Intended for internal decision-making, not user-facing validation
- Designed to be:
	- Fast
	- Safe for large matrices
	- Robust to floating-point noise

## Typical internal use cases
- Decide whether to:
	- Run count-based normalization (CPM, size factors)
	- Fit NB-GLMs
	- Warn users about inappropriate input layers
- Distinguish between:
	- counts layers
	- log1p, CPM, or other transformed layers

## See also
	•	bk.pp.normalize_cpm
	•	bk.tl.de
	•	bk.tl.de_glm
	•	bk.pp.set_raw_counts

