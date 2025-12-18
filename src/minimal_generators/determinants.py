import sympy as sy

def build_pascal(n:int)->sy.Matrix:
    """
    Build the (n+1)x(n+1) Pascal matrix
    """
    B=[]
    B.append([1]+n*[0])
    for i in range(1,n+1):
        next_row=[1]
        for j in range(1,n+1):
            next_row.append(B[i-1][j] + B[i-1][j-1])
        B.append(next_row)
    return sy.Matrix(B)

_PASCAL_MATRIX=None
_PASCAL_SIZE=-1              

def get_pascal(n:int)->sy.Matrix:
    """
    Get the cached Pascal matrix of size (n+1)x(n+1), building it if necessary.
    """
    global _PASCAL_MATRIX,_PASCAL_SIZE
    if _PASCAL_MATRIX is None or _PASCAL_SIZE < n:
        _PASCAL_MATRIX=build_pascal(n)
        _PASCAL_SIZE=n
    return _PASCAL_MATRIX

def bin_det(rows:list[int], cols:list[int], use_pascal=True)->tuple[sy.Matrix,int]:
    """
    Compute the binomial determinant of two subsets. Return the matrix and its determinant.

    Uses Pascal matrix if use_pascal is True, else computes binomial coefficients directly.
    """
    m=len(rows); n=len(cols)
    if use_pascal:
        B=get_pascal(max(rows+cols) if rows+cols else 0)
    data=[]
    for rr in rows:
        row_vals=[]
        for cc in cols:
            row_vals.append(B[rr,cc] if use_pascal else sy.binomial(rr,cc))
        data.append(row_vals)
    M=sy.Matrix(data)
    return (M,M.det()) if m==n else (M,0)