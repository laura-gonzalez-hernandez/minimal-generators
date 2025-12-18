import math
import time
import sys
import sympy as sy
from sympy import linear_eq_to_matrix

from minimal_generators.polynomials import build_monomials, build_f_polys
from minimal_generators.determinants import bin_det
from minimal_generators.utils import Colors, log_progress, log_info, log_debug

def solve_linear_system(eqs, symbols):
    """
    Solvers linear equations using an optimized Gauss-Jordan elimination approach, keeping symbolic computation.
    Method A.gauss_jordan_solve(b) solves Ax+b=0, hence we change the sign of original vector b.
    """
    exprs = list(eqs)
    A, b = linear_eq_to_matrix(exprs, list(symbols))
    sol_vec = A.gauss_jordan_solve(-b)[0]
  
    return [-sol_vec[i, 0] if hasattr(sol_vec, 'shape') and len(sol_vec.shape) > 1 else -sol_vec[i] 
            for i in range(len(symbols))]

def get_minimal_generators(n: int, singular: bool = False, emit: bool = True, progress: bool = False, debug: bool = False, lin_sys_solver: str = 'default'):
    """
    Generate minimal generators polynomials.
    
    Args:
        n: Positive integer (minimum 5)
        singular: Use singular-like format for output
        emit: Print polynomials to terminal
        progress: Show progress bar
        debug: Show debug messages
        lin_sys_solver: Linear system solver ('default' uses sy.solve, 'gauss-jordan' uses optimized solver)
    """
    a=int((n-1)*n/2)
    r=[]
    l=[] #l_r quotient
    e=[] #epsilon_r reminder
    phi1=[]
    phi2=[]
    phi3=[]

    #llistes en python comencen sempre pel 0 i no es poden inicialitzar al 1

    for k in range(1,n+1):
        r.append(a*(n-1)+1+k)
        l.append(int(r[k-1]/a))
        e.append(r[k-1]%a)
        phi1.append(l[k-1]-math.floor((e[k-1]+1)/2))
        phi2.append(e[k-1]-2*math.floor(e[k-1]/2))
        phi3.append(math.floor(e[k-1]/2))


    kappa = [min(a,c) for a,b,c in zip(phi1,phi2,phi3)]
    i = [b + abs(a-c) for a,b,c in zip(phi1,phi2,phi3)]
    c = [c-a for a,b,c in zip(phi1,phi2,phi3)]

    exp_x=[]
    exp_y=[]
    exp_z=[]

    for k in range(1,n+1):
        temp_x=[]
        temp_y=[]
        temp_z=[]
        for j in range(0,kappa[k-1]+1):
            temp_x.append(phi1[k-1]-j)
            temp_y.append(phi2[k-1]+2*j)
            temp_z.append(phi3[k-1]-j)
        exp_x.append(temp_x)
        exp_y.append(temp_y)
        exp_z.append(temp_z)

    interval_I=[]
    goods=[]
    coef=[]
    coef_h1,coef_h2,coef_h3,coef_h4=[],[],[],[]


    s1,s2,s3,s4 = [],[],[],[]
    phi_s1_1,phi_s1_2,phi_s1_3=[],[],[]
    phi_s2_1,phi_s2_2,phi_s2_3=[],[],[]
    phi_s3_1,phi_s3_2,phi_s3_3=[],[],[]
    phi_s4_1,phi_s4_2,phi_s4_3=[],[],[]

    exp_x_s1,exp_y_s1,exp_z_s1=[],[],[]
    exp_x_s2,exp_y_s2,exp_z_s2=[],[],[]
    exp_x_s3,exp_y_s3,exp_z_s3=[],[],[]
    exp_x_s4,exp_y_s4,exp_z_s4=[],[],[]
            
    if n%2 == 1:
        upper_k = n-3
    elif n%2 == 0:
        upper_k = n-4

    # Progress tracking
    tail_blocks = 3
    total_steps = max(0, upper_k) + tail_blocks
    processed = 0
    start_time = time.perf_counter() if progress else None
    def _progress(stage: str):
        if not progress: return
        bar_w = 40; pct = 0 if total_steps == 0 else processed / total_steps
        filled = int(pct * bar_w)
        bar = '█' * filled + '░' * (bar_w - filled)
        elapsed = time.perf_counter() - start_time
        log_progress(bar, processed, total_steps, pct, elapsed, stage)
    if progress:
        log_debug(f"Starting algorithm for n={n} (steps≈{total_steps})", enabled=debug)

    # Helper function to solve linear systems based on solver choice
    def _solve(eqs_or_matrix, symbols_or_b, num_vars=None):
        if lin_sys_solver == 'gauss-jordan':
            # Direct matrix solving: eqs_or_matrix is M, symbols_or_b is b
            # Solve M * x = b
            M = eqs_or_matrix
            b = symbols_or_b
            sol_vec = M.gauss_jordan_solve(b)[0]
            # Extract solution values - handle both 1D and 2D matrices
            if hasattr(sol_vec, 'shape'):
                if len(sol_vec.shape) == 1:
                    return [sol_vec[i] for i in range(sol_vec.shape[0])]
                else:
                    return [sol_vec[i, 0] for i in range(sol_vec.shape[0])]
            else:
                # If it's a list or other iterable
                return list(sol_vec)
        else:
            # Symbolic solving: eqs_or_matrix is equations, symbols_or_b is symbols
            sol = sy.solve(eqs_or_matrix, symbols_or_b)
            # Return list of values in the same order as symbols
            return [sol[sym] for sym in symbols_or_b]
            return [sol[sym] for sym in symbols_or_b]

    for k in range(1,upper_k+1):
        if k%2 == 1:
            interval_I.append(list(range(n-k-2,int((2*n-k-3)/2)+1))) #sumem +1 perquè és un range()
            goods.append([0]+list(range(int((n-k)/2),int((n-3)/2)+1)))
            temp_coef=[]
            for j in range(0,kappa[k-1]+1):
                temp_coef.append((-1)**(j)*bin_det(list(set(interval_I[k-1]) - {exp_x[k-1][j]}),goods[k-1])[1])
            coef.append(temp_coef)

            s1.append(r[k-1]+n)
            s2.append(r[k-1]+int((n-1)*n/2))
            s3.append(r[k-1]+int((n-1)*n/2)+n)
            s4.append(r[k-1]+n**2)

            phi_s1_1.append(math.floor((i[k-1]-1)/2))
            phi_s1_2.append(1)
            phi_s1_3.append(n-math.floor((i[k-1]+3)/2))

            phi_s2_1.append(int((n+i[k-1]+1)/2))
            phi_s2_2.append(0)
            phi_s2_3.append(int((n-i[k-1]-1)/2))

            phi_s3_1.append(math.floor((i[k-1]+1)/2))
            phi_s3_2.append(1)
            phi_s3_3.append(n-math.floor((i[k-1]+3)/2))

            phi_s4_1.append(math.floor((i[k-1]+3)/2))
            phi_s4_2.append(1)
            phi_s4_3.append(n-math.floor((i[k-1]+3)/2))

            kappa_s1 = [min(a,c) for a,b,c in zip(phi_s1_1,phi_s1_2,phi_s1_3)]
            kappa_s2 = [min(a,c) for a,b,c in zip(phi_s2_1,phi_s2_2,phi_s2_3)]
            kappa_s3 = [min(a,c) for a,b,c in zip(phi_s3_1,phi_s3_2,phi_s3_3)]   
            kappa_s4 = [min(a,c) for a,b,c in zip(phi_s4_1,phi_s4_2,phi_s4_3)]

            # define all exponents s1,..,s4
            temp_x=[]
            temp_y=[]
            temp_z=[]
            for j in range(0,kappa_s1[k-1]+1):
                temp_x.append(phi_s1_1[k-1]-j)
                temp_y.append(phi_s1_2[k-1]+2*j)
                temp_z.append(phi_s1_3[k-1]-j)
            exp_x_s1.append(temp_x)
            exp_y_s1.append(temp_y)
            exp_z_s1.append(temp_z)

            temp_x=[]
            temp_y=[]
            temp_z=[]
            for j in range(1,kappa_s2[k-1]+1):
                temp_x.append(phi_s2_1[k-1]-j)
                temp_y.append(phi_s2_2[k-1]+2*j)
                temp_z.append(phi_s2_3[k-1]-j)
            exp_x_s2.append(temp_x)
            exp_y_s2.append(temp_y)
            exp_z_s2.append(temp_z)

            temp_x=[]
            temp_y=[]
            temp_z=[]
            for j in range(0,kappa_s3[k-1]+1):
                temp_x.append(phi_s3_1[k-1]-j)
                temp_y.append(phi_s3_2[k-1]+2*j)
                temp_z.append(phi_s3_3[k-1]-j)
            exp_x_s3.append(temp_x)
            exp_y_s3.append(temp_y)
            exp_z_s3.append(temp_z)
            
            temp_x=[]
            temp_y=[]
            temp_z=[]
            for j in range(math.floor((i[k-1]+3)/2)-math.floor((i[k-1]-2)/2) , kappa_s4[k-1]+1):
                temp_x.append(phi_s4_1[k-1]-j)
                temp_y.append(phi_s4_2[k-1]+2*j)
                temp_z.append(phi_s4_3[k-1]-j)
            exp_x_s4.append(temp_x)
            exp_y_s4.append(temp_y)
            exp_z_s4.append(temp_z)

            # h1 coefficients
            A = sy.Matrix(bin_det(list(range(0, math.floor((i[k-1]-1)/2)+1)), list(range(0, math.floor((i[k-1]-1)/2)+1)))[0])
            C = sy.Matrix(bin_det(list(range(i[k-1], int((n-1+i[k-1])/2)+1)), list(range(1, math.floor((i[k-1]+1)/2)+1)))[0])
            o1 = sy.eye(math.floor((i[k-1]+1)/2))[:, ::-1]
            o2 = sy.eye(math.floor(int((n+1-i[k-1])/2)))[:, ::-1]
            lambda_vec = sy.Matrix(coef[k-1])
            
            if lin_sys_solver == 'gauss-jordan':
                # Direct matrix solving without creating symbols
                # Equation: A.T * o1 * lambda1_vec + C.T * o2 * lambda_vec = 0
                # Solve: (A.T * o1) * lambda1 = -(C.T * o2 * lambda_vec)
                M = A.T * o1
                b = -(C.T * o2 * lambda_vec)
                num_vars = math.floor((i[k-1]-1)/2)+1
                sol_values = _solve(M, b, num_vars)
            else:
                # Symbolic solving
                lambda1_syms = sy.symbols(f'lambda_0:{math.floor((i[k-1]-1)/2)+1}')
                lambda1_vec = sy.Matrix(lambda1_syms).reshape(len(lambda1_syms), 1)
                eqs = A.T * o1 * lambda1_vec + C.T * o2 * lambda_vec
                sol_values = _solve(eqs, lambda1_syms)

            coef_h1.append(sol_values)

            # h2 coefficients
            A = sy.Matrix(bin_det(list(range(i[k-1]+1, int((n+i[k-1]-1)/2)+1)), [0]+list(range(math.floor((i[k-1]+5)/2), math.floor(n/2)+1)))[0])
            C = sy.Matrix(bin_det(list(range(i[k-1], int((n-1+i[k-1])/2)+1)), [math.floor(n/2)])[0])
            Z = sy.zeros(C.rows,int((n-3-i[k-1])/2))
            C = C.row_join(Z)
            o1 = sy.eye(int((n-i[k-1]-1)/2))[:, ::-1]
            o2 = sy.eye(int((n+1-i[k-1])/2))[:, ::-1]
            
            if lin_sys_solver == 'gauss-jordan':
                M = A.T * o1
                b = -(C.T * o2 * lambda_vec)
                num_vars = math.floor((n-i[k-1]-1)/2)+1
                sol_values = _solve(M, b, num_vars)
            else:
                lambda2_syms = sy.symbols(f'lambda_1:{math.floor((n-i[k-1]-1)/2)+1}')
                lambda2_vec = sy.Matrix(lambda2_syms).reshape(len(lambda2_syms), 1)
                eqs = A.T * o1 * lambda2_vec + C.T * o2 * lambda_vec
                sol_values = _solve(eqs, lambda2_syms)

            coef_h2.append(sol_values)

            # h3 coefficients
            lambda2_vec = sy.Matrix(coef_h2[k-1])
            A = sy.Matrix(bin_det(list(range(0, math.floor((i[k-1]+1)/2)+1)), list(range(0, math.floor((i[k-1]+1)/2)+1)))[0])
            C1 = sy.Matrix(bin_det(list(range(i[k-1], int((n-1+i[k-1])/2)+1)), list(range(math.floor((n+2)/2),math.floor((n+2)/2)+math.floor((i[k-1]+1)/2)+1)))[0])
            C2 = sy.Matrix(bin_det(list(range(i[k-1]+1, int((n-1+i[k-1])/2)+1)), list(range(1,math.floor((i[k-1]+3)/2)+1)))[0])
            o = sy.eye(math.floor((i[k-1]+3)/2))[:, ::-1]
            o1 = sy.eye(int((n+1-i[k-1])/2))[:, ::-1] 
            o2 = sy.eye(int((n-1-i[k-1])/2))[:, ::-1]
            
            if lin_sys_solver == 'gauss-jordan':
                M = A.T * o
                b = -(C1.T * o1 * lambda_vec + C2.T * o2 * lambda2_vec)
                num_vars = math.floor((i[k-1]+1)/2)+1
                sol_values = _solve(M, b, num_vars)
            else:
                lambda3_syms = sy.symbols(f'lambda_0:{math.floor((i[k-1]+1)/2)+1}')
                lambda3_vec = sy.Matrix(lambda3_syms).reshape(len(lambda3_syms), 1)
                eqs = A.T * o * lambda3_vec + C1.T * o1 * lambda_vec + C2.T * o2 * lambda2_vec
                sol_values = _solve(eqs, lambda3_syms)

            coef_h3.append(sol_values)

            # h4 coefficients
            A = sy.Matrix(bin_det(list(range(0, math.floor((i[k-1]-2)/2)+1)), list(range(0, math.floor((i[k-1]-2)/2)+1)))[0])
            C = sy.Matrix(bin_det(list(range(i[k-1]+1, int((n-1+i[k-1])/2)+1)), list(range(math.floor((n+2)/2), int((n-1+i[k-1])/2)+1)))[0])
            o = sy.eye(math.floor((i[k-1])/2))[:, ::-1]
            o1 = sy.eye(int((n-1-i[k-1])/2))[:, ::-1]
            
            if lin_sys_solver == 'gauss-jordan':
                M = A.T * o
                b = -(C.T * o1 * lambda2_vec)
                num_vars = math.floor((i[k-1]+3)/2)+1 - (math.floor((i[k-1]+3)/2)-math.floor((i[k-1]-2)/2))
                sol_values = _solve(M, b, num_vars)
            else:
                lambda4_syms = sy.symbols(f'lambda_{math.floor((i[k-1]+3)/2)-math.floor((i[k-1]-2)/2)}:{math.floor((i[k-1]+3)/2)+1}')
                lambda4_vec = sy.Matrix(lambda4_syms).reshape(len(lambda4_syms), 1)
                eqs = A.T * o * lambda4_vec + C.T * o1 * lambda2_vec
                sol_values = _solve(eqs, lambda4_syms)

            coef_h4.append(sol_values)

        if k%2 == 0:

            interval_I.append(list(range(n-k-2,int((2*n-k-4)/2)+1))) #sumem +1 perquè és un range()
            goods.append([0]+list(range(int((n-k+1)/2),int((n-3)/2)+1)))
            temp_coef=[]
            for j in range(0,kappa[k-1]+1):
                temp_coef.append((-1)**(j)*bin_det(list(set(interval_I[k-1]) - {exp_x[k-1][j]}),goods[k-1])[1])
            coef.append(temp_coef)


            s1.append(r[k-1]+n)
            s2.append(r[k-1]+int((n-1)*n/2))
            s3.append(r[k-1]+int((n-1)*n/2)+n)
            s4.append(r[k-1]+n**2)

            phi_s1_1.append(math.floor((i[k-1]-1)/2))
            phi_s1_2.append(0)
            phi_s1_3.append(n-math.floor((i[k-1]+1)/2))

            phi_s2_1.append(int((n-1+i[k-1])/2))
            phi_s2_2.append(1)
            phi_s2_3.append(int((n-1-i[k-1])/2))

            phi_s3_1.append(math.floor((i[k-1]+1)/2))
            phi_s3_2.append(0)
            phi_s3_3.append(n-math.floor((i[k-1]+1)/2))

            phi_s4_1.append(math.floor((i[k-1]+3)/2))
            phi_s4_2.append(0)
            phi_s4_3.append(n-math.floor((i[k-1]+1)/2))

            kappa_s1 = [min(a,c) for a,b,c in zip(phi_s1_1,phi_s1_2,phi_s1_3)]
            kappa_s2 = [min(a,c) for a,b,c in zip(phi_s2_1,phi_s2_2,phi_s2_3)]
            kappa_s3 = [min(a,c) for a,b,c in zip(phi_s3_1,phi_s3_2,phi_s3_3)]   
            kappa_s4 = [min(a,c) for a,b,c in zip(phi_s4_1,phi_s4_2,phi_s4_3)]

            # define all exponents s1,..,s4
            temp_x=[]
            temp_y=[]
            temp_z=[]
            for j in range(0,kappa_s1[k-1]+1):
                temp_x.append(phi_s1_1[k-1]-j)
                temp_y.append(phi_s1_2[k-1]+2*j)
                temp_z.append(phi_s1_3[k-1]-j)
            exp_x_s1.append(temp_x)
            exp_y_s1.append(temp_y)
            exp_z_s1.append(temp_z)

            temp_x=[]
            temp_y=[]
            temp_z=[]
            for j in range(1,kappa_s2[k-1]+1):
                temp_x.append(phi_s2_1[k-1]-j)
                temp_y.append(phi_s2_2[k-1]+2*j)
                temp_z.append(phi_s2_3[k-1]-j)
            exp_x_s2.append(temp_x)
            exp_y_s2.append(temp_y)
            exp_z_s2.append(temp_z)

            temp_x=[]
            temp_y=[]
            temp_z=[]
            for j in range(0,kappa_s3[k-1]+1):
                temp_x.append(phi_s3_1[k-1]-j)
                temp_y.append(phi_s3_2[k-1]+2*j)
                temp_z.append(phi_s3_3[k-1]-j)
            exp_x_s3.append(temp_x)
            exp_y_s3.append(temp_y)
            exp_z_s3.append(temp_z)
            
            temp_x=[]
            temp_y=[]
            temp_z=[]
            for j in range(math.floor((i[k-1]+3)/2)-math.floor((i[k-1]-4)/2) , kappa_s4[k-1]+1):
                temp_x.append(phi_s4_1[k-1]-j)
                temp_y.append(phi_s4_2[k-1]+2*j)
                temp_z.append(phi_s4_3[k-1]-j)
            exp_x_s4.append(temp_x)
            exp_y_s4.append(temp_y)
            exp_z_s4.append(temp_z)

            # h1 coefficients
            A = sy.Matrix(bin_det(list(range(0, math.floor((i[k-1]-1)/2)+1)), list(range(0, math.floor((i[k-1]-1)/2)+1)))[0])
            C = sy.Matrix(bin_det(list(range(i[k-1]-1, int((n-3+i[k-1])/2)+1)), list(range(1, math.floor((i[k-1]+1)/2)+1)))[0])
            o1 = sy.eye(math.floor((i[k-1]+1)/2))[:, ::-1]
            o2 = sy.eye(math.floor(int((n+1-i[k-1])/2)))[:, ::-1]
            lambda_vec = sy.Matrix(coef[k-1])
            
            if lin_sys_solver == 'gauss-jordan':
                M = A.T * o1
                b = -(C.T * o2 * lambda_vec)
                num_vars = math.floor((i[k-1]-1)/2)+1
                sol_values = _solve(M, b, num_vars)
            else:
                lambda1_syms = sy.symbols(f'lambda_0:{math.floor((i[k-1]-1)/2)+1}')
                lambda1_vec = sy.Matrix(lambda1_syms).reshape(len(lambda1_syms), 1)
                eqs = A.T * o1 * lambda1_vec + C.T * o2 * lambda_vec
                sol_values = _solve(eqs, lambda1_syms)

            coef_h1.append(sol_values)

            # h2 coefficients
            A = sy.Matrix(bin_det(list(range(i[k-1], int((n-3+i[k-1])/2)+1)), [0]+list(range(math.floor((i[k-1]+5)/2), math.floor(n/2)+1)))[0])
            C = sy.Matrix(bin_det(list(range(i[k-1]-1, int((n-3+i[k-1])/2)+1)), [math.floor(n/2)])[0])
            Z = sy.zeros(C.rows,int((n-3-i[k-1])/2))
            C = C.row_join(Z)
            o1 = sy.eye(int((n-i[k-1]-1)/2))[:, ::-1]
            o2 = sy.eye(int((n+1-i[k-1])/2))[:, ::-1]
            
            if lin_sys_solver == 'gauss-jordan':
                M = A.T * o1
                b = -(C.T * o2 * lambda_vec)
                num_vars = math.floor((n-i[k-1]-1)/2)+1
                sol_values = _solve(M, b, num_vars)
            else:
                lambda2_syms = sy.symbols(f'lambda_1:{math.floor((n-i[k-1]-1)/2)+1}')
                lambda2_vec = sy.Matrix(lambda2_syms).reshape(len(lambda2_syms), 1)
                eqs = A.T * o1 * lambda2_vec + C.T * o2 * lambda_vec
                sol_values = _solve(eqs, lambda2_syms)

            coef_h2.append(sol_values)

            # h3 coefficients
            lambda2_vec = sy.Matrix(coef_h2[k-1])
            A = sy.Matrix(bin_det(list(range(0, math.floor((i[k-1]+1)/2)+1)), list(range(0, math.floor((i[k-1]+1)/2)+1)))[0])
            C1 = sy.Matrix(bin_det(list(range(i[k-1]-1, int((n-3+i[k-1])/2)+1)), list(range(math.floor((n+2)/2),math.floor((n+2)/2)+math.floor((i[k-1]+1)/2)+1)))[0])
            C2 = sy.Matrix(bin_det(list(range(i[k-1], int((n-3+i[k-1])/2)+1)), list(range(1,math.floor((i[k-1]+3)/2)+1)))[0])
            o = sy.eye(math.floor((i[k-1]+3)/2))[:, ::-1]
            o1 = sy.eye(int((n+1-i[k-1])/2))[:, ::-1] 
            o2 = sy.eye(int((n-1-i[k-1])/2))[:, ::-1]
            
            if lin_sys_solver == 'gauss-jordan':
                M = A.T * o
                b = -(C1.T * o1 * lambda_vec + C2.T * o2 * lambda2_vec)
                num_vars = math.floor((i[k-1]+1)/2)+1
                sol_values = _solve(M, b, num_vars)
            else:
                lambda3_syms = sy.symbols(f'lambda_0:{math.floor((i[k-1]+1)/2)+1}')
                lambda3_vec = sy.Matrix(lambda3_syms).reshape(len(lambda3_syms), 1)
                eqs = A.T * o * lambda3_vec + C1.T * o1 * lambda_vec + C2.T * o2 * lambda2_vec
                sol_values = _solve(eqs, lambda3_syms)

            coef_h3.append(sol_values)

            # h4 coefficients
            if n%2==1:
                i_min_h4=2
            elif n%2==0:
                i_min_h4=3
            if i[k-1]!=i_min_h4:
                A = sy.Matrix(bin_det(list(range(0, math.floor((i[k-1]-4)/2)+1)), list(range(0, math.floor((i[k-1]-4)/2)+1)))[0])
                C = sy.Matrix(bin_det(list(range(i[k-1], int((n-3+i[k-1])/2)+1)), list(range(math.floor((n+2)/2), int((n-3+i[k-1])/2)+1)))[0])
                o = sy.eye(math.floor((i[k-1]-2)/2))[:, ::-1]
                o1 = sy.eye(int((n-1-i[k-1])/2))[:, ::-1]
                
                if lin_sys_solver == 'gauss-jordan':
                    M = A.T * o
                    b = -(C.T * o1 * lambda2_vec)
                    num_vars = math.floor((i[k-1]+3)/2)+1 - (math.floor((i[k-1]+3)/2)-math.floor((i[k-1]-4)/2))
                    sol_values = _solve(M, b, num_vars)
                else:
                    lambda4_syms = sy.symbols(f'lambda_{math.floor((i[k-1]+3)/2)-math.floor((i[k-1]-4)/2)}:{math.floor((i[k-1]+3)/2)+1}')
                    lambda4_vec = sy.Matrix(lambda4_syms).reshape(len(lambda4_syms), 1)
                    eqs = A.T * o * lambda4_vec + C.T * o1 * lambda2_vec
                    sol_values = _solve(eqs, lambda4_syms)

                coef_h4.append(sol_values)
            
            if i[k-1]==i_min_h4:
                coef_h4.append([0])
        
        processed += 1; _progress(f"main k={k}")

    # casos especials
    if n%2 == 1:        
        k = n-2      
        interval_I.append(list(range(0,kappa[k-1]+1)))
        goods.append(list(range(0,kappa[k-1])))
        temp_coef=[]
        for j in range(0,kappa[k-1]+1):
            temp_coef.append((-1)**(j)*bin_det([int((n-1)/2)],[int((n-1-2*j)/2)])[1])
        coef.append(temp_coef)

        s1.append(r[k-1]+int((n-1)*n/2))
        s2.append(r[k-1]+int((n+1)*n/2))
        s3.append(0)
        s4.append(0)

        phi_s1_1.append(int((n+1)/2))
        phi_s1_2.append(0)
        phi_s1_3.append(int((n-1)/2))

        phi_s2_1.append(0)
        phi_s2_2.append(1)
        phi_s2_3.append(n-1)

        phi_s3_1.append(0)
        phi_s3_2.append(0)
        phi_s3_3.append(0)

        phi_s4_1.append(0)
        phi_s4_2.append(0)
        phi_s4_3.append(0)

        kappa_s1 = [min(a,c) for a,b,c in zip(phi_s1_1,phi_s1_2,phi_s1_3)]
        kappa_s2 = [min(a,c) for a,b,c in zip(phi_s2_1,phi_s2_2,phi_s2_3)]
        kappa_s3 = [min(a,c) for a,b,c in zip(phi_s2_1,phi_s2_2,phi_s2_3)]
        kappa_s4 = [min(a,c) for a,b,c in zip(phi_s2_1,phi_s2_2,phi_s2_3)]

        #### definim els exponents de s1,s2

        exp_x_s1.append([1])
        exp_y_s1.append([n-1])
        exp_z_s1.append([0])

        exp_x_s2.append([0])
        exp_y_s2.append([1])
        exp_z_s2.append([n-1])

        exp_x_s3.append([0])
        exp_y_s3.append([0])
        exp_z_s3.append([0])

        exp_x_s4.append([0])
        exp_y_s4.append([0])
        exp_z_s4.append([0])

        # h1 coefficients

        coef_h1.append([-1])


        # h2 coefficients

        coef_h2.append([1])

        # h3 h4 coef buids

        coef_h3.append([0])
        coef_h4.append([0])
        
        processed += 1; _progress("tail k=n-2")
            
        k = n-1

        interval_I.append(list(range(0,int((n-1)/2)+1))) #sumem +1 perquè és un range()
        goods.append(list(range(0,int((n-3)/2)+1)))
        temp_coef=[]
        for j in range(0,kappa[k-1]+1):
            temp_coef.append((-1)**(j)*bin_det([int((n-3)/2)],[int((n-3-2*j)/2)])[1])
        coef.append(temp_coef)

        s1.append(r[k-1]+int((n-3)*n/2))
        s2.append(r[k-1]+int((n-1)*n/2))
        s3.append(r[k-1]+(n-1)*n)
        s4.append(r[k-1]+n**2-1)

        phi_s1_1.append(n)
        phi_s1_2.append(0)
        phi_s1_3.append(0)

        phi_s2_1.append(int((n-1)/2))
        phi_s2_2.append(1)
        phi_s2_3.append(int((n-1)/2))

        phi_s3_1.append(int((n+1)/2))
        phi_s3_2.append(1)
        phi_s3_3.append(int((n-1)/2))

        phi_s4_1.append(1)
        phi_s4_2.append(0)
        phi_s4_3.append(n)

        kappa_s1 = [min(a,c) for a,b,c in zip(phi_s1_1,phi_s1_2,phi_s1_3)]
        kappa_s2 = [min(a,c) for a,b,c in zip(phi_s2_1,phi_s2_2,phi_s2_3)]
        kappa_s3 = [min(a,c) for a,b,c in zip(phi_s3_1,phi_s3_2,phi_s3_3)]   
        kappa_s4 = [min(a,c) for a,b,c in zip(phi_s4_1,phi_s4_2,phi_s4_3)]

        # define all exponents s1,..,s4

        exp_x_s1.append([n])
        exp_y_s1.append([0])
        exp_z_s1.append([0])

        temp_x=[]
        temp_y=[]
        temp_z=[]
        for j in range(1,kappa_s2[k-1]+1):
            temp_x.append(phi_s2_1[k-1]-j)
            temp_y.append(phi_s2_2[k-1]+2*j)
            temp_z.append(phi_s2_3[k-1]-j)
        exp_x_s2.append(temp_x)
        exp_y_s2.append(temp_y)
        exp_z_s2.append(temp_z)

        temp_x=[]
        temp_y=[]
        temp_z=[]
        for j in range(1,kappa_s3[k-1]+1):
            temp_x.append(phi_s3_1[k-1]-j)
            temp_y.append(phi_s3_2[k-1]+2*j)
            temp_z.append(phi_s3_3[k-1]-j)
        exp_x_s3.append(temp_x)
        exp_y_s3.append(temp_y)
        exp_z_s3.append(temp_z)


        exp_x_s4.append([0])
        exp_y_s4.append([2])
        exp_z_s4.append([n-1])

        # h1 coefficients

        coef_h1.append([-1])

        # h2 coefficients
        A = sy.Matrix(bin_det(list(range(0, int((n-3)/2)+1)), list(range(0, int((n-3)/2)+1)))[0])
        C = sy.Matrix(bin_det([n],list(range(1, int((n-1)/2)+1)))[0])
        o1 = sy.eye(int((n-1)/2))[:, ::-1]
        
        if lin_sys_solver == 'gauss-jordan':
            M = A.T * o1
            b = C.T
            sol_values = _solve(M, b)
        else:
            lambda2_syms = sy.symbols(f'lambda_1:{int((n-1)/2)+1}')
            lambda2_vec = sy.Matrix(lambda2_syms).reshape(len(lambda2_syms), 1)
            eqs = A.T * o1 * lambda2_vec - C.T
            sol_values = _solve(eqs, lambda2_syms)

        coef_h2.append(sol_values)

        # h3 coefficients
        A = sy.Matrix(bin_det(list(range(1, int((n-1)/2)+1)), [0]+list(range(2, int((n-1)/2)+1)))[0])
        C = sy.Matrix(bin_det([n],[int((n+1)/2)]+list(range(int((n+5)/2),n+1)))[0])
        o = sy.eye(int((n-1)/2))[:, ::-1]
        
        if lin_sys_solver == 'gauss-jordan':
            M = A.T * o
            b = C.T
            sol_values = _solve(M, b)
        else:
            lambda3_syms = sy.symbols(f'lambda_1:{int((n-1)/2)+1}')
            lambda3_vec = sy.Matrix(lambda3_syms).reshape(len(lambda3_syms), 1)
            eqs = A.T * o * lambda3_vec - C.T
            sol_values = _solve(eqs, lambda3_syms)

        coef_h3.append(sol_values)

        # h4 coefficients

        coef_h4.append([bin_det([n],[int((n+3)/2)])[1] - sum([bin_det([int((n+1-2*j)/2)],[1])[1] * coef_h3[k-1][j-1] for j in range(1,int((n-1)/2)+1)])])

        processed += 1; _progress("tail k=n-1")

        k = n

        interval_I.append(list(range(0,int((n-3)/2)+1))) #sumem +1 perquè és un range()
        goods.append(list(range(0,int((n-5)/2)+1)))
        temp_coef=[]
        for j in range(0,kappa[k-1]+1):
            temp_coef.append((-1)**(j)*bin_det([int((n-3)/2)],[int((n-3-2*j)/2)])[1])
        coef.append(temp_coef)

        s1.append(r[k-1]+int((n-3)*n/2))
        s2.append(r[k-1]+int((n-1)*n/2))
        s3.append(r[k-1]+(n-1)*n)

        phi_s1_1.append(n-1)
        phi_s1_2.append(1)
        phi_s1_3.append(0)

        phi_s2_1.append(int((n-1)/2))
        phi_s2_2.append(0)
        phi_s2_3.append(int((n+1)/2))

        phi_s3_1.append(int((n+1)/2))
        phi_s3_2.append(0)
        phi_s3_3.append(int((n+1)/2))

        phi_s4_1.append(0)
        phi_s4_2.append(0)
        phi_s4_3.append(0)

        kappa_s1 = [min(a,c) for a,b,c in zip(phi_s1_1,phi_s1_2,phi_s1_3)]
        kappa_s2 = [min(a,c) for a,b,c in zip(phi_s2_1,phi_s2_2,phi_s2_3)]
        kappa_s3 = [min(a,c) for a,b,c in zip(phi_s3_1,phi_s3_2,phi_s3_3)]  
        kappa_s4 = [min(a,c) for a,b,c in zip(phi_s3_1,phi_s3_2,phi_s3_3)]  

        # define all exponents s1,..,s4
        exp_x_s1.append([n-1])
        exp_y_s1.append([1])
        exp_z_s1.append([0])

        temp_x=[]
        temp_y=[]
        temp_z=[]
        for j in range(1,kappa_s2[k-1]+1):
            temp_x.append(phi_s2_1[k-1]-j)
            temp_y.append(phi_s2_2[k-1]+2*j)
            temp_z.append(phi_s2_3[k-1]-j)
        exp_x_s2.append(temp_x)
        exp_y_s2.append(temp_y)
        exp_z_s2.append(temp_z)

        temp_x=[]
        temp_y=[]
        temp_z=[]
        for j in range(2,kappa_s3[k-1]+1):
            temp_x.append(phi_s3_1[k-1]-j)
            temp_y.append(phi_s3_2[k-1]+2*j)
            temp_z.append(phi_s3_3[k-1]-j)
        exp_x_s3.append(temp_x)
        exp_y_s3.append(temp_y)
        exp_z_s3.append(temp_z)


        exp_x_s4.append([0])
        exp_y_s4.append([0])
        exp_z_s4.append([0])

        # h1 coefficients
        coef_h1.append([-1])

        # h2 coefficients
        A = sy.Matrix(bin_det(list(range(0, int((n-3)/2)+1)), list(range(0, int((n-3)/2)+1)))[0])
        C = sy.Matrix(bin_det([n-1],list(range(1, int((n-1)/2)+1)))[0])
        o1 = sy.eye(int((n-1)/2))[:, ::-1]
        
        if lin_sys_solver == 'gauss-jordan':
            M = A.T * o1
            b = C.T
            sol_values = _solve(M, b)
        else:
            lambda2_syms = sy.symbols(f'lambda_1:{int((n-1)/2)+1}')
            lambda2_vec = sy.Matrix(lambda2_syms).reshape(len(lambda2_syms), 1)
            eqs = A.T * o1 * lambda2_vec - C.T
            sol_values = _solve(eqs, lambda2_syms)

        coef_h2.append(sol_values)

        # h3 coefficients
        A = sy.Matrix(bin_det(list(range(0, int((n-3)/2)+1)),list(range(0, int((n-3)/2)+1)))[0])
        C = sy.Matrix(bin_det([n-1],list(range(int((n+1)/2),n)))[0])
        o = sy.eye(int((n-1)/2))[:, ::-1]
        
        if lin_sys_solver == 'gauss-jordan':
            M = A.T * o
            b = C.T
            sol_values = _solve(M, b)
        else:
            lambda3_syms = sy.symbols(f'lambda_2:{int((n+1)/2)+1}')
            lambda3_vec = sy.Matrix(lambda3_syms).reshape(len(lambda3_syms), 1)
            eqs = A.T * o * lambda3_vec - C.T
            sol_values = _solve(eqs, lambda3_syms)

        coef_h3.append(sol_values)

        # h4 buid

        coef_h4.append([0])
        
        processed += 1; _progress("tail k=n")

    if n%2 == 0:
        k = n-3      
        interval_I.append(list(range(1,int(n/2)+1)))
        goods.append([0]+list(range(2,int((n-2)/2)+1)))
        
        temp_coef=[]
        for j in range(0,kappa[k-1]+1):
            temp_coef.append((-1)**(j)*bin_det(list(set(interval_I[k-1]) - {exp_x[k-1][j]}),goods[k-1])[1])
        coef.append(temp_coef)

        s1.append(r[k-1]+n-1)
        s2.append(r[k-1]+int((n-1)*n/2))
        s3.append(r[k-1]+int((n+1)*n/2))
        s4.append(0)

        phi_s1_1.append(0)
        phi_s1_2.append(1)
        phi_s1_3.append(n-2)

        phi_s2_1.append(int((n+2)/2))
        phi_s2_2.append(0)
        phi_s2_3.append(int((n-2)/2))

        phi_s3_1.append(1)
        phi_s3_2.append(1)
        phi_s3_3.append(n-2)

        phi_s4_1.append(0)
        phi_s4_2.append(0)
        phi_s4_3.append(0)

        kappa_s1 = [min(a,c) for a,b,c in zip(phi_s1_1,phi_s1_2,phi_s1_3)]
        kappa_s2 = [min(a,c) for a,b,c in zip(phi_s2_1,phi_s2_2,phi_s2_3)]
        kappa_s3 = [min(a,c) for a,b,c in zip(phi_s2_1,phi_s2_2,phi_s2_3)]
        kappa_s4 = [min(a,c) for a,b,c in zip(phi_s2_1,phi_s2_2,phi_s2_3)]

    #### definim els exponents de s1,s2

        exp_x_s1.append([0])
        exp_y_s1.append([1])
        exp_z_s1.append([n-2])

        exp_x_s2.append([2])
        exp_y_s2.append([n-2])
        exp_z_s2.append([0])

        exp_x_s3.append([1,0])
        exp_y_s3.append([1,3])
        exp_z_s3.append([n-2,n-3])

        exp_x_s4.append([0])
        exp_y_s4.append([0])
        exp_z_s4.append([0])

        # h1 coefficients
        
        coef_h1.append([-sum([((-1)**j)*bin_det(list(set(interval_I[k-1]) - {exp_x[k-1][j]}),goods[k-1])[1]*bin_det([int((n-2*j)/2)],[1])[1] for j in range(0,int((n-2)/2)+1)])])


        # h2 coefficients

        coef_h2.append([-1])

        # h3 coefficients

        coef_h3.append([1,1])

        # h4 coef buids

        coef_h4.append([0])

        processed += 1; _progress("tail k=n-3")

        k=n-2
        temp_coef=[]
        for j in range(0,int((n-2)/2)+1):
            temp_coef.append((-1)**(j)*bin_det([int((n-2)/2)],[int((n-2-2*j)/2)])[1])
        coef.append(temp_coef)


        # La notació de s1s2s3s4 canvia respecre al Lemma 6.9. per facilitat

        s1.append(r[k-1]+int((n-1)*(n-2)/2))
        s2.append(r[k-1]+int((n-1)*n/2))
        s3.append(r[k-1]+n*(n-1))
        s4.append(r[k-1]+(n+1)*(n-1))

        phi_s1_1.append(n)
        phi_s1_2.append(0)
        phi_s1_3.append(0)

        phi_s2_1.append(int(n/2))
        phi_s2_2.append(1)
        phi_s2_3.append(int((n-2)/2))

        phi_s3_1.append(int((n+2)/2))
        phi_s3_2.append(1)
        phi_s3_3.append(int((n-2)/2))

        phi_s4_1.append(2)
        phi_s4_2.append(0)
        phi_s4_3.append(n-1)

        kappa_s1 = [min(a,c) for a,b,c in zip(phi_s1_1,phi_s1_2,phi_s1_3)]
        kappa_s2 = [min(a,c) for a,b,c in zip(phi_s2_1,phi_s2_2,phi_s2_3)]
        kappa_s3 = [min(a,c) for a,b,c in zip(phi_s2_1,phi_s2_2,phi_s2_3)]
        kappa_s4 = [min(a,c) for a,b,c in zip(phi_s2_1,phi_s2_2,phi_s2_3)]

    #### definim els exponents de s1,s2,s3,s4

        exp_x_s1.append([n])
        exp_y_s1.append([0])
        exp_z_s1.append([0])

        temp_x=[]
        temp_y=[]
        temp_z=[]
        for j in range(0,kappa_s2[k-1]+1):
            temp_x.append(phi_s2_1[k-1]-j)
            temp_y.append(phi_s2_2[k-1]+2*j)
            temp_z.append(phi_s2_3[k-1]-j)
        exp_x_s2.append(temp_x)
        exp_y_s2.append(temp_y)
        exp_z_s2.append(temp_z)


        temp_x=[]
        temp_y=[]
        temp_z=[]
        for j in range(2,kappa_s3[k-1]+1):
            temp_x.append(phi_s3_1[k-1]-j)
            temp_y.append(phi_s3_2[k-1]+2*j)
            temp_z.append(phi_s3_3[k-1]-j)
        exp_x_s3.append(temp_x)
        exp_y_s3.append(temp_y)
        exp_z_s3.append(temp_z)


        exp_x_s4.append([1,0])
        exp_y_s4.append([2,4])
        exp_z_s4.append([n-2,n-3])

        # h1 coefficients

        coef_h1.append([-1])


        # h2 coefficients
        A = sy.Matrix(bin_det(list(range(1, int(n/2)+1)), list(range(0, int((n-2)/2)+1)))[0])
        C = sy.Matrix(bin_det([n],list(range(1, int(n/2)+1)))[0])
        o1 = sy.eye(int(n/2))[:, ::-1]
        
        if lin_sys_solver == 'gauss-jordan':
            M = A.T * o1
            b = C.T
            sol_values = _solve(M, b)
        else:
            lambda2_syms = sy.symbols(f'lambda_0:{int((n-2)/2)+1}')
            lambda2_vec = sy.Matrix(lambda2_syms).reshape(len(lambda2_syms), 1)
            eqs = A.T * o1 * lambda2_vec - C.T
            sol_values = _solve(eqs, lambda2_syms)

        coef_h2.append(sol_values)

        # h3 coefficients
        A = sy.Matrix(bin_det(list(range(2, int((n-2)/2)+1)), [0]+list(range(3, int((n-2)/2)+1)))[0])
        C1= sy.Matrix([bin_det([n],[int((n+2)/2)])[1]-coef_h2[k-1][0]]) 
        C2 = sy.Matrix(bin_det([n],list(range(int((n+8)/2), n+1)))[0])
        C = C1.col_join(C2.T)
        o1 = sy.eye(int((n-4)/2))[:, ::-1]
        
        if lin_sys_solver == 'gauss-jordan':
            M = A.T * o1
            b = C
            sol_values = _solve(M, b)
        else:
            lambda3_syms = sy.symbols(f'lambda_2:{int((n-2)/2)+1}')
            lambda3_vec = sy.Matrix(lambda3_syms).reshape(len(lambda3_syms), 1)
            eqs = A.T * o1 * lambda3_vec - C
            sol_values = _solve(eqs, lambda3_syms)

        coef_h3.append(sol_values)

        # h4 coefficients
        aux_lambda3_vec = sy.Matrix(coef_h3[k-1])
        A = sy.Matrix(bin_det([0,1],[0,1])[0])
        C1= sy.Matrix(bin_det([n],[int((n+4)/2),int((n+6)/2)])[0])
        C2 = sy.Matrix(bin_det(list(range(2,int((n-2)/2)+1)),[1,2])[0])
        o1 = sy.eye(2)[:, ::-1]
        o2 = sy.eye(int((n-4)/2))[:, ::-1]
        
        if lin_sys_solver == 'gauss-jordan':
            M = A.T * o1
            b = C1.T - C2.T*o2*aux_lambda3_vec
            sol_values = _solve(M, b)
        else:
            lambda4_syms = sy.symbols(f'lambda_1:{2+1}')
            lambda4_vec = sy.Matrix(lambda4_syms).reshape(len(lambda4_syms), 1)
            eqs = A.T * o1 * lambda4_vec - C1.T + C2.T*o2*aux_lambda3_vec
            sol_values = _solve(eqs, lambda4_syms)

        coef_h4.append(sol_values)
            

        k = n-1

        interval_I.append(list(range(0,int((n-2)/2)+1))) #sumem +1 perquè és un range()
        temp_coef=[]
        for j in range(0,int((n-2)/2)+1):
            temp_coef.append((-1)**(j)*bin_det([int((n-2)/2)],[int((n-2-2*j)/2)])[1])
        coef.append(temp_coef)

        s1.append(r[k-1]+int((n-2)*(n-1)/2))
        s2.append(r[k-1]+int((n-1)*n/2))
        s3.append(r[k-1]+(n-1)*n)
        s4.append(r[k-1]+n*(n+1))

        phi_s1_1.append(n-1)
        phi_s1_2.append(1)
        phi_s1_3.append(0)

        phi_s2_1.append(int(n/2))
        phi_s2_2.append(0)
        phi_s2_3.append(int(n/2))

        phi_s3_1.append(int((n+2)/2))
        phi_s3_2.append(0)
        phi_s3_3.append(int(n/2))

        phi_s4_1.append(1)
        phi_s4_2.append(1)
        phi_s4_3.append(n-1)

        kappa_s1 = [min(a,c) for a,b,c in zip(phi_s1_1,phi_s1_2,phi_s1_3)]
        kappa_s2 = [min(a,c) for a,b,c in zip(phi_s2_1,phi_s2_2,phi_s2_3)]
        kappa_s3 = [min(a,c) for a,b,c in zip(phi_s3_1,phi_s3_2,phi_s3_3)]   
        kappa_s4 = [min(a,c) for a,b,c in zip(phi_s4_1,phi_s4_2,phi_s4_3)]

        # define all exponents s1,..,s4
        
        exp_x_s1.append([n-1])
        exp_y_s1.append([1])
        exp_z_s1.append([0])

        temp_x=[]
        temp_y=[]
        temp_z=[]
        for j in range(1,int(n/2)+1):
            temp_x.append(phi_s2_1[k-1]-j)
            temp_y.append(phi_s2_2[k-1]+2*j)
            temp_z.append(phi_s2_3[k-1]-j)
        exp_x_s2.append(temp_x)
        exp_y_s2.append(temp_y)
        exp_z_s2.append(temp_z)


        temp_x=[]
        temp_y=[]
        temp_z=[]
        for j in range(3,int(n/2)+1):
            temp_x.append(phi_s3_1[k-1]-j)
            temp_y.append(phi_s3_2[k-1]+2*j)
            temp_z.append(phi_s3_3[k-1]-j)
        exp_x_s3.append(temp_x)
        exp_y_s3.append(temp_y)
        exp_z_s3.append(temp_z)
        

        exp_x_s4.append([0])
        exp_y_s4.append([3])
        exp_z_s4.append([n-2])

        # h1 coefficients

        coef_h1.append([-1])

        # h2 coefficients
        A = sy.Matrix(bin_det(list(range(0, int((n-2)/2)+1)), list(range(0, int((n-2)/2)+1)))[0])
        C = sy.Matrix(bin_det([n-1],list(range(1, int(n/2)+1)))[0])
        o1 = sy.eye(int(n/2))[:, ::-1]
        
        if lin_sys_solver == 'gauss-jordan':
            M = A.T * o1
            b = C.T
            sol_values = _solve(M, b)
        else:
            lambda2_syms = sy.symbols(f'lambda_1:{int(n/2)+1}')
            lambda2_vec = sy.Matrix(lambda2_syms).reshape(len(lambda2_syms), 1)
            eqs = A.T * o1 * lambda2_vec - C.T
            sol_values = _solve(eqs, lambda2_syms)

        coef_h2.append(sol_values)
        
        # h3 coefficients
        A = sy.Matrix(bin_det(list(range(1, int((n-4)/2)+1)), [0]+list(range(2, int((n-4)/2)+1)))[0])
        C = sy.Matrix(bin_det([n-1],[int((n+2)/2)]+list(range(int((n+6)/2),n-1+1)))[0])
        o = sy.eye(int((n-4)/2))[:, ::-1]
        
        if lin_sys_solver == 'gauss-jordan':
            M = A.T * o
            b = C.T
            sol_values = _solve(M, b)
        else:
            lambda3_syms = sy.symbols(f'lambda_3:{int(n/2)+1}')
            lambda3_vec = sy.Matrix(lambda3_syms).reshape(len(lambda3_syms), 1)
            eqs = A.T * o * lambda3_vec - C.T
            sol_values = _solve(eqs, lambda3_syms)

        coef_h3.append(sol_values)

        # h4 coefficients
        coef_h4.append([bin_det([n-1],[int((n+4)/2)])[1] - sum([bin_det([int((n+2-2*j)/2)],[1])[1] * coef_h3[k-1][j-3] for j in range(3,int(n/2)+1)])]) ###duda aqui un altre cop del coef_h3[][j-1]
        
        processed += 1; _progress("tail k=n-2")
        
    g = build_monomials(coef, exp_x, exp_y, exp_z)
    h1 = build_monomials(coef_h1, exp_x_s1, exp_y_s1, exp_z_s1)
    h2 = build_monomials(coef_h2, exp_x_s2, exp_y_s2, exp_z_s2)
    h3 = build_monomials(coef_h3, exp_x_s3, exp_y_s3, exp_z_s3)
    h4 = build_monomials(coef_h4, exp_x_s4, exp_y_s4, exp_z_s4)

    f = build_f_polys(n,g,h1,h2,h3,h4)

    processed = total_steps; _progress("done")
    if progress:
        # Move to next line after progress bar
        print()
        log_debug(f"Completed algorithm for n={n}", enabled=debug)
        
    if emit:
        for idx, poly in enumerate(f):
            poly.print(singular=singular, name=f'f_{idx+1}')

    return f

# Placeholder to test new versions

def get_minimal_generators_new(n: int, singular: bool = False, emit: bool = True, progress: bool = False, debug: bool = False):
   
   return get_minimal_generators(n, singular=singular, emit=emit, progress=progress, debug=debug, lin_sys_solver='gauss-jordan')