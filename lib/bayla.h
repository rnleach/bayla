#ifndef _BAYLA_H_
#define _BAYLA_H_

#include <float.h>

#include "elk.h"
#include "magpie.h"
#include "packrat.h"

#pragma warning(push)

#define API static inline

/*---------------------------------------------------------------------------------------------------------------------------
 *                                                 Matrices and Vectors 
 *-------------------------------------------------------------------------------------------------------------------------*/
typedef struct
{
    i32 ndim;
    f64 *data;
} BayLaSquareMatrix;

#define MATRIX_IDX(row, col, ncols) ((row) * (ncols) + (col))

API BayLaSquareMatrix bayla_square_matrix_create(i32 ndim, MagAllocator *alloc);
API BayLaSquareMatrix bayla_square_matrix_invert(BayLaSquareMatrix matrix, MagAllocator *alloc, MagAllocator scratch);
API void bayla_square_matrix_solve(BayLaSquareMatrix matrix, f64 *b, MagAllocator scratch);
API void bayla_square_matrix_left_multiply(BayLaSquareMatrix matrix, f64 *input, f64 *output);
API f64 bayla_dot_product(size ndim, f64 *left, f64 *right);

typedef struct
{
    i32 ndim;
    b32 valid;
    f64 *data;
} BayLaCholeskyMatrix;

API BayLaCholeskyMatrix bayla_cholesky_decomp(BayLaSquareMatrix matrix, MagAllocator *alloc);
API void bayla_cholesky_multiply(BayLaCholeskyMatrix const chol, size const num_params, f64 const *restrict input, f64 *restrict output);

/* -------------------------------------------------------------------------------------------------------------------------
 *                                            USER PROVIDED: Model Functions
 * -----------------------------------------------------------------------------------------------------------------------*/

/* The logarithm of a normalized prior function for the parameters of this model.
 * 
 * P(Θ | Model, I) in Jaynes eq 20.1, 20.2 (page 602)
 */
typedef f64 (*BayLaLogPrior)(size num_parms, f64 const *parms, void *user_data);

/* The logarithm of a normalized likelihood function for this model and the data. 
 * 
 * The "data" is fixed and stored somewhere in the user_data. When you evaluate this function, you are evaluating the
 * likelihood of the data given a specific set of paramters to the model. Those parameters are the argument *parms to the
 * likelihood function.
 *
 * P(Data | Θ, Model, I) in Jaynes eq 20.1, 20.2 (pg 602)
 */
typedef f64 (*BayLaLogLikelihood)(size num_parms, f64 const *parms, void *user_data);

/* The location in parameter space of the max a posteriori.
 *
 * The "data" is fixed and stored somewhere in the user_data. When you evaluate this function, it fills the buffer parms_out
 * with the model parameters that maximize the posterior for this model and the data provided.
 */
typedef void (*BayLaMaxAPosteriori)(void *user_data, size num_parms, f64 *parms_out);

/* Calculates the Hessian matrix to within a constant factor at the max-a-posteriori parameter values.
 *
 * The "data" is fixed and stored somewhere in the ud. When you evaluate this function, it returns Hessian matrix (which is
 * sized num_parms x num_parms). This, along with the max_a_posteriori_parms are used to create a Gaussian approximation to
 * the posterior distribution. The Hessian is ony accurate to within a constant factor because we are lacking the value of 
 * the evidence (Z or z) which is the constant factor in the denominator of Bayes' Law. Ineed, this is one of the values we
 * are trying to estimate when we do the sampling.
 *
 * Note that this is the same as the negative of the inverse of the covariance matrix, and the result is used in 
 * bayla_mv_norm_dist_create_from_hessian.
 */
typedef BayLaSquareMatrix(*BayLaHessianMatrix)(void *ud, size nparms, f64 const *max_a_posteriori_parms, MagAllocator *alloc);

/*---------------------------------------------------------------------------------------------------------------------------
 *
 *                                                          API
 *
 *-------------------------------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------------
 *                                                       Math Tools
 *-------------------------------------------------------------------------------------------------------------------------*/

/* A type and functions for doing math in the log domain. */
typedef struct { f64 val; } BayLaLogValue;

API BayLaLogValue bayla_log_value_create(f64 log_domain_value);
API BayLaLogValue bayla_log_value_map_into_log_domain(f64 non_log_domain_value);
API f64 bayla_log_value_map_out_of_log_domain(BayLaLogValue log_val);
API BayLaLogValue bayla_log_value_add(BayLaLogValue left, BayLaLogValue right);
API BayLaLogValue bayla_log_value_subtract(BayLaLogValue left, BayLaLogValue right);
API BayLaLogValue bayla_log_value_multiply(BayLaLogValue left, BayLaLogValue right);
API BayLaLogValue bayla_log_value_divide(BayLaLogValue numerator, BayLaLogValue denominator);
API BayLaLogValue bayla_log_value_reciprocal(BayLaLogValue log_val);
API BayLaLogValue bayla_log_value_power(BayLaLogValue base, f64 exponent);
API BayLaLogValue bayla_log_values_sum(size nvals, size offset, size stride, BayLaLogValue *vals);

/* Functions that have less overflow / underflow / rounding errors than naive implementations.  */
API f64 bayla_log1mexp(f64 x);                                                 /* log(1 - exp(-x))                        */
API f64 bayla_log1pexp(f64 x);                                                 /* log(1 + exp(+x))                        */
API f64 bayla_logsumexp(f64 x, f64 y);                                         /* add two numbers in the log domain.      */
API f64 bayla_logsumexp_list(size nvals, size offset, size stride, f64 *vals); /* Sum a list of values in the log domain. */

BayLaLogValue const bayla_log_zero = { .val = -INFINITY };
BayLaLogValue const bayla_log_one = { .val = 0.0 };

/*---------------------------------------------------------------------------------------------------------------------------
 *                                                 Generate Random Numbers
 *-------------------------------------------------------------------------------------------------------------------------*/

typedef struct
{
    f64 mu;
    f64 sigma;
} BayLaNormalDistribution;

API BayLaNormalDistribution bayla_normal_distribution_create(f64 mean, f64 stddev);
API f64 bayla_normal_distribution_random_deviate(BayLaNormalDistribution *dist, ElkRandomState *state);

typedef struct
{
    i32 ndim;
    b32 valid;
    f64 *mean;
    BayLaCholeskyMatrix chol;
    BayLaSquareMatrix cov;
    BayLaSquareMatrix inv_cov;
    BayLaLogValue normalization_constant;
} BayLaMVNormDist;

API BayLaMVNormDist bayla_mv_norm_dist_create(size ndim, f64 *mean, BayLaSquareMatrix covar, MagAllocator *alloc);
API BayLaMVNormDist bayla_mv_norm_dist_create_from_hessian(size ndim, f64 *mean, BayLaSquareMatrix hessian, MagAllocator *alloc);

/* returns the probability value of the associated point, workspace needs to be at least as big as dist->ndim */
API BayLaLogValue bayla_mv_norm_dist_random_deviate(BayLaMVNormDist *dist, ElkRandomState *state, f64 *deviate, f64 *workspace);
/* returns the probability value of the input value deviate. workspace needs to be twice dist->ndim */
API BayLaLogValue bayla_mv_norm_dist_prob(BayLaMVNormDist *dist, f64 const *deviate, f64 *workspace);

/*---------------------------------------------------------------------------------------------------------------------------
 *                                              Analysis of Bayesian Models
 *-------------------------------------------------------------------------------------------------------------------------*/

/* Description of a model to be used for a Baysiean analysis. */
typedef struct
{
    
    f64 model_prior_probability;    /* Prior prob *this model* is the correct one. P(Model | I), Jaynes eq. 20.3 (pg 602) */
    size n_parameters;              /* The number of parameters in this model.                                            */

    BayLaLogPrior prior;                 /* See typedef above. P(Θ | Model, I), Jaynes eq 20.1, 20.2 (pg 602)             */
    BayLaLogLikelihood likelihood;       /* See typedef above. P(Data | Θ, Model, I), Jaynes eq 20.1, 20.2 (pg 602)       */
    BayLaMaxAPosteriori posterior_max;   /* See typedef above. argmax(Θ) P(Θ | Data, Model, I)                            */
    BayLaHessianMatrix hessian_matrix;   /* See typedef above. Used to create a Gaussian approximation.                   */

    void *user_data;                /* Any other state or data you may need for evaluating the prior, likelihood, etc.    */
} BayLaModel;

API BayLaLogValue bayla_model_evaluate(BayLaModel const *model, f64 const *parameter_vals);
/*---------------------------------------------------------------------------------------------------------------------------
 *                                                  Sampling and Statistics
 *-------------------------------------------------------------------------------------------------------------------------*/
typedef struct
{
    size n_samples;
    size ndim;
    f64 neff;
    BayLaLogValue z_evidence;
    f64 *rows;             /* n_samples x (ndim + 3) */
    b32 valid;
} BayLaSamples;

typedef struct
{
    f64 low;
    f64 high;
} BayLaCredibleInterval;

typedef struct 
{
    f64 val;
    f64 std;
} BayLaErrorValue;

/* sf is the spread factor, it means multiply the variance of the Gaussian proposal distribution to get more sampling
 * in the tails.
 */
API BayLaSamples bayla_importance_sample_gauss_approx(BayLaModel const *model, f64 sf, size n_samples, u64 seed, MagAllocator *alloc, MagAllocator scratch);
API void bayla_samples_save_csv(BayLaSamples *samples, f64 ci_p_thresh, char const *fname);
API f64 bayla_samples_calculate_ci_p_thresh(BayLaSamples *samples, f64 ci_pct, MagAllocator scratch);

API BayLaCredibleInterval bayla_samples_calculate_ci(BayLaSamples *samples, f64 ci_pct, size param_idx, MagAllocator scratch);
API BayLaErrorValue bayla_samples_estimate_evidence(BayLaSamples *samples);
API BayLaErrorValue bayla_samples_calculate_expectation(BayLaSamples *samples, f64 (*func)(size n, f64 const *params, void *ud), void *ud);

/* Since we can only know the Hessian to within a constant factor, this optomizes the scale factor to maximize the effective
 * sample size.
 */
API BayLaSamples bayla_importance_sample_gauss_approx_optimize(BayLaModel const *model, size n_samples, u64 seed, MagAllocator *alloc, MagAllocator scratch1, MagAllocator scratch2);
/*---------------------------------------------------------------------------------------------------------------------------
 *
 *
 *                                                    Implementations
 *
 *
 *-------------------------------------------------------------------------------------------------------------------------*/

/* Make sure these are always binary compatible. */
_Static_assert(sizeof(BayLaLogValue) == sizeof(f64));
_Static_assert(_Alignof(BayLaLogValue) == _Alignof(f64));

API BayLaSquareMatrix 
bayla_square_matrix_create(i32 ndim, MagAllocator *alloc)
{
    f64 *data = eco_arena_nmalloc(alloc, ndim * ndim, f64);
    Assert(data);
    return (BayLaSquareMatrix){ .ndim = ndim, .data = data };
}

API BayLaSquareMatrix 
bayla_square_matrix_invert(BayLaSquareMatrix matrix, MagAllocator *alloc, MagAllocator scratch)
{
    i32 const n = matrix.ndim;

    BayLaSquareMatrix inv = bayla_square_matrix_create(n, alloc);
    memcpy(inv.data, matrix.data, sizeof(f64) * n * n); 
    f64 *a = inv.data;

    i32 *indxc = eco_arena_nmalloc(&scratch, n, i32);
    i32 *indxr = eco_arena_nmalloc(&scratch, n, i32);
    i32 *ipiv = eco_arena_nmalloc(&scratch, n, i32);
    
    for(i32 i = 0; i < n; ++i) { Assert(ipiv[i] == 0); }

    i32 irow = 0;
    i32 icol = 0;
    for(i32 i = 0; i < n; ++i)
    {
        f64 big = 0.0;
        for(i32 j = 0; j < n; ++j)
        {
            if(ipiv[j] != 1)
            {
                for(i32 k = 0; k < n; ++k)
                {
                    if(ipiv[k] == 0)
                    {
                        if(fabs(a[MATRIX_IDX(j, k, n)]) >= big)
                        {
                            big = fabs(a[MATRIX_IDX(j, k, n)]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }

        ++(ipiv[icol]);

        if(irow != icol)
        {
            for(i32 l = 0; l < n; ++l)
            {
                f64 tmp = a[MATRIX_IDX(icol, l, n)];
                a[MATRIX_IDX(icol, l, n)] = a[MATRIX_IDX(irow, l, n)];
                a[MATRIX_IDX(irow, l, n)] = tmp;
            }
        }

        indxr[i] = irow;
        indxc[i] = icol;

        Assert(a[MATRIX_IDX(icol, icol, n)] != 0.0);

        f64 pivinv = 1.0 / a[MATRIX_IDX(icol, icol, n)];
        a[MATRIX_IDX(icol, icol, n)] = 1.0;
        for(i32 l = 0; l < n; ++l) { a[MATRIX_IDX(icol, l, n)] *= pivinv; }
        for(i32 ll = 0; ll < n; ++ll)
        {
            if(ll != icol)
            {
                f64 dum = a[MATRIX_IDX(ll, icol, n)];
                a[MATRIX_IDX(ll, icol, n)] = 0.0;
                for(i32 l = 0; l < n; ++l) { a[MATRIX_IDX(ll, l, n)] -= a[MATRIX_IDX(icol, l, n)] * dum; }
            }
        }
    }

    for(i32 l = n - 1; l >=0; --l)
    {
        if(indxr[l] != indxc[l])
        {
            for(i32 k = 0; k < n; ++k)
            {
                f64 tmp = a[MATRIX_IDX(k, indxr[l], n)];
                a[MATRIX_IDX(k, indxr[l], n)] = a[MATRIX_IDX(k, indxc[l], n)];
                a[MATRIX_IDX(k, indxc[l], n)] = tmp;
            }
        }
    }

    return inv;
}

API void 
bayla_square_matrix_solve(BayLaSquareMatrix matrix, f64 *b, MagAllocator scratch)
{
    i32 const n = matrix.ndim;

    BayLaSquareMatrix inv = bayla_square_matrix_create(n, &scratch);
    memcpy(inv.data, matrix.data, sizeof(f64) * n * n); 
    f64 *a = inv.data;

    i32 *indxc = eco_arena_nmalloc(&scratch, n, i32);
    i32 *indxr = eco_arena_nmalloc(&scratch, n, i32);
    i32 *ipiv = eco_arena_nmalloc(&scratch, n, i32);
    
    for(i32 i = 0; i < n; ++i) { Assert(ipiv[i] == 0); }

    i32 irow = 0;
    i32 icol = 0;
    for(i32 i = 0; i < n; ++i)
    {
        f64 big = 0.0;
        for(i32 j = 0; j < n; ++j)
        {
            if(ipiv[j] != 1)
            {
                for(i32 k = 0; k < n; ++k)
                {
                    if(ipiv[k] == 0)
                    {
                        if(fabs(a[MATRIX_IDX(j, k, n)]) >= big)
                        {
                            big = fabs(a[MATRIX_IDX(j, k, n)]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }

        ++(ipiv[icol]);

        if(irow != icol)
        {
            for(i32 l = 0; l < n; ++l)
            {
                f64 tmp = a[MATRIX_IDX(icol, l, n)];
                a[MATRIX_IDX(icol, l, n)] = a[MATRIX_IDX(irow, l, n)];
                a[MATRIX_IDX(irow, l, n)] = tmp;
            }
            f64 tmp = b[irow];
            b[irow] = b[icol];
            b[icol] = tmp;
        }

        indxr[i] = irow;
        indxc[i] = icol;

        Assert(a[MATRIX_IDX(icol, icol, n)] != 0.0);

        f64 pivinv = 1.0 / a[MATRIX_IDX(icol, icol, n)];
        a[MATRIX_IDX(icol, icol, n)] = 1.0;

        for(i32 l = 0; l < n; ++l) { a[MATRIX_IDX(icol, l, n)] *= pivinv; }
        b[icol] *= pivinv;

        for(i32 ll = 0; ll < n; ++ll)
        {
            if(ll != icol)
            {
                f64 dum = a[MATRIX_IDX(ll, icol, n)];
                a[MATRIX_IDX(ll, icol, n)] = 0.0;
                for(i32 l = 0; l < n; ++l) { a[MATRIX_IDX(ll, l, n)] -= a[MATRIX_IDX(icol, l, n)] * dum; }
                b[ll] -= b[icol] * dum;
            }
        }
    }

    for(i32 l = n - 1; l >=0; --l)
    {
        if(indxr[l] != indxc[l])
        {
            for(i32 k = 0; k < n; ++k)
            {
                f64 tmp = a[MATRIX_IDX(k, indxr[l], n)];
                a[MATRIX_IDX(k, indxr[l], n)] = a[MATRIX_IDX(k, indxc[l], n)];
                a[MATRIX_IDX(k, indxc[l], n)] = tmp;
            }
        }
    }
}

API void 
bayla_square_matrix_left_multiply(BayLaSquareMatrix matrix, f64 *input, f64 *output)
{
    for(i32 i = 0; i < matrix.ndim; ++i)
    {
        output[i] = input[0] * matrix.data[MATRIX_IDX(i, 0, matrix.ndim)];
        for(i32 j = 1; j < matrix.ndim; ++j)
        {
            output[i] += input[j] * matrix.data[MATRIX_IDX(i, j, matrix.ndim)];
        }
    }
}

API f64 
bayla_dot_product(size ndim, f64 *left, f64 *right)
{
    f64 val = 0.0;
    for(i32 i = 0; i < ndim; ++i)
    {
        val += left[i] * right[i];
    }

    return val;
}

API BayLaCholeskyMatrix 
bayla_cholesky_decomp(BayLaSquareMatrix matrix, MagAllocator *alloc)
{
    i32 ndim = matrix.ndim;
    f64 *data = eco_arena_nmalloc(alloc, (size)(ndim * ndim), f64);
    BayLaCholeskyMatrix chol = { .ndim = ndim, .data = data, .valid = true };
    memcpy(chol.data, matrix.data, ndim * ndim *sizeof(f64));


    for(i32 i = 0; i < ndim; ++i)
    {
        for(i32 j = i; j < ndim; ++j)
        {
            f64 sum = chol.data[MATRIX_IDX(i, j, ndim)];
            for(i32 k = i - 1; k >= 0; --k)
            {
                sum -= chol.data[MATRIX_IDX(i, k, ndim)] * chol.data[MATRIX_IDX(j, k, ndim)];
            }

            if(i == j)
            {
                if(sum <= 0.0)
                {
                    chol.valid = false; 
                    return chol;
                }
                chol.data[MATRIX_IDX(i, i, ndim)] = sqrt(sum);
            }
            else
            {
                chol.data[MATRIX_IDX(j, i, ndim)] = sum / chol.data[MATRIX_IDX(i, i, ndim)];
            }
        }
    }

    for(i32 i = 0; i < ndim; ++i)
    {
        for(i32 j = 0; j < i; ++j)
        { 
            chol.data[MATRIX_IDX(j, i, ndim)] = 0.0;
        }
    }

    return chol;
}

API void 
bayla_cholesky_multiply(BayLaCholeskyMatrix const chol, size const num_params, f64 const *restrict input, f64 *restrict output)
{
    Assert(num_params == chol.ndim);

    for(i32 i = 0; i < num_params; ++i)
    {
        output[i] = 0.0;
        for(i32 j = 0; j <= i; ++j)
        {
            output[i] += chol.data[MATRIX_IDX(i, j, chol.ndim)] * input[j];
        }
    }
}

API BayLaLogValue
bayla_log_value_create(f64 log_domain_value)
{
    return (BayLaLogValue){.val = log_domain_value };
}

API f64 
bayla_log_value_map_out_of_log_domain(BayLaLogValue log_val) 
{ 
    return exp(log_val.val); 
}

API BayLaLogValue 
bayla_log_value_map_into_log_domain(f64 non_log_domain_value)
{
    return (BayLaLogValue){.val = log(non_log_domain_value) };
}

API BayLaLogValue 
bayla_log_value_add(BayLaLogValue left, BayLaLogValue right)
{
    return (BayLaLogValue){.val = bayla_logsumexp(left.val, right.val)};
}

API BayLaLogValue 
bayla_log_value_subtract(BayLaLogValue left, BayLaLogValue right)
{
    f64 l = left.val;
    f64 r = right.val;

    if(isinf(l) && isinf(r) && r > 0.0) { return (BayLaLogValue){ .val = NAN }; }
    if(l == r) { return (BayLaLogValue){ .val = - INFINITY}; }

    f64 val = l + bayla_log1mexp(l - r);
    return (BayLaLogValue){.val = val };
}

API BayLaLogValue
bayla_log_value_multiply(BayLaLogValue left, BayLaLogValue right)
{
    return (BayLaLogValue){.val = left.val + right.val};
}

API BayLaLogValue 
bayla_log_value_divide(BayLaLogValue numerator, BayLaLogValue denominator)
{
    return (BayLaLogValue){.val = numerator.val - denominator.val};
}

API BayLaLogValue 
bayla_log_value_reciprocal(BayLaLogValue log_val)
{
    return (BayLaLogValue){.val = -log_val.val};
}

API BayLaLogValue 
bayla_log_value_power(BayLaLogValue base, f64 exponent)
{
    return (BayLaLogValue){.val = exponent * base.val};
}

/* Make sure the cast in the function below is valid. */
_Static_assert(sizeof(BayLaLogValue) == sizeof(f64), "BayLaLogValue not wrapper for f64 (double)");

API BayLaLogValue 
bayla_log_values_sum(size nvals, size offset, size stride, BayLaLogValue *vals)
{
    f64 *log_vals = (f64 *)vals;
    f64 log_sum = bayla_logsumexp_list(nvals, offset, stride, log_vals);
    return (BayLaLogValue){ .val = log_sum };
}

API f64 
bayla_log1mexp(f64 x)  /* log(1 - exp(-x)) */
{
    /* Algorithm copied from eq 7 in https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf 
     * on Sept 9, 2024.
     */
    if(x >= 0.0 && x <= ELK_LN2) { return log(-expm1(-x)); }
    else if(x > ELK_LN2) { return log1p(-exp(-x)); }
    else { return log(-1.0); /* This is a NAN and it should set the proper math flags. */ }
}

API f64
bayla_log1pexp(f64 x)  /* log(1 + exp(+x)) */
{
    /* Algorithm copied from eq 10 in https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf 
     * on Sept 9, 2024.
     */
    if(x < -37.0)     { return exp(x);        }
    else if(x < 18.0) { return log1p(exp(x)); }
    else if(x < 33.3) { return x + exp(-x);   }
    else              { return x; /* This represents a lack of enough precision to do the calculation. */ }
}

API f64 
bayla_logsumexp(f64 x, f64 y)
{
    if(x == - INFINITY && y == - INFINITY) { return - INFINITY;                }
    if(x >= y)                             { return x + bayla_log1pexp(y - x); }
    else                                   { return y + bayla_log1pexp(x - y); }
}

API f64 
bayla_logsumexp_list(size nvals, size offset, size stride, f64 *vals)
{
    /* Use a two pass algorithm, first find the max, then scale and apply addition at higher precision. */
    f64 max = -INFINITY;
    for(size i = 0; i < nvals; ++i)
    {
        size idx = i * stride + offset;
        if(vals[idx] > max) { max = vals[idx]; }
    }
    Assert(max > - INFINITY);

    ElkKahanAccumulator result = {0};
    for(size i = 0; i < nvals; ++i)
    {
        size idx = i * stride + offset;
        result = elk_kahan_accumulator_add(result, exp(vals[idx] - max));
    }

    return max + log(result.sum);
}

API BayLaNormalDistribution 
bayla_normal_distribution_create(f64 mean, f64 stddev)
{
    return (BayLaNormalDistribution){.mu = mean, .sigma = stddev};
}

API f64 
bayla_normal_distribution_random_deviate(BayLaNormalDistribution *dist, ElkRandomState *state)
{
    f64 u = 0.0, v = 0.0, q = 0.0;
    do 
    {
        u = elk_random_state_uniform_f64(state);

        v = elk_random_state_uniform_f64(state);
        v = 1.7156 * (v - 0.5);

        f64 x = u - 0.449871;

        f64 y = fabs(v) + 0.386595;

        q = x * x + y * (0.19600 * y - 0.25472 * x);

    } while(q > 0.27597 && (q > 0.27846 || (v * v > -4.0 * log(u) * u * u)));

    return dist->mu + dist->sigma * v / u;
}

API b32
bayla_mv_norm_dist_validate_internal(BayLaMVNormDist *dist)
{
    b32 result = dist->valid;

    /* Check that the normalization constant doesn't blow us up! (In the log domain,-INFINITY implies a constant of 0) */
    result &= !(isinf(dist->normalization_constant.val) || isnan(dist->normalization_constant.val));

    for(i32 i = 0; i < dist->ndim * dist->ndim; ++i)
    {
        result &= !(isinf(dist->chol.data[i]) || isnan(dist->chol.data[i]));
        result &= !(isinf(dist->cov.data[i]) || isnan(dist->cov.data[i]));
        result &= !(isinf(dist->inv_cov.data[i]) || isnan(dist->inv_cov.data[i]));
    }

    return result;
}

API BayLaMVNormDist 
bayla_mv_norm_dist_create(size ndim, f64 *mean, BayLaSquareMatrix covar, MagAllocator *alloc)
{
    MagStaticArena scratch_ = mag_static_arena_allocate_and_create(ECO_KiB(20));
    MagAllocator scratch = mag_allocator_from_static_arena(&scratch_);

    BayLaMVNormDist dist = { .ndim = ndim, .valid = true };
    f64 *dist_mean = eco_arena_nmalloc(alloc, ndim, f64);
    memcpy(dist_mean, mean, ndim * sizeof(f64));
    dist.mean = dist_mean;

    dist.chol = bayla_cholesky_decomp(covar, alloc);
    dist.valid = dist.chol.valid;
    Assert(dist.valid);

    dist.cov = bayla_square_matrix_create(ndim, alloc);
    memcpy(dist.cov.data, covar.data, sizeof(f64) * ndim * ndim);

    dist.inv_cov = bayla_square_matrix_invert(dist.cov, alloc, scratch);

    BayLaLogValue covar_det = bayla_log_one;
    for(i32 i = 0; i < dist.ndim; ++i)
    {
        BayLaLogValue diag_val = bayla_log_value_map_into_log_domain(dist.chol.data[MATRIX_IDX(i, i, dist.ndim)]);
        covar_det = bayla_log_value_multiply(covar_det, bayla_log_value_power(diag_val, 2.0));
    }
    dist.normalization_constant = bayla_log_value_reciprocal(bayla_log_value_power(covar_det, 0.5));
    dist.normalization_constant.val += (-0x1.d67f1c864beb4p-1 * ndim); /* Literal of log(1/sqrt(2*pi)) */

    dist.valid = bayla_mv_norm_dist_validate_internal(&dist);

    eco_arena_destroy(&scratch);

    return dist;
}

API BayLaMVNormDist 
bayla_mv_norm_dist_create_from_hessian(size ndim, f64 *mean, BayLaSquareMatrix hessian, MagAllocator *alloc)
{
    MagStaticArena scratch_ = mag_static_arena_allocate_and_create(ECO_KiB(20));
    MagAllocator scratch = mag_allocator_from_static_arena(&scratch_);

    BayLaMVNormDist dist = { .ndim = ndim, .valid = true };
    f64 *dist_mean = eco_arena_nmalloc(alloc, ndim, f64);
    memcpy(dist_mean, mean, ndim * sizeof(f64));
    dist.mean = dist_mean;

    dist.inv_cov = bayla_square_matrix_create(ndim, alloc);
    for(size i = 0; i < ndim * ndim; ++i) { dist.inv_cov.data[i] = -hessian.data[i]; }

    dist.cov = bayla_square_matrix_invert(dist.inv_cov, alloc, scratch);

    for(size i = 0; i < ndim; ++i) { Assert(dist.cov.data[MATRIX_IDX(i, i, ndim)] > 0.0); } /* Check positive variance. */

    dist.chol = bayla_cholesky_decomp(dist.cov, alloc);
    dist.valid = dist.chol.valid;

    BayLaLogValue covar_det = bayla_log_one;
    for(i32 i = 0; i < dist.ndim; ++i)
    {
        BayLaLogValue diag_val = bayla_log_value_map_into_log_domain(dist.chol.data[MATRIX_IDX(i, i, dist.ndim)]);
        covar_det = bayla_log_value_multiply(covar_det, bayla_log_value_power(diag_val, 2.0));
    }
    dist.normalization_constant = bayla_log_value_reciprocal(bayla_log_value_power(covar_det, 0.5));
    dist.normalization_constant.val += (-0x1.d67f1c864beb4p-1 * ndim); /* Literal of log(1/sqrt(2*pi)) */

    eco_arena_destroy(&scratch);

    return dist;
}

API BayLaLogValue 
bayla_mv_norm_dist_random_deviate(BayLaMVNormDist *dist, ElkRandomState *state, f64 *deviate, f64 *workspace)
{
    f64 log_prob = 0.0; /* log(1.0) */
    f64 u = 0.0, v = 0.0, q = 0.0;

    for(i32 i = 0; i < dist->ndim; ++i)
    {
        /* Inline normal deviate calculation. */
        do 
        {
            u = elk_random_state_uniform_f64(state);

            v = elk_random_state_uniform_f64(state);
            v = 1.7156 * (v - 0.5);

            f64 x = u - 0.449871;

            f64 y = fabs(v) + 0.386595;

            q = x * x + y * (0.19600 * y - 0.25472 * x);

        } while(q > 0.27597 && (q > 0.27846 || (v * v > -4.0 * log(u) * u * u)));

        workspace[i] = v / u;
        log_prob += workspace[i] * workspace[i];
    }

    log_prob *= -0.5;
    BayLaLogValue prob = bayla_log_value_create(log_prob);
    prob = bayla_log_value_multiply(prob, dist->normalization_constant);

    bayla_cholesky_multiply(dist->chol, dist->ndim, workspace, deviate);

    for(i32 i = 0; i < dist->ndim; ++i)
    {
        deviate[i] += dist->mean[i];
    }

    return prob;
}

API BayLaLogValue 
bayla_mv_norm_dist_prob(BayLaMVNormDist *dist, f64 const *deviate, f64 *workspace)
{
    /* Calculate the displacement from the mean. */
    f64 *displacement = workspace;
    f64 *output = &workspace[dist->ndim];
    for(i32 i = 0; i < dist->ndim; ++i)
    {
        displacement[i] = deviate[i] - dist->mean[i];
    }

    bayla_square_matrix_left_multiply(dist->inv_cov, displacement, output);
    f64 log_prob = -0.5 * bayla_dot_product(dist->ndim, output, displacement);
    BayLaLogValue prob = bayla_log_value_create(log_prob);
    prob = bayla_log_value_multiply(prob, dist->normalization_constant);
    return prob;
}

API BayLaLogValue 
bayla_model_evaluate(BayLaModel const *model, f64 const *parameter_vals)
{
    size ndim = model->n_parameters;
    void *user_data = model->user_data;
    BayLaLogPrior prior = model->prior;
    BayLaLogLikelihood likelihood = model->likelihood;

    f64 log_prob = prior(ndim, parameter_vals, user_data) + likelihood(ndim, parameter_vals, user_data);
    
    return bayla_log_value_create(log_prob);
}

API BayLaSamples 
bayla_importance_sample_gauss_approx(BayLaModel const *model, f64 sf, size n_samples, u64 seed, MagAllocator *alloc, MagAllocator scratch_)
{
    MagAllocator *scratch = &scratch_;

    /* Fixed sizes and indexes. */
    size const ndim = model->n_parameters;
    size const pidx = ndim + 0;
    size const qidx = ndim + 1;
    size const widx = ndim + 2;
    size const ncols = ndim + 3;            /* The last 3 are Q, P, and W */
    f64 *rows = eco_arena_nmalloc(alloc, ncols * n_samples, f64);
    Assert(rows);

    BayLaSamples samples =
        {
            .n_samples = n_samples,
            .ndim = ndim,
            .neff = NAN,
            .z_evidence = (BayLaLogValue){ .val = NAN},
            .rows = rows,
            .valid = true
        };

    ElkRandomState state_ = elk_random_state_create(seed);
    ElkRandomState *state = &state_;

    /* The proposal distribution. */
    f64 *prop_mean = eco_arena_nmalloc(scratch, ndim, f64);
    model->posterior_max(model->user_data, ndim, prop_mean);
    BayLaSquareMatrix prop_hessian = model->hessian_matrix(model->user_data, ndim, prop_mean, scratch);

    /* Apply the spread factor. */
    for(size i = 0; i < ndim * ndim; ++i)
    {
        prop_hessian.data[i] /= sf;
    }

    BayLaMVNormDist prop_dist = bayla_mv_norm_dist_create_from_hessian(ndim, prop_mean, prop_hessian, scratch);
    if(!prop_dist.valid)
    {
        samples.valid = prop_dist.valid;
        return samples;
    }

    f64 *workspace = eco_arena_nmalloc(scratch, ndim, f64);
    Assert(prop_dist.valid && workspace);

    BayLaLogValue *ws = eco_arena_nmalloc(scratch, n_samples, BayLaLogValue);
    BayLaLogValue *w2s = eco_arena_nmalloc(scratch, n_samples, BayLaLogValue);

    for(size s = 0; s < n_samples; ++s)
    {
        f64 *row = &rows[s * ncols];
        BayLaLogValue ld_qi = bayla_mv_norm_dist_random_deviate(&prop_dist, state, row, workspace);
        BayLaLogValue ld_pi = bayla_model_evaluate(model, row);
        row[pidx] = ld_pi.val;
        row[qidx] = ld_qi.val;
        row[widx] = row[pidx] - row[qidx];

        ws[s] = bayla_log_value_create(row[widx]);
        w2s[s] = bayla_log_value_power(ws[s], 2.0);
    }

    BayLaLogValue w_acc = bayla_log_values_sum(n_samples, 0, 1, ws);
    BayLaLogValue w2_acc = bayla_log_values_sum(n_samples, 0, 1, w2s);
    BayLaLogValue l_n_samples = bayla_log_value_map_into_log_domain(n_samples);

    /* w_acc^2 / w2_acc */
    samples.neff = bayla_log_value_map_out_of_log_domain(bayla_log_value_divide(bayla_log_value_power(w_acc, 2), w2_acc));
    samples.z_evidence = bayla_log_value_divide(w_acc, l_n_samples);

    return samples;
}

API void 
bayla_samples_save_csv(BayLaSamples *samples, f64 ci_p_thresh, char const *fname)
{
    FILE *f = fopen(fname, "wb");

    /* Sort them in order of increasing weight. */
    size const n_samples = samples->n_samples;
    size const ndim = samples->ndim;
    size const pidx = ndim + 0;
    size const qidx = ndim + 1;
    size const widx = ndim + 2;
    size const ncols = ndim + 3;

    MagStaticArena scratch = mag_static_arena_allocate_and_create(ncols * n_samples * sizeof(f64));
    f64 *row_scratch = eco_arena_nmalloc(&scratch, ncols * n_samples, f64);
    Assert(row_scratch);

    pak_radix_sort(
            samples->rows,             /* The buffer to sort.                                      */
            n_samples,                 /* The number of items in the buffer.                       */
            widx * sizeof(f64),        /* Offset into the object for the item we're sorting by.    */
            sizeof(f64) * ncols,       /* The stride, or how large a single row is.                */
            row_scratch,               /* Scratch memory for doing the search.                     */
            PAK_RADIX_SORT_F64,        /* Sorting f64s.                                            */
            PAK_SORT_ASCENDING);       /* Sort in ascending order.                                 */


    for(size r = 0; r < n_samples; ++r)
    {
        f64 *row = &samples->rows[r * ncols];
        for(size c = 0; c < ndim; ++c)
        {
            fprintf(f, "%.18e,", row[c]);
        }
        fprintf(f, "%.18e,%.18e,%.18e,%d\n", row[pidx], row[qidx], row[widx], row[pidx] >= ci_p_thresh ? 1 : 0);
    }

    fclose(f);

    mag_static_arena_destroy(&scratch);
}

API f64 
bayla_samples_calculate_ci_p_thresh(BayLaSamples *samples, f64 ci_pct, MagAllocator scratch)
{
    Assert(ci_pct > 0.0 && ci_pct <= 1.0);

    size const n_samples = samples->n_samples;
    size const ndim = samples->ndim;
    size const ncols = ndim + 3;
    size const pidx = ndim + 0;
    size const widx = ndim + 2;

    f64 total_weight = bayla_log_value_map_out_of_log_domain(samples->z_evidence) * n_samples;

    f64 *row_scratch = eco_arena_nmalloc(&scratch, ncols * n_samples, f64);
    Assert(row_scratch);

    pak_radix_sort(
            samples->rows,             /* The buffer to sort.                                      */
            n_samples,                 /* The number of items in the buffer.                       */
            pidx * sizeof(f64),        /* Offset into the object for the item we're sorting by.    */
            sizeof(f64) * ncols,       /* The stride, or how large a single row is.                */
            row_scratch,               /* Scratch memory for doing the search.                     */
            PAK_RADIX_SORT_F64,        /* Sorting f64s.                                            */
            PAK_SORT_DESCENDING);      /* Sort in descending order. We want the most weight first. */

    ElkKahanAccumulator weight = {0};
    for(size r = 0; r < n_samples; ++r)
    {
        f64 *row = &samples->rows[r * ncols];
        weight = elk_kahan_accumulator_add(weight, exp(row[widx]));
        if(weight.sum >= ci_pct * total_weight)
        {
            return row[pidx];
        }
    }

    return NAN;
}

typedef struct
{
    f64 prob;
    f64 wgt;
    f64 x;
} BayLaInternalCISortRecord;

API BayLaCredibleInterval 
bayla_samples_calculate_ci(BayLaSamples *samples, f64 ci_pct, size param_idx, MagAllocator scratch)
{
    Assert(ci_pct > 0.0 && ci_pct <= 1.0);
    Assert(param_idx >= 0 && param_idx < samples->ndim);

    size const n_samples = samples->n_samples;
    size const ndim = samples->ndim;
    size const ncols = ndim + 3;
    size const pidx = ndim + 0;
    size const widx = ndim + 2;
    size const nbins = n_samples / 10;

    f64 const total_weight = bayla_log_value_map_out_of_log_domain(samples->z_evidence) * n_samples;
    f64 const *const rows = samples->rows;

    /* Find the minimum and maximum values in our dataset. */
    f64 min_x = DBL_MAX;
    f64 max_x = -DBL_MAX;
    for(size r = 0; r < n_samples; ++r)
    {
        f64 const x = rows[r * ncols + param_idx];
        min_x = min_x < x ? min_x : x;
        max_x = max_x > x ? max_x : x;
    }

    /* Set up the bins and their values. */
    f64 dx = (max_x - min_x) / nbins;
    BayLaInternalCISortRecord *recs = eco_arena_nmalloc(&scratch, nbins, BayLaInternalCISortRecord);
    for(size b = 0; b < nbins; ++b)
    {
        recs[b].x = min_x + b * dx + 0.5 * dx;
    }

    /* Project the data probability (counts) and weights onto the parameter axis. */
    for(size r = 0; r < n_samples; ++r)
    {
        f64 const x = rows[r * ncols + param_idx];
        f64 const w = exp(rows[r * ncols + widx]);
        f64 const p = exp(rows[r * ncols + pidx]);
        size bin = (size)round((x - min_x) / dx);

        recs[bin].prob += p;
        recs[bin].wgt += w;
    }

    /* Sort in order of decreasing probability (count) */
    BayLaInternalCISortRecord *sort_scratch = eco_arena_nmalloc(&scratch, nbins, BayLaInternalCISortRecord);
    pak_radix_sort(
            recs,                  /* The buffer to sort.                                   */
            nbins,                 /* The number of items in the buffer.                    */
            0,                     /* Offset into the object for the item we're sorting by. */
            sizeof(*recs),         /* The stride, or how large a single row is.             */
            sort_scratch,          /* Scratch memory for doing the search.                  */
            PAK_RADIX_SORT_F64,    /* Sorting F64s.                                         */
            PAK_SORT_DESCENDING);  /* Sort in descending order.                             */

    BayLaCredibleInterval ci = {.low = DBL_MAX, .high = -DBL_MAX };
    ElkKahanAccumulator weight = {0};
    size bin = -1;
    while(weight.sum < ci_pct * total_weight && bin < nbins - 1)
    {
        ++bin;
        BayLaInternalCISortRecord rec = recs[bin];
        weight = elk_kahan_accumulator_add(weight, rec.wgt);
        ci.low = rec.x < ci.low ? rec.x : ci.low;
        ci.high = rec.x > ci.high ? rec.x : ci.high;
    }

    if(ci.high < ci.low)
    {
        ci.high = ci.low = NAN;
    }

    return ci;
}

API BayLaErrorValue 
bayla_samples_estimate_evidence(BayLaSamples *samples)
{
    size const n_samples = samples->n_samples;
    size const ndim = samples->ndim;
    size const ncols = ndim + 3;
    size const widx = ndim + 2;
    f64 const *rows = samples->rows;
    f64 const z = bayla_log_value_map_out_of_log_domain(samples->z_evidence);

    ElkKahanAccumulator v_acc = {0};

    for(size r = 0; r < n_samples; ++r)
    {
        f64 const *row = &rows[r * ncols];
        f64 dw = exp(row[widx]) - z;

        v_acc = elk_kahan_accumulator_add(v_acc, dw * dw);
    }

    f64 std = sqrt(v_acc.sum / (n_samples * n_samples));

    return (BayLaErrorValue){ .val = z, .std = std};
}

API BayLaErrorValue 
bayla_samples_calculate_expectation(BayLaSamples *samples, f64 (*func)(size n, f64 const *params, void *ud), void *ud)
{
    size const n_samples = samples->n_samples;
    size const ndim = samples->ndim;
    size const ncols = ndim + 3;
    size const widx = ndim + 2;
    f64 const *rows = samples->rows;
    f64 const total_weight = bayla_log_value_map_out_of_log_domain(samples->z_evidence) * n_samples;

    ElkKahanAccumulator f_acc = {0};

    for(size r = 0; r < n_samples; ++r)
    {
        f64 const *row = &rows[r * ncols];
        f64 w = exp(row[widx]);
        f64 f = func(ndim, row, ud);

        f_acc = elk_kahan_accumulator_add(f_acc, f * w);
    }

    f64 mean = f_acc.sum / total_weight;

    ElkKahanAccumulator v_acc = {0};
    for(size r = 0; r < n_samples; ++r)
    {
        f64 const *row = &rows[r * ncols];
        f64 w = row[widx];
        f64 df = func(ndim, row, ud) * w - mean;

        v_acc = elk_kahan_accumulator_add(v_acc, df * df);
    }

    f64 std = sqrt(v_acc.sum / (n_samples * n_samples));

    return (BayLaErrorValue){ .val = mean, .std = std};
}

typedef struct
{
    f64 xa;
    f64 xb;
    f64 xc;
    f64 fxa;
    f64 fxb;
    f64 fxc;
} BayLaBracketPoints;

static inline void shft2(f64 *a, f64 *b, f64 const c) { *a = *b; *b = c; }
static inline void shft3(f64 *a, f64 *b, f64 *c, f64 const d) { *a = *b; *b = *c; *c = d; }
static inline void mov3(f64 *a, f64 *b, f64 *c, f64 const d, f64 const e, f64 const f) { *a = d; *b = e; *c = f; }
static inline void swap(f64 *a, f64 *b) { f64 tmp = *a; *a = *b; *b = tmp; }
static inline f64 sign(f64 const a, f64 const b) { return b >= 0.0 ? (a >= 0.0 ? a : -a) : (a >= 0.0 ? -a : a); }
static inline f64 max(f64 const a, f64 const b) { return a >= b ? a : b; }
static inline f64 min(f64 const a, f64 const b) { return a <= b ? a : b; }

static inline f64
bayla_evaluate_samples_for_ess(f64 sf, BayLaModel const *model, size n_samples, u64 seed, MagAllocator scratch1, MagAllocator scratch2)
{
    BayLaSamples samples = bayla_importance_sample_gauss_approx(model, exp(sf), n_samples, seed, &scratch1, scratch2);

    return samples.valid ? -samples.neff : NAN;
}

/* We're looking for a maxima, but bayla_evaluate_samples_for_ess turns it to a negative, so now we're looking for a minima.
 * It's just easier that way because the algorithm was written for bracketing a minima.
 */
static BayLaBracketPoints
bayla_bracket_minima(BayLaModel const *model, size n_samples, u64 seed, MagAllocator scratch1, MagAllocator scratch2)
{
    f64 const GOLD = 1.618034;
    f64 const GLIMIT = 10.0;
    f64 const TINY = 1.0e-20;

    f64 *max_a_posteriori_parms = eco_arena_nmalloc(&scratch1, model->n_parameters, f64);
    model->posterior_max(model->user_data, model->n_parameters, max_a_posteriori_parms);
    f64 max_density = bayla_model_evaluate(model, max_a_posteriori_parms).val;

    f64 ax = max_density + log(5.0);
    f64 fa = bayla_evaluate_samples_for_ess(ax, model, n_samples, seed, scratch1, scratch2);

    f64 bx = max_density + log(0.1);;
    f64 fb = bayla_evaluate_samples_for_ess(bx, model, n_samples, seed, scratch1, scratch2);

    if(fb > fa) { swap(&ax, &bx); swap(&fa, &fb); }

    f64 cx = bx + GOLD * (bx - ax);
    f64 fc = bayla_evaluate_samples_for_ess(cx, model, n_samples, seed, scratch1, scratch2);

    f64 fu = NAN;
    while(fb > fc)
    {
        f64 r = (bx - ax) * (fb - fc);
        f64 q = (bx - cx) * (fb - fa);
        f64 u = bx - ((bx - cx) * q - (bx - ax) * r) / (2.0 * sign(max(fabs(q - r), TINY), q - r));
        f64 ulim = bx + GLIMIT * (cx - bx);

        if((bx - u) * (u - cx) > 0.0)
        {
            fu = bayla_evaluate_samples_for_ess(u, model, n_samples, seed, scratch1, scratch2);
            if(fu < fc)
            {
                ax = bx; fa = fb;
                bx = u;  fb = fu;
                return (BayLaBracketPoints){ .xa = ax, .xb = bx, .xc = cx, .fxa = fa, .fxb = fb, .fxc = fc };
            }
            else if (fu > fb)
            {
                cx = u; fc = fu;
                return (BayLaBracketPoints){ .xa = ax, .xb = bx, .xc = cx, .fxa = fa, .fxb = fb, .fxc = fc };
            }

            u = cx + GOLD * (cx - bx);
            fu = bayla_evaluate_samples_for_ess(u, model, n_samples, seed, scratch1, scratch2);
        }
        else if ((cx - u) * (u - ulim) > 0.0)
        {
            fu = bayla_evaluate_samples_for_ess(u, model, n_samples, seed, scratch1, scratch2);
            if(fu < fc)
            {
                shft3(&bx, &cx, &u, u + GOLD * (u - cx));
                shft3(&fb, &fc, &fu, bayla_evaluate_samples_for_ess(u, model, n_samples, seed, scratch1, scratch2));
            }

        }
        else if ((u - ulim) * (ulim - cx) >= 0.0)
        {
            u = ulim;
            fu = bayla_evaluate_samples_for_ess(u, model, n_samples, seed, scratch1, scratch2);
        }
        else
        {
            u = cx + GOLD * (cx - bx);
            fu = bayla_evaluate_samples_for_ess(u, model, n_samples, seed, scratch1, scratch2);
        }

        shft3(&ax, &bx, &cx, u);
        shft3(&fa, &fb, &fc, fu);
    }

    return (BayLaBracketPoints){ .xa = ax, .xb = bx, .xc = cx, .fxa = fa, .fxb = fb, .fxc = fc };
}

typedef struct
{
    f64 xmin;
    f64 fmin;
} BayLaMinimizationPair;

/* We're looking for a maxima, but bayla_evaluate_samples_for_ess turns it to a negative, so now we're looking for a minima.
 * It's just easier that way because the algorithm was written for bracketing a minima.
 */
static BayLaMinimizationPair
bayla_find_minima(BayLaModel const *model, size n_samples, u64 seed, MagAllocator scratch1, MagAllocator scratch2)
{
    BayLaBracketPoints bp = bayla_bracket_minima(model, n_samples, seed, scratch1, scratch2);
    //printf("xa = %e xb = %e xc = %e\n", bp.xa, bp.xb, bp.xc);
    //printf("fa = %e fb = %e fc = %e\n", bp.fxa, bp.fxb, bp.fxc);

    size const ITMAX = 1000;
    f64 const CGOLD = 0.3819660;
    f64 const ZEPS = 0.0; /* DBL_EPSILON * 1.0e-3; */ /* We're going to be dealing with VERY small numbers. */
    f64 const tol = 3.0e-8;

    f64 d = NAN;
    f64 e = 0.0;

    f64 a = (bp.xa < bp.xc ? bp.xa : bp.xc);
    f64 b = (bp.xa > bp.xc ? bp.xa : bp.xc);
    f64 x = bp.xb;
    f64 w = bp.xb;
    f64 v = bp.xb;
    f64 u = NAN;

    f64 fx = bayla_evaluate_samples_for_ess(x, model, n_samples, seed, scratch1, scratch2);
    f64 fw = fx;
    f64 fv = fx;
    f64 fu = NAN;

    for(size i = 0; i < ITMAX; ++i)
    {
        f64 xm = 0.5 * (a + b);
        f64 tol1 = tol * fabs(x) + ZEPS;
        f64 tol2 = 2.0 * tol1;
        if(fabs(x - xm) <= (tol2 - 0.5 *(b - a))) { return (BayLaMinimizationPair){ .xmin = exp(x), .fmin = fx }; }
        if(fabs(e) > tol1)
        {
            f64 r = (x - w) * (fx - fv);
            f64 q = (x - v) * (fx - fw);
            f64 p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if(q > 0.0) { p = -p; }
            q = fabs(q);
            f64 etemp = e;
            e = d;
            if(fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
            {
                e = x >= xm ? a - x : b - x;
                d = CGOLD * e;
            }
            else
            {
                d = p / q;
                u = x + d;
                if(u - a < tol2 || b - u < tol2) { d = sign(tol1, xm - x); }
            }
        }
        else
        {
            e = x > xm ? a - x : b - x;
            d = CGOLD * e;
        }

        u = fabs(d) >= tol1 ? x + d : x + sign(tol1, d);
        fu = bayla_evaluate_samples_for_ess(u, model, n_samples, seed, scratch1, scratch2);

        if(fu <= fx)
        {
            if(u >= x) { a = x; } else { b = x; }
            shft3(&v, &w, &x, u);
            shft3(&fv, &fw, &fx, fu);
        }
        else
        {
            if(u < x) { a = u; } else { b = u; }
            if(fu <= fw || w == x)
            {
                v = w; fv = fw;
                w = u; fw = fu;
            }
            else if(fu <= fv || v == x || v == w)
            {
                v = u; fv = fu;
            }
        }
    }

    Panic();

    return (BayLaMinimizationPair) { .xmin = NAN, .fmin = NAN };
}

API BayLaSamples 
bayla_importance_sample_gauss_approx_optimize(BayLaModel const *model, size n_samples, u64 seed, MagAllocator *alloc, MagAllocator scratch1, MagAllocator scratch2)
{
    BayLaMinimizationPair min_result = bayla_find_minima(model, n_samples, seed, scratch1, scratch2);
    //printf("min_result.xmin = %e min_result.fmin = %e\n", min_result.xmin, min_result.fmin);

    return bayla_importance_sample_gauss_approx(model, min_result.xmin, n_samples, seed, alloc, scratch1);
}

#pragma warning(pop)

#endif

