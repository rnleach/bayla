#ifndef _BAYLA_H_
#define _BAYLA_H_

#include <math.h>
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
API void bayla_square_matrix_print(BayLaSquareMatrix matrix);
API void bayla_square_matrix_print_symmetric_error(BayLaSquareMatrix matrix);
API BayLaSquareMatrix bayla_square_matrix_covariance_to_corr(BayLaSquareMatrix matrix, MagAllocator *alloc);
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

#define bayla_log_value_equal(left, right) ((left.val) == (right.val))
#define bayla_log_value_is_zero(log_val) (isinf(log_val.val) && (log_val.val) < 0.0)
#define bayla_log_value_is_nan(log_val) (isnan(log_val.val))

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

API BayLaMVNormDist bayla_mv_norm_dist_create(i32 ndim, f64 *mean, BayLaSquareMatrix covar, MagAllocator *alloc);
API BayLaMVNormDist bayla_mv_norm_dist_create_from_hessian(i32 ndim, f64 *mean, BayLaSquareMatrix hessian, MagAllocator *alloc);

/* returns the probability value of the associated point, workspace needs to be at least as big as dist->ndim */
API BayLaLogValue bayla_mv_norm_dist_random_deviate(BayLaMVNormDist const *dist, ElkRandomState *state, f64 *deviate, f64 *workspace);
/* returns the probability value of the input value deviate. workspace needs to be twice dist->ndim */
API BayLaLogValue bayla_mv_norm_dist_prob(BayLaMVNormDist *dist, f64 const *deviate, f64 *workspace);

typedef struct
{
    i32 ndim;
    f64 *mins;
    f64 *maxs;
    BayLaLogValue prob;
} BayLaMVUniformDist;

API BayLaMVUniformDist bayla_mv_uniform_dist_create(i32 ndim, f64 *mins, f64 *maxs, MagAllocator *alloc);
API BayLaLogValue bayla_mv_uniform_dist_random_deviate(BayLaMVUniformDist *dist, ElkRandomState *state, f64 *deviate);
API BayLaLogValue bayla_mv_uniform_dist_prob(BayLaMVUniformDist *dist, f64 const *deviate);

typedef struct
{
    i32 ndim;
    b32 valid;
    f64 *mean;
    f64 nu;
    BayLaCholeskyMatrix chol; /* Cholesky factor of sigma */
    f64 log_det_sigma;
    f64 log_const;
} BayLaMVStudenTDist;

API BayLaMVStudenTDist bayla_mv_studt_dist_create(i32 ndim, f64 *mean, BayLaSquareMatrix covar, f64 nu, MagAllocator *alloc);
API BayLaMVStudenTDist bayla_mv_studt_dist_create_from_hessian(i32 ndim, f64 *mean, BayLaSquareMatrix hessian, f64 nu, MagAllocator *alloc);
API BayLaLogValue bayla_mv_studt_dist_random_deviate(BayLaMVStudenTDist *dist, ElkRandomState *state, f64 *deviate);
API BayLaLogValue bayla_mv_studt_dist_prob(BayLaMVStudenTDist *dist, f64 const *deviate);

/*---------------------------------------------------------------------------------------------------------------------------
 *                                              Analysis of Bayesian Models
 *-------------------------------------------------------------------------------------------------------------------------*/

/* Description of a model to be used for a Baysiean analysis. */
typedef struct
{
    
    f64 model_prior_probability;    /* Prior prob *this model* is the correct one. P(Model | I), Jaynes eq. 20.3 (pg 602) */
    i32 n_parameters;               /* The number of parameters in this model.                                            */

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
    f64 *rows;             /* n_samples x (ndim + 3), last 3 columns are in log domain, wrap them BayLaLogValue to use */
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
API BayLaSamples bayla_importance_sample_gauss_approx(BayLaModel const *model, size n_samples, u64 seed, MagAllocator *alloc, MagAllocator scratch);
API BayLaSamples bayla_importance_sample_gauss_approx_par(BayLaModel const *model, size n_samples, u64 seed, MagAllocator *alloc, MagAllocator scratch, CoyThreadPool *pool);
API BayLaSamples bayla_importance_sample_uniform(BayLaModel const *model, size n_samples, u64 seedd, f64 *mins, f64 *maxs, MagAllocator *alloc, MagAllocator scratch);
API BayLaSamples bayla_importance_sample_studt_approx(BayLaModel const *model, f64 nu, size n_samples, u64 seed, MagAllocator *alloc, MagAllocator scratch);

API void bayla_samples_save_csv(BayLaSamples *samples, BayLaLogValue ci_p_thresh, char const *fname);
API BayLaLogValue bayla_samples_calculate_ci_p_thresh(BayLaSamples *samples, f64 ci_pct, MagAllocator scratch);

API BayLaCredibleInterval bayla_samples_calculate_ci(BayLaSamples *samples, f64 ci_pct, size param_idx, MagAllocator scratch);
API BayLaErrorValue bayla_samples_estimate_evidence(BayLaSamples *samples);
API BayLaErrorValue bayla_samples_calculate_expectation(BayLaSamples *samples, f64 (*func)(size n, f64 const *params, void *ud), void *ud);
API f64 bayla_samples_calculate_prob_predicate(BayLaSamples *samples, b32 (*func)(size n, f64 const *params, void *ud), void *ud);

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

static inline void swap(f64 *a, f64 *b) { f64 tmp = *a; *a = *b; *b = tmp; }
static inline f64 sign(f64 const a, f64 const b) { return b >= 0.0 ? (a >= 0.0 ? a : -a) : (a >= 0.0 ? -a : a); }
static inline f64 maximum(f64 const a, f64 const b) { return a >= b ? a : b; }
static inline f64 minimum(f64 const a, f64 const b) { return a <= b ? a : b; }

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

    /* Assume symmetric matrix. */
    for(i32 r = 0; r < n; ++r)
    {
        for(i32 c = 0; c < r; ++c)
        {
            a[MATRIX_IDX(r, c, n)] = a[MATRIX_IDX(c, r, n)] = 0.5 * (a[MATRIX_IDX(r, c, n)] + a[MATRIX_IDX(c, r, n)]);
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

    i32 *indxc = eco_arena_nmalloc(&scratch, n, i32); Assert(indxc);
    i32 *indxr = eco_arena_nmalloc(&scratch, n, i32); Assert(indxr);
    i32 *ipiv = eco_arena_nmalloc(&scratch, n, i32);  Assert(ipiv);
    
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

API void 
bayla_square_matrix_print(BayLaSquareMatrix matrix)
{
    printf("\n");
    i32 n = matrix.ndim;
    for(i32 r = 0; r < n; ++r)
    {
        printf("| ");
        for(i32 c = 0; c < n; ++c)
        {
            printf(" %10.2e ", matrix.data[MATRIX_IDX(r, c, n)]);
        }
        printf(" |\n");
    }
    printf("\n");
}

API void 
bayla_square_matrix_print_symmetric_error(BayLaSquareMatrix matrix)
{
    printf("\n");
    i32 n = matrix.ndim;
    for(i32 r = 0; r < n; ++r)
    {
        printf("| ");
        for(i32 c = 0; c < n; ++c)
        {
            f64 error_ratio = (matrix.data[MATRIX_IDX(r, c, n)] - matrix.data[MATRIX_IDX(c, r, n)]) / matrix.data[MATRIX_IDX(r, c, n)];
            if(isnan(error_ratio)) { error_ratio = 0.0; }
            printf(" %10.2e ", error_ratio);
        }
        printf(" |\n");
    }
    printf("\n");
}

API BayLaSquareMatrix 
bayla_square_matrix_covariance_to_corr(BayLaSquareMatrix matrix, MagAllocator *alloc)
{
    i32 const n = matrix.ndim;
    BayLaSquareMatrix corr = bayla_square_matrix_create(n, alloc);

    MagAllocator scratch = *alloc;

    f64 *sig = eco_arena_nmalloc(&scratch, n, f64);
    
    for(i32 i = 0; i < n; ++i)
    {
        sig[i] = 1.0 / sqrt(matrix.data[MATRIX_IDX(i, i, n)]);
    }

    for(i32 i = 0; i < n; ++i)
    {
        for(i32 j = 0; j < n; ++j)
        {
            corr.data[MATRIX_IDX(i, j, n)] = matrix.data[MATRIX_IDX(i, j, n)] * sig[i] * sig[j];
        }
    }

    return corr;
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
    if(max <= -INFINITY) { return -INFINITY; }
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
bayla_standard_normal_distribution_random_deviate(ElkRandomState *state)
{
    /* For internal use only. */
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

    return v / u;
}

API f64 
bayla_normal_distribution_random_deviate(BayLaNormalDistribution *dist, ElkRandomState *state)
{
    return dist->mu + dist->sigma * bayla_standard_normal_distribution_random_deviate(state);
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
bayla_mv_norm_dist_create(i32 ndim, f64 *mean, BayLaSquareMatrix covar, MagAllocator *alloc)
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
bayla_mv_norm_dist_create_from_hessian(i32 ndim, f64 *mean, BayLaSquareMatrix hessian, MagAllocator *alloc)
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
#if 0
    bayla_square_matrix_print_symmetric_error(dist.cov);
    bayla_square_matrix_print(bayla_square_matrix_covariance_to_corr(dist.cov, alloc));
#endif

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
bayla_mv_norm_dist_random_deviate(BayLaMVNormDist const *dist, ElkRandomState *state, f64 *deviate, f64 *workspace)
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

API BayLaMVUniformDist 
bayla_mv_uniform_dist_create(i32 ndim, f64 *mins, f64 *maxs, MagAllocator *alloc)
{
    BayLaMVUniformDist dist = {0};
    dist.ndim = ndim;
    
    f64 prob = 1.0;
    for(size i = 0; i < ndim; ++i)
    {
        Assert(maxs[i] > mins[i]);
        prob *= maxs[i] - mins[i];
    }

    dist.prob = bayla_log_value_reciprocal(bayla_log_value_map_into_log_domain(prob));

    dist.maxs = eco_arena_nmalloc(alloc, ndim, f64); Assert(maxs);
    memcpy(dist.maxs, maxs, ndim * sizeof(f64));
    dist.mins = eco_arena_nmalloc(alloc, ndim, f64); Assert(mins);
    memcpy(dist.mins, mins, ndim * sizeof(f64));

    return dist;
}

API BayLaLogValue 
bayla_mv_uniform_dist_random_deviate(BayLaMVUniformDist *dist, ElkRandomState *state, f64 *deviate)
{
    i32 ndim = dist->ndim;
    f64 *mins = dist->mins;
    f64 *maxs = dist->maxs;

    for(i32 i = 0; i < ndim; ++i)
    {
        f64 u = elk_random_state_uniform_f64(state);
        deviate[i] = u * (maxs[i] - mins[i]) + mins[i];
    }

    return dist->prob;
}

API BayLaLogValue 
bayla_mv_uniform_dist_prob(BayLaMVUniformDist *dist, f64 const *deviate)
{
    return dist->prob;
}

f64 
bayla_gamma_sample(f64 alpha, ElkRandomState *state)
{
    /* Marsaglia & Tsang Gamma sampler (shape=alpha, scale=1)
     *
     * Coded with assistance from Grok (xAI). Intended for internal use only.
     */

    if(alpha < 1.0) 
    {
        return bayla_gamma_sample(alpha + 1.0, state) * pow(elk_random_state_uniform_f64(state), 1.0 / alpha); 
    }

    f64 const d = alpha - 1.0 / 3.0;
    f64 const c = 1.0 / sqrt(9.0 * d);

    while(true)
    {
        f64 x, v;
        do
        {
            x = bayla_standard_normal_distribution_random_deviate(state);
            v = 1.0 + c * x;
        } while( v <= 0.0);

        v = v * v * v;
        f64 u = elk_random_state_uniform_f64(state);

        if( u < 1.0 - 0.0331 * x * x * x * x || log(u) < 0.5 * x * x + d * (1.0 - v + log(v))) 
        {
            return d * v;
        }
    }
}

API BayLaMVStudenTDist 
bayla_mv_studt_dist_create(i32 ndim, f64 *mean, BayLaSquareMatrix covar, f64 nu, MagAllocator *alloc)
{
    BayLaMVStudenTDist dist = { .ndim = ndim, .valid = true, .nu = nu };

    if(nu <= 0.0) { dist.valid = false; return dist; }

    dist.mean = eco_arena_nmalloc(alloc, ndim, f64); Assert(dist.mean);
    memcpy(dist.mean, mean, sizeof(*mean) * ndim);

    dist.chol = bayla_cholesky_decomp(covar, alloc);
    dist.valid = dist.chol.valid;
    StopIf(!dist.valid, return dist);

    f64 sum_log_diag = 0.0;
    for(i32 i = 0; i < ndim; ++i) { sum_log_diag += log(dist.chol.data[MATRIX_IDX(i, i, ndim)]); }
    dist.log_det_sigma = 2.0 * sum_log_diag;

    f64 half_nu_d = 0.5 * (nu + ndim);
    dist.log_const = lgamma(half_nu_d) - lgamma(0.5 * nu) - 0.5 * ndim * log(nu * ELK_PI) - 0.5 * dist.log_det_sigma;

    return dist;
}

API BayLaMVStudenTDist 
bayla_mv_studt_dist_create_from_hessian(i32 ndim, f64 *mean, BayLaSquareMatrix hessian, f64 nu, MagAllocator *alloc)
{
    MagStaticArena scratch_ = mag_static_arena_allocate_and_create(ECO_KiB(20));
    MagStaticArena scratch2_ = mag_static_arena_allocate_and_create(ECO_KiB(20));
    MagAllocator scratch = mag_allocator_from_static_arena(&scratch_);
    MagAllocator scratch2 = mag_allocator_from_static_arena(&scratch2_);

    BayLaSquareMatrix inv_cov = bayla_square_matrix_create(ndim, &scratch);
    for(size i = 0; i < ndim * ndim; ++i) { inv_cov.data[i] = -hessian.data[i]; }

    BayLaSquareMatrix cov = bayla_square_matrix_invert(inv_cov, &scratch, scratch2);

#if 0
    bayla_square_matrix_print_symmetric_error(dist.cov);
    bayla_square_matrix_print(bayla_square_matrix_covariance_to_corr(dist.cov, alloc));
#endif

    for(size i = 0; i < ndim; ++i) { Assert(cov.data[MATRIX_IDX(i, i, ndim)] > 0.0); } /* Check positive variance. */

    BayLaMVStudenTDist dist = bayla_mv_studt_dist_create(ndim, mean, cov, nu, alloc);

    eco_arena_destroy(&scratch);
    eco_arena_destroy(&scratch2);

    return dist;
}

API BayLaLogValue 
bayla_mv_studt_dist_random_deviate(BayLaMVStudenTDist *dist, ElkRandomState *state, f64 *deviate)
{
#define MAX_NDIM 12
    i32 const ndim = dist->ndim;
    f64 const nu = dist->nu;
    f64 const half_nu_d = 0.5 * (nu + (f64)ndim);

    f64 y_work_space[MAX_NDIM]; Assert(ndim <= MAX_NDIM);

    for(i32 i = 0; i < ndim; ++i) { y_work_space[i] = bayla_standard_normal_distribution_random_deviate(state); }
    bayla_cholesky_multiply(dist->chol, ndim, y_work_space, deviate);

    f64 g = bayla_gamma_sample(nu * 0.5, state) * (2.0 / nu);
    f64 scale = sqrt(nu / g);
    memset(y_work_space, 0, sizeof(y_work_space));
    for(i32 i = 0; i < ndim; ++i)
    {
        deviate[i] *= scale;

        f64 sum = deviate[i];
        for(i32 j = 0; j < i; ++j)
        {
            sum -= dist->chol.data[MATRIX_IDX(i, j, ndim)] * y_work_space[j];
        }
        y_work_space[i] = sum / dist->chol.data[MATRIX_IDX(i, i, ndim)];

        deviate[i] += dist->mean[i];
    }

    f64 quad = 0.0;
    for(i32 i = 0; i < ndim; ++i) { quad += y_work_space[i] * y_work_space[i]; }

    f64 logpdf = dist->log_const - half_nu_d * log(1.0 + quad / nu);

    return bayla_log_value_create(logpdf);
#undef MAX_NDIM
}

API BayLaLogValue 
bayla_mv_studt_dist_prob(BayLaMVStudenTDist *dist, f64 const *deviate)
{
#define MAX_NDIM 12
    i32 const ndim = dist->ndim;
    f64 const nu = dist->nu;
    f64 const half_nu_d = 0.5 * (nu + (f64)ndim);

    f64 z_work_space[MAX_NDIM]; Assert(ndim <= MAX_NDIM);

    for(i32 i = 0; i < ndim; ++i)
    {
        f64 sum = z_work_space[i] - dist->mean[i];
        for(i32 j = 0; j < i; ++j)
        {
            sum -= dist->chol.data[MATRIX_IDX(i, j, ndim)] * z_work_space[j];
        }
        z_work_space[i] = sum / dist->chol.data[MATRIX_IDX(i, i, ndim)];
    }

    f64 quad = 0.0;
    for(i32 i = 0; i < ndim; ++i) { quad += z_work_space[i] * z_work_space[i]; }

    f64 logpdf = dist->log_const - half_nu_d * log(1.0 + quad / nu);

    return bayla_log_value_create(logpdf);
#undef MAX_NDIM
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
bayla_importance_sample_gauss_approx(BayLaModel const *model, size n_samples, u64 seed, MagAllocator *alloc, MagAllocator scratch_)
{
    MagAllocator *scratch = &scratch_;

    /* Fixed sizes and indexes. */
    i32 const ndim = model->n_parameters;
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

    BayLaMVNormDist prop_dist = bayla_mv_norm_dist_create_from_hessian(ndim, prop_mean, prop_hessian, scratch);
    if(!prop_dist.valid)
    {
        samples.valid = prop_dist.valid;
        return samples;
    }

    f64 *workspace = eco_arena_nmalloc(scratch, ndim, f64);
    Assert(prop_dist.valid && workspace);

    BayLaLogValue *ws = eco_arena_nmalloc(scratch, n_samples, BayLaLogValue);

    for(size s = 0; s < n_samples; ++s)
    {
        f64 *row = &rows[s * ncols];
        BayLaLogValue ld_qi = bayla_mv_norm_dist_random_deviate(&prop_dist, state, row, workspace);
        BayLaLogValue ld_pi = bayla_model_evaluate(model, row);
        BayLaLogValue ld_wi = bayla_log_value_divide(ld_pi, ld_qi);

        if(
                   bayla_log_value_is_nan(ld_qi)
                || bayla_log_value_is_nan(ld_pi)
                || bayla_log_value_is_zero(ld_qi)
                || bayla_log_value_is_zero(ld_pi)
          )
        {
            --s; /* Redo this iteration. */
            continue;
        }

        row[pidx] = ld_pi.val;
        row[qidx] = ld_qi.val;
        row[widx] = ld_wi.val;

        ws[s] = ld_wi;
    }

    /* Calculate evidence */
    BayLaLogValue w_acc = bayla_log_values_sum(n_samples, 0, 1, ws);
    BayLaLogValue l_n_samples = bayla_log_value_map_into_log_domain((f64)n_samples);
    samples.z_evidence = bayla_log_value_divide(w_acc, l_n_samples);

    /* Normalize weights and square them for calculating effective sample size */
    for(size s = 0; s < n_samples; ++s)
    {
        ws[s] = bayla_log_value_divide(ws[s], w_acc);
        ws[s] = bayla_log_value_power(ws[s], 2.0);
    }

    w_acc = bayla_log_values_sum(n_samples, 0, 1, ws);
    samples.neff = bayla_log_value_map_out_of_log_domain(bayla_log_value_reciprocal(w_acc));

    if(isnan(samples.neff) || isinf(samples.neff)) { samples.valid = false; }

    return samples;
}

typedef struct
{
    BayLaModel const *model;
    BayLaMVNormDist const *prop_dist;
    f64 *workspace;
    ElkRandomState ran_state;
    size start_idx;
    size end_idx;
    f64 *rows;
    BayLaLogValue *ws;
    size n_samples;
    CoyBatchCompletion *bc;
} BayLaGaussApproxSamplerInnerArgs;

void
bayla_importance_sample_gauss_approx_par_inner(void *td)
{
    BayLaGaussApproxSamplerInnerArgs *args = td;

    i32 const ndim = args->model->n_parameters;
    size const pidx = ndim + 0;
    size const qidx = ndim + 1;
    size const widx = ndim + 2;
    size const ncols = ndim + 3;
    BayLaLogValue *ws = args->ws;
    f64 *rows = args->rows;

    for(size s = args->start_idx; s < args->end_idx; ++s)
    {
        Assert(s < args->n_samples);

        f64 *row = &rows[s * ncols];
        BayLaLogValue ld_qi = bayla_mv_norm_dist_random_deviate(args->prop_dist, &args->ran_state, row, args->workspace);
        BayLaLogValue ld_pi = bayla_model_evaluate(args->model, row);
        BayLaLogValue ld_wi = bayla_log_value_divide(ld_pi, ld_qi);

        if(
                   bayla_log_value_is_nan(ld_qi)
                || bayla_log_value_is_nan(ld_pi)
                || bayla_log_value_is_zero(ld_qi)
                || bayla_log_value_is_zero(ld_pi)
          )
        {
            --s; /* Redo this iteration. */
            continue;
        }

        row[pidx] = ld_pi.val;
        row[qidx] = ld_qi.val;
        row[widx] = ld_wi.val;

        Assert(s < args->n_samples);
        ws[s] = ld_wi;
    }

    coy_batch_completion_task_done(args->bc);
}

API BayLaSamples 
bayla_importance_sample_gauss_approx_par(BayLaModel const *model, size n_samples, u64 seed, MagAllocator *alloc, MagAllocator scratch_, CoyThreadPool *pool)
{
    MagAllocator *scratch = &scratch_;

    /* Fixed sizes and indexes. */
    i32 const ndim = model->n_parameters;
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

    /* The proposal distribution. */
    f64 *prop_mean = eco_arena_nmalloc(scratch, ndim, f64);
    model->posterior_max(model->user_data, ndim, prop_mean);
    BayLaSquareMatrix prop_hessian = model->hessian_matrix(model->user_data, ndim, prop_mean, scratch);

    BayLaMVNormDist prop_dist = bayla_mv_norm_dist_create_from_hessian(ndim, prop_mean, prop_hessian, scratch);
    if(!prop_dist.valid)
    {
        samples.valid = prop_dist.valid;
        return samples;
    }

    BayLaLogValue *ws = eco_arena_nmalloc(scratch, n_samples, BayLaLogValue);

    size const n_threads = pool->nthreads;
    BayLaGaussApproxSamplerInnerArgs *inner_args = eco_arena_nmalloc(scratch, n_threads, BayLaGaussApproxSamplerInnerArgs);
    CoyFuture *futures = eco_arena_nmalloc(scratch, n_threads, CoyFuture);
    CoyBatchCompletion bc = {0};
    coy_batch_completion_init(&bc, n_threads);

    size const block_size = n_samples / n_threads;
    size const left_overs = n_samples % n_threads;
    size start_idx = 0;
    for(size t = 0; t < n_threads; ++t)
    {
        inner_args[t].bc = &bc;
        inner_args[t].n_samples = n_samples;
        inner_args[t].model = model;
        inner_args[t].prop_dist = &prop_dist;
        inner_args[t].workspace = eco_arena_nmalloc(scratch, ndim, f64);
        Assert(inner_args[t].workspace);

        inner_args[t].ran_state = elk_random_state_create(seed + 11 * t);

        size end_idx = start_idx + block_size;
        end_idx = t < left_overs ? end_idx + 1 : end_idx;

        inner_args[t].start_idx = start_idx;
        inner_args[t].end_idx = end_idx;
        inner_args[t].rows = rows;
        inner_args[t].ws = ws;

        futures[t] = coy_future_create(bayla_importance_sample_gauss_approx_par_inner, &inner_args[t]);
        coy_threadpool_submit(pool, &futures[t]);

        start_idx = end_idx;
    }

    coy_batch_completion_wait(&bc);
    coy_batch_completion_destroy(&bc);

    /* Calculate evidence */
    BayLaLogValue w_acc = bayla_log_values_sum(n_samples, 0, 1, ws);
    BayLaLogValue l_n_samples = bayla_log_value_map_into_log_domain((f64)n_samples);
    samples.z_evidence = bayla_log_value_divide(w_acc, l_n_samples);

    /* Normalize weights and square them for calculating effective sample size */
    for(size s = 0; s < n_samples; ++s)
    {
        ws[s] = bayla_log_value_divide(ws[s], w_acc);
        ws[s] = bayla_log_value_power(ws[s], 2.0);
    }

    w_acc = bayla_log_values_sum(n_samples, 0, 1, ws);
    samples.neff = bayla_log_value_map_out_of_log_domain(bayla_log_value_reciprocal(w_acc));

    if(isnan(samples.neff) || isinf(samples.neff)) { samples.valid = false; }

    return samples;
}

API BayLaSamples 
bayla_importance_sample_uniform(BayLaModel const *model, size n_samples, u64 seed, f64 *mins, f64 *maxs, MagAllocator *alloc, MagAllocator scratch_)
{
    MagAllocator *scratch = &scratch_;

    /* Fixed sizes and indexes. */
    i32 const ndim = model->n_parameters;
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
    BayLaMVUniformDist prop_dist = bayla_mv_uniform_dist_create(ndim, mins, maxs, scratch);

    BayLaLogValue *ws = eco_arena_nmalloc(scratch, n_samples, BayLaLogValue);

    for(size s = 0; s < n_samples; ++s)
    {
        f64 *row = &rows[s * ncols];
        BayLaLogValue ld_qi = bayla_mv_uniform_dist_random_deviate(&prop_dist, state, row);
        BayLaLogValue ld_pi = bayla_model_evaluate(model, row);
        BayLaLogValue ld_wi = bayla_log_value_divide(ld_pi, ld_qi);

        if(
                   bayla_log_value_is_nan(ld_qi)
                || bayla_log_value_is_nan(ld_pi)
                || bayla_log_value_is_zero(ld_qi)
                || bayla_log_value_is_zero(ld_pi)
          )
        {
            --s; /* Redo this iteration. */
            continue;
        }

        row[pidx] = ld_pi.val;
        row[qidx] = ld_qi.val;
        row[widx] = ld_wi.val;

        ws[s] = ld_wi;
    }

    /* Calculate evidence */
    BayLaLogValue w_acc = bayla_log_values_sum(n_samples, 0, 1, ws);
    BayLaLogValue l_n_samples = bayla_log_value_map_into_log_domain((f64)n_samples);
    samples.z_evidence = bayla_log_value_divide(w_acc, l_n_samples);

    /* Normalize weights and square them for calculating effective sample size */
    for(size s = 0; s < n_samples; ++s)
    {
        ws[s] = bayla_log_value_divide(ws[s], w_acc);
        ws[s] = bayla_log_value_power(ws[s], 2.0);
    }

    w_acc = bayla_log_values_sum(n_samples, 0, 1, ws);
    samples.neff = bayla_log_value_map_out_of_log_domain(bayla_log_value_reciprocal(w_acc));

    if(isnan(samples.neff) || isinf(samples.neff)) { samples.valid = false; }

    return samples;
}

API BayLaSamples 
bayla_importance_sample_studt_approx(BayLaModel const *model, f64 nu, size n_samples, u64 seed, MagAllocator *alloc, MagAllocator scratch_)
{
    MagAllocator *scratch = &scratch_;

    /* Fixed sizes and indexes. */
    i32 const ndim = model->n_parameters;
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

    BayLaMVStudenTDist prop_dist = bayla_mv_studt_dist_create_from_hessian(ndim, prop_mean, prop_hessian, nu, scratch);
    if(!prop_dist.valid)
    {
        samples.valid = prop_dist.valid;
        return samples;
    }

    BayLaLogValue *ws = eco_arena_nmalloc(scratch, n_samples, BayLaLogValue);

    for(size s = 0; s < n_samples; ++s)
    {
        f64 *row = &rows[s * ncols];
        BayLaLogValue ld_qi = bayla_mv_studt_dist_random_deviate(&prop_dist, state, row);
        BayLaLogValue ld_pi = bayla_model_evaluate(model, row);
        BayLaLogValue ld_wi = bayla_log_value_divide(ld_pi, ld_qi);

        if(
                   bayla_log_value_is_nan(ld_qi)
                || bayla_log_value_is_nan(ld_pi)
                || bayla_log_value_is_zero(ld_qi)
                || bayla_log_value_is_zero(ld_pi)
          )
        {
            --s; /* Redo this iteration. */
            continue;
        }

        row[pidx] = ld_pi.val;
        row[qidx] = ld_qi.val;
        row[widx] = ld_wi.val;

        ws[s] = ld_wi;
    }

    /* Calculate evidence */
    BayLaLogValue w_acc = bayla_log_values_sum(n_samples, 0, 1, ws);
    BayLaLogValue l_n_samples = bayla_log_value_map_into_log_domain((f64)n_samples);
    samples.z_evidence = bayla_log_value_divide(w_acc, l_n_samples);

    /* Normalize weights and square them for calculating effective sample size */
    for(size s = 0; s < n_samples; ++s)
    {
        ws[s] = bayla_log_value_divide(ws[s], w_acc);
        ws[s] = bayla_log_value_power(ws[s], 2.0);
    }

    w_acc = bayla_log_values_sum(n_samples, 0, 1, ws);
    samples.neff = bayla_log_value_map_out_of_log_domain(bayla_log_value_reciprocal(w_acc));

    if(isnan(samples.neff) || isinf(samples.neff)) { samples.valid = false; }

    return samples;
}

API void 
bayla_samples_save_csv(BayLaSamples *samples, BayLaLogValue ci_p_thresh, char const *fname)
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


    f64 max_w = samples->rows[ncols * (n_samples - 1) + widx];

    for(size r = 0; r < n_samples; ++r)
    {
        f64 *row = &samples->rows[r * ncols];
        for(size c = 0; c < ndim; ++c)
        {
            fprintf(f, "%.18e,", row[c]);
        }
        fprintf(f, "%.18e,%.18e,%.18e,%d\n", row[pidx], row[qidx], row[widx] - max_w, row[pidx] >= ci_p_thresh.val ? 1 : 0);
    }

    fclose(f);

    mag_static_arena_destroy(&scratch);
}

API BayLaLogValue 
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
            return bayla_log_value_create(row[pidx]);
        }
    }

    return bayla_log_value_create(NAN);
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
        f64 w = exp(row[widx]);
        f64 df = func(ndim, row, ud) - mean;

        v_acc = elk_kahan_accumulator_add(v_acc, df * df * w);
    }

    f64 std = sqrt(v_acc.sum / total_weight / samples->neff);

    return (BayLaErrorValue){ .val = mean, .std = std};
}

API f64 
bayla_samples_calculate_prob_predicate(BayLaSamples *samples, b32 (*func)(size n, f64 const *params, void *ud), void *ud)
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
        b32 f = func(ndim, row, ud); /* Zero or one */

        f_acc = elk_kahan_accumulator_add(f_acc, f * w);
    }

    f64 prob = f_acc.sum / total_weight;

    return prob;
}

#pragma warning(pop)

#endif
