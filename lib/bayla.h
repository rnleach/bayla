#ifndef _BAYLA_H_
#define _BAYLA_H_

#include "elk.h"
#include "magpie.h"

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
API BayLaSquareMatrix bayla_square_matrix_invert(BayLaSquareMatrix matrix, MagAllocator *alloc, MagStaticArena scratch);
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

/* Calculates the Hessian matrix at a specified position, usually the max-a-posteriori parameter values.
 *
 * The "data" is fixed and stored somewhere in the user_data. When you evaluate this function, it returns Hessian matrix
 * (which is sized num_parms x num_parms). This, along with the max_a_posteriori_parms are used to create a Gaussian
 * approximation to the posterior distribution.
 *
 * Note that this is the same as the inverse of the covariance matrix, and the result is used in 
 * bayla_mv_norm_dist_create_from_inverse.
 */
typedef BayLaSquareMatrix (*BayLaHessianMatrix)(void *user_data, size num_parms, f64 const *max_a_posteriori_parms);

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
API BayLaLogValue bayla_log_values_sum(size nvals, BayLaLogValue *vals);

/* Functions that have less overflow / underflow / rounding errors than naive implementations.  */
API f64 bayla_log1mexp(f64 x);                       /* log(1 - exp(-x))                        */
API f64 bayla_log1pexp(f64 x);                       /* log(1 + exp(+x))                        */
API f64 bayla_logsumexp(f64 x, f64 y);               /* add two numbers in the log domain.      */
API f64 bayla_logsumexp_list(size nvals, f64 *vals); /* Sum a list of values in the log domain. */

BayLaLogValue const bayla_log_zero = { .val = -INFINITY };
BayLaLogValue const bayla_log_one = { .val = 0.0 };

/*---------------------------------------------------------------------------------------------------------------------------
 *                                            Values with an Error Information
 *-------------------------------------------------------------------------------------------------------------------------*/

/* Many values come with an error term.
 *
 * The numerical results of Monte-Carlo integration come with an error term. Also if using the Gaussian approximation near
 * a maximum likelihood value, the second derivative can be evaluated numerically to come up with a sigma value for the
 * distribution.
 *
 * This assumes the errors have a Gaussian distribution and the error term is the standard deviation.
 */
typedef struct
{
    BayLaLogValue value;
    BayLaLogValue error;
} BayLaValue;

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
API BayLaMVNormDist bayla_mv_norm_dist_create_from_inverse(size ndim, f64 *mean, BayLaSquareMatrix inverse_covar, MagAllocator *alloc);

/* returns the probability value of the associated point, workspace needs to be at least as big as dist->ndim */
API BayLaLogValue bayla_mv_norm_dist_random_deviate(BayLaMVNormDist *dist, ElkRandomState *state, f64 *deviate, f64 *workspace);
/* returns the probability value of the input value deviate. workspace needs to be twice dist->ndim */
API BayLaLogValue bayla_mv_norm_dist_prob(BayLaMVNormDist *dist, f64 *deviate, f64 *workspace);

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

    f64 *min_parameter_vals;        /* Array n_parameters long of the minimum value parameters can practically have.      */
    f64 *max_parameter_vals;        /* Array n_parameters long of the miximum value parameters can practically have.      */

    void *user_data;                /* Any other state or data you may need for evaluating the prior, likelihood, etc.    */
} BayLaModel;

API BayLaLogValue bayla_model_evaluate(BayLaModel const *model, f64 const *parameter_vals);
/*---------------------------------------------------------------------------------------------------------------------------
 *                                                  Sampling Statistics
 *-------------------------------------------------------------------------------------------------------------------------*/
API f64 bayla_effective_sample_size(size nvals, BayLaLogValue *ws);

/*---------------------------------------------------------------------------------------------------------------------------
 *                                                    Implementations
 *-------------------------------------------------------------------------------------------------------------------------*/

API BayLaSquareMatrix 
bayla_square_matrix_create(i32 ndim, MagAllocator *alloc)
{
    f64 *data = eco_arena_nmalloc(alloc, ndim * ndim, f64);
    Assert(data);
    return (BayLaSquareMatrix){ .ndim = ndim, .data = data };
}

API BayLaSquareMatrix 
bayla_square_matrix_invert(BayLaSquareMatrix matrix, MagAllocator *alloc, MagStaticArena scratch)
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
                    Assert(false);
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
bayla_log_values_sum(size nvals, BayLaLogValue *vals)
{
    f64 *log_vals = (f64 *)vals;
    f64 log_sum = bayla_logsumexp_list(nvals, log_vals);
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
bayla_logsumexp_list(size nvals, f64 *vals)
{
    /* Use a two pass algorithm, first find the max, then scale and apply addition at higher precision. */
    f64 max = -INFINITY;
    for(size i = 0; i < nvals; ++i)
    {
        if(vals[i] > max) { max = vals[i]; }
    }
    Assert(max > - INFINITY);

    ElkKahanAccumulator result = {0};
    for(size i = 0; i < nvals; ++i)
    {
        result = elk_kahan_accumulator_add(result, exp(vals[i] - max));
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
    MagStaticArena scratch = mag_static_arena_allocate_and_create(ECO_KiB(20));
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

    mag_static_arena_destroy(&scratch);

    return dist;
}

API BayLaMVNormDist 
bayla_mv_norm_dist_create_from_inverse(size ndim, f64 *mean, BayLaSquareMatrix inverse_covar, MagAllocator *alloc)
{
    MagStaticArena scratch = mag_static_arena_allocate_and_create(ECO_KiB(20));
    BayLaMVNormDist dist = { .ndim = ndim, .valid = true };
    f64 *dist_mean = eco_arena_nmalloc(alloc, ndim, f64);
    memcpy(dist_mean, mean, ndim * sizeof(f64));
    dist.mean = dist_mean;

    dist.inv_cov = bayla_square_matrix_create(ndim, alloc);
    memcpy(dist.inv_cov.data, inverse_covar.data, sizeof(f64) * ndim * ndim);

    dist.cov = bayla_square_matrix_invert(dist.inv_cov, alloc, scratch);

    dist.chol = bayla_cholesky_decomp(dist.cov, alloc);
    dist.valid = dist.chol.valid;
    Assert(dist.valid);

    BayLaLogValue covar_det = bayla_log_one;
    for(i32 i = 0; i < dist.ndim; ++i)
    {
        BayLaLogValue diag_val = bayla_log_value_map_into_log_domain(dist.chol.data[MATRIX_IDX(i, i, dist.ndim)]);
        covar_det = bayla_log_value_multiply(covar_det, bayla_log_value_power(diag_val, 2.0));
    }
    dist.normalization_constant = bayla_log_value_reciprocal(bayla_log_value_power(covar_det, 0.5));
    dist.normalization_constant.val += (-0x1.d67f1c864beb4p-1 * ndim); /* Literal of log(1/sqrt(2*pi)) */

    mag_static_arena_destroy(&scratch);

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
bayla_mv_norm_dist_prob(BayLaMVNormDist *dist, f64 *deviate, f64 *workspace)
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

API f64 
bayla_effective_sample_size(size nvals, BayLaLogValue *ws)
{
    BayLaLogValue sum = bayla_log_values_sum(nvals, ws);

    /* Square the values. */
    for(size i = 0; i < nvals; ++i) { ws[i].val *= 2.0; }

    BayLaLogValue sum_sqs = bayla_log_values_sum(nvals, ws);

    /* Undo square the values. */
    for(size i = 0; i < nvals; ++i) { ws[i].val *= 0.5; }

    BayLaLogValue neff = bayla_log_value_divide(bayla_log_value_power(sum, 2.0), sum_sqs);

    return bayla_log_value_map_out_of_log_domain(neff);
}

#pragma warning(pop)

#endif
