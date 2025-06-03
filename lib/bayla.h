#ifndef _BAYLA_H_
#define _BAYLA_H_

#include "elk.h"
#include "magpie.h"
#include "packrat.h"

#pragma warning(push)

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

/*---------------------------------------------------------------------------------------------------------------------------
 *
 *                                                          API
 *
 *-------------------------------------------------------------------------------------------------------------------------*/
#define API static inline

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

/*---------------------------------------------------------------------------------------------------------------------------
 *                                              Analysis of BayLaian Models
 *-------------------------------------------------------------------------------------------------------------------------*/

/* Description of a model to be used for a Baysiean analysis. */
typedef struct
{
    
    f64 model_prior_probability;    /* Prior prob *this model* is the correct one. P(Model | I), Jaynes eq. 20.3 (pg 602) */
    size n_parameters;              /* The number of parameters in this model.                                            */

    BayLaLogPrior prior;                       /* See typedef above. P(Θ | Model, I), Jaynes eq 20.1, 20.2 (pg 602)       */
    BayLaLogLikelihood likelihood;             /* See typedef above. P(Data | Θ, Model, I), Jaynes eq 20.1, 20.2 (pg 602) */

    f64 *min_parameter_vals;        /* Array n_parameters long of the minimum value parameters can practically have.      */
    f64 *max_parameter_vals;        /* Array n_parameters long of the miximum value parameters can practically have.      */

    void *user_data;                /* Any other state or data you may need for evaluating the prior, likelihood, etc.    */
} BayLaModel;

/* Options for the integrator. Not all integrators use all options. */
typedef enum
{
    BAYLA_INTEGRATE_INVALID, 
    BAYLA_INTEGRATE_VEGAS_PLUS,
    BAYLA_INTEGRATE_VEGAS_MISER,
    BAYLA_INTEGRATE_VEGAS, 
    BAYLA_INTEGRATE_MISER,
    BAYLA_INTEGRATE_MONTE_CARLO
} BayLaIntegratorStrategy;

/* A mesh used for importance sampling in VEGAS based methods. */
typedef struct BayLaVegasMap BayLaVegasMap;

typedef struct
{
    BayLaIntegratorStrategy strategy;
    union
    {
        /* Simple Monte Carlo Strategy Options. */
        struct 
        {
            size min_samples;          /* The minimum number of samples to try before stopping.                           */
            size max_samples;          /* Give up after this many calls to the function to be integrated.                 */
            size samples_per_batch;    /* Number of samples between checks against acceptable error.                      */
            f64 acceptable_abs_error;  /* If the absolute error falls below this, stop if reached min_samples too.        */
            f64 acceptable_prop_error; /* If the proportional error falls below this, stop if reached min_samples too.    */
            u64 seed;                  /* Random number generator seed. YOU MUST CHOOSE.                                  */
        } simple_monte_carlo;

        /* MISER Strategy Options. */
        struct
        {
            size min_points;             /* The minimum number of total samples in any sub-region.                        */
            size min_to_subdivide;       /* The minimum number of samples requred to further sub-divide a region.         */
            size total_samples;          /* The total number of samples to take.                                          */
            f64 explore_factor;          /* The fraction of points to use at each recursion to explore the variance.      */
            f64 acceptable_abs_error;    /* If the absolute error falls below this, stop if reached min_samples too.      */
            f64 acceptable_prop_error;   /* If the proportional error falls below this, stop if reached min_samples too.  */
            u64 seed;                    /* Random number generator seed. YOU MUST CHOOSE.                                */
        } miser;

        /* VEGAS Strategy Options. */
        struct
        {
            size samples_per_refinement; /* The number of samples per grid refinement iteration.                          */
            size num_grid_cells;         /* The number of grid cells for dividing up each dimension.                      */
            size min_total_samples;      /* The minium number of total samples to accept the result.                      */
            size max_total_samples;      /* The maxium number of samples to take before giving up.                        */
            f64 acceptable_abs_error;    /* If the absolute error falls below this, stop if reached min_samples too.      */
            f64 acceptable_prop_error;   /* If the proportional error falls below this, stop if reached min_samples too.  */
            f64 alpha;                   /* Used to damp the grid refinements, >=0, usually 0.5 to 2.0. 1.5 good default. */
            u64 seed;                    /* Random number generator seed. YOU MUST CHOOSE.                                */
            BayLaVegasMap *map;          /* NULL. A pre-conditioned mesh created by bayla_vegas_map_precondition.         */
        } vegas;

        /* VEGAS-MISER Strategy Options. */
        struct
        {
            size samples_per_refinement; /* The number of samples per grid refinement iteration.                          */
            size num_grid_cells;         /* The number of grid cells for dividing up each dimension.                      */
            size min_total_samples;      /* The minium number of total samples to accept the result.                      */
            size max_total_samples;      /* The maxium number of samples to take before giving up.                        */
            size min_points;             /* The minimum number of total samples in any sub-region. (MISER)                */
            size min_to_subdivide;       /* The minimum number of samples requred to further sub-divide a region. (MISER) */
            f64 explore_factor;          /* The fraction of points at each recursion used to explore the variance. (MISER)*/
            f64 acceptable_abs_error;    /* If the absolute error falls below this, stop if reached min_samples too.      */
            f64 acceptable_prop_error;   /* If the proportional error falls below this, stop if reached min_samples too.  */
            f64 alpha;                   /* Used to damp the grid refinements, >=0, usually 0.5 to 2.0. 1.5 good default. */
            u64 seed;                    /* Random number generator seed. YOU MUST CHOOSE.                                */
            BayLaVegasMap *map;          /* NULL. A pre-conditioned mesh created by bayla_vegas_map_precondition.         */
        } vegas_miser;

        /* VEGAS Plus Strategy Options. */
        struct
        {
            size samples_per_refinement; /* The number of samples per grid refinement iteration.                          */
            size num_grid_cells;         /* The number of grid cells for dividing up each dimension.                      */
            size min_total_samples;      /* The minium number of total samples to accept the result.                      */
            size max_total_samples;      /* The maxium number of samples to take before giving up.                        */
            f64 acceptable_abs_error;    /* If the absolute error falls below this, stop if reached min_samples too.      */
            f64 acceptable_prop_error;   /* If the proportional error falls below this, stop if reached min_samples too.  */
            f64 alpha;                   /* Used to damp the grid refinements, >=0, usually 0.5 to 2.0. 1.5 good default. */
            f64 beta;                    /* Used to damp the adaptive sampling. β = 1, no damping, β = 0, no adaption.    */
            u64 seed;                    /* Random number generator seed. YOU MUST CHOOSE.                                */
            BayLaVegasMap *map;          /* NULL. A pre-conditioned mesh created by bayla_vegas_map_precondition.         */
        } vegas_plus;
    };
} BayLaIntegratorOptions;

/* Parameter values and log-probability samples from a distribution. */
typedef struct
{
    size n_parameters;             /* Number of parameters per sample, or the number of columsn in the *parameters array. */
    PakArrayLedger samples_ledger; /* Ledger for parallel arrays *parameters and *probs.                                  */
    f64 *parameters;               /* elk_len(&samples_ledger) rows x n_parameters columns array.                         */
    BayLaLogValue *probs;          /* Distribution values, array parallel to parameters.                                  */
} BayLaParametersSamples;

/* Precondition a VEGAS mesh with the provided samples. */
API void bayla_vegas_map_precondition(
        BayLaModel *model,
        BayLaIntegratorOptions *opts,
        BayLaParametersSamples *samples,
        MagStaticArena *perm);

typedef struct
{
    BayLaIntegratorStrategy strategy; /* Should be passed through from the options.                                       */
    BayLaValue result;                /* The numeric value of the integral and its estimated error.                       */
    size total_samples;               /* Basically the number of function calls.                                          */
    size samples_used;                /* The number of samples used for the reportd integral value.(VEGAS and VEGAS+ only)*/
    f64 chisq;                        /* The Chi-Squared statistic. (VEGAS and VEGAS+ only)                               */
    b32 converged;                    /* true if it stopped by reaching an acceptable error or false otherwise.           */
} BayLaEvidenceResult;

/* Calculate the evidence { p(Data | Model, I) } from Jaynes eq 20.2 using Monte-Carlo integration. */
API BayLaEvidenceResult bayla_calculate_evidence(
        BayLaModel const *model,
        BayLaIntegratorOptions *opts,
        MagStaticArena scratch);

/*---------------------------------------------------------------------------------------------------------------------------
 *                                                  Sample Distributions
 *-------------------------------------------------------------------------------------------------------------------------*/

/* TODO: Notes on the implementation, e.g. NUTS or Metropolis-Hastings. */
typedef struct
{
    size n_samples;                  /* The number of samples to take.                                                    */
    size n_burn_in_samples;          /* The number of samples to take before starting to keep them.                       */
    f64 *starting_params;            /* The starting values of the params. This array should be model->n_parameters long. */
    BayLaModel *model;               /* The model posterior to sample from.                                               */
    BayLaParametersSamples *output;  /* Output, recommend using bayla_preconditioning_sampler_args_allocate_and_create.   */
} BayLaPreconditioningSamplerArgs;

API BayLaPreconditioningSamplerArgs bayla_preconditioning_sampler_args_allocate_and_create(
        BayLaModel *model,
        size n_samples,
        size n_burn_in_samples,
        f64 *starting_params,
        MagStaticArena *arena);

API void bayla_preconditioning_sample(
        BayLaPreconditioningSamplerArgs *args,
        ElkRandomState *state,
        MagStaticArena scratch);

/*---------------------------------------------------------------------------------------------------------------------------
 *                                                    Implementations
 *-------------------------------------------------------------------------------------------------------------------------*/

API BayLaLogValue bayla_log_value_create(f64 log_domain_value) { return (BayLaLogValue){.val = log_domain_value }; }
API f64 bayla_log_value_map_out_of_log_domain(BayLaLogValue log_val) { return exp(log_val.val); }

API BayLaLogValue 
bayla_log_value_map_into_log_domain(f64 non_log_domain_value)
{
    return (BayLaLogValue){.val = log(non_log_domain_value) };
}


API BayLaLogValue bayla_log_value_add(BayLaLogValue left, BayLaLogValue right)
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

API BayLaLogValue bayla_log_value_multiply(BayLaLogValue left, BayLaLogValue right)
{
    return (BayLaLogValue){.val = left.val + right.val};
}

API BayLaLogValue 
bayla_log_value_divide(BayLaLogValue numerator, BayLaLogValue denominator)
{
    return (BayLaLogValue){.val = numerator.val - denominator.val};
}

API BayLaLogValue bayla_log_value_reciprocal(BayLaLogValue log_val)
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
    /* Use a two pass algorithm */
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

        f64 y = fabs(v) + 0.3386595;

        q = x * x + y * (0.19600 * y - 0.25472 * x);

    } while(q > 0.27597 && (q > 0.27846 || (v * v > -4.0 * log(u) * u * u)));

    return dist->mu + dist->sigma * v / u;
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


API BayLaPreconditioningSamplerArgs 
bayla_preconditioning_sampler_args_allocate_and_create(
        BayLaModel *model,
        size n_samples,
        size n_burn_in_samples,
        f64 *starting_params,
        MagStaticArena *arena)
{
    size N = model->n_parameters;

    /* Need to allocate space for output, everything else should already exist. */
    BayLaParametersSamples *samples = eco_malloc(arena, BayLaParametersSamples);
    f64 *parameter_samples = eco_nmalloc(arena, N * n_samples, f64);
    BayLaLogValue *probs = eco_nmalloc(arena, n_samples, BayLaLogValue);
    Assert(probs && parameter_samples && samples);

    PakArrayLedger ledger = pak_array_ledger_create(n_samples);

    *samples = (BayLaParametersSamples)
        {
            .n_parameters = N,
            .samples_ledger = ledger,
            .parameters = parameter_samples,
            .probs = probs
        };

    return (BayLaPreconditioningSamplerArgs)
        {
            .n_samples = n_samples,
            .n_burn_in_samples = n_burn_in_samples,
            .starting_params = starting_params,
            .model = model,
            .output = samples
        };
}

API void 
bayla_preconditioning_sample(BayLaPreconditioningSamplerArgs *args, ElkRandomState *state, MagStaticArena scratch)
{
    /* This is CRUDELY based on Metropolis-Hastings algorithm. In that algorithm you would add a point to the chain again
     * if you chose not to move. But this algorithm only adds points after you have decided to move to them, but doesn't
     * add them again. This has shown in tests to produce a set of samples that is better for preconditioning a VEGAS grid.
     */

    /* unpack the arguments */
    size N = args->model->n_parameters;
    f64 const *starting_params = args->starting_params;

    BayLaParametersSamples *samples = args->output;
    PakArrayLedger *out_ledger = &samples->samples_ledger;
    f64 *out_params = samples->parameters;
    BayLaLogValue *probs = samples->probs;

    BayLaModel const *model = args->model;
    f64 *min_parms = args->model->min_parameter_vals;
    f64 *max_parms = args->model->max_parameter_vals;

    f64 *ms = eco_nmalloc(&scratch, N, f64);
    f64 *ss = eco_nmalloc(&scratch, N, f64);
    f64 *p0 = eco_nmalloc(&scratch, N, f64);
    f64 *p = eco_nmalloc(&scratch, N, f64);

    Assert(ms && ss && p0 && p);
    memcpy(ms, starting_params, sizeof(f64) * N);

    BayLaNormalDistribution norm = bayla_normal_distribution_create(0.0, 1.0);

    /* Copy the starting point into the list of parameters and calculate the first value. */
    memcpy(p, starting_params, sizeof(f64) * N);
    memcpy(p0, starting_params, sizeof(f64) * N);
    BayLaLogValue prob = bayla_model_evaluate(model, p);
    BayLaLogValue prob0 = prob;

    size n_keep = 0;
    size n_reject = 0;

    /* That should fill the 0th position, so bump up the ledger and remember the postition of the last parameters. */
    size next = pak_array_ledger_push_back_index(out_ledger);

    while(true)
    {
        if(next == PAK_COLLECTION_FULL) { break; } /* Stop sampling. */

        /* Get some values for this position. */
        do
        {
            for(size n = 0; n < N; ++n)
            {
                if(n_keep > 1)
                {
                    f64 sigma = sqrt(ss[n]/(n_keep - 1));
                    if(n_reject > 5) { sigma /= n_reject; } /* Shrink the sigma if we're on a max and can't get off. */
                    norm.sigma = sigma;
                }
                else
                {
                    f64 sigma = fabs(p0[n] / (n_reject + 1));
                    sigma = sigma < 1e-10 ? 1.0 : sigma;
                    Assert(!isnan(sigma) && !isinf(sigma));
                    norm.sigma = sigma;
                }
                do
                {
                    p[n] = p0[n] + bayla_normal_distribution_random_deviate(&norm, state);
                } while(p[n] < min_parms[n] || p[n] > max_parms[n]);
            }

            /* Calculate the probability for this position. */
            prob = bayla_model_evaluate(model, p);
        } while(isnan(prob.val));

        /* Decide whether to move to this point or not. */
        f64 log_prob_ratio = bayla_log_value_divide(prob, prob0).val;
        b32 move = true;
        if(log_prob_ratio < 0.0)
        {
            f64 ran_val = elk_random_state_uniform_f64(state);
            if(ran_val >= exp(log_prob_ratio)) { move = false; }
        }

        if(move)
        {
            if(n_keep > args->n_burn_in_samples)
            {
                /* Save output */
                f64 *dest = out_params + next * N;
                memcpy(dest, p, sizeof(f64) * N);
                probs[next] = prob;

                /* Get the next position in the output arrays. */
                next = pak_array_ledger_push_back_index(out_ledger);
            }

            /* Prepare for next step. */
            f64 *temp = p0;
            p0 = p;
            p = temp;
            prob0 = prob;
            n_reject = 0;
            n_keep++;
            for(size n = 0; n < N; ++n)
            {
                f64 mkm1 = ms[n];
                ms[n] = mkm1 + (p0[n] - mkm1) / n_keep;
                ss[n] = ss[n] + (p0[n] - mkm1) * (p0[n] - ms[n]);
            }
        }
        else
        {
            n_reject++;
        }
    }
}

API BayLaEvidenceResult
bayla_simple_monte_carlo_integrate(BayLaModel const *model, BayLaIntegratorOptions const *opts, MagStaticArena scratch)
{
#pragma warning(push)
    /* Complains about possible loss of precision when converting 'size' to 'f64'. Any value that would cause a problem
     * is unrealistically large and would cause the computer to run for far too long. (n_samples)
     */
#pragma warning(disable:4244)
    size min_samples = opts->simple_monte_carlo.min_samples;
    size max_samples = opts->simple_monte_carlo.max_samples;
    size samples_per_batch = opts->simple_monte_carlo.samples_per_batch;
    u64 seed = opts->simple_monte_carlo.seed;
    BayLaLogValue log_acceptable_abs_error = bayla_log_value_map_into_log_domain(opts->simple_monte_carlo.acceptable_abs_error);
    f64 acceptable_prop_error = opts->simple_monte_carlo.acceptable_prop_error;

    size n_params = model->n_parameters;
    f64 *min_vals = model->min_parameter_vals;
    f64 *max_vals = model->max_parameter_vals;

    ElkRandomState random_state_ = elk_random_state_create(seed);
    ElkRandomState *ran = &random_state_;

    f64 volume = 1.0;
    for(size i = 0; i < n_params; ++i) 
    {
        volume *= fabs(max_vals[i] - min_vals[i]);
    }
    BayLaLogValue log_volume = bayla_log_value_map_into_log_domain(volume);

    size n_samples = 0;
    BayLaLogValue log_sum = bayla_log_value_map_into_log_domain(0);
    BayLaLogValue log_sumsq = bayla_log_value_map_into_log_domain(0);

    BayLaLogValue log_integral = bayla_log_value_map_into_log_domain(0);
    BayLaLogValue log_integral_error = bayla_log_value_map_into_log_domain(0);

    f64 *params = eco_nmalloc(&scratch, n_params, f64);
    PanicIf(!params);

    f64 proportional_error = INFINITY;

    b32 stop = false;
    while(!stop)
    {
        for(size i = 0; i < samples_per_batch; ++i)
        {
            /* generate a random point in parameter space */
            for(size j = 0; j < n_params; ++j)
            {
                params[j] = elk_random_state_uniform_f64(ran) * (max_vals[j] - min_vals[j]) + min_vals[j];
            }

            /* evaluate the sum of the log prior & log likelihood */
            BayLaLogValue point_val = bayla_model_evaluate(model, params);
            Assert(!isnan(point_val.val));

            /*  update the sum & sumsq accumulators */
            log_sum = bayla_log_value_add(log_sum, point_val);
            log_sumsq = bayla_log_value_add(log_sumsq, bayla_log_value_power(point_val, 2.0));
        }

        n_samples += samples_per_batch;
        BayLaLogValue log_n_samples = bayla_log_value_map_into_log_domain(n_samples);

        BayLaLogValue log_mean = bayla_log_value_divide(log_sum, log_n_samples); 
        BayLaLogValue log_sqr_mean = bayla_log_value_divide(log_sumsq, log_n_samples); 

        BayLaLogValue log_mean2 = bayla_log_value_power(log_mean, 2.0);
        BayLaLogValue log_var = bayla_log_value_divide(bayla_log_value_subtract(log_sqr_mean, log_mean2), log_n_samples);

        log_integral_error = bayla_log_value_multiply(log_volume, bayla_log_value_power(log_var, 0.5));
        log_integral = bayla_log_value_multiply(log_volume, log_mean);

        proportional_error = bayla_log_value_map_out_of_log_domain(bayla_log_value_divide(log_integral_error, log_integral));

        /* Stopping conditions. */
        stop = n_samples >= max_samples || 
            log_integral_error.val < log_acceptable_abs_error.val || proportional_error <= acceptable_prop_error;

        stop = stop && n_samples > min_samples;
    }

    b32 achieved_acceptable_error = log_integral_error.val < log_acceptable_abs_error.val || proportional_error <= acceptable_prop_error;

    BayLaEvidenceResult retval = {0};
    retval.strategy = BAYLA_INTEGRATE_MONTE_CARLO;
    retval.result = (BayLaValue){.value = log_integral, .error = log_integral_error };
    retval.converged = achieved_acceptable_error;
    retval.total_samples = n_samples;

    return retval;
#pragma warning(pop)
}

/* Special type for type punning so we can map between log-domain and linear-domain in the same memory location safely. */
typedef union
{
    f64 linear_domain;
    BayLaLogValue log_domain;
} BayLaDualDomainValue;

_Static_assert(sizeof(BayLaDualDomainValue) == sizeof(f64) && sizeof(f64) == sizeof(BayLaLogValue),
        "These should be wrappers for f64 (double)");

struct BayLaVegasMap
{
    /* Configuration related. */
    f64 alpha;
    size ndim;
    size n_grid;
    f64 const *max_vals;
    f64 const *min_vals;

    /* Workspace. */
    size *n_s;
    size *iy; 
    f64 *xs;
    f64 *delta_xs;
    BayLaDualDomainValue *ds;
};

API BayLaVegasMap *
bayla_vegas_map_create(
        size ndim,
        size n_grid,
        f64 alpha,
        f64 const *max_vals,
        f64 const *min_vals,
        MagStaticArena *perm)
{
    BayLaVegasMap *map = eco_malloc(perm, BayLaVegasMap);
    
    map->ndim = ndim;
    map->n_grid = n_grid;
    map->alpha = alpha;
    map->max_vals = max_vals;
    map->min_vals = min_vals;

    map->n_s = eco_nmalloc(perm, ndim * n_grid, size);
    map->iy = eco_nmalloc(perm, ndim, size);

    map->xs = eco_nmalloc(perm, ndim * (n_grid + 1), f64);
    map->delta_xs = eco_nmalloc(perm, ndim * (n_grid + 1), f64);

    map->ds = eco_nmalloc(perm, ndim * n_grid, BayLaDualDomainValue);

    Assert(map->n_s && map->iy && map->xs && map->delta_xs && map->ds);

    /* Initialize the xs & delta_xs to uniform spacing */
    for(size d = 0; d < ndim; ++d)
    {
        f64 min = min_vals[d];
        f64 max = max_vals[d];
        f64 stepsize = (max - min) / n_grid;
        for(size g = 0; g < n_grid; ++g)
        {
            map->delta_xs[d * n_grid + g] = stepsize;
            map->xs[d * (n_grid + 1) + g] = min + g * stepsize;
        }
        map->xs[d * (n_grid + 1) + n_grid] = max;
    }

    /* Make sure these are initialized to zero so we are ready for the first iteration. */
    for(size i = 0; i < ndim * n_grid; ++i) { map->ds[i].log_domain = bayla_log_zero; }
    memset(map->n_s, 0, sizeof(*map->n_s) * ndim * n_grid);

    return map;
}

/* Don't put this in the main API because it's not intended for regular use yet, it's just intended for use in tests
 * to check the consistency of multiple runs of the post preconditioning integrator.
 */
API BayLaVegasMap *
bayla_vegas_map_deep_copy(BayLaVegasMap *src, MagStaticArena *perm)
{
    BayLaVegasMap *map = eco_malloc(perm, BayLaVegasMap);
    
    map->ndim = src->ndim;
    map->n_grid = src->n_grid;
    map->alpha = src->alpha;
    map->max_vals = src->max_vals;
    map->min_vals = src->min_vals;

    map->n_s = eco_nmalloc(perm, map->ndim * map->n_grid, size);
    memcpy(map->n_s, src->n_s, map->ndim * map->n_grid * sizeof(size));
    map->iy = eco_nmalloc(perm, map->ndim, size);
    memcpy(map->iy, src->iy, map->ndim * sizeof(size));

    map->xs = eco_nmalloc(perm, map->ndim * (map->n_grid + 1), f64);
    memcpy(map->xs, src->xs, map->ndim * (map->n_grid + 1) * sizeof(f64));
    map->delta_xs = eco_nmalloc(perm, map->ndim * (map->n_grid + 1), f64);
    memcpy(map->delta_xs, src->delta_xs, map->ndim * (map->n_grid + 1) * sizeof(f64));

    map->ds = eco_nmalloc(perm, map->ndim * map->n_grid, BayLaDualDomainValue);
    memcpy(map->ds, src->ds, map->ndim * map->n_grid * sizeof(BayLaDualDomainValue));

    Assert(map->n_s && map->iy && map->xs && map->delta_xs && map->ds);

    return map;
}

API void
bayla_vegas_map_evaluate(
        BayLaModel const *model,
        BayLaVegasMap *map,
        f64 *y,
        BayLaLogValue *out_val,
        BayLaLogValue *out_val_sq,
        MagStaticArena scratch)
{
    /* Take input y and calculate fx*jy and return that. Also accumulate necessary information in state variables for the
     * next call to refine the mesh.
     */
    f64 *x = eco_nmalloc(&scratch, map->ndim, f64);
    f64 jy = 1.0;

    /* Transform from y-space to x-space and calculate the Jacobian */
    for(size d = 0; d < map->ndim; ++d)
    {
        Assert(y[d] >= 0.0 && y[d] <= 1.0);

        /* Calculate iy, eq 10 */
        map->iy[d] = (size)floor(y[d] * map->n_grid);
        Assert(map->iy[d] >= 0 && map->iy[d] < map->n_grid);

        /* Calculate dy, eq 11 */
        f64 dy = y[d] * map->n_grid - map->iy[d];

        /* Map those coordinates from y-space to x-space, eq 9 */
        x[d] = map->xs[d * (map->n_grid + 1) + map->iy[d]] + map->delta_xs[d * map->n_grid + map->iy[d]] * dy;
        Assert(x[d] >= map->min_vals[d] && x[d] <= map->max_vals[d]);

        /* Calculate jy, eq 12 */
        jy *= map->n_grid * map->delta_xs[d * map->n_grid + map->iy[d]];
    }

    /* Calculate the log-probability. */
    BayLaLogValue log_fx = bayla_model_evaluate(model, x);
    Assert(!isnan(log_fx.val));
    BayLaLogValue log_jy = bayla_log_value_map_into_log_domain(jy);
    BayLaLogValue log_jy_fx = bayla_log_value_multiply(log_fx, log_jy);
    Assert(!isnan(log_jy_fx.val));

    /* Return (jy * fx) and its square */
    *out_val = log_jy_fx;
    *out_val_sq = bayla_log_value_power(log_jy_fx, 2.0);

    /* For each dimension (D), accumulate (jy * fx)**2 in ds[D][iy], part of eq 17 */
    for(size d = 0; d < map->ndim; ++d)
    {
        map->ds[d * map->n_grid + map->iy[d]].log_domain = 
            bayla_log_value_add(map->ds[d * map->n_grid + map->iy[d]].log_domain, *out_val_sq);
        map->n_s[d * map->n_grid + map->iy[d]]++;
    }
}

#define TINY 1.0e-100

API void
bayla_vegas_map_apply_sample(BayLaVegasMap *map, f64 *xs_input, BayLaLogValue log_prob)
{
    size const ndim = map->ndim;
    size const n_grid = map->n_grid;
    f64 const *max_vals = map->max_vals;
    f64 const *min_vals = map->min_vals;
    
    size *iy = map->iy;
    size *n_s = map->n_s;
    f64 *xs = map->xs;
    f64 *delta_xs = map->delta_xs;
    BayLaLogValue *log_ds = (BayLaLogValue*)map->ds;

    f64 jy = 1.0;

    /* Calculate the ds by converting the x-coords (parameters) into y-space.
     *
     * Need iy[d]
     *
     * Do a linear search of the xs[d*(n_grid + 1) + %] w/ % from 1 to n_grid + 1. The index of the first value 
     * greater than x, subtract 1. That's iy for the given d. Accumulate jy during this step.
     */
    for(size d = 0; d < ndim; ++d)
    {
        iy[d] = -1; /* Sentinel value to check for errors. */
        f64 param_x = xs_input[d];
        Assert(param_x >= min_vals[d] && param_x <= max_vals[d]);

        for(size ix = 1; ix <= n_grid; ++ix) /* Yes, indexes from 1 to ngrid */
        {
            f64 grid_x = xs[d * (n_grid + 1) + ix];

            if(grid_x > param_x)
            {
                iy[d] = ix - 1;
                break;
            }
        }
        Assert(iy[d] >= 0 && iy[d] < n_grid); /* Check the sentinel value. */

        jy *= n_grid * delta_xs[d * n_grid + iy[d]];
    }

    /* With jy and iy use the provided fx to calculate log_fx_jy_sq and add up log_ds and n_s. */
    BayLaLogValue log_jy = bayla_log_value_map_into_log_domain(jy);
    BayLaLogValue log_jy_fx = bayla_log_value_multiply(log_prob, log_jy);
    BayLaLogValue log_fx_jy_sq = bayla_log_value_power(log_jy_fx, 2.0);
    Assert(!isnan(log_fx_jy_sq.val));

    /* For each dimension (D), accumulate (jy * fx)**2 in ds[D][iy], part of eq 17 */
    for(size d = 0; d < ndim; ++d)
    {
        log_ds[d * n_grid + iy[d]] = bayla_log_value_add(log_ds[d * n_grid + iy[d]], log_fx_jy_sq);
        Assert(!isnan(log_ds[d * n_grid + iy[d]].val));
        n_s[d * n_grid + iy[d]]++;
    }
}

API void
bayla_vegas_map_refine(BayLaVegasMap *map, MagStaticArena scratch)
{
#pragma warning(push)
#pragma warning(disable: 4244)
    size const ndim = map->ndim;
    size const n_grid = map->n_grid;
    f64 const alpha = map->alpha;
    f64 const *max_vals = map->max_vals;
    f64 const *min_vals = map->min_vals;

    f64 *xs = map->xs;
    f64 *delta_xs = map->delta_xs;
    size *n_s = map->n_s;
    BayLaDualDomainValue *log_ds = map->ds;

    BayLaLogValue const log_tiny = bayla_log_value_map_into_log_domain(TINY);

    /* Normalize the ds, part of eq 17 */
    for(size d = 0; d < ndim; ++d)
    {
        for(size g = 0; g < n_grid; ++g)
        {
            /* MSVC complains about converting n_s[...] to from size to f64. Never going to be a problem with practically
             * possible sample sizes. Disabled warning 4244.
             */
            if(n_s[d * n_grid + g] > 0) { log_ds[d * n_grid + g].log_domain.val -= log(n_s[d * n_grid + g]); }
            if(log_ds[d * n_grid + g].log_domain.val < log_tiny.val) { log_ds[d * n_grid + g].log_domain = log_tiny; }
        }
    }

    /* Transform into the linear domain. */
    for(size i = 0; i < ndim * n_grid; ++i)
    {
        log_ds[i].linear_domain = bayla_log_value_map_out_of_log_domain(log_ds[i].log_domain);
        Assert(!isnan(log_ds[i].linear_domain));
    }

    f64 *ds = (f64*)log_ds; /* Convenience after conversion to the linear domain. */

    f64 *old_ds = eco_nmalloc(&scratch, n_grid, f64);
    f64 *old_xs = eco_nmalloc(&scratch, n_grid + 1, f64);

    PanicIf(!(old_ds && old_xs));

    for(size d = 0; d < ndim; ++d)
    {
        /* Copy the old values so I have them to reference while making the new ones. */
        memcpy(old_ds, &ds[d * n_grid + 0], sizeof(*old_ds) * n_grid);
        memcpy(old_xs, &xs[d * (n_grid + 1) + 0], sizeof(*old_xs) * (n_grid + 1));

        Assert(old_xs[0] == min_vals[d]);
        Assert(old_xs[n_grid] == max_vals[d]);

        /* ---------- Equation 18 --------- */                                                                /*   Smooth */
        f64 sum_d = 0.0;                                                                                      /*    the   */
        for(size g = 0; g < n_grid; ++g)                                                                      /*    ds    */
        {                                                                                                     /*     |    */
            sum_d += old_ds[g];                                                                               /*     |    */
        }                                                                                                     /*     |    */
        f64 recip_sum_d = 1.0 / sum_d;                                                                        /*     |    */
                                                                                                              /*     |    */
        ds[d * n_grid + 0] = recip_sum_d * (7.0 * old_ds[0] + old_ds[1]) / 8.0;                               /*     |    */
        for(size g = 1; g < n_grid - 1; ++g)                                                                  /*     |    */
        {                                                                                                     /*     |    */
            ds[d * n_grid + g] = recip_sum_d * (old_ds[g - 1] + 6.0 * old_ds[g] + old_ds[g + 1]) / 8.0;       /*     |    */
        }                                                                                                     /*     |    */
        ds[d * n_grid + n_grid - 1] = recip_sum_d * (old_ds[n_grid - 2] + 7.0 * old_ds[n_grid - 1]) / 8.0;    /*     |    */
                                                                                                              /*     |    */
        /* ---------- Equation 19 --------- */                                                                /*     |    */
        for(size g = 0; g < n_grid; ++g)                                                                      /*     |    */
        {                                                                                                     /*     |    */
            if(ds[d * n_grid + g] < TINY) { ds[d * n_grid + g] = TINY; }                                      /*     |    */
            ds[d * n_grid + g] = pow((1.0 - ds[d * n_grid + g]) / log(1.0 / ds[d * n_grid + g]), alpha);      /*     |    */
            if(ds[d * n_grid + g] < TINY) { ds[d * n_grid + g] = TINY; }                                      /*     |    */
        }                                                                                                     /*     |    */
                                                                                                              /*----------*/
        /* ---------- Equation 20 --------- */
        sum_d = 0.0;
        for(size g = 0; g < n_grid; ++g)
        {
            sum_d += ds[d * n_grid + g];
        }
        f64 dd = sum_d / n_grid;

        /* ---------- Equations 21 -------- */                                                            /* Re-partition */
        f64 xnp = max_vals[d];                                                                            /*      the     */
        size gi = 0;                                                                                      /*      xs      */
        size gj = 0;                                                                                      /*       |      */
        f64 sd = 0.0;                                                                                     /*       |      */
                                                                                                          /*       |      */
        while(++gi < n_grid) /* step 2 */                                                                 /*       |      */
        {                                                                                                 /*       |      */
            while(sd < dd)                                                                                /*       |      */
            {                                                                                             /*       |      */
                sd += ds[d * n_grid + gj++]; /* step 3 */                                                 /*       |      */
            }                                                                                             /*       |      */
                                                                                                          /*       |      */
            sd -= dd;  /* step 4 */                                                                       /*       |      */
                                                                                                          /*       |      */
            /* ---------- Equation 22 --------- */
            xs[d * (n_grid + 1) + gi] = old_xs[gj] - sd / ds[d * n_grid + gj - 1] * delta_xs[d * n_grid + gj - 1];
            Assert(!isnan(xs[d * (n_grid + 1) + gi]));
        }
        Assert(xs[d *(n_grid + 1) + 0] == min_vals[d]);                                                   /*       |      */
        Assert(xs[d *(n_grid + 1) + n_grid] == max_vals[d]);                                              /*       |      */
                                                                                                          /*       |      */
        for(size g = 0; g < n_grid - 1; ++g)                                                              /*       |      */
        {                                                                                                 /*       |      */
            delta_xs[d * n_grid + g] = xs[d * (n_grid + 1) + g + 1] - xs[d * (n_grid + 1) + g];           /*       |      */
        }                                                                                                 /*       |      */
        delta_xs[d * n_grid + n_grid - 1] = xnp - xs[d * (n_grid + 1) + n_grid - 1];                      /*--------------*/
    }

    /* Reset the accumulators for the next iteration of the grid refinement. */
    for(size i = 0; i < ndim * n_grid; ++i) { log_ds[i].log_domain = bayla_log_zero; }
    memset(n_s, 0, sizeof(*n_s) * ndim * n_grid);
#pragma warning(pop)
}

API void 
bayla_vegas_map_precondition(
        BayLaModel *model,
        BayLaIntegratorOptions *opts,
        BayLaParametersSamples *samples,
        MagStaticArena *perm)
{
    size ndim = model->n_parameters;
    size nsamples = pak_len(&samples->samples_ledger);
    Assert(ndim == samples->n_parameters);
    f64 *min_vals = model->min_parameter_vals;
    f64 *max_vals = model->max_parameter_vals;

    BayLaVegasMap **mapp = NULL;
    size n_grid = 0;
    f64 alpha = 0.0;
    switch(opts->strategy)
    {
        case BAYLA_INTEGRATE_VEGAS:
        {
            n_grid = opts->vegas.num_grid_cells;
            alpha = opts->vegas.alpha;
            mapp = &opts->vegas.map;
        } break;

        case BAYLA_INTEGRATE_VEGAS_PLUS:
        {
            n_grid = opts->vegas_plus.num_grid_cells;
            alpha = opts->vegas_plus.alpha;
            mapp = &opts->vegas_plus.map;
        } break;

        case BAYLA_INTEGRATE_VEGAS_MISER:
        {
            n_grid = opts->vegas_miser.num_grid_cells;
            alpha = opts->vegas_miser.alpha;
            mapp = &opts->vegas_miser.map;
        } break;

        default: Panic(); /* Should never get here. */
    }

    BayLaVegasMap *map = bayla_vegas_map_create(ndim, n_grid, alpha, max_vals, min_vals, perm);
    *mapp = map; /* Ties this into options. Convoluted, I know. */

    f64 *delta_xs = map->delta_xs;

    /* From here on out use the remaining space of perm as a scratch arena, but I can NEVER use perm again! */
    MagStaticArena scratch_ = *perm;
    MagStaticArena *scratch = &scratch_;
    perm = NULL; /* CANNOT USE perm AGAIN FOR PERMANENT STORAGE. */

    size *iy = map->iy;
    size *n_s = map->n_s;
    f64 *ds = (f64*)map->ds;
    BayLaLogValue *log_ds = (BayLaLogValue*)map->ds;

    Assert(iy && n_s && ds && log_ds);

    /* Some stuff to help decide when to stop. */
    size const max_trials = 1000;
    f64 const stopping_max_delta_delta_x_ratio = 1.0e-2;
    f64 *previous_delta_xs = eco_nmalloc(scratch, ndim * n_grid, f64);
    Assert(previous_delta_xs);

    f64 iter_alpha = alpha;
    for(size trials = 0; trials < max_trials; ++trials)
    {
        for(size s = 0; s < nsamples; ++s)
        {
            BayLaLogValue log_prob = samples->probs[s];
            f64 *xs_input = &samples->parameters[s * ndim];

            bayla_vegas_map_apply_sample(map, xs_input, log_prob);
        }

        memcpy(previous_delta_xs, delta_xs, sizeof(f64) * ndim * n_grid);
        bayla_vegas_map_refine(map, scratch_);

        /* Check early stopping conditions. */
        f64 max_delta_delta_x_ratio = -INFINITY;
        for(size i = 0; i < ndim * n_grid; ++i)
        {
            f64 delta_delta_x_ratio = fabs((delta_xs[i] - previous_delta_xs[i]) / previous_delta_xs[i]);
            if(max_delta_delta_x_ratio < delta_delta_x_ratio) { max_delta_delta_x_ratio = delta_delta_x_ratio; }
        }

        if(max_delta_delta_x_ratio <= stopping_max_delta_delta_x_ratio && trials > 4)
        {
            break;
        }
        else
        {
            if(max_delta_delta_x_ratio > 0.99)
            {
                iter_alpha -= 0.05 * (iter_alpha - 0.1); /* Large changes, move towards damped adjustments.              */
            }
            else
            {
                iter_alpha -= 0.05 * (iter_alpha - 1.0); /* Smaller adjustments, start accelerating the adjustment rate. */
            }
        }
    }
}


API size
bayla_integer_exp(size base, u8 exp)
{
    size result = 1;
    while(true)
    {
        if(exp & 1) { result *= base; }
        exp >>= 1;
        if(!exp) break;
        base *= base;
    }

    return result;
}

typedef struct
{
    size ndim;
    size n_samples;
    size n_strat;
    size n_samples_per_hypercube;
    size n_total;
    size n_hypercubes;
} BayLaVegasStratificationParams;

API BayLaVegasStratificationParams
bayla_vegas_stratification_params_create(size const samples_per_iter, size const ndim)
{
#pragma warning(push)
    Assert(ndim < UINT8_MAX);
#pragma warning(disable: 4244)
    size n_strat = (size)floor(pow(samples_per_iter / 2, 1.0 / ndim));
    size n_hypercubes = bayla_integer_exp(n_strat, ndim);
    size n_samples_per_hypercube = samples_per_iter / n_hypercubes;

    return (BayLaVegasStratificationParams)
        {
            .ndim = ndim,
            .n_samples = samples_per_iter,
            .n_strat = n_strat,
            .n_samples_per_hypercube = n_samples_per_hypercube,
            .n_total = n_samples_per_hypercube * n_hypercubes,
            .n_hypercubes = n_hypercubes
        };
#pragma warning(pop)
}

API void
bayla_vegas_stratification_sample_num_to_y(
        size const sn, 
        BayLaVegasStratificationParams const *strat_params,
        ElkRandomState *ran,
        f64 *y)
{

#pragma warning(push)
    Assert(strat_params->ndim <= UINT8_MAX);
#pragma warning(disable:4244)

    if(sn >= strat_params->n_total)
    {
        /* We've moved beyond stratified sampling, just sample the whole domain as in regular Monte Carlo. */
        for(size d = 0; d < strat_params->ndim; ++d)
        {
            y[d] = elk_random_state_uniform_f64(ran);
        }
    }
    else
    {
        /* We're still in the range of n_total - the number of points that can evenly be divided among the hypercubes. */

        /* Select which hypercube */
        size index = sn % strat_params->n_hypercubes;

        /* Convert list index to ndim dimensional coords. */
        for(size d = 0; d < strat_params->ndim; ++d)
        {
            size factor = bayla_integer_exp(strat_params->n_strat, strat_params->ndim - d - 1);
            size coords_of_hypercube = index / factor;
            index -= coords_of_hypercube * factor;
            f64 random = elk_random_state_uniform_f64(ran);

            y[d] = (coords_of_hypercube + random) / strat_params->n_strat;
        }
    }
#pragma warning(pop)
}

API void
bayla_vegas_accumulate_results(
        BayLaEvidenceResult *retval,
        size n_iter,
        size *n_samples,
        BayLaLogValue *log_integrals,
        BayLaLogValue *log_sigmas)
{
    BayLaLogValue log_sum_i = bayla_log_value_map_into_log_domain(0.0);
    BayLaLogValue log_sum_var_recip = bayla_log_value_map_into_log_domain(0.0);

    retval->samples_used = 0;
    for(size i = n_iter - 1; i >= 0; --i)
    {
        size num_values = n_iter - i;

        BayLaLogValue log_sigma_sq = bayla_log_value_power(log_sigmas[i], 2.0);
        log_sum_i = bayla_log_value_add(log_sum_i, bayla_log_value_divide(log_integrals[i], log_sigma_sq));
        log_sum_var_recip = bayla_log_value_add(log_sum_var_recip, bayla_log_value_reciprocal(log_sigma_sq));

        BayLaLogValue log_integral = bayla_log_value_divide(log_sum_i, log_sum_var_recip);
        BayLaLogValue log_error = bayla_log_value_reciprocal(bayla_log_value_power(log_sum_var_recip, 0.5));

        BayLaLogValue log_chi2 = bayla_log_value_map_into_log_domain(0.0);
        for(size j = n_iter - 1; j >= i; --j)
        {
            BayLaLogValue log_var = bayla_log_value_power(log_sigmas[j], 2.0);
            BayLaLogValue log_integral_i = log_integrals[j];

            BayLaLogValue temp;
            if(log_integral_i.val > log_integral.val)
            {
                temp = bayla_log_value_power(bayla_log_value_subtract(log_integral_i, log_integral), 2.0);
            }
            else
            {
                temp = bayla_log_value_power(bayla_log_value_subtract(log_integral, log_integral_i), 2.0);
            }

            temp = bayla_log_value_divide(temp, log_var);
            log_chi2 = bayla_log_value_add(log_chi2, temp);
            Assert(!isnan(log_chi2.val));
        }

        f64 chi2 = bayla_log_value_map_out_of_log_domain(log_chi2);

        if(chi2 > num_values && num_values > 1) { break; }

        retval->result = (BayLaValue){.value = log_integral, .error = log_error};
        retval->chisq = chi2;
        
        switch(retval->strategy)
        {
        case BAYLA_INTEGRATE_VEGAS:
            {
                retval->samples_used += *n_samples;
            } break;

        case BAYLA_INTEGRATE_VEGAS_MISER:
        case BAYLA_INTEGRATE_VEGAS_PLUS:
            {
                retval->samples_used += n_samples[i];
            } break;


        default: Panic();
        }
    }

    
    switch(retval->strategy)
    {
    case BAYLA_INTEGRATE_VEGAS:
        {
            retval->total_samples = n_samples[0] * n_iter; 
        } break;

    case BAYLA_INTEGRATE_VEGAS_PLUS:
        {
            retval->total_samples = 0;
            for(size i = 0; i < n_iter; ++i)
            {
                retval->total_samples += n_samples[i];
            }
        } break;

    case BAYLA_INTEGRATE_VEGAS_MISER: {} break;

    default: Panic();
    }

    Assert(retval->total_samples >= retval->samples_used);
}

API BayLaEvidenceResult
bayla_vegas_integrate(BayLaModel const *model, BayLaIntegratorOptions *opts, MagStaticArena scratch_)
{
#pragma warning(push)
    /* Complains about conversion of samples_per_refinement and n_grid to f64. */
    Assert(opts->vegas.samples_per_refinement < 9007199254740992 && opts->vegas.num_grid_cells < 9007199254740992);
#pragma warning(disable: 4244)
    size samples_per_refinement = opts->vegas.samples_per_refinement;
    BayLaLogValue log_samples = bayla_log_value_map_into_log_domain(samples_per_refinement);
    size n_grid = opts->vegas.num_grid_cells;
    size max_total_samples = opts->vegas.max_total_samples;
    size min_total_samples = opts->vegas.min_total_samples;
    BayLaLogValue acceptable_abs_error = bayla_log_value_map_into_log_domain(opts->vegas.acceptable_abs_error);
    BayLaLogValue acceptable_prop_error = bayla_log_value_map_into_log_domain(opts->vegas.acceptable_prop_error);
    f64 alpha = opts->vegas.alpha;

    u64 seed = opts->vegas.seed;

    size ndim = model->n_parameters;
    f64 *min_vals = model->min_parameter_vals;
    f64 *max_vals = model->max_parameter_vals;

    ElkRandomState random_state_ = elk_random_state_create(seed);
    ElkRandomState *ran = &random_state_;

    BayLaEvidenceResult retval = {.strategy = BAYLA_INTEGRATE_VEGAS, .converged = true};

    /* Allocate some workspace memory. */
    MagStaticArena *scratch = &scratch_;
    BayLaLogValue *log_integral = eco_nmalloc(scratch, max_total_samples / samples_per_refinement + 1, BayLaLogValue);
    BayLaLogValue *log_sigma = eco_nmalloc(scratch, max_total_samples / samples_per_refinement + 1, BayLaLogValue);

    f64 *x = eco_nmalloc(scratch, ndim, f64);
    f64 *y = eco_nmalloc(scratch, ndim, f64);

    BayLaLogValue *log_fx_jys = eco_nmalloc(scratch, samples_per_refinement, BayLaLogValue);
    BayLaLogValue *log_fx_jy_sqs = eco_nmalloc(scratch, samples_per_refinement, BayLaLogValue);

    PanicIf(!(x && y && log_fx_jys && log_fx_jy_sqs));

    BayLaVegasMap *map = NULL;
    if(opts->vegas.map)
    {
        map = opts->vegas.map;
    }
    else
    {
       map = bayla_vegas_map_create(ndim, n_grid, alpha, max_vals, min_vals, scratch);
    }

    /* Calculate quantities used repeatedly in implementing the stratified sampling. */
    BayLaVegasStratificationParams strat_params = bayla_vegas_stratification_params_create(samples_per_refinement, ndim);

    /* Iterate over several batches to refine the grid. */
    size iter = 0;
    while(true)
    {
        for(size s = 0; s < samples_per_refinement; ++ s)
        {
            /* Generate random coordinates in the 'ndim' dimensions of y-space. 
             * 
             * Stratification for stratified sampling happens in here. 
             */
            bayla_vegas_stratification_sample_num_to_y(s, &strat_params, ran, y);

            /* Evaluate the function and accumulate (jy * fx) and its square */
            bayla_vegas_map_evaluate(model, map, y, &log_fx_jys[s], &log_fx_jy_sqs[s], scratch_);
        }

        BayLaLogValue log_sum = bayla_log_values_sum(samples_per_refinement, log_fx_jys);
        log_integral[iter] = bayla_log_value_divide(log_sum, log_samples);

        BayLaLogValue log_sum_sq = bayla_log_values_sum(samples_per_refinement, log_fx_jy_sqs);
        BayLaLogValue log_mean_sq = bayla_log_value_divide(log_sum_sq, log_samples);
        BayLaLogValue log_var = bayla_log_value_subtract(log_mean_sq, bayla_log_value_power(log_integral[iter], 2.0));
        log_var = bayla_log_value_divide(log_var, bayla_log_value_map_into_log_domain(samples_per_refinement - 1));
        log_sigma[iter] = bayla_log_value_power(log_var, 0.5);

        /* At this point, we've finished an iteration of the integral, we haven't refined the grid yet, so do that if
         * needed before going back through the loop!
         */
        iter++;

        /* Update the integral values and error metrics. */
        bayla_vegas_accumulate_results(&retval, iter, &samples_per_refinement, log_integral, log_sigma);

        /* Check for stopping conditions */
        if(retval.samples_used > min_total_samples)
        {
            if(!isnan(acceptable_abs_error.val) && retval.result.error.val < acceptable_abs_error.val)
            {
                break;
            }
            else if(!isnan(acceptable_prop_error.val) && bayla_log_value_divide(retval.result.error, retval.result.value).val < acceptable_prop_error.val)
            {
                break;
            }
        }

        if(retval.total_samples >= max_total_samples)
        {
            /* If we never achieved the desired error, then say it never converged. */
            if(!isnan(acceptable_abs_error.val) || !isnan(acceptable_prop_error.val)) { retval.converged = false; }
            break;
        }

        bayla_vegas_map_refine(map, scratch_);
    } /* iter */

    return retval;
#pragma warning(pop)
}

typedef struct
{
    size n_cube_samples;
    BayLaLogValue log_jf_sum;
    BayLaLogValue log_jf_sq_sum;
} BayLaHyperCube;

typedef struct
{
    size ndim;
    size n_strat;
    size n_hypercubes;
    size n_samples_per_iter;
    f64 hypercube_volume;
    BayLaHyperCube *cubes;
} BayLaVegasPlusStratificationParams;

API BayLaVegasPlusStratificationParams
bayla_vegas_plus_stratification_params_create(size const samples_per_iter, size const ndim, MagStaticArena *perm)
{
#pragma warning(push)
    /* Complains about possible loss of precision when converting 'size' to 'f64'. Any value that would cause a problem
     * is unrealistically large and would crash the computer due to a lack of memory anyway.
     */
#pragma warning(disable:4244)
    size n_strat = (size)floor(pow(samples_per_iter / 4, 1.0 / ndim)); if(n_strat < 1) { n_strat = 1; }
    size n_hypercubes = bayla_integer_exp(n_strat, ndim);
    f64 hypercube_volume = 1.0 / n_hypercubes;

    size n_samples_per_hypercube = samples_per_iter / n_hypercubes;

    BayLaHyperCube *cubes = eco_nmalloc(perm, n_hypercubes, BayLaHyperCube);
    for(size h = 0; h < n_hypercubes; ++h)
    {
        cubes[h] = (BayLaHyperCube){ .n_cube_samples = n_samples_per_hypercube }; /*Leave accumulators zero initialized.*/
    }

    return (BayLaVegasPlusStratificationParams)
        {
            .ndim = ndim,
            .n_strat = n_strat,
            .n_hypercubes = n_hypercubes,
            .n_samples_per_iter = samples_per_iter,
            .hypercube_volume = hypercube_volume,
            .cubes = cubes
        };
#pragma warning(pop)
}

API void
bayla_vegas_plus_stratification_cube_num_to_y(
        size const hn, 
        BayLaVegasPlusStratificationParams const *strat_params,
        ElkRandomState *ran,
        f64 *y)
{

#pragma warning(push)
    /* Warns about loss of data from convertion to u8 */
    Assert(strat_params->ndim < UINT8_MAX);
#pragma warning(disable: 4244)

    size index = hn;

    /* Convert list index to ndim dimensional coords. */
    for(size d = 0; d < strat_params->ndim; ++d)
    {
        size factor = bayla_integer_exp(strat_params->n_strat, strat_params->ndim - d - 1);
        size coords_of_hypercube = index / factor;
        index -= coords_of_hypercube * factor;
        f64 random = elk_random_state_uniform_f64(ran);

        y[d] = (coords_of_hypercube + random) / strat_params->n_strat;
        Assert(y[d] >= 0.0 && y[d] <= 1.0);
    }
#pragma warning(pop)
}

API void
bayla_vegas_plus_accumulate_hypercube_integrals(
        BayLaVegasPlusStratificationParams *strat_params, 
        BayLaLogValue *log_integral, 
        BayLaLogValue *log_sigma,
        size *iteration_samples)
{
#pragma warning(push)
    /* Complains about possible loss of precision when converting 'size' to 'f64'. Any value that would cause a problem
     * is unrealistically large and would crash the computer due to a lack of memory anyway.
     */
#pragma warning(disable:4244)
    BayLaLogValue log_sigma_sq_acc = bayla_log_value_map_into_log_domain(0.0);
    BayLaLogValue log_integral_acc = bayla_log_value_map_into_log_domain(0.0);
    *iteration_samples = 0;

    for(size h = 0; h < strat_params->n_hypercubes; ++h)
    {
        BayLaHyperCube *cube = &strat_params->cubes[h];

        BayLaLogValue log_cube_n_samples = bayla_log_value_map_into_log_domain(cube->n_cube_samples);
        BayLaLogValue log_cube_integral = bayla_log_value_divide(cube->log_jf_sum, log_cube_n_samples);

        BayLaLogValue log_mean_sq = bayla_log_value_divide(cube->log_jf_sq_sum, log_cube_n_samples);
        BayLaLogValue log_sqrd_mean = bayla_log_value_power(log_cube_integral, 2.0);
        BayLaLogValue log_cube_var;
        if(log_mean_sq.val >= log_sqrd_mean.val)
        {
            log_cube_var = bayla_log_value_subtract(log_mean_sq, log_sqrd_mean);
            log_cube_var = bayla_log_value_divide(log_cube_var, log_cube_n_samples);
        }
        else
        {
            log_cube_var = bayla_log_value_map_into_log_domain(0.0);
        }
        Assert(!isnan(log_cube_var.val));

        log_integral_acc = bayla_log_value_add(log_integral_acc, log_cube_integral);
        log_sigma_sq_acc = bayla_log_value_add(log_sigma_sq_acc, bayla_log_value_divide(log_cube_var, log_cube_n_samples));
        *iteration_samples += cube->n_cube_samples;
    }

    *log_integral = bayla_log_value_divide(log_integral_acc, bayla_log_value_map_into_log_domain(strat_params->n_hypercubes));
    *log_sigma = bayla_log_value_multiply(
            bayla_log_value_power(log_sigma_sq_acc, 0.5),
            bayla_log_value_map_into_log_domain(strat_params->hypercube_volume));

#pragma warning(pop)
}

API void
bayla_vegas_plus_refine_sampling_stratification(
        BayLaVegasPlusStratificationParams *strat_params, 
        f64 beta, 
        MagStaticArena scratch)
{
#pragma warning(push)
    /* Complains about possible loss of precision when converting 'size' to 'f64'. Any value that would cause a problem
     * is unrealistically large and would crash the computer due to a lack of memory anyway.
     */
#pragma warning(disable:4244)
    size n_hypercubes = strat_params->n_hypercubes;

    BayLaDualDomainValue *ds = eco_nmalloc(&scratch, n_hypercubes, BayLaDualDomainValue);

    for(size h = 0; h < n_hypercubes; ++h)
    {
        BayLaHyperCube *cube = &strat_params->cubes[h];


        BayLaLogValue log_cube_n_samples = bayla_log_value_map_into_log_domain(cube->n_cube_samples);
        BayLaLogValue log_cube_integral = bayla_log_value_divide(cube->log_jf_sum, log_cube_n_samples);

        BayLaLogValue log_mean_sq = bayla_log_value_divide(cube->log_jf_sq_sum, log_cube_n_samples);
        BayLaLogValue log_sqrd_mean = bayla_log_value_power(log_cube_integral, 2.0);
        BayLaLogValue log_cube_var;
        if(log_mean_sq.val >= log_sqrd_mean.val)
        {
            log_cube_var = bayla_log_value_subtract(log_mean_sq, log_sqrd_mean);
            log_cube_var = bayla_log_value_divide(log_cube_var, bayla_log_value_map_into_log_domain(cube->n_cube_samples - 1));
        }
        else
        {
            log_cube_var = bayla_log_value_map_into_log_domain(0.0);
        }
        Assert(!isnan(log_cube_var.val));

        BayLaLogValue log_cube_sigma = bayla_log_value_power(
                bayla_log_value_multiply(
                    log_cube_var, bayla_log_value_map_into_log_domain(strat_params->hypercube_volume)),
                0.5);

        ds[h].log_domain = bayla_log_value_power(log_cube_sigma, beta);
        if(ds[h].log_domain.val < log(TINY) || isnan(ds[h].log_domain.val)) { ds[h].log_domain.val = log(TINY); }

        Assert(!isnan(ds[h].log_domain.val));
    }

    BayLaLogValue log_sum_ds = bayla_log_values_sum(n_hypercubes, (BayLaLogValue *)ds);
    Assert(!isnan(log_sum_ds.val));

    for(size h = 0; h < n_hypercubes; ++h)
    {
        ds[h].log_domain = bayla_log_value_divide(ds[h].log_domain, log_sum_ds);
        ds[h].linear_domain = bayla_log_value_map_out_of_log_domain(ds[h].log_domain);
        Assert(ds[h].linear_domain >= 0.0 && ds[h].linear_domain <= 1.0);
    }

    for(size h = 0; h < n_hypercubes; ++h)
    {
        size n_cube_samples = strat_params->n_samples_per_iter * ds[h].linear_domain;
        strat_params->cubes[h].n_cube_samples = n_cube_samples > 3 ? n_cube_samples : 3;
    }
#pragma warning(pop)
}


API BayLaEvidenceResult
bayla_vegas_plus_integrate(BayLaModel const *model, BayLaIntegratorOptions *opts, MagStaticArena scratch_)
{
#pragma warning(push)
    /* Complains about possible loss of precision when converting 'size' to 'f64'. Any value that would cause a problem
     * is unrealistically large and would crash the computer due to a lack of memory anyway.
     */
#pragma warning(disable:4244)
    size samples_per_refinement = opts->vegas_plus.samples_per_refinement;
    size n_grid = opts->vegas_plus.num_grid_cells;
    size max_total_samples = opts->vegas_plus.max_total_samples;
    size min_total_samples = opts->vegas_plus.min_total_samples;
    BayLaLogValue acceptable_abs_error = bayla_log_value_map_into_log_domain(opts->vegas_plus.acceptable_abs_error);
    BayLaLogValue acceptable_prop_error = bayla_log_value_map_into_log_domain(opts->vegas_plus.acceptable_prop_error);
    f64 alpha = opts->vegas_plus.alpha;
    f64 beta = opts->vegas_plus.beta;

    u64 seed = opts->vegas_plus.seed;

    size ndim = model->n_parameters;
    f64 *min_vals = model->min_parameter_vals;
    f64 *max_vals = model->max_parameter_vals;

    ElkRandomState random_state_ = elk_random_state_create(seed);
    ElkRandomState *ran = &random_state_;

    BayLaEvidenceResult retval = {.strategy = BAYLA_INTEGRATE_VEGAS_PLUS, .converged = true};

    size max_iterations = max_total_samples / samples_per_refinement * 2; /* Give some extra space here. */

    /* Allocate some workspace memory. */
    MagStaticArena *scratch = &scratch_;
    BayLaLogValue *log_integral = eco_nmalloc(scratch, max_iterations, BayLaLogValue);
    BayLaLogValue *log_sigma = eco_nmalloc(scratch, max_iterations, BayLaLogValue);
    size *samples_per_iter = eco_nmalloc(scratch, max_iterations, size);
    PakArrayLedger iter_ledger = pak_array_ledger_create(max_iterations);

    f64 *x = eco_nmalloc(scratch, ndim, f64);
    f64 *y = eco_nmalloc(scratch, ndim, f64);

    PanicIf(!(x && y));

    BayLaVegasMap *map = NULL;
    if(opts->vegas_plus.map)
    {
        map = opts->vegas_plus.map;
    }
    else
    {
       map = bayla_vegas_map_create(ndim, n_grid, alpha, max_vals, min_vals, scratch);
    }

    /* Calculate quantities used repeatedly in implementing the stratified sampling. */
    BayLaVegasPlusStratificationParams strat_params =
        bayla_vegas_plus_stratification_params_create(samples_per_refinement, ndim, scratch);

    /* Iterate over several batches to refine the grid. */
    size iter = 0;
    while(true)
    {
        for(size h = 0; h < strat_params.n_hypercubes; ++h)
        {
            strat_params.cubes[h].log_jf_sum = bayla_log_zero;
            strat_params.cubes[h].log_jf_sq_sum = bayla_log_zero;

            for(size s = 0; s < strat_params.cubes[h].n_cube_samples; ++s)
            {
                /* Generate random coordinates in the 'ndim' dimensions of y-space. 
                 * 
                 * Stratification for stratified sampling happens in here. 
                 */
                bayla_vegas_plus_stratification_cube_num_to_y(h, &strat_params, ran, y);

                BayLaLogValue log_jf;
                BayLaLogValue log_jf_sq;

                /* Evaluate the function and accumulate (jy * fx) and its square */
                bayla_vegas_map_evaluate(model, map, y, &log_jf, &log_jf_sq, scratch_);

                strat_params.cubes[h].log_jf_sum = bayla_log_value_add(strat_params.cubes[h].log_jf_sum, log_jf);
                strat_params.cubes[h].log_jf_sq_sum = bayla_log_value_add(strat_params.cubes[h].log_jf_sq_sum, log_jf_sq);
            }
        }

        size iter_idx = pak_array_ledger_push_back_index(&iter_ledger);
        bayla_vegas_plus_accumulate_hypercube_integrals(
                &strat_params, 
                &log_integral[iter_idx], 
                &log_sigma[iter_idx], 
                &samples_per_iter[iter_idx]);

        /* At this point, we've finished an iteration of the integral, we haven't refined the grid yet, so do that if
         * needed before going back through the loop!
         */
        iter++;

        bayla_vegas_accumulate_results(&retval, iter, samples_per_iter, log_integral, log_sigma);

        /* Check for stopping conditions */
        if(retval.samples_used > min_total_samples)
        {
            if(!isnan(acceptable_abs_error.val) && retval.result.error.val < acceptable_abs_error.val)
            {
                break;
            }
            else if(!isnan(acceptable_prop_error.val) && bayla_log_value_divide(retval.result.error, retval.result.value).val < acceptable_prop_error.val)
            {
                break;
            }
        }

        if(retval.total_samples >= max_total_samples || pak_array_ledger_full(&iter_ledger))
        {
            /* If we never achieved the desired error, then say it never converged. */
            if(!isnan(acceptable_abs_error.val) || !isnan(acceptable_prop_error.val)) { retval.converged = false; }
            break;
        }

        bayla_vegas_map_refine(map, scratch_);
        bayla_vegas_plus_refine_sampling_stratification(&strat_params, beta, scratch_);

    } /* iter */

    return retval;
#pragma warning(pop)
}

typedef struct
{
    BayLaModel const *model;
    size ndim;
    BayLaVegasMap *map;
} BayLaMiserInnerFunctionArgs;

API BayLaLogValue
bayla_miser_inner_function_evaluate(BayLaMiserInnerFunctionArgs args, f64 *ys, MagStaticArena scratch)
{
    if(args.map)
    {
        BayLaLogValue point_val = {0};
        BayLaLogValue p2 = {0};

        bayla_vegas_map_evaluate(args.model, args.map, ys, &point_val, &p2, scratch);

        return point_val;
    }
    else
    {
        return bayla_model_evaluate(args.model, ys);
    }
}

API BayLaValue
bayla_miser_integrate_inner(
        BayLaMiserInnerFunctionArgs args,
        f64 const *min_xs,
        f64 const *max_xs,
        size const npts,
        size const min_points,
        size const min_to_subdivide,
        size *n_samples,
        size *samples_used,
        f64 explore_factor,
        ElkRandomState *ran,
        MagStaticArena scratch)
{
#pragma warning(push)
    /* Complains about possible loss of precision when converting 'size' to 'f64'. Any value that would cause a problem
     * is unrealistically large and would crash the computer due to a lack of memory anyway.
     */
#pragma warning(disable:4244)

    BayLaLogValue const log_tiny = bayla_log_value_map_into_log_domain(TINY);
    size ndim = args.ndim;

    f64 *lengths = eco_nmalloc(&scratch, ndim, f64);
    for(size d = 0; d < ndim; ++d)
    {
        lengths[d] = (max_xs[d] - min_xs[d]);
    }

    if(npts < min_to_subdivide)
    {
        /* Too few points left to subdivide, use a simple Monte Carlo method. */

        f64 *xs = eco_nmalloc(&scratch, ndim, f64);

        BayLaLogValue log_sum = bayla_log_value_map_into_log_domain(0);
        BayLaLogValue log_sumsq = bayla_log_value_map_into_log_domain(0);

        for(size i = 0; i < npts; ++i)
        {
            /* Generate a random point inside the domain. */
            for(size d = 0; d < ndim; ++d)
            {
                xs[d] = elk_random_state_uniform_f64(ran) * lengths[d] + min_xs[d];
            }

            /* Evaluate the function value. */
            BayLaLogValue point_val = bayla_miser_inner_function_evaluate(args, xs, scratch);
            Assert(!isnan(point_val.val));

            /*  update the sum & sumsq accumulators */
            log_sum = bayla_log_value_add(log_sum, point_val);
            log_sumsq = bayla_log_value_add(log_sumsq, bayla_log_value_power(point_val, 2.0));
        }

        *n_samples += npts;
        *samples_used += npts;

        BayLaLogValue log_npts = bayla_log_value_map_into_log_domain(npts);
        BayLaLogValue log_nptsm1 = bayla_log_value_map_into_log_domain(npts - 1);

        BayLaLogValue log_mean = bayla_log_value_divide(log_sum, log_npts); 
        BayLaLogValue log_sqr_mean = bayla_log_value_divide(log_sumsq, log_npts); 

        BayLaLogValue log_mean2 = bayla_log_value_power(log_mean, 2.0);
        BayLaLogValue log_var = bayla_log_value_divide(bayla_log_value_subtract(log_sqr_mean, log_mean2), log_nptsm1);
        if(isnan(log_var.val)) log_var = log_tiny;

        return (BayLaValue){ .value = log_mean, .error = log_var };
    }
    else
    {
        /* Do some preliminary uniform sampling to decide how to bifurcate this domain. */
        f64 *xs = eco_nmalloc(&scratch, ndim, f64);
        f64 *x_mids = eco_nmalloc(&scratch, ndim, f64);
        BayLaLogValue *fmaxl = eco_nmalloc(&scratch, ndim, BayLaLogValue);
        BayLaLogValue *fmaxr = eco_nmalloc(&scratch, ndim, BayLaLogValue);
        BayLaLogValue *fminl = eco_nmalloc(&scratch, ndim, BayLaLogValue);
        BayLaLogValue *fminr = eco_nmalloc(&scratch, ndim, BayLaLogValue);
        f64 *max_xs_l = eco_nmalloc(&scratch, ndim, f64);
        f64 *min_xs_r = eco_nmalloc(&scratch, ndim, f64);

        Assert(xs && x_mids && fmaxl && fmaxr && fminl && fminr && max_xs_l && min_xs_r);

        /* Initialize the mid-points for each dimension. */
        for(size d = 0; d < ndim; ++d)
        {
            x_mids[d] = (min_xs[d] + max_xs[d]) / 2.0;
            fminl[d].val = fminr[d].val = INFINITY;
            fmaxl[d].val = fmaxr[d].val = -INFINITY;
        }

        /* Loop over the points in the sample. */
        size n_pre = (size)(npts * explore_factor);
        n_pre = n_pre > min_points ? n_pre : min_points;

        for(size n = 0; n < n_pre; ++n)
        {
            /* Generate a random point inside the domain. */
            for(size d = 0; d < ndim; ++d)
            {
                xs[d] = elk_random_state_uniform_f64(ran) * lengths[d] + min_xs[d];
            }

            /* Evaluate the function value. */
            BayLaLogValue point_val = bayla_miser_inner_function_evaluate(args, xs, scratch);
            Assert(!isnan(point_val.val));

            for(size d = 0; d < ndim; ++d)
            {
                if(xs[d] < x_mids[d])
                {
                    fminl[d] = fminl[d].val < point_val.val ? fminl[d] : point_val;
                    fmaxl[d] = fmaxl[d].val > point_val.val ? fmaxl[d] : point_val;
                }
                else
                {
                    fminr[d] = fminr[d].val < point_val.val ? fminr[d] : point_val;
                    fmaxr[d] = fmaxr[d].val > point_val.val ? fmaxr[d] : point_val;
                }
            }
        }

        /* Choose which dimension to bisect. */
        size db = -1;
        BayLaLogValue sumb = {.val = INFINITY};
        BayLaLogValue siglb = bayla_log_one;
        BayLaLogValue sigrb = bayla_log_one;
        for(size d = 0; d < ndim; ++d)
        {
            if(fmaxl[d].val > fminl[d].val && fmaxr[d].val > fminr[d].val)
            {
                BayLaLogValue lval = bayla_log_value_power(bayla_log_value_subtract(fmaxl[d], fminl[d]), 2.0/3.0);
                BayLaLogValue sigl = log_tiny.val > lval.val ? log_tiny : lval;

                BayLaLogValue rval = bayla_log_value_power(bayla_log_value_subtract(fmaxr[d], fminr[d]), 2.0/3.0);
                BayLaLogValue sigr = log_tiny.val > rval.val ? log_tiny : rval;

                BayLaLogValue sum = bayla_log_value_add(sigl, sigr);

                if(sum.val <= sumb.val)
                {
                    sumb = sum;
                    db = d;
                    siglb = sigl;
                    sigrb = sigr;
                }
            }
        }

        if(db == -1) /* min_points may be too small, assign one randomly. */
        {
            db = (size)(elk_random_state_uniform_u64(ran) % (u64)ndim);
        }

        f64 rgl = min_xs[db];
        f64 rgm = x_mids[db];
        f64 rgr = max_xs[db];
        f64 fracl = fabs((rgm - rgl) / (rgr - rgl));

        f64 siglb_linear = bayla_log_value_map_out_of_log_domain(siglb);
        f64 sigrb_linear = bayla_log_value_map_out_of_log_domain(sigrb);
        size nptl = (size)(
                min_points + (npts - n_pre - 2 * min_points) * fracl * siglb_linear /
                (fracl * siglb_linear + (1.0 - fracl) * sigrb_linear));
        size nptr = npts - n_pre - nptl;
        
        Assert(nptl >= min_points && nptr >= min_points && fracl > 0.0);

        /* Set up other domains for recursion. */
        memcpy(max_xs_l, max_xs, sizeof(f64) * ndim);
        max_xs_l[db] = x_mids[db];
        BayLaValue lave = bayla_miser_integrate_inner(
                args,
                min_xs,
                max_xs_l,
                nptl,
                min_points,
                min_to_subdivide,
                n_samples,
                samples_used,
                explore_factor,
                ran,
                scratch);

        memcpy(min_xs_r, min_xs, sizeof(f64) * ndim);
        min_xs_r[db] = x_mids[db];
        BayLaValue rave = bayla_miser_integrate_inner(
                args,
                min_xs_r,
                max_xs,
                nptr,
                min_points,
                min_to_subdivide,
                n_samples,
                samples_used,
                explore_factor,
                ran,
                scratch);

        BayLaLogValue lfracl = bayla_log_value_map_into_log_domain(fracl);
        BayLaLogValue lfracr = bayla_log_value_map_into_log_domain(1.0 - fracl);

        BayLaLogValue ave = bayla_log_value_add(
                bayla_log_value_multiply(lfracl, lave.value),
                bayla_log_value_multiply(lfracr, rave.value));

        BayLaLogValue var = bayla_log_value_add(
                bayla_log_value_multiply(bayla_log_value_power(lfracl, 2.0), lave.error),
                bayla_log_value_multiply(bayla_log_value_power(lfracr, 2.0), rave.error));

        *n_samples += n_pre;

        return (BayLaValue){ .value = ave, .error = var};
    }

    return (BayLaValue){0};
#pragma warning(pop)
}

API BayLaEvidenceResult
bayla_miser_integrate(BayLaModel const *model, BayLaIntegratorOptions *opts, MagStaticArena scratch)
{
    BayLaEvidenceResult retval = {.strategy = BAYLA_INTEGRATE_MISER};

    size min_points = opts->miser.min_points;
    size min_to_subdivide = opts->miser.min_to_subdivide;
    size total_samples = opts->miser.total_samples;
    f64 explore_factor = opts->miser.explore_factor;
    BayLaLogValue acceptable_abs_error = bayla_log_value_map_into_log_domain(opts->miser.acceptable_abs_error);
    BayLaLogValue acceptable_prop_error = bayla_log_value_map_into_log_domain(opts->miser.acceptable_prop_error);

    u64 seed = opts->miser.seed;

    size ndim = model->n_parameters;
    f64 *min_vals = model->min_parameter_vals;
    f64 *max_vals = model->max_parameter_vals;

    f64 volume = 1.0;
    for(size d = 0; d < ndim; ++d)
    {
        volume *= (max_vals[d] - min_vals[d]);
    }
    BayLaLogValue log_volume = bayla_log_value_map_into_log_domain(volume);

    BayLaMiserInnerFunctionArgs args = { .model = model, .ndim = ndim, .map = NULL };

    ElkRandomState random_state_ = elk_random_state_create(seed);
    ElkRandomState *ran = &random_state_;

    retval.result = bayla_miser_integrate_inner(
            args,
            min_vals,
            max_vals,
            total_samples,
            min_points,
            min_to_subdivide,
            &retval.total_samples,
            &retval.samples_used,
            explore_factor,
            ran,
            scratch);

    retval.result.error = bayla_log_value_multiply(log_volume, bayla_log_value_power(retval.result.error, 0.5));
    retval.result.value = bayla_log_value_multiply(log_volume, retval.result.value);

    b32 absolutely_converged = !isnan(acceptable_abs_error.val) && retval.result.error.val < acceptable_abs_error.val;
    b32 proportionaly_converged = !isnan(acceptable_prop_error.val) &&
        bayla_log_value_divide(retval.result.error, retval.result.value).val < acceptable_prop_error.val;

    if(absolutely_converged || proportionaly_converged)
    {
        retval.converged = true;
    }
    else
    {
        retval.converged = false;
    }

    return retval;
}

API BayLaEvidenceResult
bayla_vegas_miser_integrate(BayLaModel const *model, BayLaIntegratorOptions *opts, MagStaticArena scratch_)
{
    BayLaEvidenceResult retval = {.strategy = BAYLA_INTEGRATE_VEGAS_MISER, .converged = true };

    MagStaticArena *scratch = &scratch_;

    size samples_per_refinement = opts->vegas_miser.samples_per_refinement;
    size n_grid = opts->vegas_miser.num_grid_cells;
    size min_total_samples = opts->vegas_miser.min_total_samples;
    size max_total_samples = opts->vegas_miser.max_total_samples;
    size min_points = opts->vegas_miser.min_points;
    size min_to_subdivide = opts->vegas_miser.min_to_subdivide;
    f64 alpha = opts->vegas_miser.alpha;
    f64 explore_factor = opts->vegas_miser.explore_factor;
    BayLaLogValue acceptable_abs_error = bayla_log_value_map_into_log_domain(opts->vegas_miser.acceptable_abs_error);
    BayLaLogValue acceptable_prop_error = bayla_log_value_map_into_log_domain(opts->vegas_miser.acceptable_prop_error);

    u64 seed = opts->vegas_miser.seed;

    size ndim = model->n_parameters;
    f64 *min_vals = model->min_parameter_vals;
    f64 *max_vals = model->max_parameter_vals;

    ElkRandomState random_state_ = elk_random_state_create(seed);
    ElkRandomState *ran = &random_state_;

    BayLaVegasMap *map = NULL;
    if(opts->vegas_miser.map)
    {
        map = opts->vegas_miser.map;
    }
    else
    {
       map = bayla_vegas_map_create(ndim, n_grid, alpha, max_vals, min_vals, scratch);
    }

    size max_iterations = max_total_samples / samples_per_refinement * 2; /* Give some extra space here. */

    f64 *min_ys = eco_nmalloc(scratch, ndim, f64);
    f64 *max_ys = eco_nmalloc(scratch, ndim, f64);
    for(size d = 0; d < ndim; ++d)
    {
        min_ys[d] = 0.0;
        max_ys[d] = 1.0;
    }

    BayLaMiserInnerFunctionArgs args = { .model = model, .ndim = ndim, .map = map };

    /* Allocate some workspace memory. */
    BayLaLogValue *log_integral = eco_nmalloc(scratch, max_iterations, BayLaLogValue);
    BayLaLogValue *log_integral_err = eco_nmalloc(scratch, max_iterations, BayLaLogValue);
    size *samples_per_iter = eco_nmalloc(scratch, max_iterations, size);

    size iter = 0;
    while(true)
    {
        size total_samples = 0;

        BayLaValue integral = bayla_miser_integrate_inner(
                args,
                min_ys,
                max_ys,
                samples_per_refinement,
                min_points, 
                min_to_subdivide,
                &total_samples,
                &samples_per_iter[iter],
                explore_factor,
                ran,
                scratch_);

        log_integral[iter] = integral.value;
        log_integral_err[iter] = bayla_log_value_power(integral.error, 0.5);
        retval.total_samples += total_samples;

        iter++;

        bayla_vegas_accumulate_results(&retval, iter, samples_per_iter, log_integral, log_integral_err);

        /* Check stopping conditions. */
        if(retval.samples_used > min_total_samples)
        {
            if(!isnan(acceptable_abs_error.val) && retval.result.error.val < acceptable_abs_error.val)
            {
                break;
            }
            else if(!isnan(acceptable_prop_error.val) && bayla_log_value_divide(retval.result.error, retval.result.value).val < acceptable_prop_error.val)
            {
                break;
            }
        }

        if(retval.total_samples >= max_total_samples || iter >= max_iterations)
        {
            /* If we never achieved the desired error, then say it never converged. */
            if(!isnan(acceptable_abs_error.val) || !isnan(acceptable_prop_error.val)) { retval.converged = false; }
            break;
        }

        bayla_vegas_map_refine(map, scratch_);

    }

    return retval;
}

#undef TINY

API BayLaEvidenceResult 
bayla_calculate_evidence(BayLaModel const *model, BayLaIntegratorOptions *opts, MagStaticArena scratch)
{
    switch(opts->strategy)
    {
        case BAYLA_INTEGRATE_VEGAS_PLUS:  return bayla_vegas_plus_integrate(model, opts, scratch);
        case BAYLA_INTEGRATE_VEGAS_MISER: return bayla_vegas_miser_integrate(model, opts, scratch);
        case BAYLA_INTEGRATE_VEGAS:       return bayla_vegas_integrate(model, opts, scratch);
        case BAYLA_INTEGRATE_MONTE_CARLO: return bayla_simple_monte_carlo_integrate(model, opts, scratch);
        case BAYLA_INTEGRATE_MISER:       return bayla_miser_integrate(model, opts, scratch);
        default:
            {
                return (BayLaEvidenceResult)
                    { 
                        .strategy = BAYLA_INTEGRATE_INVALID,
                        .result = (BayLaValue){.value = (BayLaLogValue){.val = NAN}, .error = (BayLaLogValue){.val = NAN } },
                        .converged = false
                    };
            }
    }
}

#pragma warning(pop)

#endif
