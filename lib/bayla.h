#ifndef _BAYLA_H_
#define _BAYLA_H_

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

/* The location in parameter space of the max a posteriori.
 *
 * The "data" is fixed and stored somewhere in the user_data. When you evaluate this function, it fills the buffer parms_out
 * with the model parameters that maximize the posterior for this model and the data provided.
 */
typedef void (*BayLaMaxAPosteriori)(void *user_data, size num_parms, f64 *parms_out);

/* TODO: a function for calculating the covariance matrix.
 * The "data" is fixed and stored somewhere in the user_data. When you evaluate this function, it fills the buffer cov_matrix
 * (which is sized num_parms x num_parms) with entries in the covariance matrix. This, along with the max_a_posteriori_parms
 * are used to create a Gaussian approximation to the posterior distribution.
 */
typedef void (*BayLaCovarianceMatrix)(void *user_data, size num_parms, f64 const *max_a_posteriori_parms, f64 *cov_matrix);

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
    BayLaMaxAPosteriori posterior_max;         /* See typedef above. argmax(Θ) P(Θ | Data, Model, I)                      */
    BayLaCovarianceMatrix calc_cov_matrix;     /* See typedef above. Used to create a Gaussian approximation.             */

    f64 *min_parameter_vals;        /* Array n_parameters long of the minimum value parameters can practically have.      */
    f64 *max_parameter_vals;        /* Array n_parameters long of the miximum value parameters can practically have.      */

    void *user_data;                /* Any other state or data you may need for evaluating the prior, likelihood, etc.    */
} BayLaModel;

/*---------------------------------------------------------------------------------------------------------------------------
 *                                                    Implementations
 *-------------------------------------------------------------------------------------------------------------------------*/

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



#pragma warning(pop)

#endif
