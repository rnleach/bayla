/*------------------------------------------------- Data Gen Function ----------------------------------------------------*/

static inline f64
natural_log_approximation_by_3rd_order_polynomial(f64 x, f64 x0)
{
    /* This is a 3rd order polynomial approximation to ln around the point x0 */
    f64 d = x - x0;
    return log(x0) + d / x0 - (d * d) / (2 * x0 * x0) + (d * d * d) / ( 3.0 * x0 * x0 * x0);
}

/* Some data used by all models, this will be initialized in the test. */
#define NUM_DATA_POINTS 50
typedef struct
{
    /* For keeping track of how many times a function is called. */
    size ncalls;

    /* The original data used to compute the stats. */
    f64 *xs;
    f64 *ys;

    /* Useful statistics for calculating the probability of the data given the parameters. */
    size N;

    f64 y_sq_bar;
    f64 x_sq_bar;
    f64 x4_bar;
    f64 x6_bar;
    f64 x8_bar;
    f64 x10_bar;

    f64 y_bar;
    f64 x_bar;

    f64 xy_bar;
    f64 x3_bar;
    f64 x5_bar;
    f64 x7_bar;
    f64 x9_bar;
    f64 y_x2_bar;
    f64 y_x3_bar;
    f64 y_x4_bar;
    f64 y_x5_bar;

    f64 y_logx_bar;
    f64 logx_sq_bar;
} UserData;

f64 global_xs[NUM_DATA_POINTS];
f64 global_ys[NUM_DATA_POINTS];

/*-------------------------------------------------  Logarithmic Model  --------------------------------------------------*/
/* parameter 0 a standard deviation. */
f64 log_min_parms[1] = {0.0};
f64 log_max_parms[1] = {20.0};

static f64 log_parameters_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    /* exponential with mean 1 */
    return -parms[0] - bayla_log1mexp(log_max_parms[0]);
}

static f64 log_parameters_log_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 sigma = parms[0];
    f64 y_sq_bar = data->y_sq_bar;
    f64 y_logx_bar = data->y_logx_bar;
    f64 logx_sq_bar = data->logx_sq_bar;

    data->ncalls++;

    f64 log_prob = -N * (log(sigma) + 0.5 * log(2.0 * ELK_PI));
    log_prob -= N * (y_sq_bar - 2.0 * y_logx_bar + logx_sq_bar) / (2.0 * sigma * sigma);

    return log_prob;
}

UserData log_user_data = {0};

BayLaModel log_model = 
    {
        .model_prior_probability = 1.0 / 7.0,
        .n_parameters = 1,
        .prior = log_parameters_log_prior,
        .likelihood = log_parameters_log_likelihood,
        .min_parameter_vals = log_min_parms,
        .max_parameter_vals = log_max_parms,
        .user_data = &log_user_data
    };

/*--------------------------------------------------- Constant Model -----------------------------------------------------*/
/* parameter 0 is a constant value, parameter 1 is a standard deviation. */
f64 constant_min_parms[2] = {-2.0, 0.0};
f64 constant_max_parms[2] = { 2.0, 20.0};

static f64 constant_parameters_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    /* exponential with mean 1 for std-dev and uniform for all other parameters. */
    return log((1.0 / 4.0)) - parms[1] - bayla_log1mexp(constant_max_parms[1]);
}

static f64 constant_parameters_log_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 c = parms[0];
    f64 sigma = parms[1];
    f64 y_bar = data->y_bar;
    f64 y_sq_bar = data->y_sq_bar;

    data->ncalls++;

    f64 log_prob = -N * (log(sigma) + 0.5 * log(2.0 * ELK_PI));
    log_prob -= N * (y_sq_bar - 2.0 * y_bar * c + c * c) / (2.0 * sigma * sigma);

    return log_prob;
}

UserData constant_user_data = {0};

BayLaModel constant_model = 
    {
        .model_prior_probability = 1.0 / 7.0,
        .n_parameters = 2,
        .prior = constant_parameters_log_prior,
        .likelihood = constant_parameters_log_likelihood,
        .min_parameter_vals = constant_min_parms,
        .max_parameter_vals = constant_max_parms,
        .user_data = &constant_user_data
    };

/*---------------------------------------------------  Linear Model  -----------------------------------------------------*/
/* parameter 0 is a constant value, parameter 1 is a linear coefficent, parameter 2 is a standard deviation. */
f64 linear_min_parms[3] = {-2.0, 0.0, 0.0};
f64 linear_max_parms[3] = { 0.0, 4.0, 20.0};

static f64 linear_parameters_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    /* exponential with mean 1 for std-dev and uniform for all other parameters. */
    return log((1.0 / 2.0) * (1.0 / 4.0)) - parms[2] - bayla_log1mexp(linear_max_parms[2]);
}

static f64 linear_parameters_log_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 c = parms[0];
    f64 b = parms[1];
    f64 sigma = parms[2];
    f64 y_sq_bar = data->y_sq_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 y_bar = data->y_bar;
    f64 x_bar = data->x_bar;
    f64 xy_bar = data->xy_bar;

    data->ncalls++;

    f64 log_prob = -N * (log(sigma) + 0.5 * log(2.0 * ELK_PI));
    log_prob -= N * (y_sq_bar + 2.0 * (c * b * x_bar - b * xy_bar - c * y_bar) + b * b * x_sq_bar + c * c) / (2.0 * sigma * sigma);

    return log_prob;
}

UserData linear_user_data = {0};

BayLaModel linear_model = 
    {
        .model_prior_probability = 1.0 / 7.0,
        .n_parameters = 3,
        .prior = linear_parameters_log_prior,
        .likelihood = linear_parameters_log_likelihood,
        .min_parameter_vals = linear_min_parms,
        .max_parameter_vals = linear_max_parms,
        .user_data = &linear_user_data
    };

/*-----------------------------------------------  2nd Order Polynomial  -------------------------------------------------*/
/* param 0 is a constant, param 1 is a linear coeff, param 2 is a quadratic coeff,  param 3 is a standard deviation. */
f64 second_order_min_parms[4] = {-2.0, -4.0, -2.0, 0.0};
f64 second_order_max_parms[4] = {+0.0, +4.0, +2.0, 20.0};

static f64 second_order_parameters_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    /* exponential with mean 1 for std-dev and uniform for all other parameters. */
    return log((1.0 / 2.0) * (1.0 / 8.0) * (1.0 / 4.0)) - parms[3] - bayla_log1mexp(second_order_max_parms[3]);
}

static f64 second_order_parameters_log_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 c = parms[0];
    f64 b = parms[1];
    f64 a = parms[2];
    f64 sigma = parms[3];
    f64 y_sq_bar = data->y_sq_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x4_bar = data->x4_bar;
    f64 y_bar = data->y_bar;
    f64 x_bar = data->x_bar;
    f64 xy_bar = data->xy_bar;
    f64 x3_bar = data->x3_bar;
    f64 y_x2_bar = data->y_x2_bar;

    data->ncalls++;

    f64 log_prob = -N * (log(sigma) + 0.5 * log(2.0 * ELK_PI));

    f64 poly_term = y_sq_bar + x_sq_bar * (2.0 * a * c + b * b)
                  + 2.0 * (a * b * x3_bar + c * b * x_bar - a * y_x2_bar - b * xy_bar - c * y_bar)
                  + a * a * x4_bar + c * c;

    log_prob -= N * poly_term / (2.0 * sigma * sigma);

    return log_prob;
}

UserData second_order_user_data = {0};

BayLaModel second_order_model = 
    {
        .model_prior_probability = 1.0 / 7.0,
        .n_parameters = 4,
        .prior = second_order_parameters_log_prior,
        .likelihood = second_order_parameters_log_likelihood,
        .min_parameter_vals = second_order_min_parms,
        .max_parameter_vals = second_order_max_parms,
        .user_data = &second_order_user_data
    };

/*-----------------------------------------------  3rd Order Polynomial  -------------------------------------------------*/
/* param 0 is a constant, param 1 is a linear coeff, param 2 is a quadratic coeff, ....  param 4 is a standard deviation. */
f64 third_order_min_parms[5] = {-2.0, -4.0, -2.0, -2.0, 0.0};
f64 third_order_max_parms[5] = {+0.0, +4.0, +2.0, +2.0, 20.0};

static f64 third_order_parameters_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    /* exponential with mean 1 for std-dev and uniform for all other parameters. */
    return log((1.0 / 2.0) * (1.0 / 8.0) * (1.0 / 4.0) * (1.0 / 4.0)) - parms[4] - bayla_log1mexp(third_order_max_parms[4]);
}

static f64 third_order_parameters_log_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 c = parms[0];
    f64 b = parms[1];
    f64 a = parms[2];
    f64 d = parms[3];

    f64 sigma = parms[4];

    f64 y_sq_bar = data->y_sq_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x4_bar = data->x4_bar;
    f64 x6_bar = data->x6_bar;

    f64 y_bar = data->y_bar;
    f64 x_bar = data->x_bar;

    f64 xy_bar = data->xy_bar;
    f64 x3_bar = data->x3_bar;
    f64 x5_bar = data->x5_bar;
    f64 y_x2_bar = data->y_x2_bar;
    f64 y_x3_bar = data->y_x3_bar;

    data->ncalls++;

    f64 log_prob = -N * (log(sigma) + 0.5 * log(2.0 * ELK_PI));

    f64 poly_term = y_sq_bar + x_sq_bar * (2.0 * a * c + b * b)
                  + 2.0 *
                    ((a * b + c * d) * x3_bar + c * b * x_bar + a * d * x5_bar + b * d * x4_bar - a * y_x2_bar - b * xy_bar - c * y_bar - d * y_x3_bar)
                  + a * a * x4_bar + d * d * x6_bar + c * c;

    log_prob -= N * poly_term / (2.0 * sigma * sigma);

    return log_prob;
}

UserData third_order_user_data = {0};

BayLaModel third_order_model = 
    {
        .model_prior_probability = 1.0 / 7.0,
        .n_parameters = 5,
        .prior = third_order_parameters_log_prior,
        .likelihood = third_order_parameters_log_likelihood,
        .min_parameter_vals = third_order_min_parms,
        .max_parameter_vals = third_order_max_parms,
        .user_data = &third_order_user_data
    };

/*-----------------------------------------------  4th Order Polynomial  -------------------------------------------------*/
/* param 0 is a constant, param 1 is a linear coeff, param 2 is a quadratic coeff, ....  param 4 is a standard deviation. */
f64 fourth_order_min_parms[6] = {-2.0, -4.0, -2.0, -2.0, -2.0, 0.0};
f64 fourth_order_max_parms[6] = {+0.0, +4.0, +2.0, +2.0, +2.0, 20.0};

static f64 fourth_order_parameters_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    /* exponential with mean 1 for std-dev and uniform for all other parameters. */
    return log((1.0 / 2.0) * (1.0 / 8.0) * (1.0 / 4.0) * (1.0 / 4.0) * (1.0 / 4.0)) - parms[5] - bayla_log1mexp(fourth_order_max_parms[5]);
}

static f64 fourth_order_parameters_log_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 c = parms[0];
    f64 b = parms[1];
    f64 a = parms[2];
    f64 d = parms[3];
    f64 e = parms[4];

    f64 sigma = parms[5];

    f64 y_sq_bar = data->y_sq_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x4_bar = data->x4_bar;
    f64 x6_bar = data->x6_bar;
    f64 x8_bar = data->x8_bar;

    f64 y_bar = data->y_bar;
    f64 x_bar = data->x_bar;

    f64 xy_bar = data->xy_bar;
    f64 x3_bar = data->x3_bar;
    f64 x5_bar = data->x5_bar;
    f64 x7_bar = data->x7_bar;
    f64 y_x2_bar = data->y_x2_bar;
    f64 y_x3_bar = data->y_x3_bar;
    f64 y_x4_bar = data->y_x4_bar;

    data->ncalls++;

    f64 log_prob = -N * (log(sigma) + 0.5 * log(2.0 * ELK_PI));

    f64 poly_term = y_sq_bar + x_sq_bar * (2.0 * a * c + b * b)
                  + 2.0 *
                    (
                     (a * b + c * d) * x3_bar + c * b * x_bar + (e * b + a * d) * x5_bar + (e * c + b * d) * x4_bar
                     + e * d * x7_bar + e * a * x6_bar
                     - a * y_x2_bar - b * xy_bar - c * y_bar - d * y_x3_bar - e * y_x4_bar
                    )
                  + a * a * x4_bar + d * d * x6_bar + e * e * x8_bar + c * c;

    log_prob -= N * poly_term / (2.0 * sigma * sigma);

    return log_prob;
}

UserData fourth_order_user_data = {0};

BayLaModel fourth_order_model = 
    {
        .model_prior_probability = 1.0 / 7.0,
        .n_parameters = 6,
        .prior = fourth_order_parameters_log_prior,
        .likelihood = fourth_order_parameters_log_likelihood,
        .min_parameter_vals = fourth_order_min_parms,
        .max_parameter_vals = fourth_order_max_parms,
        .user_data = &fourth_order_user_data
    };

/*-----------------------------------------------  5th Order Polynomial  -------------------------------------------------*/
/* param 0 is a constant, param 1 is a linear coeff, param 2 is a quadratic coeff, ....  param 5 is a standard deviation. */
f64 fifth_order_min_parms[7] = {-2.0, -4.0, -2.0, -2.0, -2.0, -2.0, 0.0};
f64 fifth_order_max_parms[7] = {+0.0,  4.0, +2.0, +2.0, +2.0, +2.0, 20.0};

static f64 fifth_order_parameters_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    /* exponential with mean 1 for std-dev and uniform for all other parameters. */
    return log((1.0 / 2.0) * (1.0 / 8.0) * (1.0 / 4.0) * (1.0 / 4.0) * (1.0 / 4.0) * (1.0 / 4.0)) - parms[6] - bayla_log1mexp(fifth_order_max_parms[6]);
}

static f64 fifth_order_parameters_log_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 c = parms[0];
    f64 b = parms[1];
    f64 a = parms[2];
    f64 d = parms[3];
    f64 e = parms[4];
    f64 f = parms[5];

    f64 sigma = parms[6];

    f64 y_sq_bar = data->y_sq_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x4_bar = data->x4_bar;
    f64 x6_bar = data->x6_bar;
    f64 x8_bar = data->x8_bar;
    f64 x10_bar = data->x10_bar;

    f64 y_bar = data->y_bar;
    f64 x_bar = data->x_bar;

    f64 xy_bar = data->xy_bar;
    f64 x3_bar = data->x3_bar;
    f64 x5_bar = data->x5_bar;
    f64 x7_bar = data->x7_bar;
    f64 x9_bar = data->x9_bar;
    f64 y_x2_bar = data->y_x2_bar;
    f64 y_x3_bar = data->y_x3_bar;
    f64 y_x4_bar = data->y_x4_bar;
    f64 y_x5_bar = data->y_x5_bar;

    data->ncalls++;

    f64 log_prob = -N * (log(sigma) + 0.5 * log(2.0 * ELK_PI));

    f64 poly_term = y_sq_bar + x_sq_bar * (2.0 * a * c + b * b)
                  + 2.0 *
                    (
                     (a * b + c * d) * x3_bar + c * b * x_bar + (e * b + a * d + c * f) * x5_bar + (e * c + b * d) * x4_bar
                     + (a * f + e * d) * x7_bar + (b * f + e * a) * x6_bar + d * f * x8_bar + e * f * x9_bar
                     - a * y_x2_bar - b * xy_bar - c * y_bar - d * y_x3_bar - e * y_x4_bar - f * y_x5_bar
                    )
                  + a * a * x4_bar + d * d * x6_bar + e * e * x8_bar + f * f * x10_bar + c * c;

    log_prob -= N * poly_term / (2.0 * sigma * sigma);

    return log_prob;
}

UserData fifth_order_user_data = {0};

BayLaModel fifth_order_model = 
    {
        .model_prior_probability = 1.0 / 7.0,
        .n_parameters = 7,
        .prior = fifth_order_parameters_log_prior,
        .likelihood = fifth_order_parameters_log_likelihood,
        .min_parameter_vals = fifth_order_min_parms,
        .max_parameter_vals = fifth_order_max_parms,
        .user_data = &fifth_order_user_data
    };

/*------------------------------------------------  Initialize UserData  -------------------------------------------------*/
static inline void
initialize_global_data(void)
{
    BayLaNormalDistribution normal = bayla_normal_distribution_create(0.0, 0.25);
    ElkRandomState ran = elk_random_state_create(12);

    f64 sum_y_sq = 0.0;
    f64 sum_x_sq = 0.0;
    f64 sum_x4 = 0.0;
    f64 sum_x6 = 0.0;
    f64 sum_x8 = 0.0;
    f64 sum_x10 = 0.0;

    f64 sum_y = 0.0;
    f64 sum_x = 0.0;

    f64 sum_xy = 0.0;
    f64 sum_x3 = 0.0;
    f64 sum_x5 = 0.0;
    f64 sum_x7 = 0.0;
    f64 sum_x9 = 0.0;
    f64 sum_y_x2 = 0.0;
    f64 sum_y_x3 = 0.0;
    f64 sum_y_x4 = 0.0;
    f64 sum_y_x5 = 0.0;

    f64 sum_y_logx = 0.0;
    f64 sum_logx_sq = 0.0;

    /* Initialize the global data, with some noise! */
    for(size i = 0; i < NUM_DATA_POINTS; ++i)
    {
        // approximate log(x) from x = 0.1 to 3.0 */
        global_xs[i] = 0.1 + 2.9 * i / (NUM_DATA_POINTS - 1);
        global_ys[i] = natural_log_approximation_by_3rd_order_polynomial(global_xs[i], 0.9);
        f64 error = bayla_normal_distribution_random_deviate(&normal, &ran);
        global_ys[i] += error;

        sum_y_sq += global_ys[i] * global_ys[i];
        sum_x_sq += global_xs[i] * global_xs[i];
        sum_x4 += global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i];
        sum_x6 += global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i];
        sum_x8 += global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i];
        sum_x10 += global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i];

        sum_y += global_ys[i];
        sum_x += global_xs[i];

        sum_xy += global_xs[i] * global_ys[i];
        sum_x3 += global_xs[i] * global_xs[i] * global_xs[i];
        sum_x5 += global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i];
        sum_x7 += global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i];
        sum_x9 += global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i];
        sum_y_x2 += global_xs[i] * global_xs[i] * global_ys[i];
        sum_y_x3 += global_xs[i] * global_xs[i] * global_xs[i] * global_ys[i];
        sum_y_x4 += global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_ys[i];
        sum_y_x5 += global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_xs[i] * global_ys[i];

        f64 logx = log(global_xs[i]);
        sum_y_logx += global_ys[i] * logx;
        sum_logx_sq += logx * logx;
    }

    size N = NUM_DATA_POINTS;
    f64 y_sq_bar = sum_y_sq / NUM_DATA_POINTS;
    f64 x_sq_bar = sum_x_sq / NUM_DATA_POINTS;
    f64 x4_bar = sum_x4 / NUM_DATA_POINTS;
    f64 x6_bar = sum_x6 / NUM_DATA_POINTS;
    f64 x8_bar = sum_x8 / NUM_DATA_POINTS;
    f64 x10_bar = sum_x10 / NUM_DATA_POINTS;

    f64 y_bar = sum_y / NUM_DATA_POINTS;
    f64 x_bar = sum_x / NUM_DATA_POINTS;

    f64 xy_bar = sum_xy / NUM_DATA_POINTS;
    f64 x3_bar = sum_x3 / NUM_DATA_POINTS;
    f64 x5_bar = sum_x5 / NUM_DATA_POINTS;
    f64 x7_bar = sum_x7 / NUM_DATA_POINTS;
    f64 x9_bar = sum_x9 / NUM_DATA_POINTS;
    f64 y_x2_bar = sum_y_x2 / NUM_DATA_POINTS;
    f64 y_x3_bar = sum_y_x3 / NUM_DATA_POINTS;
    f64 y_x4_bar = sum_y_x4 / NUM_DATA_POINTS;
    f64 y_x5_bar = sum_y_x5 / NUM_DATA_POINTS;

    f64 y_logx_bar = sum_y_logx / NUM_DATA_POINTS;
    f64 logx_sq_bar = sum_logx_sq / NUM_DATA_POINTS;

    UserData const default_ud = 
        {
            .ncalls = 0,
            .xs = global_xs,
            .ys = global_ys,

            .N = N,

            .y_sq_bar = y_sq_bar,
            .x_sq_bar = x_sq_bar,
            .x4_bar = x4_bar,
            .x6_bar = x6_bar,
            .x8_bar = x8_bar,
            .x10_bar = x10_bar,

            .y_bar = y_bar,
            .x_bar = x_bar,

            .xy_bar = xy_bar,
            .x3_bar = x3_bar,
            .x5_bar = x5_bar,
            .x7_bar = x7_bar,
            .x9_bar = x9_bar,
            .y_x2_bar = y_x2_bar,
            .y_x3_bar = y_x3_bar,
            .y_x4_bar = y_x4_bar,
            .y_x5_bar = y_x5_bar,

            .y_logx_bar = y_logx_bar,
            .logx_sq_bar = logx_sq_bar
        };

    constant_user_data = default_ud;
    linear_user_data = default_ud;
    second_order_user_data = default_ud;
    third_order_user_data = default_ud;
    fourth_order_user_data = default_ud;
    fifth_order_user_data = default_ud;
    log_user_data = default_ud;
}

#undef NUM_DATA_POINTS

/*---------------------------------------------------  Test Models   -----------------------------------------------------*/

static BayLaIntegratorOptions global_opts[] =
    {
        {
            .strategy = BAYLA_INTEGRATE_MONTE_CARLO,
            .simple_monte_carlo.min_samples = 50000,
            .simple_monte_carlo.max_samples = 30000000,
            .simple_monte_carlo.samples_per_batch = 1000,
            .simple_monte_carlo.acceptable_abs_error = -1.0e-100,
            .simple_monte_carlo.acceptable_prop_error = 0.05,
            .simple_monte_carlo.seed = 43
        },

        {
            .strategy = BAYLA_INTEGRATE_MISER,
            .miser.min_points = 30,
            .miser.min_to_subdivide = 120,
            .miser.total_samples = 30000000,
            .miser.explore_factor = 0.1,
            .miser.acceptable_abs_error = -1.0e-100,
            .miser.acceptable_prop_error = 0.05,
            .miser.seed = 44
        },

        {
            .strategy = BAYLA_INTEGRATE_VEGAS,
            .vegas.min_total_samples = 50000,
            .vegas.max_total_samples = 30000000,
            .vegas.samples_per_refinement = 2000000,
            .vegas.acceptable_abs_error = -1.0e-10,
            .vegas.acceptable_prop_error = 0.05,
            .vegas.num_grid_cells = 1000,
            .vegas.alpha = 0.75,
            .vegas.seed = 45
        },

        {
            .strategy = BAYLA_INTEGRATE_VEGAS_MISER,
            .vegas_miser.min_total_samples = 50000,
            .vegas_miser.max_total_samples = 30000000,
            .vegas_miser.samples_per_refinement = 200000,
            .vegas_miser.min_points = 30,
            .vegas_miser.min_to_subdivide = 120,
            .vegas_miser.explore_factor = 0.1,
            .vegas_miser.acceptable_abs_error = -1.0e-10,
            .vegas_miser.acceptable_prop_error = 0.05,
            .vegas_miser.num_grid_cells = 1000,
            .vegas_miser.alpha = 0.75,
            .vegas_miser.seed = 46
        },

        {
            .strategy = BAYLA_INTEGRATE_VEGAS_PLUS,
            .vegas_plus.min_total_samples = 50000,
            .vegas_plus.max_total_samples = 30000000,
            .vegas_plus.samples_per_refinement = 2000000,
            .vegas_plus.acceptable_abs_error = -1.0e-10,
            .vegas_plus.acceptable_prop_error = 0.05,
            .vegas_plus.num_grid_cells = 1000,
            .vegas_plus.alpha = 0.75,
            .vegas_plus.beta = 0.75,
            .vegas_plus.seed = 47
        }
    };

typedef struct
{
    BayLaModel *model;
    BayLaIntegratorOptions opts;
    MagStaticArena scratch;
    BayLaEvidenceResult *evidence;
} EvidenceCalcThreadData;

void
evidence_thread_func(void *data)
{
    EvidenceCalcThreadData *td = data;
    *td->evidence = bayla_calculate_evidence(td->model, &td->opts, td->scratch);
}

static inline void
test_evidence_for_polynomial_degree(BayLaIntegratorOptions *opts)
{

    BayLaModel *models[] = 
        {
            &log_model, &constant_model, &linear_model,
            &second_order_model, &third_order_model, &fourth_order_model, &fifth_order_model,
        };

#define NUM_MODELS  ECO_ARRAY_SIZE(models)

    char *model_names[NUM_MODELS] = {"log", "constant", "linear", "2nd order", "3rd order", "4th order", "5th order",};
    BayLaEvidenceResult evidence[NUM_MODELS] = {0};

    initialize_global_data();

    EvidenceCalcThreadData tds[NUM_MODELS] = {0};
    CoyThread threads[NUM_MODELS] = {0};
    MagStaticArena scratches[NUM_MODELS] = {0};
    for(size m = 0; m < NUM_MODELS; ++m)
    {
        scratches[m] = mag_static_arena_allocate_and_create(ECO_MiB(32));
    }

    b32 success = true;
    for(size m = 0; m < NUM_MODELS; ++m)
    {
        tds[m].model = models[m];
        tds[m].opts = *opts;
        tds[m].scratch = mag_static_arena_borrow(&scratches[m]);
        tds[m].evidence = &evidence[m];

        success &= coy_thread_create(&threads[m], evidence_thread_func, &tds[m]);
        Assert(success);
    }

    f64 max_evidence = -INFINITY;
    size max_evidence_index = -1;
    for(size m = 0; m < NUM_MODELS; ++m)
    {
        success &= coy_thread_join(&threads[m]);
        Assert(success);
        coy_thread_destroy(&threads[m]);

        if(evidence[m].result.value.val > max_evidence)
        {
            max_evidence = evidence[m].result.value.val;
            max_evidence_index = m;
        }
    }

    size total_calls = 0;
    printf("\n%42s - %13s %10s %12s\n", "Evidence", "Total Samples", "N Used", "N-Calls");
    for(size m = 0; m < NUM_MODELS; ++m)
    {
        printf("%s%s %9s : %8.1e (%8.1e, %5.1lf%%) - %13td",
                m == max_evidence_index ? "*": " ",
                evidence[m].converged ? " ": "-",
                model_names[m],
                bayla_log_value_map_out_of_log_domain(evidence[m].result.value),
                bayla_log_value_map_out_of_log_domain(evidence[m].result.error),
                100.0 * bayla_log_value_map_out_of_log_domain(bayla_log_value_divide(evidence[m].result.error, evidence[m].result.value)),
                evidence[m].total_samples);
        
        switch(evidence[m].strategy)
        {
        case BAYLA_INTEGRATE_VEGAS:
        case BAYLA_INTEGRATE_VEGAS_PLUS:
        case BAYLA_INTEGRATE_VEGAS_MISER:
        case BAYLA_INTEGRATE_MISER:
            {
                printf(" %10td", evidence[m].samples_used);
            } break;

        case BAYLA_INTEGRATE_MONTE_CARLO:
            {
                printf("%10s", " ");
            } break;

        default: Panic();
        }

        UserData *ud = models[m]->user_data;
        printf(" %12td\n", ud->ncalls);
        total_calls += ud->ncalls;
    }
    printf("Total calls: %.0lf million\n\n\n", total_calls / 1000000.0);

    for(size m = 0; m < NUM_MODELS; ++m)
    {
#ifdef _MAG_TRACK_MEM_USAGE
        f64 pct_mem = mag_static_arena_max_ratio(&scratches[m]) * 100.0;
        b32 over_allocated = mag_static_arena_over_allocated(&scratches[m]);
        printf( "%s used %.2lf%% of scratch and scratch %s over allocated.\n",
                __func__, pct_mem, over_allocated ? "***WAS***": "was not");
#endif

        mag_static_arena_destroy(&scratches[m]);
    }

#undef NUM_MODELS

}

/*----------------------------------------------------- Consistency -------------------------------------------------------*/

#define NUM_TRIALS 30

static inline void
test_evidence_for_polynomial_degree_consistency(BayLaIntegratorOptions *opts)
{

    BayLaModel *models[] = 
        {
            &log_model, &constant_model, &linear_model,
            &second_order_model, &third_order_model, &fourth_order_model, &fifth_order_model,
        };

#define NUM_MODELS  ECO_ARRAY_SIZE(models)

    char *model_names[NUM_MODELS] = {"log", "constant", "linear", "2nd order", "3rd order", "4th order", "5th order",};
    BayLaEvidenceResult evidence[NUM_MODELS][NUM_TRIALS] = {0};

    MagStaticArena scratches[NUM_TRIALS] = {0};
    for(size cnt = 0; cnt < NUM_TRIALS; ++cnt)
    {
        scratches[cnt] = mag_static_arena_allocate_and_create(ECO_MiB(32));
    }

    initialize_global_data();

    for(size m = 0; m < NUM_MODELS; ++m)
    {

        EvidenceCalcThreadData tds[NUM_TRIALS] = {0};
        CoyThread threads[NUM_TRIALS] = {0};

        b32 success = true;
        for(size cnt = 0; cnt < NUM_TRIALS; ++cnt)
        {

            tds[cnt].model = models[m];
            tds[cnt].opts = *opts;

            u64 delta = cnt * ECO_ARRAY_SIZE(models);
            switch(tds[cnt].opts.strategy)
            {
                case BAYLA_INTEGRATE_VEGAS:       tds[cnt].opts.vegas.seed += delta;              break;
                case BAYLA_INTEGRATE_VEGAS_PLUS:  tds[cnt].opts.vegas_plus.seed += delta;         break;
                case BAYLA_INTEGRATE_VEGAS_MISER: tds[cnt].opts.vegas_miser.seed += delta;        break;
                case BAYLA_INTEGRATE_MISER:       tds[cnt].opts.miser.seed += delta;              break;
                case BAYLA_INTEGRATE_MONTE_CARLO: tds[cnt].opts.simple_monte_carlo.seed += delta; break;
                default: Panic();                                                                 break;
            }

            tds[cnt].evidence = &evidence[m][cnt];
            tds[cnt].scratch = mag_static_arena_borrow(&scratches[cnt]);

            success &= coy_thread_create(&threads[cnt], evidence_thread_func, &tds[cnt]);
            Assert(success);
        }

        for(size cnt = 0; cnt < NUM_TRIALS; ++cnt)
        {
            success &= coy_thread_join(&threads[cnt]);
            coy_thread_destroy(&threads[cnt]);
        }
        Assert(success);

        f64 min_err = INFINITY;
        f64 max_err = -INFINITY;
        f64 sum = 0.0;
        for(size cnt = 0; cnt < NUM_TRIALS; ++cnt)
        {
            f64 err = bayla_log_value_map_out_of_log_domain(evidence[m][cnt].result.error);
            min_err = min_err < err ? min_err : err;
            max_err = max_err > err ? max_err : err;

            sum += bayla_log_value_map_out_of_log_domain(evidence[m][cnt].result.value);
        }
        f64 mean = sum / NUM_TRIALS;

        f64 sum_diff_sq = 0.0;
        for(size cnt = 0; cnt < NUM_TRIALS; ++cnt)
        {
            f64 value = bayla_log_value_map_out_of_log_domain(evidence[m][cnt].result.value);
            sum_diff_sq += (value - mean) * (value - mean); 
        }
        f64 std = sqrt(sum_diff_sq / (NUM_TRIALS - 1));
        printf("%10s %8.6e ± %8.6e [%8.6e <-> %8.6e]\n", model_names[m], mean, std, min_err, max_err);
    }

#ifdef _MAG_TRACK_MEM_USAGE
    printf("\n\n");
#endif
    for(size cnt = 0; cnt < NUM_TRIALS; ++cnt)
    {

#ifdef _MAG_TRACK_MEM_USAGE

        f64 pct_mem = mag_static_arena_max_ratio(&scratches[cnt]) * 100.0;
        b32 over_allocated = mag_static_arena_over_allocated(&scratches[cnt]);
        printf( "%s used %.2lf%% of scratch and scratch %s over allocated.\n",
                __func__, pct_mem, over_allocated ? "***WAS***": "was not");

#endif

        mag_static_arena_destroy(&scratches[cnt]);
    }

#undef NUM_MODELS

}

#undef NUM_TRIALS

/*-------------------------------------------  Preconditioning Sampling  --------------------------------------------------*/

static inline void
test_sampling_preconditioning(BayLaIntegratorOptions *opts)
{

    BayLaModel *models[] = 
        {
            &log_model, &constant_model, &linear_model,
            &second_order_model, &third_order_model, &fourth_order_model, &fifth_order_model,
        };

#define NUM_MODELS ECO_ARRAY_SIZE(models)

    char *model_names[NUM_MODELS] = {"log", "constant", "linear", "2nd order", "3rd order", "4th order", "5th order",};
    BayLaEvidenceResult evidence[NUM_MODELS] = {0};

    MagStaticArena scratch = mag_static_arena_allocate_and_create(ECO_MiB(32));
    MagStaticArena perm_ = mag_static_arena_allocate_and_create(ECO_MiB(32));
    MagStaticArena *perm = &perm_;

    ElkRandomState ran = elk_random_state_create(11);
    initialize_global_data();

    f64 starting_values[NUM_MODELS] = {0};

    f64 max_evidence = -INFINITY;
    size max_evidence_index = -1;
    for(size m = 0; m < NUM_MODELS; ++m)
    {
        BayLaModel *model = models[m];

        f64 *min_vals = model->min_parameter_vals;
        f64 *max_vals = model->max_parameter_vals;
        for(size n = 0; n < model->n_parameters; ++n)
        {
            starting_values[n] = (min_vals[n] + max_vals[n]) / 2.0;
        }

        CoyProfileAnchor ap = COY_START_PROFILE_BLOCK("Precondition Sampling.");
        BayLaPreconditioningSamplerArgs args = bayla_preconditioning_sampler_args_allocate_and_create(
                model,
                1000,
                10000,
                starting_values,
                perm);

        bayla_preconditioning_sample(&args, &ran, mag_static_arena_borrow(&scratch));
        COY_END_PROFILE(ap);

        ap = COY_START_PROFILE_BLOCK("Preconditioning.");
        BayLaParametersSamples *samples = args.output;
        bayla_vegas_map_precondition(model, opts, samples, perm);
        COY_END_PROFILE(ap);

        evidence[m] = bayla_calculate_evidence(model, opts, mag_static_arena_borrow(&scratch));
        if(evidence[m].result.value.val > max_evidence)
        {
            max_evidence = evidence[m].result.value.val;
            max_evidence_index = m;
        }
    }

    size total_calls = 0;
    printf("\n%42s - %13s %10s %12s\n", "Evidence", "Total Samples", "N Used", "N-Calls");
    for(size m = 0; m < NUM_MODELS; ++m)
    {
        printf("%s%s %9s : %8.1e (%8.1e, %5.1lf%%) - %13td",
                m == max_evidence_index ? "*": " ",
                evidence[m].converged ? " ": "-",
                model_names[m],
                bayla_log_value_map_out_of_log_domain(evidence[m].result.value),
                bayla_log_value_map_out_of_log_domain(evidence[m].result.error),
                100.0 * bayla_log_value_map_out_of_log_domain(bayla_log_value_divide(evidence[m].result.error, evidence[m].result.value)),
                evidence[m].total_samples);

        switch(evidence[m].strategy)
        {
        case BAYLA_INTEGRATE_VEGAS:
        case BAYLA_INTEGRATE_VEGAS_PLUS:
        case BAYLA_INTEGRATE_VEGAS_MISER:
        case BAYLA_INTEGRATE_MISER:
            {
                printf(" %10td", evidence[m].samples_used);
            } break;

        case BAYLA_INTEGRATE_MONTE_CARLO:
            {
                printf("%10s", " ");
            } break;

        default: Panic();
        }

        UserData *ud = models[m]->user_data;
        printf(" %12td\n", ud->ncalls);
        total_calls += ud->ncalls;
    }
    printf("Total calls: %.0lf million\n", total_calls / 1000000.0);

#ifdef _MAG_TRACK_MEM_USAGE
    f64 pct_mem = mag_static_arena_max_ratio(&scratch) * 100.0;
    b32 over_allocated = mag_static_arena_over_allocated(&scratch);
    printf( "\n\n%s used %.2lf%% of scratch and scratch %s over allocated.\n",
            __func__, pct_mem, over_allocated ? "***WAS***": "was not");

    pct_mem = mag_static_arena_max_ratio(perm) * 100.0;
    over_allocated = mag_static_arena_over_allocated(perm);
    printf( "%s used %.2lf%% of perm and perm %s over allocated.\n\n\n",
            __func__, pct_mem, over_allocated ? "***WAS***": "was not");
#endif

    mag_static_arena_destroy(&scratch);

#undef NUM_MODELS

}

/*---------------------------------------------------   Pi Evidence  -----------------------------------------------------*/

f64 pi_min_parms[2] = {-1.0, -1.0};
f64 pi_max_parms[2] = { 1.0,  1.0};

static f64 pi_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    return 0;
}

static f64 pi_log_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    f64 x = parms[0];
    f64 y = parms[1];

    if(x*x + y*y <= 1.0) return 0;
    return -INFINITY;
}

BayLaModel pi_model = 
    {
        .model_prior_probability = 1.0,
        .n_parameters = 2,
        .prior = pi_log_prior,
        .likelihood = pi_log_likelihood,
        .min_parameter_vals = pi_min_parms,
        .max_parameter_vals = pi_max_parms,
        .user_data = NULL
    };

typedef struct
{
    BayLaModel *model;
    BayLaIntegratorOptions *opts;
    MagStaticArena scratch;
    BayLaEvidenceResult evidence;
    f64 true_error;
    f64 true_error_z_score;
} TestEvidencePiThreadData;

void
test_evidence_pi_thread_func(void *data)
{
    TestEvidencePiThreadData *td = data;
    td->evidence = bayla_calculate_evidence(td->model, td->opts, td->scratch);

    td->true_error = bayla_log_value_map_out_of_log_domain(td->evidence.result.value) - ELK_PI;
    td->true_error_z_score = td->true_error / bayla_log_value_map_out_of_log_domain(td->evidence.result.error);
}

static inline void
test_evidence_pi(void)
{
    CoyProfileAnchor ap = COY_START_PROFILE_FUNCTION;

    char *names[] = {"Simple Monte Carlo", "MISER", "VEGAS", "MISER-VEGAS", "VEGAS+"};

    BayLaIntegratorOptions opts[]= 
        {
            (BayLaIntegratorOptions)
            {
                .strategy = BAYLA_INTEGRATE_MONTE_CARLO,
                .simple_monte_carlo.min_samples = 10000,
                .simple_monte_carlo.max_samples = 20000000,
                .simple_monte_carlo.samples_per_batch = 10000,
                .simple_monte_carlo.acceptable_abs_error = -1.0e-10,
                .simple_monte_carlo.acceptable_prop_error = 0.0001,
                .simple_monte_carlo.seed = 43
            },

            (BayLaIntegratorOptions)
            {
                .strategy = BAYLA_INTEGRATE_MISER,
                .miser.min_points = 30,
                .miser.min_to_subdivide = 120,
                .miser.total_samples = 20000000,
                .miser.explore_factor = 0.1,
                .miser.acceptable_abs_error = -1.0e-10,
                .miser.acceptable_prop_error = 0.0001,
                .miser.seed = 44
            },

            (BayLaIntegratorOptions)
            {
                .strategy = BAYLA_INTEGRATE_VEGAS,
                .vegas.min_total_samples = 10000,
                .vegas.max_total_samples = 20000000,
                .vegas.samples_per_refinement = 10000,
                .vegas.num_grid_cells = 1000,
                .vegas.acceptable_abs_error = -1.0e-10,
                .vegas.acceptable_prop_error = 0.0001,
                .vegas.alpha = 0.75,
                .vegas.seed = 45
            },

            (BayLaIntegratorOptions)
            {
                .strategy = BAYLA_INTEGRATE_VEGAS_MISER,
                .vegas_miser.min_total_samples = 10000,
                .vegas_miser.max_total_samples = 20000000,
                .vegas_miser.samples_per_refinement = 10000,
                .vegas_miser.num_grid_cells = 1000,
                .vegas_miser.min_points = 30,
                .vegas_miser.min_to_subdivide = 120,
                .vegas_miser.explore_factor = 0.1,
                .vegas_miser.acceptable_abs_error = -1.0e-10,
                .vegas_miser.acceptable_prop_error = 0.0001,
                .vegas_miser.alpha = 0.75,
                .vegas_miser.seed = 46
            },

            (BayLaIntegratorOptions)
            {
                .strategy = BAYLA_INTEGRATE_VEGAS_PLUS,
                .vegas_plus.min_total_samples = 10000,
                .vegas_plus.max_total_samples = 20000000,
                .vegas_plus.samples_per_refinement = 10000,
                .vegas_plus.num_grid_cells = 1000,
                .vegas_plus.acceptable_abs_error = -1.0e-10,
                .vegas_plus.acceptable_prop_error = 0.0001,
                .vegas_plus.alpha = 0.75,
                .vegas_plus.beta = 0.75,
                .vegas_plus.seed = 47
            }
        };

    _Static_assert(ECO_ARRAY_SIZE(opts) == ECO_ARRAY_SIZE(names), "Mismatched Arrays.");

    TestEvidencePiThreadData tds[ECO_ARRAY_SIZE(opts)] = {0};
    CoyThread threads[ECO_ARRAY_SIZE(opts)] = {0};
    MagStaticArena scratches[ECO_ARRAY_SIZE(opts)] = {0};
    for(size method = 0; method < ECO_ARRAY_SIZE(opts); ++method)
    {
        scratches[method] = mag_static_arena_allocate_and_create(ECO_KiB(512));
    }
    
    b32 success = true;
    for(size method = 0; method < ECO_ARRAY_SIZE(opts); ++method)
    {
        tds[method].model = &pi_model;
        tds[method].opts = &opts[method];
        tds[method].scratch = mag_static_arena_borrow(&scratches[method]);

        success &= coy_thread_create(&threads[method], test_evidence_pi_thread_func, &tds[method]);
        Assert(success);
    }

    printf("\nCalculating Pi as a test of the integrators.\n\n");
    printf("%20s - %8s ± %8s [%9s] %8s [%9s] %8s [ %8s ]\n",
            "Method", "Value", "Error", "Pct Err", "True Err", "Pct True", "Samples", "Samples Used");

    for(size method = 0; method < ECO_ARRAY_SIZE(opts); ++method)
    {
        success = coy_thread_join(&threads[method]);
        Assert(success);
        coy_thread_destroy(&threads[method]);

        printf("%20s - %lf ± %lf [%8.4lf%%] %8.4lf [%8.4lf%%]",
                names[method],
                bayla_log_value_map_out_of_log_domain(tds[method].evidence.result.value), 
                bayla_log_value_map_out_of_log_domain(tds[method].evidence.result.error),
                100.0 * bayla_log_value_map_out_of_log_domain(bayla_log_value_divide(
                        tds[method].evidence.result.error, tds[method].evidence.result.value)),
                tds[method].true_error_z_score, tds[method].true_error / ELK_PI * 100.0);

        switch(tds[method].evidence.strategy)
        {
            case BAYLA_INTEGRATE_MONTE_CARLO:
            {
                printf(" %8td [ %12td ]\n", tds[method].evidence.total_samples, tds[method].evidence.total_samples);
            } break;

            case BAYLA_INTEGRATE_VEGAS_PLUS:
            case BAYLA_INTEGRATE_VEGAS:
            case BAYLA_INTEGRATE_MISER:
            case BAYLA_INTEGRATE_VEGAS_MISER:
            {
                printf(" %8td [ %12td ]\n", tds[method].evidence.total_samples, tds[method].evidence.samples_used);
            } break;

            default: Panic();
        }
    }
    printf("\n\n");

    for(size method = 0; method < ECO_ARRAY_SIZE(opts); ++method)
    {
#ifdef _MAG_TRACK_MEM_USAGE
        f64 pct_mem = mag_static_arena_max_ratio(&scratches[method]) * 100.0;
        b32 over_allocated = mag_static_arena_over_allocated(&scratches[method]);
        printf( "%s used %.2lf%% of scratch and scratch %s over allocated.\n",
                __func__, pct_mem, over_allocated ? "***WAS***": "was not");
#endif

        mag_static_arena_destroy(&scratches[method]);
    }

    COY_END_PROFILE(ap);
}

/*---------------------------------------------------   All Tests    -----------------------------------------------------*/
static inline void
all_evidence_tests(void)
{
    CoyProfileAnchor ap = {0};

    test_evidence_pi();

    printf("\nTesting Monte Carlo\n");
    ap = COY_START_PROFILE_BLOCK("Monte Carlo");
    test_evidence_for_polynomial_degree(&global_opts[0]);
    COY_END_PROFILE(ap);

    printf("\n");
    ap = COY_START_PROFILE_BLOCK("Monte Carlo Consistency");
    test_evidence_for_polynomial_degree_consistency(&global_opts[0]);
    COY_END_PROFILE(ap);

    printf("\nTesting MISER\n");
    ap = COY_START_PROFILE_BLOCK("MISER");
    test_evidence_for_polynomial_degree(&global_opts[1]);
    COY_END_PROFILE(ap);

    printf("\n");
    ap = COY_START_PROFILE_BLOCK("MISER Consistency");
    test_evidence_for_polynomial_degree_consistency(&global_opts[1]);
    COY_END_PROFILE(ap);

    printf("\nTesting VEGAS\n");
    ap = COY_START_PROFILE_BLOCK("VEGAS");
    test_evidence_for_polynomial_degree(&global_opts[2]);
    COY_END_PROFILE(ap);

    printf("\n");
    ap = COY_START_PROFILE_BLOCK("VEGAS Consistency");
    test_evidence_for_polynomial_degree_consistency(&global_opts[2]);
    COY_END_PROFILE(ap);

    printf("\nTesting Sampling Preconditioned VEGAS\n");
    ap = COY_START_PROFILE_BLOCK("Preconditioned VEGAS");
    test_sampling_preconditioning(&global_opts[2]);
    COY_END_PROFILE(ap);

    printf("\nTesting VEGAS-MISER\n");
    ap = COY_START_PROFILE_BLOCK("VEGAS-MISER");
    test_evidence_for_polynomial_degree(&global_opts[3]);
    COY_END_PROFILE(ap);

    printf("\n");
    ap = COY_START_PROFILE_BLOCK("VEGAS-MISER Consistency");
    test_evidence_for_polynomial_degree_consistency(&global_opts[3]);
    COY_END_PROFILE(ap);

    printf("\nTesting Sampling Preconditioned VEGAS-MISER\n");
    ap = COY_START_PROFILE_BLOCK("Preconditioned VEGAS-MISER");
    test_sampling_preconditioning(&global_opts[3]);
    COY_END_PROFILE(ap);

    printf("\nTesting VEGAS+\n");
    ap = COY_START_PROFILE_BLOCK("VEGAS+");
    test_evidence_for_polynomial_degree(&global_opts[4]);
    COY_END_PROFILE(ap);

    printf("\n");
    ap = COY_START_PROFILE_BLOCK("VEGAS+ Consistency");
    test_evidence_for_polynomial_degree_consistency(&global_opts[4]);
    COY_END_PROFILE(ap);

    printf("\nTesting Sampling Preconditioned VEGAS+\n");
    ap = COY_START_PROFILE_BLOCK("Preconditioned VEGAS+");
    test_sampling_preconditioning(&global_opts[4]);
    COY_END_PROFILE(ap);
}

