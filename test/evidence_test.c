/*-------------------------------------------------- Data Gen Function ----------------------------------------------------*/

static inline f64
natural_log_approximation_by_3rd_order_polynomial(f64 x, f64 x0)
{
    /* This is a 3rd order polynomial approximation to ln around the point x0 */
    f64 d = x - x0;
    return log(x0) + d / x0 - (d * d) / (2 * x0 * x0) + (d * d * d) / ( 3.0 * x0 * x0 * x0);
}

/* Some data used by all models, this will be initialized in the test. */
#define NUM_DATA_POINTS 40
typedef struct
{
    /* For keeping track of how many times a function is called. */
    size ncalls;

    /* The original data used to compute the stats. */
    f64 *xs;
    f64 *ys;

    /* Useful statistics for calculating the probability of the data given the parameters. */
    size N;

    f64 x_bar;
    f64 x_sq_bar;
    f64 x3_bar;
    f64 x4_bar;
    f64 x5_bar;
    f64 x6_bar;
    f64 x7_bar;
    f64 x8_bar;
    f64 x9_bar;
    f64 x10_bar;

    f64 y_sq_bar;

    f64 y_bar;
    f64 xy_bar;
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
static f64 const log_model_min_sigma = 1.0e-5;
static f64 const log_model_max_sigma = 15.0;

/* parameter 0 a standard deviation. */
static f64
log_model_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    f64 sigma = parms[0];
    if(sigma < log_model_min_sigma || sigma > log_model_max_sigma)
    {
        return -INFINITY;
    }

    return -log(sigma) - log(log(log_model_max_sigma / log_model_min_sigma));
}

static f64
log_model_log_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 sigma = parms[0];
    f64 y_sq_bar = data->y_sq_bar;
    f64 y_logx_bar = data->y_logx_bar;
    f64 logx_sq_bar = data->logx_sq_bar;

    data->ncalls++;

    f64 Q = logx_sq_bar - 2.0 * y_logx_bar + y_sq_bar;
    f64 log_prob =  -N * 0.5 * log(2.0 * ELK_PI) - N * sigma - 0.5 * N * Q / (sigma * sigma);

    return log_prob;
}

static void
log_model_max_a_posteriori(void *user_data, size num_params, f64 *params_out)
{
    Assert(num_params == 1);

    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 y_sq_bar = data->y_sq_bar;
    f64 y_logx_bar = data->y_logx_bar;
    f64 logx_sq_bar = data->logx_sq_bar;

    f64 Q = logx_sq_bar - 2.0 * y_logx_bar + y_sq_bar;
    params_out[0] = sqrt(N * Q / (N + 1.0));

    for(size i = 0; i < num_params; ++i)
    {
        Assert(params_out[i] >= log_model_min_sigma && params_out[i] <= log_model_max_sigma);
    }

#if 0
    printf("Log Model max-a-postiori: (σ = %e)\n", params_out[0]);
#endif
}

static BayLaSquareMatrix
log_model_2d_hessian(void *user_data, size num_params, f64 const *max_a_posteriori_params, MagAllocator *alloc)
{
    Assert(num_params == 1);

    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 sigma = max_a_posteriori_params[0];
    f64 y_sq_bar = data->y_sq_bar;
    f64 y_logx_bar = data->y_logx_bar;
    f64 logx_sq_bar = data->logx_sq_bar;

    f64 Q = logx_sq_bar - 2.0 * y_logx_bar + y_sq_bar;

    BayLaSquareMatrix hess = bayla_square_matrix_create(1, alloc);
    hess.data[0] = pow(2.0 * ELK_PI, - N / 2.0) / log(log_model_max_sigma / log_model_min_sigma);
    f64 sigma2 = sigma * sigma;
    f64 sigma4 = sigma2 * sigma2;
    hess.data[0] *= pow(sigma, -(N + 1)) * exp(- 0.5 * N * Q / sigma2) * ((N + 1) / sigma2 - (3.0 * N * Q / sigma4));

    Assert(hess.data[0] < 0.0);

    return hess;
}

UserData log_user_data = {0};

BayLaModel log_model = 
    {
        .model_prior_probability = 1.0 / 7.0,
        .n_parameters = 1,
        .prior = log_model_log_prior,
        .likelihood = log_model_log_likelihood,
        .posterior_max = log_model_max_a_posteriori,
        .hessian_matrix = log_model_2d_hessian,
        .user_data = &log_user_data
    };

/*--------------------------------------------------- Constant Model -----------------------------------------------------*/
/* parameter 0 is a constant value, parameter 1 is a standard deviation. */
f64 constant_model_min_parms[2] = {-2.0, log_model_min_sigma };
f64 constant_model_max_parms[2] = { 2.0, log_model_max_sigma };

static f64 
constant_model_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    f64 v0 = parms[0];
    f64 v0_min = constant_model_min_parms[0];
    f64 v0_max = constant_model_max_parms[0];

    f64 sigma = parms[1];
    f64 sigma_min = constant_model_min_parms[1];
    f64 sigma_max = constant_model_max_parms[1];

    if(sigma < sigma_min || sigma > sigma_max || v0 < v0_min || v0 > v0_max)
    {
        return -INFINITY;
    }

    return -log(sigma * log(sigma_max / sigma_min) * (v0_max - v0_min));
}

static f64 
constant_model_log_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    f64 v0 = parms[0];
    f64 sigma = parms[1];

    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 y_bar = data->y_bar;
    f64 y_sq_bar = data->y_sq_bar;

    data->ncalls++;

    f64 log_prob = -N * (log(sigma) + 0.5 * log(2.0 * ELK_PI));
    log_prob -= N * (y_sq_bar - 2.0 * y_bar * v0 + v0 * v0) / (2.0 * sigma * sigma);

    return log_prob;
}

static void
constant_model_max_a_posteriori(void *user_data, size num_params, f64 *params_out)
{
    Assert(num_params == 2);

    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 y_bar = data->y_bar;
    f64 y_sq_bar = data->y_sq_bar;

    params_out[0] = y_bar;
    params_out[1] = sqrt((N / (N + 1.0)) * (y_sq_bar - y_bar * y_bar));

    for(size i = 0; i < num_params; ++i)
    {
        Assert(params_out[i] >= constant_model_min_parms[i] && params_out[i] <= constant_model_max_parms[i]);
    }

#if 0
    printf("Constant Model max-a-postiori: %e (σ = %e)\n", params_out[0],  params_out[1]);
#endif
}

static BayLaSquareMatrix
constant_model_2d_hessian(void *user_data, size num_params, f64 const *max_a_posteriori_params, MagAllocator *alloc)
{
    Assert(num_params == 2);

    f64 v0 = max_a_posteriori_params[0];
    f64 v0_min = constant_model_min_parms[0];
    f64 v0_max = constant_model_max_parms[0];

    f64 sigma = max_a_posteriori_params[1];
    f64 sigma_min = constant_model_min_parms[1];
    f64 sigma_max = constant_model_max_parms[1];

    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 y_bar = data->y_bar;
    f64 y_sq_bar = data->y_sq_bar;

    f64 Q = v0 * v0 - 2.0 * v0 * y_bar + y_sq_bar;

    BayLaSquareMatrix hess = bayla_square_matrix_create(2, alloc);

    f64 prior = -log(sigma * log(sigma_max / sigma_min) * (v0_max - v0_min));
    f64 likelihood = -N * (log(sigma) + 0.5 * log(2.0 * ELK_PI));
    likelihood -= N * (y_sq_bar - 2.0 * y_bar * v0 + v0 * v0) / (2.0 * sigma * sigma);
    f64 p = exp(prior + likelihood);
    f64 sigma2 = sigma * sigma;
    f64 sigma4 = sigma2 * sigma2;

    hess.data[MATRIX_IDX(0, 0, 2)] = -N / sigma2 * p;
    hess.data[MATRIX_IDX(1, 1, 2)] = ((N + 1.0) / sigma2 - 3.0 * N * Q / sigma4) * p;

    Assert(hess.data[MATRIX_IDX(0, 0, 2)] < 0.0 && hess.data[MATRIX_IDX(1, 1, 2)] < 0.0 );

    return hess;
}

UserData constant_user_data = {0};

BayLaModel constant_model = 
    {
        .model_prior_probability = 1.0 / 7.0,
        .n_parameters = 2,
        .prior = constant_model_log_prior,
        .likelihood = constant_model_log_likelihood,
        .posterior_max = constant_model_max_a_posteriori,
        .hessian_matrix = constant_model_2d_hessian,
        .user_data = &constant_user_data
    };

/*---------------------------------------------------  Linear Model  -----------------------------------------------------*/
/* parameter 0 is a constant value, parameter 1 is a linear coefficent, parameter 2 is a standard deviation. */
f64 linear_model_min_parms[3] = {-2.0, -0.0, log_model_min_sigma };
f64 linear_model_max_parms[3] = { 2.0,  4.0, log_model_max_sigma };

static f64 
linear_model_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    f64 v0 = parms[0];
    f64 v0_min = linear_model_min_parms[0];
    f64 v0_max = linear_model_max_parms[0];

    f64 b = parms[1];
    f64 b_min = linear_model_min_parms[1];
    f64 b_max = linear_model_max_parms[1];

    f64 sigma = parms[2];
    f64 sigma_min = linear_model_min_parms[2];
    f64 sigma_max = linear_model_max_parms[2];

    if(sigma < sigma_min || sigma > sigma_max || v0 < v0_min || v0 > v0_max || b < b_min || b > b_max)
    {
        return -INFINITY;
    }

    return -log(sigma * log(sigma_max / sigma_min) * (v0_max - v0_min) * (b_max - b_min));
}

static f64 
linear_model_log_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    f64 v0 = parms[0];
    f64 b = parms[1];
    f64 sigma = parms[2];

    UserData *data = user_data;
    f64 N = (f64)data->N;

    f64 x_bar = data->x_bar;
    f64 x_sq_bar = data->x_sq_bar;

    f64 y_bar = data->y_bar;
    f64 y_sq_bar = data->y_sq_bar;

    f64 xy_bar = data->xy_bar;

    data->ncalls++;
    
    f64 Q = v0 * v0 + 2.0 * (v0 * b * x_bar - v0 * y_bar - b * xy_bar) + b * b *x_sq_bar + y_sq_bar;
    f64 log_prob = -(N + 1.0) * log(sigma) - 0.5 * N * log(2.0 * ELK_PI) - N * Q / (2.0 * sigma * sigma);

    return log_prob;
}

static void
linear_model_max_a_posteriori(void *user_data, size num_parms, f64 *parms_out)
{
    Assert(num_parms == 3);

    byte buffer2[sizeof(i32) * 2 * 3 + sizeof(f64) * 4 + _Alignof(f64)] = {0};
    MagStaticArena scratch_ = mag_static_arena_create(sizeof(buffer2) / sizeof(buffer2[0]), buffer2);
    MagAllocator scratch = mag_allocator_from_static_arena(&scratch_);

    UserData *data = user_data;
    f64 N = (f64)data->N;

    f64 x_bar = data->x_bar;
    f64 x_sq_bar = data->x_sq_bar;

    f64 y_bar = data->y_bar;
    f64 y_sq_bar = data->y_sq_bar;

    f64 xy_bar = data->xy_bar;

    f64 m_data[4] = 
        {
                 1.0,    x_bar,
               x_bar, x_sq_bar
        };
    BayLaSquareMatrix m = { .ndim = 2, .data = m_data };

    f64 input[2] = { y_bar, xy_bar };

    bayla_square_matrix_solve(m, input, scratch);
    memcpy(parms_out, input, sizeof(f64) * 2);

    f64 v0 = parms_out[0];
    f64 b = parms_out[1];

    f64 Q = v0 * v0 + 2.0 * (v0 * b * x_bar - v0 * y_bar - b * xy_bar) + b * b *x_sq_bar + y_sq_bar;
    parms_out[2] = sqrt(N / (N + 1.0) * Q);

    for(size i = 0; i < num_parms; ++i)
    {
        Assert(parms_out[i] >= linear_model_min_parms[i] && parms_out[i] <= linear_model_max_parms[i]);
    }

#if 0
    printf("Linear Model max-a-postiori: %e + %e * x (σ = %e)\n", params_out[0],  params_out[1], params_out[2]);
#endif
}

static BayLaSquareMatrix
linear_model_2d_hessian(void *user_data, size num_parms, f64 const *max_a_posteriori_parms, MagAllocator *alloc)
{
    Assert(num_parms == 3);

    f64 v0 = max_a_posteriori_parms[0];
    f64 v0_min = linear_model_min_parms[0];
    f64 v0_max = linear_model_max_parms[0];
    Assert(v0 > v0_min && v0 < v0_max);

    f64 b = max_a_posteriori_parms[1];
    f64 b_min = linear_model_min_parms[1];
    f64 b_max = linear_model_max_parms[1];
    Assert(b > b_min && b < b_max);

    f64 sigma = max_a_posteriori_parms[2];
    f64 sigma_min = linear_model_min_parms[2];
    f64 sigma_max = linear_model_max_parms[2];
    Assert(sigma > sigma_min && sigma < sigma_max);

    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 y_bar = data->y_bar;
    f64 y_sq_bar = data->y_sq_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x_bar = data->x_bar;
    f64 xy_bar = data->xy_bar;

    f64 Q = v0 * v0 + 2.0 * (v0 * b * x_bar - v0 * y_bar - b * xy_bar) + b * b *x_sq_bar + y_sq_bar;

    BayLaSquareMatrix hess = bayla_square_matrix_create(3, alloc);

    f64 prior = -log(sigma * log(sigma_max / sigma_min) * (v0_max - v0_min) * (b_max - b_min));
    f64 likelihood = -(N + 1.0) * log(sigma) - 0.5 * N * log(2.0 * ELK_PI) - N * Q / (2.0 * sigma * sigma);
    f64 p = exp(prior + likelihood);
    f64 sigma2 = sigma * sigma;
    f64 sigma4 = sigma2 * sigma2;
    f64 beta = (xy_bar - x_bar * y_bar) / (x_sq_bar - x_bar * x_bar);

    hess.data[MATRIX_IDX(0, 0, 3)] = -N / sigma2;
    hess.data[MATRIX_IDX(0, 1, 3)] = hess.data[MATRIX_IDX(1, 0, 3)] = -N * x_bar / sigma2;

    hess.data[MATRIX_IDX(1, 1, 3)] = -N *x_sq_bar / sigma2;
    hess.data[MATRIX_IDX(1, 2, 3)] = hess.data[MATRIX_IDX(2, 1, 3)] = 2.0 * N / sigma2 / sigma * (x_bar * y_bar - xy_bar + beta * (x_sq_bar - x_bar * x_bar));

    hess.data[MATRIX_IDX(2, 2, 3)] = (N +  1.0) / sigma2 - 3.0 * N * Q / sigma4;
    hess.data[MATRIX_IDX(0, 2, 3)] = hess.data[MATRIX_IDX(2, 0, 3)] = 0;

    for(size i = 0; i < 9; ++i) { hess.data[i] *= p; }

    Assert(hess.data[MATRIX_IDX(0, 0, 3)] < 0.0 && hess.data[MATRIX_IDX(1, 1, 3)] < 0.0 && hess.data[MATRIX_IDX(2, 2, 3)] < 0.0);

    return hess;
}

UserData linear_user_data = {0};

BayLaModel linear_model = 
    {
        .model_prior_probability = 1.0 / 7.0,
        .n_parameters = 3,
        .prior = linear_model_log_prior,
        .likelihood = linear_model_log_likelihood,
        .posterior_max = linear_model_max_a_posteriori,
        .hessian_matrix = linear_model_2d_hessian,
        .user_data = &linear_user_data
    };

/*-----------------------------------------------  2nd Order Polynomial  -------------------------------------------------*/
/* param 0 is a constant, param 1 is a linear coeff, param 2 is a quadratic coeff, param 3 is a standard deviation. */
f64 second_order_min_parms[4] = {-4.0, -4.0, -4.0, log_model_min_sigma };
f64 second_order_max_parms[4] = {+4.0, +4.0, +4.0, log_model_max_sigma };

static f64 
second_order_model_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    f64 v0 = parms[0];
    f64 v0_min = second_order_min_parms[0];
    f64 v0_max = second_order_max_parms[0];

    f64 b = parms[1];
    f64 b_min = second_order_min_parms[1];
    f64 b_max = second_order_max_parms[1];

    f64 c = parms[2];
    f64 c_min = second_order_min_parms[2];
    f64 c_max = second_order_max_parms[2];

    f64 sigma = parms[3];
    f64 sigma_min = second_order_min_parms[3];
    f64 sigma_max = second_order_max_parms[3];
    
    if(
            sigma < sigma_min || sigma > sigma_max 
            || v0 < v0_min || v0 > v0_max
            || b < b_min || b > b_max
            || c < c_min || c > c_max
      )
    {
        return -INFINITY;
    }

    f64 zc = (v0_max - v0_min) * (b_max - b_min) * (c_max - c_min) * log(sigma_max / sigma_min);
    return -log(zc * sigma);
}

static f64 
second_order_model_log_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    f64 v0 = parms[0];
    f64 b = parms[1];
    f64 c = parms[2];
    f64 sigma = parms[3];

    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 y_sq_bar = data->y_sq_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x4_bar = data->x4_bar;
    f64 y_bar = data->y_bar;
    f64 x_bar = data->x_bar;
    f64 xy_bar = data->xy_bar;
    f64 x3_bar = data->x3_bar;
    f64 y_x2_bar = data->y_x2_bar;

    f64 Q = v0 * v0 + b * b * x_sq_bar + c * c * x4_bar + y_sq_bar +
            2.0 * (v0 * b * x_bar + v0 * c * x_sq_bar + b * c * x3_bar - b * xy_bar - v0 * y_bar - c * y_x2_bar);

    data->ncalls++;

    return -N * 0.5 * log(ELK_2_PI) - N * log(sigma) - 0.5 * N * Q / (sigma * sigma);
}

static void
second_order_model_max_a_posteriori(void *user_data, size num_parms, f64 *parms_out)
{
    Assert(num_parms == 4);

    byte buffer2[sizeof(i32) * 3 * 3 + sizeof(f64) * 9 + _Alignof(f64)] = {0};
    MagStaticArena scratch_ = mag_static_arena_create(sizeof(buffer2) / sizeof(buffer2[0]), buffer2);
    MagAllocator scratch = mag_allocator_from_static_arena(&scratch_);

    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 y_sq_bar = data->y_sq_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x4_bar = data->x4_bar;
    f64 y_bar = data->y_bar;
    f64 x_bar = data->x_bar;
    f64 xy_bar = data->xy_bar;
    f64 x3_bar = data->x3_bar;
    f64 y_x2_bar = data->y_x2_bar;

    f64 m_data[9] = 
        {
                 1.0,    x_bar, x_sq_bar,
               x_bar, x_sq_bar,   x3_bar,
            x_sq_bar,   x3_bar,   x4_bar
        };
    BayLaSquareMatrix m = { .ndim = 3, .data = m_data };

    f64 input[3] = { y_bar, xy_bar, y_x2_bar };

    bayla_square_matrix_solve(m, input, scratch);
    memcpy(parms_out, input, sizeof(f64) * 3);

    f64 v0 = parms_out[0];
    f64 b = parms_out[1];
    f64 c = parms_out[2];

    f64 Q = v0 * v0 + b * b * x_sq_bar + c * c * x4_bar + y_sq_bar +
            2.0 * (v0 * b * x_bar + v0 * c * x_sq_bar + b * c * x3_bar - b * xy_bar - v0 * y_bar - c * y_x2_bar);

    parms_out[3] = sqrt(N * Q / (N + 1.0));

    for(size i = 0; i < num_parms; ++i)
    {
        Assert(parms_out[i] >= second_order_min_parms[i] && parms_out[i] <= second_order_max_parms[i]);
    }

#if 0
    printf("Second Order Model max-a-postiori: %e + %e * x + %e x**2 (σ = %e)\n",
            parms_out[0],  parms_out[1], parms_out[2], parms_out[3]);
#endif
}

static BayLaSquareMatrix
second_order_model_2d_hessian(void *user_data, size num_parms, f64 const *max_a_posteriori_parms, MagAllocator *alloc)
{
    Assert(num_parms == 4);

    f64 v0 = max_a_posteriori_parms[0];
    f64 b = max_a_posteriori_parms[1];
    f64 c = max_a_posteriori_parms[2];
    f64 sigma = max_a_posteriori_parms[3];

    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 y_sq_bar = data->y_sq_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x4_bar = data->x4_bar;
    f64 y_bar = data->y_bar;
    f64 x_bar = data->x_bar;
    f64 xy_bar = data->xy_bar;
    f64 x3_bar = data->x3_bar;
    f64 y_x2_bar = data->y_x2_bar;

    f64 Q = v0 * v0 + b * b * x_sq_bar + c * c * x4_bar + y_sq_bar +
            2.0 * (v0 * b * x_bar + v0 * c * x_sq_bar + b * c * x3_bar - b * xy_bar - v0 * y_bar - c * y_x2_bar);

    f64 prior = second_order_model_log_prior(num_parms, max_a_posteriori_parms, user_data);
    f64 likelihood = second_order_model_log_likelihood(num_parms, max_a_posteriori_parms, user_data);
    f64 p = exp(prior + likelihood);
    f64 sigma2 = sigma * sigma;
    f64 sigma3 = sigma2 * sigma;
    f64 sigma4 = sigma3 * sigma;

    BayLaSquareMatrix hess = bayla_square_matrix_create(4, alloc);

    hess.data[MATRIX_IDX(0, 0, 4)] = -N / sigma2;
    hess.data[MATRIX_IDX(0, 1, 4)] = hess.data[MATRIX_IDX(1, 0, 4)] = -N * x_bar / sigma2;
    hess.data[MATRIX_IDX(0, 2, 4)] = hess.data[MATRIX_IDX(2, 0, 4)] = -N * x_sq_bar / sigma2;
    hess.data[MATRIX_IDX(0, 3, 4)] = hess.data[MATRIX_IDX(3, 0, 4)] = 2 * N / sigma3 * (v0 + b * x_bar + c * x_sq_bar - y_bar);

    hess.data[MATRIX_IDX(1, 1, 4)] = -N / sigma2 * x_sq_bar;
    hess.data[MATRIX_IDX(1, 2, 4)] = hess.data[MATRIX_IDX(2, 1, 4)] = -N * x3_bar / sigma2;
    hess.data[MATRIX_IDX(1, 3, 4)] = hess.data[MATRIX_IDX(3, 1, 4)] = 2 * N / sigma3 * (v0 * x_bar + b * x_sq_bar + c * x3_bar - xy_bar);

    hess.data[MATRIX_IDX(2, 2, 4)] = -N / sigma2 * x4_bar;
    hess.data[MATRIX_IDX(2, 3, 4)] = hess.data[MATRIX_IDX(3, 2, 4)] = 2 * N / sigma3 * (v0 * x_sq_bar + b * x3_bar + c * x4_bar - y_x2_bar);

    hess.data[MATRIX_IDX(3, 3, 4)] = (N + 1.0) / sigma2 - 3.0 * N * Q / sigma4;

    for(size i = 0; i < 16; ++i) { hess.data[i] *= p; }

    Assert(hess.data[MATRIX_IDX(0, 0, 4)] < 0.0 && hess.data[MATRIX_IDX(1, 1, 4)] < 0.0 && hess.data[MATRIX_IDX(2, 2, 4)] < 0.0 && hess.data[MATRIX_IDX(3, 3, 4)] < 0.0);

    return hess;
}

UserData second_order_user_data = {0};

BayLaModel second_order_model = 
    {
        .model_prior_probability = 1.0 / 7.0,
        .n_parameters = 4,
        .prior = second_order_model_log_prior,
        .likelihood = second_order_model_log_likelihood,
        .posterior_max = second_order_model_max_a_posteriori,
        .hessian_matrix = second_order_model_2d_hessian,
        .user_data = &second_order_user_data
    };

/*-----------------------------------------------  3rd Order Polynomial  -------------------------------------------------*/
/* param 0 is constant, param 1 is linear, param 2 is quadratic , parm 3 is cubic, param 4 is a standard deviation. */
f64 third_order_min_parms[5] = {-4.0, -4.0, -4.0, -4.0, log_model_min_sigma };
f64 third_order_max_parms[5] = {+4.0, +4.0, +4.0, +4.0, log_model_max_sigma };

static f64 
third_order_model_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    f64 v0 = parms[0];
    f64 v0_min = third_order_min_parms[0];
    f64 v0_max = third_order_max_parms[0];

    f64 b = parms[1];
    f64 b_min = third_order_min_parms[1];
    f64 b_max = third_order_max_parms[1];

    f64 c = parms[2];
    f64 c_min = third_order_min_parms[2];
    f64 c_max = third_order_max_parms[2];

    f64 d = parms[3];
    f64 d_min = third_order_min_parms[3];
    f64 d_max = third_order_max_parms[3];
    
    f64 sigma = parms[4];
    f64 sigma_min = third_order_min_parms[4];
    f64 sigma_max = third_order_max_parms[4];
    
    if(
            sigma < sigma_min || sigma > sigma_max 
            || v0 < v0_min || v0 > v0_max
            || b < b_min || b > b_max
            || c < c_min || c > c_max
            || d < d_min || d > d_max
      )
    {
        return -INFINITY;
    }

    f64 zc = (v0_max - v0_min) * (b_max - b_min) * (c_max - c_min) * (d_max - d_min) * log(sigma_max / sigma_min);
    return -log(zc * sigma);
}

static f64 
third_order_model_log_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    f64 v0 = parms[0];
    f64 b = parms[1];
    f64 c = parms[2];
    f64 d = parms[3];
    f64 sigma = parms[4];

    UserData *data = user_data;
    f64 N = (f64)data->N;

    f64 x_bar = data->x_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x3_bar = data->x3_bar;
    f64 x4_bar = data->x4_bar;
    f64 x5_bar = data->x5_bar;
    f64 x6_bar = data->x6_bar;

    f64 y_bar = data->y_bar;
    f64 y_sq_bar = data->y_sq_bar;

    f64 xy_bar = data->xy_bar;
    f64 y_x2_bar = data->y_x2_bar;
    f64 y_x3_bar = data->y_x3_bar;

    f64 Q = v0 * v0 + b * b * x_sq_bar + c * c * x4_bar + d * d *x6_bar + y_sq_bar +
            2.0 * (v0 * b * x_bar + v0 * c * x_sq_bar + v0 * d * x3_bar + b * c * x3_bar + b * d * x4_bar + c * d * x5_bar
                    - v0 * y_bar - b * xy_bar - c * y_x2_bar - d * y_x3_bar);

    data->ncalls++;

    return -N * 0.5 * log(ELK_2_PI) - N * log(sigma) - 0.5 * N * Q / (sigma * sigma);
}

static void
third_order_model_max_a_posteriori(void *user_data, size num_parms, f64 *parms_out)
{
    Assert(num_parms == 5);

    byte buffer2[sizeof(i32) * 3 * 4 + sizeof(f64) * 16 + _Alignof(f64)] = {0};
    MagStaticArena scratch_ = mag_static_arena_create(sizeof(buffer2) / sizeof(buffer2[0]), buffer2);
    MagAllocator scratch = mag_allocator_from_static_arena(&scratch_);

    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 x_bar = data->x_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x3_bar = data->x3_bar;
    f64 x4_bar = data->x4_bar;
    f64 x5_bar = data->x5_bar;
    f64 x6_bar = data->x6_bar;

    f64 y_bar = data->y_bar;
    f64 y_sq_bar = data->y_sq_bar;

    f64 xy_bar = data->xy_bar;
    f64 y_x2_bar = data->y_x2_bar;
    f64 y_x3_bar = data->y_x3_bar;

    f64 m_data[16] = 
        {
                 1.0,    x_bar, x_sq_bar, x3_bar,
               x_bar, x_sq_bar,   x3_bar, x4_bar,
            x_sq_bar,   x3_bar,   x4_bar, x5_bar,
              x3_bar,   x4_bar,   x5_bar, x6_bar,
        };
    BayLaSquareMatrix m = { .ndim = 4, .data = m_data };

    f64 input[4] = { y_bar, xy_bar, y_x2_bar, y_x3_bar };

    bayla_square_matrix_solve(m, input, scratch);
    memcpy(parms_out, input, sizeof(f64) * 4);

    f64 v0 = parms_out[0];
    f64 b = parms_out[1];
    f64 c = parms_out[2];
    f64 d = parms_out[3];

    f64 Q = v0 * v0 + b * b * x_sq_bar + c * c * x4_bar + d * d *x6_bar + y_sq_bar +
            2.0 * (v0 * b * x_bar + v0 * c * x_sq_bar + v0 * d * x3_bar + b * c * x3_bar + b * d * x4_bar + c * d * x5_bar
                    - v0 * y_bar - b * xy_bar - c * y_x2_bar - d * y_x3_bar);

    parms_out[4] = sqrt(N * Q / (N + 1.0));

    for(size i = 0; i < num_parms; ++i)
    {
        Assert(parms_out[i] >= third_order_min_parms[i] && parms_out[i] <= third_order_max_parms[i]);
    }

#if 0
    printf("Third Order Model max-a-postiori: %e + %e * x + %e x**2  + %e x**3 (σ = %e)\n",
            parms_out[0],  parms_out[1], parms_out[2], parms_out[3], parms_out[4]);
#endif
}

static BayLaSquareMatrix
third_order_model_2d_hessian(void *user_data, size num_parms, f64 const *max_a_posteriori_parms, MagAllocator *alloc)
{
    Assert(num_parms == 5);

    f64 v0 = max_a_posteriori_parms[0];
    f64 b = max_a_posteriori_parms[1];
    f64 c = max_a_posteriori_parms[2];
    f64 d = max_a_posteriori_parms[3];
    f64 sigma = max_a_posteriori_parms[4];

    UserData *data = user_data;
    f64 N = (f64)data->N;

    f64 x_bar = data->x_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x3_bar = data->x3_bar;
    f64 x4_bar = data->x4_bar;
    f64 x5_bar = data->x5_bar;
    f64 x6_bar = data->x6_bar;

    f64 y_bar = data->y_bar;
    f64 y_sq_bar = data->y_sq_bar;

    f64 xy_bar = data->xy_bar;
    f64 y_x2_bar = data->y_x2_bar;
    f64 y_x3_bar = data->y_x3_bar;

    f64 Q = v0 * v0 + b * b * x_sq_bar + c * c * x4_bar + d * d *x6_bar + y_sq_bar +
            2.0 * (v0 * b * x_bar + v0 * c * x_sq_bar + v0 * d * x3_bar + b * c * x3_bar + b * d * x4_bar + c * d * x5_bar
                    - v0 * y_bar - b * xy_bar - c * y_x2_bar - d * y_x3_bar);


    f64 prior = third_order_model_log_prior(num_parms, max_a_posteriori_parms, user_data);
    f64 likelihood = third_order_model_log_likelihood(num_parms, max_a_posteriori_parms, user_data);
    f64 p = exp(prior + likelihood);
    f64 sigma2 = sigma * sigma;
    f64 sigma3 = sigma2 * sigma;
    f64 sigma4 = sigma3 * sigma;

    BayLaSquareMatrix hess = bayla_square_matrix_create(5, alloc);

    hess.data[MATRIX_IDX(0, 0, 5)] = -N / sigma2;
    hess.data[MATRIX_IDX(0, 1, 5)] = hess.data[MATRIX_IDX(1, 0, 5)] = -N * x_bar / sigma2;
    hess.data[MATRIX_IDX(0, 2, 5)] = hess.data[MATRIX_IDX(2, 0, 5)] = -N * x_sq_bar / sigma2;
    hess.data[MATRIX_IDX(0, 3, 5)] = hess.data[MATRIX_IDX(3, 0, 5)] = -N * x3_bar / sigma2;
    hess.data[MATRIX_IDX(0, 4, 5)] = hess.data[MATRIX_IDX(4, 0, 5)] = 2 * N / sigma3 * (v0 + b * x_bar + c * x_sq_bar + d * x3_bar - y_bar);

    hess.data[MATRIX_IDX(1, 1, 5)] = -N / sigma2 * x_sq_bar;
    hess.data[MATRIX_IDX(1, 2, 5)] = hess.data[MATRIX_IDX(2, 1, 5)] = -N * x3_bar / sigma2;
    hess.data[MATRIX_IDX(1, 3, 5)] = hess.data[MATRIX_IDX(3, 1, 5)] = -N * x4_bar / sigma2;
    hess.data[MATRIX_IDX(1, 4, 5)] = hess.data[MATRIX_IDX(4, 1, 5)] = 2 * N / sigma3 * (v0 * x_bar + b * x_sq_bar + c * x3_bar + d * x4_bar - xy_bar);

    hess.data[MATRIX_IDX(2, 2, 5)] = -N / sigma2 * x4_bar;
    hess.data[MATRIX_IDX(2, 3, 5)] = hess.data[MATRIX_IDX(3, 2, 5)] = -N * x5_bar / sigma2;
    hess.data[MATRIX_IDX(2, 4, 5)] = hess.data[MATRIX_IDX(4, 2, 5)] = 2 * N / sigma3 * (v0 * x_sq_bar + b * x3_bar + c * x4_bar + d * x5_bar - y_x2_bar);

    hess.data[MATRIX_IDX(3, 3, 5)] = -N / sigma2 * x6_bar;
    hess.data[MATRIX_IDX(3, 4, 5)] = hess.data[MATRIX_IDX(4, 3, 5)] = 2 * N / sigma3 * (v0 * x3_bar + b * x4_bar + c * x5_bar + d * x6_bar - y_x3_bar);

    hess.data[MATRIX_IDX(4, 4, 5)] = (N + 1.0) / sigma2 - 3.0 * N * Q / sigma4;

    for(size i = 0; i < 25; ++i) { hess.data[i] *= p; }

    Assert(hess.data[MATRIX_IDX(0, 0, 5)] < 0.0 && hess.data[MATRIX_IDX(1, 1, 5)] < 0.0 && hess.data[MATRIX_IDX(2, 2, 5)] < 0.0 && hess.data[MATRIX_IDX(3, 3, 5)] < 0.0 && hess.data[MATRIX_IDX(4, 4, 5)] < 0.0);

    return hess;
}

UserData third_order_user_data = {0};

BayLaModel third_order_model = 
    {
        .model_prior_probability = 1.0 / 7.0,
        .n_parameters = 5,
        .prior = third_order_model_log_prior,
        .likelihood = third_order_model_log_likelihood,
        .posterior_max = third_order_model_max_a_posteriori,
        .hessian_matrix = third_order_model_2d_hessian,
        .user_data = &third_order_user_data
    };

/*-----------------------------------------------  4th Order Polynomial  -------------------------------------------------*/
/* param 0 is constant, param 1 is linear, param 2 is quadratic , parm 3 is cubic, parm 4 is quartic, param 5 is a
 * standard deviation. */
f64 fourth_order_min_parms[6] = {-4.0, -4.0, -4.0, -4.0, -4.0, log_model_min_sigma };
f64 fourth_order_max_parms[6] = {+4.0, +4.0, +4.0, +4.0, +4.0, log_model_max_sigma };

static f64 
fourth_order_model_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    f64 v0 = parms[0];
    f64 v0_min = fourth_order_min_parms[0];
    f64 v0_max = fourth_order_max_parms[0];

    f64 b = parms[1];
    f64 b_min = fourth_order_min_parms[1];
    f64 b_max = fourth_order_max_parms[1];

    f64 c = parms[2];
    f64 c_min = fourth_order_min_parms[2];
    f64 c_max = fourth_order_max_parms[2];

    f64 d = parms[3];
    f64 d_min = fourth_order_min_parms[3];
    f64 d_max = fourth_order_max_parms[3];
    
    f64 e = parms[4];
    f64 e_min = fourth_order_min_parms[4];
    f64 e_max = fourth_order_max_parms[4];
    
    f64 sigma = parms[5];
    f64 sigma_min = fourth_order_min_parms[5];
    f64 sigma_max = fourth_order_max_parms[5];
    
    if(
            sigma < sigma_min || sigma > sigma_max 
            || v0 < v0_min || v0 > v0_max
            || b < b_min || b > b_max
            || c < c_min || c > c_max
            || d < d_min || d > d_max
            || e < e_min || e > e_max
      )
    {
        return -INFINITY;
    }

    f64 zc = (v0_max - v0_min) * (b_max - b_min) * (c_max - c_min) * (d_max - d_min) * (e_max - e_min) * log(sigma_max / sigma_min);
    return -log(zc * sigma);
}

static f64 
fourth_order_model_log_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    f64 v0 = parms[0];
    f64 b = parms[1];
    f64 c = parms[2];
    f64 d = parms[3];
    f64 e = parms[4];
    f64 sigma = parms[5];

    UserData *data = user_data;
    f64 N = (f64)data->N;

    f64 x_bar = data->x_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x3_bar = data->x3_bar;
    f64 x4_bar = data->x4_bar;
    f64 x5_bar = data->x5_bar;
    f64 x6_bar = data->x6_bar;
    f64 x7_bar = data->x7_bar;
    f64 x8_bar = data->x8_bar;

    f64 y_bar = data->y_bar;
    f64 y_sq_bar = data->y_sq_bar;

    f64 xy_bar = data->xy_bar;
    f64 y_x2_bar = data->y_x2_bar;
    f64 y_x3_bar = data->y_x3_bar;
    f64 y_x4_bar = data->y_x4_bar;

    f64 Q = v0 * v0 + b * b * x_sq_bar + c * c * x4_bar + d * d * x6_bar + e * e * x8_bar + y_sq_bar +
            2.0 * (v0 * b * x_bar + v0 * c * x_sq_bar + v0 * d * x3_bar + v0 * e * x4_bar
                   + b * c * x3_bar + b * d * x4_bar + b * e * x5_bar
                   + c * d * x5_bar + c * e * x6_bar 
                   + d * e * x7_bar
                   - v0 * y_bar - b * xy_bar - c * y_x2_bar - d * y_x3_bar - e * y_x4_bar);

    data->ncalls++;

    return -N * 0.5 * log(ELK_2_PI) - N * log(sigma) - 0.5 * N * Q / (sigma * sigma);
}

static void
fourth_order_model_max_a_posteriori(void *user_data, size num_parms, f64 *parms_out)
{
    Assert(num_parms == 6);

    byte buffer2[sizeof(i32) * 3 * 5 + sizeof(f64) * 25 + _Alignof(f64)] = {0};
    MagStaticArena scratch_ = mag_static_arena_create(sizeof(buffer2) / sizeof(buffer2[0]), buffer2);
    MagAllocator scratch = mag_allocator_from_static_arena(&scratch_);

    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 x_bar = data->x_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x3_bar = data->x3_bar;
    f64 x4_bar = data->x4_bar;
    f64 x5_bar = data->x5_bar;
    f64 x6_bar = data->x6_bar;
    f64 x7_bar = data->x7_bar;
    f64 x8_bar = data->x8_bar;

    f64 y_bar = data->y_bar;
    f64 y_sq_bar = data->y_sq_bar;

    f64 xy_bar = data->xy_bar;
    f64 y_x2_bar = data->y_x2_bar;
    f64 y_x3_bar = data->y_x3_bar;
    f64 y_x4_bar = data->y_x4_bar;

    f64 m_data[25] = 
        {
                 1.0,    x_bar, x_sq_bar, x3_bar, x4_bar,
               x_bar, x_sq_bar,   x3_bar, x4_bar, x5_bar,
            x_sq_bar,   x3_bar,   x4_bar, x5_bar, x5_bar,
              x3_bar,   x4_bar,   x5_bar, x6_bar, x5_bar,
              x4_bar,   x5_bar,   x6_bar, x7_bar, x8_bar
        };
    BayLaSquareMatrix m = { .ndim = 5, .data = m_data };

    f64 input[5] = { y_bar, xy_bar, y_x2_bar, y_x3_bar, y_x4_bar };

    bayla_square_matrix_solve(m, input, scratch);
    memcpy(parms_out, input, sizeof(f64) * 5);

    f64 v0 = parms_out[0];
    f64 b = parms_out[1];
    f64 c = parms_out[2];
    f64 d = parms_out[3];
    f64 e = parms_out[4];

    f64 Q = v0 * v0 + b * b * x_sq_bar + c * c * x4_bar + d * d * x6_bar + e * e * x8_bar + y_sq_bar +
            2.0 * (v0 * b * x_bar + v0 * c * x_sq_bar + v0 * d * x3_bar + v0 * e * x4_bar
                   + b * c * x3_bar + b * d * x4_bar + b * e * x5_bar
                   + c * d * x5_bar + c * e * x6_bar 
                   + d * e * x7_bar
                   - v0 * y_bar - b * xy_bar - c * y_x2_bar - d * y_x3_bar - e * y_x4_bar);

    parms_out[5] = sqrt(N * Q / (N + 1.0));

    for(size i = 0; i < num_parms; ++i)
    {
        Assert(parms_out[i] >= fourth_order_min_parms[i] && parms_out[i] <= fourth_order_max_parms[i]);
    }

#if 0
    printf("Fourth Order Model max-a-postiori: %e + %e * x + %e x**2  + %e x**3 + %e x**4 (σ = %e)\n",
            parms_out[0],  parms_out[1], parms_out[2], parms_out[3], parms_out[4], parms_out[5]);
#endif
}

static BayLaSquareMatrix
fourth_order_model_2d_hessian(void *user_data, size num_parms, f64 const *max_a_posteriori_parms, MagAllocator *alloc)
{
    Assert(num_parms == 6);

    f64 v0 = max_a_posteriori_parms[0];
    f64 b = max_a_posteriori_parms[1];
    f64 c = max_a_posteriori_parms[2];
    f64 d = max_a_posteriori_parms[3];
    f64 e = max_a_posteriori_parms[4];
    f64 sigma = max_a_posteriori_parms[5];

    UserData *data = user_data;
    f64 N = (f64)data->N;

    f64 x_bar = data->x_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x3_bar = data->x3_bar;
    f64 x4_bar = data->x4_bar;
    f64 x5_bar = data->x5_bar;
    f64 x6_bar = data->x6_bar;
    f64 x7_bar = data->x7_bar;
    f64 x8_bar = data->x8_bar;

    f64 y_bar = data->y_bar;
    f64 y_sq_bar = data->y_sq_bar;

    f64 xy_bar = data->xy_bar;
    f64 y_x2_bar = data->y_x2_bar;
    f64 y_x3_bar = data->y_x3_bar;
    f64 y_x4_bar = data->y_x4_bar;

    f64 Q = v0 * v0 + b * b * x_sq_bar + c * c * x4_bar + d * d * x6_bar + e * e * x8_bar + y_sq_bar +
            2.0 * (v0 * b * x_bar + v0 * c * x_sq_bar + v0 * d * x3_bar + v0 * e * x4_bar
                   + b * c * x3_bar + b * d * x4_bar + b * e * x5_bar
                   + c * d * x5_bar + c * e * x6_bar 
                   + d * e * x7_bar
                   - v0 * y_bar - b * xy_bar - c * y_x2_bar - d * y_x3_bar - e * y_x4_bar);

    f64 prior = fourth_order_model_log_prior(num_parms, max_a_posteriori_parms, user_data);
    f64 likelihood = fourth_order_model_log_likelihood(num_parms, max_a_posteriori_parms, user_data);
    f64 p = exp(prior + likelihood);
    f64 sigma2 = sigma * sigma;
    f64 sigma3 = sigma2 * sigma;
    f64 sigma4 = sigma3 * sigma;

    BayLaSquareMatrix hess = bayla_square_matrix_create(6, alloc);

    hess.data[MATRIX_IDX(0, 0, 6)] = -N / sigma2;
    hess.data[MATRIX_IDX(0, 1, 6)] = hess.data[MATRIX_IDX(1, 0, 6)] = -N * x_bar / sigma2;
    hess.data[MATRIX_IDX(0, 2, 6)] = hess.data[MATRIX_IDX(2, 0, 6)] = -N * x_sq_bar / sigma2;
    hess.data[MATRIX_IDX(0, 3, 6)] = hess.data[MATRIX_IDX(3, 0, 6)] = -N * x3_bar / sigma2;
    hess.data[MATRIX_IDX(0, 4, 6)] = hess.data[MATRIX_IDX(4, 0, 6)] = -N * x4_bar / sigma2;
    hess.data[MATRIX_IDX(0, 5, 6)] = hess.data[MATRIX_IDX(5, 0, 6)] = 2 * N / sigma3 * (v0 + b * x_bar + c * x_sq_bar + d * x3_bar + e * x4_bar - y_bar);

    hess.data[MATRIX_IDX(1, 1, 6)] = -N / sigma2 * x_sq_bar;
    hess.data[MATRIX_IDX(1, 2, 6)] = hess.data[MATRIX_IDX(2, 1, 6)] = -N * x3_bar / sigma2;
    hess.data[MATRIX_IDX(1, 3, 6)] = hess.data[MATRIX_IDX(3, 1, 6)] = -N * x4_bar / sigma2;
    hess.data[MATRIX_IDX(1, 4, 6)] = hess.data[MATRIX_IDX(4, 1, 6)] = -N * x5_bar / sigma2;
    hess.data[MATRIX_IDX(1, 5, 6)] = hess.data[MATRIX_IDX(5, 1, 6)] = 2 * N / sigma3 * (v0 * x_bar + b * x_sq_bar + c * x3_bar + d * x4_bar + e * x5_bar - xy_bar);

    hess.data[MATRIX_IDX(2, 2, 6)] = -N / sigma2 * x4_bar;
    hess.data[MATRIX_IDX(2, 3, 6)] = hess.data[MATRIX_IDX(3, 2, 6)] = -N * x5_bar / sigma2;
    hess.data[MATRIX_IDX(2, 4, 6)] = hess.data[MATRIX_IDX(4, 2, 6)] = -N * x6_bar / sigma2;
    hess.data[MATRIX_IDX(2, 5, 6)] = hess.data[MATRIX_IDX(5, 2, 6)] = 2 * N / sigma3 * (v0 * x_sq_bar + b * x3_bar + c * x4_bar + d * x5_bar + e * x6_bar - y_x2_bar);

    hess.data[MATRIX_IDX(3, 3, 6)] = -N / sigma2 * x6_bar;
    hess.data[MATRIX_IDX(3, 4, 6)] = hess.data[MATRIX_IDX(4, 3, 6)] = -N / sigma2 * x7_bar;
    hess.data[MATRIX_IDX(3, 5, 6)] = hess.data[MATRIX_IDX(5, 3, 6)] = 2 * N / sigma3 * (v0 * x3_bar + b * x4_bar + c * x5_bar + d * x6_bar + e * x7_bar - y_x3_bar);

    hess.data[MATRIX_IDX(4, 4, 6)] = -N / sigma2 * x8_bar;
    hess.data[MATRIX_IDX(4, 5, 6)] = hess.data[MATRIX_IDX(5, 4, 6)] =  2 * N / sigma3 * (v0 * x4_bar + b * x5_bar + c * x6_bar + d * x7_bar + e * x8_bar - y_x4_bar);

    hess.data[MATRIX_IDX(5, 5, 6)] = (N + 1.0) / sigma2 - 3.0 * N * Q / sigma4;

    for(size i = 0; i < num_parms * num_parms; ++i) { hess.data[i] *= p; }

    for(size i = 0; i < num_parms; ++i) { Assert(hess.data[MATRIX_IDX(i, i, 6)] < 0.0); }

    return hess;
}

UserData fourth_order_user_data = {0};

BayLaModel fourth_order_model = 
    {
        .model_prior_probability = 1.0 / 7.0,
        .n_parameters = 6,
        .prior = fourth_order_model_log_prior,
        .likelihood = fourth_order_model_log_likelihood,
        .posterior_max = fourth_order_model_max_a_posteriori,
        .hessian_matrix = fourth_order_model_2d_hessian,
        .user_data = &fourth_order_user_data
    };

/*-----------------------------------------------  5th Order Polynomial  -------------------------------------------------*/
/* param 0 is constant, param 1 is linear, param 2 is quadratic , parm 3 is cubic, parm 4 is quartic, parm 5 is quintic,
 * param 6 is a standard deviation. */
f64 fifth_order_min_parms[7] = {-4.0, -4.0, -4.0, -4.0, -4.0, -4.0, log_model_min_sigma };
f64 fifth_order_max_parms[7] = {+4.0, +4.0, +4.0, +4.0, +4.0, +4.0, log_model_max_sigma };

static f64 
fifth_order_model_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    f64 v0 = parms[0];
    f64 v0_min = fifth_order_min_parms[0];
    f64 v0_max = fifth_order_max_parms[0];

    f64 b = parms[1];
    f64 b_min = fifth_order_min_parms[1];
    f64 b_max = fifth_order_max_parms[1];

    f64 c = parms[2];
    f64 c_min = fifth_order_min_parms[2];
    f64 c_max = fifth_order_max_parms[2];

    f64 d = parms[3];
    f64 d_min = fifth_order_min_parms[3];
    f64 d_max = fifth_order_max_parms[3];
    
    f64 e = parms[4];
    f64 e_min = fifth_order_min_parms[4];
    f64 e_max = fifth_order_max_parms[4];
    
    f64 f = parms[5];
    f64 f_min = fifth_order_min_parms[5];
    f64 f_max = fifth_order_max_parms[5];
    
    f64 sigma = parms[6];
    f64 sigma_min = fifth_order_min_parms[6];
    f64 sigma_max = fifth_order_max_parms[6];
    
    if(
            sigma < sigma_min || sigma > sigma_max 
            || v0 < v0_min || v0 > v0_max
            || b < b_min || b > b_max
            || c < c_min || c > c_max
            || d < d_min || d > d_max
            || e < e_min || e > e_max
            || f < f_min || f > f_max
      )
    {
        return -INFINITY;
    }

    f64 zc = (v0_max - v0_min) * (b_max - b_min) * (c_max - c_min) * (d_max - d_min) * (e_max - e_min) * (f_max - f_min) * log(sigma_max / sigma_min);
    return -log(zc * sigma);
}

static f64 
fifth_order_model_log_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    f64 v0 = parms[0];
    f64 b = parms[1];
    f64 c = parms[2];
    f64 d = parms[3];
    f64 e = parms[4];
    f64 f = parms[5];
    f64 sigma = parms[6];

    UserData *data = user_data;
    f64 N = (f64)data->N;

    f64 x_bar = data->x_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x3_bar = data->x3_bar;
    f64 x4_bar = data->x4_bar;
    f64 x5_bar = data->x5_bar;
    f64 x6_bar = data->x6_bar;
    f64 x7_bar = data->x7_bar;
    f64 x8_bar = data->x8_bar;
    f64 x9_bar = data->x9_bar;
    f64 x10_bar = data->x10_bar;

    f64 y_bar = data->y_bar;
    f64 y_sq_bar = data->y_sq_bar;

    f64 xy_bar = data->xy_bar;
    f64 y_x2_bar = data->y_x2_bar;
    f64 y_x3_bar = data->y_x3_bar;
    f64 y_x4_bar = data->y_x4_bar;
    f64 y_x5_bar = data->y_x5_bar;

    f64 Q = v0 * v0 + b * b * x_sq_bar + c * c * x4_bar + d * d * x6_bar + e * e * x8_bar + f * f * x10_bar + y_sq_bar +
            2.0 * (v0 * b * x_bar + v0 * c * x_sq_bar + v0 * d * x3_bar + v0 * e * x4_bar + f * f * x5_bar
                   + b * c * x3_bar + b * d * x4_bar + b * e * x5_bar + b * f * x6_bar
                   + c * d * x5_bar + c * e * x6_bar + c * f * x7_bar
                   + d * e * x7_bar + d * f * x8_bar
                   + e * f * x9_bar
                   - v0 * y_bar - b * xy_bar - c * y_x2_bar - d * y_x3_bar - e * y_x4_bar - f * y_x5_bar);

    data->ncalls++;

    return -N * 0.5 * log(ELK_2_PI) - N * log(sigma) - 0.5 * N * Q / (sigma * sigma);
}

static void
fifth_order_model_max_a_posteriori(void *user_data, size num_parms, f64 *parms_out)
{
    Assert(num_parms == 7);

    byte buffer2[sizeof(i32) * 3 * 6 + sizeof(f64) * 36 + _Alignof(f64)] = {0};
    MagStaticArena scratch_ = mag_static_arena_create(sizeof(buffer2) / sizeof(buffer2[0]), buffer2);
    MagAllocator scratch = mag_allocator_from_static_arena(&scratch_);

    UserData *data = user_data;
    f64 N = (f64)data->N;
    f64 x_bar = data->x_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x3_bar = data->x3_bar;
    f64 x4_bar = data->x4_bar;
    f64 x5_bar = data->x5_bar;
    f64 x6_bar = data->x6_bar;
    f64 x7_bar = data->x7_bar;
    f64 x8_bar = data->x8_bar;
    f64 x9_bar = data->x9_bar;
    f64 x10_bar = data->x10_bar;

    f64 y_bar = data->y_bar;
    f64 y_sq_bar = data->y_sq_bar;

    f64 xy_bar = data->xy_bar;
    f64 y_x2_bar = data->y_x2_bar;
    f64 y_x3_bar = data->y_x3_bar;
    f64 y_x4_bar = data->y_x4_bar;
    f64 y_x5_bar = data->y_x5_bar;

    f64 m_data[36] = 
        {
                 1.0,    x_bar, x_sq_bar, x3_bar, x4_bar, x5_bar,
               x_bar, x_sq_bar,   x3_bar, x4_bar, x5_bar, x6_bar,
            x_sq_bar,   x3_bar,   x4_bar, x5_bar, x5_bar, x7_bar,
              x3_bar,   x4_bar,   x5_bar, x6_bar, x5_bar, x7_bar,
              x4_bar,   x5_bar,   x6_bar, x7_bar, x8_bar, x9_bar,
              x5_bar,   x6_bar,   x7_bar, x8_bar, x9_bar, x10_bar
        };
    BayLaSquareMatrix m = { .ndim = 6, .data = m_data };

    f64 input[6] = { y_bar, xy_bar, y_x2_bar, y_x3_bar, y_x4_bar, y_x5_bar };

    bayla_square_matrix_solve(m, input, scratch);
    memcpy(parms_out, input, sizeof(f64) * 6);

    f64 v0 = parms_out[0];
    f64 b = parms_out[1];
    f64 c = parms_out[2];
    f64 d = parms_out[3];
    f64 e = parms_out[4];
    f64 f = parms_out[5];

    f64 Q = v0 * v0 + b * b * x_sq_bar + c * c * x4_bar + d * d * x6_bar + e * e * x8_bar + f * f * x10_bar + y_sq_bar +
            2.0 * (v0 * b * x_bar + v0 * c * x_sq_bar + v0 * d * x3_bar + v0 * e * x4_bar + f * f * x5_bar
                   + b * c * x3_bar + b * d * x4_bar + b * e * x5_bar + b * f * x6_bar
                   + c * d * x5_bar + c * e * x6_bar + c * f * x7_bar
                   + d * e * x7_bar + d * f * x8_bar
                   + e * f * x9_bar
                   - v0 * y_bar - b * xy_bar - c * y_x2_bar - d * y_x3_bar - e * y_x4_bar - f * y_x5_bar);

    parms_out[6] = sqrt(N * Q / (N + 1.0));

    for(size i = 0; i < num_parms; ++i)
    {
        Assert(parms_out[i] >= fifth_order_min_parms[i] && parms_out[i] <= fifth_order_max_parms[i]);
    }

#if 0
    printf("Fifth Order Model max-a-postiori: %e + %e * x + %e x**2  + %e x**3 + %e x**4 + %e x**5 (σ = %e)\n",
            parms_out[0],  parms_out[1], parms_out[2], parms_out[3], parms_out[4], parms_out[5], parms_out[6]);
#endif
}

static BayLaSquareMatrix
fifth_order_model_2d_hessian(void *user_data, size num_parms, f64 const *max_a_posteriori_parms, MagAllocator *alloc)
{
    Assert(num_parms == 7);

    f64 v0 = max_a_posteriori_parms[0];
    f64 b = max_a_posteriori_parms[1];
    f64 c = max_a_posteriori_parms[2];
    f64 d = max_a_posteriori_parms[3];
    f64 e = max_a_posteriori_parms[4];
    f64 f = max_a_posteriori_parms[5];
    f64 sigma = max_a_posteriori_parms[6];

    UserData *data = user_data;
    f64 N = (f64)data->N;

    f64 x_bar = data->x_bar;
    f64 x_sq_bar = data->x_sq_bar;
    f64 x3_bar = data->x3_bar;
    f64 x4_bar = data->x4_bar;
    f64 x5_bar = data->x5_bar;
    f64 x6_bar = data->x6_bar;
    f64 x7_bar = data->x7_bar;
    f64 x8_bar = data->x8_bar;
    f64 x9_bar = data->x9_bar;
    f64 x10_bar = data->x10_bar;

    f64 y_bar = data->y_bar;
    f64 y_sq_bar = data->y_sq_bar;

    f64 xy_bar = data->xy_bar;
    f64 y_x2_bar = data->y_x2_bar;
    f64 y_x3_bar = data->y_x3_bar;
    f64 y_x4_bar = data->y_x4_bar;
    f64 y_x5_bar = data->y_x5_bar;

    f64 Q = v0 * v0 + b * b * x_sq_bar + c * c * x4_bar + d * d * x6_bar + e * e * x8_bar + f * f * x10_bar + y_sq_bar +
            2.0 * (v0 * b * x_bar + v0 * c * x_sq_bar + v0 * d * x3_bar + v0 * e * x4_bar + f * f * x5_bar
                   + b * c * x3_bar + b * d * x4_bar + b * e * x5_bar + b * f * x6_bar
                   + c * d * x5_bar + c * e * x6_bar + c * f * x7_bar
                   + d * e * x7_bar + d * f * x8_bar
                   + e * f * x9_bar
                   - v0 * y_bar - b * xy_bar - c * y_x2_bar - d * y_x3_bar - e * y_x4_bar - f * y_x5_bar);

    f64 prior = fifth_order_model_log_prior(num_parms, max_a_posteriori_parms, user_data);
    f64 likelihood = fifth_order_model_log_likelihood(num_parms, max_a_posteriori_parms, user_data);
    f64 p = exp(prior + likelihood);
    f64 sigma2 = sigma * sigma;
    f64 sigma3 = sigma2 * sigma;
    f64 sigma4 = sigma3 * sigma;

    BayLaSquareMatrix hess = bayla_square_matrix_create(7, alloc);

    hess.data[MATRIX_IDX(0, 0, 7)] = -N / sigma2;
    hess.data[MATRIX_IDX(0, 1, 7)] = hess.data[MATRIX_IDX(1, 0, 7)] = -N * x_bar / sigma2;
    hess.data[MATRIX_IDX(0, 2, 7)] = hess.data[MATRIX_IDX(2, 0, 7)] = -N * x_sq_bar / sigma2;
    hess.data[MATRIX_IDX(0, 3, 7)] = hess.data[MATRIX_IDX(3, 0, 7)] = -N * x3_bar / sigma2;
    hess.data[MATRIX_IDX(0, 4, 7)] = hess.data[MATRIX_IDX(4, 0, 7)] = -N * x4_bar / sigma2;
    hess.data[MATRIX_IDX(0, 5, 7)] = hess.data[MATRIX_IDX(5, 0, 7)] = -N * x5_bar / sigma2;
    hess.data[MATRIX_IDX(0, 6, 7)] = hess.data[MATRIX_IDX(6, 0, 7)] = 2 * N / sigma3 * (v0 + b * x_bar + c * x_sq_bar + d * x3_bar + e * x4_bar + f * x5_bar - y_bar);

    hess.data[MATRIX_IDX(1, 1, 7)] = -N / sigma2 * x_sq_bar;
    hess.data[MATRIX_IDX(1, 2, 7)] = hess.data[MATRIX_IDX(2, 1, 7)] = -N * x3_bar / sigma2;
    hess.data[MATRIX_IDX(1, 3, 7)] = hess.data[MATRIX_IDX(3, 1, 7)] = -N * x4_bar / sigma2;
    hess.data[MATRIX_IDX(1, 4, 7)] = hess.data[MATRIX_IDX(4, 1, 7)] = -N * x5_bar / sigma2;
    hess.data[MATRIX_IDX(1, 5, 7)] = hess.data[MATRIX_IDX(5, 1, 7)] = -N * x6_bar / sigma2;
    hess.data[MATRIX_IDX(1, 6, 7)] = hess.data[MATRIX_IDX(6, 1, 7)] = 2 * N / sigma3 * (v0 * x_bar + b * x_sq_bar + c * x3_bar + d * x4_bar + e * x5_bar + f * x6_bar - xy_bar);

    hess.data[MATRIX_IDX(2, 2, 7)] = -N / sigma2 * x4_bar;
    hess.data[MATRIX_IDX(2, 3, 7)] = hess.data[MATRIX_IDX(3, 2, 7)] = -N * x5_bar / sigma2;
    hess.data[MATRIX_IDX(2, 4, 7)] = hess.data[MATRIX_IDX(4, 2, 7)] = -N * x6_bar / sigma2;
    hess.data[MATRIX_IDX(2, 5, 7)] = hess.data[MATRIX_IDX(5, 2, 7)] = -N * x7_bar / sigma2;
    hess.data[MATRIX_IDX(2, 6, 7)] = hess.data[MATRIX_IDX(6, 2, 7)] = 2 * N / sigma3 * (v0 * x_sq_bar + b * x3_bar + c * x4_bar + d * x5_bar + e * x6_bar + f * x7_bar - y_x2_bar);

    hess.data[MATRIX_IDX(3, 3, 7)] = -N / sigma2 * x6_bar;
    hess.data[MATRIX_IDX(3, 4, 7)] = hess.data[MATRIX_IDX(4, 3, 7)] = -N / sigma2 * x7_bar;
    hess.data[MATRIX_IDX(3, 5, 7)] = hess.data[MATRIX_IDX(5, 3, 7)] = -N / sigma2 * x8_bar;
    hess.data[MATRIX_IDX(3, 6, 7)] = hess.data[MATRIX_IDX(6, 3, 7)] = 2 * N / sigma3 * (v0 * x3_bar + b * x4_bar + c * x5_bar + d * x6_bar + e * x7_bar + f * x8_bar - y_x3_bar);

    hess.data[MATRIX_IDX(4, 4, 7)] = -N / sigma2 * x8_bar;
    hess.data[MATRIX_IDX(4, 5, 7)] = hess.data[MATRIX_IDX(5, 4, 7)] = -N / sigma2 * x9_bar;
    hess.data[MATRIX_IDX(4, 6, 7)] = hess.data[MATRIX_IDX(6, 4, 7)] =  2 * N / sigma3 * (v0 * x4_bar + b * x5_bar + c * x6_bar + d * x7_bar + e * x8_bar + f * x9_bar - y_x4_bar);

    hess.data[MATRIX_IDX(5, 5, 7)] = -N / sigma2 * x10_bar;
    hess.data[MATRIX_IDX(5, 6, 7)] = hess.data[MATRIX_IDX(6, 5, 7)] =  2 * N / sigma3 * (v0 * x5_bar + b * x6_bar + c * x7_bar + d * x8_bar + e * x9_bar + f * x10_bar - y_x5_bar);

    hess.data[MATRIX_IDX(6, 6, 7)] = (N + 1.0) / sigma2 - 3.0 * N * Q / sigma4;

    for(size i = 0; i < num_parms * num_parms; ++i) { hess.data[i] *= p; }

    for(size i = 0; i < num_parms; ++i) { Assert(hess.data[MATRIX_IDX(i, i, 7)] < 0.0); }

    return hess;
}

UserData fifth_order_user_data = {0};

BayLaModel fifth_order_model = 
    {
        .model_prior_probability = 1.0 / 7.0,
        .n_parameters = 7,
        .prior = fifth_order_model_log_prior,
        .likelihood = fifth_order_model_log_likelihood,
        .posterior_max = fifth_order_model_max_a_posteriori,
        .hessian_matrix = fifth_order_model_2d_hessian,
        .user_data = &fifth_order_user_data
    };


/*------------------------------------------------  Initialize UserData  -------------------------------------------------*/
static inline void
initialize_global_data(void)
{
    BayLaNormalDistribution normal = bayla_normal_distribution_create(0.0, 0.25);
    ElkRandomState ran = elk_random_state_create(13);

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

    FILE *f = fopen("test_data.csv", "wb");

    /* Initialize the global data, with some noise! */
    for(size i = 0; i < NUM_DATA_POINTS; ++i)
    {
        // approximate log(x) from x = 0.1 to 3.0 */
        global_xs[i] = 0.1 + 2.9 * i / (NUM_DATA_POINTS - 1);
        global_ys[i] = natural_log_approximation_by_3rd_order_polynomial(global_xs[i], 0.9);
        f64 error = bayla_normal_distribution_random_deviate(&normal, &ran);
        global_ys[i] += error;

        fprintf(f, "%.16e, %.16e\n", global_xs[i], global_ys[i]);

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

    fclose(f);

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

    log_user_data = default_ud;
    constant_user_data = default_ud;
    linear_user_data = default_ud;
    second_order_user_data = default_ud;
    third_order_user_data = default_ud;
    fourth_order_user_data = default_ud;
    fifth_order_user_data = default_ud;
}

#undef NUM_DATA_POINTS

/*---------------------------------------------------  Test Models   -----------------------------------------------------*/
static inline void
test_log_model(MagAllocator alloc_, MagAllocator scratch1, MagAllocator scratch2)
{
    CoyProfileAnchor ap = {0};
    char *model_name = "Log";

    MagAllocator *alloc = &alloc_;

    ap = COY_START_PROFILE_BLOCK("  Log Model - sample");
    BayLaSamples samples = bayla_importance_sample_optimize(&log_model, 10000, 13, alloc, scratch1, scratch2);
    BayLaErrorValue z = bayla_samples_estimate_evidence(&samples);
    f64 ci_thresh = bayla_samples_calculate_ci_p_thresh(&samples, 0.68, scratch1);
    COY_END_PROFILE(ap);
    printf("%10s %4ld %9ld %6.0lf %15.2lf %11g ± %11g [%4.2lf%%]\n",
            model_name, samples.ndim, samples.n_samples, samples.neff, samples.neff / samples.n_samples,
            z.val, 3.0 * z.std, 3.0 * z.std / z.val * 100);

    bayla_samples_save_csv(&samples, ci_thresh, "log_model.csv");
}

static inline void
test_constant_model(MagAllocator alloc_, MagAllocator scratch1, MagAllocator scratch2)
{
    CoyProfileAnchor ap = {0};
    char *model_name = "Constant";

    MagAllocator *alloc = &alloc_;

    ap = COY_START_PROFILE_BLOCK("  Constant Model - sample");
    BayLaSamples samples = bayla_importance_sample_optimize(&constant_model, 10000, 13, alloc, scratch1, scratch2);
    BayLaErrorValue z = bayla_samples_estimate_evidence(&samples);
    f64 ci_thresh = bayla_samples_calculate_ci_p_thresh(&samples, 0.68, scratch1);
    COY_END_PROFILE(ap);
    printf("%10s %4ld %9ld %6.0lf %15.2lf %11g ± %11g [%4.2lf%%]\n",
            model_name, samples.ndim, samples.n_samples, samples.neff, samples.neff / samples.n_samples,
            z.val, 3.0 * z.std, 3.0 * z.std / z.val * 100);

    bayla_samples_save_csv(&samples, ci_thresh, "constant_model.csv");
}

static inline void
test_linear_model(MagAllocator alloc_, MagAllocator scratch1, MagAllocator scratch2)
{
    CoyProfileAnchor ap = {0};
    char *model_name = "Linear";

    MagAllocator *alloc = &alloc_;

    ap = COY_START_PROFILE_BLOCK("  Linear Model - sample");
    BayLaSamples samples = bayla_importance_sample_optimize(&linear_model, 10000, 13, alloc, scratch1, scratch2);
    BayLaErrorValue z = bayla_samples_estimate_evidence(&samples);
    f64 ci_thresh = bayla_samples_calculate_ci_p_thresh(&samples, 0.68, scratch1);
    COY_END_PROFILE(ap);
    printf("%10s %4ld %9ld %6.0lf %15.2lf %11g ± %11g [%4.2lf%%]\n",
            model_name, samples.ndim, samples.n_samples, samples.neff, samples.neff / samples.n_samples,
            z.val, 3.0 * z.std, 3.0 * z.std / z.val * 100);

    bayla_samples_save_csv(&samples, ci_thresh, "linear_model.csv");
}

static inline void
test_2nd_order_model(MagAllocator alloc_, MagAllocator scratch1, MagAllocator scratch2)
{
    CoyProfileAnchor ap = {0};
    char *model_name = "2nd Order";

    MagAllocator *alloc = &alloc_;

    ap = COY_START_PROFILE_BLOCK("  2nd Order Model - sample");
    BayLaSamples samples = bayla_importance_sample_optimize(&second_order_model, 10000, 13, alloc, scratch1, scratch2);
    BayLaErrorValue z = bayla_samples_estimate_evidence(&samples);
    f64 ci_thresh = bayla_samples_calculate_ci_p_thresh(&samples, 0.68, scratch1);
    COY_END_PROFILE(ap);
    printf("%10s %4ld %9ld %6.0lf %15.2lf %11g ± %11g [%4.2lf%%]\n",
            model_name, samples.ndim, samples.n_samples, samples.neff, samples.neff / samples.n_samples,
            z.val, 3.0 * z.std, 3.0 * z.std / z.val * 100);

    bayla_samples_save_csv(&samples, ci_thresh, "second_order_model.csv");
}

static inline void
test_3rd_order_model(MagAllocator alloc_, MagAllocator scratch1, MagAllocator scratch2)
{
    CoyProfileAnchor ap = {0};
    char *model_name = "3rd Order";

    MagAllocator *alloc = &alloc_;

    ap = COY_START_PROFILE_BLOCK("  3rd Order Model - sample");
    BayLaSamples samples = bayla_importance_sample_optimize(&third_order_model, 10000, 13, alloc, scratch1, scratch2);
    BayLaErrorValue z = bayla_samples_estimate_evidence(&samples);
    f64 ci_thresh = bayla_samples_calculate_ci_p_thresh(&samples, 0.68, scratch1);
    COY_END_PROFILE(ap);
    printf("%10s %4ld %9ld %6.0lf %15.2lf %11g ± %11g [%4.2lf%%]\n",
            model_name, samples.ndim, samples.n_samples, samples.neff, samples.neff / samples.n_samples,
            z.val, 3.0 * z.std, 3.0 * z.std / z.val * 100);

    bayla_samples_save_csv(&samples, ci_thresh, "third_order_model.csv");
}

static inline void
test_4th_order_model(MagAllocator alloc_, MagAllocator scratch1, MagAllocator scratch2)
{
    CoyProfileAnchor ap = {0};
    char *model_name = "4th Order";

    MagAllocator *alloc = &alloc_;

    ap = COY_START_PROFILE_BLOCK("  4th Order Model - sample");
    BayLaSamples samples = bayla_importance_sample_optimize(&fourth_order_model, 100000, 13, alloc, scratch1, scratch2);
    BayLaErrorValue z = bayla_samples_estimate_evidence(&samples);
    f64 ci_thresh = bayla_samples_calculate_ci_p_thresh(&samples, 0.68, scratch1);
    COY_END_PROFILE(ap);
    printf("%10s %4ld %9ld %6.0lf %15.2lf %11g ± %11g [%4.2lf%%]\n",
            model_name, samples.ndim, samples.n_samples, samples.neff, samples.neff / samples.n_samples,
            z.val, 3.0 * z.std, 3.0 * z.std / z.val * 100);

    bayla_samples_save_csv(&samples, ci_thresh, "fourth_order_model.csv");
}

static inline void
test_5th_order_model(MagAllocator alloc_, MagAllocator scratch1, MagAllocator scratch2)
{
    CoyProfileAnchor ap = {0};
    char *model_name = "5th Order";

    MagAllocator *alloc = &alloc_;

    ap = COY_START_PROFILE_BLOCK("  5th Order Model - sample");
    BayLaSamples samples = bayla_importance_sample_optimize(&fifth_order_model, 200000, 13, alloc, scratch1, scratch2);
    BayLaErrorValue z = bayla_samples_estimate_evidence(&samples);
    f64 ci_thresh = bayla_samples_calculate_ci_p_thresh(&samples, 0.68, scratch1);
    COY_END_PROFILE(ap);
    printf("%10s %4ld %9ld %6.0lf %15.2lf %11g ± %11g [%4.2lf%%]\n",
            model_name, samples.ndim, samples.n_samples, samples.neff, samples.neff / samples.n_samples,
            z.val, 3.0 * z.std, 3.0 * z.std / z.val * 100);

    bayla_samples_save_csv(&samples, ci_thresh, "fifth_order_model.csv");
}

/*---------------------------------------------------   All Tests    -----------------------------------------------------*/
static inline void
all_evidence_tests(void)
{
    printf("\n             Evidence Tests\n\n");

    MagAllocator alloc = mag_allocator_dyn_arena_create(ECO_MiB(2));
    MagAllocator scratch1 = mag_allocator_dyn_arena_create(ECO_MiB(2));
    MagAllocator scratch2 = mag_allocator_dyn_arena_create(ECO_MiB(2));

    CoyProfileAnchor ap_all = COY_START_PROFILE_BLOCK("Evidence Tests");
    CoyProfileAnchor ap = {0};

    ap = COY_START_PROFILE_BLOCK("Evidence Tests Init Data");
    initialize_global_data();
    COY_END_PROFILE(ap);

    printf("%-10s %-4s %-9s %-6s %-15s %-33s\n",
            "Model Name", "ndim", "n_samples", "neff", "effective ratio", "z_evidence");

    ap = COY_START_PROFILE_BLOCK("Evidence Tests Log Model");
    test_log_model(alloc, scratch1, scratch2);
    COY_END_PROFILE(ap);

    ap = COY_START_PROFILE_BLOCK("Evidence Tests Const Model");
    test_constant_model(alloc, scratch1, scratch2);
    COY_END_PROFILE(ap);

    ap = COY_START_PROFILE_BLOCK("Evidence Tests Linear Model");
    test_linear_model(alloc, scratch1, scratch2);
    COY_END_PROFILE(ap);

    ap = COY_START_PROFILE_BLOCK("Evidence Tests 2nd Order Model");
    test_2nd_order_model(alloc, scratch1, scratch2);
    COY_END_PROFILE(ap);

    ap = COY_START_PROFILE_BLOCK("Evidence Tests 3rd Order Model");
    test_3rd_order_model(alloc, scratch1, scratch2);
    COY_END_PROFILE(ap);

    ap = COY_START_PROFILE_BLOCK("Evidence Tests 4th Order Model");
    test_4th_order_model(alloc, scratch1, scratch2);
    COY_END_PROFILE(ap);

    ap = COY_START_PROFILE_BLOCK("Evidence Tests 5th Order Model");
    test_5th_order_model(alloc, scratch1, scratch2);
    COY_END_PROFILE(ap);

    eco_arena_destroy(&scratch2);
    eco_arena_destroy(&scratch1);
    eco_arena_destroy(&alloc);

    COY_END_PROFILE(ap_all);
}

