/*---------------------------------------------------   Pi Evidence  -----------------------------------------------------*/

f64 pi_min_parms[2] = {-1.0, -1.0};
f64 pi_max_parms[2] = { 1.0,  1.0};

static f64 
pi_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    return 0;
}

static f64 
pi_log_likelihood(size num_parms, f64 const *parms, void *user_data)
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
test_integrator_pi(void)
{
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
}

/*-------------------------------------------------   Known Integral  ----------------------------------------------------*/
#define KI_RANGE 50

f64 ki_min_parms[] = { -KI_RANGE, -KI_RANGE, -KI_RANGE, -KI_RANGE, -KI_RANGE, -KI_RANGE, -KI_RANGE };
f64 ki_max_parms[] = {  KI_RANGE,  KI_RANGE,  KI_RANGE,  KI_RANGE,  KI_RANGE,  KI_RANGE,  KI_RANGE };

static f64
ki_integral_true_value(size n_params)
{
    f64 result = pow(ELK_PI / 4.0, n_params / 2.0);
    for(size i = 0; i < n_params; ++i)
    {
        result *= erf(ki_max_parms[i]) - erf(ki_min_parms[i]);
    }

    return result;
}

static f64 
ki_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    return 0;
}

static f64 
ki_log_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    f64 sum_sq = 0.0;
    for(size i = 0; i < num_parms; ++i)
    {
        sum_sq += parms[i] * parms[i];
    }

    return -sum_sq;
}

BayLaModel ki_model = 
    {
        .model_prior_probability = 1.0,
        .n_parameters = ECO_ARRAY_SIZE(ki_min_parms),
        .prior = ki_log_prior,
        .likelihood = ki_log_likelihood,
        .min_parameter_vals = ki_min_parms,
        .max_parameter_vals = ki_max_parms,
        .user_data = NULL
    };

typedef struct
{
    BayLaModel model;
    BayLaIntegratorOptions opts;
    MagStaticArena scratch;
    BayLaEvidenceResult evidence;
    f64 true_error_pct;
    f64 true_error_z_score;
} TestEvidenceKiThreadData;

void
test_evidence_ki_thread_func(void *data)
{
    TestEvidenceKiThreadData *td = data;
    td->evidence = bayla_calculate_evidence(&td->model, &td->opts, td->scratch);

    size n_parms = td->model.n_parameters;
    f64 true_value = ki_integral_true_value(n_parms);
    f64 integral_estimate_value = bayla_log_value_map_out_of_log_domain(td->evidence.result.value);
    f64 integral_estimate_uncertainty = bayla_log_value_map_out_of_log_domain(td->evidence.result.error);

    td->true_error_pct = (integral_estimate_value - true_value) / true_value * 100.0;
    td->true_error_z_score = (integral_estimate_value - true_value) / integral_estimate_uncertainty;
}

static inline void
test_integrator_ki(void)
{
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

    TestEvidenceKiThreadData tds[ECO_ARRAY_SIZE(opts)] = {0};
    CoyThread threads[ECO_ARRAY_SIZE(opts)] = {0};
    MagStaticArena scratches[ECO_ARRAY_SIZE(opts)] = {0};
    for(size method = 0; method < ECO_ARRAY_SIZE(opts); ++method)
    {
        scratches[method] = mag_static_arena_allocate_and_create(ECO_KiB(512));
    }
    
    b32 success = true;
    for(size method = 0; method < ECO_ARRAY_SIZE(opts); ++method)
    {
        tds[method].model = ki_model;
        tds[method].opts = opts[method];
        tds[method].scratch = mag_static_arena_borrow(&scratches[method]);

        success &= coy_thread_create(&threads[method], test_evidence_ki_thread_func, &tds[method]);
        Assert(success);
    }

    printf("\nCalculating a known integral as a test of the integrators.\n");
    printf("%20s - %8s ± %8s [%9s] %16s [%9s] %8s [ %8s ]\n",
            "Method", "Value", "Error", "Pct Err", "True Err Z score", "Pct True", "Samples", "Samples Used");

    for(size method = 0; method < ECO_ARRAY_SIZE(opts); ++method)
    {
        success = coy_thread_join(&threads[method]);
        Assert(success);
        coy_thread_destroy(&threads[method]);

        printf("%20s - %lf ± %lf [%8.4lf%%] %16.4lf [%8.4lf%%]",
                names[method],
                bayla_log_value_map_out_of_log_domain(tds[method].evidence.result.value), 
                bayla_log_value_map_out_of_log_domain(tds[method].evidence.result.error),
                100.0 * bayla_log_value_map_out_of_log_domain(bayla_log_value_divide(
                        tds[method].evidence.result.error, tds[method].evidence.result.value)),
                tds[method].true_error_z_score, tds[method].true_error_pct);

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
}


/*---------------------------------------------------   All Tests    -----------------------------------------------------*/
static inline void
all_integrator_tests(void)
{
    test_integrator_ki();
    test_integrator_pi();
}

