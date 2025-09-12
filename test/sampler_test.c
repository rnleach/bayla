/*---------------------------------------------------------------------------------------------------------------------------
 *                                      A simple, dummy model made of a 2d Gaussian
 *-------------------------------------------------------------------------------------------------------------------------*/

typedef struct
{
    BayLaMVNormDist dist;
    f64 *workspace;
} Gauss2DLikelihoodTestArgs;

static f64 
dummy_log_prior(size num_parms, f64 const *parms, void *user_data)
{
    return 0.0;
}

static f64
gauss_2d_likelihood(size num_parms, f64 const *parms, void *user_data)
{
    Assert(num_parms == 2);
    Gauss2DLikelihoodTestArgs *args = user_data;
    BayLaMVNormDist *dist = &args->dist;
    f64 *workspace = args->workspace;

    return bayla_mv_norm_dist_prob(dist, parms, workspace).val;
}

static void
gauss_2d_max_a_posteriori(void *user_data, size num_params, f64 *parms_out)
{
    Gauss2DLikelihoodTestArgs *args = user_data;
    BayLaMVNormDist *dist = &args->dist;

    memcpy(parms_out, dist->mean, sizeof(f64) * num_params);
}

static BayLaSquareMatrix
gauss_2d_hessian(void *user_data, size num_parms, f64 const *max_post_parms, MagAllocator *alloc)
{
    Gauss2DLikelihoodTestArgs *args = user_data;
    BayLaMVNormDist *dist = &args->dist;
    size const ndim = dist->ndim;

    BayLaSquareMatrix hess = bayla_square_matrix_create(ndim, alloc);
    memcpy(hess.data, dist->inv_cov.data, ndim * ndim * sizeof(f64));
    for(size i = 0; i < ndim *ndim; ++i)
    {
        hess.data[i] *= -1.0;
    }

    return hess;
}

static inline f64
get_param0(size ndim, f64 const *params, void *ud)
{
    return params[0];
}

static inline f64
get_param1(size ndim, f64 const *params, void *ud)
{
    return params[1];
}

static void
test_simple_model(void)
{
    printf("             ...Simple Model Tests\n");

    MagAllocator alloc_ = mag_allocator_dyn_arena_create(ECO_MiB(2));
    MagAllocator *alloc = &alloc_;
    MagStaticArena scratch_ = mag_static_arena_allocate_and_create(ECO_MiB(4));
    MagAllocator scratch = mag_allocator_from_static_arena(&scratch_);

    /* The true distribution. */
    f64 true_mean[2] = { -0.3, 0.8 };
    BayLaSquareMatrix true_cov = bayla_square_matrix_create(2, alloc);
    true_cov.data[0] = 2.0;
    true_cov.data[1] = -0.0;
    true_cov.data[2] = -0.0;
    true_cov.data[3] = 0.5;

    BayLaMVNormDist true_dist = bayla_mv_norm_dist_create(2, true_mean, true_cov, alloc);
    Assert(true_dist.valid);
    f64 *workspace = eco_arena_nmalloc(alloc, true_dist.ndim * 2, f64);

    Gauss2DLikelihoodTestArgs args = { .dist = true_dist, .workspace = workspace };

    BayLaModel model = 
        {
            .model_prior_probability = NAN,
            .n_parameters = true_dist.ndim,
            .prior = dummy_log_prior,
            .likelihood = gauss_2d_likelihood,
            .posterior_max = gauss_2d_max_a_posteriori,
            .hessian_matrix = gauss_2d_hessian,
            .user_data = &args
        };

    BayLaSamples samples = bayla_importance_sample(&model, 1.0, 10000, 13, alloc, scratch);
    BayLaErrorValue z = bayla_samples_estimate_evidence(&samples);
    printf("\nndim = %2ld n_samples = %4ld neff = %4.0lf effective ratio = %4.2lf z_evidence = %g ± %g [%4.2lf%%])\n",
            samples.ndim, samples.n_samples, samples.neff, samples.neff / samples.n_samples,
            z.val, 3.0 * z.std, 3.0 * z.std / z.val * 100);

    f64 p_thresh = bayla_samples_calculate_ci_p_thresh(&samples, 0.95, scratch);
    printf("p_thresh = %lf\n", p_thresh);

    BayLaCredibleInterval x0_ci = bayla_samples_calculate_ci(&samples, 0.95, 0, scratch);
    BayLaCredibleInterval x1_ci = bayla_samples_calculate_ci(&samples, 0.95, 1, scratch);

    BayLaErrorValue x0_exp = bayla_samples_calculate_expectation(&samples, get_param0, NULL);
    BayLaErrorValue x1_exp = bayla_samples_calculate_expectation(&samples, get_param1, NULL);

    printf("\n   Credible Intervals - x0: [ %lf -> %lf -> %lf]  x1: [ %lf -> %lf -> %lf ]\n",
            x0_ci.low, x0_exp.val, x0_ci.high, x1_ci.low, x1_exp.val, x1_ci.high);
    printf("Expectation (w/error) - x0: %lf ± %lf (%0.2lf%%) x1: %lf ± %lf (%0.2lf%%)\n",
            x0_exp.val, 3.0 * x0_exp.std, 100.0 * fabs(3.0 * x0_exp.std / x0_exp.val),
            x1_exp.val, 3.0 * x1_exp.std, 100.0 * fabs(3.0 * x1_exp.std / x1_exp.val));

    bayla_samples_save_csv(&samples, p_thresh, "simple_model.csv");

    eco_arena_destroy(&scratch);
    eco_arena_destroy(alloc);
}

/*---------------------------------------------------------------------------------------------------------------------------
 *                                               Crude Test of Components
 *-------------------------------------------------------------------------------------------------------------------------*/

void
crude_sampler_tests(void)
{
    printf("             ...Crude Tests\n");

    MagAllocator alloc_ = mag_allocator_dyn_arena_create(ECO_MiB(2));
    MagAllocator *alloc = &alloc_;

    ElkRandomState state_ = elk_random_state_create(12);
    ElkRandomState *state = &state_;

    /* The true distribution. */
    f64 true_mean[2] = { -0.3, 0.8 };
    BayLaSquareMatrix true_cov = bayla_square_matrix_create(2, alloc);
    true_cov.data[0] = 2.0;
    true_cov.data[1] = 0.0;
    true_cov.data[2] = 0.0;
    true_cov.data[3] = 0.5;

    BayLaMVNormDist true_dist = bayla_mv_norm_dist_create(2, true_mean, true_cov, alloc);
    Assert(true_dist.valid);

    /* The proposal distribution. */
    f64 prop_mean[2] = { 0.0, 0.0 };
    BayLaSquareMatrix prop_cov = bayla_square_matrix_create(2, alloc);
    prop_cov.data[0] = 1.0;
    prop_cov.data[1] = 0.0;
    prop_cov.data[2] = 0.0;
    prop_cov.data[3] = 1.0;

    BayLaMVNormDist prop_dist = bayla_mv_norm_dist_create(2, prop_mean, prop_cov, alloc);
    Assert(prop_dist.valid);

    size const n_samples_max = 1000000;
    size const n_cols = 2 + 3; /* ndim + 3, 3 for p, q, and w */
    f64 const ci_pct = 0.68;

    f64 *workspace = eco_arena_nmalloc(alloc, true_dist.ndim * 2, f64);
    f64 *rows = eco_arena_nmalloc(alloc, n_samples_max * n_cols, f64);
    f64 *row_scratch = eco_arena_nmalloc(alloc, n_samples_max * n_cols, f64);

    for(size n_samples = 10; n_samples <= n_samples_max; n_samples *= 10)
    {

        ElkKahanAccumulator x1_acc = {0};
        ElkKahanAccumulator x2_acc = {0};
        ElkKahanAccumulator w_acc = {0};
        ElkKahanAccumulator w2_acc = {0};

        for(size i = 0; i < n_samples; ++i)
        {
            f64 params[2] = {0};
            f64 *row = &rows[i * n_cols];
            BayLaLogValue ld_qi = bayla_mv_norm_dist_random_deviate(&prop_dist, state, params, workspace);
            BayLaLogValue ld_pi = bayla_mv_norm_dist_prob(&true_dist, params, workspace);

            row[0] = params[0];
            row[1] = params[1];
            row[2] = bayla_log_value_map_out_of_log_domain(ld_pi);
            row[3] = bayla_log_value_map_out_of_log_domain(ld_qi);
            row[4] = row[2] / row[3];

            w_acc = elk_kahan_accumulator_add(w_acc, row[4]);
            w2_acc = elk_kahan_accumulator_add(w2_acc, row[4] * row[4]);
            x1_acc = elk_kahan_accumulator_add(x1_acc, row[4] * row[0]);
            x2_acc = elk_kahan_accumulator_add(x2_acc, row[4] * row[1]);
        }

        f64 neff = w_acc.sum * w_acc.sum / w2_acc.sum;
        f64 z = w_acc.sum / n_samples;
        f64 mean_x1 = x1_acc.sum / w_acc.sum;
        f64 mean_x2 = x2_acc.sum / w_acc.sum;

        /* Get the credible interval (or area in this case) */
        pak_radix_sort(rows, n_samples, 2 * sizeof(f64), sizeof(f64) * n_cols, row_scratch, PAK_RADIX_SORT_F64, PAK_SORT_DESCENDING);
        size n_points = 0;
        ElkKahanAccumulator total_weight = {0};
        for(size r = 0; r < n_samples; ++r)
        {
            f64 *row = &rows[r * n_cols];
            total_weight = elk_kahan_accumulator_add(total_weight, row[4]);
#if 0
            printf("%9.6lf %9.6lf %9.6lf %9.6lf %9.6lf %9.6lf\n",
                    row[0], row[1], row[2], row[3], row[4], total_weight.sum / w_acc.sum);
#endif
            if(total_weight.sum >= ci_pct * w_acc.sum)
            {
                n_points = r + 1;
                break;
            }
        }

        if(n_points == 0) { n_points = n_samples; }

        printf("n = %7ld neff = %7.0lf ratio = %5.3lf Z = %13.6e x1 = %5.2lf x2 = %5.2lf\n",
                n_samples, neff, neff / n_samples, z, mean_x1, mean_x2);

#if 0
        char fname[100] = {0};
        snprintf(fname, sizeof(fname), "samples_%ld.csv", n_samples);
        FILE *f = fopen(fname, "wb");

        for(size n = n_samples - 1; n >= 0; --n)
        {
            f64 *row = &rows[n * n_cols];
            fprintf(f, "%.18e,%.18e,%.18e,%.18e,%.18e,%d\n",
                    row[0], row[1], row[2], row[3], row[4], n < n_points ? 1 : 0);
        }

        fclose(f);
#endif

    }

    mag_allocator_destroy(alloc);
}

/*--------------------------------------------------------------------------------------------------------------------------
 *                                                 Run All Sampler Tests
 *------------------------------------------------------------------------------------------------------------------------*/
static void
all_sampler_tests(void)
{
    printf("\n             Sampler Tests\n");

#if 0
    crude_sampler_tests();
#endif

#if 1
    test_simple_model();
#endif
}

