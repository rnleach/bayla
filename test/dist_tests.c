void
distribution_tests(void)
{
    ElkRandomState state_ = elk_random_state_create(11);
    ElkRandomState *state = &state_;

    MagAllocator alloc_ = mag_allocator_dyn_arena_create(ECO_KiB(100));
    MagAllocator *alloc = &alloc_;

    f64 mean[2] = {-1230.0, +0.0123};
    f64 workspace[100] = {0};
    f64 rho = 0.75;

    BayLaSquareMatrix cov = bayla_square_matrix_create(2, alloc);

    f64 s0 = 0.0025;
    f64 s1 = 100.0;

    /* Diagonal elements. */
    cov.data[0] = s0 * s0;;
    cov.data[3] = s1 * s1;;

    cov.data[1] = rho * s0 * s1;
    cov.data[2] = cov.data[1];

    BayLaMVNormDist dist = bayla_mv_norm_dist_create(2, mean, cov, alloc);

    for(size i = 0; i < 100; ++i)
    {
        f64 deviate[2] = {0};
        BayLaLogValue lp1 = bayla_mv_norm_dist_random_deviate(&dist, state, deviate, workspace);
        BayLaLogValue lp2 = bayla_mv_norm_dist_prob(&dist, deviate, workspace);

        f64 diff = bayla_log_value_map_out_of_log_domain(lp1) - bayla_log_value_map_out_of_log_domain(lp2);
        f64 ratio = bayla_log_value_map_out_of_log_domain(bayla_log_value_divide(lp1, lp2));
#if 0
        f64 log_ratio = bayla_log_value_divide(lp1, lp2).val;

        printf("%10.6g ==? %10.6g for (%16.12g, %16.12g) ratio %e log ratio %13.6e diff %13.6e\n",
                lp1.val, lp2.val, deviate[0], deviate[1], ratio, log_ratio, diff);
#else
        Assert(fabs(diff) <= 1.0e-10);
        Assert(fabs(fabs(ratio) - 1.0) <= 1.0e-8);
#endif
    }

    mag_allocator_destroy(alloc);
}

/*--------------------------------------------------------------------------------------------------------------------------
 *                                            Run All Distribution Math Tests
 *------------------------------------------------------------------------------------------------------------------------*/
static inline void
all_distribution_math_tests(void)
{
    printf("\n             Distribution Math Tests\n");

    distribution_tests();
}
