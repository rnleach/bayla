
/*--------------------------------------------------------------------------------------------------------------------------
 *                                                 Run All Sampler Tests
 *------------------------------------------------------------------------------------------------------------------------*/

#include <float.h>
#include "../lib/packrat.h"

static void
all_sampler_tests(void)
{
    printf("\n             Sampler Tests\n");

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

        /* Get the credible x1 interval */
        pak_radix_sort(rows, n_samples, 0 * sizeof(f64), sizeof(f64) * n_cols, row_scratch, PAK_RADIX_SORT_F64, PAK_SORT_ASCENDING);
        f64 min_range = DBL_MAX;
        f64 xstart = NAN;
        f64 xend = NAN;

        size l = 0;
        size r = 1;
        f64 *lrow =  &rows[l * n_cols];
        ElkKahanAccumulator total_weight = elk_kahan_accumulator_add((ElkKahanAccumulator){0}, lrow[4]);
        while(l < n_samples)   
        {
            while(r < n_samples)
            {
                f64 *rrow = &rows[r * n_cols];
                total_weight = elk_kahan_accumulator_add(total_weight, rrow[4]);

                if(total_weight.sum >= ci_pct * w_acc.sum)
                {
                    f64 range = rrow[0] - lrow[0];
                    if(range < min_range)
                    {
                        min_range = range;
                        xstart = lrow[0];
                        xend = rrow[0];
                    }
                    ++r;
                    break;
                }

                ++r;
            }

            /* Remove the weight from the left most point. */
            while(total_weight.sum >= ci_pct * w_acc.sum && l < n_samples)
            {
                total_weight = elk_kahan_accumulator_add(total_weight, -lrow[4]);
                ++l;
                lrow =  &rows[l * n_cols];
            }

            if(r >= n_samples) { break; }
        }

        /* Get the credible interval (or area in this case) */
        pak_radix_sort(rows, n_samples, 2 * sizeof(f64), sizeof(f64) * n_cols, row_scratch, PAK_RADIX_SORT_F64, PAK_SORT_DESCENDING);
        size n_points = 0;
        total_weight = (ElkKahanAccumulator){0};
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

#if 0
        printf("n_points / n_samples = %lf\n", (f64)n_points / (f64)n_samples);
#endif

        f64 x1_min = DBL_MAX;
        f64 x1_max = -DBL_MAX;
        f64 x2_min = DBL_MAX;
        f64 x2_max = -DBL_MAX;

        for(size r = 0; r < n_points; ++r)
        {
            f64 *row = &rows[r * n_cols];
            x1_min = x1_min < row[0] ? x1_min : row[0];
            x2_min = x2_min < row[1] ? x2_min : row[1];
            x1_max = x1_max > row[0] ? x1_max : row[0];
            x2_max = x2_max > row[1] ? x2_max : row[1];
        }

        printf("n = %7ld neff = %7.0lf ratio = %5.3lf Z = %13.6e x1 = %5.2lf [%5.2lf - %5.2lf] x2 = %5.2lf [%5.2lf - %5.2lf]",
                n_samples, neff, neff / n_samples, z, mean_x1, x1_min, x1_max, mean_x2, x2_min, x2_max);
        printf(" [%5.2lf - %5.2lf]\n", xstart, xend);

#if 1
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
