#include <inttypes.h>
#include <tgmath.h>
#include <stdio.h>

/* We must have asserts working for the tests to work. */
#ifdef NDEBUG
#    undef NDEBUG
#endif

#include "../lib/coyote.h"
#include "../lib/bayla.h"

#include "log_math.c"
#include "dist_tests.c"
#include "sampler_test.c"
#include "evidence_test.c"

int
main(int argc, char *argv[])
{
    coy_profile_begin();
    printf("  ---------- Starting Tests ----------\n");

    all_log_math_tests();
    all_distribution_math_tests();
    all_sampler_tests();
    all_evidence_tests();

    printf("\n  ---------- Finished Tests ----------\n");
    coy_profile_end();

#if COY_PROFILE && 0
    printf("Total Runtime = %.3lf seconds at a frequency of %"PRIu64"\n",
            coy_global_profiler.total_elapsed, coy_global_profiler.freq);

    for(i32 i = 0; i < COY_PROFILE_NUM_BLOCKS; ++i)
    {
        CoyBlockProfiler *block = &coy_global_profiler.blocks[i];
        if(block->hit_count)
        {
            printf("%-28s Hits: %3"PRIu64" Exclusive: %6.2lf%%", block->label, block->hit_count, block->exclusive_pct);
            
            if(block->inclusive_pct != block->exclusive_pct)
            {
                printf(" Inclusive: %6.2lf%%\n", block->inclusive_pct);
            }
            else
            {
                printf("\n");
            }
        }
    }
#endif
}

COY_PROFILE_STATIC_CHECK;

