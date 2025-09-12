static inline b32
approx_equal_abs(f64 left, f64 right, f64 tolerance)
{
    return fabs(left - right) <= tolerance;
}

static inline b32
approx_equal_rel(f64 left, f64 right, f64 tolerance)
{
    f64 denom = left + (left - right) / 2.0;
    if(fabs(denom) > 0.0)
    {
        return fabs(left - right) / denom <= tolerance; 
    }
    else
    {
        return true; 
    }
}

/*--------------------------------------------------------------------------------------------------------------------------
 *                                                 Log-Domain Math Tests
 *------------------------------------------------------------------------------------------------------------------------*/
static inline void
basic_functions(void)
{
    Assert(bayla_log1mexp(0.0) == - INFINITY);
    Assert(bayla_log1mexp(INFINITY) == 0.0);
    Assert(isnan(bayla_log1mexp(-1.0e-100)));

    Assert(bayla_log1pexp(- INFINITY) == 0.0);
    Assert(bayla_log1pexp(INFINITY) == INFINITY);
    Assert(bayla_log1pexp(0.0) == ELK_LN2);

    /* All inputs and outputs are in the log domain, but the additions are in the decimal domain. */
    Assert(bayla_logsumexp(- INFINITY, - INFINITY) == - INFINITY);
    Assert(bayla_logsumexp(0.0, - INFINITY) == 0.0);
    Assert(bayla_logsumexp(1.0, - INFINITY) == 1.0);
    Assert(bayla_logsumexp(0.0, 0.0) == ELK_LN2);
    Assert(bayla_logsumexp(- INFINITY, INFINITY) == INFINITY);

    f64 vals[] = {- INFINITY, 0.0, 0.0 };
    Assert(bayla_logsumexp_list(sizeof(vals) / sizeof(vals[0]), 0, 1, vals) == ELK_LN2);
}

static inline void
wrap_and_unwrap(void)
{
    /* Wrapping a value already in the log domain should have no effect. */
    BayLaLogValue lv = bayla_log_value_create(0.0);
    Assert(lv.val == 0.0);

    /* One to zero and back should round trip. */
    lv = bayla_log_value_map_into_log_domain(1.0);
    Assert(lv.val == 0.0);
    f64 dv = bayla_log_value_map_out_of_log_domain(lv);
    Assert(dv == 1.0);

    /* Zero to negative infinity and back should round trip. */
    lv = bayla_log_value_map_into_log_domain(0.0);
    Assert(lv.val == - INFINITY);
    dv = bayla_log_value_map_out_of_log_domain(lv);
    Assert(dv == 0.0);

    /* 'e' to one and back should round trip. */
    lv = bayla_log_value_map_into_log_domain(ELK_E);
    Assert(lv.val == 1.0);
    dv = bayla_log_value_map_out_of_log_domain(lv);
    Assert(dv == ELK_E);

    /* Infinity should round trip. */
    lv = bayla_log_value_map_into_log_domain(INFINITY);
    Assert(lv.val == INFINITY);
    dv = bayla_log_value_map_out_of_log_domain(lv);
    Assert(dv == INFINITY);

    /* Anything that leads to a NAN is not a round trip. */
    lv = bayla_log_value_map_into_log_domain(-1.0e-100);
    Assert(isnan(lv.val));
    dv = bayla_log_value_map_out_of_log_domain(lv);
    Assert(isnan(dv));
}

static inline void
addition_and_subtraction(void)
{
    BayLaLogValue zero = bayla_log_value_map_into_log_domain(0.0);
    BayLaLogValue one = bayla_log_value_map_into_log_domain(1.0);
    BayLaLogValue euler = bayla_log_value_map_into_log_domain(ELK_E);
    BayLaLogValue infinity = bayla_log_value_map_into_log_domain(INFINITY);

    /* Addition. */
    Assert(bayla_log_value_add(zero, one).val == one.val);
    Assert(bayla_log_value_add(zero, euler).val == euler.val);
    Assert(bayla_log_value_add(zero, infinity).val == infinity.val);

    Assert(bayla_log_value_add(one, zero).val == one.val);
    Assert(bayla_log_value_add(euler, zero).val == euler.val);
    Assert(bayla_log_value_add(infinity, zero).val == infinity.val);

    Assert(bayla_log_value_add(infinity, zero).val == infinity.val);
    Assert(bayla_log_value_add(infinity, one).val == infinity.val);
    Assert(bayla_log_value_add(infinity, euler).val == infinity.val);

    Assert(isnan(bayla_log_value_add(infinity, infinity).val));

    Assert(bayla_log_value_add(one, one).val == bayla_log_value_map_into_log_domain(2.0).val);

    f64 left = bayla_log_value_add(euler, euler).val;
    f64 right = bayla_log_value_map_into_log_domain(2.0 * ELK_E).val;
    Assert(approx_equal_abs(left, right, 1.0e-15));
    Assert(approx_equal_rel(left, right, 1.0e-15));

    /* Subtraction. */
    Assert(bayla_log_value_subtract(zero, zero).val == zero.val);
    Assert(bayla_log_value_subtract(one, one).val == zero.val);
    Assert(bayla_log_value_subtract(euler, euler).val == zero.val);

    Assert(isnan(bayla_log_value_subtract(zero, one).val));      /* Must result in a positive value. */
    Assert(isnan(bayla_log_value_subtract(zero, euler).val));    /* Must result in a positive value. */
    Assert(isnan(bayla_log_value_subtract(zero, infinity).val)); /* Must result in a positive value. */

    Assert(bayla_log_value_subtract(one, zero).val == one.val);
    Assert(bayla_log_value_subtract(euler, zero).val == euler.val);
    Assert(bayla_log_value_subtract(infinity, zero).val == infinity.val);

    Assert(bayla_log_value_subtract(infinity, zero).val == infinity.val);
    Assert(bayla_log_value_subtract(infinity, one).val == infinity.val);
    Assert(bayla_log_value_subtract(infinity, euler).val == infinity.val);

    Assert(isnan(bayla_log_value_subtract(infinity, infinity).val));

    Assert(bayla_log_value_subtract(one, one).val == bayla_log_value_map_into_log_domain(0.0).val);
    Assert(bayla_log_value_subtract(euler, euler).val == zero.val);
}

static inline void
multiplication_and_division(void)
{
    BayLaLogValue zero = bayla_log_value_map_into_log_domain(0.0);
    BayLaLogValue half = bayla_log_value_map_into_log_domain(0.5);
    BayLaLogValue one = bayla_log_value_map_into_log_domain(1.0);
    BayLaLogValue two = bayla_log_value_map_into_log_domain(2.0);
    BayLaLogValue euler = bayla_log_value_map_into_log_domain(ELK_E);
    BayLaLogValue four = bayla_log_value_map_into_log_domain(4.0);
    BayLaLogValue infinity = bayla_log_value_map_into_log_domain(INFINITY);

    /* Multiplication. */
    Assert(bayla_log_value_multiply(zero, one).val == zero.val);
    Assert(bayla_log_value_multiply(one, zero).val == zero.val);
    Assert(bayla_log_value_multiply(zero, two).val == zero.val);
    Assert(bayla_log_value_multiply(two, zero).val == zero.val);
    Assert(bayla_log_value_multiply(zero, euler).val == zero.val);
    Assert(bayla_log_value_multiply(euler, zero).val == zero.val);

    Assert(isnan(bayla_log_value_multiply(infinity, zero).val));
    Assert(isnan(bayla_log_value_multiply(zero, infinity).val));

    Assert(bayla_log_value_multiply(one, one).val == one.val);
    Assert(bayla_log_value_multiply(euler, one).val == euler.val);
    Assert(bayla_log_value_multiply(one, euler).val == euler.val);
    Assert(bayla_log_value_multiply(two, one).val == two.val);
    Assert(bayla_log_value_multiply(one, two).val == two.val);

    Assert(bayla_log_value_multiply(one, infinity).val == infinity.val);
    Assert(bayla_log_value_multiply(infinity, one).val == infinity.val);
    Assert(bayla_log_value_multiply(two, infinity).val == infinity.val);
    Assert(bayla_log_value_multiply(infinity, two).val == infinity.val);
    Assert(bayla_log_value_multiply(euler, infinity).val == infinity.val);
    Assert(bayla_log_value_multiply(infinity, euler).val == infinity.val);

    Assert(approx_equal_abs(bayla_log_value_multiply(two, two).val, four.val, 1.0e-315));
    Assert(approx_equal_rel(bayla_log_value_multiply(two, two).val, four.val, 1.0e-315));

    /* Division */
    Assert(bayla_log_value_divide(zero, one).val == zero.val);
    Assert(bayla_log_value_divide(one, zero).val == INFINITY);
    Assert(bayla_log_value_divide(zero, two).val == zero.val);
    Assert(bayla_log_value_divide(two, zero).val == INFINITY);
    Assert(bayla_log_value_divide(zero, euler).val == zero.val);
    Assert(bayla_log_value_divide(euler, zero).val == INFINITY);

    Assert(bayla_log_value_divide(infinity, zero).val == INFINITY);
    Assert(bayla_log_value_divide(zero, infinity).val == zero.val);

    Assert(bayla_log_value_divide(one, one).val == one.val);
    Assert(bayla_log_value_divide(euler, one).val == euler.val);
    Assert(bayla_log_value_divide(one, euler).val != euler.val);
    Assert(bayla_log_value_divide(two, one).val == two.val);
    Assert(bayla_log_value_divide(one, two).val != two.val);
    Assert(bayla_log_value_divide(one, two).val == half.val);

    Assert(bayla_log_value_divide(one, infinity).val == zero.val);
    Assert(bayla_log_value_divide(infinity, one).val == infinity.val);
    Assert(bayla_log_value_divide(two, infinity).val == zero.val);
    Assert(bayla_log_value_divide(infinity, two).val == infinity.val);
    Assert(bayla_log_value_divide(euler, infinity).val == zero.val);
    Assert(bayla_log_value_divide(infinity, euler).val == infinity.val);

    Assert(approx_equal_abs(bayla_log_value_divide(one, one).val, one.val, 1.0e-315));
    Assert(approx_equal_abs(bayla_log_value_divide(two, two).val, one.val, 1.0e-315));
    Assert(approx_equal_abs(bayla_log_value_divide(euler, euler).val, one.val, 1.0e-315));

}

static inline void
reciprocal_and_power(void)
{
    BayLaLogValue zero = bayla_log_value_map_into_log_domain(0.0);
    BayLaLogValue half = bayla_log_value_map_into_log_domain(0.5);
    BayLaLogValue one = bayla_log_value_map_into_log_domain(1.0);
    BayLaLogValue two = bayla_log_value_map_into_log_domain(2.0);
    BayLaLogValue four = bayla_log_value_map_into_log_domain(4.0);
    BayLaLogValue infinity = bayla_log_value_map_into_log_domain(INFINITY);

    /* Reciprocal. */
    Assert(bayla_log_value_reciprocal(zero).val == infinity.val);
    Assert(bayla_log_value_reciprocal(one).val == one.val);
    Assert(bayla_log_value_reciprocal(two).val == half.val);
    Assert(bayla_log_value_reciprocal(infinity).val == zero.val);

    /* Powers. */
    Assert(bayla_log_value_power(half, 0.0).val == one.val);
    Assert(bayla_log_value_power(one, 0.0).val == one.val);
    Assert(bayla_log_value_power(two, 0.0).val == one.val);
    Assert(isnan(bayla_log_value_power(infinity, 0.0).val));

    Assert(bayla_log_value_power(one, 0.5).val == one.val);
    Assert(bayla_log_value_power(one, 1.0).val == one.val);
    Assert(bayla_log_value_power(one, 2.0).val == one.val);

    Assert(bayla_log_value_power(two, 2.0).val == four.val);
}

/*--------------------------------------------------------------------------------------------------------------------------
 *                                             Run All Log-Domain Math Tests
 *------------------------------------------------------------------------------------------------------------------------*/
static void
all_log_math_tests(void)
{
    printf("\n             Log-Domain Math Tests\n");

    basic_functions();
    wrap_and_unwrap();
    addition_and_subtraction();
    multiplication_and_division();
    reciprocal_and_power();
}
