#include <stdlib.h>
#include <check.h>


START_TEST (test_sum)
{
    ck_assert_int_eq (3, 3);
}
END_TEST

Suite *hello_suite( void )
{
    Suite *s = suite_create("Sum");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_sum);
    suite_add_tcase(s, tc_core);

    return s;
}

int main( void )
{
    int number_failed;
    Suite *s = hello_suite();
    SRunner *sr = srunner_create(s);
    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
