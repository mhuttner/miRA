#include <stdlib.h>
#include <check.h>
#include <stdio.h>
#include "../src/parse_sam.h"

START_TEST (test_parse_line)
{
    char sample_line[1000]="Seq23599_x1	272	scaffold_1	67	0	21M	*	0	0	GCCACCCATGCCGCATCCACA	IIIIIIIIIIIIIIIIIIIII	AS:i:-6	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:12A8	YT:Z:UU	NH:i:16	CC:Z:=	CP:i:250	HI:i:0";
    struct sam_entry test_entry;
    parse_line(&test_entry,sample_line,1000);
    ck_assert_int_eq(test_entry.flag,272);
    ck_assert(test_entry.pos == 67);
}
END_TEST

Suite *hello_suite( void )
{
    Suite *s = suite_create("Sum");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_parse_line);
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
