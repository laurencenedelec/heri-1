#!/bin/sh

test_description='log --grep/--author/--regexp-ignore-case/-S/-G'
. ./test-lib.sh

test_expect_success setup '
	>expect_nomatch &&

	>file &&
	git add file &&
	test_tick &&
	git commit -m initial &&
	git rev-parse --verify HEAD >expect_initial &&

	echo Picked >file &&
	git add file &&
	test_tick &&
	git commit --author="Another Person <another@example.com>" -m second &&
	git rev-parse --verify HEAD >expect_second
'

test_expect_success 'log --grep' '
	git log --grep=initial --format=%H >actual &&
	test_cmp expect_initial actual
'

test_expect_success 'log --grep --regexp-ignore-case' '
	git log --regexp-ignore-case --grep=InItial --format=%H >actual &&
	test_cmp expect_initial actual
'

test_expect_success 'log --grep -i' '
	git log -i --grep=InItial --format=%H >actual &&
	test_cmp expect_initial actual
'

test_expect_success 'log --author --regexp-ignore-case' '
	git log --regexp-ignore-case --author=person --format=%H >actual &&
	test_cmp expect_second actual
'

test_expect_success 'log --author -i' '
	git log -i --author=person --format=%H >actual &&
	test_cmp expect_second actual
'

test_expect_success 'log -G (nomatch)' '
	git log -Gpicked --format=%H >actual &&
	test_cmp expect_nomatch actual
'

test_expect_success 'log -G (match)' '
	git log -GPicked --format=%H >actual &&
	test_cmp expect_second actual
'

test_expect_success 'log -G --regexp-ignore-case (nomatch)' '
	git log --regexp-ignore-case -Gpickle --format=%H >actual &&
	test_cmp expect_nomatch actual
'

test_expect_success 'log -G -i (nomatch)' '
	git log -i -Gpickle --format=%H >actual &&
	test_cmp expect_nomatch actual
'

test_expect_success 'log -G --regexp-ignore-case (match)' '
	git log --regexp-ignore-case -Gpicked --format=%H >actual &&
	test_cmp expect_second actual
'

test_expect_success 'log -G -i (match)' '
	git log -i -Gpicked --format=%H >actual &&
	test_cmp expect_second actual
'

test_expect_success 'log -G --textconv (missing textconv tool)' '
	echo "* diff=test" >.gitattributes &&
	test_must_fail git -c diff.test.textconv=missing log -Gfoo &&
	rm .gitattributes
'

test_expect_success 'log -G --no-textconv (missing textconv tool)' '
	echo "* diff=test" >.gitattributes &&
	git -c diff.test.textconv=missing log -Gfoo --no-textconv >actual &&
	test_cmp expect_nomatch actual &&
	rm .gitattributes
'

test_expect_success 'log -S (nomatch)' '
	git log -Spicked --format=%H >actual &&
	test_cmp expect_nomatch actual
'

test_expect_success 'log -S (match)' '
	git log -SPicked --format=%H >actual &&
	test_cmp expect_second actual
'

test_expect_success 'log -S --regexp-ignore-case (match)' '
	git log --regexp-ignore-case -Spicked --format=%H >actual &&
	test_cmp expect_second actual
'

test_expect_success 'log -S -i (match)' '
	git log -i -Spicked --format=%H >actual &&
	test_cmp expect_second actual
'

test_expect_success 'log -S --regexp-ignore-case (nomatch)' '
	git log --regexp-ignore-case -Spickle --format=%H >actual &&
	test_cmp expect_nomatch actual
'

test_expect_success 'log -S -i (nomatch)' '
	git log -i -Spickle --format=%H >actual &&
	test_cmp expect_nomatch actual
'

test_expect_success 'log -S --textconv (missing textconv tool)' '
	echo "* diff=test" >.gitattributes &&
	test_must_fail git -c diff.test.textconv=missing log -Sfoo &&
	rm .gitattributes
'

test_expect_success 'log -S --no-textconv (missing textconv tool)' '
	echo "* diff=test" >.gitattributes &&
	git -c diff.test.textconv=missing log -Sfoo --no-textconv >actual &&
	test_cmp expect_nomatch actual &&
	rm .gitattributes
'

test_done
