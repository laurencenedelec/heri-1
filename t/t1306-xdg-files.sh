#!/bin/sh
#
# Copyright (c) 2012 Valentin Duperray, Lucien Kong, Franck Jonas,
#		     Thomas Nguy, Khoi Nguyen
#		     Grenoble INP Ensimag
#

test_description='Compatibility with $XDG_CONFIG_HOME/git/ files'

. ./test-lib.sh

test_expect_success 'read config: xdg file exists and ~/.gitconfig doesn'\''t' '
	mkdir -p .config/git &&
	echo "[alias]" >.config/git/config &&
	echo "	myalias = !echo in_config" >>.config/git/config &&
	echo in_config >expected &&
	git myalias >actual &&
	test_cmp expected actual
'


test_expect_success 'read config: xdg file exists and ~/.gitconfig exists' '
	>.gitconfig &&
	echo "[alias]" >.gitconfig &&
	echo "	myalias = !echo in_gitconfig" >>.gitconfig &&
	echo in_gitconfig >expected &&
	git myalias >actual &&
	test_cmp expected actual
'


test_expect_success 'read with --get: xdg file exists and ~/.gitconfig doesn'\''t' '
	rm .gitconfig &&
	echo "[user]" >.config/git/config &&
	echo "	name = read_config" >>.config/git/config &&
	echo read_config >expected &&
	git config --get user.name >actual &&
	test_cmp expected actual
'


test_expect_success 'read with --get: xdg file exists and ~/.gitconfig exists' '
	>.gitconfig &&
	echo "[user]" >.gitconfig &&
	echo "	name = read_gitconfig" >>.gitconfig &&
	echo read_gitconfig >expected &&
	git config --get user.name >actual &&
	test_cmp expected actual
'


test_expect_success 'read with --list: xdg file exists and ~/.gitconfig doesn'\''t' '
	rm .gitconfig &&
	echo user.name=read_config >expected &&
	git config --global --list >actual &&
	test_cmp expected actual
'


test_expect_success 'read with --list: xdg file exists and ~/.gitconfig exists' '
	>.gitconfig &&
	echo "[user]" >.gitconfig &&
	echo "	name = read_gitconfig" >>.gitconfig &&
	echo user.name=read_gitconfig >expected &&
	git config --global --list >actual &&
	test_cmp expected actual
'


test_expect_success 'Setup' '
	git init git &&
	cd git &&
	echo foo >to_be_excluded
'


test_expect_success 'Exclusion of a file in the XDG ignore file' '
	mkdir -p "$HOME"/.config/git/ &&
	echo to_be_excluded >"$HOME"/.config/git/ignore &&
	test_must_fail git add to_be_excluded
'


test_expect_success 'Exclusion in both XDG and local ignore files' '
	echo to_be_excluded >.gitignore &&
	test_must_fail git add to_be_excluded
'


test_expect_success 'Exclusion in a non-XDG global ignore file' '
	rm .gitignore &&
	echo >"$HOME"/.config/git/ignore &&
	echo to_be_excluded >"$HOME"/my_gitignore &&
	git config core.excludesfile "$HOME"/my_gitignore &&
	test_must_fail git add to_be_excluded
'


test_done
