/*
 * "git reset" builtin command
 *
 * Copyright (c) 2007 Carlos Rica
 *
 * Based on git-reset.sh, which is
 *
 * Copyright (c) 2005, 2006 Linus Torvalds and Junio C Hamano
 */
#include "builtin.h"
#include "tag.h"
#include "object.h"
#include "commit.h"
#include "run-command.h"
#include "refs.h"
#include "diff.h"
#include "diffcore.h"
#include "tree.h"
#include "branch.h"
#include "parse-options.h"
#include "unpack-trees.h"
#include "cache-tree.h"

static const char * const git_reset_usage[] = {
	N_("git reset [--mixed | --soft | --hard | --merge | --keep] [-q] [<commit>]"),
	N_("git reset [-q] <commit> [--] <paths>..."),
	N_("git reset --patch [<commit>] [--] [<paths>...]"),
	NULL
};

enum reset_type { MIXED, SOFT, HARD, MERGE, KEEP, NONE };
static const char *reset_type_names[] = {
	N_("mixed"), N_("soft"), N_("hard"), N_("merge"), N_("keep"), NULL
};

static inline int is_merge(void)
{
	return !access(git_path("MERGE_HEAD"), F_OK);
}

static int reset_index(const unsigned char *sha1, int reset_type, int quiet)
{
	int nr = 1;
	struct tree_desc desc[2];
	struct tree *tree;
	struct unpack_trees_options opts;

	memset(&opts, 0, sizeof(opts));
	opts.head_idx = 1;
	opts.src_index = &the_index;
	opts.dst_index = &the_index;
	opts.fn = oneway_merge;
	opts.merge = 1;
	if (!quiet)
		opts.verbose_update = 1;
	switch (reset_type) {
	case KEEP:
	case MERGE:
		opts.update = 1;
		break;
	case HARD:
		opts.update = 1;
		/* fallthrough */
	default:
		opts.reset = 1;
	}

	read_cache_unmerged();

	if (reset_type == KEEP) {
		unsigned char head_sha1[20];
		if (get_sha1("HEAD", head_sha1))
			return error(_("You do not have a valid HEAD."));
		if (!fill_tree_descriptor(desc, head_sha1))
			return error(_("Failed to find tree of HEAD."));
		nr++;
		opts.fn = twoway_merge;
	}

	if (!fill_tree_descriptor(desc + nr - 1, sha1))
		return error(_("Failed to find tree of %s."), sha1_to_hex(sha1));
	if (unpack_trees(nr, desc, &opts))
		return -1;

	if (reset_type == MIXED || reset_type == HARD) {
		tree = parse_tree_indirect(sha1);
		prime_cache_tree(&active_cache_tree, tree);
	}

	return 0;
}

static void print_new_head_line(struct commit *commit)
{
	const char *hex, *body;

	hex = find_unique_abbrev(commit->object.sha1, DEFAULT_ABBREV);
	printf(_("HEAD is now at %s"), hex);
	body = strstr(commit->buffer, "\n\n");
	if (body) {
		const char *eol;
		size_t len;
		body += 2;
		eol = strchr(body, '\n');
		len = eol ? eol - body : strlen(body);
		printf(" %.*s\n", (int) len, body);
	}
	else
		printf("\n");
}

static int update_index_refresh(int fd, struct lock_file *index_lock, int flags)
{
	if (!index_lock) {
		index_lock = xcalloc(1, sizeof(struct lock_file));
		fd = hold_locked_index(index_lock, 1);
	}

	refresh_index(&the_index, (flags), NULL, NULL,
		      _("Unstaged changes after reset:"));
	if (write_cache(fd, active_cache, active_nr) ||
			commit_locked_index(index_lock))
		return error ("Could not refresh index");
	return 0;
}

static void update_index_from_diff(struct diff_queue_struct *q,
		struct diff_options *opt, void *data)
{
	int i;

	for (i = 0; i < q->nr; i++) {
		struct diff_filespec *one = q->queue[i]->one;
		if (one->mode && !is_null_sha1(one->sha1)) {
			struct cache_entry *ce;
			ce = make_cache_entry(one->mode, one->sha1, one->path,
				0, 0);
			if (!ce)
				die(_("make_cache_entry failed for path '%s'"),
				    one->path);
			add_cache_entry(ce, ADD_CACHE_OK_TO_ADD |
				ADD_CACHE_OK_TO_REPLACE);
		} else
			remove_file_from_cache(one->path);
	}
}

static int read_from_tree(const char **pathspec, unsigned char *tree_sha1,
			  int refresh_flags)
{
	struct lock_file *lock = xcalloc(1, sizeof(struct lock_file));
	int index_fd;
	struct diff_options opt;

	memset(&opt, 0, sizeof(opt));
	diff_tree_setup_paths(pathspec, &opt);
	opt.output_format = DIFF_FORMAT_CALLBACK;
	opt.format_callback = update_index_from_diff;

	index_fd = hold_locked_index(lock, 1);
	read_cache();
	if (do_diff_cache(tree_sha1, &opt))
		return 1;
	diffcore_std(&opt);
	diff_flush(&opt);
	diff_tree_release_paths(&opt);

	return update_index_refresh(index_fd, lock, refresh_flags);
}

static void set_reflog_message(struct strbuf *sb, const char *action,
			       const char *rev)
{
	const char *rla = getenv("GIT_REFLOG_ACTION");

	strbuf_reset(sb);
	if (rla)
		strbuf_addf(sb, "%s: %s", rla, action);
	else if (rev)
		strbuf_addf(sb, "reset: moving to %s", rev);
	else
		strbuf_addf(sb, "reset: %s", action);
}

static void die_if_unmerged_cache(int reset_type)
{
	if (is_merge() || read_cache() < 0 || unmerged_cache())
		die(_("Cannot do a %s reset in the middle of a merge."),
		    _(reset_type_names[reset_type]));

}

static const char **parse_args(const char **argv, const char *prefix, const char **rev_ret)
{
	const char *rev = "HEAD";
	unsigned char unused[20];
	/*
	 * Possible arguments are:
	 *
	 * git reset [-opts] <rev> <paths>...
	 * git reset [-opts] <rev> -- <paths>...
	 * git reset [-opts] -- <paths>...
	 * git reset [-opts] <paths>...
	 *
	 * At this point, argv points immediately after [-opts].
	 */

	if (argv[0]) {
		if (!strcmp(argv[0], "--")) {
			argv++; /* reset to HEAD, possibly with paths */
		} else if (argv[1] && !strcmp(argv[1], "--")) {
			rev = argv[0];
			argv += 2;
		}
		/*
		 * Otherwise, argv[0] could be either <rev> or <paths> and
		 * has to be unambiguous.
		 */
		else if (!get_sha1_committish(argv[0], unused)) {
			/*
			 * Ok, argv[0] looks like a rev; it should not
			 * be a filename.
			 */
			verify_non_filename(prefix, argv[0]);
			rev = *argv++;
		} else {
			/* Otherwise we treat this as a filename */
			verify_filename(prefix, argv[0], 1);
		}
	}
	*rev_ret = rev;
	return argv[0] ? get_pathspec(prefix, argv) : NULL;
}

static int update_refs(const char *rev, const unsigned char *sha1)
{
	int update_ref_status;
	struct strbuf msg = STRBUF_INIT;
	unsigned char *orig = NULL, sha1_orig[20],
		*old_orig = NULL, sha1_old_orig[20];

	if (!get_sha1("ORIG_HEAD", sha1_old_orig))
		old_orig = sha1_old_orig;
	if (!get_sha1("HEAD", sha1_orig)) {
		orig = sha1_orig;
		set_reflog_message(&msg, "updating ORIG_HEAD", NULL);
		update_ref(msg.buf, "ORIG_HEAD", orig, old_orig, 0, MSG_ON_ERR);
	} else if (old_orig)
		delete_ref("ORIG_HEAD", old_orig, 0);
	set_reflog_message(&msg, "updating HEAD", rev);
	update_ref_status = update_ref(msg.buf, "HEAD", sha1, orig, 0, MSG_ON_ERR);
	strbuf_release(&msg);
	return update_ref_status;
}

int cmd_reset(int argc, const char **argv, const char *prefix)
{
	int reset_type = NONE, update_ref_status = 0, quiet = 0;
	int patch_mode = 0;
	const char *rev;
	unsigned char sha1[20];
	const char **pathspec = NULL;
	struct commit *commit;
	const struct option options[] = {
		OPT__QUIET(&quiet, N_("be quiet, only report errors")),
		OPT_SET_INT(0, "mixed", &reset_type,
						N_("reset HEAD and index"), MIXED),
		OPT_SET_INT(0, "soft", &reset_type, N_("reset only HEAD"), SOFT),
		OPT_SET_INT(0, "hard", &reset_type,
				N_("reset HEAD, index and working tree"), HARD),
		OPT_SET_INT(0, "merge", &reset_type,
				N_("reset HEAD, index and working tree"), MERGE),
		OPT_SET_INT(0, "keep", &reset_type,
				N_("reset HEAD but keep local changes"), KEEP),
		OPT_BOOLEAN('p', "patch", &patch_mode, N_("select hunks interactively")),
		OPT_END()
	};

	git_config(git_default_config, NULL);

	argc = parse_options(argc, argv, prefix, options, git_reset_usage,
						PARSE_OPT_KEEP_DASHDASH);
	pathspec = parse_args(argv, prefix, &rev);

	if (get_sha1_committish(rev, sha1))
		die(_("Failed to resolve '%s' as a valid ref."), rev);

	/*
	 * NOTE: As "git reset $treeish -- $path" should be usable on
	 * any tree-ish, this is not strictly correct. We are not
	 * moving the HEAD to any commit; we are merely resetting the
	 * entries in the index to that of a treeish.
	 */
	commit = lookup_commit_reference(sha1);
	if (!commit)
		die(_("Could not parse object '%s'."), rev);
	hashcpy(sha1, commit->object.sha1);

	if (patch_mode) {
		if (reset_type != NONE)
			die(_("--patch is incompatible with --{hard,mixed,soft}"));
		return run_add_interactive(rev, "--patch=reset", pathspec);
	}

	/* git reset tree [--] paths... can be used to
	 * load chosen paths from the tree into the index without
	 * affecting the working tree nor HEAD. */
	if (pathspec) {
		if (reset_type == MIXED)
			warning(_("--mixed with paths is deprecated; use 'git reset -- <paths>' instead."));
		else if (reset_type != NONE)
			die(_("Cannot do %s reset with paths."),
					_(reset_type_names[reset_type]));
	}
	if (reset_type == NONE)
		reset_type = MIXED; /* by default */

	if (reset_type != SOFT && reset_type != MIXED)
		setup_work_tree();

	if (reset_type == MIXED && is_bare_repository())
		die(_("%s reset is not allowed in a bare repository"),
		    _(reset_type_names[reset_type]));

	if (pathspec)
		return read_from_tree(pathspec, sha1,
				quiet ? REFRESH_QUIET : REFRESH_IN_PORCELAIN);

	/* Soft reset does not touch the index file nor the working tree
	 * at all, but requires them in a good order.  Other resets reset
	 * the index file to the tree object we are switching to. */
	if (reset_type == SOFT || reset_type == KEEP)
		die_if_unmerged_cache(reset_type);

	if (reset_type != SOFT) {
		struct lock_file *lock = xcalloc(1, sizeof(struct lock_file));
		int newfd = hold_locked_index(lock, 1);
		int err = reset_index(sha1, reset_type, quiet);
		if (reset_type == KEEP && !err)
			err = reset_index(sha1, MIXED, quiet);
		if (err)
			die(_("Could not reset index file to revision '%s'."), rev);
		if (write_cache(newfd, active_cache, active_nr) ||
		    commit_locked_index(lock))
			die(_("Could not write new index file."));
	}

	/* Any resets update HEAD to the head being switched to,
	 * saving the previous head in ORIG_HEAD before. */
	update_ref_status = update_refs(rev, sha1);

	if (reset_type == HARD && !update_ref_status && !quiet)
		print_new_head_line(commit);
	else if (reset_type == MIXED) /* Report what has not been updated. */
		update_index_refresh(0, NULL,
				quiet ? REFRESH_QUIET : REFRESH_IN_PORCELAIN);

	remove_branch_state();

	return update_ref_status;
}
