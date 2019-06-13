# Contributing Guidelines

This document is inspired by similar instructions from ISCE, gdal and jupyterhub. If you're reading this section, you're probably interested in contributing to MintPy. Welcome and thanks for your interest in contributing! 

These are several ways to contribute to the MintPy project:

* Submitting bug reports and feature requests
* Writing tutorials or jupyter-notebooks
* Fixing typos, code and improving documentation
* Writing code for everyone to use

If you get stuck at any point you can create an [issue on GitHub](https://github.com/insarlab/MintPy/issues) or contact us on the [user forum](https://groups.google.com/forum/#!forum/mintpy).

For more information on contributing to open source projects, [GitHub's own guide](https://guides.github.com/activities/contributing-to-open-source/)
is a great starting point if you are new to version control.

## Writing documentations ##

Documentation is written in Markdown on [GitHub Wiki](https://github.com/insarlab/MintPy/wiki). Any GitHub user can create and edit pages to use for documentation, examples, support, or anything you wish.

## Writing code ##

### 0. Setting up a development environment ###

Fork insarlab/MintPy from GitHub UI, and then

```
git clone https://github.com/my_user_name/MintPy.git
cd MintPy
git remote add upstream https://github.com/insarlab/MintPy.git
```

### 1. Working with a feature branch ###

[Here](https://thoughtbot.com/blog/git-interactive-rebase-squash-amend-rewriting-history) is a great tutorial if you are new to rewriting history with git.

```
# update to the latest upstream master
git checkout master
git fetch upstream
git rebase upstream/master
git push -f
git checkout -b my_new_feature_branch

# do work. For example:
git add my_new_file
git add my_modifid_message
git rm old_file
git commit -a 

# At end of your work, make sure history is reasonable by:
# folding non significant commits into a consistent set
# use 'fixup' for example to merge several commits together
# use 'reword' to modify commit messages
# to re-write the last 5 commits for example:
git rebase -i HEAD~5

# push your local changes to your fork on GitHub
git push

# you may need to force-push your branch with
git push -f
```

### 2. Issue a pull request from GitHub UI ###

If the pull request discussion results in changes, commit new changes to `my_user_name/my_new_feature_branch`, they will show up in the pull request in `insarlab` automatically.

## Testing ##

It's a good idea to test any changes or bugs you have fixed. We realize that we don't have a complete testing system in place yet, except for an overall testing script `test_smallbaselineApp.py`, run

```
${MINTPY_HOME}/test/test_smallbaselineApp.py
```

to see the testing result, it takes about 26 mins to finish.


## Things you should NOT do ##

(For anyone with push rights to github.com/insarlab/MintPy) Never modify a commit or the history of anything that has been committed to the `master` branch.
