# Contributing Guidelines

This document is inspired by similar instructions from GMT, ISCE, gdal and jupyterhub. If you're reading this section, you're probably interested in contributing to MintPy. Welcome and thanks for your interest in contributing! 

These are several ways to contribute to the MintPy project:

* Submitting bug reports, feature requests on [GitHub issue](https://github.com/insarlab/MintPy/issues)
* Writing/improving documentation, tutorials and jupyter-notebooks
* Fixing typos, bugs in code
* Writing code for everyone to use

If you get stuck at any point you can create an [issue on GitHub](https://github.com/insarlab/MintPy/issues) or contact us on the [user forum](https://groups.google.com/forum/#!forum/mintpy).

For more information on contributing to open source projects, [GitHub's own guide](https://guides.github.com/activities/contributing-to-open-source/)
is a great starting point if you are new to version control. Tute Costa has a great [tutorial on how to rewrite history with git rebase/squash/amend](https://thoughtbot.com/blog/git-interactive-rebase-squash-amend-rewriting-history).

## Writing documentations ##

Documentation is written in Markdown. ANY GitHub user can edit pages and/or create new pages on [GitHub Wiki](https://github.com/insarlab/MintPy/wiki) directly for documentation, examples, support, or anything you wish. Eventually, mature documents on the wiki will be moved to the `insarlab/MintPy/docs` through standard code review process to be shown in the [readthedocs](https://mintpy.readthedocs.io/en/latest/).

## Writing code ##

We follow the [git pull request workflow](https://www.asmeurer.com/git-workflow/) to make changes to our codebase. Every change made goes through a pull request, even our own, so that our [continuous integration](https://en.wikipedia.org/wiki/Continuous_integration) services have a change to check that the code is up to standards. This way, the master branch is always stable.

### General guidelines for pull requests (PRs) ###

+ **Open an issue first** describing what you want to do, except for bugs fix. If there is already an issue that matches your PR, leave a comment there instead to let us know what you plan to do. We may have easier ways to help you implement it faster. 
+ Each pull request should consist of a **small** and logical collection of changes.
+ Larger changes should be broken down into smaller components and integrated separately.
+ Bug fixes should be submitted in separate PRs.
+ Describe what your PR changes and why this is a good thing. Be as specific as you can. The PR description is how we keep track of the changes made to the project over time.
+ Do not commit changes to files that are irrelevant to your feature or bugfix (eg: `.gitignore`, IDE project files, etc).
+ Write descriptive commit messages. Chris Beams has a [guide](https://chris.beams.io/posts/git-commit/) on how to write good commit messages. Tute Costa has a great [tutorial](https://thoughtbot.com/blog/git-interactive-rebase-squash-amend-rewriting-history) on how to rewrite a nice and clean history with git rebase/squash/amend.
+ Be willing to accept criticism and work on improving your code; we don't want to break other users' code, so care must be taken not to introduce bugs.
+ Be aware that the pull request review process is not immediate, and is generally proportional to the size of the pull request.

### Code Review ###

After you've submitted a pull request, you should expect to hear at least a comment within a couple of days. We may suggest some changes or improvements or alternatives.

Some things that will increase the chance that your pull request is accepted quickly:

+ Write a good and detailed description of what the PR does.
+ Readable code is better than clever code (even with comments).
+ Write documentation for your code and leave comments explaining the _reason_ behind non-obvious things.
+ Include an example of new features in the gallery or tutorials, if possible.

Pull requests will automatically have tests run by Circle CI and Codacy. Github will show the status of these checks on the pull request. Try to get them all passing (green). If you have any trouble, leave a comment in the PR or contact us on the [user forum](https://groups.google.com/forum/#!forum/mintpy).

### Example workflow ###

This is not a git tutorial by any means. It just collects a few best practice for git usage for MintPy development. There are plenty of good resources online to help get started.

#### Set up a development environment ####

Fork insarlab/MintPy from GitHub UI, and then

```
git clone https://github.com/my_user_name/MintPy.git
cd MintPy
git remote add upstream https://github.com/insarlab/MintPy.git
```

#### Work with a feature branch ####

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

# you may need to re-synchronize against upstream/master
# if you need some bugfix or new capability that has been
# added to master since you created your branch
git fetch upstream
git rebase upstream/master

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

#### Issue a pull request from GitHub UI ####

If the pull request discussion results in changes, commit new changes to `my_user_name/my_new_feature_branch`, they will show up in the pull request in `insarlab` automatically.


## Testing ##

It's a good idea to test any changes or bugs you have fixed, in the feature branch before issuing the pull request. We realize that we don't have a complete testing system in place yet (maybe you can contribute this!), except for an overall testing script `test_smallbaselineApp.py`, run

```
${MINTPY_HOME}/test/test_smallbaselineApp.py
```

to see the testing result, it takes about 10 mins to finish.


## Things you should NOT do ##

(For anyone with push rights to github.com/insarlab/MintPy) Never modify a commit or the history of anything that has been committed to the `master` branch.
