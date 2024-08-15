# Contributing Guidelines

This document is inspired by similar instructions from numpy, GMT, ISCE, gdal and jupyterhub. If you're reading this section, you're probably interested in contributing to MintPy. Welcome and thanks for your interest in contributing!

There are several ways to contribute to the MintPy project:

* Writing or proofreading documentation (including tutorials and examples in jupyter-notebooks)
* Submitting bug reports or feature requests on [GitHub issue](https://github.com/insarlab/MintPy/issues)
* Suggesting or implementing tests
* Fixing typos or bugs in code
* Writing code for everyone to use
* Giving feedback about the projects (including giving feedback about the contribution process)

If you get stuck at any point you can open an [issue on GitHub](https://github.com/insarlab/MintPy/issues) or comment on any open issue or pull request or contact us on the [user forum](https://groups.google.com/forum/#!forum/mintpy).

For more information on contributing to open source projects, [GitHub's own guide](https://guides.github.com/activities/contributing-to-open-source/)
is a great starting point if you are new to version control.

## Development process ##

#### 1. If you are a first-time contributor: ####

+ Go to https://github.com/insarlab/MintPy.git and click "fork" button to create your own copy of the project.

+ Clone the project to your local computer and add the upstream repository:

   ```
   git clone https://github.com/your_user_name/MintPy.git
   cd MintPy
   git remote add upstream https://github.com/insarlab/MintPy.git
   ```

+ Now, `git remote -v` will show two remote repositories named:

   - `upstream`, which refers to the `insarlab` repository
   - `origin`, which refers to your personal fork

+ Setting up [`pre-commit`](https://pre-commit.com/) within `MintPy` directory:
   - Run `pre-commit install` to set up the git hook scripts, so that `pre-commit` will run automatically on `git commit`. If the `No .pre-commit-config.yaml file was found` error occurs, update your local MintPy to the latest upstream version to have this config file.


#### 2. Develop your contribution: ####

+ **Open an [issue](https://github.com/insarlab/MintPy/issues) first** if you plan to introduce a new feature or to change functionality, we may have easier ways to help you implement it. If there is already an issue that matches your idea, leave a comment there instead to let us know what you plan to do. For bug fixes, documentation updates, etc., this is generally not necessary.

+ Pull the latest changes from upstream:

   ```
   git checkout main
   git pull upstream/main
   ```

+ Create a branch for the feature you want to work on. Since the branch name will appear in the merge message, use a sensible name such as 'seasonal_fitting':

   ```
   git checkout -b seasonal_fitting
   ```

+ Work on your idea, run tests, and commit locally (`git add` and `git commit`) and/or to your fork on GitHub as you progress (`git push` in the command line or [GitHub Desktop](https://desktop.github.com/) with graphical user interface). Use a clear commit message describing the motivation for a change, the nature of a bug for bug fixes, or some details on what an enhancement does.

+ Run the following `pre-commit` commands repeatedly until all checks are passed:

   ```
   pre-commit run --all-files
   ```

+ Run the [overall test](./CONTRIBUTING.md#testing) locally.

#### 3. To submit your contribution: ####

+ Go to your fork on GitHub. The new branch will show up with a Pull Request button. Click and fill out the pull request template, and make sure the title and message are clear, concise, and self-explanatory (these descriptions are how we keep track of the changes made to the project over time). Then click the button to submit it.

#### 4. Review process: ####

We follow the [git pull request (PR) workflow](https://www.asmeurer.com/git-workflow/) to make changes to our codebase. Every change made goes through a PR, even our own, so that our [continuous integration](https://en.wikipedia.org/wiki/Continuous_integration) services have a chance to check that the code is up to standards. GitHub will show the status of these checks on the PR. Try to get them all passing (green). If you have any trouble, leave a comment in the PR.

+ Reviewers (the other developers and interested community members) will write inline and/or general comments on your PR to help you improve its implementation, documentation and style. We don't want to break the shared codebase, so care must be taken not to introduce bugs. Please don’t let the review discourage you from contributing: its only aim is to improve the quality of the project, not to criticize (we are, after all, very grateful for the time you’re donating!).

+ To update your PR, make your changes on your local repository, run tests, and only if they succeed commit and push to your fork. As soon as those changes are pushed up (to the same branch as before) the PR will update automatically. If you have no idea how to fix the test failures, you may push your changes anyway and ask for help in a PR comment.


## Divergence between `upstream/main` and your feature branch ##

If GitHub indicates that the branch of your Pull Request can no longer be merged automatically, you have to incorporate changes that have been made since you started into your branch. Our recommended way to do this is to rebase on `main`. Tute Costa has a great tutorial on how to [rewrite history with git rebase/squash/amend](https://thoughtbot.com/blog/git-interactive-rebase-squash-amend-rewriting-history).


## General guidelines for pull requests (PRs) ##

+ Each pull request should consist of a **small** and logical collection of changes.
+ Larger changes should be broken down into smaller components and integrated separately.
+ Describe what your PR changes and why this is a good thing. Be as specific as you can.
+ Do not commit changes to files that are irrelevant to your feature or bugfix (_e.g._: `.gitignore`, IDE project files, etc.).
+ Write descriptive commit messages. Chris Beams has a [guide on writing good commit messages](https://chris.beams.io/posts/git-commit/).
+ Be aware that the pull request review process is not immediate, and is generally proportional to the size of the pull request.

Some things that will increase the chance that your pull request is accepted quickly:

+ Write a good and detailed description of what the PR does.
+ Readable code is better than clever code (even with comments).
+ Write documentation for your code and leave comments explaining the _reason_ behind non-obvious things.


## Testing ##

It's a good idea to test any changes or bugs you have fixed, in the feature branch locally, before issuing/submitting the pull request. We realize that we don't have a complete testing system in place yet (maybe you can contribute this!), except for an overall testing script for `smallbaselineApp.py`:

```
${MINTPY_HOME}/tests/smallbaselineApp.py
```

It takes about 15 minutes to finish.


## Things you should NOT do ##

(For anyone with push rights to github.com/insarlab/MintPy) Never modify a commit or the history of anything that has been committed to the `main` branch.

## Looking for ideas to contribute? ##

Feel ready and look for something to try? Check our [to-do list](https://github.com/insarlab/MintPy/wiki/To-do-list) and [roadmap](https://github.com/insarlab/MintPy/wiki/version-2-roadmap) for the next major version 2.0!
