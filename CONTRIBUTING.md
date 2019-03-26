# Contributing Guidelines #  

This document is inspired by similar instructions from ISCE, gdal and jupyterhub. If you're reading this section, you're probably interested in contributing to PySAR. Welcome and thanks for your interest in contributing! 

These are several ways to contribute to the PySAR project:

  * Submitting bug reports and feature requests
  * Writing tutorials or jupyter-notebooks
  * Fixing typos, code and improving documentation
  * Writing code for everyone to use

If you get stuck at any point you can create an [issue on GitHub](https://github.com/insarlab/PySAR/issues) or contact us on the [user forum](https://groups.google.com/forum/?nomobile=true#!forum/py-sar).

For more information on contributing to open source projects, [GitHub's own guide](https://guides.github.com/activities/contributing-to-open-source/)
is a great starting point if you are new to version control.

## Writing documentations ##  

Documentation is written in Markdown on [GitHub Wiki](https://github.com/insarlab/PySAR/wiki). Any GitHub user can create and edit pages to use for documentation, examples, support, or anything you wish.

## Git workflows ##

### Setting up a development environment ##

Fork insarlab/PySAR from GitHub UI, and then

```
git clone https://github.com/insarlab/PySAR.git
cd PySAR
git remote add my_user_name https://github.com/my_user_name/PySAR.git
```

### Updating your local master against upstream master ###

```
git checkout master
git fetch origin
# Be careful: this will loose all local changes you might have done now
git reset --hard origin/master
```

### Working with a feature branch ###

[Here](https://thoughtbot.com/blog/git-interactive-rebase-squash-amend-rewriting-history) is a great tutorial if you are new to rewriting history with git.

```
git checkout master
(potentially update your local master against upstream, as described above)
git checkout -b my_new_feature_branch

# do work. For example:
git add my_new_file
git add my_modifid_message
git rm old_file
git commit -a 

# you may need to resynchronize against master if you need some bugfix
# or new capability that has been added to master since you created your
# branch
git fetch origin
git rebase origin/master

# At end of your work, make sure history is reasonable by folding non
# significant commits into a consistent set
git rebase -i master (use 'fixup' for example to merge several commits together,
and 'reword' to modify commit messages)

# or alternatively, in case there is a big number of commits and marking
# all them as 'fixup' is tedious
git fetch origin
git rebase origin/master
git reset --soft origin/master
git commit -a -m "Put here the synthetic commit message"

# push your branch
git push my_user_name my_new_feature_branch
```

From GitHub UI, issue a pull request

If the pull request discussion results in changes,
commit locally and push. To get a reasonable history, you may need to
```
git rebase -i master
```
, in which case you will have to force-push your branch with 
```
git push -f my_user_name my_new_feature_branch
```

## Testing ##

It's a good idea to test any changes or bugs you have fixed. We realize that we don't have a complete testing system in place yet, except for an overall testing script `test_pysarApp.py`, run
```bash
${PYSAR_HOME}/test/test_pysarApp.py
```
to see the testing result, it takes about 10 mins to finish.

## Things you should NOT do ##

(For anyone with push rights to github.com/insarlab/PySAR) Never modify a commit or the history of anything that has been committed to the `master` branch.
