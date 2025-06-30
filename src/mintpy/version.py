############################################################
# Program is part of MintPy                                #
# Copyright (c) 2013, Zhang Yunjun, Heresh Fattahi         #
# Author: Zhang Yunjun, Mar 2018                           #
############################################################
# Recommend import:
#   from mintpy.version import version, version_description, logo

import collections
import os
import subprocess
from importlib.metadata import PackageNotFoundError, metadata

###########################################################################
Tag = collections.namedtuple('Tag', 'version date')
release_history = (
    Tag('1.6.2', '2025-07-07'),
    Tag('1.6.1', '2024-07-31'),
    Tag('1.6.0', '2024-05-09'),
    Tag('1.5.3', '2023-11-23'),
    Tag('1.5.2', '2023-08-09'),
    Tag('1.5.1', '2023-01-03'),
    Tag('1.5.0', '2022-11-18'),
    Tag('1.4.1', '2022-08-15'),
    Tag('1.4.0', '2022-08-04'),
    Tag('1.3.3', '2022-04-14'),
    Tag('1.3.2', '2021-11-21'),
    Tag('1.3.1', '2021-08-02'),
    Tag('1.3.0', '2021-03-06'),
    Tag('1.2.3', '2020-07-14'),
    Tag('1.2.2', '2020-04-27'),
    Tag('1.2.1', '2020-03-30'),
    Tag('1.2.0', '2020-01-03'),
    Tag('1.2beta', '2019-08-08'),
    Tag('1.2alpha', '2019-07-28'),
    Tag('1.1.2', '2019-05-14'),
    Tag('1.1.1', '2019-05-08'),
    Tag('1.1.0', '2019-04-24'),
    Tag('1.0.0', '2018-11-07'),
    Tag('0.4.0', '2018-03-11'),
    Tag('0.3.0', '2017-06-03'),
    Tag('0.2.1', '2017-03-11'),
    Tag('0.2.0', '2016-07-14'),
    Tag('0.1.0', '2015-11-23'),
)

def get_version_info():
    """Grab version and its date from a git repository."""
    # go to the repository directory
    dir_orig = os.getcwd()
    os.chdir(os.path.dirname(os.path.dirname(__file__)))

    try:
        # grab from git cmd
        cmd = "git describe --tags"
        version = subprocess.check_output(cmd.split(), stderr=subprocess.DEVNULL)
        version = version.decode('utf-8').strip()[1:]
        # if there are new commits after the latest release
        if '-' in version:
            version, num_commit = version.split('-')[:2]
            version += f'.post{num_commit}'
        cmd = "git log -1 --date=short --format=%cd"
        version_date = subprocess.check_output(cmd.split(), stderr=subprocess.DEVNULL)
        version_date = version_date.decode('utf-8').strip()
    except:
        # use the latest release version/date
        version = release_history[0].version
        version_date = release_history[0].date

    # go back to the original directory
    os.chdir(dir_orig)
    return version, version_date


def get_version_info_v2():
    """Grab the version and its date from a git repository.
    Note by Yunjun on 5 Dec 2023: could not grab release/commit date info
        using importlib, thus, return None always.
    """
    try:
        package_name = os.path.basename(os.path.dirname(__file__))
        package = metadata(package_name)
        version = package['Version']
    except PackageNotFoundError:
        print('package is not installed!\n'
          'Please follow the installation instructions in the README.md.\n'
          'Or, to just get the version number, use:\n'
          '   python -m setuptools_scm')
        version = None
    return version, None


###########################################################################
version, version_date = get_version_info()
version_description = """MintPy version {v}, date {d}""".format(
    v=version,
    d=version_date,
)

# generate_from: http://patorjk.com/software/taag/
logo = r"""
___________________________________________________________

  /##      /## /##             /##     /#######
 | ###    /###|__/            | ##    | ##__  ##
 | ####  /#### /## /#######  /######  | ##  \ ## /##   /##
 | ## ##/## ##| ##| ##__  ##|_  ##_/  | #######/| ##  | ##
 | ##  ###| ##| ##| ##  \ ##  | ##    | ##____/ | ##  | ##
 | ##\  # | ##| ##| ##  | ##  | ## /##| ##      | ##  | ##
 | ## \/  | ##| ##| ##  | ##  |  ####/| ##      |  #######
 |__/     |__/|__/|__/  |__/   \___/  |__/       \____  ##
                                                 /##  | ##
                                                |  ######/
   Miami InSAR Time-series software in Python    \______/
          MintPy {v}, {d}
___________________________________________________________
""".format(v=version, d=version_date)
