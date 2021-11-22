#!/usr/bin/env python3
# grab version / date of the latest commit


import os
import subprocess
import collections


###########################################################################
Tag = collections.namedtuple('Tag', 'version date')
release_history = (
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
release_version = release_history[0].version
release_date = release_history[0].date

def get_version_info(version='v{}'.format(release_version), date=release_date):
    """Grab version and date of the latest commit from a git repository"""
    # go to the repository directory
    dir_orig = os.getcwd()
    os.chdir(os.path.dirname(os.path.dirname(__file__)))

    # grab git info into string
    try:
        cmd = "git describe --tags"
        version = subprocess.check_output(cmd.split(), stderr=subprocess.DEVNULL)
        version = version.decode('utf-8').strip()

        #if there are new commits after the latest release
        if '-' in version:
            version, num_commit = version.split('-')[:2]
            version += '-{}'.format(num_commit)

        cmd = "git log -1 --date=short --format=%cd"
        date = subprocess.check_output(cmd.split(), stderr=subprocess.DEVNULL)
        date = date.decode('utf-8').strip()
    except:
        pass

    # go back to the original directory
    os.chdir(dir_orig)
    return version, date


###########################################################################

version_num, version_date = get_version_info()
version_description = """MintPy version {v}, date {d}""".format(
    v=version_num,
    d=version_date,
)

# generate_from: http://patorjk.com/software/taag/
logo = """
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
""".format(v=version_num, d=version_date)

website = 'https://github.com/insarlab/MintPy'

description = 'Miami INsar Time-series software in PYthon'

