# grab version / date of the latest commit
import os
import subprocess
_dir_orig = os.getcwd()
os.chdir(os.path.dirname(os.path.dirname(__file__)))
try:
    release_version = subprocess.check_output(["git", "describe", "--tags"]).decode('utf-8').strip()
    release_date = subprocess.check_output(["git", "log", "-1", "--date=short", "--format=%cd"]).decode('utf-8').strip()
except:
    release_version = 'v1.1.0-dev'
    release_date = '2019-04-26'
os.chdir(_dir_orig)

# generate_from: http://patorjk.com/software/taag/
logo = """
_________________________________________________
       ____             __     __     ____  
       /    )         /    )   / |    /    )
------/____/----------\-------/__|---/___ /------
     /        /   /    \     /   |  /    |  
____/________(___/_(____/___/____|_/_____|_______
                /                           
            (_ /                            

 A Python package for InSAR time series analysis.
           PySAR v{v}, {d}
_________________________________________________
""".format(v=release_version, d=release_date)

website = 'https://yunjunz.github.io/PySAR/'

description = """PySAR version {v}, release date {d}""".format(v=release_version,
                                                               d=release_date)
