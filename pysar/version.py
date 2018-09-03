release_version = '1.0.0-dev'
release_date = '2018-08-02'

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

version_description="""PySAR  version {v}  release date {d}
Copyright (C) 2013-2018 by Zhang Yunjun, Heresh Fattahi and others
Website: {w}

PySAR comes with ABSOLUTELY NO WARRANTY. This is free software, and
you are welcome to redistribute it under certain conditions. See the
GNU General Public License v3.0 for details.""".format(v=release_version,
                                                       d=release_date,
                                                       w=website)
