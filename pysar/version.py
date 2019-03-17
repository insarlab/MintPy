release_version = '1.0.0-dev'
release_date = '2019-03-16'

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
