Hierarchy of sub-modules within MintPy. Level _N_ modules depends on level _N-1_, _N-2_, ..., 0 modules. Be cautious of introducing circular imports.

```
/mintpy
------------------ level 0 --------------------
    constants
    /defaults
        auto_path
        template
    /objects
        cluster
        colors
        euler_pole
        giant
        ramp
        sensor
    /stdproc
        multilook
    /utils
        /solvers
            l1
            l1regls
            lstl1
        arg_utils
        map
        ptime (including progress)
        utils0
    /simulation
        fractal
------------------ level 1 --------------------
    /objects
        conncomp      (objects/ramp)
        ionex         (objects/ptime, utils/map)
        stack         (utils/ptime)
    /utils
        time_func     (utils/ptime)
    /simulation
        decorrelation (utils/ptime)
        defo_model    (utils/utils0)
        variance      (utils/ptime)
------------------ level 2 --------------------
    /utils
        readfile      (objects/{stack, giant, sensor}, utils/{ptime, utils0})
        s1_utils      (objects/{stack}, utils/{ptime, time_func})
------------------ level 3 --------------------
    /objects
        resample      (utils/{utils0, ptime, readfile}, constants)
        coord         (constants, utils/{utils0, utils1, readfile})
    /simulation
        iono          (utils/{utils0, readfile})
    /utils
        writefile     (utils/{readfile})
        network       (objects/{stack, sensor}, utils/{readfile})
------------------ level 4 --------------------
    /objects
        gnss          (objects/{stack, coord}, utils/{ptime, utils1, readfile})
        stackDict     (objects/{stack}, utils/{ptime, utils0, readfile})
    /simulation
        simulation    (objects/{stack}, utils/{ptime, network}, simulation/{fractal, decorrelation, defo_model})
    /utils
        attribute     (objects/{coord}, utils/{readfile})
        utils1        (objects/{stack, ramp}, utils/{ptime, utils0, readfile, writefile})
------------------ level 5 --------------------
    /utils
        plot          (objects/{stack, coord, colors}, utils/{ptime, utils0, readfile, network, map})
        utils         (objects/{stack, coord, resample}, utils/{ptime, attribute, utils1, readfile})
        isce_utils    (constants, utils/{ptime, readfile, writefile, attribute, utils1})
------------------ level 6 --------------------
    /objects
        insar_vs_gnss (objects/{stack, giant}, utils/{readfile, gnss, plot, utils})
```
