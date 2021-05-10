Hierarchy of sub-modules within MintPy. Level _N_ modules depends on level _N-1_, _N-2_, ..., 0 modules. Be cautious of introducing circular imports.

```
/mintpy
------------------ level 0 --------------------
    /defaults
        auto_path
        template
    /objects
        cluster
        colors
        giant
        ramp
        sensor
    /utils
        /solvers
            l1
            l1regls
            lstl1
        arg_group
        ptime
        utils0
    /simulation
        fractal
------------------ level 1 --------------------
    /objects
        conncomp      (objects/ramp)
        stack         (utils/ptime)
    /simulation
        decorrelation (utils/ptime)
        defo_model    (utils/utils0)
        variance      (utils/ptime)
------------------ level 2 --------------------
    /utils
        readfile      (objects/{stack, giant})
------------------ level 3 --------------------
    /objects
        resample      (utils/{readfile, utils0, ptime})
        coord         (utils/{readfile, utils0, utils1})
    /utils
        writefile     (objects/{stack},         utils/{readfile})
        network       (objects/{stack, sensor}, utils/{readfile})
------------------ level 4 --------------------
    /objects
        gps           (objects/{stack, coord},  utils/{ptime, utils1, readfile})
        stackDict     (objects/{stack},         utils/{ptime, utils0, readfile})
    /simulation
        simulation    (objects/{stack},         utils/{ptime, network}, simulation/{fractal, decorrelation, defo_model})
    /utils
        attribute     (objects/{coord},         utils/{readfile})
        utils1        (objects/{stack, ramp},   utils/{ptime, utils0, readfile, writefile})
------------------ level 5 --------------------
    /utils
        plot          (objects/{stack, coord, colors},  utils/{ptime, utils0, readfile, network})
        utils         (objects/{stack, coord},          utils/{ptime, utils1, readfile})
        isce_utils    (utils/{ptime, readfile, writefile, utils1})
------------------ level 6 --------------------
    /objects
        insar_vs_gps  (objects/{stack, giant},          utils/{readfile, gps, plot, utils})
```
