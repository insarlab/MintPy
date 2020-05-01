Hierarchy of sub-modules within MintPy. Level N modules depends on level N-1, N-2, ..., 0 modules.

```
/mintpy
------------------ level 0 --------------------
    /defaults
        auto_path
        template
    /objects
        colors
        giant
        ramp
        sensor
        stack
    /utils
        /solvers
            l1
            l1regls
            lstl1
        ptime
        utils0
    /simulation
        fractal
------------------ level 1 --------------------
    /objects
        conncomp      (objects/ramp)
    /simulation
        decorrelation (utils/ptime)
        defo_model    (utils/utils0)
        variance      (utils/ptime)
    /utils
        readfile      (objects/{stack, giant})
------------------ level 2 --------------------
    /objects
        resample      (utils/{readfile, utils0, ptime})
        coord         (utils/{readfile, utils0, utils1})
    /utils
        writefile     (objects/{stack},         utils/{readfile})
        network       (objects/{stack, sensor}, utils/{readfile})
------------------ level 3 --------------------
    /objects
        gps           (objects/{stack, coord}, utils/{ptime, utils0, readfile})
        stackDict     (objects/{stack},        utils/{ptime, utils0, readfile})
    /simulation
        simulation    (objects/{stack},        utils/{ptime, network}, simulation/{fractal, decorrelation, defo_model})
    /utils
        utils1        (objects/{stack, ramp},  utils/{ptime, utils0, readfile, writefile})
------------------ level 4 --------------------
    /utils
        plot          (objects/{stack, coord, colors}, utils/{ptime, utils0, utils1, readfile, network})
        utils         (objects/{stack, coord},         utils/{ptime, utils0, utils1, readfile})
        isce_utils    (utils/{readfile, writefile, utils1})
------------------ level 5 --------------------
    /objects
        insar_vs_gps  (objects/{stack, giant}, utils/{readfile, gps, plot, utils})
```
