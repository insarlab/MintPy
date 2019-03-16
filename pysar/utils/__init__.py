## Sub-module dependency graph:
## Check pysar/__init__.py for more comprehensiave graph
# Level 0 modules:
# /solvers
#     l1
#     l1regls
#     lstl1
# ptime
# utils0
# Level 1 modules:
# variance  (utils/ptime)
# readfile  (objects/*)
# writefile (objects/*, utils/readfile)
# network   (objects/*, utils/readfile)
# utils1    (objects/*, utils/writefile)
#
