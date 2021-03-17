# (unreleased)
- Added an automatic fix for removable singularities in GHK equations. This feature has been moved from chaste_codegen to allow more general use.
  `load_model` has gained an optional parameter `skip_singularity_fixes` to allow users to switch the fixing off (defaults to False).
  The process looks for equations of any of the following forms, where U is a function of V:
  - `U / (exp(U) - 1.0)`
  - `U / (1.0 - exp(U))`
  - `(exp(U) - 1.0) / U`
  - `(1.0 - exp(U)) / U`  
  It replaces these with a piecewise 1e-7 either side of U==0 drawing a stright line in the region.
  For example `(V + 5)/(exp(V + 5) - 1)` becomes `((fabs(-V - 5.0000000000000000) < fabs(-4.9999999000000000 / 2 - -5.0000001000000000 / 2)) ? ((-5.0000001000000000 + 5.0) / (-1.0 + exp(-5.0000001000000000 + 5.0)) + (--5.0000001000000000 + V) * ((-4.9999999000000000 + 5.0) / (-1.0 + exp(-4.9999999000000000 + 5.0)) - (-5.0000001000000000 + 5.0) / (-1.0 + exp(-5.0000001000000000 + 5.0))) / (--5.0000001000000000 - 4.9999999000000000)) : ((5.0 + V) / (-1.0 + exp(5.0 + V))))`

  See for more details appendix B in: Johnstone, R. H. (2018). Uncertainty characterisation in action potential modelling for cardiac drug safety [PhD thesis]. University of Oxford. https://ora.ox.ac.uk/objects/uuid:0a28829c-828d-4641-bfb0-11193ef47195

# Release 0.2.3
- Fixes for sympy 1.7: 
    - Fixed the printer test (eval=false no longer exists but evaluate=false can still be used). 
    - Removed the positional arguments in Variable as these were causing issues, instead use named arguments for name and units

# Release 0.2.2
- Fixed a rendering issue with secondary trigonometric functions such as sec and acoth (PR#317). This means for example that 1 / sec(x) now renders as cos(x) instead of 1 / 1 / cos(x).

# Release 0.2.1
- Improved support for secondary trigonometric functions such as sec and acoth (PR#314).

# Release 0.2.0
- Speed of parsing CellML models has been improved by reorganising connection processing. (PR#303)
- Removed `Model.connect_variables` as it is no longer needed by the parser. (PR#303)
- It is now possible to pass `sort=False` to the `get_state_variables`, `get_derivatives`, `get_derived_quantities` and `get_variables_by_rdf` methods of `Model` to get unsorted results. This gives a small performance improvement when you are sorting the result of those calls differently or when the order is not important. (PR#302)

# Release 0.1.0
Initial release of cellmlmanip.
