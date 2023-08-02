# Release 0.3.6
- Added the ability to exclude a collection of tags from being used as variable names, by adding the `exclude_terms` parameter to  `model.get_display_name`. This can be used for example to prevent ontology annotation that can be used multiple times from appearing as the display name of a variable.

# Release 0.3.5
- Pinned Pint version to be > 0.9 and < 0.20 since 0.20 has changed the location and number of arguments of the ScaleConverter and UnitDefinition classes
- Fixed a bug in the parser where equations in a piecewise containing a boolean caused parsing errors.
  see https://github.com/ModellingWebLab/cellmlmanip/issues/350
- Added an error for duplicate unit definitions.
- Added error message when trying to connect components that do not exist.
- Added an error for duplicate component names.
- Fixed errors dealing with dimensionless units which have a multiplier or exponent, and added an error for offset units (where the offset is not 0) as those are not supported.
  see https://github.com/ModellingWebLab/cellmlmanip/issues/351
- Improved error reporting: When a unit is defined inside a component and used it now throws a ValueError indicating units in components are not supported. Previously a confusing KeyError was thrown, which suggested the unit was not found. 
- Added a better error message for inequalities <, >, <=, >= when one or both sides are a boolean. Added a warning for == and != when one side is a boolean, but not the other.
- Added better error messages for booleans appaearing in derivatives and derivative equations.

# Release 0.3.4
- Updated how substitution of functions that were changed in the parser are handled during analysis for fixing singularities, in order to make sure it workes with sympy 1.10
- Dropped support for python 3.5 as it is end of life.

# Release 0.3.3
- Minor performance upgrade for `Model.remove_fixable_singularities` using caching on fixing singularites. No functionality changes.

# Release 0.3.2
- Bug fix for parsing error created by version 6 of rdflib. This fix is backwards compatible with version 5 and does not change anything else.
 
# Release 0.3.1
- Added a better error message for unsupported unit celsius.

# Release 0.3.0
- Added a method `Model.remove_fixable_singularities` to remove fixable singularities in the model's equations.

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
