# Release 0.2.1
- Improved support for secondary trigonometric functions such as sec and acoth (PR#314).

# Release 0.2.0
- Speed of parsing CellML models has been improved by reorganising connection processing. (PR#303)
- Removed `Model.connect_variables` as it is no longer needed by the parser. (PR#303)
- It is now possible to pass `sort=False` to the `get_state_variables`, `get_derivatives`, `get_derived_quantities` and `get_variables_by_rdf` methods of `Model` to get unsorted results. This gives a small performance improvement when you are sorting the result of those calls differently or when the order is not important. (PR#302)

# Release 0.1.0
Initial release of cellmlmanip.
