# Release 0.2.0
- Speed of parsing CellML models has been improved by reorganising connection processing. (PR#303)
- Removed `Model.connect_variables` as it is no longer needed by the parser.
- It is now possible to pass `sort=False` to the `get_state_variables`, `get_derivatives`, `get_derived_quantities` and `get_variables_by_rdf methods` of `Model` to get unsorted results. This gives a small performance improvement when you are already sorting the result of those calls or when the order is not important.

# Release 0.1.0
Initial release of cellmlmanip
