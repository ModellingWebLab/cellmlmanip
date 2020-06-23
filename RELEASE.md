# Release 0.1.0
It is now possible to pass sort=False to the `get_state_variables`, `get_derivatives`, `get_derived_quantities` and `get_variables_by_rdf` methods of `Model` to retrieve unsorted results. 
This is intended to imrpove performance in situations where you are already sorting result of these calls or the ordering is not important.

# Release 0.1.0
Initial release of cellmlmanip