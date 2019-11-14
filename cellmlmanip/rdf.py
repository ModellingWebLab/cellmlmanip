"""Module for working with RDF data."""
# import logging
import rdflib


# logger = logging.getLogger(__name__)
# logger.setLevel(logging.INFO)


'''
def create_rdf_node(self, node_content=None, fragment_id=None):
    """
    Creates an RDF node from the given ``node_content``.

    The ``node_content``, if given, must either be a :class:`rdflib.term.Node`
    instance, a tuple ``(namespace, local_name)``, or a string, in which case
    it is interpreted as a literal RDF node.

    Alternatively, ``fragment_id`` may be given to refer to a ``cmeta:id``
    within the current model.

    If neither are given, a blank node is created.
    """
    # This method was adapted from PyCml

    if fragment_id:
        if node_content is not None:
            raise ValueError(
                'Node content and fragment id cannot both be given.')
        node = rdflib.URIRef(str('#' + fragment_id))

    elif node_content is not None:

        if isinstance(node_content, rdflib.term.Node):
            # Use existing node
            node = node_content

        if isinstance(node_content, tuple):
            # Create a node from the given (namespace, local_name) tuple
            ns_uri, local_name = node_content

            # Ensure namespace prefix can be appended to
            if ns_uri[-1] not in ['#', '/']:
                ns_uri = ns_uri + '#'

            ns = rdflib.Namespace(ns_uri)
            node = ns[local_name]

        elif isinstance(node_content, str):
            # Create literal node
            node = rdflib.Literal(node_content)

        else:
            raise ValueError(
                'Unable to create a node from the given node_content:'
                + str(node_content))

    else:
        # Create blank node
        node = rdflib.BNode()

    return node
'''


def create_rdf_node(namespace_uri, local_name=None):
    """Creates and returns an RDF node ``{namespace_uri}local_name``."""
    if namespace_uri is None:
        return None

    if not isinstance(namespace_uri, str) and len(namespace_uri) == 2:
        local_name = namespace_uri[1]
        namespace_uri = namespace_uri[0]

    # Ensure namespace prefix can be appended to
    if namespace_uri[-1] not in ['#', '/'] and local_name is not None:
        namespace_uri = namespace_uri + '#'

    if local_name is None:
        return rdflib.Literal(namespace_uri)
    else:
        return rdflib.Namespace(namespace_uri)[local_name]
