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


def create_rdf_node(node_content):
    """Creates and returns an RDF node

    :param symbol: node_content EITHER a string or a tuple (namespace_uri, local_name)"""
    if node_content is None:
        return None
    elif isinstance(node_content, rdflib.Literal) or isinstance(node_content, rdflib.BNode) or \
            isinstance(node_content, rdflib.URIRef):
        return node_content

    if not isinstance(node_content, str) and len(node_content) == 2:
        local_name = node_content[1]
        uri = node_content[0]
        # Ensure namespace prefix can be appended to
        if not uri.endswith('#') and not uri.endswith('/'):
            uri += '#'
        return rdflib.Namespace(uri)[local_name]
    elif isinstance(node_content, str) and node_content.startswith('#'):
        return rdflib.URIRef(node_content)
    else:
        raise NotImplementedError("Cannot create node")
