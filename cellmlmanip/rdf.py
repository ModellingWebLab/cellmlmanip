"""Module for working with RDF data."""
import rdflib


def create_rdf_node(node_content):
    """Creates and returns an RDF node.

    :param node_content: the conetnet for the node.

    The ``node_content``, if given, must either be a :class:`rdflib.term.Node` instance,
    a tuple ``(namespace, local_name)``, or a string, in which case it is interpreted
    as either a URI ref if it starts with # otherwise a literal RDF node.

    The ``node_content`` may also be `None` to return `None`, allowing easy handling of wildcard options to queries."""
    if node_content is None:
        return None
    elif isinstance(node_content, rdflib.term.Node) or isinstance(node_content, rdflib.term.Node) or \
            isinstance(node_content, rdflib.URIRef):
        return node_content

    if isinstance(node_content, tuple):
        uri, local_name = node_content
        # Ensure namespace prefix can be appended to
        if not uri.endswith('#') and not uri.endswith('/'):
            uri += '#'
        return rdflib.Namespace(uri)[local_name]
    elif isinstance(node_content, str) and node_content.startswith('#'):
        return rdflib.URIRef(node_content)
    else:
        return rdflib.Literal(node_content)
