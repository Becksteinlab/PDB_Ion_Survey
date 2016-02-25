"""
Functions for querying, downloading, and analyzing ion coordination of PDB structures
"""

import logging
import urllib2
import os.path
import sys
from collections import OrderedDict
from cStringIO import StringIO
from xml.sax.saxutils import XMLGenerator
from xml.sax.xmlreader import AttributesImpl


# this is how we identify an error message returned instead of a PDB file
error_message = 'RCSB web servers'


def _emit(key, value, content_handler, attr_prefix = '@', cdata_key = '#text', depth = 0, preprocessor = None, pretty = False, newl = '\n', indent = '\t', full_document=True):
    '''I do not know what this does.
    -------------------------
    Credit to: Martin Blech
    '''
    if preprocessor is not None:
        result = preprocessor(key, value)
        if result is None:
            return
        key, value = result
    if (not hasattr(value, '__iter__')
            or isinstance(value, basestring)
            or isinstance(value, dict)):
        value = [value]
    for index, v in enumerate(value):
        if full_document and depth == 0 and index > 0:
            raise ValueError('document with multiple roots')
        if v is None:
            v = OrderedDict()
        elif not isinstance(v, dict):
            v = unicode(v)
        if isinstance(v, basestring):
            v = OrderedDict(((cdata_key, v),))
        cdata = None
        attrs = OrderedDict()
        children = []
        for ik, iv in v.items():
            if ik == cdata_key:
                cdata = iv
                continue
            if ik.startswith(attr_prefix):
                attrs[ik[len(attr_prefix):]] = iv
                continue
            children.append((ik, iv))
        if pretty:
            content_handler.ignorableWhitespace(depth * indent)
        content_handler.startElement(key, AttributesImpl(attrs))
        if pretty and children:
            content_handler.ignorableWhitespace(newl)
        for child_key, child_value in children:
            _emit(child_key, child_value, content_handler, attr_prefix, cdata_key, depth+1, preprocessor, pretty, newl, indent)
        if cdata is not None:
            content_handler.characters(cdata)
        if pretty and children:
            content_handler.ignorableWhitespace(depth * indent)
        content_handler.endElement(key)
        if pretty and depth:
            content_handler.ignorableWhitespace(newl)


def unparse(input_dict, output = None, encoding = 'utf-8', full_document = True,
            **kwargs):
    """Emit an XML document for the given `input_dict`.
    The resulting XML document is returned as a string, but if `output` (a
    file-like object) is specified, it is written there instead.
    Dictionary keys prefixed with `attr_prefix` (default=`'@'`) are interpreted
    as XML node attributes, whereas keys equal to `cdata_key`
    (default=`'#text'`) are treated as character data.
    The `pretty` parameter (default=`False`) enables pretty-printing. In this
    mode, lines are terminated with `'\n'` and indented with `'\t'`, but this
    can be customized with the `newl` and `indent` parameters.
    -------------------------
    Credit to: Martin Blech
    """
    if full_document and len(input_dict) != 1:
        raise ValueError('Document must have exactly one root.')
    must_return = False
    if output is None:
        output = StringIO()
        must_return = True
    content_handler = XMLGenerator(output, encoding)
    if full_document:
        content_handler.startDocument()
    for key, value in input_dict.items():
        _emit(key, value, content_handler, full_document = full_document,
              **kwargs)
    if full_document:
        content_handler.endDocument()
    if must_return:
        value = output.getvalue()
        try:  # pragma no cover
            value = value.decode(encoding)
        except AttributeError:  # pragma no cover
            pass
        return value


def get_proteins(ionname, containsProtein = True, containsDna = False, containsRna = False, containsHybrid = False):
    """Searches PDB for files with a specified bound ion.
    :Arguments:
        *ionname*
            name of desired ion
        *containsProtein*
            boolean value of whether to include protein molecules in the search; default = True
        *containsDna*
            boolean value of whether to include DNA molecules in the search; default = False
        *containsRna*
            boolean value of whether to include RNA molecules in the search; default = False
        *containsHybrid*
            boolean value of whether to include DNA/RNA hybrid molecules in the search; default = False
    :Returns:
        *idlist*
            ids of all PDB files containing ions with name ionname

    -------------------------
    Based off of a function made by William Gilpin

    """
    query_params = dict()
    query_params['queryType'] = 'org.pdb.query.simple.ChemCompIdQuery'
    query_params['chemCompId'] = ionname
    scan_params = dict()
    scan_params['orgPdbQuery'] = query_params
    url = 'http://www.rcsb.org/pdb/rest/search'
    queryText = unparse(scan_params, pretty = True)
    queryText = queryText.encode()
    req = urllib2.Request(url, data = queryText)
    f = urllib2.urlopen(req)
    result = f.read()
    idlist = str(result)
    idlist = idlist.split('\n')
    query_paramsB = dict()
    query_paramsB['queryType'] = 'org.pdb.query.simple.ChainTypeQuery'
    if containsProtein:
        query_paramsB['containsProtein'] = 'Y'
    else:
        query_paramsB['containsProtein'] = 'N'
    if containsDna:
        query_paramsB['containsDna'] = 'Y'
    else:
        query_paramsB['containsDna'] = 'N'
    if containsRna:
        query_paramsB['containsRna'] = 'Y'
    else:
        query_paramsB['containsRna'] = 'N'
    if containsHybrid:
        query_paramsB['containsHybrid'] = 'Y'
    else:
        query_paramsB['containsHybrid'] = 'N'
    scan_paramsB = dict()
    scan_paramsB['orgPdbQuery'] = query_paramsB
    urlB = 'http://www.rcsb.org/pdb/rest/search'
    queryTextB = unparse(scan_paramsB, pretty = True)
    queryTextB = queryTextB.encode()
    reqB = urllib2.Request(urlB, data=queryTextB)
    fB = urllib2.urlopen(reqB)
    resultB = fB.read()
    idlistB = str(resultB)
    idlistB = idlistB.split('\n')
    idset = set(idlist)
    idsetB = set(idlistB)
    ids = idset.intersection(idsetB)
    ids = list(ids)
    return ids


def _make_logger(filename):
    """Make a logger instance that prints to the screen and writes
    to `filename`.

    """
    log = logging.getLogger('PDB_Ions')
    log.setLevel(logging.WARNING)

    if not any([isinstance(x, logging.StreamHandler) for x in log.handlers]):
        ch = logging.StreamHandler(sys.stdout)
        cf = logging.Formatter(
            '%(name)-12s: %(levelname)-8s %(message)s')
        ch.setFormatter(cf)
        log.addHandler(ch)

    if not any([isinstance(x, logging.FileHandler) for x in log.handlers]):
        fh = logging.FileHandler(filename)
        ff = logging.Formatter('%(asctime)s %(name)-12s '
                               '%(levelname)-8s %(message)s')
        fh.setFormatter(ff)
        log.addHandler(fh)

    return log


def get_pdb_file(pdb_id, path, compression = False):
    '''Get the full PDB file associated with a PDB_ID
    :Arguments:
        *pdb_id*
            four character string giving a pdb entry of interest
        *path*
            path to the desired location of the file
        *compression
            retrieve a compressed (gz) version of the file if True; default = False
    :Returns:
        *pdb_id.pdb*
            pdb file corresponding to pdb_id
    -------------------------
    Based off of a function made by William Gilpin
    '''
    # get logger instance
    log = _make_logger(os.path.join(path, 'pdb_ions_survey.log'))

    fullurl = 'http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb'
    if compression:
        fullurl += '&compression=YES'
    else:
        fullurl += '&compression=NO'
    fullurl += '&structureId=' + pdb_id

    # if file exists, don't download again 
    if os.path.isfile(os.path.join(path, pdb_id) + '.pdb'):
        return

    # creates a string containing pdb file information
    req = urllib2.Request(fullurl)
    f = urllib2.urlopen(req)
    result = f.read()
    result = result.decode('unicode_escape')

    # tests result is not error
    if error_message in result:
        log.warning('RCSB PDB does not have a .pdb file available for' 
                    'download for PDB ID {}'.format(pdb_id))
        raise KeyError("PDB ID '{}' does not have a"
                       ".pdb file available".format(pdb_id))

    f_out = open(os.path.join(path, pdb_id) + '.pdb', 'w')
    f_out.write(result)
    f_out.close()
