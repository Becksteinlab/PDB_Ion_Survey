# PDB Ion Survey 
# Copyright (c) 2016 Kacey Clark
# Published under the GPL v3
# https://github.com/Becksteinlab/PDB_Ion_Survey/

'''
Functions for querying, downloading, and analyzing ion coordination of PDB structures
'''

from __future__ import absolute_import

import os.path
import sys
import io

import logging
import urllib.request
from collections import OrderedDict
from xml.sax.saxutils import XMLGenerator
from xml.sax.xmlreader import AttributesImpl


# this is how we identify an error message returned instead of a PDB file
error_message = 'RCSB web servers'


def _emit(key, value, content_handler, attr_prefix='@', cdata_key='#text', depth=0, preprocessor=None, pretty=False, newl='\n', indent='\t', full_document=True):
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
            or isinstance(value, str)
            or isinstance(value, dict)):
        value = [value]
    for index, v in enumerate(value):
        if full_document and depth == 0 and index > 0:
            raise ValueError('document with multiple roots')
        if v is None:
            v = OrderedDict()
        elif not isinstance(v, dict):
            v = str(v)
        if isinstance(v, str):
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


def unparse(input_dict, output=None, encoding='utf-8', full_document=True,
            **kwargs):
    '''Emit an XML document for the given `input_dict`.
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
    '''
    if full_document and len(input_dict) != 1:
        raise ValueError('Document must have exactly one root.')
    must_return = False
    if output is None:
        output = io.StringIO()
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


def get_pdb_ids(ionname, containsProtein=True, containsDNA=False, containsRNA=False, containsHybrid=False):
    '''Searches PDB for files with a specified bound ion.
    :Arguments:
        *ionname*
            String name of desired ion
        *containsProtein*
            Boolean true if including protein molecules in the search; default=True
        *containsDNA*
            Boolean true if including DNA molecules in the search; default=False
        *containsRNA*
            Boolean true if includingRNA molecules in the search; default=False
        *containsHybrid*
            Boolean true if including DNA/RNA hybrid molecules in the search; default=False
    :Returns:
        *ids*
            List ids of all PDB files satisfying input conditions
    -------------------------
    Based off of a function made by William Gilpin
    '''
    query_params = dict()
    query_params['queryType'] = 'org.pdb.query.simple.ChemCompIdQuery'
    query_params['chemCompId'] = ionname
    scan_params = dict()
    scan_params['orgPdbQuery'] = query_params
    url = 'http://www.rcsb.org/pdb/rest/search'
    queryText = unparse(scan_params, pretty = True)
    queryText = queryText.encode()
    req = urllib.request.Request(url, data = queryText)
    f = urllib.request.urlopen(req)
    result = f.read()
    idlist = str(result)
    idlist = idlist.split('\\n')
    query_paramsB = dict()
    query_paramsB['queryType'] = 'org.pdb.query.simple.ChainTypeQuery'
    if containsProtein:
        query_paramsB['containsProtein'] = 'Y'
    else:
        query_paramsB['containsProtein'] = 'N'
    if containsDNA:
        query_paramsB['containsDna'] = 'Y'
    else:
        query_paramsB['containsDna'] = 'N'
    if containsRNA:
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
    queryTextB = unparse(scan_paramsB, pretty=True)
    queryTextB = queryTextB.encode()
    reqB = urllib.request.Request(urlB, data=queryTextB)
    fB = urllib.request.urlopen(reqB)
    resultB = fB.read()
    idlistB = str(resultB)
    idlistB = idlistB.split('\\n')
    try:
        idlist.remove("'")
    except:
        pass
    try:
        idlistB.remove("'")
    except:
        pass
    idset = set(idlist)
    idsetB = set(idlistB)
    ids = idset.intersection(idsetB)
    ids = list(ids)
    return ids
