"""
    add-pubmed.py -- given identifers in a csv file, add the corresponding
    papers to VIVO

    Version 0.0 MC 2014-05-12
    --  Getting started.  This code has never been run.

    To Do
    --  check tyhe date code.  Which is a date string and which is a date?
    --  check the data and resouerce proerty names for pubs
    --  add the date dictionary code to vivotools and use it from there
    --  add author list processing (argh!)
    --  run for the first time
    --  write a bunch of test cases
    --  straighten out document and rdf and add.  One per functionn, please.
"""

__author__ = "Michael Conlon"
__copyright__ = "Copyright 2014, University of Florida"
__license__ = "BSD 3-Clause license"
__version__ = "0.0"

from datetime import datetime
from vivotools import rdf_header
from vivotools import rdf_footer
from vivotools import make_concept_dictionary
from vivotools import update_pubmed
from vivotools import read_csv
from vivotools import get_vivo_uri
from vivotools import find_vivo_uri
from vivotools import document_from_pubmed
from Bio import Entrez
import vivotools as vt

class NotFound(Exception):
    pass

class TimeOut(Exception):
    pass

def get_entrez_record(pmid):
    """
    Given a pmid, use Entrez to get first record from PubMed
    """
    Entrez.email = 'mconlon@ufl.edu'

    # Get record(s) from Entrez.  Retry if Entrez does not respond

    start = 2.0
    retries = 10
    count = 0
    while True:
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
            records = Entrez.parse(handle)
            break
        except IOError:
            count = count + 1
            if count > retries:
                raise Timeout
            sleep_seconds = start**count
            print "<!-- Failed Entrez query. Count = "+str(count)+ \
                " Will sleep now for "+str(sleep_seconds)+ \
                " seconds and retry -->"
            time.sleep(sleep_seconds) # increase the wait time with each retry
    return records[0]

def make_pub_rdf(pub):
    """
    Given a pub structure, make VIVO RDF for the publication
    """
    ardf = ""
    pub_uri = get_vivo_uri()
    [add, sub] = assert_resource_property(pub_uri, "rdf:type", 
        translate_predicate("owl:Thing"))
    ardf = ardf + add
    [add, sub] = assert_resource_property(grant_uri, "rdf:type",
        translate_predicate("bibo:AcademicArticle")
    ardf = ardf + add
    properties = {'title':'rdfs:label',
                  'volume':'bibo:volume',
                  'issue':'bibo:issue',
                  'pmid':'vivo:pmid',
                  'doi':'vivo:doi',
                  'page_start':'bibo:startPage',
                  'page_end':'bibo:endPage',
                  'date_harvested':'ufVivo:dateHarvested',
                  'harvested_by':'ufVivo:harvestedBy'}
    resources = {'journal_uri':'vivo:journal',
                 'date_uri':'vivo:dateTimeValue'}
    for property in properties.keys():
        if property in pub:
            add = assert_data_property(grant_uri,
                properties[property],
                pub[property])
    for resource in resource.keys():
        if resource in pub:
            add = assert_resource_property(grant_uri,
                resources[resource],
                pub[resource])
    return [ardf, pub_uri]

def add_pubmed(pmid):
    """
    Given a pubmid identifer, return a structure containing the elements
    of the publication of interest to VIVO
    """
    ardf = ""
    record = get_entrez_record(pmid)
    if record is None:
        return ["", None]
    pub = document_from_pubmed(record)
    pub['date_harvested'] = str(datetime.now())
    pub['harvested_by'] = "Python PubMed Add "+__version__
    pub['journal_uri'] = find_vivo_uri("bibo:issn", pub['issn'])

    pub_date = pub['date']['month']+'/'+pub['date']['day']+'/'+
                                      pub['date']['year']
    if pub_date in date_dictionary:
        pub['date_uri'] = date_dictionary[pub_date]
    else:
        [add, pub_date_uri] = make_datetime_rdf(pub_date.isoformat())
        date_dictionary[pub_date] = pub_date_uri
        pub['date_uri'] = pub_date_uri
        ardf = ardf + add
    
    add = make_pub_rdf(pub_uri, pub)
    ardf = ardf + add
    return [ardf, pub_uri]

#
#  Start Here
#

srdf = rdf_header()
ardf = rdf_header()

print datetime.now, "Start"
log_file = open("add_pubs_log.txt", "w")
exc_file = open("add_pubs_exc.txt", "w")
print datetime.now(),"Making concept dictionary"
make_concept_dictionary()
print datetime.now(),"Concept dictionary has ",\
    len(vt.concept_dictionary),"entries"
add_pubs = read_csv('add_pubs.csv')

for n, row in add_pubs.items():
    pmid = row['pmid']
    author_uri = find_vivo_uri('ufVivo:ufid', row['ufid'])
    if author_uri is None:
        print >>exc_file, row['ufid'], "UFID not found in VIVO"
        continue
    pub_uri = find_vivo_uri('bibo:pmid', pmid)
    if pub_uri is not None:
        print >>exc_file, "PubMed ID", pmid, "found in VIVO at", pub_uri
        continue

    try:
        [add, pub_uri] = add_pubmed(pmid)
        ardf = ardf + add
    except NotFound:
        print >>exc_file, "PubMed ID", pmid, "not found in PubMed"
        continue
    except TimeOut:
        print >>exc_file, "Entrez timed out for ", pmid
        continue
    else:
        print >>exc_file, "Unknown error for PubMed ID", pmid
        continue

    try:
        [add,sub] = update_pubmed(pub_uri)
        ardf = ardf + add
        srdf = srdf + sub
    except:
        print "Exception in update_pubmed for", pub_uri
        continue

    print >>log_file, "Paper added for ", author_uri, "at", pub_uri


add_file = open("add.rdf","w")
ardf = ardf + rdf_footer()
print >>add_file,ardf
add_file.close()

sub_file = open("sub.rdf","w")
srdf = srdf + rdf_footer()
print >>sub_file,srdf
sub_file.close()

log_file.close()
exc_file.close()
print datetime.now(), "Finish"

