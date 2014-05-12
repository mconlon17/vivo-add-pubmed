# Add Papers to VIVO from their PubMed identifiers

## Background

Updating papers in VIVO can be a chore.  The interface is cumbersome and slow.  Many biomedical researchers have
biosketches or curricula vitae that contain their publications.  The researcher often cites a PubMed identifier for each 
publication.  Adding publications to VIVO from a PubMed identifier will greatly simplify the ask of bringing
publication collections for biomedical faculty up to date.

## Goal

Provide a simple interface for adding papers in PubMed to VIVO.  Given a person identifier and a publication identifier,
add the publication for the person to VIVO.

## Method

1. Prepare a CSV file of UFID and PubMed identifier -- could be PubMed id, PMCID, NIHMSID, DOI (to be resolved in PubMed)
1. For each row:
    1. Find the person in VIVO from the UFID.  If not found, write error to the exception file and go to the next row.
    1. Get the PubMed ID from the supplied identifier.  If none, write error to exception file and go to the next row.
	1. Find the PubMed identifier in VIVO.  If found, write error to exception file and go to next row.
	1. Get the paper information from PubMed via Entrez.  If fail, write error to exception file and go to the next row.
	1. Identify the UF author(s) -- one of the authors must match the supplied UFID author.  If no UF authors are found, write error to the exception file and go to the next row.
	1. Create a python pub object from the Entrez data.
	1. Write a text version of the python pub object to the log.
    1. Generate the RDF or the python pub object and write to the Add file.
1. Wrap-up	
