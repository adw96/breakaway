This is a followup to an update to breakaway to patch an issue occurring on mac M1 systems – because, we believe, of slightly
different linear algebra implementations on M1 systems, optim() was returning slightly different results and 
causing tests to fail. We have patched this issue by improving the robustness of our optimizations via
multiple random initializations.

In feedback to our update, the following issues were raised:

- The description field of the breakaway DESCRIPTION would be improved with references to relevant published work
- Only undirected quotation marks should be used in the DESCRIPTION
- Some .Rd files for exported functions lacked \value entries
- Some functions adjusted user options via par() 

We have, to the best of our knowledge, addressed each of these issues in this submission. The
description field of the breakaway DESCRIPTION file now includes links to papers whose methods
appear in breakaway, and no directed quotation marks are used in the DESCRIPTION. We have 
updated documentation to ensure that all exported functions have adequate documentation in .Rd
files, and additionally this submission includes vignettes providing further documentation for 
breakaway. We have amended all functions changing user options with par() so that options
are restored on exit (with the on.exit() function, as suggested).

We have run checks locally, on win-builder, and on r-hub,
all of which pass either with no notes or with a note we believe is irrelevant to this package's fitness for CRAN.

Specifically, checks on win-builder and r-hub return a note that breakaway was previously archived on CRAN. This note also lists "microbiome" as a 
possibly misspelled word in breakaway documentation, and a (working) link is flagged as possibly not working. 

Thank you!