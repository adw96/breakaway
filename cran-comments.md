This is an update to breakaway to patch an issue occurring on mac M1 systems – because, we believe, of slightly
different linear algebra implementations on M1 systems, optim() was returning slightly different results and 
causing tests to fail. We have patched this issue by improving the robustness of our optimizations via
multiple random initializations.

 We have run checks locally, on win-builder, and on r-hub,
all of which pass either with no notes or with a note we believe is irrelevant to this package's fitness for CRAN.

Specifically, checks on win-builder and r-hub return a note that breakaway was previously archived on CRAN. This note also lists "microbiome" as a 
possibly misspelled word in breakaway documentation, and a (working) link is flagged as possibly not working. 

Thank you!