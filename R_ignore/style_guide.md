# `fomes` - style guide for contributors



## Package structure

The overall package structure follows the standard structure of any R package, with the following idiosyncrasies:

- The "R_ignore" folder is setup to be ignored by R, but not by Github. This is where we can store files that are not part of a standard R package, but that we want included anyway, e.g. papers.
- If you make an "ignore" folder in your local version of the package then this will be ignored by *both* R and github. This is where you can store local files that you may want to be inside the package for organisational reasons, but not make it into anyone elses version, e.g. test scripts.
- Note, for both of the above folders, even though R claims to ignore them, in fact it copies and then deletes them when the package is built. Hence, avoid putting anything very large (e.g. whole genomes) in these folders as it will slow things down massively.


## Code structure and style

Contributors are asked to stick to the existing coding style, which follows Hadley Wickham's suggested style (http://r-pkgs.had.co.nz/r.html) with the following tweaks:

- single-line if statements, e.g. if (x < 10) y <- 5, are never OK. Always use {} over multiple lines.

Code should be structured and written in a way that is easy to pick up at a later date with minimal effort. Ease of reading takes precendence over code efficiency in most cases, as coding time is usually more costly than running time.

All R functions should be documented using the roxygen method (http://r-pkgs.had.co.nz/man.html). Both R and C++ functions should be documented clearly throughout.

All code should have comments that are helpful justifications, reminders, or explanations. Please do not provide comments that describe what the code is doing, as this should be self-evident.



<br>
<br>
### Sources
Thanks to Bob Verity & SIMPLEGEN for this style guide template
