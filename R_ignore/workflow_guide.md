# `fomes` - style guide for contributors

## Continuous integration

This package is setup for testing and continuous integration with Github Workflows. If you look at the main Github page (https://github.com/TriCEM/fomes) you should see badges showing that the current build is passing. It is crucial that certain branches pass all checks at all times (see below for details).


## Github and branching

This package uses git and Github for version control. This has the advantage that it is impossible to completely break the code as we always have saves going back through time.

**We will use the git-flow branching** pattern described here (https://nvie.com/posts/a-successful-git-branching-model/), please read this before working on the code. In simple terms:

- The master branch is the *outward-facing stable branch*. It should always represent the most recent official release, and therefore should never break. You should never work directly on this branch.
- The develop branch is the *inward-facing stable branch*. Similar to master, this branch should always pass all checks. The difference is that develop will continually move forward as new features are merged in, whereas master is frozen in time at official release points. You should never work directly on this branch.
- All large code changes should be done through feature branches. These can be small and short-lived, or more extensive. You are free to break and fix code as much as you like on feature branches. Once you are happy, these should be merged into develop using a Github pull request - not directly in the console. As long as the PR passes all checks you are free to complete the merge, or you can nominate a reviewer if you would like someone else to look over your changes. Always going through PRs in this way ensures that develop remains stable. We will periodically delete old feature branches once merged. Feature branches should be named with the following structure: `feature/mynewfeature` (_e.g._ `feature/binominal_sampling`).
- All minor code changes should be done through hotfix branches. These should only be _very_ small, short-lived branches addressing minor issues that differ from feature branches above in their scope. We encourage frequent deletion of hotfix branches once merged. Hotfix branches should be named with the following structure: `hotfix/myquickfix` (_e.g._ `hotfix/forgot_to_add_test_for_binomial_feature`).
- Version numbers will be used to keep track of releases using the three-part X.Y.Z format, where X = major change, Y = small change e.g. added feature, Z = patch/bug fix. The development version will be 0.Y.Z, the first non-development release will be 1.0.0.


## C++

This package uses C++ through the Rcpp package. Unfortunately, however, debugging and profiling tools are not yet up to scratch for Rcpp. So, where possible within C++ code, hash-defines may be used to allow us to switch between different versions of the code that compile via different methods - for example within Rcpp vs. directly in Xcode on Mac.


<br>
<br>
### Sources
Thanks to Bob Verity & SIMPLEGEN for this workflow guide template
