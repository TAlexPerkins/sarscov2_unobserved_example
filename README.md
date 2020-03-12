# sarscov2_unobserved_example

This repository contains code to reproduce results from the following thread on Twitter.

https://twitter.com/TAlexPerkins/status/1237513140578865153

A preprint describing the methods used to generate these results is still under development. This page will be updated to a link to the preprint as soon as it is available.

In the meantime, **I urge anyone using this code to do so with extreme caution**. This work has not been peer-reviewed. Furthermore, changing parameter values or other settings in ways that are biologically unrealistic (despite the user's best intentions) is very likely to generate misleading results.

The primary reason that I am making this available in advance of a preprint is because several researchers have requested this code so they can use it to rapidly perform similar analyses for other countries where SARS-CoV-2 is spreading. They will need to do some work to adapt this to their country. Unfortunately, I cannot provide support for that at this time other than to say that the main changes required involve updating data inputs on lines 12-19, 54, and 66 in code/script.R. In a nutshell, if you run script.R in R from the code folder, that will generate all of the results in the figures folder.

Please note that the results in the figures folder appear slightly different from what was on the Twitter thread. That is because we have been refining some of our assumptions and parameter settings, as this is still work in progress. Nonetheless, the conclusions from the Twitter thread still represent our best professional assessment of this situation at this time.

# License

This code is being released under the MIT License.