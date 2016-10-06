## Resubmission

Reduced size of tar.gz from 11M to 370K, as requested.

## Test environments
* local OS X install, R 3.3.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs, and 3 NOTES:

- License components with restrictions and base license permitting such:
  MIT + file LICENSE
File 'LICENSE':
  YEAR: 2014
  COPYRIGHT HOLDER: Stefan Schroedl

- running examples for arch 'i386' ... [19s]
Examples with CPU or elapsed time > 10s
          user system elapsed
plotluck 15.29   0.19   15.66
- running examples for arch 'x64' ... [26s]
Examples with CPU or elapsed time > 10s
          user system elapsed
plotluck 20.95   0.25   21.32

I have not changed the examples from last submission, and I don't get these
notes locally. I am assuming the win-builder server is busy.