# ebola-cfr

Code for case fatality risk estimation in real-time.

### Quick start guide

First, set local path in R to GitHub directory, e.g.:
`
setwd("~/Documents/GitHub/ebola-cfr/")
`
The main script to reproduce the analysis is in `scripts/main_script.R`. This calls the following R files:

> `R/data_load.R` - Script to load and format data from 1976 outbreak and early stages of 2014 epidemic

> `R/cfr_function.R` - Function to estimate real-time CFR and 95% CI, adjusting for delay from onset to death

> `R/plot_examples.R` - Functions to plot CFR estimation for 1976 outbreak and early stages of 2014 epidemic

Note these functions use the onset-to-death distribution from the 2014 Ebola epidemic, rather than 1976 as used in the below Lancet piece. Incidence data from 2014 were extrapolated from cumulative estimates at non-even intervals using a simple averaging function, so the specific daily incidence values should be interpreted with caution.

### Citations

[Kucharski AJ, Edmunds WJ. Case fatality rate for Ebola virus disease in west Africa. Lancet, 2014](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(14)61706-2/fulltext)

[Nishiura et al. Early epidemiological assessment of the virulence of emerging infectious diseases: a case study of an influenza pandemic. PLoS One, 2009](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0006852)
