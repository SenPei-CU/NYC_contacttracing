# Statistical analysis of community transmission of COVID-19 in New York City

Code and data for the statistical analysis in Pei et al. Contact tracing reveals community transmission of COVID-19 in New York City, Nature Communications (2022).

The code performs Poisson generalized linear mixed models (GLMM) with random effects and uses conditional autoregressive (CAR) priors to account for the inherent spatial-temporal autocorrelation in disease transmission data. For introduction of this method, please see: Lee, D., Rushworth, A. & Napier, G. Spatio-Temporal Areal Unit Modeling in R with Conditional Autoregressive Priors Using the CARBayesST Package. J. Stat. Softw. 84, 1–39 (2018). The code also runs Moran's I test and the Durbin-Watson test on residues. 

Details of the variables in the data table:
zip: zip code
week: week number starting from the week of 10/04/2020
localtrans: number of non-household within-zip code transmission events in this week
localtrans_1: number of non-household within-zip code transmission events 1 week after the current week
localtrans_2: number of non-household within-zip code transmission events 2 weeks after the current week
crosstrans: number of non-household cross-zip code transmission events in this week
crosstrans_1: number of non-household cross-zip code transmission events 1 week after the current week
crosstrans_2: number of non-household cross-zip code transmission events 2 weeks after the current week
logpop: log transformed population in this zip code area
logpopdensity: log transformed population density in this zip code area
logcase: log transformed weekly case numbers per capita in this zip code area
logtest: log transformed weekly test numbers per capita in this zip code area
logcumucaserate: log transformed cumulative reported cases per capita until the current week
black: percentage of Black residents
hispanic: percentage of Hispanic residents
householdincome: median household income
bachelor: percentage of residents with bachelor degree
householdsize: mean household size
totalvisitor: weekly number of POI visitors per capita in the current week
vaccoverage: percentage of fully vaccinated residents in the current week
age65plus: percentage of residents with age 65+

