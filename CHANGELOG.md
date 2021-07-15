# v19.1.2

  - Based off a meeting (15 July 2021)

## Changed

  - Changing unit of Breakthrough infection probability from RR back to percent.
    - value from template 0 ~ 100
    - `cmod`(`cm`) replaces `sigmaR`


# v19.1.1

  - Based off `covidage_v19.1.R` (12 July 2021)
  - Ref `covidage_v19.1_original.R -> covidage_v19.1.1.original.diff`

## New

 - `VOC` tab to hold "Transmissibility, Lethality, and Breakthrough infection probability"
 - `dmod_vector` capped at `1/max(pdeath_h,pdeath_ho,pdeath_hc,pdeath_hco,pdeath_icu,pdeath_icuo,pdeath_icuc,pdeath_icuco,pdeath_vent,pdeath_ventc,pdeath_vent_hc,pdeath_icu_hc,pdeath_icu_hco)`
 - `cmod_vector` capped at `1/sigmaR`

## Changed

 - Breakthrough infection probability : Unit changed from percentage (%) to Relative risk (RR)
   - `cmod` scales `sigmaR`

## Updated

 - In `covidage_v19.1.1` the post-processing part of the script is updated
    - `process_ode_outcome`
    - `multi_runs`


# v19.1.0

 - Based off `covidage_v19.1.R` (1 July 2021)
 - Objective is model variants
 - Ref `covidage_v19.1_original.R`

## New

 - SR compartment
 - QSR compartment
 - 3 x Interventions
   - Transmissibility
     - unit: RR 
     - `pmod`, `p_mod`, `pmod_vector`
     - `pm` scales `p`
   - Lethality
     - unit: RR
     - `dmod`, `d_mod`, `dmod_vector`
     - `dm` scales all `pdeath_` variables
   - Breakthrough infection probability
     - unit: %
     - `cmod`, `c_mod`, `cmod_vector`
     - `cmod` replaces `sigmaR`
 - RR: Relative risk
