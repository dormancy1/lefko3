# lefko3 6.5.0 (2025-05-02)

## USER VISIBLE CHANGES

* Functions projection3() and f_projection3() now operate more quickly, and use
  integeronly = TRUE by default.

* Main development site moved to GitHub, with links updated.

* Vignettes updated, including links, function output, and text errors.

## BUG FIXES

* Function miniMod() now correctly retains zero-inflation models.

* Empirical age-by-stage MPMs created with mpm_create() and arlefko2() no longer
  lead to fatal errors when a final age is used less than the maximum in the
  dataset.

* Function modelsearch() no longer fails on historical model runs with vital
  rate variable name vectors with 2 entries.
  
* Indexing of hstages element corrected for deVries-format historical MPMs.

# lefko3 6.4.0 (2025-03-13)

## NEW FEATURES

* Function cycle_check() created to test for stage discontinuities in matrices.

* Function-based matrix building procedures now test for NA values of model
  coefficients, and display warnings if found.
  
* Function supplemental() now allows additive offsets to transition elements in
  MPMs and projections.

## USER VISIBLE CHANGES

* Updated project URL.

* Function summary.lefkoMat() now also checks for stage discontinuities.

* Updated and standardized error messages.

* Vignette 1 updated with offset info.

## BUG FIXES

* Fixed error in projection3() resulting in shuffled starting frames under
  some conditions in which stageframes are shuffled from their original state.
  
* Fixed error in function f_projection3() coming from handling of tertiary
  individual covariate.

# lefko3 6.3.2 (2025-01-06)

## NEW FEATURES

* Function hfv_qc() now checks the number of singleton categories in all random
  variables.

## BUG FIXES

* Fixed incompatibility of lefkoMod objects created prior to version 6.3.0 with
  all matrix creation and projection functions in newer versions.

* Fixed incompatibility between integer vector inputs supplied to individual
  covariate arguments in function f_projection3().
  
* Fixed error causing incompatibility between ipm classification and certain
  combinations of minsize, maxsize, and halfbinwidths in sf_create().

# lefko3 6.3.1 (2024-08-20)

## NEW FEATURES

* Function hfv_qc() now checks the number of individuals with single transitions
  in hfv datasets.

## USER VISIBLE CHANGES

* Eliminated warnings for multiple similar transitions in mpm_create()
  when the transitions are of different types.

## BUG FIXES

* Corrected erroneous milepost messages in function modelsearch().

# lefko3 6.3.0 (2024-08-05)

## NEW FEATURES

* Function miniMod() was developed to create minimum memory vital rate model
  summaries in vrm_input format useable as substitutes for lefkoMod objects.

* Function modelsearch() can now test three annual covariates in addition to the
  three individual covariates.

* All function-based matrix estimation and projection functions can now handle
  values of three annual covariates.

* Function modelsearch() now models fixed factor interactions separately from
  size and reproductive status, allowing greater flexibility to choose global
  models.

## USER VISIBLE CHANGES

* Some help entries have been corrected and updated.

* Error messages standardized throughout much of the package.

* Function modelsearch() now outputs linear models with function calls that
  point to correctly named data subsets, all of which can be exported with the
  data_out argument.

## BUG FIXES

* Function matrix_interp() no longer gives a fatal error when run on a
  historical lefkoLTRE object with part != 2.

* Matrix creation and projection functions using vrm_input objects no longer
  produce Gaussian output when provided with negative binomial inputs.

* Function verticalize3() no longer produces errors when reproductive status,
  fecundity status, or observation status variables are not provided.

# lefko3 6.2.1 (2024-02-24)

## USER VISIBLE CHANGES

* Function matrix_interp() now also handles standard matrix and sparse matrix 
  (dgCMatrix format) inputs.

* Zero-inflated Poisson and negative binomial GLMs have been enabled again.

* Function overwrite() is now deprecated and links to function supplemental().

* Examples for functions involving function-based MPMs have been altered to
  showcase generalized linear models developed with base R.

* Function create_pm() now uses the same model parameter naming conventions as
  verticalize3() and historicalize3().

* Matrix creation and projection examples now use basic linear models rather
  than function modelsearch().
  
* Function modelsearch() now preferentially removes the most highly categorized
  random covariate when a global model fails, if using mixed modeling.
  
* Function hfv_qc() now also shows the numbers of categories in variables that
  can be used as random factors as a part of the output.
  
## BUG FIXES

* Zero-inflation estimation corrected to properly account for the mixture model
  in Poisson and negative binomial versions.

* Function summary.lefkoProj() no longer fails when summarizing appended
  projections with differing numbers of time steps.

* Function modelsearch() no longer improperly subsets data in Leslie MPM
  searches.
  
* Function modelsearch() now correctly removes independent terms to simplify
  model structure if global model development fails.
  
* Fixed structural bugs preventing matrix_interp() from working with functions
  sensitivity3() and elasticity3().

# lefko3 6.2.0 (2024-01-27)

## NEW FEATURES

* Functions projection3(), f_projection3(), slambda3(), stablestage3(),
  repvalue3(), sensitivity3(), elasticity3(), and ltre3() can now perform
  first-order Markovian temporally stochastic simulations.
  
* Functions stablestage3() and repvalue3() can now handle stochastic
  projections.
  
* Functions sensitivity3() and elasticity3() can now be forced to return output
  in either standard or sparse matrix format for single matrix inputs.

* Function append_lP() created to append projections to each other.

* Function matrix_interp() created to aid in summarizing results from analyses
  of huge matrices.
  
* Function lambda3() now estimates population growth rate and its logarithm in
  lefkoProj objects.

## USER VISIBLE CHANGES

* Some arguments have been standardized across population dynamics functions.

* Zero-inflated Poisson and negative binomial GLMs have been disabled in
  function modelsearch() pending issues in package pscl.
  
* Vignettes updated to deal with breaking changes.

## BUG FIXES

* Function f_projection3() no longer fails when non-stochastic projections are
  conducted with multiple times.
  
* Function f_projection3() no longer fails when the number of times requested is
  lower than the number of recognized times in the dataset.
  
* Functions stablestage3() and repvalue3() no longer fail when a tweights
  vector or matrix is entered.
  
* Functions sensitivity3() and elasticity3() no longer fail under some
  conditions involving sparse matrices.

# lefko3 6.1.3 (2023-11-30)

## NEW FEATURES

* Function markov_run() created to allow 1st order Markovian temporally
  stochastic simulations.
  
* Function create_lM() can now handle both regular and sparse matrices.

## USER VISIBLE CHANGES

* Some warning and error messages have been simplified.

## BUG FIXES

* Function projection3() no longer crashes when using a fixed vector of year
  inputs with lefkoMat objects holding patch-level MPMs.

# lefko3 6.1.2 (2023-11-07)

## BUG FIXES

* All MPMs can now be estimated without supplemental tables.

* Functions add_lM(), delete_lM(), and subset_lM() no longer fail with sparse
  matrices.

# lefko3 6.1.1 (2023-10-13)

## NEW FEATURES

* Function projection3() no longer adds an extra projection to matrix
  projections of mean MPMs created through the lmean() function.

## BUG FIXES

* Function density_input() no longer throws an error when using shorthand
  stage group designations with age-by-stage MPMs.
  
* Functions sensitivity3() and elasticity3() no longer produce zero division
  errors in some stochastic runs with sparse matrix output.
  
* Fixed indexing issue affecting propagation of multipliers for supplement
  tables used in function-based and raw Leslie MPMs.
  
* Fixed package documentation.

# lefko3 6.1.0 (2023-07-08)

## NEW FEATURES

* Function add_stage() has been created to add stages to existing lefkoMat
  objects.

## USER VISIBLE CHANGES

* Options repmod and fecmod were standardized to fecmod in MPM creation
  functions.
* Function hfv_qc() now only displays test results for vital rate models listed
  in the vitalrates argument.
* Functions rleslie() and fleslie() now include an entrystage column in the
  automatically generated stageframes included as ahstages object.

## BUG FIXES

* Function hfv_qc() can now handle Gaussian checks on variables with more than
  5000 data points.
* Function historicalize3() now correctly calculates first seen time, last seen
  time, observed age, and alive status in time t+1.
* Function supplemental() no longer yields an error for age-based MPM
  supplements with estage2 option set to all NAs.
* Function start_input() no longer yields errors for age-based MPMs.
* Function density_input() no longer yields errors for age-based MPMs.
* Historical matrix creation functions no longer call warnings that multiple
  rows in supplement tables show the same transition even when they do not.
* Function lambda3() no longer calls an error when an eigenvalue of 0 is
  detected.

# lefko3 6.0.5 (2023-05-03)

## USER VISIBLE CHANGES

* Function ltre3() now uses tolerance limits to determine non-zero elements
  and contributions in sLTRE and SNA-LTRE analyses.
* Stochastic LTREs are now more easily interruptible.

## BUG FIXES

* Functions summary.lefkoLTRE() and summary.lefkoElas() now properly sort
  contributions from fecundity transitions.

# lefko3 6.0.4 (2023-04-11)

## USER VISIBLE CHANGES

* Function ltre3() can now handle SNA-LTRE of large historical MPMs.

# lefko3 6.0.3 (2023-04-07)

## NEW FEATURES

* Function sf_create() now allows flexible assignment of representative stage
  sizes and associated minimum and maximum bin boundaries.

## USER VISIBLE CHANGES

* Function rlefko2() now produces an interpretable error message if an
  incorrectly developed stageframe is applied.
* Warnings added about use of large input matrices for SNA-LTRE analysis.

## BUG FIXES

* Fixed proxy transition issue in function edit_lM().

# lefko3 6.0.2 (2023-03-31)

## USER VISIBLE CHANGES

* Patch variable in labels object of dataset anthyllis has been corrected to
  match population names in the source paper.

# lefko3 6.0.1 (2023-03-30)

## NEW FEATURES

* Function create_lM() now includes an empty dataqc element in lefkoMat output.

## USER VISIBLE CHANGES

* Function edit_lM() no longer requires a dataqc element in input lefkoMat
  objects.
* Function summary.lefkoMat() correctly interprets lefkoMat objects missing a
  dataqc element.

# lefko3 6.0.0 (2023-03-27)

## NEW FEATURES

* Function summary_hfv() now includes an error checking function that searches
  hfv data frames for problems in stage assignment.
* New function mpm_create() is now the workhorse function handling all MPM
  creation. All previous matrix creation functions have been rewritten as
  wrapper functions for this new function, with arguments and handling as before
  to prevent compatibility issues.
* Function modelsearch() can now handle fixed factor individual covariates, and
  can use a different suite of independent factors for the global model of each
  vital rate.
* Function mpm_create() and all matrix creation functions now handle fixed
  factor variables as independent covariates, as well as random covariates
  and fixed quantitative variables.
* Function lmean() can produce arithmetic mean matrices for simple lists of
  square, equidimensional matrices.
* Functions f_projection3(), projection3(), stablestage3(), repvalue3(),
  sensitivity3(), elasticity3(), and ltre3() now all have override settings
  for sparse vs dense matrix encoding, allowing greater control over speed.
* Core matrix functions now allow sparse matrix output.
* Functions lmean(), lambda3(), slambda3(), stablestage3(), repvalue3(),
  sensitivity3(), elasticity3(), image3(), ltre3(), cond_hmpm(), cond_diff(),
  projection3(), summary.lefkoElas(), and summary.lefkoMat() can handle sparse
  matrix format MPMs.
* Functions stablestage3.list() and repvalue3.list() have been added.
* Function supplemental() now also handles age transitions in Leslie and
  age-by-stage MPMs, and in f_projection3() runs involving them.
* Function edit_lM() has been added to edit MPMs with external data.
* New dataset Pyrola hs been added.
* Function create_lM() now also allows direct imports from the COMPADRE and
  COMADRE databases as lefkoMat objects.

## USER-VISIBLE CHANGES

* Functions utilizing a choice of sparse vs dense matrix approaches have been
  sped up via the optimization of criteria used to determine which approach to
  take.
* Functions cond_hmpm() and cond_diff() now perform faster, with better memory
  management, and with corrected list structure under err_check mode.
* Cpp header files are now visible to other packages.
* Function hist_null() no longer attempts to remove impossible transitions.
* Overhead cut on projection3(), and related function that involve projecting
  matrices forward or backward.
* Function modelsearch() now throws an informative error if too many parameters
  are tested for exhaustive model dredging to operate.
* Function modelsearch() now simplifies certain zero-inflated models if they
  run past the limits of MuMIn's dredge() function.
* Substantial speed up to function-based matrix creation, particular when
  reducing matrices.
* All matrix creation functions now have the ability to output only U and F
  matrices ("simple" option), and to produce error-checking output ("err_check"
  option).
* Observation status variables may now be set in all raw MPM creation functions.
* NRasRep option added in raw matrix creation functions for cases where stages
  have not been assigned.
* Functions verticalize3() and historicalize3() now standardize datasets almost
  2x more quickly.
* Function modelsearch() now automatically subsets datasets down to complete
  cases of the variables used in the modeling to prevent failures.
* User-interrupt now enabled in matrix creation functions.
* Function ltre3() is now a generic function that handles lefkoMat, matrix, and
  matrix list inputs.
* Function lambda3() now handles a variety of input scenarios for the 'sparse'
  argument.
* Function overwrite() has been marked as deprecated.

## BUG FIXES

* Function hfv_qc() no longer gives an error when secondary or tertiary size
  variable names are provided.
* Functions rleslie() and fleslie() now properly create stageframe age IDs.
* Functions rlefko3() and rlefko2() now properly parse non-standard observation
  and maturity status variables for stage determination when an hfv data frame
  is provided without stages already determined.
* Corrected inaccurate age-stage indexing in object agestages resulting from
  function arlefko2().
* Corrected assignment of population matrix tags to the labels data frame in
  lefkoMat objects produced by lmean() under the "patch" option.
* Fixed fatal error when using summary.lefkoMat() with vrm_imported function-
  based MPMs.
* Fixed fatal error occurring when summary.lefkoElas() is used with a lefkoElas
  object developed from a list of matrices.
* Function sf_skeleton() added to create a basic, empty stageframe.
* Fixed bug yielding incorrect fecundity when prebreeding in fleslie().
* Eliminated bug causing rare glitches in fecunsdity estimation when stage order
  is altered by MPM creation functions.
* Fixed issue causing survival and fecundity estimation to ignore individuals
  older than max age in function rleslie() under continue = TRUE.
* Corrected issue in which survival and fecundity multipliers provided via
  supplemental() were only incorporated into function-based MPMs involving
  stages if the transitions involved proxy estimation.
* Corrected citation information.

# lefko3 5.5.0 (2022-09-14)

## NEW FEATURES

* Function modelsearch() now allows "partial" option, which shows model milepost
  messages but silences all others. Also allows "yes" and "no" in place of TRUE
  and FALSE.
* Function density_input() has been rewritten to handle all standard stage
  shorthand codes.

## USER-VISIBLE CHANGES

* NEWS file corrected to show version numbers.
* Function modelsearch() now uses a simple R-squared approach to assess
  accuracy of all models except binomial models, which is assessed as
  traditional logistic accuracy.
* Vignettes updated to reflect changes to accuracy calculation.
* End-of-line characters added to most warnings to prevent all warnings mixing
  together into single messages.

## BUG FIXES

* Corrected issue causing Linux-based matrix estimation to yield negative matrix
  element values in some instances involving truncated distributions.
* Fixed issue in density adjustment protocol used by function f_projection3()
  that would force all density adjustment onto survival probability.
* Corrected issue causing negative accuracy to be estimated for zero-truncated
  models in modelsearch().

# lefko3 5.4.2 (2022-08-08)

## USER-VISIBLE CHANGES

* NEWS file has been reformatted to conform to CRAN Markdown specifications.

# lefko3 5.4.1 (2022-07-29)

## USER-VISIBLE CHANGES

* Output to summary_hfv() now fits within screen width.

* Exponent tolerance limits and warnings added to functions projection3() and
  f_projection3(). Warnings and help file notes have also been added.

# lefko3 5.4.0 (2022-07-21)

## NEW FEATURES

* Function actualstage3() has been rewritten in C++ to handle all stage,
  age-by-stage, and age inputs. Can now remove stages as necessary, and allows
  different stages to be removed in times t1 and t2 in historical stage pair
  assessments.

* Function hist_null() can now produce historically-formatted ahistorical
  matrices that do not include impossible transitions.

* Anthyllis matrices are now also provided as a new dataset, in lefkoMat format.

## USER-VISIBLE CHANGES

* Updates and corrections to help files.

## BUG FIXES

* Lathyrus vignettes corrected to account for duplicated individual identity
  across subpopulations.

* Corrected issue in repvalue3() causing warnings in MPMs with certain rare
  conditions.

* Corrected issue in lambda3() causing rare eigen analysis failure in sparse
  matrices.

* Reference matrix warnings in ltre3() have been corrected.

# lefko3 5.3.0 (2022-06-24)

## NEW FEATURES

* Function actualstage3() now also outputs proportions of ages and age-stages
  within hfv datasets.

* Function vrm_import() created to allow IPM and other function-based model
  slope coefficients can now be imported directly without conducting linear
  model searches through modelsearch().

* Functions flefko3(), flefko2(), aflefko2(), fleslie(), and f_projection3()
  have been updated to allow the use of vrm_input objects, allowing slope
  coefficients to be imported directly.

* Function hfv_qc() developed to test hfv datasets for quality control using the
  same parameterizations and subsetting strategies used in function
  modelsearch().

## USER-VISIBLE CHANGES

* Function create_pm() can now yield paramnames objects with default model
  variable names used instead of "none" strings.

* Vignettes updated with new functions. Particularly, material for sf_distrib()
  replaced with material for hfv_qc().

* Tweaks and corrections to help files.

## BUG FIXES

* Fixed bug causing substochasticity setting 2 in functions projection3() and
  f_projection3() to yield negative survival probabilities under some
  circumstances in which survival column sums totaled more than 1.0.

* Fix bug omitting substochasticity correction in fecundity under some
  conditions when function f_projection3() is run with substoch = 2.

# lefko3 5.2.0 (2022-05-18)

## NEW FEATURES

* Function f_projection3() can now run density dependent projections in which
  the vital rates used are modified by projected density.

* Function density_vr() created to provide vital rate density dependence
  information for function f_projection3().

* Option added to functions projection3() and f_projection3() to suppress
  warnings.

* Function summary.lefkoProj() now has an option to interpret NaN values as
  positive infinity, changing the count of replicates still alive in the
  milepost section. This option, inf_alive, may also be set to FALSE to count
  NaNs as evidence of extinction.
  
* A new summary_hfv() function has been created to simplify quality control
  checks on historical vertical format datasets.

* Function verticalize3() now warns the user to check the resulting dataset if
  the dataset is censored and the blocksize option is used without the
  censorRepeat option being set to TRUE, as these cases may result in excessive
  data pruning.

* Function summary.lefkoProj() now calculates mean extinction time when
  prompted.

## USER-VISIBLE CHANGES

* Function modelsearch() now allows an "age" setting for the suite option,
  which yields the same results as "cons" but checks to make sure that an age
  variable is properly added.

* Function modelsearch() now allows a "leslie" setting for the vitalrates
  option, which provides a shorthand for c("surv", "fec") as would be used by
  fleslie().

* The documentation for functions modelsearch() and fleslie() have been
  updated to provide more information for handling Leslie MPMs.

* Added series of new error messages if expected data frame size input in
  verticalize3() is too large.

* Vignettes have been updated with new functions.

* Some error messages made more informative.

## BUG FIXES

* The stageframe developed with functions rleslie() and fleslie() has been fixed
  to include stage_id, alive, and almost_born variables, allowing
  Leslie-formatted lefkoMat objects to be passed to all lefko3 functions
  analyzing them, including lmean().

* Fixed bug yielding vector leaks from other variables when using the blocksize
  feature in verticalize3().

* Fixed incorrect fecundity summation in rlefko3(), resulting from miscounting
  offspring from individuals dying the transition in which they reproduced.

# lefko3 5.1.0 (2022-04-01)

## NEW FEATURES

* Added ability to run ordered annual matrix progressions in function
  projection3().

* Function projection3() can now run arithmetic mean matrix projections
  without any modifications required to the input lefkoMat object.

* Function plot.lefkoProj() now has an auto-title feature.

## USER-VISIBLE CHANGES

* Added cancel handling to verticalize3(), historicalize3(), slambda3(),
  sensitivity3(), elasticity3(), and ltre3().

* Added transitions to some supplemental() in function examples and vignettes,
  reflecting more exact transition specification.

## BUG FIXES

* Fixed incorrect substochastic correction in projection3() and f_projection3()
  when substoch = 2 and negative survival entries occur.

* Cleaned up minor memory issues in some functions.

* Corrected class lefkoProj structure output by function f_projection3() to
  correspond exactly to that produced by function projection3().

* Fixed bug preventing proxy multipliers from being properly integrated in
  function-based and raw MPMs.

* Corrected issue in which truncated distributions occasionally led to NA
  size probabilities when the probability of a 0 under the un-truncated
  distribution was sufficiently high.

# lefko3 5.0.1 (2022-03-08)

## USER-VISIBLE CHANGES

* Updated help files for projection functions.

## BUG FIXES

* Updated primitive class checks in some functions.

# lefko3 5.0.0 (2022-03-01)

## NEW FEATURES

* New function arlefko2() created to handle raw age-by-stage MPM estimation.

* Function aflefko2() now fixes the starting age by whether the life history
  model is pre-breeding or post-breeding.

* Functions verticalize3() and historicalize3() now handle age offsets and
  automatically offset age by 1 if set to prebreeding = TRUE.

* Further output has been added when err_check = TRUE is set in all
  function-based matrix estimators.

* New function f_projection3() runs MPM and IPM projections in which matrices
  are developed from vital rate models at each step.

## USER-VISIBLE CHANGES

* Default final_age in function aflefko2() is now the final age at time t in the
  original dataset.

* All function-based matrix estimators have been sped up.

* Population dynamics functions have been sped up.

* All population dynamics functions involving actual projection can now be
  cancelled during runs.

* Function create_lM() now produces a properly edited stageframe usable with all
  other functions taking lefkoMat objects as input.

## BUG FIXES

* Function plot.lefkoProj() now uses the correct default axis labels when
  conducting a state-space plot.

* Function hist_null() now correctly estimates the numbers of survival and
  fecundity elements.

* Eliminated bug in functions flefko3(), flefko2(), and aflefko2()
  preventing incorporation of random year effects into matrix estimation.

* Eliminated bug causing inflation of lambda under some instances of the
  Poisson distribution in function-based matrix estimation.

* Fixed erroneous mixing of QC information in lmean() output.

* Fixed issue preventing multiple replicates from being displayed in
  plot.lefkoProj under the timeseries option.

* Fixed issue causing certain impossible transitions linked to the prior stage
  to be estimated as non-zero in deVries-format function-based hMPMs.

* Fixed incorrect formatting of deVries hMPMs resulting from function
  hist_null().

* Function create_pm() updated to handle maturity status models.

# lefko3 4.2.0 (2022-01-17)

## NEW FEATURES

* Functions rleslie() and fleslie() created to estimate raw and function-based
  Leslie MPMs, respectively.

* Added a beta term setting K as hard maximum in the logistic function using
  functions density_input() and projection3().

* Accuracy check added to size and fecundity models in modelsearch().

* Function modelsearch() new runs exhaustive model selection under the
  `suite = "cons"` setting if age, density, or any individual covariate is added
  to the analysis.

* Function modelsearch() now performs accuracy checks by default, and allows an
  opt-out.

* Function actualstage3() now removes stages if desired, including stage
  NotAlive, which is now the default.

* New plot.lefkoProj() function added.

* All function-based matrix estimators now incorporate a default method for
  continuous distributions using cumulative density functions to estimate size
  transition probabilities. The midpoint rule is now only an option.

* New standalone density dependence functions ricker3(), beverton3(), usher3(),
  and logistic3() added.

* Substochasticity options for projection3() now include the option to force
  fecundity to stay non-negative.

* Function projection3() now issues warnings if density dependence yields matrix
  impossible matrix values.

* Function modelsearch() now also estimates maturity status probability in
  juveniles.

* Functions supplemental(), density_input(), and start_input() now also include
  automated handling of keywords "obs" and "nobs", which specify the use of all
  observable or unobservable stages, respectively.

## USER-VISIBLE CHANGES

* Functions rlefko2() and rlefko3() now provide warnings if individual identity
  variable name or column number is not provided.

* Function summary.lefkoMat() now provides appropriate descriptions of data
  quality control sections with missing values.

* Added spatial coordinates to cypdata and cypvert.

* Function modelsearch() is now more flexible with input errors in some key
  settings.

* Ended dependence on package stringr.

* Expanded data frame output with actualstage3() to give individual stage
  designations in times t and t-1.

* Function summary.lefkoMat() now yields warnings if any fecundity matrix
  contains a negative value, and if any matrix contains an NA value.

* Added equations in key help files.

* Corrections and updates to all vignettes.

* Function modelsearch() now tries several more parameterizations if the global
  model fails.

* Function lmean() now also incorporates modelqc and dataqc elements for quality
  control purposes.

* Function aflefko2() now treats missing minimum ages as zeroes, allowing
  estimation even when age info is missing from stageframes.

## BUG FIXES

* Matrix quality control summaries corrected when functions add_lM(),
  delete_lM(), and subset_lM() are used.

* Corrected errors in documentation regatding censoring variables used in
  rlefko3() and rlefko2().

* Corrected issue yielding fatal error in historicalize3() when spatial
  coordinates are set and reduce = TRUE.

* Corrected improper matrix formatting in some integer variables used in
  historicalize3().

* Fixed bug ignoring coordsRepeat setting in verticalize3().

* Fixed bug causing inconsistent coordinate assignment under unobserved stages
  in historicalize3().
  
* Corrected bug capable of causing the logistic function to yield negative
  matrix elements.

* Fixed improper handling of numeric entries for year and patch options in
  function modelsearch().

* Corrected issue in flefko3(), flefko2(), and aflefko2() affecting the
  incorporation of zero-inflated Poisson and negative binomial models.

* Fixed fecundity estimation issue arising in some circumstances with
  deVries-formatted historical models.

* Fixed bug leading size of 0 to be treated as zero-inflated in Poisson and
  negative binomial models.

* Added tolerance limit to binomial probability estimation in all function-based
  matrix estimators, to prevent some conditions yielding NA probabilities.
  
* Eliminated bug causing rare error in formulae creation in modelsearch().

* Eliminated bug causing unrealistic pseudo-R^2 in rare runs of modelsearch().

* Eliminated fatal error in modelsearch() occurring under rare instances when
  all attempts to build a working global model fail.

* Fixed bug incorrectly importing non-standard variable names into
  modelsearch().

* Fixed error causing density dependent projections under large, historical MPMs
  in projection3() to issue substochasticity warnings continuously regardless of
  actual substochasticity violation.

# lefko3 4.1.1 (2021-12-07)

## BUG FIXES

* Fixed C++ runtime error in function hist_null().

# lefko3 4.1.0 (2021-12-06)

## NEW FEATURES

* Included option to ignore observation status during stage assignment in
  verticalize3() and historicalize3().

* Included option to quiet all warnings in verticalize3() and historicalize3().

* New options to prevent display of overdispersion and zero-inflation tests in
  sf_distrib().

* Added hist_null() function to create historically-structured MPM to serve as a
  null model MPM assuming no individual history.

* Created function diff_lM() that creates lefkoDiff objects, which are MPM style
  objects that are actually matrix differences between two lefkoMat objects with
  equivalent dimensions.

* Added cond_diff() function to create conditional difference matrices from
  supplied lefkoDiff objects.
  
* Created actualstage3() function to calculate actual stage frequencies and
  proportions from input hfv dataset.

* Created create_pm() function to create a skeleton paramnames object for users
  building vital rate models manually.

## USER-VISIBLE CHANGES

* Introduced an integer check in function sf_distrib() that stops tests on
  non-integer variables.

* All vignettes have been updated to the latest functions and results.

## BUG FIXES

* Fixed incorrect time step calculation in summary.lefkoMat().

* Fixed incorrect P values for zero-inflation when fewer than expected 0s found
  in function sf_distrib().

* Fixed fatal issue affecting certain kinds of proxy transitions in
  function-based matrix estimation.

* Fixed incorrect stage distribution calculation in stochastic stablestage3()
  and projection3().

* Fixed fatal error in summary.lefkoProj() that occurred when summarizing
  projections with single replicates.

* Fixed fatal errors in flefko2(), flefko3(), aflefko2(), and summary.lefkoMat()
  relating to the use of user-defined vital rate functions.

* Corrected incorrect assignment of stages in flefko3(), flefko2(), and
  aflefko2() based on use of groups in supplement tables.

# lefko3 4.0.1 (2021-11-15)

## BUG FIXES

* Function aflefko2() now properly removes agestages rows when reduce = TRUE.

* Problems with mathjaxr integration fixed.

# lefko3 4.0.0 (2021-11-14)

## NEW FEATURES

* Function projection3() can now run density dependent projections using
  the Ricker, Beverton-Holt, Usher, and logistic functions.

* Functions start_input() and density_input() have been created to describe
  starting vectors and density dependence in the simplest format possible.

* Function modelsearch() now also calculates the accuracy of binomial models and
  adds the output to the QC section of the output. Function summary.lefkoMod()
  reports this output.

* Function sf_create() now handles dual- and triple-size classification.

* Functions sf_create() and supplemental() now handle the designation of core
  stage groups.

* Functions verticalize3() and historicalize3() now handle dual- and triple-size
  metric classification.

* Function modelsearch() can now test the impacts of three different size
  variables, and can also develop vital rate models for up to three separate
  response size variables.

* Function modelsearch() can now test the impacts of spatial density and stage
  group.

* Function modelsearch() now allows individual covariates to be incorporated as
  random categorical variables, as well as fixed continuous quantitative
  variables.

* Function modelsearch() now also allows the Gamma distribution for size and
  fecundity.

* Function flefko2() now handles up to three size metrics, stage groups, random
  categorical variables, and density.

* Function sf_distrib() now handles up to three size metrics, and also handles
  fecundity estimation in both times t and t+1.

* Function projection3() now allows two forms of enforced substochasticity. The
  first keeps all survival-transition probabilities within bounds, while the
  second keeps all stage survival probabilities within bounds.

## USER-VISIBLE CHANGES

* Function projection3() now creates easier-to-understand output structures
  clearly dividing numbers of individuals from stage distributions from
  reproductive values and population sizes.

* Function supplemental() now takes proxy multipliers, and gives errors for
  negative given rates and multipliers.

* Corrected some figure dimensions in vignettes.

* Corrected help text for some functions.

* Functions sf_distrib(), summary.lefkoMat(), summary.lefkoCondMat(),
  summary.lefkoMod(), and image3() no longer return a NULL to the console.

* The default for option censorRepeat in function verticalize3() has been
  changed to FALSE, in order to prevent overzealous censoring by default.

* Text output from summary functions and sf_distrib() has been altered to fit
  within standard output consoles, including typical page widths.

* Added sparse matrix override to slambda3() in case of errors incurred in
  eigen analysis.

## BUG FIXES

* Corrected bug causing code "npr" to yield propagule stages in time t+1 in
  supplement and overwrite tables.

* Fixed issue preventing estimation of size transition into observable classes
  with size 0 under the Poisson distribution in flefko3() and flefko2().

* Corrected porous stage assignment handling in historicalize3().

* Eliminated incorrect projection set assignment when running multiple patches
  or populations through function projection3().

* Eliminated bug giving warnings when running rlefko2() and rlefko3() with
  specific choices of patches or populations given as input options.

* Fixed issue causing rlefko2() and rlefko3() to estimate only some
  transition-survival probabilities and fecundities when patchcol or popcol is
  given without any respective specific patch or pop options given.

# lefko3 3.8.0 (2021-09-08)

## NEW FEATURES

* Function projection3() now includes the ability to produce replicate
  simulations.

* Function summary.lefkoProj() created to summarize the results of population
  projection using function projection3().

* Functions verticalize3() and historicalize3() now handle true radial density
  estimation, if X and Y coordinates and a valid spacing threshold are supplied.

* Functions add_lM() and delete_lM() added to add and delete matrices from
  lefkoMat objects.

* Function subset_lM() added to create new lefkoMat objects as subsets of
  other lefkoMat objects.

## USER-VISIBLE CHANGES

* Sensitivity and elasticity descriptions corrected in vignettes.

* Increased consistency in terminology used in output for stablestage3() and
  repvalue3().

* Function create_lM() now handles age-by-stage matrices.

## BUG FIXES

* Fixed issue that could lead rlefko2() and rlefko3() to stop processing
  matrices if specific years or patches are called.

* Fixed issue affecting stochastic reproductive value estimation in historical
  lefkoMat objects.

* Fixed issue affecting some deterministic reproductive value estimates,
  particularly in large historical MPMs, where a base reference value for
  standardization cannot be determined.

* Fixed issue preventing proper handling of historical matrix inputs in function
  create_lM().

# lefko3 3.7.0 (2021-08-18)

## NEW FEATURES

* Function ltre3() developed to conduct deterministic and stochastic life table
  response experiments.

* Function create_lM() developed to allow the incorporation of lists of matrices
  developed outside of package lefko3.

* Function summary.lefkoLTRE() developed to quickly summarize LTRE contributions
  according to transition type.

* Added new vignette to showcase LTRE and sLTRE analysis, and the import of
  matrices as lefkoMat objects.

## USER-VISIBLE CHANGES

* Altered some terminology in help files and vignettes to reflect most recent
  best practice in population ecology. Updated vignettes with new information,
  particularly related to LTRE and sLTRE analysis.

* Option append_mats has been added to sensitivity3.lefkoMat(),
  sensitivity3.list(), elasticity3.lefkoMat(), and elasticity3.list(). This
  option allows users to prevent the incorporation of the original matrices into
  output lefkoSens and lefkoElas objects, resulting in lower memory usage.

* Function summary.lefkoElas() can now handle age-by-stage MPMs.

* Corrected problematic Cypripedium examples involving the construction of
  matrices with stage survival probabilities rising above 1.0.

## BUG FIXES

* Fixed indexing bug in function summary.lefkoElas() that had led to incorrect
  elasticity summation.

* Fixed memory leak affecting sparse matrices in projection3() and related
  functions.

# lefko3 3.6.0 (2021-07-21)

## NEW FEATURES

* Automatic sparse matrix detection and override options have been added to all
  population dynamics functions.

## USER-VISIBLE CHANGES

* All population dynamics analysis functions have been substantially sped up
  through better memory management techniques.

# lefko3 3.5.3 (2021-07-14)

## NEW FEATURES

* Function summary.lefkoMat() now also displays column sums from all U
  (survival-transition) matrices as a quality control check.

## USER-VISIBLE CHANGES

* Package glmmTMB has been reintegrated, allowing mixed modeling of negative
  binomial, zero-inflated, and zero-truncated response distributions again.

* Error and warning messages from compiled functions provide cleaner messages.

## BUG FIXES

* Corrected incorrect handling of missing patch terms in lm, glm, vglm, and
  zeroinfl objects.

* Corrected static_cast issue that might prevent package installation in
  Solaris systems.

# lefko3 3.5.2 (2021-07-12)

## USER-VISIBLE CHANGES

* Negative binomial response, zero-truncated, and zero-inflated have been
  temporarily disabled in function modelsearch() while package glmmTMB() is
  being reprogrammed for compatibility with R 4.1.0. These distributions are
  still available through the glm approach.

* Examples have been updated for clarity and consistency.

# lefko3 3.5.1 (2021-07-07)

## NEW FEATURES

* Function supplemental() now allows entry stage proxies to be marked for raw
  historical MPMs.

* Shortcuts included covering non-propagule stages and non-reproductive mature
  stages for functions overwrite() and supplemental().

## USER-VISIBLE CHANGES

* Updated Lathyrus vignettes to correct improper handling of entry stages.

* Negative binomial response in mixed modeling has been moved to package lme4,
  and zero-truncated and zero-inflated distributions are available only through
  the glm approach in modelsearch() for the time being.

## BUG FIXES

* Corrected improper handling of prior forms of entry stages in raw historical 
  MPMs.

# lefko3 3.5.0 (2021-06-29)

## NEW FEATURES

* Historical MPM estimating functions can now estimate hMPMs in deVries format,
  in which newborns have an unique prior stage. All population dynamics
  functions have also been reworked to handle this format.

* New output has been added showing the exact order of age-stage combinations
  across all rows in estimated age-by-stage MPMs, as output from function
  aflefko2().

* Stochastic sensitivities of historical MPMs now yield ahistorical equivalents
  in all output.

* Added manual designation of ahistorical vs. historical output in functions
  sensitivity3.list() and elasticity3.list().

## USER-VISIBLE CHANGES

* Function modelsearch() can now handle invariant response terms.

* Functions flefko3() and flefko2() can now handle invariant fecundity.

* Function verticalize3() can now handle horizontal datasets without clear
  patterns in order.

* Amended overwrite() and supplemental() functions to account for possible
  prior newborn stage in deVries formatted hMPMs.

* Function supplemental() now details stage names not accounted for in the
  input stageframe in error messages.

* Some error and warning messages have been clarified.

* Many examples have been expanded to include extra conditions.

## BUG FIXES

* Function rlefko3() now handles dimension reduction properly.

* Function cond_hmpm() now produces conditional matrices from all combinations
  of population, subpopulation, and year, and also properly labels them.

* Typos and other issues corrected in vignettes.

# lefko3 3.4.0 (2021-03-31)

## NEW FEATURES

* Zero-truncated Poisson and negative binomial distributions have been added to
  function modelsearch(), and as underlying size and fecundity distributions in
  flefko3(), flefko2(), and aflefko2().

* New function image3() created to easily create matrix images for lefkoMat and
  other objects. Function acts as a wrapper for the image() function in package
  SparseM.

* New vignette showcasing the estimation and analysis of age x stage MPMs.
  
## USER-VISIBLE CHANGES

* Added err_check option to function-based matrix estimators, allowing the
  output of vital rates used in the estimation of U matrices.

* The test used to assess overdispersion in size and fecundity has been changed
  to deal more accurately with count-based variance (including tests performed
  and how data are subsetted), and to offer more choices in which tests to run.

* Function modelsearch() has been rewritten to result in a smaller installed
  package size.

* Function parasearch() has been removed pending revision of modeling
  methodology.

* Reproductive value vectors are now all standardized to the first non-zero
  value in the vector.

* Function sf_create() can now handle stage comments as input.

* Functions sensitivity3() and elasticity3() have been redesigned for increased
  speed with unusually large matrices.

* Stochastic analysis functions have been streamlined for speed and efficiency.

* Expanded all examples to include both Lathyrus and Cypripedium versions.

* Corrections and expansions to package vignettes.
  
## BUG FIXES

* Fixed issue in which loss of one of the year or patch terms in modeling of
  zero-inflated mixed models could lead to errors in function-based matrix
  estimation.

* Corrected error in estimation algorithm for size transition probability under
  negative binomial distribution, and difficulty in handling low levels of
  overdispersion in some negative binomial models.

* Fixed issue making lmean() unable to take element-wise mean matrices of
  age-by-stage MPMs produced using aflefko2().

* Fixed bug yielding erroneous reproductive values in historical matrices.

* Fixed bug yielding erroneous sensitivity elements in historically-corrected
  ahistorical sensitivity matrices.

* Fixed bug occasionally yielding 3d arrays in response to stochastic calls of
  sensitivity3() and elasticity3().

* Corrected stablestage3() and repvalue3() to properly handle age x stage MPMs.

* Corrected memory leak issue in aflefko2().

# lefko3 3.3.2 (2021-02-25)

## BUG FIXES

* Fixed Imports field related to function parasearch().

* Fixed loss of fecundity placement in age-by-stage matrices when using function
  supplemental() with function aflefko2().

* Fixed default year and patch settings in aflefko2().

* Corrected memory allocation issue in historicalize3().

# lefko3 3.3.1 (2021-02-23)

## USER-VISIBLE CHANGES

* Corrections to vignettes.

* Help file typos fixed.

# lefko3 3.3.0 (2021-02-21)

## NEW FEATURES

* Function supplemental() created to allow greater flexibility in the input of
  reproductive multipliers. All matrix creation functions can now handle its
  output.

* Stochastic sensitivity and elasticity analyses now enabled through functions
  sensitivity3() and elasticity3(), respectively.

* Functions sensitivity3() and elasticity3() now handle simple lists of A
  matrices, in addition to lefkoMat objects and simple matrices.

* Function parasearch() created to allow parallelized model building and
  selection.
  
## USER-VISIBLE CHANGES

* Censoring options in verticalize3() and historicalize3() can now handle both
  static and temporally-variable censor variables.

* Overwrite and supplemental tables can now include "prop", "immat", and "mat"
  as shorthand for suites of all propagule, immature, and mature stages.

* Function modelsearch() now gives warning if a dataset is used with NAs in
  individual ID when a mixed modeling approach is requested. It now also
  provides standard diagnostic messages at each step.

* Corrections and additions to help files.
  
## BUG FIXES

* Fixed censoring protocols in verticalize3(), historicalize3() rlefko3(), and
  rlefko2().

* Corrected distribution algorithm for zero inflation and text output in
  sf_distrib().

* Eliminated bug causing high overdispersion parameters from negative
  binomial-based size models to yield function-based matrices with NA values for
  all transitions to stages with positive size.

* Fixed issue in which mixed models with a non-zero-inflated negative binomial
  distribution fail to yield function-based matrices.

* Fixed compatibility issues in matrix creation functions with R 3.6.3.

* Fixed minor issue in summary.lefkoMat() output text.

# lefko3 3.2.0 (2021-01-03)

## NEW FEATURES

* Function cond_hmpm() added to create conditional hMPMs.

* Function slambda3() added to estimate the log stochastic growth rate.

* Function summary.lefkoElas() added to summarize elasticities by kind of
  transition.

* Added citation() data for package.

* Added NEWS section, using package lme4 NEWS as a template.
    
## USER-VISIBLE CHANGES

* Corrected inconsistent stage name variable in stage frame creation and
  manipulation functions.

* Updated and reorganized vignettes.

* Modified stageframes exported by matrix creation functions now include entry
  status variable.

* Objects of class lefkoElas now include the original A, U, and F matrices.

* Elasticity and sensitivity outputs for lefkoMat objects now include the
       original A, U, and F matrices used as input.
  
## BUG FIXES

* Matrix estimation function now create loy tables treating pop and patch as
  strings by default, eliminating possible conversion errors.

* Corrected incorrect overwrite() call in lathyrus example used in function
  flefko3() and all population dynamics analysis functions.

# lefko3 3.1.2 (2020-11-16)

## BUG FIXES

* Corrected auto-conversion of characters to factors occurring on R 3.6.3,
  which affected the creation and interpretation of the $labels element in
  atrix estimation.

# lefko3 3.1.1 (2020-11-13)

## BUG FIXES

* Fixed bug in lmean() function resulting from an implicit cast and affecting
  users operating lefko3 on Solaris systems.

# lefko3 3.1.0 (2020-11-08)

## NEW FEATURES

* Added function sfdistrib() to test whether mean = var and the level of
  zero-inflation in count variable data to be used for size and fecundity.

* Added zero-inflated Poisson and negative binomial distributions as choices for
  the underlying distribution of size and fecundity.

* Added individual and environmental covariates to function-based MPM estimation
  functions.

## USER-VISIBLE CHANGES

* Corrected, updated, and reorganized vignettes.

## BUG FIXES

* Corrected parameterization of negative binomial distribution in function
  modelsarch().

# lefko3 3.0.0 (2020-10-22)

## NEW FEATURES

* Added function aflefko2() function to estimate age x stage MPMs.

* Added function sensitivity3() to estimate sensitivity matrices of historical
  and ahistorical MPMs.

* Added Function elasticity3() to estimate elasticity matrices of historical and
  ahistorical MPMs.

## USER-VISIBLE CHANGES

* Increased consistency of output variable names across functions.

* Corrected, updated, and reorganized vignettes.

* Eliminated faulty parameterization of matrix geometric mean in lmean()
  function.

## BUG FIXES

* Corrected core kernel underlying lmean(), and redeveloped in C++ to yield
  faster results.

* Corrected parameterization of negative binomial distribution.

* Function summary.lefkoMat() now also shows which vital rate models were not
  estimated.

# lefko3 2.4.2 (2020-09-12)

## USER-VISIBLE CHANGES

* All matrix estimators sped up through lapply()-based calls to C++ kernels.

* Function lmean() now keeps the original names of populations and patches.

* Corrected incorrect reference links in vignettes.

## BUG FIXES

* Fixed memory leak in kernel powering function historicalize3().

* Fixed faulty indexing affecting all matrix estimators under certain
  conditions.

* Bugs in Cypripedium candidum vignettes and examples caused by incomplete
  overwrite tables fixed.

# lefko3 2.3.0 (2020-08-18)

## USER-VISIBLE CHANGES

* Major speed up to flefko2() and flefko3() through the incorporation of
  C++ core kernels.

## BUG FIXES

* Bug fixed that could lead to the erroneous incorporation of mature stages as
  immature in flefko3() under certain rare circumstances.

* Bug fixed that could fit proxy rates into the wrong elements in all matrix
  construction routines under certain rare circumstances.

* Bug fixed that could, under some circumstances, lead to individuals being
  treated as dead when not observed for a period of more than a year.

# lefko3 2.2.2 (2020-07-28)

## USER-VISIBLE CHANGES

* First version available on CRAN.

* Revised vignettes.

* Improved default handling of reproductive ratio and fecundity ratio.

## BUG FIXES

* Fixed calculation of observed lifespan in functions verticalize3() and
  historicalize3().

# lefko3 2.2.1 (2020-07-18)

## USER-VISIBLE CHANGES

* Removed some vignettes to reduce size.

* Added greater detail to Metadata in DESCRIPTION file

## BUG FIXES

* Added on.exit() to modelsearch() to keep user options through that routine.

* Fixed bugs in repvalue3.lefkoMat() and repvalue3.matrix().

# lefko3 2.2.0 (2020-07-08)

## USER-VISIBLE CHANGES

* All internal function documentation now complete.

* Reduced vignette sizes.

## BUG FIXES

* Fixed bugs in summary.lefkoMod().

# lefko3 2.1.0 (2020-07-06)

## NEW FEATURES

* Advanced routines developed for deterministic analysis functions, including
  repvalue3() and stablestage3().

## USER-VISIBLE CHANGES

* Cleaned up all help pages.

## BUG FIXES

* Scattered bugs eliminated.

# lefko3 2.0.0 (2020-06-27)

## NEW FEATURES

* Stageframe objects now have their own class, and stageframes now handle
  log-size based bin widths.

* Function sf_create() also now stops if vector lengths are not equal, and can
  create stageframes for IPMs.

* Mean matrix function lmean() created.

* Added size in time t to juvenile models.

* Fecundity models can now use fecundity in t+1 as a response.

* Integrated estimation function into verticalize3() and historicalize3().

## USER-VISIBLE CHANGES

* All matrix creation functions have been fully rewritten for speed (eliminated
  extraneous find() calls and long loops).

* Vignettes split up across chapter files.

## BUG FIXES

* Eliminated problems with class interpretation in modelsearch().

# lefko3 1.0.0 (2020-03-24)

## NEW FEATURES

* Lathyrus dataset added.

* Examples added to all core functions.

## USER-VISIBLE CHANGES

* All matrix estimation functions have been rewritten to include core C++
  functions, for speed up and memory efficiency.

* Added patch functionality to modelsearch(), flefko3(), and flefko2().

## BUG FIXES

* Bugs corrected in probability and rate estimation functions.

# lefko3 0.13.1 (2020-02-16)

## USER-VISIBLE CHANGES

* Reorganized flefko2() around Rcpp functions.

## BUG FIXES

* Cleaned up bugs discovered in rlefko3().

# lefko3 0.13.0 (2020-02-14)

## USER-VISIBLE CHANGES

* Re-developed core functions of rlefko3() with Rcpp and RcppArmadillo,
  substantially reducing run time.

# lefko3 0.12.1 (2020-02-11)

## USER-VISIBLE CHANGES

* Re-developed rlefko3() function per new approach to rlefko2().

# lefko3 0.12.0 (2020-02-05)

## NEW FEATURES

* Developed rlefko2() function.

## USER-VISIBLE CHANGES

* Redesigned all matrix functions to handle different populations and patches.

* Made core matrix creation kernels quiet.

* Increased standardization of inputs.

* Eliminated class lefkoMatMulti.

# lefko3 0.11.0 (2020-01-31)

## NEW FEATURES

* Created the hfvdata class for standardized, vertical demographic data frames.

* Added quality control sections to lefkoMod and lefkoMatMulti objects,
  including functions that create these objects.

* Created new S3 summary.lefkoMatMulti() function.

# lefko3 0.10.2 (2020-01-29)

## NEW FEATURES

* Added negative binomial functionality to vital rate estimators.

## USER-VISIBLE CHANGES

* Massive changes and improvements to tutorials.

# lefko3 0.10.1 (2020-01-22)

## NEW FEATURES

* Added negative binomial functionality to modelsearch().

* Added an S3 summary function for class lefkoMod.

# lefko3 0.10.0 (2020-01-19)

## NEW FEATURES

* Added Cypripedium candidum vignette.

* Data mgmt functions can now prune unused variables from output.

## USER-VISIBLE CHANGES

* Replaced all flw and frt entries with repst and fec.

# lefko3 0.9.0 (2020-01-15)

## NEW FEATURES

* Added cypdata data file.

## USER-VISIBLE CHANGES

* Formatted code and help files.

# lefko3 0.8.3 (2020-01-11)

## USER-VISIBLE CHANGES

* Eliminated year effect warnings by option.

# lefko3 0.8.2 (2020-01-10)

## NEW FEATURES

* Added the lefkoMatMulti class.

* Gave matrix creation functions the ability to create output with lists of A,
  T, and F matrices for all years requested.

# lefko3 0.8.1 (2020-01-10)

## NEW FEATURES

* Added greater flexibility to modelsearch() function to allow it to estimate
  binomial models under constant response.

## USER-VISIBLE CHANGES

* Synchronized data management function to modelsearch() function.

# lefko3 0.8.0 (2020-01-01)

## NEW FEATURES

* Added modelsearch() function to handle model building and selection workflow.

# lefko3 0.7.0 (2019-11-22)

## NEW FEATURES

* All workhorse functions are now generally debugged C++ functions.

## BUG FIXES

* Fixed overwriting and estimating portions of matrix creation functions.

# lefko3 0.6.1 (2019-11-22)

## BUG FIXES

* Fixed overwriting portion of matrix design, so that given transitions are
  properly added to historical matrices.

# lefko3 0.6.0 (2019-11-20)

## NEW FEATURES

* Fully rebuilt flefko3() and flefko2() so that fecundity is now estimated using
  binaries written in C++.

## USER-VISIBLE CHANGES

* Matrix design reworked to make everything faster.

# lefko3 0.5.1 (2019-11-19)

## NEW FEATURES

* Rebuilt flefko2() to use C++ precompiled probability density functions.

# lefko3 0.5.0 (2019-11-18)

## NEW FEATURES

* Rewrote workhorse portions of survival transition estimators as precompiled
  binaries in C++.

* Rewrote historical matrix function (flefko3()) to handle reworked pxy
  functions.

## BUG FIXES

* Corrected issues in help files.

# lefko3 0.4.4 (2019-11-11)

## NEW FEATURES

* Added a patch term to ahv2hv() function.

## USER-VISIBLE CHANGES

* Changed plantid variable name to individ.

# lefko3 0.4.3 (2019-11-08)

## NEW FEATURES

* Added density3() function.

## USER-VISIBLE CHANGES

* Sped up verticalize3() and ahv2hv() functions.

## BUG FIXES

* Corrected links throughout.

# lefko3 0.4.2 (2019-11-07)

## NEW FEATURES

* Added ability to handle character column names to verticalize_h().

# lefko3 0.4.1 (2019-11-06)

## BUG FIXES

* Corrected issues in verticalize_h function.

# lefko3 0.4.0 (2019-11-03)

## NEW FEATURES

* Added RSpectra-based population matrix projection analysis functions.

## USER-VISIBLE CHANGES

* Added GPL license.

# lefko3 0.3.1 (2019-11-02)

## BUG FIXES

* Corrected bug leading to infinite mu values in probbin for some models.

* Corrected the method of using proxy transitions from certain juvenile
  stages to mature stages in all lefko functions.

# lefko3 0.3.0 (2019-11-01)

## USER-VISIBLE CHANGES

* Increased ability to choose which unobserved stages lead to which observed
  stages.

* Standardized terms in more functions for cross-compatibility and easier
  workflow.

# lefko3 0.2.1 (2019-10-29)

## NEW FEATURES

* Added Gaussian probability calculation.

* Developed sf_create() function.

* First version of rlefko3() created.

# lefko3 0.1.1 (2019-10-26)

## NEW FEATURES

* Matrix creation functions now discriminate size and fecundity on the basis of
  probability distribution.

## USER-VISIBLE CHANGES

* Added numeric class recognition.

* Documentation improvements.

# lefko3 0.1.0 (2019-10-25)

## NEW FEATURES

* First development version. Main 2d and 3d matrix functions created.

