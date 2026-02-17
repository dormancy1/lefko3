#' Demographic Dataset of \emph{Cypripedium candidum} Population, in Horizontal
#' Format
#' 
#' A dataset containing the states and fates of \emph{Cypripedium candidum} 
#' (white lady's slipper orchids), family Orchidaceae, from a population in 
#' Illinois, USA, resulting from monitoring that occurred annually between 2004 
#' and 2009.
#' 
#' @docType data
#' 
#' @usage data(cypdata)
#' 
#' @format A data frame with 77 individuals and 29 variables. Each row 
#' corresponds to an unique individual, and each variable from \code{size.04} 
#' on refers to the state of the individual in a particular year.
#' 
#' \describe{
#'   \item{plantid}{A numeric variable giving a unique number to each 
#'   individual.}
#'   \item{patch}{A variable refering to patch within the population.}
#'   \item{X}{An X coordinate for the plant within the population.}
#'   \item{Y}{A Y coordinate for the plant within the population.}
#'   \item{censor}{A variable coding for whether the data point is valid. An
#'   entry of 1 means that it is so.}
#'   \item{Inf2.04}{Number of double inflorescences in 2004.}
#'   \item{Inf.04}{Number of inflorescences in 2004.}
#'   \item{Veg.04}{Number of stems without inflorescences in 2004.}
#'   \item{Pod.04}{Number of fruits in 2004.}
#'   \item{Inf2.05}{Number of double inflorescences in 2005.}
#'   \item{Inf.05}{Number of inflorescences in 2005.}
#'   \item{Veg.05}{Number of stems without inflorescences in 2005.}
#'   \item{Pod.05}{Number of fruits in 2005.}
#'   \item{Inf2.06}{Number of double inflorescences in 2006.}
#'   \item{Inf.06}{Number of inflorescences in 2006.}
#'   \item{Veg.06}{Number of stems without inflorescences in 2006.}
#'   \item{Pod.06}{Number of fruits in 2006.}
#'   \item{Inf2.07}{Number of double inflorescences in 2007.}
#'   \item{Inf.07}{Number of inflorescences in 2007.}
#'   \item{Veg.07}{Number of stems without inflorescences in 2007.}
#'   \item{Pod.07}{Number of fruits in 2007.}
#'   \item{Inf2.08}{Number of double inflorescences in 2008.}
#'   \item{Inf.08}{Number of inflorescences in 2008.}
#'   \item{Veg.08}{Number of stems without inflorescences in 2008.}
#'   \item{Pod.08}{Number of fruits in 2008.}
#'   \item{Inf2.09}{Number of double inflorescences in 2009.}
#'   \item{Inf.09}{Number of inflorescences in 2009.}
#'   \item{Veg.09}{Number of stems without inflorescences in 2009.}
#'   \item{Pod.09}{Number of fruits in 2009.}
#' }
#' 
#' @source Shefferson, R.P., R. Mizuta, and M.J. Hutchings. 2017. Predicting
#' evolution in response to climate change: the example of sprouting probability
#' in three dormancy-prone orchid species. \emph{Royal Society Open Science} 
#' 4(1):160647.
#' 
#' @examples 
#' data(cypdata)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
#'   "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   propstatus = propvector, immstatus = immvector, indataset = indataset,
#'   binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
#'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
#'   NRasRep = TRUE)
#' 
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#'                        
#' lambda3(cypmatrix2r)
"cypdata"

#' Demographic Dataset of \emph{Cypripedium candidum} Population, in Vertical
#' Format
#' 
#' A dataset containing the states and fates of \emph{Cypripedium candidum} 
#' (white lady's slipper orchids), family Orchidaceae, from a population in 
#' Illinois, USA, resulting from monitoring that occurred annually between 2004 
#' and 2009. Same dataset as \code{cypdata}, but arranged in an ahistorical
#' vertical format.
#' 
#' @docType data
#' 
#' @usage data(cypvert)
#' 
#' @format A data frame with 77 individuals, 322 rows, and 14 variables. Each
#' row corresponds to a specific two-year transition for a specific individual.
#' Variable codes are similar to those for \code{cypdata}, but use \code{.2} to
#' identify occasion \emph{t} and \code{.3} to identify occasion \emph{t}+1.
#' 
#' \describe{
#'   \item{plantid}{A numeric variable giving a unique number to each 
#'   individual.}
#'   \item{patch}{A variable refering to patch within the population.}
#'   \item{X}{An X coordinate for the plant within the population.}
#'   \item{Y}{A Y coordinate for the plant within the population.}
#'   \item{censor}{A variable coding for whether the data point is valid. An
#'   entry of 1 means that it is so.}
#'   \item{year2}{Year in occasion \emph{t}.}
#'   \item{Inf2.2}{Number of double inflorescences in occasion \emph{t}.}
#'   \item{Inf.2}{Number of inflorescences in occasion \emph{t}.}
#'   \item{Veg.2}{Number of stems without inflorescences in occasion \emph{t}.}
#'   \item{Pod.2}{Number of fruits in occasion \emph{t}.}
#'   \item{Inf2.3}{Number of double inflorescences in occasion \emph{t}+1.}
#'   \item{Inf.3}{Number of inflorescences in occasion \emph{t}+1.}
#'   \item{Veg.3}{Number of stems without inflorescences in occasion \emph{t}+1.}
#'   \item{Pod.3}{Number of fruits in occasion \emph{t}+1.}
#' }
#' 
#' @source Shefferson, R.P., R. Mizuta, and M.J. Hutchings. 2017. Predicting
#' evolution in response to climate change: the example of sprouting probability
#' in three dormancy-prone orchid species. \emph{Royal Society Open Science} 
#' 4(1):160647.
#' 
#' @examples 
#' data(cypvert)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
#'   "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector, 
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector, 
#'   propstatus = propvector, immstatus = immvector, indataset = indataset,
#'   binhalfwidth = binvec)
#' 
#' cypraw_v2 <- historicalize3(data = cypvert, patchidcol = "patch", 
#'   individcol = "plantid", year2col = "year2", sizea2col = "Inf2.2", 
#'   sizea3col = "Inf2.3", sizeb2col = "Inf.2", sizeb3col = "Inf.3", 
#'   sizec2col = "Veg.2", sizec3col = "Veg.3", repstra2col = "Inf2.2", 
#'   repstra3col = "Inf2.3", repstrb2col = "Inf.2", repstrb3col = "Inf.3", 
#'   feca2col = "Pod.2", feca3col = "Pod.3", repstrrel = 2, 
#'   stageassign = cypframe_raw, stagesize = "sizeadded", censorcol = "censor",
#'   censor = FALSE, NAas0 = TRUE, NRasRep = TRUE, reduce = TRUE)
#'   
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v2, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#'   
#' lambda3(cypmatrix2r)
"cypvert"

#' Demographic Dataset of \emph{Lathyrus vernus} Population
#' 
#' A dataset containing the states and fates of \emph{Lathyrus vernus} (spring
#' vetch), family Fabaceae, from a population in Sweden monitored annually
#' from 1988 to 1991 in six study plots.
#' 
#' @docType data
#' 
#' @usage data(lathyrus)
#' 
#' @format A data frame with 1119 individuals and 34 variables. Each row
#' corresponds to a unique individual, and each variable from \code{Volume88}
#' on refers to the state of the individual in a given year.
#' 
#' \describe{
#'   \item{SUBPLOT}{A variable refering to patch within the population.}
#'   \item{GENET}{A numeric variable giving a unique number to each 
#'   individual.}
#'   \item{Volume88}{Aboveground volume in cubic mm in 1988.}
#'   \item{lnVol88}{Natural logarithm of \code{Volume88}.}
#'   \item{FCODE88}{Equals 1 if flowering and 0 if not flowering in 1988.}
#'   \item{Flow88}{Number of flowers in 1988.}
#'   \item{Intactseed88}{Number of intact mature seeds produced in 1988.
#'   Not always an integer, as in some cases seed number was estimated via 
#'   linear modeling.}
#'   \item{Dead1988}{Marked as 1 if known to be dead in 1988.}
#'   \item{Dormant1988}{Marked as 1 if known to be alive but vegetatively 
#'   dormant in 1988.}
#'   \item{Missing1988}{Marked as 1 if not found in 1988.}
#'   \item{Seedling1988}{Marked as 1, 2, or 3 if observed as a seedling in year
#'   \emph{t}. Numbers refer to certainty of assignment: 1 = certain that plant
#'   is a seedling in 1988, 2 = likely that plant is a seedling in 1988,
#'   3 = probable that plant is a seedling in 1988.}
#'   \item{Volume89}{Aboveground volume in cubic mm in 1989.}
#'   \item{lnVol89}{Natural logarithm of \code{Volume89}.}
#'   \item{FCODE89}{Equals 1 if flowering and 0 if not flowering in 1989.}
#'   \item{Flow89}{Number of flowers in 1989.}
#'   \item{Intactseed89}{NZumber of intact mature seeds produced in 1989.
#'   Not always an integer, as in some cases seed number was estimated via
#'   linear modeling.}
#'   \item{Dead1989}{Marked as 1 if known to be dead in 1989.}
#'   \item{Dormant1989}{Marked as 1 if known to be alive but vegetatively 
#'   dormant in 1989.}
#'   \item{Missing1989}{Marked as 1 if not found in 1989.}
#'   \item{Seedling1989}{Marked as 1, 2, or 3 if observed as a seedling in
#'   year \emph{t}. Numbers refer to certainty of assignment: 1 = certain 
#'   that plant is a seedling in 1989, 2 = likely that plant is a seedling 
#'   in 1989, 3 = probable that plant is a seedling in 1989.}
#'   \item{Volume90}{Aboveground volume in mm<sup>3</sup> in 1990.}
#'   \item{lnVol90}{Natural logarithm of \code{Volume90}.}
#'   \item{FCODE90}{Equals 1 if flowering and 0 if not flowering in 1990.}
#'   \item{Flow90}{Number of flowers in 1990.}
#'   \item{Intactseed90}{NZumber of intact mature seeds produced in 1990.
#'   Not always an integer, as in some cases seed number was estimated via 
#'   linear modeling.}
#'   \item{Dead1990}{Marked as 1 if known to be dead in 1990.}
#'   \item{Dormant1990}{Marked as 1 if known to be alive but vegetatively 
#'   dormant in 1990.}
#'   \item{Missing1990}{Marked as 1 if not found in 1990.}
#'   \item{Seedling1990}{Marked as 1, 2, or 3 if observed as a seedling in
#'   year \emph{t}. Numbers refer to certainty of assignment: 1 = certain 
#'   that plant is a seedling in 1990, 2 = likely that plant is a seedling
#'   in 1990, 3 = probable that plant is a seedling in 1990.}
#'   \item{Volume91}{Aboveground volume in mm<sup>3</sup> in 1991.}
#'   \item{lnVol91}{Natural logarithm of \code{Volume91}.}
#'   \item{FCODE91}{Equals 1 if flowering and 0 if not flowering in 1991.}
#'   \item{Flow91}{Number of flowers in 1991.}
#'   \item{Intactseed91}{NZumber of intact mature seeds produced in 1991.
#'   Not always an integer, as in some cases seed number was estimated via
#'   linear modeling.}
#'   \item{Dead1991}{Marked as 1 if known to be dead in 1991.}
#'   \item{Dormant1991}{Marked as 1 if known to be alive but vegetatively 
#'   dormant in 1991.}
#'   \item{Missing1991}{Marked as 1 if not found in 1991.}
#'   \item{Seedling1991}{Marked as 1, 2, or 3 if observed as a seedling 
#'   in year \emph{t}. Numbers refer to certainty of assignment: 
#'   1 = certain that plant is a seedling in 1991, 2 = likely that plant 
#'   is a seedling in 1991, 3 = probable that plant is a seedling in 
#'   1991.}
#' }
#' 
#' @source Ehrlen, J. 2000. The dynamics of plant populations: does the 
#' history of individuals matter? \emph{Ecology} 81(6):1675-1684.
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET",
#'   juvcol = c("Seedling1988", "Seedling1989", "Seedling1990", "Seedling1991"),
#'   sizeacol = c("Volume88", "Volume89", "Volume90", "Volume91"),
#'   repstracol = c("FCODE88", "FCODE89", "FCODE90", "FCODE91"),
#'   fecacol = c("Intactseed88", "Intactseed89", "Intactseed90", "Intactseed91"),
#'   deadacol = c("Dead1988", "Dead1989", "Dead1990", "Dead1991"),
#'   nonobsacol = c("Dormant1988", "Dormant1989", "Dormant1990", "Dormant1991"),
#'   censorcol = c("Missing1988", "Missing1989", "Missing1990", "Missing1991"), 
#'   stageassign = lathframe, stagesize = "sizea",
#'   censorkeep = NA, censor = TRUE)
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
#'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
#'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' ehrlen3mean$A[[1]]
#' 
#' lambda3(ehrlen3mean)
"lathyrus"

#' Demographic Dataset of \emph{Pyrola japonica} and \emph{Pyrola subaphylla}
#' Populations, in Horizontal Format
#' 
#' A dataset containing the states and fates of \emph{Pyrola japonica} and
#' \emph{Pyrola subaphylla}, family Ericaceae, from populations in the vicinity
#' of Mt. Bandai, Fukushima Prefecture, Japan, resulting from monitoring that
#' occurred annually between 2015 and 2020.
#' 
#' @docType data
#' 
#' @usage data(pyrola)
#' 
#' @format A data frame with 454 individuals and 57 variables. Each row 
#' corresponds to an unique individual, and each variable from
#' \code{sprouted.2015} on refers to the state of the individual in a particular
#' year.
#' 
#' \describe{
#'   \item{species}{String denoting which of the two species the individual
#'   belongs to.}
#'   \item{population}{Integer denoting whcih population the individual belongs
#'   to. Synonymous with species in this dataset.}
#'   \item{id}{A numeric variable giving a unique number to each 
#'   individual within each species. Note that numbers are reused among the two
#'   species.}
#'   \item{sprouted.2015}{A binomial indicating whether the individual had
#'   living aboveground tissue observable in the 2015 census.}
#'   \item{lvs.num.2015}{Number of leaves in 2015.}
#'   \item{lvs.lng.2015}{Length of largest leaf in 2015.}
#'   \item{lvs.wdt.2015}{Width of largest leaf in 2015.}
#'   \item{inf.num.2015}{Number of inflorescences in 2015.}
#'   \item{inf.lng.tot.2015}{Summed inflorescence length in 2015.}
#'   \item{flo.tot.2015}{Number of flowers in 2015.}
#'   \item{frt.tot.2015}{Number of fruits in 2015.}
#'   \item{sprouted.2016}{A binomial indicating whether the individual had
#'   living aboveground tissue observable in the 2016 census.}
#'   \item{lvs.num.2016}{Number of leaves in 2016.}
#'   \item{lvs.lng.2016}{Length of largest leaf in 2016.}
#'   \item{lvs.wdt.2016}{Width of largest leaf in 2016.}
#'   \item{inf.num.2016}{Number of inflorescences in 2016.}
#'   \item{inf.lng.tot.2016}{Summed inflorescence length in 2016.}
#'   \item{flo.tot.2016}{Number of flowers in 2016.}
#'   \item{frt.tot.2016}{Number of fruits in 2016.}
#'   \item{sprouted.2017}{A binomial indicating whether the individual had
#'   living aboveground tissue observable in the 2017 census.}
#'   \item{lvs.num.2017}{Number of leaves in 2017.}
#'   \item{lvs.lng.2017}{Length of largest leaf in 2017.}
#'   \item{lvs.wdt.2017}{Width of largest leaf in 2017.}
#'   \item{inf.num.2017}{Number of inflorescences in 2017.}
#'   \item{inf.lng.tot.2017}{Summed inflorescence length in 2017.}
#'   \item{flo.tot.2017}{Number of flowers in 2017.}
#'   \item{frt.tot.2017}{Number of fruits in 2017.}
#'   \item{sprouted.2018}{A binomial indicating whether the individual had
#'   living aboveground tissue observable in the 2018 census.}
#'   \item{lvs.num.2018}{Number of leaves in 2018.}
#'   \item{lvs.lng.2018}{Length of largest leaf in 2018.}
#'   \item{lvs.wdt.2018}{Width of largest leaf in 2018.}
#'   \item{inf.num.2018}{Number of inflorescences in 2018.}
#'   \item{inf.lng.tot.2018}{Summed inflorescence length in 2018.}
#'   \item{flo.tot.2018}{Number of flowers in 2018.}
#'   \item{frt.tot.2018}{Number of fruits in 2018.}
#'   \item{sprouted.2019}{A binomial indicating whether the individual had
#'   living aboveground tissue observable in the 2019 census.}
#'   \item{lvs.num.2019}{Number of leaves in 2019.}
#'   \item{lvs.lng.2019}{Length of largest leaf in 2019.}
#'   \item{lvs.wdt.2019}{Width of largest leaf in 2019.}
#'   \item{inf.num.2019}{Number of inflorescences in 2019.}
#'   \item{inf.lng.tot.2019}{Summed inflorescence length in 2019.}
#'   \item{flo.tot.2019}{Number of flowers in 2019.}
#'   \item{frt.tot.2019}{Number of fruits in 2019.}
#'   \item{sprouted.2020}{A binomial indicating whether the individual had
#'   living aboveground tissue observable in the 2020 census.}
#'   \item{lvs.num.2020}{Number of leaves in 2020.}
#'   \item{lvs.lng.2020}{Length of largest leaf in 2020.}
#'   \item{lvs.wdt.2020}{Width of largest leaf in 2020.}
#'   \item{inf.num.2020}{Number of inflorescences in 2020.}
#'   \item{inf.lng.tot.2020}{Summed inflorescence length in 2020.}
#'   \item{flo.tot.2020}{Number of flowers in 2020.}
#'   \item{frt.tot.2020}{Number of fruits in 2020.}
#' }
#' 
#' @source Shefferson, R.P., K. Shutoh, and K. Suetsugu. \emph{In review}.
#' Vegetative dormancy and the evolution of mycoheterotrophy in sister
#' \emph{Pyrola} species. \emph{Journal of Ecology}.
#' 
#' @examples
#' 
#' data(pyrola)
#' 
#' pyrola$species <- as.factor(pyrola$species)
#' pyrola$population <- as.factor(pyrola$population)
#' jreg <- pyrola[which(pyrola$population == 1),]
#' stagevec_jp <- c("P1", "Sdl", "Dorm", "V0nr", "V1nr", "V2nr", "V3nr", "V4nr",
#'   "V0r", "V1r", "V2r", "V3r", "V4r")
#' sizeavec_jp <- c(0, 0, 0, 0, 1, 2, 3, 7, 0, 1, 2, 3, 7)
#' sizeahbin_jp <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 3.5, 0.5, 0.5, 0.5, 0.5,
#'   3.5)
#' repvec_jp <- c(0, 0, 0, 0, rep(0, 4), rep(1, 5))
#' propvec_jp <- c(1, rep(0, 12))
#' immvec_jp <- c(1, 1, rep(0, 11))
#' matvec_jp <- c(0, 0, rep(1, 11))
#' obsvec_jp <- c(0, 0, 0, rep(1, 10))
#' indata_jp <- c(0, 0, rep(1, 11))
#' comments_jp <- c("protocorm", "seedling", "dormant adult", "stump", "1lf nr",
#'   "2lf nr", "3lf nr", "4+lf nr", "0lf r", "1lf r", "2lf r", "3lf r",
#'   "4+lf r")
#' jp_frame <- sf_create(sizes = sizeavec_jp, stagenames = stagevec_jp,
#'   binhalfwidth = sizeahbin_jp, repstatus = repvec_jp, obsstatus = obsvec_jp,
#'   indataset = indata_jp, propstatus = propvec_jp, immstatus = immvec_jp,
#'   matstatus = matvec_jp, comments = comments_jp)
#'   
#' jhfv <- verticalize3(data = jreg, noyears = 6, firstyear = 2015,
#'   individcol = "id", blocksize = 8, sizeacol = "lvs.num.2015",
#'   obsacol = "sprouted.2015", repstracol = "flo.tot.2015",
#'   repstrbcol = "frt.tot.2015", fecacol = "flo.tot.2015",
#'   fecbcol = "frt.tot.2015", NAas0 = TRUE, stagesize = "sizea",
#'   stageassign = jp_frame)
#' 
#' surv_model <- glm(alive3 ~ sizea2 + as.factor(year2), data = jhfv, family = "binomial")
#' 
#' obs_data <- subset(jhfv, alive3 == 1)
#' obs_model <- glm(obsstatus3 ~ as.factor(year2), data = obs_data, family = "binomial")
#' 
#' size_data <- subset(obs_data, obsstatus3 == 1)
#' size_model <- glm(sizea3 ~ sizea2, data = size_data, family = "poisson")
#' 
#' reps_model <- glm(repstatus3 ~ sizea2, data = size_data, family = "binomial")
#' 
#' fec_data <- subset(jhfv, repstatus2 == 1)
#' fec_model <- MASS::glm.nb(fec2added ~ 1, data = fec_data)
#' 
#' mod_params <- create_pm(name_terms = TRUE)
#' mod_params$modelparams[4] <- "alive3"
#' mod_params$modelparams[5] <- "obsstatus3"
#' mod_params$modelparams[6] <- "sizea3"
#' mod_params$modelparams[9] <- "repstatus3"
#' mod_params$modelparams[11] <- "fec2added"
#' mod_params$modelparams[12] <- "sizea2"
#' mod_params$modelparams[18] <- "repstatus2"
#' 
#' jp_germ <- 0.90
#' jp_supp2 <- supplemental(stage3 = c("Sdl", "Dorm", "V0nr", "V1nr", "P1", "Sdl"), 
#'   stage2 = c("P1", "Sdl", "Sdl", "Sdl", "rep", "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA),
#'   givenrate = c(0.25, 0.35, 0.10, 0.10, NA, NA), # 0.345, 0.054
#'   multiplier = c(NA, NA, NA, NA, jp_germ * 0.5, jp_germ * 0.5),
#'   type = c(1, 1, 1, 1, 3, 3), stageframe = jp_frame, historical = FALSE)
#' 
#' jp_ahmpm <- flefko2(year = "all", stageframe = jp_frame, supplement = jp_supp2,
#'   paramnames = mod_params, surv_model = surv_model, obs_model = obs_model,
#'   size_model = size_model, repst_model = reps_model, fec_model = fec_model,
#'   data = jhfv, err_check = TRUE)
#' 
#' lambda3(jp_ahmpm)
"pyrola"

#' Matrix Set of \emph{Anthyllis vulneraria} Populations in Belgium
#' 
#' A \code{lefkoMat} object containing projection matrices developed from
#' demographic data gathered on nine \emph{Anthyllis vulneraria} populations
#' from 2003 to 2006 in southwestern Belgium.
#' 
#' @docType data
#' 
#' @usage data(anthyllis)
#' 
#' @format A \code{lefkoMat} object holding 27 matrices. The structure of the
#' object is as below:
#' 
#' \describe{
#'   \item{A}{The 27 A matrices.}
#'   \item{U}{The 27 survival-transition matrices used to develop the A
#'   matrices.}
#'   \item{F}{The 27 fecundity matrices used to develop the A matrices.}
#'   \item{hstages}{Not used, so set to \code{NA}.}
#'   \item{agestages}{Not used, so set to \code{NA}.}
#'   \item{ahstages}{The edited stageframe describing the life history of the
#'   study organism as interpreted in the original demographic study.}
#'   \item{labels}{The order of the matrices, where each population is treated
#'   as a separate patch and each matrix corresponds to a different combination
#'   of population and year in time \emph{t}.}
#'   \item{matrixqc}{A vector of integers used in the quality control section of
#'   \code{lefkoMat} summary statements.}
#'   \item{dataqc}{Currently a vector with two \code{NA} values.}
#' }
#' 
#' @source Davison, R. et al. 2010. Demographic effects of extreme weather
#' events on a short-lived calcareous grassland species: stochastic life table
#' response experiments. \emph{Journal of Ecology} 98(2):255-267.
#' 
#' @examples
#' data(anthyllis)
#' 
#' lambda3(anthyllis)
"anthyllis"

