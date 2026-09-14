#' World Bank WDI panel prepared for LHWR examples
#'
#' A country-year panel downloaded from the World Bank World Development
#' Indicators API and included with the package as an illustrative real-data
#' example for Lorentz Hyperbolic Weighted Regression. The data are real WDI
#' observations for 2010--2024, but the Lorentz-hyperboloid representation is a
#' documented demonstration of how raw economic variables can be converted into
#' the `Z` matrix required by LHWR. It should therefore be treated as a package
#' demo/template data set rather than as a definitive empirical application.
#'
#' The response and covariates are real-valued macroeconomic variables. The
#' columns `z0_lorentz`, `z1_lorentz`, and `z2_lorentz` define a
#' Lorentz-hyperboloid representation constructed from a two-component economic
#' feature map. The frozen preprocessing constants are retained explicitly:
#'
#' The index \eqn{i} identifies a country-year row. \eqn{\mathrm{GDPpc}_i}
#' denotes positive constant-price GDP per capita and its natural logarithm is
#' stored as `log_gdp_pc_const`. The first feature centers and scales that log:
#'
#' \deqn{
#' z_{i1} = \{\log(\mathrm{GDPpc}_i)-8.8071974413888245\}/1.4524958189486741,
#' }
#'
#' The second feature centers and scales trade openness, measured as a percent
#' of GDP in `trade_open`, using its separate frozen constants:
#'
#' \deqn{
#' z_{i2} = (\mathrm{trade\ openness}_i-91.8784599801564923)/58.6130974219034400,
#' }
#'
#' Both spatial features are dimensionless. The positive time coordinate is
#' then chosen to impose the unit-hyperboloid constraint:
#'
#' \deqn{
#' z_{i0} = \sqrt{1 + z_{i1}^2 + z_{i2}^2}.
#' }
#'
#' Hence each row `Z_i = (z0_lorentz, z1_lorentz, z2_lorentz)` lies on the
#' upper Lorentz hyperboloid. The coordinates satisfy
#' `-z0_lorentz^2 + z1_lorentz^2 + z2_lorentz^2 = -1` up to floating-point
#' tolerance. The location and scale constants were calculated on the
#' preprocessing reference sample before retaining this analysis panel. The
#' coordinates are not longitude-latitude coordinates; they encode an economic
#' neighborhood geometry based on GDP per capita and trade openness.
#'
#' @name wdi_lhwr
#' @format A data frame with 2262 rows and 18 columns:
#' \describe{
#'   \item{iso3}{ISO 3-letter country code.}
#'   \item{country}{Country name.}
#'   \item{region}{World Bank region.}
#'   \item{income_level}{World Bank income group.}
#'   \item{year}{Calendar year, 2010--2024.}
#'   \item{gdp_growth}{Annual GDP growth, percent.}
#'   \item{inflation}{Annual inflation, percent.}
#'   \item{fdi_gdp}{Foreign direct investment net inflows, percent of GDP.}
#'   \item{unemployment}{Total unemployment, percent of labor force.}
#'   \item{internet_users}{Individuals using the Internet, percent of population.}
#'   \item{school_secondary}{Gross secondary-school enrollment, percent.}
#'   \item{log_population}{Natural logarithm of total population.}
#'   \item{log_gdp_pc_const}{Natural logarithm of constant-price GDP per capita.}
#'   \item{trade_open}{Trade openness, percent of GDP.}
#'   \item{z0_lorentz}{Positive time coordinate of the Lorentz representation.}
#'   \item{z1_lorentz}{First economic feature coordinate on the hyperboloid.}
#'   \item{z2_lorentz}{Second economic feature coordinate on the hyperboloid.}
#'   \item{hyperboloid_error}{Numerical check of the hyperboloid constraint plus one.}
#' }
#'
#' @source [World Bank World Development Indicators](https://datacatalog.worldbank.org/search/dataset/0037712/world-development-indicators),
#' downloaded through the API in July 2026 and distributed under the
#' [Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/).
#' Indicators include `NY.GDP.MKTP.KD.ZG`, `NY.GDP.PCAP.KD`,
#' `NE.TRD.GNFS.ZS`, `BX.KLT.DINV.WD.GD.ZS`, `FP.CPI.TOTL.ZG`,
#' `SP.POP.TOTL`, `SL.UEM.TOTL.ZS`, `IT.NET.USER.ZS`, and `SE.SEC.ENRR`.
#'
#' @examples
#' data(wdi_lhwr)
#' dat <- subset(wdi_lhwr, year == max(year))
#'
#' # The package stores the constructed Lorentz coordinates directly.
#' Z <- as.matrix(dat[, c("z0_lorentz", "z1_lorentz", "z2_lorentz")])
#' check_hyperboloid(Z, stop_on_fail = TRUE)
#'
#' # Reproduce the stored coordinates with the frozen preprocessing constants.
#' z1 <- (dat$log_gdp_pc_const - 8.8071974413888245) / 1.4524958189486741
#' z2 <- (dat$trade_open - 91.8784599801564923) / 58.6130974219034400
#' Z_rebuilt <- cbind(sqrt(1 + z1^2 + z2^2), z1, z2)
#' max(abs(Z - Z_rebuilt))
#'
#' # Keep examples lightweight; use a modest grid.
#' fit_data <- data.frame(
#'   gdp_growth = dat$gdp_growth,
#'   unemployment = as.numeric(scale(dat$unemployment)),
#'   internet_users = as.numeric(scale(dat$internet_users))
#' )
#'
#' bw <- bw.lhwr(
#'   gdp_growth ~ unemployment + internet_users,
#'   data = fit_data,
#'   Z = Z,
#'   approach = "AICc",
#'   bw_grid = c(40, 60),
#'   engine = "knn"
#' )
#' fit <- lhwr(
#'   gdp_growth ~ unemployment + internet_users,
#'   data = fit_data,
#'   Z = Z,
#'   bw = bw$bw,
#'   engine = "knn",
#'   diagnostics = "aicc"
#' )
#' fit
NULL
