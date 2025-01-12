# biofilter_toc ----

test_that("biofilter_toc returns an error when water is absent or input incorrectly.", {
  expect_error(biofilter_toc(ebct = 10, ozonated = TRUE))
  expect_error(biofilter_toc(water = list(ph = 7, toc = 5), ebct = 10, ozonated = TRUE))
})

test_that("biofilter_toc returns an error or warning when arguments are input improperly or missing.", {
  water <- suppressWarnings(define_water(ph = 7, temp = 25, alk = 100, toc = 5.0, doc = 4.0, uv254 = 0.1))

  expect_error(biofilter_toc(water, ozonated = TRUE))
  expect_error(biofilter_toc(water, ebct = "4", ozonated = FALSE))

  expect_error(biofilter_toc(water, ebct = 4, ozonated = 2))
  expect_error(biofilter_toc(water, ebct = 4, ozonated = "TRUE"))
})

test_that("biofilter_toc returns an error when TOC is missing.", {
  water_no_toc <- suppressWarnings(define_water(ph = 7, temp = 15, alk = 100))
  expect_error(biofilter_toc(water_no_toc, ebct = 10, ozonated = TRUE))
})

test_that("biofilter_toc calculates correct TOC removal for non-ozonated water.", {
  water <- suppressWarnings(define_water(ph = 7, temp = 15, alk = 100, toc = 5.0, doc = 4.0, uv254 = 0.1))
  dosed_water <- suppressWarnings(biofilter_toc(water, ebct = 10, ozonated = FALSE))

  # Expect that TOC is reduced correctly using non-ozonated parameters
  expect_equal(round(dosed_water@toc, 2), 4.53) # Expected TOC after treatment
  expect_equal(round(dosed_water@doc, 2), 4.05) # Expected DOC (BDOC fraction of TOC)
})

test_that("biofilter_toc calculates correct TOC removal for ozonated water.", {
  water <- suppressWarnings(define_water(ph = 7, temp = 15, alk = 100, toc = 5.0, doc = 4.0, uv254 = 0.1))
  dosed_water <- suppressWarnings(biofilter_toc(water, ebct = 10, ozonated = TRUE))

  # Expect that TOC is reduced correctly using ozonated parameters
  expect_equal(round(dosed_water@toc, 2), 4.46) # Expected TOC after treatment
  expect_equal(round(dosed_water@doc, 2), 3.92) # Expected DOC (BDOC fraction of TOC)
})

test_that("biofilter_toc correctly handles temperatures and non-ozonated water.", {
  water <- suppressWarnings(define_water(ph = 7, temp = 45, alk = 100, toc = 5.0, doc = 4.0, uv254 = 0.1))

  # the Bad Place temperature, non-ozonated
  dosed_water_high <- suppressWarnings(biofilter_toc(water, ebct = 10, ozonated = FALSE))
  expect_equal(round(dosed_water_high@toc, 2), 4.47)

  # the Medium Place temperature, non-ozonated
  water@temp <- 19
  dosed_water_med <- suppressWarnings(biofilter_toc(water, ebct = 10, ozonated = FALSE))
  expect_equal(round(dosed_water_med@toc, 2), 4.53)

  # the Good Place temperature, non-ozonated
  water@temp <- 5
  dosed_water_low <- suppressWarnings(biofilter_toc(water, ebct = 10, ozonated = FALSE))
  expect_equal(round(dosed_water_low@toc, 2), 4.79)
})

test_that("biofilter_toc correctly handles temperatures and ozonated water.", {
  water <- suppressWarnings(define_water(ph = 7, temp = 45, alk = 100, toc = 5.0, doc = 4.0, uv254 = 0.1))

  dosed_water_high <- suppressWarnings(biofilter_toc(water, ebct = 10, ozonated = TRUE))
  expect_equal(round(dosed_water_high@toc, 2), 4.07)

  water@temp <- 19
  dosed_water_med <- suppressWarnings(biofilter_toc(water, ebct = 10, ozonated = TRUE))
  expect_equal(round(dosed_water_med@toc, 2), 4.46)

  water@temp <- 5
  dosed_water_low <- suppressWarnings(biofilter_toc(water, ebct = 10, ozonated = TRUE))
  expect_equal(round(dosed_water_low@toc, 2), 4.69)
})
