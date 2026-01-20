test_that("Rd examples run cleanly", {
  examples <- tools::testInstalledPackage("depcoeff", types = "examples")
  succeed()
})
