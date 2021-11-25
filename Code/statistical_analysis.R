# All functions dealing with significance tests

# Perform a Fisher test to compare the variances of distribution_1 and distribution_2.
F_test <- function(distribution_1, distribution_2) { # TODO
  ?var.test()
  result = var.test(distribution_1, distribution_2)
  print(result)
}
