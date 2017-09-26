library(xSom)
library(som)
context("Extended SOM")


# mit naExist und ohne
l <- 50
xDim <- 10
somSize <- 3
x <- seq(-1,1,len=l)
y <- sin((6*x-3)^2)
data <- cbind(x,y)
W <- som.init.extended(data, somSize, c(T,T))

test_that("Gleichheit zwischen som und xSom (in R)", {
  old <- som(data, xdim = somSize, ydim = somSize, alphaType="linear", neigh="gaussian"  )$code
  new <- somCheckNaR(dataSet= data, weightMatrix=W,  oldColumns=c(T,T),
                     updateParametersPerEpoch = F, sampling=0, naExist=F)
  attr(old, "class") <- NULL

  expect_equal(
    new, old
  )
})

test_that("Gleichheit zwischen som und xSom (in C++)", {
  is <- somCheckNa(dataSet= data, weightMatrix=W,  oldColumns=c(T,T), sampling=0, updateParametersPerEpoch = F, naExist = F)
  expected <- somCheckNaR(dataSet= data, weightMatrix=W,  oldColumns=c(T,T), sampling=0, updateParametersPerEpoch = F,naExist = F)

  expect_equal(
    is, expected
  )
})

test_that("Gleichheit zwischen xSom (C++) und xSom (in R), Update pro Epoche", {
  is <- somCheckNa(dataSet= data, weightMatrix=W,  oldColumns=c(T,T), sampling=0, updateParametersPerEpoch = T, naExist = F)
  expected <- somCheckNaR(dataSet= data, weightMatrix=W,  oldColumns=c(T,T), sampling=0, updateParametersPerEpoch = T, naExist = F)

  expect_equal(
    is, expected
  )
})

test_that("Gleichheit zwischen xSom (C++) und xSom (in R), Update pro Datensatz", {
  is <- somCheckNa(dataSet= data, weightMatrix=W,  oldColumns=c(T,T), sampling=0, updateParametersPerEpoch = F, naExist = F)
  expected <- somCheckNaR(dataSet= data, weightMatrix=W,  oldColumns=c(T,T), sampling=0, updateParametersPerEpoch = F, naExist = F)

  expect_equal(
    is, expected
  )
})

# Beispiel mit inverser Lernrate
test_that("Gleichheit zwischen xSom (C++) und som, inverseLernRate", {
  old <- som(data, xdim = somSize, ydim = somSize, alphaType="inverse", neigh="gaussian", inv.alp.c = c(7,7)  )$code
  new <- somCheckNa(dataSet= data, weightMatrix=W,  oldColumns=c(T,T),
                     updateParametersPerEpoch = F, sampling=0, naExist=F, learnRateReduction = -7)
  attr(old, "class") <- NULL

  expect_equal(
    new, old
  )
})

test_that("Gleichheit zwischen xSom (R) und som, inverseLernRate", {
  old <- som(data, xdim = somSize, ydim = somSize, alphaType="inverse", neigh="gaussian", inv.alp.c = c(7,7)  )$code
  new <- somCheckNaR(dataSet= data, weightMatrix=W,  oldColumns=c(T,T),
                    updateParametersPerEpoch = F, sampling=0, naExist=F, learnRateReduction = -7)
  attr(old, "class") <- NULL

  expect_equal(
    new, old
  )
})
