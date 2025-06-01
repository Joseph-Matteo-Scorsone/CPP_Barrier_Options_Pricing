# Barrier Option Pricing in C++

This project implements a pricing engine for **barrier options** using the Black-Scholes model in C++.

## Features

- Supports four types of barrier options:
  - Up-and-Out
  - Up-and-In
  - Down-and-Out
  - Down-and-In
- Also includes vanilla call and put pricing.
- Partial support for rebate payouts upon knockout.
- Uses analytical formulas for barrier option pricing (where implemented).

## Barrier Options

Barrier options are exotic options where the payoff depends on whether the underlying asset reaches a certain barrier level during the option's life.

## Parameters

Each option is defined by the following parameters:

- `S`: Spot price
- `K`: Strike price
- `H`: Barrier level
- `r`: Risk-free rate
- `q`: Dividend yield
- `sigma`: Volatility
- `T`: Time to maturity
- `OptionType`: Call or Put
- `BarrierType`: UpAndOut, DownAndOut, UpAndIn, DownAndIn
- `rebate`: Amount paid if knocked out (default is 0)

Relevant papers
- https://people.maths.ox.ac.uk/howison/barriers.pdf
- https://arxiv.org/pdf/cond-mat/0211489
