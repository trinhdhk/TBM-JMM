# TBM-JMM: Joint Mixed-Effects Model implementation for TBM brain morphometry

This repo implements a proof-of-concept joint model to investigate the longitudinal trajectory of TBM brain subcortical volume, contrasted under treatment.

## Motivation

Conventional joint models assume longitudinal processes are independent given some limited shared random effects.
They also assumes all the longitudinal processes associate with the survival process independently.

The may not be true for highly correlated features like brain region. 
This model adopts Gelman[2013]'s proposal to correct for multiple comparisions using shrinkage from mixed effects model, 
but it a joint modelling framework.

## Technical information

The model modifies some brms code, currently is quite hacky. 
This model comprises of 73 longitudinal processes fitted together by a shared set of fixed effects and process-specific random effects.

These effects are isolates and fed into the survival process, which for now is a pool logistic regression.

## TODO

Package development