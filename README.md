# SPELT: Species Pairs Evolutionary Lag Test (SPELT)

[![Build Status](https://travis-ci.org/nhcooper123/SPELT.png?branch=master)](https://travis-ci.org/nhcooper123/SPELT)

[![Build Status](https://zenodo.org/badge/4008/nhcooper123/SPELT.png)](https://github.com/nhcooper123/SPELT/releases)

This repository contains functions needed to use SPELT. 

## Installing SPELT

You can install directly from GitHub if you have the devtools package installed:

	library(devtools)
	install_github("SPELT", username = "nhcooper123")
	library(SPELT)

Totally magical!

## Using SPELT

The package has lots of internal functions but the only function you really need is SPELT:

	data(shorebird, package = "caper")
	SPELT.results <- SPELT(shorebird.tree, shorebird.data,
                         "F.Mass", "Egg.Mass", "Species")
	summary(SPELT.results)
	plot(SPELT.results)
	str(SPELT.results$data)

Check out ?SPELT in R for more details.
