# SPELT: Species Pairs Evolutionary Lag Test (SPELT)

This repository contains functions needed to use SPELT.

## Installing SPELT

You can install directly from GitHub if you have the devtools package installed:

	library(devtools)
	install_github("SPELT", username = "nhcooper123")
	library(SPELT)

Magical!

## Using SPELT

The package has lots of internal functions with help files but the only function you really need is SPELT:

	data(shorebird)
	SPELT.results <- SPELT(shorebird.tree, shorebird.data,
                         "F.Mass", "Egg.Mass", "Species")
	summary(SPELT.results)
	plot(SPELT.results)
	str(SPELT.results$data)