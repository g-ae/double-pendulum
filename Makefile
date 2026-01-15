.PHONY: optimise

optimise:
	julia optimisation_pendule.jl
	montage masse_a_x.png masse_a_y.png masse_b_x.png masse_b_y.png -tile 2x2 -geometry +0+0 -background white resultat_combine.png
	rm masse_*.png