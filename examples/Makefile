
all: ../images/pv-isotherms.svg ../images/Ts-isobars.svg

../images/%.svg: %.py
	mkdir -p ../images
	./$< -o $@
