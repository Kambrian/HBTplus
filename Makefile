include Makefile.base

synccosma: clean
	rsync -avzL $(shell pwd)/ jvbq85@cosma-a:data/HBT2/code