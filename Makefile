include Makefile.inc
	
synccosma: clean
	rsync -avzL $(shell pwd)/ jvbq85@cosma-a:data/HBT2/code