

FC  = ifort
OPT = -O0 -g -traceback -check all -check bounds -debug all 

mdi_test : module_precision.o ModelDataInstances_Class.o ModelDataInstances_Driver.o	
	${FC} ${OPT} module_precision.o ModelDataInstances_Class.o ModelDataInstances_Driver.o -o $@
     

module_precision.o : ../../main/module_precision.f90
	${FC} ${OPT} -c ../../main/module_precision.f90 -o $@

ModelDataInstances_Class.o : ../ModelDataInstances_Class.f90 module_precision.o
	${FC} ${OPT} -c ../ModelDataInstances_Class.f90 -o $@

ModelDataInstances_Driver.o : ModelDataInstances_Driver.f90 ModelDataInstances_Class.o module_precision.o
	${FC} ${OPT} -c ModelDataInstances_Driver.f90 -o $@
