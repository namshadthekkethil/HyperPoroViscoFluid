#####################################################################
## Here specify the location of the IBAMR source and the location
##IBAMR_SRC_DIR = /xlwork6/heartvalve/IBAMR-2015/sfw/IBAMR
##IBAMR_BUILD_DIR = /xlwork6/heartvalve/IBAMR-2015/sfw/IBAMR/ibamr-objs-opt
##IBAMR_SRC_DIR = /xlwork6/2101013f/LiuyangIBAMR/IBAMR
##IBAMR_BUILD_DIR = /xlwork6/2101013f/LiuyangIBAMR/IBAMR/ibamr-objs-opt
##IBAMR_SRC_DIR = /xlwork6/2101013f/LiuyangIBAMR/LAE/IBAMR
##IBAMR_BUILD_DIR = /xlwork6/2101013f/LiuyangIBAMR/LAE/IBAMR/ibamr-objs-opt
#IBAMR_SRC_DIR =/work/e645/shared/sfw/ibamr/IBAMR
#IBAMR_BUILD_DIR =/work/e645/shared/sfw/ibamr/ibamr-objs-opt
IBAMR_SRC_DIR =/work/e645/shared/sfw2/ibamr_latest/IBAMR
IBAMR_BUILD_DIR =/work/e645/shared/sfw2/ibamr_latest/ibamr-objs-opt
######################################################################
## Include variables specific to the particular IBAMR build.
include $(IBAMR_BUILD_DIR)/config/make.inc


SRC = $(wildcard *.C) $(wildcard /work/e642/e642/namshadth/source_codes/FEMLDE_7/*.C)
PDIM = 3
OBJS = $(SRC:%.C=%.o) $(IBAMR_LIB_3D) $(IBTK_LIB_3D)


main: $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -DNDIM=$(PDIM) -o main



clean:
	$(RM) main
	$(RM) *.o *.lo *.objs *.ii *.int.c *.mod
	$(RM) -r .libs

-include $(SRC:%.C=%.d)
