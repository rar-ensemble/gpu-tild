# -rdc=true   FLAG for remote cu code
CC = nvcc
CFLAGS = -std=c++11 -rdc=true -arch=sm_37 -dlto -O3 -I/usr/local/cuda-10.2/include -Wno-deprecated-gpu-targets
#CFLAGS = -std=c++11 -rdc=true -arch=sm_37 -dlto -O3 -I/usr/local/cuda-10.2/include -lineinfo -g #delete lineinfo later and -g
LIBS = -lm -lcudart -lcufft -dlto -L/usr/local/cuda-10.2/lib


#############################################################################
# nothing should be changed below here
       
SRCS = array_utils.cu bonds.cu calc_properties.cu config_utils.cu cuda_random_posits.cu \
	device_EM_integrator.cu device_array_utils.cu device_bonds.cu device_comm_utils.cu \
	device_grid_utils.cu device_utils.cu die.cu field_component.cu forces.cu \
	initialize.cu integ_utils.cu io_utils.cu main.cu pair_style.cu \
	pair_style_fieldphases.cu pair_style_gaussian.cu pbc_utils.cu read_input.cu \
	pair_style_erf.cu reduce_utils.cu update_pairstyles.cu device_GJF_integrator.cu \
	pair_style_charges.cu pair_style_gaussian_erf.cu group.cu integrator.cu \
	nlist.cu Compute.cu

OBJS = $(SRCS:.cu=.o)
   
%.o: %.cu
	$(CC) $(CFLAGS) -c $< -o $@
	 
gpu-tild: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBS)

#OBJS = $(SRCS:.cu=.o)
#
#.cu.o: .cu
#	${CC} ${CFLAGS} -c  $<
#
#gpu-tild:  ${OBJS}
#	$(CC) ${CFLAGS} -o $@ ${OBJS} $(LIBS)

clean:
	rm -f *.o
	rm -f gpu-tild
	rm -f *~

