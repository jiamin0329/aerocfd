#file name
program = ../cfd.out

#complier
fc = mpif90
fflag = -O2 -g -frecord-marker=4 -fcheck=bounds -Dmpi ${DEBUG}
##fflag = -O3 -frecord-marker=4 -Dmpi ${DEBUG}
#*********************************************************************************************#
#         -fcheck=all                           link objection files                                      #
#*********************************************************************************************#
obj_modules = ./obj/module_block.o     \
              ./obj/module_flags.o     \
              ./obj/module_glbindex.o  \
              ./obj/module_index.o     \
              ./obj/module_mpi.o       \
              ./obj/module_ns_consts.o \
              ./obj/module_output.o    \
              ./obj/module_tursa.o

obj_main = ./obj/main.o

obj_bc = ./obj/sub_bc.o                \
         ./obj/sub_bc_-10_farfield.o   \
         ./obj/sub_bc_-20_solidwall.o  \
         ./obj/sub_bc_-30_symm.o       \
         ./obj/sub_bc_-40_outflow.o    \
         ./obj/sub_bc_-50_inflow.o     \
         ./obj/sub_bc_updatebufferaverage2d.o \
         ./obj/sub_bc_updatebufferaverage.o \
         ./obj/sub_bc_updatebuffer2d.o \
         ./obj/sub_bc_updatebuffer.o \
         ./obj/sub_bc_updatebuffercoordinate2d.o \
         ./obj/sub_bc_updatebuffercoordinate.o \
         ./obj/sub_bc_updatebufferjacobian2d1.o \
         ./obj/sub_bc_updatebufferjacobian1.o \
         ./obj/sub_bc_updatebufferjacobian2d2.o \
         ./obj/sub_bc_updatebufferjacobian2.o \
         ./obj/sub_bc_updatebufferjacobian2d3.o \
         ./obj/sub_bc_updatebufferjacobian3.o
       

obj_init = ./obj/sub_init_allocatememory.o     \
           ./obj/sub_init_flowfield.o          \
           ./obj/sub_init_gridtransformation.o \
           ./obj/sub_init_indextransfer.o      \
           ./obj/sub_init_readbcinfo.o         \
           ./obj/sub_init_readcontrolfile.o    \
           ./obj/sub_init_walldist.o           \
           ./obj/sub_init_spacing.o
           
obj_inviscid = ./obj/sub_inviscid1.o                \
               ./obj/sub_inviscid_splitting.o       \
               ./obj/sub_inviscid_upwind.o
               
obj_output = ./obj/sub_output_checking.o           \
             ./obj/sub_output_flowfield.o          \
             ./obj/sub_output_forceintegration.o   \
             ./obj/sub_output_geom.o               \
             ./obj/sub_output_jacobian.o           \
             ./obj/sub_output_residual.o           \
             ./obj/sub_output_surface.o            \
             ./obj/sub_output_acousticfile.o
              
obj_time = ./obj/sub_timeadvancement_1stimplicit.o  \
           ./obj/sub_timeadvancement_finddt.o       \
           ./obj/sub_timeadvancement_jac.o

obj_pri = ./obj/sub_getpri_v.o
          
obj_vis = ./obj/sub_viscous.o  \
          ./obj/sub_secondordercentral.o

obj_tur = ./obj/sub_tur_sa.o
 
objs = ${obj_modules}  \
       ${obj_main}     \
       ${obj_bc}       \
       ${obj_init}     \
       ${obj_inviscid} \
       ${obj_time}     \
       ${obj_output}   \
       ${obj_vis}      \
       ${obj_pri}      \
       ${obj_tur}

all:	$(objs)
	$(fc) -o $(program) $(objs)
#*********************************************************************************************#
#                               end objection files                                           #
#*********************************************************************************************# 

#*********************************************************************************************#
#                             generate objective files                                        #
#*********************************************************************************************#
#main
./obj/main.o: main.f90            \
              module_block.f90    \
              module_flags.f90    \
              module_glbindex.f90 \
              module_index.f90    \
              module_mpi.f90      \
              module_ns_consts.f90\
              module_output.f90
	$(fc) -c $(fflag) -o $@ $<

#obj_bc
./obj/sub_bc.o: sub_bc.f90           \
                module_ns_consts.f90 \
                module_block.f90     \
                module_index.f90     \
                module_mpi.f90       \
                module_flags.f90
	$(fc) -c $(fflag) -o $@ $<
./obj/sub_bc_-10_farfield.o: sub_bc_-10_farfield.f90
	$(fc) -c $(fflag) -o $@ $<	
./obj/sub_bc_-20_solidwall.o: sub_bc_-20_solidwall.f90
	$(fc) -c $(fflag) -o $@ $<	
./obj/sub_bc_-30_symm.o: sub_bc_-30_symm.f90
	$(fc) -c $(fflag) -o $@ $<	
./obj/sub_bc_-40_outflow.o: sub_bc_-40_outflow.f90
	$(fc) -c $(fflag) -o $@ $<	
./obj/sub_bc_-50_inflow.o: sub_bc_-50_inflow.f90
	$(fc) -c $(fflag) -o $@ $<		
./obj/sub_bc_sub_bc_updatebufferaverage2d.o: sub_bc_sub_bc_updatebufferaverage2d.f90  \
                       module_block.f90   \
                       module_index.f90   \
                       module_mpi.f90     \
                       module_glbindex.f90
	$(fc) -c $(fflag) -o $@ $<

./obj/sub_bc_sub_bc_updatebufferaverage.o: sub_bc_sub_bc_updatebufferaverage.f90  \
                       module_block.f90   \
                       module_index.f90   \
                       module_mpi.f90     \
                       module_glbindex.f90
	$(fc) -c $(fflag) -o $@ $<

./obj/sub_bc_updatebuffer2d.o: sub_bc_updatebuffer2d.f90 \
                               module_block.f90          \
                               module_index.f90          \
                               module_mpi.f90            \
                               module_glbindex.f90
	$(fc) -c $(fflag) -o $@ $<

./obj/sub_bc_updatebuffer.o: sub_bc_updatebuffer.f90 \
                               module_block.f90          \
                               module_index.f90          \
                               module_mpi.f90            \
                               module_glbindex.f90
	$(fc) -c $(fflag) -o $@ $<

./obj/sub_bc_updatebuffercoordinate2d.o: sub_bc_updatebuffercoordinate2d.f90 \
                               module_block.f90          \
                               module_index.f90          \
                               module_mpi.f90            \
                               module_glbindex.f90
	$(fc) -c $(fflag) -o $@ $<

./obj/sub_bc_updatebuffercoordinate.o: sub_bc_updatebuffercoordinate.f90 \
                               module_block.f90          \
                               module_index.f90          \
                               module_mpi.f90            \
                               module_glbindex.f90
	$(fc) -c $(fflag) -o $@ $<
	
./obj/sub_bc_updatebufferjacobian2d1.o: sub_bc_updatebufferjacobian2d1.f90 \
                               module_block.f90          \
                               module_index.f90          \
                               module_mpi.f90            \
                               module_glbindex.f90
	$(fc) -c $(fflag) -o $@ $<

./obj/sub_bc_updatebufferjacobian1.o: sub_bc_updatebufferjacobian1.f90 \
                               module_block.f90          \
                               module_index.f90          \
                               module_mpi.f90            \
                               module_glbindex.f90
	$(fc) -c $(fflag) -o $@ $<

./obj/sub_bc_updatebufferjacobian2d2.o: sub_bc_updatebufferjacobian2d2.f90  \
                             module_block.f90         \
                             module_index.f90         \
                             module_mpi.f90           \
                             module_glbindex.f90
	$(fc) -c $(fflag) -o $@ $<	

./obj/sub_bc_updatebufferjacobian2.o: sub_bc_updatebufferjacobian2.f90  \
                             module_block.f90         \
                             module_index.f90         \
                             module_mpi.f90           \
                             module_glbindex.f90
	$(fc) -c $(fflag) -o $@ $<	

./obj/sub_bc_updatebufferjacobian2d3.o: sub_bc_updatebufferjacobian2d3.f90  \
                             module_block.f90         \
                             module_index.f90         \
                             module_mpi.f90           \
                             module_glbindex.f90
	$(fc) -c $(fflag) -o $@ $<	

./obj/sub_bc_updatebufferjacobian3.o: sub_bc_updatebufferjacobian3.f90  \
                             module_block.f90         \
                             module_index.f90         \
                             module_mpi.f90           \
                             module_glbindex.f90
	$(fc) -c $(fflag) -o $@ $<	

#obj_init 
./obj/sub_init_allocatememory.o: sub_init_allocatememory.f90 \
                                 module_block.f90            \
                                 module_flags.f90            \
                                 module_glbindex.f90         \
                                 module_index.f90            \
                                 module_mpi.f90              \
                                 module_ns_consts.f90        \
                                 module_output.f90
	$(fc) -c $(fflag) -o $@ $<	
./obj/sub_init_flowfield.o: sub_init_flowfield.f90 \
                            module_block.f90       \
                            module_flags.f90       \
                            module_index.f90       \
                            module_mpi.f90         \
                            module_ns_consts.f90 
	$(fc) -c $(fflag) -o $@ $<	
./obj/sub_init_gridtransformation.o: sub_init_gridtransformation.f90 \
                                     module_flags.f90
	$(fc) -c $(fflag) -o $@ $<	                   
./obj/sub_init_indextransfer.o: sub_init_indextransfer.f90
	$(fc) -c $(fflag) -o $@ $<	                    
./obj/sub_init_readbcinfo.o: sub_init_readbcinfo.f90  \
                             module_block.f90         \
                             module_glbindex.f90      \
                             module_index.f90         \
                             module_mpi.f90           \
                             module_output.f90
	$(fc) -c $(fflag) -o $@ $<	                     
./obj/sub_init_readcontrolfile.o: sub_init_readcontrolfile.f90 \
                                  module_flags.f90             \
                                  module_mpi.f90               \
                                  module_ns_consts.f90         \
                                  module_output.f90
	$(fc) -c $(fflag) -o $@ $<	                        
./obj/sub_init_walldist.o: sub_init_walldist.f90 \
                           module_block.f90      \
                           module_flags.f90      \
                           module_glbindex.f90   \
                           module_index.f90      \
                           module_mpi.f90
	$(fc) -c $(fflag) -o $@ $<	
./obj/sub_init_spacing.o: sub_init_spacing.f90
	$(fc) -c $(fflag) -o $@ $<		

#obj_inviscid
./obj/sub_inviscid1.o: sub_inviscid1.f90 \
                       module_flags.f90
	$(fc) -c $(fflag) -o $@ $<	 
./obj/sub_inviscid_splitting.o: sub_inviscid_splitting.f90
	$(fc) -c $(fflag) -o $@ $<	  
./obj/sub_inviscid_upwind.o: sub_inviscid_upwind.f90
	$(fc) -c $(fflag) -o $@ $<	

#obj_output
./obj/sub_output_checking.o: sub_output_checking.f90 \
                             module_block.f90        \
                             module_index.f90        \
                             module_mpi.f90
	$(fc) -c $(fflag) -o $@ $<	
./obj/sub_output_flowfield.o: sub_output_flowfield.f90 \
                              module_block.f90         \
                              module_flags.f90         \
                              module_index.f90         \
                              module_mpi.f90           \
                              module_ns_consts.f90     \
                              module_output.f90
	$(fc) -c $(fflag) -o $@ $<	
./obj/sub_output_forceintegration.o: sub_output_forceintegration.f90 \
                                     module_block.f90                \
                                     module_flags.f90                \
                                     module_index.f90                \
                                     module_mpi.f90                  \
                                     module_ns_consts.f90
	$(fc) -c $(fflag) -o $@ $<		
./obj/sub_output_geom.o: sub_output_geom.f90 \
                         module_block.f90    \
                         module_index.f90    \
                         module_mpi.f90
	$(fc) -c $(fflag) -o $@ $<		
./obj/sub_output_jacobian.o: sub_output_jacobian.f90 \
                             module_block.f90        \
                             module_index.f90        \
                             module_mpi.f90
	$(fc) -c $(fflag) -o $@ $<	
./obj/sub_output_residual.o: sub_output_residual.f90 \
                             module_block.f90        \
                             module_index.f90        \
                             module_mpi.f90          \
                             module_flags.f90
	$(fc) -c $(fflag) -o $@ $<	
	
./obj/sub_output_surface.o: sub_output_surface.f90   \
                            module_block.f90         \
                            module_flags.f90         \
                            module_index.f90         \
                            module_mpi.f90           \
                            module_ns_consts.f90     \
                            module_output.f90
	$(fc) -c $(fflag) -o $@ $<
./obj/sub_output_acousticfile.o: sub_output_acousticfile.f90   \
                            module_block.f90         \
                            module_flags.f90         \
                            module_index.f90         \
                            module_mpi.f90           \
                            module_ns_consts.f90     \
                            module_output.f90
	$(fc) -c $(fflag) -o $@ $<	

#obj_time
./obj/sub_timeadvancement_1stimplicit.o: sub_timeadvancement_1stimplicit.f90 \
                                         module_flags.f90
	$(fc) -c $(fflag) -o $@ $<	      
./obj/sub_timeadvancement_finddt.o: sub_timeadvancement_finddt.f90
	$(fc) -c $(fflag) -o $@ $<	
./obj/sub_timeadvancement_jac.o: sub_timeadvancement_jac.f90
	$(fc) -c $(fflag) -o $@ $<	
	
#obj_pri
./obj/sub_getpri_v.o: sub_getpri_v.f90
	$(fc) -c $(fflag) -o $@ $<	 
       
#obj_vis
./obj/sub_viscous.o: sub_viscous.f90 \
                     module_flags.f90
	$(fc) -c $(fflag) -o $@ $<	 

./obj/sub_secondordercentral.o: sub_secondordercentral.f90
	$(fc) -c $(fflag) -o $@ $<	 	

./obj/sub_tur_sa.o: sub_tur_sa.f90
	$(fc) -c $(fflag) -o $@ $<	 
	
#modules
./obj/%.o: %.f90                \
		   module_block.f90     \
		   module_flags.f90     \
           module_glbindex.f90  \
           module_index.f90     \
           module_mpi.f90       \
           module_ns_consts.f90 \
           module_output.f90    \
           module_tursa.f90
	$(fc) -c $(fflag) -o $@ $<
#*********************************************************************************************#
#                              end generate objective files                                   #
#*********************************************************************************************#	
clean:
	rm -f ./obj/*.o ./obj/*.mod ./obj/*.l
	rm -f       *.o       *.mod       *.l
cleanb:
	rm -f *~     
