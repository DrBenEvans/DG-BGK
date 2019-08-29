all = BGK

#FFLAGS= -fc=ifort -heap-arrays -g -traceback -check all 
#FFLAGS= -fc=ifort -heap-arrays -check all 
#FFLAGS= -fc=ifort -heap-arrays -O2
FFLAGS= -g -Wall 

OBJ = addvec.o adthem.o advnce.o alotim2.o alotim.o bcastp.o \
      bcastv.o cmminc.o compmu.o conmat.o dg_lin_convec.o edgcom.o \
      edgflxopt.o elnode.o flucon2.o getavr.o getbou.o getelc.o \
      getfla.o getflu.o getfor.o getgeo.o getinc.o getinf.o \
      getlen.o getloc.o getmac.o getmat.o getnor.o getres.o \
      getrho.o getrow.o getsid.o gteqnf.o gtinpt.o identa.o \
      identm.o ifilla.o ifillm.o ifillv.o inflow.o inicon2.o \
      inicon.o macros.o normal.o output.o ptmesh.o rfilla.o \
      rfillm.o rfillv.o spectralbgk.o vspace.o input_variable_module.o

BGK : $(OBJ)
	mpif90 $(FFLAGS) -o $@ $^

advnce.o : advnce.f90 input_variable_module.mod
	mpif90 $(FFLAGS) -c $<

ouptut.o : output.f90 input_variable_module.mod
	mpif90 $(FFLAGS) -c $<

ptmesh.o : ptmesh.f90 input_variable_module.mod getloc_module.mod
	mpif90 $(FFLAGS) -c $<

spectralbgk.o : spectralbgk.f90 dg_lin_convec_module.mod  \
	            input_variable_module.mod

advnce_module.mod : advnce.o

output_module.mod : output.o

dg_lin_convec.o : dg_lin_convec.f90 advnce_module.mod output_module.mod \
	                 ptmesh_module.mod
	mpif90 $(FFLAGS) -c $<

%.mod : %.o

%.o : %.f90 
	mpif90 $(FFLAGS) -c $<

clean:
	rm -rf BGK $(OBJ)
