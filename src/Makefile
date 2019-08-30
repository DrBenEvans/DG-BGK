all = BGK

#FFLAGS= -fc=ifort -heap-arrays -g -traceback -check all 
#FFLAGS= -fc=ifort -heap-arrays -check all 
#FFLAGS= -fc=ifort -heap-arrays -O2
#FFLAGS= -g -Wall -O#
FFLAGS= -Wall -O2 -mavx2

# getavr.o was excluded since GETAVR is not referenced anywhere in
# the code

OBJ = addvec.o adthem.o advnce.o alotim2.o alotim.o bcastp.o \
      bcastv.o cmminc.o compmu.o conmat.o dg_lin_convec.o edgcom.o \
      edgflxopt.o elnode.o flucon2.o getbou.o getelc.o \
      getfla.o getflu.o getfor.o getgeo.o getinc.o getinf.o \
      getlen.o getloc.o getmac.o getmat.o getnor.o getres.o \
      getrho.o getrow.o getsid.o gteqnf.o gtinpt.o identa.o \
      identm.o ifilla.o ifillm.o ifillv.o inflow.o inicon2.o \
      inicon.o macros.o normal.o output.o ptmesh.o rfilla.o \
      rfillm.o rfillv.o spectralbgk.o vspace.o input_variable_module.o

BGK : $(OBJ)
	mpif90 $(FFLAGS) -o $@ $^

advnce.o : advnce.f90 edgflxopt_module.mod \
	         adthem_module.mod \
	         cmminc_module.mod \
	         compmu_module.mod \
	         flucon2_module.mod \
	         getelc_module.mod \
	         getfla_module.mod \
	         getinc_module.mod \
	         getmac_module.mod \
					 getres_module.mod \
	         inflow_module.mod \
					 input_variable_module.mod macros_module.mod
	mpif90 $(FFLAGS) -c $<

edgflxopt.o : edgflxopt.f90 getbou_module.mod
	mpif90 $(FFLAGS) -c $<

inicon.o : inicon.f90 input_variable_module.mod
	mpif90 $(FFLAGS) -c $<

inicon2.o : inicon2.f90 input_variable_module.mod
	mpif90 $(FFLAGS) -c $<

ouptut.o : output.f90 input_variable_module.mod
	mpif90 $(FFLAGS) -c $<

ptmesh.o : ptmesh.f90 input_variable_module.mod getloc_module.mod
	mpif90 $(FFLAGS) -c $<

dg_lin_convec.o : dg_lin_convec.f90 advnce_module.mod getnor_module.mod \
				  inicon2_module.mod output_module.mod ptmesh_module.mod
	mpif90 $(FFLAGS) -c $<

spectralbgk.o : spectralbgk.f90 dg_lin_convec_module.mod  \
	            gtinpt_module.mod input_variable_module.mod

advnce_module.mod : advnce.o

dg_lin_convec_module.mod : dg_lin_convec.o

edgflxopt_module.mod : edgflxopt.o

adthem_module.mod : adthem.o

cmminc_module.mod : cmminc.o

compmu_module.mod : compmu.o

flucon2_module.mod : flucon2.o

getbou_module.mod : getbou.o

getelc_module.mod : getelc.o

getfla_module.mod : getfla.o

getinc_module.mod : getinc.o

getloc_module.mod : getloc.o

getmac_module.mod : getmac.o

getres_module.mod : getres.o

getnor_module.mod : getnor.o

gtinpt_module.mod : gtinpt.o

inflow_module.mod : inflow.o

inicon_module.mod : inicon.o

inicon2_module.mod : inicon2.o

input_variable_module.mod : input_variable_module.o

macros_module.mod : macros.o

output_module.mod : output.o

ptmesh_module.mod : ptmesh.o


%.mod : %.o

%.o : %.f90 
	mpif90 $(FFLAGS) -c $<

clean:
	rm -rf BGK $(OBJ) *.mod
