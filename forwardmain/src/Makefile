TOP = ../..

include $(TOP)/buildoptions.mk

CFLAGS = $(GLOBAL_CFLAGS) $(MPI_INCDIR) $(METIS_INCDIR) $(LESLIB_DEFS) $(LESLIB_INCDIR) $(SOLVERIO_INCDIR) $(FLOWSOLVER_INCDIR)
CXXFLAGS  = $(GLOBAL_CXXFLAGS) $(MPI_INCDIR) $(METIS_INCDIR) $(LESLIB_DEFS) $(LESLIB_INCDIR) $(SOLVERIO_INCDIR) $(FLOWSOLVER_INCDIR) $(VTK_INCDIR) 
FFLAGS   = $(GLOBAL_FFLAGS) $(MPI_INCDIR) $(METIS_INCDIR)

HDRS  =   

CXXSRCS = main.cxx ../../estimation/src/PhGlobalArrayTransfer.cxx

CSRCS   = 

FSRCS   = 

SRCS  = $(FSRCS) $(CSRCS) $(CXXSRCS)

OBJS  = $(FSRCS:.f90=.$(OBJECTEXT)) $(CSRCS:.c=.$(OBJECTEXT)) $(CXXSRCS:.cxx=.$(OBJECTEXT))

$(TARGET_FLOWSOLVER): $(OBJS)
	$(LINK_EXE)$(TARGET_FLOWSOLVER) $(GLOBAL_LFLAGS) $(OBJS) $(FLOWSOLVER_LIBS) $(SOLVERIO_LIBS) $(METIS_LIBS) $(LESLIB_LIBS) $(MEMLS_LIBS) $(VTK_LIBS) $(FLOWSOLVER_LIBS) $(MPI_LIBS) $(F90_LIBS) $(CXX_LIBS)
					
	mv $(TARGET_FLOWSOLVER) $(TOP)/bin	

.PHONY : clean
clean : 
	-$(RM) $(OBJS) $(TARGET_FLOWSOLVER) *.mod *.exe *.pdb *.manifest
	-$(RM) $(TOP)/bin/$(TARGET_FLOWSOLVER)*
