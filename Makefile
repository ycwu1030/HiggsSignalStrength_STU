PACKAGE = HiggsSignalStrengthSTU

SrcSuf := cpp
IncSuf := h
ObjSuf := o
DepSuf := d
DllSuf := so

MAINDIR := $(shell pwd)
SrcDir := $(MAINDIR)
IncDir := $(MAINDIR)
ObjDir := $(MAINDIR)
DepDir := $(ObjDir)

INCDIRS = $(shell pwd)
EXTRAHDR = 

INCLUDES      = $(addprefix -I, $(INCDIRS) ) 
CXX           = clang++
LD            = clang++
SOFLAGS       = -shared
SOFLAGS       = -dynamiclib -single_module -install_name $(MAINDIR)/

ifdef DEBUG
  CXXFLAGS    = -g -O0 
  LDFLAGS     = -g -O0
  DEFINES     =
else
  CXXFLAGS    = -O3
  LDFLAGS     = -O3
  DEFINES     = 
#  DEFINES     = -DNDEBUG
endif

DEFINES      += -DLINUXVERS -DHAS_SSTREAM
ifdef VERBOSE
DEFINES      += -DVERBOSE
endif
ifdef TESTCODE
DEFINES      += -DTESTCODE
endif

CXXFLAGS     += $(DEFINES) $(INCLUDES)
LIBS = 
LDFLAGS = 

SRC = $(wildcard $(SrcDir)/*.$(SrcSuf))
HEAD = $(wildcard $(IncDir)/*.$(IncSuf))
OBJ = $(patsubst $(SrcDir)/%.$(SrcSuf),$(ObjDir)/%.$(ObjSuf),$(SRC))
DEP = $(patsubst $(SrcDir)/%.$(SrcSuf),$(DepDir)/%.$(DepSuf),$(SRC))
USERLIB = lib$(PACKAGE).$(DllSuf)

all:$(USERLIB) 

$(USERLIB):$(OBJ) $(SRC) $(HEAD)
	@echo "Linking Package $@..."
	@$(LD) -fPIC $(SOFLAGS)$@ $(LDFLAGS) -o $@ $(OBJ) $(LIBS)
	@echo "--->$@ Done!"

$(ObjDir)/%.$(ObjSuf):$(SrcDir)/%.$(SrcSuf)
	@echo "Compile $<..."
	@$(CXX) -fPIC $(CXXFLAGS) -o $@ -c $<

#-------------------------------------------------------------------------------
clean:
	@rm -rf $(OBJ) $(USERLIB)
	@echo "--->Clean Done"
