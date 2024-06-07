CPPFLAGS += -I./include -I$(PACS_ROOT) -I$(PACS_ROOT)/include
LDLIBS += -L$(PACS_ROOT)/lib -lpacs -lmuparser -lMesh1D
