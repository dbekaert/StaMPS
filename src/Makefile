#################################################################
# "make [all]"		all exe's				#
# "make prog.o"		compile only prog.o			#
# "make prog"		compile only prog			#
# "make clean"		remove junk caused by makes		#
# "make cleaner"	remove all caused by makes 		#
# "make uninstall"	remove all caused by makes (incl. inst)	#
# "make install"	installs executables in INSTALL_DIR 	#
#                       					#
#%// Andy Hooper 08-May_2006					#
# ============================================================= #
# 10/09 AH: -lXp -lXext added to XLIBS                          #
# 11/09 AH: display removed from defaults                       #
# ============================================================= #
#################################################################

### X11 libraries
XLIBS = -L/usr/X11R6/lib -lXm -lXt -lX11 -lXp -lXext
XINCS = -I/usr/X11R6/include

### X11 libraries for Mac OS X Tiger or Leopard with IST-Inc 
### OpenMotif (/usr/OpenMotif)
### IST-INC: http://www.ist-inc.com/DOWNLOADS/motif_download.html
#XLIBS = -L/usr/X11R6/lib -lX11 -L/usr/OpenMotif/lib -lXm -lXt -lXp -lXext
#XINCS = -I/usr/X11R6/include -I/usr/OpenMotif/include
#
### X11 libraries for Mac OS X Tiger or Leopard with Macports.org
### OpenMotif (installed in /opt/local/)
#XLIBS = -L/usr/X11R6/lib -lX11 -L/opt/local/lib -lXm -lXt -lXp -lXext
#XINCS = -I/usr/X11R6/include -I/opt/local/include


### The shell we're using ###
SHELL	=	/bin/sh

### Install directory 
INSTALL_DIR =	../bin

### compiler
CC 	=	g++

### Compilation flags
CFLAGS	=	-O3 -std=c++11


##################################################################
### THERE SHOULD BE NOTHING YOU WANT TO CHANGE BELOW THIS LINE ###
##################################################################

### The programs we want
PROGS	=	calamp \
		selpsc_patch \
		selsbc_patch \
		cpxsum \
		pscphase \
		psclonlat \
		pscdem 

# usedin install:
SCRIPTS =	

### Compilation
default:	all
all:		$(PROGS) dummy

# dummy, since else makefile want to recompile the last util.
dummy:

# the utilities.
selsbc_patch_new: selsbc_patch_new.o
		  $(CC) $(CFLAGS) $@.o -o $@
selpsc_patch:	selpsc_patch.o
		$(CC) $(CFLAGS) $@.o -o $@
selsbc_patch:	selsbc_patch.o
		$(CC) $(CFLAGS) $@.o -o $@
calamp:		calamp.o
		$(CC) $(CFLAGS) $@.o -o $@
cpxsum:		cpxsum.o
		$(CC) $(CFLAGS) $@.o -o $@
pscphase:	pscphase.o
		$(CC) $(CFLAGS) $@.o -o $@
psclonlat:	psclonlat.o
		$(CC) $(CFLAGS) $@.o -o $@
pscdem:		pscdem.o
		$(CC) $(CFLAGS) $@.o -o $@
dismph:		dismph.o XInfo.o CDisp.o CGetData.o CDispComp.o bytescale.o
		$(CC) $(CFLAGS) $? -o $@ ${XLIBS}
.cpp.o: 
	 	${CC} -c ${XINCS} $*.cpp

### Install in INSTALL_DIR by linking/copying executables ###
install:	$(PROGS)
		$(MAKE) definstall; 


### Use symbolic links at our system, copy to /usr/local/bin on other.
definstall:	$(PROGS)
		dir=$(INSTALL_DIR); \
        	if test ! -d $$dir; then \
		  echo "Sorry, dir $(INSTALL_DIR) does not exist, exiting..."; exit; fi; \
		echo Installing in directory: $$dir; \
		list='$(PROGS) $(SCRIPTS)'; for p in $$list; do \
		  echo "Installing (copy): $$p"; \
		  cp -f $$p $(INSTALL_DIR)/$$p; \
		done
		$(MAKE) cleaner

### Helpers ###
clean:	
		@rm -f *.o *dummy* *.bak
		@echo "* Removed junk."
### some reason under cygwin progs are not removed without ".exe"???
CYGWIN_EXE = 	$(PROGS:=.exe)
cleaner:	clean
		@rm -f $(PROGS) $(CYGWIN_EXE) a.out
		@echo "* Removed executables in source dir:  `pwd`."
uninstall:	cleaner
		dir=$(INSTALL_DIR); \
        	if test -d $$dir; then \
		  cd $(INSTALL_DIR); rm -f $(PROGS) $(CYGWIN_EXE); \
		  echo "* Removed executables in install dir: $(INSTALL_DIR)."; fi;

### How to make object files ###
.cc.o:		
	$(CC) $(CFLAGS) -c -o $(@) $<

### EOF.

