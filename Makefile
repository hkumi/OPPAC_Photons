# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.3

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /mnt/misc/sw/x86_64/Debian/8/cmake/3.3.1/bin/cmake

# The command to remove a file.
RM = /mnt/misc/sw/x86_64/Debian/8/cmake/3.3.1/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /user/cortesi/g4/oppac_0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /user/cortesi/g4/oppac_0

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/mnt/misc/sw/x86_64/Debian/8/cmake/3.3.1/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target list_install_components
list_install_components:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Available install components are: \"Unspecified\""
.PHONY : list_install_components

# Special rule for the target list_install_components
list_install_components/fast: list_install_components

.PHONY : list_install_components/fast

# Special rule for the target install
install: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/mnt/misc/sw/x86_64/Debian/8/cmake/3.3.1/bin/cmake -P cmake_install.cmake
.PHONY : install

# Special rule for the target install
install/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/mnt/misc/sw/x86_64/Debian/8/cmake/3.3.1/bin/cmake -P cmake_install.cmake
.PHONY : install/fast

# Special rule for the target install/local
install/local: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/mnt/misc/sw/x86_64/Debian/8/cmake/3.3.1/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local

# Special rule for the target install/local
install/local/fast: install/local

.PHONY : install/local/fast

# Special rule for the target install/strip
install/strip: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/mnt/misc/sw/x86_64/Debian/8/cmake/3.3.1/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip

# Special rule for the target install/strip
install/strip/fast: install/strip

.PHONY : install/strip/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/mnt/misc/sw/x86_64/Debian/8/cmake/3.3.1/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /user/cortesi/g4/oppac_0/CMakeFiles /user/cortesi/g4/oppac_0/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /user/cortesi/g4/oppac_0/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named exe_scinti

# Build rule for target.
exe_scinti: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 exe_scinti
.PHONY : exe_scinti

# fast build rule for target.
exe_scinti/fast:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/build
.PHONY : exe_scinti/fast

#=============================================================================
# Target rules for targets named scinti

# Build rule for target.
scinti: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 scinti
.PHONY : scinti

# fast build rule for target.
scinti/fast:
	$(MAKE) -f CMakeFiles/scinti.dir/build.make CMakeFiles/scinti.dir/build
.PHONY : scinti/fast

scinti.o: scinti.cc.o

.PHONY : scinti.o

# target to build an object file
scinti.cc.o:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/scinti.cc.o
.PHONY : scinti.cc.o

scinti.i: scinti.cc.i

.PHONY : scinti.i

# target to preprocess a source file
scinti.cc.i:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/scinti.cc.i
.PHONY : scinti.cc.i

scinti.s: scinti.cc.s

.PHONY : scinti.s

# target to generate assembly for a file
scinti.cc.s:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/scinti.cc.s
.PHONY : scinti.cc.s

src/DC.o: src/DC.cc.o

.PHONY : src/DC.o

# target to build an object file
src/DC.cc.o:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/DC.cc.o
.PHONY : src/DC.cc.o

src/DC.i: src/DC.cc.i

.PHONY : src/DC.i

# target to preprocess a source file
src/DC.cc.i:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/DC.cc.i
.PHONY : src/DC.cc.i

src/DC.s: src/DC.cc.s

.PHONY : src/DC.s

# target to generate assembly for a file
src/DC.cc.s:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/DC.cc.s
.PHONY : src/DC.cc.s

src/PG.o: src/PG.cc.o

.PHONY : src/PG.o

# target to build an object file
src/PG.cc.o:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/PG.cc.o
.PHONY : src/PG.cc.o

src/PG.i: src/PG.cc.i

.PHONY : src/PG.i

# target to preprocess a source file
src/PG.cc.i:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/PG.cc.i
.PHONY : src/PG.cc.i

src/PG.s: src/PG.cc.s

.PHONY : src/PG.s

# target to generate assembly for a file
src/PG.cc.s:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/PG.cc.s
.PHONY : src/PG.cc.s

src/PL.o: src/PL.cc.o

.PHONY : src/PL.o

# target to build an object file
src/PL.cc.o:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/PL.cc.o
.PHONY : src/PL.cc.o

src/PL.i: src/PL.cc.i

.PHONY : src/PL.i

# target to preprocess a source file
src/PL.cc.i:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/PL.cc.i
.PHONY : src/PL.cc.i

src/PL.s: src/PL.cc.s

.PHONY : src/PL.s

# target to generate assembly for a file
src/PL.cc.s:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/PL.cc.s
.PHONY : src/PL.cc.s

src/Run.o: src/Run.cc.o

.PHONY : src/Run.o

# target to build an object file
src/Run.cc.o:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/Run.cc.o
.PHONY : src/Run.cc.o

src/Run.i: src/Run.cc.i

.PHONY : src/Run.i

# target to preprocess a source file
src/Run.cc.i:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/Run.cc.i
.PHONY : src/Run.cc.i

src/Run.s: src/Run.cc.s

.PHONY : src/Run.s

# target to generate assembly for a file
src/Run.cc.s:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/Run.cc.s
.PHONY : src/Run.cc.s

src/RunAction.o: src/RunAction.cc.o

.PHONY : src/RunAction.o

# target to build an object file
src/RunAction.cc.o:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/RunAction.cc.o
.PHONY : src/RunAction.cc.o

src/RunAction.i: src/RunAction.cc.i

.PHONY : src/RunAction.i

# target to preprocess a source file
src/RunAction.cc.i:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/RunAction.cc.i
.PHONY : src/RunAction.cc.i

src/RunAction.s: src/RunAction.cc.s

.PHONY : src/RunAction.s

# target to generate assembly for a file
src/RunAction.cc.s:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/RunAction.cc.s
.PHONY : src/RunAction.cc.s

src/StepAction.o: src/StepAction.cc.o

.PHONY : src/StepAction.o

# target to build an object file
src/StepAction.cc.o:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/StepAction.cc.o
.PHONY : src/StepAction.cc.o

src/StepAction.i: src/StepAction.cc.i

.PHONY : src/StepAction.i

# target to preprocess a source file
src/StepAction.cc.i:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/StepAction.cc.i
.PHONY : src/StepAction.cc.i

src/StepAction.s: src/StepAction.cc.s

.PHONY : src/StepAction.s

# target to generate assembly for a file
src/StepAction.cc.s:
	$(MAKE) -f CMakeFiles/exe_scinti.dir/build.make CMakeFiles/exe_scinti.dir/src/StepAction.cc.s
.PHONY : src/StepAction.cc.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... list_install_components"
	@echo "... install"
	@echo "... install/local"
	@echo "... install/strip"
	@echo "... edit_cache"
	@echo "... scinti"
	@echo "... exe_scinti"
	@echo "... scinti.o"
	@echo "... scinti.i"
	@echo "... scinti.s"
	@echo "... src/DC.o"
	@echo "... src/DC.i"
	@echo "... src/DC.s"
	@echo "... src/PG.o"
	@echo "... src/PG.i"
	@echo "... src/PG.s"
	@echo "... src/PL.o"
	@echo "... src/PL.i"
	@echo "... src/PL.s"
	@echo "... src/Run.o"
	@echo "... src/Run.i"
	@echo "... src/Run.s"
	@echo "... src/RunAction.o"
	@echo "... src/RunAction.i"
	@echo "... src/RunAction.s"
	@echo "... src/StepAction.o"
	@echo "... src/StepAction.i"
	@echo "... src/StepAction.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

