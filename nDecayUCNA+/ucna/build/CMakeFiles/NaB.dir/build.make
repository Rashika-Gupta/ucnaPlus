# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.6

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /usr/bin/cmake3

# The command to remove a file.
RM = /usr/bin/cmake3 -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build

# Include any dependencies generated for this target.
include CMakeFiles/NaB.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/NaB.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/NaB.dir/flags.make

CMakeFiles/NaB.dir/newNab.cc.o: CMakeFiles/NaB.dir/flags.make
CMakeFiles/NaB.dir/newNab.cc.o: ../newNab.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/NaB.dir/newNab.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NaB.dir/newNab.cc.o -c /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/newNab.cc

CMakeFiles/NaB.dir/newNab.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NaB.dir/newNab.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/newNab.cc > CMakeFiles/NaB.dir/newNab.cc.i

CMakeFiles/NaB.dir/newNab.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NaB.dir/newNab.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/newNab.cc -o CMakeFiles/NaB.dir/newNab.cc.s

CMakeFiles/NaB.dir/newNab.cc.o.requires:

.PHONY : CMakeFiles/NaB.dir/newNab.cc.o.requires

CMakeFiles/NaB.dir/newNab.cc.o.provides: CMakeFiles/NaB.dir/newNab.cc.o.requires
	$(MAKE) -f CMakeFiles/NaB.dir/build.make CMakeFiles/NaB.dir/newNab.cc.o.provides.build
.PHONY : CMakeFiles/NaB.dir/newNab.cc.o.provides

CMakeFiles/NaB.dir/newNab.cc.o.provides.build: CMakeFiles/NaB.dir/newNab.cc.o


CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.o: CMakeFiles/NaB.dir/flags.make
CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.o: ../src/UCNBAnalysisManager.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.o -c /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBAnalysisManager.cc

CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBAnalysisManager.cc > CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.i

CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBAnalysisManager.cc -o CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.s

CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.o.requires:

.PHONY : CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.o.requires

CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.o.provides: CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.o.requires
	$(MAKE) -f CMakeFiles/NaB.dir/build.make CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.o.provides.build
.PHONY : CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.o.provides

CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.o.provides.build: CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.o


CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.o: CMakeFiles/NaB.dir/flags.make
CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.o: ../src/UCNBDetectorConstruction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.o -c /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBDetectorConstruction.cc

CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBDetectorConstruction.cc > CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.i

CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBDetectorConstruction.cc -o CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.s

CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.o.requires:

.PHONY : CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.o.requires

CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.o.provides: CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.o.requires
	$(MAKE) -f CMakeFiles/NaB.dir/build.make CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.o.provides.build
.PHONY : CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.o.provides

CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.o.provides.build: CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.o


CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.o: CMakeFiles/NaB.dir/flags.make
CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.o: ../src/UCNBDetectorMessenger.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.o -c /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBDetectorMessenger.cc

CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBDetectorMessenger.cc > CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.i

CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBDetectorMessenger.cc -o CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.s

CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.o.requires:

.PHONY : CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.o.requires

CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.o.provides: CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.o.requires
	$(MAKE) -f CMakeFiles/NaB.dir/build.make CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.o.provides.build
.PHONY : CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.o.provides

CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.o.provides.build: CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.o


CMakeFiles/NaB.dir/src/UCNBEventAction.cc.o: CMakeFiles/NaB.dir/flags.make
CMakeFiles/NaB.dir/src/UCNBEventAction.cc.o: ../src/UCNBEventAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/NaB.dir/src/UCNBEventAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NaB.dir/src/UCNBEventAction.cc.o -c /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBEventAction.cc

CMakeFiles/NaB.dir/src/UCNBEventAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NaB.dir/src/UCNBEventAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBEventAction.cc > CMakeFiles/NaB.dir/src/UCNBEventAction.cc.i

CMakeFiles/NaB.dir/src/UCNBEventAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NaB.dir/src/UCNBEventAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBEventAction.cc -o CMakeFiles/NaB.dir/src/UCNBEventAction.cc.s

CMakeFiles/NaB.dir/src/UCNBEventAction.cc.o.requires:

.PHONY : CMakeFiles/NaB.dir/src/UCNBEventAction.cc.o.requires

CMakeFiles/NaB.dir/src/UCNBEventAction.cc.o.provides: CMakeFiles/NaB.dir/src/UCNBEventAction.cc.o.requires
	$(MAKE) -f CMakeFiles/NaB.dir/build.make CMakeFiles/NaB.dir/src/UCNBEventAction.cc.o.provides.build
.PHONY : CMakeFiles/NaB.dir/src/UCNBEventAction.cc.o.provides

CMakeFiles/NaB.dir/src/UCNBEventAction.cc.o.provides.build: CMakeFiles/NaB.dir/src/UCNBEventAction.cc.o


CMakeFiles/NaB.dir/src/UCNBField.cc.o: CMakeFiles/NaB.dir/flags.make
CMakeFiles/NaB.dir/src/UCNBField.cc.o: ../src/UCNBField.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/NaB.dir/src/UCNBField.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NaB.dir/src/UCNBField.cc.o -c /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBField.cc

CMakeFiles/NaB.dir/src/UCNBField.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NaB.dir/src/UCNBField.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBField.cc > CMakeFiles/NaB.dir/src/UCNBField.cc.i

CMakeFiles/NaB.dir/src/UCNBField.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NaB.dir/src/UCNBField.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBField.cc -o CMakeFiles/NaB.dir/src/UCNBField.cc.s

CMakeFiles/NaB.dir/src/UCNBField.cc.o.requires:

.PHONY : CMakeFiles/NaB.dir/src/UCNBField.cc.o.requires

CMakeFiles/NaB.dir/src/UCNBField.cc.o.provides: CMakeFiles/NaB.dir/src/UCNBField.cc.o.requires
	$(MAKE) -f CMakeFiles/NaB.dir/build.make CMakeFiles/NaB.dir/src/UCNBField.cc.o.provides.build
.PHONY : CMakeFiles/NaB.dir/src/UCNBField.cc.o.provides

CMakeFiles/NaB.dir/src/UCNBField.cc.o.provides.build: CMakeFiles/NaB.dir/src/UCNBField.cc.o


CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.o: CMakeFiles/NaB.dir/flags.make
CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.o: ../src/UCNBPhysicsList.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.o -c /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBPhysicsList.cc

CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBPhysicsList.cc > CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.i

CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBPhysicsList.cc -o CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.s

CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.o.requires:

.PHONY : CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.o.requires

CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.o.provides: CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.o.requires
	$(MAKE) -f CMakeFiles/NaB.dir/build.make CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.o.provides.build
.PHONY : CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.o.provides

CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.o.provides.build: CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.o


CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.o: CMakeFiles/NaB.dir/flags.make
CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.o: ../src/UCNBPrimaryGeneratorAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.o -c /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBPrimaryGeneratorAction.cc

CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBPrimaryGeneratorAction.cc > CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.i

CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBPrimaryGeneratorAction.cc -o CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.s

CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.o.requires:

.PHONY : CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.o.requires

CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.o.provides: CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.o.requires
	$(MAKE) -f CMakeFiles/NaB.dir/build.make CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.o.provides.build
.PHONY : CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.o.provides

CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.o.provides.build: CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.o


CMakeFiles/NaB.dir/src/UCNBRunAction.cc.o: CMakeFiles/NaB.dir/flags.make
CMakeFiles/NaB.dir/src/UCNBRunAction.cc.o: ../src/UCNBRunAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/NaB.dir/src/UCNBRunAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NaB.dir/src/UCNBRunAction.cc.o -c /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBRunAction.cc

CMakeFiles/NaB.dir/src/UCNBRunAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NaB.dir/src/UCNBRunAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBRunAction.cc > CMakeFiles/NaB.dir/src/UCNBRunAction.cc.i

CMakeFiles/NaB.dir/src/UCNBRunAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NaB.dir/src/UCNBRunAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBRunAction.cc -o CMakeFiles/NaB.dir/src/UCNBRunAction.cc.s

CMakeFiles/NaB.dir/src/UCNBRunAction.cc.o.requires:

.PHONY : CMakeFiles/NaB.dir/src/UCNBRunAction.cc.o.requires

CMakeFiles/NaB.dir/src/UCNBRunAction.cc.o.provides: CMakeFiles/NaB.dir/src/UCNBRunAction.cc.o.requires
	$(MAKE) -f CMakeFiles/NaB.dir/build.make CMakeFiles/NaB.dir/src/UCNBRunAction.cc.o.provides.build
.PHONY : CMakeFiles/NaB.dir/src/UCNBRunAction.cc.o.provides

CMakeFiles/NaB.dir/src/UCNBRunAction.cc.o.provides.build: CMakeFiles/NaB.dir/src/UCNBRunAction.cc.o


CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.o: CMakeFiles/NaB.dir/flags.make
CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.o: ../src/UCNBStackingAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.o -c /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBStackingAction.cc

CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBStackingAction.cc > CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.i

CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBStackingAction.cc -o CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.s

CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.o.requires:

.PHONY : CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.o.requires

CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.o.provides: CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.o.requires
	$(MAKE) -f CMakeFiles/NaB.dir/build.make CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.o.provides.build
.PHONY : CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.o.provides

CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.o.provides.build: CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.o


CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.o: CMakeFiles/NaB.dir/flags.make
CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.o: ../src/UCNBSteppingAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.o -c /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBSteppingAction.cc

CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBSteppingAction.cc > CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.i

CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBSteppingAction.cc -o CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.s

CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.o.requires:

.PHONY : CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.o.requires

CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.o.provides: CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.o.requires
	$(MAKE) -f CMakeFiles/NaB.dir/build.make CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.o.provides.build
.PHONY : CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.o.provides

CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.o.provides.build: CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.o


CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.o: CMakeFiles/NaB.dir/flags.make
CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.o: ../src/UCNBSteppingVerbose.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.o -c /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBSteppingVerbose.cc

CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBSteppingVerbose.cc > CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.i

CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/src/UCNBSteppingVerbose.cc -o CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.s

CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.o.requires:

.PHONY : CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.o.requires

CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.o.provides: CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.o.requires
	$(MAKE) -f CMakeFiles/NaB.dir/build.make CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.o.provides.build
.PHONY : CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.o.provides

CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.o.provides.build: CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.o


# Object files for target NaB
NaB_OBJECTS = \
"CMakeFiles/NaB.dir/newNab.cc.o" \
"CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.o" \
"CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.o" \
"CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.o" \
"CMakeFiles/NaB.dir/src/UCNBEventAction.cc.o" \
"CMakeFiles/NaB.dir/src/UCNBField.cc.o" \
"CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.o" \
"CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.o" \
"CMakeFiles/NaB.dir/src/UCNBRunAction.cc.o" \
"CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.o" \
"CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.o" \
"CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.o"

# External object files for target NaB
NaB_EXTERNAL_OBJECTS =

NaB: CMakeFiles/NaB.dir/newNab.cc.o
NaB: CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.o
NaB: CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.o
NaB: CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.o
NaB: CMakeFiles/NaB.dir/src/UCNBEventAction.cc.o
NaB: CMakeFiles/NaB.dir/src/UCNBField.cc.o
NaB: CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.o
NaB: CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.o
NaB: CMakeFiles/NaB.dir/src/UCNBRunAction.cc.o
NaB: CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.o
NaB: CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.o
NaB: CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.o
NaB: CMakeFiles/NaB.dir/build.make
NaB: /scif/apps/geant4/lib64/libG4Tree.so
NaB: /scif/apps/geant4/lib64/libG4GMocren.so
NaB: /scif/apps/geant4/lib64/libG4visHepRep.so
NaB: /scif/apps/geant4/lib64/libG4RayTracer.so
NaB: /scif/apps/geant4/lib64/libG4VRML.so
NaB: /scif/apps/geant4/lib64/libG4OpenGL.so
NaB: /scif/apps/geant4/lib64/libG4gl2ps.so
NaB: /scif/apps/geant4/lib64/libG4interfaces.so
NaB: /scif/apps/geant4/lib64/libG4persistency.so
NaB: /scif/apps/geant4/lib64/libG4analysis.so
NaB: /scif/apps/geant4/lib64/libG4error_propagation.so
NaB: /scif/apps/geant4/lib64/libG4readout.so
NaB: /scif/apps/geant4/lib64/libG4physicslists.so
NaB: /scif/apps/geant4/lib64/libG4parmodels.so
NaB: /scif/apps/geant4/lib64/libG4FR.so
NaB: /scif/apps/geant4/lib64/libG4vis_management.so
NaB: /scif/apps/geant4/lib64/libG4modeling.so
NaB: /usr/lib64/libXm.so
NaB: /usr/lib64/libSM.so
NaB: /usr/lib64/libICE.so
NaB: /usr/lib64/libX11.so
NaB: /usr/lib64/libXext.so
NaB: /usr/lib64/libXt.so
NaB: /usr/lib64/libXmu.so
NaB: /usr/lib64/libGLU.so
NaB: /usr/lib64/libGL.so
NaB: /usr/lib64/libQtOpenGL.so
NaB: /usr/lib64/libQtGui.so
NaB: /usr/lib64/libQtGui_debug.so
NaB: /usr/lib64/libQtCore.so
NaB: /usr/lib64/libQtCore_debug.so
NaB: /usr/lib64/libxerces-c.so
NaB: /scif/apps/geant4/lib64/libG4run.so
NaB: /scif/apps/geant4/lib64/libG4event.so
NaB: /scif/apps/geant4/lib64/libG4tracking.so
NaB: /scif/apps/geant4/lib64/libG4processes.so
NaB: /scif/apps/geant4/lib64/libG4zlib.so
NaB: /usr/lib64/libexpat.so
NaB: /scif/apps/geant4/lib64/libG4digits_hits.so
NaB: /scif/apps/geant4/lib64/libG4track.so
NaB: /scif/apps/geant4/lib64/libG4particles.so
NaB: /scif/apps/geant4/lib64/libG4geometry.so
NaB: /scif/apps/geant4/lib64/libG4materials.so
NaB: /scif/apps/geant4/lib64/libG4graphics_reps.so
NaB: /scif/apps/geant4/lib64/libG4intercoms.so
NaB: /scif/apps/geant4/lib64/libG4global.so
NaB: /scif/apps/clhep/lib/libCLHEP.so
NaB: CMakeFiles/NaB.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Linking CXX executable NaB"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/NaB.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/NaB.dir/build: NaB

.PHONY : CMakeFiles/NaB.dir/build

CMakeFiles/NaB.dir/requires: CMakeFiles/NaB.dir/newNab.cc.o.requires
CMakeFiles/NaB.dir/requires: CMakeFiles/NaB.dir/src/UCNBAnalysisManager.cc.o.requires
CMakeFiles/NaB.dir/requires: CMakeFiles/NaB.dir/src/UCNBDetectorConstruction.cc.o.requires
CMakeFiles/NaB.dir/requires: CMakeFiles/NaB.dir/src/UCNBDetectorMessenger.cc.o.requires
CMakeFiles/NaB.dir/requires: CMakeFiles/NaB.dir/src/UCNBEventAction.cc.o.requires
CMakeFiles/NaB.dir/requires: CMakeFiles/NaB.dir/src/UCNBField.cc.o.requires
CMakeFiles/NaB.dir/requires: CMakeFiles/NaB.dir/src/UCNBPhysicsList.cc.o.requires
CMakeFiles/NaB.dir/requires: CMakeFiles/NaB.dir/src/UCNBPrimaryGeneratorAction.cc.o.requires
CMakeFiles/NaB.dir/requires: CMakeFiles/NaB.dir/src/UCNBRunAction.cc.o.requires
CMakeFiles/NaB.dir/requires: CMakeFiles/NaB.dir/src/UCNBStackingAction.cc.o.requires
CMakeFiles/NaB.dir/requires: CMakeFiles/NaB.dir/src/UCNBSteppingAction.cc.o.requires
CMakeFiles/NaB.dir/requires: CMakeFiles/NaB.dir/src/UCNBSteppingVerbose.cc.o.requires

.PHONY : CMakeFiles/NaB.dir/requires

CMakeFiles/NaB.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/NaB.dir/cmake_clean.cmake
.PHONY : CMakeFiles/NaB.dir/clean

CMakeFiles/NaB.dir/depend:
	cd /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build /mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/build/CMakeFiles/NaB.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/NaB.dir/depend

