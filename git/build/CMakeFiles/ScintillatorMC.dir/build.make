# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build

# Include any dependencies generated for this target.
include CMakeFiles/ScintillatorMC.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ScintillatorMC.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ScintillatorMC.dir/flags.make

CMakeFiles/ScintillatorMC.dir/ScintillatorMC.cc.o: CMakeFiles/ScintillatorMC.dir/flags.make
CMakeFiles/ScintillatorMC.dir/ScintillatorMC.cc.o: ../ScintillatorMC.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ScintillatorMC.dir/ScintillatorMC.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ScintillatorMC.dir/ScintillatorMC.cc.o -c /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/ScintillatorMC.cc

CMakeFiles/ScintillatorMC.dir/ScintillatorMC.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ScintillatorMC.dir/ScintillatorMC.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/ScintillatorMC.cc > CMakeFiles/ScintillatorMC.dir/ScintillatorMC.cc.i

CMakeFiles/ScintillatorMC.dir/ScintillatorMC.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ScintillatorMC.dir/ScintillatorMC.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/ScintillatorMC.cc -o CMakeFiles/ScintillatorMC.dir/ScintillatorMC.cc.s

CMakeFiles/ScintillatorMC.dir/src/Analysis.cc.o: CMakeFiles/ScintillatorMC.dir/flags.make
CMakeFiles/ScintillatorMC.dir/src/Analysis.cc.o: ../src/Analysis.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/ScintillatorMC.dir/src/Analysis.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ScintillatorMC.dir/src/Analysis.cc.o -c /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/Analysis.cc

CMakeFiles/ScintillatorMC.dir/src/Analysis.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ScintillatorMC.dir/src/Analysis.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/Analysis.cc > CMakeFiles/ScintillatorMC.dir/src/Analysis.cc.i

CMakeFiles/ScintillatorMC.dir/src/Analysis.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ScintillatorMC.dir/src/Analysis.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/Analysis.cc -o CMakeFiles/ScintillatorMC.dir/src/Analysis.cc.s

CMakeFiles/ScintillatorMC.dir/src/DetectorConstruction.cc.o: CMakeFiles/ScintillatorMC.dir/flags.make
CMakeFiles/ScintillatorMC.dir/src/DetectorConstruction.cc.o: ../src/DetectorConstruction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/ScintillatorMC.dir/src/DetectorConstruction.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ScintillatorMC.dir/src/DetectorConstruction.cc.o -c /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/DetectorConstruction.cc

CMakeFiles/ScintillatorMC.dir/src/DetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ScintillatorMC.dir/src/DetectorConstruction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/DetectorConstruction.cc > CMakeFiles/ScintillatorMC.dir/src/DetectorConstruction.cc.i

CMakeFiles/ScintillatorMC.dir/src/DetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ScintillatorMC.dir/src/DetectorConstruction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/DetectorConstruction.cc -o CMakeFiles/ScintillatorMC.dir/src/DetectorConstruction.cc.s

CMakeFiles/ScintillatorMC.dir/src/HadrontherapyStepMax.cc.o: CMakeFiles/ScintillatorMC.dir/flags.make
CMakeFiles/ScintillatorMC.dir/src/HadrontherapyStepMax.cc.o: ../src/HadrontherapyStepMax.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/ScintillatorMC.dir/src/HadrontherapyStepMax.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ScintillatorMC.dir/src/HadrontherapyStepMax.cc.o -c /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/HadrontherapyStepMax.cc

CMakeFiles/ScintillatorMC.dir/src/HadrontherapyStepMax.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ScintillatorMC.dir/src/HadrontherapyStepMax.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/HadrontherapyStepMax.cc > CMakeFiles/ScintillatorMC.dir/src/HadrontherapyStepMax.cc.i

CMakeFiles/ScintillatorMC.dir/src/HadrontherapyStepMax.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ScintillatorMC.dir/src/HadrontherapyStepMax.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/HadrontherapyStepMax.cc -o CMakeFiles/ScintillatorMC.dir/src/HadrontherapyStepMax.cc.s

CMakeFiles/ScintillatorMC.dir/src/LocalIonIonInelasticPhysic.cc.o: CMakeFiles/ScintillatorMC.dir/flags.make
CMakeFiles/ScintillatorMC.dir/src/LocalIonIonInelasticPhysic.cc.o: ../src/LocalIonIonInelasticPhysic.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/ScintillatorMC.dir/src/LocalIonIonInelasticPhysic.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ScintillatorMC.dir/src/LocalIonIonInelasticPhysic.cc.o -c /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/LocalIonIonInelasticPhysic.cc

CMakeFiles/ScintillatorMC.dir/src/LocalIonIonInelasticPhysic.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ScintillatorMC.dir/src/LocalIonIonInelasticPhysic.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/LocalIonIonInelasticPhysic.cc > CMakeFiles/ScintillatorMC.dir/src/LocalIonIonInelasticPhysic.cc.i

CMakeFiles/ScintillatorMC.dir/src/LocalIonIonInelasticPhysic.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ScintillatorMC.dir/src/LocalIonIonInelasticPhysic.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/LocalIonIonInelasticPhysic.cc -o CMakeFiles/ScintillatorMC.dir/src/LocalIonIonInelasticPhysic.cc.s

CMakeFiles/ScintillatorMC.dir/src/OrganicMaterial.cc.o: CMakeFiles/ScintillatorMC.dir/flags.make
CMakeFiles/ScintillatorMC.dir/src/OrganicMaterial.cc.o: ../src/OrganicMaterial.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/ScintillatorMC.dir/src/OrganicMaterial.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ScintillatorMC.dir/src/OrganicMaterial.cc.o -c /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/OrganicMaterial.cc

CMakeFiles/ScintillatorMC.dir/src/OrganicMaterial.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ScintillatorMC.dir/src/OrganicMaterial.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/OrganicMaterial.cc > CMakeFiles/ScintillatorMC.dir/src/OrganicMaterial.cc.i

CMakeFiles/ScintillatorMC.dir/src/OrganicMaterial.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ScintillatorMC.dir/src/OrganicMaterial.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/OrganicMaterial.cc -o CMakeFiles/ScintillatorMC.dir/src/OrganicMaterial.cc.s

CMakeFiles/ScintillatorMC.dir/src/PhysicsList.cc.o: CMakeFiles/ScintillatorMC.dir/flags.make
CMakeFiles/ScintillatorMC.dir/src/PhysicsList.cc.o: ../src/PhysicsList.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/ScintillatorMC.dir/src/PhysicsList.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ScintillatorMC.dir/src/PhysicsList.cc.o -c /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/PhysicsList.cc

CMakeFiles/ScintillatorMC.dir/src/PhysicsList.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ScintillatorMC.dir/src/PhysicsList.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/PhysicsList.cc > CMakeFiles/ScintillatorMC.dir/src/PhysicsList.cc.i

CMakeFiles/ScintillatorMC.dir/src/PhysicsList.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ScintillatorMC.dir/src/PhysicsList.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/PhysicsList.cc -o CMakeFiles/ScintillatorMC.dir/src/PhysicsList.cc.s

CMakeFiles/ScintillatorMC.dir/src/PrimaryGeneratorAction.cc.o: CMakeFiles/ScintillatorMC.dir/flags.make
CMakeFiles/ScintillatorMC.dir/src/PrimaryGeneratorAction.cc.o: ../src/PrimaryGeneratorAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/ScintillatorMC.dir/src/PrimaryGeneratorAction.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ScintillatorMC.dir/src/PrimaryGeneratorAction.cc.o -c /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/PrimaryGeneratorAction.cc

CMakeFiles/ScintillatorMC.dir/src/PrimaryGeneratorAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ScintillatorMC.dir/src/PrimaryGeneratorAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/PrimaryGeneratorAction.cc > CMakeFiles/ScintillatorMC.dir/src/PrimaryGeneratorAction.cc.i

CMakeFiles/ScintillatorMC.dir/src/PrimaryGeneratorAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ScintillatorMC.dir/src/PrimaryGeneratorAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/PrimaryGeneratorAction.cc -o CMakeFiles/ScintillatorMC.dir/src/PrimaryGeneratorAction.cc.s

CMakeFiles/ScintillatorMC.dir/src/SensitiveDetector.cc.o: CMakeFiles/ScintillatorMC.dir/flags.make
CMakeFiles/ScintillatorMC.dir/src/SensitiveDetector.cc.o: ../src/SensitiveDetector.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/ScintillatorMC.dir/src/SensitiveDetector.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ScintillatorMC.dir/src/SensitiveDetector.cc.o -c /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/SensitiveDetector.cc

CMakeFiles/ScintillatorMC.dir/src/SensitiveDetector.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ScintillatorMC.dir/src/SensitiveDetector.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/SensitiveDetector.cc > CMakeFiles/ScintillatorMC.dir/src/SensitiveDetector.cc.i

CMakeFiles/ScintillatorMC.dir/src/SensitiveDetector.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ScintillatorMC.dir/src/SensitiveDetector.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/SensitiveDetector.cc -o CMakeFiles/ScintillatorMC.dir/src/SensitiveDetector.cc.s

CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorScintillator.cc.o: CMakeFiles/ScintillatorMC.dir/flags.make
CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorScintillator.cc.o: ../src/SensitiveDetectorScintillator.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorScintillator.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorScintillator.cc.o -c /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/SensitiveDetectorScintillator.cc

CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorScintillator.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorScintillator.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/SensitiveDetectorScintillator.cc > CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorScintillator.cc.i

CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorScintillator.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorScintillator.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/SensitiveDetectorScintillator.cc -o CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorScintillator.cc.s

CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorTracker.cc.o: CMakeFiles/ScintillatorMC.dir/flags.make
CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorTracker.cc.o: ../src/SensitiveDetectorTracker.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorTracker.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorTracker.cc.o -c /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/SensitiveDetectorTracker.cc

CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorTracker.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorTracker.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/SensitiveDetectorTracker.cc > CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorTracker.cc.i

CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorTracker.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorTracker.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/SensitiveDetectorTracker.cc -o CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorTracker.cc.s

CMakeFiles/ScintillatorMC.dir/src/SteppingAction.cc.o: CMakeFiles/ScintillatorMC.dir/flags.make
CMakeFiles/ScintillatorMC.dir/src/SteppingAction.cc.o: ../src/SteppingAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/ScintillatorMC.dir/src/SteppingAction.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ScintillatorMC.dir/src/SteppingAction.cc.o -c /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/SteppingAction.cc

CMakeFiles/ScintillatorMC.dir/src/SteppingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ScintillatorMC.dir/src/SteppingAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/SteppingAction.cc > CMakeFiles/ScintillatorMC.dir/src/SteppingAction.cc.i

CMakeFiles/ScintillatorMC.dir/src/SteppingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ScintillatorMC.dir/src/SteppingAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/SteppingAction.cc -o CMakeFiles/ScintillatorMC.dir/src/SteppingAction.cc.s

CMakeFiles/ScintillatorMC.dir/src/pCTconfig.cc.o: CMakeFiles/ScintillatorMC.dir/flags.make
CMakeFiles/ScintillatorMC.dir/src/pCTconfig.cc.o: ../src/pCTconfig.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/ScintillatorMC.dir/src/pCTconfig.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ScintillatorMC.dir/src/pCTconfig.cc.o -c /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/pCTconfig.cc

CMakeFiles/ScintillatorMC.dir/src/pCTconfig.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ScintillatorMC.dir/src/pCTconfig.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/pCTconfig.cc > CMakeFiles/ScintillatorMC.dir/src/pCTconfig.cc.i

CMakeFiles/ScintillatorMC.dir/src/pCTconfig.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ScintillatorMC.dir/src/pCTconfig.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/src/pCTconfig.cc -o CMakeFiles/ScintillatorMC.dir/src/pCTconfig.cc.s

# Object files for target ScintillatorMC
ScintillatorMC_OBJECTS = \
"CMakeFiles/ScintillatorMC.dir/ScintillatorMC.cc.o" \
"CMakeFiles/ScintillatorMC.dir/src/Analysis.cc.o" \
"CMakeFiles/ScintillatorMC.dir/src/DetectorConstruction.cc.o" \
"CMakeFiles/ScintillatorMC.dir/src/HadrontherapyStepMax.cc.o" \
"CMakeFiles/ScintillatorMC.dir/src/LocalIonIonInelasticPhysic.cc.o" \
"CMakeFiles/ScintillatorMC.dir/src/OrganicMaterial.cc.o" \
"CMakeFiles/ScintillatorMC.dir/src/PhysicsList.cc.o" \
"CMakeFiles/ScintillatorMC.dir/src/PrimaryGeneratorAction.cc.o" \
"CMakeFiles/ScintillatorMC.dir/src/SensitiveDetector.cc.o" \
"CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorScintillator.cc.o" \
"CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorTracker.cc.o" \
"CMakeFiles/ScintillatorMC.dir/src/SteppingAction.cc.o" \
"CMakeFiles/ScintillatorMC.dir/src/pCTconfig.cc.o"

# External object files for target ScintillatorMC
ScintillatorMC_EXTERNAL_OBJECTS =

ScintillatorMC: CMakeFiles/ScintillatorMC.dir/ScintillatorMC.cc.o
ScintillatorMC: CMakeFiles/ScintillatorMC.dir/src/Analysis.cc.o
ScintillatorMC: CMakeFiles/ScintillatorMC.dir/src/DetectorConstruction.cc.o
ScintillatorMC: CMakeFiles/ScintillatorMC.dir/src/HadrontherapyStepMax.cc.o
ScintillatorMC: CMakeFiles/ScintillatorMC.dir/src/LocalIonIonInelasticPhysic.cc.o
ScintillatorMC: CMakeFiles/ScintillatorMC.dir/src/OrganicMaterial.cc.o
ScintillatorMC: CMakeFiles/ScintillatorMC.dir/src/PhysicsList.cc.o
ScintillatorMC: CMakeFiles/ScintillatorMC.dir/src/PrimaryGeneratorAction.cc.o
ScintillatorMC: CMakeFiles/ScintillatorMC.dir/src/SensitiveDetector.cc.o
ScintillatorMC: CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorScintillator.cc.o
ScintillatorMC: CMakeFiles/ScintillatorMC.dir/src/SensitiveDetectorTracker.cc.o
ScintillatorMC: CMakeFiles/ScintillatorMC.dir/src/SteppingAction.cc.o
ScintillatorMC: CMakeFiles/ScintillatorMC.dir/src/pCTconfig.cc.o
ScintillatorMC: CMakeFiles/ScintillatorMC.dir/build.make
ScintillatorMC: /home/ryan/G4/install/lib/libG4Tree.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4GMocren.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4visHepRep.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4RayTracer.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4VRML.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4OpenGL.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4gl2ps.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4interfaces.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4persistency.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4error_propagation.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4readout.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4physicslists.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4parmodels.so
ScintillatorMC: /home/ryan/root/install/lib/libCore.so
ScintillatorMC: /home/ryan/root/install/lib/libImt.so
ScintillatorMC: /home/ryan/root/install/lib/libRIO.so
ScintillatorMC: /home/ryan/root/install/lib/libNet.so
ScintillatorMC: /home/ryan/root/install/lib/libHist.so
ScintillatorMC: /home/ryan/root/install/lib/libGraf.so
ScintillatorMC: /home/ryan/root/install/lib/libGraf3d.so
ScintillatorMC: /home/ryan/root/install/lib/libGpad.so
ScintillatorMC: /home/ryan/root/install/lib/libROOTDataFrame.so
ScintillatorMC: /home/ryan/root/install/lib/libTree.so
ScintillatorMC: /home/ryan/root/install/lib/libTreePlayer.so
ScintillatorMC: /home/ryan/root/install/lib/libRint.so
ScintillatorMC: /home/ryan/root/install/lib/libPostscript.so
ScintillatorMC: /home/ryan/root/install/lib/libMatrix.so
ScintillatorMC: /home/ryan/root/install/lib/libPhysics.so
ScintillatorMC: /home/ryan/root/install/lib/libMathCore.so
ScintillatorMC: /home/ryan/root/install/lib/libThread.so
ScintillatorMC: /home/ryan/root/install/lib/libMultiProc.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4FR.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4vis_management.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4modeling.so
ScintillatorMC: /usr/lib/x86_64-linux-gnu/libGLU.so
ScintillatorMC: /usr/lib/x86_64-linux-gnu/libGL.so
ScintillatorMC: /usr/lib/x86_64-linux-gnu/libXmu.so
ScintillatorMC: /usr/lib/x86_64-linux-gnu/libXext.so
ScintillatorMC: /usr/lib/x86_64-linux-gnu/libXt.so
ScintillatorMC: /usr/lib/x86_64-linux-gnu/libSM.so
ScintillatorMC: /usr/lib/x86_64-linux-gnu/libICE.so
ScintillatorMC: /usr/lib/x86_64-linux-gnu/libX11.so
ScintillatorMC: /usr/lib/x86_64-linux-gnu/libQt5OpenGL.so.5.12.8
ScintillatorMC: /usr/lib/x86_64-linux-gnu/libQt5PrintSupport.so.5.12.8
ScintillatorMC: /usr/lib/x86_64-linux-gnu/libQt5Widgets.so.5.12.8
ScintillatorMC: /usr/lib/x86_64-linux-gnu/libQt5Gui.so.5.12.8
ScintillatorMC: /usr/lib/x86_64-linux-gnu/libQt5Core.so.5.12.8
ScintillatorMC: /home/ryan/G4/install/lib/libG4run.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4event.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4tracking.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4processes.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4analysis.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4zlib.so
ScintillatorMC: /usr/lib/x86_64-linux-gnu/libexpat.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4digits_hits.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4track.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4particles.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4geometry.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4materials.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4graphics_reps.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4intercoms.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4global.so
ScintillatorMC: /home/ryan/G4/install/lib/libG4clhep.so
ScintillatorMC: CMakeFiles/ScintillatorMC.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Linking CXX executable ScintillatorMC"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ScintillatorMC.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ScintillatorMC.dir/build: ScintillatorMC

.PHONY : CMakeFiles/ScintillatorMC.dir/build

CMakeFiles/ScintillatorMC.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ScintillatorMC.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ScintillatorMC.dir/clean

CMakeFiles/ScintillatorMC.dir/depend:
	cd /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build /home/ryan/Documents/ScintillatorMC/Full_scint/ScintillatorMC/build/CMakeFiles/ScintillatorMC.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ScintillatorMC.dir/depend
