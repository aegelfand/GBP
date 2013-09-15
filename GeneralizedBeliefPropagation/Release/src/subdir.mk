################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/JoinGraph.cpp \
../src/LogFunction.cpp \
../src/Main.cpp \
../src/Partition.cpp \
../src/RegionGraph.cpp \
../src/mersenne.cpp 

OBJS += \
./src/JoinGraph.o \
./src/LogFunction.o \
./src/Main.o \
./src/Partition.o \
./src/RegionGraph.o \
./src/mersenne.o 

CPP_DEPS += \
./src/JoinGraph.d \
./src/LogFunction.d \
./src/Main.d \
./src/Partition.d \
./src/RegionGraph.d \
./src/mersenne.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


