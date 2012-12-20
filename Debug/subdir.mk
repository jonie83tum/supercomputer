################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../check_allreduce.c \
../check_indexed.c \
../compute_solution.c \
../finalization.c \
../gccg.c \
../initialization.c \
../initialization_v2.c \
../sendrecv.c \
../test_functions.c \
../util_read_files.c \
../util_write_files.c 

OBJS += \
./check_allreduce.o \
./check_indexed.o \
./compute_solution.o \
./finalization.o \
./gccg.o \
./initialization.o \
./initialization_v2.o \
./sendrecv.o \
./test_functions.o \
./util_read_files.o \
./util_write_files.o 

C_DEPS += \
./check_allreduce.d \
./check_indexed.d \
./compute_solution.d \
./finalization.d \
./gccg.d \
./initialization.d \
./initialization_v2.d \
./sendrecv.d \
./test_functions.d \
./util_read_files.d \
./util_write_files.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -I/usr/local/include -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


