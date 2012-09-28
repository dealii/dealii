#
# Check for various system features:
#


CHECK_INCLUDE_FILE("sys/stat.h" HAVE_SYS_STAT_H)

CHECK_INCLUDE_FILE("sys/syscall.h" HAVE_SYS_SYSCALL_H)

CHECK_INCLUDE_FILE("sys/times.h" HAVE_SYS_TIMES_H)

CHECK_INCLUDE_FILE("sys/types.h" HAVE_SYS_TYPES_H)

CHECK_INCLUDE_FILE("unistd.h" HAVE_UNISTD_H)
