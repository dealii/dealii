@echo off
REM
REM Copyright 2005-2010 Intel Corporation.  All Rights Reserved.
REM
REM This file is part of Threading Building Blocks.
REM
REM Threading Building Blocks is free software; you can redistribute it
REM and/or modify it under the terms of the GNU General Public License
REM version 2 as published by the Free Software Foundation.
REM
REM Threading Building Blocks is distributed in the hope that it will be
REM useful, but WITHOUT ANY WARRANTY; without even the implied warranty
REM of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
REM GNU General Public License for more details.
REM
REM You should have received a copy of the GNU General Public License
REM along with Threading Building Blocks; if not, write to the Free Software
REM Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
REM
REM As a special exception, you may use this file as part of a free software
REM library without restriction.  Specifically, if other files instantiate
REM templates or use macros or inline functions from this file, or you compile
REM this file and link it with other files to produce an executable, this
REM file does not by itself cause the resulting executable to be covered by
REM the GNU General Public License.  This exception does not however
REM invalidate any other reasons why the executable file might be covered by
REM the GNU General Public License.
REM
setlocal
for %%D in ("%tbb_root%") do set actual_root=%%~fD
set fslash_root=%actual_root:\=/%
set bin_dir=%CD%
set fslash_bin_dir=%bin_dir:\=/%
set _INCLUDE=INCLUDE& set _LIB=LIB
if not x%UNIXMODE%==x set _INCLUDE=CPATH& set _LIB=LIBRARY_PATH

if exist tbbvars.bat goto skipbat
echo Generating local tbbvars.bat
echo @echo off>tbbvars.bat
echo SET TBB30_INSTALL_DIR=%actual_root%>>tbbvars.bat
echo SET TBB_ARCH_PLATFORM=%arch%\%runtime%>>tbbvars.bat
echo SET TBB_TARGET_ARCH=%arch%>>tbbvars.bat
echo SET %_INCLUDE%=%%TBB30_INSTALL_DIR%%\include;%%%_INCLUDE%%%>>tbbvars.bat
echo SET %_LIB%=%bin_dir%;%%%_LIB%%%>>tbbvars.bat
echo SET PATH=%bin_dir%;%%PATH%%>>tbbvars.bat
if not x%UNIXMODE%==x echo SET LD_LIBRARY_PATH=%bin_dir%;%%LD_LIBRARY_PATH%%>>tbbvars.bat
:skipbat

if exist tbbvars.sh goto skipsh
echo Generating local tbbvars.sh
echo #!/bin/sh>tbbvars.sh
echo export TBB30_INSTALL_DIR="%fslash_root%">>tbbvars.sh
echo export TBB_ARCH_PLATFORM="%arch%\%runtime%">>tbbvars.sh
echo export TBB_TARGET_ARCH="%arch%">>tbbvars.sh
echo export %_INCLUDE%="${TBB30_INSTALL_DIR}/include;$%_INCLUDE%">>tbbvars.sh
echo export %_LIB%="%fslash_bin_dir%;$%_LIB%">>tbbvars.sh
echo export PATH="%fslash_bin_dir%;$PATH">>tbbvars.sh
if not x%UNIXMODE%==x echo export LD_LIBRARY_PATH="%fslash_bin_dir%;$LD_LIBRARY_PATH">>tbbvars.sh
:skipsh

if exist tbbvars.csh goto skipcsh
echo Generating local tbbvars.csh
echo #!/bin/csh>tbbvars.csh
echo setenv TBB30_INSTALL_DIR "%actual_root%">>tbbvars.csh
echo setenv TBB_ARCH_PLATFORM "%arch%\%runtime%">>tbbvars.csh
echo setenv TBB_TARGET_ARCH "%arch%">>tbbvars.csh
echo setenv %_INCLUDE% "${TBB30_INSTALL_DIR}\include;$%_INCLUDE%">>tbbvars.csh
echo setenv %_LIB% "%bin_dir%;$%_LIB%">>tbbvars.csh
echo setenv PATH "%bin_dir%;$PATH">>tbbvars.csh
if not x%UNIXMODE%==x echo setenv LD_LIBRARY_PATH "%bin_dir%;$LD_LIBRARY_PATH">>tbbvars.csh
:skipcsh

endlocal
exit
