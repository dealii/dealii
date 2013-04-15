@echo off
REM
REM Copyright 2005-2013 Intel Corporation.  All Rights Reserved.
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

:: Getting parameters
if ("%1") == ("") goto error
if ("%2") == ("") goto error
if ("%3") == ("") goto error
set arch=%1
if ("%2") == ("debug") set postfix=_debug
set output_dir=%3

:: Optional 4th parameter to set install root
if ("%4") NEQ ("") set TBBROOT=%4
:: Actually we can set install root by ourselves
if ("%TBBROOT%") == ("") set TBBROOT=%~d0%~p0..\..\

:: Getting vs folders in case vc_mt binaries are not provided
if ("%VS90COMNTOOLS%")  NEQ ("") set vc_dir=vc9
if ("%VS100COMNTOOLS%") NEQ ("") set vc_dir=vc10
if ("%VS110COMNTOOLS%") NEQ ("") set vc_dir=vc11

:: Are we standalone/oss or inside compiler?
if exist "%TBBROOT%\bin\%arch%\vc9\tbb%postfix%.dll" set interim_path=bin\%arch%
if exist "%TBBROOT%\..\redist\%arch%\tbb\vc9\tbb%postfix%.dll" set interim_path=..\redist\%arch%\tbb
if ("%interim_path%") == ("") goto error

:: Do we provide vc_mt binaries?
if exist "%TBBROOT%\%interim_path%\vc_mt\tbb%postfix%.dll" set vc_dir=vc_mt
if ("%vc_dir%") == ("") goto error

:: We know everything we wanted and there are no errors
:: Copying binaries

copy "%TBBROOT%\%interim_path%\%vc_dir%\tbb%postfix%.dll" "%output_dir%"
copy "%TBBROOT%\%interim_path%\%vc_dir%\tbb%postfix%.pdb" "%output_dir%"
copy "%TBBROOT%\%interim_path%\%vc_dir%\tbbmalloc%postfix%.dll" "%output_dir%"
copy "%TBBROOT%\%interim_path%\%vc_dir%\tbbmalloc%postfix%.pdb" "%output_dir%"
if exist "%TBBROOT%\%interim_path%\%vc_dir%\tbb_preview%postfix%.dll" copy "%TBBROOT%\%interim_path%\%vc_dir%\tbb_preview%postfix%.dll" "%output_dir%"
if exist "%TBBROOT%\%interim_path%\%vc_dir%\tbb_preview%postfix%.pdb" copy "%TBBROOT%\%interim_path%\%vc_dir%\tbb_preview%postfix%.pdb" "%output_dir%"

goto end
:error
echo Error occurred in libraries copying during post-build step.
exit /B 1

:end
exit /B 0

