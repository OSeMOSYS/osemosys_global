
@echo off
setlocal EnableDelayedExpansion

REM the for-loop is for (1 start position, 1 how many steps to jump each time, until which number)
for /L %%i in (1,1,3) do (
	
	REM lpName is the name of the output from glpsol and solName is the output from CPLEX
	SET lpName=matrix_%%i.lp
	SET solName=solution_%%i.sol
	REM this command starts glpsol, updpate the OSeMOSYS_model and data to your corresponding files.
	
		START /B CMD /C CALL glpsol -m OSeMOSYS_model.txt -d data_%%i.txt --wlp matrix_%%i.lp --check
	
	REM checks if the lp file is in the folder
	call :waitfile "!lpName!"
	REM lp file is found, but wait 180 sek to make sure the file is built before starting CPLEX
	ECHO Found file !lpName!
	TIMEOUT /T 180 >nul

	REM make a new empty file for mycplexcommands
	break>mycplexcommands
	REM echo writes each line to mycplexcommands that I want to execute in CPLEX
	(
	echo read matrix_%%i.lp
	echo optimize
	echo write
	echo solution_%%i.sol
	echo quit
	) > mycplexcommands

	REM executes the cplex script written above
	cplex < mycplexcommands

	REM checks if the sol file is in the folder otherwise it waits
	call :waitfilesol "!solName!"
	REM sol file is found, wait 180 sek to make sure the file is built before starting python
	ECHO Found file !solName!
	TIMEOUT /T 180 >nul

	REM the .sol file is input to transform python script
	python transform_31072013.py solution_%%i.sol solution_%%i.txt
	REM the text file is sorted to alphabetical sort
	sort/+1<solution_%%i.txt>solution_%%i_sorted.txt

	REM These commands delete lp and sol file which are large. Comment out is you don't want to delete them.
	del matrix_%%i.lp
	del solution_%%i.sol
	del solution_%%i.txt

)

ECHO All done!
pause
exit /b 0

REM this command is active as long as glpsol is building the lp file (will be printed every 120 sek)
:waitfile
    ECHO waitforlpfile: %~1
:waitfile2
	ECHO Waiting...
    TIMEOUT /T 120 >nul
    IF not EXIST "%~1" goto waitfile2
exit /b 0

REM this command is active as long as glpsol is building the sol file (will be printed every 120 sek)
:waitfilesol
    ECHO waitforsolfile: %~1
:waitfilesol2
	ECHO Waiting...
    TIMEOUT /T 120 >nul
    IF not EXIST "%~1" goto waitfilesol2
exit /b 0

