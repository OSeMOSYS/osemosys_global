# Taco Niet
# Date:  September 26, 2016 - 2018
# Script to run OSeMOSYS with various parameters using full nodes on a SLURM scheduler managed HPC system.

#!/usr/bin/python
import subprocess # for running glpsol
from shutil import copy2 # Allow us to make the needed copies of the parameter file
from math import sqrt
import zipfile, os, time, io

if __name__ == '__main__':
    modelruns = []
    tailcommands = []

	# Setup task list here.  This can be created by a series of nested for loops.
    param3List = ['10.6', '15.9', '21.2']
    param2List = ['0000', '0050', '0100']
    param1List = ['35', '65', '100']    

    NumRunsTotal = len(param1List)*len(param2List)*len(param3List)
    
    NumRunstoGroup = 600  # Number of runs to group for each qsub command
    NumberofCores = 32
    WallTime = '06:00:00' # Walltime in hh:mm:ss format
    ZipWallTime = '08:00:00' # Walltime in hh:mm:ss format for zipping files
    
    for param1 in param1List:
        for param2 in param2List:
            for param3 in param3List:
                print('Setting up '+param1+' '+ param2+' '+param3+'.')
                copy2('params_base.dat', param3 +'_'+ param1 +'_'+ param2 + 'params.dat')
                f = open(param3 +'_'+ param1 +'_'+ param2 + 'params.dat', 'a') # Open file for appending...
                f.write('\n')

                # Write Parameter 1 into data file...
                f.write('param	Parameter1 default 0	[REG,*,2020]	:=\n')
                f.write('tparam3  '+str(float(param1)/1000)+'\n')
                f.write('tEXAMPLE	1.000   ;\n\n')

                # Set the param2 size:
                f.write('param TechnologyMaxCapacity default 99999999999 := REG DAM 2020 '+param2+';\n\n')
                
                # Set the param3 residual capacity:
                f.write('param	ResidualCapacity	default	0	:=[REG,*,*]	:   2020	:=\n')
                f.write('tEXAMPLE	0.24    tparam3	'+param3+';\n\n')               
               
                # Specify output files
                # NOTE:  This needs to be adjusted based on which version of OSeMOSYS you are using...  See  OutputFileAsParameter branch for one example.
                f.write('param SelectedResults := \''+ param3 +'_'+ param1 +'_'+ param2 + 'SelectedResults.csv\';\n')
                                    
                # End the file and close it.
                f.write('end;\n')
                f.close()
                
                # append file name to list of files to process

                modelruns.append('glpsol -m OSeMOSYS.mod -d ' + param3 +'_'+ param1 +'_'+ param2 + '_'+ 'params.dat --log '  + param3 +'_'+ param1 +'_'+ param2 + 'GLPScreen.txt -o ' + param3 +'_'+ param1 +'_'+ param2 + 'GLPOutput.txt')

    # Setup qsub batch job scripts
    # First open the files, overwritting previous versions, and write the header material
    
    NumJobs = int(len(modelruns)/(NumRunstoGroup*NumberofCores))+1
        
    print('NumJobs = '+str(NumJobs))
    
    # Setup job files to request full nodes:
    for count in range(1, NumJobs+1):    
        with open("RunGroup"+str(int(count)) + ".sh", "w") as fpbs:
            fpbs.write("#!/bin/bash\n")
            fpbs.write("#SBATCH --time="+WallTime+"\n")
            
            fpbs.write("#SBATCH --job-name=RunGroup"+str(int(count))+"\n")
            fpbs.write("#SBATCH --output=%x-%j.out\n")
            fpbs.write("#SBATCH --mail-user=tniet@bcit.ca\n")
            fpbs.write("#SBATCH --mail-type=ALL\n")
            
            fpbs.write("#SBATCH --nodes=1\n")
            fpbs.write("#SBATCH --ntasks-per-node="+str(NumberofCores)+"\n")

            fpbs.write("# Script to run job group "+str(int(count))+".\n")

            fpbs.write("\n")
            fpbs.write("module load python/3.6.3\n")

            fpbs.write("\n")
            
            fpbs.write("echo \"Current working directory is `pwd`\"\n")
            fpbs.write("echo \"Running on hostname `hostname`\"\n")
        # Create the command file for the run...
        with open("RunCommands"+str(int(count)) + ".sh", "w") as fpbs:
            for countjobs in range(((count-1)*NumRunstoGroup*NumberofCores + 1), min(((count)*NumRunstoGroup*NumberofCores)+1, len(modelruns))):
                fpbs.write ("echo \"Starting run at: `date`\"")
                fpbs.write (" && $HOME/project/glpk-4.61/examples/"+modelruns[countjobs])
                fpbs.write (" && " + tailcommands[countjobs]+'\n')
    
    # And finally print the footer material into the file
    for count in range(1, NumJobs+1):    
        with open("RunGroup"+str(int(count)) + ".sh", "a") as fpbs:
            fpbs.write("parallel < RunCommands"+str(int(count)) + ".sh\n")
            fpbs.write ("echo \"Runs for group "+str(int(count))+"finished with exit code $? at: `date`\"\n")
            fpbs.write("#End of file\n")
    
    # Submit jobs to job queue with sbatch
    for count in range(1, NumJobs+1):
        filetosubmit = "RunGroup"+str(int(count)) + ".sh"
        cmd = ["sbatch"]
        cmd.append(filetosubmit)
        print ("Submitting Job "+ str(count) +" with command: "+str(cmd))
        status = subprocess.check_output(cmd)
        print (status.decode('UTF-8'))
        #print(status.decode('UTF-8').split(" ")[3])
        jobnum = int(status.decode('UTF-8').split(" ")[3])

    # And we're done...
    print("RunCedar completed.")